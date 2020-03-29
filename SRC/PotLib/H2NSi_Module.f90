!===========================================================================
!===========================================================================
!This file is part of ModelLib.
!
!    ModelLib is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ModelLib is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ModelLib.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2016  David LAUVERGNAT, Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================

!> @brief Module which makes the initialization, calculation of the H2NSi potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_H2NSi

  USE mod_NumParameters

  IMPLICIT NONE

!> @brief Derived type in which the H2NSi parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE H2NSiPot_t

     PRIVATE

     real(kind=Rkind), allocatable, public :: Qref(:)

     integer                       :: ndim,nb_func
     real(kind=Rkind), allocatable :: F(:)
     integer, allocatable          :: tab_func(:,:)

     integer :: option    = -1
     logical :: PubliUnit = .FALSE.
 
  END TYPE H2NSiPot_t
 
  PRIVATE eval_H2NSiPot1,dealloc_H2NSi
 
  CONTAINS
!> @brief Subroutine which makes the initialization of the H2NSi parameters.
!!
!! @param H2NSiPot          TYPE(H2NSiPot_t):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.

  SUBROUTINE dealloc_H2NSi(H2NSiPot)
    USE mod_Lib
    IMPLICIT NONE

    TYPE (H2NSiPot_t),          intent(inout) :: H2NSiPot

    H2NSiPot%option    = -1
    H2NSiPot%PubliUnit = .FALSE.

    IF (allocated(H2NSiPot%Qref))     deallocate(H2NSiPot%Qref)
    IF (allocated(H2NSiPot%F))        deallocate(H2NSiPot%F)
    IF (allocated(H2NSiPot%tab_func)) deallocate(H2NSiPot%tab_func)

  END SUBROUTINE dealloc_H2NSi
  SUBROUTINE Init_H2NSi(H2NSiPot,ndim,option,nio,read_param,PubliUnit)
    USE mod_Lib
    IMPLICIT NONE

    TYPE (H2NSiPot_t),          intent(inout) :: H2NSiPot
    integer,                     intent(in)    :: option,ndim
    integer,           optional, intent(in)    :: nio
    logical,           optional, intent(in)    :: read_param,PubliUnit

    logical :: read_param_loc
    integer :: nio_fit,nb_func,nb_columns,i,j,k
    character (len=50) :: name_Qref


    character (len=:), allocatable  :: FileName


    CALL dealloc_H2NSi(H2NSiPot)

    IF (present(PubliUnit)) H2NSiPot%PubliUnit = PubliUnit

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_H2NSiPot '
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present '
       write(out_unitp,*) ' => impossible to read the input file '
       STOP 'STOP in Init_H2NSiPot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    H2NSiPot%option = option

    IF (ndim /= 6) THEN
       write(out_unitp,*) ' ERROR in Init_H2NSiPot '
       write(out_unitp,*) ' ndim /= 6. it is not possible for this potential'
       STOP 'STOP in Init_H2NSiPot: ndim MUST equal to 6'
    END IF
    H2NSiPot%ndim = 6


    IF (H2NSiPot%option < 1 .OR. H2NSiPot%option > 2) H2NSiPot%option = 1

    IF (read_param_loc) THEN
      STOP 'Init_H2NSi: nothing to read'

    ELSE

      SELECT CASE (H2NSiPot%option)
      CASE (1)

        FileName = make_FileName('InternalData/H2NSi/h2nsif12a.txt')
        CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
        allocate(H2NSiPot%Qref(H2NSiPot%ndim))

        DO i=1,H2NSiPot%ndim
          read(nio_fit,*) name_Qref,H2NSiPot%Qref(i)
        END DO
        H2NSiPot%Qref(4:6) = H2NSiPot%Qref(4:6) / (180._Rkind/pi)
        H2NSiPot%Qref(6)   = pi

        read(nio_fit,*) H2NSiPot%nb_func
        allocate(H2NSiPot%F(H2NSiPot%nb_func))
        allocate(H2NSiPot%tab_func(H2NSiPot%ndim,H2NSiPot%nb_func))


        DO i=1,H2NSiPot%nb_func
           read(nio_fit,*) name_Qref,H2NSiPot%tab_func(:,i),H2NSiPot%F(i)
           H2NSiPot%F(i) = H2NSiPot%F(i) * (180._Rkind/pi)**sum(H2NSiPot%tab_func(4:H2NSiPot%ndim,i))
        END DO
        close(nio_fit)
        deallocate(FileName)

      CASE (2)

        FileName = make_FileName('InternalData/H2NSi/h2nsicc.pot')
        CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
        read(nio_fit,*) H2NSiPot%nb_func
        allocate(H2NSiPot%Qref(H2NSiPot%ndim))
        allocate(H2NSiPot%F(H2NSiPot%nb_func))
        allocate(H2NSiPot%tab_func(H2NSiPot%ndim,H2NSiPot%nb_func))

         k = 0
         DO
          nb_columns = min(6,H2NSiPot%nb_func-k)
          !write(6,*) k+1,k+nb_columns,nb_columns
          IF (nb_columns == 0) EXIT
           read(nio_fit,11) (H2NSiPot%tab_func(1:H2NSiPot%ndim,j),H2NSiPot%F(j),j=k+1,k+nb_columns)
 11        format(6i1,f15.8,5(2x,6i1,f15.8))
           k = k + nb_columns
         END DO
         read(nio_fit,*) H2NSiPot%Qref(:)
         H2NSiPot%Qref(6) = pi ! Qref(6) must be changed to be compatible with a z-matrix.
         close(nio_fit)
         deallocate(FileName)


      CASE Default

          write(out_unitp,*) ' ERROR in Init_H2NSi '
          write(out_unitp,*) ' This option is not possible. option: ',H2NSiPot%option
          write(out_unitp,*) ' Its value MUST be 1 or 2'

          STOP

      END SELECT
    END IF

    IF (H2NSiPot%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    END IF

  END SUBROUTINE Init_H2NSi
!> @brief Subroutine wich prints the H2NSi parameters.
!!
!! @param H2NSiPot         TYPE(H2NSiPot_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_H2NSi(H2NSiPot,nio)
    
    TYPE(H2NSiPot_t), intent(in) :: H2NSiPot
    integer         , intent(in) :: nio

    write(nio,*) 'H2NSi current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '   Q(2) \ Q(3)                         '
    write(nio,*) '         N-------------Si              '
    write(nio,*) '   Q(4) / Q(5)    Q(1)                 '
    write(nio,*) '       /                               '
    write(nio,*) '      H   +Q(6)                        '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Minimum:                             '
    write(nio,*) '  Q(1) = rSiN (Bohr)                   '
    write(nio,*) '  Q(2) = rH1  (Bohr)                   '
    write(nio,*) '  Q(3) = aH1  (Radian)                 '
    write(nio,*) '  Q(4) = rH2  (Bohr)                   '
    write(nio,*) '  Q(5) = aH3  (Radian)                 '
    write(nio,*) '  Q(6) = phi  (Radian)                 '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) 'Ref:'
    write(nio,*) '  D. Lauvergnat, M.L. Senent, L. Jutier, and M. Hochlaf, JCP 135, 074301 (2011).'
    write(nio,*) '      https://doi.org/10.1063/1.3624563'
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',H2NSiPot%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',H2NSiPot%option
    write(nio,*)
    write(nio,*) '---------------------------------------'

    SELECT CASE (H2NSiPot%option)

    CASE (1)
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Level: RCCSD(T)-F12a/cc-pVTZ-F12     '
    write(nio,*) '                                       '
    write(nio,*) '  Minimum:                             '
    write(nio,*) '  Q(1) = rSiN = 3.2262192852 (Bohr)    '
    write(nio,*) '  Q(2) = rH1  = 1.9084802130 (Bohr)    '
    write(nio,*) '  Q(3) = aH1  = 124.2468848  (°)       '
    write(nio,*) '  Q(4) = rH2  = 1.9084802130 (Bohr)    '
    write(nio,*) '  Q(5) = aH3  = 124.2468848  (°)       '
    write(nio,*) '  Q(6) = phi  = 180.         (°)       '
    write(nio,*) '                                       '
    write(nio,*) '  V = -344.92729937 Hartree            '
    write(nio,*) ' grad(:) =[ 0.0000000000, 0.0000000000,'
    write(nio,*) '            0.0000000000, 0.0000000000,'
    write(nio,*) '            0.0000000000, 0.0000000000]'
    write(nio,*) '                                       '
    write(nio,*) '  The angles are convert in radian.    '

    CASE (2)
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Level: RCCSD(T)/Aug-cc-pV5Z          '
    write(nio,*) '                                       '
    write(nio,*) '  Reference geometry:                  '
    write(nio,*) '  Q(1) = rSiN = 3.20         (Bohr)    '
    write(nio,*) '  Q(2) = rH1  = 1.90         (Bohr)    '
    write(nio,*) '  Q(3) = aH1  = 2.1816615395 (Radian)  '
    write(nio,*) '  Q(4) = rH2  = 1.90         (Bohr)    '
    write(nio,*) '  Q(5) = aH3  = 2.1816615395 (Radian)  '
    write(nio,*) '  Q(6) = phi  = pi           (Radian)  '
    write(nio,*) '                                       '
    write(nio,*) '  V = -344.92196350 Hartree            '
    write(nio,*) ' grad(:) =[-0.0077950700,-0.0035043800,'
    write(nio,*) '            0.0024132300,-0.0035043800,'
    write(nio,*) '            0.0024132300,0.0000000000] '
    write(nio,*) '                                       '
    write(nio,*) ' Unpublished potential                 '

    CONTINUE

    CASE Default
        write(out_unitp,*) ' ERROR in write_H2NSiPot '
        write(out_unitp,*) ' This option is not possible. option: ',H2NSiPot%option
        write(out_unitp,*) ' Its value MUST be 1 or 2'

        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end H2NSi current parameters'

  END SUBROUTINE Write_H2NSi
!> @brief Subroutine wich prints the H2NSi parameters.
!!
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write0_H2NSi(nio)

    integer         , intent(in) :: nio

    write(nio,*) 'H2NSi default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end H2NSi default parameters'

  END SUBROUTINE Write0_H2NSi
  SUBROUTINE get_Q0_H2NSi(Q0,H2NSiPot,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (H2NSiPot_t),          intent(in)    :: H2NSiPot
    integer,                     intent(in)    :: option

    IF (size(Q0) /= 6) THEN
      write(out_unitp,*) ' ERROR in get_Q0_H2NSi '
      write(out_unitp,*) ' The size of Q0 is not ndim=6: '
      write(out_unitp,*) ' size(Q0)',size(Q0)
      STOP
    END IF

    SELECT CASE (option)
    CASE (0) ! ref
      Q0(:) = H2NSiPot%Qref([3,1,4,2,5,6])
    CASE Default ! ref
      Q0(:) = H2NSiPot%Qref([3,1,4,2,5,6])
    END SELECT

  END SUBROUTINE get_Q0_H2NSi
!> @brief Subroutine wich calculates the H2NSi potential (for the 3 models) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param H2NSiPot         TYPE(H2NSiPot_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_H2NSiPot(Mat_OF_PotDia,dnQ,H2NSiPot,nderiv)
    USE mod_dnS

    TYPE(H2NSiPot_t) , intent(in)   :: H2NSiPot
    TYPE (dnS_t),       intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),       intent(in)    :: dnQ(:) !
    integer          , intent(in)    :: nderiv

    SELECT CASE (H2NSiPot%option)

    CASE (1,2)
      CALL eval_H2NSiPot1(Mat_OF_PotDia,dnQ,H2NSiPot)

    CASE Default
        write(out_unitp,*) ' ERROR in eval_H2NSiPot '
        write(out_unitp,*) ' This option is not possible. option: ',H2NSiPot%option
        write(out_unitp,*) ' Its value MUST be 1 or 2'

        STOP
    END SELECT

  END SUBROUTINE eval_H2NSiPot

!> @brief Subroutine wich calculates the H2NSi potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param H2NSiPot         TYPE(H2NSiPot_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).

  SUBROUTINE eval_H2NSiPot1(Mat_OF_PotDia,dnQ,H2NSiPot)
    !Unpublished model potential
    USE mod_dnS

    TYPE (dnS_t),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)    :: dnQ(:)
    TYPE(H2NSiPot_t) , intent(in)    :: H2NSiPot


    TYPE (dnS_t)        :: DQ(6,6)
    TYPE (dnS_t)        :: Vtemp
    integer            :: i,j

    !write(6,*) ' sub eval_H2NSiPot1' ; flush(6)

      ! Warning, the coordinate ordering in the potential data (from the file) is different from the z-matrix one.
      DQ(:,1) = dnQ([2,4,1,3,5,6]) - H2NSiPot%Qref(:)
      DO j=2,size(DQ,dim=2)
        DQ(:,j) = DQ(:,j-1) * DQ(:,1)
      END DO


      Mat_OF_PotDia(1,1) = ZERO
      Vtemp = Mat_OF_PotDia(1,1) ! to have a correct initialization

      DO j=1,H2NSiPot%nb_func
        Vtemp = H2NSiPot%F(j)
        DO i=1,H2NSiPot%ndim
          !Vtemp = Vtemp * DQ(i,1)**H2NSiPot%tab_func(i,j)
          IF (H2NSiPot%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,H2NSiPot%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO


!-----------------------------------------------------------------------!

   CALL QML_dealloc_dnS(Vtemp)
   CALL QML_dealloc_dnS(DQ)

   !write(6,*) ' end eval_H2NSiPot1' ; flush(6)

  END SUBROUTINE eval_H2NSiPot1

END MODULE mod_H2NSi
