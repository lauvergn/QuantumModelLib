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

!> @brief Module which makes the initialization, calculation of the HNNHp potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_HNNHp

  USE mod_NumParameters

  IMPLICIT NONE

!> @brief Derived type in which the HNNHp parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE Param_HNNHp

     PRIVATE

     real(kind=Rkind), allocatable, public :: Qref(:)

     integer                       :: ndim,nb_func
     real(kind=Rkind), allocatable :: F(:)
     integer, allocatable          :: tab_func(:,:)

     integer :: option    = -1
     logical :: PubliUnit = .FALSE.
 
  END TYPE Param_HNNHp
 
  PRIVATE eval_HNNHpPot1,dealloc_HNNHp
 
  CONTAINS
!> @brief Subroutine which makes the initialization of the HNNHp parameters.
!!
!! @param Para_HNNHp          TYPE(Param_HNNHp):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.

  SUBROUTINE dealloc_HNNHp(Para_HNNHp)
    USE mod_Lib
    IMPLICIT NONE

    TYPE (Param_HNNHp),          intent(inout) :: Para_HNNHp

    Para_HNNHp%option    = -1
    Para_HNNHp%PubliUnit = .FALSE.

    IF (allocated(Para_HNNHp%Qref))     deallocate(Para_HNNHp%Qref)
    IF (allocated(Para_HNNHp%F))        deallocate(Para_HNNHp%F)
    IF (allocated(Para_HNNHp%tab_func)) deallocate(Para_HNNHp%tab_func)

  END SUBROUTINE dealloc_HNNHp
  SUBROUTINE Init_HNNHp(Para_HNNHp,ndim,option,nio,read_param,PubliUnit)
    USE mod_Lib
    IMPLICIT NONE

    TYPE (Param_HNNHp),          intent(inout) :: Para_HNNHp
    integer,                     intent(in)    :: option,ndim
    integer,           optional, intent(in)    :: nio
    logical,           optional, intent(in)    :: read_param,PubliUnit

    logical :: read_param_loc
    integer :: nio_fit,nb_func,nb_columns,i,j,k
    character (len=50) :: name_Qref


    character (len=:), allocatable  :: FileName


    CALL dealloc_HNNHp(Para_HNNHp)

    IF (present(PubliUnit)) Para_HNNHp%PubliUnit = PubliUnit

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_HNNHpPot '
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present '
       write(out_unitp,*) ' => impossible to read the input file '
       STOP 'STOP in Init_HNNHpPot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    Para_HNNHp%option = option

    IF (ndim /= 6) THEN
       write(out_unitp,*) ' ERROR in Init_HNNHpPot '
       write(out_unitp,*) ' ndim /= 6. it is not possible for this potential'
       STOP 'STOP in Init_HNNHpPot: ndim MUST equal to 6'
    END IF
    Para_HNNHp%ndim = 6


    IF (Para_HNNHp%option < 1 .OR. Para_HNNHp%option > 1) Para_HNNHp%option = 1

    IF (read_param_loc) THEN
      STOP 'Init_HNNHp: nothing to read'

    ELSE

      SELECT CASE (Para_HNNHp%option)
      CASE (1)

        FileName = make_FileName('InternalData/HNNHp/n2h2pcc5.txt')
        CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
        read(nio_fit,*) Para_HNNHp%nb_func
        !write(6,*) 'nb_func',Para_HNNHp%nb_func

        allocate(Para_HNNHp%Qref(Para_HNNHp%ndim))
        allocate(Para_HNNHp%F(Para_HNNHp%nb_func))
        allocate(Para_HNNHp%tab_func(Para_HNNHp%ndim,Para_HNNHp%nb_func))

         k = 0
         DO
          nb_columns = min(6,Para_HNNHp%nb_func-k)
          !write(6,*) k+1,k+nb_columns,nb_columns
          IF (nb_columns == 0) EXIT
           read(nio_fit,11) (Para_HNNHp%tab_func(1:Para_HNNHp%ndim,j),Para_HNNHp%F(j),j=k+1,k+nb_columns)
 11        format(6i1,f15.8,5(2x,6i1,f15.8))
           k = k + nb_columns
         END DO
         read(nio_fit,*) Para_HNNHp%Qref(:)
         !write(6,*) 'Qref',Para_HNNHp%Qref

         close(nio_fit)
         deallocate(FileName)


      CASE Default

          write(out_unitp,*) ' ERROR in Init_HNNHp '
          write(out_unitp,*) ' This option is not possible. option: ',Para_HNNHp%option
          write(out_unitp,*) ' Its value MUST be 1'

          STOP

      END SELECT
    END IF

    IF (Para_HNNHp%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    END IF

  END SUBROUTINE Init_HNNHp
!> @brief Subroutine wich prints the HNNHp parameters.
!!
!! @param Para_HNNHp         TYPE(Param_HNNHp):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_HNNHp(Para_HNNHp,nio)
    
    TYPE(Param_HNNHp), intent(in) :: Para_HNNHp
    integer         , intent(in) :: nio

    write(nio,*) 'trans-HNNH+ current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '   Q(2) \ Q(3)                         '
    write(nio,*) '         N------------N                '
    write(nio,*) '            Q(1)   Q(5) \  Q(4)        '
    write(nio,*) '                          \            '
    write(nio,*) '                            H    +Q(6) '
    write(nio,*) '                                       '
    write(nio,*) '  Q(1) = dNN (Bohr)                    '
    write(nio,*) '  Q(2) = dNH (Bohr)                    '
    write(nio,*) '  Q(3) = aHNN (Radian)                 '
    write(nio,*) '  Q(4) = dNH (Bohr)                    '
    write(nio,*) '  Q(5) = aHNN (Radian)                 '
    write(nio,*) '  Q(6) = Torsion (HNNH) (Radian)       '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) 'Ref:'
    write(nio,*) '  D. Lauvergnat and M. Hochlaf, J. Chem. Phys. 130, 224312 (2009).'
    write(nio,*) '      https://doi.org/10.1063/1.3154141'
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',Para_HNNHp%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',Para_HNNHp%option
    write(nio,*)
    write(nio,*) '---------------------------------------'

    SELECT CASE (Para_HNNHp%option)

    CASE (1)
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Level: RCCSD(T)/aug-cc-pV5Z*         '
    write(nio,*) '    spdfg for N and spdf for H         '
    write(nio,*) '                                       '
    write(nio,*) '  Reference Trans-geometry:            '
    write(nio,*) '  Q(1) = rNN  = 2.2           (Bohr)   '
    write(nio,*) '  Q(2) = rH1  = 2.0           (Bohr)   '
    write(nio,*) '  Q(3) = aH1  = 2.1816615395  (rad)    '
    write(nio,*) '  Q(4) = rH2  = 2.0           (Bohr)   '
    write(nio,*) '  Q(5) = aH3  = 2.1816615395  (rad)    '
    write(nio,*) '  Q(6) = phi  = 3.1415926169  (rad)    '
    write(nio,*) '                                       '
    write(nio,*) '  V = -110.16613531 Hartree            '
    write(nio,*) ' grad(:) =[-0.0002097800, 0.0133780500,'
    write(nio,*) '            0.0054179000, 0.0133780500,'
    write(nio,*) '            0.0054179000, 0.0000000000]'
    write(nio,*) '                                       '

    CONTINUE

    CASE Default
        write(out_unitp,*) ' ERROR in write_HNNHpPot '
        write(out_unitp,*) ' This option is not possible. option: ',Para_HNNHp%option
        write(out_unitp,*) ' Its value MUST be 1'

        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end HNNHp current parameters'

  END SUBROUTINE Write_HNNHp
!> @brief Subroutine wich prints the HNNHp parameters.
!!
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write0_HNNHp(nio)

    integer         , intent(in) :: nio

    write(nio,*) 'HNNHp default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end HNNHp default parameters'

  END SUBROUTINE Write0_HNNHp
  SUBROUTINE get_Q0_HNNHp(Q0,Para_HNNHp,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (Param_HNNHp),          intent(in)    :: Para_HNNHp
    integer,                     intent(in)    :: option

    IF (size(Q0) /= 6) THEN
      write(out_unitp,*) ' ERROR in get_Q0_HNNHp '
      write(out_unitp,*) ' The size of Q0 is not ndim=6: '
      write(out_unitp,*) ' size(Q0)',size(Q0)
      STOP
    END IF

    SELECT CASE (option)
    CASE (0) ! ref
      Q0(:) = Para_HNNHp%Qref
    CASE Default ! ref
      Q0(:) = Para_HNNHp%Qref
    END SELECT

  END SUBROUTINE get_Q0_HNNHp


!> @brief Subroutine wich calculates the HNNHp potential (for the 3 models) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_HNNHp         TYPE(Param_HNNHp):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_HNNHpPot(Mat_OF_PotDia,dnQ,Para_HNNHp,nderiv)
    USE mod_dnSca

    TYPE(Param_HNNHp) , intent(in)   :: Para_HNNHp
    TYPE(dnSca),       intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE(dnSca),       intent(in)    :: dnQ(:) !
    integer          , intent(in)    :: nderiv

    SELECT CASE (Para_HNNHp%option)

    CASE (1)
      CALL eval_HNNHpPot1(Mat_OF_PotDia,dnQ,Para_HNNHp)

    CASE Default
        write(out_unitp,*) ' ERROR in eval_HNNHpPot '
        write(out_unitp,*) ' This option is not possible. option: ',Para_HNNHp%option
        write(out_unitp,*) ' Its value MUST be 1'

        STOP
    END SELECT

  END SUBROUTINE eval_HNNHpPot

!> @brief Subroutine wich calculates the HNNHp potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_HNNHp         TYPE(Param_HNNHp):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).

  SUBROUTINE eval_HNNHpPot1(Mat_OF_PotDia,dnQ,Para_HNNHp)
    USE mod_dnSca

    TYPE(dnSca),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE(dnSca),        intent(in)    :: dnQ(:)
    TYPE(Param_HNNHp) , intent(in)    :: Para_HNNHp


    TYPE(dnSca)        :: DQ(6,6)
    TYPE(dnSca)        :: Vtemp
    integer            :: i,j

    !write(6,*) ' sub eval_HNNHpPot1' ; flush(6)

      ! Warning, the coordinate ordering in the potential data (from the file) is different from the z-matrix one.
      DQ(:,1) = dnQ([2,1,4,3,5,6]) - Para_HNNHp%Qref(:)
      DO j=2,size(DQ,dim=2)
        DQ(1:5,j) = DQ(1:5,j-1) * DQ(1:5,1)
      END DO
      DO i=size(DQ,dim=2),1,-1
        DQ(6,i) = cos( DQ(6,1) * real(i,kind=Rkind) )
      END DO


      Mat_OF_PotDia(1,1) = ZERO
      Vtemp = Mat_OF_PotDia(1,1) ! to have a correct initialization

      DO j=1,Para_HNNHp%nb_func
        Vtemp = Para_HNNHp%F(j)
        DO i=1,Para_HNNHp%ndim
          IF (Para_HNNHp%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,Para_HNNHp%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO


!-----------------------------------------------------------------------!

   CALL dealloc_dnSca(Vtemp)
   CALL dealloc_dnSca(DQ)

   !write(6,*) ' end eval_HNNHpPot1' ; flush(6)

  END SUBROUTINE eval_HNNHpPot1

  ! here we suppose that the atom ordering: N1-N2-H1-H2
  ! the bounds are N1-N2, N1-H1, n2-H2
  SUBROUTINE Cart_TO_Q_HNNHp(dnX,dnQ,Para_HNNHp,nderiv)
    USE mod_dnSca

    TYPE(dnSca),       intent(in)    :: dnX(:,:)
    TYPE(dnSca),       intent(inout) :: dnQ(:)
    TYPE(Param_HNNHp) , intent(in)   :: Para_HNNHp
    integer          , intent(in)    :: nderiv


    ! local vector
    TYPE(dnSca)    :: VecNN(3),VecN1H1(3),VecN2H2(3)


    VecNN(:)   = dnX(:,1)-dnX(:,2)
    VecN1H1(:) = dnX(:,1)-dnX(:,3)
    VecN2H2(:) = dnX(:,2)-dnX(:,4)

    dnQ(1) = sqrt(dot_product(VecNN,VecNN))
    dnQ(2) = sqrt(dot_product(VecN1H1,VecN1H1))
    dnQ(4) = sqrt(dot_product(VecN2H2,VecN2H2))

    dnQ(3) = acos(dot_product(VecNN,VecN1H1)/(dnQ(1)*dnQ(2)))
    dnQ(5) = acos(dot_product(VecNN,VecN2H2)/(dnQ(1)*dnQ(4)))



  END SUBROUTINE Cart_TO_Q_HNNHp

END MODULE mod_HNNHp
