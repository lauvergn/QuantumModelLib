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

!> @brief Module which makes the initialization, calculation of the H2SiN potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_H2SiN

  USE mod_NumParameters

  IMPLICIT NONE

!> @brief Derived type in which the H2SiN parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE Param_H2SiN

     PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

     integer                       :: ndim,nb_func
     real(kind=Rkind), allocatable :: F(:)
     integer, allocatable          :: tab_func(:,:)

     integer :: option    = 1
     logical :: PubliUnit = .FALSE.
 
  END TYPE Param_H2SiN
 
  PRIVATE eval_H2SiNPot1
 
  CONTAINS
!> @brief Subroutine which makes the initialization of the H2SiN parameters.
!!
!! @param Para_H2SiN          TYPE(Param_H2SiN):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.

  SUBROUTINE Init_H2SiN(Para_H2SiN,ndim,option,nio,read_param,PubliUnit)
    USE mod_Lib
    IMPLICIT NONE

    TYPE (Param_H2SiN),          intent(inout) :: Para_H2SiN
    integer,                     intent(in)    :: option,ndim
    integer,           optional, intent(in)    :: nio
    logical,           optional, intent(in)    :: read_param,PubliUnit

    logical :: read_param_loc
    integer :: nio_fit,nb_func,nb_columns,j,k


    character (len=:), allocatable  :: FileName

    IF (present(PubliUnit)) Para_H2SiN%PubliUnit = PubliUnit

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_H2SiNPot '
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present '
       write(out_unitp,*) ' => impossible to read the input file '
       STOP 'STOP in Init_H2SiNPot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    Para_H2SiN%option = option

    IF (ndim /= 6) THEN
       write(out_unitp,*) ' ERROR in Init_H2SiNPot '
       write(out_unitp,*) ' ndim /= 6. it is not possible for this potential'
       STOP 'STOP in Init_H2SiNPot: ndim MUST equal to 6'
    END IF
    Para_H2SiN%ndim = 6


    IF (Para_H2SiN%option < 1 .OR. Para_H2SiN%option > 2) Para_H2SiN%option = 1

    IF (read_param_loc) THEN
      !CALL Read_H2SiNPot(Para_H2SiN,nio)
      STOP 'Init_H2SiN: nothing to read'

    ELSE

      SELECT CASE (Para_H2SiN%option)
      CASE (1)

        FileName = make_FileName('InternalData/H2SiN/h2sinf12a.pot')
        CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
        read(nio_fit,*) Para_H2SiN%nb_func
        allocate(Para_H2SiN%Qref(Para_H2SiN%ndim))
        allocate(Para_H2SiN%F(Para_H2SiN%nb_func))
        allocate(Para_H2SiN%tab_func(Para_H2SiN%ndim,Para_H2SiN%nb_func))

         k = 0
         DO
          nb_columns = min(6,Para_H2SiN%nb_func-k)
          !write(6,*) k+1,k+nb_columns,nb_columns
          IF (nb_columns == 0) EXIT
           read(nio_fit,11) (Para_H2SiN%tab_func(1:Para_H2SiN%ndim,j),Para_H2SiN%F(j),j=k+1,k+nb_columns)
 11        format(6i1,f15.8,5(2x,6i1,f15.8))
           k = k + nb_columns
         END DO
         read(nio_fit,*) Para_H2SiN%Qref(:)
         Para_H2SiN%Qref(6) = pi ! Qref(6) must be changed to be compatible with a z-matrix.
         close(nio_fit)
         deallocate(FileName)

      CASE Default

          write(out_unitp,*) ' ERROR in Init_H2SiN '
          write(out_unitp,*) ' This option is not possible. option: ',Para_H2SiN%option
          write(out_unitp,*) ' Its value MUST be 1'

          STOP
      END SELECT
    END IF

    IF (Para_H2SiN%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    END IF

  END SUBROUTINE Init_H2SiN
!> @brief Subroutine wich prints the H2SiN parameters.
!!
!! @param Para_H2SiN         TYPE(Param_H2SiN):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_H2SiN(Para_H2SiN,nio)
    
    TYPE(Param_H2SiN), intent(in) :: Para_H2SiN
    integer         , intent(in) :: nio

    write(nio,*) 'H2SiN current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '   Q(2) \ Q(3)                         '
    write(nio,*) '         Si------------N               '
    write(nio,*) '   Q(4) / Q(5)    Q(1)                 '
    write(nio,*) '       /                               '
    write(nio,*) '      H   +Q(6)                        '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Minimum:                             '
    write(nio,*) '  Q(1) = rSiN = 3.1103305087 (Bohr)    '
    write(nio,*) '  Q(2) = rH1  = 2.7870835493 (Bohr)    '
    write(nio,*) '  Q(3) = aH1  = 2.1336938250 (Radian)  '
    write(nio,*) '  Q(4) = rH2  = 2.7870835493 (Bohr)    '
    write(nio,*) '  Q(5) = aH3  = 2.1336938250 (Radian)  '
    write(nio,*) '  Q(6) = phi  =  pi (Radian)           '
    write(nio,*) '                                       '
    write(nio,*) '  V = -344.8478588900 Hartree          '
    write(nio,*) ' grad(:) =[-0.0000000073,-0.0000000118,'
    write(nio,*) '            0.0000000291,-0.0000000010,'
    write(nio,*) '            0.0000000479,0.0000000000] '
    write(nio,*) '                                       '
    write(nio,*) 'Ref:'
    write(nio,*) '  D. Lauvergnat, M.L. Senent, L. Jutier, and M. Hochlaf, JCP 135, 074301 (2011).'
    write(nio,*) '      https://doi.org/10.1063/1.3624563'
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',Para_H2SiN%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',Para_H2SiN%option
    write(nio,*)

    SELECT CASE (Para_H2SiN%option)

    CASE (1)

    CONTINUE

    CASE Default
        write(out_unitp,*) ' ERROR in write_H2SiNPot '
        write(out_unitp,*) ' This option is not possible. option: ',Para_H2SiN%option
        write(out_unitp,*) ' Its value MUST be 1'

        STOP
    END SELECT

    write(nio,*)
    write(nio,*) 'end H2SiN current parameters'

  END SUBROUTINE Write_H2SiN
!> @brief Subroutine wich prints the H2SiN parameters.
!!
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write0_H2SiN(nio)

    integer         , intent(in) :: nio

    write(nio,*) 'H2SiN default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end H2SiN default parameters'

  END SUBROUTINE Write0_H2SiN
!> @brief Subroutine wich calculates the H2SiN potential (for the 3 models) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_H2SiN         TYPE(Param_H2SiN):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_H2SiNPot(Mat_OF_PotDia,dnQ,Para_H2SiN,nderiv)
    USE mod_dnSca

    TYPE(Param_H2SiN) , intent(in)   :: Para_H2SiN
    TYPE(dnSca),       intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE(dnSca),       intent(in)    :: dnQ(:) !
    integer          , intent(in)    :: nderiv

    SELECT CASE (Para_H2SiN%option)

    CASE (1)
      CALL eval_H2SiNPot1(Mat_OF_PotDia,dnQ,Para_H2SiN,nderiv)

    CASE Default
        write(out_unitp,*) ' ERROR in eval_H2SiNPot '
        write(out_unitp,*) ' This option is not possible. option: ',Para_H2SiN%option
        write(out_unitp,*) ' Its value MUST be 1'

        STOP
    END SELECT

  END SUBROUTINE eval_H2SiNPot

!> @brief Subroutine wich calculates the H2SiN potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_H2SiN         TYPE(Param_H2SiN):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).

  SUBROUTINE eval_H2SiNPot1(Mat_OF_PotDia,dnQ,Para_H2SiN,nderiv)
    !Not Published model potential 
    USE mod_dnSca

    TYPE(dnSca),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE(dnSca),        intent(in)    :: dnQ(:)
    TYPE(Param_H2SiN) , intent(in)    :: Para_H2SiN
    integer,            intent(in)    :: nderiv


    TYPE(dnSca)        :: DQ(6,4)
    TYPE(dnSca)        :: Vtemp
    integer            :: i,j

    !write(6,*) ' sub eval_H2SiNPot1' ; flush(6)

      ! Warning, the coordinate ordering in the potential data (from the file) is different from the z-matrix one.
      DQ(:,1) = dnQ([2,4,1,3,5,6]) - Para_H2SiN%Qref(:)
      DO j=2,4
        DQ(:,j) = DQ(:,j-1) * DQ(:,1)
      END DO


      Mat_OF_PotDia(1,1) = ZERO
      Vtemp = Mat_OF_PotDia(1,1) ! to have a correct initialization

      DO j=1,Para_H2SiN%nb_func
        Vtemp = Para_H2SiN%F(j)
        DO i=1,Para_H2SiN%ndim
          !Vtemp = Vtemp * DQ(i,1)**Para_H2SiN%tab_func(i,j)
          IF (Para_H2SiN%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,Para_H2SiN%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO


!-----------------------------------------------------------------------!

   CALL dealloc_dnSca(Vtemp)
   CALL dealloc_dnSca(DQ)

   !write(6,*) ' end eval_H2SiNPot1' ; flush(6)

  END SUBROUTINE eval_H2SiNPot1

END MODULE mod_H2SiN
