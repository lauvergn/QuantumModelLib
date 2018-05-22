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

!> @brief Module which makes the initialization, calculation of the Buckingham potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
MODULE mod_BuckPot
  USE mod_NumParameters
  IMPLICIT NONE

!> @brief Derived type in which the Buckingham parameters are set-up.
!> @brief V(R) = A.Exp(-B*r)-C/r^6
!> @brief Default parameters for Ar-Ar
!> @brief Reference: R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283. doi:10.1098/rspa.1938.0173.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param A,B,C              real: Buckingham parameters
  TYPE Param_Buck
     private
     real (kind=Rkind) :: A   = 387.63744459726228783977_Rkind  ! for Ar2 (eq 27)
     real (kind=Rkind) :: B   =   1.93837805257707347985_Rkind  ! for Ar2 (eq 27)
     real (kind=Rkind) :: C   = 106.54483566475760255666_Rkind  ! for Ar2 (eq 27)
     ! A= 1.69 10^-8 erg      = 387.63744459726228783977 Hartree
     ! B= 1/0.273 A^-1        = 1.93837805257707347985 bohr^-1
     ! C= 102 10^-12 erg A^-6 = 106.54483566475760255666 Hartree bohr^-6
  END TYPE Param_Buck

  PRIVATE  Read_BuckPot

CONTAINS
!> @brief Subroutine which makes the initialization of the Buckingham parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Buck          TYPE(Param_Buck):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param A,B,C              real (optional):    parameters
  SUBROUTINE Init_BuckPot(Para_Buck,read_param,nio,A,B,C)
    TYPE (Param_Buck),           intent(inout)   :: Para_Buck
    real (kind=Rkind), optional, intent(in)      :: A,B,C

    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param


    logical :: read_param_loc

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_BuckPot'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'ERROR in Init_BuckPot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    !Default for Ar-Ar
    Para_Buck = Param_Buck(387.63744459726228783977_Rkind,1.93837805257707347985_Rkind,106.54483566475760255666_Rkind)

    IF (read_param_loc) THEN
      CALL Read_BuckPot(Para_Buck,nio)
    ELSE
      IF (present(A))   Para_Buck%A = A
      IF (present(B))   Para_Buck%B = B
      IF (present(C))   Para_Buck%C = C
    END IF

  END SUBROUTINE Init_BuckPot

!> @brief Subroutine wich reads the Buckingham parameters with a namelist.
!!   This can be called only from the "Init_BuckPot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Buck          TYPE(Param_Buck):    derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to read the parameters.
  SUBROUTINE Read_BuckPot(Para_Buck,nio)
    TYPE (Param_Buck), intent(inout) :: Para_Buck
    integer, intent(in) :: nio

    real (kind=Rkind) :: A,B,C
    integer           :: err_read
    namelist /Buck/ A,B,C

    A = Para_Buck%A
    B = Para_Buck%B
    C = Para_Buck%C

    read(nio,nml=Buck,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_BuckPot'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Buck" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_BuckPot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_BuckPot'
      write(out_unitp,*) ' Some parameter names of the namelist "Buck" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Buck)
      STOP ' ERROR in Read_BuckPot'
    END IF

    !write(out_unitp,nml=Buck)

    Para_Buck = Param_Buck(A,B,C)


  END SUBROUTINE Read_BuckPot
!> @brief Subroutine wich prints the Buckingham parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Buck          TYPE(Param_Buck):    derived type with the Buckingham parameters.
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write_BuckPot(Para_Buck,nio)
    TYPE (Param_Buck), intent(in) :: Para_Buck
    integer, intent(in) :: nio

    write(nio,*) 'Buckingham parameters:'
    write(nio,*) ' For Ar-Ar, eq 27 of reference:'

    write(nio,*) 'R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283.'
    write(nio,*) 'Value at: R=7.0 Bohr'
    write(nio,*) 'V        = -0.000409 Hartree'
    write(nio,*) 'gradient = -0.000186'
    write(nio,*) 'hessian  =  0.001088'

    write(nio,*)
    write(nio,*) 'Current parameters:'
    write(nio,*) '    V(R) = A.Exp(-B.r) - C/r^6'
    write(nio,*) '  A:   ',Para_Buck%A
    write(nio,*) '  B:   ',Para_Buck%B
    write(nio,*) '  B:   ',Para_Buck%C
    write(nio,*) 'end Buckingham parameters'

  END SUBROUTINE Write_BuckPot

!> @brief Subroutine wich calculates the Buckingham potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_Buck          TYPE(Param_Buck):    derived type with the Buckingham parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE Eval_BuckPot(PotVal,r,Para_Buck,nderiv)
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_Buck), intent(in)     :: Para_Buck
    TYPE(dnMatPot),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: r
    integer,           intent(in)     :: nderiv

    !local variables (derived type). They have to be deallocated
    TYPE(dnSca)     :: dnPot,dnR

    !write(out_unitp,*) 'BEGINNING in Eval_BuckPot'
    !flush(out_unitp)

    IF ( Check_NotAlloc_dnMatPot(PotVal,nderiv) ) THEN
      CALL alloc_dnMatPot(PotVal,nsurf=1,ndim=1,nderiv=nderiv)
    END IF

    dnR     = init_dnSca(r,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives

    dnPot = dnBuck(dnR,Para_Buck)

    !transfert the 1D-potential and its derivatives (dnPot) to the matrix form (PotVal)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=1)


    CALL dealloc_dnSca(dnPot)
    CALL dealloc_dnSca(dnR)

    !write(out_unitp,*) 'Buckingham PotVal at',r
    !CALL Write_dnMatPot(PotVal)
    !write(out_unitp,*) 'END in Eval_BuckPot'
    !flush(out_unitp)

  END SUBROUTINE Eval_BuckPot

  SUBROUTINE Eval_BuckPot_old(PotVal,r,Para_Buck,nderiv)
    USE mod_dnMatPot

    TYPE (Param_Buck), intent(in)    :: Para_Buck
    TYPE(dnMatPot), intent(inout)    :: PotVal
    real (kind=Rkind),intent(in)     :: r
    integer, intent(in)              :: nderiv


    ! Potential calculation
    PotVal%d0(1,1) = Para_Buck%A*exp(-Para_Buck%B*r)-Para_Buck%C/r**6

    ! gradient calculation
    if (nderiv >= 1) then
      PotVal%d1(1,1,1) = -Para_Buck%A*Para_Buck%B*exp(-Para_Buck%B*r) + &
                          SIX*Para_Buck%C/r**7
    endif

    ! Hessian calculation
    if (nderiv >= 2) then
      PotVal%d2(1,1,1,1) = Para_Buck%A*Para_Buck%B**2*exp(-Para_Buck%B*r) - &
                          42_Rkind*Para_Buck%C/r**8
    endif

  END SUBROUTINE Eval_BuckPot_old

!> @brief Function wich calculates the Buckingham potential with derivatives up to the 2d order is required.
!> @brief V(R) = A.Exp(-B*r)-C/r^6
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return dnBuck            TYPE(dnSca):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE(dnSca):           derived type with the value of "r" and,if required, its derivatives.
!! @param Para_Buck          TYPE(Param_Buck):    derived type with the Buckingham parameters.
  FUNCTION dnBuck(dnR,Para_Buck)
    USE mod_dnMatPot
    USE mod_dnSca


    TYPE(dnSca)                          :: dnBuck
    TYPE(dnSca),          intent(in)     :: dnR

    TYPE (Param_Buck), intent(in)      :: Para_Buck


    !write(out_unitp,*) 'BEGINNING in dnBuck'
    !write(out_unitp,*) 'dnR'
    !CALL write_dnSca(dnR)

    dnBuck = Para_Buck%A * exp(-Para_Buck%B*dnR) - Para_Buck%C * dnR**(-6)

    !write(out_unitp,*) 'Buckingham at',get_d0_FROM_dnSca(dnR)
    !CALL Write_dnSca(dnBuck)
    !write(out_unitp,*) 'END in dnBuck'
    !flush(out_unitp)

  END FUNCTION dnBuck

END MODULE mod_BuckPot
