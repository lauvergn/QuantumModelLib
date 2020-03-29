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
  TYPE BuckPot_t
     private
     real (kind=Rkind) :: A   = 387.63744459726228783977_Rkind  ! for Ar2 (eq 27)
     real (kind=Rkind) :: B   =   1.93837805257707347985_Rkind  ! for Ar2 (eq 27)
     real (kind=Rkind) :: C   = 106.54483566475760255666_Rkind  ! for Ar2 (eq 27)
     ! A= 1.69 10^-8 erg      = 387.63744459726228783977 Hartree
     ! B= 1/0.273 A^-1        = 1.93837805257707347985 bohr^-1
     ! C= 102 10^-12 erg A^-6 = 106.54483566475760255666 Hartree bohr^-6
  END TYPE BuckPot_t

  PRIVATE  Read_BuckPot

CONTAINS
!> @brief Subroutine which makes the initialization of the Buckingham parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param BuckPot          TYPE(BuckPot_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param A,B,C              real (optional):    parameters
  SUBROUTINE Init_BuckPot(BuckPot,read_param,nio,A,B,C)
    TYPE (BuckPot_t),           intent(inout)   :: BuckPot
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
    BuckPot = BuckPot_t(387.63744459726228783977_Rkind,1.93837805257707347985_Rkind,106.54483566475760255666_Rkind)

    IF (read_param_loc) THEN
      CALL Read_BuckPot(BuckPot,nio)
    ELSE
      IF (present(A))   BuckPot%A = A
      IF (present(B))   BuckPot%B = B
      IF (present(C))   BuckPot%C = C
    END IF

  END SUBROUTINE Init_BuckPot

!> @brief Subroutine wich reads the Buckingham parameters with a namelist.
!!   This can be called only from the "Init_BuckPot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param BuckPot          TYPE(BuckPot_t):    derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to read the parameters.
  SUBROUTINE Read_BuckPot(BuckPot,nio)
    TYPE (BuckPot_t), intent(inout) :: BuckPot
    integer, intent(in) :: nio

    real (kind=Rkind) :: A,B,C
    integer           :: err_read
    namelist /Buck/ A,B,C

    A = BuckPot%A
    B = BuckPot%B
    C = BuckPot%C

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

    BuckPot = BuckPot_t(A,B,C)


  END SUBROUTINE Read_BuckPot
!> @brief Subroutine wich prints the Buckingham parameters from his publication.
!!
!> @author David Lauvergnat
!! @date 20/07/2019
!!
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write0_BuckPot(nio)
    integer, intent(in) :: nio

    write(nio,*) 'Buckingham parameters:'
    write(nio,*) ' For Ar-Ar, eq 27 of reference:'
    write(nio,*) '  R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283.'
    write(nio,*)
    write(nio,*) 'Buckingham(r) = A.Exp(-B*r)-C/r^6'
    write(nio,*)
    write(nio,*) '   A= 1.69 10^-8 erg      = 387.63744459726228783977 Hartree'
    write(nio,*) '   B= 1/0.273 A^-1        = 1.93837805257707347985   bohr^-1'
    write(nio,*) '   C= 102 10^-12 erg A^-6 = 106.54483566475760255666 Hartree bohr^-6'
    write(nio,*)
    write(nio,*) 'Value at: R=7.0 Bohr'
    write(nio,*) 'V        = -0.000409 Hartree'
    write(nio,*) 'gradient = -0.000186'
    write(nio,*) 'hessian  =  0.001088'
    write(nio,*)
    write(nio,*) 'end Buckingham parameters'

  END SUBROUTINE Write0_BuckPot
!> @brief Subroutine wich prints the Buckingham parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param BuckPot          TYPE(BuckPot_t):    derived type with the Buckingham parameters.
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write_BuckPot(BuckPot,nio)
    TYPE (BuckPot_t), intent(in) :: BuckPot
    integer, intent(in) :: nio

    write(nio,*) 'Buckingham current parameters:'
    write(nio,*)
    write(nio,*) '    V(R) = A.Exp(-B.r) - C/r^6'
    write(nio,*) '  A:   ',BuckPot%A
    write(nio,*) '  B:   ',BuckPot%B
    write(nio,*) '  C:   ',BuckPot%C
    write(nio,*)
    write(nio,*) 'end Buckingham current parameters'

  END SUBROUTINE Write_BuckPot

  SUBROUTINE get_Q0_Buck(R0,BuckPot)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: R0
    TYPE (BuckPot_t),           intent(in)    :: BuckPot

    real (kind=Rkind) :: Rt1,Rt2
    integer           :: i

    Rt1 = TEN/BuckPot%B
    DO i=1,100
      !Rt2 = (BuckPot%B*BuckPot%A/(SIX*BuckPot%C)*exp(-BuckPot%B*Rt1))**(-ONE/SEVEN)
      Rt2 = -ONE/BuckPot%B* log(SIX*BuckPot%C/(BuckPot%B*BuckPot%A)*Rt1**(-7))
      !write(6,*) i,RT2
      IF (abs(Rt1-Rt2) < ONETENTH**10) EXIT
      Rt1 = Rt2
    END DO
    R0 = Rt1

  END SUBROUTINE get_Q0_Buck

!> @brief Subroutine wich calculates the Buckingham potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnS_t):         Matrix of dnS with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnR                TYPE (dnS_t):         derived type wich contain the value for which the potential is calculated: dnR%d0
!! @param BuckPot          TYPE(BuckPot_t):    derived type with the Buckingham parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE Eval_BuckPot(Mat_OF_PotDia,dnR,BuckPot,nderiv)
    USE mod_dnS

    TYPE (BuckPot_t), intent(in)     :: BuckPot
    TYPE (dnS_t),     intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),     intent(in)     :: dnR
    integer,          intent(in)     :: nderiv

    !write(out_unitp,*) 'BEGINNING in Eval_BuckPot'

    Mat_OF_PotDia(1,1) = dnBuck(dnR,BuckPot)

    !write(out_unitp,*) 'END in Eval_BuckPot'
    !flush(out_unitp)

  END SUBROUTINE Eval_BuckPot

!> @brief Function wich calculates the Buckingham potential with derivatives up to the 2d order is required.
!> @brief V(R) = A.Exp(-B*r)-C/r^6
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return dnBuck            TYPE (dnS_t):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE (dnS_t):           derived type with the value of "r" and,if required, its derivatives.
!! @param BuckPot          TYPE(BuckPot_t):    derived type with the Buckingham parameters.
  FUNCTION dnBuck(dnR,BuckPot)
    USE mod_dnMat
    USE mod_dnS


    TYPE (dnS_t)                          :: dnBuck
    TYPE (dnS_t),          intent(in)     :: dnR

    TYPE (BuckPot_t),    intent(in)     :: BuckPot


    !write(out_unitp,*) 'BEGINNING in dnBuck'
    !write(out_unitp,*) 'dnR'
    !CALL write_dnS(dnR)

    dnBuck = BuckPot%A * exp(-BuckPot%B*dnR) - BuckPot%C * dnR**(-6)

    !write(out_unitp,*) 'Buckingham at',get_d0_FROM_dnS(dnR)
    !CALL Write_dnS(dnBuck)
    !write(out_unitp,*) 'END in dnBuck'
    !flush(out_unitp)

  END FUNCTION dnBuck

END MODULE mod_BuckPot
