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
!    Copyright 2016 David Lauvergnat [1]
!      with contributions of:
!        Félix MOUHAT [2]
!        Liang LIANG [3]
!        Emanuele MARSILI [1,4]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France
![3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France
![4]: Durham University, Durham, UK
!* Originally, it has been developed during the Quantum-Dynamics E-CAM project :
!     https://www.e-cam2020.eu/quantum-dynamics
!
!===========================================================================
!===========================================================================
!> @brief Module which makes the initialization, calculation of the Buckingham potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
MODULE QML_Buck_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Buckingham parameters are set-up.
!> @brief V(R) = A.Exp(-B*r)-C/r^6
!> @brief Default parameters for Ar-Ar
!> @brief Reference: R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283. doi:10.1098/rspa.1938.0173.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param A,B,C              real: Buckingham parameters
  TYPE, EXTENDS (QML_Empty_t) :: QML_Buck_t
     private
     real (kind=Rkind) :: A   = 387.63744459726228783977_Rkind  ! for Ar2 (eq 27)
     real (kind=Rkind) :: B   =   1.93837805257707347985_Rkind  ! for Ar2 (eq 27)
     real (kind=Rkind) :: C   = 106.54483566475760255666_Rkind  ! for Ar2 (eq 27)
     ! A= 1.69 10^-8 erg      = 387.63744459726228783977 Hartree
     ! B= 1/0.273 A^-1        = 1.93837805257707347985 bohr^-1
     ! C= 102 10^-12 erg A^-6 = 106.54483566475760255666 Hartree bohr^-6
     real (kind=Rkind) :: mu  = 36423.484024390622_Rkind        ! au
  CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_Buck
    PROCEDURE :: Write_QModel    => Write_QML_Buck
    PROCEDURE :: Write0_QModel   => Write0_QML_Buck
  END TYPE QML_Buck_t

  PUBLIC :: QML_Buck_t,Init_QML_Buck,Init0_QML_Buck,Write_QML_Buck,QML_dnBuck

CONTAINS
!> @brief Subroutine which makes the initialization of the Buckingham parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel          TYPE(QML_Buck_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param A,B,C              real (optional):    parameters
  FUNCTION Init_QML_Buck(QModel_in,read_param,nio_param_file,A,B,C) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_Buck_t)                           :: QModel

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: A,B,C

    !local variables
    real (kind=Rkind) :: R0

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Buck'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    !Default for Ar-Ar
    CALL Init0_QML_Buck(QModel,A=387.63744459726228783977_Rkind,       &
                                 B=1.93837805257707347985_Rkind,        &
                                 C=106.54483566475760255666_Rkind,      &
                                 model_name='buck')

    IF (read_param) THEN
      CALL Read_QML_Buck(QModel,nio_param_file)
    ELSE
      IF (present(A))   QModel%A = A
      IF (present(B))   QModel%B = B
      IF (present(C))   QModel%C = C
    END IF

    IF (debug) write(out_unitp,*) 'init Q0 of Buck'
    CALL get_Q0_QML_Buck(R0,QModel) ! no analytical value
    QModel%Q0 = [R0]

    IF (debug) write(out_unitp,*) 'init d0GGdef of Buck'
    QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])

    IF (debug) THEN
      write(out_unitp,*) 'A,B,C: ',QModel%A,QModel%B,QModel%C
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_Buck
  SUBROUTINE Init0_QML_Buck(QModel,A,B,C,model_name)
  IMPLICIT NONE

    TYPE (QML_Buck_t),            intent(inout)   :: QModel
    real (kind=Rkind), optional,   intent(in)      :: A,B,C
    character (len=*), optional,   intent(in)      :: model_name

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init0_QML_Buck'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    QModel%In_a_Model = .TRUE.

    QModel%ndim       = 1
    QModel%nsurf      = 1
    QModel%pot_name   = 'Buck'
    IF (present(model_name)) QModel%pot_name = model_name

    IF (debug) write(out_unitp,*) 'init Buck parameters (A,B,C), if present'

    IF (present(A))   QModel%A   = A
    IF (present(B))   QModel%B   = B
    IF (present(C))   QModel%C   = C

    IF (debug) THEN
      write(out_unitp,*) 'A,B,C: ',QModel%A,QModel%B,QModel%C
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Init0_QML_Buck
!> @brief Subroutine wich reads the Buckingham parameters with a namelist.
!!   This can be called only from the "Init_QML_Buck" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel          TYPE(QML_Buck_t):    derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to read the parameters.
  SUBROUTINE Read_QML_Buck(QModel,nio)
    TYPE (QML_Buck_t), intent(inout) :: QModel
    integer,            intent(in)    :: nio

    real (kind=Rkind) :: A,B,C
    integer           :: err_read
    namelist /Buck/ A,B,C

    A = QModel%A
    B = QModel%B
    C = QModel%C

    read(nio,nml=Buck,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_QML_Buck'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Buck" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_QML_Buck'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_QML_Buck'
      write(out_unitp,*) ' Some parameter names of the namelist "Buck" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Buck)
      STOP ' ERROR in Read_QML_Buck'
    END IF

    !write(out_unitp,nml=Buck)
    QModel%A = A
    QModel%B = B
    QModel%C = C

  END SUBROUTINE Read_QML_Buck
!> @brief Subroutine wich prints the Buckingham parameters from his publication.
!!
!> @author David Lauvergnat
!! @date 20/07/2019
!!
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write0_QML_Buck(QModel,nio)
    CLASS (QML_Buck_t), intent(in) :: QModel
    integer,             intent(in) :: nio

    write(nio,*) 'Buckingham default parameters:'
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

  END SUBROUTINE Write0_QML_Buck
!> @brief Subroutine wich prints the Buckingham parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel          TYPE(QML_Buck_t):    derived type with the Buckingham parameters.
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write_QML_Buck(QModel,nio)
    CLASS (QML_Buck_t), intent(in) :: QModel
    integer,             intent(in) :: nio

    write(nio,*) 'Buckingham current parameters:'
    write(nio,*)
    write(nio,*) '    V(R) = A.Exp(-B.r) - C/r^6'
    write(nio,*) '  A:   ',QModel%A
    write(nio,*) '  B:   ',QModel%B
    write(nio,*) '  C:   ',QModel%C
    write(nio,*)
    write(nio,*) 'end Buckingham current parameters'

  END SUBROUTINE Write_QML_Buck

  SUBROUTINE get_Q0_QML_Buck(R0,QModel)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: R0
    TYPE (QML_Buck_t),           intent(in)    :: QModel

    real (kind=Rkind) :: Rt1,Rt2
    integer           :: i

    Rt1 = TEN/QModel%B
    DO i=1,100
      !Rt2 = (QModel%B*QModel%A/(SIX*QModel%C)*exp(-QModel%B*Rt1))**(-ONE/SEVEN)
      Rt2 = -ONE/QModel%B* log(SIX*QModel%C/(QModel%B*QModel%A)*Rt1**(-7))
      !write(6,*) i,RT2
      IF (abs(Rt1-Rt2) < ONETENTH**10) EXIT
      Rt1 = Rt2
    END DO
    R0 = Rt1

  END SUBROUTINE get_Q0_QML_Buck

!> @brief Subroutine wich calculates the Buckingham potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnS_t):         Matrix of dnS with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnR                TYPE (dnS_t):         derived type wich contain the value for which the potential is calculated: dnR%d0
!! @param QModel          TYPE(QML_Buck_t):    derived type with the Buckingham parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Buck(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS(QML_Buck_t),  intent(in)     :: QModel
    TYPE(dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         intent(in)     :: dnQ(:)
    integer,             intent(in)     :: nderiv

    !write(out_unitp,*) 'BEGINNING in EvalPot_QML_Buck'

    Mat_OF_PotDia(1,1) = QML_dnBuck(dnQ(1),QModel)

    !write(out_unitp,*) 'END in EvalPot_QML_Buck'
    !flush(out_unitp)

  END SUBROUTINE EvalPot_QML_Buck

!> @brief Function wich calculates the Buckingham potential with derivatives up to the 2d order is required.
!> @brief V(R) = A.Exp(-B*r)-C/r^6
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return QML_dnBuck            TYPE (dnS_t):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE (dnS_t):           derived type with the value of "r" and,if required, its derivatives.
!! @param BuckPot          TYPE(QML_Buck_t):    derived type with the Buckingham parameters.
  FUNCTION QML_dnBuck(dnR,BuckPot)
    USE QMLdnSVM_dnMat_m
    USE QMLdnSVM_dnS_m


    TYPE (dnS_t)                          :: QML_dnBuck
    TYPE (dnS_t),          intent(in)     :: dnR

    TYPE (QML_Buck_t),    intent(in)     :: BuckPot


    !write(out_unitp,*) 'BEGINNING in QML_dnBuck'
    !write(out_unitp,*) 'dnR'
    !CALL QML_Write_dnS(dnR)

    QML_dnBuck = BuckPot%A * exp(-BuckPot%B*dnR) - BuckPot%C * dnR**(-6)

    !write(out_unitp,*) 'Buckingham at',get_d0_FROM_dnS(dnR)
    !CALL QML_Write_dnS(QML_dnBuck)
    !write(out_unitp,*) 'END in QML_dnBuck'
    !flush(out_unitp)

  END FUNCTION QML_dnBuck

END MODULE QML_Buck_m
