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
!> @brief Module which makes the initialization, calculation of the sigmoid function (tanh) (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
MODULE QML_Sigmoid_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the sigmoid parameters are set-up.
!> @brief s(x) = A*0.5*(1+e*tanh( (x-B)/C )) (e=1 or -1)
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param A,B,C,e        real: sigmoid parameters
  TYPE, EXTENDS (QML_Empty_t) :: QML_Sigmoid_t
     PRIVATE
     real (kind=Rkind) :: A   = ONE
     real (kind=Rkind) :: B   = ZERO ! equivalent to Req
     real (kind=Rkind) :: C   = ONE
     real (kind=Rkind) :: e   = ONE
  CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_Sigmoid
    PROCEDURE :: Write_QModel    => Write_QML_Sigmoid
    PROCEDURE :: Write0_QModel   => Write0_QML_Sigmoid
  END TYPE QML_Sigmoid_t

  PUBLIC :: QML_Sigmoid_t,Init_QML_Sigmoid,Init0_QML_Sigmoid,Write_QML_Sigmoid,QML_dnSigmoid
CONTAINS

!> @brief Subroutine which makes the initialization of the Sigmoid parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel             TYPE(QML_Sigmoid_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param A,B,C,e            real (optional):    sigmoid parameters
  FUNCTION Init_QML_Sigmoid(QModel_in,read_param,nio_param_file,A,B,C,e) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_Sigmoid_t)                        :: QModel

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: A,B,C,e

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Sigmoid'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) '  read_param:     ',read_param
      write(out_unitp,*) '  nio_param_file: ',nio_param_file
      flush(out_unitp)
    END IF


    !QModel_loc%QML_Empty_t = Init_QML_Empty(QModel_in) ! it does not work with nagfor
    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    IF (debug) write(out_unitp,*) 'init default Sigmoid parameters'
    CALL Init0_QML_Sigmoid(QModel,A=ONE,B=ZERO,C=ONE,e=ONE,model_name='sigmoid')

    IF (read_param) THEN
      CALL Read_QML_Sigmoid(QModel,nio_param_file)
    ELSE
      IF (debug) write(out_unitp,*) 'init Sigmoid parameters (A,B,C,e), if present'
      IF (present(A)) QModel%A = A
      IF (present(B)) QModel%B = B
      IF (present(C)) QModel%C = C
      IF (present(e)) QModel%e = e
    END IF

    IF (debug) write(out_unitp,*) 'init Q0 of Sigmoid'
    QModel%Q0 = [QModel%B]

    IF (debug) write(out_unitp,*) 'init d0GGdef of Sigmoid'
    QModel%d0GGdef = reshape([ONE],shape=[1,1])

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_Sigmoid
  SUBROUTINE Init0_QML_Sigmoid(QModel,A,B,C,e,model_name)
  IMPLICIT NONE

    TYPE (QML_Sigmoid_t),           intent(inout)   :: QModel
    real (kind=Rkind),   optional,   intent(in)      :: A,B,C,e
    character (len=*),   optional,   intent(in)      :: model_name

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init0_QML_Sigmoid'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    QModel%ndim     = 1
    QModel%nsurf    = 1

    QModel%pot_name = 'sigmoid'
    IF (present(model_name)) QModel%pot_name = model_name

    IF (debug) write(out_unitp,*) 'init Sigmoid parameters (A,B,C,e), if present'

    IF (present(A))   QModel%A = A
    IF (present(B))   QModel%B = B
    IF (present(C))   QModel%C = C
    IF (present(e))   QModel%e = e

    IF (debug) THEN
      write(out_unitp,*) ',A,B,C,e: ',QModel%A,QModel%B,QModel%C,QModel%e
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Init0_QML_Sigmoid


!> @brief Subroutine wich reads the Sigmoid parameters with a namelist.
!!   This can be called only from the "Init_Sigmoid" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel       TYPE(QML_Sigmoid_t):   derived type in which the parameters are set-up.
!! @param nio                integer            :   file unit to read the parameters.
  SUBROUTINE Read_QML_Sigmoid(QModel,nio)
    TYPE (QML_Sigmoid_t), intent(inout) :: QModel
    integer, intent(in) :: nio

    real (kind=Rkind) :: A,B,C,e
    integer           :: err_read


    namelist /Sigmoid/ A,B,C,e

    A = QModel%A
    B = QModel%B
    C = QModel%C
    e = QModel%e

    read(nio,nml=Sigmoid,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_QML_Sigmoid'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Sigmoid" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_QML_Sigmoid'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_QML_Sigmoid'
      write(out_unitp,*) ' Some parameter names of the namelist "Sigmoid" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Sigmoid)
      STOP ' ERROR in Read_QML_Sigmoid'
    END IF
    !write(out_unitp,nml=Sigmoid)

    QModel%A = A
    QModel%B = B
    QModel%C = C
    QModel%e =e

  END SUBROUTINE Read_QML_Sigmoid

!> @brief Subroutine wich prints the Sigmoid parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel             CLASS(QML_Sigmoid_t):   derived type with the Sigmoid parameters.
!! @param nio                integer            :   file unit to print the parameters.
  SUBROUTINE Write_QML_Sigmoid(QModel,nio)
    CLASS(QML_Sigmoid_t),    intent(in) :: QModel
    integer,                  intent(in) :: nio

    write(nio,*) 'Sigmoid current parameters'
    CALL Write_QML_Empty(QModel%QML_Empty_t,nio)
    write(nio,*) '  s(x) = A*0.5*(1+e*tanh( (x-B)/C )) (e=1 or -1)'
    write(nio,*) '  A:   ',QModel%A
    write(nio,*) '  B:   ',QModel%B
    write(nio,*) '  C:   ',QModel%C
    write(nio,*) '  e:   ',QModel%e
    write(nio,*) 'end Sigmoid current parameters'

  END SUBROUTINE Write_QML_Sigmoid
  SUBROUTINE Write0_QML_Sigmoid(QModel,nio)
    CLASS(QML_Sigmoid_t),    intent(in) :: QModel
    integer,                  intent(in) :: nio

    write(nio,*) 'Sigmoid default parameters'
    write(nio,*) '  s(x) = A*0.5*(1+e*tanh( (x-B)/C )) (e=1 or -1)'
    write(nio,*) '  A:   1.'
    write(nio,*) '  B:   0.'
    write(nio,*) '  C:   1.'
    write(nio,*) '  e:   1.'
    write(nio,*) 'end Sigmoid default parameters'

  END SUBROUTINE Write0_QML_Sigmoid
!> @brief Subroutine wich calculates the Sigmoid potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             CLASS(QML_Sigmoid_t): derived type with the Sigmoid parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Sigmoid(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE QMLdnSVM_dnS_m

    CLASS (QML_Sigmoid_t), intent(in)     :: QModel
    TYPE (dnS_t),           intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),           intent(in)     :: dnQ(:)
    integer,                intent(in)     :: nderiv


    !write(out_unitp,*) 'BEGINNING in EvalPot_QML_Sigmoid'
    !flush(out_unitp)

    Mat_OF_PotDia(1,1) = QML_dnSigmoid(dnQ(1),QModel)


    !write(out_unitp,*) 'END in EvalPot_QML_Sigmoid'
    !flush(out_unitp)

  END SUBROUTINE EvalPot_QML_Sigmoid

!> @brief Function wich calculates the Sigmoid potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return QML_dnSigmoid         TYPE (dnS_t):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE (dnS_t):           derived type with the value of "r" and,if required, its derivatives.
!! @param SigmoidPot       TYPE(QML_Sigmoid_t): derived type with the Sigmoid parameters.
  FUNCTION QML_dnSigmoid(dnR,SigmoidPot) !A*0.5*(1+e*tanh( (x-B)/C )) (e=1 or -1)
    USE QMLdnSVM_dnMat_m
    USE QMLdnSVM_dnS_m

    TYPE (dnS_t)                          :: QML_dnSigmoid
    TYPE (dnS_t),          intent(in)     :: dnR

    TYPE (QML_Sigmoid_t), intent(in)   :: SigmoidPot

    !local variable
    TYPE (dnS_t)     :: dnbeta

    !write(out_unitp,*) 'BEGINNING in Sigmoid'
    !write(out_unitp,*) 'dnR'
    !CALL QML_Write_dnS(dnR)

    dnbeta  = tanh((dnR-SigmoidPot%B)/SigmoidPot%C)
    !write(out_unitp,*) 'dnbeta'
    !CALL QML_Write_dnS(dnbeta)
    QML_dnSigmoid = (HALF*SigmoidPot%A) * (ONE+SigmoidPot%e*dnbeta)

     CALL QML_dealloc_dnS(dnbeta)

    !write(out_unitp,*) 'Sigmoid at',get_d0_FROM_dnS(dnR)
    !CALL QML_Write_dnS(QML_dnSigmoid)
    !write(out_unitp,*) 'END in Sigmoid'
    !flush(out_unitp)

  END FUNCTION QML_dnSigmoid
END MODULE QML_Sigmoid_m
