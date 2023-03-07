!===========================================================================
!===========================================================================
!This file is part of QuantumModelLib (QML).
!===============================================================================
! MIT License
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!    Copyright (c) 2022 David Lauvergnat [1]
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
  USE QDUtil_NumParameters_m, out_unit => out_unit
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
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) '  read_param:     ',read_param
      write(out_unit,*) '  nio_param_file: ',nio_param_file
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    IF (debug) write(out_unit,*) 'init default Sigmoid parameters'
    CALL Init0_QML_Sigmoid(QModel,A=ONE,B=ZERO,C=ONE,e=ONE,model_name='sigmoid')

    IF (read_param) THEN
      CALL Read_QML_Sigmoid(QModel,nio_param_file)
    ELSE
      IF (debug) write(out_unit,*) 'init Sigmoid parameters (A,B,C,e), if present'
      IF (present(A)) QModel%A = A
      IF (present(B)) QModel%B = B
      IF (present(C)) QModel%C = C
      IF (present(e)) QModel%e = e
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of Sigmoid'
    QModel%Q0 = [QModel%B]

    IF (debug) write(out_unit,*) 'init d0GGdef of Sigmoid'
    QModel%d0GGdef = reshape([ONE],shape=[1,1])

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
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
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%In_a_Model = .TRUE.
    QModel%ndim       = 1
    QModel%nsurf      = 1

    QModel%pot_name = 'sigmoid'
    IF (present(model_name)) QModel%pot_name = model_name

    IF (debug) write(out_unit,*) 'init Sigmoid parameters (A,B,C,e), if present'

    IF (present(A))   QModel%A = A
    IF (present(B))   QModel%B = B
    IF (present(C))   QModel%C = C
    IF (present(e))   QModel%e = e

    IF (debug) THEN
      write(out_unit,*) ',A,B,C,e: ',QModel%A,QModel%B,QModel%C,QModel%e
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
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
    IMPLICIT NONE

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
      write(out_unit,*) ' ERROR in Read_QML_Sigmoid'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "Sigmoid" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Sigmoid'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Sigmoid'
      write(out_unit,*) ' Some parameter names of the namelist "Sigmoid" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=Sigmoid)
      STOP ' ERROR in Read_QML_Sigmoid'
    END IF
    !write(out_unit,nml=Sigmoid)

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
    IMPLICIT NONE

    CLASS(QML_Sigmoid_t),    intent(in) :: QModel
    integer,                  intent(in) :: nio

    write(nio,*) 'Sigmoid current parameters'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*) '  s(x) = A*0.5*(1+e*tanh( (x-B)/C )) (e=1 or -1)'
    write(nio,*) '  A:   ',QModel%A
    write(nio,*) '  B:   ',QModel%B
    write(nio,*) '  C:   ',QModel%C
    write(nio,*) '  e:   ',QModel%e
    write(nio,*) 'end Sigmoid current parameters'

  END SUBROUTINE Write_QML_Sigmoid
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
    USE ADdnSVM_m, ONLY :  dnS_t
    IMPLICIT NONE

    CLASS (QML_Sigmoid_t),  intent(in)     :: QModel
    TYPE (dnS_t),           intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),           intent(in)     :: dnQ(:)
    integer,                intent(in)     :: nderiv


    !write(out_unit,*) 'BEGINNING in EvalPot_QML_Sigmoid'
    !flush(out_unit)

    Mat_OF_PotDia(1,1) = QML_dnSigmoid(dnQ(1),QModel)


    !write(out_unit,*) 'END in EvalPot_QML_Sigmoid'
    !flush(out_unit)

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
    USE ADdnSVM_m, ONLY :  dnS_t,tanh,dealloc_dnS,Write_dnS
    IMPLICIT NONE

    TYPE (dnS_t)                          :: QML_dnSigmoid
    TYPE (dnS_t),          intent(in)     :: dnR

    TYPE (QML_Sigmoid_t), intent(in)   :: SigmoidPot

    !local variable
    TYPE (dnS_t)     :: dnbeta

    !write(out_unit,*) 'BEGINNING in Sigmoid'
    !write(out_unit,*) 'dnR'
    !CALL Write_dnS(dnR)

    dnbeta  = tanh((dnR-SigmoidPot%B)/SigmoidPot%C)
    !write(out_unit,*) 'dnbeta'
    !CALL Write_dnS(dnbeta)
    QML_dnSigmoid = (HALF*SigmoidPot%A) * (ONE+SigmoidPot%e*dnbeta)

     CALL dealloc_dnS(dnbeta)

    !write(out_unit,*) 'Sigmoid at',get_d0_FROM_dnS(dnR)
    !CALL Write_dnS(QML_dnSigmoid)
    !write(out_unit,*) 'END in Sigmoid'
    !flush(out_unit)

  END FUNCTION QML_dnSigmoid
END MODULE QML_Sigmoid_m
