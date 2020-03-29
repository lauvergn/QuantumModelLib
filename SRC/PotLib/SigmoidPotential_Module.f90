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

!> @brief Module which makes the initialization, calculation of the sigmoid function (tanh) (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
MODULE mod_SigmoidPot
  USE mod_NumParameters


  IMPLICIT NONE

!> @brief Derived type in which the sigmoid parameters are set-up.
!> @brief s(x) = A*0.5*(1+e*tanh( (x-B)/C )) (e=1 or -1)
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param A,B,C,e        real: sigmoid parameters
  TYPE SigmoidPot_t
     PRIVATE
     real (kind=Rkind) :: A   = ONE
     real (kind=Rkind) :: B   = ZERO
     real (kind=Rkind) :: C   = ONE
     real (kind=Rkind) :: e   = ONE
  END TYPE SigmoidPot_t

  PRIVATE  Read_SigmoidPot
CONTAINS

!> @brief Subroutine which makes the initialization of the Sigmoid parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param SigmoidPot       TYPE(SigmoidPot_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param A,B,C,e            real (optional):    sigmoid parameters
  SUBROUTINE Init_SigmoidPot(SigmoidPot,nio,read_param,A,B,C,e)
    TYPE (SigmoidPot_t),        intent(inout)   :: SigmoidPot
    real (kind=Rkind), optional, intent(in)      :: A,B,C,e
    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param


    logical :: read_param_loc

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_SigmoidPot'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'STOP in Init_SigmoidPot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    SigmoidPot = SigmoidPot_t(ONE,ZERO,ONE,ONE)

    IF (read_param) THEN
      CALL Read_SigmoidPot(SigmoidPot,nio)
    ELSE

      IF (present(A)) SigmoidPot%A = A
      IF (present(B)) SigmoidPot%B = B
      IF (present(C)) SigmoidPot%C = C
      IF (present(e)) SigmoidPot%e = e

    END IF

  END SUBROUTINE Init_SigmoidPot

!> @brief Subroutine wich reads the Sigmoid parameters with a namelist.
!!   This can be called only from the "Init_Sigmoid" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param SigmoidPot       TYPE(SigmoidPot_t):   derived type in which the parameters are set-up.
!! @param nio                integer            :   file unit to read the parameters.
  SUBROUTINE Read_SigmoidPot(SigmoidPot,nio)
    TYPE (SigmoidPot_t), intent(inout) :: SigmoidPot
    integer, intent(in) :: nio

    real (kind=Rkind) :: A,B,C,e
    integer           :: err_read


    namelist /Sigmoid/ A,B,C,e

    A = SigmoidPot%A
    B = SigmoidPot%B
    C = SigmoidPot%C
    e = SigmoidPot%e

    read(nio,nml=Sigmoid,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_SigmoidPot'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Sigmoid" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_SigmoidPot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_SigmoidPot'
      write(out_unitp,*) ' Some parameter names of the namelist "Sigmoid" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Sigmoid)
      STOP ' ERROR in Read_SigmoidPot'
    END IF
    !write(out_unitp,nml=Sigmoid)

    SigmoidPot = SigmoidPot_t(A,B,C,e)

  END SUBROUTINE Read_SigmoidPot

!> @brief Subroutine wich prints the Sigmoid parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param SigmoidPot       TYPE(SigmoidPot_t):   derived type with the Sigmoid parameters.
!! @param nio                integer            :   file unit to print the parameters.
  SUBROUTINE Write_SigmoidPot(SigmoidPot,nio)
    TYPE (SigmoidPot_t), intent(in) :: SigmoidPot
    integer, intent(in) :: nio

    write(nio,*) 'Sigmoid curent parameters'
    write(nio,*) '  s(x) = A*0.5*(1+e*tanh( (x-B)/C )) (e=1 or -1)'
    write(nio,*) '  A:   ',SigmoidPot%A
    write(nio,*) '  B:   ',SigmoidPot%B
    write(nio,*) '  C:   ',SigmoidPot%C
    write(nio,*) '  e:   ',SigmoidPot%e
    write(nio,*) 'end Sigmoid curent parameters'

  END SUBROUTINE Write_SigmoidPot
  SUBROUTINE get_Q0_Sigmoid(R0,SigmoidPot)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: R0
    TYPE (SigmoidPot_t),        intent(in)    :: SigmoidPot

    R0 = SigmoidPot%B

  END SUBROUTINE get_Q0_Sigmoid

!> @brief Subroutine wich calculates the Sigmoid potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param SigmoidPot       TYPE(SigmoidPot_t): derived type with the Sigmoid parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE Eval_SigmoidPot(Mat_OF_PotDia,dnR,SigmoidPot,nderiv)
    USE mod_dnS

    TYPE (SigmoidPot_t), intent(in)     :: SigmoidPot
    TYPE (dnS_t),          intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),          intent(in)     :: dnR
    integer, intent(in)                  :: nderiv


    !write(out_unitp,*) 'BEGINNING in Eval_SigmoidPot'
    !flush(out_unitp)


    Mat_OF_PotDia(1,1) = dnSigmoid(dnR,SigmoidPot)


    !write(out_unitp,*) 'END in Eval_SigmoidPot'
    !flush(out_unitp)

  END SUBROUTINE Eval_SigmoidPot

!> @brief Function wich calculates the Sigmoid potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return dnSigmoid         TYPE (dnS_t):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE (dnS_t):           derived type with the value of "r" and,if required, its derivatives.
!! @param SigmoidPot       TYPE(SigmoidPot_t): derived type with the Sigmoid parameters.
  FUNCTION dnSigmoid(dnR,SigmoidPot) !A*0.5*(1+e*tanh( (x-B)/C )) (e=1 or -1)
    USE mod_dnMat
    USE mod_dnS

    TYPE (dnS_t)                          :: dnSigmoid
    TYPE (dnS_t),          intent(in)     :: dnR

    TYPE (SigmoidPot_t), intent(in)   :: SigmoidPot

    !local variable
    TYPE (dnS_t)     :: dnbeta

    !write(out_unitp,*) 'BEGINNING in Sigmoid'
    !write(out_unitp,*) 'dnR'
    !CALL QML_Write_dnS(dnR)

    dnbeta  = tanh((dnR-SigmoidPot%B)/SigmoidPot%C)
    !write(out_unitp,*) 'dnbeta'
    !CALL QML_Write_dnS(dnbeta)
    dnSigmoid = (HALF*SigmoidPot%A) * (ONE+SigmoidPot%e*dnbeta)

     CALL QML_dealloc_dnS(dnbeta)

    !write(out_unitp,*) 'Sigmoid at',get_d0_FROM_dnS(dnR)
    !CALL QML_Write_dnS(dnSigmoid)
    !write(out_unitp,*) 'END in Sigmoid'
    !flush(out_unitp)

  END FUNCTION dnSigmoid
END MODULE mod_SigmoidPot
