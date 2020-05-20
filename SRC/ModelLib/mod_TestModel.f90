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

!> @brief Module which makes the initialization, calculation of the test potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_TestModel
  USE mod_NumParameters
  USE mod_EmptyModel
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the test parameters are set-up.
  TYPE, EXTENDS (EmptyModel_t) ::  TestModel_t

   PRIVATE

   CONTAINS
    PROCEDURE :: Eval_QModel_Pot => eval_TestPot
    PROCEDURE :: Write_QModel    => Write_TestModel
    PROCEDURE :: Write0_QModel   => Write0_TestModel
  END TYPE TestModel_t
 
  PUBLIC :: TestModel_t,Init_TestModel
 
  CONTAINS
!> @brief Function which makes the initialization of the test parameters.
!!
!! @param QModel             TYPE(TestModel_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(EmptyModel_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_TestModel(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (TestModel_t)                           :: QModel ! RESULT

    TYPE(EmptyModel_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_TestModel'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_EmptyModel(QModel%EmptyModel_t,QModel_in)

    QModel%nsurf    = 2
    QModel%ndim     = 1
    QModel%pot_name = 'test'


    IF (debug) write(out_unitp,*) 'init Q0 of test'
    QModel%Q0 = [ZERO,ZERO,ZERO]

    IF (debug) write(out_unitp,*) 'init d0GGdef of test'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_TestModel
!> @brief Subroutine wich prints the current TestModel parameters.
!!
!! @param QModel            CLASS(TestModel_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_TestModel(QModel,nio)
  IMPLICIT NONE

    CLASS(TestModel_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

  END SUBROUTINE Write_TestModel
!> @brief Subroutine wich prints the default TestModel parameters.
!!
!! @param QModel            CLASS(TestModel_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_TestModel(QModel,nio)
  IMPLICIT NONE

    CLASS(TestModel_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'test default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end test default parameters'


  END SUBROUTINE Write0_TestModel

!> @brief Subroutine wich calculates the test potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(TestModel_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TestPot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(TestModel_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    Mat_OF_PotDia(1,1) =  ONE+dnQ(1)**1
    Mat_OF_PotDia(2,2) =  -Mat_OF_PotDia(1,1)
    Mat_OF_PotDia(1,2) =  ONE
    Mat_OF_PotDia(2,1) =  ONE

!Mat_OF_PotDia(1,1) = ONE                + dnQ(1)+dnQ(1)**2+dnQ(1)**3+dnQ(1)**4+dnQ(1)**5
!Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + dnQ(2)+dnQ(2)**2+dnQ(2)**3+dnQ(2)**4+dnQ(2)**5
!Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + dnQ(3)+dnQ(3)**2+dnQ(3)**3+dnQ(3)**4+dnQ(3)**5
!
!Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + dnQ(1)*dnQ(2)+dnQ(1)*dnQ(3)+dnQ(2)*dnQ(3)
!
!Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + dnQ(1)**2*dnQ(2)+TWO*dnQ(1)**2*dnQ(3)+THREE*dnQ(2)**2*dnQ(3)
!Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + HALF*dnQ(1)*dnQ(2)**2+FIVE*dnQ(1)*dnQ(3)**2+FOUR*dnQ(2)*dnQ(3)**2
!
!!Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + dnQ(1)*dnQ(2)*dnQ(3)

  END SUBROUTINE eval_TestPot

END MODULE mod_TestModel
