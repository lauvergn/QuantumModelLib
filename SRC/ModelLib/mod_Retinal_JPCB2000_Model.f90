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

!> @brief Module which makes the initialization, calculation of the Retinal_JPCB2000 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_Retinal_JPCB2000_Model
  USE mod_NumParameters
  USE mod_EmptyModel
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Retinal_JPCB2000 parameters are set-up.
  TYPE, EXTENDS (EmptyModel_t) ::  Retinal_JPCB2000_Model_t

   PRIVATE

   !The parameters of the model are in eV:
   !  m􏰂=4.84E-4,E1=􏰌2.48,W0=􏰌3.6,W1=1.09,w=0.19,kappa=􏰌0.1,and lambda=􏰌0.19.

   real(kind=Rkind) :: m      = ONE/4.84e-4_Rkind
   real(kind=Rkind) :: E1     = 2.48_Rkind
   real(kind=Rkind) :: W0     = 3.6_Rkind
   real(kind=Rkind) :: W1     = 1.09_Rkind
   real(kind=Rkind) :: w      = 0.19_Rkind
   real(kind=Rkind) :: kappa  = 0.1_Rkind
   real(kind=Rkind) :: lambda = 0.19_Rkind


   CONTAINS
    PROCEDURE :: Eval_QModel_Pot => eval_Retinal_JPCB2000_Pot
    PROCEDURE :: Write_QModel    => Write_Retinal_JPCB2000_Model
    PROCEDURE :: Write0_QModel   => Write0_Retinal_JPCB2000_Model
  END TYPE Retinal_JPCB2000_Model_t
 
  PUBLIC :: Retinal_JPCB2000_Model_t,Init_Retinal_JPCB2000_Model
 
  CONTAINS
!> @brief Function which makes the initialization of the Retinal_JPCB2000 parameters.
!!
!! @param QModel             TYPE(Retinal_JPCB2000_Model_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(EmptyModel_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_Retinal_JPCB2000_Model(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (Retinal_JPCB2000_Model_t)              :: QModel ! RESULT

    TYPE(EmptyModel_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    real (kind=Rkind) :: auTOeV  = 27.211384_Rkind


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_Retinal_JPCB2000_Model'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_EmptyModel(QModel%EmptyModel_t,QModel_in)

    QModel%nsurf    = 2
    QModel%ndim     = 2
    QModel%pot_name = 'retinal_jpcb2000'

    IF (.NOT. QModel%PubliUnit) THEN
      QModel%m      = QModel%m      * auTOeV   ! 1/m has to be divided by auTOeV
      QModel%E1     = QModel%E1     / auTOeV
      QModel%W0     = QModel%W0     / auTOeV
      QModel%W1     = QModel%W1     / auTOeV
      QModel%w      = QModel%w      / auTOeV
      QModel%kappa  = QModel%kappa  / auTOeV
      QModel%lambda = QModel%lambda / auTOeV
    END IF


    IF (debug) write(out_unitp,*) 'init Q0 of Retinal_JPCB2000'
    QModel%Q0 = [ZERO,ZERO]

    IF (debug) write(out_unitp,*) 'init d0GGdef of Retinal_JPCB2000'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%m
    QModel%d0GGdef(2,2) = QModel%w


    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_Retinal_JPCB2000_Model
!> @brief Subroutine wich prints the current Retinal_JPCB2000_Model parameters.
!!
!! @param QModel            CLASS(Retinal_JPCB2000_Model_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_Retinal_JPCB2000_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(Retinal_JPCB2000_Model_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*)
    write(nio,*) 'Retinal_JPCB2000 current parameters'
    write(nio,*) '   Parameters of the model from:'
    write(nio,*) '  S. Hahn, G. Stock / Chemical Physics 259 (2000) 297-312'
    write(nio,*) '     doi: 10.1016/S0301-0104(00)00201-9'
    write(nio,*)
    IF (QModel%PubliUnit) THEN
      write(nio,*) ' Unit in eV (as in the publication)'
    ELSE
      write(nio,*) ' Unit in au (atomic units)'
    END IF
    write(nio,*)
    write(nio,*) '   1/m    = ',ONE/QModel%m
    write(nio,*) '   m      = ',QModel%m
    write(nio,*) '   E1     = ',QModel%E1
    write(nio,*) '   W0     = ',QModel%W0
    write(nio,*) '   W1     = ',QModel%W1
    write(nio,*) '   w      = ',QModel%w
    write(nio,*) '   kappa  = ',QModel%kappa
    write(nio,*) '   lambda = ',QModel%lambda
    write(nio,*)
    write(nio,*) 'end Retinal_JPCB2000 default parameters'
    write(nio,*)

  END SUBROUTINE Write_Retinal_JPCB2000_Model
!> @brief Subroutine wich prints the default Retinal_JPCB2000_Model parameters.
!!
!! @param QModel            CLASS(Retinal_JPCB2000_Model_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_Retinal_JPCB2000_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(Retinal_JPCB2000_Model_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*)
    write(nio,*) 'Retinal_JPCB2000 default parameters'
    write(nio,*) '   Parameters (in eV) of the model from:'
    write(nio,*) '  S. Hahn, G. Stock / Chemical Physics 259 (2000) 297-312'
    write(nio,*) '     doi: 10.1016/S0301-0104(00)00201-9'
    write(nio,*)
    write(nio,*) '   1/m    = 4.84e-4'
    write(nio,*) '   E1     = 2.48'
    write(nio,*) '   W0     = 3.6'
    write(nio,*) '   W1     = 1.09'
    write(nio,*) '   w      = 0.19'
    write(nio,*) '   kappa  = 0.1'
    write(nio,*) '   lambda = 0.19'
    write(nio,*)
    write(nio,*) 'end Retinal_JPCB2000 default parameters'
    write(nio,*)

  END SUBROUTINE Write0_Retinal_JPCB2000_Model

!> @brief Subroutine wich calculates the Retinal_JPCB2000 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(Retinal_JPCB2000_Model_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_Retinal_JPCB2000_Pot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(Retinal_JPCB2000_Model_t), intent(in)    :: QModel
    TYPE (dnS_t),                    intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                    intent(in)    :: dnQ(:)
    integer,                         intent(in)    :: nderiv

    Mat_OF_PotDia(1,1) =             HALF*QModel%W0*(ONE-cos(dnQ(1))) + &
                           HALF*QModel%W*dnQ(2)**2

    Mat_OF_PotDia(2,2) = QModel%E1 - HALF*QModel%W1*(ONE-cos(dnQ(1))) + &
                           HALF*QModel%W*dnQ(2)**2 + QModel%kappa*dnQ(2)

    Mat_OF_PotDia(1,2) = QModel%lambda*dnQ(2)
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


  END SUBROUTINE eval_Retinal_JPCB2000_Pot

END MODULE mod_Retinal_JPCB2000_Model
