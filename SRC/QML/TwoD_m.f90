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

!> @brief Module which makes the initialization, calculation of the TwoD potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_TwoD_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the TwoD parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_TwoD_t
   PRIVATE

   real(kind=Rkind)     :: KX    = 0.02_Rkind
   real(kind=Rkind)     :: KY    = 0.1_Rkind
   real(kind=Rkind)     :: DELTA = 0.01_Rkind
   real(kind=Rkind)     :: X1    = 6._Rkind
   real(kind=Rkind)     :: X2    = 2._Rkind
   real(kind=Rkind)     :: X3    = 31._Rkind/8.0_Rkind
   real(kind=Rkind)     :: GAMMA = 0.01_Rkind
   real(kind=Rkind)     :: ALPHA = 3._Rkind
   real(kind=Rkind)     :: BETA  = 1.5_Rkind

   real (kind=Rkind)    :: muX  = 20000._Rkind
   real (kind=Rkind)    :: muY  = 6667._Rkind

   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_TwoD
    PROCEDURE :: Write_QModel    => Write_QML_TwoD
    PROCEDURE :: Write0_QModel   => Write0_QML_TwoD
  END TYPE QML_TwoD_t

  PUBLIC :: QML_TwoD_t,Init_QML_TwoD

  CONTAINS
!> @brief Subroutine which makes the initialization of the TwoD parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 12/07/2019
!!
!! @param QModel             TYPE(QML_TwoD_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_TwoD(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_TwoD_t)                          :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_TwoD'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 2
    QModel%ndim     = 2
    QModel%pot_name = 'twod'


    IF (debug) write(out_unitp,*) 'init Q0 of TwoD'
    SELECT CASE (QModel%option)
    CASE (1) ! minimum of V(1,1)
      QModel%Q0 = [QModel%X1,ZERO]
    CASE (2)  ! minimum of V(2,2)
      QModel%Q0 = [QModel%X2,ZERO]
    CASE Default ! minimum of V(1,1)
      QModel%Q0 = [QModel%X1,ZERO]
    END SELECT

    IF (debug) write(out_unitp,*) 'init d0GGdef of TwoD'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%muX
    QModel%d0GGdef(2,2) = ONE/QModel%muY

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_TwoD
!> @brief Subroutine wich prints the current QML_TwoD parameters.
!!
!! @param QModel            CLASS(QML_TwoD_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_TwoD(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_TwoD_t),  intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'TwoD current parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) '---- WARNNING ---------------------------'
    write(nio,*) 'The parameters are different from the published ones: '
    write(nio,*) ' A. Ferretti, G. Granucci, A. Lami, M. Persico, G. Villani, ...'
    write(nio,*) '  .... J. Chem. Phys. 104, 5517 (1996); https://doi.org/10.1063/1.471791'
    write(nio,*) 'with the X=Q(1), Y=Q(2) in bohr.'
    write(nio,*) '     and the energy in Hartree.'


    write(nio,*) 'PubliUnit: ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) 'Current parameters:'
    write(nio,*) '-----------------------------------------'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'KX    =',QModel%KX
    write(nio,*) 'KY    =',QModel%KY
    write(nio,*) 'DELTA =',QModel%DELTA
    write(nio,*) 'X1    =',QModel%X1
    write(nio,*) 'X2    =',QModel%X2
    write(nio,*) 'X3    =',QModel%X3
    write(nio,*) 'GAMMA =',QModel%GAMMA
    write(nio,*) 'ALPHA =',QModel%ALPHA
    write(nio,*) 'BETA  =',QModel%BETA
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'MuX =',QModel%MuX
    write(nio,*) 'MuY =',QModel%MuY
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Q0  =',QModel%Q0
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'end TwoD parameters'

  END SUBROUTINE Write_QML_TwoD
!> @brief Subroutine wich prints the default QML_TwoD parameters.
!!
!! @param QModel            CLASS(QML_TwoD_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_TwoD(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_TwoD_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'TwoD default parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) '---- WARNNING ---------------------------'
    write(nio,*) 'The parameters are different from the published ones: '
    write(nio,*) ' A. Ferretti, G. Granucci, A. Lami, M. Persico, G. Villani, ...'
    write(nio,*) '  .... J. Chem. Phys. 104, 5517 (1996); https://doi.org/10.1063/1.471791'
    write(nio,*) 'with the X=Q(1), Y=Q(2) in bohr.'
    write(nio,*) '     and the energy in Hartree.'
    write(nio,*)
    write(nio,*) 'Diabatic Potential Values (in Hartree) at: R=3.875 bohr and theta=0.5 bohr'
    write(nio,*) '1        0.05765625  0.00343645'
    write(nio,*) '2        0.00343645  0.05765625'


    write(nio,*) 'PubliUnit: ',QModel%PubliUnit

    write(nio,*)
    write(nio,*) 'Default parameters:'
    write(nio,*) 'KX    = 0.02'
    write(nio,*) 'KY    = 0.1'
    write(nio,*) 'DELTA =0.01'
    write(nio,*) 'X1    = 6.'
    write(nio,*) 'X2    = 2.'
    write(nio,*) 'X3    = 31./8.'
    write(nio,*) 'GAMMA = 0.01'
    write(nio,*) 'ALPHA = 3.'
    write(nio,*) 'BETA  = 1.5'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'MuX = 20000.'
    write(nio,*) 'MuY = 6667.'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Q0  = [6.,0.]'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'end TwoD default parameters'

  END SUBROUTINE Write0_QML_TwoD

!> @brief Subroutine wich calculates the TwoD potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_TwoD_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_TwoD(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS(QML_TwoD_t),  intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv


   !Hel(1,1,x,y)=0.5d0*KX*(R(1,x)-X1)**2 + 0.5d0*KY*(R(2,y))**2
   Mat_OF_PotDia(1,1) = HALF * QModel%KX*(dnQ(1)-QModel%X1)**2 + HALF * QModel%KY*(dnQ(2))**2

   !Hel(2,2,x,y)=0.5d0*KX*(R(1,x)-X2)**2 + 0.5d0*KY*(R(2,y))**2 +DELTA
   Mat_OF_PotDia(2,2) = HALF * QModel%KX*(dnQ(1)-QModel%X2)**2 + HALF * QModel%KY*(dnQ(2))**2 + QModel%DELTA

   !Hel(1,2,x,y)=GAMMA*R(2,y)*dexp(-ALPHA*(R(1,x)-X3)**2)*dexp(-BETA*(R(2,y))**2)
   Mat_OF_PotDia(1,2) = QModel%GAMMA * dnQ(2) * exp(-QModel%ALPHA*(dnQ(1)-QModel%X3)**2) * exp(-QModel%BETA*(dnQ(2))**2)
   Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)

  END SUBROUTINE EvalPot_QML_TwoD

END MODULE QML_TwoD_m
