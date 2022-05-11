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
MODULE QML_Retinal_JPCB2000_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Retinal_JPCB2000 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_Retinal_JPCB2000_t

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

   real(kind=Rkind) :: wi(25)       = real([   0.0,   0.0,  792.8,  842.8,  866.2,  &
                                             882.4, 970.3,   976.0,  997.0, 1017.1, &
                                            1089.6, 1189.0, 1214.7, 1238.1, 1267.9, &
                                            1317.0, 1359.0, 1389.0, 1428.4, 1434.9, &
                                            1451.8, 1572.8, 1612.1, 1629.2, 1659.1],kind=Rkind)
   real(kind=Rkind) :: ciwi_inv(25) = real([ 0.0,   0.0,    0.175,  0.2,    0.175, &
                                             0.225, 0.55,   0.3,    0.33,   0.45,  &
                                             0.125,  0.175, 0.44,   0.5,    0.475, &
                                             0.238,  0.25,  0.25,   0.25,   0.225, &
                                             0.225,  0.25,  0.225,  0.125,  0.225],kind=Rkind)
   real(kind=Rkind) :: ci(25)       = ZERO ! it will be initialized after

   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_Retinal_JPCB2000
    PROCEDURE :: Write_QModel    => Write_QML_Retinal_JPCB2000
    PROCEDURE :: Write0_QModel   => Write0_QML_Retinal_JPCB2000
  END TYPE QML_Retinal_JPCB2000_t

  PUBLIC :: QML_Retinal_JPCB2000_t,Init_QML_Retinal_JPCB2000

  CONTAINS
!> @brief Function which makes the initialization of the Retinal_JPCB2000 parameters.
!!
!! @param QModel             TYPE(QML_Retinal_JPCB2000_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_Retinal_JPCB2000(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_Retinal_JPCB2000_t)              :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    real (kind=Rkind), parameter :: auTOeV      = 27.211384_Rkind
    real (kind=Rkind), parameter :: auTOcm_inv  = 219474.631443_Rkind
    integer :: i


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Retinal_JPCB2000'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 2

    IF (QModel%ndim == 0) THEN
      ! it means that ndim was not present in the CALL Init_Model().
      ! => ndim is set to the default value (2)
      QModel%ndim  = 2
      QModel%pot_name = 'retinal_jpcb2000'
    END IF
    IF (QModel%ndim  > 2) QModel%pot_name = 'retinal_cp2000'

    !The value of QModel%ndim must be in [2:25]
    IF (QModel%ndim < 2 .OR. QModel%ndim > 25) THEN
       write(out_unitp,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unitp)
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' ndim range MUST be [2:25]. ndim: ',QModel%ndim
       STOP 'ERROR in Init_QML_Retinal_JPCB2000: ndim range MUST be [2:25]'
    END IF
    IF (debug) write(out_unitp,*) 'ndim',QModel%ndim



    IF (.NOT. QModel%PubliUnit) THEN
      QModel%m      = QModel%m      * auTOeV   ! 1/m has to be divided by auTOeV
      QModel%E1     = QModel%E1     / auTOeV
      QModel%W0     = QModel%W0     / auTOeV
      QModel%W1     = QModel%W1     / auTOeV
      QModel%w      = QModel%w      / auTOeV
      QModel%kappa  = QModel%kappa  / auTOeV
      QModel%lambda = QModel%lambda / auTOeV

      QModel%wi     = QModel%wi     / auTOcm_inv

    ELSE ! if PubliUnit the wi must be convert in eV
      QModel%wi     = QModel%wi     / auTOcm_inv*auTOeV
    END IF

    QModel%ci       = QModel%ciwi_inv * QModel%wi

    IF (debug) write(out_unitp,*) 'init Q0 of Retinal_JPCB2000'
    allocate(QModel%Q0(QModel%ndim))
    QModel%Q0 = ZERO

    IF (debug) write(out_unitp,*) 'init d0GGdef of Retinal_JPCB2000'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    QModel%d0GGdef(1,1) = ONE/QModel%m
    QModel%d0GGdef(2,2) = QModel%w

    DO i=3,QModel%ndim
      QModel%d0GGdef(i,i) = QModel%wi(i)
    END DO

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_Retinal_JPCB2000
!> @brief Subroutine wich prints the current QML_Retinal_JPCB2000 parameters.
!!
!! @param QModel            CLASS(QML_Retinal_JPCB2000_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_Retinal_JPCB2000(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Retinal_JPCB2000_t),   intent(in) :: QModel
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

  END SUBROUTINE Write_QML_Retinal_JPCB2000
!> @brief Subroutine wich prints the default QML_Retinal_JPCB2000 parameters.
!!
!! @param QModel            CLASS(QML_Retinal_JPCB2000_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_Retinal_JPCB2000(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Retinal_JPCB2000_t),   intent(in) :: QModel
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

  END SUBROUTINE Write0_QML_Retinal_JPCB2000

!> @brief Subroutine wich calculates the Retinal_JPCB2000 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_Retinal_JPCB2000_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Retinal_JPCB2000(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_Retinal_JPCB2000_t),   intent(in)    :: QModel
    TYPE (dnS_t),                    intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                    intent(in)    :: dnQ(:)
    integer,                         intent(in)    :: nderiv

    TYPE (dnS_t) :: VBath
    integer      :: i

    IF (QModel%ndim > 2) THEN
      VBath = HALF*sum( QModel%wi(3:QModel%ndim)*dnQ(3:QModel%ndim)**2 )

      Mat_OF_PotDia(1,1) = HALF*QModel%W0*(ONE-cos(dnQ(1))) +                   &
                           HALF*QModel%W*dnQ(2)**2 + VBath

      Mat_OF_PotDia(2,2) = QModel%E1 - HALF*QModel%W1*(ONE-cos(dnQ(1))) +       &
                           HALF*QModel%W*dnQ(2)**2 + QModel%kappa*dnQ(2) +      &
                        VBath + sum( QModel%ci(3:QModel%ndim)*dnQ(3:QModel%ndim) )

      Mat_OF_PotDia(1,2) = QModel%lambda*dnQ(2)
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)

      CALL dealloc_dnS(VBath)

    ELSE
      Mat_OF_PotDia(1,1) = HALF*QModel%W0*(ONE-cos(dnQ(1))) +                   &
                           HALF*QModel%W*dnQ(2)**2

      Mat_OF_PotDia(2,2) = QModel%E1 - HALF*QModel%W1*(ONE-cos(dnQ(1))) +       &
                           HALF*QModel%W*dnQ(2)**2 + QModel%kappa*dnQ(2)

      Mat_OF_PotDia(1,2) = QModel%lambda*dnQ(2)
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)
    END IF


  END SUBROUTINE EvalPot_QML_Retinal_JPCB2000

END MODULE QML_Retinal_JPCB2000_m
