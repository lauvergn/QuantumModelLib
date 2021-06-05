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
!
!> @brief Module which makes the initialization, calculation of the PSB3 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_PSB3_Model
  USE QML_NumParameters_m
  USE mod_EmptyModel
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the PSB3 parameters are set-up.
  TYPE, EXTENDS (EmptyModel_t) ::  PSB3_Model_t
    PRIVATE

    real (kind=Rkind) :: blamin   = 0.0912615_Rkind
    real (kind=Rkind) :: blaTSdir = 0.025079167_Rkind
    real (kind=Rkind) :: deepth   = 2000_Rkind

    real (kind=Rkind) :: kf1 = 3733.5_Rkind
    real (kind=Rkind) :: d2  = 54.633_Rkind
    real (kind=Rkind) :: d3  = 3.8008_Rkind
    real (kind=Rkind) :: kf4 = 1097.19_Rkind

    real (kind=Rkind) :: c1  = 437.068_Rkind
    real (kind=Rkind) :: c2  = 16.7317_Rkind
    real (kind=Rkind) :: c3  = 7.35468_Rkind
    real (kind=Rkind) :: c4  = 88.517_Rkind
    real (kind=Rkind) :: c5  = 5.95221_Rkind

    real (kind=Rkind) :: h1  = 155.749_Rkind

    real (kind=Rkind) :: k1  = 24.0358_Rkind

    !---Additional parameter for the yet unpublished model-------!
    real (kind=Rkind) :: d1  = 2571.3_Rkind
    real (kind=Rkind) :: d4  = 594.97_Rkind

    real (kind=Rkind) :: hd1  = 126.22_Rkind
    real (kind=Rkind) :: hd2  = 75.003_Rkind
    real (kind=Rkind) :: hc1  = 86.464_Rkind
    real (kind=Rkind) :: hc2  = 43.461_Rkind

    real (kind=Rkind) :: k3  = 1.1347_Rkind
    !-----------------------------------------------------!

    ! Warning the parameters are given as in the publication.
    !   Therefore, the BLA(=Q(1)) is in Angstrom and the energy is in kcal.mol^-1.
   CONTAINS
    PROCEDURE :: Eval_QModel_Pot => eval_PSB3_Pot
    PROCEDURE :: Write_QModel    => Write_PSB3_Model
    PROCEDURE :: Write0_QModel   => Write0_PSB3_Model
  END TYPE PSB3_Model_t
 
  PUBLIC :: PSB3_Model_t,Init_PSB3_Model
 
  CONTAINS
!> @brief Function which makes the initialization of the PSB3 parameters.
!!
!! @param QModel             TYPE(PSB3_Model_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(EmptyModel_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_PSB3_Model(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (PSB3_Model_t)                          :: QModel ! RESULT

    TYPE(EmptyModel_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_PSB3_Model'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_EmptyModel(QModel%EmptyModel_t,QModel_in)

    QModel%nsurf    = 2
    QModel%ndim     = 3
    QModel%pot_name = 'psb3'

    IF (QModel%option < 1 .OR. QModel%option > 2) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1) ! unpublished model

      QModel%blamin   = 0.0912615_Rkind
      QModel%blaTSdir = 0.025079167_Rkind
      QModel%deepth   = 2000_Rkind

      QModel%kf1      = 3733.5_Rkind
      QModel%d2       = 54.633_Rkind
      QModel%d3       = 3.8008_Rkind
      QModel%kf4      = 1097.19_Rkind

      QModel%c1       = 437.068_Rkind
      QModel%c2       = 16.7317_Rkind
      QModel%c3       = 7.35468_Rkind
      QModel%c4       = 88.517_Rkind
      QModel%c5       = 5.95221_Rkind

      QModel%h1       = 155.749_Rkind
      QModel%k1       = 24.0358_Rkind

    CASE (2) ! Published model

      QModel%blamin   = 0.0912615_Rkind
      QModel%d1       = 2571.3_Rkind
      QModel%d2       = 54.633_Rkind
      QModel%d3       = 3.8008_Rkind
      QModel%d4       = 594.97_Rkind

      QModel%c1       = 437.068_Rkind
      QModel%c2       = 16.7317_Rkind
      QModel%c3       = 7.35468_Rkind
      QModel%c4       = 88.517_Rkind
      QModel%c5       = 5.95221_Rkind

      QModel%hd1      = 126.22_Rkind
      QModel%hd2      = 75.003_Rkind
      QModel%hc1      = 86.464_Rkind
      QModel%hc2      = 43.461_Rkind

      QModel%k1       = 13.699_Rkind
      QModel%k3       = 1.1347_Rkind

    CASE Default
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' This option is not possible. option:',QModel%option
        write(out_unitp,*) ' Its value MUST be 1 or 2'
        STOP
    END SELECT


    IF (debug) write(out_unitp,*) 'init Q0 of PSB3'
    QModel%Q0 = [0.172459_Rkind, -3.14_Rkind, ZERO]
    IF (debug) write(out_unitp,*) 'init d0GGdef of PSB3'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
    QModel%d0GGdef(1,:) = [0.00007981_Rkind, ZERO,             ZERO            ]
    QModel%d0GGdef(2,:) = [ZERO,             0.00002599_Rkind, 0.00004025_Rkind]
    QModel%d0GGdef(3,:) = [ZERO,             0.00004025_Rkind, 0.00040375_Rkind]

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_PSB3_Model
!> @brief Subroutine wich prints the PSB3_Model parameters.
!!
!! @param QModel            CLASS(PSB3_Model_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_PSB3_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(PSB3_Model_t),   intent(in) :: QModel
    integer,               intent(in) :: nio

    write(nio,*) 'PSB3 current parameters'
    write(nio,*)
    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  adiabatic:      ',QModel%adiabatic
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)

    SELECT CASE (QModel%option)

    CASE (1)

    write(nio,*) '  ... for first Model:'
    write(nio,*) '  blamin   :      ',QModel%blamin
    write(nio,*) '  blaTSdir :      ',QModel%blaTSdir
    write(nio,*) '  deepth   :      ',QModel%deepth
    write(nio,*) '  kf1      :      ',QModel%kf1
    write(nio,*) '  d2       :      ',QModel%d2
    write(nio,*) '  d3       :      ',QModel%d3
    write(nio,*) '  kf4      :      ',QModel%kf4
    write(nio,*) '  c1       :      ',QModel%c1
    write(nio,*) '  c2       :      ',QModel%c2
    write(nio,*) '  c3       :      ',QModel%c3
    write(nio,*) '  c4       :      ',QModel%c4
    write(nio,*) '  c5       :      ',QModel%c5
    write(nio,*) '  h1       :      ',QModel%h1
    write(nio,*) '  k1       :      ',QModel%k1

    CASE (2)

    write(nio,*) '  ... for second Model:'
    write(nio,*) '  blamin   :      ',QModel%blamin
    write(nio,*) '  d1       :      ',QModel%d1
    write(nio,*) '  d2       :      ',QModel%d2
    write(nio,*) '  d3       :      ',QModel%d3
    write(nio,*) '  d4       :      ',QModel%d4
    write(nio,*) '  c1       :      ',QModel%c1
    write(nio,*) '  c2       :      ',QModel%c2
    write(nio,*) '  c3       :      ',QModel%c3
    write(nio,*) '  c4       :      ',QModel%c4
    write(nio,*) '  c5       :      ',QModel%c5
    write(nio,*) '  hd1      :      ',QModel%hd1
    write(nio,*) '  hd2      :      ',QModel%hd2
    write(nio,*) '  hc1      :      ',QModel%hc1
    write(nio,*) '  hc2      :      ',QModel%hc2
    write(nio,*) '  k1       :      ',QModel%k1
    write(nio,*) '  k3       :      ',QModel%k3

    CASE Default
        write(out_unitp,*) ' ERROR in Write_PSB3_Model '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1 or 2 '
        STOP
    END SELECT

    write(nio,*)
    write(nio,*) 'end PSB3 current parameters'

  END SUBROUTINE Write_PSB3_Model
  SUBROUTINE Write0_PSB3_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(PSB3_Model_t),   intent(in) :: QModel
    integer,               intent(in) :: nio

    write(nio,*) 'PSB3 default parameters'
    write(nio,*)
    write(nio,*) ' Warning the parameters are given as in the publication.'
    write(nio,*) '  Therefore, the BLA(=Q(1)) is in Angstrom and the energy is in kcal.mol^-1.'
    write(nio,*)
    write(nio,*) 'end PSB3 default parameters'


  END SUBROUTINE Write0_PSB3_Model

!> @brief Subroutine wich calculates the PSB3 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(PSB3_Model_t):  derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_PSB3_Pot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE QML_dnS_m
  IMPLICIT NONE

    CLASS(PSB3_Model_t),  intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    SELECT CASE (QModel%option)

    CASE (1)
      CALL eval_PSB3Pot1(Mat_OF_PotDia,dnQ,QModel,nderiv)

    CASE (2)
      CALL eval_PSB3Pot2(Mat_OF_PotDia,dnQ,QModel,nderiv)

    CASE Default
        write(out_unitp,*) ' ERROR in eval_PSB3_Pot'
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1 or 2 '

        STOP
    END SELECT

  END SUBROUTINE eval_PSB3_Pot

!> @brief Subroutine wich calculates the PSB3 potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param PSB3Pot          TYPE(PSB3Pot_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_PSB3Pot1(Mat_OF_PotDia,dnQ,PSB3Pot,nderiv)
    !Unpublished model potential (yet)
    USE QML_dnS_m
    IMPLICIT NONE

    TYPE (dnS_t),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)    :: dnQ(:) ! BLA, Tors, HOOP
    TYPE(PSB3_Model_t),  intent(in)    :: PSB3Pot
    integer,             intent(in)    :: nderiv


    !local variables (derived type). They have to be deallocated
    TYPE (dnS_t)        :: dnPot,BLA,Tors,HOOP
    TYPE (dnS_t)        :: MorseBLAP,Morsemin,Overlap
    TYPE (dnS_t)        :: Hdir2D,Hct2D
    real (kind=Rkind)   :: d1,d4

    real (kind=Rkind),   parameter     :: EnergyConv = 627.509_Rkind
    real (kind=Rkind),   parameter     :: LenghtConv = 0.52917721067121_Rkind



    BLA      = dnQ(1)
    Tors     = dnQ(2)
    HOOP     = dnQ(3)

    !If PubliUnit is False conversion must be done, the potential is expressed in Angstrom
    !and it requires the proper conversion into Bhor
    IF(.NOT. PSB3Pot%PubliUnit) THEN
       BLA = BLA * LenghtConv
    END IF

!-----------------------------------------------------------------------!
    d1 = SQRT(PSB3Pot%kf1/(TWO * PSB3Pot%deepth))
    d4 = SQRT(PSB3Pot%kf4/(TWO * PSB3Pot%deepth))

    MorseBLAP = PSB3Pot%deepth * (-ONE + EXP(-d1 * (BLA - PSB3Pot%blaTSdir)))**2
    Morsemin  = PSB3Pot%deepth * (-ONE + EXP(-d4 * (BLA - PSB3Pot%blamin  )))**2

    Overlap   = Tors - HOOP/TWO

    Hdir2D = sin(Overlap)**2 * (MorseBLAP + PSB3Pot%d2) + PSB3Pot%d3 * cos(Overlap / TWO)**2 + Morsemin * Cos(Overlap)**2
    Hct2D  = (ONE + PSB3Pot%c5 * sin(Overlap)**2) * (PSB3Pot%c1 * BLA**2 + &
             PSB3Pot%c2 * BLA + PSB3Pot%c3) + PSB3Pot%c4 * cos(Overlap)**2
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
! V11 matrix element

   dnPot =  Hdir2D + PSB3Pot%h1 * sin(HOOP/FOUR)**2

   IF(.NOT. PSB3Pot%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(1,1) = dnPot
!-----------------------------------------------------------------------!
! V22 matrix element

   dnPot =  Hct2D + PSB3Pot%h1 * sin(HOOP/FOUR)**2

   IF(.NOT. PSB3Pot%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(2,2) = dnPot

!-----------------------------------------------------------------------!
!V12 = V21

   dnPot = PSB3Pot%k1 * sin((Overlap) * TWO)

   IF(.NOT. PSB3Pot%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(1,2) = dnPot
   Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


!-----------------------------------------------------------------------!
   CALL QML_dealloc_dnS(BLA)
   CALL QML_dealloc_dnS(Tors)
   CALL QML_dealloc_dnS(HOOP)

   CALL QML_dealloc_dnS(Overlap)
   CALL QML_dealloc_dnS(MorseBLAP)
   CALL QML_dealloc_dnS(Morsemin)
   CALL QML_dealloc_dnS(Hdir2D)
   CALL QML_dealloc_dnS(Hct2D)

   CALL QML_dealloc_dnS(dnPot)

  END SUBROUTINE eval_PSB3Pot1

!> @brief Subroutine wich calculates the PSB3 potential (Published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param PSB3Pot          TYPE(PSB3Pot_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_PSB3Pot2(Mat_OF_PotDia,dnQ,PSB3Pot,nderiv) !Second PSB3's potential
  ! Published potential
  USE QML_dnS_m
  IMPLICIT NONE

    TYPE (PSB3_Model_t), intent(in)     :: PSB3Pot
    TYPE (dnS_t),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)     :: dnQ(:) ! BLA, Tors, HOOP
    integer,             intent(in)     :: nderiv

    !local variables (derived type). They have to be deallocated
    TYPE (dnS_t)       :: dnPot,BLA,Tors,HOOP
    TYPE (dnS_t)       :: Overlap
    TYPE (dnS_t)       :: Hdir2D,Hct2D
    real (kind=Rkind)  :: d1,d4

    real (kind=Rkind),   parameter      :: EnergyConv = 627.509_Rkind
    real (kind=Rkind),   parameter      :: LenghtConv = 0.52917721067121_Rkind



    BLA      = dnQ(1)
    Tors     = dnQ(2)
    HOOP     = dnQ(3)


    !If PubliUnit is False conversion must be done, the potential is expressed in Angstrom
    !and it requires the proper conversion into Bhor
    IF(.NOT. PSB3Pot%PubliUnit) THEN
      BLA = BLA * LenghtConv    ! to set up the correct derivatives with respect to (R*LenghtConv)
    END IF

!-----------------------------------------------------------------------!
    Overlap   = Tors - HOOP/TWO

    Hdir2D = SIN(Tors)**2 * (PSB3Pot%d1 * BLA**2 + PSB3Pot%d2) + &
             PSB3Pot%d3 * COS(Tors / TWO)**2 + PSB3Pot%d4 * (BLA - PSB3Pot%blamin)**2
    Hct2D  = (ONE + PSB3Pot%c5 * SIN(Tors)**2) *                   &
             (PSB3Pot%c1 * BLA**2 + PSB3Pot%c2 * BLA + PSB3Pot%c3) + PSB3Pot%c4 * COS(Tors)**2
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
! V11 matrix element

   dnPot =  Hdir2D + PSB3Pot%hd1 * sin(HOOP/FOUR)**2 - PSB3Pot%hd2 * Sin(HOOP/FOUR) * sin(Tors * TWO)

   IF(.NOT. PSB3Pot%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(1,1) = dnPot

!-----------------------------------------------------------------------!
! V22 matrix element

   dnPot =  Hct2D  + PSB3Pot%hc1 * Sin(HOOP/FOUR)**2 + PSB3Pot%hc2 * Sin(HOOP/FOUR) * sin(Tors * TWO)

   IF(.NOT. PSB3Pot%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(2,2) = dnPot

!-----------------------------------------------------------------------!
! V12 = V21

   dnPot = (ONE + PSB3Pot%k3 * SIN(Tors)**2) *  PSB3Pot%k1 * SIN((Overlap) * TWO)

   IF(.NOT. PSB3Pot%PubliUnit) dnPot = dnPot/EnergyConv

   Mat_OF_PotDia(1,2) = dnPot
   Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


!-----------------------------------------------------------------------!
   CALL QML_dealloc_dnS(BLA)
   CALL QML_dealloc_dnS(Tors)
   CALL QML_dealloc_dnS(HOOP)

   CALL QML_dealloc_dnS(Overlap)
   CALL QML_dealloc_dnS(Hdir2D)
   CALL QML_dealloc_dnS(Hct2D)

   CALL QML_dealloc_dnS(dnPot)

  END SUBROUTINE eval_PSB3Pot2

END MODULE mod_PSB3_Model
