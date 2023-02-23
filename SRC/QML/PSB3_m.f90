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
!> @brief Module which makes the initialization, calculation of the PSB3 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_PSB3_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the PSB3 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_PSB3_t
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
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_PSB3
    PROCEDURE :: Write_QModel    => Write_QML_PSB3
  END TYPE QML_PSB3_t

  PUBLIC :: QML_PSB3_t,Init_QML_PSB3

  CONTAINS
!> @brief Function which makes the initialization of the PSB3 parameters.
!!
!! @param QModel             TYPE(QML_PSB3_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_PSB3(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_PSB3_t)                          :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_PSB3'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 2
    QModel%ndim     = 3
    QModel%pot_name = 'psb3'

    IF (QModel%option < 1 .OR. QModel%option > 3) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1,3) ! published model (2020)

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
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' This option is not possible. option:',QModel%option
        write(out_unit,*) ' Its value MUST be 1, 3 or 2'
        STOP
    END SELECT


    IF (debug) write(out_unit,*) 'init Q0 of PSB3'
    QModel%Q0 = [0.172459_Rkind, -3.14_Rkind, ZERO]
    IF (debug) write(out_unit,*) 'init d0GGdef of PSB3'
    QModel%d0GGdef      = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,:) = [0.00007981_Rkind, ZERO,             ZERO            ]
    QModel%d0GGdef(2,:) = [ZERO,             0.00002599_Rkind, 0.00004025_Rkind]
    QModel%d0GGdef(3,:) = [ZERO,             0.00004025_Rkind, 0.00040375_Rkind]

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_PSB3
!> @brief Subroutine wich prints the QML_PSB3 parameters.
!!
!! @param QModel            CLASS(QML_PSB3_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_PSB3(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_PSB3_t),   intent(in) :: QModel
    integer,               intent(in) :: nio

    write(nio,*) 'PSB3 default parameters'
    write(nio,*)
    write(nio,*) ' Warning the parameters are given as in the publication.'
    write(nio,*) '  Therefore, the BLA(=Q(1)) is in Angstrom and the energy is in kcal.mol^-1.'
    write(nio,*)
    write(nio,*) 'end PSB3 default parameters'
    write(nio,*)
    write(nio,*) 'PSB3 current parameters'
    write(nio,*)
    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  adiabatic:      ',QModel%adiabatic
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)

    SELECT CASE (QModel%option)

    CASE (1,3)

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

    IF (QModel%option == 3) write(nio,*) 'Q(3) = Q(3)+2.Pi'

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
        write(out_unit,*) ' ERROR in Write_QML_PSB3 '
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1 or 2 '
        STOP
    END SELECT

    write(nio,*)
    write(nio,*) 'end PSB3 current parameters'

  END SUBROUTINE Write_QML_PSB3
!> @brief Subroutine wich calculates the PSB3 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_PSB3_t):  derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_PSB3(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_PSB3_t),    intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    SELECT CASE (QModel%option)

    CASE (1,3)
      CALL EvalPot1_QML_PSB3(Mat_OF_PotDia,dnQ,QModel,nderiv)

    CASE (2)
      CALL EvalPot2_QML_PSB3(Mat_OF_PotDia,dnQ,QModel,nderiv)

    CASE Default
        write(out_unit,*) ' ERROR in EvalPot_QML_PSB3'
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1,3 or 2 '

        STOP
    END SELECT

  END SUBROUTINE EvalPot_QML_PSB3

!> @brief Subroutine wich calculates the PSB3 potential with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param PSB3Pot          TYPE(PSB3Pot_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot1_QML_PSB3(Mat_OF_PotDia,dnQ,PSB3Pot,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)    :: dnQ(:) ! BLA, Tors, HOOP
    TYPE(QML_PSB3_t),    intent(in)    :: PSB3Pot
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
    IF (PSB3Pot%option == 3) THEN
      HOOP     = dnQ(3)+TWO*Pi
    ELSE
      HOOP     = dnQ(3)
    END IF
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
   CALL dealloc_dnS(BLA)
   CALL dealloc_dnS(Tors)
   CALL dealloc_dnS(HOOP)

   CALL dealloc_dnS(Overlap)
   CALL dealloc_dnS(MorseBLAP)
   CALL dealloc_dnS(Morsemin)
   CALL dealloc_dnS(Hdir2D)
   CALL dealloc_dnS(Hct2D)

   CALL dealloc_dnS(dnPot)

  END SUBROUTINE EvalPot1_QML_PSB3

!> @brief Subroutine wich calculates the PSB3 potential (Published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param PSB3Pot          TYPE(PSB3Pot_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot2_QML_PSB3(Mat_OF_PotDia,dnQ,PSB3Pot,nderiv) !Second PSB3's potential
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (QML_PSB3_t), intent(in)     :: PSB3Pot
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
   CALL dealloc_dnS(BLA)
   CALL dealloc_dnS(Tors)
   CALL dealloc_dnS(HOOP)

   CALL dealloc_dnS(Overlap)
   CALL dealloc_dnS(Hdir2D)
   CALL dealloc_dnS(Hct2D)

   CALL dealloc_dnS(dnPot)

  END SUBROUTINE EvalPot2_QML_PSB3

END MODULE QML_PSB3_m
