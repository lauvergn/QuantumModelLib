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
!
!> @brief Module which makes the initialization, calculation of the PSB3 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 30/09/2021
!!
MODULE QML_NO3_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

  TYPE no3minus_pot_t
    !common /no3minus_pot/e0ref,e0rho,p,pst,iref
    real (kind=Rkind) :: e0ref,e0rho
    integer           :: iref
    integer           :: n = 337
    real (kind=Rkind) :: p(337)
    integer           :: npst = 100
    integer           :: pst(2,100)
  END TYPE no3minus_pot_t
  TYPE tmc_t 
    !common /tmc/a,npoly
    integer           :: npoly(2)
    real (kind=Rkind) :: a(21)  !copy of the first elements of p array
  END TYPE tmc_t
  TYPE transformcoordblock_t
    !common /transformcoordblock/ phi_ref,le_ref,beta_ref
    real (kind=Rkind) :: phi_ref,le_ref,beta_ref
  END TYPE transformcoordblock_t

  TYPE acoord_vjt_vee_t
    real (kind=Rkind) :: pa(8), pb(4)

    real (kind=Rkind) :: va(6), vb(6)
    real (kind=Rkind) :: wa(9), wb(9)
    real (kind=Rkind) :: za(9), zb(9)

    real (kind=Rkind) :: vee(28),  wee(49), zee(49)
  END TYPE acoord_vjt_vee_t
    !> @brief Derived type in which the NO3 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_NO3_t
    PRIVATE

    TYPE (no3minus_pot_t)         :: no3minus_pot
    TYPE (tmc_t)                  :: tmc
    TYPE (transformcoordblock_t)  :: transformcoordblock

  CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_NO3
    PROCEDURE :: Write_QModel     => Write_QML_NO3
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_QML_NO3
  END TYPE QML_NO3_t

  PUBLIC :: QML_NO3_t,Init_QML_NO3

  CONTAINS
  ! adapted for QML 06/10/2021
  ! A Viel 2015 03 16
  ! update reference geometry for coordinates transformation
  ! update parameter set 337 parameters
  ! for NO3 E" surface
  ! two coupled surface case
  ! from file tmc-ang-prec-37um2b_extra_b5_20150316.par

!> @brief Function which makes the initialization of the NO3 parameters.
!!
!! @param QModel             TYPE(QML_NO3_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_NO3(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_NO3_t)                             :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    integer :: i,j

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_NO3'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      Write(out_unit,*)'# dated:30/09/2021'
      Write(out_unit,*)'# From cartesian to symmetrized internar coord'
      Write(out_unit,*)'# Scaled umbrella     '
      Write(out_unit,*)'# NO3 E'' JCP2017 case '
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf       = 2
    QModel%pot_name    = 'NO3'
    QModel%no_ana_der  = .TRUE.

    QModel%ndimCart    = 12
    QModel%ndimQ       = 12 ! always in CC


    IF (QModel%option /= 1) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1) ! unpublished model
      IF (QModel%Cart_TO_Q) THEN
        QModel%ndim       = QModel%ndimCart
      ELSE
        QModel%ndim       = QModel%ndimQ
      END IF
      !reference geometry
      QModel%transformcoordblock%phi_ref  = TWO/THREE*pi  ! 120°
      QModel%transformcoordblock%le_ref   = 2.344419_Rkind
      !for completness reno=2.344419d0,reoo=4.0606528d0,rad=0.017453293d0)
      QModel%transformcoordblock%beta_ref = pi*HALF
      Write(out_unit,*)'# ref geometry:'
      Write(out_unit,*)'# phi_ref: ',QModel%transformcoordblock%phi_ref
      Write(out_unit,*)'# le_ref:  ',QModel%transformcoordblock%le_ref
      Write(out_unit,*)'# beta_ref:',QModel%transformcoordblock%beta_ref


      ! all parameters
      QModel%no3minus_pot%iref = 2
      QModel%no3minus_pot%n    = 337
      QModel%no3minus_pot%npst = 30

      ! tmc-ang-prec-37um2b_extra_b5_20150316.par
      ! including more extra points, 16.03.2015
      !                               Fitting RMS:      0.85624381E-03   (   187.9 cm-1)
      !                               Merit function:   0.85624381E-03
      !                               Reduced RMS:      0.30159465E-03   (    66.2 cm-1) for data below  -279.500000
      !                               Sum of point weights:  0.9458
      !                               Number of weighted points:    2404
      QModel%no3minus_pot%p(  1)=  -0.123440260E+00_Rkind
      QModel%no3minus_pot%p(  2)=   0.231012750E+01_Rkind
      QModel%no3minus_pot%p(  3)=  -0.131815690E+01_Rkind
      QModel%no3minus_pot%p(  4)=   0.505372680E+01_Rkind
      QModel%no3minus_pot%p(  5)=  -0.202745100E+02_Rkind
      QModel%no3minus_pot%p(  6)=   0.123400003E+00_Rkind
      QModel%no3minus_pot%p(  7)=   0.556995241E+00_Rkind
      QModel%no3minus_pot%p(  8)=   0.505947070E+01_Rkind
      QModel%no3minus_pot%p(  9)=  -0.840552473E+02_Rkind
      QModel%no3minus_pot%p( 10)=   0.303348770E+01_Rkind
      QModel%no3minus_pot%p( 11)=   0.915482720E+00_Rkind
      QModel%no3minus_pot%p( 12)=   0.597447070E+01_Rkind
      QModel%no3minus_pot%p( 13)=   0.129964090E+01_Rkind
      QModel%no3minus_pot%p( 14)=   0.923982010E+01_Rkind
      QModel%no3minus_pot%p( 15)=   0.123400003E+00_Rkind
      QModel%no3minus_pot%p( 16)=   0.537126740E+01_Rkind
      QModel%no3minus_pot%p( 17)=  -0.528341500E+01_Rkind
      QModel%no3minus_pot%p( 18)=   0.662076740E+02_Rkind
      QModel%no3minus_pot%p( 19)=   0.833679770E+01_Rkind
      QModel%no3minus_pot%p( 20)=   0.286309010E+02_Rkind
      QModel%no3minus_pot%p( 21)=   0.123400003E+00_Rkind
      QModel%no3minus_pot%p( 22)=  -0.304431280E+01_Rkind
      QModel%no3minus_pot%p( 23)=   0.748191810E+01_Rkind
      QModel%no3minus_pot%p( 24)=  -0.448249980E+01_Rkind
      QModel%no3minus_pot%p( 25)=   0.206854010E+02_Rkind
      QModel%no3minus_pot%p( 26)=  -0.172197980E+02_Rkind
      QModel%no3minus_pot%p( 27)=  -0.642527480E+02_Rkind
      QModel%no3minus_pot%p( 28)=   0.565738800E+02_Rkind
      QModel%no3minus_pot%p( 29)=   0.303284550E+02_Rkind
      QModel%no3minus_pot%p( 30)=   0.302292120E+02_Rkind
      QModel%no3minus_pot%p( 31)=   0.792999760E+02_Rkind
      QModel%no3minus_pot%p( 32)=  -0.711077690E+02_Rkind
      QModel%no3minus_pot%p( 33)=   0.575018080E+02_Rkind
      QModel%no3minus_pot%p( 34)=  -0.347775870E+02_Rkind
      QModel%no3minus_pot%p( 35)=   0.930777310E+01_Rkind
      QModel%no3minus_pot%p( 36)=  -0.145823480E+02_Rkind
      QModel%no3minus_pot%p( 37)=   0.520607240E+03_Rkind
      QModel%no3minus_pot%p( 38)=  -0.123662040E+04_Rkind
      QModel%no3minus_pot%p( 39)=   0.125896900E+04_Rkind
      QModel%no3minus_pot%p( 40)=  -0.681959190E+03_Rkind
      QModel%no3minus_pot%p( 41)=   0.181199110E+03_Rkind
      QModel%no3minus_pot%p( 42)=  -0.419560080E+03_Rkind
      QModel%no3minus_pot%p( 43)=   0.105702100E+04_Rkind
      QModel%no3minus_pot%p( 44)=  -0.607165880E+03_Rkind
      QModel%no3minus_pot%p( 45)=  -0.640175490E+03_Rkind
      QModel%no3minus_pot%p( 46)=   0.701414450E+03_Rkind
      QModel%no3minus_pot%p( 47)=   0.237988350E+03_Rkind
      QModel%no3minus_pot%p( 48)=  -0.595329960E+03_Rkind
      QModel%no3minus_pot%p( 49)=   0.202039810E+03_Rkind
      QModel%no3minus_pot%p( 50)=  -0.872841948E+01_Rkind
      QModel%no3minus_pot%p( 51)=  -0.409836740E+00_Rkind
      QModel%no3minus_pot%p( 52)=   0.801472910E+01_Rkind
      QModel%no3minus_pot%p( 53)=   0.526269750E+02_Rkind
      QModel%no3minus_pot%p( 54)=   0.227930150E+02_Rkind
      QModel%no3minus_pot%p( 55)=  -0.349545650E+02_Rkind
      QModel%no3minus_pot%p( 56)=   0.380465090E+02_Rkind
      QModel%no3minus_pot%p( 57)=  -0.744954770E+02_Rkind
      QModel%no3minus_pot%p( 58)=   0.959469230E+01_Rkind
      QModel%no3minus_pot%p( 59)=   0.123400003E+00_Rkind
      QModel%no3minus_pot%p( 60)=  -0.531160440E+01_Rkind
      QModel%no3minus_pot%p( 61)=  -0.872510030E+01_Rkind
      QModel%no3minus_pot%p( 62)=   0.254766890E+02_Rkind
      QModel%no3minus_pot%p( 63)=   0.545690110E+03_Rkind
      QModel%no3minus_pot%p( 64)=  -0.234308810E+03_Rkind
      QModel%no3minus_pot%p( 65)=  -0.676611530E+01_Rkind
      QModel%no3minus_pot%p( 66)=  -0.131862020E+04_Rkind
      QModel%no3minus_pot%p( 67)=   0.385114210E+00_Rkind
      QModel%no3minus_pot%p( 68)=  -0.162871150E+03_Rkind
      QModel%no3minus_pot%p( 69)=   0.123400003E+00_Rkind
      QModel%no3minus_pot%p( 70)=  -0.487582430E-01_Rkind
      QModel%no3minus_pot%p( 71)=  -0.493363570E-01_Rkind
      QModel%no3minus_pot%p( 72)=   0.493413080E+00_Rkind
      QModel%no3minus_pot%p( 73)=  -0.786792930E+00_Rkind
      QModel%no3minus_pot%p( 74)=   0.148580440E+01_Rkind
      QModel%no3minus_pot%p( 75)=   0.275701800E+01_Rkind
      QModel%no3minus_pot%p( 76)=  -0.235189530E+01_Rkind
      QModel%no3minus_pot%p( 77)=   0.147278820E+02_Rkind
      QModel%no3minus_pot%p( 78)=   0.123400003E+00_Rkind
      QModel%no3minus_pot%p( 79)=  -0.238514350E+00_Rkind
      QModel%no3minus_pot%p( 80)=   0.297004000E+00_Rkind
      QModel%no3minus_pot%p( 81)=  -0.120273480E+00_Rkind
      QModel%no3minus_pot%p( 82)=   0.191271270E+01_Rkind
      QModel%no3minus_pot%p( 83)=  -0.226871020E+01_Rkind
      QModel%no3minus_pot%p( 84)=  -0.114691680E+02_Rkind
      QModel%no3minus_pot%p( 85)=  -0.897240320E+02_Rkind
      QModel%no3minus_pot%p( 86)=   0.147408570E+02_Rkind
      QModel%no3minus_pot%p( 87)=   0.123400003E+00_Rkind
      QModel%no3minus_pot%p( 88)=  -0.624627310E+00_Rkind
      QModel%no3minus_pot%p( 89)=   0.984178850E-02_Rkind
      QModel%no3minus_pot%p( 90)=  -0.212520820E+01_Rkind
      QModel%no3minus_pot%p( 91)=  -0.127773090E+01_Rkind
      QModel%no3minus_pot%p( 92)=   0.408838780E+01_Rkind
      QModel%no3minus_pot%p( 93)=   0.310835790E+01_Rkind
      QModel%no3minus_pot%p( 94)=  -0.787746700E+01_Rkind
      QModel%no3minus_pot%p( 95)=   0.572379100E+01_Rkind
      QModel%no3minus_pot%p( 96)=  -0.883120710E-01_Rkind
      QModel%no3minus_pot%p( 97)=  -0.845179740E+01_Rkind
      QModel%no3minus_pot%p( 98)=  -0.923696170E+01_Rkind
      QModel%no3minus_pot%p( 99)=   0.150396360E+02_Rkind
      QModel%no3minus_pot%p(100)=   0.129763030E+02_Rkind
      QModel%no3minus_pot%p(101)=  -0.127820180E+02_Rkind
      QModel%no3minus_pot%p(102)=   0.724479130E+01_Rkind
      QModel%no3minus_pot%p(103)=  -0.226340140E+02_Rkind
      QModel%no3minus_pot%p(104)=   0.193995670E+02_Rkind
      QModel%no3minus_pot%p(105)=  -0.101918030E+02_Rkind
      QModel%no3minus_pot%p(106)=   0.120124630E+03_Rkind
      QModel%no3minus_pot%p(107)=  -0.195995920E+02_Rkind
      QModel%no3minus_pot%p(108)=  -0.821429710E+01_Rkind
      QModel%no3minus_pot%p(109)=  -0.146887710E+02_Rkind
      QModel%no3minus_pot%p(110)=  -0.731600850E+02_Rkind
      QModel%no3minus_pot%p(111)=   0.673921420E+02_Rkind
      QModel%no3minus_pot%p(112)=   0.150943880E+02_Rkind
      QModel%no3minus_pot%p(113)=  -0.867078970E+01_Rkind
      QModel%no3minus_pot%p(114)=   0.147406530E+03_Rkind
      QModel%no3minus_pot%p(115)=  -0.260872930E+02_Rkind
      QModel%no3minus_pot%p(116)=  -0.840659750E+03_Rkind
      QModel%no3minus_pot%p(117)=   0.504515730E+03_Rkind
      QModel%no3minus_pot%p(118)=   0.621689150E+03_Rkind
      QModel%no3minus_pot%p(119)=  -0.524583500E+03_Rkind
      QModel%no3minus_pot%p(120)=  -0.382318100E+02_Rkind
      QModel%no3minus_pot%p(121)=   0.201076520E+03_Rkind
      QModel%no3minus_pot%p(122)=  -0.733776780E+02_Rkind
      QModel%no3minus_pot%p(123)=  -0.296733560E+02_Rkind
      QModel%no3minus_pot%p(124)=  -0.211305860E+03_Rkind
      QModel%no3minus_pot%p(125)=  -0.382490340E+03_Rkind
             !p(126)=   0.626682910E+03_Rkind ! to fix a symmetry pbs.
             !p(127)=   0.648222640E+03_Rkind ! we take the average. DML 05/10/2021
      QModel%no3minus_pot%p(126)=   (0.626682910E+03_Rkind+0.648222640E+03_Rkind)/TWO
      QModel%no3minus_pot%p(127)=   (0.626682910E+03_Rkind+0.648222640E+03_Rkind)/TWO
      QModel%no3minus_pot%p(128)=   0.510189060E+02_Rkind
      QModel%no3minus_pot%p(129)=  -0.300347270E+03_Rkind
      QModel%no3minus_pot%p(130)=  -0.131928350E+04_Rkind
      QModel%no3minus_pot%p(131)=  -0.902324530E+01_Rkind
      QModel%no3minus_pot%p(132)=  -0.121951920E+01_Rkind
      QModel%no3minus_pot%p(133)=   0.617623300E+03_Rkind
      QModel%no3minus_pot%p(134)=   0.307972350E+03_Rkind
      QModel%no3minus_pot%p(135)=  -0.228997370E+02_Rkind
      QModel%no3minus_pot%p(136)=  -0.163914390E+03_Rkind
      QModel%no3minus_pot%p(137)=  -0.169780300E+00_Rkind
      QModel%no3minus_pot%p(138)=   0.125500700E+01_Rkind
      QModel%no3minus_pot%p(139)=   0.482992570E+00_Rkind
      QModel%no3minus_pot%p(140)=   0.317988550E+01_Rkind
      QModel%no3minus_pot%p(141)=  -0.436770420E+01_Rkind
      QModel%no3minus_pot%p(142)=   0.356941810E+00_Rkind
      QModel%no3minus_pot%p(143)=   0.900703840E+01_Rkind
      QModel%no3minus_pot%p(144)=  -0.282621440E+02_Rkind
      QModel%no3minus_pot%p(145)=  -0.137130480E+02_Rkind
      QModel%no3minus_pot%p(146)=   0.604692600E+01_Rkind
      QModel%no3minus_pot%p(147)=  -0.276681540E+02_Rkind
      QModel%no3minus_pot%p(148)=  -0.564288980E+02_Rkind
      QModel%no3minus_pot%p(149)=   0.279994500E+02_Rkind
      QModel%no3minus_pot%p(150)=  -0.750562690E+02_Rkind
      QModel%no3minus_pot%p(151)=  -0.802644670E+02_Rkind
      QModel%no3minus_pot%p(152)=   0.511596150E+02_Rkind
      QModel%no3minus_pot%p(153)=  -0.274657630E+02_Rkind
      QModel%no3minus_pot%p(154)=   0.885353980E+02_Rkind
      QModel%no3minus_pot%p(155)=   0.203854730E+00_Rkind
      QModel%no3minus_pot%p(156)=   0.343485190E+01_Rkind
      QModel%no3minus_pot%p(157)=   0.115761950E+01_Rkind
      QModel%no3minus_pot%p(158)=  -0.108960670E+02_Rkind
      QModel%no3minus_pot%p(159)=   0.810512110E+01_Rkind
      QModel%no3minus_pot%p(160)=  -0.149700950E+01_Rkind
      QModel%no3minus_pot%p(161)=  -0.632516450E+02_Rkind
      QModel%no3minus_pot%p(162)=  -0.182788980E+03_Rkind
      QModel%no3minus_pot%p(163)=  -0.607637350E+02_Rkind
      QModel%no3minus_pot%p(164)=   0.796587110E+02_Rkind
      QModel%no3minus_pot%p(165)=  -0.325214470E+02_Rkind
      QModel%no3minus_pot%p(166)=   0.246171810E+03_Rkind
      QModel%no3minus_pot%p(167)=  -0.139468350E+02_Rkind
      QModel%no3minus_pot%p(168)=   0.659294940E+03_Rkind
      QModel%no3minus_pot%p(169)=  -0.142032340E+01_Rkind
      QModel%no3minus_pot%p(170)=   0.418184160E+03_Rkind
      QModel%no3minus_pot%p(171)=   0.151015790E+00_Rkind
      QModel%no3minus_pot%p(172)=  -0.796880880E+01_Rkind
      QModel%no3minus_pot%p(173)=   0.217108089E+00_Rkind
      QModel%no3minus_pot%p(174)=  -0.232760231E+02_Rkind
      QModel%no3minus_pot%p(175)=   0.269121775E+03_Rkind
      QModel%no3minus_pot%p(176)=  -0.365601790E+03_Rkind
      QModel%no3minus_pot%p(177)=   0.245177337E+02_Rkind
      QModel%no3minus_pot%p(178)=   0.737181669E+02_Rkind
      QModel%no3minus_pot%p(179)=   0.230068180E+04_Rkind
      QModel%no3minus_pot%p(180)=  -0.919786868E+03_Rkind
      QModel%no3minus_pot%p(181)=   0.203758567E+01_Rkind
      QModel%no3minus_pot%p(182)=   0.125347943E+02_Rkind
      QModel%no3minus_pot%p(183)=  -0.240018368E+02_Rkind
      QModel%no3minus_pot%p(184)=   0.315636912E+02_Rkind
      QModel%no3minus_pot%p(185)=  -0.586391234E+03_Rkind
      QModel%no3minus_pot%p(186)=  -0.541519483E+02_Rkind
      QModel%no3minus_pot%p(187)=   0.541543516E+02_Rkind
      QModel%no3minus_pot%p(188)=  -0.110047897E+01_Rkind
      QModel%no3minus_pot%p(189)=  -0.984270831E+01_Rkind
      QModel%no3minus_pot%p(190)=   0.176431382E+02_Rkind
      QModel%no3minus_pot%p(191)=   0.619129772E+03_Rkind
      QModel%no3minus_pot%p(192)=  -0.309888310E+02_Rkind
      QModel%no3minus_pot%p(193)=   0.516577868E+03_Rkind
      QModel%no3minus_pot%p(194)=   0.763453420E+03_Rkind
      QModel%no3minus_pot%p(195)=   0.476767245E+00_Rkind
      QModel%no3minus_pot%p(196)=   0.396602201E+00_Rkind
      QModel%no3minus_pot%p(197)=  -0.157742464E+02_Rkind
      QModel%no3minus_pot%p(198)=  -0.345412065E+02_Rkind
      QModel%no3minus_pot%p(199)=   0.106966811E+02_Rkind
      QModel%no3minus_pot%p(200)=   0.685088030E+03_Rkind
      QModel%no3minus_pot%p(201)=   0.673178507E+03_Rkind
      QModel%no3minus_pot%p(202)=  -0.106366340E+04_Rkind
      QModel%no3minus_pot%p(203)=  -0.332883674E+04_Rkind
      QModel%no3minus_pot%p(204)=   0.481034455E+04_Rkind
      QModel%no3minus_pot%p(205)=  -0.189700034E+02_Rkind
      QModel%no3minus_pot%p(206)=   0.152614812E+03_Rkind
      QModel%no3minus_pot%p(207)=  -0.709471913E+03_Rkind
      QModel%no3minus_pot%p(208)=   0.126748116E+03_Rkind
      QModel%no3minus_pot%p(209)=  -0.240266010E+03_Rkind
      QModel%no3minus_pot%p(210)=   0.148598985E+04_Rkind
      QModel%no3minus_pot%p(211)=   0.115713133E+04_Rkind
      QModel%no3minus_pot%p(212)=  -0.897094470E+03_Rkind
      QModel%no3minus_pot%p(213)=  -0.834863905E+03_Rkind
      QModel%no3minus_pot%p(214)=   0.545372351E+02_Rkind
      QModel%no3minus_pot%p(215)=   0.317678018E+03_Rkind
      QModel%no3minus_pot%p(216)=   0.899708039E+02_Rkind
      QModel%no3minus_pot%p(217)=   0.461995866E+03_Rkind
      QModel%no3minus_pot%p(218)=  -0.657237461E+03_Rkind
      QModel%no3minus_pot%p(219)=  -0.129355471E+04_Rkind
      QModel%no3minus_pot%p(220)=  -0.139661058E+03_Rkind
      QModel%no3minus_pot%p(221)=  -0.189084040E+03_Rkind
      QModel%no3minus_pot%p(222)=  -0.209609912E+02_Rkind
      QModel%no3minus_pot%p(223)=   0.312119590E+02_Rkind
      QModel%no3minus_pot%p(224)=  -0.875470382E+02_Rkind
      QModel%no3minus_pot%p(225)=  -0.578320555E+02_Rkind
      QModel%no3minus_pot%p(226)=   0.819108852E+01_Rkind
      QModel%no3minus_pot%p(227)=   0.120109292E+03_Rkind
      QModel%no3minus_pot%p(228)=   0.207725902E+03_Rkind
      QModel%no3minus_pot%p(229)=  -0.220092259E+03_Rkind
      QModel%no3minus_pot%p(230)=  -0.219373924E+03_Rkind
      QModel%no3minus_pot%p(231)=  -0.104377843E+04_Rkind
      QModel%no3minus_pot%p(232)=   0.188456908E+03_Rkind
      QModel%no3minus_pot%p(233)=  -0.187745233E+02_Rkind
      QModel%no3minus_pot%p(234)=   0.115501296E+03_Rkind
      QModel%no3minus_pot%p(235)=  -0.680341930E+02_Rkind
      QModel%no3minus_pot%p(236)=   0.135039895E+03_Rkind
      QModel%no3minus_pot%p(237)=  -0.286175591E+03_Rkind
      QModel%no3minus_pot%p(238)=  -0.150033749E+02_Rkind
      QModel%no3minus_pot%p(239)=   0.642889917E+03_Rkind
      QModel%no3minus_pot%p(240)=   0.496437690E+01_Rkind
      QModel%no3minus_pot%p(241)=   0.579415330E+01_Rkind
      QModel%no3minus_pot%p(242)=   0.819178020E+00_Rkind
      QModel%no3minus_pot%p(243)=  -0.641777110E+01_Rkind
      QModel%no3minus_pot%p(244)=   0.162797280E+03_Rkind
      QModel%no3minus_pot%p(245)=  -0.180520690E+03_Rkind
      QModel%no3minus_pot%p(246)=  -0.536904030E+03_Rkind
      QModel%no3minus_pot%p(247)=   0.604811650E+03_Rkind
      QModel%no3minus_pot%p(248)=   0.284787840E+03_Rkind
      QModel%no3minus_pot%p(249)=  -0.135967270E+03_Rkind
      QModel%no3minus_pot%p(250)=   0.153868400E+02_Rkind
      QModel%no3minus_pot%p(251)=  -0.240188550E+03_Rkind
      QModel%no3minus_pot%p(252)=   0.314345890E+03_Rkind
      QModel%no3minus_pot%p(253)=   0.664042750E+03_Rkind
      QModel%no3minus_pot%p(254)=  -0.120296430E+04_Rkind
      QModel%no3minus_pot%p(255)=   0.296121720E+03_Rkind
      QModel%no3minus_pot%p(256)=  -0.123934530E+03_Rkind
      QModel%no3minus_pot%p(257)=  -0.652689590E+02_Rkind
      QModel%no3minus_pot%p(258)=   0.147837820E+04_Rkind
      QModel%no3minus_pot%p(259)=  -0.532670070E+04_Rkind
      QModel%no3minus_pot%p(260)=  -0.202987460E+04_Rkind
      QModel%no3minus_pot%p(261)=   0.272007950E+04_Rkind
      QModel%no3minus_pot%p(262)=   0.454035570E+04_Rkind
      QModel%no3minus_pot%p(263)=  -0.174507240E+04_Rkind
      QModel%no3minus_pot%p(264)=   0.259475430E+03_Rkind
      QModel%no3minus_pot%p(265)=   0.640285120E+03_Rkind
      QModel%no3minus_pot%p(266)=  -0.489340180E+01_Rkind
      QModel%no3minus_pot%p(267)=  -0.180109170E+01_Rkind
      QModel%no3minus_pot%p(268)=   0.214753300E+02_Rkind
      QModel%no3minus_pot%p(269)=  -0.103364920E+02_Rkind
      QModel%no3minus_pot%p(270)=   0.336812050E+01_Rkind
      QModel%no3minus_pot%p(271)=   0.331311770E+01_Rkind
      QModel%no3minus_pot%p(272)=  -0.184333060E+02_Rkind
      QModel%no3minus_pot%p(273)=   0.363236360E+02_Rkind
      QModel%no3minus_pot%p(274)=  -0.733841910E+02_Rkind
      QModel%no3minus_pot%p(275)=   0.733114270E+02_Rkind
      QModel%no3minus_pot%p(276)=  -0.103140980E+03_Rkind
      QModel%no3minus_pot%p(277)=   0.643779480E+02_Rkind
      QModel%no3minus_pot%p(278)=  -0.300323700E+03_Rkind
      QModel%no3minus_pot%p(279)=   0.116331260E+02_Rkind
      QModel%no3minus_pot%p(280)=   0.333795200E+03_Rkind
      QModel%no3minus_pot%p(281)=   0.151003270E+02_Rkind
      QModel%no3minus_pot%p(282)=  -0.294676480E+02_Rkind
      QModel%no3minus_pot%p(283)=   0.137832210E+02_Rkind
      QModel%no3minus_pot%p(284)=  -0.150401300E+02_Rkind
      QModel%no3minus_pot%p(285)=   0.128243380E+03_Rkind
      QModel%no3minus_pot%p(286)=   0.515344760E+03_Rkind
      QModel%no3minus_pot%p(287)=  -0.104975170E+04_Rkind
      QModel%no3minus_pot%p(288)=   0.696065160E+03_Rkind
      QModel%no3minus_pot%p(289)=  -0.850704770E+02_Rkind
      QModel%no3minus_pot%p(290)=   0.552233050E+03_Rkind
      QModel%no3minus_pot%p(291)=  -0.126727650E+03_Rkind
      QModel%no3minus_pot%p(292)=   0.553415220E+03_Rkind
      QModel%no3minus_pot%p(293)=  -0.532309150E+02_Rkind
      QModel%no3minus_pot%p(294)=  -0.112746260E+04_Rkind
      QModel%no3minus_pot%p(295)=   0.186326010E+03_Rkind
      QModel%no3minus_pot%p(296)=  -0.364935110E+03_Rkind
      QModel%no3minus_pot%p(297)=   0.207343810E+03_Rkind
      QModel%no3minus_pot%p(298)=  -0.293492870E+03_Rkind
      QModel%no3minus_pot%p(299)=  -0.922632000E+01_Rkind
      QModel%no3minus_pot%p(300)=  -0.199198270E+03_Rkind
      QModel%no3minus_pot%p(301)=   0.302964050E+03_Rkind
      QModel%no3minus_pot%p(302)=  -0.278387770E+03_Rkind
      QModel%no3minus_pot%p(303)=   0.140546710E+03_Rkind
      QModel%no3minus_pot%p(304)=   0.529782740E+03_Rkind
      QModel%no3minus_pot%p(305)=  -0.273028890E+03_Rkind
      QModel%no3minus_pot%p(306)=  -0.149698510E+03_Rkind
      QModel%no3minus_pot%p(307)=  -0.142095790E+04_Rkind
      QModel%no3minus_pot%p(308)=   0.748850490E+03_Rkind
      QModel%no3minus_pot%p(309)=   0.925607750E+03_Rkind
      QModel%no3minus_pot%p(310)=  -0.460069010E+03_Rkind
      QModel%no3minus_pot%p(311)=  -0.206325460E+03_Rkind
      QModel%no3minus_pot%p(312)=   0.857738680E-02_Rkind
      ! Vertical Energy
      QModel%no3minus_pot%p(313)=  -0.279544420E+03_Rkind
      QModel%no3minus_pot%p(313)=   ZERO               !-0.279544420E+03_Rkind
      QModel%no3minus_pot%p(314)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(315)=   0.109055120E+00_Rkind
      QModel%no3minus_pot%p(316)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(317)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(318)=   0.696349180E-01_Rkind
      QModel%no3minus_pot%p(319)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(320)=   0.281557390E+00_Rkind
      QModel%no3minus_pot%p(321)=  -0.827685040E-01_Rkind
      QModel%no3minus_pot%p(322)=   0.156327730E-01_Rkind
      QModel%no3minus_pot%p(323)=  -0.416376780E-03_Rkind
      QModel%no3minus_pot%p(324)=   0.630990090E-04_Rkind
      QModel%no3minus_pot%p(325)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(326)=   0.750000000E+00_Rkind
      QModel%no3minus_pot%p(327)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(328)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(329)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(330)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(331)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(332)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(333)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(334)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(335)=   0.000000000E+00_Rkind
      QModel%no3minus_pot%p(336)=   0.953536010E+01_Rkind
      QModel%no3minus_pot%p(337)=   0.103373110E+01_Rkind

      ! Vertical energy
      !    e0ref=  -0.279402540D+03
      QModel%no3minus_pot%e0ref=  -0.279402540E+03_Rkind +0.279544420E+03_Rkind
      QModel%no3minus_pot%e0rho=   0.800000000E-01_Rkind

      QModel%no3minus_pot%pst(1,  1)=  1
      QModel%no3minus_pot%pst(2,  1)=  6
      QModel%no3minus_pot%pst(1,  2)=  7
      QModel%no3minus_pot%pst(2,  2)=  3
      QModel%no3minus_pot%pst(1,  3)= 10
      QModel%no3minus_pot%pst(2,  3)=  6
      QModel%no3minus_pot%pst(1,  4)= 16
      QModel%no3minus_pot%pst(2,  4)=  6
      QModel%no3minus_pot%pst(1,  5)= 22
      QModel%no3minus_pot%pst(2,  5)= 28
      QModel%no3minus_pot%pst(1,  6)= 50
      QModel%no3minus_pot%pst(2,  6)= 10
      QModel%no3minus_pot%pst(1,  7)= 60
      QModel%no3minus_pot%pst(2,  7)= 10
      QModel%no3minus_pot%pst(1,  8)= 70
      QModel%no3minus_pot%pst(2,  8)=  9
      QModel%no3minus_pot%pst(1,  9)= 79
      QModel%no3minus_pot%pst(2,  9)=  9
      QModel%no3minus_pot%pst(1, 10)= 88
      QModel%no3minus_pot%pst(2, 10)= 49
      QModel%no3minus_pot%pst(1, 11)=137
      QModel%no3minus_pot%pst(2, 11)= 18
      QModel%no3minus_pot%pst(1, 12)=155
      QModel%no3minus_pot%pst(2, 12)= 18
      QModel%no3minus_pot%pst(1, 13)=173
      QModel%no3minus_pot%pst(2, 13)=  4
      QModel%no3minus_pot%pst(1, 14)=177
      QModel%no3minus_pot%pst(2, 14)=  4
      QModel%no3minus_pot%pst(1, 15)=181
      QModel%no3minus_pot%pst(2, 15)=  7
      QModel%no3minus_pot%pst(1, 16)=188
      QModel%no3minus_pot%pst(2, 16)=  7
      QModel%no3minus_pot%pst(1, 17)=195
      QModel%no3minus_pot%pst(2, 17)=  2
      QModel%no3minus_pot%pst(1, 18)=197
      QModel%no3minus_pot%pst(2, 18)=  8
      QModel%no3minus_pot%pst(1, 19)=205
      QModel%no3minus_pot%pst(2, 19)= 15
      QModel%no3minus_pot%pst(1, 20)=220
      QModel%no3minus_pot%pst(2, 20)=  3
      QModel%no3minus_pot%pst(1, 21)=223
      QModel%no3minus_pot%pst(2, 21)=  3
      QModel%no3minus_pot%pst(1, 22)=226
      QModel%no3minus_pot%pst(2, 22)=  7
      QModel%no3minus_pot%pst(1, 23)=233
      QModel%no3minus_pot%pst(2, 23)=  7
      QModel%no3minus_pot%pst(1, 24)=240
      QModel%no3minus_pot%pst(2, 24)= 26
      QModel%no3minus_pot%pst(1, 25)=266
      QModel%no3minus_pot%pst(2, 25)= 46
      QModel%no3minus_pot%pst(1, 26)=312
      QModel%no3minus_pot%pst(2, 26)=  1
      QModel%no3minus_pot%pst(1, 27)=313
      QModel%no3minus_pot%pst(2, 27)=  1
      QModel%no3minus_pot%pst(1, 28)=314
      QModel%no3minus_pot%pst(2, 28)= 11
      QModel%no3minus_pot%pst(1, 29)=325
      QModel%no3minus_pot%pst(2, 29)= 11
      QModel%no3minus_pot%pst(1, 30)=336
      QModel%no3minus_pot%pst(2, 30)=  2

      Write(out_unit,*)'# NO3 E" 2 states 16/03/2015'
      Write(out_unit,*)'# tmc-ang-prec-37um2b_extra_b5_20150316.par  '
      Write(out_unit,*)'# switching method with total energy and polynomial'
      Write(out_unit,*)'# iref ',QModel%no3minus_pot%iref
      Write(out_unit,*)'# iref is not used'
      Write(out_unit,*)'# npar ',QModel%no3minus_pot%n
      !do i=1,size(QModel%no3minus_pot%pst,dim=2)
      !  Write(out_unit,2)i,QModel%no3minus_pot%pst(1,i),i,QModel%no3minus_pot%pst(2,i)
      !enddo
2     format('# pst(1,',i3,')=',i5,' pst(2,',i3,')=',i5)

      !do i=1,size(QModel%no3minus_pot%p)
      !  Write(out_unit,1)i,QModel%no3minus_pot%p(i)
      !enddo
1    format('# p(',i3,')=',f16.10)

      Write(out_unit,*)'# e0ref au and ev',QModel%no3minus_pot%e0ref,QModel%no3minus_pot%e0ref*27.21_Rkind
      Write(out_unit,*)'# e0rho ',QModel%no3minus_pot%e0rho


      ! tmc parameter
      ! now in p array directly
      !      do i=1,11
      !        a(1,i)=ZERO
      !        a(2,i)=ZERO
      !      enddo
      !      do i=1,11
      !        a(1,i) = p(pst(1,28)+i-1)    ! NO
      !        a(2,i) = p(pst(1,29)+i-1)    ! OO
      !      enddo
      ! copy of first p array into a -> tmc common
      ! only reason (not a good one): not to put the
      ! no3minus_pot/ common (parameter) into the
      ! coordinate transformation code.
      ! 20 is arbitrary number
      !ICI
      j=1
      do i=QModel%no3minus_pot%pst(1,28),QModel%no3minus_pot%pst(1,28)+20
        QModel%tmc%a(j)=QModel%no3minus_pot%p(i)
        !Write(out_unit,*)'ici',i,QModel%p(i)
        j=j+1
      enddo
      QModel%tmc%a(2)=abs(QModel%tmc%a(2))  ! DML: it was in QML_NO3_ff subroutine
      QModel%tmc%a(5)=abs(QModel%tmc%a(5))  ! DML: it was in QML_NO3_ff subroutine
      QModel%tmc%a(4)=abs(QModel%tmc%a(4))  ! DML: it was in QML_NO3_ff subroutine

      QModel%tmc%npoly(:)=[5,5]
      Write(out_unit,*)'# TMC     npoly(1) : ',QModel%tmc%npoly(1)
      Write(out_unit,*)'# TMC     npoly(2) : ',QModel%tmc%npoly(2)

    CASE Default
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' This option is not possible. option:',QModel%option
        write(out_unit,*) ' Its value MUST be 1'
        STOP
    END SELECT


    IF (debug) write(out_unit,*) 'init Q0 of NO3'
    QModel%Q0 = [ZERO,ZERO,ZERO,ZERO,ZERO,ZERO]
    IF (debug) write(out_unit,*) 'init d0GGdef of NO3'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_NO3
!> @brief Subroutine wich prints the QML_NO3 parameters.
!!
!! @param QModel            CLASS(QML_NO3_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_NO3(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_NO3_t),   intent(in) :: QModel
    integer,               intent(in) :: nio
 ! A Viel 2013.09.05
    ! rs.f put at the end
    !
    ! A Viel 2012.07.13
    ! Version with iref=1 hard coded
    !
    ! ATTENTION limited to 32 threads (mx=32)
    !
    ! A Viel 2012 06.13
    ! generic functions needed for jt problems
    ! NO3 case
    ! TAKEN from genetic suite (W. Eisfeld)
    ! version -prec-
    write(nio,*) 'NO3 current parameters'
    write(nio,*)
    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  adiabatic:      ',QModel%adiabatic
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)

    SELECT CASE (QModel%option)

    CASE (1)


    CASE Default
        write(out_unit,*) ' ERROR in Write_QML_NO3 '
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1'
        STOP
    END SELECT

    write(nio,*)
    write(nio,*) 'end NO3 current parameters'

  END SUBROUTINE Write_QML_NO3

  SUBROUTINE Cart_TO_Q_QML_NO3(QModel,dnX,dnQ,nderiv)
    USE QDUtil_m,         ONLY : TO_string
    USE ADdnSVM_m
    IMPLICIT NONE
  
      CLASS(QML_NO3_t),        intent(in)    :: QModel
      TYPE (dnS_t),            intent(in)    :: dnX(:,:)
      TYPE (dnS_t),            intent(inout) :: dnQ(:)
      integer,                 intent(in)    :: nderiv
  
      integer                   :: i,j,ij
  
      !----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Cart_TO_Q_QML_NO3'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'size(dnQ)',size(dnQ)
        write(out_unit,*) 'dnQ:'
        DO j=1,size(dnQ,dim=1)
          CALL Write_dnS(dnQ(j),out_unit,info='dnQ('// TO_string(j) // ')')
        END DO
        write(out_unit,*) 'shape dnX',shape(dnX)
        write(out_unit,*) 'dnX'
        DO i=1,size(dnX,dim=2)
        DO j=1,size(dnX,dim=1)
          CALL Write_dnS(dnX(j,i),out_unit)
        END DO
        END DO
        flush(out_unit)
      END IF

      ij = 0
      DO i=1,size(dnX,dim=2)
      DO j=1,size(dnX,dim=1)
        ij = ij + 1
        dnQ(ij) = dnX(j,i)
      END DO
      END DO

      IF (debug) THEN
        DO j=1,size(dnQ,dim=1)
          CALL Write_dnS(dnQ(j),out_unit,info='dnQ('// TO_string(j) // ')')
        END DO
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
  END SUBROUTINE Cart_TO_Q_QML_NO3
!> @brief Subroutine wich calculates the NO3 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_NO3_t):  derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_NO3(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_NO3_t),     intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    SELECT CASE (QModel%option)

    CASE (1)
      CALL EvalPot1_QML_NO3(Mat_OF_PotDia,dnQ,QModel,nderiv)

    CASE Default
        write(out_unit,*) ' ERROR in EvalPot_QML_NO3'
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1'
        STOP
    END SELECT

  END SUBROUTINE EvalPot_QML_NO3

!> @brief Subroutine wich calculates the NO3 potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param NO3Pot          TYPE(NO3Pot_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot1_QML_NO3(Mat_OF_PotDia,dnQ,QModel,nderiv)
    !Unpublished model potential (yet)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)    :: dnQ(:)
    TYPE(QML_NO3_t),     intent(in)    :: QModel
    integer,             intent(in)    :: nderiv


    !local variables (derived type). They have to be deallocated
    TYPE (dnS_t)        :: dnPot,BLA,Tors,HOOP
    TYPE (dnS_t)        :: MorseBLAP,Morsemin,Overlap
    TYPE (dnS_t)        :: Hdir2D,Hct2D
    real (kind=Rkind)   :: d1,d4

    real (kind=Rkind) :: XNO3(0:3,3)
    real (kind=Rkind) :: vdia(2,2)

    XNO3(0,:) = get_d0(dnQ( 1: 3))
    XNO3(1,:) = get_d0(dnQ( 4: 6))
    XNO3(2,:) = get_d0(dnQ( 7: 9))
    XNO3(3,:) = get_d0(dnQ(10:12))

    ! test geometry
    ! 0 is N
    !XNO3(0,1)=2.9506658729288038e-003_Rkind
    !XNO3(0,2)=ZERO()
    !XNO3(0,3)=-3.7754294953142864e-003_Rkind

    ! 1,2,3 are the 3 Oxygen atoms
    !XNO3(1,1)=9.5060742943492427e-004_Rkind
    !XNO3(1,2)=ZERO
    !XNO3(1,3)=2.3817433236872598_Rkind

    !XNO3(2,1)=2.0639330934921434_Rkind
    !XNO3(2,2)=ZERO
    !XNO3(2,3)=-1.1939442518313312_Rkind

    !XNO3(3,1)=-2.0674669214991690_Rkind
    !XNO3(3,2)=ZERO
    !XNO3(3,3)=-1.1844937951551615_Rkind

    CALL QML_NO3_potential(XNO3,vdia,QModel)
    Mat_OF_PotDia(1,1) = vdia(1,1)
    Mat_OF_PotDia(1,2) = vdia(1,2)
    Mat_OF_PotDia(2,1) = vdia(2,1)
    Mat_OF_PotDia(2,2) = vdia(2,2)


  END SUBROUTINE EvalPot1_QML_NO3


  ! A Viel 30.09.2021
! 'clean' version of the NO3 E'' potential used in  JCP 146, 034303 (2017)
!
! return 2x2 diabatic matrix (atomic unit)
! input cartesian coordinates N, O, O, O
!
! subroutine init_pot_para : initialisation
! subroutine potential :
!            trans_coord(x,qinter)  coordinates transformation
!            diabatic matrix elements
!-----------------------------------------------------------
! deux test geometries: expected output
! cart           0   2.9506658729288038E-003   0.0000000000000000       -3.7754294953142864E-003
! cart           1   9.5060742943492427E-004   0.0000000000000000        2.3817433236872598
! cart           2   2.0639330934921434        0.0000000000000000       -1.1939442518313312
! cart           3  -2.0674669214991690        0.0000000000000000       -1.1844937951551615
! vdia 1  -2.0938061190543672E-003   3.4218375116754394E-004
! vdia 2   3.4218375116754394E-004  -1.2190023654554047E-003
! adiacm  -485.42325517890794       -241.65416575993908

! cart           0 -0.28276917777022903      -0.41734715646392218       -1.0997467054058241E-002
! cart           1  -9.1098922338089527E-002 -0.13445551772839726        2.2271508290107942
! cart           2   2.2747337906099050      -0.13445551772839726      -0.81179512186018443
! cart           3  -1.9360788283195443       0.63428610976387545       -1.4057277504727415
! vdia 1   5.3234458438400425E-002  -5.0114745937633876E-003
! vdia 2  -5.0114745937633876E-003   4.8620099394501577E-002
!-----------------------------------------------------------
  subroutine QML_NO3_potential(x,vdia,QModel)
    !$      use omp_lib
    USE QDUtil_m
    implicit none
          real (kind=Rkind),   intent(in)         :: x(0:3,3)  ! DML erreur x(0-3,3)
          real (kind=Rkind),   intent(inout)      :: vdia(2,2)
          TYPE(QML_NO3_t),     intent(in), target :: QModel

          integer itd ,i,itemp
          real (kind=Rkind) :: v
          real (kind=Rkind) :: qinter(7),e,damp
          real (kind=Rkind) :: eref,w,z,damplow
          real (kind=Rkind) :: polynom,polynom2,vjt6
          real (kind=Rkind) :: vcoup_ea,vcoup_ea2,vcoup_aa,vcoup_ee
          real (kind=Rkind) :: vcoup_eea2,vcoup_eaa2,vcoup_eea1
          real (kind=Rkind) :: dum
          real (kind=Rkind) :: eloworder(2,2),vloworder(2),uloworder(2,2)
          real (kind=Rkind) :: fv1(2),fv2(2)

          real (kind=Rkind) :: etest(2,2),vtest(2),utest(2,2)

          integer,           pointer :: iref,pst(:,:)
          real (kind=Rkind), pointer :: p(:),e0ref,e0rho

    
          integer mx
          TYPE(acoord_vjt_vee_t), target  :: avv

          real(kind=Rkind), pointer :: pa(:), pb(:)
          real(kind=Rkind), pointer :: va(:), vb(:)
          real(kind=Rkind), pointer :: wa(:), wb(:)
          real(kind=Rkind), pointer :: za(:), zb(:)
          real(kind=Rkind), pointer :: vee(:),  wee(:), zee(:)

          pa => avv%pa
          pb => avv%pb
          va => avv%va
          vb => avv%vb
          wa => avv%wa
          wb => avv%wb
          za => avv%za
          zb => avv%zb

          vee => avv%vee
          wee => avv%wee
          zee => avv%zee

          iref  => QModel%no3minus_pot%iref
          pst   => QModel%no3minus_pot%pst
          p     => QModel%no3minus_pot%p
          e0ref => QModel%no3minus_pot%e0ref
          e0rho => QModel%no3minus_pot%e0rho

          call QML_NO3_trans_coord(x,qinter,QModel%transformcoordblock,QModel%tmc)
          ! write(out_unit,'(7g20.10)')(qinter(i),i=1,7)

    ! PES diabatic matrix computation
    
    
    ! vertical excitation
    !      e=p(pst(1,27))  !in low order term eref
           e=ZERO
    !..   damping factor for umbrella motion, using hyperradius:
    ! qinter(7): (no1+no2+no3)/3
          damp=exp(-abs(p(pst(1,26))*qinter(7)))       !damping of pure umbrella mode
    
    ! precompute terms
          CALL QML_NO3_vwzprec(qinter,avv)

          ! reference low order surface -> in eloworder
    
    ! computing first order matrix (iref=0 in Wofgang Eisfeld code)
          itemp=0
    ! vertical excitation
          eref=p(pst(1,27))
    !
    ! symmetric stretching
    !..   a1 symmetric stretch:
    !      call QML_NO3_Epolynom(dum,itd, p(pst(1,1)), pst(2,1), itemp)
          dum=pa(1) * p(pst(1,1))+pa(2) * p(pst(1,1)+1)
    !       write(out_unit,*)'dum low',dum
          eref=eref+dum
    
    ! umbrella mode
    !..   a2 umbrella:
    !      call QML_NO3_Epolynom2(dum,itd, p(pst(1,2)), pst(2,2), itemp)
          dum=pb(1) * p(pst(1,2))
    !       write(out_unit,*)'dum low',dum
          eref=eref+dum*damp
    
    ! e stretching mode:
    ! e bending mode:
    !      CALL QML_NO3_Evjt6a(dum,itd, p(pst(1,3)), pst(2,3), itemp)
          dum= p(pst(1,3)) * va(1)
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum
    !      call Evjt6b(dum,itd, p(pst(1,4)), pst(2,4), itemp)
          dum= p(pst(1,4)) * vb(1)
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum
    
    !..   mode-mode coupling of a1-e
    !      call Evcoup_e1a(dum,itd, p(pst(1,6)), pst(2,6), itemp)  ! e-a coupling 1. e mode
          dum= (p(pst(1,6))*pa(1)*va(1))
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum
    !      call Evcoup_e2a(dum,itd, p(pst(1,7)), pst(2,7), itemp)  ! e-a coupling 2. e mode
          dum=(p(pst(1,7))*pa(1)*vb(1))
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum
    
    !..   mode-mode coupling of a2-e
    !      call Evcoup_e1a2(dum,itd, p(pst(1,13)), pst(2,13), itemp)     !  coupling with asym. str.
          dum=p(pst(1,13))*pb(1)*va(1)
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum*damp
    !      call Evcoup_e2a2(dum,itd, p(pst(1,14)), pst(2,14), itemp)     !  coupling with asym. bend
          dum=p(pst(1,14))*pb(1)*vb(1)
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum*damp
    
    
    !..   mode-mode coupling of a1-a2
    !      call Evcoup_aa(dum,itd,p(pst(1,17)),pst(2,17), itemp)         !  coupling of a1 and a2
          dum=p(pst(1,17)) * pa(1) * pb(1)
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum*damp
    
    
    !..   mode-mode coupling e-e
    !      call Evcoup_ee(dum,itd, p(pst(1,5)), pst(2,5), itemp)        ! e-e coupling
          dum=p(pst(1,5)) * vee(1)
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum
    
    
    !..   3-mode coupling e-e-a2
    !      call Evcoup_eea2(dum,itd,p(pst(1,18)),pst(2,18),itemp)        !  e-e-a2 coupling
          dum=p(pst(1,18)) * pb(1) * vee(1)
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum*damp
    
    
    !..   3-mode coupling e-a1-a2
    !      call Evcoup_e1aa2(dum,itd,p(pst(1,20)),pst(2,20),itemp)        !  e-a1-a2 coupling (asym. stretch)
          dum=(p(pst(1,20)) * pa(1) * pb(1) * va(1))
    !       write(out_unit,*)'dum low',dum
    !      eref = eref + dum BUG detected 30.09.2015
          eref = eref + dum*damp
    !      call Evcoup_e2aa2(dum,itd,p(pst(1,21)),pst(2,21),itemp)        !  e-a1-a2 coupling (asym. bend)
          dum=(p(pst(1,21)) * pa(1) * pb(1) * vb(1))
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum*damp
    
    
    !..   3-mode coupling e-e-a1
    !      call Evcoup_eea1(dum,itd,p(pst(1,24)),pst(2,24),itemp)        !  e-e-a1 coupling
          dum=p(pst(1,24)) * pa(1) * vee(1)
    !       write(out_unit,*)'dum low',dum
          eref = eref + dum
    !
    !..   dynamic damping depending on lowest reference state energy:
    !     we use switching function with origin at 3 eV above NO2 + O asymptote and
    !     damp out the correction matrix within +/- 1 eV
    
    ! low order diagonal part
          eloworder(1,1)=eref
          eloworder(2,2)=eref
    
    !..   add JT coupling for 1. e mode
    !      call Ewjt6a(dum,itd, p(pst(1,8)), pst(2,8), itemp)
          dum=p(pst(1,8)) * wa(1)
          w = dum
    
    !      call Ezjt6a(dum,itd, p(pst(1,8)), pst(2,8), itemp)
          dum=p(pst(1,8)) * za(1)
          z = dum
    
    !..   add JT coupling for 2. e mode
    !      call Ewjt6b(dum, itd, p(pst(1, 9)), pst(2, 9), itemp)
          dum=p(pst(1,9)) * wb(1)
          w=w+dum
    !      call Ezjt6b(dum, itd, p(pst(1, 9)), pst(2, 9), itemp)
          dum=p(pst(1,9)) * zb(1)
          z=z+dum
    
    !..   add JT mode-mode coupling between two e modes:
    !      call Ewcoup_ee(dum,itd,p(pst(1,10)), pst(2,10), itemp)
          dum=p(pst(1,10)) * wee(1)
          w=w+dum
    !      call Ezcoup_ee(dum,itd,p(pst(1,10)), pst(2,10), itemp)
          dum=p(pst(1,10)) * zee(1)
          z=z+dum
    !..   add JT mode-mode coupling between a1 and e modes:
    !      call Ewcoup_e1a(dum,itd, p(pst(1,11)), pst(2,11), itemp)
          dum=p(pst(1,11))*pa(1)*wa(1)
          w=w+dum
    !      call Ezcoup_e1a(dum,itd, p(pst(1,11)), pst(2,11), itemp)
          dum=p(pst(1,11))*pa(1)*za(1)
          z=z+dum
    
    
    !      call Ewcoup_e2a(dum,itd, p(pst(1,12)), pst(2,12), itemp)
          dum=p(pst(1,12))*pa(1)*wb(1)
          w=w+dum
    !      call Ezcoup_e2a(dum,itd, p(pst(1,12)), pst(2,12), itemp)
          dum=p(pst(1,12))*pa(1)*zb(1)
          z=z+dum
    
    !--------------------------------------------------------------------------------
    
    !..   add JT mode-mode coupling between a2 and e modes:
    !      damp=exp(-abs(p(pst(1,26))*r(1))) !already computed
    !      call Ewcoup_e1a2(dum,itd,p(pst(1,15)),pst(2,15), itemp)
          dum=p(pst(1,15)) * pb(1) * wa(1)
          w=w+dum*damp
    !      call Ezcoup_e1a2(dum,itd,p(pst(1,15)),pst(2,15), itemp)
          dum=p(pst(1,15)) * pb(1) * za(1)
          z=z+dum*damp
    
    !      call Ewcoup_e2a2(dum,itd,p(pst(1,16)),pst(2,16), itemp)
          dum=p(pst(1,16)) * pb(1) * wb(1)
          w=w+dum*damp
    !      call Ezcoup_e2a2(dum,itd,p(pst(1,16)),pst(2,16), itemp)
          dum=p(pst(1,16)) * pb(1) * zb(1)
          z=z+dum*damp
    
    
    !..   add 3-mode coupling between e-e-a2 modes:
    !      call Ewcoup_eea2(dum,itd,p(pst(1,19)),pst(2,19), itemp)
          dum=p(pst(1,19)) * pb(1) * wee(1)
          w=w+dum*damp
    !      call Ezcoup_eea2(dum,itd,p(pst(1,19)),pst(2,19), itemp)
          dum=p(pst(1,19)) * pb(1) * zee(1)
          z=z+dum*damp
    
    !..   add 3-mode coupling between e-a1-a2 modes:
    !      call Ewcoup_e1aa2(dum,itd,p(pst(1,22)),pst(2,22), itemp)
          dum=p(pst(1,22)) * pa(1) * pb(1) *wa(1)
          w=w+dum*damp
    !      call Ezcoup_e1aa2(dum,itd,p(pst(1,22)),pst(2,22), itemp)
          dum=p(pst(1,22)) * pa(1) * pb(1) *za(1)
          z=z+dum*damp
    !      call Ewcoup_e2aa2(dum,itd,p(pst(1,23)),pst(2,23), itemp)
          dum=p(pst(1,23)) * pa(1) * pb(1) *wb(1)
          w=w+dum*damp
    !      call Ezcoup_e2aa2(dum,itd,p(pst(1,23)),pst(2,23), itemp)
          dum=p(pst(1,23)) * pa(1) * pb(1) *zb(1)
          z=z+dum*damp
    
    !ICI      endif PLANAR logical
    !..   add 3-mode coupling between e-e-a1 modes:
    !      call Ewcoup_eea1(dum,itd,p(pst(1,25)),pst(2,25), itemp)
          dum=p(pst(1,25)) * pa(1) * wee(1)
          w=w+dum     !*damp
    !      call Ezcoup_eea1(dum,itd,p(pst(1,25)),pst(2,25), itemp)
          dum=p(pst(1,25)) * pa(1) * zee(1)
          z=z+dum     !*damp
    
          eloworder(1,1)=eloworder(1,1) + w
          eloworder(2,2)=eloworder(2,2) - w
    
          eloworder(1,2)= z
          eloworder(2,1)= z
    
          CALL diagonalization(eloworder,vloworder,uloworder)
    
    !.. damping high order
    ! from Eisfeld 09.03.2015
           dum=(vloworder(1)-e0ref)/e0rho        !shift energy to onset
           if (dum < ZERO) then
             damplow=ONE
           elseif (dum > ONE) then
             damplow=ZERO
           else
            damplow=ONE - dum**2 * (TWO - dum**2)
           endif
    
    !
          IREF=1    ! Check if used in potential_functions
    ! corrections to low order matrix
    !--------------------------------------------------
    ! common part of e(1,1) and e(2,2) : e and w
    !--------------------------------------------------
    
    ! symmetric stretching
    !..   a1 symmetric stretch:
          call QML_NO3_Epolynom(dum,itd, QModel%no3minus_pot%p(pst(1,1):), pst(2,1), iref,avv)
          e=dum
    
    ! umbrella mode
    !..   a2 umbrella:
          call QML_NO3_Epolynom2(dum,itd, QModel%no3minus_pot%p(pst(1,2):), pst(2,2), iref,avv)
          e=e+dum*damp
    
    ! e stretching mode:
    !..   e stretching/bending
          CALL QML_NO3_Evjt6a(dum,itd, QModel%no3minus_pot%p(pst(1,3):), pst(2,3), iref,avv)
          e = e + dum
          call QML_NO3_Evjt6b(dum,itd, QModel%no3minus_pot%p(pst(1,4):), pst(2,4), iref,avv)
          e = e + dum
    
    !..   mode-mode coupling of a1-e
          call QML_NO3_Evcoup_e1a(dum,itd, QModel%no3minus_pot%p(pst(1,6):), pst(2,6), iref,avv)  ! e-a coupling 1. e mode
          e = e + dum
          call QML_NO3_Evcoup_e2a(dum,itd, QModel%no3minus_pot%p(pst(1,7):), pst(2,7), iref,avv)  ! e-a coupling 2. e mode
          e = e + dum
    
    !..   mode-mode coupling of a2-e
          call QML_NO3_Evcoup_e1a2(dum,itd, QModel%no3minus_pot%p(pst(1,13):), pst(2,13), iref,avv) !  coupling with asym. str.
          e = e + dum*damp
          call QML_NO3_Evcoup_e2a2(dum,itd, QModel%no3minus_pot%p(pst(1,14):), pst(2,14), iref,avv) !  coupling with asym. bend
          e = e + dum*damp
    
    !..   mode-mode coupling of a1-a2
          call QML_NO3_Evcoup_aa(dum,itd,QModel%no3minus_pot%p(pst(1,17):),pst(2,17), iref,avv)    !  coupling of a1 and a2
          e = e + dum*damp
    
    !..   mode-mode coupling e-e
          call QML_NO3_Evcoup_ee(dum,itd, QModel%no3minus_pot%p(pst(1,5):), pst(2,5), iref,avv)        ! e-e coupling
          e = e + dum
    
    !..   3-mode coupling e-e-a2
          call QML_NO3_Evcoup_eea2(dum,itd,QModel%no3minus_pot%p(pst(1,18):),pst(2,18),iref,avv)        !  e-e-a2 coupling
          e = e + dum*damp
    
    
    !..   3-mode coupling e-a1-a2
          call QML_NO3_Evcoup_e1aa2(dum,itd,QModel%no3minus_pot%p(pst(1,20):),pst(2,20),iref,avv)  !  e-a1-a2 coupling (asym. stretch)
          e = e + dum*damp
    !      e = e + dum bug detected on 30.09.2015
          call QML_NO3_Evcoup_e2aa2(dum,itd,QModel%no3minus_pot%p(pst(1,21):),pst(2,21),iref,avv)  !  e-a1-a2 coupling (asym. bend)
          e = e + dum*damp
    
    
    !..   3-mode coupling e-e-a1
          call QML_NO3_Evcoup_eea1(dum,itd,QModel%no3minus_pot%p(pst(1,24):),pst(2,24),iref,avv)        !  e-e-a1 coupling
          e = e + dum
    
          w=ZERO
    !..   add JT coupling for 1. e mode
          call QML_NO3_Ewjt6a(dum,itd, QModel%no3minus_pot%p(pst(1,8):), pst(2,8), iref,avv)
          w = dum
    
    !..   add JT coupling for 2. e mode
          call QML_NO3_Ewjt6b(dum, itd, QModel%no3minus_pot%p(pst(1, 9):), pst(2, 9), iref,avv)
          w=w+dum
    
    !..   add JT mode-mode coupling between two e modes:
          call QML_NO3_Ewcoup_ee(dum,itd,QModel%no3minus_pot%p(pst(1,10):), pst(2,10), iref,avv)
          w=w+dum
    
    !..   add JT mode-mode coupling between a1 and e modes:
          call QML_NO3_Ewcoup_e1a(dum,itd, QModel%no3minus_pot%p(pst(1,11):), pst(2,11), iref,avv)
          w=w+dum
          call QML_NO3_Ewcoup_e2a(dum,itd, QModel%no3minus_pot%p(pst(1,12):), pst(2,12), iref,avv)
          w=w+dum
    
    !..   add JT mode-mode coupling between a2 and e modes:
          call QML_NO3_Ewcoup_e1a2(dum,itd,QModel%no3minus_pot%p(pst(1,15):),pst(2,15), iref,avv)
          w=w+dum*damp
          call QML_NO3_Ewcoup_e2a2(dum,itd,QModel%no3minus_pot%p(pst(1,16):),pst(2,16), iref,avv)
          w=w+dum*damp
    
    !..   add 3-mode coupling between e-e-a2 modes:
          call QML_NO3_Ewcoup_eea2(dum,itd,QModel%no3minus_pot%p(pst(1,19):),pst(2,19), iref,avv)
          w=w+dum*damp
    
    !..   add 3-mode coupling between e-a1-a2 modes:
          call QML_NO3_Ewcoup_e1aa2(dum,itd,QModel%no3minus_pot%p(pst(1,22):),pst(2,22), iref,avv)
          w=w+dum*damp
          call QML_NO3_Ewcoup_e2aa2(dum,itd,QModel%no3minus_pot%p(pst(1,23):),pst(2,23), iref,avv)
          w=w+dum*damp
    
    !..   add 3-mode coupling between e-e-a1 modes:
          call QML_NO3_Ewcoup_eea1(dum,itd,QModel%no3minus_pot%p(pst(1,25):),pst(2,25), iref,avv)
          w=w+dum     !*damp
    !
    !--------------------------------------------------
    ! off diagonal term  e(1,2) : z
          z=ZERO
    !..   add JT coupling for 1. e mode
          call QML_NO3_Ezjt6a(dum,itd, QModel%no3minus_pot%p(pst(1,8):), pst(2,8), iref,avv)
          z = dum
    
    !..   add JT coupling for 2. e mode
          call QML_NO3_Ezjt6b(dum, itd, QModel%no3minus_pot%p(pst(1, 9):), pst(2, 9), iref,avv)
          z=z+dum
    
    !..   add JT mode-mode coupling between two e modes:
          call QML_NO3_Ezcoup_ee(dum,itd,QModel%no3minus_pot%p(pst(1,10):), pst(2,10), iref,avv)
          z=z+dum
    
    !..   add JT mode-mode coupling between a1 and e modes:
          call QML_NO3_Ezcoup_e1a(dum,itd, QModel%no3minus_pot%p(pst(1,11):), pst(2,11), iref,avv)
          z=z+dum
          call QML_NO3_Ezcoup_e2a(dum,itd, QModel%no3minus_pot%p(pst(1,12):), pst(2,12), iref,avv)
          z=z+dum
    
    !..   add JT mode-mode coupling between a2 and e modes:
          call QML_NO3_Ezcoup_e1a2(dum,itd,QModel%no3minus_pot%p(pst(1,15):),pst(2,15), iref,avv)
          z=z+dum*damp
          call QML_NO3_Ezcoup_e2a2(dum,itd,QModel%no3minus_pot%p(pst(1,16):),pst(2,16), iref,avv)
          z=z+dum*damp
    
    !..   add 3-mode coupling between e-e-a2 modes:
          call QML_NO3_Ezcoup_eea2(dum,itd,QModel%no3minus_pot%p(pst(1,19):),pst(2,19), iref,avv)
          z=z+dum*damp
    
    !..   add 3-mode coupling between e-a1-a2 modes:
          call QML_NO3_Ezcoup_e1aa2(dum,itd,QModel%no3minus_pot%p(pst(1,22):),pst(2,22), iref,avv)
          z=z+dum*damp
          call QML_NO3_Ezcoup_e2aa2(dum,itd,QModel%no3minus_pot%p(pst(1,23):),pst(2,23), iref,avv)
          z=z+dum*damp
    
    !..   add 3-mode coupling between e-e-a1 modes:
          call QML_NO3_Ezcoup_eea1(dum,itd,QModel%no3minus_pot%p(pst(1,25):),pst(2,25), iref,avv)
          z=z+dum     !*damp
    
    !
    ! update v
          vdia(1,1)=eloworder(1,1)+damplow*(e+w)
          vdia(2,2)=eloworder(2,2)+damplow*(e-w)
          vdia(1,2)=eloworder(1,2)+damplow*z
          vdia(2,1)=eloworder(1,2)+damplow*z
  End Subroutine QML_NO3_potential
    !
    !------------------------------------------------------------------------------
    ! returns exponent of tunable Morse coordinate
    ! exponent is polynomial * gaussian (skewed)
    ! ilabel = 1 or 2 selects the parameters a and sfac to be used
    !      function f(x,ilabel)
    !      function f(x,ii)
  subroutine QML_NO3_ff(x,ii,f,tmc)
    implicit none
    TYPE (tmc_t),                 intent(in)         :: tmc
    real (kind=Rkind),            intent(in)         :: x
    integer,                      intent(in)         :: ii
    real (kind=Rkind),            intent(inout)      :: f

    integer i, pn, nmax
    real (kind=Rkind) :: r, gaus, poly, skew

    
    ! a(1):        position of equilibrium
    ! a(2):        constant of exponent
    ! a(3):        constant for skewing the gaussian
    ! a(4):        tuning for skewing the gaussian
    ! a(5):        Gaussian exponent
    ! a(6):        Shift of Gaussian maximum
    ! a(7)...:     polynomial coefficients
    ! a(8+n)...:   coefficients of Morse Power series
    
    ! Tunable Morse function
    ! Power series in Tunable Morse coordinates of order m
    ! exponent is polynomial of order npoly * gaussian + switching function
    
    !       write(out_unit,*) 'x:', x
    !       write(out_unit,'(10f12.6)') (a(i), i=1,6)
    !       write(out_unit,'(10f12.6)') (a(i), i=7,11)
    !       write(out_unit,'(10f12.6)') (a(i), i=12,14)
    !      stop
    
    !.....set r  r-r_e
          r=x-tmc%a(1)
    
          !a(2)=abs(a(2)) ! done in the initialization
          !a(5)=abs(a(5))
    
    !..   set up skewing function:
    !      skew=HALF * a(ii,3) *(tanh(a(ii,4)*r) + ONE)
          skew=HALF * tmc%a(3) *(tanh(tmc%a(4)*(r-tmc%a(6))) + ONE)
    
    !..   set up gaussian function:
    !      gaus=exp(-abs(a(ii,5))*(r+a(ii,6))**2)
          gaus=exp(-tmc%a(5)*r**2)
    
    !..   set up power series:
          poly=ZERO
          do i=0,tmc%npoly(ii)-1
            poly=poly+tmc%a(7+i)*r**i
          enddo
    
    !..   set up full exponent function:
          f=tmc%a(2)  + skew + gaus*poly
    !       write(out_unit,*) ii, a(ii,2)
    
    end subroutine QML_NO3_ff
    
    !---------------------------------------------------------------------------------
    !     function to generate polynomial of any order
    
    !      function polynom(n, p, j, iref,avv)
  subroutine QML_NO3_Epolynom(polynom,n, p, j, iref,avv)
    implicit none
    real (kind=Rkind),            intent(inout) :: polynom
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    polynom = sum( avv%pa(3:j) * p(3:j) )
    
  end subroutine QML_NO3_Epolynom
    
    !---------------------------------------------------------------------------------
    !     function to generate polynomial of even power
    
    !      function polynom2(n, p, j, iref,avv)
  subroutine QML_NO3_Epolynom2(polynom2,n, p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: polynom2
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref

    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    polynom2 = sum( avv%pb(2:j) * p(2:j) )
    
  end subroutine QML_NO3_Epolynom2
    
  !---------------------------------------------------------------------------------
  !     function to generate V Jahn-Teller matrix elements up to 6th order
  subroutine QML_NO3_Evjt6a(vjt6a,n, p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: vjt6a
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv
 
    vjt6a = sum(p(2:j)*avv%va(2:j))

  end subroutine QML_NO3_Evjt6a
  !---------------------------------------------------------------------------------
  !     function to generate V Jahn-Teller matrix elements up to 6th order
  subroutine QML_NO3_Evjt6b(vjt6b, n, p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: vjt6b
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    vjt6b = sum(p(2:j)*avv%vb(2:j))

  end subroutine QML_NO3_Evjt6b
    
  !---------------------------------------------------------------------------
  !     function to generate W Jahn-Teller matrix elements up to 6th order
  subroutine QML_NO3_Ewjt6a(wjt6a, n, p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: wjt6a
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    wjt6a = sum(p(2:j)*avv%wa(2:j))

  end subroutine QML_NO3_Ewjt6a
    
  !---------------------------------------------------------------------------
  !     function to generate W Jahn-Teller matrix elements up to 6th order
  subroutine QML_NO3_Ewjt6b(wjt6b, n, p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: wjt6b
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    wjt6b = sum(p(2:j)*avv%wb(2:j))

  end subroutine QML_NO3_Ewjt6b

  !---------------------------------------------------------------------------
  !     function to generate Z Jahn-Teller matrix elements up to 6th order
  subroutine QML_NO3_Ezjt6a(zjt6a, n, p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: zjt6a
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    zjt6a = sum(p(2:j)*avv%za(2:j))

  end subroutine QML_NO3_Ezjt6a
    
  !---------------------------------------------------------------------------
  !     function to generate Z Jahn-Teller matrix elements up to 6th order
  subroutine QML_NO3_Ezjt6b(zjt6b, n, p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: zjt6b
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    zjt6b = sum(p(2:j)*avv%zb(2:j))

  end subroutine QML_NO3_Ezjt6b

  !------------------------------------------------------------------------------
  !      function vcoup_ee(n, par, j, iref,avv)
  subroutine QML_NO3_Evcoup_ee(vcoup_ee, n, p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: vcoup_ee
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    vcoup_ee = sum(p(2:j)*avv%vee(2:j))

  end subroutine QML_NO3_Evcoup_ee
    
  !------------------------------------------------------------------------------
  !      function wcoup_ee(n,par, j, iref,avv)
  subroutine QML_NO3_Ewcoup_ee(wcoup_ee, n,p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: wcoup_ee
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    wcoup_ee = sum(p(2:j)*avv%wee(2:j))

  end subroutine QML_NO3_Ewcoup_ee
    
  !------------------------------------------------------------------------------
  !      function zcoup_ee(n,par, j, iref,avv)
  subroutine QML_NO3_Ezcoup_ee(zcoup_ee, n,p, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout) :: zcoup_ee
    real (kind=Rkind),            intent(in)    :: p(:)
    integer,                      intent(in)    :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in)    :: avv

    zcoup_ee = sum(p(2:j)*avv%zee(2:j))

  end subroutine QML_NO3_Ezcoup_ee
    
  !------------------------------------------------------------------------------
  !     function to generate V matrix elements for the coupling between e and a modes up to fourth order
  subroutine QML_NO3_Evcoup_e1a(vcoup_e1a, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_e1a
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=10
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),va(:)

    pa => avv%pa
    va => avv%va

    p(:)   = ZERO
    p(2:j) = par(2:j) 

    !     3rd order
    vcoup_e1a =  p(1)*pa(1)*va(1)
    
    !     4th order
    vcoup_e1a = vcoup_e1a + p(2)*pa(1)*va(2)+p(3)*pa(2)*va(1)
    
    !     5th order
    vcoup_e1a = vcoup_e1a + p(4)*pa(1)*va(3)+p(5)*pa(2)*va(2) + p(6)*pa(3)*va(1)
    
    !     6th order
    vcoup_e1a= vcoup_e1a + p(7)*pa(1)*va(4) + p(8)*pa(2)*va(3) + p(9)*pa(3)*va(2) + p(10)*pa(4)*va(1)
    
  end subroutine QML_NO3_Evcoup_e1a
    
  !------------------------------------------------------------------------------
  !     function to generate V matrix elements for the coupling between e and a modes up to fourth order
  subroutine QML_NO3_Evcoup_e2a(vcoup_e2a, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_e2a
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=10
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),vb(:)

    pa => avv%pa
    vb => avv%vb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     3rd order
    vcoup_e2a =  p(1)*pa(1)*vb(1)
    
    !     4th order
    vcoup_e2a = vcoup_e2a + p(2)*pa(1)*vb(2)+p(3)*pa(2)*vb(1)
    
    !     5th order
    vcoup_e2a=vcoup_e2a + p(4)*pa(1)*vb(3) + p(5)*pa(2)*vb(2) + p(6)*pa(3)*vb(1)
    
    !     6th order
    vcoup_e2a= vcoup_e2a + p(7)*pa(1)*vb(4) + p(8)*pa(2)*vb(3) + p(9)*pa(3)*vb(2) + p(10)*pa(4)*vb(1)
    
  end subroutine QML_NO3_Evcoup_e2a
    
  !-------------------------------------------------------------------------------------
  !     function to generate W matrix elements for the coupling between e and a modes up to sixth order    
  subroutine QML_NO3_Ewcoup_e1a(wcoup_e1a, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: wcoup_e1a
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=18
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),wa(:)

    pa => avv%pa
    wa => avv%wa

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     2nd order
    wcoup_e1a = p(1)*pa(1)*wa(1)
    
    !     3rd order
    wcoup_e1a = wcoup_e1a + p(2)*pa(1)*wa(2) + p(3)*pa(2)*wa(1)
    
    !     4th order
    wcoup_e1a = wcoup_e1a + p(4)*pa(1)*wa(3) + p(5)*pa(2)*wa(2) + p(6)*pa(3)*wa(1)
    
    !     5th order
    wcoup_e1a = wcoup_e1a + p(7)*pa(4)*wa(1) + p(8)*pa(3)*wa(2) +           &
          p(9)*pa(2)*wa(3) + p(10)*pa(1)*wa(4) + p(11)*pa(1)*wa(5)
    
    !     6th order
    wcoup_e1a = wcoup_e1a +                                                       &
         &    p(12)*pa(5)*wa(1)+p(13)*pa(4)*wa(2)                               &
         &    + p(14)*pa(3)*wa(3)                                                   &
         &    + p(15)*pa(2)*wa(4)+p(16)*pa(2)*wa(5)                             &
         &    + pa(1)*(p(17)*wa(6)+p(18)*wa(7))
    
  end subroutine QML_NO3_Ewcoup_e1a
    
  !-------------------------------------------------------------------------------------
  !     function to generate W matrix elements for the coupling between e and a modes up to sixth order    
  subroutine QML_NO3_Ewcoup_e2a(wcoup_e2a, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: wcoup_e2a
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=18
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),wb(:)

    pa => avv%pa
    wb => avv%wb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     2nd order
    wcoup_e2a = p(1)*pa(1)*wb(1)
    
    !     3rd order
    wcoup_e2a=wcoup_e2a +(p(2)*pa(1)*wb(2)+p(3)*pa(2)*wb(1))   !/ fact(3)
    
    !     4th order
    wcoup_e2a=wcoup_e2a +(p(4)*pa(1)*wb(3)+p(5)*pa(2)*wb(2)               &
         & + p(6)*pa(3)*wb(1))   !/ fact(4)
    
    !     5th order
    wcoup_e2a = wcoup_e2a                                                         &
         & + (p(7)*pa(4)*wb(1) + p(8)*pa(3)*wb(2)                               &
         & + p(9)*pa(2)*wb(3)+p(10)*pa(1)*wb(4)                                 &
         & + p(11)*pa(1)*wb(5))  !/ fact(5)
    
    !     6th order
    wcoup_e2a = wcoup_e2a +                                                       &
         &    p(12)*pa(5)*wb(1)+p(13)*pa(4)*wb(2)                               &
         &   +p(14)*pa(3)*wb(3)                                                     &
         &    + p(15)*pa(2)*wb(4)+p(16)*pa(2)*wb(5)                             &
         &    + pa(1)*(p(17)*wb(6)+p(18)*wb(7))
    
  end subroutine QML_NO3_Ewcoup_e2a
    
  !-------------------------------------------------------------------------------------
  !     function to generate Z matrix elements for the coupling betzeen e and a modes up to sixth order
  subroutine QML_NO3_Ezcoup_e1a(zcoup_e1a, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: zcoup_e1a
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=18
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),za(:)

    pa => avv%pa
    za => avv%za

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     2nd order
    zcoup_e1a = p(1)*pa(1)*za(1)
    
    !     3rd order
    zcoup_e1a=zcoup_e1a+ (p(2)*pa(1)*za(2) + p(3)*pa(2)*za(1))   !/ fact(3)
    
    !     4th order
    zcoup_e1a=zcoup_e1a + (p(4)*pa(1)*za(3) + p(5)*pa(2)*za(2)            &
         & + p(6)*pa(3)*za(1))   !/ fact(4)
    
    !     5th order
    zcoup_e1a=zcoup_e1a + (p(7)*pa(4)*za(1) + p(8)*pa(3)*za(2)            &
         & +p(9)*pa(2)*za(3)+p(10)*pa(1)*za(4)                                  &
         & +p(11)*pa(1)*za(5))  !/ fact(5)
    
    !     6th order
    zcoup_e1a = zcoup_e1a +                                                       &
         &    p(12)*pa(5)*za(1)+p(13)*pa(4)*za(2)                               &
         &   +p(14)*pa(3)*za(3)                                                     &
         &    + p(15)*pa(2)*za(4)+p(16)*pa(2)*za(5)                             &
         &    + pa(1)*(p(17)*za(6)+p(18)*za(7))
    
  end subroutine QML_NO3_Ezcoup_e1a

  !-------------------------------------------------------------------------------------
  !     function to generate Z matrix elements for the coupling betzeen e and a modes up to sixth order
  subroutine QML_NO3_Ezcoup_e2a(zcoup_e2a, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: zcoup_e2a
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=18
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),zb(:)

    pa => avv%pa
    zb => avv%zb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     2nd order
    zcoup_e2a = p(1)*pa(1)*zb(1)
    
    !     3rd order
    zcoup_e2a=zcoup_e2a +(p(2)*pa(1)*zb(2) + p(3)*pa(2)*zb(1))   !/ fact(3)
    
    !     4th order
    zcoup_e2a=zcoup_e2a + (p(4)*pa(1)*zb(3) + p(5)*pa(2)*zb(2) + p(6)*pa(3)*zb(1))   !/ fact(4)
    
    !     5th order
    zcoup_e2a=zcoup_e2a + (p(7)*pa(4)*zb(1) + p(8)*pa(3)*zb(2)            &
         & +p(9)*pa(2)*zb(3)+p(10)*pa(1)*zb(4) +p(11)*pa(1)*zb(5))  !/ fact(5)
    
    !     6th order
    zcoup_e2a = zcoup_e2a +                                                       &
         &    p(12)*pa(5)*zb(1)+p(13)*pa(4)*zb(2)                               &
         &   +p(14)*pa(3)*zb(3)                                                     &
         &    + p(15)*pa(2)*zb(4)+p(16)*pa(2)*zb(5)                             &
         &    + pa(1)*(p(17)*zb(6)+p(18)*zb(7))
    
  end subroutine QML_NO3_Ezcoup_e2a

  !------------------------------------------------------------------------------
  !     function to generate V matrix elements for the coupling between a1 and a2 modes up to sixth order
  !     only even orders in a2
  subroutine QML_NO3_Evcoup_aa(vcoup_aa, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_aa
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=6
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),pb(:)

    pa => avv%pa
    pb => avv%pb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !..   third order
    vcoup_aa = p(1) * pa(1) * pb(1)
    
    !..   fourth order
    vcoup_aa = vcoup_aa + p(2) * pa(2) * pb(1)
    
    !..   fifth order
    vcoup_aa = vcoup_aa + p(3) * pa(1) * pb(2) + p(4) * pa(3) * pb(1)
    
    !..   sixths order
    vcoup_aa = vcoup_aa + p(5) * pa(2) * pb(2) + p(6) * pa(4) * pb(1)
    
  end subroutine QML_NO3_Evcoup_aa

  !------------------------------------------------------------------------------
  !     function to generate V matrix elements for the coupling between e and a modes up to sixth order
  !     only even orders in a
  subroutine QML_NO3_Evcoup_e1a2(vcoup_e1a2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_e1a2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=9
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: va(:),pb(:)

    va => avv%va
    pb => avv%pb

    p(:)   = ZERO
    p(2:j) = par(2:j)

 
    !     4th order
    vcoup_e1a2 = p(1)*pb(1)*va(1)
    
    !     5th order
    vcoup_e1a2 = vcoup_e1a2 + p(2)*pb(1)*va(2) 
    
    !     6th order
    vcoup_e1a2 = vcoup_e1a2 + p(3)*pb(1)*va(3) + p(4)*pb(2)*va(1)
    
    !     7th order
    vcoup_e1a2 = vcoup_e1a2 + p(5)*pb(1)*va(4) + p(6)*pb(2)*va(2)
    
    !     8th order
    vcoup_e1a2 = vcoup_e1a2 + p(7)*pb(1)*va(5) + p(8)*pb(1)*va(6) + p(9)*pb(2)*va(3)
    
  end subroutine QML_NO3_Evcoup_e1a2
    
  !------------------------------------------------------------------------------
  !     function to generate V matrix elements for the coupling between e and a modes up to sixth order
  !     only even orders in a
  subroutine QML_NO3_Evcoup_e2a2(vcoup_e2a2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_e2a2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=9
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: vb(:),pb(:)

    vb => avv%vb
    pb => avv%pb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     4th order
          vcoup_e2a2 = p(1)*pb(1)*vb(1)
    
    !     5th order
          vcoup_e2a2 = vcoup_e2a2 + p(2)*pb(1)*vb(2)                                &
         & ! /fact(5)
    
    !     6th order
          vcoup_e2a2 = vcoup_e2a2 + p(3)*pb(1)*vb(3)                                &
         &                    + p(4)*pb(2)*vb(1)
    
    !     7th order
          vcoup_e2a2 = vcoup_e2a2 + p(5)*pb(1)*vb(4)                                &
         &                      + p(6)*pb(2)*vb(2)
    
    !     8th order
          vcoup_e2a2 = vcoup_e2a2 + p(7)*pb(1)*vb(5)                                &
         &                      + p(8)*pb(1)*vb(6)                                  &
         &                      + p(9)*pb(2)*vb(3)
    
  end subroutine QML_NO3_Evcoup_e2a2
    
  !-------------------------------------------------------------------------------------
  !     function to generate W matrix elements for the coupling between e and a modes up to sixth order
  !     only even orders in a
  subroutine QML_NO3_Ewcoup_e1a2(wcoup_e1a2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: wcoup_e1a2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=14
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: wa(:),pb(:)

    wa => avv%wa
    pb => avv%pb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     3rd order
          wcoup_e1a2 = p(1) * pb(1) * wa(1)
    
    !     4th order
          wcoup_e1a2 = wcoup_e1a2 + p(2) * pb(1) * wa(2)
    
    !     5th order
          wcoup_e1a2 = wcoup_e1a2 + p(3) * pb(2) * wa(1)                            &
         &                    + p(4) * pb(1) * wa(3)
    
    !     6th order
          wcoup_e1a2 = wcoup_e1a2 + p(5) * pb(2) * wa(2)                            &
         &                    + p(6) * pb(1) * wa(4)                                &
         &                    + p(7) * pb(1) * wa(5)
    
    !     7th order
          wcoup_e1a2 = wcoup_e1a2 + p(8) * pb(1) * wa(6)                            &
         &                      + p(9) * pb(1) * wa(7)                              &
         &                      + p(10)* pb(2) * wa(3)
    
    !     8th order
          wcoup_e1a2 = wcoup_e1a2 + p(11) * pb(1) * wa(8)                           &
         &                      + p(12) * pb(1) * wa(9)                             &
         &                      + p(13) * pb(2) * wa(4)                             &
         &                      + p(14) * pb(2) * wa(5)
    
  end subroutine QML_NO3_Ewcoup_e1a2
    
  !-------------------------------------------------------------------------------------
  !     function to generate W matrix elements for the coupling between e and a modes up to sixth order
  !     only even orders in a
  subroutine QML_NO3_Ewcoup_e2a2(wcoup_e2a2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: wcoup_e2a2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=14
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: wb(:),pb(:)

    wb => avv%wb
    pb => avv%pb

    p(:)   = ZERO
    p(2:j) = par(2:j)

 
    !     3rd order
          wcoup_e2a2 = p(1) * pb(1) * wb(1)
    
    !     4th order
          wcoup_e2a2 = wcoup_e2a2 + p(2) * pb(1) * wb(2)
    
    !     5th order
          wcoup_e2a2 = wcoup_e2a2 + p(3) * pb(2) * wb(1)                            &
         &                    + p(4) * pb(1) * wb(3)
    
    !     6th order
          wcoup_e2a2 = wcoup_e2a2 + p(5) * pb(2) * wb(2)                            &
         &                    + p(6) * pb(1) * wb(4)                                &
         &                    + p(7) * pb(1) * wb(5)
    
    !     7th order
          wcoup_e2a2 = wcoup_e2a2 + p(8) * pb(1) * wb(6)                            &
         &                      + p(9) * pb(1) * wb(7)                              &
         &                      + p(10)* pb(2) * wb(3)
    
    !     8th order
          wcoup_e2a2 = wcoup_e2a2 + p(11) * pb(1) * wb(8)                           &
         &                      + p(12) * pb(1) * wb(9)                             &
         &                      + p(13) * pb(2) * wb(4)                             &
         &                      + p(14) * pb(2) * wb(5)
    
  end subroutine QML_NO3_Ewcoup_e2a2
    
  !-------------------------------------------------------------------------------------
  !     function to generate Z matrix elements for the coupling betzeen e and a modes up to sixth order
  !     only even orders in a
  subroutine QML_NO3_Ezcoup_e1a2(zcoup_e1a2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: zcoup_e1a2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=14
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pb(:),za(:)

    pb => avv%pb
    za => avv%za

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     3rd order
          zcoup_e1a2 = p(1) * pb(1) * za(1)
    
    !     4th order
          zcoup_e1a2 = zcoup_e1a2 + p(2) * pb(1) * za(2)
    
    !     5th order
          zcoup_e1a2 = zcoup_e1a2 + p(3) * pb(2) * za(1)                            &
         &                    + p(4) * pb(1) * za(3)
    
    !     6th order
          zcoup_e1a2 = zcoup_e1a2 + p(5) * pb(2) * za(2)                            &
         &                    + p(6) * pb(1) * za(4)                                &
         &                    + p(7) * pb(1) * za(5)
    
    !     7th order
          zcoup_e1a2 = zcoup_e1a2 + p(8) * pb(1) * za(6)                            &
         &                      + p(9) * pb(1) * za(7)                              &
         &                      + p(10)* pb(2) * za(3)
    
    !     8th order
          zcoup_e1a2 = zcoup_e1a2 + p(11) * pb(1) * za(8)                           &
         &                      + p(12) * pb(1) * za(9)                             &
         &                      + p(13) * pb(2) * za(4)                             &
         &                      + p(14) * pb(2) * za(5)
    
  end subroutine QML_NO3_Ezcoup_e1a2
    
  !-------------------------------------------------------------------------------------
  !     function to generate Z matrix elements for the coupling betzeen e and a modes up to sixth order
  !     only even orders in a
  subroutine QML_NO3_Ezcoup_e2a2(zcoup_e2a2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: zcoup_e2a2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=14
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pb(:),zb(:)

    pb => avv%pb
    zb => avv%zb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     3rd order
          zcoup_e2a2 = p(1) * pb(1) * zb(1)
    
    !     4th order
          zcoup_e2a2 = zcoup_e2a2 + p(2) * pb(1) * zb(2)
    
    !     5th order
          zcoup_e2a2 = zcoup_e2a2 + p(3) * pb(2) * zb(1)                            &
         &                    + p(4) * pb(1) * zb(3)
    
    !     6th order
          zcoup_e2a2 = zcoup_e2a2 + p(5) * pb(2) * zb(2)                            &
         &                    + p(6) * pb(1) * zb(4)                                &
         &                    + p(7) * pb(1) * zb(5)
    
    !     7th order
          zcoup_e2a2 = zcoup_e2a2 + p(8) * pb(1) * zb(6)                            &
         &                      + p(9) * pb(1) * zb(7)                              &
         &                      + p(10)* pb(2) * zb(3)
    
    !     8th order
          zcoup_e2a2 = zcoup_e2a2 + p(11) * pb(1) * zb(8)                           &
         &                      + p(12) * pb(1) * zb(9)                             &
         &                      + p(13) * pb(2) * zb(4)                             &
         &                      + p(14) * pb(2) * zb(5)
    
  end subroutine QML_NO3_Ezcoup_e2a2

  !------------------------------------------------------------------------------
  !  3-mode coupling of a1 + e + e
  !      function vcoup_eea1(n, par, j, iref,avv)
  subroutine QML_NO3_Evcoup_eea1(vcoup_eea1, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_eea1
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=26
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),vee(:)

    pa  => avv%pa
    vee => avv%vee

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     3rd order
          vcoup_eea1 = p(1) * pa(1) * vee(1)
    
    !     4th order
          vcoup_eea1=vcoup_eea1 + pa(1) * (p(2)*vee(2) + p(3)*vee(3))
          vcoup_eea1=vcoup_eea1 + p(4) * pa(2) * vee(1)
    
    !     5th order
          vcoup_eea1 = vcoup_eea1 + pa(1)*(p(5)*vee(4) + p(6)*vee(5)              &
         &                        + p(7)*vee(6) + p(8)*vee(7))                      &
         & +pa(2)*(p(9)*vee(2)+ p(10)*vee(3)) + pa(3)*p(11)*vee(1)
    
    !     6th order
          vcoup_eea1=vcoup_eea1 + pa(2) * (p(12)*vee(4) + p(13)*vee(5)            &
         &                        + p(14)*vee(6) + p(15)*vee(7))                    &
         & +pa(3)*(p(16)*vee(2)+p(17)*vee(3)) + pa(4)*p(18)*vee(1)
    
    !................................................................................
    ! order not logical because added later!
          vcoup_eea1 = vcoup_eea1 + pa(1)*(p(19)*vee(8)+p(20)*vee(9)              &
         & +p(21)*vee(10) +p(22)*vee(11)+p(23)*vee(12)+p(24)*vee(13)            &
         & +p(25)*vee(14) +p(26)*vee(15) )
    
  end subroutine QML_NO3_Evcoup_eea1
    
  !------------------------------------------------------------------------------
  !  3-mode coupling
  !      function wcoup_eea1(n,par, j, iref,avv)
  subroutine QML_NO3_Ewcoup_eea1(wcoup_eea1, n,par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: wcoup_eea1
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=46
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),wee(:)

    pa  => avv%pa
    wee => avv%wee

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     3rd order
          wcoup_eea1 = p(1) * pa(1) * wee(1)
    
    !     4th order
          wcoup_eea1 = wcoup_eea1 + pa(1) * (p(2)*wee(2)+p(3)*wee(3)              &
         &                        + p(4)*wee(4) + p(5)*wee(5))
          wcoup_eea1 = wcoup_eea1 + pa(2) * p(6)*wee(1)
    
    !     5th order
          wcoup_eea1 = wcoup_eea1 + pa(1) * (p(7)*wee(6)+p(8)*wee(7)              &
         &             + p(9)*wee(8) + p(10)*wee(9) + p(11)*wee(10)               &
         &             + p(12)*wee(11) + p(13)*wee(12) + p(14)*wee(13)            &
         &             + p(15)*wee(14))
          wcoup_eea1 = wcoup_eea1 + pa(2) * (p(16)*wee(2)+p(17)*wee(3)            &
         & + p(18)*wee(4)  + p(19)*wee(5))
          wcoup_eea1 = wcoup_eea1 + pa(3) * p(20)*wee(1)
    
    !     6th order
          wcoup_eea1 = wcoup_eea1 + pa(2)*(p(21)*wee(6) + p(22)*wee(7)            &
         &             + p(23)*wee(8) + p(24)*wee(9) + p(25)*wee(10)              &
         &             + p(26)*wee(11) + p(27)*wee(12) + p(28)*wee(13)            &
         &             + p(29)*wee(14))
          wcoup_eea1 = wcoup_eea1 + pa(3) * (p(30)*wee(2)+p(31)*wee(3)            &
         &                        + p(32)*wee(4) + p(33)*wee(5))
          wcoup_eea1 = wcoup_eea1 + pa(4) * p(34)*wee(1)
    
    !................................................................................
    ! order not logical because added later!
          wcoup_eea1 = wcoup_eea1 + pa(1)*(p(35)*wee(15)+p(36)*wee(16)            &
         &             + p(37)*wee(17) + p(38)*wee(18) + p(39)*wee(19)            &
         &             + p(40)*wee(20) + p(41)*wee(21) + p(42)*wee(22)            &
         &             + p(43)*wee(23) + p(44)*wee(24)+ p(45)*wee(25)             &
         &             + p(46)*wee(26) )
    
  end subroutine QML_NO3_Ewcoup_eea1
    
    !------------------------------------------------------------------------------
    !  3-mode coupling
    !      function zcoup_eea1(n,par, j, iref,avv)
  subroutine QML_NO3_Ezcoup_eea1(zcoup_eea1, n,par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: zcoup_eea1
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=46
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),zee(:)

    pa  => avv%pa
    zee => avv%zee

    p(:)   = ZERO
    p(2:j) = par(2:j)
 
    !     3rd order
          zcoup_eea1 = p(1) * pa(1) * zee(1)
    
    !     4th order
          zcoup_eea1 = zcoup_eea1 + pa(1) * (p(2)*zee(2)+p(3)*zee(3)              &
         &                        + p(4)*zee(4) + p(5)*zee(5))
          zcoup_eea1 = zcoup_eea1 + pa(2) * p(6)*zee(1)
    
    !     5th order
          zcoup_eea1 = zcoup_eea1 + pa(1) * (p(7)*zee(6)+p(8)*zee(7)              &
         &             + p(9)*zee(8) + p(10)*zee(9) + p(11)*zee(10)               &
         &             + p(12)*zee(11) + p(13)*zee(12) + p(14)*zee(13)            &
         &             + p(15)*zee(14))
          zcoup_eea1 = zcoup_eea1 + pa(2) * (p(16)*zee(2)+p(17)*zee(3)            &
         & + p(18)*zee(4)  + p(19)*zee(5))
          zcoup_eea1 = zcoup_eea1 + pa(3) * p(20)*zee(1)
    
    !     6th order
          zcoup_eea1 = zcoup_eea1 + pa(2)*(p(21)*zee(6) + p(22)*zee(7)            &
         &             + p(23)*zee(8) + p(24)*zee(9) + p(25)*zee(10)              &
         &             + p(26)*zee(11) + p(27)*zee(12) + p(28)*zee(13)            &
         &             + p(29)*zee(14))
          zcoup_eea1 = zcoup_eea1 + pa(3) * (p(30)*zee(2)+p(31)*zee(3)            &
         &                        + p(32)*zee(4) + p(33)*zee(5))
          zcoup_eea1 = zcoup_eea1 + pa(4) * p(34)*zee(1)
    
    !................................................................................
    ! order not logical because added later!
          zcoup_eea1 = zcoup_eea1 + pa(1)*(p(35)*zee(15)+p(36)*zee(16)            &
         &             + p(37)*zee(17) + p(38)*zee(18) + p(39)*zee(19)            &
         &             + p(40)*zee(20) + p(41)*zee(21) + p(42)*zee(22)            &
         &             + p(43)*zee(23) + p(44)*zee(24)+ p(45)*zee(25)             &
         &             + p(46)*zee(26) )
    
  end subroutine QML_NO3_Ezcoup_eea1
    
    
    !------------------------------------------------------------------------------
    !  3-mode coupling
    !      function vcoup_eea2(n, par, j, iref,avv)
  subroutine QML_NO3_Evcoup_eea2(vcoup_eea2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_eea2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=8
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pb(:),vee(:)

    pb  => avv%pb
    vee => avv%vee

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     4th order
          vcoup_eea2 = p(1) * pb(1) * vee(1)
    
    !     5th order
          vcoup_eea2 = vcoup_eea2 + pb(1)*(vee(2)*p(2)+vee(3)*p(3))
    
    !     6th order
          vcoup_eea2 = vcoup_eea2 + p(4) * pb(2) * vee(1)                           &
         & + pb(1) * (vee(4) * p(5) + vee(5) * p(6)                               &
         & + vee(6) * p(7) +  vee(7) * p(8))
    
  end subroutine QML_NO3_Evcoup_eea2
    
    !------------------------------------------------------------------------------
    !  3-mode coupling
    !      function wcoup_eea2(n,par, j, iref,avv)
  subroutine QML_NO3_Ewcoup_eea2(wcoup_eea2, n,par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: wcoup_eea2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=15
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pb(:),wee(:)

    pb  => avv%pb
    wee => avv%wee

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     4th order
          wcoup_eea2 = p(1) * pb(1) * wee(1)
    
    !     5th order
          wcoup_eea2 = wcoup_eea2 + pb(1) * (wee(2)*p(2)+wee(3)*p(3)        &
         & + wee(4) * p(4) + wee(5) * p(5))
    
    !     6th order
          wcoup_eea2 = wcoup_eea2 + p(6) * pb(2) * wee(1)                           &
         & + pb(1)*(wee(6)*p( 7) + wee(7) * p( 8) + wee(8) * p( 9)              &
         & + wee(9) * p(10) + wee(10) * p(11) + wee(11) * p(12)                   &
         & + wee(12) * p(13) + wee(13) * p(14) + wee(14) * p(15)  )
    
  end subroutine QML_NO3_Ewcoup_eea2
    
    !------------------------------------------------------------------------------
    !  3-mode coupling
    !      function zcoup_eea2(n,par, j, iref,avv)
  subroutine QML_NO3_Ezcoup_eea2(zcoup_eea2, n,par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: zcoup_eea2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=15
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pb(:),zee(:)

    pb  => avv%pb
    zee => avv%zee

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     4th order
          zcoup_eea2 = p(1) * pb(1) * zee(1)
    
    !     5th order
          zcoup_eea2 = zcoup_eea2 + pb(1) * (zee(2)*p(2)+zee(3)*p(3)              &
         & + zee(4) * p(4) + zee(5) * p(5))
    
    !     6th order
          zcoup_eea2 = zcoup_eea2 + p(6) * pb(2) * zee(1)                           &
         & + pb(1)*(zee(6)*p( 7) + zee(7) * p( 8) + zee(8) * p( 9)              &
         & + zee(9) * p(10) + zee(10) * p(11) + zee(11) * p(12)                   &
         & + zee(12) * p(13) + zee(13) * p(14) + zee(14) * p(15)  )
    
  end subroutine QML_NO3_Ezcoup_eea2
    
  !------------------------------------------------------------------------------
  !     function to generate V matrix elements for the coupling between e and a modes up to fourth order
  subroutine QML_NO3_Evcoup_e1aa2(vcoup_e1aa2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_e1aa2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=3
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),pb(:),va(:)

    pa => avv%pa
    pb => avv%pb
    va => avv%va

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     5th order
    vcoup_e1aa2 =  (p(1) * pa(1) * pb(1) * va(1))
    
    !     6th order
    vcoup_e1aa2 = vcoup_e1aa2 + (p(2)*pa(1)*pb(1)*va(2) + p(3)*pa(2)*pb(1)*va(1))
    
  end subroutine QML_NO3_Evcoup_e1aa2
    
    !------------------------------------------------------------------------------
    !     function to generate V matrix elements for the coupling between e and a modes up to fourth order
    
    !      function vcoup_e2aa2(n, par, j, iref,avv)
  subroutine QML_NO3_Evcoup_e2aa2(vcoup_e2aa2, n, par, j, iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: vcoup_e2aa2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=3
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),pb(:),vb(:)

    pa => avv%pa
    pb => avv%pb
    vb => avv%vb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     5th order
    vcoup_e2aa2 =  (p(1) * pa(1) * pb(1) * vb(1))
    
    !     6th order
    vcoup_e2aa2 = vcoup_e2aa2 + (p(2)*pa(1)*pb(1)*vb(2) + p(3)*pa(2)*pb(1)*vb(1))
    
  end subroutine QML_NO3_Evcoup_e2aa2
    
    !-------------------------------------------------------------------------------------
    !     function to generate W matrix elements for the coupling between e and a modes up to sixth order
    !     6 parameters
    
    !      function wcoup_e1aa2(n, par, j,iref,avv)
  subroutine QML_NO3_Ewcoup_e1aa2(wcoup_e1aa2, n, par, j,iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: wcoup_e1aa2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=7
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),pb(:),wa(:)

    pa => avv%pa
    pb => avv%pb
    wa => avv%wa

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     4th order
    wcoup_e1aa2 = p(1) * pa(1) * pb(1) * wa(1)
    
    !     5th order
    wcoup_e1aa2 = wcoup_e1aa2 + (p(2)*pa(1)*pb(1)*wa(2)                     &
         &  + p(3)*pa(2)*pb(1)*wa(1))
    
    !     6th order
    wcoup_e1aa2 = wcoup_e1aa2 + p(4)*pa(1)*pb(1)*wa(3)                      &
         &  + p(5)*pa(2)*pb(1)*wa(2)                                              &
         &  + p(6)*pa(3)*pb(1)*wa(1) + p(7)*pa(1)*pb(2)*wa(1)
    
    
  end subroutine QML_NO3_Ewcoup_e1aa2
    
    !-------------------------------------------------------------------------------------
    !     function to generate W matrix elements for the coupling between e and a modes up to sixth order
    !     6 parameters
    
    !      function wcoup_e2aa2(n, par, j,iref,avv)
  subroutine QML_NO3_Ewcoup_e2aa2(wcoup_e2aa2, n, par, j,iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: wcoup_e2aa2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=7
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),pb(:),wb(:)

    pa => avv%pa
    pb => avv%pb
    wb => avv%wb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     4th order
    wcoup_e2aa2 = p(1) * pa(1) * pb(1) * wb(1)
    
    !     5th order
    wcoup_e2aa2 = wcoup_e2aa2 + (p(2)*pa(1)*pb(1)*wb(2)                     &
         &  + p(3)*pa(2)*pb(1)*wb(1))
    
    !     6th order
    wcoup_e2aa2 = wcoup_e2aa2 + p(4)*pa(1)*pb(1)*wb(3)                      &
         &  + p(5)*pa(2)*pb(1)*wb(2)                                              &
         &  + p(6)*pa(3)*pb(1)*wb(1) + p(7)*pa(1)*pb(2)*wb(1)
    
    
  end subroutine QML_NO3_Ewcoup_e2aa2
    
    !-------------------------------------------------------------------------------------
    !     function to generate Z matrix elements for the coupling betzeen e and a modes up to sixth order
    !     6 parameters
    
    !      function zcoup_e1aa2(n, par, j,iref,avv)
  subroutine QML_NO3_Ezcoup_e1aa2(zcoup_e1aa2, n, par, j,iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: zcoup_e1aa2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=7
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),pb(:),za(:)

    pa => avv%pa
    pb => avv%pb
    za => avv%za

    p(:)   = ZERO
    p(2:j) = par(2:j)
 
    !     4th order
    zcoup_e1aa2 = p(1) * pa(1) * pb(1) * za(1)
    
    !     5th order
    zcoup_e1aa2 = zcoup_e1aa2 + (p(2)*pa(1)*pb(1)*za(2)                     &
         &  + p(3)*pa(2)*pb(1)*za(1))
    
    !     6th order
    zcoup_e1aa2 = zcoup_e1aa2 + p(4)*pa(1)*pb(1)*za(3)                      &
         &  + p(5)*pa(2)*pb(1)*za(2)                                              &
         &  + p(6)*pa(3)*pb(1)*za(1) + p(7)*pa(1)*pb(2)*za(1)
    
    
  end subroutine QML_NO3_Ezcoup_e1aa2
    
    !-------------------------------------------------------------------------------------
    !     function to generate W matrix elements for the coupling betzeen e and a modes up to sixth order
    !     6 parameters
    
    !      function zcoup_e2aa2(n, par, j,iref,avv)
  subroutine QML_NO3_Ezcoup_e2aa2(zcoup_e2aa2, n, par, j,iref,avv)
    implicit none

    real (kind=Rkind),            intent(inout)      :: zcoup_e2aa2
    real (kind=Rkind),            intent(in)         :: par(:)
    integer,                      intent(in)         :: n,j,iref
    TYPE(acoord_vjt_vee_t),       intent(in), target :: avv

    integer, parameter  :: lnx=7
    real (kind=Rkind)   :: p(lnx)
    real (kind=Rkind), pointer :: pa(:),pb(:),zb(:)

    pa => avv%pa
    pb => avv%pb
    zb => avv%zb

    p(:)   = ZERO
    p(2:j) = par(2:j)

    !     4th order
    zcoup_e2aa2 = p(1) * pa(1) * pb(1) * zb(1)
    
    !     5th order
    zcoup_e2aa2 = zcoup_e2aa2 + (p(2)*pa(1)*pb(1)*zb(2)                     &
         &  + p(3)*pa(2)*pb(1)*zb(1))
    
    !     6th order
    zcoup_e2aa2 = zcoup_e2aa2 + p(4)*pa(1)*pb(1)*zb(3)                      &
         &  + p(5)*pa(2)*pb(1)*zb(2)                                              &
         &  + p(6)*pa(3)*pb(1)*zb(1) + p(7)*pa(1)*pb(2)*zb(1)
    
  end subroutine QML_NO3_Ezcoup_e2aa2
  subroutine QML_NO3_vwzprec(q,avv)
    implicit none

    real(kind=Rkind),       intent(in)             :: q(6)
    TYPE(acoord_vjt_vee_t), intent(inout), target  :: avv
    
    integer :: j
    real(kind=Rkind) :: a, b, x3, x4, y3, y4

    real(kind=Rkind), pointer :: pa(:), pb(:)
    real(kind=Rkind), pointer :: va(:), vb(:)
    real(kind=Rkind), pointer :: wa(:), wb(:)
    real(kind=Rkind), pointer :: za(:), zb(:)
    real(kind=Rkind), pointer :: vee(:),  wee(:), zee(:)

    pa => avv%pa
    pb => avv%pb
    va => avv%va
    vb => avv%vb
    wa => avv%wa
    wb => avv%wb
    za => avv%za
    zb => avv%zb

    vee => avv%vee
    wee => avv%wee
    zee => avv%zee

    a=q(1)
    b=q(2)
    x3=q(3)
    y3=q(4)
    x4=q(5)
    y4=q(6)
    
    ! compute powers of a1 mode:
    pa(1) = a
    pa(2) = a**2
    pa(3) = a**3
    pa(4) = a**4
    pa(5) = a**5
    pa(6) = a**6
    pa(7) = a**7
    pa(8) = a**8

    !     compute powers of a2 mode:
    pb(1) = b**2
    pb(2) = b**4
    pb(3) = b**6
    pb(4) = b**8

!..   compute vjt terms for first e set:
    va(1) = x3**2 + y3**2
    va(2) = TWO * x3**3 - SIX * x3 * y3**2
    va(3) = x3**4 + TWO * x3**2 * y3**2 + y3**4
    va(4) = TWO * x3**5 - FOUR * x3**3 * y3**2 - SIX * x3 * y3**4
    va(5) = TWO * x3**6 - 30._Rkind*x3**4*y3**2 + 30._Rkind*x3**2*y3**4 -TWO * y3**6
    va(6) = x3**6 + THREE * x3**4 * y3**2 + THREE * x3**2 * y3**4 + y3**6

!..   compute vjt terms for second e set:
    vb(1) = x4**2 + y4**2
    vb(2) = TWO * x4**3 - SIX * x4 * y4**2
    vb(3) = x4**4 + TWO * x4**2 * y4**2 + y4**4
    vb(4) = TWO * x4**5 - FOUR * x4**3 * y4**2 - SIX * x4 * y4**4
    vb(5) = TWO * x4**6 - 30._Rkind*x4**4*y4**2 + 30._Rkind*x4**2*y4**4 -TWO * y4**6
    vb(6) = x4**6 + THREE * x4**4 * y4**2 + THREE * x4**2 * y4**4 + y4**6

!..   compute wjt terms for first e set:
    wa(1) = x3
    wa(2) = x3**2 - y3**2
    wa(3) = x3**3 + x3*y3**2
    wa(4) = x3**4 - SIX * x3**2 * y3**2 + y3**4
    wa(5) = x3**4 - y3**4
    wa(6) = x3**5 - TEN * x3**3 * y3**2 + FIVE * x3 * y3**4
    wa(7) = x3**5 + TWO * x3**3 * y3**2 + x3 * y3**4
    wa(8) = x3**6 - FIVE * x3**4 * y3**2 - FIVE * x3**2*y3**4+y3**6
    wa(9) = x3**6 + x3**4 * y3**2 - x3**2 * y3**4 - y3**6

!..   compute wjt terms for second e set:
    wb(1) = x4
    wb(2) = x4**2 - y4**2
    wb(3) = x4**3 + x4*y4**2
    wb(4) = x4**4 - SIX * x4**2 * y4**2 + y4**4
    wb(5) = x4**4 - y4**4
    wb(6) = x4**5 - TEN * x4**3 * y4**2 + FIVE * x4 * y4**4
    wb(7) = x4**5 + TWO * x4**3 * y4**2 + x4 * y4**4
    wb(8) = x4**6 - FIVE * x4**4 * y4**2 - FIVE * x4**2*y4**4+y4**6
    wb(9) = x4**6 + x4**4 * y4**2 - x4**2 * y4**4 - y4**6

!..   compute zjt terms for first e set:
    za(1) = y3
    za(2) = -TWO*x3*y3
    za(3) = y3**3 + y3*x3**2
    za(4) = FOUR*x3**3*y3 - FOUR*x3*y3**3
    za(5) = -TWO*x3**3*y3 - TWO*x3*y3**3
    za(6) = -FIVE*x3**4*y3 + TEN*x3**2*y3**3 - y3**5
    za(7) =  x3**4*y3 + TWO*x3**2*y3**3 + y3**5
    za(8) = FOUR*x3**5*y3 - FOUR*x3*y3**5
    za(9) = -TWO*x3**5*y3 - FOUR*x3**3*y3**3 - TWO*x3*y3**5

!..   compute zjt terms for second e set:
    zb(1) = y4
    zb(2) = -TWO*x4*y4
    zb(3) = y4**3 + y4*x4**2
    zb(4) = FOUR*x4**3*y4 - FOUR*x4*y4**3
    zb(5) = -TWO*x4**3*y4 - TWO*x4*y4**3
    zb(6) = -FIVE*x4**4*y4 + TEN*x4**2*y4**3 - y4**5
    zb(7) =  x4**4*y4 + TWO*x4**2*y4**3 + y4**5
    zb(8) = FOUR*x4**5*y4 - FOUR*x4*y4**5
    zb(9) = -TWO*x4**5*y4 - FOUR*x4**3*y4**3 - TWO*x4*y4**5

!..   compute Vcoup_ee terms:
    vee(1) =  TWO*x3*x4 + TWO*y3*y4
    vee(2) = TWO*x3*x4**2 - TWO*x3*y4**2 - FOUR*x4*y3*y4
    vee(3) = TWO*x3**2*x4 - TWO*x4*y3**2 - FOUR*x3*y3*y4
    vee(4) = TWO*x3**2*x4**2 - TWO*x3**2*y4**2                                &
   &          - TWO*x4**2*y3**2 + TWO*y3**2*y4**2                               &
   &          + EIGHT*x3*x4*y3*y4
    vee(5) = TWO*x3**3*x4 + TWO*x3**2*y3*y4                                   &
   &          + TWO*x3*x4*y3**2 + TWO*y3**3*y4
    vee(6) = TWO*x4**3*x3 + TWO*x3*x4*y4**2                                   &
   &         + TWO*x4**2*y3*y4 + TWO*y4**3*y3
    vee(7) = x3**2*x4**2 + x3**2*y4**2                                          &
   &          + x4**2*y3**2 + y3**2*y4**2

    vee(8) = x3*y4**4-4*x4*y3*y4**3-6*x3*x4**2*y4**2+4*x4**3*y3*y4              &
   &          +x3*x4**4
    vee(9) = (y4**2+x4**2)*(x3*y4**2+2*x4*y3*y4-x3*x4**2)
    vee(10) = (y3**2+x3**2)*(x3*y4**2+2*x4*y3*y4-x3*x4**2)
    vee(11) = x4*(y3**2+x3**2)*(THREE*y4**2-x4**2)
    vee(12) = (TWO*x3*y3*y4+x4*y3**2-x3**2*x4)*(y4**2+x4**2)
    vee(13) = x3*(THREE*y3**2-x3**2)*(y4**2+x4**2)
    vee(14) = FOUR*x3*y3**3*y4-FOUR*x3**3*y3*y4-x4*y3**4                        &
   &          +SIX*x3**2*x4*y3**2-x3**4*x4
    vee(15) = (y3**2+x3**2)*(TWO*x3*y3*y4+x4*y3**2-x3**2*x4)


    vee(16) = (y3*y4**5-FIVE*x3*x4*y4**4-TEN*x4**2*y3*y4**3                   &
   &          +TEN*x3*x4**3*y4**2+FIVE*x4**4*y3*y4-x3*x4**5)
    vee(17) = (y3*y4**2-x3*y4**2-TWO*x4*y3*y4-TWO*x3*x4*y4                    &
   &          -x4**2*y3+x3*x4**2)*(y3*y4**2+x3*y4**2+TWO*x4*y3*y4                &
   &          -TWO*x3*x4*y4-x4**2*y3-x3*x4**2)
    vee(18) = (y3*y4-x3*x4)*(y3**2*y4**2-THREE*x3**2*y4**2                       &
   &          -EIGHT*x3*x4*y3*y4-THREE*x4**2*y3**2+x3**2*x4**2)
    vee(19) = (y3**2*y4-TWO*x3*y3*y4-x3**2*y4-x4*y3**2                         &
   &          -TWO*x3*x4*y3+x3**2*x4)*(y3**2*y4+TWO*x3*y3*y4                    &
   &          -x3**2*y4+x4*y3**2-TWO*x3*x4*y3-x3**2*x4)
    vee(20) = (y3**5*y4-TEN*x3**2*y3**3*y4+FIVE*x3**4*y3*y4                   &
   &          -FIVE*x3*x4*y3**4+TEN*x3**3*x4*y3**2-x3**5*x4)

    vee(21) = (y3*y4+x3*x4)*(y4**2+x4**2)**2
    vee(22) = (y3*y4-x3*y4+x4*y3+x3*x4)*(y3*y4+x3*y4                            &
   &          -x4*y3+x3*x4)*(y4**2+x4**2)
    vee(23) = (y3**2+x3**2)*(y4**2+x4**2)**2
    vee(24) = (y3*y4+x3*x4)*(y3**2*y4**2-THREE*x3**2*y4**2                       &
   &          +EIGHT*x3*x4*y3*y4-THREE*x4**2*y3**2+x3**2*x4**2)
    vee(25) = (y3**2+x3**2)*(y3*y4+x3*x4)*(y4**2+x4**2)
    vee(26) = (y3**2+x3**2)*(y3*y4-x3*y4+x4*y3+x3*x4)                           &
   &          *(y3*y4+x3*y4-x4*y3+x3*x4)
    vee(27) = (y3**2+x3**2)**2*(y4**2+x4**2)
    vee(28) = (y3**2+x3**2)**2*(y3*y4+x3*x4)

!..   compute Wcoup_ee terms:
    wee(1) = x3 * x4 - y3 * y4
    wee(2) = x3**2 * x4 + x4 * y3**2
    wee(3) = x3 * y4**2 + x3 * x4**2
    wee(4) = TWO * x3 * y3 * y4 + x3**2 * x4 - x4 * y3**2
    wee(5) = x3 * x4**2 + TWO * x4 * y3*y4 - x3 * y4**2
    wee(6) = x3**3 * x4 - THREE * x3**2 *y3*y4                                   &
   &            - THREE*x3*x4*y3**2 + y3**3*y4
    wee(7) = x3**2*x4**2 - x3**2*y4**2 - x4**2*y3**2 + y3**2*y4**2              &
   & - FOUR*x3*x4*y3*y4
    wee(8) =x4**3*x3 - THREE*x3*x4*y4**2 - THREE*x4**2*y3*y4+y4**3*y3
!               wee44(9) =x3**3*x4 - THREE*x3**2*y3*y4 -THREE*x3*x4*y3**2-y3**3*y4
    wee(9) =x3**3*x4 + THREE*x3**2*y3*y4 - THREE*x3*x4*y3**2-y3**3*y4
!               wee45(10) =x4**3*x3 - THREE*x3*x4*y4**2 - THREE*x4**2*y3*y4-y4**3*y3
    wee(10) =x4**3*x3 - THREE*x3*x4*y4**2 + THREE*x4**2*y3*y4-y4**3*y3
    wee(11) = x3**3*x4 - x3**2*y3*y4 + x3*x4*y3**2 - y3**3*y4
    wee(12) = x3**2*x4**2 + x3**2*y4**2 - x4**2*y3**2 - y3**2*y4**2
    wee(13) = x3**2*x4**2 - x3**2*y4**2 + x4**2*y3**2 - y3**2*y4**2
    wee(14) = x3*x4**3 + x3*x4*y4**2 - x4**2*y3*y4 - y4**3*y3

    wee(15) = (x3*y4**4+FOUR*x4*y3*y4**3-SIX*x3*x4**2*y4**2                    &
   &          -FOUR*x4**3*y3*y4+x3*x4**4)
    wee(16) = (TWO*x3*y3*y4**3+THREE*x4*y3**2*y4**2                             &
   & -THREE*x3**2*x4*y4**2-SIX*x3*x4**2*y3*y4-x4**3*y3**2+x3**2*x4**3)
    wee(17) = (THREE*x3*y3**2*y4**2-x3**3*y4**2+TWO*x4*y3**3*y4                 &
   &          -SIX*x3**2*x4*y3*y4-THREE*x3*x4**2*y3**2+x3**3*x4**2)
    wee(18) = (FOUR*x3*y3**3*y4-FOUR*x3**3*y3*y4+x4*y3**4                       &
   &          -SIX*x3**2*x4*y3**2+x3**4*x4)

    wee(19) = -(y4**2+x4**2)*(x3*y4**2-TWO*x4*y3*y4-x3*x4**2)
    wee(20) = x3*(y4**2+x4**2)**2
    wee(21) = -(TWO*x3*y3*y4**3-THREE*x4*y3**2*y4**2+THREE*x3**2                 &
   &          *x4*y4**2-SIX*x3*x4**2*y3*y4+x4**3*y3**2-x3**2*x4**3)
    wee(22) = x4*(y3**2+x3**2)*(y4**2+x4**2)
    wee(23) = -(y3**2+x3**2)*(x3*y4**2-TWO*x4*y3*y4-x3*x4**2)
    wee(24) = x3*(y3**2+x3**2)*(y4**2+x4**2)
    wee(25) = x4*(y3**2+x3**2)**2
    wee(26) = (y3**2+x3**2)*(TWO*x3*y3*y4-x4*y3**2+x3**2*x4)

    wee(27) = (y3*y4**5+FIVE*x3*x4*y4**4-TEN*x4**2*y3*y4**3               &
   &    -TEN*x3*x4**3*y4**2+FIVE*x4**4*y3*y4+x3*x4**5)
    wee(28) = (y4**2+x4**2)*(y3*y4**3-THREE*x3*x4*y4**2                      &
   &    -THREE*x4**2*y3*y4+x3*x4**3)
    wee(29) = (y3**2+x3**2)*(y4**2-TWO*x4*y4-x4**2)                        &
   &    *(y4**2+TWO*x4*y4-x4**2)
    wee(30) = (y3*y4-x3*y4-x4*y3-x3*x4)*(y3*y4+x3*y4                        &
   &    +x4*y3-x3*x4)*(y4**2+x4**2)
    wee(31) = (y3**2+x3**2)*(y3*y4**3-THREE*x3*x4*y4**2                      &
   &    -THREE*x4**2*y3*y4+x3*x4**3)
    wee(32) = (y3**3*y4-THREE*x3**2*y3*y4-THREE*x3*x4*y3**2+x3**3*x4)         &
   &    *(y4**2+x4**2)
    wee(33) = (y3**2+x3**2)*(y3*y4-x3*y4-x4*y3-x3*x4)                       &
   &    *(y3*y4+x3*y4+x4*y3-x3*x4)
    wee(34) = (y3**2-TWO*x3*y3-x3**2)*(y3**2+TWO*x3*y3-x3**2)             &
   &    *(y4**2+x4**2)
    wee(35) = (y3**2+x3**2)*(y3**3*y4-THREE*x3**2*y3*y4                      &
   &    -THREE*x3*x4*y3**2+x3**3*x4)
    wee(36)= (y3**5*y4-TEN*x3**2*y3**3*y4+FIVE*x3**4*y3*y4                &
   &    +FIVE*x3*x4*y3**4-TEN*x3**3*x4*y3**2+x3**5*x4)

    wee(37) = (y3*y4-x3*x4)*(y4**2+x4**2)**2
    wee(38) = (y4**2+x4**2)*(y3*y4**3+THREE*x3*x4*y4**2                      &
   &    -THREE*x4**2*y3*y4-x3*x4**3)
    wee(39) = (y3**2+x3**2)*(y4-x4)*(y4+x4)*(y4**2+x4**2)
    wee(40) = (y3-x3)*(y3+x3)*(y4**2+x4**2)**2
    wee(41) = (y3*y4**2-x3*y4**2+TWO*x4*y3*y4+TWO*x3*x4*y4                &
   &    -x4**2*y3+x3*x4**2)*(y3*y4**2+x3*y4**2-TWO*x4*y3*y4                  &
   &    +TWO*x3*x4*y4-x4**2*y3-x3*x4**2)
    wee(42) = (y3**3*y4-THREE*x3**2*y3*y4+THREE*x3*x4*y3**2-x3**3*x4)         &
   &    *(y4**2+x4**2)
    wee(43) = (y3**2+x3**2)*(y3*y4-x3*x4)*(y4**2+x4**2)
    wee(44) = (y3**2+x3**2)*(y3*y4**3+THREE*x3*x4*y4**2                      &
   &    -THREE*x4**2*y3*y4-x3*x4**3)
    wee(45) = (y3**2*y4-TWO*x3*y3*y4-x3**2*y4+x4*y3**2                     &
   &    +TWO*x3*x4*y3-x3**2*x4)*(y3**2*y4+TWO*x3*y3*y4-x3**2*y4             &
   &    -x4*y3**2+TWO*x3*x4*y3+x3**2*x4)
    wee(46) = (y3-x3)*(y3+x3)*(y3**2+x3**2)*(y4**2+x4**2)
    wee(47) = (y3**2+x3**2)**2*(y4-x4)*(y4+x4)
    wee(48) = (y3**2+x3**2)*(y3**3*y4-THREE*x3**2*y3*y4                      &
   &    +THREE*x3*x4*y3**2-x3**3*x4)
    wee(49) = (y3**2+x3**2)**2*(y3*y4-x3*x4)

!..   compute Zcoup_ee terms:
    zee(1) = -x3*y4 - x4*y3

    zee(2) = y3**2*y4 + x3**2*y4
    zee(3) = y3*y4**2 + x4**2*y3
    zee(4) = y3**2*y4 + TWO*x3*x4*y3 - x3**2*y4
    zee(5) = y4**2*y3 + TWO*x3*x4*y4 - x4**2*y3

    zee(6) = x3**3*y4 + THREE*x3**2*x4*y3                                        &
   &          - THREE*x3*y3**2*y4 - x4*y3**3
    zee(7) = TWO*x3**2*x4*y4 + TWO*x3*x4**2*y3 - TWO*x3*y3*y4**2             &
   & - TWO*x4*y3**2*y4
    zee(8) =x4**3*y3 +THREE*x3*x4**2*y4 - THREE*x4*y3*y4**2 - x3*y4**3
    zee(9) =x3**3*y4 -THREE*x3**2*x4*y3 - THREE*x3*y3**2*y4 + x4*y3**3
    zee(10)=x4**3*y3 -THREE*x3*x4**2*y4 - THREE*x4*y3*y4**2 + x3*y4**3
    zee(11) = - x3**3*y4 - x3**2*x4*y3 - x3*y3**2*y4 - x4*y3**3
    zee(12) = - TWO*x3*x4**2*y3 - TWO*x3*y3*y4**2
    zee(13) = - TWO*x3**2*x4*y4 - TWO*x4*y3**2*y4
    zee(14) = - y4**3*x3 - x3*x4**2*y4 - x4*y3*y4**2 - y3*x4**3

    zee(15) = -(y3*y4**4-FOUR*x3*x4*y4**3-SIX*x4**2*y3*y4**2                                   &
   &           +FOUR*x3*x4**3*y4+x4**4*y3)
    zee(16) = -(y3**2*y4**3-x3**2*y4**3-SIX*x3*x4*y3*y4**2                                     &
   &     -THREE*x4**2*y3**2*y4+THREE*x3**2*x4**2*y4+TWO*x3*x4**3*y3)
    zee(17) = -(y3**3*y4**2-THREE*x3**2*y3*y4**2-SIX*x3*x4*y3**2*y4                             &
   &           +TWO*x3**3*x4*y4-x4**2*y3**3+THREE*x3**2*x4**2*y3)
    zee(18)=-(y3**4*y4-SIX*x3**2*y3**2*y4+x3**4*y4-FOUR*x3*x4*y3**3                          &
   &           +FOUR*x3**3*x4*y3)

    zee(19) = (y4**2+x4**2)*(y3*y4**2+TWO*x3*x4*y4-x4**2*y3)
    zee(20) = y3*(y4**2+x4**2)**2
    zee(21) = (y3**2*y4**3-x3**2*y4**3+SIX*x3*x4*y3*y4**2                                      &
   &     -THREE*x4**2*y3**2*y4+THREE*x3**2*x4**2*y4-TWO*x3*x4**3*y3)
    zee(22) = (y3**2+x3**2)*y4*(y4**2+x4**2)
    zee(23) = (y3**2+x3**2)*(y3*y4**2+TWO*x3*x4*y4-x4**2*y3)
    zee(24) = y3*(y3**2+x3**2)*(y4**2+x4**2)
    zee(25) = (y3**2+x3**2)**2*y4
    zee(26) = (y3**2+x3**2)*(y3**2*y4-x3**2*y4+TWO*x3*x4*y3)

    zee(27) = (x3*y4**5-FIVE*x4*y3*y4**4-TEN*x3*x4**2*y4**3               &
   &   +TEN*x4**3*y3*y4**2+FIVE*x3*x4**4*y4-x4**5*y3)
    zee(28) = -(y4**2+x4**2)*(x3*y4**3+THREE*x4*y3*y4**2                     &
   &   -THREE*x3*x4**2*y4-x4**3*y3)
    zee(29) = -FOUR*x4*(y3**2+x3**2)*y4*(y4-x4)*(y4+x4)
    zee(30) = -TWO*(x3*y4+x4*y3)*(y3*y4-x3*x4)*(y4**2+x4**2)
    zee(31) = -(y3**2+x3**2)*(x3*y4**3+THREE*x4*y3*y4**2                     &
   &   -THREE*x3*x4**2*y4-x4**3*y3)
    zee(32) = -(THREE*x3*y3**2*y4-x3**3*y4+x4*y3**3                          &
   &   -THREE*x3**2*x4*y3)*(y4**2+x4**2)
    zee(33) = -TWO*(y3**2+x3**2)*(x3*y4+x4*y3)*(y3*y4-x3*x4)
    zee(34) = -FOUR*x3*y3*(y3-x3)*(y3+x3)*(y4**2+x4**2)
    zee(35) = -(y3**2+x3**2)*(THREE*x3*y3**2*y4-x3**3*y4+x4*y3**3            &
   &   -THREE*x3**2*x4*y3)
    zee(36)= -(FIVE*x3*y3**4*y4-TEN*x3**3*y3**2*y4+x3**5*y4              &
   &   -x4*y3**5+TEN*x3**2*x4*y3**3-FIVE*x3**4*x4*y3)


    zee(37) = (x3*y4+x4*y3)*(y4**2+x4**2)**2
    zee(38) = -(y4**2+x4**2)*(x3*y4**3-THREE*x4*y3*y4**2                    &
   &   -THREE*x3*x4**2*y4+x4**3*y3)
    zee(40) = TWO*x3*y3*(y4**2+x4**2)**2                 !permutation 39/40 DML 5/10/2021
    zee(39) = TWO*x4*(y3**2+x3**2)*y4*(y4**2+x4**2)      !permutation 39/40 DML 5/10/2021
    zee(41) = -TWO*(x3*y4**2-TWO*x4*y3*y4-x3*x4**2)                      &
   &   *(y3*y4**2+TWO*x3*x4*y4-x4**2*y3)
    zee(42) = (THREE*x3*y3**2*y4-x3**3*y4-x4*y3**3+THREE*x3**2*x4*y3)        &
   &   *(y4**2+x4**2)
    zee(43) = (y3**2+x3**2)*(x3*y4+x4*y3)*(y4**2+x4**2)
    zee(44) = -(y3**2+x3**2)*(x3*y4**3-THREE*x4*y3*y4**2                    &
   &   -THREE*x3*x4**2*y4+x4**3*y3)
    zee(45) = TWO*(TWO*x3*y3*y4-x4*y3**2+x3**2*x4)                       &
   &   *(y3**2*y4-x3**2*y4+TWO*x3*x4*y3)
    zee(46) = TWO*x3*y3*(y3**2+x3**2)*(y4**2+x4**2)
    zee(47) = TWO*x4*(y3**2+x3**2)**2*y4
    zee(48) = (y3**2+x3**2)*(THREE*x3*y3**2*y4-x3**3*y4-x4*y3**3            &
   &   +THREE*x3**2*x4*y3)
    zee(49) = (y3**2+x3**2)**2*(x3*y4+x4*y3)
    
  end subroutine QML_NO3_vwzprec
    
    !---------------------------------------------------
    ! A Viel 30.09.2021
    ! from cartesian coordinates to internal of NO3E" JCP20217 potential
    ! N (0,1-3) O (2-4,1-3)   1-3 : x, y, z
  subroutine QML_NO3_trans_coord(x,qinter,transformcoordblock,tmc)
    USE QDUtil_m
    implicit none

    TYPE (transformcoordblock_t), intent(in)         :: transformcoordblock
    TYPE (tmc_t),                 intent(in)         :: tmc
    real (kind=Rkind),            intent(inout)      :: qinter(:)
    real (kind=Rkind),            intent(in)         :: x(0:3,3)

    integer i, j, k

    real (kind=Rkind) ::   xvec(3,0:2),xvecinvers(3,3)

    real (kind=Rkind) ::   xintern(0:6) !beta, 3 distances, 3 angles
    real (kind=Rkind) ::   tnorm,det,t(3)
    real (kind=Rkind) :: f  ! function for morse evaluation
    real (kind=Rkind) :: rnodist,rnoprod
    
    real (kind=Rkind) ::   phi_ref,le_ref,beta_ref

    real (kind=Rkind), parameter :: sq2 = ONE/sqrt(TWO)
    real (kind=Rkind), parameter :: sq3 = ONE/sqrt(THREE)
    real (kind=Rkind), parameter :: sq6 = ONE/sqrt(SIX)

    logical, parameter :: debug = .FALSE.

    if (debug) write(out_unit,*) 'in QML_NO3_trans_coord'  ; flush(6)

    phi_ref  = transformcoordblock%phi_ref
    le_ref   = transformcoordblock%le_ref
    beta_ref = transformcoordblock%beta_ref
    !
    ! write cartesian
          !do i=0,3
          !  write(out_unit,*)'cart',i,(x(i,j),j=1,3)
          !enddo
    
    ! copy of internal.f subroutine
    ! NH vectors matrix (3x3)
    do i=0,2
      xvec(:,i) = x(i+1,:)-x(0,:)
      if (debug) write(out_unit,*) 'xvec',i,xvec(:,i)
    enddo
    ! no distances
    do i=0,2
      xintern(i+1)=sqrt(dot_product(xvec(:,i),xvec(:,i)))
    enddo
    if (debug) write(out_unit,*) 'xintern',xintern

    rnodist = sum(xintern(1:3))
    rnoprod = product(xintern(1:3))
    if (debug) write(out_unit,*) 'rnodist,rnoprod',rnodist,rnoprod


    ! normalization
    do i=0,2
      xvec(:,i) = xvec(:,i)/xintern(i+1)
      if (debug) write(out_unit,*) 'xvec',i,xvec(:,i)
    enddo

    ! determinant
    !det=QML_NO3_det3t3(xvec)
    det = Det_OF(xvec)

    if (abs(det) < ONETENTH**13) then
      ! planar case
      call QML_NO3_planarnew(xvec,xintern(4),xintern(0))
    else
      ! Pyramdal case
      ! trisector
      !call QML_NO3_mat_invers(xvec,xvecinvers,det)
      xvecinvers = inv_OF_Mat_TO(xvec)
      t(1)=sum(xvecinvers(:,1))
      t(2)=sum(xvecinvers(:,2))
      t(3)=sum(xvecinvers(:,3))
      tnorm=sqrt(t(1)**2+t(2)**2+t(3)**2)
      xintern(0)=acos(ONE/tnorm)
      ! projected angles
      call QML_NO3_proangnew(xvec,t,tnorm,xintern(0),xintern(4))
    endif
    !
    ! xintern(0)=trisector
    ! xintern(1)=ditance no1
    ! xintern(2)=ditance no2
    ! xintern(3)=ditance no3
    ! xintern(4)=angle o2no3
    ! xintern(5)=angle o3no1
    ! xintern(6)=angle o1no2
    ! 13.06.2012 order coherent with surface definition
    !         write(out_unit,*)xintern(0),'trisec'   !printPES
    !         write(out_unit,*)xintern(1),'no1   '
    !         write(out_unit,*)xintern(2),'no2   '
    !         write(out_unit,*)xintern(3),'no3   '
    !         write(out_unit,*)xintern(4),'o2no3 '
    !         write(out_unit,*)xintern(5),'o3no1 '
    !         write(out_unit,*)xintern(6),'o1no2 '
    !
    !...    now transform angles so they properly dissociate:
    !        a1=(a1)/((no2+reno)*(no3+reno))
    !        a2=(a2)/((no1+reno)*(no3+reno))
    !        a3=(a3)/((no1+reno)*(no2+reno))
    xintern(4)=xintern(4)/(xintern(2)*xintern(3))
    xintern(5)=xintern(5)/(xintern(1)*xintern(3))
    xintern(6)=xintern(6)/(xintern(1)*xintern(2))
    !
    ! Displace NO distances
    xintern(1:3) = xintern(1:3) -le_ref
    xintern(0)   = xintern(0)   -beta_ref ! planar = 0 for Eisfeld PES
    !       write(out_unit,*)xintern(0),"umbrella"  ! printPES
    !
    ! compute morse or tmc morse
    do i=1,3
      call QML_NO3_ff(xintern(i),1,f,tmc)
      xintern(i)= ONE - exp(-f*xintern(i))
    enddo

    ! symmetrization
    qinter(1)=(xintern(1)+xintern(2)+xintern(3))*sq3
    qinter(2)=xintern(0)/rnoprod*le_ref**3
    qinter(3)=(TWO*xintern(1)-xintern(2)-xintern(3))*sq6
    qinter(4)=(xintern(2)-xintern(3))*sq2
    qinter(5)=(TWO*xintern(4)-xintern(5)-xintern(6))*sq6
    qinter(6)=(xintern(5)-xintern(6))*sq2
    ! radial variable for damping in the trisector variable
    qinter(7)=rnodist/THREE

    if (debug) then
      write(out_unit,'(7g20.10)') qinter
      write(out_unit,*) 'end QML_NO3_trans_coord'
      flush(out_unit)
    END IF

  End subroutine QML_NO3_trans_coord
    !***********************************************************
    !
    !   calculation of a 3*3 determinant
    !
    !***********************************************************
    
          function QML_NO3_det3t3(M)
          implicit none
          integer n
          parameter (n=3)
          real (kind=Rkind) :: QML_NO3_det3t3, M(n,n)
    
           QML_NO3_det3t3 = M(1,1)*M(2,2)*M(3,3)       &
         &         +M(1,2)*M(2,3)*M(3,1)       &
         &         +M(1,3)*M(2,1)*M(3,2)       &
         &         -M(1,3)*M(2,2)*M(3,1)       &
         &         -M(1,1)*M(2,3)*M(3,2)       &
         &         -M(1,2)*M(2,1)*M(3,3)
    
          end function QML_NO3_det3t3
    
    !***********************************************************
    !
    !   makes the inverse of a 3*3 matrix
    !
    !***********************************************************
  subroutine QML_NO3_mat_invers(M,Im,det)
          implicit none
          integer i,j,ma
          parameter (ma=3)
          real (kind=Rkind) :: M(ma,ma), A(ma,ma), Im(ma,ma)
          real (kind=Rkind) :: det
    
    !.....makes the inverse of a 3*3 matrix
    !     ---------------------------------
    !      det=QML_NO3_det3t3(M)
    !      if (abs(det) < ONETENTH**10) then
    !         write(out_unit,*) 'determinant zero! matrix cannot be inverted!'
    !        stop
    !      endif
    
    !.....calculation of determinats of the submatrices
    
          A(1,1)= M(2,2)*M(3,3)-M(2,3)*M(3,2)
          A(1,2)=-M(2,1)*M(3,3)+M(2,3)*M(3,1)
          A(1,3)= M(2,1)*M(3,2)-M(2,2)*M(3,1)
    
          A(2,1)=-M(1,2)*M(3,3)+M(1,3)*M(3,2)
          A(2,2)= M(1,1)*M(3,3)-M(3,1)*M(1,3)
          A(2,3)=-M(1,1)*M(3,2)+M(1,2)*M(3,1)
    
          A(3,1)= M(1,2)*M(2,3)-M(1,3)*M(2,2)
          A(3,2)=-M(1,1)*M(2,3)+M(1,3)*M(2,1)
          A(3,3)= M(1,1)*M(2,2)-M(1,2)*M(2,1)
    
    !.....calculation of the inverse (Im) of the matrix
    
          do i=1,ma
           do j=1,ma
             Im(j,i)=(ONE/det)*A(i,j)
           enddo
          enddo
    
  end subroutine QML_NO3_mat_invers
    
    
    !********************************************************************
    !
    ! calculates projected angles in plane -> trisector is normal vector
    ! of this plane (for nonplanar systems)
    !
    ! rotate trisector on z axis (2 rotations) -> projected angles are
    ! angles of projections in xy plane (use xy components of position
    ! vectors of rotated system)
    !
    !********************************************************************
  subroutine QML_NO3_proangnew(M,t,tnorm,beta,phi)
          implicit none
          integer i,j, n, d
          parameter (n=3, d=2)
          real (kind=Rkind) :: M(n,*), Q(n,n), S(n,n)
          real (kind=Rkind) :: t(n), t_z(n), phi(n)
          real (kind=Rkind) :: v2(3), v3(3), v4(3), z(3)
          real (kind=Rkind) :: tnorm
          real (kind=Rkind) :: epsi, alpha, delta, beta, beta2
          real (kind=Rkind) :: sia, l, u, ws, betr, rad
          real (kind=Rkind) :: dum
    
          do i=1,3
          do j=1,3
               Q(j,i)=ZERO
               S(j,i)=ZERO
          enddo
             phi(i)=ZERO
             v2(i)=ZERO
             v3(i)=ZERO
             v4(i)=ZERO
          enddo
    
    !..   unit vector along z axis
          do i=1,3
              z(i)=ZERO
          enddo
          z(3)=ONE
    
    !..   calculate angle of trisector with z axis
          alpha=t(3)/tnorm
          sia=sqrt(ONE-alpha**2)
          alpha=acos(alpha)
          if (alpha < ZERO) sia=-sia
    
    
    !      alpha=angle(z,t,n)
    !      sia=sin(alpha)
    
          if (abs(sia) < ONETENTH**10) then
            do i=1,3
              do j=1,3
                S(i,j)=M(i,j)
              enddo
            enddo
    
          else
    
    !..   use polar coordinates; calculate angle (delta) in xy plane
    
            l=tnorm*sia
    !        delta=acos(t(1)/l)
            l=t(1)/l
            l=dmax1(-ONE,l)
            l=dmin1(ONE,l)
            delta=acos(l)
    
            if (delta > (pi*HALF)) then
              delta=pi-delta
            endif
    
            delta=delta*sign(ONE,t(1))*sign(ONE,t(2))
            epsi=-delta
    !
    !..     rotate trisector in xz plane (rot. around z axis)
    
            call QML_NO3_rotate(M,n,Q,n,epsi)
    
            if (alpha > (pi*HALF)) then
              alpha=pi-alpha
            endif
    
            call QML_NO3_trise(Q,beta2,t_z)
    
    
            alpha=alpha*sign(ONE,t_z(1))*sign(ONE,t_z(3))
    
    !..     rotate trisector on z axis (rot. around y axis)
    
            call QML_NO3_rotate(Q,n,S,d,alpha)
    
          endif
    
    !..   calculate projected angles phi(1-3)
    
          do i=1,3
            if (i == 3) then
              v2(i)=ZERO
              v3(i)=ZERO
              v4(i)=ZERO
            else
              v2(i)=S(i,1)
              v3(i)=S(i,2)
              v4(i)=S(i,3)
            endif
          enddo
    
    
    
    !      phi(1)=angle(v3,v4,n)
    !      phi(2)=angle(v2,v4,n)
    !      phi(3)=angle(v2,v3,n)
           dum=ONE/sqrt(v2(1)**2+v2(2)**2)
           v2(1)=v2(1)*dum
           v2(2)=v2(2)*dum
           dum=ONE/sqrt(v3(1)**2+v3(2)**2)
           v3(1)=v3(1)*dum
           v3(2)=v3(2)*dum
           dum=ONE/sqrt(v4(1)**2+v4(2)**2)
           v4(1)=v4(1)*dum
           v4(2)=v4(2)*dum
           phi(1)=v3(1)*v4(1)+v3(2)*v4(2)
           phi(2)=v2(1)*v4(1)+v2(2)*v4(2)
           phi(3)=v2(1)*v3(1)+v2(2)*v3(2)
          do i=1,3
          if (abs(phi(i)-ONE) < ONETENTH**10) then
            phi(i)=ZERO
          elseif (abs(phi(i)+ONE) < ONETENTH**10) then
            phi(i)=pi
          else
            phi(i)=acos(phi(i))
          endif
          enddo
    
    
    !..   scalp gives wrong angle -> in this angle lies third vector!!
    !..   to prevent this:
    
          ws=(phi(1)+phi(2)+phi(3))-(2*pi)
    
          if (abs(ws) > ONETENTH**10) then
            if (abs(phi(1)+phi(2)-phi(3)) < ONETENTH**10) then
               phi(3)=(TWO*pi)-phi(3)
            else if (abs(phi(2)+phi(3)-phi(1)) < ONETENTH**10) then
               phi(1)=(TWO*pi)-phi(1)
            else if (abs(phi(1)+phi(3)-phi(2)) < ONETENTH**10) then
               phi(2)=(TWO*pi)-phi(2)
            endif
          endif
    
    !       write(out_unit,*)
    !       write(out_unit,*) 'matrix S:'
    !      do i=1,3
    !          write(out_unit,*) (S(i,j), j=1,3)
    !      enddo
    !       write(out_unit,*)
    
  end subroutine QML_NO3_proangnew
    !----------------------------------------------------------------
    !***********************************************************
    !
    !     rotation of a matrix
    !     n defines rotational axis (x: n=1; y: n=2; z: n=3)
    !     dimensions: 3*u matrix rotated
    !
    !***********************************************************
    
  subroutine QML_NO3_rotate(m1,u,m2,n,theta)
          implicit none
          integer i,j,n,max,p,u
          parameter (max=3)
          real (kind=Rkind) :: m1(max,*), m2(max,*), rotmat(max,max)
          real (kind=Rkind) :: theta, s, c
    
    !      call init_m(rotmat,max,max)
           do i=1,3
           do j=1,3
           rotmat(j,i)=ZERO
           enddo
           enddo
    
          if ((n < 1).or.(n > 3)) then
             write(out_unit,*) 'rotational axis not specified! n is not 1,2 or 3'
            stop
          endif
    
          s = sin(theta)
          c = cos(theta)
    
          do i=1,3
            if (i.ne.n) then
              do j=1,i-1
                if (j.ne.n) then
                  rotmat(i,j)=s
                  rotmat(j,i)=-s
                endif
              enddo
            endif
            rotmat(i,i)=c
          enddo
    
          rotmat(n,n)=ONE
    
          call QML_NO3_matmult(rotmat,max,max,m1,u,m2)
    
    !      do i=1,max
    !           write(out_unit,100) (rotmat(i,j),j=1,3)
    !      enddo
    !       write(out_unit,*)
    !
    !       write(out_unit,*) 'matrix rotated:'
    !      do i=1,max
    !           write(out_unit,100) (m2(i,j),j=1,3)
    !      enddo
    !       write(out_unit,*)
    100   format(3e15.6)
    
  end subroutine QML_NO3_rotate
    
    !----------------------------------------------------------------
    !***********************************************************
    !
    ! matrix multiplication of two matrices
    ! sizes: M1 -> n1*n2; M2 -> n2*n3; EM -> n1*n3
    !
    !***********************************************************
  subroutine QML_NO3_matmult(M1,n1,n2,M2,n3,EM)
          implicit none
          integer i,j,k,n1,n2,n3
          real (kind=Rkind) :: M1(n1,*), M2(n2,*), EM(n1,*)
    
    
          do i=1,n1
             do j=1,n3
               EM(i,j)=ZERO
                do k=1,n2
                  EM(i,j)=EM(i,j)+(M1(i,k)*M2(k,j))
                enddo
             enddo
          enddo
    
  end subroutine QML_NO3_matmult

  !****************************************************************
  !
  ! calculation of trisector and pyramidalization angle beta
  ! for non planar systems
  !
  !****************************************************************
  subroutine QML_NO3_trise(M,beta,t)
          implicit none
          Integer i,k,j,n
          parameter (n=3)
          real (kind=Rkind) :: M(n,n), Im(n,n)
          real (kind=Rkind) :: t(3)
          real (kind=Rkind) :: beta, betrag
          real (kind=Rkind) :: norm,det
    
          do i=1,3
             t(i)=ZERO
          enddo
    
    !.....calculate trisector angle beta
    
          det=QML_NO3_det3t3(M)
          call QML_NO3_mat_invers(M,Im,det)
    
    !.....calculate trisector t
    !..   transpose of Im because x,y,z values on columns!
    !..   (M_transpose)^-1 = (M^-1)_transpose
    
          do i=1,3
            t(1)=t(1)+Im(i,1)
            t(2)=t(2)+Im(i,2)
            t(3)=t(3)+Im(i,3)
          enddo
    !      betrag=norm(t,3)
           betrag=sqrt(t(1)**2+t(2)**2+t(3)**2)
    
    !..   calculate angle beta
    
          beta=acos(ONE/betrag)
    
  end subroutine QML_NO3_trise
    
    !----------------------------------------------------------------
    !********************************************************
    !
    ! treat case when molecule is planar (determinant of
    ! NH bonds is zero)
    !
    !********************************************************
  subroutine QML_NO3_planarnew(M,phi,beta)
    implicit none
          integer i,j,k,n
          parameter (n=3)
          real (kind=Rkind) :: M(n,*), phi(n)
          real (kind=Rkind) :: v2(n), v3(n), v4(n)
          real (kind=Rkind) :: scalp, angle, ws, beta
          real (kind=Rkind) :: s1, s2, s3,dum
    
    !.....columns of matrix M as vectors
    
          do i=1,3
            v2(i)=M(i,1)
            v3(i)=M(i,2)
            v4(i)=M(i,3)
          enddo
    
    !      s1=scalp(v2,v3,n)
    !      s2=scalp(v2,v4,n)
    !      s3=scalp(v3,v4,n)
           s1=v2(1)*v3(1)+v2(2)*v3(2)+v2(3)*v3(3)
           s2=v2(1)*v4(1)+v2(2)*v4(2)+v2(3)*v4(3)
           s3=v3(1)*v4(1)+v3(2)*v4(2)+v3(3)*v4(3)
    
    !.....check if atoms collinear (on one line)
    
          if ((abs(s1-ONE) < ONETENTH**07).and.(abs(s2-ONE)       &
         &       < ONETENTH**07).and.(abs(s3-ONE) < ONETENTH**07)) then
            beta=ZERO
    
             write(out_unit,*)
             write(out_unit,*) 'all atoms collinear in same direction!!'
             write(out_unit,*) 'coordinates not defined!!'
             write(out_unit,*)
    
            stop
          else
    
    !.....planar case and atoms not collinear => pyramidalization angle
    !..   is 90°, projected angles = bonding angles
    
            beta=pi*HALF
    
    !        phi(1)=angle(v3,v4,3)
    !        phi(2)=angle(v2,v4,3)
    !        phi(3)=angle(v2,v3,3)
           dum=ONE/sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
           v2(1)=v2(1)*dum
           v2(2)=v2(2)*dum
           v2(3)=v2(3)*dum
           dum=ONE/sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
           v3(1)=v3(1)*dum
           v3(2)=v3(2)*dum
           v3(3)=v3(3)*dum
           dum=ONE/sqrt(v4(1)**2+v4(2)**2+v4(3)**2)
           v4(1)=v4(1)*dum
           v4(2)=v4(2)*dum
           v4(3)=v4(3)*dum
           phi(1)=v3(1)*v4(1)+v3(2)*v4(2)+v3(3)*v4(3)
           phi(2)=v2(1)*v4(1)+v2(2)*v4(2)+v2(3)*v4(3)
           phi(3)=v2(1)*v3(1)+v2(2)*v3(2)+v2(3)*v3(3)
          do i=1,3
          if (abs(phi(i)-ONE) < ONETENTH**10) then
            phi(i)=ZERO
          else
            phi(i)=acos(phi(i))
          endif
          enddo
          endif
    
    !..   scalp gives wrong angle -> in this angle lies third vector!!
    !..   sum of angles is not 360°
    !..   to prevent this:
    
          ws=(phi(1)+phi(2)+phi(3))-(TWO*pi)
    
          if (abs(ws) > ONETENTH**6) then
            if (abs(phi(1)+phi(2)-phi(3)) < ONETENTH**10) then
               phi(3)=(TWO*pi)-phi(3)
            else if (abs(phi(2)+phi(3)-phi(1)) < ONETENTH**10) then
               phi(1)=(TWO*pi)-phi(1)
            else if (abs(phi(1)+phi(3)-phi(2)) < ONETENTH**10) then
               phi(2)=(TWO*pi)-phi(2)
            endif
          endif
    
  end subroutine QML_NO3_planarnew
    
END MODULE QML_NO3_m
