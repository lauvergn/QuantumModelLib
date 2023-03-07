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

!> @brief Derived type in which the NO3 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_NO3_t
    PRIVATE

    integer           :: iref,n,npst

    real (kind=Rkind) :: p(337),e0ref,e0rho
    integer           :: pst(2,100)

    integer           :: npoly(2)
    real (kind=Rkind) :: a(21)  !copy of the first elements of p array

    real (kind=Rkind) :: phi_ref,le_ref,beta_ref
    real (kind=Rkind) :: pi,sq2,sq3,sq6

    !common /no3minus_pot/e0ref,e0rho,p,pst,iref

    !common /tmc/a,npoly

    !common /transformcoordblock/ &
    !   phi_ref,le_ref,beta_ref   &
    !   ,pi,sq2,sq3,sq6

   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_NO3
    PROCEDURE :: Write_QModel    => Write_QML_NO3
    PROCEDURE :: Write0_QModel   => Write0_QML_NO3
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

    TYPE (QML_NO3_t)                          :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
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

    QModel%nsurf    = 2
    QModel%ndim     = 6
    QModel%pot_name = 'no3'

    IF (QModel%option /= 1) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1) ! unpublished model

      ! all parameters
      QModel%iref = 2
      QModel%n    = 337
      QModel%npst = 30

      !pi = FOUR*atan(ONE)
      QModel%sq2 = ONE/sqrt(TWO)
      QModel%sq3 = ONE/sqrt(THREE)
      QModel%sq6 = ONE/sqrt(SIX)

      !reference geometry
      QModel%phi_ref  = TWO/THREE*pi  ! 120°
      QModel%le_ref   = 2.344419_Rkind
      !for completness reno=2.344419d0,reoo=4.0606528d0,rad=0.017453293d0)
      QModel%beta_ref = pi*HALF
      Write(out_unit,*)'# ref geometry:'
      Write(out_unit,*)'# phi_ref: ',QModel%phi_ref
      Write(out_unit,*)'# le_ref:  ',QModel%le_ref
      Write(out_unit,*)'# beta_ref:',QModel%beta_ref


      ! tmc-ang-prec-37um2b_extra_b5_20150316.par
      ! including more extra points, 16.03.2015
      !                               Fitting RMS:      0.85624381E-03   (   187.9 cm-1)
      !                               Merit function:   0.85624381E-03
      !                               Reduced RMS:      0.30159465E-03   (    66.2 cm-1) for data below  -279.500000
      !                               Sum of point weights:  0.9458
      !                               Number of weighted points:    2404
      QModel%p(  1)=  -0.123440260E+00_Rkind
      QModel%p(  2)=   0.231012750E+01_Rkind
      QModel%p(  3)=  -0.131815690E+01_Rkind
      QModel%p(  4)=   0.505372680E+01_Rkind
      QModel%p(  5)=  -0.202745100E+02_Rkind
      QModel%p(  6)=   0.123400003E+00_Rkind
      QModel%p(  7)=   0.556995241E+00_Rkind
      QModel%p(  8)=   0.505947070E+01_Rkind
      QModel%p(  9)=  -0.840552473E+02_Rkind
      QModel%p( 10)=   0.303348770E+01_Rkind
      QModel%p( 11)=   0.915482720E+00_Rkind
      QModel%p( 12)=   0.597447070E+01_Rkind
      QModel%p( 13)=   0.129964090E+01_Rkind
      QModel%p( 14)=   0.923982010E+01_Rkind
      QModel%p( 15)=   0.123400003E+00_Rkind
      QModel%p( 16)=   0.537126740E+01_Rkind
      QModel%p( 17)=  -0.528341500E+01_Rkind
      QModel%p( 18)=   0.662076740E+02_Rkind
      QModel%p( 19)=   0.833679770E+01_Rkind
      QModel%p( 20)=   0.286309010E+02_Rkind
      QModel%p( 21)=   0.123400003E+00_Rkind
      QModel%p( 22)=  -0.304431280E+01_Rkind
      QModel%p( 23)=   0.748191810E+01_Rkind
      QModel%p( 24)=  -0.448249980E+01_Rkind
      QModel%p( 25)=   0.206854010E+02_Rkind
      QModel%p( 26)=  -0.172197980E+02_Rkind
      QModel%p( 27)=  -0.642527480E+02_Rkind
      QModel%p( 28)=   0.565738800E+02_Rkind
      QModel%p( 29)=   0.303284550E+02_Rkind
      QModel%p( 30)=   0.302292120E+02_Rkind
      QModel%p( 31)=   0.792999760E+02_Rkind
      QModel%p( 32)=  -0.711077690E+02_Rkind
      QModel%p( 33)=   0.575018080E+02_Rkind
      QModel%p( 34)=  -0.347775870E+02_Rkind
      QModel%p( 35)=   0.930777310E+01_Rkind
      QModel%p( 36)=  -0.145823480E+02_Rkind
      QModel%p( 37)=   0.520607240E+03_Rkind
      QModel%p( 38)=  -0.123662040E+04_Rkind
      QModel%p( 39)=   0.125896900E+04_Rkind
      QModel%p( 40)=  -0.681959190E+03_Rkind
      QModel%p( 41)=   0.181199110E+03_Rkind
      QModel%p( 42)=  -0.419560080E+03_Rkind
      QModel%p( 43)=   0.105702100E+04_Rkind
      QModel%p( 44)=  -0.607165880E+03_Rkind
      QModel%p( 45)=  -0.640175490E+03_Rkind
      QModel%p( 46)=   0.701414450E+03_Rkind
      QModel%p( 47)=   0.237988350E+03_Rkind
      QModel%p( 48)=  -0.595329960E+03_Rkind
      QModel%p( 49)=   0.202039810E+03_Rkind
      QModel%p( 50)=  -0.872841948E+01_Rkind
      QModel%p( 51)=  -0.409836740E+00_Rkind
      QModel%p( 52)=   0.801472910E+01_Rkind
      QModel%p( 53)=   0.526269750E+02_Rkind
      QModel%p( 54)=   0.227930150E+02_Rkind
      QModel%p( 55)=  -0.349545650E+02_Rkind
      QModel%p( 56)=   0.380465090E+02_Rkind
      QModel%p( 57)=  -0.744954770E+02_Rkind
      QModel%p( 58)=   0.959469230E+01_Rkind
      QModel%p( 59)=   0.123400003E+00_Rkind
      QModel%p( 60)=  -0.531160440E+01_Rkind
      QModel%p( 61)=  -0.872510030E+01_Rkind
      QModel%p( 62)=   0.254766890E+02_Rkind
      QModel%p( 63)=   0.545690110E+03_Rkind
      QModel%p( 64)=  -0.234308810E+03_Rkind
      QModel%p( 65)=  -0.676611530E+01_Rkind
      QModel%p( 66)=  -0.131862020E+04_Rkind
      QModel%p( 67)=   0.385114210E+00_Rkind
      QModel%p( 68)=  -0.162871150E+03_Rkind
      QModel%p( 69)=   0.123400003E+00_Rkind
      QModel%p( 70)=  -0.487582430E-01_Rkind
      QModel%p( 71)=  -0.493363570E-01_Rkind
      QModel%p( 72)=   0.493413080E+00_Rkind
      QModel%p( 73)=  -0.786792930E+00_Rkind
      QModel%p( 74)=   0.148580440E+01_Rkind
      QModel%p( 75)=   0.275701800E+01_Rkind
      QModel%p( 76)=  -0.235189530E+01_Rkind
      QModel%p( 77)=   0.147278820E+02_Rkind
      QModel%p( 78)=   0.123400003E+00_Rkind
      QModel%p( 79)=  -0.238514350E+00_Rkind
      QModel%p( 80)=   0.297004000E+00_Rkind
      QModel%p( 81)=  -0.120273480E+00_Rkind
      QModel%p( 82)=   0.191271270E+01_Rkind
      QModel%p( 83)=  -0.226871020E+01_Rkind
      QModel%p( 84)=  -0.114691680E+02_Rkind
      QModel%p( 85)=  -0.897240320E+02_Rkind
      QModel%p( 86)=   0.147408570E+02_Rkind
      QModel%p( 87)=   0.123400003E+00_Rkind
      QModel%p( 88)=  -0.624627310E+00_Rkind
      QModel%p( 89)=   0.984178850E-02_Rkind
      QModel%p( 90)=  -0.212520820E+01_Rkind
      QModel%p( 91)=  -0.127773090E+01_Rkind
      QModel%p( 92)=   0.408838780E+01_Rkind
      QModel%p( 93)=   0.310835790E+01_Rkind
      QModel%p( 94)=  -0.787746700E+01_Rkind
      QModel%p( 95)=   0.572379100E+01_Rkind
      QModel%p( 96)=  -0.883120710E-01_Rkind
      QModel%p( 97)=  -0.845179740E+01_Rkind
      QModel%p( 98)=  -0.923696170E+01_Rkind
      QModel%p( 99)=   0.150396360E+02_Rkind
      QModel%p(100)=   0.129763030E+02_Rkind
      QModel%p(101)=  -0.127820180E+02_Rkind
      QModel%p(102)=   0.724479130E+01_Rkind
      QModel%p(103)=  -0.226340140E+02_Rkind
      QModel%p(104)=   0.193995670E+02_Rkind
      QModel%p(105)=  -0.101918030E+02_Rkind
      QModel%p(106)=   0.120124630E+03_Rkind
      QModel%p(107)=  -0.195995920E+02_Rkind
      QModel%p(108)=  -0.821429710E+01_Rkind
      QModel%p(109)=  -0.146887710E+02_Rkind
      QModel%p(110)=  -0.731600850E+02_Rkind
      QModel%p(111)=   0.673921420E+02_Rkind
      QModel%p(112)=   0.150943880E+02_Rkind
      QModel%p(113)=  -0.867078970E+01_Rkind
      QModel%p(114)=   0.147406530E+03_Rkind
      QModel%p(115)=  -0.260872930E+02_Rkind
      QModel%p(116)=  -0.840659750E+03_Rkind
      QModel%p(117)=   0.504515730E+03_Rkind
      QModel%p(118)=   0.621689150E+03_Rkind
      QModel%p(119)=  -0.524583500E+03_Rkind
      QModel%p(120)=  -0.382318100E+02_Rkind
      QModel%p(121)=   0.201076520E+03_Rkind
      QModel%p(122)=  -0.733776780E+02_Rkind
      QModel%p(123)=  -0.296733560E+02_Rkind
      QModel%p(124)=  -0.211305860E+03_Rkind
      QModel%p(125)=  -0.382490340E+03_Rkind
             !p(126)=   0.626682910E+03_Rkind ! to fix a symmetry pbs.
             !p(127)=   0.648222640E+03_Rkind ! we take the average. DML 05/10/2021
      QModel%p(126)=   (0.626682910E+03_Rkind+0.648222640E+03_Rkind)/TWO
      QModel%p(127)=   (0.626682910E+03_Rkind+0.648222640E+03_Rkind)/TWO
      QModel%p(128)=   0.510189060E+02_Rkind
      QModel%p(129)=  -0.300347270E+03_Rkind
      QModel%p(130)=  -0.131928350E+04_Rkind
      QModel%p(131)=  -0.902324530E+01_Rkind
      QModel%p(132)=  -0.121951920E+01_Rkind
      QModel%p(133)=   0.617623300E+03_Rkind
      QModel%p(134)=   0.307972350E+03_Rkind
      QModel%p(135)=  -0.228997370E+02_Rkind
      QModel%p(136)=  -0.163914390E+03_Rkind
      QModel%p(137)=  -0.169780300E+00_Rkind
      QModel%p(138)=   0.125500700E+01_Rkind
      QModel%p(139)=   0.482992570E+00_Rkind
      QModel%p(140)=   0.317988550E+01_Rkind
      QModel%p(141)=  -0.436770420E+01_Rkind
      QModel%p(142)=   0.356941810E+00_Rkind
      QModel%p(143)=   0.900703840E+01_Rkind
      QModel%p(144)=  -0.282621440E+02_Rkind
      QModel%p(145)=  -0.137130480E+02_Rkind
      QModel%p(146)=   0.604692600E+01_Rkind
      QModel%p(147)=  -0.276681540E+02_Rkind
      QModel%p(148)=  -0.564288980E+02_Rkind
      QModel%p(149)=   0.279994500E+02_Rkind
      QModel%p(150)=  -0.750562690E+02_Rkind
      QModel%p(151)=  -0.802644670E+02_Rkind
      QModel%p(152)=   0.511596150E+02_Rkind
      QModel%p(153)=  -0.274657630E+02_Rkind
      QModel%p(154)=   0.885353980E+02_Rkind
      QModel%p(155)=   0.203854730E+00_Rkind
      QModel%p(156)=   0.343485190E+01_Rkind
      QModel%p(157)=   0.115761950E+01_Rkind
      QModel%p(158)=  -0.108960670E+02_Rkind
      QModel%p(159)=   0.810512110E+01_Rkind
      QModel%p(160)=  -0.149700950E+01_Rkind
      QModel%p(161)=  -0.632516450E+02_Rkind
      QModel%p(162)=  -0.182788980E+03_Rkind
      QModel%p(163)=  -0.607637350E+02_Rkind
      QModel%p(164)=   0.796587110E+02_Rkind
      QModel%p(165)=  -0.325214470E+02_Rkind
      QModel%p(166)=   0.246171810E+03_Rkind
      QModel%p(167)=  -0.139468350E+02_Rkind
      QModel%p(168)=   0.659294940E+03_Rkind
      QModel%p(169)=  -0.142032340E+01_Rkind
      QModel%p(170)=   0.418184160E+03_Rkind
      QModel%p(171)=   0.151015790E+00_Rkind
      QModel%p(172)=  -0.796880880E+01_Rkind
      QModel%p(173)=   0.217108089E+00_Rkind
      QModel%p(174)=  -0.232760231E+02_Rkind
      QModel%p(175)=   0.269121775E+03_Rkind
      QModel%p(176)=  -0.365601790E+03_Rkind
      QModel%p(177)=   0.245177337E+02_Rkind
      QModel%p(178)=   0.737181669E+02_Rkind
      QModel%p(179)=   0.230068180E+04_Rkind
      QModel%p(180)=  -0.919786868E+03_Rkind
      QModel%p(181)=   0.203758567E+01_Rkind
      QModel%p(182)=   0.125347943E+02_Rkind
      QModel%p(183)=  -0.240018368E+02_Rkind
      QModel%p(184)=   0.315636912E+02_Rkind
      QModel%p(185)=  -0.586391234E+03_Rkind
      QModel%p(186)=  -0.541519483E+02_Rkind
      QModel%p(187)=   0.541543516E+02_Rkind
      QModel%p(188)=  -0.110047897E+01_Rkind
      QModel%p(189)=  -0.984270831E+01_Rkind
      QModel%p(190)=   0.176431382E+02_Rkind
      QModel%p(191)=   0.619129772E+03_Rkind
      QModel%p(192)=  -0.309888310E+02_Rkind
      QModel%p(193)=   0.516577868E+03_Rkind
      QModel%p(194)=   0.763453420E+03_Rkind
      QModel%p(195)=   0.476767245E+00_Rkind
      QModel%p(196)=   0.396602201E+00_Rkind
      QModel%p(197)=  -0.157742464E+02_Rkind
      QModel%p(198)=  -0.345412065E+02_Rkind
      QModel%p(199)=   0.106966811E+02_Rkind
      QModel%p(200)=   0.685088030E+03_Rkind
      QModel%p(201)=   0.673178507E+03_Rkind
      QModel%p(202)=  -0.106366340E+04_Rkind
      QModel%p(203)=  -0.332883674E+04_Rkind
      QModel%p(204)=   0.481034455E+04_Rkind
      QModel%p(205)=  -0.189700034E+02_Rkind
      QModel%p(206)=   0.152614812E+03_Rkind
      QModel%p(207)=  -0.709471913E+03_Rkind
      QModel%p(208)=   0.126748116E+03_Rkind
      QModel%p(209)=  -0.240266010E+03_Rkind
      QModel%p(210)=   0.148598985E+04_Rkind
      QModel%p(211)=   0.115713133E+04_Rkind
      QModel%p(212)=  -0.897094470E+03_Rkind
      QModel%p(213)=  -0.834863905E+03_Rkind
      QModel%p(214)=   0.545372351E+02_Rkind
      QModel%p(215)=   0.317678018E+03_Rkind
      QModel%p(216)=   0.899708039E+02_Rkind
      QModel%p(217)=   0.461995866E+03_Rkind
      QModel%p(218)=  -0.657237461E+03_Rkind
      QModel%p(219)=  -0.129355471E+04_Rkind
      QModel%p(220)=  -0.139661058E+03_Rkind
      QModel%p(221)=  -0.189084040E+03_Rkind
      QModel%p(222)=  -0.209609912E+02_Rkind
      QModel%p(223)=   0.312119590E+02_Rkind
      QModel%p(224)=  -0.875470382E+02_Rkind
      QModel%p(225)=  -0.578320555E+02_Rkind
      QModel%p(226)=   0.819108852E+01_Rkind
      QModel%p(227)=   0.120109292E+03_Rkind
      QModel%p(228)=   0.207725902E+03_Rkind
      QModel%p(229)=  -0.220092259E+03_Rkind
      QModel%p(230)=  -0.219373924E+03_Rkind
      QModel%p(231)=  -0.104377843E+04_Rkind
      QModel%p(232)=   0.188456908E+03_Rkind
      QModel%p(233)=  -0.187745233E+02_Rkind
      QModel%p(234)=   0.115501296E+03_Rkind
      QModel%p(235)=  -0.680341930E+02_Rkind
      QModel%p(236)=   0.135039895E+03_Rkind
      QModel%p(237)=  -0.286175591E+03_Rkind
      QModel%p(238)=  -0.150033749E+02_Rkind
      QModel%p(239)=   0.642889917E+03_Rkind
      QModel%p(240)=   0.496437690E+01_Rkind
      QModel%p(241)=   0.579415330E+01_Rkind
      QModel%p(242)=   0.819178020E+00_Rkind
      QModel%p(243)=  -0.641777110E+01_Rkind
      QModel%p(244)=   0.162797280E+03_Rkind
      QModel%p(245)=  -0.180520690E+03_Rkind
      QModel%p(246)=  -0.536904030E+03_Rkind
      QModel%p(247)=   0.604811650E+03_Rkind
      QModel%p(248)=   0.284787840E+03_Rkind
      QModel%p(249)=  -0.135967270E+03_Rkind
      QModel%p(250)=   0.153868400E+02_Rkind
      QModel%p(251)=  -0.240188550E+03_Rkind
      QModel%p(252)=   0.314345890E+03_Rkind
      QModel%p(253)=   0.664042750E+03_Rkind
      QModel%p(254)=  -0.120296430E+04_Rkind
      QModel%p(255)=   0.296121720E+03_Rkind
      QModel%p(256)=  -0.123934530E+03_Rkind
      QModel%p(257)=  -0.652689590E+02_Rkind
      QModel%p(258)=   0.147837820E+04_Rkind
      QModel%p(259)=  -0.532670070E+04_Rkind
      QModel%p(260)=  -0.202987460E+04_Rkind
      QModel%p(261)=   0.272007950E+04_Rkind
      QModel%p(262)=   0.454035570E+04_Rkind
      QModel%p(263)=  -0.174507240E+04_Rkind
      QModel%p(264)=   0.259475430E+03_Rkind
      QModel%p(265)=   0.640285120E+03_Rkind
      QModel%p(266)=  -0.489340180E+01_Rkind
      QModel%p(267)=  -0.180109170E+01_Rkind
      QModel%p(268)=   0.214753300E+02_Rkind
      QModel%p(269)=  -0.103364920E+02_Rkind
      QModel%p(270)=   0.336812050E+01_Rkind
      QModel%p(271)=   0.331311770E+01_Rkind
      QModel%p(272)=  -0.184333060E+02_Rkind
      QModel%p(273)=   0.363236360E+02_Rkind
      QModel%p(274)=  -0.733841910E+02_Rkind
      QModel%p(275)=   0.733114270E+02_Rkind
      QModel%p(276)=  -0.103140980E+03_Rkind
      QModel%p(277)=   0.643779480E+02_Rkind
      QModel%p(278)=  -0.300323700E+03_Rkind
      QModel%p(279)=   0.116331260E+02_Rkind
      QModel%p(280)=   0.333795200E+03_Rkind
      QModel%p(281)=   0.151003270E+02_Rkind
      QModel%p(282)=  -0.294676480E+02_Rkind
      QModel%p(283)=   0.137832210E+02_Rkind
      QModel%p(284)=  -0.150401300E+02_Rkind
      QModel%p(285)=   0.128243380E+03_Rkind
      QModel%p(286)=   0.515344760E+03_Rkind
      QModel%p(287)=  -0.104975170E+04_Rkind
      QModel%p(288)=   0.696065160E+03_Rkind
      QModel%p(289)=  -0.850704770E+02_Rkind
      QModel%p(290)=   0.552233050E+03_Rkind
      QModel%p(291)=  -0.126727650E+03_Rkind
      QModel%p(292)=   0.553415220E+03_Rkind
      QModel%p(293)=  -0.532309150E+02_Rkind
      QModel%p(294)=  -0.112746260E+04_Rkind
      QModel%p(295)=   0.186326010E+03_Rkind
      QModel%p(296)=  -0.364935110E+03_Rkind
      QModel%p(297)=   0.207343810E+03_Rkind
      QModel%p(298)=  -0.293492870E+03_Rkind
      QModel%p(299)=  -0.922632000E+01_Rkind
      QModel%p(300)=  -0.199198270E+03_Rkind
      QModel%p(301)=   0.302964050E+03_Rkind
      QModel%p(302)=  -0.278387770E+03_Rkind
      QModel%p(303)=   0.140546710E+03_Rkind
      QModel%p(304)=   0.529782740E+03_Rkind
      QModel%p(305)=  -0.273028890E+03_Rkind
      QModel%p(306)=  -0.149698510E+03_Rkind
      QModel%p(307)=  -0.142095790E+04_Rkind
      QModel%p(308)=   0.748850490E+03_Rkind
      QModel%p(309)=   0.925607750E+03_Rkind
      QModel%p(310)=  -0.460069010E+03_Rkind
      QModel%p(311)=  -0.206325460E+03_Rkind
      QModel%p(312)=   0.857738680E-02_Rkind
      ! Vertical Energy
      QModel%p(313)=  -0.279544420E+03_Rkind
      QModel%p(313)=   ZERO               !-0.279544420E+03_Rkind
      QModel%p(314)=   0.000000000E+00_Rkind
      QModel%p(315)=   0.109055120E+00_Rkind
      QModel%p(316)=   0.000000000E+00_Rkind
      QModel%p(317)=   0.000000000E+00_Rkind
      QModel%p(318)=   0.696349180E-01_Rkind
      QModel%p(319)=   0.000000000E+00_Rkind
      QModel%p(320)=   0.281557390E+00_Rkind
      QModel%p(321)=  -0.827685040E-01_Rkind
      QModel%p(322)=   0.156327730E-01_Rkind
      QModel%p(323)=  -0.416376780E-03_Rkind
      QModel%p(324)=   0.630990090E-04_Rkind
      QModel%p(325)=   0.000000000E+00_Rkind
      QModel%p(326)=   0.750000000E+00_Rkind
      QModel%p(327)=   0.000000000E+00_Rkind
      QModel%p(328)=   0.000000000E+00_Rkind
      QModel%p(329)=   0.000000000E+00_Rkind
      QModel%p(330)=   0.000000000E+00_Rkind
      QModel%p(331)=   0.000000000E+00_Rkind
      QModel%p(332)=   0.000000000E+00_Rkind
      QModel%p(333)=   0.000000000E+00_Rkind
      QModel%p(334)=   0.000000000E+00_Rkind
      QModel%p(335)=   0.000000000E+00_Rkind
      QModel%p(336)=   0.953536010E+01_Rkind
      QModel%p(337)=   0.103373110E+01_Rkind

             ! Vertical energy
             !    e0ref=  -0.279402540D+03
            QModel%e0ref=  -0.279402540E+03_Rkind +0.279544420E+03_Rkind
            QModel%e0rho=   0.800000000E-01_Rkind
            QModel%pst(1,  1)=  1
            QModel%pst(2,  1)=  6
            QModel%pst(1,  2)=  7
            QModel%pst(2,  2)=  3
            QModel%pst(1,  3)= 10
            QModel%pst(2,  3)=  6
            QModel%pst(1,  4)= 16
            QModel%pst(2,  4)=  6
            QModel%pst(1,  5)= 22
            QModel%pst(2,  5)= 28
            QModel%pst(1,  6)= 50
            QModel%pst(2,  6)= 10
            QModel%pst(1,  7)= 60
            QModel%pst(2,  7)= 10
            QModel%pst(1,  8)= 70
            QModel%pst(2,  8)=  9
            QModel%pst(1,  9)= 79
            QModel%pst(2,  9)=  9
            QModel%pst(1, 10)= 88
            QModel%pst(2, 10)= 49
            QModel%pst(1, 11)=137
            QModel%pst(2, 11)= 18
            QModel%pst(1, 12)=155
            QModel%pst(2, 12)= 18
            QModel%pst(1, 13)=173
            QModel%pst(2, 13)=  4
            QModel%pst(1, 14)=177
            QModel%pst(2, 14)=  4
            QModel%pst(1, 15)=181
            QModel%pst(2, 15)=  7
            QModel%pst(1, 16)=188
            QModel%pst(2, 16)=  7
            QModel%pst(1, 17)=195
            QModel%pst(2, 17)=  2
            QModel%pst(1, 18)=197
            QModel%pst(2, 18)=  8
            QModel%pst(1, 19)=205
            QModel%pst(2, 19)= 15
            QModel%pst(1, 20)=220
            QModel%pst(2, 20)=  3
            QModel%pst(1, 21)=223
            QModel%pst(2, 21)=  3
            QModel%pst(1, 22)=226
            QModel%pst(2, 22)=  7
            QModel%pst(1, 23)=233
            QModel%pst(2, 23)=  7
            QModel%pst(1, 24)=240
            QModel%pst(2, 24)= 26
            QModel%pst(1, 25)=266
            QModel%pst(2, 25)= 46
            QModel%pst(1, 26)=312
            QModel%pst(2, 26)=  1
            QModel%pst(1, 27)=313
            QModel%pst(2, 27)=  1
            QModel%pst(1, 28)=314
            QModel%pst(2, 28)= 11
            QModel%pst(1, 29)=325
            QModel%pst(2, 29)= 11
            QModel%pst(1, 30)=336
            QModel%pst(2, 30)=  2

                   Write(out_unit,*)'# NO3 E" 2 states 16/03/2015'
                   Write(out_unit,*)'# tmc-ang-prec-37um2b_extra_b5_20150316.par  '
                   Write(out_unit,*)'# switching method with total energy and polynomial'
                   Write(out_unit,*)'# iref ',QModel%iref
                   Write(out_unit,*)'# iref is not used'
                   Write(out_unit,*)'# npar ',QModel%n
                   do i=1,QModel%npst
                     Write(out_unit,2)i,QModel%pst(1,i),i,QModel%pst(2,i)
                   enddo
                   2    format('# pst(1,',i3,')=',i5,' pst(2,',i3,')=',i5)

                   do i=1,QModel%n
                     Write(out_unit,1)i,QModel%p(i)
                   enddo
                   1    format('# p(',i3,')=',f16.10)

                   Write(out_unit,*)'# e0ref au and ev', QModel%e0ref,QModel%e0ref*27.21d0
                   Write(out_unit,*)'# e0rho ',QModel%e0rho


             ! tmc parameter
             ! now in p array directly
             !      do i=1,11
             !        a(1,i)=0.d0
             !        a(2,i)=0.d0
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
                   do i=QModel%pst(1,28),QModel%pst(1,28)+20
                     QModel%a(j)=QModel%p(i)
                     !Write(out_unit,*)'ici',i,QModel%p(i)
                     j=j+1
                   enddo
                   QModel%a(2)=abs(QModel%a(2))  ! DML: it was in ff subroutine
                   QModel%a(5)=abs(QModel%a(5))  ! DML: it was in ff subroutine
                   QModel%a(4)=abs(QModel%a(4))  ! DML: it was in ff subroutine


                   QModel%npoly(1)=5
                   QModel%npoly(2)=5

                   Write(out_unit,*)'# TMC     npoly(1) : ',QModel%npoly(1)
                   Write(out_unit,*)'# TMC     npoly(2) : ',QModel%npoly(2)

             ! A Viel 2013.09.05
             ! rs.f put at the end
             !
             ! A Viel 2012.07.13
             ! Version with iref=1 hard coded
             !
             ! ATTENTION limited to 32 threads (mx=32)

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
  SUBROUTINE Write0_QML_NO3(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_NO3_t),   intent(in) :: QModel
    integer,               intent(in) :: nio

    write(nio,*) 'NO3 default parameters'
    write(nio,*)
    write(nio,*) 'end NO3 default parameters'


  END SUBROUTINE Write0_QML_NO3

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

    CLASS(QML_NO3_t),  intent(in)    :: QModel
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
  SUBROUTINE EvalPot1_QML_NO3(Mat_OF_PotDia,dnQ,NO3Pot,nderiv)
    !Unpublished model potential (yet)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)    :: dnQ(:) ! BLA, Tors, HOOP
    TYPE(QML_NO3_t),     intent(in)    :: NO3Pot
    integer,             intent(in)    :: nderiv


    !local variables (derived type). They have to be deallocated
    TYPE (dnS_t)        :: dnPot,BLA,Tors,HOOP
    TYPE (dnS_t)        :: MorseBLAP,Morsemin,Overlap
    TYPE (dnS_t)        :: Hdir2D,Hct2D
    real (kind=Rkind)   :: d1,d4

    real (kind=Rkind) :: XNO3(0:3,3)
    real (kind=Rkind) :: vdia(2,2)

    real (kind=Rkind),   parameter     :: EnergyConv = 627.509_Rkind
    real (kind=Rkind),   parameter     :: LenghtConv = 0.52917721067121_Rkind


    ! test geometry
    ! 0 is N
    XNO3(0,1)=2.9506658729288038e-003_Rkind
    XNO3(0,2)=ZERO
    XNO3(0,3)=-3.7754294953142864e-003_Rkind

    ! 1,2,3 are the 3 Oxygen atoms
    XNO3(1,1)=9.5060742943492427e-004_Rkind
    XNO3(1,2)=ZERO
    XNO3(1,3)=2.3817433236872598_Rkind

    XNO3(2,1)=2.0639330934921434_Rkind
    XNO3(2,2)=ZERO
    XNO3(2,3)=-1.1939442518313312_Rkind

    XNO3(3,1)=-2.0674669214991690_Rkind
    XNO3(3,2)=ZERO
    XNO3(3,3)=-1.1844937951551615_Rkind

stop
    !CALL potential_NO3(XNO3,vdia,NO3Pot)

    Mat_OF_PotDia(1,1) = vdia(1,1)
    Mat_OF_PotDia(1,2) = vdia(1,2)
    Mat_OF_PotDia(2,1) = vdia(2,1)
    Mat_OF_PotDia(2,2) = vdia(2,2)


  END SUBROUTINE EvalPot1_QML_NO3


END MODULE QML_NO3_m
