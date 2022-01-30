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

!> @brief Module which makes the initialization, calculation of the H3 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_H3_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

  ! in the block data, some variable names are not the same as in the SLTH subroutine
  ! => Here, we use the subroutine variable names.

  ! in the block data:
  !       COMMON/VCOM/C,A,A1,F, FN,FB1,FB2,FB3, XN1,XN2,XN3,XN4, B1,B2,B3, G1,G2,G3, D1,D2,D3,D4,XL1,XL2

  ! in the SLTH subroutine:
        !COMMON/VCOM/C,A,A1,F, FNS,F1,F2,F3,   AN1,AN2,AN3,AN4,  B1,B2,B3, W1,W2,W3, D1,D2,D3,D4,XL1,XL2

 ! B3 is a local variable

    real(kind=Rkind), parameter :: C6  = 6.89992032_Rkind
    real(kind=Rkind), parameter :: C8  = 219.9997304_Rkind

    real(kind=Rkind), parameter :: C   = -1.2148730613_Rkind
    real(kind=Rkind), parameter :: A   = -1.514663474_Rkind
    real(kind=Rkind), parameter :: A1  = -1.46_Rkind
    real(kind=Rkind), parameter :: F   = 2.088442_Rkind

 ! different names: FN,XN1,XN2,XN3,XN4 => FNS,AN1,AN2,AN3,AN4
  !       DATA FN,XN1,XN2,XN3,XN4/.0035,.0012646477,-.0001585792,
  !    1 .0000079707,-.0000001151/
    real(kind=Rkind), parameter :: FNS = .0035_Rkind
    real(kind=Rkind), parameter :: AN1 = .0012646477_Rkind
    real(kind=Rkind), parameter :: AN2 = -.0001585792_Rkind
    real(kind=Rkind), parameter :: AN3 = .0000079707_Rkind
    real(kind=Rkind), parameter :: AN4 = -.0000001151_Rkind

 ! different names: FB1,FB2,FB3  => F1,F2,F3
    real(kind=Rkind), parameter :: F1  = .52_Rkind
    real(kind=Rkind), parameter :: B1  = 3.0231771503_Rkind
    real(kind=Rkind), parameter :: B2  = -1.08935219_Rkind

 ! different names: G1,G2,G3  =>  W1,W2,W3
    real(kind=Rkind), parameter :: F2  = .052_Rkind
    real(kind=Rkind), parameter :: W1  = 1.7732141742_Rkind
    real(kind=Rkind), parameter :: W2  = -2.0979468223_Rkind
    real(kind=Rkind), parameter :: W3  = -3.9788502171_Rkind

    real(kind=Rkind), parameter :: D1 = .4908116374_Rkind
    real(kind=Rkind), parameter :: D2 = -.8718696387_Rkind
    real(kind=Rkind), parameter :: D3 = .1612118092_Rkind
    real(kind=Rkind), parameter :: D4 = -.12737311045_Rkind

    real(kind=Rkind), parameter :: F3  = .79_Rkind
    real(kind=Rkind), parameter :: XL2 = .9877930913_Rkind
    real(kind=Rkind), parameter :: XL1 = -13.3599568553_Rkind

    real(kind=Rkind), parameter :: RKW(87)= [                                                  &
      .400000000_Rkind,     .4500000000_Rkind,    .5000000000_Rkind,    .550000000_Rkind,      &
      .600000000_Rkind,     .6500000000_Rkind,    .7000000000_Rkind,    .750000000_Rkind,      &
      .800000000_Rkind,     .9000000000_Rkind,    .1000000000E+01_Rkind,.110000000E+01_Rkind,  &
      .1200000000E+01_Rkind,.1300000000E+01_Rkind,.1350000000E+01_Rkind,.1390000000E+01_Rkind, &
      .1400000000E+01_Rkind,.1401000010E+01_Rkind,.1401099990E+01_Rkind,.1410000000E+01_Rkind, &
      .1450000000E+01_Rkind,.1500000000E+01_Rkind,.1600000000E+01_Rkind,.1700000000E+01_Rkind, &
      .1800000000E+01_Rkind,.1900000000E+01_Rkind,.2000000000E+01_Rkind,.2100000000E+01_Rkind, &
      .2200000000E+01_Rkind,.2300000000E+01_Rkind,.2400000000E+01_Rkind,.2500000000E+01_Rkind, &
      .2600000000E+01_Rkind,.2700000000E+01_Rkind,.2800000000E+01_Rkind,.2900000000E+01_Rkind, &
      .3000000000E+01_Rkind,.3100000000E+01_Rkind,.3200000000E+01_Rkind,.3300000000E+01_Rkind, &
      .3400000000E+01_Rkind,.3500000000E+01_Rkind,.3600000000E+01_Rkind,.3700000000E+01_Rkind, &
      .3800000000E+01_Rkind,.3900000000E+01_Rkind,.4000000000E+01_Rkind,.4100000000E+01_Rkind, &
      .4200000000E+01_Rkind,.4300000000E+01_Rkind,.4400000000E+01_Rkind,.4500000000E+01_Rkind, &
      .4600000000E+01_Rkind,.4700000000E+01_Rkind,.4800000000E+01_Rkind,.4900000000E+01_Rkind, &
      .5000000000E+01_Rkind,.5100000000E+01_Rkind,.5200000000E+01_Rkind,.5300000000E+01_Rkind, &
      .5400000000E+01_Rkind,.5500000000E+01_Rkind,.5600000000E+01_Rkind,.5700000000E+01_Rkind, &
      .5800000000E+01_Rkind,.5900000000E+01_Rkind,.6000000000E+01_Rkind,.6100000000E+01_Rkind, &
      .6200000000E+01_Rkind,.6300000000E+01_Rkind,.6400000000E+01_Rkind,.6500000000E+01_Rkind, &
      .6600000000E+01_Rkind,.6700000000E+01_Rkind,.6800000000E+01_Rkind,.6900000000E+01_Rkind, &
      .7000000000E+01_Rkind,.7200000000E+01_Rkind,.7400000000E+01_Rkind,.7600000000E+01_Rkind, &
      .7800000000E+01_Rkind,.8000000000E+01_Rkind,.8250000000E+01_Rkind,.8500000000E+01_Rkind, &
      .9000000000E+01_Rkind,.9500000000E+01_Rkind,.1000000000E+02_Rkind]

    real(kind=Rkind), parameter :: EKW(87) = [                                                  &
      .879796188_Rkind,      .649071056_Rkind,     .473372447_Rkind,    .337228924_Rkind,       &
      .230365628_Rkind,      .145638432_Rkind,     .779738117E-01_Rkind,.236642733E-01_Rkind,   &
      -.200555771E-01_Rkind,-.836421044E-01_Rkind,-.124538356_Rkind,    -.150056027_Rkind,      &
      -.164934012_Rkind,    -.172345701_Rkind,    -.173962500_Rkind,    -.174451499_Rkind,      &
      -.174474200_Rkind,    -.174474400_Rkind,    -.174474400_Rkind,    -.174459699_Rkind,      &
      -.174055600_Rkind,    -.172853502_Rkind,    -.168579707_Rkind,    -.162456813_Rkind,      &
      -.155066822_Rkind,    -.146849432_Rkind,    -.138131041_Rkind,    -.129156051_Rkind,      &
      -.120123163_Rkind,    -.111172372_Rkind,    -.102412583_Rkind,    -.939271927E-01_Rkind,  &
      -.857809026E-01_Rkind,-.780163108E-01_Rkind,-.706699181E-01_Rkind,-.637640270E-01_Rkind,  &
      -.573117349E-01_Rkind,-.513184414E-01_Rkind,-.457831464E-01_Rkind,-.407002530E-01_Rkind,  &
      -.360577581E-01_Rkind,-.318401624E-01_Rkind,-.280271683E-01_Rkind,-.245977718E-01_Rkind,  &
      -.215296753E-01_Rkind,-.187966785E-01_Rkind,-.163688812E-01_Rkind,-.142246837E-01_Rkind,  &
      -.123370858E-01_Rkind,-.106809878E-01_Rkind,-.923028934E-02_Rkind,-.796819096E-02_Rkind,  &
      -.687029215E-02_Rkind,-.591779314E-02_Rkind,-.509229414E-02_Rkind,-.437819496E-02_Rkind,  &
      -.376259562E-02_Rkind,-.323089623E-02_Rkind,-.277399691E-02_Rkind,-.237999732E-02_Rkind,  &
      -.204229767E-02_Rkind,-.175209799E-02_Rkind,-.150299828E-02_Rkind,-.128989853E-02_Rkind,  &
      -.110689874E-02_Rkind,-.949798920E-03_Rkind,-.814999069E-03_Rkind,-.700199190E-03_Rkind,  &
      -.602999302E-03_Rkind,-.516199400E-03_Rkind,-.446599479E-03_Rkind,-.386399548E-03_Rkind,  &
      -.332799617E-03_Rkind,-.290599668E-03_Rkind,-.246599722E-03_Rkind,-.215399753E-03_Rkind,  &
      -.188899784E-03_Rkind,-.143399836E-03_Rkind,-.108599875E-03_Rkind,-.867998994E-04_Rkind,  &
      -.681999214E-04_Rkind,-.527999393E-04_Rkind,-.403999540E-04_Rkind,-.313999636E-04_Rkind,  &
      -.184999787E-04_Rkind,-.120999861E-04_Rkind,-.909998949E-05_Rkind]

    real(kind=Rkind), parameter :: WKW(87) = [                                                  &
      .308019605E+02_Rkind,  .214419954E+02_Rkind,.154937452E+02_Rkind,.115151545E+02_Rkind,    &
      .871827707E+01_Rkind,  .673831756E+01_Rkind,.527864661E+01_Rkind,.419929947E+01_Rkind,    &
      .333940643E+01_Rkind,  .219403463E+01_Rkind,.149861953E+01_Rkind,.103863661E+01_Rkind,    &
      .730647471_Rkind,      .518552387_Rkind,    .441110777_Rkind,    .383461006_Rkind,        &
      .373946396_Rkind,      .358559402_Rkind,    .372215569_Rkind,    .356670198_Rkind,        &
      .312744133_Rkind,      .261523038_Rkind,    .180817537_Rkind,    .124665543_Rkind,        &
      .807794104E-01_Rkind,  .486562494E-01_Rkind,.251952492E-01_Rkind,.452257820E-02_Rkind,    &
      -.854560161E-02_Rkind,-.196001146E-01_Rkind,-.276538076E-01_Rkind,-.344244662E-01_Rkind,  &
      -.381080935E-01_Rkind,-.421628973E-01_Rkind,-.441600287E-01_Rkind,-.454966841E-01_Rkind,  &
      -.460129217E-01_Rkind,-.458513118E-01_Rkind,-.453815149E-01_Rkind,-.440623159E-01_Rkind,  &
      -.426089183E-01_Rkind,-.404417185E-01_Rkind,-.383839285E-01_Rkind,-.361823035E-01_Rkind,  &
      -.336666088E-01_Rkind,-.302110314E-01_Rkind,-.286090554E-01_Rkind,-.255125522E-01_Rkind,  &
      -.233005599E-01_Rkind,-.201850499E-01_Rkind,-.191990995E-01_Rkind,-.161784216E-01_Rkind,  &
      -.146071006E-01_Rkind,-.126330766E-01_Rkind,-.110605069E-01_Rkind,-.996481997E-02_Rkind,  &
      -.818014482E-02_Rkind,-.765454189E-02_Rkind,-.608163613E-02_Rkind,-.575887028E-02_Rkind,  &
      -.466284400E-02_Rkind,-.408972107E-02_Rkind,-.363824334E-02_Rkind,-.295728079E-02_Rkind,  &
      -.259261281E-02_Rkind,-.221225014E-02_Rkind,-.193837141E-02_Rkind,-.203425060E-02_Rkind,  &
      -.484614204E-03_Rkind,-.226728547E-02_Rkind,-.766232140E-03_Rkind,-.307779418E-03_Rkind,  &
      -.196264565E-02_Rkind,+.131836977E-02_Rkind,-.223083472E-02_Rkind,-.750220030E-04_Rkind,  &
      -.289074004E-03_Rkind,-.220265690E-03_Rkind,-.434861384E-03_Rkind,+.971346041E-05_Rkind,  &
      -.839919101E-04_Rkind,-.153745275E-03_Rkind,-.369227366E-04_Rkind,-.249634065E-04_Rkind,  &
      -.290482724E-04_Rkind,-.148433244E-04_Rkind,+.682166282E-05_Rkind]

  !       BLOCK DATA
	! IMPLICIT REAL*8 (A-H,O-Z)
  !       COMMON/POTCOM/C6,C8,RKW(87),EKW(87),WKW(87)
  !       COMMON/VCOM/C,A,A1,F,FN,FB1,FB2,FB3,XN1,XN2,XN3,XN4,B1,B2,
  !    1  B3,G1,G2,G3,D1,D2,D3,D4,XL1,XL2
  !       DATA C6,C8/6.89992032,219.9997304/
  !       DATA C,A,A1,F/-1.2148730613,-1.514663474,-1.46,2.088442/
  !       DATA FN,XN1,XN2,XN3,XN4/.0035,.0012646477,-.0001585792,
  !    1 .0000079707,-.0000001151/
  !       DATA FB1,B1,B2/.52,3.0231771503,-1.08935219/
  !       DATA FB2,G1,G2,G3/.052,1.7732141742,-2.0979468223,-3.9788502171/
  !       DATA D1,D2,D3,D4/.4908116374,-.8718696387,.1612118092,-.127373
  !    1  11045/
  !       DATA FB3,XL2,XL1/.79,.9877930913,-13.3599568553/
  !       DATA (RKW(I),I=1,40)/
  !    1.400000000E+00,.4500000000E+00,.5000000000E+00,.550000000E+00,
  !    2.600000000E+00,.6500000000E+00,.7000000000E+00,.750000000E+00,
  !    3.800000000E+00,.9000000000E+00,.1000000000E+01,.110000000E+01,
  !    4.1200000000E+01,.1300000000E+01,.1350000000E+01,.1390000000E+01,
  !    5.1400000000E+01,.1401000010E+01,.1401099990E+01,.1410000000E+01,
  !    6.1450000000E+01,.1500000000E+01,.1600000000E+01,.1700000000E+01,
  !    7.1800000000E+01,.1900000000E+01,.2000000000E+01,.2100000000E+01,
  !    8.2200000000E+01,.2300000000E+01,.2400000000E+01,.2500000000E+01,
  !    9.2600000000E+01,.2700000000E+01,.2800000000E+01,.2900000000E+01,
  !    9.3000000000E+01,.3100000000E+01,.3200000000E+01,.3300000000E+01/
  !       DATA(RKW(I),I=41,87)/
  !    1.3400000000E+01,.3500000000E+01,.3600000000E+01,.3700000000E+01,
  !    2.3800000000E+01,.3900000000E+01,.4000000000E+01,.4100000000E+01,
  !    3.4200000000E+01,.4300000000E+01,.4400000000E+01,.4500000000E+01,
  !    4.4600000000E+01,.4700000000E+01,.4800000000E+01,.4900000000E+01,
  !    5.5000000000E+01,.5100000000E+01,.5200000000E+01,.5300000000E+01,
  !    6.5400000000E+01,.5500000000E+01,.5600000000E+01,.5700000000E+01,
  !    7.5800000000E+01,.5900000000E+01,.6000000000E+01,.6100000000E+01,
  !    8.6200000000E+01,.6300000000E+01,.6400000000E+01,.6500000000E+01,
  !    9.6600000000E+01,.6700000000E+01,.6800000000E+01,.6900000000E+01,
  !    9.7000000000E+01,.7200000000E+01,.7400000000E+01,.7600000000E+01,
  !    1.7800000000E+01,.8000000000E+01,.8250000000E+01,.8500000000E+01,
  !    2.9000000000E+01,.9500000000E+01,.1000000000E+02/
  !      DATA(EKW(I),I=1,40)/
  !    1.879796188E+00,.649071056E+00,.473372447E+00,.337228924E+00,
  !    2.230365628E+00,.145638432E+00,.779738117E-01,.236642733E-01,
  !    3-.200555771E-01,-.836421044E-01,-.124538356E+00,-.150056027E+00,
  !    4-.164934012E+00,-.172345701E+00,-.173962500E+00,-.174451499E+00,
  !    5-.174474200E+00,-.174474400E+00,-.174474400E+00,-.174459699E+00,
  !    6-.174055600E+00,-.172853502E+00,-.168579707E+00,-.162456813E+00,
  !    7-.155066822E+00,-.146849432E+00,-.138131041E+00,-.129156051E+00,
  !    8-.120123163E+00,-.111172372E+00,-.102412583E+00,-.939271927E-01,
  !    9-.857809026E-01,-.780163108E-01,-.706699181E-01,-.637640270E-01,
  !    9-.573117349E-01,-.513184414E-01,-.457831464E-01,-.407002530E-01/
  !     DATA(EKW(I),I=41,87)/
  !    1-.360577581E-01,-.318401624E-01,-.280271683E-01,-.245977718E-01,
  !    2-.215296753E-01,-.187966785E-01,-.163688812E-01,-.142246837E-01,
  !    3-.123370858E-01,-.106809878E-01,-.923028934E-02,-.796819096E-02,
  !    4-.687029215E-02,-.591779314E-02,-.509229414E-02,-.437819496E-02,
  !    5-.376259562E-02,-.323089623E-02,-.277399691E-02,-.237999732E-02,
  !    6-.204229767E-02,-.175209799E-02,-.150299828E-02,-.128989853E-02,
  !    7-.110689874E-02,-.949798920E-03,-.814999069E-03,-.700199190E-03,
  !    8-.602999302E-03,-.516199400E-03,-.446599479E-03,-.386399548E-03,
  !    9-.332799617E-03,-.290599668E-03,-.246599722E-03,-.215399753E-03,
  !    9-.188899784E-03,-.143399836E-03,-.108599875E-03,-.867998994E-04,
  !    1-.681999214E-04,-.527999393E-04,-.403999540E-04,-.313999636E-04,
  !    2-.184999787E-04,-.120999861E-04,-.909998949E-05/
  !      DATA(WKW(I),I=1,40)/
  !    1.308019605E+02,.214419954E+02,.154937452E+02,.115151545E+02,
  !    2.871827707E+01,.673831756E+01,.527864661E+01,.419929947E+01,
  !    3.333940643E+01,.219403463E+01,.149861953E+01,.103863661E+01,
  !    4.730647471E+00,.518552387E+00,.441110777E+00,.383461006E+00,
  !    5.373946396E+00,.358559402E+00,.372215569E+00,.356670198E+00,
  !    6.312744133E+00,.261523038E+00,.180817537E+00,.124665543E+00,
  !    7.807794104E-01,.486562494E-01,.251952492E-01,.452257820E-02,
  !    8-.854560161E-02,-.196001146E-01,-.276538076E-01,-.344244662E-01,
  !    9-.381080935E-01,-.421628973E-01,-.441600287E-01,-.454966841E-01,
  !    9-.460129217E-01,-.458513118E-01,-.453815149E-01,-.440623159E-01/
  !     DATA(WKW(I),I=41,87)/
  !    1-.426089183E-01,-.404417185E-01,-.383839285E-01,-.361823035E-01,
  !    2-.336666088E-01,-.302110314E-01,-.286090554E-01,-.255125522E-01,
  !    3-.233005599E-01,-.201850499E-01,-.191990995E-01,-.161784216E-01,
  !    4-.146071006E-01,-.126330766E-01,-.110605069E-01,-.996481997E-02,
  !    5-.818014482E-02,-.765454189E-02,-.608163613E-02,-.575887028E-02,
  !    6-.466284400E-02,-.408972107E-02,-.363824334E-02,-.295728079E-02,
  !    7-.259261281E-02,-.221225014E-02,-.193837141E-02,-.203425060E-02,
  !    8-.484614204E-03,-.226728547E-02,-.766232140E-03,-.307779418E-03,
  !    9-.196264565E-02,+.131836977E-02,-.223083472E-02,-.750220030E-04,
  !    9-.289074004E-03,-.220265690E-03,-.434861384E-03,+.971346041E-05,
  !    1-.839919101E-04,-.153745275E-03,-.369227366E-04,-.249634065E-04,
  !    2-.290482724E-04,-.148433244E-04,+.682166282E-05/
  !       END

!> @brief Derived type in which the H3 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_H3_t

   PRIVATE

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_H3
    PROCEDURE :: Write_QModel     => Write_QML_H3
    PROCEDURE :: Write0_QModel    => Write0_QML_H3
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_QML_H3
    PROCEDURE :: Eval_QModel_Func => EvalFunc_QML_H3
  END TYPE QML_H3_t

  PUBLIC :: QML_H3_t,Init_QML_H3


  CONTAINS
!> @brief Function which makes the initialization of the H3 parameters.
!!
!! @param QModel             TYPE(QML_H3_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_H3(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_H3_t)                           :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H3'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf      = 1
    QModel%ndimCart   = 9

    IF (QModel%option == -1) QModel%option = 0

    SELECT CASE(QModel%option)
    CASE (0) !

      QModel%ndimQ      = 3
      IF (QModel%Cart_TO_Q) THEN
        QModel%ndim       = QModel%ndimCart
      ELSE
        QModel%ndim       = QModel%ndimQ
      END IF

      QModel%pot_name   = 'H3_LSTH'
      QModel%no_ana_der = .TRUE.
    CASE (1)
      QModel%ndimQ    = 1
      QModel%ndim     = 1

      QModel%ndimFunc = 1
      QModel%nb_Func  = 4 ! V, R1op,R2opt,R3opt

      QModel%pot_name   = 'H3_LSTH_IRC'
      QModel%no_ana_der = .FALSE.
    CASE Default
       write(out_unitp,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unitp)
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' option: ',QModel%option
       write(out_unitp,*) ' the possible option are: 0 (3D) or 1 (1D-IRC)'
       STOP 'ERROR in Init_QML_H3: wrong option'
    END SELECT


    IF (debug) write(out_unitp,*) 'init Q0 of H3 (H3 minimum)'
    QModel%Q0 = [2.806_Rkind,2.271_Rkind,2.271_Rkind]

    IF (debug) write(out_unitp,*) 'init d0GGdef of H3'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      !CALL Write_QML_H3(QModel,nio=out_unitp)
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_H3
!> @brief Subroutine wich prints the current QML_H3 parameters.
!!
!! @param QModel            CLASS(QML_H3_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_H3(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_H3_t),   intent(in) :: QModel
    integer,                   intent(in) :: nio

    write(nio,*) 'H3 IV current parameters'

    CALL Write_QML_Empty(QModel%QML_Empty_t,nio)

    write(nio,*) 'end H3 current parameters'

  END SUBROUTINE Write_QML_H3
!> @brief Subroutine wich prints the default QML_H3 parameters.
!!
!! @param QModel            CLASS(QML_H3_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_H3(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_H3_t),   intent(in) :: QModel
    integer,                   intent(in) :: nio

    write(nio,*) 'H3 IV default parameters'
    write(nio,*)
    write(nio,*) 'end H3 IV default parameters'


  END SUBROUTINE Write0_QML_H3

!> @brief Subroutine wich calculates the H3 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_H3_t):    derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H3(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS(QML_H3_t),       intent(in)    :: QModel
    TYPE (dnS_t),          intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),          intent(in)    :: dnQ(:)
    integer,               intent(in)    :: nderiv

    real(kind=Rkind) :: V,Q(3)
    TYPE (dnS_t)    :: Func(QModel%nb_Func)


    SELECT CASE(QModel%option)
    CASE (0) ! 3D
      Q(:) = QML_get_d0_FROM_dnS(dnQ)
      CALL QML_LSTH(Q,V)
      CALL QML_set_dnS(Mat_OF_PotDia(1,1),d0=V)

    CASE (1) ! IRC
      CALL EvalFunc_QML_H3(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)

    END SELECT


  END SUBROUTINE EvalPot_QML_H3

  SUBROUTINE Cart_TO_Q_QML_H3(QModel,dnX,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS(QML_H3_t),       intent(in)    :: QModel
    TYPE (dnS_t),            intent(in)    :: dnX(:,:)
    TYPE (dnS_t),            intent(inout) :: dnQ(:)
    integer,                 intent(in)    :: nderiv

    ! local vector
    integer         :: i,j
    TYPE (dnS_t)    :: Vec12(3),Vec13(3),Vec23(3)


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Cart_TO_Q_QML_H3'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'dnX'
      DO i=1,size(dnX,dim=2)
      DO j=1,size(dnX,dim=1)
        CALL QML_Write_dnS(dnX(j,i),out_unitp)
      END DO
      END DO
      flush(out_unitp)
    END IF

    Vec23(:) = dnX(:,3)-dnX(:,2)
    Vec12(:) = dnX(:,2)-dnX(:,1)
    Vec13(:) = dnX(:,3)-dnX(:,1)
    IF (debug) write(out_unitp,*) 'Cart_TO_Q_QML_H3 vect done'

    dnQ(1) = sqrt(dot_product(Vec23,Vec23))
    dnQ(2) = sqrt(dot_product(Vec12,Vec12))
    dnQ(3) = sqrt(dot_product(Vec13,Vec13))

    CALL QML_dealloc_dnS(Vec23)
    CALL QML_dealloc_dnS(Vec12)
    CALL QML_dealloc_dnS(Vec13)

    IF (debug) THEN
      CALL QML_Write_dnS(dnQ(1),out_unitp,info='dnQ(1)')
      CALL QML_Write_dnS(dnQ(2),out_unitp,info='dnQ(2)')
      CALL QML_Write_dnS(dnQ(3),out_unitp,info='dnQ(3)')
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF
  END SUBROUTINE Cart_TO_Q_QML_H3

  SUBROUTINE QML_LSTH(X,VXD)
  USE QMLLib_NumParameters_m
  IMPLICIT NONE

!       tag='LSTH PES of H2+H. [P. Siegbahn, B. Liu, D.G. Truhlar and C.J.Horowitz, JCP 68, 2457(1978)'
!       X refers to the three hh distances.

    real (kind=Rkind) :: X(3),VXD

    real (kind=Rkind) :: S1(3),S2(3),S3(3)
    real (kind=Rkind) :: EF1,EF2,EF3,X21,X22,X23,T1,T2,T3
    real (kind=Rkind) :: XQ1,XQ2,XQ3,XJ1,XJ2,XJ3,XQ,XJ,ELOND
    real (kind=Rkind) :: WNT,WN,WN2,WN3,WN4,WN5,R,R2,R3
    real (kind=Rkind) :: EXNS,ENS,E,EB1,EB2,EB4,EB1T,EB3T,EB4A,EB4B,EQ,EXF1,EXF2,EXF3,RI,VX
    real (kind=Rkind) :: WB,WB2,WB3,WB4

    real (kind=Rkind) :: COS0,COS1,COS2,COS3
    real (kind=Rkind), parameter :: autTOeV_SLTH = 627.510_Rkind/23.06_Rkind
    real (kind=Rkind), parameter :: V0 = 4.7466_Rkind/autTOeV_SLTH

!  E    LLONDON
!       ***************************************************************
!       TYPE 2
2       FORMAT('ENTERING SLTH')
        EF1=EXP(F*X(1))
        EF2=EXP(F*X(2))
        EF3=EXP(F*X(3))
        X21=X(1)*X(1)
        X22=X(2)*X(2)
        X23=X(3)*X(3)
        T1=C*(A+X(1)+A1*X21)/EF1
        T2=C*(A+X(2)+A1*X22)/EF2
        T3=C*(A+X(3)+A1*X23)/EF3
        CALL QML_VH2(X,S1,S2,S3)
        XQ1=S1(1)+T1
        XQ2=S2(1)+T2
        XQ3=S3(1)+T3
        XJ1=S1(1)-T1
        XJ2=S2(1)-T2
        XJ3=S3(1)-T3
        XQ=(XQ1+XQ2+XQ3)/TWO
        XJ=SQRT(((XJ1-XJ2)**2+(XJ2-XJ3)**2+(XJ3-XJ1)**2)/EIGHT)
        ELOND=XQ-XJ
!       ENS
        WNT=(X(1)-X(2))*(X(2)-X(3))*(X(3)-X(1))
        WN=ABS(WNT)
        WN2=WN*WN
        WN3=WN2*WN
        WN4=WN3*WN
        WN5=WN4*WN
        R=X(1)+X(2)+X(3)
        R2=R*R
        R3=R2*R
!       EXNS=EXP(FNS*R3)
!       ENS=(AN1*WN2+AN2*WN3+AN3*WN4+AN4*WN5)/EXNS
        EXNS=EXP(-FNS*R3)
        ENS=(AN1*WN2+AN2*WN3+AN3*WN4+AN4*WN5)*EXNS
!CCCCC
!CCCCC  NONLINEAR CORRECTIONS
!CCCCC
        COS0=(X21+X22+X23)/TWO
        COS1=(X21-COS0)/X(2)/X(3)
        COS2=(X22-COS0)/X(1)/X(3)
        COS3=(X23-COS0)/X(1)/X(2)
        WB=COS1+COS2+COS3+ONE
        WB2=WB*WB
        WB3=WB2*WB
        WB4=WB3*WB
        EXF1=EXP(-F1*R)
        EXF2=EXP(-F2*R2)
        EXF3=EXP(-F3*R)
        EB1T=(B1+B2*R)*EXF1
        EB3T=(XL1+XL2*R2)*EXF3
        EB1=WB*(EB1T+EB3T)
        EB2=(WB2*W1+WB3*W2+WB4*W3)*EXF2
        EQ=(X(1)-X(2))**2+(X(2)-X(3))**2+(X(3)-X(1))**2
        RI=ONE/X(1)+ONE/X(2)+ONE/X(3)
        EB4A=WB*D1*EXF1+WB2*D2*EXF2
        EB4B=D3*EXF1+D4*EXF2
        EB4=EB4A*RI+EB4B*WB*EQ
        E=ELOND+ENS+EB1+EB2+EB4
        VXD = E ! in Hartree
        !E=E*627.510
        !VX=E/23.06+4.7466
        !VXD=DBLE(VX)
        !VXD = (E+V0)*autTOeV_SLTH
  END SUBROUTINE QML_LSTH

  SUBROUTINE QML_VH2(X,S1,S2,S3)
  USE QMLLib_NumParameters_m
  IMPLICIT NONE
        real(kind=Rkind), intent(in)    :: X(3)
        real(kind=Rkind), intent(inout) :: S1(3),S2(3),S3(3)

        !COMMON/POTCOM/C6,C8,RKW(87),EKW(87),WKW(87)
        IF (X(1) > TEN) THEN
          CALL QML_VBIGR(X(1),S1)
        ELSE
          CALL QML_SPLID2(87,RKW,EKW,WKW,1,X(1),S1)
        END IF
        IF (X(2) > TEN) THEN
          CALL QML_VBIGR(X(2),S2)
        ELSE
          CALL QML_SPLID2(87,RKW,EKW,WKW,1,X(2),S2)
        END IF

        IF(X(3) > TEN) THEN
          CALL QML_VBIGR(X(3),S3)
        ELSE
          CALL QML_SPLID2(87,RKW,EKW,WKW,1,X(3),S3)
        END IF
  END SUBROUTINE QML_VH2
!       *************************************************************
  SUBROUTINE QML_VBIGR(X,S)
  USE QMLLib_NumParameters_m
  IMPLICIT NONE
        real(kind=Rkind), intent(in)    ::  X
        real(kind=Rkind), intent(inout) ::  S(3)

        real(kind=Rkind)   ::  X2,X3,X6,C8A

        !COMMON/POTCOM/C6,C8,RKW(87),EKW(87),WKW(87)

        X2=X*X
        X3=X2*X
        X6=X3*X3
        C8A=C8/X2
        S(1)=-(C6+C8A)/X6
        S(2)=(C6*SIX + C8A*EIGHT)/X6/X

  END SUBROUTINE QML_VBIGR

  SUBROUTINE QML_SPLID2(N,X,F,W,IJ,Y,TAB)
  USE QMLLib_NumParameters_m
  IMPLICIT NONE

        integer,          intent(in)    :: N,IJ
        real(kind=Rkind), intent(in)    :: X(N)
        real(kind=Rkind), intent(in)    :: Y,F(N),W(N)
        real(kind=Rkind), intent(inout) :: TAB(3)

        INTEGER          :: I,K,MI,KI
        real(kind=Rkind) :: FLK,AA,BB,CC

!         IF(Y-X(1))10,10,20
! 10      I=1
!         GO TO 30
! 20      IF(Y-X(N))15,40,40
! 40      I=N-1
!         GO TO 30
! 15      I=0
!         DO K=1,N
!           IF(X(K) > Y) EXIT
!           I=I+1
!         END DO
! 30      MI=(I-1)*IJ+1

        IF(Y-X(1)> ZERO )THEN ! label 20
          IF (Y-X(N) < ZERO) THEN ! label 15
            I=0
            DO K=1,N
              IF(X(K) > Y) EXIT
              I=I+1
            END DO
          ELSE ! label 40
            I=N-1
          END IF
        ELSE ! <= 0 label 10 + goto 30
          I = 1
        END IF
        MI=(I-1)*IJ+1

        KI=MI+IJ
        FLK=X(I+1)-X(I)
        AA=(W(MI)*(X(I+1)-Y)**3 + W(KI)*(Y-X(I))**3)/(SIX*FLK)
        BB=(F(KI)/FLK-W(KI)*FLK/SIX)*(Y-X(I))
        CC=(F(MI)/FLK-FLK*W(MI)/SIX)*(X(I+1)-Y)
        TAB(1)=AA+BB+CC
        AA=(W(KI)*(Y-X(I))**2-W(MI)*(X(I+1)-Y)**2)/(TWO*FLK)
        BB=(F(KI)-F(MI))/FLK
        CC=FLK*(W(MI)-W(KI))/SIX
        TAB(2)=AA+BB+CC
        TAB(3)=(W(MI)*(X(I+1)-Y)+W(KI)*(Y-X(I)))/FLK
  END SUBROUTINE QML_SPLID2

  ! for IRC
  SUBROUTINE EvalFunc_QML_H3(QModel,Func,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  USE QMLdnSVM_dnPoly_m
  IMPLICIT NONE

    CLASS(QML_H3_t),      intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    TYPE (dnS_t)    :: s,ts
    integer         :: i
    integer, parameter :: max_deg = 16
    TYPE (dnS_t)       :: tab_Pl(0:max_deg)


    s  = dnQ(1)
    ts = tanh(s)

    DO i=0,max_deg
      tab_Pl(i) = QML_dnLegendre0(ts,i,ReNorm=.FALSE.)
    END DO

    ! potential
    Func(1) = -0.17447440045043028_Rkind                   +                    &
              (0.01085897589385337_Rkind)     * tab_Pl(0)  +                    &
              (-0.00966749217922542_Rkind)    * tab_Pl(2)  +                    &
              (-0.000608958364909583_Rkind)   * tab_Pl(4)  +                    &
              (-0.0005463286353325903_Rkind)  * tab_Pl(6)  +                    &
              (-7.674511378388051e-05_Rkind)  * tab_Pl(8)  +                    &
              (-1.2703773915173453e-06_Rkind) * tab_Pl(10) +                    &
              (1.758821635592407e-05_Rkind)   * tab_Pl(12) +                    &
              (1.705154432909527e-05_Rkind)   * tab_Pl(14) +                    &
              (-3.748171912010247e-05_Rkind)  * tab_Pl(16)
   ! R1eq
   Func(2) = 1.4010444_Rkind * QML_dnSigmoid_H3(-s,ONE)           +             &
             QML_dnSigmoid_H3(s,ONE)*(s+1.6803643117748104_Rkind) +             &
            (0.135717682405616_Rkind)      * tab_Pl(0)            +             &
            (-0.0009647928498475703_Rkind) * tab_Pl(1)            +             &
            (-0.14960721942361688_Rkind)   * tab_Pl(2)            +             &
            (-0.004622944008638028_Rkind)  * tab_Pl(3)            +             &
            (0.013899404906010182_Rkind)   * tab_Pl(4)            +             &
            (0.004063799184701145_Rkind)   * tab_Pl(5)            +             &
            (-0.001052823208516841_Rkind)  * tab_Pl(6)            +             &
            (0.0006881469085736082_Rkind)  * tab_Pl(7)            +             &
            (0.0009920029080177583_Rkind)  * tab_Pl(8)            +             &
            (0.00014255050945269853_Rkind) * tab_Pl(9)            +             &
            (2.495474739587226e-05_Rkind)  * tab_Pl(10)           +             &
            (0.0004543931972376581_Rkind)  * tab_Pl(11)

   ! R2eq(s) = R1eq(-s)
   Func(3) = 1.4010444_Rkind * QML_dnSigmoid_H3(s,ONE)              +           &
             QML_dnSigmoid_H3(-s,ONE)*(-s+1.6803643117748104_Rkind) +           &
            (0.135717682405616_Rkind)      * tab_Pl(0)              +           &
            (0.0009647928498475703_Rkind)  * tab_Pl(1)              +           &
            (-0.14960721942361688_Rkind)   * tab_Pl(2)              +           &
            (0.004622944008638028_Rkind)   * tab_Pl(3)              +           &
            (0.013899404906010182_Rkind)   * tab_Pl(4)              +           &
            (-0.004063799184701145_Rkind)  * tab_Pl(5)              +           &
            (-0.001052823208516841_Rkind)  * tab_Pl(6)              +           &
            (-0.0006881469085736082_Rkind) * tab_Pl(7)              +           &
            (0.0009920029080177583_Rkind)  * tab_Pl(8)              +           &
            (-0.00014255050945269853_Rkind)* tab_Pl(9)              +           &
            (2.495474739587226e-05_Rkind)  * tab_Pl(10)             +           &
            (-0.0004543931972376581_Rkind) * tab_Pl(11)

    ! R3eq(s) = R1eq(s)+R2eq(s) (!linear HHH)
    Func(4) = Func(2) + Func(3)

    DO i=0,max_deg
      CALL QML_dealloc_dnS(tab_Pl(i))
    END DO
    CALL QML_dealloc_dnS(s)
    CALL QML_dealloc_dnS(ts)

  END SUBROUTINE EvalFunc_QML_H3
  FUNCTION QML_dnSigmoid_H3(x,a)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    TYPE (dnS_t)                        :: QML_dnSigmoid_H3

    TYPE (dnS_t),         intent(in)    :: x
    real(kind=Rkind),     INTENT(IN)    :: a

    QML_dnSigmoid_H3 = (ONE+tanh(a*x))/TWO

  END FUNCTION QML_dnSigmoid_H3
END MODULE QML_H3_m
