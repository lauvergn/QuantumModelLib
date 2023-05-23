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
!> @brief Module which makes the initialization, calculation of the H3 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_H3_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
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
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_H3
    PROCEDURE :: Write_QModel     => Write_QML_H3
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_QML_H3
    PROCEDURE :: EvalFunc_QModel => EvalFunc_QML_H3
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
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_H3_t)                              :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    real(kind=Rkind), parameter :: mH = 1837.1526464003414_Rkind

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H3'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf      = 1
    QModel%ndimCart   = 9

    allocate(QModel%masses(QModel%ndimCart/3))
    QModel%masses(:) = [mH,mH,mH]

    IF (QModel%option == -1) QModel%option = 0

    SELECT CASE(QModel%option)
    CASE (0,10,20) !
      IF (QModel%ndim == 2) THEN
        QModel%ndimQ      = 2
      ELSE
        QModel%ndimQ      = 3
      END IF
      IF (QModel%Cart_TO_Q) THEN
        QModel%ndim       = QModel%ndimCart
      ELSE
        QModel%ndim       = QModel%ndimQ
      END IF

      ! so that, we can get the irc function with the 3D model
      QModel%ndimFunc = 1
      QModel%nb_Func  = 6 ! V, R1op,R2opt,R3opt,RPH_a0,RPH_a1

      QModel%pot_name   = 'H3_LSTH'
      QModel%no_ana_der = .TRUE.
    CASE (1,11,21)
      QModel%ndimQ    = 1
      QModel%ndim     = 1

      QModel%ndimFunc = 1
      QModel%nb_Func  = 6 ! V, R1op,R2opt,R3opt,RPH_a0,RPH_a1

      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 2
      QModel%IndexFunc_Grad = 0
      QModel%IndexFunc_Hess = 0

      QModel%pot_name   = 'H3_LSTH_IRC'
      QModel%no_ana_der = .FALSE.

    CASE (31)
      QModel%ndimQ    = 1
      QModel%ndim     = 1

      QModel%ndimFunc = 1
      QModel%nb_Func  = 3 ! V, Hessian, Rho 

      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 3
      QModel%IndexFunc_Grad = 0
      QModel%IndexFunc_Hess = 2

      QModel%pot_name   = 'h3_jo'
      QModel%no_ana_der = .FALSE.

    CASE (3) ! 2d avec le rho fitté de 31
      IF (QModel%ndim == 2) THEN
        QModel%ndimQ      = 2
      ELSE
        QModel%ndimQ      = 3
      END IF
      IF (QModel%Cart_TO_Q) THEN
        QModel%ndim       = QModel%ndimCart
      ELSE
        QModel%ndim       = QModel%ndimQ
      END IF

      QModel%ndimFunc = 1
      QModel%nb_Func  = 3 ! V, Hessian, Rho 

      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 3
      QModel%IndexFunc_Grad = 0
      QModel%IndexFunc_Hess = 2

      QModel%pot_name   = 'H3_LSTH'
      QModel%no_ana_der = .TRUE.

    CASE (41) !CASE 31 avec transfo de coordonnée sinh(s)
      QModel%ndimQ    = 1
      QModel%ndim     = 1

      QModel%ndimFunc = 1
      QModel%nb_Func  = 3 ! V, Hessian, Rho 

      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 3
      QModel%IndexFunc_Grad = 0
      QModel%IndexFunc_Hess = 2

      QModel%pot_name   = 'h3_jo'
      QModel%no_ana_der = .FALSE.

    CASE (4) ! 2d avec le rho fitté de 41 (sinh(s))
        IF (QModel%ndim == 2) THEN
          QModel%ndimQ      = 2
        ELSE
          QModel%ndimQ      = 3
        END IF
        IF (QModel%Cart_TO_Q) THEN
          QModel%ndim       = QModel%ndimCart
        ELSE
          QModel%ndim       = QModel%ndimQ
        END IF
  
        QModel%ndimFunc = 1
        QModel%nb_Func  = 3 ! V, Hessian, Rho 
  
        QModel%IndexFunc_Ene  = 1
        QModel%IndexFunc_Qop  = 3
        QModel%IndexFunc_Grad = 0
        QModel%IndexFunc_Hess = 2
  
        QModel%pot_name   = 'H3_LSTH'
        QModel%no_ana_der = .TRUE.

    CASE (51) !CASE 51 avec transfo 2D
      QModel%ndimQ    = 1
      QModel%ndim     = 1

      QModel%ndimFunc = 1
      QModel%nb_Func  = 3 ! V, Hessian, Rho 

      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 3
      QModel%IndexFunc_Grad = 0
      QModel%IndexFunc_Hess = 2

      QModel%pot_name   = 'h3_jo'
      QModel%no_ana_der = .FALSE.

    CASE (5) ! 2d avec le rho fitté de 51 (transfo2D)
        IF (QModel%ndim == 2) THEN
          QModel%ndimQ      = 2
        ELSE
          QModel%ndimQ      = 3
        END IF
        IF (QModel%Cart_TO_Q) THEN
          QModel%ndim       = QModel%ndimCart
        ELSE
          QModel%ndim       = QModel%ndimQ
        END IF
  
        QModel%ndimFunc = 1
        QModel%nb_Func  = 3 ! V, Hessian, Rho 
  
        QModel%IndexFunc_Ene  = 1
        QModel%IndexFunc_Qop  = 3
        QModel%IndexFunc_Grad = 0
        QModel%IndexFunc_Hess = 2
  
        QModel%pot_name   = 'H3_LSTH'
        QModel%no_ana_der = .TRUE.

    CASE (61) !CASE 61 VO + ZPE avec transfo 2D
      QModel%ndimQ    = 1
      QModel%ndim     = 1

      QModel%ndimFunc = 1
      QModel%nb_Func  = 4 ! V, Hessian, Rho, ZPE 

      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 3
      QModel%IndexFunc_Grad = 0
      QModel%IndexFunc_Hess = 2

      QModel%pot_name   = 'h3_jo'
      QModel%no_ana_der = .FALSE.

    CASE (72) !CASE 72 VO + 0.5*V''(fit le long de rho puis s)*drho^2
      IF (QModel%ndim == 2) THEN
         QModel%ndimQ      = 2
      ELSE
        QModel%ndimQ      = 3
      END IF
      IF (QModel%Cart_TO_Q) THEN
        QModel%ndim       = QModel%ndimCart
      ELSE
        QModel%ndim       = QModel%ndimQ
      END IF

      QModel%ndimQ    = 2
      QModel%ndim     = 2


      QModel%ndimFunc = 1
      QModel%nb_Func  = 4 ! V, 0.5*Hessian, rhoOpt, V'  

      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 3
      QModel%IndexFunc_Grad = 0
      QModel%IndexFunc_Hess = 2

      QModel%pot_name   = 'h3_LSTH' !LSTH ou jo ?
      QModel%no_ana_der = .FALSE.

    CASE (74) !CASE 74 VO + 0.5*V''*dro^2 + 1/6*V'''*dro^3+1/24*V''''*dro^4(fit le long de rho puis s)
      IF (QModel%ndim == 2) THEN
         QModel%ndimQ      = 2
      ELSE
        QModel%ndimQ      = 3
      END IF
      IF (QModel%Cart_TO_Q) THEN
        QModel%ndim       = QModel%ndimCart
      ELSE
        QModel%ndim       = QModel%ndimQ
      END IF

      QModel%ndimQ    = 2
      QModel%ndim     = 2


      QModel%ndimFunc = 1
      QModel%nb_Func  = 6 ! V, Hessian, rhoOpt, V', V''', V''''  

      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 3
      QModel%IndexFunc_Grad = 0
      QModel%IndexFunc_Hess = 2

      QModel%pot_name   = 'h3_LSTH' !LSTH ou jo ?
      QModel%no_ana_der = .FALSE.


    CASE Default
       write(out_unit,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unit)
       write(out_unit,*) ' ERROR in ',name_sub
       write(out_unit,*) ' option: ',QModel%option
       write(out_unit,*) ' the possible option are: 0,10,20 (3D) or 1,11,21 (1D-IRC)'
       STOP 'ERROR in Init_QML_H3: wrong option'
    END SELECT

    IF (QModel%AbInitio) QModel%no_ana_der = .FALSE.


    IF (debug) write(out_unit,*) 'init Q0 of H3 (H3 minimum)'
    QModel%Q0 = [2.806_Rkind,2.271_Rkind,2.271_Rkind]

    IF (debug) write(out_unit,*) 'init d0GGdef of H3'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)

    IF (debug) THEN
      CALL Write_QML_H3(QModel,nio=out_unit)
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_H3
!> @brief Subroutine wich prints the current QML_H3 parameters.
!!
!! @param QModel            CLASS(QML_H3_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_H3(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_H3_t),   intent(in) :: QModel
    integer,           intent(in) :: nio

    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) 'H3 LSTH current parameters'
    write(nio,*)
    write(nio,*) ' Units: Energy in Hartree and distances in bohr'
    write(nio,*) 'refs: '
    write(nio,*) ' P. Siegbahn, B. Liu,  J. Chem. Phys. 68, 2457(1978).'
    write(nio,*) ' D.G. Truhlar and C.J. Horowitz, J. Chem. Phys. 68, 2466 (1978); https://doi.org/10.1063/1.436019'
    write(nio,*)
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*)

    SELECT CASE(QModel%option)
    CASE (0) ! 2D or 3D
      IF (QModel%ndim == 2) write(nio,*) 'Linear 2D-H3 LSTH model, with 2 distances'
      IF (QModel%ndim == 3) write(nio,*) '3D-H3 LSTH model, with 3 distances'
      IF (QModel%ndim == 9) write(nio,*) 'Cartessian H3 LSTH model'

      write(nio,*) 'Can be used with the first 1D-IRC H3 LSTH model (from polar tranformation)'

    CASE (10) ! 2D or 3D
      IF (QModel%ndim == 2) write(nio,*) 'Linear 2D-H3 LSTH model, with 2 distances'
      IF (QModel%ndim == 3) write(nio,*) '3D-H3 LSTH model, with 3 distances'
      IF (QModel%ndim == 9) write(nio,*) 'Cartessian LSTH LST model'

      write(nio,*) 'Can be used with the second 1D-IRC H3 LSTH model (from sum and difference)'

    CASE (20) ! 2D or 3D
      IF (QModel%ndim == 2) write(nio,*) 'Linear 2D-H3 LSTH model, with 2 distances'
      IF (QModel%ndim == 3) write(nio,*) '3D-H3 LSTH model, with 3 distances'
      IF (QModel%ndim == 9) write(nio,*) 'Cartessian LSTH LST model'

      write(nio,*) 'Can be used with the third (correct?) 1D-IRC H3 LSTH model'
    CASE (1) ! IRC
      write(nio,*) 'First 1D-IRC H3 LSTH model (from polar transformation)'
    CASE (11) ! IRC
      write(nio,*) 'Second 1D-IRC H3 LSTH model (from sum and difference)'
    CASE (21) ! IRC
      write(nio,*) 'Third (correct?) 1D-IRC H3 LSTH model'

    CASE (31) 
            write(nio,*) 'Case 31 for Jo'
    CASE (3) ! 2D 
        IF (QModel%ndim == 2) write(nio,*) 'Linear 2D-H3 LSTH model, with 2 distances, FOR Jo case 3'
        IF (QModel%ndim == 3) write(nio,*) '3D-H3 LSTH model, with 3 distances'
        IF (QModel%ndim == 9) write(nio,*) 'Cartessian H3 LSTH model'
    CASE (41)
            write(nio,*) 'Case 41 for Jo : 31 + transfo ts=sinh(s)'
    CASE (4) ! 2D 
          IF (QModel%ndim == 2) write(nio,*) 'Linear 2D-H3 LSTH model, with 2 distances, FOR Jo case 4 ts=sinh(s)'
          IF (QModel%ndim == 3) write(nio,*) '3D-H3 LSTH model, with 3 distances'
          IF (QModel%ndim == 9) write(nio,*) 'Cartessian H3 LSTH model'

    CASE (51)
            write(nio,*) 'Case 51 for Jo : 51 + transfo 2D rho.s'
    CASE (5) ! 2D 
          IF (QModel%ndim == 2) write(nio,*) 'Linear 2D-H3 LSTH model, with 2 distances, FOR Jo case 5 transfo 2D rho.s'
          IF (QModel%ndim == 3) write(nio,*) '3D-H3 LSTH model, with 3 distances'
          IF (QModel%ndim == 9) write(nio,*) 'Cartessian H3 LSTH model'

    CASE (61)
            write(nio,*) 'Case 61 for Jo : 61 V0+ZPE 1D'

    CASE (72)
            write(nio,*) 'Case 72 for Jo : 72 V0+harmonique(fit)'

    CASE (74)
            write(nio,*) 'Case 74 for Jo : 74 V0+ harmo+ anharmonique ordre 4(fit)'
    END SELECT
    write(nio,*)

    write(nio,*)
    write(nio,*) '-------------------------------------------------------------'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*)
    write(nio,*) 'end H3 LSTH current parameters'
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) '-------------------------------------------------------------'

  END SUBROUTINE Write_QML_H3

!> @brief Subroutine wich calculates the H3 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_H3_t):    derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H3(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_H3_t),       intent(in)    :: QModel
    TYPE (dnS_t),          intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),          intent(in)    :: dnQ(:)
    integer,               intent(in)    :: nderiv

    real(kind=Rkind) :: V,Q(3)
    TYPE (dnS_t)    :: Func(QModel%nb_Func)
    TYPE (dnS_t)    :: deltaRho

    !write(*,*) "In EvalPot", "with option", QModel%option
    !write(*,*) "dnQ", get_d0(dnQ) 
    SELECT CASE(QModel%option)
    CASE (0,10,20,3,4,5) ! 2D or 3D
      IF (size(dnQ) == 2) THEN
        Q(1:2) = get_d0(dnQ)
        Q(3)   = Q(1) + Q(2)
      ELSE
        Q(:) = get_d0(dnQ)
      END IF
      CALL QML_LSTH_refactoring(Q,V)
      CALL set_dnS(Mat_OF_PotDia(1,1),d0=V)

    CASE (1) ! IRC
      CALL EvalFunc_QML_H3_v1(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)

    CASE (11) ! IRC
      CALL EvalFunc_QML_H3_v11(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)

    CASE (21) ! IRC
      CALL EvalFunc_QML_H3_v11(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)

    CASE (31)
      CALL EvalFunc_QML_H3_v31(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)

    CASE (41) 
      CALL EvalFunc_QML_H3_v41(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)

    CASE (51) 
      CALL EvalFunc_QML_H3_v51(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)

    CASE (61) 
      CALL EvalFunc_QML_H3_v61(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1) + Func(4)

    CASE (72) 
      CALL EvalFunc_QML_H3_v72(QModel,Func,[dnQ(1)],nderiv)
      deltaRho=(dnQ(2)-Func(3))
      !Mat_OF_PotDia(1,1) = Func(1) + Func(4)*(dnQ(2)-Func(3)) + HALF*Func(2)*(dnQ(2)-Func(3))*(dnQ(2)-Func(3))
      !Mat_OF_PotDia(1,1) = Func(1) + Func(4)*(dnQ(2)-Func(3)) + Func(2)*(dnQ(2)-Func(3))*(dnQ(2)-Func(3))
      Mat_OF_PotDia(1,1) = Func(1) + Func(4)*deltaRho + Func(2)*deltaRho*deltaRho

    CASE (74) 
      CALL EvalFunc_QML_H3_v74(QModel,Func,[dnQ(1)],nderiv)
      deltaRho=(dnQ(2)-Func(3))
      !Mat_OF_PotDia(1,1) = Func(1) + Func(4)*(dnQ(2)-Func(3)) + HALF*Func(2)*(dnQ(2)-Func(3))*(dnQ(2)-Func(3))  &
      !+ SIXTH*Func(5)*(dnQ(2)-Func(3))*(dnQ(2)-Func(3))*(dnQ(2)-Func(3)) + QUARTER*THIRD*HALF*Func(6) &
      !* (dnQ(2)-Func(3))*(dnQ(2)-Func(3))*(dnQ(2)-Func(3))*(dnQ(2)-Func(3))
      Mat_OF_PotDia(1,1) = Func(1) + Func(4)*deltaRho + Func(2)*deltaRho*deltaRho  &
      + Func(5)*deltaRho*deltaRho*deltaRho + Func(6)*deltaRho*deltaRho*deltaRho*deltaRho
      !Mat_OF_PotDia(1,1) = Func(1) + Func(4)*(dnQ(2)-Func(3)) + Func(2)*((dnQ(2)-Func(3))*(dnQ(2)-Func(3)))  &
      !+ Func(5)*((dnQ(2)-Func(3))*(dnQ(2)-Func(3))*(dnQ(2)-Func(3))) + Func(6) &
      !* ((dnQ(2)-Func(3))*(dnQ(2)-Func(3))*(dnQ(2)-Func(3))*(dnQ(2)-Func(3)))
    END SELECT
    !write(*,*) "Mat_pot", get_d0(Mat_OF_PotDia) 
  END SUBROUTINE EvalPot_QML_H3

  SUBROUTINE Cart_TO_Q_QML_H3(QModel,dnX,dnQ,nderiv)
    USE QDUtil_m,         ONLY : TO_string
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H3_t),         intent(in)    :: QModel
    TYPE (dnS_t),            intent(in)    :: dnX(:,:)
    TYPE (dnS_t),            intent(inout) :: dnQ(:)
    integer,                 intent(in)    :: nderiv

    ! local vectors
    integer                   :: i,j
    real (kind=Rkind)         :: sm1,sm2,sm3
    TYPE (dnS_t), allocatable :: Vec12(:),Vec13(:),Vec23(:)

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Cart_TO_Q_QML_H3'
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

    allocate(Vec23(3))
    allocate(Vec12(3))
    allocate(Vec13(3))
    IF (QModel%MassWeighted) THEN
      ! here dnX are mass weighted, we have to un-mass weighted the vectors
      sm1 = ONE/sqrt(QModel%masses(1))
      sm2 = ONE/sqrt(QModel%masses(2))
      sm3 = ONE/sqrt(QModel%masses(3))

      Vec23(:) = dnX(:,3)*sm3-dnX(:,2)*sm2
      Vec12(:) = dnX(:,2)*sm2-dnX(:,1)*sm1
      Vec13(:) = dnX(:,3)*sm3-dnX(:,1)*sm1
    ELSE
      Vec23(:) = dnX(:,3)-dnX(:,2)
      Vec12(:) = dnX(:,2)-dnX(:,1)
      Vec13(:) = dnX(:,3)-dnX(:,1)
    END IF

    IF (debug) THEN
      write(out_unit,*) 'Cart_TO_Q_QML_H3 vect done'
      flush(out_unit)
      DO j=1,size(Vec23,dim=1)
        CALL Write_dnS(Vec23(j),out_unit,info='Vec23')
      END DO
      DO j=1,size(Vec12,dim=1)
        CALL Write_dnS(Vec23(j),out_unit,info='Vec12')
      END DO
      DO j=1,size(Vec23,dim=1)
        CALL Write_dnS(Vec13(j),out_unit,info='Vec13')
      END DO
      flush(out_unit)
    END IF

    dnQ(1) = sqrt(dot_product(Vec23,Vec23))
    dnQ(2) = sqrt(dot_product(Vec12,Vec12))
    dnQ(3) = sqrt(dot_product(Vec13,Vec13))

    CALL dealloc_dnS(Vec23)
    CALL dealloc_dnS(Vec12)
    CALL dealloc_dnS(Vec13)

    IF (debug) THEN
      CALL Write_dnS(dnQ(1),out_unit,info='dnQ(1)')
      CALL Write_dnS(dnQ(2),out_unit,info='dnQ(2)')
      CALL Write_dnS(dnQ(3),out_unit,info='dnQ(3)')
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
  END SUBROUTINE Cart_TO_Q_QML_H3

  SUBROUTINE QML_LSTH_refactoring(X,VXD)
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
        EF1=EXP(-F*X(1))   ! change in EXP(-F*X(:)) instead of EXP(F*X(:)) to avoid a possible overflow
        EF2=EXP(-F*X(2))
        EF3=EXP(-F*X(3))
        X21=X(1)*X(1)
        X22=X(2)*X(2)
        X23=X(3)*X(3)
        T1=C*(A+X(1)+A1*X21)*EF1 ! change the / by * to avoid a possible overflow (see above)
        T2=C*(A+X(2)+A1*X22)*EF2
        T3=C*(A+X(3)+A1*X23)*EF3
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
  END SUBROUTINE QML_LSTH_refactoring

  SUBROUTINE QML_LSTH(X,VXD)
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

  SUBROUTINE QML_SPLID2(N,X,FF,W,IJ,Y,TAB)
    IMPLICIT NONE

        integer,          intent(in)    :: N,IJ
        real(kind=Rkind), intent(in)    :: X(N)
        real(kind=Rkind), intent(in)    :: Y,FF(N),W(N)
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
        BB=(FF(KI)/FLK-W(KI)*FLK/SIX)*(Y-X(I))
        CC=(FF(MI)/FLK-FLK*W(MI)/SIX)*(X(I+1)-Y)
        TAB(1)=AA+BB+CC
        AA=(W(KI)*(Y-X(I))**2-W(MI)*(X(I+1)-Y)**2)/(TWO*FLK)
        BB=(FF(KI)-FF(MI))/FLK
        CC=FLK*(W(MI)-W(KI))/SIX
        TAB(2)=AA+BB+CC
        TAB(3)=(W(MI)*(X(I+1)-Y)+W(KI)*(Y-X(I)))/FLK
  END SUBROUTINE QML_SPLID2

  SUBROUTINE EvalFunc_QML_H3(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H3_t),      intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    SELECT CASE(QModel%option)
    CASE (0,1) ! IRC
      CALL EvalFunc_QML_H3_v1(QModel,Func,dnQ,nderiv)
    CASE (10,11) ! IRC
      CALL EvalFunc_QML_H3_v11(QModel,Func,dnQ,nderiv)
    CASE (20,21) ! IRC
      CALL EvalFunc_QML_H3_v21(QModel,Func,dnQ,nderiv)
    CASE (31,3) 
      CALL EvalFunc_QML_H3_v31(QModel,Func,dnQ,nderiv)
    CASE (41,4) 
      CALL EvalFunc_QML_H3_v41(QModel,Func,dnQ,nderiv)
    CASE (51,5) 
      CALL EvalFunc_QML_H3_v51(QModel,Func,dnQ,nderiv)
    CASE (61) 
      CALL EvalFunc_QML_H3_v61(QModel,Func,dnQ,nderiv)
    CASE (72) 
      CALL EvalFunc_QML_H3_v72(QModel,Func,dnQ,nderiv)
    CASE (74) 
      CALL EvalFunc_QML_H3_v74(QModel,Func,dnQ,nderiv)
    END SELECT
  END SUBROUTINE EvalFunc_QML_H3

  ! for IRC, RPH
  SUBROUTINE EvalFunc_QML_H3_v11(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H3_t),      intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    TYPE (dnS_t)                  :: s,ts,am,ap
    integer                       :: i
    integer,           parameter  :: max_deg = 20
    TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
    real (kind=Rkind), parameter  :: betaQ = 1.4_Rkind


    s  = dnQ(1)
    ts = tanh(s)

    DO i=0,max_deg
      tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
    END DO

    ! potential
    Func(1) = -0.17447440045043028_Rkind +                                      &
              (0.010858975893853413_Rkind)    * tab_Pl(0)  +                    &
              (-0.009667492179225417_Rkind)   * tab_Pl(2)  +                    &
              (-0.0006089583649096594_Rkind)  * tab_Pl(4)  +                    &
              (-0.0005463286353325458_Rkind)  * tab_Pl(6)  +                    &
              (-7.674511378387464e-05_Rkind)  * tab_Pl(8)  +                    &
              (-1.2703773915447383e-06_Rkind) * tab_Pl(10) +                    &
              (1.7588216355925958e-05_Rkind)  * tab_Pl(12) +                    &
              (1.7051544329081514e-05_Rkind)  * tab_Pl(14) +                    &
              (-3.7481719120090565e-05_Rkind) * tab_Pl(16)

   ! R1eq
   Func(2) = 1.4010444_Rkind * QML_dnSigmoid_H3(-s,betaQ)         +             &
           QML_dnSigmoid_H3(s,betaQ)*(s+1.6803643117748104_Rkind) +             &
            (0.1357176824058676_Rkind)     * tab_Pl(0)            +             &
            (-0.0009647928495994389_Rkind) * tab_Pl(1)            +             &
            (-0.1496072194236861_Rkind)    * tab_Pl(2)            +             &
            (-0.004622944008685905_Rkind)  * tab_Pl(3)            +             &
            (0.013899404906043201_Rkind)   * tab_Pl(4)            +             &
            (0.00406379918470803_Rkind)    * tab_Pl(5)            +             &
            (-0.0010528232085257284_Rkind) * tab_Pl(6)            +             &
            (0.0006881469085679271_Rkind)  * tab_Pl(7)            +             &
            (0.0009920029080101385_Rkind)  * tab_Pl(8)            +             &
            (0.00014255050945681492_Rkind) * tab_Pl(9)            +             &
            (2.4954747402158953e-05_Rkind) * tab_Pl(10)           +             &
            (0.0004543931972375909_Rkind)  * tab_Pl(11)

   ! R2eq(s) = R1eq(-s)
   Func(3) = 1.4010444_Rkind * QML_dnSigmoid_H3(s,betaQ)          +             &
         QML_dnSigmoid_H3(-s,betaQ)*(-s+1.6803643117748104_Rkind) +             &
            (0.1357176824058676_Rkind)     * tab_Pl(0)            +             &
            (0.0009647928495994389_Rkind)  * tab_Pl(1)            +             &
            (-0.1496072194236861_Rkind)    * tab_Pl(2)            +             &
            (0.004622944008685905_Rkind)   * tab_Pl(3)            +             &
            (0.013899404906043201_Rkind)   * tab_Pl(4)            +             &
            (-0.00406379918470803_Rkind)   * tab_Pl(5)            +             &
            (-0.0010528232085257284_Rkind) * tab_Pl(6)            +             &
            (-0.0006881469085679271_Rkind) * tab_Pl(7)            +             &
            (0.0009920029080101385_Rkind)  * tab_Pl(8)            +             &
            (-0.00014255050945681492_Rkind)* tab_Pl(9)            +             &
            (2.4954747402158953e-05_Rkind) * tab_Pl(10)           +             &
            (-0.0004543931972375909_Rkind) * tab_Pl(11)

    ! R3eq(s) = R1eq(s)+R2eq(s) (!linear HHH)
    Func(4) = Func(2) + Func(3)


    ! RPH parameters ap=(a(1)+a(2))/2
    ap  =                  &
    (0.10654140540118079_Rkind) * tab_Pl(0) +  &
    (-0.08156956206363981_Rkind) * tab_Pl(2) +  &
    (0.05165961900675448_Rkind) * tab_Pl(4) +  &
    (-0.02712759178768974_Rkind) * tab_Pl(6) +  &
    (0.01297477538446188_Rkind) * tab_Pl(8) +  &
    (-0.006231119449146584_Rkind) * tab_Pl(10) +  &
    (0.0029194223124278493_Rkind) * tab_Pl(12) +  &
    (-0.0015729890889698303_Rkind) * tab_Pl(14) +  &
    (0.0007856122821166629_Rkind) * tab_Pl(16) +  &
    (-0.0004478849431914204_Rkind) * tab_Pl(18) +  &
    (0.0002493501993969134_Rkind) * tab_Pl(20)

    ! RPH paramter: ap=(a(1)-a(2))/2
    am  =                  &
    (-0.24947389963365996_Rkind) * tab_Pl(1) +  &
    (0.11071787203324072_Rkind) * tab_Pl(3) +  &
    (-0.051987344951231605_Rkind) * tab_Pl(5) +  &
    (0.02429237299028829_Rkind) * tab_Pl(7) +  &
    (-0.01126303087651565_Rkind) * tab_Pl(9) +  &
    (0.005631218234846266_Rkind) * tab_Pl(11) +  &
    (-0.003014340083025705_Rkind) * tab_Pl(13) +  &
    (0.0014056596059369348_Rkind) * tab_Pl(15) +  &
    (-0.000660824039222194_Rkind) * tab_Pl(17) +  &
    (-7.516971220985713e-05_Rkind) * tab_Pl(19)

    Func(5) = ap+am
    Func(6) = ap-am

    DO i=0,max_deg
      CALL dealloc_dnS(tab_Pl(i))
    END DO
    CALL dealloc_dnS(s)
    CALL dealloc_dnS(ts)

  END SUBROUTINE EvalFunc_QML_H3_v11


  SUBROUTINE EvalFunc_QML_H3_v1(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H3_t),      intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    TYPE (dnS_t)                  :: s,ts,dnR,dnTh
    integer                       :: i
    integer,           parameter  :: max_deg = 29
    TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
    real (kind=Rkind), parameter  :: betaQ = 1.4_Rkind


    s  = dnQ(1)
    ts = tanh(s)

    DO i=0,max_deg
      tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
    END DO

    ! potential
    Func(1) = -0.17447440045043028_Rkind +                                      &
              (0.010858975893853413_Rkind)    * tab_Pl(0)  +                    &
              (-0.009667492179225417_Rkind)   * tab_Pl(2)  +                    &
              (-0.0006089583649096594_Rkind)  * tab_Pl(4)  +                    &
              (-0.0005463286353325458_Rkind)  * tab_Pl(6)  +                    &
              (-7.674511378387464e-05_Rkind)  * tab_Pl(8)  +                    &
              (-1.2703773915447383e-06_Rkind) * tab_Pl(10) +                    &
              (1.7588216355925958e-05_Rkind)  * tab_Pl(12) +                    &
              (1.7051544329081514e-05_Rkind)  * tab_Pl(14) +                    &
              (-3.7481719120090565e-05_Rkind) * tab_Pl(16)

   ! R1eq
   Func(2) = 1.4010444_Rkind * QML_dnSigmoid_H3(-s,betaQ)         +             &
           QML_dnSigmoid_H3(s,betaQ)*(s+1.6803643117748104_Rkind) +             &
            (0.1357176824058676_Rkind)     * tab_Pl(0)            +             &
            (-0.0009647928495994389_Rkind) * tab_Pl(1)            +             &
            (-0.1496072194236861_Rkind)    * tab_Pl(2)            +             &
            (-0.004622944008685905_Rkind)  * tab_Pl(3)            +             &
            (0.013899404906043201_Rkind)   * tab_Pl(4)            +             &
            (0.00406379918470803_Rkind)    * tab_Pl(5)            +             &
            (-0.0010528232085257284_Rkind) * tab_Pl(6)            +             &
            (0.0006881469085679271_Rkind)  * tab_Pl(7)            +             &
            (0.0009920029080101385_Rkind)  * tab_Pl(8)            +             &
            (0.00014255050945681492_Rkind) * tab_Pl(9)            +             &
            (2.4954747402158953e-05_Rkind) * tab_Pl(10)           +             &
            (0.0004543931972375909_Rkind)  * tab_Pl(11)

   ! R2eq(s) = R1eq(-s)
   Func(3) = 1.4010444_Rkind * QML_dnSigmoid_H3(s,betaQ)          +             &
         QML_dnSigmoid_H3(-s,betaQ)*(-s+1.6803643117748104_Rkind) +             &
            (0.1357176824058676_Rkind)     * tab_Pl(0)            +             &
            (0.0009647928495994389_Rkind)  * tab_Pl(1)            +             &
            (-0.1496072194236861_Rkind)    * tab_Pl(2)            +             &
            (0.004622944008685905_Rkind)   * tab_Pl(3)            +             &
            (0.013899404906043201_Rkind)   * tab_Pl(4)            +             &
            (-0.00406379918470803_Rkind)   * tab_Pl(5)            +             &
            (-0.0010528232085257284_Rkind) * tab_Pl(6)            +             &
            (-0.0006881469085679271_Rkind) * tab_Pl(7)            +             &
            (0.0009920029080101385_Rkind)  * tab_Pl(8)            +             &
            (-0.00014255050945681492_Rkind)* tab_Pl(9)            +             &
            (2.4954747402158953e-05_Rkind) * tab_Pl(10)           +             &
            (-0.0004543931972375909_Rkind) * tab_Pl(11)

    ! R3eq(s) = R1eq(s)+R2eq(s) (!linear HHH)
    Func(4) = Func(2) + Func(3)


    ! RPH parameters r
    dnR      =                                                                  &
              (0.26766607524797303_Rkind)     * tab_Pl(0)  +                    &
              (6.806052032746861e-09_Rkind)   * tab_Pl(1)  +                    &
              (-0.0015535754772908871_Rkind)  * tab_Pl(2)  +                    &
              (2.0009552260491068e-08_Rkind)  * tab_Pl(3)  +                    &
              (-0.012582288956034828_Rkind)   * tab_Pl(4)  +                    &
              (-2.8073095955518356e-08_Rkind) * tab_Pl(5)  +                    &
              (0.010502145645054565_Rkind)    * tab_Pl(6)  +                    &
              (2.5402220408683394e-08_Rkind)  * tab_Pl(7)  +                    &
              (-0.00678496091874973_Rkind)    * tab_Pl(8)  +                    &
              (-1.847484282317999e-08_Rkind)  * tab_Pl(9)  +                    &
              (0.003474989262273117_Rkind)    * tab_Pl(10) +                    &
              (1.2017866799932617e-08_Rkind)  * tab_Pl(11) +                    &
              (-0.0019395131094161956_Rkind)  * tab_Pl(12) +                    &
              (-7.233204033268473e-09_Rkind)  * tab_Pl(13) +                    &
              (0.0013262986972957342_Rkind)   * tab_Pl(14) +                    &
              (4.2016095723736e-09_Rkind)     * tab_Pl(15) +                    &
              (-0.0007299784927963528_Rkind)  * tab_Pl(16) +                    &
              (-2.401414821263817e-09_Rkind)  * tab_Pl(17) +                    &
              (0.00044078321280719524_Rkind)  * tab_Pl(18) +                    &
              (1.3108656875294227e-09_Rkind)  * tab_Pl(19) +                    &
              (0.00022352966903015168_Rkind)  * tab_Pl(20) +                    &
              (-7.518732295064561e-10_Rkind)  * tab_Pl(21) +                    &
              (-0.00036549376491577264_Rkind) * tab_Pl(22) +                    &
              (4.1612329783873094e-10_Rkind)  * tab_Pl(23) +                    &
              (0.0004459035924032867_Rkind)   * tab_Pl(24) +                    &
              (-2.6223332018163205e-10_Rkind) * tab_Pl(25) +                    &
              (2.273672442589093e-05_Rkind)   * tab_Pl(26) +                    &
              (1.447938227854266e-10_Rkind)   * tab_Pl(27) +                    &
              (-1.855501791553842e-05_Rkind)  * tab_Pl(28) +                    &
              (-9.176168686383644e-11_Rkind)  * tab_Pl(29)

    ! RPH paramter: th
    dnTh     =                                                                  &
              (0.7853980834870691_Rkind)      * tab_Pl(0)  +                    &
              (1.6556773456797165_Rkind)      * tab_Pl(1)  +                    &
              (3.4626851455522983e-07_Rkind)  * tab_Pl(2)  +                    &
              (-0.5731854583939033_Rkind)     * tab_Pl(3)  +                    &
              (-2.8771835705821473e-07_Rkind) * tab_Pl(4)  +                    &
              (0.2350158022234596_Rkind)      * tab_Pl(5)  +                    &
              (1.989230244253878e-07_Rkind)   * tab_Pl(6)  +                    &
              (-0.10219520909830185_Rkind)    * tab_Pl(7)  +                    &
              (-1.2621208480613884e-07_Rkind) * tab_Pl(8)  +                    &
              (0.048468715568805845_Rkind)    * tab_Pl(9)  +                    &
              (7.696291823409225e-08_Rkind)   * tab_Pl(10) +                    &
              (-0.024524927534999522_Rkind)   * tab_Pl(11) +                    &
              (-4.58667842666333e-08_Rkind)   * tab_Pl(12) +                    &
              (0.014023814755791649_Rkind)    * tab_Pl(13) +                    &
              (2.773916912783092e-08_Rkind)   * tab_Pl(14) +                    &
              (-0.006958147875138725_Rkind)   * tab_Pl(15) +                    &
              (-1.639868076367661e-08_Rkind)  * tab_Pl(16) +                    &
              (0.003770984041691324_Rkind)    * tab_Pl(17) +                    &
              (9.721641212735887e-09_Rkind)   * tab_Pl(18) +                    &
              (-0.0019178090260791802_Rkind)  * tab_Pl(19) +                    &
              (-5.691839285240476e-09_Rkind)  * tab_Pl(20) +                    &
              (0.001074145733356513_Rkind)    * tab_Pl(21) +                    &
              (3.3623843187612974e-09_Rkind)  * tab_Pl(22) +                    &
              (-0.0005259206048766796_Rkind)  * tab_Pl(23) +                    &
              (-1.9704799248725513e-09_Rkind) * tab_Pl(24) +                    &
              (0.0003092457483074583_Rkind)   * tab_Pl(25) +                    &
              (1.2005965897009222e-09_Rkind)  * tab_Pl(26) +                    &
              (-0.00013246994648197604_Rkind) * tab_Pl(27) +                    &
              (-4.627433514270627e-10_Rkind)  * tab_Pl(28) +                    &
              (9.114113633803292e-05_Rkind)   * tab_Pl(29)

    Func(5) = dnR * cos(dnTh)
    Func(6) = dnR * sin(dnTh)

    DO i=0,max_deg
      CALL dealloc_dnS(tab_Pl(i))
    END DO
    CALL dealloc_dnS(s)
    CALL dealloc_dnS(ts)

  END SUBROUTINE EvalFunc_QML_H3_v1

 ! for IRC, RPH (IRC with the correct data ?)
  SUBROUTINE EvalFunc_QML_H3_v21(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE
  
      CLASS(QML_H3_t),      intent(in)    :: QModel
      TYPE (dnS_t),         intent(inout) :: Func(:)
      TYPE (dnS_t),         intent(in)    :: dnQ(:)
      integer,              intent(in)    :: nderiv
  
      TYPE (dnS_t)                  :: s,ts,am,ap
      integer                       :: i
      integer,           parameter  :: max_deg = 20
      TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
      real (kind=Rkind), parameter  :: betaQ = 1.4_Rkind
  
  
      s  = dnQ(1)
      ts = tanh(s)
  
      DO i=0,max_deg
        tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
      END DO
  
      ! potential
      Func(1)  = -0.17447440045043028_Rkind +                  &
      (0.008932139934203165_Rkind) * tab_Pl(0) +  &
      (-0.010702285888249915_Rkind) * tab_Pl(2) +  &
      (0.0020654276667500387_Rkind) * tab_Pl(4) +  &
      (-0.0008489800485404381_Rkind) * tab_Pl(6) +  &
      (0.0006273408592540359_Rkind) * tab_Pl(8) +  &
      (-0.00022735810861969512_Rkind) * tab_Pl(10) +  &
      (0.0001325111874514195_Rkind) * tab_Pl(12) +  &
      (-0.00010543827157114539_Rkind) * tab_Pl(14) +  &
      (2.1653267283543284e-05_Rkind) * tab_Pl(16) +  &
      (-3.974237016186115e-05_Rkind) * tab_Pl(18) +  &
      (-3.794957260182394e-06_Rkind) * tab_Pl(20)
  
     ! R1eq
      Func(2) = 1.4010444_Rkind * QML_dnSigmoid_H3(-s,ONE) +                      &
            QML_dnSigmoid_H3(s,ONE) *                                             &
              (1.2249168019329821_Rkind * s + 1.8405491180735511_Rkind) +         &
            (0.14390507395498883_Rkind) * tab_Pl(0) +  &
            (0.07653458032053517_Rkind) * tab_Pl(1) +  &
            (-0.05206829932317228_Rkind) * tab_Pl(2) +  &
            (-0.10987189662022206_Rkind) * tab_Pl(3) +  &
            (-0.08378050418420241_Rkind) * tab_Pl(4) +  &
            (0.04798298767730978_Rkind) * tab_Pl(5) +  &
            (0.006719137603736121_Rkind) * tab_Pl(6) +  &
            (-0.021578130275728205_Rkind) * tab_Pl(7) +  &
            (-0.007702494895032822_Rkind) * tab_Pl(8) +  &
            (0.009200844085327998_Rkind) * tab_Pl(9) +  &
            (-0.003084477761700385_Rkind) * tab_Pl(10) +  &
            (-0.0025349840362887985_Rkind) * tab_Pl(11) +  &
            (0.0008514056605674219_Rkind) * tab_Pl(12) +  &
            (-0.00018447867084684657_Rkind) * tab_Pl(13) +  &
            (-0.0026926646885450938_Rkind) * tab_Pl(14) +  &
            (0.0012472189535898588_Rkind) * tab_Pl(15) +  &
            (0.0009512718361847219_Rkind) * tab_Pl(16) +  &
            (-0.0012918715574876767_Rkind) * tab_Pl(17) +  &
            (-0.0011677786096285978_Rkind) * tab_Pl(18) +  &
            (0.0007329160177868408_Rkind) * tab_Pl(19)
  
     ! R2eq(s) = R1eq(-s)
      Func(3) = 1.4010444_Rkind * QML_dnSigmoid_H3(s,ONE) +                       &
            QML_dnSigmoid_H3(-s,ONE) *                                            &
              (-1.2249168019329821_Rkind * s + 1.8405491180735511_Rkind) +        &
            (0.14390507395498883_Rkind) * tab_Pl(0) +  &
            (-0.07653458032053517_Rkind) * tab_Pl(1) +  &
            (-0.05206829932317228_Rkind) * tab_Pl(2) +  &
            (0.10987189662022206_Rkind) * tab_Pl(3) +  &
            (-0.08378050418420241_Rkind) * tab_Pl(4) +  &
            (-0.04798298767730978_Rkind) * tab_Pl(5) +  &
            (0.006719137603736121_Rkind) * tab_Pl(6) +  &
            (0.021578130275728205_Rkind) * tab_Pl(7) +  &
            (-0.007702494895032822_Rkind) * tab_Pl(8) +  &
            (-0.009200844085327998_Rkind) * tab_Pl(9) +  &
            (-0.003084477761700385_Rkind) * tab_Pl(10) +  &
            (0.0025349840362887985_Rkind) * tab_Pl(11) +  &
            (0.0008514056605674219_Rkind) * tab_Pl(12) +  &
            (0.00018447867084684657_Rkind) * tab_Pl(13) +  &
            (-0.0026926646885450938_Rkind) * tab_Pl(14) +  &
            (-0.0012472189535898588_Rkind) * tab_Pl(15) +  &
            (0.0009512718361847219_Rkind) * tab_Pl(16) +  &
            (0.0012918715574876767_Rkind) * tab_Pl(17) +  &
            (-0.0011677786096285978_Rkind) * tab_Pl(18) +  &
            (-0.0007329160177868408_Rkind) * tab_Pl(19)  
      ! R3eq(s) = R1eq(s)+R2eq(s) (!linear HHH)
      Func(4) = Func(2) + Func(3)
  
  
      ! RPH parameters ap=(a(1)+a(2))/2 wrong
      ap  =                  &
      (0.10654140540118079_Rkind) * tab_Pl(0) +  &
      (-0.08156956206363981_Rkind) * tab_Pl(2) +  &
      (0.05165961900675448_Rkind) * tab_Pl(4) +  &
      (-0.02712759178768974_Rkind) * tab_Pl(6) +  &
      (0.01297477538446188_Rkind) * tab_Pl(8) +  &
      (-0.006231119449146584_Rkind) * tab_Pl(10) +  &
      (0.0029194223124278493_Rkind) * tab_Pl(12) +  &
      (-0.0015729890889698303_Rkind) * tab_Pl(14) +  &
      (0.0007856122821166629_Rkind) * tab_Pl(16) +  &
      (-0.0004478849431914204_Rkind) * tab_Pl(18) +  &
      (0.0002493501993969134_Rkind) * tab_Pl(20)
  
      ! RPH paramter: ap=(a(1)-a(2))/2 wrong
      am  =                  &
      (-0.24947389963365996_Rkind) * tab_Pl(1) +  &
      (0.11071787203324072_Rkind) * tab_Pl(3) +  &
      (-0.051987344951231605_Rkind) * tab_Pl(5) +  &
      (0.02429237299028829_Rkind) * tab_Pl(7) +  &
      (-0.01126303087651565_Rkind) * tab_Pl(9) +  &
      (0.005631218234846266_Rkind) * tab_Pl(11) +  &
      (-0.003014340083025705_Rkind) * tab_Pl(13) +  &
      (0.0014056596059369348_Rkind) * tab_Pl(15) +  &
      (-0.000660824039222194_Rkind) * tab_Pl(17) +  &
      (-7.516971220985713e-05_Rkind) * tab_Pl(19)
  
      Func(5) = ap+am
      Func(6) = ap-am
  
      DO i=0,max_deg
        CALL dealloc_dnS(tab_Pl(i))
      END DO
      CALL dealloc_dnS(s)
      CALL dealloc_dnS(ts)
  
    END SUBROUTINE EvalFunc_QML_H3_v21

     SUBROUTINE EvalFunc_QML_H3_v31(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

      CLASS(QML_H3_t),      intent(in)    :: QModel
      TYPE (dnS_t),         intent(inout) :: Func(:)
      TYPE (dnS_t),         intent(in)    :: dnQ(:)
      integer,              intent(in)    :: nderiv

      TYPE (dnS_t)                  :: s,ts,am,ap
      integer                       :: i
      integer,           parameter  :: max_deg = 20
      TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
      real (kind=Rkind), parameter  :: betaQ = 1._Rkind


      s  = dnQ(1)
      ts = tanh(s*betaQ)

      DO i=0,max_deg
        tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
      END DO
     
      !  V
      Func(1)= &
      (-0.16597446466905263_Rkind) * tab_Pl(0) +  &
      (-0.012383613848784496_Rkind) * tab_Pl(2) +  &
      (0.003326606346209517_Rkind) * tab_Pl(4) +  &
      (0.0008779335935679601_Rkind) * tab_Pl(6) +  &
      (-0.00024065622197537675_Rkind) * tab_Pl(8) +  &
      (-0.0001241975119023965_Rkind) * tab_Pl(10) +  &
      (-1.3361183199880677e-05_Rkind) * tab_Pl(12) +  &
      (8.34191640793267e-05_Rkind) * tab_Pl(14) +  &
      (-3.0024892518874162e-06_Rkind) * tab_Pl(16)
      
      ! hessian
      Func(2) = &
      (0.5160182951754864_Rkind) * tab_Pl(0) +  &
      (-0.17971203676655206_Rkind) * tab_Pl(2) +  &
      (0.065778695179622_Rkind) * tab_Pl(4) +  &
      (-0.056073060962388566_Rkind) * tab_Pl(6) +  &
      (0.004030092812838212_Rkind) * tab_Pl(8) +  &
      (0.021678208287918613_Rkind) * tab_Pl(10) +  &
      (0.00825140847920465_Rkind) * tab_Pl(12) +  &
      (-0.007827161286743232_Rkind) * tab_Pl(14) +  &
      (-0.005497405564310822_Rkind) * tab_Pl(16)        

      ! rho
      Func(3) = &
      (1.2933688878722176_Rkind) * tab_Pl(0) +  &
      (0.08694736684022922_Rkind) * tab_Pl(2) +  &
      (-0.016117447456099533_Rkind) * tab_Pl(4) +  &
      (0.011110605473470134_Rkind) * tab_Pl(6) +  &
      (0.009578713881238925_Rkind) * tab_Pl(8) +  &
      (0.0015319141117551551_Rkind) * tab_Pl(10) +  &
      (-0.001136293339165834_Rkind) * tab_Pl(12) +  &
      (0.0007092933469349128_Rkind) * tab_Pl(14) +  &
      (0.007085765600105608_Rkind) * tab_Pl(16)
                                     
      write(out_unit,*) " hi universe we're in the EvalFunc 31 !!"

      DO i=0,max_deg
        CALL dealloc_dnS(tab_Pl(i))
      END DO
      CALL dealloc_dnS(s)
      CALL dealloc_dnS(ts)

    END SUBROUTINE EvalFunc_QML_H3_v31

     SUBROUTINE EvalFunc_QML_H3_v41(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

      CLASS(QML_H3_t),      intent(in)    :: QModel
      TYPE (dnS_t),         intent(inout) :: Func(:)
      TYPE (dnS_t),         intent(in)    :: dnQ(:)
      integer,              intent(in)    :: nderiv

      TYPE (dnS_t)                  :: s,ts,am,ap
      integer                       :: i
      integer,           parameter  :: max_deg = 20
      TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
      real (kind=Rkind), parameter  :: betaQ = 1._Rkind


      !s  = dnQ(1)
      s  = sinh(dnQ(1))
      ts = tanh(s*betaQ)
      DO i=0,max_deg
        tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
      END DO
     
      !  V
      Func(1)= &
      (-0.16597446466905263_Rkind) * tab_Pl(0) +  &
      (-0.012383613848784496_Rkind) * tab_Pl(2) +  &
      (0.003326606346209517_Rkind) * tab_Pl(4) +  &
      (0.0008779335935679601_Rkind) * tab_Pl(6) +  &
      (-0.00024065622197537675_Rkind) * tab_Pl(8) +  &
      (-0.0001241975119023965_Rkind) * tab_Pl(10) +  &
      (-1.3361183199880677e-05_Rkind) * tab_Pl(12) +  &
      (8.34191640793267e-05_Rkind) * tab_Pl(14) +  &
      (-3.0024892518874162e-06_Rkind) * tab_Pl(16)
      
      ! hessian
      Func(2) = &
      (0.5160182951754864_Rkind) * tab_Pl(0) +  &
      (-0.17971203676655206_Rkind) * tab_Pl(2) +  &
      (0.065778695179622_Rkind) * tab_Pl(4) +  &
      (-0.056073060962388566_Rkind) * tab_Pl(6) +  &
      (0.004030092812838212_Rkind) * tab_Pl(8) +  &
      (0.021678208287918613_Rkind) * tab_Pl(10) +  &
      (0.00825140847920465_Rkind) * tab_Pl(12) +  &
      (-0.007827161286743232_Rkind) * tab_Pl(14) +  &
      (-0.005497405564310822_Rkind) * tab_Pl(16)        

      ! rho
      Func(3) = &
      (1.2933688878722176_Rkind) * tab_Pl(0) +  &
      (0.08694736684022922_Rkind) * tab_Pl(2) +  &
      (-0.016117447456099533_Rkind) * tab_Pl(4) +  &
      (0.011110605473470134_Rkind) * tab_Pl(6) +  &
      (0.009578713881238925_Rkind) * tab_Pl(8) +  &
      (0.0015319141117551551_Rkind) * tab_Pl(10) +  &
      (-0.001136293339165834_Rkind) * tab_Pl(12) +  &
      (0.0007092933469349128_Rkind) * tab_Pl(14) +  &
      (0.007085765600105608_Rkind) * tab_Pl(16)
                                     
!      write(out_unit,*) 'new s,s,V',get_d0(dnQ(1)),get_d0(s),get_d0(Func(1))

      DO i=0,max_deg
        CALL dealloc_dnS(tab_Pl(i))
      END DO
      CALL dealloc_dnS(s)
      CALL dealloc_dnS(ts)
    END SUBROUTINE EvalFunc_QML_H3_v41


     SUBROUTINE EvalFunc_QML_H3_v51(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

      CLASS(QML_H3_t),      intent(in)    :: QModel
      TYPE (dnS_t),         intent(inout) :: Func(:)
      TYPE (dnS_t),         intent(in)    :: dnQ(:)
      integer,              intent(in)    :: nderiv

      TYPE (dnS_t)                  :: s,ts,am,ap
      integer                       :: i
      integer,           parameter  :: max_deg = 32
      TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
      real (kind=Rkind), parameter  :: betaQ = 1._Rkind


      s  = dnQ(1)
      ts  = sin(atan(s))
      !ts = tanh(s*betaQ)
      DO i=0,max_deg
        tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
      END DO
     
      !  V
      Func(1)=  &
      (-0.1620718837989824_Rkind) * tab_Pl(0) +  & 
      (-0.008993773246295272_Rkind) * tab_Pl(2) +  &
      (-0.004330516349166968_Rkind) * tab_Pl(4) +  &
      (-0.0010620585935713837_Rkind) * tab_Pl(6) +  &
      (0.0003858633034138086_Rkind) * tab_Pl(8) +  &
      (0.0007669172511520921_Rkind) * tab_Pl(10) +  &
      (0.0006466391218323583_Rkind) * tab_Pl(12) +  &
      (0.0003753917050891673_Rkind) * tab_Pl(14) +  &
      (0.0001217153671666159_Rkind) * tab_Pl(16) +  &
      (-3.806281720474181e-05_Rkind) * tab_Pl(18) +  &
      (-0.00010291012057228436_Rkind) * tab_Pl(20) +  &
      (-0.00010013347559712096_Rkind) * tab_Pl(22) +  &
      (-6.607022342401301e-05_Rkind) * tab_Pl(24) +  &
      (-2.5997737071039032e-05_Rkind) * tab_Pl(26) +  &
      (1.564175665510225e-06_Rkind) * tab_Pl(28) +  &
      (1.5410394805893913e-05_Rkind) * tab_Pl(30) +  &
      (3.4267515659467252e-06_Rkind) * tab_Pl(32)

      
      ! hessian
      Func(2) =  &
      (0.6075899355887548_Rkind) * tab_Pl(0) +  & 
      (-0.10567883781305588_Rkind) * tab_Pl(2) +  &
      (-0.05714381041852741_Rkind) * tab_Pl(4) +  &
      (-0.02935995581332255_Rkind) * tab_Pl(6) +  &
      (-0.016991066224729384_Rkind) * tab_Pl(8) +  &
      (-0.014537719234377633_Rkind) * tab_Pl(10) +  &
      (-0.013177254444038779_Rkind) * tab_Pl(12) +  &
      (-0.0060939317509078016_Rkind) * tab_Pl(14) +  &
      (-0.0003379189027574508_Rkind) * tab_Pl(16) +  &
      (-0.003984506386244121_Rkind) * tab_Pl(18) +  &
      (0.0025798395807923933_Rkind) * tab_Pl(20) +  &
      (0.0015280532748144886_Rkind) * tab_Pl(22) +  &
      (0.0025480533183245015_Rkind) * tab_Pl(24) +  &
      (-0.0003994112348223169_Rkind) * tab_Pl(26) +  &
      (0.0010596234077735716_Rkind) * tab_Pl(28) +  &
      (-0.0007077615662376131_Rkind) * tab_Pl(30) +  &
      (0.0025738127648113544_Rkind) * tab_Pl(32)

 

      ! rho
      Func(3) = &
      (1.261243931645815_Rkind) * tab_Pl(0) +  &
      (0.05498541052261184_Rkind) * tab_Pl(2) +  &
      (0.03446819696272064_Rkind) * tab_Pl(4) +  &
      (0.020291272930553106_Rkind) * tab_Pl(6) +  &
      (0.011855010123711414_Rkind) * tab_Pl(8) +  &
      (0.007536688407402705_Rkind) * tab_Pl(10) +  &
      (0.004940560392031361_Rkind) * tab_Pl(12) +  &
      (0.0034254508799222323_Rkind) * tab_Pl(14) +  &
      (0.002319413680131285_Rkind) * tab_Pl(16) +  &
      (0.0010470828509148065_Rkind) * tab_Pl(18) +  &
      (4.281256394175054e-05_Rkind) * tab_Pl(20) +  &
      (-0.00030975651066805104_Rkind) * tab_Pl(22) +  &
      (-0.0005399637278206589_Rkind) * tab_Pl(24) +  &
      (-0.00041279694640769006_Rkind) * tab_Pl(26) +  &
      (-0.00017611772257669623_Rkind) * tab_Pl(28) +  &
      (-2.5934623118451503e-05_Rkind) * tab_Pl(30) +  &
      (0.00036253207816267066_Rkind) * tab_Pl(32)

      DO i=0,max_deg
        CALL dealloc_dnS(tab_Pl(i))
      END DO
      CALL dealloc_dnS(s)
      CALL dealloc_dnS(ts)
    END SUBROUTINE EvalFunc_QML_H3_v51


     SUBROUTINE EvalFunc_QML_H3_v61(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

      CLASS(QML_H3_t),      intent(in)    :: QModel
      TYPE (dnS_t),         intent(inout) :: Func(:)
      TYPE (dnS_t),         intent(in)    :: dnQ(:)
      integer,              intent(in)    :: nderiv

      TYPE (dnS_t)                  :: s,ts,am,ap
      integer                       :: i
      integer,           parameter  :: max_deg = 32
      TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
      real (kind=Rkind), parameter  :: betaQ = 1._Rkind


      s  = dnQ(1)
      ts  = sin(atan(s))
      !ts = tanh(s*betaQ)
      DO i=0,max_deg
        tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
      END DO
     
      !  V
      Func(1)=  &
      (-0.1620718837989824_Rkind) * tab_Pl(0) +  & 
      (-0.008993773246295272_Rkind) * tab_Pl(2) +  &
      (-0.004330516349166968_Rkind) * tab_Pl(4) +  &
      (-0.0010620585935713837_Rkind) * tab_Pl(6) +  &
      (0.0003858633034138086_Rkind) * tab_Pl(8) +  &
      (0.0007669172511520921_Rkind) * tab_Pl(10) +  &
      (0.0006466391218323583_Rkind) * tab_Pl(12) +  &
      (0.0003753917050891673_Rkind) * tab_Pl(14) +  &
      (0.0001217153671666159_Rkind) * tab_Pl(16) +  &
      (-3.806281720474181e-05_Rkind) * tab_Pl(18) +  &
      (-0.00010291012057228436_Rkind) * tab_Pl(20) +  &
      (-0.00010013347559712096_Rkind) * tab_Pl(22) +  &
      (-6.607022342401301e-05_Rkind) * tab_Pl(24) +  &
      (-2.5997737071039032e-05_Rkind) * tab_Pl(26) +  &
      (1.564175665510225e-06_Rkind) * tab_Pl(28) +  &
      (1.5410394805893913e-05_Rkind) * tab_Pl(30) +  &
      (3.4267515659467252e-06_Rkind) * tab_Pl(32)

      
      ! hessian
      Func(2) =  &
      (0.6075899355887548_Rkind) * tab_Pl(0) +  & 
      (-0.10567883781305588_Rkind) * tab_Pl(2) +  &
      (-0.05714381041852741_Rkind) * tab_Pl(4) +  &
      (-0.02935995581332255_Rkind) * tab_Pl(6) +  &
      (-0.016991066224729384_Rkind) * tab_Pl(8) +  &
      (-0.014537719234377633_Rkind) * tab_Pl(10) +  &
      (-0.013177254444038779_Rkind) * tab_Pl(12) +  &
      (-0.0060939317509078016_Rkind) * tab_Pl(14) +  &
      (-0.0003379189027574508_Rkind) * tab_Pl(16) +  &
      (-0.003984506386244121_Rkind) * tab_Pl(18) +  &
      (0.0025798395807923933_Rkind) * tab_Pl(20) +  &
      (0.0015280532748144886_Rkind) * tab_Pl(22) +  &
      (0.0025480533183245015_Rkind) * tab_Pl(24) +  &
      (-0.0003994112348223169_Rkind) * tab_Pl(26) +  &
      (0.0010596234077735716_Rkind) * tab_Pl(28) +  &
      (-0.0007077615662376131_Rkind) * tab_Pl(30) +  &
      (0.0025738127648113544_Rkind) * tab_Pl(32)

 

      ! rho
      Func(3) = &
      (1.261243931645815_Rkind) * tab_Pl(0) +  &
      (0.05498541052261184_Rkind) * tab_Pl(2) +  &
      (0.03446819696272064_Rkind) * tab_Pl(4) +  &
      (0.020291272930553106_Rkind) * tab_Pl(6) +  &
      (0.011855010123711414_Rkind) * tab_Pl(8) +  &
      (0.007536688407402705_Rkind) * tab_Pl(10) +  &
      (0.004940560392031361_Rkind) * tab_Pl(12) +  &
      (0.0034254508799222323_Rkind) * tab_Pl(14) +  &
      (0.002319413680131285_Rkind) * tab_Pl(16) +  &
      (0.0010470828509148065_Rkind) * tab_Pl(18) +  &
      (4.281256394175054e-05_Rkind) * tab_Pl(20) +  &
      (-0.00030975651066805104_Rkind) * tab_Pl(22) +  &
      (-0.0005399637278206589_Rkind) * tab_Pl(24) +  &
      (-0.00041279694640769006_Rkind) * tab_Pl(26) +  &
      (-0.00017611772257669623_Rkind) * tab_Pl(28) +  &
      (-2.5934623118451503e-05_Rkind) * tab_Pl(30) +  &
      (0.00036253207816267066_Rkind) * tab_Pl(32)

      !Fitted Frequency along Rho
      Func(4) = &
      (0.006005922154491725_Rkind) * tab_Pl(0) +  &
      (0.0031574819729319747_Rkind) * tab_Pl(2) +  &
      (0.0008492844551779132_Rkind) * tab_Pl(4) +  &
      (0.0002242725684159113_Rkind) * tab_Pl(6) +  &
      (1.4799836843707549e-05_Rkind) * tab_Pl(8) +  &
      (-7.148382422659551e-05_Rkind) * tab_Pl(10) +  &
      (-9.724998095625515e-05_Rkind) * tab_Pl(12) +  &
      (-6.622114649532386e-05_Rkind) * tab_Pl(14) +  &
      (-3.2641975806064845e-05_Rkind) * tab_Pl(16) +  &
      (-3.581797770238322e-05_Rkind) * tab_Pl(18) +  &
      (-2.6512344100050384e-06_Rkind) * tab_Pl(20) +  &
      (4.311136006173978e-06_Rkind) * tab_Pl(22) +  &
      (4.528040827611103e-06_Rkind) * tab_Pl(24) +  &
      (4.431197572981137e-07_Rkind) * tab_Pl(26) +  &
      (-8.354627843217928e-06_Rkind) * tab_Pl(28) +  &
      (3.2796344679545e-06_Rkind) * tab_Pl(30) +  &
      (7.753347481333872e-05_Rkind) * tab_Pl(32)

                                     

      DO i=0,max_deg
        CALL dealloc_dnS(tab_Pl(i))
      END DO
      CALL dealloc_dnS(s)
      CALL dealloc_dnS(ts)
    END SUBROUTINE EvalFunc_QML_H3_v61

     SUBROUTINE EvalFunc_QML_H3_v72(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

      CLASS(QML_H3_t),      intent(in)    :: QModel
      TYPE (dnS_t),         intent(inout) :: Func(:)
      TYPE (dnS_t),         intent(in)    :: dnQ(:)
      integer,              intent(in)    :: nderiv

      TYPE (dnS_t)                  :: s,ts,am,ap
      integer                       :: i
      integer,           parameter  :: max_deg = 32
      TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
      real (kind=Rkind), parameter  :: betaQ = 1._Rkind


      s  = dnQ(1)
      ts  = sin(atan(s))
      !ts = tanh(s*betaQ)
      DO i=0,max_deg
        tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
      END DO

      ! Fit polynomial en c4 + c3*drho + c2*drho^2
      !  c4( =V0) fitté atour de rhoOPt à différentes val de s, puis fitté le long de s
      ! maxdiff =3.47995241222e-05
      Func(1)=  &
      (-0.1620716487102793_Rkind) * tab_Pl(0) +  &
      (-0.008992725629195533_Rkind) * tab_Pl(2) +  &
      (-0.00432915791758284_Rkind) * tab_Pl(4) +  &
      (-0.001061111040246674_Rkind) * tab_Pl(6) +  &
      (0.0003857228005224359_Rkind) * tab_Pl(8) +  &
      (0.0007652895174396595_Rkind) * tab_Pl(10) +  &
      (0.0006436145988118415_Rkind) * tab_Pl(12) +  &
      (0.00037157203206966737_Rkind) * tab_Pl(14) +  &
      (0.00011819866129507563_Rkind) * tab_Pl(16) +  &
      (-4.02663761147522e-05_Rkind) * tab_Pl(18) +  &
      (-0.00010232354419382081_Rkind) * tab_Pl(20) +  &
      (-9.760506805819436e-05_Rkind) * tab_Pl(22) +  &
      (-5.885890983333011e-05_Rkind) * tab_Pl(24) +  &
      (-2.1436826791342742e-05_Rkind) * tab_Pl(26) +  &
      (1.4784209082551226e-05_Rkind) * tab_Pl(28) +  &
      (2.4150163942799644e-05_Rkind) * tab_Pl(30) +  &
      (-3.501962307059586e-05_Rkind) * tab_Pl(32)


      
      ! c2 (c2*drho^2) fitté atour de rhoOPt à différentes val de s, puis fitté le long de s
      ! maxdiff = 2.75102563607e-05
      Func(2) =  &
      (0.30375204657676236_Rkind) * tab_Pl(0) +  &
      (-0.05261249923316819_Rkind) * tab_Pl(2) +  &
      (-0.028245088307616156_Rkind) * tab_Pl(4) +  &
      (-0.01376309967876779_Rkind) * tab_Pl(6) +  &
      (-0.008615637365472464_Rkind) * tab_Pl(8) +  &
      (-0.007122038317701137_Rkind) * tab_Pl(10) +  &
      (-0.005141765724950869_Rkind) * tab_Pl(12) +  &
      (-0.003418104118376876_Rkind) * tab_Pl(14) +  &
      (-0.0010447833949072116_Rkind) * tab_Pl(16) +  &
      (0.00040790797785053524_Rkind) * tab_Pl(18) +  &
      (0.0009037612950454548_Rkind) * tab_Pl(20) +  &
      (0.0005452967481115013_Rkind) * tab_Pl(22) +  &
      (0.000282087543543148_Rkind) * tab_Pl(24) +  &
      (-0.0002715479728813458_Rkind) * tab_Pl(26) +  &
      (-0.0007667274135405311_Rkind) * tab_Pl(28) +  &
      (-0.0007199692637246045_Rkind) * tab_Pl(30) +  &
      (0.0006760031365286_Rkind) * tab_Pl(32)

 

      ! rho
      Func(3) = &
      (1.261243931645815_Rkind) * tab_Pl(0) +  &
      (0.05498541052261184_Rkind) * tab_Pl(2) +  &
      (0.03446819696272064_Rkind) * tab_Pl(4) +  &
      (0.020291272930553106_Rkind) * tab_Pl(6) +  &
      (0.011855010123711414_Rkind) * tab_Pl(8) +  &
      (0.007536688407402705_Rkind) * tab_Pl(10) +  &
      (0.004940560392031361_Rkind) * tab_Pl(12) +  &
      (0.0034254508799222323_Rkind) * tab_Pl(14) +  &
      (0.002319413680131285_Rkind) * tab_Pl(16) +  &
      (0.0010470828509148065_Rkind) * tab_Pl(18) +  &
      (4.281256394175054e-05_Rkind) * tab_Pl(20) +  &
      (-0.00030975651066805104_Rkind) * tab_Pl(22) +  &
      (-0.0005399637278206589_Rkind) * tab_Pl(24) +  &
      (-0.00041279694640769006_Rkind) * tab_Pl(26) +  &
      (-0.00017611772257669623_Rkind) * tab_Pl(28) +  &
      (-2.5934623118451503e-05_Rkind) * tab_Pl(30) +  &
      (0.00036253207816267066_Rkind) * tab_Pl(32)
      
      ! c3(c3*drho) fitté atour de rhoOPt à différentes val de s, puis fitté le long de s 
      ! maxdiff = 5.57117016741e-06
      Func(4) = &   
      (0.00012054894006886368_Rkind) * tab_Pl(0) +  &
      (-4.159861752693137e-05_Rkind) * tab_Pl(2) +  &
      (-3.088431820462534e-06_Rkind) * tab_Pl(4) +  &
      (-1.5430069424391806e-05_Rkind) * tab_Pl(6) +  &
      (-8.389126757210443e-05_Rkind) * tab_Pl(8) +  &
      (9.425785955142528e-05_Rkind) * tab_Pl(10) +  &
      (-3.255469762512486e-05_Rkind) * tab_Pl(12) +  &
      (-3.114281736068148e-05_Rkind) * tab_Pl(14) +  &
      (6.793677577444186e-05_Rkind) * tab_Pl(16) +  &
      (-3.75084784501627e-06_Rkind) * tab_Pl(18) +  &
      (-5.64014489843293e-05_Rkind) * tab_Pl(20) +  &
      (4.147527388974017e-05_Rkind) * tab_Pl(22) +  &
      (-1.5965985196154143e-05_Rkind) * tab_Pl(24) +  &
      (-7.97292699073988e-06_Rkind) * tab_Pl(26) +  &
      (1.0429200371081538e-05_Rkind) * tab_Pl(28) +  &
      (-7.53503532862055e-06_Rkind) * tab_Pl(30) +  &
      (6.714016001365854e-06_Rkind) * tab_Pl(32)

      
      DO i=0,max_deg
        CALL dealloc_dnS(tab_Pl(i))
      END DO
      CALL dealloc_dnS(s)
      CALL dealloc_dnS(ts)
  END SUBROUTINE EvalFunc_QML_H3_v72

  SUBROUTINE EvalFunc_QML_H3_v74(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

      CLASS(QML_H3_t),      intent(in)    :: QModel
      TYPE (dnS_t),         intent(inout) :: Func(:)
      TYPE (dnS_t),         intent(in)    :: dnQ(:)
      integer,              intent(in)    :: nderiv

      TYPE (dnS_t)                  :: s,ts,am,ap
      integer                       :: i
      integer,           parameter  :: max_deg = 48
      TYPE (dnS_t)                  :: tab_Pl(0:max_deg)
      real (kind=Rkind), parameter  :: betaQ = 1._Rkind


      s  = dnQ(1)
      ts  = sin(atan(s))
      !ts = tanh(s*betaQ)
      DO i=0,max_deg
        tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
      END DO
     
      !  c4(V0) fitté atour de rhoOPt à différentes val de s, puis fitté le long de s
      ! maxdiff = 2.63656814312e-06
      Func(1)=  &!ces coeff correspondent aux pts enlevé entre 6.6 et 3.6 et 3.6 et 2.8
      (-0.16207139294104606_Rkind) * tab_Pl(0) +  &
      (-0.008989917909474423_Rkind) * tab_Pl(2) +  &
      (-0.004329368522365633_Rkind) * tab_Pl(4) +  &
      (-0.0010482901505167747_Rkind) * tab_Pl(6) +  &
      (0.000360202857892648_Rkind) * tab_Pl(8) +  &
      (0.0008368367197993513_Rkind) * tab_Pl(10) +  &
      (0.0004787805555940436_Rkind) * tab_Pl(12) +  &
      (0.0007176755497951032_Rkind) * tab_Pl(14) +  &
      (-0.0005593458069133195_Rkind) * tab_Pl(16) +  &
      (0.0011473839801117314_Rkind) * tab_Pl(18) +  &
      (-0.002011414626981134_Rkind) * tab_Pl(20) +  &
      (0.0026394004582187968_Rkind) * tab_Pl(22) +  &
      (-0.003547876846881965_Rkind) * tab_Pl(24) +  &
      (0.0037935717520674265_Rkind) * tab_Pl(26) +  &
      (-0.003370977874011358_Rkind) * tab_Pl(28) +  &
      (0.0019903812545488146_Rkind) * tab_Pl(30) +  &
      (0.0002660231202153635_Rkind) * tab_Pl(32) +  &
      (-0.0027567787885957656_Rkind) * tab_Pl(34) +  &
      (0.004760341218380286_Rkind) * tab_Pl(36) +  &
      (-0.005627991509208083_Rkind) * tab_Pl(38) +  &
      (0.005123336481187769_Rkind) * tab_Pl(40) +  &
      (-0.0037251294865690305_Rkind) * tab_Pl(42) +  &
      (0.002083326971132794_Rkind) * tab_Pl(44) +  &
      (-0.0008403379506349807_Rkind) * tab_Pl(46) +  &
      (0.00020809646457610622_Rkind) * tab_Pl(48)

      !(-0.1620716487102793_Rkind) * tab_Pl(0) +  &
      !(-0.008992725629195533_Rkind) * tab_Pl(2) +  &
      !(-0.00432915791758284_Rkind) * tab_Pl(4) +  &
      !(-0.001061111040246674_Rkind) * tab_Pl(6) +  &
      !(0.0003857228005224359_Rkind) * tab_Pl(8) +  &
      !(0.0007652895174396595_Rkind) * tab_Pl(10) +  &
      !(0.0006436145988118415_Rkind) * tab_Pl(12) +  &
      !(0.00037157203206966737_Rkind) * tab_Pl(14) +  &
      !(0.00011819866129507563_Rkind) * tab_Pl(16) +  &
      !(-4.02663761147522e-05_Rkind) * tab_Pl(18) +  &
      !(-0.00010232354419382081_Rkind) * tab_Pl(20) +  &
      !(-9.760506805819436e-05_Rkind) * tab_Pl(22) +  &
      !(-5.885890983333011e-05_Rkind) * tab_Pl(24) +  &
      !(-2.1436826791342742e-05_Rkind) * tab_Pl(26) +  &
      !(1.4784209082551226e-05_Rkind) * tab_Pl(28) +  &
      !(2.4150163942799644e-05_Rkind) * tab_Pl(30) +  &
      !(-3.501962307059586e-05_Rkind) * tab_Pl(32)


      
      ! c2 (c2*drho^2) fitté atour de rhoOPt à différentes val de s, puis fitté le long de s
      ! maxdiff = 3.29572217026e-05
      Func(2) =  & !ces coeff correspondent aux pts enlevé entre 6.6 et 3.6 et 3.6 et 2.8
      (0.30374950443390847_Rkind) * tab_Pl(0) +  &
      (-0.05263777329872866_Rkind) * tab_Pl(2) +  &
      (-0.02824499169155245_Rkind) * tab_Pl(4) +  &
      (-0.013869244641591425_Rkind) * tab_Pl(6) +  &
      (-0.008403214399978586_Rkind) * tab_Pl(8) +  &
      (-0.007696945284070528_Rkind) * tab_Pl(10) +  &
      (-0.003806082903400655_Rkind) * tab_Pl(12) +  &
      (-0.006166450639817748_Rkind) * tab_Pl(14) +  &
      (0.004309716438982306_Rkind) * tab_Pl(16) +  &
      (-0.008859755489454838_Rkind) * tab_Pl(18) +  &
      (0.015668516816827203_Rkind) * tab_Pl(20) +  &
      (-0.020386609731335353_Rkind) * tab_Pl(22) +  &
      (0.026726056241208166_Rkind) * tab_Pl(24) +  &
      (-0.0288510778775844_Rkind) * tab_Pl(26) +  &
      (0.024431057911522034_Rkind) * tab_Pl(28) +  &
      (-0.014898755884266347_Rkind) * tab_Pl(30) +  &
      (-0.0023882469597527427_Rkind) * tab_Pl(32) +  &
      (0.020524707114045054_Rkind) * tab_Pl(34) +  &
      (-0.035159777872531876_Rkind) * tab_Pl(36) +  &
      (0.0419505860292803_Rkind) * tab_Pl(38) +  &
      (-0.038291605815803_Rkind) * tab_Pl(40) +  &
      (0.028667463162531014_Rkind) * tab_Pl(42) +  &
      (-0.01659056531057078_Rkind) * tab_Pl(44) +  &
      (0.007217862809646366_Rkind) * tab_Pl(46) +  &
      (-0.00203367463365633_Rkind) * tab_Pl(48)

      !Coeff associé au fichier 'recap_optRhoEqui' dans dossier Fitpy
      !(0.30375204657676236_Rkind) * tab_Pl(0) +  &
      !(-0.05261249923316819_Rkind) * tab_Pl(2) +  &
      !(-0.028245088307616156_Rkind) * tab_Pl(4) +  &
      !(-0.01376309967876779_Rkind) * tab_Pl(6) +  &
      !(-0.008615637365472464_Rkind) * tab_Pl(8) +  &
      !(-0.007122038317701137_Rkind) * tab_Pl(10) +  &
      !(-0.005141765724950869_Rkind) * tab_Pl(12) +  &
      !(-0.003418104118376876_Rkind) * tab_Pl(14) +  &
      !(-0.0010447833949072116_Rkind) * tab_Pl(16) +  &
      !(0.00040790797785053524_Rkind) * tab_Pl(18) +  &
      !(0.0009037612950454548_Rkind) * tab_Pl(20) +  &
      !(0.0005452967481115013_Rkind) * tab_Pl(22) +  &
      !(0.000282087543543148_Rkind) * tab_Pl(24) +  &
      !(-0.0002715479728813458_Rkind) * tab_Pl(26) +  &
      !(-0.0007667274135405311_Rkind) * tab_Pl(28) +  &
      !(-0.0007199692637246045_Rkind) * tab_Pl(30) +  &
      !(0.0006760031365286_Rkind) * tab_Pl(32)


 

      ! rho
      Func(3) = &
      (1.261243931645815_Rkind) * tab_Pl(0) +  &
      (0.05498541052261184_Rkind) * tab_Pl(2) +  &
      (0.03446819696272064_Rkind) * tab_Pl(4) +  &
      (0.020291272930553106_Rkind) * tab_Pl(6) +  &
      (0.011855010123711414_Rkind) * tab_Pl(8) +  &
      (0.007536688407402705_Rkind) * tab_Pl(10) +  &
      (0.004940560392031361_Rkind) * tab_Pl(12) +  &
      (0.0034254508799222323_Rkind) * tab_Pl(14) +  &
      (0.002319413680131285_Rkind) * tab_Pl(16) +  &
      (0.0010470828509148065_Rkind) * tab_Pl(18) +  &
      (4.281256394175054e-05_Rkind) * tab_Pl(20) +  &
      (-0.00030975651066805104_Rkind) * tab_Pl(22) +  &
      (-0.0005399637278206589_Rkind) * tab_Pl(24) +  &
      (-0.00041279694640769006_Rkind) * tab_Pl(26) +  &
      (-0.00017611772257669623_Rkind) * tab_Pl(28) +  &
      (-2.5934623118451503e-05_Rkind) * tab_Pl(30) +  &
      (0.00036253207816267066_Rkind) * tab_Pl(32)
      
      ! c3(c3*drho) fitté atour de rhoOPt à différentes val de s, puis fitté le long de s 
      ! maxdiff = 7.1548287231e-06 
      Func(4) = &  !ces coeff correspondent aux pts enlevé entre 6.6 et 3.6 et 3.6 et 2.8 
      (0.00012063674082565115_Rkind) * tab_Pl(0) +  &
      (-4.1425739977470036e-05_Rkind) * tab_Pl(2) +  &
      (-2.4239116984059995e-06_Rkind) * tab_Pl(4) +  &
      (-1.6650686583068975e-05_Rkind) * tab_Pl(6) +  &
      (-8.079451360298356e-05_Rkind) * tab_Pl(8) +  &
      (8.429523860703813e-05_Rkind) * tab_Pl(10) +  &
      (-1.2203672188874432e-05_Rkind) * tab_Pl(12) +  &
      (-7.800676003009566e-05_Rkind) * tab_Pl(14) +  &
      (0.00015804831530465938_Rkind) * tab_Pl(16) +  &
      (-0.00016490891576258253_Rkind) * tab_Pl(18) +  &
      (0.00020743883827847572_Rkind) * tab_Pl(20) +  &
      (-0.0003393323849648792_Rkind) * tab_Pl(22) +  &
      (0.0004795101489175562_Rkind) * tab_Pl(24) +  &
      (-0.0005551463058609534_Rkind) * tab_Pl(26) +  &
      (0.0005039296106884514_Rkind) * tab_Pl(28) +  &
      (-0.00030089936753135966_Rkind) * tab_Pl(30) +  &
      (-4.1211848421952196e-05_Rkind) * tab_Pl(32) +  &
      (0.0004070924592144858_Rkind) * tab_Pl(34) +  &
      (-0.0006790178634567731_Rkind) * tab_Pl(36) +  &
      (0.0008075681707767094_Rkind) * tab_Pl(38) +  &
      (-0.0007497832421263899_Rkind) * tab_Pl(40) +  &
      (0.0005213860035859344_Rkind) * tab_Pl(42) +  &
      (-0.0002640804121683683_Rkind) * tab_Pl(44) +  &
      (9.866549309015685e-05_Rkind) * tab_Pl(46) +  &
      (-2.131966131723614e-05_Rkind) * tab_Pl(48)

      !Coeff associé au fichier 'recap_optRhoEqui' dans dossier Fitpy
      !(0.00012054894006886368_Rkind) * tab_Pl(0) +  &
      !(-4.159861752693137e-05_Rkind) * tab_Pl(2) +  &
      !(-3.088431820462534e-06_Rkind) * tab_Pl(4) +  &
      !(-1.5430069424391806e-05_Rkind) * tab_Pl(6) +  &
      !(-8.389126757210443e-05_Rkind) * tab_Pl(8) +  &
      !(9.425785955142528e-05_Rkind) * tab_Pl(10) +  &
      !(-3.255469762512486e-05_Rkind) * tab_Pl(12) +  &
      !(-3.114281736068148e-05_Rkind) * tab_Pl(14) +  &
      !(6.793677577444186e-05_Rkind) * tab_Pl(16) +  &
      !(-3.75084784501627e-06_Rkind) * tab_Pl(18) +  &
      !(-5.64014489843293e-05_Rkind) * tab_Pl(20) +  &
      !(4.147527388974017e-05_Rkind) * tab_Pl(22) +  &
      !(-1.5965985196154143e-05_Rkind) * tab_Pl(24) +  &
      !(-7.97292699073988e-06_Rkind) * tab_Pl(26) +  &
      !(1.0429200371081538e-05_Rkind) * tab_Pl(28) +  &
      !(-7.53503532862055e-06_Rkind) * tab_Pl(30) +  &
      !(6.714016001365854e-06_Rkind) * tab_Pl(32)

      
      ! c1 (c1*drho^3) fitté atour de rhoOPt à différentes val de s, puis fitté le long de s
      ! maxdiff = 5.49896226448e-05
      Func(5) = & !ces coeff correspondent aux pts enlevé entre 6.6 et 3.6 et 3.6 et 2.8
      (-0.42808970081556635_Rkind) * tab_Pl(0) +  &
      (0.11247166878455189_Rkind) * tab_Pl(2) +  &
      (0.059746640846453425_Rkind) * tab_Pl(4) +  &
      (0.019129889700503465_Rkind) * tab_Pl(6) +  &
      (0.0035237126723394446_Rkind) * tab_Pl(8) +  &
      (0.0074005177924327675_Rkind) * tab_Pl(10) +  &
      (-0.0033182465799101793_Rkind) * tab_Pl(12) +  &
      (0.023520802321938636_Rkind) * tab_Pl(14) +  &
      (-0.031316730641233526_Rkind) * tab_Pl(16) +  &
      (0.06314698041546835_Rkind) * tab_Pl(18) +  &
      (-0.10244257501198704_Rkind) * tab_Pl(20) +  &
      (0.14247360279431864_Rkind) * tab_Pl(22) +  &
      (-0.18657356209919096_Rkind) * tab_Pl(24) +  &
      (0.20197102516549023_Rkind) * tab_Pl(26) +  &
      (-0.18011147007817652_Rkind) * tab_Pl(28) +  &
      (0.10668386666389848_Rkind) * tab_Pl(30) +  &
      (0.014103648389836868_Rkind) * tab_Pl(32) +  &
      (-0.14637879317205393_Rkind) * tab_Pl(34) +  &
      (0.25445905188916274_Rkind) * tab_Pl(36) +  &
      (-0.2996700529026105_Rkind) * tab_Pl(38) +  &
      (0.272544782790872_Rkind) * tab_Pl(40) +  &
      (-0.19654051234395137_Rkind) * tab_Pl(42) +  &
      (0.10833980693031284_Rkind) * tab_Pl(44) +  &
      (-0.04322740125039924_Rkind) * tab_Pl(46) +  &
      (0.009767270370872093_Rkind) * tab_Pl(48)

      !Coeff associé au fichier 'recap_optRhoEqui' dans dossier Fitpy
      !(-0.4280712417313688_Rkind) * tab_Pl(0) +  &
      !(0.11247534169282206_Rkind) * tab_Pl(2) +  &
      !(0.05998037605796082_Rkind) * tab_Pl(4) +  &
      !(0.01870873415739227_Rkind) * tab_Pl(6) +  &
      !(0.005027009034584681_Rkind) * tab_Pl(8) +  &
      !(0.003776335643769335_Rkind) * tab_Pl(10) +  &
      !(0.005189361844932756_Rkind) * tab_Pl(12) +  &
      !(0.005285749904536103_Rkind) * tab_Pl(14) +  &
      !(0.0038340861120990787_Rkind) * tab_Pl(16) +  &
      !(0.0004693807417369628_Rkind) * tab_Pl(18) +  &
      !(-0.002231515868601842_Rkind) * tab_Pl(20) +  &
      !(-0.002355656089486632_Rkind) * tab_Pl(22) +  &
      !(-0.0018142662659441826_Rkind) * tab_Pl(24) +  &
      !(-0.0007192492715343131_Rkind) * tab_Pl(26) +  &
      !(0.0004425287789617937_Rkind) * tab_Pl(28) +  &
      !(0.001596012400134642_Rkind) * tab_Pl(30) +  &
      !(0.0007501493129836831_Rkind) * tab_Pl(32)


      ! c0 (c0*drho^4) fitté atour de rhoOPt à différentes val de s, puis fitté le long de s
      ! maxdiff = 0.000715618043238
      Func(6) = & !ces coeff correspondent aux pts enlevé entre 6.6 et 3.6 et 3.6 et 2.8
      (0.3878984339121209_Rkind) * tab_Pl(0) +  &
      (-0.11984756856262713_Rkind) * tab_Pl(2) +  &
      (-0.062360429222014756_Rkind) * tab_Pl(4) +  &
      (-0.018791720873183277_Rkind) * tab_Pl(6) +  &
      (0.006886309376733128_Rkind) * tab_Pl(8) +  &
      (-0.011861079630599196_Rkind) * tab_Pl(10) +  &
      (0.02838129356224232_Rkind) * tab_Pl(12) +  &
      (-0.06706018831054988_Rkind) * tab_Pl(14) +  &
      (0.12618831678511594_Rkind) * tab_Pl(16) +  &
      (-0.2488281924518596_Rkind) * tab_Pl(18) +  &
      (0.3956917615053644_Rkind) * tab_Pl(20) +  &
      (-0.5703720557053324_Rkind) * tab_Pl(22) +  &
      (0.7358984543634043_Rkind) * tab_Pl(24) +  &
      (-0.8107843762174871_Rkind) * tab_Pl(26) +  &
      (0.7313663688956455_Rkind) * tab_Pl(28) +  &
      (-0.4306584398099331_Rkind) * tab_Pl(30) +  &
      (-0.04926423367854221_Rkind) * tab_Pl(32) +  &
      (0.5899530871873381_Rkind) * tab_Pl(34) +  &
      (-1.0295825892712576_Rkind) * tab_Pl(36) +  &
      (1.2079830149064956_Rkind) * tab_Pl(38) +  &
      (-1.0932292143867703_Rkind) * tab_Pl(40) +  &
      (0.7751899409504499_Rkind) * tab_Pl(42) +  &
      (-0.41669810982940203_Rkind) * tab_Pl(44) +  &
      (0.1559804684815059_Rkind) * tab_Pl(46) +  &
      (-0.030382992175241014_Rkind) * tab_Pl(48)

      !Coeff associé au fichier 'recap_optRhoEqui' dans dossier Fitpy
      !(0.38784157838415956_Rkind) * tab_Pl(0) +  &
      !(-0.11980773514140389_Rkind) * tab_Pl(2) +  &
      !(-0.06312819738591222_Rkind) * tab_Pl(4) +  &
      !(-0.017092165956158818_Rkind) * tab_Pl(6) +  &
      !(0.0013762595357292883_Rkind) * tab_Pl(8) +  &
      !(0.0019457694456141176_Rkind) * tab_Pl(10) +  &
      !(-0.004058733658769555_Rkind) * tab_Pl(12) +  &
      !(0.0028543497034913407_Rkind) * tab_Pl(14) +  &
      !(-0.010147768764494407_Rkind) * tab_Pl(16) +  &
      !(-0.004763434006431677_Rkind) * tab_Pl(18) +  &
      !(0.0017120023223853517_Rkind) * tab_Pl(20) +  &
      !(0.002508091645298231_Rkind) * tab_Pl(22) +  &
      !(-2.904687912089893e-07_Rkind) * tab_Pl(24) +  &
      !(0.002880671741730294_Rkind) * tab_Pl(26) +  &
      !(0.003464168323350602_Rkind) * tab_Pl(28) +  &
      !(0.0011655181631145695_Rkind) * tab_Pl(30) +  &
      !(-0.005870081331958774_Rkind) * tab_Pl(32)

      DO i=0,max_deg
        CALL dealloc_dnS(tab_Pl(i))
      END DO
      CALL dealloc_dnS(s)
      CALL dealloc_dnS(ts)
  END SUBROUTINE EvalFunc_QML_H3_v74

  FUNCTION QML_dnSigmoid_H3(x,sc)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t)                        :: QML_dnSigmoid_H3

    TYPE (dnS_t),         intent(in)    :: x
    real(kind=Rkind),     INTENT(IN)    :: sc

    QML_dnSigmoid_H3 = (ONE+tanh(sc*x))/TWO

  END FUNCTION QML_dnSigmoid_H3
END MODULE QML_H3_m
