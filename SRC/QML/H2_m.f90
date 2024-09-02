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

!> @brief Module which makes the initialization, calculation of the H2 potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 18/10/2021
!!
MODULE QML_H2_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the H2 parameters are set-up.
!> @brief Default parameters for H-H
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param mu             real: Reduced mass of HH (in au)
  TYPE, EXTENDS (QML_Empty_t) :: QML_H2_t
     PRIVATE
     real (kind=Rkind), PUBLIC      :: mu  = 1837.1526464003414_Rkind/TWO !< Reduced mass of H2 (in au)
     real (kind=Rkind)              :: Req = 1.4_Rkind
     real (kind=Rkind)              :: E0  = -1.17404460_Rkind ! CCSD(T)-F12B/VTZ-F12
     real (kind=Rkind), allocatable :: TaylorPot(:)

  CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_H2
    PROCEDURE :: Write_QModel    => Write_QML_H2
  END TYPE QML_H2_t

  PUBLIC :: QML_H2_t,Init_QML_H2,Write_QML_H2,QML_dnH2

CONTAINS
!> @brief Subroutine which makes the initialization of the H2 parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param H2Pot              TYPE(QML_H2_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_H2(QModel_in,read_param,nio_param_file) RESULT(QModel)
    IMPLICIT NONE

    !TYPE (QML_H2_t), allocatable                 :: QModel
    TYPE (QML_H2_t)                              :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H2'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) '  read_param:     ',read_param
      write(out_unit,*) '  nio_param_file: ',nio_param_file
      flush(out_unit)
    END IF

    !allocate(QML_H2_t :: QModel)
    QModel%QML_Empty_t = QModel_in

    QModel%pot_name = 'h2'
    QModel%ndim     = 1
    QModel%nsurf    = 1


    IF (debug) write(out_unit,*) 'init default H2 parameters (HF)'
    QModel%E0  = -1.17404460_Rkind
    QModel%Req = 1.4_Rkind

    IF (read_param) THEN
      CALL Read_QML_H2(QModel%Req,QModel%E0,nio_param_file)
    ELSE
      IF (debug) write(out_unit,*) 'init H2 parameters, if present'
    END IF

    SELECT CASE (QModel%option)
    CASE (1)
      write(out_unit,*) 'Fourth order expansion in (R-Req)'

      QModel%TaylorPot = [-0.000498716666667_Rkind,0.185668083333_Rkind,-0.219403333333_Rkind,0.182041666667_Rkind]

      IF (debug) write(out_unit,*) 'init Q0 of H2'
      QModel%Q0 = [QModel%req]

      IF (debug) write(out_unit,*) 'init d0GGdef of H2'
      QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])
      QModel%no_ana_der = .FALSE.

    CASE(2) ! The coordinate x is defined as x = Req/R
      write(out_unit,*) 'Fourth order expansion in (x-1) with x = Req/R'

      ! with DQ=+/-0.1 and +/- 0.2
      !QModel%TaylorPot = [-0.000753325_Rkind,0.363249875_Rkind,-0.1412875_Rkind,0.0165625_Rkind]

      ! with DQ=+/-0.05 and +/- 0.1
      QModel%TaylorPot = [0.00077465_Rkind,0.363216166667_Rkind,-0.14342_Rkind,0.019933333338_Rkind]

      QModel%Q0 = [ONE]
      ! 1/mu d2/dR^2  =>   1/mu [x" d/dx + x'^2 d2/dx^2]
      !   => G(x) = 1/mu * x'^2 = 1/mu * (-Req/R^2)^2 = 1/mu * (x^2/Req)^2
      !  and xeq = 1.
      QModel%d0GGdef = reshape([ONE/QModel%req**2 * ONE/QModel%mu],shape=[1,1])
      QModel%no_ana_der = .FALSE.
    CASE(3) ! H2 from H+H2 LSTH potential
      QModel%E0  = ZERO
      QModel%Req = 1.4_Rkind

      QModel%Q0 = [QModel%req]

      IF (debug) write(out_unit,*) 'init d0GGdef of H2'
      QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])
      QModel%no_ana_der = .TRUE.

    CASE default
      write(out_unit,*) 'ERROR in Init_QML_H2'
      write(out_unit,*) 'Possible option: 1,2,3'
      write(out_unit,*) 'Actual value, option=',QModel%option
      STOP 'ERROR in Init_QML_H2: wrong option value'
    END SELECT

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_H2

!> @brief Subroutine wich reads the H2 parameters with a namelist.
!!   This can be called only from the "Init_QML_H2" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_H2_t):   derived type in which the parameters are set-up.
!! @param nio                integer          :   file unit to read the parameters.
  SUBROUTINE Read_QML_H2(Req_inout,E0_inout,nio)
    IMPLICIT NONE

    integer,           intent(in)    :: nio
    real (kind=Rkind), intent(inout) :: Req_inout,E0_inout

    !local variable
    real (kind=Rkind) :: Req,E0
    integer           :: err_read

    namelist /H2/ Req,E0

    E0  = E0_inout
    Req = Req_inout
    write(out_unit,*) 'read H2 namelist ...'

    read(nio,nml=H2,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_H2'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "H2" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_H2'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_H2'
      write(out_unit,*) ' Some parameter names of the namelist "H2" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=H2)
      STOP ' ERROR in Read_QML_H2'
    END IF
    !write(out_unit,nml=H2)

    E0_inout   = E0
    Req_inout  = Req

  END SUBROUTINE Read_QML_H2
!! === README ==
!! H2 potential
!! pot_name  = 'H2'
!! ndim      = 1
!! nsurf     = 1
!!
!! options, (1) Talyor expansion: $V(R) = \sum_i a_i \cdot (R-Req)^{i-1}$
!!     Level: CCSD(T)-F12B/VTZ-F12 (with molpro 2010)
!! options, (2) Talyor expansion: $V(R) = \sum_i a_i \cdot x^{i-1}$ with $x=Req/R-1$
!!     Level: CCSD(T)-F12B/VTZ-F12 (with molpro 2010)
!! options(3) extract for the H+H2 LSTH potential
!!
!! reduced mass      = 1837.1526464003414/2 au
!! === END README ==
!> @brief Subroutine wich prints the H2 current parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_H2_t):   derived type with the H2 parameters.
!! @param nio            integer          :   file unit to print the parameters.
  SUBROUTINE Write_QML_H2(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_H2_t),    intent(in) :: QModel
    integer,            intent(in) :: nio


    write(nio,*) 'H2 current parameters:'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) '    V(R) = E0 + Sum_i a_i(R-Req)^i'
    write(nio,*) '  E0:            ',QModel%E0
    write(nio,*) '  Req:           ',QModel%Req
    write(nio,*) '  TaylorPot(:)   ',QModel%TaylorPot(:)
    write(nio,*)
    write(nio,*) 'end H2 current parameters'

  END SUBROUTINE Write_QML_H2

!> @brief Subroutine wich calculates the H2 potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      Potential with derivatives,.
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_H2_t):   derived type with the H2 parameters.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H2(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2_t),     intent(in)     :: QModel
    TYPE(dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         intent(in)     :: dnQ(:)
    integer,             intent(in)     :: nderiv

    integer :: i
    real (kind=Rkind) :: R,S1(3)
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_H2'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) ' nderiv:',nderiv
      write(out_unit,*) ' Q(:):',(get_d0(dnQ(i)),i=1,size(dnQ))
    END IF

    SELECT CASE (QModel%option)
    CASE (1,2)
      Mat_OF_PotDia(1,1) = QML_dnH2(dnQ(1),QModel)
    CASE(3)
      R = get_d0(dnQ(1))
      CALL QMLH2_VH2(R,S1)
      Mat_OF_PotDia(1,1) = S1(1)
    END SELECT

    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
  END SUBROUTINE EvalPot_QML_H2

!> @brief Function wich calculates the H2 potential with derivatives.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return QML_dnH2      TYPE (dnS_t):     derived type with a value (pot), if required, its derivatives.
!! @param dnQ            TYPE (dnS_t):     derived type with the value of "r" and,if required, its derivatives.
!! @param QModel         TYPE(QML_H2_t):   derived type with the H2 parameters.
  FUNCTION QML_dnH2(dnQ,QModel)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t)                     :: QML_dnH2

    TYPE (dnS_t),      intent(in)    :: dnQ
    TYPE (QML_H2_t),   intent(in)    :: QModel

    !local variable
    TYPE (dnS_t)        :: dnDeltaQ
    integer             :: i
    !logical, parameter  :: debug = .TRUE.
    logical, parameter  :: debug = .FALSE.

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING in QML_dnH2'
      CALL Write_dnS(dnQ,info='dnQ')
    END IF

    SELECT CASE (QModel%option)
    CASE (1)
      dnDeltaQ  = dnQ-QModel%Req
    CASE(2)
      dnDeltaQ  = QModel%Req/dnQ-ONE
    END SELECT
    IF (debug) CALL Write_dnS(dnDeltaQ,info='dnDeltaQ')

    QML_dnH2 = QModel%E0
    DO i=1,size(QModel%TaylorPot)
      QML_dnH2 = QML_dnH2 + QModel%TaylorPot(i)*dnDeltaQ**i
    END DO

    CALL dealloc_dnS(dnDeltaQ)

    IF (debug) THEN
      write(out_unit,*) 'H2 at',get_d0(dnQ)
      CALL Write_dnS(QML_dnH2)
      write(out_unit,*) 'END in QML_dnH2'
      flush(out_unit)
    END IF

  END FUNCTION QML_dnH2


  SUBROUTINE QMLH2_VH2(X,S1)
    IMPLICIT NONE
        real(kind=Rkind), intent(in)    :: X
        real(kind=Rkind), intent(inout) :: S1(3)


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
  
        !COMMON/POTCOM/C6,C8,RKW(87),EKW(87),WKW(87)
        IF (X > TEN) THEN
          CALL QMLH2_VBIGR(X,S1)
        ELSE
          CALL QMLH2_SPLID2(87,RKW,EKW,WKW,1,X,S1)
        END IF

  END SUBROUTINE QMLH2_VH2
!       *************************************************************
  SUBROUTINE QMLH2_VBIGR(X,S)
    IMPLICIT NONE
        real(kind=Rkind), intent(in)    ::  X
        real(kind=Rkind), intent(inout) ::  S(3)

        real(kind=Rkind)   ::  X2,X3,X6,C8A

        !COMMON/POTCOM/C6,C8,RKW(87),EKW(87),WKW(87)

        real(kind=Rkind), parameter :: C6  = 6.89992032_Rkind
        real(kind=Rkind), parameter :: C8  = 219.9997304_Rkind
    
        X2=X*X
        X3=X2*X
        X6=X3*X3
        C8A=C8/X2
        S(1)=-(C6+C8A)/X6
        S(2)=(C6*SIX + C8A*EIGHT)/X6/X

  END SUBROUTINE QMLH2_VBIGR

  SUBROUTINE QMLH2_SPLID2(N,X,FF,W,IJ,Y,TAB)
    IMPLICIT NONE

        integer,          intent(in)    :: N,IJ
        real(kind=Rkind), intent(in)    :: X(N)
        real(kind=Rkind), intent(in)    :: Y,FF(N),W(N)
        real(kind=Rkind), intent(inout) :: TAB(3)

        INTEGER          :: I,K,MI,KI
        real(kind=Rkind) :: FLK,AA,BB,CC

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
  END SUBROUTINE QMLH2_SPLID2
END MODULE QML_H2_m
