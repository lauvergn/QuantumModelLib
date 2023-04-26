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

!> @brief Module which makes the initialization, calculation of the CNH_Murrell potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_CNH_Murrell_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the CNH_Murrell parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_CNH_Murrell_t

   PRIVATE

   CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_CNH_Murrell
    PROCEDURE :: Write_QModel     => Write_QML_CNH_Murrell
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_QML_CNH_Murrell
    PROCEDURE :: EvalFunc_QModel  => EvalFunc_QML_CNH_Murrell
  END TYPE QML_CNH_Murrell_t

  PUBLIC :: QML_CNH_Murrell_t,Init_QML_CNH_Murrell


  CONTAINS
!> @brief Function which makes the initialization of the CNH_Murrell parameters.
!!
!! @param QModel             TYPE(QML_CNH_Murrell_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_CNH_Murrell(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    !TYPE (QML_CNH_Murrell_t), allocatable        :: QModel ! RESULT
    TYPE (QML_CNH_Murrell_t)                     :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_CNH_Murrell'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    !allocate(QML_CNH_Murrell_t :: QModel)

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf      = 1
    QModel%ndimCart   = 9

    QModel%ndimQ      = 3

    IF (QModel%Cart_TO_Q)    QModel%option = 0
    IF (QModel%option == -1) QModel%option = 0

    SELECT CASE(QModel%option)
    CASE (0,1,11) !

      IF (QModel%Cart_TO_Q) THEN
        QModel%ndim       = QModel%ndimCart
      ELSE
        QModel%ndim       = QModel%ndimQ
      END IF

      ! so that, we can get the MEP function with the 3D model
      QModel%ndimFunc       = 1
      QModel%nb_Func        = 1 + 2 + 4 ! V,Req,Hess
      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 2
      QModel%IndexFunc_Hess = 4

      QModel%pot_name       = 'CNH_Murrell'
    CASE (2,21)
      QModel%ndimQ          = 1
      QModel%ndim           = 1

      ! so that, we can get the MEP function with the 3D model
      QModel%ndimFunc       = 1
      QModel%nb_Func        = 1 + 2 + 4 ! V,Req,Hess
      QModel%IndexFunc_Ene  = 1
      QModel%IndexFunc_Qop  = 2
      QModel%IndexFunc_Hess = 4

      QModel%pot_name       = 'CNH_Murrell_MEP'
    CASE Default
       write(out_unit,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unit)
       write(out_unit,*) ' ERROR in ',name_sub
       write(out_unit,*) ' option: ',QModel%option
       write(out_unit,*) ' the possible option are: 0, 1, 11 (3D) or 2,21 (1D-MEP)'
       STOP 'ERROR in Init_QML_CNH_Murrell: wrong option'
    END SELECT


    IF (debug) write(out_unit,*) 'init Q0 of CNH_Murrell (CNH_Murrell minimum)'
    QModel%Q0 = [ZERO,3.18722_Rkind,2.17926_Rkind]

    IF (debug) write(out_unit,*) 'init d0GGdef of CNH_Murrell'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)


    IF (debug) THEN
      CALL Write_QML_CNH_Murrell(QModel,nio=out_unit)
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_CNH_Murrell
!> @brief Subroutine wich prints the current QML_CNH_Murrell parameters.
!!
!! @param QModel            CLASS(QML_CNH_Murrell_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_CNH_Murrell(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_CNH_Murrell_t),  intent(in) :: QModel
    integer,                   intent(in) :: nio

    write(nio,*) 'CNH_Murrell current parameters'
    write(nio,*) 'ref: J. N. Murrell, S. Carter, and L. O. Halonen, ...'
    write(nio,*) '  ... J. Mol. Spectrosc. 93,307 (1982).'
    write(nio,*) 'https://doi.org/10.1016/0022-2852(82)90170-9'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*) 'Atomic order : C,N,H'
    SELECT CASE(QModel%option)
    CASE (0) ! 3D (three distances)
      write(nio,*) '3D with 3 distances'
    CASE (1) ! 3D Jacobi
      write(nio,*) '3D With Jacobi coordinates: Theta, R_H-CN, R_CN'
    CASE (11) ! 3D Jacobi
      write(nio,*) '3D With Jacobi coordinates: cos(Theta), R_H-CN, R_CN'
    CASE (2) ! MEP
      write(nio,*) '1D-MEP: Theta'
    CASE (21) ! MEP
      write(nio,*) '1D-MEP: cos(Theta)'
    END SELECT

    write(nio,*) 'end CNH_Murrell current parameters'

  END SUBROUTINE Write_QML_CNH_Murrell


!> @brief Subroutine wich calculates the CNH_Murrell potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_CNH_Murrell_t):    derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_CNH_Murrell(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_CNH_Murrell_t),  intent(in)    :: QModel
    TYPE (dnS_t),              intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),              intent(in)    :: dnQ(:)
    integer,                   intent(in)    :: nderiv

    TYPE (dnS_t)    :: Func(QModel%nb_Func)


    SELECT CASE(QModel%option)
    CASE (0) ! 3D (three distances)
      CALL QML_EvalPot3D_Murrell(Mat_OF_PotDia,dnQ,.FALSE.,.FALSE.)
    CASE (1) ! 3D Jacobi: Theta (H-CN), RH-CN, RCN
      CALL QML_EvalPot3D_Murrell(Mat_OF_PotDia,dnQ,.TRUE.,.FALSE.)
    CASE (11) ! 3D Jacobi : cos(Theta) (H-CN), RH-CN, RCN
      CALL QML_EvalPot3D_Murrell(Mat_OF_PotDia,dnQ,.TRUE.,.TRUE.)
    CASE (2,21) ! MEP : Theta (H-CN) or cos(Theta) (H-CN)
      CALL EvalFunc_QML_CNH_Murrell(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)
    END SELECT


  END SUBROUTINE EvalPot_QML_CNH_Murrell

  SUBROUTINE QML_EvalPot3D_Murrell(Mat_OF_PotDia,dnQ,Jacobi,cosTh)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),            intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),            intent(in)    :: dnQ(:)
    logical,                 intent(in)    :: Jacobi,cosTh

    real(kind=Rkind), parameter :: mB = 12._Rkind
    real(kind=Rkind), parameter :: mC = 14.003074_Rkind
    real(kind=Rkind), parameter :: mH = 1.007825_Rkind

    real(kind=Rkind), parameter :: autoA=0.52917715_Rkind
    real(kind=Rkind), parameter :: autoeV=27.21183_Rkind

    TYPE (dnS_t)  :: gR,pR
    TYPE (dnS_t)  :: T,RN,RC,R1,R2,R3

    TYPE (dnS_t)  :: Z1,Z12,V1
    TYPE (dnS_t)  :: Z2,Z22,V2
    TYPE (dnS_t)  :: Z3,Z32,V3

    TYPE (dnS_t)  :: S1,S2,S3,S12,S22,S32,POLY
    TYPE (dnS_t)  :: E1,E2,E3,HY1,HY2,HY3,SWITCH,V123

    IF (jacobi) THEN
      pR = dnQ(2)
      gR = dnQ(3)

      IF (cosTh) THEN
        T = dnQ(1)
      ELSE
        T = cos(dnQ(1))
      END IF
      RN=mC/(mB+mC)
      RC=ONE-RN
      R1=sqrt(pR**2+(RN*gR)**2-TWO*pR*gR*RN*T)
      R2=gR
      R3=sqrt(pR**2+(RC*gR)**2+TWO*pR*gR*RC*T)
    ELSE
      R1 = dnQ(1)
      R2 = dnQ(2)
      R3 = dnQ(3)
    END IF
    !write(out_unit,*) 'R1,R2,R3',get_d0(R1),get_d0(R2),get_d0(R3)

    R1=R1*autoA
    R2=R2*autoA
    R3=R3*autoA

    !.....CH
    Z1  = R1-1.0823_Rkind
    Z12 = Z1*Z1
    V1  = -2.8521_Rkind*(ONE+5.5297_Rkind*Z1+8.7166_Rkind*Z12+5.3082_Rkind*Z1*Z12)*&
          exp(-5.5297_Rkind*Z1)
    !.....CN
    Z2  = R2-1.1718_Rkind
    Z22 = Z2*Z2
    V2  = -7.9282_Rkind*(ONE+5.2448_Rkind*Z2+7.3416_Rkind*Z22+4.9785_Rkind*Z2*Z22)*&
          exp(-5.2448_Rkind*Z2)
    !.....NH
    Z3  = R3-1.0370_Rkind
    Z32 = Z3*Z3
    V3  = -3.9938_Rkind*(ONE+3.0704_Rkind*Z3)*exp(-3.0704_Rkind*Z3)

    !.....THREE BODY TERMS
    Z1 = R1-1.9607_Rkind
    Z2 = R2-2.2794_Rkind
    Z3 = R3-1.8687_Rkind
    S1 = 0.4436_Rkind*Z1+0.6091_Rkind*Z2+0.6575_Rkind*Z3
    S2 = -.8941_Rkind*Z1+0.2498_Rkind*Z2+0.3718_Rkind*Z3
    S3 = 0.0622_Rkind*Z1-0.7527_Rkind*Z2+0.6554_Rkind*Z3

    S12  = S1*S1
    S22  = S2*S2
    S32  = S3*S3
    POLY = -3.0578_Rkind*(ONE+1.9076_Rkind*S1-0.5008_Rkind*S2-                  &
          0.0149_Rkind*S3+0.6695_Rkind*S12-                                     &
          1.3535_Rkind*S22-1.0501_Rkind*S32+                                    &
          0.2698_Rkind*S1*S2-1.1120_Rkind*S1*S3+                                &
          1.9310_Rkind*S2*S3-0.0877_Rkind*S1*S12+                               &
          0.0044_Rkind*S2*S22+0.0700_Rkind*S3*S32+                              &
          0.0898_Rkind*S12*S2-1.0186_Rkind*S1*S22-                              &
          0.0911_Rkind*S12*S3+                                                  &
          0.0017_Rkind*S1*S32+0.4567_Rkind*S22*S3-                              &
          0.8840_Rkind*S2*S32+                                                  &
          0.3333_Rkind*S1*S2*S3-0.0367_Rkind*S12*S12+                           &
          0.4821_Rkind*S22*S22+                                                 &
          0.2564_Rkind*S32*S32-0.0017_Rkind*S12*S1*S2-                          &
          0.2278_Rkind*S12*S22-                                                 &
          0.1287_Rkind*S1*S2*S22+0.1759_Rkind*S1*S12*S3-                        &
          0.0399_Rkind*S12*S32-                                                 &
          0.1447_Rkind*S1*S3*S32-0.3147_Rkind*S2*S22*S3+                        &
          0.1233_Rkind*S22*S32+                                                 &
          0.3161_Rkind*S2*S3*S32+0.0919_Rkind*S12*S2*S3-                        &
          0.0954_Rkind*S1*S22*S3+                                               &
          0.1778_Rkind*S1*S2*S32-0.1892_Rkind*S22*S22*S2)


    E1     = exp(3.9742_Rkind*Z1/TWO)
    E2     = exp(4.3688_Rkind*Z2/TWO)
    E3     = exp(1.5176_Rkind*Z3/TWO)
    HY1    = (E1-ONE/E1)/(E1+ONE/E1)
    HY2    = (E2-ONE/E2)/(E2+ONE/E2)
    HY3    = (E3-ONE/E3)/(E3+ONE/E3)
    SWITCH = (ONE-HY1)*(ONE-HY2)*(ONE-HY3)
    V123   = SWITCH*POLY

    Mat_OF_PotDia(1,1) = (V1+V2+V3+V123)/autoeV

  END SUBROUTINE QML_EvalPot3D_Murrell

  ! here we assume the C,N,H atomic order
  SUBROUTINE Cart_TO_Q_QML_CNH_Murrell(QModel,dnX,dnQ,nderiv)
    USE QDUtil_m,         ONLY : TO_string
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_CNH_Murrell_t), intent(in)    :: QModel
    TYPE (dnS_t),             intent(in)    :: dnX(:,:)
    TYPE (dnS_t),             intent(inout) :: dnQ(:)
    integer,                  intent(in)    :: nderiv

    ! local vector
    integer         :: i,j
    TYPE (dnS_t), allocatable :: Vec12(:),Vec13(:),Vec23(:)

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Cart_TO_Q_QML_CNH_Murrell'
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

    Vec23(:) = dnX(:,3)-dnX(:,2) ! NH
    Vec12(:) = dnX(:,2)-dnX(:,1) ! CN
    Vec13(:) = dnX(:,3)-dnX(:,1) ! CH

    IF (debug) THEN
      write(out_unit,*) 'Cart_TO_Q_QML_CNH_Murrell vect done'
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
    dnQ(1) = sqrt(dot_product(Vec13,Vec13)) ! CH
    dnQ(2) = sqrt(dot_product(Vec12,Vec12)) ! CN
    dnQ(3) = sqrt(dot_product(Vec23,Vec23)) ! NH

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
  END SUBROUTINE Cart_TO_Q_QML_CNH_Murrell


  ! for MEP
  SUBROUTINE EvalFunc_QML_CNH_Murrell(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_CNH_Murrell_t), intent(in)    :: QModel
    TYPE (dnS_t),             intent(inout) :: Func(:)
    TYPE (dnS_t),             intent(in)    :: dnQ(:)
    integer,                  intent(in)    :: nderiv


    TYPE (dnS_t)                  :: dnCth
    integer                       :: i
    integer,           parameter  :: max_deg = 29
    TYPE (dnS_t)                  :: tab_Pl(0:max_deg)


    IF (QModel%option == 11 .OR. QModel%option == 21) THEN
      dnCth = dnQ(1)
    ELSE
      dnCth = cos(dnQ(1))
    END IF

    DO i=0,max_deg
      tab_Pl(i) = dnLegendre0(dnCth,i,ReNorm=.TRUE.)
    END DO

    ! potential: MEP
    Func(1) = (-0.655285009791138906_Rkind)     * tab_Pl( 0) +    &
              (0.368866311567775168E-02_Rkind)  * tab_Pl( 1) +    &
              (-0.124410758695529707E-01_Rkind) * tab_Pl( 2) +    &
              (-0.990511349333085718E-02_Rkind) * tab_Pl( 3) +    &
              (-0.500503540822320580E-02_Rkind) * tab_Pl( 4) +    &
              (0.270621059491787551E-02_Rkind)  * tab_Pl( 5) +    &
              (0.141103658983222240E-02_Rkind)  * tab_Pl( 6) +    &
              (-0.437949934617522837E-03_Rkind) * tab_Pl( 7) +    &
              (-0.252460909929458637E-03_Rkind) * tab_Pl( 8) +    &
              (-0.358625641419419249E-04_Rkind) * tab_Pl( 9) +    &
              (0.976619474825319656E-04_Rkind)  * tab_Pl(10) +    &
              (0.217366800231287497E-04_Rkind)  * tab_Pl(11) +    &
              (-0.268212617015022991E-04_Rkind) * tab_Pl(12) +    &
              (-0.342047081164582388E-05_Rkind) * tab_Pl(13) +    &
              (-0.362810339599638992E-06_Rkind) * tab_Pl(14) +    &
              (0.259881013065570441E-05_Rkind)  * tab_Pl(15) +    &
              (0.257930636674634998E-05_Rkind)  * tab_Pl(16) +    &
              (-0.158479604070947891E-05_Rkind) * tab_Pl(17) +    &
              (-0.920544265357076132E-06_Rkind) * tab_Pl(18) +    &
              (0.338251259556195095E-06_Rkind)  * tab_Pl(19) +    &
              (0.285425243520912436E-06_Rkind)  * tab_Pl(20) +    &
              (0.425409677567662173E-07_Rkind)  * tab_Pl(21) +    &
              (-0.101054782382119211E-06_Rkind) * tab_Pl(22) +    &
              (-0.654484114758141540E-07_Rkind) * tab_Pl(23) +    &
              (0.335637142899206936E-07_Rkind)  * tab_Pl(24) +    &
              (0.463079565321430222E-07_Rkind)  * tab_Pl(25) +    &
              (-0.336153084693118319E-07_Rkind) * tab_Pl(26) +    &
              (-0.158611305955054091E-07_Rkind) * tab_Pl(27) +    &
              (0.134961732975956485E-07_Rkind)  * tab_Pl(28) +    &
              (0.106816325366728205E-07_Rkind)  * tab_Pl(29)

    ! Rop : R_H-CN
    Func(2) = (3.35875511812137484_Rkind) * tab_Pl(0) +    &
              (0.148871082021365564_Rkind) * tab_Pl(1) +    &
              (0.481094554409106090_Rkind) * tab_Pl(2) +    &
              (-0.253922825825884796E-01_Rkind) * tab_Pl(3) +    &
              (-0.537271857241523043E-01_Rkind) * tab_Pl(4) +    &
              (0.287137955899733916E-02_Rkind) * tab_Pl(5) +    &
              (0.883287278072170806E-02_Rkind) * tab_Pl(6) +    &
              (0.341914213510483619E-02_Rkind) * tab_Pl(7) +    &
              (-0.143761468388556584E-02_Rkind) * tab_Pl(8) +    &
              (-0.155642376002690548E-02_Rkind) * tab_Pl(9) +    &
              (-0.243938464552128774E-03_Rkind) * tab_Pl(10) +    &
              (0.324046391033408054E-03_Rkind) * tab_Pl(11) +    &
              (0.452748654304650116E-03_Rkind) * tab_Pl(12) +    &
              (-0.832793760270137800E-04_Rkind) * tab_Pl(13) +    &
              (-0.208508529872799653E-03_Rkind) * tab_Pl(14) +    &
              (0.188747686492221380E-04_Rkind) * tab_Pl(15) +    &
              (0.460225477732536980E-04_Rkind) * tab_Pl(16) +    &
              (0.115136940873498727E-04_Rkind) * tab_Pl(17) +    &
              (-0.227801221698345384E-05_Rkind) * tab_Pl(18) +    &
              (-0.113574404267734137E-04_Rkind) * tab_Pl(19) +    &
              (-0.349582880383952149E-05_Rkind) * tab_Pl(20) +    &
              (0.404005631088145106E-05_Rkind) * tab_Pl(21) +    &
              (0.187086530648967789E-05_Rkind) * tab_Pl(22) +    &
              (-0.474547142158161413E-06_Rkind) * tab_Pl(23) +    &
              (-0.112906549967473481E-05_Rkind) * tab_Pl(24) +    &
              (-0.537314466513270506E-06_Rkind) * tab_Pl(25) +    &
              (0.702706581417238486E-06_Rkind) * tab_Pl(26) +    &
              (-0.107428406076785338E-05_Rkind) * tab_Pl(27) +    &
              (0.348530750491854912E-07_Rkind) * tab_Pl(28) +    &
              (0.108653261459408565E-05_Rkind) * tab_Pl(29)

    ! Rop : R_CN
    Func(3) = (3.05567618651521178_Rkind) * tab_Pl(0) +    &
              (-0.387552561077800402E-02_Rkind) * tab_Pl(1) +    &
              (0.362104655719323598E-01_Rkind) * tab_Pl(2) +    &
              (-0.530589035832196559E-02_Rkind) * tab_Pl(3) +    &
              (-0.150883035992536076E-01_Rkind) * tab_Pl(4) +    &
              (0.180163772289059442E-02_Rkind) * tab_Pl(5) +    &
              (0.149572376873839718E-02_Rkind) * tab_Pl(6) +    &
              (0.805635074701905652E-04_Rkind) * tab_Pl(7) +    &
              (0.176212224375385838E-03_Rkind) * tab_Pl(8) +    &
              (-0.210362333581846055E-03_Rkind) * tab_Pl(9) +    &
              (-0.169109695229127790E-03_Rkind) * tab_Pl(10) +    &
              (0.803399445418478573E-04_Rkind) * tab_Pl(11) +    &
              (0.921372890564473051E-04_Rkind) * tab_Pl(12) +    &
              (-0.210460594038160419E-04_Rkind) * tab_Pl(13) +    &
              (-0.324489980716477540E-04_Rkind) * tab_Pl(14) +    &
              (0.117169011613050838E-05_Rkind) * tab_Pl(15) +    &
              (0.474614394156325088E-05_Rkind) * tab_Pl(16) +    &
              (0.399296936512522564E-05_Rkind) * tab_Pl(17) +    &
              (0.153363120618120103E-05_Rkind) * tab_Pl(18) +    &
              (-0.353877162048743810E-05_Rkind) * tab_Pl(19) +    &
              (-0.238404368608901790E-05_Rkind) * tab_Pl(20) +    &
              (0.543490026077864873E-06_Rkind) * tab_Pl(21) +    &
              (0.108342236176666158E-05_Rkind) * tab_Pl(22) +    &
              (0.343696074723912407E-06_Rkind) * tab_Pl(23) +    &
              (-0.812400205583611045E-06_Rkind) * tab_Pl(24) +    &
              (-0.128163223845606709E-05_Rkind) * tab_Pl(25) +    &
              (-0.646090897300889352E-06_Rkind) * tab_Pl(26) +    &
              (0.327319857948809944E-06_Rkind) * tab_Pl(27) +    &
              (0.812961909547953328E-06_Rkind) * tab_Pl(28) +    &
              (0.415922247152272469E-06_Rkind) * tab_Pl(29)

    ! hess: R_H-CN,R_H-CN
    Func(4) = (0.490102788711937298_Rkind) * tab_Pl(0) +    &
              (-0.607405214287308134E-01_Rkind) * tab_Pl(1) +    &
              (0.471897394697577272E-01_Rkind) * tab_Pl(2) +    &
              (0.215199509086547482E-01_Rkind) * tab_Pl(3) +    &
              (0.187589053704660357E-01_Rkind) * tab_Pl(4) +    &
              (-0.730896227164504544E-02_Rkind) * tab_Pl(5) +    &
              (-0.495680012901110034E-02_Rkind) * tab_Pl(6) +    &
              (0.545601888552436865E-03_Rkind) * tab_Pl(7) +    &
              (0.288205091824102807E-02_Rkind) * tab_Pl(8) +    &
              (-0.150455202101582481E-03_Rkind) * tab_Pl(9) +    &
              (-0.135822812297228617E-02_Rkind) * tab_Pl(10) +    &
              (0.378620952994913158E-03_Rkind) * tab_Pl(11) +    &
              (0.223172238114486974E-03_Rkind) * tab_Pl(12) +    &
              (-0.113639874719230033E-03_Rkind) * tab_Pl(13) +    &
              (0.479670914337290878E-04_Rkind) * tab_Pl(14) +    &
              (-0.295643176360794466E-04_Rkind) * tab_Pl(15) +    &
              (-0.220891025839080444E-04_Rkind) * tab_Pl(16) +    &
              (0.233380028833382715E-04_Rkind) * tab_Pl(17) +    &
              (0.226740013135044743E-05_Rkind) * tab_Pl(18) +    &
              (-0.383700217098819257E-05_Rkind) * tab_Pl(19) +    &
              (-0.176925895920402338E-05_Rkind) * tab_Pl(20) +    &
              (0.126527708534225592E-06_Rkind) * tab_Pl(21) +    &
              (0.102875267686912131E-05_Rkind) * tab_Pl(22) +    &
              (-0.244257932870602493E-06_Rkind) * tab_Pl(23) +    &
              (-0.132903576796384672E-06_Rkind) * tab_Pl(24) +    &
              (0.318903079193596794E-07_Rkind) * tab_Pl(25) +    &
              (0.270895958135604389E-06_Rkind) * tab_Pl(26) +    &
              (-0.334319014937836363E-07_Rkind) * tab_Pl(27) +    &
              (-0.215243670177822377E-06_Rkind) * tab_Pl(28) +    &
              (-0.532921517866474822E-07_Rkind) * tab_Pl(29)

    ! hess: R_H-CN,R_CN
    Func(5) = (-0.269621628391861212_Rkind) * tab_Pl(0) +    &
              (0.298216690078651966E-01_Rkind) * tab_Pl(1) +    &
              (-0.252263406635114523E-02_Rkind) * tab_Pl(2) +    &
              (-0.189365831326591162E-01_Rkind) * tab_Pl(3) +    &
              (-0.135862653973030569E-01_Rkind) * tab_Pl(4) +    &
              (0.180882447027921715E-02_Rkind) * tab_Pl(5) +    &
              (-0.807456875985611278E-02_Rkind) * tab_Pl(6) +    &
              (0.162126376060913813E-02_Rkind) * tab_Pl(7) +    &
              (0.246532720643025506E-02_Rkind) * tab_Pl(8) +    &
              (-0.529028555072996428E-03_Rkind) * tab_Pl(9) +    &
              (-0.604375862043439872E-03_Rkind) * tab_Pl(10) +    &
              (0.450771375745010218E-04_Rkind) * tab_Pl(11) +    &
              (0.158382700747678866E-03_Rkind) * tab_Pl(12) +    &
              (0.838098940644294889E-05_Rkind) * tab_Pl(13) +    &
              (0.157629766319828802E-05_Rkind) * tab_Pl(14) +    &
              (-0.182030790249838681E-04_Rkind) * tab_Pl(15) +    &
              (-0.271668080516251952E-04_Rkind) * tab_Pl(16) +    &
              (0.156823664488067170E-04_Rkind) * tab_Pl(17) +    &
              (0.129700582797636253E-04_Rkind) * tab_Pl(18) +    &
              (-0.588984022955290420E-05_Rkind) * tab_Pl(19) +    &
              (-0.329310498897404528E-05_Rkind) * tab_Pl(20) +    &
              (0.382968507344196263E-06_Rkind) * tab_Pl(21) +    &
              (0.662203964149338334E-06_Rkind) * tab_Pl(22) +    &
              (0.537619208727332427E-06_Rkind) * tab_Pl(23) +    &
              (0.503885020561694568E-07_Rkind) * tab_Pl(24) +    &
              (-0.117926135430141456E-06_Rkind) * tab_Pl(25) +    &
              (-0.950694993110524940E-07_Rkind) * tab_Pl(26) +    &
              (0.661655056484692365E-07_Rkind) * tab_Pl(27) +    &
              (-0.444581410057473655E-07_Rkind) * tab_Pl(28) +    &
              (-0.327143067519841403E-07_Rkind) * tab_Pl(29)

    Func(6) =   Func(5)      ! hess: R_CN,R_H-CN

    ! hess: R_CN,R_CN
    Func(7) = (1.89775051256421889_Rkind) * tab_Pl(0) +    &
              (0.189212178373149228E-01_Rkind) * tab_Pl(1) +    &
              (-0.177346711887173797_Rkind) * tab_Pl(2) +    &
              (0.285390255705121020E-01_Rkind) * tab_Pl(3) +    &
              (0.117421295103266571_Rkind) * tab_Pl(4) +    &
              (-0.108113416642430672E-01_Rkind) * tab_Pl(5) +    &
              (-0.181273718435335900E-01_Rkind) * tab_Pl(6) +    &
              (0.219170948917988715E-02_Rkind) * tab_Pl(7) +    &
              (0.254074210398737418E-02_Rkind) * tab_Pl(8) +    &
              (0.579675850294103122E-03_Rkind) * tab_Pl(9) +    &
              (0.285040611946511661E-03_Rkind) * tab_Pl(10) +    &
              (-0.573959392146183443E-03_Rkind) * tab_Pl(11) +    &
              (-0.437665242169073482E-03_Rkind) * tab_Pl(12) +    &
              (0.220390748629632348E-03_Rkind) * tab_Pl(13) +    &
              (0.218654827462475313E-03_Rkind) * tab_Pl(14) +    &
              (-0.476404921173668880E-04_Rkind) * tab_Pl(15) +    &
              (-0.682533881162194865E-04_Rkind) * tab_Pl(16) +    &
              (-0.472076754440356353E-05_Rkind) * tab_Pl(17) +    &
              (0.917501271955926315E-05_Rkind) * tab_Pl(18) +    &
              (0.142547973671411078E-04_Rkind) * tab_Pl(19) +    &
              (0.801385332200406864E-05_Rkind) * tab_Pl(20) +    &
              (-0.431806540655166862E-05_Rkind) * tab_Pl(21) +    &
              (-0.549017839362386120E-05_Rkind) * tab_Pl(22) +    &
              (-0.380986883790360773E-06_Rkind) * tab_Pl(23) +    &
              (0.409272286953666955E-05_Rkind) * tab_Pl(24) +    &
              (0.514630181802097401E-05_Rkind) * tab_Pl(25) +    &
              (0.268744213230514938E-05_Rkind) * tab_Pl(26) +    &
              (-0.133719501498734565E-05_Rkind) * tab_Pl(27) +    &
              (-0.354094525621361928E-05_Rkind) * tab_Pl(28) +    &
              (-0.191216423100044473E-05_Rkind) * tab_Pl(29)

    DO i=0,max_deg
      CALL dealloc_dnS(tab_Pl(i))
    END DO
    CALL dealloc_dnS(dnCth)

  END SUBROUTINE EvalFunc_QML_CNH_Murrell

END MODULE QML_CNH_Murrell_m
