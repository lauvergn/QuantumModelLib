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

!> @brief Module which makes the initialization, calculation of the HCN_Murrell potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_HCN_Murrell_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the HCN_Murrell parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_HCN_Murrell_t

   PRIVATE

   CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_HCN_Murrell
    PROCEDURE :: Write_QModel     => Write_QML_HCN_Murrell
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_QML_HCN_Murrell
    PROCEDURE :: Eval_QModel_Func => EvalFunc_QML_HCN_Murrell
  END TYPE QML_HCN_Murrell_t

  PUBLIC :: QML_HCN_Murrell_t,Init_QML_HCN_Murrell


  CONTAINS
!> @brief Function which makes the initialization of the HCN_Murrell parameters.
!!
!! @param QModel             TYPE(QML_HCN_Murrell_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_HCN_Murrell(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_HCN_Murrell_t)                     :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_HCN_Murrell'
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

    QModel%ndimQ      = 3

    IF (QModel%Cart_TO_Q)    QModel%option = 0
    IF (QModel%option == -1) QModel%option = 0

    SELECT CASE(QModel%option)
    CASE (0,1) !

      IF (QModel%Cart_TO_Q) THEN
        QModel%ndim       = QModel%ndimCart
      ELSE
        QModel%ndim       = QModel%ndimQ
      END IF

      ! so that, we can get the irc function with the 3D model
      QModel%ndimFunc = 1
      QModel%nb_Func  = -1 ! V, R1op,R2opt,R3opt,RPH_r,RPH_th

      QModel%pot_name   = 'HCN_Murrell'
    CASE (2)
      QModel%ndimQ    = 1
      QModel%ndim     = 1

      QModel%ndimFunc = 1
      QModel%nb_Func  = -1 ! V, R1op,R2opt,R3opt,RPH_r,RPH_th

      QModel%pot_name   = 'HCN_Murrell_IRC'
      QModel%no_ana_der = .FALSE.
    CASE Default
       write(out_unitp,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unitp)
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' option: ',QModel%option
       write(out_unitp,*) ' the possible option are: 0, 1 (3D) or 1 (1D-IRC)'
       STOP 'ERROR in Init_QML_HCN_Murrell: wrong option'
    END SELECT


    IF (debug) write(out_unitp,*) 'init Q0 of HCN_Murrell (HCN_Murrell minimum)'
    QModel%Q0 = [ZERO,3.18722_Rkind,2.17926_Rkind]

    IF (debug) write(out_unitp,*) 'init d0GGdef of HCN_Murrell'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      !CALL Write_QML_HCN_Murrell(QModel,nio=out_unitp)
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_HCN_Murrell
!> @brief Subroutine wich prints the current QML_HCN_Murrell parameters.
!!
!! @param QModel            CLASS(QML_HCN_Murrell_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_HCN_Murrell(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HCN_Murrell_t),  intent(in) :: QModel
    integer,                   intent(in) :: nio

    write(nio,*) 'HCN_Murrell current parameters'

    CALL Write_QML_Empty(QModel%QML_Empty_t,nio)
    write(nio,*) 'Atomic order : C,N,H'
    SELECT CASE(QModel%option)
    CASE (0) ! 3D (three distances)
      write(nio,*) '3D with 3 distances'
    CASE (1) ! 3D Jacobi
      write(nio,*) '3D With Jacobi coordinates'
    CASE (2) ! IRC
      write(nio,*) '1D-IRC + RPH'
    END SELECT

    write(nio,*) 'end HCN_Murrell current parameters'

  END SUBROUTINE Write_QML_HCN_Murrell


!> @brief Subroutine wich calculates the HCN_Murrell potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_HCN_Murrell_t):    derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_HCN_Murrell(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HCN_Murrell_t),  intent(in)    :: QModel
    TYPE (dnS_t),              intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),              intent(in)    :: dnQ(:)
    integer,                   intent(in)    :: nderiv

    TYPE (dnS_t)    :: Func(QModel%nb_Func)


    SELECT CASE(QModel%option)
    CASE (0) ! 3D (three distances)
      CALL QML_EvalPot3D_Murrell(Mat_OF_PotDia,dnQ,.FALSE.)
    CASE (1) ! 3D Jacobi
      CALL QML_EvalPot3D_Murrell(Mat_OF_PotDia,dnQ,.TRUE.)

    CASE (2) ! IRC
      CALL EvalFunc_QML_HCN_Murrell(QModel,Func,dnQ,nderiv)
      Mat_OF_PotDia(1,1) = Func(1)
    END SELECT


  END SUBROUTINE EvalPot_QML_HCN_Murrell

  SUBROUTINE QML_EvalPot3D_Murrell(Mat_OF_PotDia,dnQ,Jacobi)
  USE ADdnSVM_m
  IMPLICIT NONE

    TYPE (dnS_t),            intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),            intent(in)    :: dnQ(:)
    logical,                 intent(in)    :: Jacobi

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

      T=cos(dnQ(1))
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
    !write(out_unitp,*) 'R1,R2,R3',get_d0(R1),get_d0(R2),get_d0(R3)

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
  SUBROUTINE Cart_TO_Q_QML_HCN_Murrell(QModel,dnX,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HCN_Murrell_t), intent(in)    :: QModel
    TYPE (dnS_t),             intent(in)    :: dnX(:,:)
    TYPE (dnS_t),             intent(inout) :: dnQ(:)
    integer,                  intent(in)    :: nderiv

    ! local vector
    integer         :: i,j
    TYPE (dnS_t), allocatable :: Vec12(:),Vec13(:),Vec23(:)

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Cart_TO_Q_QML_HCN_Murrell'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'size(dnQ)',size(dnQ)
      write(out_unitp,*) 'dnQ:'
      DO j=1,size(dnQ,dim=1)
        CALL Write_dnS(dnQ(j),out_unitp,info='dnQ('// int_to_char(j) // ')')
      END DO
      write(out_unitp,*) 'shape dnX',shape(dnX)
      write(out_unitp,*) 'dnX'
      DO i=1,size(dnX,dim=2)
      DO j=1,size(dnX,dim=1)
        CALL Write_dnS(dnX(j,i),out_unitp)
      END DO
      END DO
      flush(out_unitp)
    END IF

    allocate(Vec23(3))
    allocate(Vec12(3))
    allocate(Vec13(3))

    Vec23(:) = dnX(:,3)-dnX(:,2) ! NH
    Vec12(:) = dnX(:,2)-dnX(:,1) ! CN
    Vec13(:) = dnX(:,3)-dnX(:,1) ! CH

    IF (debug) THEN
      write(out_unitp,*) 'Cart_TO_Q_QML_HCN_Murrell vect done'
      flush(out_unitp)
      DO j=1,size(Vec23,dim=1)
        CALL Write_dnS(Vec23(j),out_unitp,info='Vec23')
      END DO
      DO j=1,size(Vec12,dim=1)
        CALL Write_dnS(Vec23(j),out_unitp,info='Vec12')
      END DO
      DO j=1,size(Vec23,dim=1)
        CALL Write_dnS(Vec13(j),out_unitp,info='Vec13')
      END DO
      flush(out_unitp)
    END IF
    dnQ(1) = sqrt(dot_product(Vec13,Vec13)) ! CH
    dnQ(2) = sqrt(dot_product(Vec12,Vec12)) ! CN
    dnQ(3) = sqrt(dot_product(Vec23,Vec23)) ! NH

    CALL dealloc_dnS(Vec23)
    CALL dealloc_dnS(Vec12)
    CALL dealloc_dnS(Vec13)

    IF (debug) THEN
      CALL Write_dnS(dnQ(1),out_unitp,info='dnQ(1)')
      CALL Write_dnS(dnQ(2),out_unitp,info='dnQ(2)')
      CALL Write_dnS(dnQ(3),out_unitp,info='dnQ(3)')
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF
  END SUBROUTINE Cart_TO_Q_QML_HCN_Murrell


  ! for IRC, RPH
  SUBROUTINE EvalFunc_QML_HCN_Murrell(QModel,Func,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HCN_Murrell_t), intent(in)    :: QModel
    TYPE (dnS_t),             intent(inout) :: Func(:)
    TYPE (dnS_t),             intent(in)    :: dnQ(:)
    integer,                  intent(in)    :: nderiv

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
   Func(2) = (0.1357176824058676_Rkind)     * tab_Pl(0)            +             &
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
   Func(3) = (0.1357176824058676_Rkind)     * tab_Pl(0)            +             &
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

  END SUBROUTINE EvalFunc_QML_HCN_Murrell

END MODULE QML_HCN_Murrell_m
