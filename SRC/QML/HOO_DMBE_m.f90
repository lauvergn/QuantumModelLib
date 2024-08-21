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
!> @brief Module which makes the initialization, calculation of the HOO_DMBE potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_HOO_DMBE_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the HOO_DMBE parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_HOO_DMBE_t

   PRIVATE

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_HOO_DMBE
    PROCEDURE :: Write_QModel     => Write_QML_HOO_DMBE
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_QML_HOO_DMBE
  END TYPE QML_HOO_DMBE_t

  PUBLIC :: QML_HOO_DMBE_t,Init_QML_HOO_DMBE


    REAL (kind=Rkind), parameter :: C(52) = [                                                &
     .49040645E+01_Rkind, -.86748216E+01_Rkind,  .50555792E+01_Rkind,  .42941301E+01_Rkind,  &
    -.41874792E+01_Rkind,  .13461379_Rkind,     -.99064922_Rkind,      .13358488E+01_Rkind,  &
     .13495231E+01_Rkind, -.18529696_Rkind,     -.23534213E+02_Rkind,  .24289930E+02_Rkind,  &
    -.50209026E+01_Rkind, -.10365484E+02_Rkind,  .46692224E+01_Rkind, -.14747138E+01_Rkind,  &
     .23119718E+01_Rkind, -.18247842E+01_Rkind, -.28472166_Rkind,      .51036509_Rkind,      &
     .19124083_Rkind,      .45405729E+01_Rkind,  .11087611_Rkind,     -.19990481_Rkind,      &
    -.37356178_Rkind,      .46142042E-01_Rkind, -.20565580_Rkind,     -.27015963_Rkind,      &
     .34085281_Rkind,      .28321162_Rkind,     -.11558481_Rkind,     -.29448886_Rkind,      &
    -.52932488_Rkind,      .58159523E-01_Rkind, -.48649560E-02_Rkind,  .11949167E-01_Rkind,  &
     .21409804E-01_Rkind, -.20620608E-02_Rkind,  .30177088E-01_Rkind,  .27880291E-01_Rkind,  &
     .88458711E-02_Rkind,  .13137410E-01_Rkind, -.24705619E-01_Rkind,  -.31085889E-01_Rkind, &
     .34317857E-02_Rkind,  .52593878E-01_Rkind,  .79500714E-01_Rkind,  -.79782216E-02_Rkind, &
     .31164575E-01_Rkind, -.28737598E-01_Rkind,  .98201698_Rkind,       .62000000_Rkind]

    REAL (kind=Rkind), parameter  :: RMOO   = 2.2818_Rkind
    REAL (kind=Rkind), parameter  :: R0OO   = 5.661693_Rkind
    REAL (kind=Rkind), parameter  :: RMOH   = 1.8344_Rkind
    REAL (kind=Rkind), parameter  :: R0OH   = 6.294894_Rkind

    REAL (kind=Rkind), parameter  :: COO(10)= [ZERO,ZERO,ZERO,ZERO,ZERO,        &
                        15.40_Rkind,ZERO,235.219943_Rkind,ZERO,4066.23929_Rkind]
    REAL (kind=Rkind), parameter  :: COH(10)= [ZERO,ZERO,ZERO,ZERO,ZERO,TEN,    &
                                    ZERO,180.447673_Rkind,ZERO,3685.25842_Rkind]

    REAL (kind=Rkind), parameter  :: C4   = -0.92921_Rkind
    REAL (kind=Rkind), parameter  :: C5   = -1.79000_Rkind

    REAL (kind=Rkind), parameter  :: RK0OO(10)= [ZERO,ZERO,ZERO,ZERO,ZERO,      &
                  -.27847758_Rkind,ZERO,-.46815641_Rkind,ZERO,-1.20506384_Rkind]
    REAL (kind=Rkind), parameter  :: RK1OO(10)= [ZERO,ZERO,ZERO,                &
                   3.35224980_Rkind,3.35224980_Rkind,0.95273753_Rkind,ZERO,     &
                   0.94148408_Rkind,ZERO,0.72379129_Rkind]
    REAL (kind=Rkind), parameter  :: RK0OH(10)= [ZERO,ZERO,ZERO,ZERO,ZERO,      &
                    0.02465005_Rkind,ZERO,0.05036950_Rkind,ZERO,0.06294371_Rkind]
    REAL (kind=Rkind), parameter  :: RK1OH(10)= [ZERO,ZERO,ZERO,                &
                    2.54532760_Rkind,2.54532760_Rkind,0.68758097_Rkind,ZERO,    &
                    0.82542359_Rkind,ZERO,0.94034225_Rkind]

    REAL (kind=Rkind), parameter  :: ADAMP(10) = [ZERO,ZERO,ZERO,               &
                    5.0079875_Rkind,3.8428294_Rkind,3.0951333_Rkind,ZERO,       &
                    2.1999000_Rkind,ZERO,1.6880714_Rkind]
    REAL (kind=Rkind), parameter  :: BDAMP(10) = [ZERO,ZERO,ZERO,               &
                             10.6645006_Rkind,9.6758155_Rkind,8.7787895_Rkind,  &
                             ZERO,7.2265123_Rkind,ZERO,5.9487108_Rkind]

    REAL (kind=Rkind), parameter :: R10 = 2.5143000_Rkind
    REAL (kind=Rkind), parameter :: R20 = 2.6469057_Rkind
    REAL (kind=Rkind), parameter :: R30 = 2.6469057_Rkind

!   BLOCK DATA HO2DAT_HOO_DMBE4
!       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!       COMMON/COEFF_HOO_DMBE4/C(52)
!       COMMON/DISPC_HOO_DMBE4/COO(10),COH(10)
!       COMMON/DIATDI_HOO_DMBE4/R0OO,RMOO,R0OH,RMOH
!       COMMON/RKVAL_HOO_DMBE4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
!       COMMON/POLAR_HOO_DMBE4/C4,C5
!       COMMON/DAMPC_HOO_DMBE4/ADAMP(10),BDAMP(10)
!       COMMON/REFGEO_HOO_DMBE4/R10,R20,R30
! !     ***************************************************************
!       DATA C/ &
!        .49040645E+01_Rkind, -.86748216E+01_Rkind,  .50555792E+01_Rkind,  .42941301E+01_Rkind,&
!       -.41874792E+01_Rkind,  .13461379_Rkind,     -.99064922_Rkind,      .13358488E+01_Rkind,&
!        .13495231E+01_Rkind, -.18529696_Rkind,     -.23534213E+02_Rkind,  .24289930E+02_Rkind,&
!       -.50209026E+01_Rkind, -.10365484E+02_Rkind,  .46692224E+01_Rkind, -.14747138E+01_Rkind,&
!        .23119718E+01_Rkind, -.18247842E+01_Rkind, -.28472166_Rkind,      .51036509_Rkind,&
!        .19124083_Rkind,      .45405729E+01_Rkind,  .11087611_Rkind,      -.19990481_Rkind,&
!       -.37356178_Rkind,      .46142042E-01_Rkind, -.20565580_Rkind,      -.27015963_Rkind,&
!        .34085281_Rkind,      .28321162_Rkind,     -.11558481_Rkind,      -.29448886_Rkind,&
!       -.52932488_Rkind,      .58159523E-01_Rkind, -.48649560E-02_Rkind,   .11949167E-01_Rkind,&
!        .21409804E-01_Rkind, -.20620608E-02_Rkind,  .30177088E-01_Rkind,   .27880291E-01_Rkind,&
!        .88458711E-02_Rkind,  .13137410E-01_Rkind, -.24705619E-01_Rkind,   -.31085889E-01_Rkind,&
!        .34317857E-02_Rkind,  .52593878E-01_Rkind,  .79500714E-01_Rkind,    -.79782216E-02_Rkind,&
!        .31164575E-01_Rkind, -.28737598E-01_Rkind,  .98201698_Rkind,         .62000000_Rkind/
!     DATA R0OO,RMOO,R0OH,RMOH/5.661693_Rkind,2.2818_Rkind,6.294894_Rkind,1.8344_Rkind/
!     DATA COO/ZERO,ZERO,ZERO,ZERO,0._Rkind,15.40_Rkind,ZERO,235.219943_Rkind,&
!              ZERO,4066.23929_Rkind/
!     DATA COH/ZERO,ZERO,ZERO,ZERO,0._Rkind,10.00_Rkind,ZERO,180.447673_Rkind,&
!               ZERO,3685.25842_Rkind/
!     DATA C4,C5/-0.92921_Rkind,-1.79000_Rkind/
!     DATA RK0OO/ZERO,ZERO,ZERO,ZERO,ZERO,-.27847758_Rkind,ZERO,&
!                 -.46815641_Rkind,ZERO,-1.20506384_Rkind/
!     DATA RK1OO/ZERO,ZERO,ZERO,3.35224980_Rkind,3.35224980_Rkind,&
!                 0.95273753_Rkind,ZERO,0.94148408_Rkind,ZERO,0.72379129_Rkind/
!     DATA RK0OH/ZERO,ZERO,ZERO,ZERO,ZERO,0.02465005_Rkind,ZERO,&
!                 0.05036950_Rkind,ZERO,0.06294371_Rkind/
!     DATA RK1OH/ZERO,ZERO,ZERO,2.54532760_Rkind,2.54532760_Rkind,&
!                 0.68758097_Rkind,ZERO,0.82542359_Rkind,ZERO,0.94034225_Rkind/
!     DATA ADAMP/ZERO,ZERO,ZERO,5.0079875_Rkind,3.8428294_Rkind,3.0951333_Rkind,&
!                 ZERO,2.1999000_Rkind,ZERO,1.6880714_Rkind/
!     DATA BDAMP/ZERO,ZERO,0._Rkind,10.6645006_Rkind,9.6758155_Rkind,8.7787895_Rkind,&
!                 ZERO,7.2265123_Rkind,ZERO,5.9487108_Rkind/
!     DATA R10,R20,R30/2.5143000_Rkind,2.6469057_Rkind,2.6469057_Rkind/
!   END

  CONTAINS
!> @brief Function which makes the initialization of the HOO_DMBE parameters.
!!
!! @param QModel             TYPE(QML_HOO_DMBE_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_HOO_DMBE(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_HOO_DMBE_t)                           :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_HOO_DMBE'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf      = 1
    QModel%ndimQ      = 3
    QModel%ndimCart   = 9

    IF (QModel%Cart_TO_Q) THEN
      QModel%ndim       = QModel%ndimCart
    ELSE
      QModel%ndim       = QModel%ndimQ
    END IF

    QModel%pot_name   = 'hoo_dmbe'
    QModel%no_ana_der = .TRUE.


    IF (debug) write(out_unit,*) 'init Q0 of HOO_DMBE (HOO minimum)'
    QModel%Q0 = [2.806_Rkind,2.271_Rkind,2.271_Rkind]

    IF (debug) write(out_unit,*) 'init d0GGdef of HOO_DMBE'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)

    IF (debug) THEN
      CALL Write_QML_HOO_DMBE(QModel,nio=out_unit)
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_HOO_DMBE
!> @brief Subroutine wich prints the current QML_HOO_DMBE parameters.
!!
!! @param QModel            CLASS(QML_HOO_DMBE_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_HOO_DMBE(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HOO_DMBE_t),   intent(in) :: QModel
    integer,                   intent(in) :: nio

    write(nio,*) 'HOO_DMBE IV default parameters'

    write(nio,*) 'Ref: '
    write(nio,*) '  M. R. Pastrana, L. A. M. Quintales, J. Brandão and A. J. C. Varandas'
    write(nio,*) '  JCP, 1990, 94, 8073-8080, doi: 10.1021/j100384a019'
    write(nio,*)
    write(nio,*) ' Values from Table VII of the JCP paper:'
    write(nio,*)
    write(nio,*) '  OH...O '
    write(nio,*) '    R1,R2,R3= [5.663,3.821,1.842] bohr'
    write(nio,*) '    E = -0.1738 Hartree'
    write(nio,*) '    E = -0.1737745714 Hartree (QML)'
    write(nio,*)
    write(nio,*) '  H...O-O '
    write(nio,*) '    R1,R2,R3= [2.282,7.547,9.829] bohr'
    write(nio,*) '    E = -0.1916 Hartree'
    write(nio,*) '    E = -0.1916168740 Hartree (QML)'
    write(nio,*)
    write(nio,*) '  HO2 '
    write(nio,*) '    R1,R2,R3= [2.806,2.271,2.271] bohr'
    write(nio,*) '    E = -0.2141 Hartree'
    write(nio,*) '    E = -0.2140801239 Hartree (QML)'

    write(nio,*)
    write(nio,*) 'end HOO_DMBE IV default parameters'

    write(nio,*) 'HOO_DMBE IV current parameters'

    CALL QModel%QML_Empty_t%Write_QModel(nio)

    write(nio,*) 'end HOO_DMBE IV current parameters'

  END SUBROUTINE Write_QML_HOO_DMBE
!> @brief Subroutine wich calculates the HOO_DMBE potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_HOO_DMBE_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         Potential with derivatives,.
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to secify the derivative order:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_HOO_DMBE(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HOO_DMBE_t), intent(in)    :: QModel
    TYPE (dnS_t),            intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),            intent(in)    :: dnQ(:)
    integer,                 intent(in)    :: nderiv

    real(kind=Rkind) :: V,Q(3)

    Q(:) = get_d0(dnQ)

    CALL HOO_DMBE4_pes(Q,V)

    CALL set_dnS(Mat_OF_PotDia(1,1),d0=V)

  END SUBROUTINE EvalPot_QML_HOO_DMBE

  ! here we suppose that the atom ordering: H1-O2-O3
  SUBROUTINE Cart_TO_Q_QML_HOO_DMBE(QModel,dnX,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HOO_DMBE_t), intent(in)    :: QModel
    TYPE (dnS_t),            intent(in)    :: dnX(:,:)
    TYPE (dnS_t),            intent(inout) :: dnQ(:)
    integer,                 intent(in)    :: nderiv

    ! local vector
    integer         :: i,j
    TYPE (dnS_t)    :: VecOO(3),VecHO2(3),VecHO3(3)


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Cart_TO_Q_QML_HOO_DMBE'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'dnX'
      DO i=1,size(dnX,dim=2)
      DO j=1,size(dnX,dim=1)
        CALL Write_dnS(dnX(j,i),out_unit)
      END DO
      END DO
      flush(out_unit)
    END IF

    VecOO(:)  = dnX(:,3)-dnX(:,2)
    VecHO2(:) = dnX(:,2)-dnX(:,1)
    VecHO3(:) = dnX(:,3)-dnX(:,1)
    IF (debug) write(out_unit,*) 'Cart_TO_Q_QML_HOO_DMBE vect done'

    dnQ(1) = sqrt(dot_product(VecOO,VecOO))
    dnQ(2) = sqrt(dot_product(VecHO2,VecHO2))
    dnQ(3) = sqrt(dot_product(VecHO3,VecHO3))

    CALL dealloc_dnS(VecOO)
    CALL dealloc_dnS(VecHO2)
    CALL dealloc_dnS(VecHO3)

    IF (debug) THEN
      CALL Write_dnS(dnQ(1),out_unit,info='dnQ(1)')
      CALL Write_dnS(dnQ(2),out_unit,info='dnQ(2)')
      CALL Write_dnS(dnQ(3),out_unit,info='dnQ(3)')
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
  END SUBROUTINE Cart_TO_Q_QML_HOO_DMBE


  SUBROUTINE HOO_DMBE4_pes(X,F)
  ! This is the DMBE IV potential energy surface for H + O2
  IMPLICIT NONE

    REAL (kind=Rkind), intent(in)    :: X(3)
    REAL (kind=Rkind), intent(inout) :: F

    REAL (kind=Rkind) :: R1,R2,R3,Q1,Q2,Q3

    R1=X(1) ! O-O distance
    R2=X(2) ! H-O(2) distance
    R3=X(3) ! H-O(3) distance

    Q1=ONE/SQRT(THREE)*(R1+R2+R3)
    Q2=ONE/SQRT(TWO)*(R2-R3)
    Q3=ONE/SQRT(SIX)*(TWO*R1-R2-R3)

    F = VOO_HOO_DMBE4(R1) + VOH_HOO_DMBE4(R2) + VOH_HOO_DMBE4(R3) +             &
        THREBQ_HOO_DMBE4(Q1,Q2,Q3,R1,R2,R3) +                                   &
        EXDIS_HOO_DMBE4(R1,R2,R3)  + ELECT_HOO_DMBE4(R1,R2,R3)

  END SUBROUTINE HOO_DMBE4_pes

  FUNCTION THREBQ_HOO_DMBE4(Q1,Q2,Q3,R1,R2,R3)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: THREBQ_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: Q1,Q2,Q3,R1,R2,R3


    REAL (kind=Rkind) :: POLQ,DECAY1,DECAY2,DECAY3
    REAL (kind=Rkind) :: Q12,Q13,Q14,Q15,Q16,Q22,Q32
    REAL (kind=Rkind) :: TQ1,TQ2,TQ3,TQ12,TQ13,TQ22
    REAL (kind=Rkind) :: S1,S2,S3


    Q12=Q1*Q1
    Q13=Q12*Q1
    Q14=Q13*Q1
    Q15=Q14*Q1
    Q16=Q15*Q1
    Q22=Q2*Q2
    Q32=Q3*Q3
    TQ1=Q22+Q32
    TQ2=Q32-THREE*Q22
    TQ3=Q22-Q32
    TQ12=TQ1*TQ1
    TQ13=TQ12*TQ1
    TQ22=TQ2*TQ2
    S1=R1-R10
    S2=R2-R20
    S3=R3-R30

    POLQ=C(1)*Q1+C(2)*Q12+C(3)*TQ1+C(4)*Q13+C(5)*Q1*TQ1+&
      C(6)*Q3*TQ2+C(7)*Q14+C(8)*Q12*TQ1+C(9)*TQ1**2+C(10)*Q1*Q3*TQ2+&
      C(11)*Q3+C(12)*Q1*Q3+C(13)*TQ3+C(14)*Q12*Q3+C(15)*Q1*TQ3+&
      C(16)*Q3*TQ1+C(17)*Q13*Q3+C(18)*Q12*TQ3+C(19)*Q1*Q3*TQ1+&
      C(20)*Q32*TQ2+C(21)*TQ1*TQ3+C(22)+C(23)*Q15+C(24)*Q13*TQ1+&
      C(25)*Q1*TQ12+C(26)*Q12*Q3*TQ2+C(27)*Q3*TQ1*TQ2+C(28)*Q14*Q3+&
      C(29)*Q13*TQ3+C(30)*Q12*Q3*TQ1+C(31)*Q1*Q32*TQ2+C(32)*Q1*TQ1*TQ3+&
      C(33)*Q3*TQ12+C(34)*Q3*TQ2*TQ3+C(35)*Q16+C(36)*Q14*TQ1+&
      C(37)*Q12*TQ12+C(38)*Q13*Q3*TQ2+C(39)*Q1*Q3*TQ1*TQ2+C(40)*TQ13+&
      C(41)*Q32*TQ22+C(42)*Q15*Q3+C(43)*Q14*TQ3+C(44)*Q13*Q3*TQ1+&
      C(45)*Q12*Q32*TQ2+C(46)*Q12*TQ1*TQ3+C(47)*Q1*Q3*TQ12+&
      C(48)*Q1*Q3*TQ2*TQ3+C(49)*Q32*TQ1*TQ2+C(50)*TQ12*TQ3

    DECAY1=ONE-TANH(C(51)*S1)
    DECAY2=ONE-TANH(C(52)*S2)
    DECAY3=ONE-TANH(C(52)*S3)

    THREBQ_HOO_DMBE4 = POLQ*DECAY1*DECAY2*DECAY3

  END FUNCTION THREBQ_HOO_DMBE4
  FUNCTION VOH_HOO_DMBE4(R)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: VOH_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R

    VOH_HOO_DMBE4 = EHFOH_HOO_DMBE4(R) + DISOH_HOO_DMBE4(R)

  END FUNCTION VOH_HOO_DMBE4

  FUNCTION EHFOH_HOO_DMBE4(R)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: EHFOH_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R


    REAL (kind=Rkind)             :: X,R2,R3
    REAL (kind=Rkind), parameter  :: D      = 0.13825385_Rkind
    REAL (kind=Rkind), parameter  :: ASV(4) = [2.6564788_Rkind,1.7450528_Rkind, &
                                               0.71014391_Rkind,2.5453276_Rkind]

    X=R-RMOH
    R2=X*X
    R3=R2*X

    EHFOH_HOO_DMBE4 = -D*(ONE+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*EXP(-ASV(4)*X)

  END FUNCTION EHFOH_HOO_DMBE4

  FUNCTION DISOH_HOO_DMBE4(R)
  IMPLICIT NONE

    REAL (kind=Rkind)             :: DISOH_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R

    DISOH_HOO_DMBE4 = DISP_HOO_DMBE4(R,COH(6),COH(8),COH(10),R0OH,RMOH)

  END FUNCTION DISOH_HOO_DMBE4

  FUNCTION VOO_HOO_DMBE4(R)
  IMPLICIT NONE

    REAL (kind=Rkind)             :: VOO_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R

    !REAL (kind=Rkind)             :: EHFOO_HOO_DMBE4,DISOO_HOO_DMBE4 ! functions

    VOO_HOO_DMBE4 = EHFOO_HOO_DMBE4(R)+DISOO_HOO_DMBE4(R)

  END FUNCTION VOO_HOO_DMBE4

  FUNCTION EHFOO_HOO_DMBE4(R)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: EHFOO_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R


    REAL (kind=Rkind)             :: X,R2,R3
    REAL (kind=Rkind), parameter  :: D      = 0.14291202_Rkind
    REAL (kind=Rkind), parameter  :: ASV(4) = [3.6445906_Rkind,3.9281238_Rkind, &
                                               2.0986689_Rkind,3.3522498_Rkind]

    X  = R-RMOO
    R2 = X*X
    R3 = R2*X

    EHFOO_HOO_DMBE4 = -D*(ONE+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*EXP(-ASV(4)*X)

  END FUNCTION EHFOO_HOO_DMBE4

  FUNCTION DISOO_HOO_DMBE4(R)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: DISOO_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R


    DISOO_HOO_DMBE4 = DISP_HOO_DMBE4(R,COO(6),COO(8),COO(10),R0OO,RMOO)

  END FUNCTION DISOO_HOO_DMBE4

  FUNCTION CEF_HOO_DMBE4(CAS,RK01,RK11,RK02,RK12,RE1,RE2,R1,R2)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: CEF_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: CAS
    REAL (kind=Rkind), intent(in) :: RK01,RK11,RK02,RK12
    REAL (kind=Rkind), intent(in) :: RE1,RE2
    REAL (kind=Rkind), intent(in) :: R1,R2

    CEF_HOO_DMBE4 = HALF*CAS*((ONE-RK01*EXP(-RK11*(R1-RE1)))*TANH(RK12*R2)+&
                    (ONE-RK02*EXP(-RK12*(R2-RE2)))*TANH(RK11*R1))

  END FUNCTION CEF_HOO_DMBE4

  FUNCTION EXDIS_HOO_DMBE4(R1,R2,R3)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: EXDIS_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R1,R2,R3


    REAL (kind=Rkind)  :: CEFOO(10),CEFOH2(10),CEFOH3(10),CEDOO(10),CEDOH2(10),CEDOH3(10)
    integer            :: i

    DO i=6,10,2
      CEFOO(i)=CEF_HOO_DMBE4(COO(i),RK0OH(i),RK1OH(i),RK0OH(i),RK1OH(i),RMOH,RMOH,R2,R3)
      CEDOO(i)=CEFOO(i)-COO(i)
      CEFOH2(i)=CEF_HOO_DMBE4(COH(i),RK0OO(i),RK1OO(i),RK0OH(i),RK1OH(i),RMOO,RMOH,R1,R3)
      CEDOH2(i)=CEFOH2(i)-COH(i)
      CEFOH3(i)=CEF_HOO_DMBE4(COH(i),RK0OO(i),RK1OO(i),RK0OH(i),RK1OH(i),RMOO,RMOH,R1,R2)
      CEDOH3(i)=CEFOH3(i)-COH(i)
    END DO

    EXDIS_HOO_DMBE4 = DISP_HOO_DMBE4(R1,CEDOO(6),CEDOO(8),CEDOO(10),R0OO,RMOO)    + &
                      DISP_HOO_DMBE4(R2,CEDOH2(6),CEDOH2(8),CEDOH2(10),R0OH,RMOH) + &
                      DISP_HOO_DMBE4(R3,CEDOH3(6),CEDOH3(8),CEDOH3(10),R0OH,RMOH)

  END FUNCTION EXDIS_HOO_DMBE4

  FUNCTION ELECT_HOO_DMBE4(R1,R2,R3)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: ELECT_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R1,R2,R3

    REAL (kind=Rkind) :: C4OHR2,C5OHR2,C4OHR3,C5OHR3,C4OO,C5OO,TERM4,TERM5
    REAL (kind=Rkind) :: C42,C43,C52,C53,R23,R33,R34,R14,R15,R25,R35,R24
    REAL (kind=Rkind) :: RMQ,RMQ5,RMR3,RMR2,RMR33,RMR23
    REAL (kind=Rkind) :: TAO,TAH2,TAH3,EX3,EX2,R3E3,R2E2,CRE43,CRE42,CRE53,CRE52
    REAL (kind=Rkind) :: RROH2,RROH3,RROO

    C42=C4
    C43=C4
    C52=C5
    C53=C5
    R23=R2**3
    R24=R23*R2
    R33=R3**3
    R34=R33*R3
    R14=R1**4
    R15=R14*R1
    R25=R24*R2
    R35=R34*R3
    RMQ=RMOH**4
    RMQ5=0.50_Rkind/RMQ
    RMR3=RMQ5*R34
    RMR2=RMQ5*R24
    RMR33=RMQ5*R33
    RMR23=RMQ5*R23
    TAO=TANH(RK1OO(4)*R1)
    TAH2=TANH(RK1OH(4)*R2)
    TAH3=TANH(RK1OH(4)*R3)
    EX3=EXP(-RK1OH(4)*(R3-RMOH))
    EX2=EXP(-RK1OH(4)*(R2-RMOH))
    R3E3=RMR3*EX3
    R2E2=RMR2*EX2
    CRE43=C4*R3E3
    CRE42=C4*R2E2
    CRE53=C5*R3E3
    CRE52=C5*R2E2
    C4OHR2=CRE43*TAO
    C5OHR2=CRE53*TAO
    C4OHR3=CRE42*TAO
    C5OHR3=CRE52*TAO
    C4OO=CRE43*TAH2+CRE42*TAH3
    C5OO=CRE53*TAH2+CRE52*TAH3
    RROH2=TWO*R2/(RMOH+2.5_Rkind*R0OH)
    RROH3=TWO*R3/(RMOH+2.5_Rkind*R0OH)
    RROO=TWO*R1/(RMOO+2.5_Rkind*R0OO)

    TERM4 = C4OO/R14*(ONE-EXP(-ADAMP(4)*RROO-BDAMP(4)*RROO**2))**4+             &
            C4OHR2/R24*(ONE-EXP(-ADAMP(4)*RROH2-BDAMP(4)*RROH2**2))**4+         &
            C4OHR3/R34*(ONE-EXP(-ADAMP(4)*RROH3-BDAMP(4)*RROH3**2))**4

    TERM5 = C5OO/R15*(ONE-EXP(-ADAMP(5)*RROO-BDAMP(5)*RROO**2))**5+             &
            C5OHR2/R25*(ONE-EXP(-ADAMP(5)*RROH2-BDAMP(5)*RROH2**2))**5+         &
            C5OHR3/R35*(ONE-EXP(-ADAMP(5)*RROH3-BDAMP(5)*RROH3**2))**5

    ELECT_HOO_DMBE4 = TERM4 + TERM5

  END FUNCTION ELECT_HOO_DMBE4
  FUNCTION DISP_HOO_DMBE4(R,C6,C8,C10,R0,RM)
    IMPLICIT NONE

    REAL (kind=Rkind)             :: DISP_HOO_DMBE4
    REAL (kind=Rkind), intent(in) :: R,C6,C8,C10,R0,RM

    REAL (kind=Rkind) :: RR,R6,R8,R10,D6,D8,D10


      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=TWO*R/(RM+2.5_Rkind*R0)
      D6=(ONE-EXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6
      D8=(ONE-EXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8
      D10=(ONE-EXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10

      DISP_HOO_DMBE4 = -C6/R6*D6-C8/R8*D8-C10/R10*D10

  END FUNCTION DISP_HOO_DMBE4
END MODULE QML_HOO_DMBE_m
