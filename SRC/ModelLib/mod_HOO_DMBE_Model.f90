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

!> @brief Module which makes the initialization, calculation of the HOO_DMBE potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_HOO_DMBE_Model
  USE mod_NumParameters
  USE mod_EmptyModel
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the HOO_DMBE parameters are set-up.
  TYPE, EXTENDS (EmptyModel_t) ::  HOO_DMBE_Model_t

   PRIVATE

   CONTAINS
    PROCEDURE :: Eval_QModel_Pot  => eval_HOO_DMBE_Pot
    PROCEDURE :: Write_QModel     => Write_HOO_DMBE_Model
    PROCEDURE :: Write0_QModel    => Write0_HOO_DMBE_Model
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_HOO_DMBE_Model
  END TYPE HOO_DMBE_Model_t

  PUBLIC :: HOO_DMBE_Model_t,Init_HOO_DMBE_Model

  CONTAINS
!> @brief Function which makes the initialization of the HOO_DMBE parameters.
!!
!! @param QModel             TYPE(HOO_DMBE_Model_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(EmptyModel_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_HOO_DMBE_Model(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (HOO_DMBE_Model_t)                           :: QModel ! RESULT

    TYPE(EmptyModel_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_HOO_DMBE_Model'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_EmptyModel(QModel%EmptyModel_t,QModel_in)

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


    IF (debug) write(out_unitp,*) 'init Q0 of HOO_DMBE (HOO minimum)'
    QModel%Q0 = [2.806_Rkind,2.271_Rkind,2.271_Rkind]

    IF (debug) write(out_unitp,*) 'init d0GGdef of HOO_DMBE'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      !CALL Write_HOO_DMBE_Model(QModel,nio=out_unitp)
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_HOO_DMBE_Model
!> @brief Subroutine wich prints the current HOO_DMBE_Model parameters.
!!
!! @param QModel            CLASS(HOO_DMBE_Model_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_HOO_DMBE_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(HOO_DMBE_Model_t),   intent(in) :: QModel
    integer,                   intent(in) :: nio

    write(nio,*) 'HOO_DMBE IV current parameters'

    CALL Write_EmptyModel(QModel%EmptyModel_t,nio)

    write(nio,*) 'end HOO_DMBE IV current parameters'

  END SUBROUTINE Write_HOO_DMBE_Model
!> @brief Subroutine wich prints the default HOO_DMBE_Model parameters.
!!
!! @param QModel            CLASS(HOO_DMBE_Model_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_HOO_DMBE_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(HOO_DMBE_Model_t),   intent(in) :: QModel
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


  END SUBROUTINE Write0_HOO_DMBE_Model

!> @brief Subroutine wich calculates the HOO_DMBE potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(HOO_DMBE_Model_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_HOO_DMBE_Pot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(HOO_DMBE_Model_t), intent(in)    :: QModel
    TYPE (dnS_t),            intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),            intent(in)    :: dnQ(:)
    integer,                 intent(in)    :: nderiv

    real(kind=Rkind) :: V,Q(3)

    Q(:) = QML_get_d0_FROM_dnS(dnQ)

    CALL HOO_DMBE4_pes(Q,V)

    CALL QML_set_dnS(Mat_OF_PotDia(1,1),d0=V)

  END SUBROUTINE eval_HOO_DMBE_Pot

  ! here we suppose that the atom ordering: H1-O2-O3
  SUBROUTINE Cart_TO_Q_HOO_DMBE_Model(QModel,dnX,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(HOO_DMBE_Model_t), intent(in)    :: QModel
    TYPE (dnS_t),            intent(in)    :: dnX(:,:)
    TYPE (dnS_t),            intent(inout) :: dnQ(:)
    integer,                 intent(in)    :: nderiv

    ! local vector
    integer         :: i,j
    TYPE (dnS_t)    :: VecOO(3),VecHO2(3),VecHO3(3)


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Cart_TO_Q_HOO_DMBE_Model'
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

    VecOO(:)  = dnX(:,3)-dnX(:,2)
    VecHO2(:) = dnX(:,2)-dnX(:,1)
    VecHO3(:) = dnX(:,3)-dnX(:,1)
    IF (debug) write(out_unitp,*) 'Cart_TO_Q_HOO_DMBE_Model vect done'

    dnQ(1) = sqrt(dot_product(VecOO,VecOO))
    dnQ(2) = sqrt(dot_product(VecHO2,VecHO2))
    dnQ(3) = sqrt(dot_product(VecHO3,VecHO3))

    CALL QML_dealloc_dnS(VecOO)
    CALL QML_dealloc_dnS(VecHO2)
    CALL QML_dealloc_dnS(VecHO3)

    IF (debug) THEN
      CALL QML_Write_dnS(dnQ(1),out_unitp,info='dnQ(1)')
      CALL QML_Write_dnS(dnQ(2),out_unitp,info='dnQ(2)')
      CALL QML_Write_dnS(dnQ(3),out_unitp,info='dnQ(3)')
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF
  END SUBROUTINE Cart_TO_Q_HOO_DMBE_Model

END MODULE mod_HOO_DMBE_Model
  SUBROUTINE HOO_DMBE4_pes(X,F)
  USE mod_NumParameters
  ! This is the DMBE IV potential energy surface for H + O2
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(3)
      COMMON/COEFF_HOO_DMBE4/C(52)
      COMMON/THRBOD_HOO_DMBE4/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3
      COMMON/REFGEO_HOO_DMBE4/R10,R20,R30

      R1=X(1) ! O-O distance
      R2=X(2) ! H-O(2) distance
      R3=X(3) ! H-O(3) distance

      Q1=ONE/SQRT(THREE)*(R1+R2+R3)
      Q2=ONE/SQRT(TWO)*(R2-R3)
      Q3=ONE/SQRT(SIX)*(TWO*R1-R2-R3)

      F=VOO_HOO_DMBE4(R1)+VOH_HOO_DMBE4(R2)+VOH_HOO_DMBE4(R3)+ &
        THREBQ_HOO_DMBE4(Q1,Q2,Q3)+ &
        EXDIS_HOO_DMBE4(R1,R2,R3)+ELECT_HOO_DMBE4(R1,R2,R3)

  END
  FUNCTION THREBQ_HOO_DMBE4(Q1,Q2,Q3)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/COEFF_HOO_DMBE4/C(52)
      COMMON/THRBOD_HOO_DMBE4/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3
      COMMON/REFGEO_HOO_DMBE4/R10,R20,R30
!     ****************************************************************
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

      THREBQ_HOO_DMBE4=POLQ*DECAY1*DECAY2*DECAY3

  END
  FUNCTION VOH_HOO_DMBE4(R)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      VOH_HOO_DMBE4=EHFOH_HOO_DMBE4(R)+DISOH_HOO_DMBE4(R)

  END
  FUNCTION EHFOH_HOO_DMBE4(R)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DIATDI_HOO_DMBE4/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DATA D,ASV/0.13825385_Rkind,2.6564788_Rkind,1.7450528_Rkind,0.71014391_Rkind,2.5453276_Rkind/
!     ****************************************************************
      X=R-RMOH
      R2=X*X
      R3=R2*X

      EHFOH_HOO_DMBE4=-D*(ONE+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*EXP(-ASV(4)*X)

  END
  FUNCTION DISOH_HOO_DMBE4(R)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_HOO_DMBE4/COO(10),COH(10)
      COMMON/DIATDI_HOO_DMBE4/R0OO,RMOO,R0OH,RMOH

      DISOH_HOO_DMBE4=DISP_HOO_DMBE4(R,COH(6),COH(8),COH(10),R0OH,RMOH)

  END
  FUNCTION VOO_HOO_DMBE4(R)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      VOO_HOO_DMBE4=EHFOO_HOO_DMBE4(R)+DISOO_HOO_DMBE4(R)

  END
  FUNCTION EHFOO_HOO_DMBE4(R)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DIATDI_HOO_DMBE4/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DATA D,ASV/0.14291202_Rkind,3.6445906_Rkind,3.9281238_Rkind,2.0986689_Rkind,3.3522498_Rkind/
!     ****************************************************************
      X=R-RMOO
      R2=X*X
      R3=R2*X

      EHFOO_HOO_DMBE4=-D*(ONE+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*EXP(-ASV(4)*X)

  END
  FUNCTION DISOO_HOO_DMBE4(R)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_HOO_DMBE4/COO(10),COH(10)
      COMMON/DIATDI_HOO_DMBE4/R0OO,RMOO,R0OH,RMOH

      DISOO_HOO_DMBE4=DISP_HOO_DMBE4(R,COO(6),COO(8),COO(10),R0OO,RMOO)

  END
  FUNCTION CEF_HOO_DMBE4(CAS,RK01,RK11,RK02,RK12,RE1,RE2,R1,R2)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      CEF_HOO_DMBE4=0.5_Rkind*CAS*((ONE-RK01*EXP(-RK11*(R1-RE1)))*TANH(RK12*R2)+&
       (ONE-RK02*EXP(-RK12*(R2-RE2)))*TANH(RK11*R1))

  END
  FUNCTION EXDIS_HOO_DMBE4 (R1,R2,R3)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_HOO_DMBE4/COO(10),COH(10)
      COMMON/RKVAL_HOO_DMBE4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DISCO_HOO_DMBE4/CEFOO(10),CEFOH2(10),CEFOH3(10),CEDOO(10),CEDOH2(10)
      COMMON/DISCO2_HOO_DMBE4/CEDOH3(10)
      COMMON/DIATDI_HOO_DMBE4/R0OO,RMOO,R0OH,RMOH

      DO IN=6,10,2
        CEFOO(IN)=CEF_HOO_DMBE4(COO(IN),RK0OH(IN),RK1OH(IN),RK0OH(IN),RK1OH(IN),RMOH,RMOH,R2,R3)
        CEDOO(IN)=CEFOO(IN)-COO(IN)
        CEFOH2(IN)=CEF_HOO_DMBE4(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN),RMOO,RMOH,R1,R3)
        CEDOH2(IN)=CEFOH2(IN)-COH(IN)
        CEFOH3(IN)=CEF_HOO_DMBE4(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN),RMOO,RMOH,R1,R2)
        CEDOH3(IN)=CEFOH3(IN)-COH(IN)
      END DO

      EXDIS_HOO_DMBE4=DISP_HOO_DMBE4(R1,CEDOO(6),CEDOO(8),CEDOO(10),R0OO,RMOO) +&
          DISP_HOO_DMBE4(R2,CEDOH2(6),CEDOH2(8),CEDOH2(10),R0OH,RMOH) +         &
          DISP_HOO_DMBE4(R3,CEDOH3(6),CEDOH3(8),CEDOH3(10),R0OH,RMOH)

  END
  FUNCTION ELECT_HOO_DMBE4(R1,R2,R3)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/POLAR_HOO_DMBE4/C4,C5
      COMMON/RKVAL_HOO_DMBE4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DIATDI_HOO_DMBE4/R0OO,RMOO,R0OH,RMOH
      COMMON/DAMPC_HOO_DMBE4/ADAMP(10),BDAMP(10)
      COMMON/WELECT_HOO_DMBE4/C4OHR2,C5OHR2,C4OHR3,C5OHR3,C4OO,C5OO,TERM4,TERM5

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

    TERM4=C4OO/R14*(ONE-EXP(-ADAMP(4)*RROO-BDAMP(4)*RROO**2))**4+         &
          C4OHR2/R24*(ONE-EXP(-ADAMP(4)*RROH2-BDAMP(4)*RROH2**2))**4+     &
          C4OHR3/R34*(ONE-EXP(-ADAMP(4)*RROH3-BDAMP(4)*RROH3**2))**4

    TERM5=C5OO/R15*(ONE-EXP(-ADAMP(5)*RROO-BDAMP(5)*RROO**2))**5+         &
          C5OHR2/R25*(ONE-EXP(-ADAMP(5)*RROH2-BDAMP(5)*RROH2**2))**5+     &
          C5OHR3/R35*(ONE-EXP(-ADAMP(5)*RROH3-BDAMP(5)*RROH3**2))**5

    ELECT_HOO_DMBE4=TERM4+TERM5

  END
  FUNCTION DISP_HOO_DMBE4(R,C6,C8,C10,R0,RM)
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DAMPC_HOO_DMBE4/ADAMP(10),BDAMP(10)

      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=TWO*R/(RM+2.5_Rkind*R0)
      D6=(ONE-EXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6
      D8=(ONE-EXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8
      D10=(ONE-EXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10

      DISP_HOO_DMBE4=-C6/R6*D6-C8/R8*D8-C10/R10*D10

  END
  BLOCK DATA HO2DAT_HOO_DMBE4
  USE mod_NumParameters
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/COEFF_HOO_DMBE4/C(52)
      COMMON/DISPC_HOO_DMBE4/COO(10),COH(10)
      COMMON/DIATDI_HOO_DMBE4/R0OO,RMOO,R0OH,RMOH
      COMMON/RKVAL_HOO_DMBE4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/POLAR_HOO_DMBE4/C4,C5
      COMMON/DAMPC_HOO_DMBE4/ADAMP(10),BDAMP(10)
      COMMON/REFGEO_HOO_DMBE4/R10,R20,R30
!     ***************************************************************
      DATA C/ &
       .49040645E+01_Rkind, -.86748216E+01_Rkind,  .50555792E+01_Rkind,  .42941301E+01_Rkind,&
      -.41874792E+01_Rkind,  .13461379_Rkind,     -.99064922_Rkind,      .13358488E+01_Rkind,&
       .13495231E+01_Rkind, -.18529696_Rkind,     -.23534213E+02_Rkind,  .24289930E+02_Rkind,&
      -.50209026E+01_Rkind, -.10365484E+02_Rkind,  .46692224E+01_Rkind, -.14747138E+01_Rkind,&
       .23119718E+01_Rkind, -.18247842E+01_Rkind, -.28472166_Rkind,      .51036509_Rkind,&
       .19124083_Rkind,      .45405729E+01_Rkind,  .11087611_Rkind,      -.19990481_Rkind,&
      -.37356178_Rkind,      .46142042E-01_Rkind, -.20565580_Rkind,      -.27015963_Rkind,&
       .34085281_Rkind,      .28321162_Rkind,     -.11558481_Rkind,      -.29448886_Rkind,&
      -.52932488_Rkind,      .58159523E-01_Rkind, -.48649560E-02_Rkind,   .11949167E-01_Rkind,&
       .21409804E-01_Rkind, -.20620608E-02_Rkind,  .30177088E-01_Rkind,   .27880291E-01_Rkind,&
       .88458711E-02_Rkind,  .13137410E-01_Rkind, -.24705619E-01_Rkind,   -.31085889E-01_Rkind,&
       .34317857E-02_Rkind,  .52593878E-01_Rkind,  .79500714E-01_Rkind,    -.79782216E-02_Rkind,&
       .31164575E-01_Rkind, -.28737598E-01_Rkind,  .98201698_Rkind,         .62000000_Rkind/
    DATA R0OO,RMOO,R0OH,RMOH/5.661693_Rkind,2.2818_Rkind,6.294894_Rkind,1.8344_Rkind/
    DATA COO/0.0_Rkind,0.0_Rkind,0.0_Rkind,0.0_Rkind,0._Rkind,15.40_Rkind,0.0_Rkind,235.219943_Rkind,&
             0.0_Rkind,4066.23929_Rkind/
    DATA COH/0.0_Rkind,0.0_Rkind,0.0_Rkind,0.0_Rkind,0._Rkind,10.00_Rkind,0.0_Rkind,180.447673_Rkind,&
              0.0_Rkind,3685.25842_Rkind/
    DATA C4,C5/-0.92921_Rkind,-1.79000_Rkind/
    DATA RK0OO/0.0_Rkind,0.0_Rkind,0.0_Rkind,0.0_Rkind,0.0_Rkind,-.27847758_Rkind,0.0_Rkind,&
                -.46815641_Rkind,0.0_Rkind,-1.20506384_Rkind/
    DATA RK1OO/0.0_Rkind,0.0_Rkind,0.0_Rkind,3.35224980_Rkind,3.35224980_Rkind,&
                0.95273753_Rkind,0.0_Rkind,0.94148408_Rkind,0.0_Rkind,0.72379129_Rkind/
    DATA RK0OH/0.0_Rkind,0.0_Rkind,0.0_Rkind,0.0_Rkind,0.0_Rkind,0.02465005_Rkind,0.0_Rkind,&
                0.05036950_Rkind,0.0_Rkind,0.06294371_Rkind/
    DATA RK1OH/0.0_Rkind,0.0_Rkind,0.0_Rkind,2.54532760_Rkind,2.54532760_Rkind,&
                0.68758097_Rkind,0.0_Rkind,0.82542359_Rkind,0.0_Rkind,0.94034225_Rkind/
    DATA ADAMP/0.0_Rkind,0.0_Rkind,0.0_Rkind,5.0079875_Rkind,3.8428294_Rkind,3.0951333_Rkind,&
                0.0_Rkind,2.1999000_Rkind,0.0_Rkind,1.6880714_Rkind/
    DATA BDAMP/0.0_Rkind,0.0_Rkind,0._Rkind,10.6645006_Rkind,9.6758155_Rkind,8.7787895_Rkind,&
                0.0_Rkind,7.2265123_Rkind,0.0_Rkind,5.9487108_Rkind/
    DATA R10,R20,R30/2.5143000_Rkind,2.6469057_Rkind,2.6469057_Rkind/
  END
