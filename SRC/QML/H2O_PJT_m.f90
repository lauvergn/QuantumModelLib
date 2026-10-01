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

!> @brief Module which makes the initialization, calculation of the H2O_PJT potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 30/09/2026
!!
MODULE QML_H2O_PJT_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the H2O_PJT parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE, EXTENDS (QML_Empty_t) ::  QML_H2O_PJT_t

     PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

  CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_H2O_PJT
    PROCEDURE :: Write_QModel     => Write_QML_H2O_PJT
    PROCEDURE :: RefValues_QModel => RefValues_QML_H2O_PJT
  END TYPE QML_H2O_PJT_t

  PUBLIC :: QML_H2O_PJT_t,Init_QML_H2O_PJT

  CONTAINS
!> @brief Subroutine which makes the initialization of the H2O_PJT parameters.
!!
!! @param H2O_PJTPot          TYPE(QML_H2O_PJT_t):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_H2O_PJT(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE (QML_H2O_PJT_t)                         :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !local variable
    integer                     :: err_read,nio_fit,i,j,k,idum

    real(kind=Rkind), parameter :: TOANG = 0.5291772_Rkind ! from the PJT routine
    real(kind=Rkind), parameter :: TORAD = 3.141592654_Rkind/180._Rkind
    real(kind=Rkind), parameter :: PILOC = 3.141592654_Rkind

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H2O_PJT'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf    = 1
    QModel%pot_name = 'H2O_PJT'
    QModel%ndim     = 3


    IF (QModel%option /= 2) QModel%option = 2

    SELECT CASE (QModel%option)
    CASE (1)

      QModel%d0GGdef = Identity_Mat(QModel%ndim)

      QModel%Qref = [0.95792059_Rkind/TOANG,0.95792059_Rkind/TOANG,(180._Rkind-75.50035308_Rkind)*TORAD]

      IF (QModel%PubliUnit) THEN
        write(out_unit,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad], Energy: [Hartree]'
      ELSE
        write(out_unit,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad]:, Energy: [Hartree]'
      END IF
    CASE (2)

      QModel%d0GGdef = Identity_Mat(QModel%ndim)

      QModel%Qref = [0.95792059_Rkind/TOANG,0.95792059_Rkind/TOANG,(180._Rkind-75.50035308_Rkind)*TORAD]

      IF (QModel%PubliUnit) THEN
        write(out_unit,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad], Energy: [Hartree]'
      ELSE
        write(out_unit,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad]:, Energy: [Hartree]'
      END IF
    CASE Default

      write(out_unit,*) ' ERROR in Init_QML_H2O_PJT '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 2'
      STOP 'ERROR in Init_QML_H2O_PJT: wrong option'

    END SELECT


    IF (debug) write(out_unit,*) 'init Q0 of H2O_PJT'
    allocate(QModel%Q0(QModel%ndim))
    CALL get_Q0_QML_H2O_PJT(QModel%Q0,QModel,option=0)
    IF (debug) write(out_unit,*) 'QModel%Q0',QModel%Q0

    IF (debug) write(out_unit,*) 'init d0GGdef of H2O_PJT'
    flush(out_unit)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_H2O_PJT
!> @brief Subroutine wich prints the QML_H2O_PJT parameters.
!!
!! @param QModel            CLASS(QML_H2O_PJT_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_QML_H2O_PJT(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_H2O_PJT_t), intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'H2O_PJT current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '    R2  \ a                            '
    write(nio,*) '         O------------H               '
    write(nio,*) '              R1                       '
    write(nio,*) '                                       '
    write(nio,*) '  Coordinates (option 1):          '
    write(nio,*) '  Q(1) = R1         (Bohr)             '
    write(nio,*) '  Q(2) = R2         (Bohr)             '
    write(nio,*) '  Q(3) = a          (Radian)           '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) ' Water potential, from:                '
    write(nio,*) ' Polyansky, Jensen and Tennyson        '
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)
    write(nio,*) '---------------------------------------'

    SELECT CASE (QModel%option)

    CASE (1)
      write(nio,*) '                                       '
      write(nio,*) '  NOT YET                              '
      write(nio,*) '                                       '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = R1   = 1.8107934895553828 (Bohr)'
      write(nio,*) '  Q(2) = R2   = 1.8107934895553828 (Bohr)'
      write(nio,*) '  Q(3) = a    = 1.8230843491285216 (Rad) '
      write(nio,*) '                                       '
      write(nio,*) '  V = 0.0          Hartree             '
      write(nio,*) ' grad(:) =[0.0,0.0,0.0]                '
      write(nio,*) ' hess    =[0.4726717042854987, 0.00 ,0.00'
      write(nio,*) '           0.00, 0.4726717042854987, 0.00'
      write(nio,*) '           0.00, 0.00, 0.1085242579142319]'
      write(nio,*) 'From, Polyansky, Jensen and Tennyson   '
      write(nio,*) 'J Chem Phys 101, 7651 (1994)           '
      write(nio,*) '                                       '
    CASE (2)
      write(nio,*) 'From, Polyansky, Jensen and Tennyson   '
      write(nio,*) ' J. Chem. Phys., 105, 6490-6497 (1996) '
      write(nio,*) 'Update from J Chem Phys 101, 7651 (1994)'
    CASE Default
        write(out_unit,*) ' ERROR in write_QModel '
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1'
        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end H2O_PJT current parameters'

  END SUBROUTINE Write_QML_H2O_PJT

  SUBROUTINE get_Q0_QML_H2O_PJT(Q0,QModel,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (QML_H2O_PJT_t),        intent(in)    :: QModel
    integer,                     intent(in)    :: option

    IF (size(Q0) /= 3) THEN
      write(out_unit,*) ' ERROR in get_Q0_QML_H2O_PJT '
      write(out_unit,*) ' The size of Q0 is not ndim=3: '
      write(out_unit,*) ' size(Q0)',size(Q0)
      STOP 'ERROR in get_Q0_QML_H2O_PJT: wrong Q0 size'
    END IF


    SELECT CASE (QModel%option)
    CASE (1) ! R1,R2,a
      Q0(:) = QModel%Qref

    CASE (2) ! R1,R2,a
      Q0(:) = QModel%Qref

    CASE Default
      write(out_unit,*) ' ERROR in get_Q0_QML_H2O_PJT '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1,2'
      STOP 'ERROR in get_Q0_QML_H2O_PJT: wrong option'
    END SELECT

  END SUBROUTINE get_Q0_QML_H2O_PJT
!> @brief Subroutine wich calculates the H2O_PJT potential (unpublished model) with derivatives.
!!
!! @param PotVal             TYPE (dnMat_t):      Potential with derivatives,.
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_H2O_PJT_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H2O_PJT(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2O_PJT_t), intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable :: dnQsym(:)


    SELECT CASE (QModel%option)
    !CASE (1) ! R1,R2,a
      !CALL EvalPot1_QML_H2O_PJT(QModel,Mat_OF_PotDia,dnQ,nderiv)
    CASE (2) ! R1,R2,a
      CALL EvalPot2_QML_H2O_PJT(QModel,Mat_OF_PotDia,dnQ,nderiv)
    CASE Default
      write(out_unit,*) ' ERROR in EvalPot_QML_H2O_PJT '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 2'
      STOP 'ERROR in EvalPot_QML_H2O_PJT: wrong option'
    END SELECT


  END SUBROUTINE EvalPot_QML_H2O_PJT

  SUBROUTINE EvalPot2_QML_H2O_PJT(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2O_PJT_t), intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv

    real(kind=Rkind), parameter :: PILOC  = 3.141592654_Rkind ! from the PJT routine
    real(kind=Rkind), parameter :: TOANG  = 0.5291772_Rkind   ! from the PJT routine
    real(kind=Rkind), parameter :: TORAD  = PILOC/180._Rkind  ! from the PJT routine
    real(kind=Rkind), parameter :: CMTOAU = 219474.624_Rkind  ! from the PJT routine

    real(kind=Rkind), parameter ::  X1      =     1.0_Rkind
    real(kind=Rkind), parameter ::  RHO1    =    75.50035308_Rkind
    real(kind=Rkind), parameter ::  FA1     =      .00000000_Rkind
    real(kind=Rkind), parameter ::  FA2     = 18902.44193433_Rkind
    real(kind=Rkind), parameter ::  FA3     =  1893.99788146_Rkind
    real(kind=Rkind), parameter ::  FA4     =  4096.73443772_Rkind
    real(kind=Rkind), parameter ::  FA5     = -1959.60113289_Rkind
    real(kind=Rkind), parameter ::  FA6     =  4484.15893388_Rkind
    real(kind=Rkind), parameter ::  FA7     =  4044.55388819_Rkind
    real(kind=Rkind), parameter ::  FA8     = -4771.45043545_Rkind
    real(kind=Rkind), parameter ::  FA9     =     0.00000000_Rkind
    real(kind=Rkind), parameter ::  FA10    =     0.00000000_Rkind
    real(kind=Rkind), parameter ::  RZ      =      .95792059_Rkind
    real(kind=Rkind), parameter ::  A       =     2.22600000_Rkind
    real(kind=Rkind), parameter ::  F1A1    = -6152.40141181_Rkind
    real(kind=Rkind), parameter ::  F2A1    = -2902.13912267_Rkind
    real(kind=Rkind), parameter ::  F3A1    = -5732.68460689_Rkind
    real(kind=Rkind), parameter ::  F4A1    =   953.88760833_Rkind
    real(kind=Rkind), parameter ::  F11     = 42909.88869093_Rkind
    real(kind=Rkind), parameter ::  F1A11   = -2767.19197173_Rkind
    real(kind=Rkind), parameter ::  F2A11   = -3394.24705517_Rkind
    real(kind=Rkind), parameter ::  F3A11   =      .00000000_Rkind
    real(kind=Rkind), parameter ::  F13     = -1031.93055205_Rkind
    real(kind=Rkind), parameter ::  F1A13   =  6023.83435258_Rkind
    real(kind=Rkind), parameter ::  F2A13   =      .00000000_Rkind
    real(kind=Rkind), parameter ::  F3A13   =      .00000000_Rkind
    real(kind=Rkind), parameter ::  F111    =      .00000000_Rkind
    real(kind=Rkind), parameter ::  F1A111  =   124.23529382_Rkind
    real(kind=Rkind), parameter ::  F2A111  = -1282.50661226_Rkind
    real(kind=Rkind), parameter ::  F113    = -1146.49109522_Rkind
    real(kind=Rkind), parameter ::  F1A113  =  9884.41685141_Rkind
    real(kind=Rkind), parameter ::  F2A113  =  3040.34021836_Rkind
    real(kind=Rkind), parameter ::  F1111   =  2040.96745268_Rkind
    real(kind=Rkind), parameter ::  FA1111  =      .00000000_Rkind
    real(kind=Rkind), parameter ::  F1113   =  -422.03394198_Rkind
    real(kind=Rkind), parameter ::  FA1113  = -7238.09979404_Rkind
    real(kind=Rkind), parameter ::  F1133   =      .00000000_Rkind
    real(kind=Rkind), parameter ::  FA1133  =      .00000000_Rkind
    real(kind=Rkind), parameter ::  F11111  = -4969.24544932_Rkind
    real(kind=Rkind), parameter ::  f111111 =  8108.49652354_Rkind
    real(kind=Rkind), parameter ::  F71     =    90.00000000_Rkind

    real(kind=Rkind), parameter ::  Fa11    = 0.0
    real(kind=Rkind), parameter ::  F1a3    = F1a1
    real(kind=Rkind), parameter ::  F2a3    = F2a1
    real(kind=Rkind), parameter ::  F3a3    = F3a1
    real(kind=Rkind), parameter ::  F4a3    = F4a1
    real(kind=Rkind), parameter ::  F33     = F11
    real(kind=Rkind), parameter ::  F1a33   = F1a11
    real(kind=Rkind), parameter ::  F2a33   = F2a11
    real(kind=Rkind), parameter ::  F333    = F111
    real(kind=Rkind), parameter ::  F1a333  = F1a111
    real(kind=Rkind), parameter ::  F2a333  = F2a111
    real(kind=Rkind), parameter ::  F133    = F113
    real(kind=Rkind), parameter ::  F1a133  = F1a113
    real(kind=Rkind), parameter ::  F2a133  = F2a113
    real(kind=Rkind), parameter ::  F3333   = F1111
    real(kind=Rkind), parameter ::  Fa3333  = Fa1111
    real(kind=Rkind), parameter ::  F1333   = F1113
    real(kind=Rkind), parameter ::  Fa1333  = Fa1113
    real(kind=Rkind), parameter ::  F33333  = F11111
    real(kind=Rkind), parameter ::  F333333 = F111111
    real(kind=Rkind), parameter ::  F73     = F71

    ! RZ  = OH equilibrium value
    ! RHO = equilibrium value of pi - bond angle(THETA)

    real(kind=Rkind), parameter :: c1     = 50._Rkind
    real(kind=Rkind), parameter :: c2     = 10.0_Rkind
    real(kind=Rkind), parameter :: beta1  = 22.0_Rkind
    real(kind=Rkind), parameter :: beta2  = 13.5_Rkind
    real(kind=Rkind), parameter :: gammas = 0.05_Rkind
    real(kind=Rkind), parameter :: gammaa = 0.10_Rkind
    real(kind=Rkind), parameter :: delta  = 0.85_Rkind
    real(kind=Rkind), parameter :: rhh0   = 1.40_Rkind
    real(kind=Rkind), parameter :: RHO    = RHO1*TORAD

    ! modification by Choi & Light, J. Chem. Phys., 97, 7031 (1992).
    real(kind=Rkind), parameter :: sqrt2  = sqrt(TWO)
    real(kind=Rkind), parameter :: xmup1  = sqrt2/THREE+HALF
    real(kind=Rkind), parameter :: xmum1  = xmup1-x1
    TYPE (dnS_t) :: term, r1, r2, rhh, rbig, rlit
    TYPE (dnS_t) :: alpha, alpha1, alpha2, drhh, DOLEG

    TYPE (dnS_t) :: Q1,Q2,THETA
    TYPE (dnS_t) :: DR,DS,Y1,Y3,CORO
    TYPE (dnS_t) :: V,V0
    TYPE (dnS_t) :: FE1,FE3,FE11,FE33,FE13,FE111,FE333,FE113,FE133,FE1111,FE3333,FE1113
    TYPE (dnS_t) :: FE1333,FE1133,FE11111,FE33333,FE111111,FE333333,FE71,FE73

    Q1    = dnQ(1)
    Q2    = dnQ(2)
    THETA = dnQ(3)

    ! fa11=0.0
    ! f1a3=f1a1
    ! f2a3=f2a1
    ! f3a3=f3a1
    ! f4a3=f4a1
    ! f33=f11
    ! f1a33=f1a11
    ! f2a33=f2a11
    ! f333=f111
    ! f1a333=f1a111
    ! f2a333=f2a111
    ! f133=f113
    ! f1a133=f1a113
    ! f2a133=f2a113
    ! f3333=f1111
    ! fa3333=fa1111
    ! f1333=f1113
    ! fa1333=fa1113
    ! f33333=f11111
    ! f333333 =f111111
    ! f73     =f71

    ! Find value for DR and DS
    DR = TOANG*Q1 - RZ
    DS = TOANG*Q2 - RZ

    ! Transform to Morse coordinates
    Y1 = X1 - EXP(-A * DR)
    Y3 = X1 - EXP(-A * DS)

    ! transform to Jensens angular coordinate
    CORO = COS(THETA) + COS(RHO)

    ! Now for the potential
    V0=(FA2+FA3*CORO+FA4*CORO**2+FA6*CORO**4+FA7*CORO**5)*CORO**2
    V0=V0+(FA8*CORO**6+FA5*CORO**3+FA9*CORO**7+FA10*CORO**8 )*CORO**2
    V0=V0+(                                    FA11*CORO**9 )*CORO**2
    FE1= F1A1*CORO+F2A1*CORO**2+F3A1*CORO**3+F4A1*CORO**4
    FE3= F1A3*CORO+F2A3*CORO**2+F3A3*CORO**3+F4A3*CORO**4
    FE11= F11+F1A11*CORO+F2A11*CORO**2
    FE33= F33+F1A33*CORO+F2A33*CORO**2
    FE13= F13+F1A13*CORO
    FE111= F111+F1A111*CORO+F2A111*CORO**2
    FE333= F333+F1A333*CORO+F2A333*CORO**2
    FE113= F113+F1A113*CORO+F2A113*CORO**2
    FE133= F133+F1A133*CORO+F2A133*CORO**2
    FE1111= F1111+FA1111*CORO
    FE3333= F3333+FA3333*CORO
    FE1113= F1113+FA1113*CORO
    FE1333= F1333+FA1333*CORO
    FE1133=       FA1133*CORO
    FE11111=F11111
    FE33333=F33333
    FE111111=F111111
    FE333333=F333333
    FE71    =F71
    FE73    =F73
    V   = V0 +  FE1*Y1+FE3*Y3                             &
             +  FE11*Y1**2+FE33*Y3**2+FE13*Y1*Y3          &
             +  FE111*Y1**3+FE333*Y3**3+FE113*Y1**2*Y3    &
             +  FE133*Y1*Y3**2                            &
             +  FE1111*Y1**4+FE3333*Y3**4+FE1113*Y1**3*Y3 &
             +  FE1333*Y1*Y3**3+FE1133*Y1**2*Y3**2        &
             +  FE11111*Y1**5+FE33333*Y3**5               &
             +  FE111111*Y1**6+FE333333*Y3**6             &
             +  FE71    *Y1**7+FE73    *Y3**7
    ! modification by Choi & Light, J. Chem. Phys., 97, 7031 (1992).
    term = TWO*xmum1*xmup1*q1*q2*cos(theta)
    r1   = toang*sqrt((xmup1*q1)**2+(xmum1*q2)**2-term)
    r2   = toang*sqrt((xmum1*q1)**2+(xmup1*q2)**2-term)
    rhh  = sqrt(q1**2+q2**2-TWO*q1*q2*cos(theta))
    rbig = (r1+r2)/sqrt2
    rlit = (r1-r2)/sqrt2

    alpha  = (x1-tanh(gammas*rbig**2))*(x1-tanh(gammaa*rlit**2))
    alpha1 = beta1*alpha
    alpha2 = beta2*alpha
    drhh   = toang*(rhh-delta*rhh0)
    DOLEG  = (1.4500_Rkind-THETA)
    v = v + c1*exp(-alpha1*drhh) + c2*exp(-alpha2*drhh)

    ! Convert to Hartree
    Mat_OF_PotDia(1,1) = V/CMTOAU

   !write(out_unit,*) ' end EvalPot2_QML_H2O_PJT' ; flush(6)

  END SUBROUTINE EvalPot2_QML_H2O_PJT

  SUBROUTINE RefValues_QML_H2O_PJT(QModel,err,nderiv,Q0,dnMatV,d0GGdef,option)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2O_PJT_t), intent(in)              :: QModel
    integer,              intent(inout)           :: err
    integer,              intent(in)              :: nderiv

    real (kind=Rkind),    intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),       intent(inout), optional :: dnMatV
    real (kind=Rkind),    intent(inout), optional :: d0GGdef(:,:)
    integer,              intent(in),    optional :: option

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_H2O_PJT'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      flush(out_unit)
    END IF

    IF (.NOT. QModel%Init) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'The model is not initialized!'
      err = -1
      RETURN
    ELSE
      err = 0
    END IF

    SELECT CASE (option)
    !CASE (1) ! first version of the PJT potential
    !  CONTINUE
    CASE (2) ! second version of the PJT potential
      IF (present(Q0))      CALL RefValues_QML_H2O_PJT_2(QModel,err,nderiv=nderiv,Q0=Q0)
      IF (present(dnMatV))  CALL RefValues_QML_H2O_PJT_2(QModel,err,nderiv=nderiv,dnMatV=dnMatV)
      IF (present(d0GGdef)) CALL RefValues_QML_H2O_PJT_2(QModel,err,nderiv=nderiv,d0GGdef=d0GGdef)
    CASE Default
      STOP 'ERROR in RefValues_QML_H2O_PJT: wrong option. Possible values: 2'
    END SELECT


    IF (debug) THEN
      write(out_unit,*) 'present Q0 dnMatV d0GGdef',present(Q0),present(dnMatV),present(d0GGdef)
      IF (present(Q0))      write(out_unit,*) 'Q0',Q0
      IF (present(dnMatV))  THEN
        write(out_unit,*) 'dnMatV is present'
        write(out_unit,*) 'dnMatV is allocated',(.NOT. Check_NotAlloc_dnMat(dnMatV,nderiv))
        CALL write_dnMat(dnMatV,info='dnMatV')
      END IF
      IF (present(d0GGdef)) write(out_unit,*) 'd0GGdef',d0GGdef
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_H2O_PJT
  SUBROUTINE RefValues_QML_H2O_PJT_2(QModel,err,Q0,dnMatV,d0GGdef,nderiv)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2O_PJT_t), intent(in)              :: QModel

    integer,           intent(inout)           :: err

    integer,           intent(in)              :: nderiv
    real (kind=Rkind), intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),    intent(inout), optional :: dnMatV

    real (kind=Rkind), intent(inout), optional :: d0GGdef(:,:)

    real (kind=Rkind), allocatable :: d0(:,:),d1(:,:,:),d2(:,:,:,:),d3(:,:,:,:,:),V(:)
    integer        :: i,n0,n1,n2

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_H2O_PJT_2'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      flush(out_unit)
    END IF

    IF (.NOT. QModel%Init) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'The model is not initialized!'
      err = -1
      RETURN
    ELSE
      err = 0
    END IF

    IF (present(Q0)) THEN
      IF (size(Q0) /= QModel%ndim) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'incompatible Q0 size:'
        write(out_unit,*) 'size(Q0), ndimQ:',size(Q0),QModel%ndim
        err = 1
        Q0(:) = HUGE(ONE)
        RETURN
      END IF
      Q0(:) = [1.8_Rkind,1.7_Rkind,1.8_Rkind]
    END IF

    n0 = QModel%nsurf*QModel%nsurf
    n1 = QModel%nsurf*QModel%nsurf*QModel%ndim
    n2 = QModel%nsurf*QModel%nsurf*QModel%ndim**2

    ! 3.9241287205339535E-003 ! from the original pot
    ! 3.9241223130162426E-003_Rkind, from QML (error 6e-9, due to the conversion single to double precision of the parameters)
    V  = [3.9241223130162426E-003_Rkind,                                                               &
         -5.8454569221278469E-003_Rkind,-7.3625140821373239E-002_Rkind,-7.8925822708984916E-003_Rkind, &
          0.56343652633495600_Rkind,    -4.7395997875263312E-003_Rkind, 3.6047904034701532E-002_Rkind, &
         -4.7395997875263312E-003_Rkind, 0.79594606285205305_Rkind,     3.2817176363485974E-002_Rkind, &
          3.6047904034701532E-002_Rkind, 3.2817176363485974E-002_Rkind, 0.17020160680264132_Rkind]
    
    IF (present(dnMatV)) THEN
      err = 0
      IF (nderiv >= 0) THEN ! no derivative

        d0 = reshape(V(1:n0),shape=[QModel%nsurf,QModel%nsurf])
      END IF

      IF (nderiv >= 1) THEN ! 1st order derivatives
        d1 = reshape(V(1+n0:n0+n1),shape=[QModel%nsurf,QModel%nsurf,QModel%ndim])
      END IF

      IF (nderiv >= 2) THEN ! 2d order derivatives
        d2 = reshape(V(1+n0+n1:n0+n1+n2),shape=[QModel%nsurf,QModel%nsurf,QModel%ndim,QModel%ndim])
      END IF
      SELECT CASE (nderiv)
      CASE(0)
        CALL set_dnMat(dnMatV,d0=d0)
      CASE(1)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1)
      CASE(2)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1,d2=d2)
      CASE Default
        STOP 'ERROR in RefValues_QML_H2O_PJT_2: nderiv MUST < 3'
      END SELECT

    END IF

    IF (present(d0GGdef)) d0GGdef = Identity_Mat(QModel%ndim)


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_H2O_PJT_2
END MODULE QML_H2O_PJT_m
