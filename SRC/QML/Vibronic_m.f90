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
!!
!> @author David Lauvergnat
!! @date 07/02/2023
!!
MODULE QML_Vibronic_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE ADdnSVM_m, ONLY : dnS_t
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

  TYPE QML_Diabatic_t
    real (kind=Rkind), allocatable :: Qref(:)
    TYPE (dnS_t)                   :: Ene
  END TYPE QML_Diabatic_t

  TYPE, EXTENDS (QML_Empty_t) ::  QML_Vibronic_t
    PRIVATE

    TYPE (QML_Diabatic_t), allocatable :: Diab(:,:)
    character (len=:),     allocatable :: Vibronic_name

     ! Diab(i,j), the diabatic potentials (i=j) or couplings (i/=j), are
     !   expanded as a taylor expension at second order around a reference geometry Diab(i,j)%Qref:
     !      Ediab(i,j)(Q) = Diab(i,j)%Ene%d0 + 
     !            Sum_k Diab(i,j)%Ene%d1(k) * (Q(k)-Diab(i,j)%Qref(k)) + 
     !    1/2*Sum_kl Diab(i,j)%Ene%d2(k,l) * (Q(k)-Diab(i,j)%Qref(k)) * (Q(l)-Diab(i,j)%Qref(l))

  CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_Vibronic
    PROCEDURE :: Write_QModel    => Write_QML_Vibronic
  END TYPE QML_Vibronic_t

  PUBLIC :: QML_Vibronic_t, Init_QML_Vibronic

  CONTAINS
!> @brief Subroutine which makes the initialization of a vibronic models.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 07/02/2023
!!
!! @param QModel             TYPE(QML_Vibronic_t):   derived type in which the parameters are set-up.
!! @param PubliUnit          logical (optional):   when PubliUnit=.TRUE., the units (Angstrom and eV) are used. Default (atomic unit).
  FUNCTION Init_QML_Vibronic(QModel_in,read_param,nio_param_file,Vibronic_name) RESULT(QModel)
    USE QDUtil_m,         ONLY : Read_Vec, Read_Mat, Identity_Mat, TO_string
    USE ADdnSVM_m,        ONLY : dnS_t, set_dnS, get_nderiv

    IMPLICIT NONE

    TYPE (QML_Vibronic_t)                             :: QModel

    TYPE(QML_Empty_t),           intent(in)           :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)           :: nio_param_file
    logical,                     intent(in)           :: read_param
    character (len=*),           intent(in), optional :: Vibronic_name

    integer :: i,j
    real (kind=Rkind), allocatable :: d1(:),d2(:,:)


    QModel%QML_Empty_t = QModel_in
    IF (present(Vibronic_name)) QModel%Vibronic_name = Vibronic_name

    QModel%pot_name      = 'Vibronic'

    IF (read_param) THEN
      STOP 'STOP in Init_QML_Vibronic: read data not yet!'
    ELSE IF (allocated(QModel%Vibronic_name)) THEN
      CALL Internal_QML_Vibronic(QModel)
    ELSE ! default model
      QModel%ndim          = 1
      QModel%nsurf         = 2
      QModel%Vibronic_name = 'nsurf' // TO_string(QModel%nsurf) // '_ndim' // TO_string(QModel%ndim)

      allocate(QModel%Diab(QModel%nsurf,QModel%nsurf))
      allocate(d1(QModel%ndim))
      d1 = ZERO
      d2 = Identity_Mat(QModel%ndim)
      ! quadratique potential
      DO i=1,QModel%nsurf
        allocate(QModel%Diab(i,i)%Qref(QModel%ndim))
        QModel%Diab(i,i)%Qref = ZERO
        CALL set_dnS(QModel%Diab(i,i)%Ene,d0=real(i-1,Rkind),d1=d1,d2=d2)
      END DO

      ! linear coupling
      d1 = ZERO
      d1(1) = ONE
      DO i=1,QModel%nsurf
      DO j=i+1,QModel%nsurf
        allocate(QModel%Diab(j,i)%Qref(QModel%ndim))
        QModel%Diab(j,i)%Qref = ZERO
        CALL set_dnS(QModel%Diab(j,i)%Ene,d0=ZERO,d1=d1)
        QModel%Diab(i,j)%Ene  = QModel%Diab(j,i)%Ene
        QModel%Diab(i,j)%Qref = QModel%Diab(j,i)%Qref
      END DO
      END DO

    END IF
  END FUNCTION Init_QML_Vibronic

  SUBROUTINE Internal_QML_Vibronic(QModel)
    USE QDUtil_m,         ONLY : TO_string
    USE ADdnSVM_m,        ONLY : dnS_t,Write_dnS
    IMPLICIT NONE

    CLASS (QML_Vibronic_t), intent(inout)   :: QModel

    SELECT CASE (QModel%Vibronic_name)
    CASE ('twod_rjdi2014')
      CALL Internal_QML_RJDI2014(QModel)
    CASE ('7dm23_jgbl2023_test')
      CALL Internal_QML_7DM23_JGBL2023(QModel)
    CASE ('so2')
      CALL Internal_QML_SO2(QModel)
    CASE ('xxx')
      CALL Internal_QML_XXX(QModel)
    CASE Default
      STOP 'no default in Internal_QML_Vibronic'
    END SELECT

  END SUBROUTINE Internal_QML_Vibronic
  SUBROUTINE Internal_QML_RJDI2014(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    USE ADdnSVM_m,        ONLY : dnS_t, set_dnS
    IMPLICIT NONE

    CLASS (QML_Vibronic_t), intent(inout)   :: QModel

    real(kind=Rkind), parameter     :: w1    = 9.557e-3_Rkind
    real(kind=Rkind), parameter     :: w2    = 3.3515e-3_Rkind
    real(kind=Rkind), parameter     :: Delta = 0.01984_Rkind
    real(kind=Rkind), parameter     :: a     = 20.07_Rkind
    real(kind=Rkind), parameter     :: c     = 6.127e-4_Rkind
 
    real (kind=Rkind), parameter    :: muX  = 1._Rkind
    real (kind=Rkind), parameter    :: muY  = 1._Rkind

    integer :: i,j
    real (kind=Rkind), allocatable :: d1(:),d2(:,:)


    QModel%ndim  = 2
    QModel%nsurf = 2

    allocate(QModel%Diab(QModel%nsurf,QModel%nsurf))

    !V11: 0.5*w1**2 * (X+a/2)**2 + 0.5d0*w2**2 * Y**2 + Delta/2
    i=1
    QModel%Diab(i,i)%Qref = [a/2,ZERO]
    d1 = [ZERO,ZERO]
    d2 = Identity_Mat(QModel%ndim)
    d2(1,1) = w1**2
    d2(2,2) = w2**2
    CALL set_dnS(QModel%Diab(i,i)%Ene,d0=Delta/2,d1=d1,d2=d2)

    !V22: 0.5*w1**2 * (X-a/2)**2 + 0.5d0*w2**2 * Y**2 - Delta/2
    i=2
    QModel%Diab(i,i)%Qref = [-a/2,ZERO]
    d1 = [ZERO,ZERO]
    d2 = Identity_Mat(QModel%ndim)
    d2(1,1) = w1**2
    d2(2,2) = w2**2
    CALL set_dnS(QModel%Diab(i,i)%Ene,d0=-Delta/2,d1=d1,d2=d2)

    ! V12: c * Y (linear coupling)
    i=2 ; j=1
    QModel%Diab(j,i)%Qref = [ZERO,ZERO]
    CALL set_dnS(QModel%Diab(j,i)%Ene,d0=ZERO,d1=[ZERO,c])
    ! V21
    QModel%Diab(i,j)%Ene  = QModel%Diab(j,i)%Ene
    QModel%Diab(i,j)%Qref = QModel%Diab(j,i)%Qref

  END SUBROUTINE Internal_QML_RJDI2014
  SUBROUTINE Internal_QML_7DM23_JGBL2023(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    USE ADdnSVM_m,        ONLY : dnS_t, set_dnS
    IMPLICIT NONE

    CLASS (QML_Vibronic_t), intent(inout)   :: QModel



    integer :: i,j,k
    real (kind=Rkind), allocatable :: d1(:),ks(:),d2(:,:)
    real (kind=Rkind), parameter     :: E1 = -1153.0899777201398_Rkind
    real (kind=Rkind), parameter     :: E2 = -1153.0794579197113_Rkind


    QModel%ndim  = 7
    QModel%nsurf = 2

    allocate(QModel%Diab(QModel%nsurf,QModel%nsurf))
    allocate(d1(QModel%ndim))
    allocate(ks(QModel%ndim))
    allocate(d2(QModel%ndim,QModel%ndim))

    ks(:) = [0.0000443070,0.0000568669,0.0000574396,0.0000577292,0.0001151293,0.0001028125,0.0001021209]
    !V11:
    i=1

    allocate(QModel%Diab(i,i)%Qref(QModel%ndim))
    QModel%Diab(i,i)%Qref(:) = ZERO

    d1(:) = [-0.0001922727,0.0003186424,-0.0001563349,0.0003059985,-0.0009354068,0.0002594575,-0.0001718730] ! gradient

    d2(:,:) = ZERO ! hessian
    DO k=1,QModel%ndim
      d2(k,k) = ks(k)
    END DO
 
    CALL set_dnS(QModel%Diab(i,i)%Ene,d0=E1,d1=d1,d2=d2)

    !V22: 
    i=2
    allocate(QModel%Diab(i,i)%Qref(QModel%ndim))
    QModel%Diab(i,i)%Qref(:) = ZERO

    d1(:) = [-0.0003068858,-0.0003058996,-0.0001800164,0.0000235402,-0.0000709817,-0.0002218948,-0.0009107390] ! gradient

    d2(:,:) = ZERO ! hessian
    DO k=1,QModel%ndim
      d2(k,k) = ks(k)
    END DO
    CALL set_dnS(QModel%Diab(i,i)%Ene,d0=E2,d1=d1,d2=d2)

    ! V12: 1 * Q(1) (linear coupling)
    i=2 ; j=1
    allocate(QModel%Diab(j,i)%Qref(QModel%ndim))
    QModel%Diab(j,i)%Qref(:) = ZERO

    !d1(:) = ZERO ! gradient
    !d1(1) = ONE
    d1(:) = [-0.0001124937,0.0001443638,0.0000963052,0.0000066299,-0.0000932210,0.0001905968,0.0001902977]


    CALL set_dnS(QModel%Diab(j,i)%Ene,d0=ZERO,d1=d1)

    ! V21
    QModel%Diab(i,j)%Ene  = QModel%Diab(j,i)%Ene
    QModel%Diab(i,j)%Qref = QModel%Diab(j,i)%Qref

    ! metric tensor
    QModel%d0GGdef      = Identity_Mat(QModel%ndim)

  END SUBROUTINE Internal_QML_7DM23_JGBL2023
  SUBROUTINE Internal_QML_SO2(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat, TO_string
    USE ADdnSVM_m,        ONLY : dnS_t, set_dnS, write_dnS
    IMPLICIT NONE

    CLASS (QML_Vibronic_t), intent(inout)   :: QModel

    real(kind=Rkind)  :: epsilon(13),SOC(13,13)
    real(kind=Rkind)  :: w(3)
    real(kind=Rkind)  :: kappa(3,13),lambda(3,13,13)


    integer :: i,j,ic
    real (kind=Rkind) :: d0
    real (kind=Rkind), allocatable :: d1(:),d2(:,:)

    epsilon(:) = [ZERO,0.1553277839_Rkind,0.1688652931_Rkind,0.3088375445_Rkind, &
                  0.1231069775_Rkind,0.1231069775_Rkind,0.1231069775_Rkind,      &
                  0.1547471701_Rkind,0.1547471701_Rkind,0.1547471701_Rkind,      &
                  0.1600894729_Rkind,0.1600894729_Rkind,0.1600894729_Rkind]
    w(:)    = [0.0023635993_Rkind, 0.0053089059_Rkind, 0.0064027455_Rkind]

    kappa(:,:)    = ZERO
    ! the 4 singlet contributions
    kappa(1:2, 1) = [-6.48289e-04_Rkind,-9.15365e-04_Rkind]
    kappa(1:2, 2) = [1.06541e-03_Rkind, 1.08185e-02_Rkind]
    kappa(1:2, 3) = [-7.92795e-03_Rkind, 1.63364e-02_Rkind]
    kappa(1:2, 4) = [-4.80026e-03_Rkind, 1.70393e-02_Rkind]
    ! the 3 singlet contributions (3 components)
    kappa(1:2, 5) = [ 2.60893e-03_Rkind, 6.95032e-03_Rkind]
    kappa(1:2, 6) = [ 2.60893e-03_Rkind, 6.95032e-03_Rkind]
    kappa(1:2, 7) = [ 2.60893e-03_Rkind, 6.95032e-03_Rkind]
    kappa(1:2, 8) = [-5.24793e-03_Rkind, 1.71110e-02_Rkind]
    kappa(1:2, 9) = [-5.24793e-03_Rkind, 1.71110e-02_Rkind]
    kappa(1:2,10) = [-5.24793e-03_Rkind, 1.71110e-02_Rkind]
    kappa(1:2,11) = [-7.98845e-03_Rkind, 1.62033e-02_Rkind]
    kappa(1:2,12) = [-7.98845e-03_Rkind, 1.62033e-02_Rkind]
    kappa(1:2,13) = [-7.98845e-03_Rkind, 1.62033e-02_Rkind]

    lambda(:,:,:) = ZERO
    ! 4 singlets
    lambda(3,1,4) = 1.65920e-02_Rkind
    lambda(3,4,1) = 1.65920e-02_Rkind
    lambda(3,2,3) = -5.55909e-03_Rkind
    lambda(3,3,2) = -5.55909e-03_Rkind
    lambda(3,1,3) = -3.66103e-03_Rkind
    lambda(3,3,1) = -3.66103e-03_Rkind

    ! triplets          5, 6, 7   8, 9,10  11,12,13
    ! corresponding ?: 1x,1y,1z  2x,2y,2z  3x,3y,3z
    lambda(3,5,11) = -3.66103e-03_Rkind ! x component
    lambda(3,6,12) = -3.66103e-03_Rkind ! y component
    lambda(3,7,13) = -3.66103e-03_Rkind ! z component
    lambda(3,11,5) = -3.66103e-03_Rkind ! x component
    lambda(3,12,6) = -3.66103e-03_Rkind ! y component
    lambda(3,13,7) = -3.66103e-03_Rkind ! z component
  

    SOC(:, 1) = [ZERO,ZERO,ZERO,ZERO,ZERO,   &
                -1.3669005e-07_Rkind,-2.0899909e-04_Rkind,-6.6636399e-04_Rkind,ZERO,1.3669005e-07_Rkind, &
                ZERO,-2.0503508e-06_Rkind,4.5563350e-08_Rkind]
    SOC(:, 2) = [ZERO,ZERO,ZERO,ZERO,ZERO,   &
                1.2270210e-04_Rkind,-4.5563350e-08_Rkind,ZERO,-9.1126700e-08_Rkind,-ZERO,&
                -9.1126700e-08_Rkind,-4.5563350e-08_Rkind,-2.5014279e-04_Rkind]
    SOC(:, 3) = [ZERO,ZERO,ZERO,ZERO,ZERO,   &
                4.5563350e-08_Rkind,ZERO,-4.5563350e-08_Rkind,-2.0070656e-04_Rkind,-1.3669005e-07_Rkind,&
                2.3451456e-04_Rkind,ZERO,4.1007015e-07_Rkind]
    SOC(:, 4) = [ZERO,ZERO,ZERO,ZERO,1.3669005e-07_Rkind,   &
                -1.0023937e-06_Rkind,-1.6612397e-04_Rkind,-5.7469053e-04_Rkind,2.7793644e-06_Rkind,9.5683035e-07_Rkind, &
                2.4148575e-06_Rkind,-2.0047874e-06_Rkind,2.2781675e-07_Rkind]

    SOC(:, 5) = [ZERO,ZERO,ZERO,1.3669005e-07_Rkind,ZERO,&
                ZERO,ZERO,ZERO,9.1126700e-08_Rkind,2.4080230e-04_Rkind,&
                ZERO,ZERO,ZERO]
    SOC(:, 6) = [-1.3669005e-07_Rkind,1.2270210e-04_Rkind,4.5563350e-08_Rkind,-1.0023937e-06_Rkind,ZERO, &
                ZERO,ZERO,-9.1126700e-08_Rkind,ZERO,ZERO,&
                ZERO,ZERO,2.1915971e-04_Rkind]
    SOC(:, 7) = [-2.0899909e-04_Rkind,-4.5563350e-08_Rkind,ZERO,-1.6612397e-04_Rkind,ZERO,&
                ZERO,ZERO,-2.4080230e-04_Rkind,ZERO,ZERO,&
                ZERO,-2.1915971e-04_Rkind,ZERO]
    SOC(:, 8) = [-6.6636399e-04_Rkind,ZERO,-4.5563350e-08_Rkind,-5.7469053e-04_Rkind,ZERO,&
                -9.1126700e-08_Rkind,-2.4080230e-04_Rkind,ZERO,ZERO,ZERO,  &
                ZERO,-9.3860501e-05_Rkind,-ZERO]
    SOC(:, 9) = [ZERO,-9.1126700e-08_Rkind,-2.0070656e-04_Rkind,2.7793644e-06_Rkind,9.1126700e-08_Rkind, &
                ZERO,ZERO,ZERO,ZERO,ZERO,&
                9.3860501e-05_Rkind,ZERO,2.2781675e-07_Rkind]
    SOC(:,10) = [ 1.3669005e-07_Rkind,-ZERO,-1.3669005e-07_Rkind,9.5683035e-07_Rkind,2.4080230e-04_Rkind,&
                  ZERO,ZERO,ZERO,ZERO,ZERO,&
                  ZERO,-2.2781675e-07_Rkind,ZERO]
    SOC(:,11) = [ZERO,-9.1126700e-08_Rkind,2.3451456e-04_Rkind,2.4148575e-06_Rkind,ZERO,&
                ZERO,ZERO,ZERO,9.3860501e-05_Rkind,ZERO, &
                ZERO,ZERO,ZERO]
    SOC(:,12) = [-2.0503508e-06_Rkind,-4.5563350e-08_Rkind,ZERO,-2.0047874e-06_Rkind,ZERO,&
                ZERO,-2.1915971e-04_Rkind,-9.3860501e-05_Rkind,ZERO,-2.2781675e-07_Rkind,&
                ZERO,ZERO,ZERO]
    SOC(:,13) = [4.5563350e-08_Rkind,-2.5014279e-04_Rkind,4.1007015e-07_Rkind,2.2781675e-07_Rkind,ZERO,&
                2.1915971e-04_Rkind,ZERO,-ZERO,2.2781675e-07_Rkind,ZERO,&
                ZERO,ZERO,ZERO]
  
  
    QModel%nsurf    = 13 ! 4 singlets and 3 triplets (3 components each)
    QModel%ndim     = 3

    allocate(QModel%Diab(QModel%nsurf,QModel%nsurf))

    DO i=1,QModel%nsurf
      QModel%Diab(i,i)%Qref = [ZERO,ZERO,ZERO]
      d0 = epsilon(i) + SOC(i,i)
      d1 = kappa(:,i)
      d2 = Identity_Mat(QModel%ndim)
      DO ic=1,QModel%ndim
        d2(ic,ic) = w(ic)
      END DO
      CALL set_dnS(QModel%Diab(i,i)%Ene,d0=d0,d1=d1,d2=d2)
      !CALL write_dnS(QModel%Diab(i,i)%Ene,info=('diabtic state #' // TO_string(i)))
    END DO

    DO i=1,QModel%nsurf
    DO j=1,QModel%nsurf
      QModel%Diab(j,i)%Qref = [ZERO,ZERO,ZERO]
      d0 = SOC(j,i)
      d1 = lambda(:,j,i)
      CALL set_dnS(QModel%Diab(j,i)%Ene,d0=d0,d1=d1)
    END DO
    END DO

  END SUBROUTINE Internal_QML_SO2
  SUBROUTINE Internal_QML_XXX(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    USE ADdnSVM_m,        ONLY : dnS_t, set_dnS
    IMPLICIT NONE

    CLASS (QML_Vibronic_t), intent(inout)   :: QModel



    integer :: i,j,k
    real (kind=Rkind), allocatable :: d1(:),d2(:,:)


    QModel%ndim  = 7
    QModel%nsurf = 2

    allocate(QModel%Diab(QModel%nsurf,QModel%nsurf))
    allocate(d1(QModel%ndim))
    allocate(d2(QModel%ndim,QModel%ndim))
    !V11:
    i=1

    allocate(QModel%Diab(i,i)%Qref(QModel%ndim))
    QModel%Diab(i,i)%Qref(:) = ZERO

    d1(:) = ZERO ! gradient

    d2(:,:) = ZERO ! hessian
    DO k=1,QModel%ndim
      d2(k,k) = ONE
    END DO
 
    CALL set_dnS(QModel%Diab(i,i)%Ene,d0=ZERO,d1=d1,d2=d2)

    !V22: 
    i=2
    allocate(QModel%Diab(i,i)%Qref(QModel%ndim))
    QModel%Diab(i,i)%Qref(:) = ZERO

    d1(:) = ZERO ! gradient

    d2(:,:) = ZERO ! hessian
    DO k=1,QModel%ndim
      d2(k,k) = ONE
    END DO
    CALL set_dnS(QModel%Diab(i,i)%Ene,d0=ONE,d1=d1,d2=d2)

    ! V12: 1 * Q(1) (linear coupling)
    i=2 ; j=1
    allocate(QModel%Diab(j,i)%Qref(QModel%ndim))
    QModel%Diab(j,i)%Qref(:) = ZERO

    !d1(:) = ZERO ! gradient
    !d1(1) = ONE
    d1(:) = [ONE,ZERO,ZERO,ZERO,ZERO,ZERO,ZERO]


    CALL set_dnS(QModel%Diab(j,i)%Ene,d0=ZERO,d1=d1)

    ! V21
    QModel%Diab(i,j)%Ene  = QModel%Diab(j,i)%Ene
    QModel%Diab(i,j)%Qref = QModel%Diab(j,i)%Qref

    QModel%d0GGdef      = Identity_Mat(QModel%ndim)

  END SUBROUTINE Internal_QML_XXX
!> @brief Subroutine wich prints the vibronic parameters.
!!
!> @author David Lauvergnat
!! @date 07/02/2023
!!
!! @param QModel        TYPE(QML_Vibronic_t):  derived type for the vibronic model.
!! @param nio           integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_Vibronic(QModel,nio)
    USE QDUtil_m,         ONLY : TO_string
    USE ADdnSVM_m,        ONLY : dnS_t,Write_dnS
    IMPLICIT NONE

    CLASS (QML_Vibronic_t), intent(in)   :: QModel
    integer,                intent(in)   :: nio

    integer :: i,j

    write(nio,*) 'Vibronics parameters'
    write(nio,*) '  Vibronic name: ',QModel%Vibronic_name
    write(nio,*) '  ndim:  ',QModel%ndim
    write(nio,*) '  nsurf: ',QModel%nsurf
    SELECT CASE (QModel%Vibronic_name)
    CASE ('twod_rjdi2014')
      write(nio,*) 'Published model from: '
      write(nio,*) ' Ilya G. Ryabinkin, Loïc Joubert-Doriol, and Artur F. Izmaylov, ...'
      write(nio,*) '  .... J. Chem. Phys. 140, 214116 (2014); https://doi.org/10.1063/1.4881147'
    CASE ('so2')
      write(nio,*) 'Published model from: '
      write(nio,*) ' Felix Plasser, Sandra Gómez, Maximilian F. S. J. Menger, Sebastian Mai and Leticia González, ...'
      write(nio,*) '  .... Phys. Chem. Chem. Phys., 2019, 21, 57--69; https://doi.org/10.1039/c8cp05662e'
    END SELECT
    !IF (print_level > -1) THEN
      DO i=1,QModel%nsurf
        write(nio,*) 'Qref_' // TO_String(i) // '-' // TO_String(i),QModel%Diab(i,i)%Qref
        CALL Write_dnS(QModel%Diab(i,i)%Ene,nio=nio,info='Diab_' // TO_String(i) // '-' // TO_String(i))
      END DO
      DO i=1,QModel%nsurf
      DO j=i+1,QModel%nsurf
        write(nio,*) 'Qref_' // TO_String(i) // '-' // TO_String(j),QModel%Diab(j,i)%Qref
        CALL Write_dnS(QModel%Diab(j,i)%Ene,nio=nio,info='Coupling_' // TO_String(j) // '-' // TO_String(i))
      END DO
      END DO
    !END IF
    write(nio,*) 'END Vibronics parameters'

  END SUBROUTINE Write_QML_Vibronic

!> @brief Subroutine wich calculates the Phenol potential with derivatives up to the 3d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      Potential with derivatives,.
!! @param Q                  real:                table of two values for which the potential is calculated (R,theta)
!! @param Para_Phenol        TYPE(Param_Phenol):  derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Vibronic(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS (QML_Vibronic_t), intent(in)     :: QModel
    TYPE (dnS_t),           intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),           intent(in)     :: dnQ(:)
    integer,                intent(in)     :: nderiv

    integer :: i,j
    TYPE (dnS_t),  allocatable    :: DeltaQ(:)

    DO i=1,QModel%nsurf
      DeltaQ = dnQ - QModel%Diab(i,i)%Qref
      !write(out_unit,*) i,i,'DeltaQ',get_d0(DeltaQ)
      Mat_OF_PotDia(i,i) = TO_Taylor(QModel%Diab(i,i)%Ene,DeltaQ)
      DO j=i+1,QModel%nsurf
        DeltaQ = dnQ - QModel%Diab(j,i)%Qref
        !write(out_unit,*) i,j,'DeltaQ',get_d0(DeltaQ)
        Mat_OF_PotDia(j,i) = TO_Taylor(QModel%Diab(j,i)%Ene,DeltaQ)
        Mat_OF_PotDia(i,j) = Mat_OF_PotDia(j,i)
      END DO
    END DO

  END SUBROUTINE EvalPot_QML_Vibronic


END MODULE QML_Vibronic_m
