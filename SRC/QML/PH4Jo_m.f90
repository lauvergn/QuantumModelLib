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
!> @brief Module which makes the initialization, calculation of the PH4 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 08/04/2021
!!
MODULE QML_PH4Jo_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE


  integer, parameter :: ndim      = 1 ! ????
  integer, parameter :: max_ndim  = 9 ! dimension of the system, here 9

  integer, parameter :: max_fit   = 1 + (max_ndim-1) + (max_ndim-1) + (max_ndim-1)**2 + 5
  integer, parameter :: max_nn    = 20
                                            !Qopt  1  2  3  4  5  6  7  8  9
  integer, parameter :: listQop_fit3(max_ndim) = [-1, 2, 3,-1, 5,-1, 7,-1,-1]

  character (len=*), parameter :: base_fit3_Ene1_fileName='InternalData/PH4/fit3/interEneMP2'
  character (len=*), parameter :: base_fit3_Ene2_fileName='InternalData/PH4/fit3/interEneCCSDT-F12'
  character (len=*), parameter :: base_fit3_Qopt_fileName='InternalData/PH4/fit3/interQ_'
  character (len=*), parameter :: base_fit3_grad_fileName='InternalData/PH4/fit3/interGrad_'
  character (len=*), parameter :: base_fit3_hess_fileName='InternalData/PH4/fit3/interHess_'
  character (len=*), parameter :: base_fit3_AnHar_fileName='InternalData/PH4/fit3/inter_anharCCSDT-F12_'

!> @brief Derived type in which the PH4 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_PH4Jo_t

    PRIVATE

    real (kind=Rkind)  :: F(max_nn,max_fit)   = ZERO
    integer            :: nn(0:ndim,max_fit)  = 0
    integer            :: nt(max_fit)         = 0
    integer            :: largest_nn          = 0
    integer            :: ifunc_Ene           = 0
    integer            :: ifunc_Qopt          = 0
    integer            :: ifunc_Grad          = 0
    integer            :: ifunc_Hess          = 0
    integer            :: ifunc_AnHar         = 0


    real (kind=Rkind)  :: a(max_fit)          = ZERO
    real (kind=Rkind)  :: b(max_fit)          = ZERO

    logical            :: file_exist(max_fit) = .FALSE.

    integer, allocatable :: iQopt_TO_ifunc(:)
    integer, allocatable :: iQgrad_TO_ifunc(:)
    integer, allocatable :: iQjQHess_TO_ifunc(:,:)
    integer, allocatable :: AnHar_TO_ifunc(:)

   CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_PH4Jo
    PROCEDURE :: EvalFunc_QModel  => EvalFunc_QML_PH4Jo
    PROCEDURE :: Write_QModel     => Write_QML_PH4Jo
  END TYPE QML_PH4Jo_t

  PUBLIC :: QML_PH4Jo_t,Init_QML_PH4Jo

  CONTAINS
!> @brief Function which makes the initialization of the PH4 parameters.
!!
!! @param QModel             TYPE(QML_PH4Jo_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_PH4Jo(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat, TO_string, file_open2
    USE QMLLib_UtilLib_m, ONLY : make_QMLInternalFileName
    IMPLICIT NONE

    TYPE (QML_PH4Jo_t)                             :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    integer :: i,j,ifunc
    logical :: exist,read_ab
    character (len=:), allocatable  :: FileName

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_PH4Jo'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    IF (debug) write(out_unit,*) 'option',QModel%option

    IF (QModel%ndim == 0) THEN
      ! it means that ndim was not present in the CALL Init_Model().
      ! => ndim is set to the default value (1)
      QModel%ndim  = 1
    END IF

    !The value of QModel%ndim must be 1, 2 or 9
    IF (QModel%ndim /= 1 .AND. QModel%ndim /= 9 .AND. QModel%ndim /= 2) THEN
       write(out_unit,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unit)
       write(out_unit,*) ' ERROR in ',name_sub
       write(out_unit,*) ' ndim MUST equal to 1 or 2 or 9. ndim: ',QModel%ndim
       STOP 'ERROR in Init_QML_PH4Jo: ndim MUST equal to 1 or 2 or 9'
    END IF

    IF (debug) write(out_unit,*) 'init d0GGdef of PH4'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

    SELECT CASE (QModel%option)
    CASE (7)
        QModel%nsurf    = 1
        QModel%ndim     = 1
        QModel%pot_name = 'ph4_1D'
        QModel%ndimFunc = 1
        QModel%nb_Func  = 53
    CASE (8,18,80,81,82,83,84)
        QModel%nsurf    = 1
        QModel%ndim     = 2
        QModel%pot_name = 'ph4_2D'
        QModel%ndimFunc = 1
        QModel%nb_Func  = 53
    CASE (6,9)
        QModel%nsurf    = 1
        QModel%ndim     = 9
        QModel%pot_name = 'ph4_9D'
        QModel%ndimFunc = 1
        QModel%nb_Func   = 53
    CASE default
        STOP 'STOP in Init_QML_PH4Jo: this option is not mpossible.'
    END SELECT

  END FUNCTION Init_QML_PH4Jo
!> @brief Subroutine wich prints the current QML_PH4 parameters.
!!
!! @param QModel            CLASS(QML_PH4Jo_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_PH4Jo(QModel,nio)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_PH4Jo_t),     intent(in) :: QModel
    integer,              intent(in) :: nio

    TYPE (dnS_t)    :: sN,rhoN
    TYPE (dnS_t)    :: sO,rhoO

    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) 'Model potential of PH4'
    write(nio,*)
    write(nio,*) ' Units: Energy in Hartree and distances in bohr'
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*)
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) '|||||||||||||||| option ',QModel%option,' of PH4 ||||||||||||||||'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) '-------------------------------------------------------------'
    write(nio,*) '----------------- 1D-function order -------------------------'
    write(nio,*) 'index:  1    => V0(s)'
    write(nio,*) 'index:  2-9  => Q^i_opt(s)'
    write(nio,*) 'index: 10-17 => Grad^i_opt(s)'
    write(nio,*) 'index: 18-25 => Hess^ii(s)'
    write(nio,*) 'index: 26-53 => Hess^ij(s) (j>i)'
    write(nio,*) '-------------------------------------------------------------'
 
    SELECT CASE(QModel%option)
    CASE (6)
        write(nio,*) '-------------------- 9D model -------------------------------'
        write(nio,*) '-------------------- V=V0+grad(i)*dqi+1/2Hess(i,i)*dqi*dqi --'
    CASE (7)
        write(nio,*) '-------------------- 1D model: s ----------------------------'
        write(nio,*) '-------------------- V=V0 -----------------------------------'
    CASE (8)
        write(nio,*) '-------------------- 2D model: s,rho ------------------------'
        write(nio,*) '------------- V=V0+grad(rho)*drho+1/2Hess(rho,rho)*drho**2 --'
    CASE (18)
        write(nio,*) '-------------------- 2D model: s,rho ------------------------'
        write(nio,*) '------------- V=V0+grad(rho)*drho+1/2Hess(rho,rho)*drho**2 --'
        write(nio,*) '------------- with new coordinates (add a transfo) ----------'
      sO   = ZERO
      rhoO = 0.89_Rkind
      CALL srho_old_TO_srho_new(sO,rhoO,sN,rhoN)
      write(nio,*) 'coord old => new',get_d0(sO),get_d0(rhoO),get_d0(sN),get_d0(rhoN)
      CALL srho_new_TO_srho_old(sN,rhoN,sO,rhoO)
      write(nio,*) 'coord new => old',get_d0(sN),get_d0(rhoN),get_d0(sO),get_d0(rhoO)
      write(nio,*)
    CASE (80)
        write(nio,*) '-------------------- 2D model: s,rho ------------------------'
        write(nio,*) '-------------------- V=V0 (without hessian contributions) ---'
    CASE (81)
        write(nio,*) '-------------------- 2D model: s,rho ------------------------'
        write(nio,*) '-------------------- V=V0 + rho_bottel(s) + hess_bottle(s) --'
    CASE (82)
        write(nio,*) '-------------------- 2D model: s,rho ------------------------'
        write(nio,*) '-------------------- V=V0 + rho_opt_PH4 + hess_bottle -------'
    CASE (83)
        write(nio,*) '-------------------- 2D model: s,rho ------------------------'
        write(nio,*) '-------------------- V=V0 + rho_opt_bottle + hess_PH4 -------'
    CASE (84)
        write(nio,*) '-------------------- 2D model: s,rho ------------------------'
        write(nio,*) '-------------------- V=V0_bottle + rho_opt_bottle + hess_bottle'
    CASE (9)
        write(nio,*) '-------------------- 9D model -------------------------------'
        write(nio,*) '-------------------- V=V0+grad(i)*dqi+Hess(i,j)*dqi*dqj -----'

    CASE DEFAULT
        write(nio,*) '------------- WARNING: No option default --------------------'
    END SELECT
    write(nio,*)

  END SUBROUTINE Write_QML_PH4Jo
!> @brief Subroutine wich calculates the PH4 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_PH4Jo_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_PH4Jo(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_PH4Jo_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv


    integer, parameter          :: nonAct =1, nonGrad =13, nonHess =13, nqact =8
    TYPE (dnS_t), allocatable   :: dnPoly(:)
    TYPE (dnS_t)    :: Func(QModel%nb_Func), Gradient(QModel%nb_Func-nonGrad)
    TYPE (dnS_t)    :: Hess(QModel%nb_Func-nonHess), miHess(nqact, nqact)
    TYPE (dnS_t)    :: dnDQ(QModel%ndim), vh,Rm,tRm, V0, V1, Somme, Somme2, Dq(QModel%ndim-nonAct)
    TYPE (dnS_t)    :: dnQFunc(QModel%ndimFunc)
    integer         :: i1,i2,i,ifunc,j,k

    TYPE (dnS_t)    :: sO,rhoO

    ! for the bottleneck model
    real(kind=Rkind), parameter :: alpha    = 7.3_Rkind
    real(kind=Rkind), parameter :: V0bottle = 0.00551239856_Rkind
    real(kind=Rkind), parameter :: b        = 0.1_Rkind
    real(kind=Rkind), parameter :: a        = 1.5_Rkind ! bohr^-1
    real(kind=Rkind), parameter :: k0       = 1783.31376308_Rkind*0.01_Rkind**2
    real(kind=Rkind), parameter :: qref     = 1._Rkind

    !write(out_unit,*) 'coucou pot PH4 '

    SELECT CASE(QModel%option)
    CASE(6) ! Model with V0 +gradQi+HessQiQi  
      dnQFunc = [dnQ(2)] ! s
      CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
      V0 =Func(1)
      Gradient = ZERO
      Hess     = ZERO
      somme    = ZERO
      somme2   = ZERO
      Dq(:) = ZERO
      Dq(1) = (dnQ(1)-Func(2))
      Dq(2) = (dnQ(3)-Func(3))
      Dq(3) = (dnQ(4)-Func(4))
      Dq(4) = (dnQ(5)-Func(5))
      Dq(5) = (dnQ(6)-Func(6))
      Dq(6) = (dnQ(7)-Func(7))
      Dq(7) = (dnQ(8)-Func(8))
      Dq(8) = (dnQ(9)-Func(9))
      DO j=10,17
         Gradient(j-9) = Func(j)
      END DO
      DO k=18,25
         Hess(k-17) = Func(k)
      END DO
      DO i=1,8
         somme = somme + (Dq(i))*Gradient(i)
         somme2 = somme2 +  HALF*(Dq(i)*Dq(i))*Hess(i)  
      END DO
      Mat_OF_PotDia(1,1) = V0 + somme + somme2

    CASE(7) ! Model with V0 only  
      dnQFunc = [dnQ(1)] ! s ???
      CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
      V0 =Func(1)

      Mat_OF_PotDia(1,1) = V0 

    CASE(80) ! 2D-Model with V0
      dnQFunc = [dnQ(2)] ! s
      CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
      V0 =Func(1)
      Mat_OF_PotDia(1,1) = V0

    CASE(81) ! Model with V0 + rho_opt_bottle + hess_bottle  
        dnQFunc = [dnQ(2)] ! s
        CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
        V0 =Func(1)
        !Dq(1)   = (dnQ(1)-Func(2))
        !Hess(1) = Func(18)

        Dq(1)   = (dnQ(1)-qref)
        hess(1) = k0*(ONE + b / cosh(a * dnQ(2))**2 )

        Mat_OF_PotDia(1,1) = V0 + HALF*hess(1)*Dq(1)**2

    CASE(82) ! Model with V0 + rho_opt_PH4 + hess_bottle  
         dnQFunc = [dnQ(2)] ! s
        CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
        V0 =Func(1)
        Dq(1)   = (dnQ(1)-Func(2))
        !Hess(1) = Func(18)

        !Dq(1)   = (dnQ(1)-qref)
        hess(1) = k0*(ONE + b / cosh(a * dnQ(2))**2 )

        Mat_OF_PotDia(1,1) = V0 + HALF*hess(1)*Dq(1)**2
    CASE(83) ! Model with V0 + rho_opt_bottle + hess_PH4
        dnQFunc = [dnQ(2)] ! s
       CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
       V0 =Func(1)
       !Dq(1)   = (dnQ(1)-Func(2))
       Hess(1) = Func(18)

       Dq(1)   = (dnQ(1)-qref)
       !hess(1) = k0*(ONE + b / cosh(a * dnQ(2))**2 )

       Mat_OF_PotDia(1,1) = V0 + HALF*hess(1)*Dq(1)**2
    CASE(84) ! Model with V0_bottle + rho_opt_bottle + hess_bottle
       V0      = V0bottle * ( (1-alpha)/(1+exp(-TWO*a*dnQ(2)))+ (HALF*(1+sqrt(alpha))/cosh(a * dnQ(2)))**2)
       Dq(1)   = (dnQ(1)-qref)
       hess(1) = k0*(ONE + b / cosh(a * dnQ(2))**2 )

       Mat_OF_PotDia(1,1) = V0 + HALF*hess(1)*Dq(1)**2
    CASE(8) ! Model with V0 + grad(rhoOpt)*DeltarhoOpt + HALF*Hess(rhoOpt)*DeltarhoOpt**2   
      dnQFunc = [dnQ(2)] ! s
      CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
      V0 =Func(1)
      Gradient = ZERO
      Hess = ZERO
      somme = ZERO
      somme2 = ZERO
      Dq(:) = ZERO
      Dq(1) = (dnQ(1)-Func(2))

      DO j=10,17
         Gradient(j-9) = Func(j)
      END DO
      DO k=18,25
         Hess(k-17) = Func(k)
      END DO
      somme = Dq(1)*Gradient(1) + HALF*(Dq(1)*Dq(1))*Hess(1)  

      Mat_OF_PotDia(1,1) = V0 + somme 
    CASE(18) ! Model with V0 + grad(rhoOpt)*DeltarhoOpt + HALF*Hess(rhoOpt)*DeltarhoOpt**2   
        CALL srho_new_TO_srho_old(dnQ(2),dnQ(1),sO,rhoO)
        dnQFunc = [sO] ! s
        CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
        V0      = Func(1)
        Dq(1)   = (rhoO-Func(2))
        Gradient(1) = Func(10)
        Hess(1) = Func(18)
  
        Mat_OF_PotDia(1,1) = V0 + Dq(1)*Gradient(1) + HALF*(Dq(1)*Dq(1))*Hess(1) 
    CASE(9) ! Model with V0 +grad*di+Hess*dQi*dQj  
      dnQFunc = [dnQ(2)] ! s
      CALL EvalFunc_QML_PH4Jo(QModel,Func,dnQFunc,nderiv)
      V0 =Func(1)
      V1 =ZERO
      Gradient = ZERO
      miHess(:,:) = ZERO
      somme = ZERO
      somme2 = ZERO
      Dq(:) = ZERO
      Dq(1) = (dnQ(1)-Func(2))
      Dq(2) = (dnQ(3)-Func(3))
      Dq(3) = (dnQ(4)-Func(4))
      Dq(4) = (dnQ(5)-Func(5))
      Dq(5) = (dnQ(6)-Func(6))
      Dq(6) = (dnQ(7)-Func(7))
      Dq(7) = (dnQ(8)-Func(8))
      Dq(8) = (dnQ(9)-Func(9))

      DO j=10,17
         Gradient(j-9) = Func(j)
      END DO
      DO k=18,25
         miHess(k-17,k-17) = Func(k)
      END DO

!############ Filling the half offdiagonal elements of the hessian the others are equal to 0 ############
!######### The coresponding coordinates are in the order : rho, Rph, z(0), A, zz(0), D, tAH2, pi #########
!################################### cf Jo's M2 internship notebook ###################################

      miHess(2,1) = Func(53)
      
      miHess(3,1) = Func(52); miHess(3,2) = Func(51)
      
      miHess(4,1) = Func(49); miHess(4,2) = Func(48)
      miHess(4,3) = Func(50); miHess(5,1) = Func(47); miHess(5,2) = Func(46); miHess(5,3) = Func(45); miHess(5,4) = Func(44)
      
      miHess(6,1) = Func(43); miHess(6,2) = Func(40); miHess(6,3) = Func(37); miHess(6,4) = Func(34); miHess(6,5) = Func(31)
      
      miHess(7,1) = Func(42); miHess(7,2) = Func(39); miHess(7,3) = Func(36); miHess(7,4) = Func(33); miHess(7,5) = Func(30)
      miHess(7,6) = Func(31);
      
      miHess(8,1) = Func(41); miHess(8,2) = Func(38); miHess(8,3) = Func(35); miHess(8,4) = Func(32); miHess(8,5) = Func(29);
      miHess(8,6) = Func(27); miHess(8,7) = Func(26)

      DO i=1,8
         DO j=1,8
            IF(j <= i) THEN
               V1 = V1  + (Dq(i))*Gradient(i) +  HALF*(Dq(i)*Dq(j))*miHess(i,j)
            END IF
         END DO  
      END DO
      Mat_OF_PotDia(1,1) = V0 + V1
    CASE default
        STOP 'ERROR in EvalPot_QML_PH4Jo: wrong option'
    END SELECT
    CALL dealloc_dnS(dnDQ)

  END SUBROUTINE EvalPot_QML_PH4Jo

  SUBROUTINE EvalFunc_QML_PH4Jo(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_PH4Jo_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    TYPE (dnS_t)    :: s, ts 
    integer, parameter     :: max_deg = 52
    TYPE (dnS_t)           :: tab_Pl(0:max_deg)
    integer         :: i,i1,i2,ifunc
  
    s= dnQ(1)
    !write(out_unit,*) 's EvalFunc =', get_d0(dnQ(1))
    ts= sin(atan(s))
    DO i=0,max_deg
        tab_Pl(i) = dnLegendre0(ts,i,ReNorm=.FALSE.)
    END DO

!V0, grille ts = -0.99,-0.98...0.99, maxdiff= 1.87537231682e-05, npoly=20
    Func(1) = &
    (-343.1648425550447_Rkind) * tab_Pl(0) +  &
(-0.02360558201626247_Rkind) * tab_Pl(1) +  &
(-0.009808350981168643_Rkind) * tab_Pl(2) +  &
(0.008228731540991713_Rkind) * tab_Pl(3) +  &
(0.00017440093734658762_Rkind) * tab_Pl(4) +  &
(0.0008980191074552608_Rkind) * tab_Pl(5) +  &
(0.000862562810530883_Rkind) * tab_Pl(6) +  &
(-0.0005829711697572901_Rkind) * tab_Pl(7) +  &
(0.0003657406988670201_Rkind) * tab_Pl(8) +  &
(-0.0006014390836917339_Rkind) * tab_Pl(9) +  &
(0.00018547138077772536_Rkind) * tab_Pl(10) +  &
(-0.0002472615890990517_Rkind) * tab_Pl(11) +  &
(4.1364788438304396e-05_Rkind) * tab_Pl(12) +  &
(-3.454881279697795e-05_Rkind) * tab_Pl(13) +  &
(-2.585447869911675e-05_Rkind) * tab_Pl(14) +  &
(3.762768966799136e-05_Rkind) * tab_Pl(15) +  &
(-4.3720678971352355e-05_Rkind) * tab_Pl(16) +  &
(4.361267600135664e-05_Rkind) * tab_Pl(17) +  &
(-3.296858570874348e-05_Rkind) * tab_Pl(18) +  &
(2.978968580294176e-05_Rkind) * tab_Pl(19) +  &
(-1.6513395295385247e-05_Rkind) * tab_Pl(20)
 
!################# OPTIMAL VALUES OF THE COORDINATES ON THE MEP #################

!rho_opt, grille ts = -0.99,-0.98...0.99, maxdiff= 2.39125686359e-05, npoly=20
    Func(2) = &
(0.8894713089773137_Rkind) * tab_Pl(0) +  &
(-0.003924888793412194_Rkind) * tab_Pl(1) +  &
(0.06058955355753642_Rkind) * tab_Pl(2) +  &
(0.0013818962617743753_Rkind) * tab_Pl(3) +  &
(0.03391272639351013_Rkind) * tab_Pl(4) +  &
(0.005043797287396934_Rkind) * tab_Pl(5) +  &
(0.010416568945493613_Rkind) * tab_Pl(6) +  &
(-0.0009757158681710542_Rkind) * tab_Pl(7) +  &
(0.0038195349569326146_Rkind) * tab_Pl(8) +  &
(-0.001371170516489801_Rkind) * tab_Pl(9) +  &
(0.0011299691093062517_Rkind) * tab_Pl(10) +  &
(-0.0005190388330162536_Rkind) * tab_Pl(11) +  &
(0.0004091308244872469_Rkind) * tab_Pl(12) +  &
(-1.2480237888512805e-05_Rkind) * tab_Pl(13) +  &
(0.00013881683410524206_Rkind) * tab_Pl(14) +  &
(0.00015092366204668198_Rkind) * tab_Pl(15) +  &
(8.680212399877544e-07_Rkind) * tab_Pl(16) +  &
(0.00010122005011163127_Rkind) * tab_Pl(17) +  &
(-3.9416353715322884e-05_Rkind) * tab_Pl(18) +  &
(4.899923260962432e-05_Rkind) * tab_Pl(19) +  &
(-2.626504792056346e-05_Rkind) * tab_Pl(20)

!Rph_opt, grille ts = -0.99,-0.98...0.99, maxdiff= 3.49860764475e-05, npoly=20
    Func(3) = &
(2.6716813763435274_Rkind) * tab_Pl(0) +  &
(0.002893059245846_Rkind) * tab_Pl(1) +  &
(-0.0008283749005720888_Rkind) * tab_Pl(2) +  &
(0.0006317091510585448_Rkind) * tab_Pl(3) +  &
(-0.00026466010307466307_Rkind) * tab_Pl(4) +  &
(-0.0002610551400225393_Rkind) * tab_Pl(5) +  &
(2.0613010447660722e-05_Rkind) * tab_Pl(6) +  &
(-0.00010914759261107084_Rkind) * tab_Pl(7) +  &
(6.914001126858114e-05_Rkind) * tab_Pl(8) +  &
(-3.768606261289415e-05_Rkind) * tab_Pl(9) +  &
(4.8556448271212604e-05_Rkind) * tab_Pl(10) +  &
(-3.3227139665329485e-06_Rkind) * tab_Pl(11) +  &
(1.3217159388108674e-05_Rkind) * tab_Pl(12) +  &
(1.101402649448329e-05_Rkind) * tab_Pl(13) +  &
(-1.0240577887012914e-05_Rkind) * tab_Pl(14) +  &
(6.7392926633381256e-06_Rkind) * tab_Pl(15) +  &
(-9.954739418150701e-06_Rkind) * tab_Pl(16) +  &
(4.331997684645042e-06_Rkind) * tab_Pl(17) +  &
(-9.585678823293377e-06_Rkind) * tab_Pl(18) +  &
(8.226546079795832e-08_Rkind) * tab_Pl(19) +  &
(-8.743243535823865e-06_Rkind) * tab_Pl(20)

    Func(4) = 0.0_Rkind ! sym

!A_opt, grille ts = -0.99,-0.98...0.99, maxdiff= 4.29436342259e-05, npoly=20
    Func(5) = &
(0.8077272635307288_Rkind) * tab_Pl(0) +  &
(-0.008352588424887657_Rkind) * tab_Pl(1) +  &
(0.0029892616010482233_Rkind) * tab_Pl(2) +  &
(0.0005399789340536514_Rkind) * tab_Pl(3) +  &
(-0.0006817579157627922_Rkind) * tab_Pl(4) +  &
(0.0006132559875251533_Rkind) * tab_Pl(5) +  &
(2.2203888135675074e-05_Rkind) * tab_Pl(6) +  &
(-7.25194165215201e-05_Rkind) * tab_Pl(7) +  &
(-8.400200984198919e-05_Rkind) * tab_Pl(8) +  &
(1.2380081332694868e-05_Rkind) * tab_Pl(9) +  &
(-3.9941561647210174e-05_Rkind) * tab_Pl(10) +  &
(-2.5703632918380975e-05_Rkind) * tab_Pl(11) +  &
(1.2548627772949678e-05_Rkind) * tab_Pl(12) +  &
(-2.249364016854002e-05_Rkind) * tab_Pl(13) +  &
(1.1811199947951431e-05_Rkind) * tab_Pl(14) +  &
(2.8808074238000615e-07_Rkind) * tab_Pl(15) +  &
(1.1269624718307157e-05_Rkind) * tab_Pl(16) +  &
(-1.5597529899306575e-06_Rkind) * tab_Pl(17) +  &
(1.0109815978324593e-05_Rkind) * tab_Pl(18) +  &
(6.311195514938651e-06_Rkind) * tab_Pl(19) +  &
(1.1554241035030076e-06_Rkind) * tab_Pl(20)

    Func(6) = 0.0_Rkind ! sym

!D_opt, grille ts = -0.99,-0.98...0.89-->1.0, maxdiff= 0.00145097187801, npoly=30, ! dernier pts placé à D=1.12 (arbitraire) 
    Func(7) = &
(1.533767968893748_Rkind) * tab_Pl(0) +  &
(-0.08180754815100122_Rkind) * tab_Pl(1) +  &
(0.20519988645062265_Rkind) * tab_Pl(2) +  &
(0.3220289831847085_Rkind) * tab_Pl(3) +  &
(0.35827940517077694_Rkind) * tab_Pl(4) +  &
(0.3908760118236028_Rkind) * tab_Pl(5) +  &
(0.3857257247346453_Rkind) * tab_Pl(6) +  &
(0.341953135090798_Rkind) * tab_Pl(7) +  &
(0.26501315300326916_Rkind) * tab_Pl(8) +  &
(0.1709165799208731_Rkind) * tab_Pl(9) +  &
(0.0640555857904282_Rkind) * tab_Pl(10) +  &
(-0.04302654670232631_Rkind) * tab_Pl(11) +  &
(-0.13994464772554482_Rkind) * tab_Pl(12) +  &
(-0.21964585527385916_Rkind) * tab_Pl(13) +  &
(-0.2763813223049032_Rkind) * tab_Pl(14) +  &
(-0.30925148998718743_Rkind) * tab_Pl(15) +  &
(-0.3183890541549751_Rkind) * tab_Pl(16) +  &
(-0.30749822606034394_Rkind) * tab_Pl(17) +  &
(-0.28092551803505084_Rkind) * tab_Pl(18) +  &
(-0.24394707001100002_Rkind) * tab_Pl(19) +  &
(-0.2019272861272134_Rkind) * tab_Pl(20) +  &
(-0.15917185524411118_Rkind) * tab_Pl(21) +  &
(-0.11955045688213226_Rkind) * tab_Pl(22) +  &
(-0.08495502959639412_Rkind) * tab_Pl(23) +  &
(-0.05706964246177972_Rkind) * tab_Pl(24) +  &
(-0.035721315376102986_Rkind) * tab_Pl(25) +  &
(-0.020681747286873678_Rkind) * tab_Pl(26) +  &
(-0.01076795855216272_Rkind) * tab_Pl(27) +  &
(-0.004917417430352599_Rkind) * tab_Pl(28) +  &
(-0.001807456470637829_Rkind) * tab_Pl(29) +  &
(-0.0004289935921470803_Rkind) * tab_Pl(30)

    Func(8) = 0.0_Rkind !tAH2   
    Func(9) = 3.1415926535897931_Rkind ! pi

!############################ Gradients ######################

!gradro, grille ts = -0.99,-0.98...0.99, maxdiff= 3.41837809738e-05, npoly=20
    Func(10) = &
(-3.45015714992917e-06_Rkind) * tab_Pl(0) +  &
(3.7090666717660546e-06_Rkind) * tab_Pl(1) +  &
(-6.319399726499579e-06_Rkind) * tab_Pl(2) +  &
(3.0253312513887784e-06_Rkind) * tab_Pl(3) +  &
(-3.540035070842396e-06_Rkind) * tab_Pl(4) +  &
(-6.812207534332742e-06_Rkind) * tab_Pl(5) +  &
(-8.119503585021632e-06_Rkind) * tab_Pl(6) +  &
(-1.522831061681485e-05_Rkind) * tab_Pl(7) +  &
(1.0676885602931986e-05_Rkind) * tab_Pl(8) +  &
(1.6301320261820365e-06_Rkind) * tab_Pl(9) +  &
(1.6907983367590193e-06_Rkind) * tab_Pl(10) +  &
(-6.191515673259921e-06_Rkind) * tab_Pl(11) +  &
(-1.2882882693569997e-05_Rkind) * tab_Pl(12) +  &
(-2.201722295804488e-06_Rkind) * tab_Pl(13) +  &
(4.927173900647309e-06_Rkind) * tab_Pl(14) +  &
(1.28732256877055e-05_Rkind) * tab_Pl(15) +  &
(-9.134104280480697e-06_Rkind) * tab_Pl(16) +  &
(-8.663035801061593e-06_Rkind) * tab_Pl(17) +  &
(-7.054309854615688e-08_Rkind) * tab_Pl(18) +  &
(-2.678572316084898e-06_Rkind) * tab_Pl(19) +  &
(1.3246185725023105e-07_Rkind) * tab_Pl(20)

!gradRph, grille ts = -0.99,-0.98...0.99, maxdiff= 1.47292375807e-05, npoly=20
    Func(11) = &
(1.9622838127584655e-06_Rkind) * tab_Pl(0) +  &
(3.951486614681554e-06_Rkind) * tab_Pl(1) +  &
(-3.5596072659631525e-06_Rkind) * tab_Pl(2) +  &
(-1.1912929991519755e-06_Rkind) * tab_Pl(3) +  &
(-1.4481522792466443e-06_Rkind) * tab_Pl(4) +  &
(-4.368735425640163e-07_Rkind) * tab_Pl(5) +  &
(5.826657826674407e-06_Rkind) * tab_Pl(6) +  &
(-2.5657819846798235e-06_Rkind) * tab_Pl(7) +  &
(6.402784901440222e-07_Rkind) * tab_Pl(8) +  &
(3.156703863059984e-08_Rkind) * tab_Pl(9) +  &
(-2.592648699988574e-07_Rkind) * tab_Pl(10) +  &
(2.4892208574906205e-06_Rkind) * tab_Pl(11) +  &
(-1.0135828513542424e-06_Rkind) * tab_Pl(12) +  &
(2.6075928064966495e-06_Rkind) * tab_Pl(13) +  &
(-5.343960602152863e-06_Rkind) * tab_Pl(14) +  &
(7.493943563990112e-07_Rkind) * tab_Pl(15) +  &
(-8.024063570020884e-07_Rkind) * tab_Pl(16) +  &
(1.9404518078781234e-06_Rkind) * tab_Pl(17) +  &
(1.4566267805327955e-06_Rkind) * tab_Pl(18) +  &
(2.1024798499887706e-06_Rkind) * tab_Pl(19) +  &
(5.379407334129039e-07_Rkind) * tab_Pl(20)


!gradz, grille ts = -0.99,-0.98...0.99, maxdiff=0.0, il est nul, npoly=20
    Func(12) = &
(0.0_Rkind) * tab_Pl(0)

!gradA, grille ts = -0.99,-0.98...0.99, maxdiff= 3.06971478016e-05, npoly=20
    Func(13) = &
(2.5828232167454165e-06_Rkind) * tab_Pl(0) +  &
(-6.567648941827709e-06_Rkind) * tab_Pl(1) +  &
(-1.4344846623269361e-06_Rkind) * tab_Pl(2) +  &
(1.1285431890430704e-05_Rkind) * tab_Pl(3) +  &
(-6.108240161703317e-06_Rkind) * tab_Pl(4) +  &
(-5.553988787112458e-06_Rkind) * tab_Pl(5) +  &
(2.387901675234784e-06_Rkind) * tab_Pl(6) +  &
(1.0965774322336397e-07_Rkind) * tab_Pl(7) +  &
(5.315841513323174e-06_Rkind) * tab_Pl(8) +  &
(-2.495657466349989e-06_Rkind) * tab_Pl(9) +  &
(-3.5132554794631532e-06_Rkind) * tab_Pl(10) +  &
(1.876402530383104e-06_Rkind) * tab_Pl(11) +  &
(1.3034170295953487e-06_Rkind) * tab_Pl(12) +  &
(-2.975032515700928e-06_Rkind) * tab_Pl(13) +  &
(-8.90746782654235e-07_Rkind) * tab_Pl(14) +  &
(5.792899657383682e-06_Rkind) * tab_Pl(15) +  &
(-2.532080740967144e-06_Rkind) * tab_Pl(16) +  &
(-2.4498552180850536e-06_Rkind) * tab_Pl(17) +  &
(-1.115850491857519e-06_Rkind) * tab_Pl(18) +  &
(6.218935006248005e-07_Rkind) * tab_Pl(19) +  &
(-3.4835220118801224e-06_Rkind) * tab_Pl(20)


!gradzz, grille ts = -0.99,-0.98...0.99, maxdiff=0.0, il est nul, npoly=20
    Func(14) = &
(0.0_Rkind) * tab_Pl(0)

!gradD, grille ts = -0.99,-0.98...0.99, maxdiff= 1.19673106115e-05, npoly=20
    Func(15) = &
(-1.632563180200473e-06_Rkind) * tab_Pl(0) +  &
(-7.025219565883e-07_Rkind) * tab_Pl(1) +  &
(-2.5593270531392556e-06_Rkind) * tab_Pl(2) +  &
(-3.7094249761884745e-06_Rkind) * tab_Pl(3) +  &
(6.530885489348272e-08_Rkind) * tab_Pl(4) +  &
(-5.498959454443389e-06_Rkind) * tab_Pl(5) +  &
(2.3344418400899365e-06_Rkind) * tab_Pl(6) +  &
(-2.937163533320034e-06_Rkind) * tab_Pl(7) +  &
(6.313707416381471e-06_Rkind) * tab_Pl(8) +  &
(-2.633147856553801e-06_Rkind) * tab_Pl(9) +  &
(1.0454725790837134e-05_Rkind) * tab_Pl(10) +  &
(-1.3531540866644168e-06_Rkind) * tab_Pl(11) +  &
(8.115570549551766e-06_Rkind) * tab_Pl(12) +  &
(2.224988803435241e-07_Rkind) * tab_Pl(13) +  &
(4.892637896323156e-06_Rkind) * tab_Pl(14) +  &
(-3.2794238336879473e-07_Rkind) * tab_Pl(15) +  &
(7.423359966486157e-07_Rkind) * tab_Pl(16) +  &
(-1.7084200892002648e-06_Rkind) * tab_Pl(17) +  &
(-4.512311956065796e-06_Rkind) * tab_Pl(18) +  &
(-1.8963041049498412e-06_Rkind) * tab_Pl(19) +  &
(-3.898371342679767e-06_Rkind) * tab_Pl(20)



!gradtAH2, grille ts = -0.99,-0.98...0.99, maxdiff= 5.59798471928e-06, npoly=20
    Func(16) = &
(0.0020751735585111336_Rkind) * tab_Pl(0) +  &
(-0.000598058304796521_Rkind) * tab_Pl(1) +  &
(-0.0029550638520216017_Rkind) * tab_Pl(2) +  &
(0.0010289985496287181_Rkind) * tab_Pl(3) +  &
(0.0007943490030897431_Rkind) * tab_Pl(4) +  &
(-0.0004422027668465322_Rkind) * tab_Pl(5) +  &
(9.766623334471966e-05_Rkind) * tab_Pl(6) +  &
(3.0359325724504228e-05_Rkind) * tab_Pl(7) +  &
(7.739060789978698e-06_Rkind) * tab_Pl(8) +  &
(-4.6319096579039554e-05_Rkind) * tab_Pl(9) +  &
(4.2182849890666805e-06_Rkind) * tab_Pl(10) +  &
(1.0820476549851587e-05_Rkind) * tab_Pl(11) +  &
(-2.0217142140798873e-05_Rkind) * tab_Pl(12) +  &
(1.4607273387152568e-05_Rkind) * tab_Pl(13) +  &
(-1.0667322977154771e-05_Rkind) * tab_Pl(14) +  &
(1.4464254129178251e-05_Rkind) * tab_Pl(15) +  &
(-4.177583038948286e-06_Rkind) * tab_Pl(16) +  &
(9.159276234433603e-06_Rkind) * tab_Pl(17) +  &
(-2.150724115332347e-07_Rkind) * tab_Pl(18) +  &
(3.401644746299729e-06_Rkind) * tab_Pl(19) +  &
(5.051485021596769e-07_Rkind) * tab_Pl(20)

!gradpi, grille ts = -0.99,-0.98...0.99, maxdiff= 0.0, il est nul, npoly=20
    Func(17) = &
(0.0_Rkind) * tab_Pl(0)

!############################ DIAGONAL TERMS OF THE HESSIAN MATRIX ######################

    !HessRo**2, grille ts = -0.99,-0.98...0.99, maxdiff= 0.000285664243025, npoly=20
    Func(18) = &
(1.799104581984879_Rkind) * tab_Pl(0) +  &
(-0.7619517944051413_Rkind) * tab_Pl(1) +  &
(-0.18217275020316742_Rkind) * tab_Pl(2) +  &
(0.25098660787657573_Rkind) * tab_Pl(3) +  &
(-0.34165739186435945_Rkind) * tab_Pl(4) +  &
(0.05645605163867279_Rkind) * tab_Pl(5) +  &
(-0.025767846060808043_Rkind) * tab_Pl(6) +  &
(-0.004279092696059913_Rkind) * tab_Pl(7) +  &
(-0.03445489681212831_Rkind) * tab_Pl(8) +  &
(0.022433601596826606_Rkind) * tab_Pl(9) +  &
(-0.02054284065027388_Rkind) * tab_Pl(10) +  &
(0.0054469204806154425_Rkind) * tab_Pl(11) +  &
(-0.003438747870988446_Rkind) * tab_Pl(12) +  &
(-0.000620350775432839_Rkind) * tab_Pl(13) +  &
(-0.00011012121822739766_Rkind) * tab_Pl(14) +  &
(-0.0007895496708164337_Rkind) * tab_Pl(15) +  &
(0.0008955709001055649_Rkind) * tab_Pl(16) +  &
(-0.00033847315821839465_Rkind) * tab_Pl(17) +  &
(0.0004640475656483518_Rkind) * tab_Pl(18) +  &
(-0.00024588262865054715_Rkind) * tab_Pl(19) +  &
(6.710639973198874e-05_Rkind) * tab_Pl(20)

    !HessRph**2, grille ts = -0.99,-0.98...0.99, maxdiff= 3.85403602237e-05, npoly=20
    Func(19) =  &
(0.45041242654706115_Rkind) * tab_Pl(0) +  &
(-0.0040431106268413505_Rkind) * tab_Pl(1) +  &
(0.0015484425097386219_Rkind) * tab_Pl(2) +  &
(-0.0005096647123239589_Rkind) * tab_Pl(3) +  &
(0.00015681639316428268_Rkind) * tab_Pl(4) +  &
(0.00019917966005773076_Rkind) * tab_Pl(5) +  &
(-7.476958792145235e-05_Rkind) * tab_Pl(6) +  &
(0.00010998885451747423_Rkind) * tab_Pl(7) +  &
(-6.230076154744126e-05_Rkind) * tab_Pl(8) +  &
(5.3218166977085293e-05_Rkind) * tab_Pl(9) +  &
(-4.8064477605248427e-05_Rkind) * tab_Pl(10) +  &
(2.0355061330251923e-05_Rkind) * tab_Pl(11) +  &
(-1.7030248090445517e-05_Rkind) * tab_Pl(12) +  &
(-4.191488089193228e-06_Rkind) * tab_Pl(13) +  &
(9.46366259964564e-06_Rkind) * tab_Pl(14) +  &
(-4.961288670825114e-06_Rkind) * tab_Pl(15) +  &
(8.984257699641388e-06_Rkind) * tab_Pl(16) +  &
(-5.327227842755125e-06_Rkind) * tab_Pl(17) +  &
(9.293289330968835e-06_Rkind) * tab_Pl(18) +  &
(-1.5959398363165459e-06_Rkind) * tab_Pl(19) +  &
(9.869713963056842e-06_Rkind) * tab_Pl(20)

    !Hessz**2, grille ts = -0.99,-0.98...0.99, maxdiff= 4.06286704381e-05, npoly=20
    Func(20) = &
(0.45271703317552603_Rkind) * tab_Pl(0) +  &
(-0.003109416154766369_Rkind) * tab_Pl(1) +  &
(0.000838508530714479_Rkind) * tab_Pl(2) +  &
(-0.000475724666186965_Rkind) * tab_Pl(3) +  &
(0.00019200173225036646_Rkind) * tab_Pl(4) +  &
(0.00015252554945385425_Rkind) * tab_Pl(5) +  &
(4.136661422908762e-06_Rkind) * tab_Pl(6) +  &
(9.542065398798985e-05_Rkind) * tab_Pl(7) +  &
(-3.3053490586592255e-05_Rkind) * tab_Pl(8) +  &
(4.936800190874369e-05_Rkind) * tab_Pl(9) +  &
(-4.391854480436164e-05_Rkind) * tab_Pl(10) +  &
(1.737069219777189e-05_Rkind) * tab_Pl(11) +  &
(-2.1308524406476963e-05_Rkind) * tab_Pl(12) +  &
(-5.054033358861408e-06_Rkind) * tab_Pl(13) +  &
(3.1818413849422487e-06_Rkind) * tab_Pl(14) +  &
(-4.45801456155866e-06_Rkind) * tab_Pl(15) +  &
(4.06630246692692e-06_Rkind) * tab_Pl(16) +  &
(-5.852245224742946e-06_Rkind) * tab_Pl(17) +  &
(6.82796954194428e-06_Rkind) * tab_Pl(18) +  &
(-1.4202159286639745e-06_Rkind) * tab_Pl(19) +  &
(8.857944457677366e-06_Rkind) * tab_Pl(20)

    !HessA**2, grille ts = -0.99,-0.98...0.99, maxdiff= 1.74560096718e-05, npoly=20
    Func(21) = &
(0.6779485553151433_Rkind) * tab_Pl(0) +  &
(0.01806118497233541_Rkind) * tab_Pl(1) +  &
(-0.0008476376400036711_Rkind) * tab_Pl(2) +  &
(-0.001825788589649897_Rkind) * tab_Pl(3) +  &
(0.000616958941705111_Rkind) * tab_Pl(4) +  &
(-0.0015611755286730339_Rkind) * tab_Pl(5) +  &
(-0.00045168788499405193_Rkind) * tab_Pl(6) +  &
(0.0001555577257159602_Rkind) * tab_Pl(7) +  &
(5.514607789813697e-05_Rkind) * tab_Pl(8) +  &
(7.029221693941506e-05_Rkind) * tab_Pl(9) +  &
(9.614286216241614e-05_Rkind) * tab_Pl(10) +  &
(0.00016100354253445394_Rkind) * tab_Pl(11) +  &
(3.718564728605394e-05_Rkind) * tab_Pl(12) +  &
(7.228004025757565e-05_Rkind) * tab_Pl(13) +  &
(7.211894144963612e-06_Rkind) * tab_Pl(14) +  &
(3.919870136341541e-06_Rkind) * tab_Pl(15) +  &
(-3.24768351000465e-05_Rkind) * tab_Pl(16) +  &
(-3.251434105002917e-05_Rkind) * tab_Pl(17) +  &
(-3.815540068344877e-05_Rkind) * tab_Pl(18) +  &
(-3.333905182842828e-05_Rkind) * tab_Pl(19) +  &
(-2.0779802703414682e-05_Rkind) * tab_Pl(20)


    !Hesszz**2, grille ts = -0.99,-0.98...0.99, maxdiff= 3.30556369793e-05, npoly=20
    Func(22) = &
(0.07738304437134337_Rkind) * tab_Pl(0) +  &
(-0.10465031099283306_Rkind) * tab_Pl(1) +  &
(0.020036938058923814_Rkind) * tab_Pl(2) +  &
(0.007902702614552323_Rkind) * tab_Pl(3) +  &
(-0.005796533303431934_Rkind) * tab_Pl(4) +  &
(0.005047688467336888_Rkind) * tab_Pl(5) +  &
(0.0017814673875459024_Rkind) * tab_Pl(6) +  &
(-0.0004317783260883068_Rkind) * tab_Pl(7) +  &
(-0.001094017585603884_Rkind) * tab_Pl(8) +  &
(0.00018274504364294727_Rkind) * tab_Pl(9) +  &
(-0.0003808090556277894_Rkind) * tab_Pl(10) +  &
(-9.246041571076469e-05_Rkind) * tab_Pl(11) +  &
(3.525985148278144e-06_Rkind) * tab_Pl(12) +  &
(1.8771571282793925e-05_Rkind) * tab_Pl(13) +  &
(5.6785090416819694e-05_Rkind) * tab_Pl(14) +  &
(-3.5132278600914584e-05_Rkind) * tab_Pl(15) +  &
(4.4624723399222785e-05_Rkind) * tab_Pl(16) +  &
(-3.32922142256894e-05_Rkind) * tab_Pl(17) +  &
(2.5506063802254586e-05_Rkind) * tab_Pl(18) +  &
(-2.0079728866285352e-05_Rkind) * tab_Pl(19) +  &
(1.7743108529267465e-05_Rkind) * tab_Pl(20)

    !HessD**2, grille ts = -0.99,-0.98...0.99, maxdiff= 2.28654271172e-05, npoly=20
    Func(23) = &
(0.06543979040717575_Rkind) * tab_Pl(0) +  &
(-0.08160323908861015_Rkind) * tab_Pl(1) +  &
(0.00804988706806203_Rkind) * tab_Pl(2) +  &
(0.006272505234968165_Rkind) * tab_Pl(3) +  &
(-0.002900958353445971_Rkind) * tab_Pl(4) +  &
(0.004094527322391887_Rkind) * tab_Pl(5) +  &
(0.0022502500073827685_Rkind) * tab_Pl(6) +  &
(-0.00032639128256260253_Rkind) * tab_Pl(7) +  &
(-0.0008755610897258378_Rkind) * tab_Pl(8) +  &
(2.3493013648110894e-05_Rkind) * tab_Pl(9) +  &
(-0.00041114966648902193_Rkind) * tab_Pl(10) +  &
(-0.00010259095784701788_Rkind) * tab_Pl(11) +  &
(-5.711774672856952e-05_Rkind) * tab_Pl(12) +  &
(4.9980102036121976e-05_Rkind) * tab_Pl(13) +  &
(2.630323825831114e-05_Rkind) * tab_Pl(14) +  &
(-5.044775741869458e-06_Rkind) * tab_Pl(15) +  &
(2.345928959612586e-05_Rkind) * tab_Pl(16) +  &
(-2.0089187785845613e-05_Rkind) * tab_Pl(17) +  &
(8.430687855244108e-06_Rkind) * tab_Pl(18) +  &
(-2.0909274963087583e-05_Rkind) * tab_Pl(19) +  &
(7.08060514827649e-06_Rkind) * tab_Pl(20)
 


    !HesstAH2**2, grille ts = -0.99,-0.98...0.99, maxdiff= 1.7319247564e-05, npoly=20
    Func(24) = &
(0.00977114326642278_Rkind) * tab_Pl(0) +  &
(-0.003599894265318701_Rkind) * tab_Pl(1) +  &
(-0.011701781631389838_Rkind) * tab_Pl(2) +  &
(0.004275187370510517_Rkind) * tab_Pl(3) +  &
(0.0007257496354772547_Rkind) * tab_Pl(4) +  &
(-1.3607799116940379e-05_Rkind) * tab_Pl(5) +  &
(0.0007296778684732192_Rkind) * tab_Pl(6) +  &
(-0.00015484158189192014_Rkind) * tab_Pl(7) +  &
(0.0003455676301793171_Rkind) * tab_Pl(8) +  &
(-0.0004031059565433909_Rkind) * tab_Pl(9) +  &
(0.00018679114095039538_Rkind) * tab_Pl(10) +  &
(-0.00017211989312388999_Rkind) * tab_Pl(11) +  &
(4.8098453378384165e-05_Rkind) * tab_Pl(12) +  &
(-4.984170651038166e-05_Rkind) * tab_Pl(13) +  &
(-1.0886796280145449e-06_Rkind) * tab_Pl(14) +  &
(1.4919946977188191e-05_Rkind) * tab_Pl(15) +  &
(-1.432371002603214e-05_Rkind) * tab_Pl(16) +  &
(3.834975914110647e-05_Rkind) * tab_Pl(17) +  &
(-7.897200885241417e-06_Rkind) * tab_Pl(18) +  &
(3.380551974991322e-05_Rkind) * tab_Pl(19) +  &
(-5.647734143788724e-06_Rkind) * tab_Pl(20)

    !Hesspi**2, grille ts = -0.99,-0.98...0.99, maxdiff= 1.36163773673e-05, npoly=20
    Func(25) =  &
(0.010454362408658087_Rkind) * tab_Pl(0) +  &
(-0.0037048456188888707_Rkind) * tab_Pl(1) +  &
(-0.012749364392705478_Rkind) * tab_Pl(2) +  &
(0.004569031855441761_Rkind) * tab_Pl(3) +  &
(0.0010549447888864776_Rkind) * tab_Pl(4) +  &
(-0.0002782483567478416_Rkind) * tab_Pl(5) +  &
(0.0008356510992181551_Rkind) * tab_Pl(6) +  &
(-0.00012089325679123623_Rkind) * tab_Pl(7) +  &
(0.00029162134876976247_Rkind) * tab_Pl(8) +  &
(-0.00036953189715434494_Rkind) * tab_Pl(9) +  &
(0.00017673758254304204_Rkind) * tab_Pl(10) +  &
(-0.00015761374049872652_Rkind) * tab_Pl(11) +  &
(3.254051498697846e-05_Rkind) * tab_Pl(12) +  &
(-3.7939048801410664e-05_Rkind) * tab_Pl(13) +  &
(-1.1351319139635261e-05_Rkind) * tab_Pl(14) +  &
(1.7300784548391476e-05_Rkind) * tab_Pl(15) +  &
(-1.7103238864378383e-05_Rkind) * tab_Pl(16) +  &
(3.727032695033268e-05_Rkind) * tab_Pl(17) +  &
(-5.858486521326748e-06_Rkind) * tab_Pl(18) +  &
(3.173037941496075e-05_Rkind) * tab_Pl(19) +  &
(-2.0897699449896404e-06_Rkind) * tab_Pl(20)

!############################ OFF DIAGONAL TERMS OF THE HESSIAN MATRIX ######################

!Hesspita, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 2.6339888214715572e-09
    Func(26) =  &
(3.3097491509931774e-11_Rkind) * tab_Pl(0) +  &
(4.2975989319490426e-10_Rkind) * tab_Pl(1) +  &
(-6.375939653535891e-10_Rkind) * tab_Pl(2) +  &
(-4.1938878810642023e-10_Rkind) * tab_Pl(3) +  &
(5.304418959554867e-10_Rkind) * tab_Pl(4) +  &
(2.1906532893317103e-10_Rkind) * tab_Pl(5) +  &
(1.1085771698402418e-09_Rkind) * tab_Pl(6) +  &
(-1.0494892092957112e-09_Rkind) * tab_Pl(7) +  &
(-1.0596451456964153e-09_Rkind) * tab_Pl(8) +  &
(5.073586187207524e-10_Rkind) * tab_Pl(9) +  &
(3.927089637233884e-10_Rkind) * tab_Pl(10) +  &
(1.2763483379139436e-09_Rkind) * tab_Pl(11) +  &
(-4.213902507542741e-10_Rkind) * tab_Pl(12) +  &
(-1.4916628573639149e-09_Rkind) * tab_Pl(13) +  &
(4.593860357440681e-10_Rkind) * tab_Pl(14) +  &
(4.946203544982609e-10_Rkind) * tab_Pl(15) +  &
(8.441359617186701e-10_Rkind) * tab_Pl(16) +  &
(-2.1938345157933267e-10_Rkind) * tab_Pl(17) +  &
(-1.24567718463281e-09_Rkind) * tab_Pl(18) +  &
(9.55761185576788e-12_Rkind) * tab_Pl(19) +  &
(1.2167122742824445e-09_Rkind) * tab_Pl(20)
 
!HesspiD, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.502197092647489e-08
    Func(27) =  &
(-4.5525724577101553e-10_Rkind) * tab_Pl(0) +  &
(3.4096936489145727e-09_Rkind) * tab_Pl(1) +  &
(2.3769588377227906e-09_Rkind) * tab_Pl(2) +  &
(1.2048242704435133e-09_Rkind) * tab_Pl(3) +  &
(-3.108447413279234e-09_Rkind) * tab_Pl(4) +  &
(-1.1147653065042505e-08_Rkind) * tab_Pl(5) +  &
(-4.8566933730458705e-09_Rkind) * tab_Pl(6) +  &
(4.218028136586465e-09_Rkind) * tab_Pl(7) +  &
(6.670603801803451e-09_Rkind) * tab_Pl(8) +  &
(6.690343004821456e-09_Rkind) * tab_Pl(9) +  &
(-2.332271631546395e-09_Rkind) * tab_Pl(10) +  &
(-4.460725508209964e-09_Rkind) * tab_Pl(11) +  &
(1.8131314652550294e-10_Rkind) * tab_Pl(12) +  &
(8.766326964988432e-10_Rkind) * tab_Pl(13) +  &
(-3.5193484139794733e-09_Rkind) * tab_Pl(14) +  &
(-3.4915143194191607e-09_Rkind) * tab_Pl(15) +  &
(-3.5941402570596886e-09_Rkind) * tab_Pl(16) +  &
(4.75590561596884e-09_Rkind) * tab_Pl(17) +  &
(5.134653442941441e-09_Rkind) * tab_Pl(18) +  &
(-2.1039998526162417e-09_Rkind) * tab_Pl(19) +  &
(-7.981799569983433e-09_Rkind) * tab_Pl(20)
 
!HessDta, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.2629994880191461e-05 
    Func(28) =  &
(0.0011834844143964233_Rkind) * tab_Pl(0) +  &
(0.0022071343797770843_Rkind) * tab_Pl(1) +  &
(-0.0028462992780473725_Rkind) * tab_Pl(2) +  &
(-0.0029663649828559_Rkind) * tab_Pl(3) +  &
(0.0019368056950353888_Rkind) * tab_Pl(4) +  &
(0.000644121600539268_Rkind) * tab_Pl(5) +  &
(-0.00012526327679801961_Rkind) * tab_Pl(6) +  &
(-2.2736210599167337e-05_Rkind) * tab_Pl(7) +  &
(6.132108340908813e-05_Rkind) * tab_Pl(8) +  &
(2.0174097643413158e-05_Rkind) * tab_Pl(9) +  &
(-0.00016560077854785703_Rkind) * tab_Pl(10) +  &
(9.983540146744227e-05_Rkind) * tab_Pl(11) +  &
(-5.7705808271949185e-05_Rkind) * tab_Pl(12) +  &
(5.306448449566367e-05_Rkind) * tab_Pl(13) +  &
(-3.524144312057126e-05_Rkind) * tab_Pl(14) +  &
(1.7176319553741273e-05_Rkind) * tab_Pl(15) +  &
(-1.7666954147179086e-05_Rkind) * tab_Pl(16) +  &
(-6.07781860021474e-06_Rkind) * tab_Pl(17) +  &
(-3.740542008922121e-06_Rkind) * tab_Pl(18) +  &
(-7.747014754051466e-06_Rkind) * tab_Pl(19) +  &
(5.053683202371419e-06_Rkind) * tab_Pl(20)
 
!HessZzPi, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.5976394687253262e-05
    Func(29) =  &
(0.001583041965354985_Rkind) * tab_Pl(0) +  &
(-0.0029235428786131084_Rkind) * tab_Pl(1) +  &
(-0.00039440034872019643_Rkind) * tab_Pl(2) +  &
(0.003832940748165007_Rkind) * tab_Pl(3) +  &
(-0.0018457044433028242_Rkind) * tab_Pl(4) +  &
(-0.0006175233689728172_Rkind) * tab_Pl(5) +  &
(0.0003997629920085231_Rkind) * tab_Pl(6) +  &
(-9.377590700961736e-05_Rkind) * tab_Pl(7) +  &
(4.084504158778548e-05_Rkind) * tab_Pl(8) +  &
(-6.940252682056513e-05_Rkind) * tab_Pl(9) +  &
(0.00016764708826312987_Rkind) * tab_Pl(10) +  &
(-0.0001265327031400613_Rkind) * tab_Pl(11) +  &
(7.335509386568646e-05_Rkind) * tab_Pl(12) +  &
(-5.9149467589467565e-05_Rkind) * tab_Pl(13) +  &
(4.050779868861533e-05_Rkind) * tab_Pl(14) +  &
(-1.2715053967680411e-05_Rkind) * tab_Pl(15) +  &
(1.7284806231732362e-05_Rkind) * tab_Pl(16) +  &
(1.1007408300688593e-05_Rkind) * tab_Pl(17) +  &
(-1.023515787794369e-06_Rkind) * tab_Pl(18) +  &
(1.1077748781528652e-05_Rkind) * tab_Pl(19) +  &
(-9.27696511661586e-06_Rkind) * tab_Pl(20)
 
!HessZzta, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.644056497415068e-08
    Func(30) =  &
(1.6627058092450162e-09_Rkind) * tab_Pl(0) +  &
(-6.877392441652261e-10_Rkind) * tab_Pl(1) +  &
(2.0954777027606735e-09_Rkind) * tab_Pl(2) +  &
(-3.3379262406339418e-09_Rkind) * tab_Pl(3) +  &
(-3.579250549194597e-09_Rkind) * tab_Pl(4) +  &
(1.7781595764083365e-10_Rkind) * tab_Pl(5) +  &
(3.797514781461969e-09_Rkind) * tab_Pl(6) +  &
(6.858585088229106e-09_Rkind) * tab_Pl(7) +  &
(4.5249256706441396e-09_Rkind) * tab_Pl(8) +  &
(-1.9761895387667936e-09_Rkind) * tab_Pl(9) +  &
(-1.1121057528913378e-08_Rkind) * tab_Pl(10) +  &
(-2.944791716630141e-09_Rkind) * tab_Pl(11) +  &
(-5.738369246827375e-09_Rkind) * tab_Pl(12) +  &
(1.0818941688875703e-08_Rkind) * tab_Pl(13) +  &
(3.938008195082549e-09_Rkind) * tab_Pl(14) +  &
(5.636118797340932e-09_Rkind) * tab_Pl(15) +  &
(-2.931771996895781e-09_Rkind) * tab_Pl(16) +  &
(-2.511276066765326e-09_Rkind) * tab_Pl(17) +  &
(-1.2431134832421502e-09_Rkind) * tab_Pl(18) +  &
(3.020218406443747e-10_Rkind) * tab_Pl(19) +  &
(5.659261200949319e-09_Rkind) * tab_Pl(20)

!HessZzD, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 3.912865898919908e-08
    Func(31) =  & 
(4.381825533213257e-09_Rkind) * tab_Pl(0) +  &
(-5.9977254175405505e-09_Rkind) * tab_Pl(1) +  &
(6.3260021960720234e-09_Rkind) * tab_Pl(2) +  &
(4.662095206396266e-09_Rkind) * tab_Pl(3) +  &
(-2.8327854615879553e-09_Rkind) * tab_Pl(4) +  &
(1.3836435071010647e-08_Rkind) * tab_Pl(5) +  &
(3.1754641217621965e-09_Rkind) * tab_Pl(6) +  &
(-9.874102074243807e-09_Rkind) * tab_Pl(7) +  &
(-3.2933545702192735e-10_Rkind) * tab_Pl(8) +  &
(-8.487286852740891e-09_Rkind) * tab_Pl(9) +  &
(-1.1527877071869253e-08_Rkind) * tab_Pl(10) +  &
(1.637552278856249e-08_Rkind) * tab_Pl(11) +  &
(-1.1192066800165707e-08_Rkind) * tab_Pl(12) +  &
(1.0227910044665707e-08_Rkind) * tab_Pl(13) +  &
(9.153308230505552e-09_Rkind) * tab_Pl(14) +  &
(8.610406623631362e-09_Rkind) * tab_Pl(15) +  &
(1.0132146423180588e-08_Rkind) * tab_Pl(16) +  &
(1.5001862451363032e-08_Rkind) * tab_Pl(17) +  &
(-4.225975938626148e-09_Rkind) * tab_Pl(18) +  &
(1.026170403977936e-08_Rkind) * tab_Pl(19) +  &
(1.0790642498178312e-08_Rkind) * tab_Pl(20)
 
!HessAnpi, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 5.1006870896143016e-09
    Func(32) =  &
(-1.6535115908997164e-09_Rkind) * tab_Pl(0) +  &
(7.251748612795148e-10_Rkind) * tab_Pl(1) +  &
(-9.880946813797354e-10_Rkind) * tab_Pl(2) +  &
(2.6418112890781285e-09_Rkind) * tab_Pl(3) +  &
(2.9391723526845086e-09_Rkind) * tab_Pl(4) +  &
(-2.8596291088581247e-09_Rkind) * tab_Pl(5) +  &
(1.6391424319834902e-09_Rkind) * tab_Pl(6) +  &
(-7.721922227432013e-10_Rkind) * tab_Pl(7) +  &
(4.631033937099313e-10_Rkind) * tab_Pl(8) +  &
(1.9123153785322747e-10_Rkind) * tab_Pl(9) +  &
(-4.175901352189042e-09_Rkind) * tab_Pl(10) +  &
(2.3243972378774033e-10_Rkind) * tab_Pl(11) +  &
(2.1184489816409773e-09_Rkind) * tab_Pl(12) +  &
(4.0927061165427566e-09_Rkind) * tab_Pl(13) +  &
(6.39161880721275e-10_Rkind) * tab_Pl(14) +  &
(-2.364394283052851e-09_Rkind) * tab_Pl(15) +  &
(-2.7140259656308163e-09_Rkind) * tab_Pl(16) +  &
(2.098861588502198e-09_Rkind) * tab_Pl(17) +  &
(4.364726958271089e-09_Rkind) * tab_Pl(18) +  &
(2.5865328123617384e-09_Rkind) * tab_Pl(19) +  &
(-1.6499488045231147e-09_Rkind) * tab_Pl(20)
 
!HessAnta, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.3609906215631457e-05
    Func(33) =  &
(0.0008112932192716962_Rkind) * tab_Pl(0) +  &
(1.2627788671339986e-06_Rkind) * tab_Pl(1) +  &
(-0.0011971309944768561_Rkind) * tab_Pl(2) +  &
(-5.804920258538941e-05_Rkind) * tab_Pl(3) +  &
(0.0004659215594493744_Rkind) * tab_Pl(4) +  &
(2.0802051311704433e-05_Rkind) * tab_Pl(5) +  &
(-7.534941692858729e-05_Rkind) * tab_Pl(6) +  &
(8.780734157316942e-05_Rkind) * tab_Pl(7) +  &
(2.17905277848781e-06_Rkind) * tab_Pl(8) +  &
(-3.5681919079347085e-05_Rkind) * tab_Pl(9) +  &
(-1.1091997426841942e-05_Rkind) * tab_Pl(10) +  &
(-6.312877423294159e-06_Rkind) * tab_Pl(11) +  &
(-6.067919189498752e-06_Rkind) * tab_Pl(12) +  &
(-6.429908197382552e-06_Rkind) * tab_Pl(13) +  &
(3.154990859649583e-06_Rkind) * tab_Pl(14) +  &
(9.736034009098133e-07_Rkind) * tab_Pl(15) +  &
(6.718541652026072e-06_Rkind) * tab_Pl(16) +  &
(4.841132475813797e-06_Rkind) * tab_Pl(17) +  &
(8.727157208315401e-06_Rkind) * tab_Pl(18) +  &
(5.397945240941394e-06_Rkind) * tab_Pl(19) +  &
(6.523199073411057e-06_Rkind) * tab_Pl(20)
 
!HessAnD, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.063275495825775e-05
    Func(34) =  &
(-0.01679280273123769_Rkind) * tab_Pl(0) +  &
(0.02388201908733393_Rkind) * tab_Pl(1) +  &
(-0.006541662458693422_Rkind) * tab_Pl(2) +  &
(-0.0021322079344829005_Rkind) * tab_Pl(3) +  &
(0.0025897839518964983_Rkind) * tab_Pl(4) +  &
(-0.000668143394635906_Rkind) * tab_Pl(5) +  &
(-0.00029109149310124044_Rkind) * tab_Pl(6) +  &
(-4.668979109832323e-05_Rkind) * tab_Pl(7) +  &
(0.0001186608804818445_Rkind) * tab_Pl(8) +  &
(-0.000139184959282619_Rkind) * tab_Pl(9) +  &
(-4.816101965434525e-06_Rkind) * tab_Pl(10) +  &
(3.516012998676312e-05_Rkind) * tab_Pl(11) +  &
(-4.126178987591031e-05_Rkind) * tab_Pl(12) +  &
(1.1661391201337226e-05_Rkind) * tab_Pl(13) +  &
(-3.6906953803773507e-05_Rkind) * tab_Pl(14) +  &
(1.1780887873489234e-05_Rkind) * tab_Pl(15) +  &
(-2.0131618549274535e-05_Rkind) * tab_Pl(16) +  &
(9.678463332866107e-06_Rkind) * tab_Pl(17) +  &
(-1.8948857851863556e-06_Rkind) * tab_Pl(18) +  &
(6.143320785920671e-06_Rkind) * tab_Pl(19) +  &
(2.397394713188703e-06_Rkind) * tab_Pl(20)
 
!HessZpi, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 6.114277905450242e-06
    Func(35) =  &
(0.0012325666583301218_Rkind) * tab_Pl(0) +  &
(-3.7463563121389624e-05_Rkind) * tab_Pl(1) +  &
(-0.0015349524867922222_Rkind) * tab_Pl(2) +  &
(0.00017355049065349436_Rkind) * tab_Pl(3) +  &
(0.00010361860308125614_Rkind) * tab_Pl(4) +  &
(-0.00013147762613819431_Rkind) * tab_Pl(5) +  &
(0.00018703414689092776_Rkind) * tab_Pl(6) +  &
(-2.19884022580464e-05_Rkind) * tab_Pl(7) +  &
(2.803014646082093e-05_Rkind) * tab_Pl(8) +  &
(-4.613822473514448e-06_Rkind) * tab_Pl(9) +  &
(-3.5266917362710227e-06_Rkind) * tab_Pl(10) +  &
(5.226642113097011e-06_Rkind) * tab_Pl(11) +  &
(-7.773040703651397e-06_Rkind) * tab_Pl(12) +  &
(1.1252048384718246e-05_Rkind) * tab_Pl(13) +  &
(-6.566796011791039e-06_Rkind) * tab_Pl(14) +  &
(7.156928374355649e-06_Rkind) * tab_Pl(15) +  &
(-6.1438763910631905e-06_Rkind) * tab_Pl(16) +  &
(1.178808026799229e-06_Rkind) * tab_Pl(17) +  &
(-5.503829275541059e-06_Rkind) * tab_Pl(18) +  &
(-1.730074637238755e-06_Rkind) * tab_Pl(19) +  &
(-3.099961110094051e-06_Rkind) * tab_Pl(20)
 
!HessZta, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 2.7772730053805758e-09
    Func(36) =  &
(2.81372366747926e-10_Rkind) * tab_Pl(0) +  &
(1.0222025909188988e-09_Rkind) * tab_Pl(1) +  &
(-4.6019282247668157e-10_Rkind) * tab_Pl(2) +  &
(-1.1591539706973424e-09_Rkind) * tab_Pl(3) +  &
(-1.835382500417838e-10_Rkind) * tab_Pl(4) +  &
(-3.387707035608997e-10_Rkind) * tab_Pl(5) +  &
(3.4916314349127836e-11_Rkind) * tab_Pl(6) +  &
(-7.429276069303993e-10_Rkind) * tab_Pl(7) +  &
(-8.233797641480066e-10_Rkind) * tab_Pl(8) +  &
(7.260117589312264e-10_Rkind) * tab_Pl(9) +  &
(3.156239307076389e-09_Rkind) * tab_Pl(10) +  &
(4.382980402662963e-10_Rkind) * tab_Pl(11) +  &
(-3.852399085672866e-10_Rkind) * tab_Pl(12) +  &
(-2.7836604927768294e-09_Rkind) * tab_Pl(13) +  &
(-1.3924543160771366e-09_Rkind) * tab_Pl(14) +  &
(1.0785999729284924e-09_Rkind) * tab_Pl(15) +  &
(1.5045379485926124e-09_Rkind) * tab_Pl(16) +  &
(-1.1188739978239344e-09_Rkind) * tab_Pl(17) +  &
(-9.912401764763419e-10_Rkind) * tab_Pl(18) +  &
(-1.26714087756372e-09_Rkind) * tab_Pl(19) +  &
(8.179980255588148e-10_Rkind) * tab_Pl(20)
 
!HessZD, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 2.5697652486691395e-08
    Func(37) =  &
(1.0666026678914382e-09_Rkind) * tab_Pl(0) +  &
(6.395422377215013e-11_Rkind) * tab_Pl(1) +  &
(1.0094195294228541e-08_Rkind) * tab_Pl(2) +  &
(1.4648136287471984e-08_Rkind) * tab_Pl(3) +  &
(5.855280704751842e-09_Rkind) * tab_Pl(4) +  &
(1.926899008642167e-09_Rkind) * tab_Pl(5) +  &
(-1.1215690998846404e-08_Rkind) * tab_Pl(6) +  &
(-1.5211279885112256e-08_Rkind) * tab_Pl(7) +  &
(-6.940653184016761e-09_Rkind) * tab_Pl(8) +  &
(-4.1415734032332104e-09_Rkind) * tab_Pl(9) +  &
(-2.2794844150315065e-09_Rkind) * tab_Pl(10) +  &
(9.129505440481763e-09_Rkind) * tab_Pl(11) +  &
(1.0364395884482917e-09_Rkind) * tab_Pl(12) +  &
(4.43288380556854e-09_Rkind) * tab_Pl(13) +  &
(2.6956410765408596e-09_Rkind) * tab_Pl(14) +  &
(2.4408546526401712e-09_Rkind) * tab_Pl(15) +  &
(1.1678151806115635e-08_Rkind) * tab_Pl(16) +  &
(2.2648554050947156e-08_Rkind) * tab_Pl(17) +  &
(1.2784377239579821e-08_Rkind) * tab_Pl(18) +  &
(8.412358183911222e-09_Rkind) * tab_Pl(19) +  &
(-4.387744349433441e-09_Rkind) * tab_Pl(20)
 
!HessRhpi, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.2632970271427035e-09 
    Func(38) =  &
(-4.4408775572114225e-10_Rkind) * tab_Pl(0) +  &
(1.0166856702493142e-09_Rkind) * tab_Pl(1) +  &
(-7.698286716233555e-10_Rkind) * tab_Pl(2) +  &
(-4.2356265287507896e-10_Rkind) * tab_Pl(3) +  &
(8.513897916153987e-10_Rkind) * tab_Pl(4) +  &
(-1.713883454609742e-09_Rkind) * tab_Pl(5) +  &
(1.060997129966488e-09_Rkind) * tab_Pl(6) +  &
(4.4340410730675984e-10_Rkind) * tab_Pl(7) +  &
(-4.440455085200101e-10_Rkind) * tab_Pl(8) +  &
(8.09131299205951e-10_Rkind) * tab_Pl(9) +  &
(-8.221263284964996e-10_Rkind) * tab_Pl(10) +  &
(-6.940479254073789e-10_Rkind) * tab_Pl(11) +  &
(4.158456317984518e-10_Rkind) * tab_Pl(12) +  &
(-2.1257114846750365e-10_Rkind) * tab_Pl(13) +  &
(-1.9361985471383632e-11_Rkind) * tab_Pl(14) +  &
(2.7573119932118623e-10_Rkind) * tab_Pl(15) +  &
(-4.598501956271064e-10_Rkind) * tab_Pl(16) +  &
(-4.1848477679471833e-10_Rkind) * tab_Pl(17) +  &
(-5.744515061827558e-11_Rkind) * tab_Pl(18) +  &
(-8.970912328681817e-10_Rkind) * tab_Pl(19) +  &
(-3.865337156382796e-10_Rkind) * tab_Pl(20)
 
!HessRhta, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.817986182391554e-06 
    Func(39) =  &
(-0.0009299848567277572_Rkind) * tab_Pl(0) +  &
(8.405416553784717e-05_Rkind) * tab_Pl(1) +  &
(0.0010548478339741127_Rkind) * tab_Pl(2) +  &
(-0.0001463635840200269_Rkind) * tab_Pl(3) +  &
(3.7158571369797476e-05_Rkind) * tab_Pl(4) +  &
(5.284118260107424e-05_Rkind) * tab_Pl(5) +  &
(-0.00013381658710947506_Rkind) * tab_Pl(6) +  &
(5.049305985259327e-06_Rkind) * tab_Pl(7) +  &
(-2.4161828215619866e-05_Rkind) * tab_Pl(8) +  &
(4.7088667938235965e-06_Rkind) * tab_Pl(9) +  &
(-8.641937599722844e-06_Rkind) * tab_Pl(10) +  &
(5.220703745407529e-06_Rkind) * tab_Pl(11) +  &
(1.0343582046753133e-07_Rkind) * tab_Pl(12) +  &
(-7.645991331870151e-07_Rkind) * tab_Pl(13) +  &
(1.8768275368852836e-06_Rkind) * tab_Pl(14) +  &
(-1.4638645513377496e-06_Rkind) * tab_Pl(15) +  &
(2.1578670972073563e-06_Rkind) * tab_Pl(16) +  &
(-7.258129850235891e-07_Rkind) * tab_Pl(17) +  &
(1.9007413049951436e-06_Rkind) * tab_Pl(18) +  &
(1.0010533567225492e-07_Rkind) * tab_Pl(19) +  &
(1.1072091422676509e-06_Rkind) * tab_Pl(20)
 
!HessRhD, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.3402721008915852e-05
    Func(40) =  &
(0.003866673371577077_Rkind) * tab_Pl(0) +  &
(-0.008792135368519416_Rkind) * tab_Pl(1) +  &
(0.005494667470024973_Rkind) * tab_Pl(2) +  &
(0.0007468260095592898_Rkind) * tab_Pl(3) +  &
(-0.0011873267077107822_Rkind) * tab_Pl(4) +  &
(0.0004706266364632608_Rkind) * tab_Pl(5) +  &
(-0.0004814504547327046_Rkind) * tab_Pl(6) +  &
(-0.00023738648757052033_Rkind) * tab_Pl(7) +  &
(-7.641579760915283e-05_Rkind) * tab_Pl(8) +  &
(6.082412852090778e-05_Rkind) * tab_Pl(9) +  &
(3.199313206123406e-05_Rkind) * tab_Pl(10) +  &
(3.4313332461288075e-05_Rkind) * tab_Pl(11) +  &
(6.585397135930844e-05_Rkind) * tab_Pl(12) +  &
(1.8184280675050758e-05_Rkind) * tab_Pl(13) +  &
(4.19209045919912e-05_Rkind) * tab_Pl(14) +  &
(9.15427206263082e-06_Rkind) * tab_Pl(15) +  &
(2.455738067435636e-05_Rkind) * tab_Pl(16) +  &
(2.486606555861679e-06_Rkind) * tab_Pl(17) +  &
(6.18852204554642e-06_Rkind) * tab_Pl(18) +  &
(-3.4808295484476012e-06_Rkind) * tab_Pl(19) +  &
(-4.370844419921678e-06_Rkind) * tab_Pl(20)
 
!HessRhopi, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 8.224288536505827e-09 
    Func(41) =  &
(5.459093430355546e-10_Rkind) * tab_Pl(0) +  &
(-3.714479215928717e-09_Rkind) * tab_Pl(1) +  &
(1.5059757001526248e-09_Rkind) * tab_Pl(2) +  &
(1.7208929832726858e-09_Rkind) * tab_Pl(3) +  &
(-2.4752962851919963e-09_Rkind) * tab_Pl(4) +  &
(5.173692344263941e-09_Rkind) * tab_Pl(5) +  &
(-2.2977875994736843e-09_Rkind) * tab_Pl(6) +  &
(-5.632398367989677e-10_Rkind) * tab_Pl(7) +  &
(2.6973032268194727e-09_Rkind) * tab_Pl(8) +  &
(-6.578352173413469e-09_Rkind) * tab_Pl(9) +  &
(4.22699495163663e-10_Rkind) * tab_Pl(10) +  &
(3.8020155805257645e-10_Rkind) * tab_Pl(11) +  &
(-1.1023113166398077e-09_Rkind) * tab_Pl(12) +  &
(5.277643635642196e-09_Rkind) * tab_Pl(13) +  &
(-2.030854211483892e-09_Rkind) * tab_Pl(14) +  &
(-2.6404168017574475e-09_Rkind) * tab_Pl(15) +  &
(1.1939002196513865e-09_Rkind) * tab_Pl(16) +  &
(-2.2230095954810098e-09_Rkind) * tab_Pl(17) +  &
(1.4386835997964666e-09_Rkind) * tab_Pl(18) +  &
(-3.4872407785416656e-10_Rkind) * tab_Pl(19) +  &
(-3.1357462955602856e-09_Rkind) * tab_Pl(20)
 
!HessRhota, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.5364169535359398e-05
    Func(42) =  &
(-0.007101558282085527_Rkind) * tab_Pl(0) +  &
(0.002225830522349768_Rkind) * tab_Pl(1) +  &
(0.009842374732185681_Rkind) * tab_Pl(2) +  &
(-0.002517415076921176_Rkind) * tab_Pl(3) +  &
(-0.002604656812741194_Rkind) * tab_Pl(4) +  &
(-0.0007328677145053199_Rkind) * tab_Pl(5) +  &
(0.00045623907257729985_Rkind) * tab_Pl(6) +  &
(0.0008528232620815053_Rkind) * tab_Pl(7) +  &
(-0.0007966722896381915_Rkind) * tab_Pl(8) +  &
(0.0003615780052140993_Rkind) * tab_Pl(9) +  &
(7.374032847956842e-05_Rkind) * tab_Pl(10) +  &
(-0.00011192989895974917_Rkind) * tab_Pl(11) +  &
(7.416603051971318e-05_Rkind) * tab_Pl(12) +  &
(-1.9497574077774826e-05_Rkind) * tab_Pl(13) +  &
(3.5594843379340484e-05_Rkind) * tab_Pl(14) +  &
(-6.295321493845545e-05_Rkind) * tab_Pl(15) +  &
(1.98456557573669e-05_Rkind) * tab_Pl(16) +  &
(-5.544922897565762e-05_Rkind) * tab_Pl(17) +  &
(8.845790847643446e-06_Rkind) * tab_Pl(18) +  &
(-2.515659477517537e-05_Rkind) * tab_Pl(19) +  &
(1.1748018861864634e-05_Rkind) * tab_Pl(20)
 
!HessRhoD grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 5.0890907880806905e-05 
    Func(43) =  &
(0.006249547181647878_Rkind) * tab_Pl(0) +  &
(-0.0204108052471768_Rkind) * tab_Pl(1) +  &
(0.017159688551770174_Rkind) * tab_Pl(2) +  &
(0.0014505705050275944_Rkind) * tab_Pl(3) +  &
(-0.004442874264716959_Rkind) * tab_Pl(4) +  &
(0.0006690220580563662_Rkind) * tab_Pl(5) +  &
(-0.0021272227197485522_Rkind) * tab_Pl(6) +  &
(0.0018518373107978968_Rkind) * tab_Pl(7) +  &
(-0.0008435860585145145_Rkind) * tab_Pl(8) +  &
(0.00027440291052640106_Rkind) * tab_Pl(9) +  &
(0.0003216880341582153_Rkind) * tab_Pl(10) +  &
(-0.000141381227889017_Rkind) * tab_Pl(11) +  &
(3.6188103799356554e-05_Rkind) * tab_Pl(12) +  &
(-6.082402949664792e-06_Rkind) * tab_Pl(13) +  &
(-3.1096191205753564e-05_Rkind) * tab_Pl(14) +  &
(-8.055769091220619e-05_Rkind) * tab_Pl(15) +  &
(-7.62279926489148e-06_Rkind) * tab_Pl(16) +  &
(-4.13464160896169e-05_Rkind) * tab_Pl(17) +  &
(1.7981439863562586e-05_Rkind) * tab_Pl(18) +  &
(-6.880724261316996e-06_Rkind) * tab_Pl(19) +  &
(1.5355775101155088e-05_Rkind) * tab_Pl(20)
 
!HessZzAn, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 2.1728102824379366e-08
    Func(44) =  &
(2.3781541437559926e-09_Rkind) * tab_Pl(0) +  &
(-6.369587492285533e-09_Rkind) * tab_Pl(1) +  &
(4.419920974189995e-09_Rkind) * tab_Pl(2) +  &
(5.063385243737963e-09_Rkind) * tab_Pl(3) +  &
(-4.357756159567804e-09_Rkind) * tab_Pl(4) +  &
(9.446446300630645e-09_Rkind) * tab_Pl(5) +  &
(-1.2698242634458288e-08_Rkind) * tab_Pl(6) +  &
(-1.2447959421601986e-08_Rkind) * tab_Pl(7) +  &
(-2.3458474685340356e-09_Rkind) * tab_Pl(8) +  &
(-1.854838374682251e-09_Rkind) * tab_Pl(9) +  &
(1.5634133763599562e-08_Rkind) * tab_Pl(10) +  &
(1.1451793583488786e-08_Rkind) * tab_Pl(11) +  &
(1.630687921869664e-10_Rkind) * tab_Pl(12) +  &
(-1.0732339288524848e-08_Rkind) * tab_Pl(13) +  &
(-4.341934690829229e-09_Rkind) * tab_Pl(14) +  &
(-6.679024919647917e-09_Rkind) * tab_Pl(15) +  &
(6.806485172792141e-09_Rkind) * tab_Pl(16) +  &
(3.1674810631087222e-09_Rkind) * tab_Pl(17) +  &
(-5.3467516434929e-09_Rkind) * tab_Pl(18) +  &
(-5.720092928581802e-09_Rkind) * tab_Pl(19) +  &
(-5.040013878358486e-09_Rkind) * tab_Pl(20)
 
!HessZz_z, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff : 1.3713522114190067e-05
    Func(45) =  &
(0.0015673116956567664_Rkind) * tab_Pl(0) +  &
(-0.005738139389109408_Rkind) * tab_Pl(1) +  &
(0.0047403555567041_Rkind) * tab_Pl(2) +  &
(0.0007908246849767056_Rkind) * tab_Pl(3) +  &
(-0.0010478316039068325_Rkind) * tab_Pl(4) +  &
(0.00029411814512235294_Rkind) * tab_Pl(5) +  &
(-0.0005323302275464585_Rkind) * tab_Pl(6) +  &
(-0.0002376742441327883_Rkind) * tab_Pl(7) +  &
(-4.064668873431057e-05_Rkind) * tab_Pl(8) +  &
(4.003650182441933e-05_Rkind) * tab_Pl(9) +  &
(5.164545378708012e-05_Rkind) * tab_Pl(10) +  &
(2.8278385465589863e-05_Rkind) * tab_Pl(11) +  &
(6.288475983491751e-05_Rkind) * tab_Pl(12) +  &
(7.932186064425872e-06_Rkind) * tab_Pl(13) +  &
(3.20983682320824e-05_Rkind) * tab_Pl(14) +  &
(1.752832316871233e-06_Rkind) * tab_Pl(15) +  &
(1.3910458231890923e-05_Rkind) * tab_Pl(16) +  &
(-2.8190288382324557e-06_Rkind) * tab_Pl(17) +  &
(-1.0133312315431874e-06_Rkind) * tab_Pl(18) +  &
(-6.31824446116855e-06_Rkind) * tab_Pl(19) +  &
(-7.891345063451792e-06_Rkind) * tab_Pl(20)
 
!HessZzRh, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.2491618938924514e-08
    Func(46) =  &
(1.1248592369473597e-09_Rkind) * tab_Pl(0) +  &
(-4.386102618685795e-09_Rkind) * tab_Pl(1) +  &
(3.309120351035958e-09_Rkind) * tab_Pl(2) +  &
(4.424612048614941e-09_Rkind) * tab_Pl(3) +  &
(-2.8203308780901323e-09_Rkind) * tab_Pl(4) +  &
(3.4304539382760456e-09_Rkind) * tab_Pl(5) +  &
(-8.36987199877581e-09_Rkind) * tab_Pl(6) +  &
(-9.040779057543994e-09_Rkind) * tab_Pl(7) +  &
(1.2882297539228339e-09_Rkind) * tab_Pl(8) +  &
(1.6493378404868602e-09_Rkind) * tab_Pl(9) +  &
(8.807427255691443e-09_Rkind) * tab_Pl(10) +  &
(1.3275189213081374e-08_Rkind) * tab_Pl(11) +  &
(-1.741252608270208e-09_Rkind) * tab_Pl(12) +  &
(-1.7969103343173652e-09_Rkind) * tab_Pl(13) +  &
(-8.523937203965875e-09_Rkind) * tab_Pl(14) +  &
(-7.974962571878466e-09_Rkind) * tab_Pl(15) +  &
(-1.6919700593984273e-09_Rkind) * tab_Pl(16) +  &
(3.974235083202602e-09_Rkind) * tab_Pl(17) +  &
(-9.48003074526905e-10_Rkind) * tab_Pl(18) +  &
(1.7306105360702442e-09_Rkind) * tab_Pl(19) +  &
(-2.086194212486486e-09_Rkind) * tab_Pl(20)
 
!HessZzRho, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 4.239814851940848e-08 
    Func(47) =  &
(-2.645261510215641e-09_Rkind) * tab_Pl(0) +  &
(1.8338550841035727e-08_Rkind) * tab_Pl(1) +  &
(-7.833174584935245e-09_Rkind) * tab_Pl(2) +  &
(-1.3433740495008297e-08_Rkind) * tab_Pl(3) +  &
(2.8282420333406643e-09_Rkind) * tab_Pl(4) +  &
(-2.0237708757212996e-08_Rkind) * tab_Pl(5) +  &
(8.539253060651236e-09_Rkind) * tab_Pl(6) +  &
(2.3041587122999563e-08_Rkind) * tab_Pl(7) +  &
(-4.8713632758526386e-09_Rkind) * tab_Pl(8) +  &
(4.523884897037846e-09_Rkind) * tab_Pl(9) +  &
(6.957920148295666e-09_Rkind) * tab_Pl(10) +  &
(-3.404024182014667e-08_Rkind) * tab_Pl(11) +  &
(3.52685207982774e-08_Rkind) * tab_Pl(12) +  &
(-1.1336252566527538e-08_Rkind) * tab_Pl(13) +  &
(1.4448274587544914e-08_Rkind) * tab_Pl(14) +  &
(-5.047099860312962e-09_Rkind) * tab_Pl(15) +  &
(-1.5481038903754423e-08_Rkind) * tab_Pl(16) +  &
(-2.6370264354874428e-08_Rkind) * tab_Pl(17) +  &
(1.1267323811718775e-08_Rkind) * tab_Pl(18) +  &
(5.778802966875766e-10_Rkind) * tab_Pl(19) +  &
(2.35616988046824e-08_Rkind) * tab_Pl(20)
 
!HessAnRh, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.1303973689202884e-05
    Func(48) =  &
(0.03648773136552473_Rkind) * tab_Pl(0) +  &
(0.0015298358205206003_Rkind) * tab_Pl(1) +  &
(-0.0015019964560929223_Rkind) * tab_Pl(2) +  &
(-0.00046972973795924195_Rkind) * tab_Pl(3) +  &
(0.00043328017520241337_Rkind) * tab_Pl(4) +  &
(-5.2931341975699005e-05_Rkind) * tab_Pl(5) +  &
(3.1254323922588414e-05_Rkind) * tab_Pl(6) +  &
(6.490363877455166e-05_Rkind) * tab_Pl(7) +  &
(7.282292216605888e-06_Rkind) * tab_Pl(8) +  &
(5.799224026116917e-06_Rkind) * tab_Pl(9) +  &
(-2.138102053229143e-06_Rkind) * tab_Pl(10) +  &
(7.2481443854401095e-06_Rkind) * tab_Pl(11) +  &
(1.0500458875060114e-06_Rkind) * tab_Pl(12) +  &
(7.280672346436011e-06_Rkind) * tab_Pl(13) +  &
(2.988295216138689e-06_Rkind) * tab_Pl(14) +  &
(-2.0351861741299624e-07_Rkind) * tab_Pl(15) +  &
(-1.6697064976867324e-06_Rkind) * tab_Pl(16) +  &
(-3.8133295264609436e-06_Rkind) * tab_Pl(17) +  &
(-6.992388053500056e-06_Rkind) * tab_Pl(18) +  &
(-7.75815970144272e-06_Rkind) * tab_Pl(19) +  &
(-6.014220585955108e-06_Rkind) * tab_Pl(20)
 
!HessAnRho, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 3.375163022319924e-05 
    Func(49) =  &
(0.009079622508368114_Rkind) * tab_Pl(0) +  &
(-0.008470698600657078_Rkind) * tab_Pl(1) +  &
(-0.0048718028225831845_Rkind) * tab_Pl(2) +  &
(0.0039651366070961045_Rkind) * tab_Pl(3) +  &
(4.3679907213885205e-05_Rkind) * tab_Pl(4) +  &
(-1.5991737977609906e-05_Rkind) * tab_Pl(5) +  &
(0.0008899514791406559_Rkind) * tab_Pl(6) +  &
(2.3985171966467807e-05_Rkind) * tab_Pl(7) +  &
(-0.0005546209557771005_Rkind) * tab_Pl(8) +  &
(-0.00010672262119454037_Rkind) * tab_Pl(9) +  &
(4.589418980707196e-05_Rkind) * tab_Pl(10) +  &
(-0.00011796433889484753_Rkind) * tab_Pl(11) +  &
(1.0229779564883014e-05_Rkind) * tab_Pl(12) +  &
(5.641187342757613e-05_Rkind) * tab_Pl(13) +  &
(9.75671332567052e-06_Rkind) * tab_Pl(14) +  &
(3.918484668450556e-05_Rkind) * tab_Pl(15) +  &
(1.2191265404227115e-05_Rkind) * tab_Pl(16) +  &
(3.130685021187893e-05_Rkind) * tab_Pl(17) +  &
(4.9349770646407355e-06_Rkind) * tab_Pl(18) +  &
(1.5063632858713411e-05_Rkind) * tab_Pl(19) +  &
(-3.435140973814467e-07_Rkind) * tab_Pl(20)
 
!HesszAn, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 1.2591154414422572e-08
    Func(50) =  &
(4.433090298754876e-10_Rkind) * tab_Pl(0) +  &
(-2.5806779759020707e-09_Rkind) * tab_Pl(1) +  &
(3.2969704329834975e-09_Rkind) * tab_Pl(2) +  &
(2.045079070301699e-09_Rkind) * tab_Pl(3) +  &
(-3.262770350239897e-09_Rkind) * tab_Pl(4) +  &
(2.2725419767846892e-10_Rkind) * tab_Pl(5) +  &
(-5.330412283610979e-09_Rkind) * tab_Pl(6) +  &
(3.734241119931934e-09_Rkind) * tab_Pl(7) +  &
(4.48823416468277e-09_Rkind) * tab_Pl(8) +  &
(-1.1050136939653337e-09_Rkind) * tab_Pl(9) +  &
(-2.0823809490488886e-09_Rkind) * tab_Pl(10) +  &
(-4.6216266459538056e-09_Rkind) * tab_Pl(11) +  &
(2.303730845295319e-09_Rkind) * tab_Pl(12) +  &
(4.454000125468102e-09_Rkind) * tab_Pl(13) +  &
(8.149334622690221e-10_Rkind) * tab_Pl(14) +  &
(-4.427779040971365e-09_Rkind) * tab_Pl(15) +  &
(-6.897261298009876e-09_Rkind) * tab_Pl(16) +  &
(2.350607447999494e-10_Rkind) * tab_Pl(17) +  &
(3.4150196189839695e-09_Rkind) * tab_Pl(18) +  &
(3.2192850245358144e-09_Rkind) * tab_Pl(19) +  &
(-7.984728841912268e-10_Rkind) * tab_Pl(20)
 
!HesszRh, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 3.8540360223882075e-05
    Func(51) =  &
(0.4504124265470616_Rkind) * tab_Pl(0) +  &
(-0.004043110626842922_Rkind) * tab_Pl(1) +  &
(0.00154844250974417_Rkind) * tab_Pl(2) +  &
(-0.0005096647123275679_Rkind) * tab_Pl(3) +  &
(0.00015681639317392213_Rkind) * tab_Pl(4) +  &
(0.00019917966005239527_Rkind) * tab_Pl(5) +  &
(-7.476958790844816e-05_Rkind) * tab_Pl(6) +  &
(0.00010998885451078226_Rkind) * tab_Pl(7) +  &
(-6.230076153203528e-05_Rkind) * tab_Pl(8) +  &
(5.32181669694862e-05_Rkind) * tab_Pl(9) +  &
(-4.806447758822834e-05_Rkind) * tab_Pl(10) +  &
(2.035506132224415e-05_Rkind) * tab_Pl(11) +  &
(-1.703024807303361e-05_Rkind) * tab_Pl(12) +  &
(-4.191488097110382e-06_Rkind) * tab_Pl(13) +  &
(9.463662616465133e-06_Rkind) * tab_Pl(14) +  &
(-4.9612886781792085e-06_Rkind) * tab_Pl(15) +  &
(8.98425771497299e-06_Rkind) * tab_Pl(16) +  &
(-5.327227849143387e-06_Rkind) * tab_Pl(17) +  &
(9.293289343861444e-06_Rkind) * tab_Pl(18) +  &
(-1.595939841430081e-06_Rkind) * tab_Pl(19) +  &
(9.869713973530872e-06_Rkind) * tab_Pl(20)
 
!HesszRho, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 0.0002856642430317802 
    Func(52) =  &
(1.7991045819848797_Rkind) * tab_Pl(0) +  &
(-0.7619517944051455_Rkind) * tab_Pl(1) +  &
(-0.18217275020314527_Rkind) * tab_Pl(2) +  &
(0.2509866078765607_Rkind) * tab_Pl(3) +  &
(-0.3416573918643211_Rkind) * tab_Pl(4) +  &
(0.0564560516386514_Rkind) * tab_Pl(5) +  &
(-0.025767846060756924_Rkind) * tab_Pl(6) +  &
(-0.004279092696086496_Rkind) * tab_Pl(7) +  &
(-0.03445489681206721_Rkind) * tab_Pl(8) +  &
(0.022433601596796856_Rkind) * tab_Pl(9) +  &
(-0.020542840650207207_Rkind) * tab_Pl(10) +  &
(0.005446920480584016_Rkind) * tab_Pl(11) +  &
(-0.0034387478709192107_Rkind) * tab_Pl(12) +  &
(-0.0006203507754641026_Rkind) * tab_Pl(13) +  &
(-0.00011012121816119454_Rkind) * tab_Pl(14) +  &
(-0.0007895496708456194_Rkind) * tab_Pl(15) +  &
(0.0008955709001659035_Rkind) * tab_Pl(16) +  &
(-0.00033847315824347496_Rkind) * tab_Pl(17) +  &
(0.0004640475657002285_Rkind) * tab_Pl(18) +  &
(-0.00024588262867177393_Rkind) * tab_Pl(19) +  &
(6.710639976490186e-05_Rkind) * tab_Pl(20)
 
!HessRhRho, grille ts = -0.99,-0.98...0.99, npoly=20, maxdiff: 2.5410056007119938e-05 
    Func(53) =  &
(-0.002945157297878704_Rkind) * tab_Pl(0) +  &
(0.0004681335014800168_Rkind) * tab_Pl(1) +  &
(0.0035953727999425405_Rkind) * tab_Pl(2) +  &
(0.00011525801545831526_Rkind) * tab_Pl(3) +  &
(-0.0010772520380865124_Rkind) * tab_Pl(4) +  &
(-5.598891008794806e-06_Rkind) * tab_Pl(5) +  &
(-0.0002691706031355797_Rkind) * tab_Pl(6) +  &
(-2.7643297261951077e-07_Rkind) * tab_Pl(7) +  &
(3.785074269605117e-05_Rkind) * tab_Pl(8) +  &
(6.082481148137796e-05_Rkind) * tab_Pl(9) +  &
(6.942731563933316e-05_Rkind) * tab_Pl(10) +  &
(-4.9652177842906056e-05_Rkind) * tab_Pl(11) +  &
(3.2128885375496755e-05_Rkind) * tab_Pl(12) +  &
(-3.47528797513614e-05_Rkind) * tab_Pl(13) +  &
(1.5703010361108163e-05_Rkind) * tab_Pl(14) +  &
(-5.691287396188627e-06_Rkind) * tab_Pl(15) +  &
(1.3645972038115859e-05_Rkind) * tab_Pl(16) +  &
(3.4769437956417124e-06_Rkind) * tab_Pl(17) +  &
(4.810674012632833e-06_Rkind) * tab_Pl(18) +  &
(-2.752598852829914e-07_Rkind) * tab_Pl(19) +  &
(-5.0723780382159494e-06_Rkind) * tab_Pl(20)

! ################ END OF FUNC #####################

    DO i=0,max_deg
      CALL dealloc_dnS(tab_Pl(i))
    END DO
    CALL dealloc_dnS(s)
    CALL dealloc_dnS(ts)

  END SUBROUTINE EvalFunc_QML_PH4Jo
  SUBROUTINE srho_new_TO_srho_old(sN,rhoN,sO,rhoO)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),         intent(in)    :: sN,rhoN
    TYPE (dnS_t),         intent(inout) :: sO,rhoO

    TYPE (dnS_t)                        :: th,R1,R2


    th = atan(HALF*sN/rhoN)*HALF+PI/FOUR
    R1 = rhoN / cos(th) *0.3748687959_Rkind ! R1 with the scalling
    R2 = rhoN / sin(th) *0.717978_Rkind     ! R2 with scalling

    th   = atan2(R1,R2)
    rhoO = R1*R2 / sqrt(R1*R1 + R2*R2)
    sO   = TWO*rhoO * tan(TWO*th - PI/TWO)

  END SUBROUTINE srho_new_TO_srho_old
  SUBROUTINE srho_old_TO_srho_new(sO,rhoO,sN,rhoN)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),         intent(in)    :: sO,rhoO
    TYPE (dnS_t),         intent(inout) :: sN,rhoN

    TYPE (dnS_t)                        :: th,R1,R2

    th = atan(HALF*sO/rhoO)*HALF+PI/FOUR
    R1 = rhoO / cos(th) /0.3748687959_Rkind ! new R1 without the scalling
    R2 = rhoO / sin(th) /0.717978_Rkind     ! new R2 without the scalling

    th   = atan2(R1,R2)
    rhoN = R1*R2 / sqrt(R1*R1 + R2*R2)
    sN   = TWO*rhoN * tan(TWO*th - PI/TWO)

  END SUBROUTINE srho_old_TO_srho_new
END MODULE QML_PH4Jo_m
