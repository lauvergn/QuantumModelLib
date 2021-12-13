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

MODULE Opt_m
!$ USE omp_lib
  USE QMLLib_NumParameters_m

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Set_Opt_param,QML_Opt

  TYPE :: Opt_t

    integer                           :: Max_it       = -1 ! it will be set-up after

    integer                           :: nb_neg       = 0 ! 0=>minimum, 1=>TS, 2=>top ...
    integer                           :: i_surf       = 0 ! on which surface the optimization is performed (default 1)

    integer                           :: hessian_type = 0 ! 1=> analytical hessian
    real (kind=Rkind)                 :: max_grad     = 0.000450_Rkind
    real (kind=Rkind)                 :: RMS_grad     = 0.000300_Rkind
    real (kind=Rkind)                 :: max_disp     = 0.001800_Rkind
    real (kind=Rkind)                 :: RMS_disp     = 0.001200_Rkind
  END TYPE Opt_t

CONTAINS

  SUBROUTINE Set_Opt_param(Opt_param,QModel,read_zmt)
  USE Model_m
  IMPLICIT NONE

    TYPE (Opt_t),       intent(inout)            :: Opt_param
    TYPE (Model_t),     intent(in)               :: QModel
    logical,            intent(in),    optional  :: read_zmt

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Set_Opt_param'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(read_zmt)) write(out_unitp,*) '   read_zmt',read_zmt
    flush(out_unitp)
  END IF

  CALL check_alloc_QM(QModel,name_sub)

  Opt_param = Opt_t(Max_it=20+QModel%ndim,nb_neg=0,i_surf=1,hessian_type=1)

  IF (debug) THEN
    CALL Write_Opt_param(Opt_param)
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

  END SUBROUTINE Set_Opt_param
  SUBROUTINE Write_Opt_param(Opt_param)
  IMPLICIT NONE

    TYPE (Opt_t),       intent(in)            :: Opt_param

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Write_Opt_param'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    write(out_unitp,*) ' BEGINNING ',name_sub

    write(out_unitp,*) ' Maxt_it      ',Opt_param%Max_it
    write(out_unitp,*) ' i_surf       ',Opt_param%i_surf
    write(out_unitp,*) ' nb_neg       ',Opt_param%nb_neg
    write(out_unitp,*) ' hessian_type ',Opt_param%hessian_type

    write(out_unitp,*) ' max_grad     ',Opt_param%max_grad
    write(out_unitp,*) ' RMS_grad     ',Opt_param%RMS_grad
    write(out_unitp,*) ' max_disp     ',Opt_param%max_disp
    write(out_unitp,*) ' RMS_disp     ',Opt_param%RMS_disp

    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)

  END SUBROUTINE Write_Opt_param
  SUBROUTINE QML_Opt(Q,QModel,Opt_param,Q0)
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q(:)
    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (Opt_t),       intent(in)               :: Opt_param

    real (kind=Rkind),  intent(in),    optional  :: Q0(:)


    TYPE (dnMat_t)                  :: PotVal
    integer                         :: it,iq
    real (kind=Rkind), allocatable  :: Qit(:)
    real (kind=Rkind), allocatable  :: mDQit(:)   ! -DelatQ
    real (kind=Rkind), allocatable  :: hess(:,:),grad(:)

    real (kind=Rkind)               :: max_grad
    real (kind=Rkind)               :: RMS_grad
    real (kind=Rkind)               :: max_disp
    real (kind=Rkind)               :: RMS_disp
    logical                         :: conv

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_Opt'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(Q0)) write(out_unitp,*) '   Q0',Q0
    CALL Write_Opt_param(Opt_param)
    flush(out_unitp)
  END IF

  IF (Opt_param%Max_it < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Opt_param is not initialized'
    STOP 'ERROR in QML_Opt: Opt_param is not initialized'
  END IF

  allocate(Qit(QModel%ndim))
  allocate(mDQit(QModel%ndim))
  allocate(grad(QModel%ndim))
  allocate(hess(QModel%ndim,QModel%ndim))

  IF (present(Q0)) THEN
    Qit(:) = Q0
  ELSE
    CALL get_Q0_Model(Qit,QModel,0)
  END IF

  write(out_unitp) '=================================================='
  DO it=0,Opt_param%Max_it

    CALL Eval_Pot(QModel,Qit,PotVal,nderiv=0)

    grad(:)   = PotVal%d1(Opt_param%i_surf,Opt_param%i_surf,:)
    hess(:,:) = PotVal%d2(Opt_param%i_surf,Opt_param%i_surf,:,:)

    CALL Linear_Sys(hess,grad,mDQit,QModel%ndim)

    max_grad = maxval(abs(grad))
    RMS_grad = sqrt(dot_product(grad,grad)/QModel%ndim)
    max_disp = maxval(abs(mDQit))
    RMS_disp = sqrt(dot_product(mDQit,mDQit)/QModel%ndim)


    write(out_unitp) '--------------------------------------------------'
    write(out_unitp) 'it,E',it,PotVal%d0(Opt_param%i_surf,Opt_param%i_surf)

    conv = (max_grad <= Opt_param%max_grad)
    write(out_unitp) 'max_grad,treshold',max_grad,Opt_param%max_grad,conv
    conv = (RMS_grad <= Opt_param%RMS_grad)
    write(out_unitp) 'RMS_grad,treshold',RMS_grad,Opt_param%RMS_grad,conv
    conv = (max_disp <= Opt_param%max_disp)
    write(out_unitp) 'max_disp,treshold',max_disp,Opt_param%max_disp,conv
    conv = (RMS_disp <= Opt_param%RMS_disp)
    write(out_unitp) 'RMS_disp,treshold',RMS_disp,Opt_param%RMS_disp,conv

    DO iq=1,QModel%ndim
      write(out_unitp) 'iq,Q(iq),grad(iq),DelatQ(iq)',iq,Qit(iq),grad(iq),-mDQit(iq)
    END DO

    conv = (max_grad <= Opt_param%max_grad) .AND.                               &
           (RMS_grad <= Opt_param%RMS_grad) .AND.                               &
           (RMS_disp <= Opt_param%RMS_disp) .AND.                               &
           (RMS_disp <= Opt_param%RMS_disp)

    Qit(:) = Qit-mDQit

    IF (conv) EXIT
  END DO
  IF (Opt_param%Max_it > 0) THEN
    write(out_unitp) 'Geometry optimization is converged?',conv
    write(out_unitp) 'Optimized geometry:'
    DO iq=1,QModel%ndim
      write(out_unitp) 'iq,Q(iq),',iq,Qit(iq)
    END DO
    Q(:) = Qit(:)
  ELSE
    write(out_unitp) 'No optimization (Max_it=0)'
    Q(:) = Q0(:)
  END IF
  write(out_unitp) '=================================================='

  deallocate(Qit)
  deallocate(mDQit)
  deallocate(grad)
  deallocate(hess)

  IF (debug) THEN
    write(out_unitp,*) '   it',it
    write(out_unitp,*) '   Q',Q
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_Opt

END MODULE Opt_m
