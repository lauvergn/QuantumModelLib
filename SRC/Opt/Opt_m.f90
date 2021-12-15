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
  PUBLIC :: Opt_t,Init_QML_Opt,QML_Opt,Write_QML_Opt

  TYPE :: Opt_t

    integer                           :: Max_it       = -1 ! it will be set-up after

    integer                           :: nb_neg       = 0 ! 0=>minimum, 1=>TS, 2=>top ...
    integer                           :: i_surf       = 0 ! on which surface the optimization is performed (default 1)

    integer                           :: hessian_type = 0 ! 1=> analytical hessian

    real (kind=Rkind)                 :: Thresh_max_grad     = 0.000450_Rkind
    real (kind=Rkind)                 :: Thresh_RMS_grad     = 0.000300_Rkind
    real (kind=Rkind)                 :: Thresh_max_disp     = 0.001800_Rkind
    real (kind=Rkind)                 :: Thresh_RMS_disp     = 0.001200_Rkind
    real (kind=Rkind)                 :: Largest_disp        = 0.5_Rkind

  END TYPE Opt_t

CONTAINS

  SUBROUTINE Init_QML_Opt(Opt_param,QModel,                                     &
                          read_param,param_file_name,nio_param_file)

  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

    TYPE (Opt_t),       intent(inout)            :: Opt_param
    TYPE (Model_t),     intent(in)               :: QModel
    logical,            intent(in),    optional  :: read_param
    integer,            intent(in),    optional  :: nio_param_file
    character (len=*),  intent(in),    optional  :: param_file_name


    integer                        :: icv,nb_neg,i_surf,Max_it,hessian_type
    logical                        :: TS
    real(kind=Rkind)               :: Largest_disp
    character (len=Name_longlen)   :: hessian_method

    integer                        :: err_read,nio_loc
    logical                        :: read_param_loc
    character (len=:), allocatable :: param_file_name_loc

    namelist /opt/ icv,nb_neg,i_surf,TS,hessian_method,Largest_disp,max_it

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Opt'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(read_param)) write(out_unitp,*) '   read_param',read_param
    flush(out_unitp)
  END IF

  CALL check_alloc_QM(QModel,name_sub)

  IF (present(read_param)) THEN
    read_param_loc = read_param
  ELSE
    read_param_loc = .FALSE.
  END IF
  IF (present(nio_param_file)) THEN
    nio_loc = nio_param_file
  ELSE
    IF (present(param_file_name)) THEN
      IF (len_trim(param_file_name) == 0) THEN
        param_file_name_loc = strdup("input.dat")
        nio_loc = 99
      ELSE
        param_file_name_loc = strdup(param_file_name)
        nio_loc = 99
      END IF
    ELSE
    nio_loc = in_unitp
    END IF
  END IF



  icv             = -1 ! to be able to change the convergence criteria
  Max_it          = -1
  Largest_disp    = 0.5_Rkind

  nb_neg          = -1 ! to be able to change the minimum, TS, ... optimization (see TS)
  TS              = .FALSE.

  i_surf          = 1
  hessian_method  = 'analytical'


  IF (read_param_loc) THEN
    IF (nio_loc /= in_unitp) THEN
      open(unit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
    END IF
    IF (allocated(param_file_name_loc)) deallocate(param_file_name_loc)

    read(nio_loc,nml=opt,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Opt" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Init_QML_Opt'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Some parameter names of the namelist "Opt" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Opt)
      STOP ' ERROR in Init_QML_Opt'
    END IF

  END IF

  IF (icv < 0)      icv    = 0
  IF (max_it < 0)   Max_it = (10+QModel%ndim)*(icv+2)
  Largest_disp             = abs(Largest_disp)

  IF (i_surf < 0 .OR. i_surf > QModel%nsurf) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' i_surf',i_surf
    write(out_unitp,*) ' i_surf is out-of-range ([1,',int_TO_char(QModel%nsurf),'])'
    write(out_unitp,*) ' check your data!'
    write(out_unitp,*)
    STOP ' ERROR in Set_Opt_param: i_surf is out-of-range'
  END IF
  IF (TS .AND. nb_neg /= 1 .AND. nb_neg /= -1) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' TS=.TRUE. and nb_neg /= 1',TS,nb_neg
    write(out_unitp,*) ' TS and nb_neg are not compatible'
    write(out_unitp,*) ' check your data!'
    write(out_unitp,*)
    STOP ' ERROR in Init_QML_Opt: TS and nb_neg are not compatible'
  END IF
  IF (TS)         nb_neg = 1
  IF (nb_neg < 0) nb_neg = 0

  CALL string_uppercase_TO_lowercase(hessian_method)
  SELECT CASE (hessian_method)
  CASE ('analytical','ana')
    hessian_type = 1
  CASE Default
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong hessian_method'
    write(out_unitp,*) ' The posibilities: "analytical" or "ana"'
    write(out_unitp,*) ' check your data!'
    write(out_unitp,*)
    STOP ' ERROR in Init_QML_Opt: Wrong hessian_method'
  END SELECT

  Opt_param = Opt_t(Max_it=Max_it,nb_neg=nb_neg,i_surf=i_surf,                  &
                    hessian_type=hessian_type,Largest_disp=Largest_disp)

  Opt_param%Thresh_max_grad = Opt_param%Thresh_max_grad/TEN**icv
  Opt_param%Thresh_RMS_grad = Opt_param%Thresh_RMS_grad/TEN**icv
  Opt_param%Thresh_max_disp = Opt_param%Thresh_max_disp/TEN**icv
  Opt_param%Thresh_RMS_disp = Opt_param%Thresh_RMS_disp/TEN**icv

  IF (debug) THEN
    CALL Write_QML_Opt(Opt_param)
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

  END SUBROUTINE Init_QML_Opt
  SUBROUTINE Write_QML_Opt(Opt_param)
  IMPLICIT NONE

    TYPE (Opt_t),       intent(in)            :: Opt_param

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Write_QML_Opt'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    write(out_unitp,*) ' BEGINNING ',name_sub

    write(out_unitp,*) ' Maxt_it         ',Opt_param%Max_it
    write(out_unitp,*) ' i_surf          ',Opt_param%i_surf
    write(out_unitp,*) ' nb_neg          ',Opt_param%nb_neg
    write(out_unitp,*) ' hessian_type    ',Opt_param%hessian_type

    write(out_unitp,*) ' Thresh_max_grad ',Opt_param%Thresh_max_grad
    write(out_unitp,*) ' Thresh_RMS_grad ',Opt_param%Thresh_RMS_grad
    write(out_unitp,*) ' Thresh_max_disp ',Opt_param%Thresh_max_disp
    write(out_unitp,*) ' Thresh_RMS_disp ',Opt_param%Thresh_RMS_disp

    write(out_unitp,*) ' Largest_disp    ',Opt_param%Largest_disp

    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)

  END SUBROUTINE Write_QML_Opt
  SUBROUTINE QML_Opt(Q,QModel,Opt_param,Q0)
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q(:)
    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (Opt_t),       intent(in)               :: Opt_param

    real (kind=Rkind),  intent(in),    optional  :: Q0(:)


    TYPE (dnMat_t)                  :: PotVal
    integer                         :: it,iq,i
    real (kind=Rkind), allocatable  :: Qit(:)
    real (kind=Rkind), allocatable  :: mDQit(:)   ! -DelatQ
    real (kind=Rkind), allocatable  :: Thess(:,:),hess(:,:),grad(:)
    real (kind=Rkind), allocatable  :: diag(:),Vec(:,:),tVec(:,:)

    real (kind=Rkind)               :: max_grad,RMS_grad
    real (kind=Rkind)               :: max_disp,RMS_disp,norm_disp
    logical                         :: conv

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_Opt'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(Q0)) write(out_unitp,*) '   Q0',Q0
    CALL Write_QML_Opt(Opt_param)
    CALL Write_Model(QModel)
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
  allocate(vec(QModel%ndim,QModel%ndim))
  allocate(diag(QModel%ndim))

  IF (present(Q0)) THEN
    Qit(:) = Q0
  ELSE
    CALL get_Q0_Model(Qit,QModel,0)
  END IF

  write(out_unitp,*) '=================================================='
  DO it=0,Opt_param%Max_it

    CALL Eval_Pot(QModel,Qit,PotVal,nderiv=2)

    grad(:)   = PotVal%d1(Opt_param%i_surf,Opt_param%i_surf,:)
    hess(:,:) = PotVal%d2(Opt_param%i_surf,Opt_param%i_surf,:,:)

    CALL diagonalization(hess,diag,Vec,QModel%ndim)
    write(out_unitp,*) 'diag',diag

    tvec = transpose(vec)
    IF (Opt_param%nb_neg == 0) THEN
      DO i=1,QModel%ndim
        Vec(:,i) = Vec(:,i) * abs(diag(i))
      END DO
      hess = matmul(Vec,tVec)
    ELSE IF (count(diag < ZERO) /= Opt_param%nb_neg) THEN
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) '    Wrong number of negative hessian eigenvalues!'
      write(out_unitp,*) '    Expected: ',Opt_param%nb_neg
      write(out_unitp,*) '    it has: ',count(diag < ZERO)
      STOP 'ERROR in QML_Opt: Wrong number of negative hessian eigenvalues'
    END IF


    !write(out_unitp,*) 'hess',hess
    !write(out_unitp,*) 'hess?',matmul(Vec,tVec)

    CALL Linear_Sys(hess,grad,mDQit,QModel%ndim)



    max_grad = maxval(abs(grad))
    RMS_grad = sqrt(dot_product(grad,grad)/QModel%ndim)
    max_disp = maxval(abs(mDQit))
    RMS_disp = sqrt(dot_product(mDQit,mDQit)/QModel%ndim)

    write(out_unitp,*) '--------------------------------------------------'
    write(out_unitp,*) 'it,E',it,PotVal%d0(Opt_param%i_surf,Opt_param%i_surf)

    conv = (max_grad <= Opt_param%Thresh_max_grad)
    write(out_unitp,*) 'max_grad,treshold',max_grad,Opt_param%Thresh_max_grad,conv
    conv = (RMS_grad <= Opt_param%Thresh_RMS_grad)
    write(out_unitp,*) 'RMS_grad,treshold',RMS_grad,Opt_param%Thresh_RMS_grad,conv
    conv = (max_disp <= Opt_param%Thresh_max_disp)
    write(out_unitp,*) 'max_disp,treshold',max_disp,Opt_param%Thresh_max_disp,conv
    conv = (RMS_disp <= Opt_param%Thresh_RMS_disp)
    write(out_unitp,*) 'RMS_disp,treshold',RMS_disp,Opt_param%Thresh_RMS_disp,conv

    norm_disp = sqrt(dot_product(mDQit,mDQit))
    IF (norm_disp > Opt_param%Largest_disp) THEN
      write(out_unitp,*) ' The displacements are too large.'
      write(out_unitp,*) ' The displacements:',-mDQit
      write(out_unitp,*) ' Its norm:',norm_disp
      write(out_unitp,*) '  => They are scaled by ',Opt_param%Largest_disp/norm_disp
      mDQit = mDQit * Opt_param%Largest_disp/norm_disp
    END IF

    DO iq=1,QModel%ndim
      write(out_unitp,*) 'iq,Q(iq),grad(iq),DelatQ(iq)',iq,Qit(iq),grad(iq),-mDQit(iq)
    END DO

    conv = (max_grad <= Opt_param%Thresh_max_grad) .AND.                               &
           (RMS_grad <= Opt_param%Thresh_RMS_grad) .AND.                               &
           (max_disp <= Opt_param%Thresh_max_disp) .AND.                               &
           (RMS_disp <= Opt_param%Thresh_RMS_disp)

    Qit(:) = Qit-mDQit

    IF (conv) EXIT
  END DO
  IF (Opt_param%Max_it > 0) THEN
    write(out_unitp,*) 'Geometry optimization is converged?',conv
    write(out_unitp,*) 'Optimized geometry:'
    DO iq=1,QModel%ndim
      write(out_unitp,*) 'iq,Q(iq),',iq,Qit(iq)
    END DO
    Q(:) = Qit(:)
  ELSE
    write(out_unitp,*) 'No optimization (Max_it=0)'
    Q(:) = Q0(:)
  END IF
  write(out_unitp,*) '=================================================='

  deallocate(Qit)
  deallocate(mDQit)
  deallocate(grad)
  deallocate(hess)
  deallocate(vec)
  deallocate(tvec)
  deallocate(diag)

  IF (debug) THEN
    write(out_unitp,*) '   it',it
    write(out_unitp,*) '   Q',Q
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_Opt

END MODULE Opt_m
