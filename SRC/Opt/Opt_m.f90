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
MODULE Opt_m
  !$ USE omp_lib
  USE QDUtil_NumParameters_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: QML_Opt_t,Init_QML_Opt,QML_Opt,Write_QML_Opt

  TYPE :: QML_Opt_t

    integer                           :: Max_it       = -1 ! it will be set-up after

    integer                           :: nb_neg       = 0 ! 0=>minimum, 1=>TS, 2=>top ...
    integer                           :: i_surf       = 1 ! on which surface the optimization is performed (default 1)

    integer, allocatable              :: list_act(:)      ! default all coordinates


    integer                           :: hessian_type = 1 ! 1=> analytical hessian

    real (kind=Rkind)                 :: Thresh_max_grad     = 0.000450_Rkind
    real (kind=Rkind)                 :: Thresh_RMS_grad     = 0.000300_Rkind
    real (kind=Rkind)                 :: Thresh_max_disp     = 0.001800_Rkind
    real (kind=Rkind)                 :: Thresh_RMS_disp     = 0.001200_Rkind
    real (kind=Rkind)                 :: Largest_disp        = 0.5_Rkind

  END TYPE QML_Opt_t

CONTAINS

  SUBROUTINE Init_QML_Opt(Opt_param,QModel,                                     &
                          read_param,param_file_name,nio_param_file,icv,list_act)

    USE QDUtil_m,         ONLY : TO_string
    USE Model_m
    IMPLICIT NONE

    TYPE (QML_Opt_t),   intent(inout)            :: Opt_param
    TYPE (Model_t),     intent(in)               :: QModel
    logical,            intent(in),    optional  :: read_param
    integer,            intent(in),    optional  :: nio_param_file
    character (len=*),  intent(in),    optional  :: param_file_name
    integer,            intent(in),    optional  :: icv
    integer,            intent(in),    optional  :: list_act(:)


    real (kind=Rkind), parameter      :: Thresh_max_grad     = 0.000450_Rkind
    real (kind=Rkind), parameter      :: Thresh_RMS_grad     = 0.000300_Rkind
    real (kind=Rkind), parameter      :: Thresh_max_disp     = 0.001800_Rkind
    real (kind=Rkind), parameter      :: Thresh_RMS_disp     = 0.001200_Rkind
    real (kind=Rkind), parameter      :: Largest_disp        = 0.5_Rkind


    integer                        :: icv_loc,nb_neg,i_surf,Max_it,hessian_type

    integer                        :: err_read,nio_loc,i
    logical                        :: read_param_loc
    character (len=:), allocatable :: param_file_name_loc

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Opt'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    write(out_unit,*) '   read_param      present?',present(read_param)
    write(out_unit,*) '   nio_param_file  present?',present(nio_param_file)
    write(out_unit,*) '   param_file_name present?',present(param_file_name)
    IF (present(param_file_name)) write(out_unit,*) '   param_file_name ',param_file_name
    flush(out_unit)
  END IF

  CALL check_alloc_QM(QModel,name_sub)

  IF (present(list_act)) THEN
    Opt_param%list_act = list_act
  ELSE IF (allocated(Opt_param%list_act)) THEN
    deallocate(Opt_param%list_act)
  END IF


  IF (present(icv)) THEN
    icv_loc = icv
  ELSE
    icv_loc = -1
  END IF

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
        param_file_name_loc = trim("input.dat")
        nio_loc = 99
      ELSE
        param_file_name_loc = trim(param_file_name)
        nio_loc = 99
      END IF
    ELSE
    nio_loc = in_unit
    END IF
  END IF
  IF (debug) THEN
    write(out_unit,*) '   read_param      ',read_param_loc
    write(out_unit,*) '   nio             ',nio_loc
    write(out_unit,*) '   allo param_file_name ',allocated(param_file_name_loc)
    IF (allocated(param_file_name_loc)) write(out_unit,*) '   param_file_name ',trim(param_file_name_loc)
    flush(out_unit)
  END IF

  IF (read_param_loc) THEN
    IF (nio_loc /= in_unit .AND. allocated(param_file_name_loc)) THEN
      open(unit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
    END IF
    IF (allocated(param_file_name_loc)) deallocate(param_file_name_loc)

    CALL Read_QML_Opt(Opt_param,QModel,nio_loc,icv_loc)

  ELSE
    Opt_param%Largest_disp  = abs(Largest_disp)
  END IF

  IF (icv_loc < 0)          icv_loc = 0

  IF (Opt_param%max_it < 0)   Opt_param%Max_it = (10+QModel%ndim)*(icv_loc+2)

  IF (Opt_param%i_surf < 1 .OR. Opt_param%i_surf > QModel%nsurf) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' i_surf',Opt_param%i_surf
    write(out_unit,*) ' i_surf is out-of-range ([1,',TO_string(QModel%nsurf),'])'
    write(out_unit,*) ' check your data!'
    write(out_unit,*)
    STOP ' ERROR in Init_QML_Opt: i_surf is out-of-range'
  END IF
  IF (Opt_param%nb_neg < 0) Opt_param%nb_neg = 0



  Opt_param%Thresh_max_grad = Thresh_max_grad/TEN**icv_loc
  Opt_param%Thresh_RMS_grad = Thresh_RMS_grad/TEN**icv_loc
  Opt_param%Thresh_max_disp = Thresh_max_disp/TEN**icv_loc
  Opt_param%Thresh_RMS_disp = Thresh_RMS_disp/TEN**icv_loc

  IF (allocated(Opt_param%list_act)) THEN
    Opt_param%list_act        = pack(Opt_param%list_act,mask=(Opt_param%list_act /= 0))

    IF (count(Opt_param%list_act /=0 ) == 0) THEN
      Opt_param%list_act = [(i,i=1,QModel%ndim)]
    END IF
  ELSE
    Opt_param%list_act = [(i,i=1,QModel%ndim)]
  END IF

  IF (debug) THEN
    CALL Write_QML_Opt(Opt_param)
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE Init_QML_Opt

  SUBROUTINE Read_QML_Opt(Opt_param,QModel,nio_param_file,icv_inout)
    USE QDUtil_m,         ONLY : TO_lowercase, TO_string
    USE Model_m
    IMPLICIT NONE

    TYPE (QML_Opt_t),   intent(inout)   :: Opt_param
    TYPE (Model_t),     intent(in)      :: QModel
    integer,            intent(in)      :: nio_param_file
    integer,            intent(inout)   :: icv_inout


    integer                        :: icv,nb_neg,i_surf,Max_it,hessian_type
    logical                        :: TS
    real(kind=Rkind)               :: Largest_disp
    character (len=Name_longlen)   :: hessian_method

    integer                        :: err_read,i
    integer,           allocatable :: list_act(:)      ! default all coordinates

    namelist /opt/ Max_it,nb_neg,i_surf,list_act,Largest_disp,icv,TS,hessian_method

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Read_QML_Opt'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) '   nio_param_file  ',nio_param_file
      flush(out_unit)
    END IF

    allocate(list_act(QModel%ndim))
    list_act(:) = 0
    IF (allocated(Opt_param%list_act)) THEN
      list_act(:) = Opt_param%list_act(1:size(Opt_param%list_act))
    END IF
  
    icv             = icv_inout ! to be able to change the convergence criteria
    Max_it          = -1
    Largest_disp    = 0.5_Rkind
  
    nb_neg          = -1 ! to be able to change the minimum, TS, ... optimization (see TS)
    TS              = .FALSE.
  
    i_surf          = 1
    hessian_method  = 'analytical'
  
  
    read(nio_param_file,nml=opt,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "Opt" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Opt'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Some parameter names of the namelist "Opt" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=Opt)
      STOP ' ERROR in Read_QML_Opt'
    END IF
  
    IF (TS .AND. nb_neg /= 1 .AND. nb_neg /= -1) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' TS=.TRUE. and nb_neg /= 1',TS,nb_neg
      write(out_unit,*) ' TS and nb_neg are not compatible'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Opt: TS and nb_neg are not compatible'
    END IF
    IF (TS)         nb_neg = 1
    IF (nb_neg < 0) nb_neg = 0
  
    hessian_method = TO_lowercase(hessian_method)
    SELECT CASE (hessian_method)
    CASE ('analytical','ana')
      hessian_type = 1
    CASE Default
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Wrong hessian_method'
      write(out_unit,*) ' The posibilities: "analytical" or "ana"'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Opt: Wrong hessian_method'
    END SELECT
  
    Opt_param = QML_Opt_t(Max_it=Max_it,nb_neg=nb_neg,i_surf=i_surf,              &
                          list_act=list_act,hessian_type=hessian_type,            &
                          Largest_disp=Largest_disp)
  
    icv_inout = icv
  
    IF (debug) THEN
      write(out_unit,*) ' icv_inout ',icv_inout
      CALL Write_QML_Opt(Opt_param)
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Read_QML_Opt

  SUBROUTINE Write_QML_Opt(Opt_param)
    IMPLICIT NONE

    TYPE (QML_Opt_t),       intent(in)            :: Opt_param

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Write_QML_Opt'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    write(out_unit,*) ' BEGINNING ',name_sub

    write(out_unit,*) ' Maxt_it         ',Opt_param%Max_it
    write(out_unit,*) ' i_surf          ',Opt_param%i_surf
    write(out_unit,*) ' nb_neg          ',Opt_param%nb_neg
    write(out_unit,*) ' hessian_type    ',Opt_param%hessian_type
    IF (allocated(Opt_param%list_act)) THEN
      write(out_unit,*) ' list_act        ',Opt_param%list_act
    END IF
    write(out_unit,*) ' Thresh_max_grad ',Opt_param%Thresh_max_grad
    write(out_unit,*) ' Thresh_RMS_grad ',Opt_param%Thresh_RMS_grad
    write(out_unit,*) ' Thresh_max_disp ',Opt_param%Thresh_max_disp
    write(out_unit,*) ' Thresh_RMS_disp ',Opt_param%Thresh_RMS_disp

    write(out_unit,*) ' Largest_disp    ',Opt_param%Largest_disp

    write(out_unit,*) ' END ',name_sub
    flush(out_unit)

  END SUBROUTINE Write_QML_Opt
  SUBROUTINE QML_Opt(Q,QModel,Opt_param,Q0)
    USE QDUtil_m,         ONLY : Identity_Mat, diagonalization, LinearSys_Solve
    USE ADdnSVM_m
    USE Model_m
    IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q(:)
    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_Opt_t),   intent(in)               :: Opt_param

    real (kind=Rkind),  intent(in),    optional  :: Q0(:)


    TYPE (dnMat_t)                  :: PotVal
    integer                         :: it,iq,i,nb_act
    real (kind=Rkind), allocatable  :: Qit(:),Qit_act(:),Q0_loc(:)
    real (kind=Rkind), allocatable  :: mDQit(:)   ! -DelatQ
    real (kind=Rkind), allocatable  :: Thess(:,:),hess(:,:),grad(:)
    real (kind=Rkind), allocatable  :: diag(:),Vec(:,:),tVec(:,:)

    real (kind=Rkind)               :: max_grad,RMS_grad
    real (kind=Rkind)               :: max_disp,RMS_disp,norm_disp
    logical                         :: conv

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_Opt'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      IF (present(Q0)) write(out_unit,*) '   Q0',Q0
      CALL Write_QML_Opt(Opt_param)
      CALL Write_Model(QModel)
      flush(out_unit)
    ELSE
      write(out_unit,*) '=================================================='
      write(out_unit,*) '=================================================='
      write(out_unit,*) '=== Optimization on the "',QModel%QM%pot_name,'" model.'
      write(out_unit,*) '=== model option:',QModel%QM%option
      write(out_unit,*) '=================================================='
      write(out_unit,*) '=================================================='
      flush(out_unit)
    END IF

    IF (Opt_param%Max_it < 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Opt_param is not initialized'
      STOP 'ERROR in QML_Opt: Opt_param is not initialized'
    END IF
  
    nb_act = count(Opt_param%list_act>0)
  
    allocate(Qit_act(nb_act))
    allocate(mDQit(nb_act))
    allocate(grad(nb_act))
    allocate(hess(nb_act,nb_act))
    allocate(vec(nb_act,nb_act))
    allocate(diag(nb_act))
  
    allocate(Q0_loc(QModel%ndim))
    IF (present(Q0)) THEN
      Q0_loc(:) = Q0
    ELSE
      CALL get_Q0_Model(Q0_loc,QModel,0)
    END IF
  
    allocate(Qit(QModel%ndim))
    Qit(:)     = Q0_loc
    Qit_act(:) = Qit(Opt_param%list_act)
  
  
    write(out_unit,*) '=================================================='
    DO it=0,Opt_param%Max_it
  
      CALL Eval_Pot(QModel,Qit,PotVal,nderiv=2)
      IF (debug) CALL Write_dnMat(PotVal,nio=out_unit)
  
      grad = PotVal%d1(Opt_param%i_surf,Opt_param%i_surf,Opt_param%list_act)
      hess = PotVal%d2(Opt_param%i_surf,Opt_param%i_surf,Opt_param%list_act,Opt_param%list_act)
  
      CALL diagonalization(hess,diag,Vec,nb_act)
      write(out_unit,*) 'diag',diag
  
      tvec = transpose(vec)
      IF (Opt_param%nb_neg == 0) THEN
        DO i=1,nb_act
          Vec(:,i) = Vec(:,i) * abs(diag(i))
        END DO
        hess = matmul(Vec,tVec)
      ELSE IF (count(diag < ZERO) /= Opt_param%nb_neg) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) '    Wrong number of negative hessian eigenvalues!'
        write(out_unit,*) '    Expected: ',Opt_param%nb_neg
        write(out_unit,*) '    it has: ',count(diag < ZERO)
        STOP 'ERROR in QML_Opt: Wrong number of negative hessian eigenvalues'
      END IF
  
  
      !write(out_unit,*) 'hess',hess
      !write(out_unit,*) 'hess?',matmul(Vec,tVec)
  
      mDQit = LinearSys_Solve(hess,grad)
  
  
      max_grad = maxval(abs(grad))
      RMS_grad = sqrt(dot_product(grad,grad)/nb_act)
      max_disp = maxval(abs(mDQit))
      RMS_disp = sqrt(dot_product(mDQit,mDQit)/nb_act)
  
      write(out_unit,*) '--------------------------------------------------'
      write(out_unit,*) 'it,E',it,PotVal%d0(Opt_param%i_surf,Opt_param%i_surf)
  
      conv = (max_grad <= Opt_param%Thresh_max_grad)
      write(out_unit,*) 'max_grad,treshold',max_grad,Opt_param%Thresh_max_grad,conv
      conv = (RMS_grad <= Opt_param%Thresh_RMS_grad)
      write(out_unit,*) 'RMS_grad,treshold',RMS_grad,Opt_param%Thresh_RMS_grad,conv
      conv = (max_disp <= Opt_param%Thresh_max_disp)
      write(out_unit,*) 'max_disp,treshold',max_disp,Opt_param%Thresh_max_disp,conv
      conv = (RMS_disp <= Opt_param%Thresh_RMS_disp)
      write(out_unit,*) 'RMS_disp,treshold',RMS_disp,Opt_param%Thresh_RMS_disp,conv
  
      norm_disp = sqrt(dot_product(mDQit,mDQit))
      IF (norm_disp > Opt_param%Largest_disp) THEN
        write(out_unit,*) ' The displacements are too large.'
        write(out_unit,*) ' The displacements:',-mDQit
        write(out_unit,*) ' Its norm:',norm_disp
        write(out_unit,*) '  => They are scaled by ',Opt_param%Largest_disp/norm_disp
        mDQit = mDQit * Opt_param%Largest_disp/norm_disp
      END IF
  
      DO iq=1,nb_act
        write(out_unit,*) 'iq,Q(iq),grad(iq),DelatQ(iq)',iq,                     &
                                                  Qit_act(iq),grad(iq),-mDQit(iq)
      END DO
  
      conv = (max_grad <= Opt_param%Thresh_max_grad) .AND.                               &
             (RMS_grad <= Opt_param%Thresh_RMS_grad) .AND.                               &
             (max_disp <= Opt_param%Thresh_max_disp) .AND.                               &
             (RMS_disp <= Opt_param%Thresh_RMS_disp)
  
      Qit_act(:) = Qit_act-mDQit
      CALL Qact_TO_Q(Qit_act,Qit,Opt_param%list_act)
  
      IF (conv) EXIT
    END DO
    IF (Opt_param%Max_it > 0) THEN
      write(out_unit,*) 'Geometry optimization is converged?',conv
      write(out_unit,*) 'Optimized geometry:'
      DO iq=1,QModel%ndim
        write(out_unit,*) 'iq,Q(iq),',iq,Qit(iq)
      END DO
      Q(:) = Qit(:)
    ELSE
      write(out_unit,*) 'No optimization (Max_it=0)'
      Q(:) = Q0_loc(:)
    END IF
    write(out_unit,*) '=================================================='
    flush(out_unit)
  
    deallocate(Qit_act)
    deallocate(Qit)
    deallocate(Q0_loc)
    deallocate(mDQit)
    deallocate(grad)
    deallocate(hess)
    deallocate(vec)
    deallocate(tvec)
    deallocate(diag)
  
    IF (debug) THEN
      write(out_unit,*) '   it',it
      write(out_unit,*) '   Q',Q
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    ELSE
      write(out_unit,*) '=================================================='
      write(out_unit,*) '=================================================='
      write(out_unit,*) '=== End of the optimization'
      write(out_unit,*) '=================================================='
      write(out_unit,*) '=================================================='
      flush(out_unit)
    END IF

  END SUBROUTINE QML_Opt

END MODULE Opt_m
