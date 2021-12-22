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

MODULE IRC_m
!$ USE omp_lib
  USE QMLLib_NumParameters_m
  USE Opt_m

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: QML_IRC_t,Init_QML_IRC,QML_IRC,Write_QML_IRC

  TYPE, EXTENDS (QML_Opt_t) :: QML_IRC_t

    integer                           :: IRC_Max_it       = -1 ! it will be set-up after
    real (kind=Rkind)                 :: Delta_s          = ONETENTH**2

    real (kind=Rkind)                 :: Ene_TS           = ZERO
    real (kind=Rkind), allocatable    :: Grad_QactTS(:)
    real (kind=Rkind), allocatable    :: QTS(:)
    real (kind=Rkind), allocatable    :: QactTS(:)
    real (kind=Rkind), allocatable    :: EigenVec_QactTS(:)

  END TYPE QML_IRC_t

  TYPE :: QML_IRC_at_s

    real (kind=Rkind)                    :: s
    real (kind=Rkind)                    :: Ene
    real (kind=Rkind), allocatable       :: Qact(:),DeltaQact(:),Grad_Qact(:)

    real (kind=Rkind), allocatable       :: QTS(:)
    real (kind=Rkind), allocatable       :: DeltaQact_TS(:)

  END TYPE QML_IRC_at_s

CONTAINS

  SUBROUTINE Init_QML_IRC(IRC_p,QModel,                                         &
                          read_param,param_file_name,nio_param_file)

  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

    TYPE (QML_IRC_t),   intent(inout)            :: IRC_p
    TYPE (Model_t),     intent(in)               :: QModel
    logical,            intent(in),    optional  :: read_param
    integer,            intent(in),    optional  :: nio_param_file
    character (len=*),  intent(in),    optional  :: param_file_name



    integer                        :: err_read,nio_loc,i
    logical                        :: read_param_loc
    character (len=:), allocatable :: param_file_name_loc

    integer                        :: Max_it
    real (kind=Rkind)              :: Delta_s

    namelist /IRC/ max_it,Delta_s

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_IRC'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    write(out_unitp,*) '   read_param      present?',present(read_param)
    write(out_unitp,*) '   nio_param_file  present?',present(nio_param_file)
    write(out_unitp,*) '   param_file_name present?',present(param_file_name)
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

  IF (debug) THEN
    write(out_unitp,*) '   read_param      ',read_param_loc
    write(out_unitp,*) '   nio             ',nio_loc
    write(out_unitp,*) '   param_file_name ',strdup(param_file_name_loc)
    flush(out_unitp)
  END IF


  Max_it          = -1
  Delta_s         = ONETENTH**2

  IF (read_param_loc) THEN
    IF (nio_loc /= in_unitp) THEN
      open(unit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
    END IF


    CALL Init_QML_Opt(IRC_p%QML_Opt_t,QModel,read_param=.TRUE.,nio_param_file=nio_loc)

    IF (allocated(param_file_name_loc)) deallocate(param_file_name_loc)


    read(nio_loc,nml=IRC,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "IRC" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Init_QML_IRC'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Some parameter names of the namelist "IRC" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=IRC)
      STOP ' ERROR in Init_QML_IRC'
    END IF

  END IF

  IF (max_it < 0)   Max_it = 10

  IRC_p = QML_IRC_t(IRC_Max_it=Max_it,Delta_s=Delta_s,QML_Opt_t=IRC_p%QML_Opt_t)

  IF (debug) THEN
    CALL Write_QML_IRC(IRC_p)
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

  END SUBROUTINE Init_QML_IRC
  SUBROUTINE Write_QML_IRC(IRC_p)
  IMPLICIT NONE

    TYPE (QML_IRC_t),       intent(in)            :: IRC_p

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Write_QML_IRC'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    write(out_unitp,*) ' BEGINNING ',name_sub

    CALL Write_QML_Opt(IRC_p%QML_Opt_t)

    write(out_unitp,*) ' IRC_Maxt_it     ',IRC_p%IRC_Max_it
    write(out_unitp,*) ' Delta_s         ',IRC_p%Delta_s

    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)

  END SUBROUTINE Write_QML_IRC
  SUBROUTINE QML_IRC(Q,QModel,IRC_p,Q0)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q(:)
    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(inout)            :: IRC_p

    real (kind=Rkind),  intent(in),    optional  :: Q0(:)


    TYPE (QML_Opt_t)                :: Opt_p
    integer                         :: it
    real (kind=Rkind), allocatable  :: QactOld(:),QactNew(:)
    real (kind=Rkind)               :: s,Ene_AT_s
    logical                         :: forward

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(Q0)) write(out_unitp,*) '   Q0',Q0
    CALL Write_QML_IRC(IRC_p)
    CALL Write_Model(QModel)
    flush(out_unitp)
  END IF

  IF (IRC_p%Max_it < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' IRC_p is not initialized'
    STOP 'ERROR in QML_IRC: IRC_p is not initialized'
  END IF

  ! initialization of IRC_at_s%QTS
  allocate(IRC_p%QTS(QModel%ndim))
  IF (present(Q0)) THEN
    IRC_p%QTS(:) = Q0
  ELSE
    CALL get_Q0_Model(IRC_p%QTS,QModel,0)
  END IF

  CALL QML_Opt(IRC_p%QTS,QModel,IRC_p%QML_Opt_t,Q0=IRC_p%QTS)
  !CALL Write_RVec(IRC_at_s%QTS,    out_unitp,3,name_info='Qit')

  !first point + check if the geometry is a TS (s=0)
  CALL QML_IRC_at_TS(IRC_p,QModel)


  s       = ZERO
  QactOld = IRC_p%QactTS
  allocate(QactNew(size(QactOld)))
  forward = .TRUE.

  DO it=0,IRC_p%IRC_Max_it-1

    CALL QML_IRC_mEuler(s,QactOld,QactNew,Ene_AT_s,QModel,IRC_p,m=100,forward=forward)

    write(out_unitp,*) 's,Qact,E',s,QactOld,Ene_AT_s
    flush(out_unitp)

    s       = s    +         IRC_p%Delta_s
    QactOld = QactNew

  END DO
  CALL QML_IRC_fcn(s,QactOld,QactNew,Ene_AT_s,QModel,IRC_p,forward)
  write(out_unitp,*) 's,Qact,E',s,QactOld,Ene_AT_s
  flush(out_unitp)

  s       = ZERO
  QactOld = IRC_p%QactTS
  forward = .FALSE.

  DO it=0,IRC_p%IRC_Max_it-1

    CALL QML_IRC_mEuler(s,QactOld,QactNew,Ene_AT_s,QModel,IRC_p,m=100,forward=forward)

    write(out_unitp,*) 's,Qact,E',s,QactOld,Ene_AT_s
    flush(out_unitp)

    s       = s    -         IRC_p%Delta_s
    QactOld = QactNew

  END DO
  CALL QML_IRC_fcn(s,QactOld,QactNew,Ene_AT_s,QModel,IRC_p,forward)
  write(out_unitp,*) 's,Qact,E',s,QactOld,Ene_AT_s
  flush(out_unitp)

STOP

  IF (debug) THEN
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

  END SUBROUTINE QML_IRC

  SUBROUTINE QML_IRC_mEuler(s,QactOld,QactNew,Ene_AT_s,QModel,IRC_p,m,forward)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(in)               :: QactOld(:)
    real (kind=Rkind),  intent(inout)            :: QactNew(:)
    real (kind=Rkind),  intent(inout)            :: Ene_AT_s

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p
    logical,            intent(in)               :: forward
    integer,            intent(in)               :: m


    integer                         :: it
    real (kind=Rkind), allocatable  :: Qact(:),dQact(:)
    real (kind=Rkind)               :: s_loc,Ene_loc,Delta_s

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_mEuler'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    write(out_unitp,*) '   s      ',s
    write(out_unitp,*) '   QactOld',QactOld
    flush(out_unitp)
  END IF

  IF (m < 1) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' m < 1',m
    write(out_unitp,*) ' m MUST be larger than 0'
    STOP 'ERROR in QML_IRC_mEuler: m < 1'
  END IF

  s_loc = s
  Qact = QactOld
  allocate(dQact(size(QactOld)))
  Delta_s = IRC_p%Delta_s/m

  IF (forward) THEN
    DO it=1,m

      CALL QML_IRC_fcn(s_loc,Qact,dQact,Ene_loc,QModel,IRC_p,forward=.TRUE.)

      IF (it == 1) Ene_AT_s = Ene_loc

      IF (debug) write(out_unitp,*) 's,Qact,E',s_loc,Qact,Ene_loc
      IF (debug) flush(out_unitp)

      s_loc = s_loc +         Delta_s
      Qact  = Qact  + dQact * Delta_s

    END DO
  ELSE
    DO it=1,m

      CALL QML_IRC_fcn(s_loc,Qact,dQact,Ene_loc,QModel,IRC_p,forward=.FALSE.)

      IF (it == 1) Ene_AT_s = Ene_loc

      IF (debug) write(out_unitp,*) 's,Qact,E',s_loc,Qact,Ene_loc
      IF (debug) flush(out_unitp)

      s_loc = s_loc -         Delta_s
      Qact  = Qact  + dQact * Delta_s

    END DO
  END IF

  QactNew = Qact


  IF (debug) THEN
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_IRC_mEuler

SUBROUTINE QML_IRC_fcn(s,Qact,dQact,Ene_AT_s,QModel,IRC_p,forward)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(in)               :: Qact(:)
    real (kind=Rkind),  intent(inout)            :: dQact(:)
    real (kind=Rkind),  intent(inout)            :: Ene_AT_s

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p
    logical,            intent(in)               :: forward

    TYPE (dnMat_t)                  :: PotVal
    real (kind=Rkind), allocatable  :: Qit(:)

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_fcn'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    flush(out_unitp)
  END IF

  IF (s == ZERO) THEN
    IF (forward) THEN
       dQact =  IRC_p%EigenVec_QactTS
    ELSE
       dQact = -IRC_p%EigenVec_QactTS
    END IF
    Ene_AT_s = IRC_p%Ene_TS
  ELSE
    Qit = IRC_p%QTS

    CALL Qact_TO_Q(Qact,Qit,IRC_p%list_act)

    !CALL Write_RVec(Qit,    out_unitp,3,name_info='Qit')
    CALL Eval_Pot(QModel,Qit,PotVal,nderiv=1)

    Ene_AT_s = PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)

    dQact    = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act)
    dQact    = -dQact/norm2(dQact)

    deallocate(Qit)
  END IF

  IF (debug) THEN
    write(out_unitp,*) 's,Ene    ',s,Ene_AT_s
    write(out_unitp,*) 'Qact     ',Qact
    write(out_unitp,*) 'dQact    ',dQact
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_IRC_fcn

  SUBROUTINE QML_IRC_at_TS(IRC_p,QModel)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

    TYPE (QML_IRC_t),   intent(inout)            :: IRC_p
    TYPE (Model_t),     intent(inout)            :: QModel

    TYPE (dnMat_t)                  :: PotVal
    integer                         :: nb_act
    real (kind=Rkind), allocatable  :: hess(:,:)
    real (kind=Rkind), allocatable  :: diag(:),Vec(:,:)

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_at_TS'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    flush(out_unitp)
  END IF

  IF (IRC_p%Max_it < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' IRC_p is not initialized'
    STOP 'ERROR in QML_IRC_at_TS: IRC_p is not initialized'
  END IF

  nb_act = size(IRC_p%list_act)

  allocate(hess(nb_act,nb_act))
  allocate(vec(nb_act,nb_act))
  allocate(diag(nb_act))


  !first check if the geometry in Q0 is a TS
  CALL Eval_Pot(QModel,IRC_p%QTS,PotVal,nderiv=2)
  hess              = PotVal%d2(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act,IRC_p%list_act)
  IRC_p%Grad_QactTS = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act)

  CALL diagonalization(hess,diag,Vec,nb_act,sort=1)
  write(out_unitp,*) 'grad',IRC_p%Grad_QactTS
  write(out_unitp,*) 'diag',diag
  write(out_unitp,*) 'TS',(count(diag < ZERO) == 1)
  IF ((count(diag < ZERO) /= 1)) STOP 'STOP in QML_IRC_at_TS: Not a TS at s=0'

  !s=0 (TS)
  IRC_p%Ene_TS          = PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)
  IRC_p%QactTS          = IRC_p%QTS(IRC_p%list_act)
  IRC_p%EigenVec_QactTS = Vec(:,1)

  deallocate(hess)
  deallocate(vec)
  deallocate(diag)

  IF (debug) THEN
    write(out_unitp,*) 's,Ene          ',ZERO,IRC_p%Ene_TS
    write(out_unitp,*) 'Qact           ',IRC_p%QactTS
    write(out_unitp,*) 'EigenVec_QactTS',IRC_p%EigenVec_QactTS
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_IRC_at_TS
  SUBROUTINE QML_IRC_v2(Q,QModel,IRC_p,Q0)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q(:)
    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(inout)            :: IRC_p

    real (kind=Rkind),  intent(in),    optional  :: Q0(:)


    TYPE (QML_Opt_t)                :: Opt_p
    integer                         :: it
    real (kind=Rkind), allocatable  :: Qact(:),dQact(:)
    real (kind=Rkind)               :: s,Ene_AT_s

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_v2'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(Q0)) write(out_unitp,*) '   Q0',Q0
    CALL Write_QML_IRC(IRC_p)
    CALL Write_Model(QModel)
    flush(out_unitp)
  END IF

  IF (IRC_p%Max_it < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' IRC_p is not initialized'
    STOP 'ERROR in QML_IRC_v2: IRC_p is not initialized'
  END IF

  ! initialization of IRC_at_s%QTS
  allocate(IRC_p%QTS(QModel%ndim))
  IF (present(Q0)) THEN
    IRC_p%QTS(:) = Q0
  ELSE
    CALL get_Q0_Model(IRC_p%QTS,QModel,0)
  END IF

  CALL QML_Opt(IRC_p%QTS,QModel,IRC_p%QML_Opt_t,Q0=IRC_p%QTS)
  !CALL Write_RVec(IRC_at_s%QTS,    out_unitp,3,name_info='Qit')

  !first point + check if the geometry is a TS (s=0)
  CALL QML_IRC_at_TS(IRC_p,QModel)
  s    = ZERO
  Qact = IRC_p%QactTS
  allocate(dQact(size(Qact)))

  DO it=0,IRC_p%IRC_Max_it

    CALL QML_IRC_fcn(s,Qact,dQact,Ene_AT_s,QModel,IRC_p,forward=.TRUE.)

    write(out_unitp,*) 's,Qact,E',s,Qact,Ene_AT_s
    flush(out_unitp)

    s    = s    +         IRC_p%Delta_s
    Qact = Qact + dQact * IRC_p%Delta_s

  END DO

  s    = ZERO
  Qact = IRC_p%QactTS

  DO it=0,IRC_p%IRC_Max_it

    CALL QML_IRC_fcn(s,Qact,dQact,Ene_AT_s,QModel,IRC_p,forward=.FALSE.)

    write(out_unitp,*) 's,Qact,E',s,Qact,Ene_AT_s
    flush(out_unitp)

    s    = s    -         IRC_p%Delta_s
    Qact = Qact + dQact * IRC_p%Delta_s

  END DO

STOP

  IF (debug) THEN
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_IRC_v2
  SUBROUTINE QML_IRC_v1(Q,QModel,IRC_p,Q0)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q(:)
    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p

    real (kind=Rkind),  intent(in),    optional  :: Q0(:)



    TYPE(QML_IRC_at_s)              :: IRC_at_s,IRC_at_s_old,IRC_at_TS

    TYPE (QML_Opt_t)                :: Opt_p
    TYPE (dnMat_t)                  :: PotVal
    integer                         :: it

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_v1'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(Q0)) write(out_unitp,*) '   Q0',Q0
    CALL Write_QML_IRC(IRC_p)
    CALL Write_Model(QModel)
    flush(out_unitp)
  END IF

  IF (IRC_p%Max_it < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' IRC_p is not initialized'
    STOP 'ERROR in QML_IRC: IRC_p is not initialized'
  END IF

  ! initialization of IRC_at_s%QTS
  allocate(IRC_at_TS%QTS(QModel%ndim))
  IF (present(Q0)) THEN
    IRC_at_TS%QTS(:) = Q0
  ELSE
    CALL get_Q0_Model(IRC_at_TS%QTS,QModel,0)
  END IF

  CALL QML_Opt(IRC_at_TS%QTS,QModel,IRC_p%QML_Opt_t,Q0=IRC_at_TS%QTS)
  !CALL Write_RVec(IRC_at_s%QTS,    out_unitp,3,name_info='Qit')

  !first point + check if the geometry is a TS (s=0)
  CALL QML_IRC_FirstPoint(IRC_at_TS,QModel,IRC_p)
  write(out_unitp,*) 's,Qit,|grad|,E',IRC_at_TS%s,IRC_at_TS%Qact,               &
                      norm2(IRC_at_TS%Grad_Qact),IRC_at_TS%Ene

  IRC_at_s     = IRC_at_TS
  IRC_at_s_old = IRC_at_s

  DO it=1,IRC_p%IRC_Max_it

    CALL QML_IRC_NewPoint(IRC_at_s,IRC_at_s_old,QModel,IRC_p,forward=.TRUE.)
    write(out_unitp,*) 's,Qit,|grad|,E',IRC_at_s%s,IRC_at_s%Qact,               &
                      norm2(IRC_at_s%Grad_Qact),IRC_at_s%Ene
    flush(out_unitp)

    IRC_at_s_old = IRC_at_s

  END DO

  IRC_at_s     = IRC_at_TS
  IRC_at_s%DeltaQact = -IRC_at_s%DeltaQact_TS
  IRC_at_s_old = IRC_at_s

  DO it=1,IRC_p%IRC_Max_it

    CALL QML_IRC_NewPoint(IRC_at_s,IRC_at_s_old,QModel,IRC_p,forward=.FALSE.)
    write(out_unitp,*) 's,Qit,|grad|,E',IRC_at_s%s,IRC_at_s%Qact,               &
                      norm2(IRC_at_s%Grad_Qact),IRC_at_s%Ene
    flush(out_unitp)

    IRC_at_s_old = IRC_at_s

  END DO

  IF (debug) THEN
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_IRC_v1

  SUBROUTINE QML_IRC_v0(Q,QModel,IRC_p,Q0)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q(:)
    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p

    real (kind=Rkind),  intent(in),    optional  :: Q0(:)

    TYPE (QML_Opt_t)                :: Opt_p
    TYPE (dnMat_t)                  :: PotVal
    integer                         :: it,iq,i,nb_act
    real (kind=Rkind), allocatable  :: Qit(:),QTS(:)
    !real (kind=Rkind), allocatable  :: mDQit(:)   ! -DelatQ
    real (kind=Rkind), allocatable  :: Qit_act(:),QTS_act(:)
    real (kind=Rkind), allocatable  :: mDQit_act(:)   ! -DelatQ
    real (kind=Rkind), allocatable  :: Thess(:,:),hess(:,:),grad(:)
    real (kind=Rkind), allocatable  :: diag(:),Vec(:,:),tVec(:,:)

    real (kind=Rkind)               :: max_grad,RMS_grad,s,Grad_Vec_Sign
    real (kind=Rkind)               :: max_disp,RMS_disp,norm_disp
    logical                         :: conv

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_v0'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(Q0)) write(out_unitp,*) '   Q0',Q0
    CALL Write_QML_IRC(IRC_p)
    CALL Write_Model(QModel)
    flush(out_unitp)
  END IF

  IF (IRC_p%Max_it < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' IRC_p is not initialized'
    STOP 'ERROR in QML_IRC_v0: IRC_p is not initialized'
  END IF

  nb_act = size(IRC_p%list_act)

  allocate(QTS_act(nb_act))
  allocate(Qit_act(nb_act))
  allocate(mDQit_act(nb_act))

  allocate(QTS(QModel%ndim))
  allocate(Qit(QModel%ndim))

  allocate(grad(nb_act))
  allocate(hess(nb_act,nb_act))
  allocate(vec(nb_act,nb_act))
  allocate(diag(nb_act))


  IF (present(Q0)) THEN
    QTS(:) = Q0
  ELSE
    CALL get_Q0_Model(QTS,QModel,0)
  END IF

  CALL QML_Opt(QTS,QModel,IRC_p%QML_Opt_t,Q0=QTS)
  !CALL Write_RVec(QTS,    out_unitp,3,name_info='Qit')


  !first check if the geometry in Q0 is a TS
  CALL Eval_Pot(QModel,QTS,PotVal,nderiv=2)
  grad   = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act)
  hess   = PotVal%d2(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act,IRC_p%list_act)

  CALL diagonalization(hess,diag,Vec,nb_act,sort=1)
  write(out_unitp,*) 'grad',grad
  write(out_unitp,*) 'diag',diag
  write(out_unitp,*) 'TS',(count(diag < ZERO) == 1)
  IF ((count(diag < ZERO) /= 1)) STOP 'STOP in QML_IRC: Not a TS'

  !s=0 (TS)
  s = ZERO
  Qit_act(:) = QTS(IRC_p%list_act)
  write(out_unitp,*) 's,Qit,|grad|,E',s,Qit_act,norm2(grad),PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)

  ! forward along the TS vector
  ! new point after the ts
  Qit(:) = QTS
  s = s + IRC_p%Delta_s
  mDQit_act(:) = IRC_p%Delta_s * Vec(:,1)
  Qit_act(:)   = Qit_act + mDQit_act
  CALL Qact_TO_Q(Qit_act,Qit,IRC_p%list_act)

  DO it=1,IRC_p%IRC_Max_it

    !CALL Write_RVec(Qit_act,out_unitp,3,name_info='Qit_act')
    !CALL Write_RVec(Qit,    out_unitp,3,name_info='Qit')
    CALL Eval_Pot(QModel,Qit,PotVal,nderiv=2)
    grad   = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act)
    IF (it == 1) Grad_Vec_Sign = sign(one,dot_product(Vec(:,1),grad))
    !CALL Write_RVec(grad,out_unitp,3,name_info='grad')
    hess   = PotVal%d2(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act,IRC_p%list_act)
    write(out_unitp,*) 's,Qit,|grad|,E',s,Qit_act,norm2(grad),PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)

    s = s + IRC_p%Delta_s
    mDQit_act(:) = Grad_Vec_Sign * IRC_p%Delta_s * grad/norm2(grad)
    Qit_act(:)   = Qit_act + mDQit_act
    !CALL Write_RVec(mDQit_act,out_unitp,3,name_info='mDQit_act')
    CALL Qact_TO_Q(Qit_act,Qit,IRC_p%list_act)

  END DO

  ! backward along the TS vector
  ! new point after the ts
  Qit(:)       = QTS
  Qit_act(:)   = QTS(IRC_p%list_act)
  s            = -IRC_p%Delta_s
  mDQit_act(:) = -IRC_p%Delta_s * Vec(:,1)
  Qit_act(:)   = Qit_act + mDQit_act
  CALL Qact_TO_Q(Qit_act,Qit,IRC_p%list_act)

  DO it=1,IRC_p%IRC_Max_it

    !CALL Write_RVec(Qit_act,out_unitp,3,name_info='Qit_act')
    !CALL Write_RVec(Qit,    out_unitp,3,name_info='Qit')
    CALL Eval_Pot(QModel,Qit,PotVal,nderiv=2)
    grad   = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act)
    IF (it == 1) Grad_Vec_Sign = -sign(one,dot_product(Vec(:,1),grad))
    !CALL Write_RVec(grad,out_unitp,3,name_info='grad')
    hess   = PotVal%d2(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act,IRC_p%list_act)
    write(out_unitp,*) 's,Qit,|grad|,E',s,Qit_act,norm2(grad),PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)

    s = s - IRC_p%Delta_s
    mDQit_act(:) = Grad_Vec_Sign * IRC_p%Delta_s * grad/norm2(grad)
    Qit_act(:)   = Qit_act + mDQit_act
    !CALL Write_RVec(mDQit_act,out_unitp,3,name_info='mDQit_act')
    CALL Qact_TO_Q(Qit_act,Qit,IRC_p%list_act)

  END DO

  STOP


  IF (debug) THEN
    write(out_unitp,*) '   it',it
    write(out_unitp,*) '   Q',Q
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_IRC_v0

  SUBROUTINE QML_IRC_FirstPoint(IRC_at_s,QModel,IRC_p)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

    TYPE(QML_IRC_at_s), intent(inout)            :: IRC_at_s

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p

    TYPE (dnMat_t)                  :: PotVal
    integer                         :: nb_act
    real (kind=Rkind), allocatable  :: hess(:,:)
    real (kind=Rkind), allocatable  :: diag(:),Vec(:,:)

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_FirstPoint'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    flush(out_unitp)
  END IF

  IF (IRC_p%Max_it < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' IRC_p is not initialized'
    STOP 'ERROR in QML_IRC: IRC_p is not initialized'
  END IF

  nb_act = size(IRC_p%list_act)

  allocate(hess(nb_act,nb_act))
  allocate(vec(nb_act,nb_act))
  allocate(diag(nb_act))


  !first check if the geometry in Q0 is a TS
  CALL Eval_Pot(QModel,IRC_at_s%QTS,PotVal,nderiv=2)
  hess   = PotVal%d2(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act,IRC_p%list_act)
  IRC_at_s%Grad_Qact    = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act)

  CALL diagonalization(hess,diag,Vec,nb_act,sort=1)
  write(out_unitp,*) 'grad',IRC_at_s%Grad_Qact
  write(out_unitp,*) 'diag',diag
  write(out_unitp,*) 'TS',(count(diag < ZERO) == 1)
  IF ((count(diag < ZERO) /= 1)) STOP 'STOP in QML_IRC_NewPoint: Not a TS at s=0'

  !s=0 (TS)
  IRC_at_s%s            = ZERO
  IRC_at_s%Ene          = PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)

  IRC_at_s%Qact         = IRC_at_s%QTS(IRC_p%list_act)
  IRC_at_s%DeltaQact    = Vec(:,1)

  IRC_at_s%DeltaQact_TS = Vec(:,1)

  deallocate(hess)
  deallocate(vec)
  deallocate(diag)

  IF (debug) THEN
    write(out_unitp,*) 's,Ene    ',IRC_at_s%s,IRC_at_s%Ene
    write(out_unitp,*) 'Qact     ',IRC_at_s%Qact
    write(out_unitp,*) 'DeltaQact',IRC_at_s%DeltaQact
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_IRC_FirstPoint
SUBROUTINE QML_IRC_NewPoint(IRC_at_s,IRC_at_s_old,QModel,IRC_p,forward)
  USE QMLLib_UtilLib_m
  USE QMLdnSVM_dnMat_m
  USE QMLLib_Matrix_m
  USE QMLLib_diago_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

    TYPE(QML_IRC_at_s), intent(inout)            :: IRC_at_s

    TYPE(QML_IRC_at_s), intent(in)               :: IRC_at_s_old

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p

    logical,            intent(in)               :: forward

    TYPE (dnMat_t)                  :: PotVal
    real (kind=Rkind), allocatable  :: Qit(:)

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_NewPoint'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    flush(out_unitp)
  END IF

  IF (IRC_p%Max_it < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' IRC_p is not initialized'
    STOP 'ERROR in QML_IRC: IRC_p is not initialized'
  END IF

  allocate(Qit(QModel%ndim))


  Qit(:) = IRC_at_s%QTS

  IF (forward) THEN
    IRC_at_s%s      = IRC_at_s_old%s    + IRC_p%Delta_s
  ELSE
    IRC_at_s%s      = IRC_at_s_old%s    - IRC_p%Delta_s
  END IF

  IRC_at_s%Qact   = IRC_at_s_old%Qact + IRC_p%Delta_s * IRC_at_s_old%DeltaQact
  CALL Qact_TO_Q(IRC_at_s%Qact,Qit,IRC_p%list_act)

  !CALL Write_RVec(Qit,    out_unitp,3,name_info='Qit')
  CALL Eval_Pot(QModel,Qit,PotVal,nderiv=1)
  IRC_at_s%Ene         = PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)

  IRC_at_s%Grad_Qact   = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_act)
  IRC_at_s%DeltaQact   = -IRC_at_s%Grad_Qact/norm2(IRC_at_s%Grad_Qact)


  deallocate(Qit)


  IF (debug) THEN
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

END SUBROUTINE QML_IRC_NewPoint

END MODULE IRC_m
