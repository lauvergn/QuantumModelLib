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
MODULE IRC_m
  !$ USE omp_lib
  USE QDUtil_NumParameters_m
  USE Opt_m

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: QML_IRC_t,Init_QML_IRC,QML_IRC,Write_QML_IRC

  TYPE, EXTENDS (QML_Opt_t) :: QML_IRC_t

    integer                           :: IRC_Max_it       = -1 ! it will be set-up after
    real (kind=Rkind)                 :: Delta_s          = ONETENTH**2
    character (len=:), allocatable    :: Direction
    character (len=:), allocatable    :: Method

    integer                           :: order2           = 8  ! used for micro-iteration and with Method='BS'
    character (len=:), allocatable    :: Method2               ! used only when Method='BS'
    integer                           :: m0_BS            = 0

    real (kind=Rkind)                 :: Ene_TS           = ZERO
    real (kind=Rkind), allocatable    :: Grad_QactTS(:)
    real (kind=Rkind), allocatable    :: QTS(:)
    real (kind=Rkind), allocatable    :: QactTS(:)
    real (kind=Rkind), allocatable    :: EigenVec_QactTS(:)
    integer,           allocatable    :: list_actIRC(:)      ! default all coordinates




  END TYPE QML_IRC_t

  TYPE :: QML_IRC_at_s

    real (kind=Rkind)                    :: s
    real (kind=Rkind)                    :: Ene
    real (kind=Rkind), allocatable       :: Qact(:),Grad_Qact(:)

  END TYPE QML_IRC_at_s

  integer                           :: nb_PotEval       = 0


CONTAINS

  SUBROUTINE Init_QML_IRC(IRC_p,QModel,                                         &
                          read_param,param_file_name,nio_param_file)

  USE QMLLib_UtilLib_m
  USE QDUtil_m,         ONLY : TO_lowercase, TO_string
  USE Model_m
  IMPLICIT NONE

    TYPE (QML_IRC_t),   intent(inout)            :: IRC_p
    TYPE (Model_t),     intent(inout)            :: QModel
    logical,            intent(in),    optional  :: read_param
    integer,            intent(in),    optional  :: nio_param_file
    character (len=*),  intent(in),    optional  :: param_file_name



    integer                        :: err_read,nio_loc,i,iact,nb_act
    logical                        :: read_param_loc,MassWeighted
    character (len=:), allocatable :: param_file_name_loc

    integer                        :: Max_it,order2,m0_BS
    real (kind=Rkind)              :: Delta_s
    character (len=Name_longlen)   :: Method,Method2,Direction
    integer,           allocatable :: list_actIRC(:)      ! default the same for opt

    namelist /IRC/ list_actIRC,max_it,Delta_s,Direction,Method,Method2,order2,m0_BS,MassWeighted

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_IRC'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    write(out_unit,*) '   read_param      present?',present(read_param)
    write(out_unit,*) '   nio_param_file  present?',present(nio_param_file)
    write(out_unit,*) '   param_file_name present?',present(param_file_name)
    flush(out_unit)
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
    write(out_unit,*) '   param_file_name ',trim(param_file_name_loc)
    flush(out_unit)
  END IF

  allocate(list_actIRC(QModel%ndim))
  list_actIRC(:) = 0
  
  Max_it          = -1
  Delta_s         = ONETENTH**2
  Direction       = 'both'
  Method          = 'BS'
  Method2         = 'ModMidPoint'
  MassWeighted    = QModel%QM%MassWeighted
  order2          = 8
  m0_BS           = 0

  IF (read_param_loc) THEN
    IF (nio_loc /= in_unit) THEN
      open(unit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
    END IF


    CALL Init_QML_Opt(IRC_p%QML_Opt_t,QModel,read_param=.TRUE.,nio_param_file=nio_loc)

    IF (allocated(param_file_name_loc)) deallocate(param_file_name_loc)


    read(nio_loc,nml=IRC,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "IRC" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Init_QML_IRC'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Some parameter names of the namelist "IRC" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=IRC)
      STOP ' ERROR in Init_QML_IRC'
    END IF

  END IF

  IF (max_it < 0)   Max_it = 10

  method2   = TO_lowercase(method2)
  method    = TO_lowercase(method)
  Direction = TO_lowercase(Direction)
  IF (Direction /= 'both' .AND. Direction /= 'forward' .AND. Direction /= 'reverse') THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' Direction keyword is wrong',trim(Direction)
    write(out_unit,*) ' The possibles values: "both", "forward", "reverse"'
    write(out_unit,*) ' check your data!'
    STOP ' ERROR in Init_QML_IRC'
  END IF

  nb_act = count(list_actIRC > 0)
  IF (nb_act == 0) THEN
    IRC_p%list_actIRC = IRC_p%list_act ! we use the opt ones
  ELSE
    allocate(IRC_p%list_actIRC(nb_act))
    iact = 0
    DO i=1,size(list_actIRC)
      IF (list_actIRC(i) > 0) THEN 
        iact = iact + 1
        IRC_p%list_actIRC(iact) = list_actIRC(i)
      END IF
    END DO
  END IF

  ! the initialisation with QML_Opt_t=IRC_p%QML_Opt_t does not work with ifort 18
  IRC_p%IRC_Max_it    = Max_it
  IRC_p%Delta_s       = Delta_s
  IRC_p%Method        = trim(method)
  IRC_p%Method2       = trim(method2)
  IRC_p%Direction     = trim(Direction)
  IRC_p%order2        = order2
  IRC_p%m0_BS         = m0_BS

  QModel%QM%MassWeighted = MassWeighted
  nb_PotEval    = 0

  CALL Write_QML_IRC(IRC_p)

  IF (debug) THEN
    CALL Write_QML_IRC(IRC_p)
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
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

    write(out_unit,*) ' BEGINNING ',name_sub

    CALL Write_QML_Opt(IRC_p%QML_Opt_t)
    write(out_unit,*) ' list_actIRC     ',IRC_p%list_actIRC
    write(out_unit,*) ' IRC_Maxt_it     ',IRC_p%IRC_Max_it
    write(out_unit,*) ' Delta_s         ',IRC_p%Delta_s
    write(out_unit,*) ' Direction       ',IRC_p%Direction
    write(out_unit,*) ' Method          ',IRC_p%Method
    write(out_unit,*) ' Method2         ',IRC_p%Method2
    write(out_unit,*) ' order2          ',IRC_p%order2
    write(out_unit,*) ' m0_BS           ',IRC_p%m0_BS

    write(out_unit,*) ' nb_PotEval      ',nb_PotEval

    write(out_unit,*) ' END ',name_sub
    flush(out_unit)

  END SUBROUTINE Write_QML_IRC
  SUBROUTINE QML_IRC(Q,QModel,IRC_p,Q0)
    USE Model_m
    USE Opt_m
    IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q(:)
    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(inout)            :: IRC_p

    real (kind=Rkind),  intent(in),    optional  :: Q0(:)


    TYPE (QML_Opt_t)                :: Opt_p
    integer                         :: it
    real (kind=Rkind), allocatable  :: QactOld(:),QactNew(:),grad(:),dQactds_AT_s(:)
    real (kind=Rkind)               :: s,Ene_AT_s
    real (kind=Rkind)               :: forward

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    IF (present(Q0)) write(out_unit,*) '   Q0',Q0
    CALL Write_QML_IRC(IRC_p)
    CALL Write_Model(QModel)
    flush(out_unit)
  END IF

  IF (IRC_p%Max_it < 0) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' IRC_p is not initialized'
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

  write(out_unit,*) ' nb_PotEval      ',nb_PotEval
  !first point + check if the geometry is a TS (s=0)
  CALL QML_IRC_at_TS(IRC_p,QModel)
  write(out_unit,*) ' nb_PotEval      ',nb_PotEval


  s       = ZERO
  QactOld = IRC_p%QactTS
  CALL QML_Recenter_CenterOfMass(s,QactOld,QModel,IRC_p)

  CALL QML_CenterOfMass(s,QactOld,QModel,IRC_p)

  allocate(QactNew(size(QactOld)))
  allocate(grad(size(QactOld)))
  allocate(dQactds_AT_s(size(QactOld)))
  forward = ONE

  IF (IRC_p%Direction == 'both' .OR. IRC_p%Direction == 'forward') THEN
    DO it=0,IRC_p%IRC_Max_it-1

      CALL QML_IRC_ODE(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,                 &
                       QModel,IRC_p,forward=forward,                            &
                       Method=IRC_p%Method,order=IRC_p%order2,grad_AT_s=grad)

      write(out_unit,*) ' nb_PotEval      ',nb_PotEval
      CALL QML_write_IRC_Res(QactOld,s,grad,Ene_AT_s,dQactds_AT_s,QModel,IRC_p)

      s       = s + forward*IRC_p%Delta_s
      QactOld = QactNew

      CALL QML_CenterOfMass(s,QactOld,QModel,IRC_p)

    END DO
    nb_PotEval = nb_PotEval + 1
    CALL QML_IRC_fcn(s,QactOld,dQactds_AT_s,Ene_AT_s,QModel,IRC_p,forward,grad)

    write(out_unit,*) ' nb_PotEval      ',nb_PotEval
    CALL QML_write_IRC_Res(QactOld,s,grad,Ene_AT_s,dQactds_AT_s,QModel,IRC_p)
  END IF

  s       = ZERO
  QactOld = IRC_p%QactTS
  forward = -ONE

  IF (IRC_p%Direction == 'both' .OR. IRC_p%Direction == 'reverse') THEN
    DO it=0,IRC_p%IRC_Max_it-1

      CALL QML_IRC_ODE(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,                 &
                       QModel,IRC_p,forward=forward,                            &
                       Method=IRC_p%Method,order=IRC_p%order2,grad_AT_s=grad)

      write(out_unit,*) ' nb_PotEval      ',nb_PotEval
      CALL QML_write_IRC_Res(QactOld,s,grad,Ene_AT_s,dQactds_AT_s,QModel,IRC_p)

      s       = s + forward*IRC_p%Delta_s
      QactOld = QactNew

      CALL QML_CenterOfMass(s,QactOld,QModel,IRC_p)

    END DO
    nb_PotEval = nb_PotEval + 1
    CALL QML_IRC_fcn(s,QactOld,dQactds_AT_s,Ene_AT_s,QModel,IRC_p,forward,grad)

    write(out_unit,*) ' nb_PotEval      ',nb_PotEval
    CALL QML_write_IRC_Res(QactOld,s,grad,Ene_AT_s,dQactds_AT_s,QModel,IRC_p)
  END IF

  IF (debug) THEN
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE QML_IRC

  SUBROUTINE QML_write_IRC_Res(Qact,s,grad,Ene,dQactds,QModel,IRC_p)
    USE Model_m
    IMPLICIT NONE

    real (kind=Rkind), allocatable, intent(in)  :: Qact(:),grad(:),dQactds(:)
    real (kind=Rkind),              intent(in)  :: s,Ene
    TYPE (Model_t),                 intent(in)  :: QModel
    TYPE (QML_IRC_t),               intent(in)  :: IRC_p



    real (kind=Rkind)               :: u,g_mtu
    real (kind=Rkind), allocatable  :: grad_loc(:),Qact_loc(:)

    IF (QModel%QM%Cart_TO_Q .AND. QModel%QM%MassWeighted) THEN
      Qact_loc = Qact
      grad_loc = grad
      CALL QML_UnMassWeighted_Qact(s,Qact_loc,QModel,IRC_p)
      write(out_unit,*) 's,QxyzMW,|grad|,E',s,Qact,norm2(grad),Ene
      write(out_unit,*) 's,Qxyz,|grad|,E',s,Qact_loc,norm2(grad),Ene
    ELSE
      write(out_unit,*) 's,Qact,|grad|,E',s,Qact,norm2(grad),Ene
    END IF
    write(out_unit,*) 's,grad',s,grad
    write(out_unit,*) 's,g_mt',s,dot_product(dQactds,dQactds)
    !IF (Print_extral) write(out_unit,*) 's,g_mt',s,dot_product(grad,grad)
    flush(out_unit)

  END SUBROUTINE QML_write_IRC_Res

  RECURSIVE SUBROUTINE QML_IRC_ODE(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,     &
                                   QModel,IRC_p,forward,                        &
                                   grad_AT_s,Method,order)
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(in)               :: QactOld(:)
    real (kind=Rkind),  intent(inout)            :: QactNew(:)
    real (kind=Rkind),  intent(inout)            :: dQactds_AT_s(:)
    real (kind=Rkind),  intent(inout), optional  :: grad_AT_s(:)
    real (kind=Rkind),  intent(inout)            :: Ene_AT_s

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p
    real (kind=Rkind),  intent(in)               :: forward
    integer,            intent(in),    optional  :: order
    character (len=*),  intent(in),    optional  :: Method


    character (len=:), allocatable  :: Method_loc
    integer                         :: order_loc


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_ODE'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    write(out_unit,*) '   s      ',s
    write(out_unit,*) '   QactOld',QactOld
    write(out_unit,*) '   forward',forward
    flush(out_unit)
  END IF

  IF (present(order)) THEN
    order_loc = order
  ELSE
    order_loc = IRC_p%order2
  END IF

  IF (present(Method)) THEN
    Method_loc = Method
  ELSE
    Method_loc = IRC_p%Method
  END IF

  IF (present(grad_AT_s)) THEN
    SELECT CASE (Method_loc)
    CASE('midpoint','modmidpoint')
      CALL QML_IRC_ModMidPoint(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,         &
                               QModel,IRC_p,order_loc,forward,grad_AT_s)
    CASE('euler')
      CALL QML_IRC_mEuler(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,              &
                          QModel,IRC_p,order_loc,forward,grad_AT_s)
    CASE('bs','bulirsch-stoer')
      CALL QML_IRC_BS(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,                  &
                      QModel,IRC_p,order_loc,forward,grad_AT_s)
    CASE Default
      CALL QML_IRC_mEuler(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,              &
                          QModel,IRC_p,order_loc,forward,grad_AT_s)
    END SELECT
  ELSE
    SELECT CASE (Method_loc)
    CASE('midpoint','modmidpoint')
      CALL QML_IRC_ModMidPoint(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,         &
                               QModel,IRC_p,order_loc,forward)
    CASE('euler')
      CALL QML_IRC_mEuler(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,              &
                          QModel,IRC_p,order_loc,forward)
    CASE('bs','bulirsch-stoer')
      STOP 'ERROR in QML_IRC_ODE: bulirsch-stoer method sould be called with grad_AT_s(:)'
    CASE Default
      CALL QML_IRC_mEuler(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,              &
                          QModel,IRC_p,order_loc,forward)
    END SELECT
  END IF

  IF (debug) THEN
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

END SUBROUTINE QML_IRC_ODE

  SUBROUTINE QML_IRC_BS(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,                &
                        QModel,IRC_p,order,forward,grad_AT_s)

    USE Model_m
    IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(in)               :: QactOld(:)
    real (kind=Rkind),  intent(inout)            :: QactNew(:)
    real (kind=Rkind),  intent(inout)            :: grad_AT_s(:)
    real (kind=Rkind),  intent(inout)            :: Ene_AT_s
    real (kind=Rkind),  intent(inout)            :: dQactds_AT_s(:)

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p
    real (kind=Rkind),  intent(in)               :: forward
    integer,            intent(in)               :: order


    integer                         :: nb_act
    real (kind=Rkind)               :: Ene_loc
    real (kind=Rkind), allocatable  :: dQactds_loc(:)

    real (kind=Rkind), allocatable :: yt0(:,:),yt1(:,:)
    real (kind=Rkind), allocatable :: yerr(:)
    real(kind=Rkind)               :: x,err0,err1
    integer                        :: i,j,m

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_BS'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    write(out_unit,*) '   s      ',s
    write(out_unit,*) '   QactOld',QactOld
    write(out_unit,*) '   forward',forward
    write(out_unit,*) '   order  ',order
    write(out_unit,*) '   m0_BS  ',IRC_p%m0_BS
    flush(out_unit)
  END IF

  IF (order < 1) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' order < 1. order:',order
    write(out_unit,*) ' order MUST be larger than 0'
    STOP 'ERROR in QML_IRC_BS: order < 1.'
  END IF
  nb_act = size(QactOld)
  allocate(yt0(nb_act,0:0))
  allocate(yerr(nb_act))
  allocate(dQactds_loc(nb_act))

  yerr     = ZERO
  m        = QML_BS_m(0,IRC_p%m0_BS)
  !yt0(:,0) = QactOld

  CALL QML_IRC_ODE(s,QactOld,yt0(:,0),Ene_AT_s,dQactds_AT_s,QModel,IRC_p,forward,&
                   Method=IRC_p%Method2,order=m,grad_AT_s=grad_AT_s)

  err0 = huge(ONE)

  DO j=1,order
    allocate(yt1(nb_act,0:j))
    m = QML_BS_m(j,IRC_p%m0_BS)

    CALL QML_IRC_ODE(s,QactOld,yt1(:,0),Ene_loc,dQactds_loc,QModel,IRC_p,forward,&
                     Method=IRC_p%Method2,order=m)

    !extrapolation
    DO i=1,j
      !x = real(m0+2*j+2,kind=Rkind)/real(2*(j-i)+2,kind=Rkind)
      x = real(m,kind=Rkind)/real(m-QML_BS_m(i-1,IRC_p%m0_BS),kind=Rkind)
      x = ONE/(x**2-ONE)
      yerr(:)  = yt1(:,i-1) - yt0(:,i-1)
      yt1(:,i) = yt1(:,i-1) + yerr(:)*x
    END DO

    yerr(:) = yt1(:,j) - yt0(:,j-1)
    err1    = norm2(yerr)

    deallocate(yt0)
    allocate(yt0(nb_act,0:j))
    DO i=lbound(yt1,dim=2),ubound(yt1,dim=2)
      yt0(:,i) = yt1(:,i)
    END DO
    deallocate(yt1)

    IF (err1 < 1.d-10) EXIT

    err0 = err1
  END DO
  write(out_unit,*) 'end QML_IRC_BS',min(j,order),err1

  QactNew(:) = yt0(:,min(j,order))


  IF (debug) THEN
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

END SUBROUTINE QML_IRC_BS

SUBROUTINE QML_IRC_mEuler(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,            &
                            QModel,IRC_p,m,forward,grad)
    USE Model_m
    IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(in)               :: QactOld(:)
    real (kind=Rkind),  intent(inout)            :: QactNew(:)
    real (kind=Rkind),  intent(inout)            :: dQactds_AT_s(:)
    real (kind=Rkind),  intent(inout), optional  :: grad(:)
    real (kind=Rkind),  intent(inout)            :: Ene_AT_s

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p
    real (kind=Rkind),  intent(in)               :: forward
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
    write(out_unit,*) ' BEGINNING ',name_sub
    write(out_unit,*) '   s      ',s
    write(out_unit,*) '   QactOld',QactOld
    flush(out_unit)
  END IF

  IF (m < 1) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' m < 1',m
    write(out_unit,*) ' m MUST be larger than 0'
    STOP 'ERROR in QML_IRC_mEuler: m < 1'
  END IF

  s_loc = s
  Qact  = QactOld
  allocate(dQact(size(QactOld)))
  Delta_s = IRC_p%Delta_s/m

  DO it=1,m

    IF (it == 1) THEN
      IF (present(grad)) THEN
        CALL QML_IRC_fcn(s_loc,Qact,dQact,Ene_loc,QModel,IRC_p,forward,grad)
      ELSE
        CALL QML_IRC_fcn(s_loc,Qact,dQact,Ene_loc,QModel,IRC_p,forward)
      END IF
      Ene_AT_s     = Ene_loc
      dQactds_AT_s = dQact
    ELSE
      CALL QML_IRC_fcn(s_loc,Qact,dQact,Ene_loc,QModel,IRC_p,forward)
    END IF
    nb_PotEval = nb_PotEval + 1
    IF (debug) write(out_unit,*) 's,Qact,E',s_loc,Qact,Ene_loc
    IF (debug) flush(out_unit)

    s_loc = s_loc + forward * Delta_s
    Qact  = Qact  +   dQact * Delta_s
  END DO

  QactNew = Qact


  IF (debug) THEN
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE QML_IRC_mEuler
  SUBROUTINE QML_IRC_ModMidPoint(s,QactOld,QactNew,Ene_AT_s,dQactds_AT_s,       &
                                 QModel,IRC_p,m,forward,grad)
    USE Model_m
    IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(in)               :: QactOld(:)
    real (kind=Rkind),  intent(inout)            :: QactNew(:)
    real (kind=Rkind),  intent(inout), optional  :: grad(:)
    real (kind=Rkind),  intent(inout)            :: dQactds_AT_s(:)
    real (kind=Rkind),  intent(inout)            :: Ene_AT_s

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p
    real (kind=Rkind),  intent(in)               :: forward
    integer,            intent(in)               :: m


    integer                         :: it
    real (kind=Rkind), allocatable  :: Qact(:),dQact(:)
    real (kind=Rkind)               :: s_loc,Ene_loc,Delta_s
    real (kind=Rkind), allocatable  :: zkm(:),zk(:),zkp(:)

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_ModMidPoint'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    write(out_unit,*) '   s      ',s
    write(out_unit,*) '   QactOld',QactOld
    flush(out_unit)
  END IF

  IF (m < 1) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' m < 1',m
    write(out_unit,*) ' m MUST be larger than 0'
    STOP 'ERROR in QML_IRC_ModMidPoint: m < 1'
  END IF

  s_loc = s
  Qact  = QactOld
  allocate(dQact(size(QactOld)))
  Delta_s = IRC_p%Delta_s/m

  zkm = QactOld
  IF (present(grad)) THEN
    CALL QML_IRC_fcn(s_loc,Qact,dQact,Ene_AT_s,QModel,IRC_p,forward,grad)
  ELSE
    CALL QML_IRC_fcn(s_loc,Qact,dQact,Ene_AT_s,QModel,IRC_p,forward)
  END IF
  dQactds_AT_s = dQact
  nb_PotEval   = nb_PotEval + 1
  zk  = QactOld + dQact * Delta_s


  DO it=1,m-1
    CALL QML_IRC_fcn(s_loc,zk,dQact,Ene_loc,QModel,IRC_p,forward)
    nb_PotEval = nb_PotEval + 1

    zkp = zkm + TWO*Delta_s * dQact
    zkm = zk
    zk  = zkp
  END DO

  CALL QML_IRC_fcn(s_loc,zk,zkp,Ene_loc,QModel,IRC_p,forward)
  nb_PotEval = nb_PotEval + 1

  dQact = Delta_s * zkp

  zkp = zk + zkm

  zk =  zkp + dQact
  QactNew = HALF * zk


  IF (debug) THEN
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE QML_IRC_ModMidPoint

  SUBROUTINE QML_IRC_fcn(s,Qact,dQact,Ene_AT_s,QModel,IRC_p,forward,Grad)
  USE ADdnSVM_m, ONLY : dnMat_t
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(in)               :: Qact(:)
    real (kind=Rkind),  intent(inout)            :: dQact(:)
    real (kind=Rkind),  intent(inout)            :: Ene_AT_s

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p
    real (kind=Rkind),  intent(in)               :: forward
    real (kind=Rkind),  intent(inout), optional  :: Grad(:)

    TYPE (dnMat_t)                  :: PotVal
    real (kind=Rkind), allocatable  :: Qit(:)
    integer                         :: i,iat

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_IRC_fcn'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    write(out_unit,*) 's,forward  ',s,forward

    flush(out_unit)
  END IF

  IF (s == ZERO) THEN
    dQact    =  forward * IRC_p%EigenVec_QactTS
    Ene_AT_s = IRC_p%Ene_TS
    IF (present(grad)) grad = IRC_p%Grad_QactTS
  ELSE
    Qit = IRC_p%QTS

    CALL Qact_TO_Q(Qact,Qit,IRC_p%list_actIRC)

    CALL Eval_Pot(QModel,Qit,PotVal,nderiv=1)

    Ene_AT_s = PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)

    dQact    = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_actIRC)

    CALL QML_Recenter_CenterOfMass(s,dQact,QModel,IRC_p)

    !write(out_unit,*) 'grad at s',s,norm2(dQact)
    !write(out_unit,*) 'grad/s at s',s,norm2(dQact)/abs(s)
    IF (norm2(dQact) < IRC_p%Thresh_RMS_grad*TEN) write(out_unit,*) 'WARNING small grad at s',s
    IF (present(grad)) grad = dQact



    dQact    = -dQact/norm2(dQact)

    deallocate(Qit)
  END IF

  IF (debug) THEN
    write(out_unit,*) 's,Ene    ',s,Ene_AT_s
    write(out_unit,*) 'Qact     ',Qact
    write(out_unit,*) 'dQact    ',dQact
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

END SUBROUTINE QML_IRC_fcn
SUBROUTINE QML_MassWeighted_Qact(s,Qact,QModel,IRC_p)
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(inout)            :: Qact(:)

    TYPE (Model_t),     intent(in)               :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p

    integer           :: i,iat

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_MassWeighted_Qact'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    flush(out_unit)
  END IF

  IF (QModel%QM%Cart_TO_Q .AND. QModel%QM%MassWeighted) THEN
    DO i=1,size(Qact)
      iat = (IRC_p%list_actIRC(i)+2)/3
      Qact(i) = Qact(i) * sqrt(QModel%QM%masses(iat))
    END DO
  END IF

  IF (debug) THEN
    write(out_unit,*) 'Qact     ',Qact
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

END SUBROUTINE QML_MassWeighted_Qact
SUBROUTINE QML_UnMassWeighted_Qact(s,Qact,QModel,IRC_p)
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(inout)            :: Qact(:)

    TYPE (Model_t),     intent(in)               :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p

    integer           :: i,iat

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_UnMassWeighted_Qact'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    flush(out_unit)
  END IF

  IF (QModel%QM%Cart_TO_Q .AND. QModel%QM%MassWeighted) THEN
    DO i=1,size(Qact)
      iat = (IRC_p%list_actIRC(i)+2)/3
      Qact(i) = Qact(i) / sqrt(QModel%QM%masses(iat))
    END DO
  END IF

  IF (debug) THEN
    write(out_unit,*) 'Qact     ',Qact
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

END SUBROUTINE QML_UnMassWeighted_Qact
SUBROUTINE QML_CenterOfMass(s,Qact,QModel,IRC_p)
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(in)               :: Qact(:)

    TYPE (Model_t),     intent(in)               :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p

    real (kind=Rkind)                :: QCOM(3)
    real (kind=Rkind), allocatable   :: Qact_loc(:)
    integer                          :: i,iat
    real (kind=Rkind)                :: Mtot

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_CenterOfMass'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    flush(out_unit)
  END IF

  IF (QModel%QM%Cart_TO_Q) THEN

    Mtot = sum(QModel%QM%masses)
    QCOM(:) = ZERO

    IF (QModel%QM%MassWeighted) THEN
      Qact_loc = Qact
      CALL QML_UnMassWeighted_Qact(s,Qact_loc,QModel,IRC_p)
      DO i=1,size(Qact)
        iat = (IRC_p%list_actIRC(i)+2)/3
        SELECT CASE (mod(IRC_p%list_actIRC(i),3))
        CASE (1) ! x
          QCOM(1) = QCOM(1) + Qact_loc(i)*QModel%QM%masses(iat)/Mtot
        CASE (2) ! y
          QCOM(2) = QCOM(2) + Qact_loc(i)*QModel%QM%masses(iat)/Mtot
        CASE (0) ! x
          QCOM(3) = QCOM(3) + Qact_loc(i)*QModel%QM%masses(iat)/Mtot
        END SELECT
      END DO
      CALL QML_MassWeighted_Qact(s,Qact_loc,QModel,IRC_p)
    ELSE
      DO i=1,size(Qact)
        iat = (IRC_p%list_actIRC(i)+2)/3
        SELECT CASE (mod(IRC_p%list_actIRC(i),3))
        CASE (1) ! x
          QCOM(1) = QCOM(1) + Qact(i)*QModel%QM%masses(iat)/Mtot
        CASE (2) ! y
          QCOM(2) = QCOM(2) + Qact(i)*QModel%QM%masses(iat)/Mtot
        CASE (0) ! x
          QCOM(3) = QCOM(3) + Qact(i)*QModel%QM%masses(iat)/Mtot
        END SELECT
      END DO
    END IF
  END IF

  write(out_unit,*) s,'QCOM    ',QCOM

  IF (debug) THEN
    write(out_unit,*) 'Qact     ',Qact
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

END SUBROUTINE QML_CenterOfMass

SUBROUTINE QML_Recenter_CenterOfMass(s,Qact,QModel,IRC_p)
  USE Model_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)               :: s
    real (kind=Rkind),  intent(inout)            :: Qact(:)

    TYPE (Model_t),     intent(in)               :: QModel
    TYPE (QML_IRC_t),   intent(in)               :: IRC_p

    real (kind=Rkind) :: QCOM(3)
    integer           :: i,iat
    real (kind=Rkind) :: Mtot

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_Recenter_CenterOfMass'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unit,*) ' BEGINNING ',name_sub
    flush(out_unit)
  END IF

  IF (QModel%QM%Cart_TO_Q) THEN

    Mtot = sum(QModel%QM%masses)
    QCOM(:) = ZERO

    IF (QModel%QM%MassWeighted) CALL QML_UnMassWeighted_Qact(s,Qact,QModel,IRC_p)

    DO i=1,size(Qact)
      iat = (IRC_p%list_actIRC(i)+2)/3
      SELECT CASE (mod(IRC_p%list_actIRC(i),3))
      CASE (1) ! x
        QCOM(1) = QCOM(1) + Qact(i)*QModel%QM%masses(iat)/Mtot
      CASE (2) ! y
        QCOM(2) = QCOM(2) + Qact(i)*QModel%QM%masses(iat)/Mtot
      CASE (0) ! x
        QCOM(3) = QCOM(3) + Qact(i)*QModel%QM%masses(iat)/Mtot
      END SELECT
    END DO

    DO i=1,size(Qact)
      iat = (IRC_p%list_actIRC(i)+2)/3
      SELECT CASE (mod(IRC_p%list_actIRC(i),3))
      CASE (1) ! x
        Qact(i) = Qact(i) - QCOM(1)
      CASE (2) ! y
        Qact(i) = Qact(i) - QCOM(2)
      CASE (0) ! x
        Qact(i) = Qact(i) - QCOM(3)
      END SELECT
    END DO

    IF (QModel%QM%MassWeighted) CALL QML_MassWeighted_Qact(s,Qact,QModel,IRC_p)

  END IF


  IF (debug) THEN
    write(out_unit,*) 'Qact     ',Qact
    write(out_unit,*) s,'QCOM    ',QCOM
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

END SUBROUTINE QML_Recenter_CenterOfMass

  SUBROUTINE QML_IRC_at_TS(IRC_p,QModel)
    USE QDUtil_m,         ONLY : diagonalization
    USE ADdnSVM_m
    USE Model_m
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
      write(out_unit,*) ' BEGINNING ',name_sub
      flush(out_unit)
    END IF
  
    IF (IRC_p%Max_it < 0) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' IRC_p is not initialized'
      STOP 'ERROR in QML_IRC_at_TS: IRC_p is not initialized'
    END IF
  
    nb_act = size(IRC_p%list_actIRC)
  
    allocate(hess(nb_act,nb_act))
    allocate(vec(nb_act,nb_act))
    allocate(diag(nb_act))
  
  
    !first check if the geometry in Q0 is a TS
    CALL Eval_Pot(QModel,IRC_p%QTS,PotVal,nderiv=2)
    hess              = PotVal%d2(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_actIRC,IRC_p%list_actIRC)
    IRC_p%Grad_QactTS = PotVal%d1(IRC_p%i_surf,IRC_p%i_surf,IRC_p%list_actIRC)
  
    CALL diagonalization(hess,diag,Vec,nb_act,sort=1)
    write(out_unit,*) 'grad',IRC_p%Grad_QactTS
    write(out_unit,*) 'diag',diag
    write(out_unit,*) 'TS?',(count(diag < ZERO) == 1)
    IF ((count(diag < ZERO) /= 1)) STOP 'STOP in QML_IRC_at_TS: Not a TS at s=0'
  
    !s=0 (TS)
    IRC_p%Ene_TS          = PotVal%d0(IRC_p%i_surf,IRC_p%i_surf)
    IRC_p%QactTS          = IRC_p%QTS(IRC_p%list_actIRC)
    IRC_p%EigenVec_QactTS = Vec(:,1)
  
    deallocate(hess)
    deallocate(vec)
    deallocate(diag)

    IF (debug) THEN
      write(out_unit,*) 's,Ene          ',ZERO,IRC_p%Ene_TS
      write(out_unit,*) 'Qact           ',IRC_p%QactTS
      write(out_unit,*) 'EigenVec_QactTS',IRC_p%EigenVec_QactTS
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE QML_IRC_at_TS
  FUNCTION QML_BS_m(j,m0) RESULT(m)
    IMPLICIT NONE

    integer :: m
    integer, intent(in) :: j,m0
    m = m0+2*j+2
  END FUNCTION QML_BS_m
END MODULE IRC_m
