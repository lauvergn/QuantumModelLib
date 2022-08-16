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

MODULE Model_m
!$ USE omp_lib
  USE QMLLib_NumParameters_m

  !USE ADdnSVM_m

  USE QML_Empty_m

  USE QML_Template_m
  USE QML_Test_m
  USE QML_Morse_m
  USE QML_Poly1D_m
  USE QML_H2_m

  USE QML_HenonHeiles_m
  USE QML_Tully_m

  USE QML_PSB3_m
  USE QML_Retinal_JPCB2000_m

  USE QML_HONO_m
  USE QML_HNNHp_m
  USE QML_H2SiN_m
  USE QML_H2NSi_m
  USE QML_ClH2p_m
  USE QML_ClH2p_Botschwina_m
  USE QML_HNO3_m
  USE QML_NO3_m
  USE QML_CH5_m
  USE QML_PH4_m

  USE QML_HOO_DMBE_m
  USE QML_H3_m
  USE QML_HCN_Murrell_m

  USE QML_OneDSOC_1S1T_m
  USE QML_OneDSOC_2S1T_m

  USE QML_LinearHBond_m
  USE QML_TwoD_MullerBrown_m
  USE QML_Buck_m
  USE QML_Phenol_m
  USE QML_Sigmoid_m
  USE QML_TwoD_m
  USE QML_TwoD_RJDI2014_m

  USE AdiaChannels_Basis_m

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Model_t,Init_Model,Eval_Pot,Eval_Func,Eval_tab_HMatVibAdia,Eval_dnHVib_ana
  PUBLIC :: check_alloc_QM,check_Init_QModel
  PUBLIC :: Write0_Model,Write_Model,Write_QdnV_FOR_Model
  PUBLIC :: calc_pot,calc_grad,calc_hess,calc_pot_grad,calc_pot_grad_hess
  PUBLIC :: Check_analytical_numerical_derivatives
  PUBLIC :: Eval_pot_ON_Grid,get_Q0_Model,Qact_TO_Q
  PUBLIC :: set_step_epsi_Model

  TYPE :: Model_t
    ! add nsurf and ndim to avoid crash when using the driver without initialization
    ! at the intialization, the variables are set-up to the correct values and are
    !   identical to QM%nsurf and QM%ndim ones respectively.
    integer                           :: nsurf       = 0
    integer                           :: ndim        = 0
    CLASS (QML_Empty_t),  allocatable :: QM
    TYPE (QML_Basis_t),   allocatable :: Basis ! Basis for the adiabatic separation between coordinates
    logical                           :: opt         = .FALSE.
    logical                           :: irc         = .FALSE.
  END TYPE Model_t

  !real (kind=Rkind)                     :: step = ONETENTH**4 ! model TWOD => 0.4e-7 (nderiv=2)
  real (kind=Rkind)                     :: step = ONETENTH**3 ! model TWOD => 0.6e-9 (nderiv=2)
  !real (kind=Rkind)                     :: step = ONETENTH**2 ! model TWOD => 0.1e-7 (nderiv=2)

  real (kind=Rkind)                     :: epsi = ONETENTH**10

#if defined(__QML_VER)
      character (len=Name_len) :: QML_version = __QML_VER
#else
      character (len=Name_len) :: QML_version = "unknown: -D__QML_VER=?"
#endif

#if defined(__QMLPATH)
      character (len=Line_len) :: QML_path   =                         &
       __QMLPATH
#else
      character (len=Line_len) :: QML_path   = '~/QuantumModelLib'
#endif

#if defined(__COMPILE_DATE)
      character (len=Line_len) :: compile_date = __COMPILE_DATE
#else
      character (len=Line_len) :: compile_date = "unknown: -D__COMPILE_DATE=?"
#endif

#if defined(__COMPILE_HOST)
      character (len=Line_len) :: compile_host = __COMPILE_HOST
#else
      character (len=Line_len) :: compile_host = "unknown: -D__COMPILE_HOST=?"
#endif
#if defined(__COMPILER)
      character (len=Line_len) :: compiler = __COMPILER
#else
      character (len=Line_len) :: compiler = "unknown: -D__COMPILER=?"
#endif
#if defined(__COMPILER_VER)
      character (len=Line_len) :: compiler_ver = __COMPILER_VER
#else
      character (len=Line_len) :: compiler_ver = "unknown: -D__COMPILER_VER=?"
#endif
#if defined(__COMPILER_OPT)
      character (len=Line_len) :: compiler_opt = &
      __COMPILER_OPT
#else
      character (len=Line_len) :: compiler_opt = "unknown: -D__COMPILER_OPT=?"
#endif
#if defined(__COMPILER_LIBS)
      character (len=Line_len) :: compiler_libs = __COMPILER_LIBS
#else
      character (len=Line_len) :: compiler_libs = "unknown: -D__COMPILER_LIBS=?"
#endif



  TYPE(Model_t), PUBLIC  :: QuantumModel

CONTAINS


  SUBROUTINE Read_Model(QModel_inout,nio,read_nml1,opt1,IRC1)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    TYPE (QML_Empty_t), intent(inout)  :: QModel_inout ! variable to transfer info to the init
    integer,             intent(in)    :: nio
    logical,             intent(inout) :: read_nml1,opt1,IRC1

    ! local variable
    integer, parameter :: max_act = 10
    integer :: ndim,nsurf,nderiv,option,printlevel,nb_Channels
    logical :: adiabatic,numeric,PubliUnit,read_nml
    logical :: Vib_adia,print_EigenVec_Grid,print_EigenVec_Basis
    logical :: opt,IRC
    logical :: Phase_checking,Phase_Following
    logical :: Cart_TO_Q

    character (len=20) :: pot_name
    integer :: err_read,nb_act
    integer :: list_act(max_act)

    ! Namelists for input file
    namelist /potential/ ndim,nsurf,pot_name,numeric,adiabatic,option,PubliUnit,&
                         Phase_Checking,Phase_Following,Cart_TO_Q,              &
                         read_nml,printlevel,Vib_adia,nb_Channels,list_act,     &
                         print_EigenVec_Grid,print_EigenVec_Basis,opt,IRC

!    ! Default values defined
    printlevel      = 0
    ndim            = QModel_inout%ndim
    nsurf           = QModel_inout%nsurf
    adiabatic       = QModel_inout%adiabatic

    Phase_Checking  = QModel_inout%Phase_Checking
    Phase_Following = QModel_inout%Phase_Following
    Cart_TO_Q       = QModel_inout%Cart_TO_Q

    Vib_adia        = QModel_inout%Vib_adia
    nb_Channels     = 0
    list_act(:)     = 0

    print_EigenVec_Grid  = .FALSE.
    print_EigenVec_Basis = .FALSE.

    pot_name    = 'morse'
    numeric     = .FALSE.
    PubliUnit   = .FALSE.
    read_nml    = read_nml1 ! if T, read the namelist in PotLib (HenonHeiles ....)
    option      = -1 ! no option

    !extra actions
    opt         = .FALSE.
    IRC         = .FALSE.

    write(out_unitp,*) 'Reading input file . . .'
    read(nio,nml=potential,IOSTAT=err_read)

    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_Model'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "potential" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_Model: End-of-file or End-of-record while reading the namelist'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_Model'
      write(out_unitp,*) ' Some parameter names of the namelist "potential" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=potential)
      STOP ' ERROR in Read_Model: wrong parameter(s) in the namelist'
    END IF

    !write(out_unitp,nml=potential)

    read_nml1                         = read_nml
    opt1                              = opt
    IRC1                              = IRC

    QModel_inout%option               = option
    print_level                       = printlevel ! from them module QMLLib_NumParameters_m.f90
    QModel_inout%ndim                 = ndim
    QModel_inout%nsurf                = nsurf
    QModel_inout%adiabatic            = adiabatic
    QModel_inout%numeric              = numeric
    QModel_inout%Phase_Checking       = (adiabatic .AND. Phase_Checking)
    QModel_inout%Phase_Following      = (adiabatic .AND. Phase_Following)
    QModel_inout%Cart_TO_Q            = Cart_TO_Q

    QModel_inout%pot_name             = trim(pot_name)
    QModel_inout%PubliUnit            = PubliUnit

    QModel_inout%print_EigenVec_Grid  = print_EigenVec_Grid
    QModel_inout%print_EigenVec_Basis = print_EigenVec_Basis


    IF (Vib_adia) THEN
      QModel_inout%nb_Channels    = nb_Channels
      QModel_inout%Vib_adia       = Vib_adia

      IF (nb_Channels == 0) THEN
        write(out_unitp,*) ' ERROR in Read_Model'
        write(out_unitp,*) ' Vib_adia=t and nb_Channels = 0'
        write(out_unitp,*) ' You have to define "nb_Channels" in the namelist.'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*)
        STOP ' ERROR in Read_Model: define "nb_Channels" in the namelist'
      END IF

      nb_act = count(list_act /= 0)
      IF (nb_act == 0) THEN
        write(out_unitp,*) ' ERROR in Read_Model'
        write(out_unitp,*) ' Vib_adia=t and nb_act = 0'
        write(out_unitp,*) ' You have to define "list_act(:)" in the namelist.'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*)
        STOP ' ERROR in Read_Model: define "list_act(:)" in the namelist.'
      END IF
      QModel_inout%list_act = list_act(1:nb_act)

      IF (count(QModel_inout%list_act == 0) /= 0) THEN
        write(out_unitp,*) ' ERROR in Read_Model'
        write(out_unitp,*) ' list_act(:) is wrong.'
        write(out_unitp,*) ' list_act(1:nb_act) has some 0 :',list_act(1:nb_act)
        write(out_unitp,*) ' You have to define in the namelist list_act(:)'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*)
        STOP ' ERROR in Read_Model: "list_act(:)" with 0 in the namelist.'
      END IF
    END IF


  END SUBROUTINE Read_Model

  SUBROUTINE Init_Model(QModel,pot_name,ndim,nsurf,adiabatic,Cart_TO_Q,         &
                        read_param,param_file_name,nio_param_file,              &
                        option,PubliUnit,Print_init,Vib_adia,                   &
                        Phase_Following,Phase_checking)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    TYPE (Model_t),      intent(inout)            :: QModel

    character (len=*),   intent(in),    optional :: pot_name
    integer,             intent(in),    optional :: ndim,nsurf
    logical,             intent(in),    optional :: adiabatic
    logical,             intent(in),    optional :: Cart_TO_Q

    logical,             intent(in),    optional :: read_param
    integer,             intent(in),    optional :: nio_param_file
    character (len=*),   intent(in),    optional :: param_file_name

    integer,             intent(in),    optional :: option
    logical,             intent(in),    optional :: PubliUnit
    logical,             intent(in),    optional :: Print_init

    logical,             intent(in),    optional :: Vib_adia
    logical,             intent(in),    optional :: Phase_Following
    logical,             intent(in),    optional :: Phase_checking


    ! local variables
    TYPE(QML_Empty_t)             :: QModel_in ! variable to transfer info to the init
    integer                        :: i,nio_loc,i_inact,nb_inact
    logical                        :: read_param_loc,read_nml,Print_init_loc
    logical,           allocatable :: list_Q(:)
    character (len=:), allocatable :: param_file_name_loc,pot_name_loc
    real (kind=Rkind), allocatable :: Q0(:)

    Print_init_loc = .TRUE.
    IF (present(Print_init)) Print_init_loc = Print_init

    IF (Print_init_loc) THEN
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== QML: Quantum Model Lib (E-CAM) ==============='
      write(out_unitp,*) '== QML version:       ',QML_version
      write(out_unitp,*) '== QML path:          ',trim(adjustl(QML_path))
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) '== Compiled on       "',trim(compile_host), '" the ',trim(compile_date)
      write(out_unitp,*) '== Compiler:         ',trim(compiler)
      write(out_unitp,*) '== Compiler version: ',trim(compiler_ver)
      write(out_unitp,*) '== Compiler options: ',trim(compiler_opt)
      write(out_unitp,*) '== Compiler libs:     ',trim(compiler_libs)
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) 'QML is under GNU LGPL3 license and '
      write(out_unitp,*) '  is written by David Lauvergnat [1]'
      write(out_unitp,*) '  with contributions of'
      write(out_unitp,*) '     Félix MOUHAT [2]'
      write(out_unitp,*) '     Liang LIANG [3]'
      write(out_unitp,*) '     Emanuele MARSILI [1,4]'
      write(out_unitp,*)
      write(out_unitp,*) '[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
      write(out_unitp,*) '[2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France'
      write(out_unitp,*) '[3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France'
      write(out_unitp,*) '[4]: Durham University, Durham, UK'
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== Initialization of the Model =================='
    END IF

    ! set the "File_path" in the Lib_module.f90
    File_path = trim(adjustl(QML_path))

    ! test the QML_path variable (it enables to test is the QML directory has been moved)
    CALL check_QML_Path()


    IF (allocated(QModel%QM)) deallocate(QModel%QM)

    IF (present(ndim)) THEN
      QModel_in%ndim      = ndim
    ELSE
      QModel_in%ndim      = 0
    END IF
    IF (present(nsurf)) THEN
      QModel_in%nsurf     = nsurf
    ELSE
      QModel_in%nsurf     = 0
    END IF

    IF (present(adiabatic)) THEN
      QModel_in%adiabatic = adiabatic
    ELSE
      QModel_in%adiabatic = .TRUE.
    END IF

    IF (present(Cart_TO_Q)) THEN
      QModel_in%Cart_TO_Q = Cart_TO_Q
    ELSE
      QModel_in%Cart_TO_Q = .FALSE.
    END IF

    IF (present(option)) THEN
      QModel_in%option = option
    ELSE
      QModel_in%option = -1
    END IF
    !IF (Print_init_loc) write(out_unitp,*) 'option: ',QModel_in%option

    IF (present(Vib_adia)) THEN
      QModel_in%Vib_adia = Vib_adia
    ELSE
      QModel_in%Vib_adia = .FALSE.
    END IF

    IF (present(Phase_Following)) THEN
      QModel_in%Phase_Following = Phase_Following
    ELSE
      QModel_in%Phase_Following = .TRUE.
    END IF

    IF (present(Phase_checking)) THEN
      QModel_in%Phase_checking = Phase_checking
    ELSE
      QModel_in%Phase_checking = .TRUE.
    END IF

    IF (present(PubliUnit)) THEN
      QModel_in%PubliUnit      = PubliUnit
    ELSE
      QModel_in%PubliUnit      = .FALSE.
    END IF


    IF (present(read_param)) THEN
      read_param_loc = read_param
    ELSE
      read_param_loc = .FALSE.
    END IF
    IF (QModel_in%Vib_adia) read_param_loc = .TRUE.


    IF (present(nio_param_file)) THEN
      nio_loc = nio_param_file
    ELSE
      IF (present(param_file_name)) THEN
        IF (len_trim(param_file_name) == 0) THEN
          param_file_name_loc = strdup("input.dat")
          nio_loc = 99 ! this value is not used
        ELSE
          param_file_name_loc = strdup(param_file_name)
          nio_loc = 99 ! this value is not used
        END IF
      ELSE
      nio_loc = in_unitp
      END IF
    END IF

    IF (present(pot_name)) THEN
      pot_name_loc  = strdup(pot_name)
      CALL string_uppercase_TO_lowercase(pot_name_loc)
      read_param_loc = (read_param_loc .OR.  pot_name_loc == 'read_model')
    ELSE
      IF (.NOT. read_param_loc) THEN
        write(out_unitp,*) 'ERROR in Init_Model'
        write(out_unitp,*) ' pot_name is not present and read_param=F'
        STOP 'ERROR in Init_Model: pot_name is not present and read_param=F'
      END IF
    END IF

    read_nml       = .FALSE.
    IF (read_param_loc) THEN
      IF (nio_loc /= in_unitp) THEN
        open(newunit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
      END IF
      CALL Read_Model(QModel_in,nio_loc,read_nml,QModel%opt,QModel%IRC)
      pot_name_loc = QModel_in%pot_name

      IF (allocated(param_file_name_loc)) deallocate(param_file_name_loc)
    END IF

    IF (QModel%opt) write(out_unitp,*) ' Geometry optimization will be performed'
    IF (QModel%IRC) write(out_unitp,*) ' IRC will be performed'


    QModel_in%Init = .TRUE.

    IF (QModel_in%adiabatic) THEN
      IF (Print_init_loc) THEN
        write(out_unitp,*) 'Adiabatic potential . . .'
        write(out_unitp,*) 'Phase_Checking',QModel_in%Phase_Checking
        write(out_unitp,*) 'Phase_Following',QModel_in%Phase_Following
      END IF
    ELSE
      IF (Print_init_loc) write(out_unitp,*) 'Non-adiabatic potential . . .'
    END IF

    IF (QModel_in%numeric .AND. Print_init_loc) THEN
      write(out_unitp,*) 'You have decided to perform a numeric checking of the analytic formulas.'
    END IF

    CALL string_uppercase_TO_lowercase(pot_name_loc)
    IF (Print_init_loc) write(out_unitp,*) 'pot_name_loc: ',pot_name_loc

    SELECT CASE (pot_name_loc)
    CASE ('morse')
      !! === README ==
      !! Morse potential: V(R) = D*(1-exp(-a*(r-Req))**2
      !! pot_name  = 'Morse'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 1744.60504565084306291455 au
      !! remark: Default parameters for H-F
      !! === END README ==

      allocate(QML_Morse_t :: QModel%QM)
      QModel%QM = Init_QML_Morse(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('poly1d')
      !! === README ==
      !! Polynomial potential: V(R) = sum_i coef(i) * (r-Req)**i
      !! pot_name  = 'Poly1D'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 1744.60504565084306291455 au
      !! remark: Default parameters for H-F
      !! === END README ==

      allocate(QML_Poly1D_t :: QModel%QM)
      QModel%QM = Init_QML_Poly1D(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('h2')
      !! === README ==
      !! H2 potential: V(R) = Sum_i a_i * (R-Req)**(i-1)
      !! pot_name  = 'H2'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 1837.1526464003414/2 au
      !! Level: CCSD(T)-F12B/VTZ-F12 (with molpro 2010)
      !! === END README ==

      allocate(QML_H2_t :: QModel%QM)
      QModel%QM = Init_QML_H2(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('sigmoid')
      !! sigmoid function: A * 1/2(1+e*tanh((x-B)/C))  remark: e=+/-1
      allocate(QML_Sigmoid_t :: QModel%QM)
      QModel%QM = Init_QML_Sigmoid(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('buck')
      !! === README ==
      !! Buckingham potential: V(R) = A*exp(-B*R)-C/R^6
      !! pot_name  = 'buck'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 36423.484024390622 au
      !! remark: default parameters for Ar2
      !! ref:  R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283. doi:10.1098/rspa.1938.0173
      !! === END README ==
      allocate(QML_Buck_t :: QModel%QM)
      QModel%QM = Init_QML_Buck(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('hbond','linearhbond ')
      !! === README ==
      !! LinearHBond potential: Morse1(QQ/2+q,param1)+Morse2(QQ/2-q,param2)+Eref2+Buckingham(QQ)
      !! pot_name  = 'HBond'
      !! ndim      = 2   (QQ,q)
      !! nsurf     = 1
      !! reduced masses      = (/ 29156.946380706224, 1837.1526464003414 /) au
      !! remark:
      !!    A--------------H-----X--------------------B
      !!     <--------------QQ----------------------->
      !!                    <-q->
      !! ref: Dana Codruta Marinica, Marie-Pierre Gaigeot, Daniel Borgis,
      !!    Chemical Physics Letters 423 (2006) 390–394
      !!    DOI: 10.1016/j.cplett.2006.04.007
      !! ref:  Eq 3.79 of J. Beutier, thesis.
      !!
      !! remark: when option=2 is selected, a contribution is added along QQ:
      !!    Dm.exp(-betam.(QQ-QQm)) + Dp.exp(+betap.(QQ-QQp))
      !!
      !! === END README ==
      allocate(QML_LinearHBond_t :: QModel%QM)
      QModel%QM = Init_QML_LinearHBond(QModel_in,read_param=read_nml,           &
                                       nio_param_file=nio_loc)

    CASE ('2d_mb','2d_mullerbrown','mullerbrown')
      !! === README ==
      !! 2D Müller-Brown potential:
      !! pot_name  = '2d_mb'
      !! ndim      = 2   (x,y)
      !! nsurf     = 1
      !! reduced masses      = [1000.,1000.]
      !!
      !! ref:   Klaus Müller and Leo D. Brown, ...
      !!        ... Theoret. Chim. Acta (Berl.) 53, 75-93 (1979)
      !!        https://doi.org/10.1007/BF00547608.
      !!
      !! remark: the option enables one to select among three minima and two TS.
      !!         default (option=1), the first minimum (A in the reference)
      !!
      !! === END README ==
      allocate(QML_TwoD_MullerBrown_t :: QModel%QM)
      QModel%QM = Init_QML_TwoD_MullerBrown(QModel_in,read_param=read_nml,      &
                                            nio_param_file=nio_loc)

    CASE ('henonheiles')
      !! === README ==
      !! HenonHeiles potential
      !! pot_name  = 'HenonHeiles'
      !! ndim      > 1
      !! nsurf     = 1
      !! reduced masses(:)      = ONE au
      !! ref:  parameters taken from M. Nest, H.-D. Meyer, J. Chem. Phys. 117 (2002) 10499. doi:10.1063/1.1521129
      !! === END README ==
      allocate(QML_HenonHeiles_t :: QModel%QM)
      QModel%QM = Init_QML_HenonHeiles(QModel_in,                               &
                                       read_param=read_nml,nio_param_file=nio_loc)

    CASE ('tully')
      !! === README ==
      !! Tully potential: three options
      !! pot_name  = 'Tully'
      !! ndim      = 1
      !! nsurf     = 2
      !! reduced mass      = 2000. au
      !! remark: three options are possible (option = 1,2,3)
      !! ref:  Tully, J. Chem. Phys. V93, pp15, 1990
      !! === END README ==
      allocate(QML_Tully_t :: QModel%QM)
      QModel%QM = Init_QML_Tully(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('1dsoc','1dsoc_1s1t')
      !! === README ==
      !! Spin Orbit coupling model
      !! pot_name  = '1DSOC_1S1T'
      !! ndim      = 1
      !! nsurf     = 4 or 2
      !! reduced mass      = 20000. au
      !! remarks: 1 singlet and 3 triplet components                           => nsurf     = 4
      !!  or      1 singlet and 1 linear combibation of the triplet components => nsurf     = 2
      !! ref: Giovanni Granucci, Maurizio Persico, and Gloria Spighi, J. Chem. Phys. V137, p22A501 (2012)
      !! === END README ==
      allocate(QML_OneDSOC_1S1T_t :: QModel%QM)
      QModel%QM = Init_QML_OneDSOC_1S1T(QModel_in,                              &
                                        read_param=read_nml,nio_param_file=nio_loc)

    CASE ('1dsoc_2s1t')
      !! === README ==
      !! Spin Orbit coupling model
      !! pot_name  = '1DSOC_2S1T'
      !! ndim      = 1
      !! nsurf     = 4
      !! reduced mass      = 20000. au
      !! remark: 2 singlets and 1 triplet (2 linear combinations of the triplet components are not included) => nsurf     = 4
      !! ref: Giovanni Granucci, Maurizio Persico, and Gloria Spighi, J. Chem. Phys. V137, p22A501 (2012)
      !! === END README ==
      allocate(QML_OneDSOC_2S1T_t :: QModel%QM)
      QModel%QM = Init_QML_OneDSOC_2S1T(QModel_in,                              &
                                        read_param=read_nml,nio_param_file=nio_loc)

    CASE ('phenol')
      !! === README ==
      !! Phenol model
      !! pot_name  = 'Phenol'
      !! ndim      = 2 (R=rOH, th=OH-torsion)
      !! nsurf     = 3
      !! Diagonal Metric Tensor(:)      = (/ 0.0005786177, 0.0002550307 /) au
      !! remark:
      !! ref: Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218.
      !! === END README ==
      allocate(QML_Phenol_t :: QModel%QM)
      QModel%QM = Init_QML_Phenol(QModel_in,read_param=read_nml,nio_param_file=nio_loc)


    CASE ('twod')
      !! === README ==
      !! 2D model
      !! pot_name  = 'TwoD'
      !! ndim      = 2 (X,Y)
      !! nsurf     = 2
      !! Reduced masses(:)      = (/ 20000., 6667. /) au
      !! remark: The parameter values have been modified
      !! ref: A. Ferretti, G. Granucci, A. Lami, M. Persico, G. Villani, J. Chem. Phys. 104, 5517 (1996); https://doi.org/10.1063/1.471791
      !! === END README ==
      allocate(QML_TwoD_t :: QModel%QM)
      QModel%QM = Init_QML_TwoD(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('twod_rjdi2014')
      !! === README ==
      !! 2D model
      !! pot_name  = 'TwoD_RJDI2014'
      !! ndim      = 2 (X,Y)
      !! nsurf     = 2
      !! Reduced masses(:)      = [1. , 1.] au
      !! ref:  Ilya G. Ryabinkin, Loïc Joubert-Doriol, and Artur F. Izmaylov, ...
      !!       ... J. Chem. Phys. 140, 214116 (2014); https://doi.org/10.1063/1.4881147
      !! === END README ==
      allocate(QML_TwoD_RJDI2014_t :: QModel%QM)
      QModel%QM = Init_QML_TwoD_RJDI2014(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('psb3')
      !! === README ==
      !! Model for the photo-isomerization of the penta-2,4-dieniminium (PSB3) cation.
      !! pot_name  = 'PSB3'
      !! ndim      = 3
      !! nsurf     = 2
      !! remarks: two options are possible (option = 1,2)
      !! The default is option=1 (ref2).
      !! The parameters for option=2 come from the following reference.
      !! ref1: E. Marsili, M. H. Farag, X. Yang, L. De Vico, and M. Olivucci, JPCA, 123, 1710–1719 (2019).
      !!         https://doi.org/10.1021/acs.jpca.8b10010
      !! ref2: 1 E. Marsili, M. Olivucci, D. Lauvergnat, and F. Agostini, JCTC 16, 6032 (2020).
      !!        https://pubs.acs.org/doi/10.1021/acs.jctc.0c00679
      !! === END README ==
      allocate(QML_PSB3_t :: QModel%QM)
      QModel%QM = Init_QML_PSB3(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('retinal_jpcb2000','retinal_cp2000')
      !! === README ==
      !! Model for the photo-isomerization of retinal.
      !! pot_name  = 'Retinal_JPCB2000'
      !! ndim      = 2 or up to 2+23
      !! nsurf     = 2
      !! ref:  S. Hahn, G. Stock / Chemical Physics 259 (2000) 297-312.
      !!              doi: 10.1016/S0301-0104(00)00201-9
      !! === END README ==

      allocate(QML_Retinal_JPCB2000_t :: QModel%QM)
      QModel%QM = Init_QML_Retinal_JPCB2000(QModel_in,                          &
                                            read_param=read_nml,nio_param_file=nio_loc)

    CASE ('hono')
      allocate(QML_HONO_t :: QModel%QM)
      QModel%QM = Init_QML_HONO(QModel_in,read_param=read_nml,nio_param_file=nio_loc)
      !! === README ==
      !! Model for the HONO.
      !! pot_name  = 'HONO'
      !! ndim      = 6
      !! nsurf     = 1
      !! ref1:  F. Richter, M. Hochlaf, P. Rosmus, F. Gatti, and H.-D. Meyer,
      !!        J. Chem. Phys. 120, 1306 (2004).
      !!       doi: 10.1063/1.1632471
      !! ref2: F. Richter, F. Gatti, C. Léonard, F. Le Quéré, and H.-D. Meyer,
      !!       J. Chem. Phys. 127, 164315 (2007)
      !!       doi: 10.1063/1.2784553
      !! === END README ==

    CASE ('hno3')
      allocate(QML_HNO3_t :: QModel%QM)
      QModel%QM = Init_QML_HNO3(QModel_in,read_param=read_nml,  &
                                  nio_param_file=nio_loc)
    CASE ('no3')
      allocate(QML_NO3_t :: QModel%QM)
      QModel%QM = Init_QML_NO3(QModel_in,read_param=read_nml,  &
                               nio_param_file=nio_loc)
    CASE ('ch5')
      !! === README ==
      !! H + CH4 -> H-H + CH3 potential
      !!    Quadratic potential along the reaction path'
      !!    Reaction coordinate: R- = 1/2(RCH - RHH)'
      !!    Optimal coordinates along the path at CCSD(T)-F12/cc-pVTZ-F12'
      !!    V0 along the path at CCSD(T)-F12/cc-pVTZ-F12'
      !!    Hessian along the path at MP2/cc-pVDZ'
      !! pot_name  = 'CH5'
      !! ndim      = 12 or 1
      !! nsurf     = 1
      !! option = 4 (default) or 5
      !! === END README ==
      allocate(QML_CH5_t :: QModel%QM)
      QModel%QM = Init_QML_CH5(QModel_in,read_param=read_nml,  &
                                 nio_param_file=nio_loc)

    CASE ('ph4')
      !! === README ==
      !! H + PH3 -> H-H + PH2 potential
      !!    Quadratic potential along the reaction path'
      !!    Reaction coordinate: R- = 1/2(RPH - RHH)'
      !!    Optimal coordinates along the path at MP2/cc-pVTZ'
      !!    V0 along the path at CCSD(T)-F12/cc-pVTZ-F12 (option 4) or ...'
      !!      ... MP2/cc-pVTZ'
      !!    Hessian and gradient along the path at MP2/cc-pVTZ'
      !! pot_name  = 'PH4'
      !! ndim      = 9 or 1
      !! nsurf     = 1
      !! option = 4 (default)
      !! === END README ==
      allocate(QML_PH4_t :: QModel%QM)
      QModel%QM = Init_QML_PH4(QModel_in,read_param=read_nml,  &
                                 nio_param_file=nio_loc)

    CASE ('hnnhp')
      allocate(QML_HNNHp_t :: QModel%QM)
      QModel%QM = Init_QML_HNNHp(QModel_in,read_param=read_nml, &
                                   nio_param_file=nio_loc)

    CASE ('h2sin')
      allocate(QML_H2SiN_t :: QModel%QM)
      QModel%QM = Init_QML_H2SiN(QModel_in,read_param=read_nml, &
                                   nio_param_file=nio_loc)

    CASE ('h2nsi')
      allocate(QML_H2NSi_t :: QModel%QM)
      QModel%QM = Init_QML_H2NSi(QModel_in,read_param=read_nml, &
                                   nio_param_file=nio_loc)

    CASE ('hoo_dmbe')
      !! === README ==
      !! HOO potential: DMBE IV of Varandas group
      !! pot_name  = 'HOO_DMBE'
      !! ndim      = 3   (R1=dOO,R2=dHO1,R3=dHO2)
      !! nsurf     = 1
      !! ref:    M. R. Pastrana, L. A. M. Quintales, J. Brandão and A. J. C. Varandas'
      !!         JCP, 1990, 94, 8073-8080, doi: 10.1021/j100384a019.
      !! === END README ==
      allocate(QML_HOO_DMBE_t :: QModel%QM)
      QModel%QM = Init_QML_HOO_DMBE(QModel_in,read_param=read_nml,      &
                                      nio_param_file=nio_loc)

    CASE ('h3_lsth','h3')
      !! === README ==
      !! H3 potential:
      !! pot_name  = 'H3'
      !! option    = 0,1,10,11 (LSTH)
      !! ndim      = 3   (the 3 H-H distances)
      !! nsurf     = 1
      !! refs (option=0):
      !! P. Siegbahn, B. Liu,  J. Chem. Phys. 68, 2457(1978).
      !! D.G. Truhlar and C.J. Horowitz, J. Chem. Phys. 68, 2466 (1978); https://doi.org/10.1063/1.436019
      !! options  0 and 10 : 3D model with IRC functions (1 potential + parameters)
      !! options  0 and  1 : first IRC funtions fitted in polar representation (alpha)
      !! options 10 and 11 : second IRC funtions fitted with the sum and the difference (alpha)
      !! === END README ==
      allocate(QML_H3_t :: QModel%QM)
      QModel%QM = Init_QML_H3(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('hcn_murrell')
      !! === README ==
      !! H3 potential:
      !! pot_name  = 'HCN_Murrell'
      !! option    = 0 (3D-3distances), 1 (3D-Jacobi), 3 IRC+RPH
      !! ndim      = 3
      !! nsurf     = 1
      !! J. N. Murrell, S. Carter and L. O. Halonene, J. Mol. Spectrosc. vo93 p307 1982
      !!  doi: https://doi.org/10.1016/0022-2852(82)90170-9
      !! === END README ==
      allocate(QML_HCN_Murrell_t :: QModel%QM)
      QModel%QM = Init_QML_HCN_Murrell(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('clh2+','clh2p')
      !! === README ==
      !! ClH2+ potential:
      !! pot_name  = 'ClH2+'
      !! option    = 1 to 6 (default 6)
      !! ndim      = 3
      !! nsurf     = 1
      !! remark, 6 options are possible:
      !!    options = 1,3,5: the coordinates are [angle, R+, R-] (in bohr and radian)
      !!    options = 2,4,6: the coordinates are [R1, R2, angle] (in bohr and radian)
      !!
      !!    options = 1,2: B3LYP/cc-pVTZ 1st version (do not use)
      !!    options = 3,4: B3LYP/cc-pVTZ 2d  version
      !!    options = 5,6: CCSD(T)-F12b/cc-pVTZ-F12
      !! refs: unpublished
      !! === END README ==
      allocate(QML_ClH2p_t :: QModel%QM)
      QModel%QM = Init_QML_ClH2p(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('clh2+_botschwina','clh2p_botschwina')
      !! === README ==
      !! ClH2+ potential:
      !! pot_name  = 'ClH2+_Botschwina'
      !! option    = 1,2 (default 1)
      !! ndim      = 3
      !! nsurf     = 1
      !! remark, 2 options are possible:
      !!    options = 1: CEPA-1
      !!    options = 2: CEPA-1 corrected
      !! ref: Peter Botschwina, J. Chem. Soc., Faraday Trans. 2, 1988, 84(9), 1263-1276'
      !!      DOI: 10.1039/F29888401263
      !! === END README ==
      allocate(QML_ClH2p_Botschwina_t :: QModel%QM)
      QModel%QM = Init_QML_ClH2p_Botschwina(QModel_in,read_param=read_nml,      &
                                            nio_param_file=nio_loc)

    CASE ('template')
      !! 3D-potential with 1 surface
      allocate(QML_Template_t :: QModel%QM)
      QModel%QM = Init_QML_Template(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('test')
      !! test-potential
      allocate(QML_Template_t :: QModel%QM)
      QModel%QM = Init_QML_Test(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE DEFAULT
        write(out_unitp,*) ' ERROR in Init_Model'
        write(out_unitp,*) ' This model/potential is unknown. pot_name: ',pot_name_loc
        STOP 'STOP in Init_Model: Other potentials have to be done'
    END SELECT

    IF (present(ndim)) THEN
      IF (ndim > QModel%QM%ndim) THEN
          write(out_unitp,*) ' ERROR in Init_Model'
          write(out_unitp,*) ' ndim is present and ...'
          write(out_unitp,*) ' its value is larger than QModel%QM%ndim'
          write(out_unitp,*) ' ndim,QModel%QM%ndim',ndim,QModel%QM%ndim
          write(out_unitp,*) ' check your data!'
          STOP 'STOP in Init_Model: wrong ndim'
      END IF
      IF (ndim < QModel%QM%ndim  .AND. ndim > 0) THEN
          write(out_unitp,*) ' WARNING in Init_Model'
          write(out_unitp,*) ' ndim is present and ...'
          write(out_unitp,*) ' its value is smaller than QModel%QM%ndim'
          write(out_unitp,*) ' ndim,QModel%QM%ndim',ndim,QModel%QM%ndim
          write(out_unitp,*) ' => We assume that ...'
          write(out_unitp,*) ' ... all variables (QModel%QM%ndim) will be given!'
      END IF
    END IF
    IF (present(nsurf)) THEN
      IF (nsurf /= QModel%QM%nsurf .AND. nsurf > 0) THEN
          write(out_unitp,*) ' ERROR in Init_Model'
          write(out_unitp,*) ' nsurf is present and ...'
          write(out_unitp,*) ' its value is not equal to QModel%QM%nsurf'
          write(out_unitp,*) ' nsurf,QModel%QM%nsurf',nsurf,QModel%QM%nsurf
          write(out_unitp,*) ' check your data!'
          STOP 'STOP in Init_Model: wrong nsurf'
        END IF
    END IF

    QModel%ndim  = QModel%QM%ndim
    QModel%nsurf = QModel%QM%nsurf

    ! special feature when Vib_adia = .TRUE.
    IF (QModel%QM%Vib_adia) THEN

      ! check the value of list_act (>0 and <= ndim)
      IF (any(QModel%QM%list_act < 1)           .OR. &
          any(QModel%QM%list_act > QModel%QM%ndim)) THEN
        write(out_unitp,*) ' ERROR in Init_Model'
        write(out_unitp,*) ' Some values of list_act(:) are out of range.'
        write(out_unitp,*) '   list_act(:): ',QModel%QM%list_act(:)
        write(out_unitp,*) '   range = [1,',QModel%QM%ndim,']'
        write(out_unitp,*) ' check your data!'
        STOP 'STOP in Init_Model: wrong list_act'
      END IF

      allocate(list_Q(QModel%QM%ndim))
      list_Q(:) = .FALSE.
      DO i=1,size(QModel%QM%list_act)
        list_Q(QModel%QM%list_act(i)) = .TRUE.
      END DO
      i_inact  = 0
      nb_inact = QModel%QM%ndim - size(QModel%QM%list_act)
      allocate(QModel%QM%list_inact(nb_inact))

      DO i=1,QModel%QM%ndim
        IF (.NOT. list_Q(i)) THEN
          i_inact = i_inact + 1
          QModel%QM%list_inact(i_inact) = i
          list_Q(i) = .TRUE.
        END IF
      END DO

      IF (count(list_Q) /= QModel%QM%ndim) THEN
        write(out_unitp,*) ' ERROR in Init_Model'
        write(out_unitp,*) ' Some coordinate indexes are missing in ...'
        write(out_unitp,*) ' ... list_act(:):   ',QModel%QM%list_act(:)
        write(out_unitp,*) ' and list_inact(:): ',QModel%QM%list_inact(:)
        write(out_unitp,*) ' check your data!'
        STOP 'STOP in Init_Model: wrong list_act'
      END IF

      write(out_unitp,*) 'Vib_adia   ',QModel%QM%Vib_adia
      write(out_unitp,*) 'nb_Channels',QModel%QM%nb_Channels
      write(out_unitp,*) 'list_act   ',QModel%QM%list_act
      write(out_unitp,*) 'list_inact ',QModel%QM%list_inact

      QModel%ndim  = size(QModel%QM%list_act)
      QModel%nsurf = QModel%QM%nb_Channels

      allocate(QModel%Basis)
      CALL Read_Basis(QModel%Basis,nio_loc)

    END IF

    IF (Print_init_loc) THEN
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Quantum Model'
      CALL Write_Model(QModel,nio=out_unitp)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '================================================='
      flush(out_unitp)
    END IF

    IF (read_param_loc .AND. nio_loc /= in_unitp) THEN
       close(unit=nio_loc)
    END IF

  END SUBROUTINE Init_Model
  SUBROUTINE set_step_epsi_Model(step_in,epsi_in)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(in), optional  :: step_in,epsi_in


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='set_step_epsi_Model'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    IF (present(step_in)) THEN
      write(out_unitp,*) ' WARNING: step has been changed.'
      write(out_unitp,*) ' Old and new values',step,step_in
      step = step_in
    END IF

    IF (present(epsi_in)) THEN
      write(out_unitp,*) ' WARNING: epsi has been changed.'
      write(out_unitp,*) ' Old and new values',epsi,epsi_in
      epsi = epsi_in
    END IF

    IF (debug) THEN
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE set_step_epsi_Model

  SUBROUTINE get_Q0_Model(Q0,QModel,option)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q0(:)
    TYPE (Model_t),     intent(in)               :: QModel
    integer,            intent(in)               :: option


    integer :: err_Q0
!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='get_Q0_Model'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL check_alloc_QM(QModel,name_sub)

    IF (size(Q0) /= QModel%QM%ndim) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The size of Q0 is not QModel%QM%ndim: '
      write(out_unitp,*) ' size(Q0)',size(Q0)
      write(out_unitp,*) ' ndim',QModel%QM%ndim
      STOP 'STOP in get_Q0_Model: Wrong Q0 size'
    END IF

    Q0(:) = ZERO

    CALL get_Q0_QModel(QModel%QM,Q0,err_Q0)
    IF (err_Q0 /= 0) THEN
      CALL Write_Model(QModel,out_unitp)
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Q0 is not set-up in the model'
      STOP 'STOP Q0 is not set-up in the model'
    END IF

    IF (debug) THEN
      CALL Write_RVec(Q0,out_unitp,nbcol1=5,name_info='Q0: ')
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE get_Q0_Model

  SUBROUTINE check_QML_Path()
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

  character (len=:), allocatable :: FileName
  logical :: file_exist


  FileName = make_FileName('InternalData/Test_QML_Path.txt')

  inquire(file=FileName,exist=file_exist)

  IF (.NOT. file_exist) THEN
    write(out_unitp,*) 'ERROR: the QML directory path is wrong !!'
    write(out_unitp,*) ' FileName: ',FileName
    write(out_unitp,*) ' QML_path: ',trim(adjustl(QML_path))
    write(out_unitp,*) ' Probably, the QML directory has been moved'
    write(out_unitp,*) ' Recompile again QML.'
    STOP 'Wrong QML_path'
  END IF

  END SUBROUTINE check_QML_Path

  ! check if the QM [CLASS(QML_Empty_t)] is allocated
  SUBROUTINE check_alloc_QM(QModel,name_sub_in)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    TYPE (Model_t),     intent(in)     :: QModel
    character (len=*),  intent(in)     :: name_sub_in

    IF ( .NOT. allocated(QModel%QM)) THEN
      write(out_unitp,*) ' ERROR in check_alloc_QM'
      write(out_unitp,*) ' QM is not allocated in QModel.'
      write(out_unitp,*) '  check_alloc_QM is called from ',name_sub_in
      write(out_unitp,*) '  You MUST initialize the model with:'
      write(out_unitp,*) '    CALL init_Model(...) in Model_m.f90'
      write(out_unitp,*) ' or'
      write(out_unitp,*) '    CALL sub_Init_Qmodel(...) in Model_driver.f90'
      STOP 'STOP in check_alloc_QM: QM is not allocated in QModel.'
    END IF

  END SUBROUTINE check_alloc_QM

  ! check if the check_Init_QModel [TYPE(Model_t)] is initialized
  FUNCTION check_Init_QModel(QModel)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    logical                             :: check_Init_QModel
    TYPE (Model_t),     intent(in)      :: QModel

    check_Init_QModel = allocated(QModel%QM)
    IF (allocated(QModel%QM)) THEN
      check_Init_QModel = QModel%QM%init
    END IF

  END FUNCTION check_Init_QModel

  SUBROUTINE Eval_tab_HMatVibAdia(QModel,Qact,tab_MatH)
  USE ADdnSVM_m
  IMPLICIT NONE

  TYPE (Model_t),                 intent(inout)            :: QModel
  real (kind=Rkind),              intent(in)               :: Qact(:)
  real (kind=Rkind), allocatable, intent(inout)            :: tab_MatH(:,:,:)



  integer                        :: nb_terms,nsurf

  integer                        :: ia,ja,iterm,ii,ji,i,iq,ib,jb,nb,nq

  TYPE (dnMat_t)                 :: PotVal_dia,PotVal,Vec,NAC,Mat_diag
  LOGICAL                        :: PF ! phase_following
  LOGICAL                        :: PC ! phase checking
  real (kind=Rkind), allocatable :: Q(:)
  real (kind=Rkind), allocatable :: d0GGdef_ii(:,:),d0GGdef_aa(:,:),d0GGdef_ai(:,:),d0GGdef_ia(:,:)

  integer                        :: ndim_act,ndim_inact

  TYPE (dnS_t), allocatable      :: dnV(:),dnHB(:)
  TYPE (dnS_t)                   :: dnVfull,dnHij
  TYPE (dnMat_t)                 :: dnH ! derivative of the Hamiltonian


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_tab_HMatVibAdia'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    flush(out_unitp)
  END IF

  CALL check_alloc_QM(QModel,name_sub)

  PF = QModel%QM%Phase_Following
  PC = QModel%QM%Phase_Checking
  CALL Eval_dnHVib_ana(QModel,Qact,PotVal_dia,nderiv=2)
  IF (.NOT. allocated(QModel%QM%Vec0)) allocate(QModel%QM%Vec0)
  CALL dia_TO_adia(PotVal_dia,PotVal,Vec,QModel%QM%Vec0,NAC,PF,PC,nderiv=2)

  !CALL Write_dnMat(PotVal,nio=out_unitp,info='PotVal (adia)')

  !Mat_diag = matmul(transpose(Vec),matmul(PotVal_dia,Vec))
  !CALL Write_dnMat(Mat_diag,nio=out_unitp,info='Mat_diag')

  !write(out_unitp,*) 'nsurf,ndim',QModel%nsurf,QModel%ndim
  nsurf    = QModel%nsurf
  nb_terms = (QModel%ndim + 1)*(QModel%ndim + 2)/2
  IF (.NOT. allocated(tab_MatH)) THEN
    allocate(tab_MatH(nsurf,nsurf,nb_terms))
  END IF

  !write(out_unitp,*) Qact,'NAC1',NAC%d1(1:nsurf,1:nsurf,:)


  IF ( any([nsurf,nsurf,nb_terms] /= shape(tab_MatH)) ) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' The shape of tab_MatH is wrong.'
    write(out_unitp,*) '    shape(tab_MatH):',shape(tab_MatH)
    write(out_unitp,*) '  It MUST be:      [',nsurf,nsurf,nb_terms,']'
    write(out_unitp,*) '  Check your data or the code!'
    STOP 'STOP in Eval_tab_HMatVibAdia: The shape of tab_MatH is wrong.'
  END IF
  tab_MatH(:,:,:) = ZERO

  ndim_act    = size(QModel%QM%list_act)
  ndim_inact  = size(QModel%QM%list_inact)

  ! Veff
  d0GGdef_aa = QModel%QM%d0GGdef(QModel%QM%list_act,QModel%QM%list_act)
  iterm = 1
  tab_MatH(:,:,iterm) = PotVal%d0(1:nsurf,1:nsurf)
  DO ia=1,ndim_act
  DO ja=1,ndim_act
    tab_MatH(:,:,iterm) = tab_MatH(:,:,iterm) - HALF*d0GGdef_aa(ja,ia) *        &
       NAC%d2(1:nsurf,1:nsurf,ja,ia)
  END DO
  END DO

  ! F2^jaia
  d0GGdef_aa = QModel%QM%d0GGdef(QModel%QM%list_act,QModel%QM%list_act)
  DO ia=1,ndim_act
  DO ja=ia,ndim_act
    iterm = iterm + 1
    IF (ja == ia) THEN
      tab_MatH(:,:,iterm) = -HALF*d0GGdef_aa(ja,ia) * NAC%d0(1:nsurf,1:nsurf)
    ELSE
      tab_MatH(:,:,iterm) = -d0GGdef_aa(ja,ia) * NAC%d0(1:nsurf,1:nsurf)
    END IF
  END DO
  END DO

  ! F1^ia
  d0GGdef_aa = QModel%QM%d0GGdef(QModel%QM%list_act,QModel%QM%list_act)
  DO ia=1,ndim_act
    iterm = iterm + 1
    tab_MatH(:,:,iterm) = ZERO
    DO ja=1,ndim_act
      tab_MatH(:,:,iterm) = tab_MatH(:,:,iterm) - d0GGdef_aa(ja,ia) * &
           NAC%d1(1:nsurf,1:nsurf,ja)
    END DO
    !the contribution d0GGdef_ai is missing
  END DO

  IF (debug) THEN
    DO iterm=1,nb_terms
      write(out_unitp,*) iterm
      CALL Write_RMat(tab_MatH(:,:,iterm),nio=out_unitp,nbcol1=5)
    END DO
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

  END SUBROUTINE Eval_tab_HMatVibAdia

  SUBROUTINE Eval_Pot(QModel,Q,PotVal,nderiv,NAC,Vec,numeric)
  USE ADdnSVM_m
  IMPLICIT NONE

    TYPE (Model_t),     intent(inout)            :: QModel
    TYPE (dnMat_t),     intent(inout)            :: PotVal
    real (kind=Rkind),  intent(in)               :: Q(:)
    integer,            intent(in),    optional  :: nderiv
    TYPE (dnMat_t),     intent(inout), optional  :: NAC,Vec
    logical,            intent(in),    optional  :: numeric

    ! local variables
    integer                    :: i,nderiv_loc
    TYPE (dnMat_t)             :: Vec_loc,NAC_loc,PotVal_dia,PotVal_loc
    logical                    :: numeric_loc,adia_loc
    logical                    :: PF ! phase_following
    logical                    :: PC ! phase_checking

    !real (kind=Rkind), allocatable :: G(:,:)

    integer :: numeric_option = 3   ! 0 old (up to 2d derivatives
                                    ! 3 version up to 3d derivatives less points than 4
                                    ! 4 version up to 3d derivatives more points than 3

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) ' BEGINNING ',name_sub
    IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
    flush(out_unitp)
  END IF

  CALL check_alloc_QM(QModel,name_sub)
  PF = QModel%QM%Phase_Following
  PC = QModel%QM%Phase_Checking
  IF (debug) THEN
    write(out_unitp,*) '  QModel%QM%numeric         ',QModel%QM%numeric
    write(out_unitp,*) '  QModel%QM%adiabatic       ',QModel%QM%adiabatic
    write(out_unitp,*) '  QModel%QM%Vib_adia        ',QModel%QM%Vib_adia
    write(out_unitp,*) '  QModel%QM%Phase_Following ',QModel%QM%Phase_Following
    write(out_unitp,*) '  QModel%QM%Phase_Checking  ',QModel%QM%Phase_Checking
    flush(out_unitp)
  END IF

  IF (present(nderiv)) THEN
    nderiv_loc = max(0,nderiv)
    nderiv_loc = min(3,nderiv_loc)
  ELSE
    nderiv_loc = 0
  END IF

  IF (present(numeric)) THEN
    numeric_loc = (numeric  .OR. QModel%QM%no_ana_der)
  ELSE
    numeric_loc = (QModel%QM%numeric .OR. QModel%QM%no_ana_der)
  END IF
  numeric_loc = (numeric_loc .AND. nderiv_loc > 0)


  IF (QModel%QM%Vib_adia) THEN
    IF (present(Vec)) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Vib_adia=t and Vec is present'
      write(out_unitp,*) ' This is not possible yet !'
      STOP 'ERROR in Eval_Pot: Vib_adia=t is not compatible with Vec'
    END IF

    CALL Eval_dnHVib_ana(QModel,Q,PotVal_dia,nderiv_loc)

    !write(out_unitp,*) 'PotVal (Vib_dia)'
    !CALL Write_dnMat(PotVal_dia,nio=out_unitp)

    IF (.NOT. allocated(QModel%QM%Vec0)) allocate(QModel%QM%Vec0)
    CALL dia_TO_adia(PotVal_dia,PotVal_loc,Vec_loc,QModel%QM%Vec0,NAC_loc,      &
                     PF,PC,nderiv_loc,type_diag=1)

    CALL submatrix_dnMat2_TO_dnMat1(PotVal,PotVal_loc,lb=1,ub=QModel%QM%nb_Channels)

    IF (present(Vec)) THEN
      ! it needs dnMat as a rectangular matrix !!
      CALL submatrix_dnMat2_TO_dnMat1(Vec,Vec_loc,lb=1,ub=QModel%QM%nb_Channels)
    END IF

    IF (present(NAC)) THEN
      CALL submatrix_dnMat2_TO_dnMat1(NAC,NAC_loc,lb=1,ub=QModel%QM%nb_Channels)
    END IF

    ! print the Vec%d0 if required
    CALL Write_Vec(Q,Vec_loc,QModel,nio=out_unitp)

    CALL dealloc_dnMat(NAC_loc)
    CALL dealloc_dnMat(Vec_loc)
    CALL dealloc_dnMat(PotVal_loc)
    CALL dealloc_dnMat(PotVal_dia)

  ELSE

    adia_loc = (QModel%QM%adiabatic .AND. QModel%QM%nsurf > 1)

    IF (numeric_loc) THEN  ! numerical
      IF (.NOT. adia_loc) THEN
         SELECT CASE (numeric_option)
         CASE (0)
           CALL Eval_Pot_Numeric_dia_old(QModel,Q,PotVal,nderiv_loc)
         CASE (3)
           CALL Eval_Pot_Numeric_dia_v3(QModel,Q,PotVal,nderiv_loc)
         CASE (4)
           CALL Eval_Pot_Numeric_dia_v4(QModel,Q,PotVal,nderiv_loc)
         CASE Default
           CALL Eval_Pot_Numeric_dia_old(QModel,Q,PotVal,nderiv_loc)
         END SELECT
      ELSE

        IF (present(Vec)) THEN
          IF (present(NAC)) THEN
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,      &
                                                 Vec,NAC,numeric_option)
          ELSE
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,      &
                                             Vec,NAC_loc,numeric_option)
            CALL dealloc_dnMat(NAC_loc)
          END IF
        ELSE
          IF (present(NAC)) THEN
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,      &
                                             Vec_loc,NAC,numeric_option)
          ELSE
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,      &
                                         Vec_loc,NAC_loc,numeric_option)
            CALL dealloc_dnMat(NAC_loc)
          END IF
          CALL dealloc_dnMat(Vec_loc)
        END IF
      END IF
    ELSE ! analytical calculation
      IF (present(Vec)) THEN
        IF (present(NAC)) THEN
          CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv_loc,Vec=Vec,Nac=NAC)
        ELSE
          CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv_loc,Vec=Vec)
        END IF
      ELSE
        IF (present(NAC)) THEN
          CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv_loc,Nac=NAC)
        ELSE
          CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv_loc)
        END IF
      END IF

    END IF
  END IF

  IF (debug) THEN
    IF ( QModel%QM%adiabatic) THEN
      write(out_unitp,*) 'PotVal (adia)'
    ELSE
      write(out_unitp,*) 'PotVal (dia)'
    END IF
    CALL Write_dnMat(PotVal,nio=out_unitp)
    write(out_unitp,*) ' END ',name_sub
    flush(out_unitp)
  END IF

  END SUBROUTINE Eval_Pot

  SUBROUTINE Eval_Pot_ana(QModel,Q,PotVal,nderiv,NAC,Vec)
  USE ADdnSVM_m, ONLY : dnS_t,alloc_dnS,dealloc_dnS,Variable, &
     dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat,Write_dnMat
  IMPLICIT NONE

    TYPE (Model_t),        intent(inout)            :: QModel

    TYPE (dnMat_t),        intent(inout)            :: PotVal
    real (kind=Rkind),     intent(in)               :: Q(:)
    integer,               intent(in)               :: nderiv
    TYPE (dnMat_t),        intent(inout), optional  :: NAC,Vec

    ! local variables
    integer                     :: i,j,ij,id,nat
    TYPE (dnMat_t)              :: PotVal_dia,Vec_loc,NAC_loc
    TYPE (dnS_t), allocatable   :: dnQ(:)
    TYPE (dnS_t), allocatable   :: dnX(:,:)
    TYPE (dnS_t), allocatable   :: Mat_OF_PotDia(:,:)
    logical                     :: PF,PC ! phase_following,Phase_Checking

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot_ana'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      write(out_unitp,*) '   nderiv:       ',nderiv
      write(out_unitp,*) '   present(NAC): ',present(NAC)
      write(out_unitp,*) '   present(Vec): ',present(Vec)

      flush(out_unitp)
    END IF

    CALL check_alloc_QM(QModel,name_sub)

    PF = QModel%QM%Phase_Following
    PC = QModel%QM%Phase_Checking
    IF (debug) write(out_unitp,*) '   adiabatic ',QModel%QM%adiabatic


    IF ( Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,     &
                           nderiv=nderiv)
    END IF
    PotVal = ZERO
    IF (debug) write(out_unitp,*) '   init PotVal  ' ; flush(out_unitp)

    ! allocate Mat_OF_PotDia
    allocate(Mat_OF_PotDia(QModel%QM%nsurf,QModel%QM%nsurf))
    DO j=1,size(Mat_OF_PotDia(1,:))
    DO i=1,size(Mat_OF_PotDia(:,1))
        CALL alloc_dnS(Mat_OF_PotDia(i,j),QModel%QM%ndim,nderiv)
    END DO
    END DO
    IF (debug) write(out_unitp,*) '   alloc Mat_OF_PotDia  ' ; flush(out_unitp)

    ! intialization of the dnQ(:)
    IF (QModel%QM%Cart_TO_Q) THEN
      !in Q(:) we have the cartesian coordinates
      allocate(dnQ(QModel%QM%ndimQ))
      nat = int(QModel%QM%ndim/3)
      allocate(dnX(3,nat))
      ij = 0
      DO i=1,nat
      DO j=1,3
        ij = ij + 1
        dnX(j,i) = Variable(Q(ij),nVar=QModel%QM%ndim,nderiv=nderiv,iVar=ij) ! to set up the derivatives
      END DO
      END DO

      CALL QModel%QM%Cart_TO_Q_QModel(dnX,dnQ,nderiv=nderiv)

      CALL dealloc_dnS(dnX)
      deallocate(dnX)

    ELSE
      allocate(dnQ(QModel%QM%ndim))
      DO i=1,QModel%QM%ndim
        dnQ(i) = Variable(Q(i),nVar=QModel%QM%ndim,nderiv=nderiv,iVar=i) ! to set up the derivatives
      END DO
    END IF
    IF (debug) write(out_unitp,*) '   init dnQ(:)  ' ; flush(out_unitp)

    CALL QModel%QM%EvalPot_QModel(Mat_OF_PotDia,dnQ,nderiv=nderiv)

    PotVal = Mat_OF_PotDia ! transfert the potential and its derivatives to the matrix form (PotVal)

    ! deallocation
    DO i=1,size(dnQ)
      CALL dealloc_dnS(dnQ(i))
    END DO
    deallocate(dnQ)

    DO j=1,size(Mat_OF_PotDia(1,:))
    DO i=1,size(Mat_OF_PotDia(:,1))
      CALL dealloc_dnS(Mat_OF_PotDia(i,j))
    END DO
    END DO
    deallocate(Mat_OF_PotDia)
    ! end deallocation

    IF ( QModel%QM%adiabatic .AND. QModel%QM%nsurf > 1) THEN
      IF (debug) THEN
        write(out_unitp,*) 'PotVal (dia)'
        CALL Write_dnMat(PotVal,nio=out_unitp)
        flush(out_unitp)
      END IF
      IF (.NOT. allocated(QModel%QM%Vec0)) allocate(QModel%QM%Vec0)

      PotVal_dia = PotVal

      IF (present(Vec)) THEN
        IF (present(NAC)) THEN
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec,QModel%QM%Vec0,NAC,PF,PC,nderiv)
        ELSE
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec,QModel%QM%Vec0,NAC_loc,PF,PC,nderiv)
          CALL dealloc_dnMat(NAC_loc)
        END IF
      ELSE
        IF (present(NAC)) THEN
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec_loc,QModel%QM%Vec0,NAC,PF,PC,nderiv)
        ELSE
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec_loc,QModel%QM%Vec0,NAC_loc,PF,PC,nderiv)
          CALL dealloc_dnMat(NAC_loc)
        END IF
        CALL dealloc_dnMat(Vec_loc)
      END IF
      CALL dealloc_dnMat(PotVal_dia)
    END IF


    IF (debug) THEN
      IF ( QModel%QM%adiabatic) THEN
        write(out_unitp,*) 'PotVal (adia)'
      ELSE
        write(out_unitp,*) 'PotVal (dia)'
      END IF
      CALL Write_dnMat(PotVal,nio=out_unitp)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Eval_Pot_ana

  SUBROUTINE Eval_Pot_Numeric_dia_v4(QModel,Q,PotVal,nderiv)
  USE QMLLib_FiniteDiff_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j,k,ip,jp,kp
    integer                            :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_dia_v4')


    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,nderiv=0)
    PotVal = ZERO

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0)
    CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=4)

    IF (nderiv >= 1) THEN ! along ONE coordinates (first derivatives and higher)

      ! Numeric evaluation of forces
      DO i=1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                          indQ=[i],indDQ=ind1DQ,option=4)
        END DO

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=1,QModel%QM%ndim
        IF (i == j) CYCLE

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                        indQ=[i,j],indDQ=ind2DQ,option=4)
        END DO

      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN ! 3d derivatives:d3/dQidQjdQk

      DO i=1,QModel%QM%ndim
      DO j=1,QModel%QM%ndim
      IF (i == j) CYCLE
      DO k=1,QModel%QM%ndim
        IF (i == k .OR. j == k) CYCLE

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                     indQ=[i,j,k],indDQ=ind3DQ,option=4)
        END DO

      END DO
      END DO
      END DO
    END IF

    CALL FiniteDiff_Finalize_dnMat(PotVal,step)

    deallocate(Q_loc)
    CALL dealloc_dnMat(PotVal_loc0)

  END SUBROUTINE Eval_Pot_Numeric_dia_v4
  SUBROUTINE Eval_Pot_Numeric_dia_v3(QModel,Q,PotVal,nderiv)
  USE QMLLib_FiniteDiff_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j,k,ip,jp,kp

    integer                            :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_dia_v3')


    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO


    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0)
    PotVal%d0 = PotVal_loc0%d0

    CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=3)


    IF (nderiv >= 1) THEN ! 1st derivatives

      DO i=1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                          indQ=[i],indDQ=ind1DQ,option=3)
        END DO

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                        indQ=[i,j],indDQ=ind2DQ,option=3)
        END DO

        CALL FiniteDiff3_SymPerm_OF_dnMat(PotVal,indQ=[i,j])

      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN ! 3d derivatives: d3/dQidQidQj

      ! d3/dQidQjdQk
      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim
      DO k=j+1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                      indQ=[i,j,k],indDQ=ind3DQ,option=3)
        END DO

        CALL FiniteDiff3_SymPerm_OF_dnMat(PotVal,indQ=[i,j,k])

      END DO
      END DO
      END DO
    END IF

    CALL FiniteDiff_Finalize_dnMat(PotVal,step)

    deallocate(Q_loc)
    CALL dealloc_dnMat(PotVal_loc0)

  END SUBROUTINE Eval_Pot_Numeric_dia_v3

  SUBROUTINE Eval_Pot_Numeric_dia_old(QModel,Q,PotVal,nderiv)
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_dia_old')


    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0)
    PotVal%d0 = PotVal_loc0%d0


    IF (nderiv >= 1) THEN ! 1st derivatives

      ! Numeric evaluation of forces
      DO i=1,QModel%QM%ndim

        Q_loc(i) = Q(i) + step        ! q+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0) ! Ep(q+dq)
        PotVal%d1(:,:,i) = PotVal_loc0%d0

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = PotVal_loc0%d0
        END IF

        Q_loc(i) = Q(i) - step        ! q-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0) ! Ep(q-dq)
        PotVal%d1(:,:,i) = (PotVal%d1(:,:,i)-PotVal_loc0%d0)/(TWO*step)


        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = (PotVal%d2(:,:,i,i) + PotVal_loc0%d0 - TWO*PotVal%d0)/ &
                                step**2
        END IF

        Q_loc(i) = Q(i)

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) + PotVal_loc0%d0

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0

        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i)/(FOUR*step**2)
        PotVal%d2(:,:,i,j) = PotVal%d2(:,:,j,i)

        Q_loc(i) = Q(i)
        Q_loc(j) = Q(j)
      END DO
      END DO
    END IF

    deallocate(Q_loc)
    CALL dealloc_dnMat(PotVal_loc0)

  END SUBROUTINE Eval_Pot_Numeric_dia_old

  SUBROUTINE Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv,Vec,NAC,option)
  USE ADdnSVM_m, ONLY : dnMat_t
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv
    TYPE (dnMat_t),    intent(inout)  :: Vec,NAC
    integer,           intent(in)     :: option


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot_Numeric_adia'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      write(out_unitp,*) '   nderiv',nderiv
      write(out_unitp,*) '   option',option
      flush(out_unitp)
    END IF

    SELECT CASE (option)
    CASE (0)
      CALL Eval_Pot_Numeric_adia_old(QModel,Q,PotVal,nderiv,Vec,NAC)
    CASE (3)
      CALL Eval_Pot_Numeric_adia_v3(QModel,Q,PotVal,nderiv,Vec,NAC)
    CASE (4)
      STOP 'Eval_Pot_Numeric_adia: option=4, not yet'
    !  CALL Eval_Pot_Numeric_adia_v4(QModel,Q,PotVal,nderiv,Vec,NAC)
    CASE Default
      CALL Eval_Pot_Numeric_adia_old(QModel,Q,PotVal,nderiv,Vec,NAC)
    END SELECT

    IF (debug) THEN
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF
  END SUBROUTINE Eval_Pot_Numeric_adia
  SUBROUTINE Eval_Pot_Numeric_adia_old(QModel,Q,PotVal,nderiv,Vec,NAC)
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv
    TYPE (dnMat_t),    intent(inout)  :: Vec,NAC

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0,Vec_loc0
    integer                            :: i,j
    real (kind=Rkind), allocatable     :: tVec(:,:)

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_adia_old')


    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                           nderiv=nderiv)
    END IF
    PotVal = ZERO

    IF (Check_NotAlloc_dnMat(Vec,nderiv) ) THEN
      CALL alloc_dnMat(Vec,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    Vec = ZERO

    IF (Check_NotAlloc_dnMat(NAC,nderiv) ) THEN
      CALL alloc_dnMat(NAC,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    NAC = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,nderiv=0)
    CALL alloc_dnMat(Vec_loc0,   nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0,vec=Vec_loc0)

    PotVal%d0 = PotVal_loc0%d0
    Vec%d0    = Vec_loc0%d0
    CALL Init_IdMat(NAC%d0,QModel%QM%nsurf)

    allocate(tVec(QModel%QM%nsurf,QModel%QM%nsurf))
    tVec(:,:)      = transpose(Vec%d0)

    IF (nderiv >= 1) THEN ! 1st derivatives

      ! Numeric evaluation of forces
      DO i=1,QModel%QM%ndim

        Q_loc(i) = Q(i) + step        ! q+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0) ! Ep(q+dq)
        PotVal%d1(:,:,i) = PotVal_loc0%d0
        Vec%d1(:,:,i)    = Vec_loc0%d0

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = PotVal_loc0%d0
          Vec%d2(:,:,i,i)    = Vec_loc0%d0
        END IF

        Q_loc(i) = Q(i) - step        ! q-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0) ! Ep(q-dq)
        PotVal%d1(:,:,i) = (PotVal%d1(:,:,i)-PotVal_loc0%d0)/(TWO*step)
        Vec%d1(:,:,i)    = (Vec%d1(:,:,i)-Vec_loc0%d0)/(TWO*step)

        NAC%d1(:,:,i)   = matmul(tVec,Vec%d1(:,:,i))

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = (PotVal%d2(:,:,i,i) + PotVal_loc0%d0 - TWO*PotVal%d0)/ &
                                step**2
          Vec%d2(:,:,i,i)    = (Vec%d2(:,:,i,i) + Vec_loc0%d0 - TWO*Vec%d0)/ &
                                step**2
        END IF

        Q_loc(i) = Q(i)

      END DO
    END IF


    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        PotVal%d2(:,:,j,i) = PotVal_loc0%d0
        Vec%d2(:,:,j,i)    = Vec_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) + PotVal_loc0%d0
        Vec%d2(:,:,j,i)    = Vec%d2(:,:,j,i)    + Vec_loc0%d0

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0
        Vec%d2(:,:,j,i)   = Vec%d2(:,:,j,i)     - Vec_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0
        Vec%d2(:,:,j,i)    = Vec%d2(:,:,j,i)    - Vec_loc0%d0

        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i)/(FOUR*step**2)
        PotVal%d2(:,:,i,j) = PotVal%d2(:,:,j,i)
        Vec%d2(:,:,j,i)    = Vec%d2(:,:,j,i)/(FOUR*step**2)
        Vec%d2(:,:,i,j)    = Vec%d2(:,:,j,i)

        Q_loc(i) = Q(i)
        Q_loc(j) = Q(j)
      END DO
      END DO
    END IF


    deallocate(tVec)
    deallocate(Q_loc)
    CALL dealloc_dnMat(PotVal_loc0)
    CALL dealloc_dnMat(Vec_loc0)

  END SUBROUTINE Eval_Pot_Numeric_adia_old
  SUBROUTINE Eval_Pot_Numeric_adia_v3(QModel,Q,PotVal,nderiv,Vec,NAC)
  USE QMLLib_FiniteDiff_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv
    TYPE (dnMat_t),    intent(inout)  :: Vec,NAC

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0,Vec_loc0

    integer                            :: i,j,k,ip,jp,kp

    integer                            :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)

     real (kind=Rkind), allocatable     :: tVec(:,:)

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_adia_v3')

    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                           nderiv=nderiv)
    END IF
    PotVal = ZERO

    IF (Check_NotAlloc_dnMat(Vec,nderiv) ) THEN
      CALL alloc_dnMat(Vec,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    Vec = ZERO

    IF (Check_NotAlloc_dnMat(NAC,nderiv) ) THEN
      CALL alloc_dnMat(NAC,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    NAC = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,nderiv=0)
    CALL alloc_dnMat(Vec_loc0,   nsurf=QModel%QM%nsurf,nVar=QModel%QM%ndim,nderiv=0)
    !write(out_unitp,*) 'coucou1 Eval_Pot_Numeric_adia_v3' ; flush(6)


    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0,vec=Vec_loc0)
    !write(out_unitp,*) 'coucou1.1 Eval_Pot_Numeric_adia_v3' ; flush(6)



    CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=3)
    CALL FiniteDiff_AddMat_TO_dnMat(Vec,   Vec_loc0%d0,   option=3)
    !write(out_unitp,*) 'coucou1.2 Eval_Pot_Numeric_adia_v3' ; flush(6)

    CALL Init_IdMat(NAC%d0,QModel%QM%nsurf)
    !write(out_unitp,*) 'coucou1.3 Eval_Pot_Numeric_adia_v3' ; flush(6)

    allocate(tVec(QModel%QM%nsurf,QModel%QM%nsurf))
    tVec(:,:)      = transpose(Vec%d0)
    !write(out_unitp,*) 'coucou2 0-order Eval_Pot_Numeric_adia_v3' ; flush(6)

    IF (nderiv >= 1) THEN ! 1st derivatives

      DO i=1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                          indQ=[i],indDQ=ind1DQ,option=3)
          CALL FiniteDiff_AddMat_TO_dnMat(Vec,Vec_loc0%d0,              &
                                          indQ=[i],indDQ=ind1DQ,option=3)
        END DO

        NAC%d1(:,:,i)                      = matmul(tVec,Vec%d1(:,:,i))
        IF (nderiv >= 2) NAC%d2(:,:,i,i)   = matmul(tVec,Vec%d2(:,:,i,i))
        IF (nderiv >= 3) NAC%d3(:,:,i,i,i) = matmul(tVec,Vec%d3(:,:,i,i,i))


      END DO
    END IF
    !write(out_unitp,*) 'coucou2 1-order Eval_Pot_Numeric_adia_v3' ; flush(6)

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                        indQ=[i,j],indDQ=ind2DQ,option=3)
          CALL FiniteDiff_AddMat_TO_dnMat(Vec,Vec_loc0%d0,              &
                                        indQ=[i,j],indDQ=ind2DQ,option=3)
        END DO

        NAC%d2(:,:,j,i)    = matmul(tVec,Vec%d2(:,:,j,i))
        IF (nderiv >= 3) THEN
          NAC%d3(:,:,i,i,j)    = matmul(tVec,Vec%d3(:,:,i,i,j))
          NAC%d3(:,:,j,j,i)    = matmul(tVec,Vec%d3(:,:,j,j,i))
        END IF

        CALL FiniteDiff3_SymPerm_OF_dnMat(PotVal,indQ=[i,j])
        CALL FiniteDiff3_SymPerm_OF_dnMat(Vec,indQ=[i,j])
        CALL FiniteDiff3_SymPerm_OF_dnMat(NAC,indQ=[i,j])


      END DO
      END DO
    END IF
    !write(out_unitp,*) 'coucou2 2-order Eval_Pot_Numeric_adia_v3' ; flush(6)

    IF (nderiv >= 3) THEN ! 3d derivatives: d3/dQidQidQj

      ! d3/dQidQjdQk
      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim
      DO k=j+1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                      indQ=[i,j,k],indDQ=ind3DQ,option=3)
          CALL FiniteDiff_AddMat_TO_dnMat(Vec,Vec_loc0%d0,              &
                                      indQ=[i,j,k],indDQ=ind3DQ,option=3)
        END DO

        NAC%d3(:,:,k,j,i)    = matmul(tVec,Vec%d3(:,:,k,j,i))

        CALL FiniteDiff3_SymPerm_OF_dnMat(PotVal,indQ=[i,j,k])
        CALL FiniteDiff3_SymPerm_OF_dnMat(Vec,indQ=[i,j,k])
        CALL FiniteDiff3_SymPerm_OF_dnMat(NAC,indQ=[i,j,k])

      END DO
      END DO
      END DO
    END IF
    !write(out_unitp,*) 'coucou2 3-order Eval_Pot_Numeric_adia_v3' ; flush(6)

    CALL FiniteDiff_Finalize_dnMat(PotVal,step)
    CALL FiniteDiff_Finalize_dnMat(Vec,step)
    CALL FiniteDiff_Finalize_dnMat(NAC,step)
    !write(out_unitp,*) 'coucou3 Eval_Pot_Numeric_adia_v3' ; flush(6)

    deallocate(tVec)
    deallocate(Q_loc)
    CALL dealloc_dnMat(PotVal_loc0)
    CALL dealloc_dnMat(Vec_loc0)
    !write(out_unitp,*) 'coucouf Eval_Pot_Numeric_adia_v3' ; flush(6)

  END SUBROUTINE Eval_Pot_Numeric_adia_v3
  SUBROUTINE dia_TO_adia(PotVal_dia,PotVal_adia,Vec,Vec0,NAC,Phase_Following,   &
                         Phase_checking,nderiv,type_diag)
    USE QMLLib_diago_m
    USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,Write_dnMat,DIAG_dnMat,  &
      Check_NotAlloc_dnMat,get_nsurf,get_nVar
    IMPLICIT NONE

    TYPE (dnMat_t), intent(in)               :: PotVal_dia
    TYPE (dnMat_t), intent(inout)            :: PotVal_adia,Vec,Vec0,NAC

    logical, intent(in)                      :: Phase_Following,Phase_checking

    integer, intent(in), optional            :: nderiv
    integer, intent(in), optional            :: type_diag

    ! local variable
    integer                        :: i,j,k,id,jd,kd,nderiv_loc,ndim,nsurf
    real (kind=Rkind)              :: ai,aj,aii,aij,aji,ajj,th,cc,ss,DEne
    real (kind=Rkind), allocatable :: Eig(:),tVec(:,:),Vdum(:),Vi(:)

    TYPE (dnMat_t)                :: PotVal_dia_onadia



    !test DIAG_dnMat
    TYPE (dnMat_t)              :: dnVec,dnDiag,dnMat
    integer                     :: type_diag_loc = 2    ! tred+tql
    !integer                     :: type_diag_loc = 1 ! jacobi

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='dia_TO_adia'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (present(type_diag))  type_diag_loc = type_diag

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
      !n = size(PotVal_dia%d0,dim=1)
      write(out_unitp,*) 'max Val odd-even',maxval(abs(PotVal_dia%d0(1::2,2::2)))
      write(out_unitp,*) 'max Val even-odd',maxval(abs(PotVal_dia%d0(2::2,1::2)))
      write(out_unitp,*) 'type_diag_loc',type_diag_loc
      write(out_unitp,*) 'Phase_Checking',Phase_Checking
      write(out_unitp,*) 'Phase_Following',Phase_Following
      flush(out_unitp)
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(3,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF

    IF ( Check_NotAlloc_dnMat(PotVal_dia,nderiv_loc) ) THEN
      write(out_unitp,*) ' The diabatic potential MUST be allocated!'
      CALL Write_dnMat(PotVal_dia)
      STOP 'PotVal_dia%dn NOT allocated in "dia_TO_adia"'
    END IF
    IF (debug) THEN
      write(out_unitp,*) 'PotVal_dia'
      CALL Write_dnMat(PotVal_dia,nio=out_unitp)
      flush(out_unitp)
    END IF

    nsurf = get_nsurf(PotVal_dia)
    ndim  = get_nVar(PotVal_dia)


    IF (Check_NotAlloc_dnMat(Vec0,nderiv=0)) THEN
       !$OMP CRITICAL (CRIT_dia_TO_adia)
       CALL alloc_dnMat(Vec0,nsurf=nsurf,nVar=ndim,nderiv=0)

       allocate(Eig(nsurf))

       CALL diagonalization(PotVal_dia%d0,Eig,Vec0%d0,nsurf,sort=1,phase=.TRUE.,type_diag=type_diag_loc)
       !write(out_unitp,*) 'Eig (full diag)',Eig(1:2)
       !CALL diagonalization(PotVal_dia%d0,Eig,Vec0%d0,n=2,sort=1,phase=.TRUE.,type_diag=5)

       deallocate(Eig)

       IF (debug) write(out_unitp,*) 'init Vec0 done'

       !$OMP END CRITICAL (CRIT_dia_TO_adia)
    END IF

    IF (Phase_checking) THEN
      CALL DIAG_dnMat(dnMat=PotVal_dia,dnMatDiag=PotVal_adia,                &
                         dnVec=Vec,dnVecProj=NAC,dnVec0=Vec0,type_diag=type_diag_loc)

      IF (Phase_Following) Vec0%d0 = Vec%d0
    ELSE
      CALL DIAG_dnMat(dnMat=PotVal_dia,dnMatDiag=PotVal_adia,                &
                         dnVec=Vec,dnVecProj=NAC,type_diag=type_diag_loc)
    END IF

    IF (debug) THEN

      write(out_unitp,*) 'Eig',(PotVal_adia%d0(i,i),i=1,nsurf)

      write(out_unitp,*) 'PotVal_adia',PotVal_adia%d0(1:2,1:2)
      CALL Write_dnMat(PotVal_adia,nio=out_unitp)

      write(out_unitp,*) 'Vec'
      CALL Write_dnMat(Vec,nio=out_unitp)

      write(out_unitp,*) 'NAC'
      CALL Write_dnMat(NAC,nio=out_unitp)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE dia_TO_adia

  SUBROUTINE Eval_dnHVib_ana(QModel,Qact,dnH,nderiv)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE ADdnSVM_m, ONLY : dnMat_t,dnS_t,alloc_dnMat,dnMat_TO_dnS,dot_product,       &
    ReduceDerivatives_dnS2_TO_dnS1,Write_dnMat,Write_dnS,ReduceDerivatives_dnS2_TO_dnS1, &
    dnS_TO_dnMat,SYM_dnMat
  USE AdiaChannels_Basis_m
  IMPLICIT NONE

  real (kind=Rkind),              intent(in)    :: Qact(:)
  TYPE (Model_t),                 intent(inout) :: QModel
  TYPE (dnMat_t),                 intent(inout) :: dnH ! derivative of the Hamiltonian
  integer,                        intent(in)    :: nderiv


  integer                        :: ii,ji,i,iq,ib,jb,nb,nq

  TYPE (dnMat_t)                 :: PotVal
  real (kind=Rkind), allocatable :: Q(:),d0GGdef(:,:)
  integer                        :: ndim_act

  TYPE (dnS_t), allocatable      :: dnV(:),dnHB(:)
  TYPE (dnS_t)                   :: dnVfull,dnHij


  ndim_act = size(QModel%QM%list_act)

  nb = QModel%Basis%nb
  nq = QModel%Basis%nq


  CALL alloc_dnMat(dnH,nsurf=nb,nVar=ndim_act,nderiv=nderiv,name_var='dnH')

  allocate(dnV(nq))
  allocate(dnHB(nq))

  allocate(Q(QModel%QM%ndim))
  DO i=1,size(QModel%QM%list_act)
    Q(QModel%QM%list_act(i)) = Qact(i)
  END DO

  !Buid H
  DO iq=1,QModel%Basis%nq
    DO ii=1,size(QModel%QM%list_inact)
      Q(QModel%QM%list_inact(ii)) = QModel%Basis%x(iq)
    END DO

    CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv=nderiv)
    CALL dnMat_TO_dnS(PotVal,dnVfull,i=1,j=1)
    CALL ReduceDerivatives_dnS2_TO_dnS1(dnV(iq),dnVfull,QModel%QM%list_act)
  END DO

  !CALL Write_dnMat(PotVal,6,info='PotVal')
  !CALL Write_dnS(dnV(nq),6,info='dnV',all_type=.TRUE.)

  d0GGdef = QModel%QM%d0GGdef(QModel%QM%list_inact,QModel%QM%list_inact)
  DO ib=1,nb
    ! H B(:,ib)>
    DO iq=1,nq
      dnHB(iq) = -HALF*d0GGdef(1,1)*QModel%Basis%d2gb(iq,ib,1,1) + &
                 dnV(iq)*QModel%Basis%d0gb(iq,ib)
      dnHB(iq) = dnHB(iq) * QModel%Basis%w(iq)
    END DO
    !CALL Write_dnS(dnHB(1),6,info='dnHB',all_type=.TRUE.)
    !write(out_unitp,*) 'coucou dnHB: done',ib ; flush(6)
    DO jb=1,nb
      IF (QModel%Basis%tab_symab(ib) == QModel%Basis%tab_symab(jb)) THEN
        dnHij = dot_product(QModel%Basis%d0gb(:,jb),dnHB(:))
      ELSE
        dnHij = ZERO
      END IF
      CALL dnS_TO_dnMat(dnHij,dnH,jb,ib)
    END DO
  END DO

  dnH = SYM_dnMat(dnH)

  !CALL Write_RMat(dnH%d0,6,5,name_info='H')

  END SUBROUTINE Eval_dnHVib_ana

  SUBROUTINE Eval_Func(QModel,Q,Func,nderiv)
  USE ADdnSVM_m, ONLY : dnS_t,dealloc_dnS,alloc_dnS,Variable,Write_dnS
  IMPLICIT NONE

    TYPE (Model_t),                 intent(inout)            :: QModel

    TYPE (dnS_t),     allocatable,  intent(inout)            :: Func(:)
    real (kind=Rkind),              intent(in)               :: Q(:)
    integer,                        intent(in)               :: nderiv

    ! local variables
    integer                     :: i
    TYPE (dnS_t), allocatable   :: dnQ(:)


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Func'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      write(out_unitp,*) '   nderiv    ',nderiv
      flush(out_unitp)
    END IF

    CALL check_alloc_QM(QModel,name_sub)

    IF (QModel%QM%nb_Func > 0) THEN
      IF (allocated(Func)) THEN
        DO i=1,size(Func)
          CALL dealloc_dnS(Func(i))
        END DO
        deallocate(Func)
      END IF
      allocate(Func(QModel%QM%nb_Func))
      DO i=1,size(Func)
        CALL alloc_dnS(Func(i),QModel%QM%ndimFunc,nderiv)
      END DO

      allocate(dnQ(QModel%QM%ndimFunc))
      DO i=1,QModel%QM%ndimFunc
        dnQ(i) = Variable(Q(i),nVar=QModel%QM%ndimFunc,nderiv=nderiv,iVar=i) ! to set up the derivatives
      END DO

      CALL QModel%QM%Eval_QModel_Func(Func,dnQ,nderiv=nderiv)

      ! deallocation
      DO i=1,size(dnQ)
        CALL dealloc_dnS(dnQ(i))
      END DO
      deallocate(dnQ)
      ! end deallocation

    END IF

    IF (debug) THEN
      write(out_unitp,*) 'Func',size(Func)
      DO i=1,size(Func)
        CALL Write_dnS(Func(i),nio=out_unitp,all_type=.TRUE.)
      END DO
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Eval_Func

  SUBROUTINE Write_Vec(Q,Vec,QModel,nio)
  USE QMLLib_UtilLib_m
  USE AdiaChannels_Basis_m
  USE ADdnSVM_m, ONLY : dnMat_t,Check_NotAlloc_dnMat
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)              :: Q(:)
    TYPE (dnMat_t),     intent(in)              :: Vec
    TYPE (Model_t),     intent(in)              :: QModel
    integer,            intent(in), optional    :: nio

    integer :: i,nio_loc
    real (kind=Rkind), allocatable :: G(:,:)

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    CALL check_alloc_QM(QModel,'Write_Vec')

    IF (.NOT. QModel%QM%print_EigenVec_Basis .AND. .NOT. QModel%QM%print_EigenVec_Grid) RETURN

    IF (Check_NotAlloc_dnMat(Vec,nderiv=0)) THEN
        write(nio_loc,*) '-----------------------------------------------'
        write(nio_loc,*) 'Vec%d0 cannot be printed, it is not allocated'
        write(nio_loc,*) '-----------------------------------------------'
    ELSE

      IF (QModel%QM%print_EigenVec_Basis) THEN
        write(nio_loc,*) '-----------------------------------------------'
        DO i=1,QModel%Basis%nb
          write(nio_loc,*) 'wfb',Q,i,Vec%d0(i,1:QModel%QM%nb_Channels)
        END DO
      END IF

      IF (QModel%QM%print_EigenVec_Grid .AND. QModel%QM%Vib_adia) THEN
        write(nio_loc,*) '-----------------------------------------------'
        allocate(G(QModel%Basis%nq,QModel%QM%nb_Channels))
        DO i=1,QModel%QM%nb_Channels
          CALL BasisTOGrid_Basis(G(:,i),Vec%d0(:,i),QModel%Basis)
        END DO
        DO i=1,QModel%Basis%nq
          write(nio_loc,*) 'wfg',Q,QModel%Basis%x(i),G(i,:)
        END DO

        deallocate(G)

      END IF
      write(nio_loc,*) '-----------------------------------------------'
      flush(nio_loc)
    END IF

  END SUBROUTINE Write_Vec


  SUBROUTINE Write_Model(QModel,nio)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    TYPE(Model_t),      intent(in)              :: QModel
    integer,            intent(in), optional    :: nio

    integer :: nio_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    CALL check_alloc_QM(QModel,'Write_Model')

    IF (nio_loc /= out_unitp) THEN
      open(nio_loc,file=trim(adjustl(QModel%QM%pot_name))//'.out',form='formatted')
    END IF


    write(nio_loc,*) '-----------------------------------------------'
    write(nio_loc,*) 'Output file for potential library'

    CALL QModel%QM%Write_QModel(nio=nio_loc)
    write(nio_loc,*)
    IF (allocated(QModel%QM%d0GGdef)) CALL Write_RMat(QModel%QM%d0GGdef,&
                               nio_loc,5,name_info='d0GGdef')
    write(nio_loc,*)
    IF (allocated(QModel%QM%Q0)) CALL Write_RVec(QModel%QM%Q0,          &
                               nio_loc,5,name_info='Q0')
    write(nio_loc,*)
    write(nio_loc,*) '-----------------------------------------------'
    write(nio_loc,*) 'Extra action(s):'
    write(nio_loc,*) 'opt',QModel%opt
    write(nio_loc,*) 'irc',QModel%irc
    write(nio_loc,*) '-----------------------------------------------'
    flush(nio_loc)


  END SUBROUTINE Write_Model
  SUBROUTINE Write0_Model(QModel,nio)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    TYPE(Model_t),    intent(in)              :: QModel
    integer,          intent(in), optional    :: nio

    integer :: nio_loc


    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    CALL check_alloc_QM(QModel,'Write0_Model')

    IF (nio_loc /= out_unitp) THEN
      open(nio_loc,file=trim(adjustl(QModel%QM%pot_name))//'.out',form='formatted')
    END IF


    CALL QModel%QM%Write0_QModel(nio=nio_loc)


     IF (nio_loc /= out_unitp) THEN
      close(nio_loc)
    END IF


  END SUBROUTINE Write0_Model
  SUBROUTINE Write_QdnV_FOR_Model(Q,PotVal,QModel,Vec,NAC,info,name_file)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    TYPE (Model_t),    intent(in)           :: QModel
    TYPE (dnMat_t),    intent(in)           :: PotVal
    real (kind=Rkind), intent(in)           :: Q(:)
    TYPE (dnMat_t),    intent(in), optional :: Vec ! for non adiabatic couplings
    TYPE (dnMat_t),    intent(in), optional :: NAC ! for non adiabatic couplings
    character(len=*),  intent(in), optional :: info
    character(len=*),  intent(in), optional :: name_file

    integer :: nio_loc,err_io

    IF (present(name_file)) THEN
      CALL file_open2(trim(adjustl(name_file)),                                 &
                      nio_loc,lformatted=.TRUE.,append=.TRUE.,err_file=err_io)
    ELSE
      CALL file_open2(trim(adjustl(QModel%QM%pot_name))//'.txt',                &
                      nio_loc,lformatted=.TRUE.,append=.TRUE.,err_file=err_io)
    END IF

    IF (err_io /= 0) THEN
      write(out_unitp,*) 'ERROR in Write_QdnV_FOR_Model'
      write(out_unitp,*) ' Impossible to open the file "',                      &
                          trim(adjustl(QModel%QM%pot_name))//'.txt','"'
      STOP 'Impossible to open the file'
    END IF


    write(nio_loc,'(a)',advance='no') 'TEST output: '
    IF (present(info)) THEN
      write(nio_loc,'(a)',advance='no') info
    END IF
    IF (QModel%QM%adiabatic) THEN
      write(nio_loc,'(a)') ' Adiabatic'
    ELSE
      write(nio_loc,'(a)') ' Diabatic'
    END IF

    write(nio_loc,*) 'Q'
    write(nio_loc,*) size(Q)
    write(nio_loc,*) Q

    IF (allocated(PotVal%d0)) THEN
      write(nio_loc,*) 'V'
      write(nio_loc,*) size(PotVal%d0)
      write(nio_loc,*) PotVal%d0
    END IF
    IF (allocated(PotVal%d1)) THEN
      write(nio_loc,*) 'Grad'
      write(nio_loc,*) size(PotVal%d1)
      write(nio_loc,*) PotVal%d1
    END IF
    IF (allocated(PotVal%d2)) THEN
      write(nio_loc,*) 'Hess'
      write(nio_loc,*) size(PotVal%d2)
      write(nio_loc,*) PotVal%d2
    END IF

    IF (present(Vec)) THEN
    IF (allocated(Vec%d0)) THEN
      write(nio_loc,*) 'Vec'
      write(nio_loc,*) size(Vec%d0)
      write(nio_loc,*) Vec%d0
    END IF
    IF (allocated(Vec%d1)) THEN
      write(nio_loc,*) 'd1Vec'
      write(nio_loc,*) size(Vec%d1)
      write(nio_loc,*) Vec%d1
    END IF
    IF (allocated(Vec%d2)) THEN
      write(nio_loc,*) 'd2Vec'
      write(nio_loc,*) size(Vec%d2)
      write(nio_loc,*) Vec%d2
    END IF
    END IF

    IF (present(NAC)) THEN
    IF (allocated(NAC%d1)) THEN
      write(nio_loc,*) 'NAC'
      write(nio_loc,*) size(NAC%d1)
      write(nio_loc,*) NAC%d1
    END IF
    END IF

    IF (allocated(QModel%QM%d0GGdef)) THEN
      write(nio_loc,*) 'd0GGdef'
      write(nio_loc,*) size(QModel%QM%d0GGdef)
      write(nio_loc,*) QModel%QM%d0GGdef
    END IF

    write(nio_loc,'(a)',advance='no') 'END_TEST output: '
    IF (present(info)) THEN
      write(nio_loc,'(a)',advance='no') info
    END IF
    IF (QModel%QM%adiabatic) THEN
      write(nio_loc,'(a)') ' Adiabatic'
    ELSE
      write(nio_loc,'(a)') ' Diabatic'
    END IF

    close(nio_loc)

  END SUBROUTINE Write_QdnV_FOR_Model
  SUBROUTINE Check_analytical_numerical_derivatives(QModel,Q,nderiv)
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Write_dnMat, &
    get_maxval_OF_dnMat, operator(-)
  IMPLICIT NONE

    TYPE (Model_t),       intent(inout)   :: QModel
    real (kind=Rkind),    intent(in)      :: Q(:)
    integer,              intent(in)      :: nderiv

    ! local variables
    TYPE (dnMat_t)            :: Mat_diff
    TYPE (dnMat_t)            :: PotVal_ana,PotVal_num
    TYPE (dnMat_t)            :: NAC_ana,NAC_num
    TYPE (dnMat_t)            :: Vec_ana,Vec_num

    real (kind=Rkind)         :: MaxMat,MaxDiffMat

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Check_analytical_numerical_derivatives'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (QModel%QM%no_ana_der) RETURN

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      write(out_unitp,*) '   nderiv    ',nderiv
      flush(out_unitp)
    END IF

    CALL check_alloc_QM(QModel,name_sub)


    CALL alloc_dnMat(PotVal_ana,nsurf=QModel%QM%nsurf,              &
                         nVar=QModel%QM%ndim,nderiv=nderiv)

    CALL alloc_dnMat(PotVal_num,nsurf=QModel%QM%nsurf,              &
                         nVar=QModel%QM%ndim,nderiv=nderiv)


    IF (QModel%QM%adiabatic .AND. QModel%QM%nsurf > 1) THEN
      CALL Eval_Pot(QModel,Q,PotVal_ana,nderiv,NAC_ana,Vec_ana,numeric=.FALSE.)
    ELSE
      CALL Eval_Pot(QModel,Q,PotVal_ana,nderiv,numeric=.FALSE.)
    END IF

    IF (debug) THEN
      write(out_unitp,*)   'PotVal_ana'
      CALL Write_dnMat(PotVal_ana,nio=out_unitp)
      flush(out_unitp)
    END IF

    IF (QModel%QM%adiabatic .AND. QModel%QM%nsurf > 1) THEN
      CALL Eval_Pot(QModel,Q,PotVal_num,nderiv,NAC_num,Vec_num,numeric=.TRUE.)
    ELSE
      CALL Eval_Pot(QModel,Q,PotVal_num,nderiv,numeric=.TRUE.)
    END IF
    IF (debug) THEN
      write(out_unitp,*)   'PotVal_num'
      CALL Write_dnMat(PotVal_num,nio=out_unitp)
      flush(out_unitp)
    END IF


    MaxMat      = get_maxval_OF_dnMat(PotVal_ana)
    IF (MaxMat < ONETENTH**6) MaxMat = ONE
    Mat_diff    = PotVal_num - PotVal_ana
    MaxDiffMat  = get_maxval_OF_dnMat(Mat_diff)

    write(out_unitp,'(3a,e9.2)') 'With ',QModel%QM%pot_name,                    &
               ': max of the relative Potential diff:',MaxDiffMat/MaxMat
    write(out_unitp,'(3a,l9)')   'With ',QModel%QM%pot_name,                    &
     ': Potential diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

    IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
      write(out_unitp,*)   'Potential diff (ana-numer)'
      CALL Write_dnMat(Mat_diff,nio=out_unitp)
    END IF

    IF (QModel%QM%adiabatic .AND. QModel%QM%nsurf > 1) THEN

      MaxMat      = get_maxval_OF_dnMat(NAC_ana)
      IF (MaxMat < ONETENTH**6) MaxMat = ONE
      Mat_diff    = NAC_num - NAC_ana
      MaxDiffMat  = get_maxval_OF_dnMat(Mat_diff)

      write(out_unitp,'(3a,e9.2)') 'With ',QModel%QM%pot_name,                  &
                 ': max of the relative NAC diff:',MaxDiffMat/MaxMat
      write(out_unitp,'(3a,l9)')   'With ',QModel%QM%pot_name,                  &
       ': NAC diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

      IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
        write(out_unitp,*)   'NAC diff (ana-numer)'
        CALL Write_dnMat(Mat_diff,nio=out_unitp)
        write(out_unitp,*)   'NAC_ana'
        CALL Write_dnMat(NAC_ana,nio=out_unitp)
        write(out_unitp,*)   'NAC_num'
        CALL Write_dnMat(NAC_num,nio=out_unitp)
      END IF

      MaxMat      = get_maxval_OF_dnMat(Vec_ana)
      IF (MaxMat < ONETENTH**6) MaxMat = ONE
      Mat_diff    = Vec_num - Vec_ana
      MaxDiffMat  = get_maxval_OF_dnMat(Mat_diff)

      write(out_unitp,'(3a,e9.2)') 'With ',QModel%QM%pot_name,            &
                 ': max of the relative Vec diff:',MaxDiffMat/MaxMat
      write(out_unitp,'(3a,l9)')   'With ',QModel%QM%pot_name,            &
       ': Vec diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

      IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
        write(out_unitp,*)   'Vec diff (ana-numer)'
        CALL Write_dnMat(Mat_diff,nio=out_unitp)
        write(out_unitp,*)   'Vec_ana'
        CALL Write_dnMat(Vec_ana,nio=out_unitp)
        write(out_unitp,*)   'Vec_num'
        CALL Write_dnMat(Vec_num,nio=out_unitp)
      END IF

    END IF

    CALL dealloc_dnMat(PotVal_ana)
    CALL dealloc_dnMat(PotVal_num)
    CALL dealloc_dnMat(NAC_ana)
    CALL dealloc_dnMat(NAC_num)
    CALL dealloc_dnMat(Vec_ana)
    CALL dealloc_dnMat(Vec_num)
    CALL dealloc_dnMat(Mat_diff)

    IF (debug) THEN
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Check_analytical_numerical_derivatives

  SUBROUTINE Eval_pot_ON_Grid(QModel,Qmin,Qmax,nb_points,nderiv,grid_file)
  USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),               intent(inout)   :: QModel
    real (kind=Rkind),            intent(in)      :: Qmin(:),Qmax(:)
    integer, optional,            intent(in)      :: nb_points,nderiv
    character (len=*), optional,  intent(in)      :: grid_file


    ! local variables
    integer                        :: unit_grid_file
    integer                        :: i,iq,jq,i1,i2,nb_points_loc,nderiv_loc,ndim_loc
    integer,           allocatable :: i_Q(:)
    real (kind=Rkind), allocatable :: dQ(:),Q(:)
    TYPE (dnMat_t)                 :: PotVal,NAC


    CALL check_alloc_QM(QModel,'Eval_pot_ON_Grid')


    IF (size(Qmin) /= QModel%QM%ndim .OR. size(Qmax) /= QModel%QM%ndim) THEN
       write(out_unitp,*) ' ERROR in Eval_pot_ON_Grid'
       write(out_unitp,*) ' The size of Qmin or Qmax are different from QModel%QM%ndim'
       write(out_unitp,*) '  size(Qmin)    ',size(Qmin)
       write(out_unitp,*) '  size(Qmax)    ',size(Qmax)
       write(out_unitp,*) '  QModel%QM%ndim',QModel%QM%ndim
       write(out_unitp,*) ' => Check the fortran'
       STOP 'ERROR in Eval_pot_ON_Grid: problem with QModel%QM%ndim'
    END IF

    IF (present(grid_file)) THEN
      IF (len_trim(grid_file) == 0) THEN
        unit_grid_file = out_unitp
      ELSE
        unit_grid_file = 99
        open(unit=unit_grid_file,file=trim(grid_file) )
      END IF
    ELSE
      unit_grid_file = out_unitp
    END IF

    nb_points_loc = 100
    IF (present(nb_points)) nb_points_loc = nb_points
    nb_points_loc = max(nb_points_loc,2)

    IF (present(nderiv)) THEN
      nderiv_loc = nderiv
    ELSE
      nderiv_loc = 0
    END IF

    allocate(dQ(QModel%QM%ndim))
    allocate(Q(QModel%QM%ndim))
    allocate(i_Q(QModel%QM%ndim))

    dQ(:)       = (Qmax-Qmin) / real(nb_points_loc-1,kind=Rkind)
    ndim_loc    = 0
    i_Q(:)      = 0
    DO i=1,QModel%QM%ndim
      IF (dQ(i) /= ZERO) THEN
        ndim_loc = ndim_loc + 1
        i_Q(ndim_loc) = i
      END IF
    END DO
    write(out_unitp,*) 'Grid. File name: "',trim(grid_file),'"'
    !write(out_unitp,*) 'Coordinates indices, i_Q: ',i_Q(1:ndim_loc)
    !write(out_unitp,*) 'QModel%QM%ndim',QModel%QM%ndim
    !write(out_unitp,*) 'QModel%QM%numeric',QModel%QM%numeric
    !write(out_unitp,*) 'ndim for the grid',ndim_loc


    Q(:) = Qmin

    IF (ndim_loc == 1) THEN
      i1 = i_Q(1)
      DO iq=0,nb_points_loc-1
        Q(i1) = Qmin(i1) + dQ(i1)*real(iq,kind=Rkind)
        IF (QModel%QM%nsurf > 1 .AND. QModel%QM%adiabatic) THEN
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=max(1,nderiv_loc),NAC=NAC)

          IF (nderiv_loc == 0) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%QM%nsurf),NAC%d1
          ELSE IF (nderiv_loc == 1) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%QM%nsurf),(PotVal%d1(i,i,:),i=1,QModel%QM%nsurf)
          ELSE
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%QM%nsurf),(PotVal%d1(i,i,:),i=1,QModel%QM%nsurf), &
                            (PotVal%d2(i,i,:,:),i=1,QModel%QM%nsurf)
          END IF
          flush(unit_grid_file)

        ELSE
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv_loc)

          IF (nderiv_loc == 0) THEN
            write(unit_grid_file,*) Q(i1),PotVal%d0
          ELSE IF (nderiv_loc == 1) THEN
            write(unit_grid_file,*) Q(i1),PotVal%d0,PotVal%d1
          ELSE
            write(unit_grid_file,*) Q(i1),PotVal%d0,PotVal%d1,PotVal%d2
          END IF
          flush(unit_grid_file)

        END IF
      END DO
    ELSE IF (ndim_loc == 2) THEN
      i1 = i_Q(1)
      i2 = i_Q(2)
      DO iq=0,nb_points_loc-1
      DO jq=0,nb_points_loc-1
        Q(i1) = Qmin(i1) + dQ(i1)*real(iq,kind=Rkind)
        Q(i2) = Qmin(i2) + dQ(i2)*real(jq,kind=Rkind)

        IF (QModel%QM%nsurf > 1 .AND. QModel%QM%adiabatic) THEN
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=0,NAC=NAC)
          write(unit_grid_file,*) Q(i1),Q(i2),(PotVal%d0(i,i),i=1,QModel%QM%nsurf)
        ELSE
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv_loc)
          write(unit_grid_file,*) Q(i1),Q(i2),PotVal%d0
        END IF
        flush(unit_grid_file)

      END DO
      write(unit_grid_file,*)
      END DO


    END IF

    CALL dealloc_dnMat(PotVal)
    CALL dealloc_dnMat(NAC)
    deallocate(dQ)
    deallocate(Q)
    deallocate(i_Q)

    IF (unit_grid_file /= out_unitp) THEN
      close(unit_grid_file)
    END IF


  END SUBROUTINE Eval_pot_ON_Grid


  SUBROUTINE calc_pot(V,QModel,Q)
  USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),         intent(inout)   :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated

    TYPE (dnMat_t)                         :: PotVal

    CALL check_alloc_QM(QModel,'calc_pot')

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=0)

    V = PotVal%d0

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot
  SUBROUTINE calc_pot_grad(V,g,QModel,Q)
  USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),         intent(inout)   :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(QModel,'calc_pot_grad')

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=1)

    V = PotVal%d0
    g = PotVal%d1

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot_grad
  SUBROUTINE calc_grad(g,QModel,Q)
  USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),         intent(inout)   :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(QModel,'calc_grad')

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=1)

    g = PotVal%d1

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_grad
  SUBROUTINE calc_pot_grad_hess(V,g,h,QModel,Q)
  USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),         intent(inout)   :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(QModel,'calc_pot_grad_hess')

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)

    V = PotVal%d0
    g = PotVal%d1
    h = PotVal%d2

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot_grad_hess
  SUBROUTINE calc_hess(h,QModel,Q)
  USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
  IMPLICIT NONE

    TYPE (Model_t),     intent(inout)     :: QModel
    real (kind=Rkind),  intent(in)        :: Q(:)
    real (kind=Rkind),  intent(inout)     :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(QModel,'calc_hess')

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)

    h = PotVal%d2

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_hess
END MODULE Model_m
