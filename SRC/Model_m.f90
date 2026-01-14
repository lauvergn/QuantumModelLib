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
MODULE Model_m
  !$ USE omp_lib
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  USE AdiaChannels_Basis_m, ONLY : QML_Basis_t
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Model_t,Init_Model,Eval_Pot,Eval_Func,Eval_tab_HMatVibAdia,Eval_dnHVib_ana
  PUBLIC :: check_alloc_QM,check_Init_QModel,dealloc_Model
  PUBLIC :: check_alloc_d0GGdef
  PUBLIC :: Write_Model
  PUBLIC :: calc_pot,calc_grad,calc_hess,calc_pot_grad,calc_pot_grad_hess
  PUBLIC :: Check_analytical_numerical_derivatives
  PUBLIC :: Eval_pot_ON_Grid,get_Q0_Model,get_d0GGdef_Model,Qact_TO_Q
  PUBLIC :: Set_step_epsi_Model
  PUBLIC :: Write_QdnV_FOR_Model,Test_QdnV_FOR_Model,Test_QVG_FOR_Model

  TYPE :: QML_t
    CLASS (QML_Empty_t),  allocatable :: QM
  END TYPE QML_t
  TYPE :: Model_t
    ! Add nsurf and ndim to avoid crash when using the driver without initialization
    ! At the intialization, the variables are set-up to the correct values and are
    !   identical to QM%nsurf and QM%ndim ones respectively.
    integer                           :: nsurf       = 0
    integer                           :: ndim        = 0
    integer                           :: ipot        = -1
    integer                           :: icap        = -1
    integer                           :: idipx       = -1
    integer                           :: idipy       = -1
    integer                           :: idipz       = -1
    TYPE (QML_t),         allocatable :: tab_Op(:)
    CLASS (QML_Empty_t),  allocatable :: QM
    TYPE (QML_Basis_t),   allocatable :: Basis ! Basis for the adiabatic separation between coordinates
    logical                           :: opt         = .FALSE.
    logical                           :: irc         = .FALSE.
  END TYPE Model_t



  !real (kind=Rkind)                     :: step = ONETENTH**4 ! model TWOD => 0.4e-7 (nderiv=2)
  real (kind=Rkind)                     :: step = ONETENTH**3 ! model TWOD => 0.6e-9 (nderiv=2)
  !real (kind=Rkind)                     :: step = ONETENTH**2 ! model TWOD => 0.1e-7 (nderiv=2)

  real (kind=Rkind)                     :: epsi = ONETENTH**10

  character (len=*), parameter :: QML_version =                         &
#if defined(__QML_VER)
      __QML_VER
#else
      'unknown: -D__QML_VER=?'
#endif

  character (len=*), parameter :: compile_date =                          &
#if defined(__COMPILE_DATE)
      __COMPILE_DATE
#else
      'unknown: -D__COMPILE_DATE=?'
#endif

  character (len=*), parameter :: compile_host =                          &
#if defined(__COMPILE_HOST)
      __COMPILE_HOST
#else
      "unknown: -D__COMPILE_HOST=?"
#endif

  logical, private :: Print_Version_done = .FALSE.

  TYPE(Model_t), PUBLIC  :: QuantumModel

CONTAINS

  SUBROUTINE version_QML(Print_Version)
    USE iso_fortran_env
    USE QDUtil_m
    USE QMLLib_UtilLib_m, ONLY : QML_path,check_QML_Path
    IMPLICIT NONE

    logical,             intent(in)    :: Print_Version

    IF (Print_Version) THEN
      CALL check_QML_Path()
      Print_Version_done = .TRUE.
      write(out_unit,*) '================================================='
      write(out_unit,*) '================================================='
      write(out_unit,*) '== QML: Quantum Model Lib (E-CAM) ==============='
      write(out_unit,*) '== QML version:       ',QML_version
      write(out_unit,*) '== QML path:          ',QML_path
      write(out_unit,*) '-------------------------------------------------'
      write(out_unit,*) '== Compiled on       "',compile_host, '" the ',compile_date
      write(out_unit,*) '== Compiler:         ',compiler_version()
      write(out_unit,*) '== Compiler options: ',compiler_options()
      write(out_unit,*) '-------------------------------------------------'
      write(out_unit,*) 'QuantumModelLib* is a free software under the MIT Licence.'
      write(out_unit,*) '  Copyright (c) 2022 David Lauvergnat [1]'
      write(out_unit,*) '  with contributions of:'
      write(out_unit,*) '    Félix MOUHAT [2]'
      write(out_unit,*) '    Liang LIANG [3]'
      write(out_unit,*) '    Emanuele MARSILI [1,4]'
      write(out_unit,*) '    Evaristo Villaseco Arribas [5]'
      write(out_unit,*) 
      write(out_unit,*) '  [1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
      write(out_unit,*) '  [2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France'
      write(out_unit,*) '  [3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France'
      write(out_unit,*) '  [4]: Durham University, Durham, UK'
      write(out_unit,*) '  [5]: Department of Physics, Rutgers University, Newark, New Jersey 07102, USA'
      write(out_unit,*) '    * Originally, it has been developed during the Quantum-Dynamics E-CAM project :'
      write(out_unit,*) '       https://www.e-cam2020.eu/quantum-dynamics'
      write(out_unit,*) '================================================='
      write(out_unit,*) '============ WARNING ============================'
      write(out_unit,*) ' From the version 23.3:'
      write(out_unit,*) '   the "Phase_Following" and "Phase_Checking" default are .FALSE.'
      write(out_unit,*) '================================================='
    END IF
    
  END SUBROUTINE version_QML

  SUBROUTINE Read_Model(QModel_inout,nio,read_nml1,opt1,IRC1)
    IMPLICIT NONE

    TYPE (QML_Empty_t),  intent(inout) :: QModel_inout ! variable to transfer info to the init
    integer,             intent(in)    :: nio
    logical,             intent(inout) :: read_nml1,opt1,IRC1

    ! local variable
    integer, parameter :: max_act = 10
    integer, parameter :: max_Op = 3

    integer :: ndim,nsurf,nderiv,option,printlevel,nb_Channels
    logical :: adiabatic,numeric,PubliUnit,read_nml
    logical :: Vib_adia,print_EigenVec_Grid,print_EigenVec_Basis
    logical :: opt,IRC
    logical :: Phase_checking,Phase_Following
    logical :: Cart_TO_Q,AbInitio,MassWeighted

    character (len=100) :: pot_name
    integer :: err_read,nb_act
    integer :: list_act(max_act)
    integer :: list_Op(max_Op)

    ! Namelists for input file
    namelist /potential/ ndim,nsurf,pot_name,numeric,adiabatic,option,PubliUnit,&
                         Phase_Checking,Phase_Following,                        &
                         Cart_TO_Q,MassWeighted,AbInitio,                       &
                         list_Op,                                               &
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
    MassWeighted    = QModel_inout%MassWeighted
    AbInitio        = .FALSE.
    list_Op(:)      = -1 ! 0: potential, then other scalar operators

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

    write(out_unit,*) 'Reading input file . . .'
    read(nio,nml=potential,IOSTAT=err_read)

    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_Model'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "potential" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_Model: End-of-file or End-of-record while reading the namelist'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_Model'
      write(out_unit,*) ' Some parameter names of the namelist "potential" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=potential)
      STOP ' ERROR in Read_Model: wrong parameter(s) in the namelist'
    END IF

    !write(out_unit,nml=potential)

    read_nml1                         = read_nml
    opt1                              = opt
    IRC1                              = IRC

    QModel_inout%option               = option
    CALL set_print_level(printlevel) ! from the module QDUtil lib
    QModel_inout%ndim                 = ndim
    QModel_inout%nsurf                = nsurf
    QModel_inout%adiabatic            = adiabatic
    QModel_inout%numeric              = numeric
    QModel_inout%Phase_Checking       = (adiabatic .AND. Phase_Checking)
    QModel_inout%Phase_Following      = (adiabatic .AND. Phase_Following)
    QModel_inout%Cart_TO_Q            = Cart_TO_Q
    QModel_inout%AbInitio             = AbInitio
    QModel_inout%MassWeighted         = (Cart_TO_Q .AND. MassWeighted)

    QModel_inout%pot_name             = trim(pot_name)
    QModel_inout%PubliUnit            = PubliUnit

    QModel_inout%print_EigenVec_Grid  = print_EigenVec_Grid
    QModel_inout%print_EigenVec_Basis = print_EigenVec_Basis


    IF (Vib_adia) THEN
      QModel_inout%nb_Channels    = nb_Channels
      QModel_inout%Vib_adia       = Vib_adia

      IF (nb_Channels == 0) THEN
        write(out_unit,*) ' ERROR in Read_Model'
        write(out_unit,*) ' Vib_adia=t and nb_Channels = 0'
        write(out_unit,*) ' You have to define "nb_Channels" in the namelist.'
        write(out_unit,*) ' check your data!'
        write(out_unit,*)
        STOP ' ERROR in Read_Model: define "nb_Channels" in the namelist'
      END IF

      nb_act = count(list_act /= 0)
      IF (nb_act == 0) THEN
        write(out_unit,*) ' ERROR in Read_Model'
        write(out_unit,*) ' Vib_adia=t and nb_act = 0'
        write(out_unit,*) ' You have to define "list_act(:)" in the namelist.'
        write(out_unit,*) ' check your data!'
        write(out_unit,*)
        STOP ' ERROR in Read_Model: define "list_act(:)" in the namelist.'
      END IF
      QModel_inout%list_act = list_act(1:nb_act)

      IF (count(QModel_inout%list_act == 0) /= 0) THEN
        write(out_unit,*) ' ERROR in Read_Model'
        write(out_unit,*) ' list_act(:) is wrong.'
        write(out_unit,*) ' list_act(1:nb_act) has some 0 :',list_act(1:nb_act)
        write(out_unit,*) ' You have to define in the namelist list_act(:)'
        write(out_unit,*) ' check your data!'
        write(out_unit,*)
        STOP ' ERROR in Read_Model: "list_act(:)" with 0 in the namelist.'
      END IF
    END IF

  END SUBROUTINE Read_Model

  SUBROUTINE Init_Model(QModel,pot_name,ndim,nsurf,adiabatic,Cart_TO_Q,         &
                        read_param,param_file_name,nio_param_file,              &
                        option,PubliUnit,Print_init,Vib_adia,                   &
                        Phase_Following,Phase_checking)

  USE QDUtil_m,         ONLY : TO_lowercase
  USE QMLLib_UtilLib_m

  USE QML_Empty_m

  USE QML_Template_m
  USE QML_Test_m
  USE QML_ExtModel_m

  USE QML_Morse_m
  USE QML_Poly1D_m
  USE QML_H2_m

  USE QML_HenonHeiles_m
  USE QML_DoubleWell_m
  USE QML_Tully_m

  USE QML_PSB3_m
  USE QML_Retinal_JPCB2000_m

  USE QML_HONO_m
  USE QML_HNNHp_m
  USE QML_H2SiN_m
  USE QML_H2NSi_m
  USE QML_CHFClBr_m

  USE QML_H2O_m
  USE QML_H2_H2On_m

  USE QML_ClH2p_m
  USE QML_ClH2p_Botschwina_m
  USE QML_Bottleneck_m
  USE QML_HNO3_m
  USE QML_NO3_m
  USE QML_CH5_m
  USE QML_PH4Jo_m
  USE QML_PH4_m

  USE QML_HOO_DMBE_m
  USE QML_H3_m
  USE QML_CNH_Murrell_m

  USE QML_OneDSOC_1S1T_m
  USE QML_OneDSOC_2S1T_m

  USE QML_LinearHBond_m
  USE QML_TwoD_MullerBrown_m
  USE QML_Buck_m
  USE QML_Phenol_m
  USE QML_Sigmoid_m

  USE QML_TwoD_m
  USE QML_TwoD_RJDI2014_m
  USE QML_TwoD_Valahu2022_m
  USE QML_Vibronic_m
  USE QML_Uracil_m
  USE QML_fulvene_m
  USE QML_dmabn_m

  USE AdiaChannels_Basis_m

  USE QML_OneD_Photons_m
  USE QML_OneD_Photons2_m
  IMPLICIT NONE

    TYPE (Model_t),      intent(inout)           :: QModel

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
    TYPE(QML_Empty_t)              :: QModel_in ! to transfer info into QModel
    integer                        :: i,nio_loc,i_inact,nb_inact
    logical                        :: read_param_loc,read_nml,Print_init_loc
    logical,           allocatable :: list_Q(:)
    character (len=:), allocatable :: param_file_name_loc,pot_name_loc,tab_pot_name(:)
    real (kind=Rkind), allocatable :: Q0(:)

    Print_init_loc = .TRUE.
    IF (present(Print_init)) Print_init_loc = Print_init

    CALL version_QML(Print_init_loc)

    IF (Print_init_loc) THEN
      write(out_unit,*) '================================================='
      write(out_unit,*) '== Initialization of the Model =================='
    END IF

    ! test the QML_path variable (it enables to test is the QML directory has been moved)
    CALL check_QML_Path()

    CALL dealloc_Model(QModel)

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

    IF (present(Vib_adia)) THEN
      QModel_in%Vib_adia = Vib_adia
    ELSE
      QModel_in%Vib_adia = .FALSE.
    END IF

    IF (present(Phase_Following)) THEN
      QModel_in%Phase_Following = Phase_Following
    ELSE
      QModel_in%Phase_Following = .FALSE.
    END IF

    IF (present(Phase_checking)) THEN
      QModel_in%Phase_checking = Phase_checking
    ELSE
      QModel_in%Phase_checking = .FALSE.
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
          param_file_name_loc = trim("input.dat")
          nio_loc = 99 ! this value is not used
        ELSE
          param_file_name_loc = trim(param_file_name)
          nio_loc = 99 ! this value is not used
        END IF
      ELSE
      nio_loc = in_unit
      END IF
    END IF

    IF (present(pot_name)) THEN
      pot_name_loc   = TO_lowercase(trim(pot_name))
      read_param_loc = (read_param_loc .OR.  pot_name_loc == 'read_model')
    ELSE
      IF (.NOT. read_param_loc) THEN
        write(out_unit,*) 'ERROR in Init_Model'
        write(out_unit,*) ' pot_name is not present and read_param=F'
        STOP 'ERROR in Init_Model: pot_name is not present and read_param=F'
      END IF
    END IF

    read_nml       = .FALSE.
    IF (read_param_loc) THEN
      IF (nio_loc /= in_unit) THEN
        open(newunit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
      END IF
      CALL Read_Model(QModel_in,nio_loc,read_nml,QModel%opt,QModel%IRC)
      pot_name_loc = QModel_in%pot_name

      IF (allocated(param_file_name_loc)) deallocate(param_file_name_loc)
    END IF

    IF (QModel%opt) write(out_unit,*) ' Geometry optimization will be performed'
    IF (QModel%IRC) write(out_unit,*) ' IRC will be performed'


    QModel_in%Init = .TRUE.

    IF (QModel_in%adiabatic) THEN
      IF (Print_init_loc) THEN
        write(out_unit,*) 'Adiabatic potential . . .'
        write(out_unit,*) 'Phase_Checking',QModel_in%Phase_Checking
        write(out_unit,*) 'Phase_Following',QModel_in%Phase_Following
      END IF
    ELSE
      IF (Print_init_loc) write(out_unit,*) 'Non-adiabatic potential . . .'
    END IF

    IF (QModel_in%numeric .AND. Print_init_loc) THEN
      write(out_unit,*) 'You have decided to perform a numeric checking of the analytic formulas.'
    END IF

    pot_name_loc = TO_lowercase(pot_name_loc)
    CALL Pot_Name_Analysis(pot_name_loc,tab_pot_name)
    IF (size(tab_pot_name) < 1) STOP 'ERROR in Pot_Name_Analysis'
    IF (Print_init_loc) THEN
      write(out_unit,*) 'pot_name_loc: ',pot_name_loc
      IF (allocated(tab_pot_name)) THEN
        DO i=1,size(tab_pot_name)
          write(out_unit,*) 'tab_pot_name(i): ',i,tab_pot_name(i)
        END DO
      END IF
    END IF

    SELECT CASE (trim(tab_pot_name(1)))
    CASE ('morse')
      allocate(QML_Morse_t :: QModel%QM)
      QModel%QM = Init_QML_Morse(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('poly1d')
      allocate(QML_Poly1D_t :: QModel%QM)
      QModel%QM = Init_QML_Poly1D(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('h2')
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
      !! remark: three options are possible (option = 1,2,3)
      !!    option 1: usual HenonHeiles (default)
      !!    option 2: quadratic contribution + morse potentials and tanh contributions
      !!    option 3: quadratic contribution + tanh contributions for the anharmonic part
      !! reduced masses      = ONE au
      !! ref:  parameters taken from M. Nest, H.-D. Meyer, J. Chem. Phys. 117 (2002) 10499. doi:10.1063/1.1521129
      !! === END README ==
      allocate(QML_HenonHeiles_t :: QModel%QM)
      QModel%QM = Init_QML_HenonHeiles(QModel_in,                               &
                                       read_param=read_nml,nio_param_file=nio_loc)

    CASE ('doublewell')
      !! === README ==
      !! DoubleWell potential
      !! pot_name  = 'DoubleWell'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced masses      = ONE au
      !! ref:  Nancy Makri & William H. Miller, The Journal of Chemical Physics 87, 5781 (1987); doi: 10.1063/1.453501
      !! === END README ==
      allocate(QML_DoubleWell_t :: QModel%QM)
      QModel%QM = Init_QML_DoubleWell(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

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
      !! Diagonal Metric Tensor      = [ 0.0005786177, 0.0002550307 ] au
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
      !! Reduced masses$      = [ 20000., 6667. ] au
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
      !! Reduced masses      = [1. , 1.] au
      !! ref:  Ilya G. Ryabinkin, Loïc Joubert-Doriol, and Artur F. Izmaylov, ...
      !!       ... J. Chem. Phys. 140, 214116 (2014); https://doi.org/10.1063/1.4881147
      !! === END README ==
      allocate(QML_TwoD_RJDI2014_t :: QModel%QM)
      QModel%QM = Init_QML_TwoD_RJDI2014(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('twod_valahu2022')
      !! === README ==
      !! 2D model
      !! pot_name  = 'TwoD_Valahu2022'
      !! ndim      = 2 (X,Y)
      !! nsurf     = 2
      !! ref:  xxxxx
      !! === END README ==
      allocate(QML_TwoD_Valahu2022_t :: QModel%QM)
      QModel%QM = Init_QML_TwoD_Valahu2022(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

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
      QModel%QM = Init_QML_Retinal_JPCB2000(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('vibronic')
      allocate(QML_Vibronic_t :: QModel%QM)
      IF (size(tab_pot_name) > 1) THEN
        QModel%QM = Init_QML_Vibronic(QModel_in,read_param=read_nml,nio_param_file=nio_loc, &
                                      Vibronic_name=trim(tab_pot_name(2)))
      ELSE
        QModel%QM = Init_QML_Vibronic(QModel_in,read_param=read_nml,nio_param_file=nio_loc)
      END IF

    CASE ('uracyl','uracil')
      !! === README ==
      !! Model for the photo-dissociation of the uracil cation.
      !! pot_name  = 'Uracil'
      !! ndim      = 36
      !! nsurf     = 4
      !! ref1: Mariana Assmann, Horst Köppel, and Spiridoula Matsika,
      !!       J. Phys. Chem. A 2015, 119, 866−875, 
      !!       DOI: 10.1021/jp512221x
      !! ref2: Patricia Vindel Zandbergen, Spiridoula Matsika, and Neepa T. Maitra
      !!       J. Phys. Chem. Lett. 2022, 13, 7, 1785–1790
      !!       https://doi.org/10.1021/acs.jpclett.1c04132
      !! === END README ==

      allocate(QML_Uracil_t :: QModel%QM)
      QModel%QM = Init_QML_Uracil(QModel_in,read_param=read_nml,nio_param_file=nio_loc)


    CASE ('fulvene')
      allocate(QML_fulvene_t :: QModel%QM)
      QModel%QM = Init_QML_fulvene(QModel_in,read_param=read_nml,nio_param_file=nio_loc)
    CASE ('dmabn')
      allocate(QML_dmabn_t :: QModel%QM)
      QModel%QM = Init_QML_dmabn(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

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
      QModel%QM = Init_QML_HNO3(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('no3')
      allocate(QML_NO3_t :: QModel%QM)
      QModel%QM = Init_QML_NO3(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

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
      QModel%QM = Init_QML_CH5(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('ph4old')
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
      QModel%QM = Init_QML_PH4(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('ph4','ph4jo')
      allocate(QML_PH4jo_t :: QModel%QM)
      QModel%QM = Init_QML_PH4Jo(QModel_in,read_param=read_nml,nio_param_file=nio_loc)
    CASE ('hnnhp')
      allocate(QML_HNNHp_t :: QModel%QM)
      QModel%QM = Init_QML_HNNHp(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('h2sin')
      allocate(QML_H2SiN_t :: QModel%QM)
      QModel%QM = Init_QML_H2SiN(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('h2nsi')
      allocate(QML_H2NSi_t :: QModel%QM)
      QModel%QM = Init_QML_H2NSi(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

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
      QModel%QM = Init_QML_HOO_DMBE(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('h3_lsth','h3')
      !! === README ==
      !! H3 potential:
      !! pot_name  = 'H3'
      !! option    = 0,1,10,11 (LSTH)
      !! ndim      = 3   (the 3 H-H distances)
      !! nsurf     = 1
      !! Units: Energy in Hartree and distances in bohr.
      !! refs (option=0):
      !! P. Siegbahn, B. Liu,  J. Chem. Phys. 68, 2457(1978).
      !! D.G. Truhlar and C.J. Horowitz, J. Chem. Phys. 68, 2466 (1978); https://doi.org/10.1063/1.436019
      !! options  0 and 10 : 3D model with IRC functions (1 potential + parameters)
      !! options  0 and  1 : first IRC funtions fitted in polar representation (alpha)
      !! options 10 and 11 : second IRC funtions fitted with the sum and the difference (alpha)
      !! === END README ==
      allocate(QML_H3_t :: QModel%QM)
      QModel%QM = Init_QML_H3(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('hcn_murrell','cnh_murrell')
      !! === README ==
      !! CNH or HCN potential:
      !! pot_name  = 'CNH_Murrell'
      !! option    = 0 (3D-3distances, default), 1,11 (3D-Jacobi), 2,21 (1D-Jacobi MEP)
      !! ndim      = 3
      !! nsurf     = 1
      !! remarks: 
      !!   - Atomic order: C, N, H
      !!   - Cart_TO_Q is possible
      !!   - The options 11 and 21, the third coordinate is cos(theta)
      !! ref: J. N. Murrell, S. Carter and L. O. Halonene, J. Mol. Spectrosc. vo93 p307 1982
      !!  doi: https://doi.org/10.1016/0022-2852(82)90170-9
      !! === END README ==
      allocate(QML_CNH_Murrell_t :: QModel%QM)
      QModel%QM = Init_QML_CNH_Murrell(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('h2o','water')
      !! === README ==
      !! H2O potential:
      !! pot_name  = 'H2O'
      !! option    = 1(default 1)
      !! ndim      = 3
      !! nsurf     = 1
      !! refs: Quadratic model potential for H2O; TIPS force constants taken from:  
      !!       Dang and Pettitt, J. Chem. Phys. 91 (1987)
      !! === END README ==
      allocate(QML_H2O_t :: QModel%QM)
      QModel%QM = Init_QML_H2O(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('h2_h2on','h2_clathrate')
      !! === README ==
      !! H2_H2On potential:
      !! pot_name  = 'H2_H2On'
      !! option    = 1 (default 1)
      !! ndim      = 6
      !! nsurf     = 1
      !! === END README ==
      allocate(QML_H2_H2On_t :: QModel%QM)
      QModel%QM = Init_QML_H2_H2On(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('bottleneck','eckart')
      !! === README ==
      !!Bottleneck potential: 1D Eckart Barrier + quadratic contributions'
      !! pot_name  = 'Bottleneck' or 'Eckart'
      !! option    = 1, 2 (default 2)
      !! ndim      >= 1
      !! nsurf     = 1
      !!
      !! ref (option 1): Trahan, Wyatt and Poirier, J Chem Phys 122, 164104 (2005)'
      !!   Multidimensional quantum trajectories: Applications of the derivative propagation method.'
      !! ref (option 2): Dupuy, Lauvergnat and Scribano, CPL 787, 139241 (2022)'
      !!   Smolyak representations with absorbing boundary conditions ...'
      !!       for reaction path Hamiltonian model of reactive scattering.'
      !!   DOI: 10.1016/j.cplett.2021.139241'
      !! === END README ==
      allocate(QML_Bottleneck_t :: QModel%QM)
      QModel%QM = Init_QML_Bottleneck(QModel_in,read_param=read_nml,nio_param_file=nio_loc)
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
      QModel%QM = Init_QML_ClH2p_Botschwina(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('chfclbr')
      !! === README ==
      !! CHFClBr QFF potential:
      !! pot_name  = 'CHFClBr'
      !! option    = no option (yet)
      !! ndim      = 9
      !! nsurf     = 1
      !! refs: unpublished
      !! === END README ==
      allocate(QML_CHFClBr_t :: QModel%QM)
      QModel%QM = Init_QML_CHFClBr(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('template')
      !! 3D-potential with 1 surface
      allocate(QML_Template_t :: QModel%QM)
      QModel%QM = Init_QML_Template(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('test')
      !! test-potential
      allocate(QML_Test_t :: QModel%QM)
      QModel%QM = Init_QML_Test(QModel_in,read_param=read_nml,nio_param_file=nio_loc)
    CASE ('extmodel')
      !! external-potential
      allocate(QML_ExtModel_t :: QModel%QM)
      QModel%QM = Init_QML_ExtModel(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE ('oned_photons')
      allocate(QML_OneD_Photons_t :: QModel%QM)
      QModel%QM = Init_QML_OneD_Photons(QModel_in,read_param=read_nml,nio_param_file=nio_loc)
    CASE ('oned_photons2')
      allocate(QML_OneD_Photons2_t :: QModel%QM)
      QModel%QM = Init_QML_OneD_Photons2(QModel_in,read_param=read_nml,nio_param_file=nio_loc)

    CASE DEFAULT
        write(out_unit,*) ' ERROR in Init_Model'
        write(out_unit,*) ' This model/potential is unknown. pot_name: ',pot_name_loc
        STOP 'STOP in Init_Model: Other potentials have to be done'
    END SELECT

    IF (present(ndim)) THEN
      IF (ndim > QModel%QM%ndim) THEN
          write(out_unit,*) ' ERROR in Init_Model'
          write(out_unit,*) ' ndim is present and ...'
          write(out_unit,*) ' its value is larger than QModel%QM%ndim'
          write(out_unit,*) ' ndim,QModel%QM%ndim',ndim,QModel%QM%ndim
          write(out_unit,*) ' check your data!'
          STOP 'STOP in Init_Model: wrong ndim'
      END IF
      IF (ndim < QModel%QM%ndim  .AND. ndim > 0) THEN
          write(out_unit,*) ' WARNING in Init_Model'
          write(out_unit,*) ' ndim is present and ...'
          write(out_unit,*) ' its value is smaller than QModel%QM%ndim'
          write(out_unit,*) ' ndim,QModel%QM%ndim',ndim,QModel%QM%ndim
          write(out_unit,*) ' => We assume that ...'
          write(out_unit,*) ' ... all variables (QModel%QM%ndim) will be given!'
      END IF
    END IF
    IF (present(nsurf)) THEN
      IF (nsurf /= QModel%QM%nsurf .AND. nsurf > 0) THEN
          write(out_unit,*) ' ERROR in Init_Model'
          write(out_unit,*) ' nsurf is present and ...'
          write(out_unit,*) ' its value is not equal to QModel%QM%nsurf'
          write(out_unit,*) ' nsurf,QModel%QM%nsurf',nsurf,QModel%QM%nsurf
          write(out_unit,*) ' check your data!'
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
        write(out_unit,*) ' ERROR in Init_Model'
        write(out_unit,*) ' Some values of list_act(:) are out of range.'
        write(out_unit,*) '   list_act(:): ',QModel%QM%list_act(:)
        write(out_unit,*) '   range = [1,',QModel%QM%ndim,']'
        write(out_unit,*) ' check your data!'
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
        write(out_unit,*) ' ERROR in Init_Model'
        write(out_unit,*) ' Some coordinate indexes are missing in ...'
        write(out_unit,*) ' ... list_act(:):   ',QModel%QM%list_act(:)
        write(out_unit,*) ' and list_inact(:): ',QModel%QM%list_inact(:)
        write(out_unit,*) ' check your data!'
        STOP 'STOP in Init_Model: wrong list_act'
      END IF

      write(out_unit,*) 'Vib_adia   ',QModel%QM%Vib_adia
      write(out_unit,*) 'nb_Channels',QModel%QM%nb_Channels
      write(out_unit,*) 'list_act   ',QModel%QM%list_act
      write(out_unit,*) 'list_inact ',QModel%QM%list_inact

      QModel%ndim  = size(QModel%QM%list_act)
      QModel%nsurf = QModel%QM%nb_Channels

      allocate(QModel%Basis)
      CALL Read_Basis(QModel%Basis,nio_loc)

    END IF

    IF (Print_init_loc) THEN
      write(out_unit,*) '================================================='
      write(out_unit,*) ' Quantum Model'
      CALL Write_Model(QModel,nio=out_unit)
      write(out_unit,*) '================================================='
      write(out_unit,*) '================================================='
      flush(out_unit)
    END IF

    IF (read_param_loc .AND. nio_loc /= in_unit) THEN
       close(unit=nio_loc)
    END IF

  END SUBROUTINE Init_Model
  SUBROUTINE dealloc_Model(Model)
    USE AdiaChannels_Basis_m, ONLY : QML_Basis_t,dealloc_Basis
    IMPLICIT NONE

    TYPE(Model_t),      intent(inout)           :: Model

    IF (allocated(Model%QM))    deallocate(Model%QM)
    IF (allocated(Model%Basis)) CALL dealloc_Basis(Model%Basis)
    Model%nsurf       = 0
    Model%ndim        = 0
    Model%opt         = .FALSE.
    Model%irc         = .FALSE.

  END SUBROUTINE dealloc_Model
  SUBROUTINE Set_step_epsi_Model(step_in,epsi_in)
    IMPLICIT NONE

    real (kind=Rkind),  intent(in), optional  :: step_in,epsi_in


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Set_step_epsi_Model'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      flush(out_unit)
    END IF

    IF (present(step_in)) THEN
      write(out_unit,*) ' WARNING: step has been changed.'
      write(out_unit,*) ' Old and new values',step,step_in
      step = step_in
    END IF

    IF (present(epsi_in)) THEN
      write(out_unit,*) ' WARNING: epsi has been changed.'
      write(out_unit,*) ' Old and new values',epsi,epsi_in
      epsi = epsi_in
    END IF

    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Set_step_epsi_Model

  SUBROUTINE get_Q0_Model(Q0,Model,option)
  USE QDUtil_m,         ONLY : Identity_Mat, TO_string, Write_Vec
  USE QML_Empty_m
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q0(:)
    TYPE (Model_t),     intent(in)               :: Model
    integer,            intent(in)               :: option


    integer :: err_Q0
!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='get_Q0_Model'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      flush(out_unit)
    END IF

    CALL check_alloc_QM(Model,name_sub)

    IF (size(Q0) /= Model%QM%ndim) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' The size of Q0 is not Model%QM%ndim: '
      write(out_unit,*) ' size(Q0)',size(Q0)
      write(out_unit,*) ' ndim',Model%QM%ndim
      STOP 'STOP in get_Q0_Model: Wrong Q0 size'
    END IF

    Q0(:) = ZERO

    CALL get_Q0_QModel(Model%QM,Q0,err_Q0)
    IF (err_Q0 /= 0) THEN
      CALL Write_Model(Model,out_unit)
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Q0 is not set-up in the model'
      STOP 'STOP Q0 is not set-up in the model'
    END IF

    IF (debug) THEN
      CALL Write_Vec(Q0,out_unit,nbcol=5,info='Q0: ')
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE get_Q0_Model
  FUNCTION get_d0GGdef_Model(Model) RESULT(d0GGdef)
    IMPLICIT NONE

    real (kind=Rkind),   allocatable              :: d0GGdef(:,:)
    TYPE (Model_t),                  intent(in)   :: Model

    d0GGdef = Model%QM%get_d0GGdef_QModel()

  END FUNCTION get_d0GGdef_Model

  ! check if the QM [CLASS(QML_Empty_t)] is allocated
  SUBROUTINE check_alloc_QM(QModel,name_sub_in)
    IMPLICIT NONE

    TYPE (Model_t),     intent(in)     :: QModel
    character (len=*),  intent(in)     :: name_sub_in

    IF ( .NOT. allocated(QModel%QM)) THEN
      write(out_unit,*) ' ERROR in check_alloc_QM'
      write(out_unit,*) ' QM is not allocated in QModel.'
      write(out_unit,*) '  check_alloc_QM is called from ',name_sub_in
      write(out_unit,*) '  You MUST initialize the model with:'
      write(out_unit,*) '    CALL init_Model(...) in Model_m.f90'
      write(out_unit,*) ' or'
      write(out_unit,*) '    CALL sub_Init_Qmodel(...) in Model_driver.f90'
      STOP 'ERROR in check_alloc_QM: QM is not allocated in QModel.'
    END IF

  END SUBROUTINE check_alloc_QM

  ! check if the check_Init_QModel [TYPE(Model_t)] is initialized
  FUNCTION check_Init_QModel(Model)
    IMPLICIT NONE

    logical                             :: check_Init_QModel
    TYPE (Model_t),     intent(in)      :: Model

    check_Init_QModel = allocated(Model%QM)
    IF (allocated(Model%QM)) THEN
      check_Init_QModel = Model%QM%init
    END IF

  END FUNCTION check_Init_QModel
  FUNCTION check_alloc_d0GGdef(QModel) RESULT(alloc)
    IMPLICIT NONE

    logical                           :: alloc
    CLASS(QML_Empty_t), intent(in)    :: QModel

    alloc = allocated(QModel%d0GGdef)

  END FUNCTION check_alloc_d0GGdef
  SUBROUTINE Eval_tab_HMatVibAdia(Model,Qact,tab_MatH)
    USE QDUtil_m,         ONLY : Write_Mat
    USE ADdnSVM_m
    IMPLICIT NONE

  TYPE (Model_t),                 intent(inout)            :: Model
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
    write(out_unit,*) ' BEGINNING ',name_sub
    flush(out_unit)
  END IF

  CALL check_alloc_QM(Model,name_sub)

  PF = Model%QM%Phase_Following
  PC = Model%QM%Phase_Checking
  CALL Eval_dnHVib_ana(Model,Qact,PotVal_dia,nderiv=2)
  IF (.NOT. allocated(Model%QM%Vec0)) allocate(Model%QM%Vec0)
  CALL dia_TO_adia(PotVal_dia,PotVal,Vec,Model%QM%Vec0,NAC,PF,PC,nderiv=2)

  !CALL Write_dnMat(PotVal,nio=out_unit,info='PotVal (adia)')

  !Mat_diag = matmul(transpose(Vec),matmul(PotVal_dia,Vec))
  !CALL Write_dnMat(Mat_diag,nio=out_unit,info='Mat_diag')

  !write(out_unit,*) 'nsurf,ndim',Model%nsurf,Model%ndim
  nsurf    = Model%nsurf
  nb_terms = (Model%ndim + 1)*(Model%ndim + 2)/2
  IF (.NOT. allocated(tab_MatH)) THEN
    allocate(tab_MatH(nsurf,nsurf,nb_terms))
  END IF

  !write(out_unit,*) Qact,'NAC1',NAC%d1(1:nsurf,1:nsurf,:)


  IF ( any([nsurf,nsurf,nb_terms] /= shape(tab_MatH)) ) THEN
    write(out_unit,*) ' ERROR in ',name_sub
    write(out_unit,*) ' The shape of tab_MatH is wrong.'
    write(out_unit,*) '    shape(tab_MatH):',shape(tab_MatH)
    write(out_unit,*) '  It MUST be:      [',nsurf,nsurf,nb_terms,']'
    write(out_unit,*) '  Check your data or the code!'
    STOP 'ERROR in Eval_tab_HMatVibAdia: The shape of tab_MatH is wrong.'
  END IF
  tab_MatH(:,:,:) = ZERO

  ndim_act    = size(Model%QM%list_act)
  ndim_inact  = size(Model%QM%list_inact)

  ! Veff
  d0GGdef_aa = Model%QM%d0GGdef(Model%QM%list_act,Model%QM%list_act)
  iterm = 1
  tab_MatH(:,:,iterm) = PotVal%d0(1:nsurf,1:nsurf)
  DO ia=1,ndim_act
  DO ja=1,ndim_act
    tab_MatH(:,:,iterm) = tab_MatH(:,:,iterm) - HALF*d0GGdef_aa(ja,ia) *        &
       NAC%d2(1:nsurf,1:nsurf,ja,ia)
  END DO
  END DO

  ! F2^jaia
  d0GGdef_aa = Model%QM%d0GGdef(Model%QM%list_act,Model%QM%list_act)
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
  d0GGdef_aa = Model%QM%d0GGdef(Model%QM%list_act,Model%QM%list_act)
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
      write(out_unit,*) iterm
      CALL Write_Mat(tab_MatH(:,:,iterm),nio=out_unit,nbcol=5)
    END DO
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF

  END SUBROUTINE Eval_tab_HMatVibAdia

  SUBROUTINE Eval_Pot(Model,Q,PotVal,nderiv,NAC,Vec,numeric,PotVal_dia,Vec0)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (Model_t),     intent(inout)            :: Model
    TYPE (dnMat_t),     intent(inout)            :: PotVal
    real (kind=Rkind),  intent(in)               :: Q(:)
    integer,            intent(in),    optional  :: nderiv
    TYPE (dnMat_t),     intent(inout), optional  :: NAC,Vec
    TYPE (dnMat_t),     intent(inout), optional  :: Vec0
    logical,            intent(in),    optional  :: numeric
    TYPE (dnMat_t),     intent(inout), optional  :: PotVal_dia

    ! local variables
    integer                    :: i,nderiv_loc
    TYPE (dnMat_t)             :: Vec_loc,NAC_loc,PotVal_dia_loc,PotVal_loc
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
    write(out_unit,*) ' BEGINNING ',name_sub
    IF (present(nderiv)) write(out_unit,*) '   nderiv',nderiv
    write(out_unit,*) '   present(NAC): ',present(NAC)
    write(out_unit,*) '   present(Vec): ',present(Vec)
    write(out_unit,*) '   present(Vec0):',present(Vec0)
    IF (present(Vec0)) CALL Write_dnMat(Vec0,nio=out_unit)
    flush(out_unit)
  END IF
  !write(out_unit,*) 'in Eval_Pot Q ',Q

  CALL check_alloc_QM(Model,name_sub)
  PF = Model%QM%Phase_Following
  PC = Model%QM%Phase_Checking
  IF (debug) THEN
    write(out_unit,*) '  Model%QM%numeric         ',Model%QM%numeric
    write(out_unit,*) '  Model%QM%adiabatic       ',Model%QM%adiabatic
    write(out_unit,*) '  Model%QM%Vib_adia        ',Model%QM%Vib_adia
    write(out_unit,*) '  Model%QM%Phase_Following ',Model%QM%Phase_Following
    write(out_unit,*) '  Model%QM%Phase_Checking  ',Model%QM%Phase_Checking
    flush(out_unit)
  END IF

  IF (present(nderiv)) THEN
    nderiv_loc = max(0,nderiv)
    nderiv_loc = min(3,nderiv_loc)
  ELSE
    nderiv_loc = 0
  END IF

  IF (present(numeric)) THEN
    numeric_loc = (numeric  .OR. Model%QM%no_ana_der)
  ELSE
    numeric_loc = (Model%QM%numeric .OR. Model%QM%no_ana_der)
  END IF
  numeric_loc = (numeric_loc .AND. nderiv_loc > 0)


  IF (Model%QM%Vib_adia) THEN
    IF (present(Vec)) THEN
      write(out_unit,*) ' ERROR in ',name_sub
      write(out_unit,*) ' Vib_adia=t and Vec is present'
      write(out_unit,*) ' This is not possible yet !'
      STOP 'ERROR in Eval_Pot: Vib_adia=t is not compatible with Vec'
    END IF

    CALL Eval_dnHVib_ana(Model,Q,PotVal_dia_loc,nderiv_loc)

    !write(out_unit,*) 'PotVal (Vib_dia)'
    !CALL Write_dnMat(PotVal_dia_loc,nio=out_unit)

    IF (present(Vec0)) THEN
      CALL dia_TO_adia(PotVal_dia_loc,PotVal_loc,Vec_loc,Vec0,NAC_loc,      &
                       PF,PC,nderiv_loc,type_diag=1)
    ELSE IF (.NOT. allocated(Model%QM%Vec0)) THEN 
      allocate(Model%QM%Vec0)
      CALL dia_TO_adia(PotVal_dia_loc,PotVal_loc,Vec_loc,Model%QM%Vec0,NAC_loc,      &
                     PF,PC,nderiv_loc,type_diag=1)
    END IF
    CALL submatrix_dnMat2_TO_dnMat1(PotVal,PotVal_loc,lb=1,ub=Model%QM%nb_Channels)

    IF (present(Vec)) THEN
      ! it needs dnMat as a rectangular matrix !!
      CALL submatrix_dnMat2_TO_dnMat1(Vec,Vec_loc,lb=1,ub=Model%QM%nb_Channels)
    END IF

    IF (present(NAC)) THEN
      CALL submatrix_dnMat2_TO_dnMat1(NAC,NAC_loc,lb=1,ub=Model%QM%nb_Channels)
    END IF

    ! print the Vec%d0 if required
    CALL Write_QML_EigenVec(Q,Vec_loc,Model,nio=out_unit)

    IF (present(PotVal_dia)) PotVal_dia = PotVal_dia_loc

    CALL dealloc_dnMat(NAC_loc)
    CALL dealloc_dnMat(Vec_loc)
    CALL dealloc_dnMat(PotVal_loc)
    CALL dealloc_dnMat(PotVal_dia_loc)

  ELSE

    adia_loc = (Model%QM%adiabatic .AND. Model%QM%nsurf > 1)

    IF (numeric_loc) THEN  ! numerical
      IF (.NOT. adia_loc) THEN
         SELECT CASE (numeric_option)
         CASE (0)
           CALL Eval_Pot_Numeric_dia_old(Model,Q,PotVal,nderiv_loc)
         CASE (3)
           CALL Eval_Pot_Numeric_dia_v3(Model,Q,PotVal,nderiv_loc)
         CASE (4)
           CALL Eval_Pot_Numeric_dia_v4(Model,Q,PotVal,nderiv_loc)
         CASE Default
           CALL Eval_Pot_Numeric_dia_old(Model,Q,PotVal,nderiv_loc)
         END SELECT
      ELSE

        IF (present(Vec)) THEN
          IF (present(NAC)) THEN
            CALL Eval_Pot_Numeric_adia(Model,Q,PotVal,nderiv_loc,      &
                                                 Vec,NAC,numeric_option)
          ELSE
            CALL Eval_Pot_Numeric_adia(Model,Q,PotVal,nderiv_loc,      &
                                             Vec,NAC_loc,numeric_option)
            CALL dealloc_dnMat(NAC_loc)
          END IF
        ELSE
          IF (present(NAC)) THEN
            CALL Eval_Pot_Numeric_adia(Model,Q,PotVal,nderiv_loc,      &
                                             Vec_loc,NAC,numeric_option)
          ELSE
            CALL Eval_Pot_Numeric_adia(Model,Q,PotVal,nderiv_loc,      &
                                         Vec_loc,NAC_loc,numeric_option)
            CALL dealloc_dnMat(NAC_loc)
          END IF
          CALL dealloc_dnMat(Vec_loc)
        END IF
      END IF
    ELSE ! analytical calculation
      IF (present(Vec0)) THEN
        IF (present(Vec)) THEN
          IF (present(NAC)) THEN
            CALL Eval_Pot_ana(Model,Q,PotVal,nderiv_loc,Vec=Vec,Nac=NAC,PotVal_dia=PotVal_dia_loc,Vec0=Vec0)
          ELSE
            CALL Eval_Pot_ana(Model,Q,PotVal,nderiv_loc,Vec=Vec,PotVal_dia=PotVal_dia_loc,Vec0=Vec0)
          END IF
        ELSE
          IF (present(NAC)) THEN
            CALL Eval_Pot_ana(Model,Q,PotVal,nderiv_loc,Nac=NAC,PotVal_dia=PotVal_dia_loc,Vec0=Vec0)
          ELSE
            CALL Eval_Pot_ana(Model,Q,PotVal,nderiv_loc,PotVal_dia=PotVal_dia_loc,Vec0=Vec0)
          END IF
        END IF
        IF (present(PotVal_dia)) PotVal_dia = PotVal_dia_loc
      ELSE
        IF (present(Vec)) THEN
          IF (present(NAC)) THEN
            CALL Eval_Pot_ana(Model,Q,PotVal,nderiv_loc,Vec=Vec,Nac=NAC,PotVal_dia=PotVal_dia_loc)
          ELSE
            CALL Eval_Pot_ana(Model,Q,PotVal,nderiv_loc,Vec=Vec,PotVal_dia=PotVal_dia_loc)
          END IF
        ELSE
          IF (present(NAC)) THEN
            CALL Eval_Pot_ana(Model,Q,PotVal,nderiv_loc,Nac=NAC,PotVal_dia=PotVal_dia_loc)
          ELSE
            CALL Eval_Pot_ana(Model,Q,PotVal,nderiv_loc,PotVal_dia=PotVal_dia_loc)
          END IF
        END IF
        IF (present(PotVal_dia)) PotVal_dia = PotVal_dia_loc
      END IF
    END IF
  END IF

  IF (debug) THEN
    IF ( Model%QM%adiabatic) THEN
      write(out_unit,*) 'PotVal (adia)'
    ELSE
      write(out_unit,*) 'PotVal (dia)'
    END IF
    CALL Write_dnMat(PotVal,nio=out_unit)
    IF (present(Vec0)) THEN
      write(out_unit,*) 'Vec0'
      CALL Write_dnMat(Vec0,nio=out_unit)
    END IF
    write(out_unit,*) ' END ',name_sub
    flush(out_unit)
  END IF
  !CALL Write_dnMat(PotVal,nio=out_unit,info='in Eval_Pot: PotVal')

  END SUBROUTINE Eval_Pot

  SUBROUTINE Eval_Pot_ana(Model,Q,PotVal,nderiv,NAC,Vec,PotVal_dia,Vec0)
    USE ADdnSVM_m, ONLY : dnS_t,alloc_dnS,dealloc_dnS,Variable,                &
             dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat,Write_dnMat
    IMPLICIT NONE

    TYPE (Model_t),        intent(inout)            :: Model

    TYPE (dnMat_t),        intent(inout)            :: PotVal
    real (kind=Rkind),     intent(in)               :: Q(:)
    integer,               intent(in)               :: nderiv
    TYPE (dnMat_t),        intent(inout), optional  :: NAC,Vec,PotVal_dia
    TYPE (dnMat_t),        intent(inout), optional  :: Vec0

    ! local variables
    integer                     :: i,j,ij,id,nat
    TYPE (dnMat_t)              :: PotVal_dia_loc,Vec_loc,NAC_loc
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
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) '   nderiv:       ',nderiv
      write(out_unit,*) '   present(NAC): ',present(NAC)
      write(out_unit,*) '   present(Vec): ',present(Vec)
      write(out_unit,*) '   present(Vec0):',present(Vec0)
      IF (present(Vec0)) CALL Write_dnMat(Vec0,nio=out_unit)
      flush(out_unit)
    END IF

    CALL check_alloc_QM(Model,name_sub)

    PF = Model%QM%Phase_Following
    PC = Model%QM%Phase_Checking
    IF (debug) write(out_unit,*) '   adiabatic ',Model%QM%adiabatic


    IF ( Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,nderiv=nderiv)
    END IF
    PotVal = ZERO
    IF (debug) write(out_unit,*) '   init PotVal  ' ; flush(out_unit)

    ! allocate Mat_OF_PotDia
    allocate(Mat_OF_PotDia(Model%QM%nsurf,Model%QM%nsurf))
    DO j=1,size(Mat_OF_PotDia(1,:))
    DO i=1,size(Mat_OF_PotDia(:,1))
        CALL alloc_dnS(Mat_OF_PotDia(i,j),Model%QM%ndim,nderiv)
    END DO
    END DO
    IF (debug) write(out_unit,*) '   alloc Mat_OF_PotDia  ' ; flush(out_unit)

    ! intialization of the dnQ(:)
    IF (Model%QM%Cart_TO_Q) THEN
      !in Q(:) we have the cartesian coordinates
      allocate(dnQ(Model%QM%ndimQ))
      nat = int(Model%QM%ndim/3)
      allocate(dnX(3,nat))
      ij = 0
      DO i=1,nat
      DO j=1,3
        ij = ij + 1
        dnX(j,i) = Variable(Q(ij),nVar=Model%QM%ndim,nderiv=nderiv,iVar=ij) ! to set up the derivatives
      END DO
      END DO


      IF (Model%QM%AbInitio) THEN
        CALL Model%QM%EvalPotAbInitio_QModel(Mat_OF_PotDia,dnX,nderiv=nderiv)
      ELSE
        CALL Model%QM%Cart_TO_Q_QModel(dnX,dnQ,nderiv=nderiv)
      END IF

      CALL dealloc_dnS(dnX)
      deallocate(dnX)

    ELSE
      allocate(dnQ(Model%QM%ndim))
      DO i=1,Model%QM%ndim
        dnQ(i) = Variable(Q(i),nVar=Model%QM%ndim,nderiv=nderiv,iVar=i) ! to set up the derivatives
      END DO
    END IF
    IF (debug) write(out_unit,*) '   init dnQ(:)  ' ; flush(out_unit)

    IF (Model%QM%AbInitio) THEN
      IF (debug) write(out_unit,*) 'PotVal already done'
    ELSE
      CALL Model%QM%EvalPot_QModel(Mat_OF_PotDia,dnQ,nderiv=nderiv)
      IF (debug) write(out_unit,*) ' PotVal done' ; flush(out_unit)
    END IF

    PotVal = Mat_OF_PotDia ! transfert the potential and its derivatives to the matrix form (PotVal)
    IF (debug) write(out_unit,*) ' transfert done' ; flush(out_unit)

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

    IF ( Model%QM%adiabatic .AND. Model%QM%nsurf > 1) THEN
      IF (debug) THEN
        write(out_unit,*) 'PotVal (dia)'
        CALL Write_dnMat(PotVal,nio=out_unit)
        flush(out_unit)
      END IF

      PotVal_dia_loc = PotVal
      IF (present(PotVal_dia)) PotVal_dia = PotVal_dia_loc
      IF (debug) write(out_unit,*) ' save pot dia  done' ; flush(out_unit)

      IF (present(Vec0)) THEN
        IF (present(Vec)) THEN
          IF (present(NAC)) THEN
            CALL dia_TO_adia(PotVal_dia_loc,PotVal,Vec,Vec0,NAC,PF,PC,nderiv)
          ELSE
            CALL dia_TO_adia(PotVal_dia_loc,PotVal,Vec,Vec0,NAC_loc,PF,PC,nderiv)
            CALL dealloc_dnMat(NAC_loc)
          END IF
        ELSE
          IF (present(NAC)) THEN
            CALL dia_TO_adia(PotVal_dia_loc,PotVal,Vec_loc,Vec0,NAC,PF,PC,nderiv)
          ELSE
            CALL dia_TO_adia(PotVal_dia_loc,PotVal,Vec_loc,Vec0,NAC_loc,PF,PC,nderiv)
            CALL dealloc_dnMat(NAC_loc)
          END IF
          CALL dealloc_dnMat(Vec_loc)
        END IF
      ELSE
        IF (.NOT. allocated(Model%QM%Vec0)) allocate(Model%QM%Vec0)
  
        IF (present(Vec)) THEN
          IF (present(NAC)) THEN
            CALL dia_TO_adia(PotVal_dia_loc,PotVal,Vec,Model%QM%Vec0,NAC,PF,PC,nderiv)
          ELSE
            CALL dia_TO_adia(PotVal_dia_loc,PotVal,Vec,Model%QM%Vec0,NAC_loc,PF,PC,nderiv)
            CALL dealloc_dnMat(NAC_loc)
          END IF
        ELSE
          IF (present(NAC)) THEN
            CALL dia_TO_adia(PotVal_dia_loc,PotVal,Vec_loc,Model%QM%Vec0,NAC,PF,PC,nderiv)
          ELSE
            CALL dia_TO_adia(PotVal_dia_loc,PotVal,Vec_loc,Model%QM%Vec0,NAC_loc,PF,PC,nderiv)
            CALL dealloc_dnMat(NAC_loc)
          END IF
          CALL dealloc_dnMat(Vec_loc)
        END IF
      END IF
      IF (debug) write(out_unit,*) ' dia => adia  done' ; flush(out_unit)
      CALL dealloc_dnMat(PotVal_dia_loc)
    END IF


    IF (debug) THEN
      IF ( Model%QM%adiabatic) THEN
        write(out_unit,*) 'PotVal (adia)'
      ELSE
        write(out_unit,*) 'PotVal (dia)'
      END IF
      CALL Write_dnMat(PotVal,nio=out_unit)
      IF (present(Vec0)) THEN
        write(out_unit,*) 'Vec0'
        CALL Write_dnMat(Vec0,nio=out_unit)
      END IF
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Eval_Pot_ana

  SUBROUTINE Eval_Pot_Numeric_dia_v4(Model,Q,PotVal,nderiv)
    USE QMLLib_FiniteDiff_m
    USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: Model
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j,k,ip,jp,kp
    integer                            :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)

    CALL check_alloc_QM(Model,'Eval_Pot_Numeric_dia_v4')


    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    allocate(Q_loc(Model%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,nderiv=0)
    PotVal = ZERO

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(Model,Q,PotVal_loc0,nderiv=0)
    CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=4)

    IF (nderiv >= 1) THEN ! along ONE coordinates (first derivatives and higher)

      ! Numeric evaluation of forces
      DO i=1,Model%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                          indQ=[i],indDQ=ind1DQ,option=4)
        END DO

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,Model%QM%ndim
      DO j=1,Model%QM%ndim
        IF (i == j) CYCLE

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                        indQ=[i,j],indDQ=ind2DQ,option=4)
        END DO

      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN ! 3d derivatives:d3/dQidQjdQk

      DO i=1,Model%QM%ndim
      DO j=1,Model%QM%ndim
      IF (i == j) CYCLE
      DO k=1,Model%QM%ndim
        IF (i == k .OR. j == k) CYCLE

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
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
  SUBROUTINE Eval_Pot_Numeric_dia_v3(Model,Q,PotVal,nderiv)
    USE QMLLib_FiniteDiff_m
    USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: Model
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j,k,ip,jp,kp

    integer                            :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)

    CALL check_alloc_QM(Model,'Eval_Pot_Numeric_dia_v3')


    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO


    allocate(Q_loc(Model%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(Model,Q,PotVal_loc0,nderiv=0)
    PotVal%d0 = PotVal_loc0%d0

    CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=3)


    IF (nderiv >= 1) THEN ! 1st derivatives

      DO i=1,Model%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                          indQ=[i],indDQ=ind1DQ,option=3)
        END DO

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,Model%QM%ndim
      DO j=i+1,Model%QM%ndim

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
          CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,        &
                                        indQ=[i,j],indDQ=ind2DQ,option=3)
        END DO

        CALL FiniteDiff3_SymPerm_OF_dnMat(PotVal,indQ=[i,j])

      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN ! 3d derivatives: d3/dQidQidQj

      ! d3/dQidQjdQk
      DO i=1,Model%QM%ndim
      DO j=i+1,Model%QM%ndim
      DO k=j+1,Model%QM%ndim

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
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

  SUBROUTINE Eval_Pot_Numeric_dia_old(Model,Q,PotVal,nderiv)
    USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: Model
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j

    CALL check_alloc_QM(Model,'Eval_Pot_Numeric_dia_old')


    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    allocate(Q_loc(Model%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(Model,Q,PotVal_loc0,nderiv=0)
    PotVal%d0 = PotVal_loc0%d0


    IF (nderiv >= 1) THEN ! 1st derivatives

      ! Numeric evaluation of forces
      DO i=1,Model%QM%ndim

        Q_loc(i) = Q(i) + step        ! q+dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0) ! Ep(q+dq)
        PotVal%d1(:,:,i) = PotVal_loc0%d0

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = PotVal_loc0%d0
        END IF

        Q_loc(i) = Q(i) - step        ! q-dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0) ! Ep(q-dq)
        PotVal%d1(:,:,i) = (PotVal%d1(:,:,i)-PotVal_loc0%d0)/(TWO*step)


        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = (PotVal%d2(:,:,i,i) + PotVal_loc0%d0 - TWO*PotVal%d0)/ &
                                step**2
        END IF

        Q_loc(i) = Q(i)

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,Model%QM%ndim
      DO j=i+1,Model%QM%ndim

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) + PotVal_loc0%d0

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0)
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

  SUBROUTINE Eval_Pot_Numeric_adia(Model,Q,PotVal,nderiv,Vec,NAC,option)
    USE ADdnSVM_m, ONLY : dnMat_t
    IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: Model
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
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) '   nderiv',nderiv
      write(out_unit,*) '   option',option
      flush(out_unit)
    END IF

    SELECT CASE (option)
    CASE (0)
      CALL Eval_Pot_Numeric_adia_old(Model,Q,PotVal,nderiv,Vec,NAC)
    CASE (3)
      CALL Eval_Pot_Numeric_adia_v3(Model,Q,PotVal,nderiv,Vec,NAC)
    CASE (4)
      STOP 'Eval_Pot_Numeric_adia: option=4, not yet'
    !  CALL Eval_Pot_Numeric_adia_v4(Model,Q,PotVal,nderiv,Vec,NAC)
    CASE Default
      CALL Eval_Pot_Numeric_adia_old(Model,Q,PotVal,nderiv,Vec,NAC)
    END SELECT

    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Eval_Pot_Numeric_adia
  SUBROUTINE Eval_Pot_Numeric_adia_old(Model,Q,PotVal,nderiv,Vec,NAC)
    USE QDUtil_m,  ONLY : Identity_Mat
    USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: Model
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv
    TYPE (dnMat_t),    intent(inout)  :: Vec,NAC

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0,Vec_loc0
    integer                            :: i,j
    real (kind=Rkind), allocatable     :: tVec(:,:)

    CALL check_alloc_QM(Model,'Eval_Pot_Numeric_adia_old')


    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                           nderiv=nderiv)
    END IF
    PotVal = ZERO

    IF (Check_NotAlloc_dnMat(Vec,nderiv) ) THEN
      CALL alloc_dnMat(Vec,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                          nderiv=nderiv)
    END IF
    Vec = ZERO

    IF (Check_NotAlloc_dnMat(NAC,nderiv) ) THEN
      CALL alloc_dnMat(NAC,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                          nderiv=nderiv)
    END IF
    NAC = ZERO

    allocate(Q_loc(Model%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,nderiv=0)
    CALL alloc_dnMat(Vec_loc0,   nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(Model,Q,PotVal_loc0,nderiv=0,vec=Vec_loc0)

    PotVal%d0      =  PotVal_loc0%d0
    Vec%d0         = Vec_loc0%d0
    NAC%d0         = Identity_Mat(Model%QM%nsurf)

    allocate(tVec(Model%QM%nsurf,Model%QM%nsurf))
    tVec(:,:)      = transpose(Vec%d0)

    IF (nderiv >= 1) THEN ! 1st derivatives

      ! Numeric evaluation of forces
      DO i=1,Model%QM%ndim

        Q_loc(i) = Q(i) + step        ! q+dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0) ! Ep(q+dq)
        CALL Change_EigenVecPhase(Vec=Vec_loc0%d0,Vec0=Vec%d0)

        PotVal%d1(:,:,i) = PotVal_loc0%d0
        Vec%d1(:,:,i)    = Vec_loc0%d0

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = PotVal_loc0%d0
          Vec%d2(:,:,i,i)    = Vec_loc0%d0
        END IF

        Q_loc(i) = Q(i) - step        ! q-dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0) ! Ep(q-dq)
        CALL Change_EigenVecPhase(Vec=Vec_loc0%d0,Vec0=Vec%d0)

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

      DO i=1,Model%QM%ndim
      DO j=i+1,Model%QM%ndim

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        CALL Change_EigenVecPhase(Vec_loc0%d0,Vec_loc0%d0)

        PotVal%d2(:,:,j,i) = PotVal_loc0%d0
        Vec%d2(:,:,j,i)    = Vec_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        CALL Change_EigenVecPhase(Vec=Vec_loc0%d0,Vec0=Vec%d0)

        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) + PotVal_loc0%d0
        Vec%d2(:,:,j,i)    = Vec%d2(:,:,j,i)    + Vec_loc0%d0

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        CALL Change_EigenVecPhase(Vec=Vec_loc0%d0,Vec0=Vec%d0)

        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0
        Vec%d2(:,:,j,i)   = Vec%d2(:,:,j,i)     - Vec_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        CALL Change_EigenVecPhase(Vec=Vec_loc0%d0,Vec0=Vec%d0)

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
  SUBROUTINE Eval_Pot_Numeric_adia_v3(Model,Q,PotVal,nderiv,Vec,NAC)
    USE QDUtil_m,         ONLY : Identity_Mat
    USE QMLLib_FiniteDiff_m
    USE ADdnSVM_m,        ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Check_NotAlloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: Model
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

    CALL check_alloc_QM(Model,'Eval_Pot_Numeric_adia_v3')

    IF (Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                           nderiv=nderiv)
    END IF
    PotVal = ZERO

    IF (Check_NotAlloc_dnMat(Vec,nderiv) ) THEN
      CALL alloc_dnMat(Vec,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                          nderiv=nderiv)
    END IF
    Vec = ZERO

    IF (Check_NotAlloc_dnMat(NAC,nderiv) ) THEN
      CALL alloc_dnMat(NAC,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,&
                          nderiv=nderiv)
    END IF
    NAC = ZERO

    allocate(Q_loc(Model%QM%ndim))
    Q_loc(:) = Q
    CALL alloc_dnMat(PotVal_loc0,nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,nderiv=0)
    CALL alloc_dnMat(Vec_loc0,   nsurf=Model%QM%nsurf,nVar=Model%QM%ndim,nderiv=0)
    !write(out_unit,*) 'coucou1 Eval_Pot_Numeric_adia_v3' ; flush(6)


    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(Model,Q,PotVal_loc0,nderiv=0,vec=Vec_loc0)
    !write(out_unit,*) 'coucou1.1 Eval_Pot_Numeric_adia_v3' ; flush(6)



    CALL FiniteDiff_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=3)
    CALL FiniteDiff_AddMat_TO_dnMat(Vec,   Vec_loc0%d0,   option=3)
    !write(out_unit,*) 'coucou1.2 Eval_Pot_Numeric_adia_v3' ; flush(6)

    NAC%d0 = Identity_Mat(Model%QM%nsurf)
    !write(out_unit,*) 'coucou1.3 Eval_Pot_Numeric_adia_v3' ; flush(6)

    allocate(tVec(Model%QM%nsurf,Model%QM%nsurf))
    tVec(:,:)      = transpose(Vec%d0)
    !write(out_unit,*) 'coucou2 0-order Eval_Pot_Numeric_adia_v3' ; flush(6)

    IF (nderiv >= 1) THEN ! 1st derivatives

      DO i=1,Model%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL Change_EigenVecPhase(Vec=Vec_loc0%d0,Vec0=Vec%d0)

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
    !write(out_unit,*) 'coucou2 1-order Eval_Pot_Numeric_adia_v3' ; flush(6)

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,Model%QM%ndim
      DO j=i+1,Model%QM%ndim

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL Change_EigenVecPhase(Vec=Vec_loc0%d0,Vec0=Vec%d0)

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
    !write(out_unit,*) 'coucou2 2-order Eval_Pot_Numeric_adia_v3' ; flush(6)

    IF (nderiv >= 3) THEN ! 3d derivatives: d3/dQidQidQj

      ! d3/dQidQjdQk
      DO i=1,Model%QM%ndim
      DO j=i+1,Model%QM%ndim
      DO k=j+1,Model%QM%ndim

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(Model,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL Change_EigenVecPhase(Vec=Vec_loc0%d0,Vec0=Vec%d0)

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
    !write(out_unit,*) 'coucou2 3-order Eval_Pot_Numeric_adia_v3' ; flush(6)

    CALL FiniteDiff_Finalize_dnMat(PotVal,step)
    CALL FiniteDiff_Finalize_dnMat(Vec,step)
    CALL FiniteDiff_Finalize_dnMat(NAC,step)
    !write(out_unit,*) 'coucou3 Eval_Pot_Numeric_adia_v3' ; flush(6)

    deallocate(tVec)
    deallocate(Q_loc)
    CALL dealloc_dnMat(PotVal_loc0)
    CALL dealloc_dnMat(Vec_loc0)
    !write(out_unit,*) 'coucouf Eval_Pot_Numeric_adia_v3' ; flush(6)

  END SUBROUTINE Eval_Pot_Numeric_adia_v3
  SUBROUTINE Change_EigenVecPhase(Vec,Vec0)
    USE QDUtil_m
    IMPLICIT NONE

    real(kind=Rkind), intent(in)    :: Vec0(:,:)
    real(kind=Rkind), intent(inout) :: Vec(:,:)

    integer :: i
    real(kind=Rkind) :: Sii

    DO i=1,size(Vec0,dim=2)
      Sii = dot_product(Vec0(:,i),Vec(:,i))
      IF (Sii < 0) Vec(:,i) = -Vec(:,i)
    END DO

  END SUBROUTINE Change_EigenVecPhase
  SUBROUTINE dia_TO_adia(PotVal_dia,PotVal_adia,Vec,Vec0,NAC,Phase_Following,   &
                         Phase_checking,nderiv,type_diag)

    USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,Write_dnMat,DIAG_dnMat,  &
      Check_NotAlloc_dnMat,get_nsurf,get_nVar
    USE QDUtil_m,         ONLY : diagonalization

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
      write(out_unit,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unit,*) '   nderiv',nderiv
      !n = size(PotVal_dia%d0,dim=1)
      write(out_unit,*) 'type_diag_loc',type_diag_loc
      write(out_unit,*) 'Phase_Checking',Phase_Checking
      write(out_unit,*) 'Phase_Following',Phase_Following
      flush(out_unit)
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(3,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF

    IF ( Check_NotAlloc_dnMat(PotVal_dia,nderiv_loc) ) THEN
      write(out_unit,*) ' The diabatic potential MUST be allocated!'
      CALL Write_dnMat(PotVal_dia)
      STOP 'PotVal_dia%dn NOT allocated in "dia_TO_adia"'
    END IF
    IF (debug) THEN
      write(out_unit,*) 'PotVal_dia'
      CALL Write_dnMat(PotVal_dia,nio=out_unit)
      flush(out_unit)
    END IF

    nsurf = get_nsurf(PotVal_dia)
    ndim  = get_nVar(PotVal_dia)


    IF (Check_NotAlloc_dnMat(Vec0,nderiv=0)) THEN
       !$OMP CRITICAL (CRIT_dia_TO_adia)
       CALL alloc_dnMat(Vec0,nsurf=nsurf,nVar=ndim,nderiv=0)

       allocate(Eig(nsurf))

       CALL diagonalization(PotVal_dia%d0,Eig,Vec0%d0,nsurf,sort=1,phase=.TRUE.,diago_type=type_diag_loc)

       deallocate(Eig)

       IF (debug) write(out_unit,*) 'init Vec0 done'

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

      write(out_unit,*) 'Eig',(PotVal_adia%d0(i,i),i=1,nsurf)

      write(out_unit,*) 'PotVal_adia'
      CALL Write_dnMat(PotVal_adia,nio=out_unit)

      write(out_unit,*) 'Vec'
      CALL Write_dnMat(Vec,nio=out_unit)

      write(out_unit,*) 'Vec0'
      CALL Write_dnMat(Vec0,nio=out_unit)

      write(out_unit,*) 'NAC'
      CALL Write_dnMat(NAC,nio=out_unit)
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE dia_TO_adia

  SUBROUTINE Eval_dnHVib_ana(Model,Qact,dnH,nderiv)
    USE QDUtil_m,  ONLY : Write_Mat
    USE ADdnSVM_m, ONLY : dnMat_t,dnS_t,alloc_dnMat,dnMat_TO_dnS,dot_product,  &
        ReduceDerivatives_dnS2_TO_dnS1,Write_dnMat,Write_dnS,dnS_TO_dnMat,SYM_dnMat
    USE AdiaChannels_Basis_m
    IMPLICIT NONE

  real (kind=Rkind),              intent(in)    :: Qact(:)
  TYPE (Model_t),                 intent(inout) :: Model
  TYPE (dnMat_t),                 intent(inout) :: dnH ! derivative of the Hamiltonian
  integer,                        intent(in)    :: nderiv


  integer                        :: ii,ji,i,iq,ib,jb,nb,nq

  TYPE (dnMat_t)                 :: PotVal
  real (kind=Rkind), allocatable :: Q(:),d0GGdef(:,:)
  integer                        :: ndim_act

  TYPE (dnS_t), allocatable      :: dnV(:),dnHB(:)
  TYPE (dnS_t)                   :: dnVfull,dnHij


  ndim_act = size(Model%QM%list_act)

  nb = Model%Basis%nb
  nq = Model%Basis%nq


  CALL alloc_dnMat(dnH,nsurf=nb,nVar=ndim_act,nderiv=nderiv,name_var='dnH')

  allocate(dnV(nq))
  allocate(dnHB(nq))

  allocate(Q(Model%QM%ndim))
  DO i=1,size(Model%QM%list_act)
    Q(Model%QM%list_act(i)) = Qact(i)
  END DO

  !Buid H
  DO iq=1,Model%Basis%nq
    DO ii=1,size(Model%QM%list_inact)
      Q(Model%QM%list_inact(ii)) = Model%Basis%x(iq)
    END DO

    CALL Eval_Pot_ana(Model,Q,PotVal,nderiv=nderiv)
    CALL dnMat_TO_dnS(PotVal,dnVfull,i=1,j=1)
    CALL ReduceDerivatives_dnS2_TO_dnS1(dnV(iq),dnVfull,Model%QM%list_act)
  END DO

  !CALL Write_dnMat(PotVal,6,info='PotVal')
  !CALL Write_dnS(dnV(nq),6,info='dnV',all_type=.TRUE.)

  d0GGdef = Model%QM%d0GGdef(Model%QM%list_inact,Model%QM%list_inact)
  DO ib=1,nb
    ! H B(:,ib)>
    DO iq=1,nq
      dnHB(iq) = -HALF*d0GGdef(1,1)*Model%Basis%d2gb(iq,ib,1,1) + &
                 dnV(iq)*Model%Basis%d0gb(iq,ib)
      dnHB(iq) = dnHB(iq) * Model%Basis%w(iq)
    END DO
    !CALL Write_dnS(dnHB(1),6,info='dnHB',all_type=.TRUE.)
    !write(out_unit,*) 'coucou dnHB: done',ib ; flush(6)
    DO jb=1,nb
      IF (Model%Basis%tab_symab(ib) == Model%Basis%tab_symab(jb)) THEN
        dnHij = dot_product(Model%Basis%d0gb(:,jb),dnHB(:))
      ELSE
        dnHij = ZERO
      END IF
      CALL dnS_TO_dnMat(dnHij,dnH,jb,ib)
    END DO
  END DO

  dnH = SYM_dnMat(dnH)

  !CALL Write_Mat(dnH%d0,6,5,name_info='H')

  END SUBROUTINE Eval_dnHVib_ana

  SUBROUTINE Eval_Func(Model,Q,Func,nderiv)
  USE ADdnSVM_m, ONLY : dnS_t,dealloc_dnS,alloc_dnS,Variable,Write_dnS
  IMPLICIT NONE

    TYPE (Model_t),                 intent(inout)            :: Model

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
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) '   nderiv    ',nderiv
      flush(out_unit)
    END IF

    CALL check_alloc_QM(Model,name_sub)

    IF (Model%QM%nb_Func > 0) THEN
      IF (allocated(Func)) THEN
        DO i=1,size(Func)
          CALL dealloc_dnS(Func(i))
        END DO
        deallocate(Func)
      END IF
      allocate(Func(Model%QM%nb_Func))
      DO i=1,size(Func)
        CALL alloc_dnS(Func(i),Model%QM%ndimFunc,nderiv)
      END DO

      allocate(dnQ(Model%QM%ndimFunc))
      DO i=1,Model%QM%ndimFunc
        dnQ(i) = Variable(Q(i),nVar=Model%QM%ndimFunc,nderiv=nderiv,iVar=i) ! to set up the derivatives
      END DO

      CALL Model%QM%EvalFunc_QModel(Func,dnQ,nderiv=nderiv)

      ! deallocation
      DO i=1,size(dnQ)
        CALL dealloc_dnS(dnQ(i))
      END DO
      deallocate(dnQ)
      ! end deallocation

    END IF

    IF (debug) THEN
      write(out_unit,*) 'Func',size(Func)
      DO i=1,size(Func)
        CALL Write_dnS(Func(i),nio=out_unit,all_type=.TRUE.)
      END DO
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Eval_Func

  SUBROUTINE Write_QML_EigenVec(Q,Vec,Model,nio)
  USE AdiaChannels_Basis_m
  USE ADdnSVM_m, ONLY : dnMat_t,Check_NotAlloc_dnMat
  IMPLICIT NONE

    real (kind=Rkind),  intent(in)              :: Q(:)
    TYPE (dnMat_t),     intent(in)              :: Vec
    TYPE (Model_t),     intent(in)              :: Model
    integer,            intent(in), optional    :: nio

    integer :: i,nio_loc
    real (kind=Rkind), allocatable :: G(:,:)

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unit
    END IF

    CALL check_alloc_QM(Model,'Write_QML_EigenVec')

    IF (.NOT. Model%QM%print_EigenVec_Basis .AND. .NOT. Model%QM%print_EigenVec_Grid) RETURN

    IF (Check_NotAlloc_dnMat(Vec,nderiv=0)) THEN
        write(nio_loc,*) '-----------------------------------------------'
        write(nio_loc,*) 'Vec%d0 cannot be printed, it is not allocated'
        write(nio_loc,*) '-----------------------------------------------'
    ELSE

      IF (Model%QM%print_EigenVec_Basis) THEN
        write(nio_loc,*) '-----------------------------------------------'
        DO i=1,Model%Basis%nb
          write(nio_loc,*) 'wfb',Q,i,Vec%d0(i,1:Model%QM%nb_Channels)
        END DO
      END IF

      IF (Model%QM%print_EigenVec_Grid .AND. Model%QM%Vib_adia) THEN
        write(nio_loc,*) '-----------------------------------------------'
        allocate(G(Model%Basis%nq,Model%QM%nb_Channels))
        DO i=1,Model%QM%nb_Channels
          CALL BasisTOGrid_Basis(G(:,i),Vec%d0(:,i),Model%Basis)
        END DO
        DO i=1,Model%Basis%nq
          write(nio_loc,*) 'wfg',Q,Model%Basis%x(i),G(i,:)
        END DO

        deallocate(G)

      END IF
      write(nio_loc,*) '-----------------------------------------------'
      flush(nio_loc)
    END IF

  END SUBROUTINE Write_QML_EigenVec

  SUBROUTINE Write_Model(Model,nio)
    USE QDUtil_m,         ONLY : Write_Mat, Write_Vec
    IMPLICIT NONE

    TYPE(Model_t),      intent(in)              :: Model
    integer,            intent(in), optional    :: nio

    integer :: nio_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unit
    END IF

    CALL check_alloc_QM(Model,'Write_Model')

    IF (nio_loc /= out_unit) THEN
      open(nio_loc,file=trim(adjustl(Model%QM%pot_name))//'.out',form='formatted')
    END IF


    write(nio_loc,*) '-----------------------------------------------'
    write(nio_loc,*) 'Output file for potential library'

    CALL Model%QM%Write_QModel(nio=nio_loc)
    write(nio_loc,*)
    IF (allocated(Model%QM%d0GGdef)) THEN 
      CALL Write_Mat(Model%QM%d0GGdef,nio_loc,5,info='d0GGdef')
    END IF
    write(nio_loc,*)
    IF (allocated(Model%QM%Q0)) THEN
      CALL Write_Vec(Model%QM%Q0,nio_loc,5,info='Q0')
    END IF
    write(nio_loc,*)
    write(nio_loc,*) '-----------------------------------------------'
    write(nio_loc,*) 'Extra action(s):'
    write(nio_loc,*) 'opt',Model%opt
    write(nio_loc,*) 'irc',Model%irc
    write(nio_loc,*) '-----------------------------------------------'
    flush(nio_loc)

  END SUBROUTINE Write_Model
  SUBROUTINE Write_QdnV_FOR_Model(Q,PotVal,Model,Vec,NAC,info,name_file)
    USE QDUtil_m, ONLY : file_open2
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (Model_t),    intent(in)           :: Model
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
      CALL file_open2(trim(adjustl(Model%QM%pot_name))//'.txt',                &
                      nio_loc,lformatted=.TRUE.,append=.TRUE.,err_file=err_io)
    END IF

    IF (err_io /= 0) THEN
      write(out_unit,*) 'ERROR in Write_QdnV_FOR_Model'
      write(out_unit,*) ' Impossible to open the file "',                      &
                          trim(adjustl(Model%QM%pot_name))//'.txt','"'
      STOP 'Impossible to open the file'
    END IF


    write(nio_loc,'(a)',advance='no') 'TEST output: '
    IF (present(info)) THEN
      write(nio_loc,'(a)',advance='no') info
    END IF
    IF (Model%QM%adiabatic) THEN
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

    IF (allocated(Model%QM%d0GGdef)) THEN
      write(nio_loc,*) 'd0GGdef'
      write(nio_loc,*) size(Model%QM%d0GGdef)
      write(nio_loc,*) Model%QM%d0GGdef
    END IF

    write(nio_loc,'(a)',advance='no') 'END_TEST output: '
    IF (present(info)) THEN
      write(nio_loc,'(a)',advance='no') info
    END IF
    IF (Model%QM%adiabatic) THEN
      write(nio_loc,'(a)') ' Adiabatic'
    ELSE
      write(nio_loc,'(a)') ' Diabatic'
    END IF

    close(nio_loc)

  END SUBROUTINE Write_QdnV_FOR_Model
  SUBROUTINE Test_QdnV_FOR_Model(Q,PotVal,Model,Vec,NAC,info,name_file,test_file_path,test_var,last_test)
    USE QDUtil_m, ONLY : file_open2, TO_string, Write_Mat
    USE QDUtil_Test_m
    USE QMLLib_UtilLib_m
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (Model_t),    intent(in)              :: Model
    TYPE (dnMat_t),    intent(in)              :: PotVal
    real (kind=Rkind), intent(in)              :: Q(:)
    TYPE (dnMat_t),    intent(in),    optional :: Vec ! for non adiabatic couplings
    TYPE (dnMat_t),    intent(in),    optional :: NAC ! for non adiabatic couplings
    character(len=*),  intent(in),    optional :: info
    character(len=*),  intent(in),    optional :: name_file,test_file_path
    TYPE (test_t),     intent(inout)           :: test_var
    logical,           intent(in),    optional :: last_test

    integer :: write_file_unit,read_file_unit,err_io,nsize
    logical :: read_file_open,read_file_exist,last_test_loc
    character(len=:),  allocatable :: write_file_name,read_file_name,test_file_path_loc
    real (kind=Rkind), allocatable :: Qref(:)
    TYPE (dnMat_t)                 :: PotValref,Vecref
    real (kind=Rkind), allocatable :: d1NAC(:)
    real (kind=Rkind), allocatable :: Gdef(:)
    character(len=:),  allocatable :: info_loc

    real (kind=Rkind), parameter   :: ZeroTresh    = ONETENTH**9
    logical :: res_test
    character(len=Name_longlen) ::dum_name

    ! open file to write
    IF (present(test_file_path)) THEN
      test_file_path_loc = trim(adjustl(test_file_path))
    ELSE
      test_file_path_loc = 'RES_files/'
    END IF
    IF (present(name_file)) THEN
      write_file_name = trim(adjustl(name_file))
    ELSE
      write_file_name = trim(adjustl(Model%QM%pot_name)) // '.txt'
    END IF
    CALL file_open2((test_file_path_loc // write_file_name),write_file_unit,lformatted=.TRUE.,append=.TRUE.,err_file=err_io)
    IF (err_io /= 0) THEN
      write(out_unit,*) 'ERROR in Test_QdnV_FOR_Model'
      write(out_unit,*) ' Impossible to open the file: ',(test_file_path_loc // write_file_name)
      STOP 'ERROR in Test_QdnV_FOR_Model: Impossible to open the file'
    END IF
    
    ! check if the file to be read is open, then get the unit. Otherwise, open it
    read_file_name = 'RES_old/' // write_file_name
    inquire(FILE=read_file_name, OPENED=read_file_open, EXIST=read_file_exist, NUMBER=read_file_unit)
    IF (read_file_exist .AND. .NOT. read_file_open) THEN
      CALL file_open2(read_file_name,read_file_unit,lformatted=.TRUE.,err_file=err_io)
      IF (err_io /= 0) THEN
        write(out_unit,*) 'ERROR in Test_QdnV_FOR_Model'
        write(out_unit,*) ' Impossible to open the file: ',read_file_name
        STOP 'ERROR in Test_QdnV_FOR_Model: Impossible to open the file'
      END IF
      read_file_open = .TRUE.
    END IF
    write(out_unit,*) ' Read old file: ',read_file_name, ', open?',read_file_open
    IF (present(last_test)) THEN
      last_test_loc = last_test
    ELSE
      last_test_loc = .FALSE.
    END IF

    IF (present(info)) THEN
      info_loc = info
    ELSE
      info_loc = ''
    END IF

    ! write, read and test
    IF (Model%QM%adiabatic) THEN
      write(write_file_unit,*) 'TEST output: ',info_loc,' Adiabatic'
    ELSE
      write(write_file_unit,*) 'TEST output: ',info_loc,' Diabatic'
    END IF

    IF (read_file_open) THEN
      read(read_file_unit,*,IOSTAT=err_io)
      IF (err_io /=0) THEN
        close(read_file_unit)
        read_file_open = .FALSE.
        CALL Logical_Test(test_var,test1=res_test,info=info_loc // ' The old file is incomplete !' )
        CALL Finalize_Intermediate_Test(test_var,info=info_loc)
      END IF
    END IF

    write(write_file_unit,*) 'Q'
    write(write_file_unit,*) size(Q)
    write(write_file_unit,*) Q

    IF (read_file_open) THEN
      Qref = Read_alloc_Vect(read_file_unit,err_io)
      ! For testing the model, Q, PotVal, G
      res_test = (err_io == 0)
      IF (res_test) res_test = ( size(Q) == size(Qref) )
      IF (res_test) res_test = all(abs(Q-Qref) < ZeroTresh)
      CALL Logical_Test(test_var,test1=res_test,info=info_loc // ': Q == Qref:   T ? ' // TO_string(res_test) )
      IF (.NOT. res_test) THEN
        IF (err_io /= 0) write(out_unit,*) 'Problem while reading the old file: ',read_file_name
        write(out_unit,*) info_loc // ', Q:   ',Q
        IF (err_io == 0) write(out_unit,*) info_loc // ', Qref:',Qref
        IF (err_io == 0 .AND. size(Q) == size(Qref) ) write(out_unit,*) info_loc // ', Diff:',Q-Qref
      END IF
    END IF

    IF (allocated(PotVal%d0)) THEN
      write(write_file_unit,*) 'V'
      write(write_file_unit,*) size(PotVal%d0)
      write(write_file_unit,*) PotVal%d0
    END IF
    IF (allocated(PotVal%d1)) THEN
      write(write_file_unit,*) 'Grad'
      write(write_file_unit,*) size(PotVal%d1)
      write(write_file_unit,*) PotVal%d1
    END IF
    IF (allocated(PotVal%d2)) THEN
      write(write_file_unit,*) 'Hess'
      write(write_file_unit,*) size(PotVal%d2)
      write(write_file_unit,*) PotVal%d2
    END IF

    PotValref = PotVal
    IF (read_file_open .AND. allocated(PotVal%d0)) THEN
      PotValref%d0 = reshape(Read_alloc_Vect(read_file_unit,err_io),shape=shape(PotValref%d0))
    END IF
    IF (read_file_open .AND. allocated(PotVal%d1)) THEN
      PotValref%d1 = reshape(Read_alloc_Vect(read_file_unit,err_io),shape=shape(PotValref%d1))
    END IF
    IF (read_file_open .AND. allocated(PotVal%d2)) THEN
      PotValref%d2 = reshape(Read_alloc_Vect(read_file_unit,err_io),shape=shape(PotValref%d2))
    END IF
    IF (read_file_open) THEN
      res_test = (err_io == 0)
      IF (res_test) res_test = all(abs(get_Flatten(PotVal)-get_Flatten(PotValref)) < ZeroTresh)
      CALL Logical_Test(test_var,test1=res_test,info=info_loc // ': PotVal == PotValref:   T ? ' // TO_string(res_test) )
      IF (.NOT. res_test) THEN
        IF (err_io /= 0) write(out_unit,*) 'Problem while reading the old file: ',read_file_name
        CALL Write_dnMat(PotVal,          nio=out_unit,info=info_loc // ', PotVal')
        IF (err_io == 0) CALL Write_dnMat(PotValref,       nio=out_unit,info=info_loc // ', PotValref')
        IF (err_io == 0) CALL Write_dnMat(PotVal-PotValref,nio=out_unit,info=info_loc // ', diff')
      END IF
    END IF
  

    IF (present(Vec)) THEN
      IF (allocated(Vec%d0)) THEN
        write(write_file_unit,*) 'Vec'
        write(write_file_unit,*) size(Vec%d0)
        write(write_file_unit,*) Vec%d0
      END IF
      IF (allocated(Vec%d1)) THEN
        write(write_file_unit,*) 'd1Vec'
        write(write_file_unit,*) size(Vec%d1)
        write(write_file_unit,*) Vec%d1
      END IF
      IF (allocated(Vec%d2)) THEN
        write(write_file_unit,*) 'd2Vec'
        write(write_file_unit,*) size(Vec%d2)
        write(write_file_unit,*) Vec%d2
      END IF

      Vecref = Vec
      IF (read_file_open .AND. allocated(Vec%d0)) THEN
        Vecref%d0 = reshape(Read_alloc_Vect(read_file_unit,err_io),shape=shape(Vecref%d0))
      END IF
      IF (read_file_open .AND. allocated(Vec%d1)) THEN
        Vecref%d1 = reshape(Read_alloc_Vect(read_file_unit,err_io),shape=shape(Vecref%d1))
      END IF
      IF (read_file_open .AND. allocated(Vec%d2)) THEN
        Vecref%d2 = reshape(Read_alloc_Vect(read_file_unit,err_io),shape=shape(Vecref%d2))
      END IF
      IF (read_file_open) THEN
        res_test = (err_io == 0)
        IF (res_test) res_test = all(abs(get_Flatten(Vec)-get_Flatten(Vecref)) < ZeroTresh)
        CALL Logical_Test(test_var,test1=res_test,info=info_loc // ': Vec == Vecref:   T ? ' // TO_string(res_test) )
        IF (.NOT. res_test) THEN
          IF (err_io /= 0) write(out_unit,*) 'Problem while reading the old file: ',read_file_name
          CALL Write_dnMat(Vec,       nio=out_unit,info=info_loc // ', Vec')
          IF (err_io == 0) CALL Write_dnMat(Vecref,    nio=out_unit,info=info_loc // ', Vecref')
          IF (err_io == 0) CALL Write_dnMat(Vec-Vecref,nio=out_unit,info=info_loc // ', diff')
        END IF
      END IF

    END IF

    IF (present(NAC)) THEN
      IF (allocated(NAC%d1)) THEN
        write(write_file_unit,*) 'NAC'
        write(write_file_unit,*) size(NAC%d1)
        write(write_file_unit,*) NAC%d1
      END IF
      IF (read_file_open .AND. allocated(NAC%d1)) THEN
        d1NAC = Read_alloc_Vect(read_file_unit,err_io)
        res_test = (err_io == 0)
        IF (res_test) res_test = all(abs(reshape(NAC%d1,shape=[size(NAC%d1)])-d1NAC) < ZeroTresh)
        CALL Logical_Test(test_var,test1=res_test,info=info_loc // ': NAC == NACref:   T ? ' // TO_string(res_test) )
        IF (.NOT. res_test) THEN
          IF (err_io /= 0) write(out_unit,*) 'Problem while reading the old file: ',read_file_name
          write(out_unit,*) info_loc // ', NAC:    ',NAC%d1
          IF (err_io == 0) write(out_unit,*) info_loc // ', QNACref:',d1NAC
          IF (err_io == 0) write(out_unit,*) info_loc // ', diff:   ',reshape(NAC%d1,shape=[size(NAC%d1)])-d1NAC
        END IF
      END IF
    END IF

    IF (allocated(Model%QM%d0GGdef)) THEN
      write(write_file_unit,*) 'd0GGdef'
      write(write_file_unit,*) size(Model%QM%d0GGdef)
      write(write_file_unit,*) Model%QM%d0GGdef
    END IF
    IF (read_file_open .AND. allocated(Model%QM%d0GGdef)) THEN
      Gdef = Read_alloc_Vect(read_file_unit,err_io)
      res_test = (err_io == 0)
      IF (res_test) res_test = all(abs(reshape(Model%QM%d0GGdef,shape=[size(Model%QM%d0GGdef)])-Gdef) < ZeroTresh)
      CALL Logical_Test(test_var,test1=res_test,info=info_loc // ': d0GGdef == Gref:   T ? ' // TO_string(res_test) )
      IF (.NOT. res_test) THEN
        IF (err_io /= 0) write(out_unit,*) 'Problem while reading the old file: ',read_file_name
        CALL Write_Mat(Model%QM%d0GGdef,                            &
            nio=out_unit,nbcol=5,info=info_loc // ', d0GGdef')
        IF (err_io == 0) CALL Write_Mat(reshape(Gdef,shape=shape(Model%QM%d0GGdef)), &
            nio=out_unit,nbcol=5,info=info_loc // ', Gref')
        IF (err_io == 0) CALL Write_Mat(Model%QM%d0GGdef-reshape(Gdef,shape=shape(Model%QM%d0GGdef)), &
            nio=out_unit,nbcol=5,info=info_loc // ', diff')
      END IF
    END IF

    write(write_file_unit,'(a)',advance='no') 'END_TEST output: '
    IF (read_file_open) THEN
      DO
        read(read_file_unit,*) dum_name
        IF (dum_name == 'END_TEST') EXIT
      END DO
    END IF

    close(write_file_unit)
    IF (last_test_loc .AND. read_file_open) THEN
      close(read_file_unit)
      CALL Finalize_Intermediate_Test(test_var,info=info_loc)
    END IF

  END SUBROUTINE Test_QdnV_FOR_Model
  SUBROUTINE Check_analytical_numerical_derivatives(Model,Q,nderiv,test_var,AnaNum_Test)
    USE QDUtil_m,  ONLY : TO_string
    USE QDUtil_Test_m
    USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat,Write_dnMat,get_maxval_OF_dnMat, operator(-)
    IMPLICIT NONE

    TYPE (Model_t),       intent(inout)           :: Model
    real (kind=Rkind),    intent(in)              :: Q(:)
    integer,              intent(in)              :: nderiv
    TYPE (test_t),        intent(inout), optional :: test_var
    logical,              intent(inout), optional :: AnaNum_Test

    ! local variables
    TYPE (dnMat_t)            :: Mat_diff
    TYPE (dnMat_t)            :: PotVal_ana,PotVal_num
    TYPE (dnMat_t)            :: NAC_ana,NAC_num
    TYPE (dnMat_t)            :: Vec_ana,Vec_num

    real (kind=Rkind)         :: MaxMat,MaxDiffMat
    logical                   :: AnaNum_Test_loc

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Check_analytical_numerical_derivatives'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (Model%QM%no_ana_der) RETURN

    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) '   nderiv    ',nderiv
      flush(out_unit)
    END IF

    CALL check_alloc_QM(Model,name_sub)

    CALL alloc_dnMat(PotVal_ana,nsurf=Model%QM%nsurf,              &
                         nVar=Model%QM%ndim,nderiv=nderiv)

    CALL alloc_dnMat(PotVal_num,nsurf=Model%QM%nsurf,              &
                         nVar=Model%QM%ndim,nderiv=nderiv)

    IF (Model%QM%adiabatic .AND. Model%QM%nsurf > 1) THEN
      CALL Eval_Pot(Model,Q,PotVal_ana,nderiv,NAC_ana,Vec_ana,numeric=.FALSE.)
    ELSE
      CALL Eval_Pot(Model,Q,PotVal_ana,nderiv,numeric=.FALSE.)
    END IF

    IF (debug) THEN
      write(out_unit,*)   'PotVal_ana'
      CALL Write_dnMat(PotVal_ana,nio=out_unit)
      flush(out_unit)
    END IF

    IF (Model%QM%adiabatic .AND. Model%QM%nsurf > 1) THEN
      CALL Eval_Pot(Model,Q,PotVal_num,nderiv,NAC_num,Vec_num,numeric=.TRUE.)
    ELSE
      CALL Eval_Pot(Model,Q,PotVal_num,nderiv,numeric=.TRUE.)
    END IF
    IF (debug) THEN
      write(out_unit,*)   'PotVal_num'
      CALL Write_dnMat(PotVal_num,nio=out_unit)
      flush(out_unit)
    END IF


    MaxMat      = get_maxval_OF_dnMat(PotVal_ana)
    IF (MaxMat < ONETENTH**6) MaxMat = ONE
    Mat_diff    = PotVal_num - PotVal_ana
    MaxDiffMat  = get_maxval_OF_dnMat(Mat_diff)

    write(out_unit,'(3a,e9.2)') 'With ',Model%QM%pot_name,                    &
               ': max of the relative Potential diff:',MaxDiffMat/MaxMat
    write(out_unit,'(3a,l9)')   'With ',Model%QM%pot_name,                    &
     ': Potential diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

    AnaNum_Test_loc = (MaxDiffMat/MaxMat <= step)

     IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
      write(out_unit,*)   'Potential diff (ana-numer)'
      CALL Write_dnMat(Mat_diff,nio=out_unit)
    END IF

    IF (Model%QM%adiabatic .AND. Model%QM%nsurf > 1) THEN

      MaxMat      = get_maxval_OF_dnMat(NAC_ana)
      IF (MaxMat < ONETENTH**6) MaxMat = ONE
      Mat_diff    = NAC_num - NAC_ana
      MaxDiffMat  = get_maxval_OF_dnMat(Mat_diff)

      write(out_unit,'(3a,e9.2)') 'With ',Model%QM%pot_name,                  &
                 ': max of the relative NAC diff:',MaxDiffMat/MaxMat
      write(out_unit,'(3a,l9)')   'With ',Model%QM%pot_name,                  &
       ': NAC diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

      AnaNum_Test_loc = AnaNum_Test_loc .AND. (MaxDiffMat/MaxMat <= step)


      IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
        write(out_unit,*)   'NAC diff (ana-numer)'
        CALL Write_dnMat(Mat_diff,nio=out_unit)
        write(out_unit,*)   'NAC_ana'
        CALL Write_dnMat(NAC_ana,nio=out_unit)
        write(out_unit,*)   'NAC_num'
        CALL Write_dnMat(NAC_num,nio=out_unit)
      END IF

      MaxMat      = get_maxval_OF_dnMat(Vec_ana)
      IF (MaxMat < ONETENTH**6) MaxMat = ONE
      Mat_diff    = Vec_num - Vec_ana
      MaxDiffMat  = get_maxval_OF_dnMat(Mat_diff)

      write(out_unit,'(3a,e9.2)') 'With ',Model%QM%pot_name,            &
                 ': max of the relative Vec diff:',MaxDiffMat/MaxMat
      write(out_unit,'(3a,l9)')   'With ',Model%QM%pot_name,            &
       ': Vec diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

      AnaNum_Test_loc = AnaNum_Test_loc .AND. (MaxDiffMat/MaxMat <= step)

      IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
        write(out_unit,*)   'Vec diff (ana-numer)'
        CALL Write_dnMat(Mat_diff,nio=out_unit)
        write(out_unit,*)   'Vec_ana'
        CALL Write_dnMat(Vec_ana,nio=out_unit)
        write(out_unit,*)   'Vec_num'
        CALL Write_dnMat(Vec_num,nio=out_unit)
      END IF

    END IF

    CALL dealloc_dnMat(PotVal_ana)
    CALL dealloc_dnMat(PotVal_num)
    CALL dealloc_dnMat(NAC_ana)
    CALL dealloc_dnMat(NAC_num)
    CALL dealloc_dnMat(Vec_ana)
    CALL dealloc_dnMat(Vec_num)
    CALL dealloc_dnMat(Mat_diff)

    IF (present(test_var)) THEN
      CALL Logical_Test(test_var,test1=AnaNum_Test_loc, &
        info=Model%QM%pot_name // ': ana == numer:   T ? ' // TO_string(AnaNum_Test_loc) )
    END IF

    IF (present(AnaNum_Test)) AnaNum_Test = AnaNum_Test_loc

    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Check_analytical_numerical_derivatives
  SUBROUTINE Test_QVG_FOR_Model(Model,Q,test_var,nderiv,option)
    USE QDUtil_m, ONLY : file_open2, TO_string, Write_Mat
    USE QDUtil_Test_m
    USE QMLLib_UtilLib_m
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (Model_t),    intent(inout)           :: Model
    real (kind=Rkind), intent(in)              :: Q(:)
    TYPE (test_t),     intent(inout)           :: test_var
    integer,           intent(in)              :: nderiv,option


    real (kind=Rkind), allocatable :: Qref(:)
    real (kind=Rkind), allocatable :: G(:,:),Gref(:,:) ! metric tensor
    integer                        :: ndim,nsurf,i,err
    logical                        :: Lerr
    TYPE (dnMat_t)                 :: PotVal
    TYPE (dnMat_t)                 :: PotValref
    TYPE (dnMat_t)                 :: dnErr

    character(len=:), allocatable  :: fmt


    Lerr = check_Init_QModel(Model)
    IF (.NOT. Lerr) THEN
      write(out_unit,*) 'The model in not initialized'
      CALL Logical_Test(test_var,test1=Lerr,info='The model in not initialized')
      RETURN
    END IF


    CALL Eval_Pot(Model,Q,PotVal,nderiv=nderiv)
    fmt = '(a,' // TO_string(Model%ndim) // 'f12.6)'
    write(out_unit,fmt) 'Q',Q(:)
    write(out_unit,*) 'Energy'
    CALL Write_dnMat(PotVal,nio=out_unit)
    flush(out_unit)

    ! For testing the model
    allocate(Qref(Model%ndim))
    allocate(Gref(Model%ndim,Model%ndim))
    CALL Model%QM%RefValues_QModel(err,nderiv=nderiv,Q0=Qref,d0GGdef=Gref,dnMatV=PotValref,option=option)

    write(out_unit,*) 'Reference Energy' ; flush(out_unit)
    CALL Write_dnMat(PotValref,nio=out_unit)
    flush(out_unit)
    dnErr = PotValref-PotVal
    Lerr  = Check_dnMat_IS_ZERO(dnErr)

    CALL Logical_Test(test_var,test1=Lerr,info='dnMatV')
    IF (.NOT. Lerr) CALL Write_dnMat(dnErr,nio=out_unit,info='dnErr')

    Lerr = all(abs(Q-Qref) < epsi)
    CALL Logical_Test(test_var,test1=Lerr,info='Q(:)')
    IF (.NOT. Lerr) Write(out_unit,*) 'Q-Qref',Q-Qref

    G = get_d0GGdef_Model(Model=Model)
    Lerr = all(abs(G-Gref) < epsi)
    CALL Logical_Test(test_var,test1=Lerr,info='G (metrix tensor)')
    IF (.NOT. Lerr) Write(out_unit,*) 'G-Gref',G-Gref

    deallocate(Qref)
    deallocate(Gref)

  END SUBROUTINE Test_QVG_FOR_Model
  SUBROUTINE Eval_pot_ON_Grid(Model,Qmin,Qmax,nb_points,nderiv,grid_file)
    USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),               intent(inout)   :: Model
    real (kind=Rkind),            intent(in)      :: Qmin(:),Qmax(:)
    integer, optional,            intent(in)      :: nb_points,nderiv
    character (len=*), optional,  intent(in)      :: grid_file


    ! local variables
    integer                        :: unit_grid_file
    integer                        :: i,iq,jq,i1,i2,nb_points_loc,nderiv_loc,ndim_loc
    integer,           allocatable :: i_Q(:)
    real (kind=Rkind), allocatable :: dQ(:),Q(:)
    TYPE (dnMat_t)                 :: PotVal,NAC


    CALL check_alloc_QM(Model,'Eval_pot_ON_Grid')


    IF (size(Qmin) /= Model%QM%ndim .OR. size(Qmax) /= Model%QM%ndim) THEN
       write(out_unit,*) ' ERROR in Eval_pot_ON_Grid'
       write(out_unit,*) ' The size of Qmin or Qmax are different from Model%QM%ndim'
       write(out_unit,*) '  size(Qmin)    ',size(Qmin)
       write(out_unit,*) '  size(Qmax)    ',size(Qmax)
       write(out_unit,*) '  Model%QM%ndim',Model%QM%ndim
       write(out_unit,*) ' => Check the fortran'
       STOP 'ERROR in Eval_pot_ON_Grid: problem with Model%QM%ndim'
    END IF

    IF (present(grid_file)) THEN
      IF (len_trim(grid_file) == 0) THEN
        unit_grid_file = out_unit
      ELSE
        unit_grid_file = 99 ! not used
        open(newunit=unit_grid_file,file=trim(grid_file) )
      END IF
    ELSE
      unit_grid_file = out_unit
    END IF

    nb_points_loc = 100
    IF (present(nb_points)) nb_points_loc = nb_points
    nb_points_loc = max(nb_points_loc,2)

    IF (present(nderiv)) THEN
      nderiv_loc = nderiv
    ELSE
      nderiv_loc = 0
    END IF

    allocate(dQ(Model%QM%ndim))
    allocate(Q(Model%QM%ndim))
    allocate(i_Q(Model%QM%ndim))

    dQ(:)       = (Qmax-Qmin) / real(nb_points_loc-1,kind=Rkind)
    ndim_loc    = 0
    i_Q(:)      = 0
    DO i=1,Model%QM%ndim
      IF (dQ(i) /= ZERO) THEN
        ndim_loc = ndim_loc + 1
        i_Q(ndim_loc) = i
      END IF
    END DO
    write(out_unit,*) 'Grid. File name: "',trim(grid_file),'"'
    !write(out_unit,*) 'Coordinates indices, i_Q: ',i_Q(1:ndim_loc)
    !write(out_unit,*) 'Model%QM%ndim',Model%QM%ndim
    !write(out_unit,*) 'Model%QM%numeric',Model%QM%numeric
    !write(out_unit,*) 'ndim for the grid',ndim_loc


    Q(:) = Qmin

    IF (ndim_loc == 1) THEN
      i1 = i_Q(1)
      DO iq=0,nb_points_loc-1
        Q(i1) = Qmin(i1) + dQ(i1)*real(iq,kind=Rkind)
        IF (Model%QM%nsurf > 1 .AND. Model%QM%adiabatic) THEN
          CALL Eval_Pot(Model,Q,PotVal,nderiv=max(1,nderiv_loc),NAC=NAC)

          IF (nderiv_loc == 0) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,Model%QM%nsurf),NAC%d1
          ELSE IF (nderiv_loc == 1) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,Model%QM%nsurf),(PotVal%d1(i,i,:),i=1,Model%QM%nsurf)
          ELSE
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,Model%QM%nsurf),(PotVal%d1(i,i,:),i=1,Model%QM%nsurf), &
                            (PotVal%d2(i,i,:,:),i=1,Model%QM%nsurf)
          END IF
          flush(unit_grid_file)

        ELSE
          CALL Eval_Pot(Model,Q,PotVal,nderiv=nderiv_loc)

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

        IF (Model%QM%nsurf > 1 .AND. Model%QM%adiabatic) THEN
          CALL Eval_Pot(Model,Q,PotVal,nderiv=0,NAC=NAC)
          write(unit_grid_file,*) Q(i1),Q(i2),(PotVal%d0(i,i),i=1,Model%QM%nsurf)
        ELSE
          CALL Eval_Pot(Model,Q,PotVal,nderiv=nderiv_loc)
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

    IF (unit_grid_file /= out_unit) THEN
      close(unit_grid_file)
    END IF


  END SUBROUTINE Eval_pot_ON_Grid


  SUBROUTINE calc_pot(V,Model,Q)
    USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),         intent(inout)   :: Model
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated

    TYPE (dnMat_t)                         :: PotVal

    CALL check_alloc_QM(Model,'calc_pot')

    CALL Eval_Pot(Model,Q,PotVal,nderiv=0)

    V = PotVal%d0

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot
  SUBROUTINE calc_pot_grad(V,g,Model,Q)
    USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),         intent(inout)   :: Model
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(Model,'calc_pot_grad')

    CALL Eval_Pot(Model,Q,PotVal,nderiv=1)

    V = PotVal%d0
    g = PotVal%d1

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot_grad
  SUBROUTINE calc_grad(g,Model,Q)
    USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),         intent(inout)   :: Model
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(Model,'calc_grad')

    CALL Eval_Pot(Model,Q,PotVal,nderiv=1)

    g = PotVal%d1

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_grad
  SUBROUTINE calc_pot_grad_hess(V,g,h,Model,Q)
    USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),         intent(inout)   :: Model
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(Model,'calc_pot_grad_hess')

    CALL Eval_Pot(Model,Q,PotVal,nderiv=2)

    V = PotVal%d0
    g = PotVal%d1
    h = PotVal%d2

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot_grad_hess
  SUBROUTINE calc_hess(h,Model,Q)
    USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
    IMPLICIT NONE

    TYPE (Model_t),     intent(inout)     :: Model
    real (kind=Rkind),  intent(in)        :: Q(:)
    real (kind=Rkind),  intent(inout)     :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(Model,'calc_hess')

    CALL Eval_Pot(Model,Q,PotVal,nderiv=2)

    h = PotVal%d2

    CALL dealloc_dnMat(PotVal)


  END SUBROUTINE calc_hess
END MODULE Model_m
