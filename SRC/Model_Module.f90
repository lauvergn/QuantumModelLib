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
!    Copyright 2016 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================
MODULE mod_Model
!$ USE omp_lib
  USE mod_NumParameters
  USE mod_dnMat, ONLY: dnMat_t,QML_alloc_dnMat,QML_dealloc_dnMat,       &
                       QML_Check_NotAlloc_dnMat,QML_get_maxval_OF_dnMat,&
                       QML_Write_dnMat,QML_get_nsurf_FROM_dnMat,        &
                       QML_get_ndim_FROM_dnMat,QML_dnMat2_minus_dnMat1, &
                       QML_sub_dnS_TO_dnMat
  USE mod_MorsePot
  USE mod_HenonHeilesPot
  USE mod_TullyPot

  USE mod_PSB3Pot
  USE mod_HONO
  USE mod_HNNHp
  USE mod_H2SiN
  USE mod_H2NSi

  USE mod_OneDSOC_Model
  USE mod_OneDSOC_2S1T_Model

  USE mod_LinearHBondPot
  USE mod_BuckPot
  USE mod_PhenolPot
  USE mod_SigmoidPot
  USE mod_TwoDPot

  USE mod_TemplatePot

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: QModel_t,Init_Model,Eval_Pot,Write0_Model,Write_Model,Write_QdnV_FOR_Model
  PUBLIC :: calc_pot,calc_grad,calc_hess,calc_pot_grad,calc_pot_grad_hess
  PUBLIC :: Check_analytical_numerical_derivatives
  PUBLIC :: Eval_pot_ON_Grid,get_Q0_Model

  TYPE QModel_t
    integer :: nsurf       = 1
    integer :: ndim        = 1
    logical :: numeric     = .FALSE.
    logical :: adiabatic   = .TRUE.
    logical :: PubliUnit   = .FALSE. ! when PubliUnit=.TRUE., the units of a reference (publi ...) are used. Default (atomic unit)

    character (len=:), allocatable :: pot_name

    real (kind=Rkind), allocatable :: d0GGdef(:,:)

    TYPE (dnMat_t)          :: Vec0 ! to get the correct phase of the adiatic couplings

    ! list of potentials ....
    TYPE (MorsePot_t)       :: MorsePot
    TYPE (BuckPot_t)        :: BuckPot
    TYPE (SigmoidPot_t)     :: SigmoidPot

    TYPE (LinearHBondPot_t) :: LinearHBondPot
    TYPE (HenonHeilesPot_t) :: HenonHeilesPot
    TYPE (TullyPot_t)       :: TullyPot
    TYPE (OneDSOCPot_t)       :: OneDSOCPot
    TYPE (OneDSOC_2S1T_Pot_t)  :: OneDSOC_2S1T_Pot

    TYPE (PhenolPot_t)      :: PhenolPot
    TYPE (TwoDPot_t)        :: Para_TwoD
    TYPE (PSB3Pot_t)        :: PSB3Pot

    TYPE (HONOPot_t)        :: HONOPot
    TYPE (HNNHpPot_t)       :: HNNHpPot
    TYPE (H2SiNPot_t)       :: H2SiNPot
    TYPE (H2NSiPot_t)       :: H2NSiPot

    TYPE (TemplatePot_t)    :: TemplatePot

  END TYPE QModel_t

  real (kind=Rkind)                     :: step = ONETENTH**4
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

  TYPE(QModel_t), PUBLIC  :: QuantumModel

CONTAINS


  SUBROUTINE Read_Model(QModel,nio,option1,read_nml1)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (QModel_t), intent(inout) :: QModel
    integer,            intent(in)    :: nio
    integer,            intent(inout) :: option1
    logical,            intent(inout) :: read_nml1

    integer :: ndim,nsurf,nderiv,option
    logical :: adiabatic,numeric,PubliUnit,read_nml
    character (len=20) :: pot_name
    integer :: err_read

    ! Namelists for input file
    namelist /potential/ ndim,nsurf,pot_name,numeric,adiabatic,option,PubliUnit,read_nml


    ! Default values defined
    ndim      = QModel%ndim
    nsurf     = QModel%nsurf
    adiabatic = QModel%adiabatic
    pot_name  = 'morse'
    numeric   = .FALSE.
    PubliUnit = .FALSE.
    read_nml  = .TRUE. ! if T, read the namelist in PotLib (HenonHeiles ....)
    option    = -1 ! no option


    write(out_unitp,*) 'Reading input file . . .'
    read(nio,nml=potential,IOSTAT=err_read)

    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_Model'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "potential" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_Model'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_Model'
      write(out_unitp,*) ' Some parameter names of the namelist "potential" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=potential)
      STOP ' ERROR in Read_Model'
    END IF

    !write(out_unitp,nml=potential)

    option1              = option
    read_nml1            = read_nml
    QModel%ndim      = ndim
    QModel%nsurf     = nsurf
    QModel%adiabatic = adiabatic
    QModel%numeric   = numeric
    QModel%pot_name  = strdup(pot_name)
    QModel%PubliUnit = PubliUnit

  END SUBROUTINE Read_Model
 
  SUBROUTINE Init_IdMat(Mat,ndim)
  IMPLICIT NONE

  integer,                        intent(in)    :: ndim
  real (kind=Rkind), allocatable, intent(inout) :: mat(:,:)

  integer :: i

    IF (allocated(mat)) deallocate(mat)

    allocate(mat(ndim,ndim))
    mat(:,:) = ZERO
    DO i=1,ndim
      mat(i,i) = ONE
    END DO

  END SUBROUTINE Init_IdMat


  SUBROUTINE Init_Model(QModel,pot_name,ndim,nsurf,adiabatic,       &
                        read_param,param_file_name,nio_param_file,      &
                        option,PubliUnit)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (QModel_t),  intent(inout)        :: QModel
    character (len=*),   intent(in), optional :: pot_name
    integer,             intent(in), optional :: ndim,nsurf
    logical,             intent(in), optional :: adiabatic

    logical,             intent(in), optional :: read_param
    integer,             intent(in), optional :: nio_param_file
    character (len=*),   intent(in), optional :: param_file_name

    integer,             intent(in), optional :: option
    logical,             intent(in), optional :: PubliUnit

    integer ::i, option_loc,nio_loc
    logical :: read_param_loc,adiabatic_loc,read_nml
    character (len=:), allocatable :: param_file_name_loc
    real (kind=Rkind), allocatable :: Q0(:)

    write(out_unitp,*) '================================================='
    write(out_unitp,*) '================================================='
    write(out_unitp,*) '== QML: Quantum Model Lib (E-CAM) ==============='
    write(out_unitp,*) '== QML version: ',QML_version
    write(out_unitp,*) '== QML path:    ',trim(adjustl(QML_path))
    write(out_unitp,*) '== Initialization of the Model =================='

    ! set the "File_path" in the Lib_module.f90
    File_path = trim(adjustl(QML_path))


    IF (present(adiabatic)) THEN
      adiabatic_loc = adiabatic
    ELSE
      adiabatic_loc = .TRUE.
    END IF

    IF (present(option)) THEN
      option_loc = option
    ELSE
      option_loc = -1
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

    IF (present(pot_name)) THEN
      QModel%pot_name  = strdup(pot_name)
      CALL string_uppercase_TO_lowercase(QModel%pot_name)
    ELSE
      IF (.NOT. read_param_loc) THEN
        write(out_unitp,*) 'ERROR in Init_Model'
        write(out_unitp,*) ' pot_name is not present and read_param=F'
        STOP ' pot_name is not present and read_param=F'
      END IF
    END IF


    IF (present(ndim)) THEN
      QModel%ndim      = ndim
    ELSE
      QModel%ndim      = 1
    END IF
    IF (present(nsurf)) THEN
      QModel%nsurf     = nsurf
    ELSE
      QModel%nsurf     = 1
    END IF

    IF (present(PubliUnit)) THEN
      QModel%PubliUnit      = PubliUnit
    ELSE
      QModel%PubliUnit      = .FALSE.
    END IF

    CALL QML_dealloc_dnMat(QModel%Vec0)

    read_nml = .FALSE.
    IF (read_param_loc .OR.  QModel%pot_name == 'read_model') THEN

      IF (nio_loc /= in_unitp) THEN
        open(unit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
      END IF
      CALL Read_Model(QModel,nio_loc,option_loc,read_nml)

      IF (allocated(param_file_name_loc)) deallocate(param_file_name_loc)
    END IF

    IF (QModel%adiabatic) THEN
      write(out_unitp,*) 'Adiabatic potential . . .'
    ELSE
      write(out_unitp,*) 'Non-adiabatic potential . . .'
    END IF

    IF (QModel%numeric) write(out_unitp,*) 'You have decided to perform a numeric checking of the analytic formulas.'

    read_param_loc = (read_param_loc .AND. read_nml) ! this enables to not read the next namelist when read_param_loc=t

    CALL string_uppercase_TO_lowercase(QModel%pot_name)
    write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name

    SELECT CASE (QModel%pot_name)
    CASE ('morse')
      !! === README ==
      !! Morse potential: V(R) = D*(1-exp(-a*(r-Req))**2
      !! pot_name  = 'Morse'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 1744.60504565084306291455 au
      !! remark: Default parameters for H-F
      !! === END README ==
      QModel%ndim      = 1
      QModel%nsurf     = 1

      CALL Init_MorsePot(QModel%MorsePot,nio=nio_loc,read_param=read_param_loc)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
      QModel%d0GGdef(1,1) = ONE/QModel%MorsePot%mu

    CASE ('sigmoid')
      !! sigmoid function: A * 1/2(1+e*tanh((x-B)/C))  remark: e=+/-1
      QModel%nsurf     = 1
      QModel%ndim      = 1

      CALL Init_SigmoidPot(QModel%SigmoidPot,nio=nio_loc,read_param=read_param_loc)
      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

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
      QModel%ndim      = 1
      QModel%nsurf     = 1

      CALL Init_BuckPot(QModel%BuckPot,nio=nio_loc,read_param=read_param_loc)
      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    CASE ('hbond')
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
      !! ref:  Eq 3.79 of J. Beutier, thesis.
      !! === END README ==
      QModel%ndim      = 2
      QModel%nsurf     = 1

      CALL Init_LinearHBondPot(QModel%LinearHBondPot,             &
                               nio=nio_loc,read_param=read_param_loc,   &
                               PubliUnit=QModel%PubliUnit)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
      QModel%d0GGdef(1,1) = ONE/QModel%LinearHBondPot%muQQ
      QModel%d0GGdef(2,2) = ONE/QModel%LinearHBondPot%muq

    CASE ('henonheiles')
      !! === README ==
      !! HenonHeiles potential
      !! pot_name  = 'HenonHeiles'
      !! ndim      > 1
      !! nsurf     = 1
      !! reduced masses(:)      = ONE au
      !! ref:  parameters taken from M. Nest, H.-D. Meyer, J. Chem. Phys. 117 (2002) 10499. doi:10.1063/1.1521129
      !! === END README ==
      QModel%nsurf     = 1

      CALL Init_HenonHeilesPot(QModel%HenonHeilesPot,ndim=QModel%ndim, &
                                 nio=nio_loc,read_param=read_param_loc)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


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
      QModel%nsurf     = 2
      QModel%ndim      = 1

      write(out_unitp,*) 'option_loc',option_loc
      CALL Init_TullyPot(QModel%TullyPot, option=option_loc,      &
                           nio=nio_loc,read_param=read_param_loc)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
      QModel%d0GGdef(1,1) = ONE/QModel%TullyPot%mu


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
      IF (.NOT. present(nsurf)) QModel%nsurf = 4
      QModel%ndim      = 1

      CALL Init_OneDSOC(QModel%OneDSOCPot,option=option_loc,nsurf=QModel%nsurf,  &
                        nio=nio_loc,read_param=read_param_loc)
      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
      QModel%d0GGdef(1,1) = ONE/QModel%OneDSOCPot%mu

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
      QModel%nsurf     = 4
      QModel%ndim      = 1

      CALL Init_OneDSOC_2S1T(QModel%OneDSOC_2S1T_Pot,option=option_loc,&
                             nio=nio_loc,read_param=read_param_loc)
      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
      QModel%d0GGdef(1,1) = ONE/QModel%OneDSOC_2S1T_Pot%mu


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
      QModel%nsurf     = 3
      QModel%ndim      = 2

      CALL Init_PhenolPot(QModel%PhenolPot,PubliUnit=QModel%PubliUnit)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
      ! The metric tensor of Tnum with rigid_type=100 from B3LYP/6-31G** of the ground stat
      QModel%d0GGdef(1,1) = QModel%PhenolPot%G_RR
      QModel%d0GGdef(2,2) = QModel%PhenolPot%G_ThTh


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
      QModel%nsurf     = 2
      QModel%ndim      = 2

      CALL Init_TwoDPot(QModel%Para_TwoD,PubliUnit=QModel%PubliUnit)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
      ! The metric tensor from the publication
      QModel%d0GGdef(1,1) = ONE/QModel%Para_TwoD%muX
      QModel%d0GGdef(2,2) = ONE/QModel%Para_TwoD%muY

    CASE ('psb3')
      !! === README ==
      !! Model for the photo-isomerization of the penta-2,4-dieniminium (PSB3) cation.
      !! pot_name  = 'PSB3'
      !! ndim      = 3
      !! nsurf     = 2
      !! remarks: two options are possible (option = 1,2)
      !! The default is option=1 (unpublished).
      !! The parameters for option=2 come from the following reference.
      !! ref: E. Marsili, M. H. Farag, X. Yang, L. De Vico, and M. Olivucci, JPCA, 123, 1710–1719 (2019). https://doi.org/10.1021/acs.jpca.8b10010
      !! === END README ==
      QModel%nsurf     = 2
      QModel%ndim      = 3

      CALL Init_PSB3Pot(QModel%PSB3Pot,option=option_loc,      &
                           nio=nio_loc,PubliUnit=QModel%PubliUnit)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    CASE ('hono')
      QModel%nsurf     = 1
      QModel%ndim      = 6

      CALL Init_HONO(QModel%HONOPot,option=option_loc,            &
                             nio=nio_loc,PubliUnit=QModel%PubliUnit)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    CASE ('hnnhp')
      QModel%nsurf     = 1
      QModel%ndim      = 6

      CALL Init_HNNHp(QModel%HNNHpPot,QModel%ndim,option=option_loc,&
                             nio=nio_loc,PubliUnit=QModel%PubliUnit)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    CASE ('h2sin')
      QModel%nsurf     = 1
      QModel%ndim      = 6

      CALL Init_H2SiN(QModel%H2SiNPot,QModel%ndim,option=option_loc,&
                      nio=nio_loc,PubliUnit=QModel%PubliUnit)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    CASE ('h2nsi')
      QModel%nsurf     = 1
      QModel%ndim      = 6

      CALL Init_H2NSi(QModel%H2NSiPot,QModel%ndim,option=option_loc,&
                      nio=nio_loc,PubliUnit=QModel%PubliUnit)

      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    CASE ('template')
      !! 3D-potential with 1 surface
      QModel%nsurf     = 1
      QModel%ndim      = 3

      CALL Init_TemplatePot(QModel%TemplatePot)
      CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
      QModel%d0GGdef = QModel%d0GGdef * 2000._Rkind

    CASE DEFAULT
        write(out_unitp,*) ' ERROR in Init_Model'
        write(out_unitp,*) ' This model/potential is unknown. QModel%pot_name: ',QModel%pot_name
        STOP 'STOP in Init_Model: Other potentials have to be done'
    END SELECT

    IF (read_param_loc .AND. nio_loc /= in_unitp) THEN
       close(nio_loc)
    END IF

    write(out_unitp,*) '================================================='
    write(out_unitp,*) ' ref geometry or minimum (option=0)'
    allocate(Q0(QModel%ndim))

    CALL get_Q0_Model(Q0,QModel,option=0)
    CALL Write_RVec(Q0,out_unitp,nbcol1=6,name_info='Q0: ')

    deallocate(Q0)

    !CALL Write_Model(QModel)
    write(out_unitp,*) '================================================='
    write(out_unitp,*) '================================================='


  END SUBROUTINE Init_Model

  SUBROUTINE get_Q0_Model(Q0,QModel,option)
  USE mod_Lib
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q0(:)
    TYPE (QModel_t), intent(in)               :: QModel
    integer,            intent(in)               :: option


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='get_Q0_Model'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    IF (size(Q0) /= QModel%ndim) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The size of Q0 is not QModel%ndim: '
      write(out_unitp,*) ' size(Q0)',size(Q0)
      write(out_unitp,*) ' ndim',QModel%ndim
      STOP 'STOP in get_Q0_Model: Wrong Q0 size'
    END IF

    Q0(:) = ZERO

    SELECT CASE (QModel%pot_name)
    CASE ('morse')
      CALL get_Q0_Morse(Q0(1),QModel%MorsePot)

    CASE ('buck')
      CALL get_Q0_Buck(Q0(1),QModel%BuckPot)

    CASE ('sigmoid')
      CALL get_Q0_Sigmoid(Q0(1),QModel%SigmoidPot)

    CASE ('hbond')
      CONTINUE

    CASE ('henonheiles')
      Q0(:) = ZERO

    CASE ('tully')
      CALL get_Q0_Tully(Q0(1),QModel%TullyPot)

    CASE ('1dsoc','1dsoc_1s1t')
      CALL get_Q0_OneDSOC(Q0(1),QModel%OneDSOCPot)

    CASE ('1dsoc_2s1t')
      CALL get_Q0_OneDSOC_2S1T(Q0(1),QModel%OneDSOC_2S1T_Pot)

    CASE ('phenol')
      CONTINUE

    CASE ('twod')
      CALL get_Q0_TwoD(Q0,QModel%Para_TwoD,option)

    CASE ('psb3')
      CONTINUE

    CASE ('hono')
      CALL get_Q0_HONO(Q0,QModel%HONOPot,option)

    CASE ('hnnhp')
      CALL get_Q0_HNNHp(Q0,QModel%HNNHpPot,option)

    CASE ('h2sin')
      CALL get_Q0_H2SiN(Q0,QModel%H2SiNPot,option)

    CASE ('h2nsi')
      CALL get_Q0_H2NSi(Q0,QModel%H2NSiPot,option)

    CASE ('template')
      CONTINUE

    CASE DEFAULT
      write(out_unitp,*) ' ERROR in name_sub'
      write(out_unitp,*) ' This model/potential is unknown. QModel%pot_name: ',QModel%pot_name
      STOP 'STOP in get_Q0_Model: Other models have to be done'
    END SELECT

    IF (debug) THEN
      CALL Write_RVec(Q0,out_unitp,nbcol1=5,name_info='Q0: ')
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE get_Q0_Model


  SUBROUTINE Eval_Pot(QModel,Q,PotVal,nderiv,NAC,Vec,numeric)
  USE mod_dnS
  !USE mod_diago
  IMPLICIT NONE

    TYPE (QModel_t),  intent(inout)            :: QModel
    TYPE (dnMat_t),     intent(inout)            :: PotVal
    real (kind=Rkind),   intent(in)               :: Q(:)
    integer,             intent(in),    optional  :: nderiv
    TYPE (dnMat_t),     intent(inout), optional  :: NAC,Vec
    logical,             intent(in),    optional  :: numeric

    ! local variables
    integer                    :: nderiv_loc
    TYPE (dnMat_t)            :: Vec_loc,NAC_loc
    logical                    :: numeric_loc,adia_loc

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
      write(out_unitp,*) '  numeric   ',QModel%numeric
      write(out_unitp,*) '  adiabatic ',QModel%adiabatic
      flush(out_unitp)
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(2,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF

    IF (present(numeric)) THEN
      numeric_loc = (numeric .AND. nderiv_loc > 0)
    ELSE
      numeric_loc = (QModel%numeric .AND. nderiv_loc > 0)
    END IF

    adia_loc = (QModel%adiabatic .AND. QModel%nsurf > 1)

    IF (numeric_loc) THEN  ! numerical
      IF (adia_loc) THEN
         CALL Eval_Pot_Numeric_dia(QModel,Q,PotVal,nderiv_loc)
      ELSE
        IF (present(Vec)) THEN
          IF (present(NAC)) THEN
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,Vec,NAC)
          ELSE
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,Vec,NAC_loc)
            CALL QML_dealloc_dnMat(NAC_loc)
          END IF
        ELSE
          IF (present(NAC)) THEN
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,Vec_loc,NAC)
          ELSE
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,Vec_loc,NAC_loc)
            CALL QML_dealloc_dnMat(NAC_loc)
          END IF
          CALL QML_dealloc_dnMat(Vec_loc)
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



    IF (debug) THEN
      IF ( QModel%adiabatic) write(out_unitp,*) 'PotVal (adia)'
      IF ( .NOT. QModel%adiabatic) write(out_unitp,*) 'PotVal (dia)'
      CALL QML_Write_dnMat(PotVal,6)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF


  END SUBROUTINE Eval_Pot

  SUBROUTINE Eval_Pot_ana(QModel,Q,PotVal,nderiv,NAC,Vec)
  USE mod_dnS
  IMPLICIT NONE

    TYPE (QModel_t), intent(inout)            :: QModel
    TYPE (dnMat_t),    intent(inout)            :: PotVal
    real (kind=Rkind),  intent(in)               :: Q(:)
    integer,            intent(in)               :: nderiv
    TYPE (dnMat_t),    intent(inout), optional  :: NAC,Vec

    ! local variables
    integer                    :: i,j,id,nderiv_loc
    TYPE (dnMat_t)            :: PotVal_dia,Vec_loc,NAC_loc
    TYPE (dnS_t), allocatable   :: dnQ(:)
    TYPE (dnS_t), allocatable   :: Mat_OF_PotDia(:,:)

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot_ana'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      write(out_unitp,*) '   nderiv    ',nderiv
      write(out_unitp,*) '   adiabatic ',QModel%adiabatic
      flush(out_unitp)
    END IF

    nderiv_loc = nderiv


    IF ( QML_Check_NotAlloc_dnMat(PotVal,nderiv_loc) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%nsurf,ndim=QModel%ndim,&
                          nderiv=nderiv_loc)
    END IF
    PotVal = ZERO


    ! allocate Mat_OF_PotDia
    allocate(Mat_OF_PotDia(QModel%nsurf,QModel%nsurf))
    DO j=1,size(Mat_OF_PotDia(1,:))
    DO i=1,size(Mat_OF_PotDia(:,1))
        CALL QML_alloc_dnS(Mat_OF_PotDia(i,j),QModel%ndim,nderiv_loc)
    END DO
    END DO

    ! intialization of the dnQ(:)
    allocate(dnQ(QModel%ndim))
    DO i=1,QModel%ndim
      dnQ(i) = QML_init_dnS(Q(i),ndim=QModel%ndim,nderiv=nderiv_loc,iQ=i) ! to set up the derivatives
    END DO

    SELECT CASE (QModel%pot_name)
    CASE ('morse')
      CALL Eval_MorsePot(Mat_OF_PotDia,dnQ(1),QModel%MorsePot,nderiv_loc)

    CASE ('buck')
      CALL Eval_BuckPot(Mat_OF_PotDia,dnQ(1),QModel%BuckPot,nderiv_loc)

    CASE ('sigmoid')
      CALL Eval_SigmoidPot(Mat_OF_PotDia,dnQ(1),QModel%SigmoidPot,nderiv_loc)

    CASE ('hbond')
      CALL Eval_LinearHBondPot(Mat_OF_PotDia,dnQ,QModel%LinearHBondPot,nderiv_loc)

    CASE ('henonheiles') ! Q(:) is used insted of dnQ(:)
      CALL Eval_HenonHeilesPot(Mat_OF_PotDia,Q,QModel%HenonHeilesPot,nderiv_loc)

    CASE ('tully')
      CALL Eval_TullyPot(Mat_OF_PotDia,dnQ(1),QModel%TullyPot,nderiv_loc)

    CASE ('1dsoc','1dsoc_1s1t')
      CALL Eval_OneDSOC(Mat_OF_PotDia,dnQ(1),QModel%OneDSOCPot,nderiv_loc)

    CASE ('1dsoc_2s1t')
      CALL Eval_OneDSOC_2S1T(Mat_OF_PotDia,dnQ(1),QModel%OneDSOC_2S1T_Pot,nderiv_loc)

    CASE ('phenol')
      CALL Eval_PhenolPot(Mat_OF_PotDia,dnQ,QModel%PhenolPot,nderiv_loc)

    CASE ('twod')
      CALL Eval_TwoDPot(Mat_OF_PotDia,dnQ,QModel%Para_TwoD,nderiv_loc)

    CASE ('psb3')
      CALL Eval_PSB3Pot(Mat_OF_PotDia,dnQ,QModel%PSB3Pot,nderiv_loc)

    CASE ('hono')
      CALL Eval_HONOPot(Mat_OF_PotDia,dnQ,QModel%HONOPot,nderiv_loc)

    CASE ('hnnhp')
      CALL Eval_HNNHpPot(Mat_OF_PotDia,dnQ,QModel%HNNHpPot,nderiv_loc)

    CASE ('h2sin')
      CALL Eval_H2SiNPot(Mat_OF_PotDia,dnQ,QModel%H2SiNPot,nderiv_loc)

    CASE ('h2nsi')
      CALL Eval_H2NSiPot(Mat_OF_PotDia,dnQ,QModel%H2NSiPot,nderiv_loc)

    CASE ('template')
      CALL Eval_TemplatePot(Mat_OF_PotDia,dnQ,QModel%TemplatePot,nderiv_loc)

    CASE DEFAULT
      write(out_unitp,*) ' ERROR in Eval_Pot'
      write(out_unitp,*) ' This model/potential is unknown. QModel%pot_name: ',QModel%pot_name
      STOP 'STOP in Eval_Pot: Other potentials have to be done'
    END SELECT

    PotVal = Mat_OF_PotDia ! transfert the potential and its derivatives to the matrix form (PotVal)

    ! deallocation
    DO i=1,size(dnQ)
      CALL QML_dealloc_dnS(dnQ(i))
    END DO
    deallocate(dnQ)

    DO j=1,size(Mat_OF_PotDia(1,:))
    DO i=1,size(Mat_OF_PotDia(:,1))
      CALL QML_dealloc_dnS(Mat_OF_PotDia(i,j))
    END DO
    END DO
    deallocate(Mat_OF_PotDia)
    ! end deallocation

    IF ( QModel%adiabatic .AND. QModel%nsurf > 1) THEN
      IF (debug) THEN
        write(out_unitp,*) 'PotVal (dia)'
        CALL QML_Write_dnMat(PotVal,6)
        flush(out_unitp)
      END IF

      PotVal_dia = PotVal

      IF (present(Vec)) THEN
        IF (present(NAC)) THEN
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec,QModel%Vec0,NAC,nderiv)
        ELSE
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec,QModel%Vec0,NAC_loc,nderiv)
          CALL QML_dealloc_dnMat(NAC_loc)
        END IF
      ELSE
        IF (present(NAC)) THEN
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec_loc,QModel%Vec0,NAC,nderiv)
        ELSE
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec_loc,QModel%Vec0,NAC_loc,nderiv)
          CALL QML_dealloc_dnMat(NAC_loc)
        END IF
        CALL QML_dealloc_dnMat(Vec_loc)
      END IF
      CALL QML_dealloc_dnMat(PotVal_dia)

    END IF


    IF (debug) THEN
      IF ( QModel%adiabatic) write(out_unitp,*) 'PotVal (adia)'
      IF ( .NOT. QModel%adiabatic) write(out_unitp,*) 'PotVal (dia)'
      CALL QML_Write_dnMat(PotVal,6)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF


  END SUBROUTINE Eval_Pot_ana

  SUBROUTINE Eval_Pot_Numeric_dia(QModel,Q,PotVal,nderiv)
  IMPLICIT NONE

    TYPE (QModel_t), intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind) , intent(in)     :: Q(:)
    integer,            intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                    :: PotVal_loc0
    integer                            :: i,j

    !write(out_unitp,*) 'coucou0 Eval_Pot_Numeric' ; flush(out_unitp)

    IF (QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%nsurf,ndim=QModel%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    allocate(Q_loc(QModel%ndim))
    Q_loc(:) = Q
    CALL QML_alloc_dnMat(PotVal_loc0,QModel%nsurf,QModel%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0)
    PotVal%d0 = PotVal_loc0%d0


    IF (nderiv >= 1) THEN ! 1st derivatives

      ! Numeric evaluation of forces
      DO i=1,QModel%ndim

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

      DO i=1,QModel%ndim
      DO j=i+1,QModel%ndim

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
    CALL QML_dealloc_dnMat(PotVal_loc0)

  END SUBROUTINE Eval_Pot_Numeric_dia
  SUBROUTINE Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv,Vec,NAC)
  IMPLICIT NONE

    TYPE (QModel_t), intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind) , intent(in)     :: Q(:)
    integer,            intent(in)     :: nderiv
    TYPE (dnMat_t),    intent(inout)  :: Vec,NAC

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                    :: PotVal_loc0,Vec_loc0
    integer                            :: i,j
    real (kind=Rkind), allocatable     :: tVec(:,:)

    !write(out_unitp,*) 'coucou0 Eval_Pot_Numeric' ; flush(out_unitp)

    IF (QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%nsurf,ndim=QModel%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    IF (QML_Check_NotAlloc_dnMat(Vec,nderiv) ) THEN
      CALL QML_alloc_dnMat(Vec,nsurf=QModel%nsurf,ndim=QModel%ndim,&
                          nderiv=nderiv)
    END IF
    Vec = ZERO

    IF (QML_Check_NotAlloc_dnMat(NAC,nderiv) ) THEN
      CALL QML_alloc_dnMat(NAC,nsurf=QModel%nsurf,ndim=QModel%ndim,&
                          nderiv=nderiv)
    END IF
    NAC = ZERO

    allocate(Q_loc(QModel%ndim))
    Q_loc(:) = Q
    CALL QML_alloc_dnMat(PotVal_loc0,QModel%nsurf,QModel%ndim,nderiv=0)
    CALL QML_alloc_dnMat(Vec_loc0,   QModel%nsurf,QModel%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0,vec=Vec_loc0)

    PotVal%d0 = PotVal_loc0%d0
    Vec%d0    = Vec_loc0%d0
    CALL Init_IdMat(NAC%d0,QModel%nsurf)

    allocate(tVec(QModel%nsurf,QModel%nsurf))
    tVec(:,:)      = transpose(Vec%d0)


    IF (nderiv >= 1) THEN ! 1st derivatives

      ! Numeric evaluation of forces
      DO i=1,QModel%ndim

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

      DO i=1,QModel%ndim
      DO j=i+1,QModel%ndim

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
    CALL QML_dealloc_dnMat(PotVal_loc0)
    CALL QML_dealloc_dnMat(Vec_loc0)

  END SUBROUTINE Eval_Pot_Numeric_adia
  SUBROUTINE dia_TO_adia(PotVal_dia,PotVal_adia,Vec,Vec0,NAC,nderiv)
    USE mod_diago
    IMPLICIT NONE

    TYPE (dnMat_t), intent(in)               :: PotVal_dia
    TYPE (dnMat_t), intent(inout)            :: PotVal_adia,Vec,Vec0,NAC

    integer, intent(in), optional             :: nderiv

    ! local variable
    integer                        :: i,j,k,id,jd,nderiv_loc,ndim,nsurf
    real (kind=Rkind)              :: ai,aj,aii,aij,aji,ajj,th,cc,ss
    real (kind=Rkind), allocatable :: Eig(:),tVec(:,:),Vdum(:),Vi(:)

    TYPE (dnMat_t)                :: PotVal_dia_onadia

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='dia_TO_adia'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(2,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF


    IF ( QML_Check_NotAlloc_dnMat(PotVal_dia,nderiv_loc) ) THEN
      write(out_unitp,*) ' The diabatic potential MUST be allocated!'
      CALL QML_Write_dnMat(PotVal_dia)
      STOP 'PotVal_dia%dn NOT allocated in "dia_TO_adia"'
    END IF
    IF (debug) THEN
      write(out_unitp,*) 'PotVal_dia'
      CALL QML_Write_dnMat(PotVal_dia,6)
      flush(out_unitp)
    END IF

    nsurf = QML_get_nsurf_FROM_dnMat(PotVal_dia)
    ndim  = QML_get_ndim_FROM_dnMat(PotVal_dia)

    IF ( QML_Check_NotAlloc_dnMat(PotVal_adia,nderiv_loc) ) THEN
      CALL QML_alloc_dnMat(PotVal_adia,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    END IF
    PotVal_adia = ZERO


    IF ( QML_Check_NotAlloc_dnMat(Vec,nderiv_loc) ) THEN
      CALL QML_alloc_dnMat(Vec,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    END IF
    Vec = ZERO

    IF ( QML_Check_NotAlloc_dnMat(NAC,nderiv_loc) ) THEN
      CALL QML_alloc_dnMat(NAC,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    END IF
    NAC = ZERO


    ! local variables
    CALL QML_alloc_dnMat(PotVal_dia_onadia,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    PotVal_dia_onadia = ZERO


    allocate(Eig(nsurf))
    allocate(tVec(nsurf,nsurf))
    CALL diagonalization(PotVal_dia%d0,Eig,Vec%d0,nsurf)
    IF (QML_Check_NotAlloc_dnMat(Vec0,nderiv=0)) THEN
       !$OMP CRITICAL (CRIT_dia_TO_adia)
       IF (debug) write(out_unitp,*) 'init Vec0'
       CALL QML_alloc_dnMat(Vec0,nsurf=nsurf,ndim=ndim,nderiv=0)
       Vec0%d0 = Vec%d0
       !$OMP END CRITICAL (CRIT_dia_TO_adia)
    ELSE ! change the phase if required
       IF (debug) write(out_unitp,*) 'Change phase?'
       flush(out_unitp)

       DO i=1,nsurf
         IF (dot_product(Vec0%d0(:,i),Vec%d0(:,i)) < ZERO) Vec%d0(:,i) = -Vec%d0(:,i)
       END DO

       IF (debug) THEN
         write(out_unitp,*) 'Vec before rotation'
         CALL QML_Write_dnMat(Vec,6)
       END IF
       !For degenerated eigenvectors (works only with 2 vectors)
       DO i=1,nsurf-1
        IF ( abs(Eig(i)-Eig(i+1)) < epsi) THEN
          j = i+1
          IF (debug) write(out_unitp,*) 'degenerated vectors',i,j

          aii = dot_product(Vec0%d0(:,i),Vec%d0(:,i))
          aji = dot_product(Vec0%d0(:,j),Vec%d0(:,i))
          aij = dot_product(Vec0%d0(:,i),Vec%d0(:,j))
          ajj = dot_product(Vec0%d0(:,j),Vec%d0(:,j))

          th = ( atan2(aij,ajj) -atan2(aji,aii) ) * HALF
          IF (debug) write(out_unitp,*) 'theta',th

          cc = cos(th)
          ss = sin(th)

          DO k=1,nsurf
           ai = Vec%d0(k,i)
           aj = Vec%d0(k,j)

           Vec%d0(k,i) =  cc * ai + ss * aj
           Vec%d0(k,j) = -ss * ai + cc * aj

          END DO

        END IF
       END DO

       IF (debug) THEN
         write(out_unitp,*) 'Vec after rotation'
         CALL QML_Write_dnMat(Vec,6)
       END IF

    END IF
    tVec = transpose(Vec%d0)

    IF (debug) write(out_unitp,*) 'Eig',Eig

    ! transformation of PotVal_dia on the adiabatic basis (Vec)
    !PotVal_dia_onadia%d0 = matmul(tVec,matmul(PotVal_dia%d0,Vec%d0))
    DO i=1,nsurf
      PotVal_dia_onadia%d0(i,i) = Eig(i)
    END DO
    IF (nderiv_loc > 0) THEN
      DO id=1,ndim
        PotVal_dia_onadia%d1(:,:,id) = matmul(tVec,matmul(PotVal_dia%d1(:,:,id),Vec%d0))
      END DO
    END IF
    IF (nderiv_loc > 1) THEN
      DO id=1,ndim
      DO jd=1,ndim
        PotVal_dia_onadia%d2(:,:,id,jd) = matmul(tVec,matmul(PotVal_dia%d2(:,:,id,jd),Vec%d0)  )
      END DO
      END DO
    END IF
    deallocate(Eig)
    deallocate(tVec)
    IF (debug) THEN
      write(out_unitp,*) "< Psi_j I dnHdia I Psi_i>"
      CALL QML_Write_dnMat(PotVal_dia_onadia,6)
      flush(out_unitp)
    END IF


    ! project the eigenvector derivatives on the eigenvectors => NAC

    ! no derivative
    PotVal_adia%d0 = ZERO
    NAC%d0         = ZERO
    DO i=1,nsurf
      PotVal_adia%d0(i,i) = PotVal_dia_onadia%d0(i,i)
      NAC%d0(i,i)         = ONE
    END DO

    ! 1st order derivatives
    IF (nderiv_loc > 0) THEN

      ! eigenvalue derivatives
      PotVal_adia%d1 = ZERO
      DO id=1,ndim
      DO i=1,nsurf
        PotVal_adia%d1(i,i,id) = PotVal_dia_onadia%d1(i,i,id)
      END DO
      END DO

      ! eigenvector derivatives projected on the eigenvectors
      DO id=1,ndim
      DO i=1,nsurf ! I Psi_i' >
      DO j=1,nsurf ! projection on < Psi_j I
        IF (j /= i) THEN
          NAC%d1(j,i,id) = - PotVal_dia_onadia%d1(j,i,id)/              &
                            ( PotVal_adia%d0(j,j) - PotVal_adia%d0(i,i) )
        ELSE
          NAC%d1(j,i,id) = ZERO
        END IF
      END DO
      END DO
      END DO

    END IF


    ! 2d order derivatives
    IF (nderiv_loc > 1) THEN
      PotVal_adia%d2 = ZERO
      allocate(Vdum(nsurf))
      ! eigenvalue derivatives
      DO id=1,ndim
      DO jd=1,ndim
        DO i=1,nsurf
          Vdum = matmul(PotVal_dia_onadia%d2(:,:,id,jd),NAC%d0(:,i))    + &
                 matmul(PotVal_dia_onadia%d1(:,:,jd),   NAC%d1(:,i,id)) + &
                 matmul(PotVal_dia_onadia%d1(:,:,id),   NAC%d1(:,i,jd))

          !write(out_unitp,*) 'Dum',id,jd,i,':',Vdum
          PotVal_adia%d2(i,i,id,jd) = dot_product(NAC%d0(:,i),Vdum)
        END DO

      END DO
      END DO
      deallocate(Vdum)
    END IF


    CALL QML_dealloc_dnMat(PotVal_dia_onadia)

    IF (debug) THEN
      write(out_unitp,*) 'PotVal_adia'
      CALL QML_Write_dnMat(PotVal_adia,6)

      write(out_unitp,*) 'Vec'
      CALL QML_Write_dnMat(Vec,6)

      write(out_unitp,*) 'NAC'
      CALL QML_Write_dnMat(NAC,6)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE dia_TO_adia

  SUBROUTINE Write_Model(QModel,nio)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (QModel_t), intent(in)              :: QModel
    integer,            intent(in), optional    :: nio

    integer :: nio_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = 6
    END IF

    IF (nio_loc /= 6) THEN
      open(nio_loc,file=trim(adjustl(QModel%pot_name))//'.out',form='formatted')
    END IF


    write(nio_loc,*) '-----------------------------------------------'
    write(nio_loc,*) 'Output file for potential library'
    write(nio_loc,*)
    write(nio_loc,*) 'Potential parameters are written just below'
    write(nio_loc,*)
    write(nio_loc,*) 'nsurf:     ',QModel%nsurf
    write(nio_loc,*) 'ndim:      ',QModel%ndim
    write(nio_loc,*) 'numeric:   ',QModel%numeric
    write(nio_loc,*) 'adiabatic: ',QModel%adiabatic
    write(nio_loc,*)

    IF (allocated(QModel%d0GGdef)) THEN
      write(nio_loc,*) 'Deformation metric tensor (~ 1/Mii)'
      CALL Write_RMat(QModel%d0GGdef,nio_loc,nbcol1=5)
    END IF

    SELECT CASE (QModel%pot_name)
    CASE ('morse')
      CALL Write_MorsePot(QModel%MorsePot,nio=nio_loc)
    CASE ('sigmoid')
      CALL Write_SigmoidPot(QModel%SigmoidPot,nio=nio_loc)
    CASE ('buck')
      CALL Write_BuckPot(QModel%BuckPot,nio=nio_loc)
   CASE ('hbond')
      CALL Write_LinearHBondPot(QModel%LinearHBondPot,nio=nio_loc)
    CASE ('henonheiles')
        CALL Write_HenonHeilesPot(QModel%HenonHeilesPot,nio=nio_loc)
    CASE ('tully')
        CALL Write_TullyPot(QModel%TullyPot,nio=nio_loc)
    CASE ('1dsoc','1dsoc_1s1t')
        CALL Write_OneDSOC(QModel%OneDSOCPot,nio=nio_loc)
    CASE ('1dsoc_2s1t')
        CALL Write_OneDSOC_2S1T(QModel%OneDSOC_2S1T_Pot,nio=nio_loc)
    CASE ('phenol')
        CALL Write_PhenolPot(QModel%PhenolPot,nio=nio_loc)
    CASE ('psb3')
        CALL Write_PSB3Pot(QModel%PSB3Pot,nio=nio_loc)
    CASE ('hono')
        CALL Write_HONO(QModel%HONOPot,nio=nio_loc)
    CASE ('hnnhp')
        CALL Write_HNNHp(QModel%HNNHpPot,nio=nio_loc)
    CASE ('h2sin')
        CALL Write_H2SiN(QModel%H2SiNPot,nio=nio_loc)
    CASE ('h2nsi')
        CALL Write_H2NSi(QModel%H2NSiPot,nio=nio_loc)
    CASE ('template')
        CALL  Write_TemplatePot(QModel%TemplatePot,nio=nio_loc)
    CASE DEFAULT
        write(nio_loc,*) 'WARNING in Write_Model: Other potentials have to be done'
    END SELECT
 
     IF (nio_loc /= 6) THEN
      close(nio_loc)
    END IF
    write(nio_loc,*) '-----------------------------------------------'
    flush(nio_loc)


  END SUBROUTINE Write_Model
  SUBROUTINE Write0_Model(QModel,pot_name,nio)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (QModel_t), intent(in), optional    :: QModel
    character (len=*),  intent(in), optional    :: pot_name
    integer,            intent(in), optional    :: nio

    character (len=:), allocatable :: pot_name_loc
    integer :: nio_loc

    IF (present(pot_name) .AND. present(QModel)) THEN
      write(out_unitp,*) 'ERROR in Write0_Model'
      write(out_unitp,*) ' pot_name and QModel are both present'
      write(out_unitp,*) ' ONE only MUST be present! CHECK the source'
      STOP 'pot_name and QModel are both present'
    END IF
    IF (.NOT. present(pot_name) .AND. .NOT. present(QModel)) THEN
      write(out_unitp,*) 'ERROR in Write0_Model'
      write(out_unitp,*) ' pot_name and QModel are both absent'
      write(out_unitp,*) ' ONE only MUST be present! CHECK the source'
      STOP ' pot_name and QModel are both absent'
    END IF

    IF (present(pot_name)) THEN
      pot_name_loc = strdup(pot_name)
    ELSE
      pot_name_loc = strdup(QModel%pot_name)
    END IF


    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = 6
    END IF

    IF (nio_loc /= 6) THEN
      open(nio_loc,file=trim(adjustl(QModel%pot_name))//'.out',form='formatted')
    END IF



    write(nio_loc,*) 'QUANTUM MODEL default parameters'
    flush(nio_loc)


    write(nio_loc,*)
    write(nio_loc,*) 'Potential parameters are written just below'
    write(nio_loc,*)
    IF (present(QModel)) THEN
      write(nio_loc,*) 'nsurf:     ',QModel%nsurf
      write(nio_loc,*) 'ndim:      ',QModel%ndim
      write(nio_loc,*) 'numeric:   ',QModel%numeric
      write(nio_loc,*) 'adiabatic: ',QModel%adiabatic
      write(nio_loc,*)

      IF (allocated(QModel%d0GGdef)) THEN
        write(nio_loc,*) 'Deformation metric tensor (~ 1/Mii)'
        CALL Write_RMat(QModel%d0GGdef,nio_loc,nbcol1=5)
      END IF
    END IF

    CALL string_uppercase_TO_lowercase(pot_name_loc)
    SELECT CASE (pot_name_loc)
    CASE ('morse')
      CALL Write0_MorsePot(nio=nio_loc)
    CASE ('sigmoid')
      CONTINUE
    CASE ('buck')
      CALL Write0_BuckPot(nio=nio_loc)
   CASE ('hbond')
      CALL Write0_LinearHBondPot(nio=nio_loc)
    CASE ('henonheiles')
        CALL Write0_HenonHeilesPot(nio=nio_loc)
    CASE ('tully')
        CALL Write0_TullyPot(nio=nio_loc)
    CASE ('1dsoc','1dsoc_1s1t')
        CALL Write0_OneDSOC(nio=nio_loc)
    CASE ('1dsoc_2s1t')
        CALL Write0_OneDSOC_2S1T(nio=nio_loc)
    CASE ('phenol')
        CALL Write0_PhenolPot(nio=nio_loc)
    CASE ('psb3')
        CALL Write0_PSB3Pot(nio=nio_loc)
    CASE ('hono')
        CALL Write0_HONO(nio=nio_loc)
    CASE ('h2sin')
        CALL Write0_H2SiN(nio=nio_loc)
    CASE ('h2nsi')
        CALL Write0_H2NSi(nio=nio_loc)
    CASE ('template')
        CALL  Write0_TemplatePot(nio=nio_loc)
    CASE DEFAULT
        write(nio_loc,*) 'WARNING in Write0_Model: Other potentials have to be done'
    END SELECT

     IF (nio_loc /= 6) THEN
      close(nio_loc)
    END IF

    write(nio_loc,*) 'END QUANTUM MODEL default parameters'
    flush(nio_loc)

    IF (allocated(pot_name_loc)) deallocate(pot_name_loc)


  END SUBROUTINE Write0_Model
  SUBROUTINE Write_QdnV_FOR_Model(Q,PotVal,QModel,Vec,NAC,info)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (QModel_t), intent(in)           :: QModel
    TYPE (dnMat_t),    intent(in)           :: PotVal
    real (kind=Rkind),  intent(in)           :: Q(:)
    TYPE (dnMat_t),    intent(in), optional :: Vec ! for non adiabatic couplings
    TYPE (dnMat_t),    intent(in), optional :: NAC ! for non adiabatic couplings
    character(len=*),   intent(in), optional :: info

    integer :: nio_loc,err_io

    CALL file_open2(trim(adjustl(QModel%pot_name))//'.txt',nio_loc,lformatted=.TRUE.,append=.TRUE.,err_file=err_io)
    IF (err_io /= 0) THEN
      write(out_unitp,*) 'ERROR in Write_QdnV_FOR_Model'
      write(out_unitp,*) ' Impossible to open the file "',trim(adjustl(QModel%pot_name))//'.txt','"'
      STOP 'Impossible to open the file'
    END IF
    !nio_loc = 99
    !open(nio_loc,file=trim(adjustl(QModel%pot_name))//'.xyz',form='formatted',POSITION='append')

    write(nio_loc,'(a)',advance='no') 'TEST output: '
    IF (present(info)) THEN
      write(nio_loc,'(a)',advance='no') info
    END IF
    IF (QModel%adiabatic) THEN
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

    IF (allocated(QModel%d0GGdef)) THEN
      write(nio_loc,*) 'd0GGdef'
      write(nio_loc,*) size(QModel%d0GGdef)
      write(nio_loc,*) QModel%d0GGdef
    END IF

    write(nio_loc,'(a)',advance='no') 'END_TEST output: '
    IF (present(info)) THEN
      write(nio_loc,'(a)',advance='no') info
    END IF
    IF (QModel%adiabatic) THEN
      write(nio_loc,'(a)') ' Adiabatic'
    ELSE
      write(nio_loc,'(a)') ' Diabatic'
    END IF

    close(nio_loc)

  END SUBROUTINE Write_QdnV_FOR_Model
  SUBROUTINE Check_analytical_numerical_derivatives(QModel,Q,nderiv)
  IMPLICIT NONE

    TYPE (QModel_t), intent(inout)            :: QModel
    real (kind=Rkind),dimension(:),intent(in)    :: Q
    integer, intent(in)                          :: nderiv

    TYPE (dnMat_t)           :: PotVal_ana,PotVal_num,PotVal_diff
    logical                   :: numeric_save
    real (kind=Rkind)         :: MaxPot,MaxDiffPot


      CALL QML_alloc_dnMat(PotVal_ana,nsurf=QModel%nsurf,ndim=QModel%ndim,&
                          nderiv=nderiv)

      CALL QML_alloc_dnMat(PotVal_num,nsurf=QModel%nsurf,ndim=QModel%ndim,&
                          nderiv=nderiv)

      numeric_save = QModel%numeric


      QModel%numeric = .FALSE.
      !write(out_unitp,*) 'coucou'
      CALL Eval_Pot(QModel,Q,PotVal_ana,nderiv)
      !write(out_unitp,*)   'PotVal_ana'
      !CALL QML_Write_dnMat(PotVal_ana,nio=out_unitp)
      !flush(out_unitp)
      QModel%numeric = .TRUE.
      CALL Eval_Pot(QModel,Q,PotVal_num,nderiv)
      !write(out_unitp,*)   'PotVal_num'
      !CALL QML_Write_dnMat(PotVal_num,nio=out_unitp)
      !flush(out_unitp)

      MaxPot     = QML_get_maxval_OF_dnMat(PotVal_ana)


      PotVal_diff = QML_dnMat2_MINUS_dnMat1(PotVal_num,PotVal_ana)
      MaxDiffPot  = QML_get_maxval_OF_dnMat(PotVal_diff)

      write(out_unitp,'(a,e9.2)') 'max of the relative Potential diff:',MaxDiffPot/MaxPot
      write(out_unitp,'(a,l9)')   'Potential diff (numer-ana), ZERO?  ',(MaxDiffPot/MaxPot <= step)

      IF (MaxDiffPot/MaxPot > step) THEN
        write(out_unitp,*)   'Potential diff (ana-numer)'
        CALL QML_Write_dnMat(PotVal_diff,nio=6)
      END IF

      CALL QML_dealloc_dnMat(PotVal_ana)
      CALL QML_dealloc_dnMat(PotVal_num)
      CALL QML_dealloc_dnMat(PotVal_diff)

      QModel%numeric = numeric_save

  END SUBROUTINE Check_analytical_numerical_derivatives

  SUBROUTINE Eval_pot_ON_Grid(QModel,Qmin,Qmax,nb_points,nderiv,grid_file)
  IMPLICIT NONE

    TYPE (QModel_t),           intent(inout)   :: QModel
    real (kind=Rkind),            intent(in)      :: Qmin(:),Qmax(:)
    integer, optional,            intent(in)      :: nb_points,nderiv
    character (len=*), optional,  intent(in)      :: grid_file

    integer           :: unit_grid_file

    integer                        :: i,iq,jq,i1,i2,nb_points_loc,nderiv_loc,ndim_loc
    integer, allocatable           :: i_Q(:)
    real (kind=Rkind), allocatable :: dQ(:),Q(:)
    TYPE (dnMat_t)                :: PotVal,NAC


    IF (size(Qmin) /= QModel%ndim .OR. size(Qmax) /= QModel%ndim) THEN
       write(out_unitp,*) ' ERROR in Eval_pot_ON_Grid'
       write(out_unitp,*) ' The size of Qmin or Qmax are different from QModel%ndim',size(Qmin),size(Qmax),QModel%ndim
       write(out_unitp,*) ' => Check the fortran'
       STOP 'ERROR in Eval_pot_ON_Grid: problem with QModel%ndim'
    END IF

    IF (present(grid_file)) THEN
      IF (len_trim(grid_file) == 0) THEN
        unit_grid_file = 6
      ELSE
        unit_grid_file = 99
        open(unit=unit_grid_file,file=trim(grid_file) )
      END IF
    ELSE
      unit_grid_file = 6
    END IF

    nb_points_loc = 100
    IF (present(nb_points)) nb_points_loc = nb_points
    nb_points_loc = max(nb_points_loc,2)

    IF (present(nderiv)) THEN
      nderiv_loc = nderiv
    ELSE
      nderiv_loc = 0
    END IF

    allocate(dQ(QModel%ndim))
    allocate(Q(QModel%ndim))
    allocate(i_Q(QModel%ndim))

    dQ(:)       = (Qmax-Qmin) / real(nb_points_loc-1,kind=Rkind)
    ndim_loc    = 0
    i_Q(:)      = 0
    DO i=1,QModel%ndim
      IF (dQ(i) /= ZERO) THEN
        ndim_loc = ndim_loc + 1
        i_Q(ndim_loc) = i
      END IF
    END DO
    write(out_unitp,*) 'QModel%ndim',QModel%ndim
    write(out_unitp,*) 'QModel%numeric',QModel%numeric
    write(out_unitp,*) 'ndim for the grid',ndim_loc
    write(out_unitp,*) 'i_Q',i_Q(1:ndim_loc)


    Q(:) = Qmin

    IF (ndim_loc == 1) THEN
      i1 = i_Q(1)
      DO iq=0,nb_points_loc-1
        Q(i1) = Qmin(i1) + dQ(i1)*real(iq,kind=Rkind)
        IF (QModel%nsurf > 1 .AND. QModel%adiabatic) THEN
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=max(1,nderiv_loc),NAC=NAC)

          IF (nderiv_loc == 0) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%nsurf),NAC%d1
          ELSE IF (nderiv_loc == 1) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%nsurf),(PotVal%d1(i,i,:),i=1,QModel%nsurf)
          ELSE
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%nsurf),(PotVal%d1(i,i,:),i=1,QModel%nsurf), &
                            (PotVal%d2(i,i,:,:),i=1,QModel%nsurf)
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

        IF (QModel%nsurf > 1 .AND. QModel%adiabatic) THEN
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=0,NAC=NAC)
          write(unit_grid_file,*) Q(i1),Q(i2),(PotVal%d0(i,i),i=1,QModel%nsurf)
        ELSE
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv_loc)
          write(unit_grid_file,*) Q(i1),Q(i2),PotVal%d0
        END IF
        flush(unit_grid_file)

      END DO
      write(unit_grid_file,*)
      END DO


    END IF

    CALL QML_dealloc_dnMat(PotVal)
    CALL QML_dealloc_dnMat(NAC)
    deallocate(dQ)
    deallocate(Q)
    deallocate(i_Q)

    IF (unit_grid_file /= 6) THEN
      close(unit_grid_file)
    END IF


  END SUBROUTINE Eval_pot_ON_Grid


  SUBROUTINE calc_pot(V,QModel,Q)
  IMPLICIT NONE

    TYPE (QModel_t),       intent(inout) :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated

    TYPE (dnMat_t)                         :: PotVal



    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%nsurf,ndim=QModel%ndim,nderiv=0)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=0)

    V = PotVal%d0

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot
  SUBROUTINE calc_pot_grad(V,g,QModel,Q)
  IMPLICIT NONE

    TYPE (QModel_t),   intent(inout)     :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMat_t)           :: PotVal



    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%nsurf,ndim=QModel%ndim,nderiv=1)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=1)

    V = PotVal%d0
    g = PotVal%d1

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot_grad
  SUBROUTINE calc_grad(g,QModel,Q)
  IMPLICIT NONE

    TYPE (QModel_t),   intent(inout)     :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMat_t)           :: PotVal



    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%nsurf,ndim=QModel%ndim,nderiv=1)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=1)

    g = PotVal%d1

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_grad
  SUBROUTINE calc_pot_grad_hess(V,g,h,QModel,Q)
  IMPLICIT NONE

    TYPE (QModel_t),   intent(inout)     :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMat_t)           :: PotVal



    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%nsurf,ndim=QModel%ndim,nderiv=2)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)

    V = PotVal%d0
    g = PotVal%d1
    h = PotVal%d2

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot_grad_hess
  SUBROUTINE calc_hess(h,QModel,Q)
  IMPLICIT NONE

    TYPE (QModel_t),   intent(inout)   :: QModel
    real (kind=Rkind),  intent(in)        :: Q(:)
    real (kind=Rkind),  intent(inout)     :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMat_t)           :: PotVal



    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%nsurf,ndim=QModel%ndim,nderiv=2)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)

    h = PotVal%d2

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_hess
END MODULE mod_Model
