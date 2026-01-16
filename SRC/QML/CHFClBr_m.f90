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

!> @brief Module which makes the initialization, calculation of the CHFClBr potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 30/11/2021
!!
MODULE QML_CHFClBr_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

  integer, parameter :: max_order = 100

!> @brief Derived type in which the CHFClBr parameters are set-up.
!> @brief CHFClBr(R) = Sum_i C_i * r-Req)**i
!> @brief Default parameters for H-F
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param norder      integer: order of the expansion (up to 4). The default is 4
  TYPE, EXTENDS (QML_Empty_t) :: QML_CHFClBr_t ! V(R) = Sum_i C_i * r-Req)**i
     PRIVATE
     integer          :: norder = 4
  CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_CHFClBr
    PROCEDURE :: Write_QModel     => Write_QML_CHFClBr
    PROCEDURE :: RefValues_QModel => RefValues_QML_CHFClBr
  END TYPE QML_CHFClBr_t

  PUBLIC :: QML_CHFClBr_t,Init_QML_CHFClBr,Write_QML_CHFClBr

CONTAINS
!> @brief Subroutine which makes the initialization of the CHFClBr parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param CHFClBrPot         TYPE(QML_CHFClBr_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_CHFClBr(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_CHFClBr_t)                         :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    integer :: i
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_CHFClBr'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) '  read_param:     ',read_param
      write(out_unit,*) '  nio_param_file: ',nio_param_file
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%pot_name = 'CHFClBr'
    QModel%ndim     = 9
    QModel%nsurf    = 1


    IF (debug) write(out_unit,*) 'init default CHFClBr parameters'
    QModel%norder = 4

    IF (read_param) THEN
      CALL Read_QML_CHFClBr(QModel,nio_param_file)
    END IF

    IF (QModel%norder > 4 .OR. QModel%norder < 2) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'Wrong norder value:',QModel%norder
      write(out_unit,*) 'Possible values: 2, 3, 4'
      STOP 'ERROR in Init_QML_CHFClBr: Wrong norder value'
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of CHFClBr'
    QModel%Q0 = [(ZERO,i=1,QModel%ndim)]

    IF (debug) write(out_unit,*) 'init d0GGdef of CHFClBr'
    QModel%masses = [952.84844038464450_Rkind,681.50589213793319_Rkind,505.79570923770643_Rkind, &
                     319.83919070617901_Rkind,270.72255106228903_Rkind,198.46622168004251_Rkind, &
                     179.17740535161565_Rkind,165.50265366523007_Rkind,66.985713015455389_Rkind]
    QModel%d0GGdef = Identity_Mat(QModel%ndim)
    DO i=1,QModel%ndim
      QModel%d0GGdef(i,i) = ONE/QModel%masses(i)
    END DO

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_CHFClBr

!> @brief Subroutine wich reads the CHFClBr parameters with a namelist.
!!   This can be called only from the "Init_QML_CHFClBr" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_CHFClBr_t):   derived type in which the parameters are set-up.
!! @param nio            integer           :   file unit to read the parameters.
  SUBROUTINE Read_QML_CHFClBr(QModel,nio)
    IMPLICIT NONE

    TYPE (QML_CHFClBr_t), intent(inout) :: QModel
    integer,              intent(in)    :: nio

    !local variable
    integer                       :: err_read
    integer                       :: norder

    namelist /CHFClBr/ norder

    write(out_unit,*) 'read CHFClBr namelist ...'
    flush(out_unit)

    read(nio,nml=CHFClBr,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_CHFClBr'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "CHFClBr" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_CHFClBr'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_CHFClBr'
      write(out_unit,*) ' Some parameter names of the namelist "CHFClBr" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=CHFClBr)
      STOP ' ERROR in Read_QML_CHFClBr'
    END IF
    !write(out_unit,nml=CHFClBr)
    Qmodel%norder      = norder

  END SUBROUTINE Read_QML_CHFClBr
!! === README ==
!! Polynomial potential: $V(R) = \sum_i coef(i) \cdot (r-Req)^i$
!! pot_name  = 'CHFClBr'
!! ndim      = 1
!! nsurf     = 1
!! reduced mass      = 1744.60504565084306291455 au
!! remark: Default parameters for H-F
!! === END README ==
!> @brief Subroutine wich prints the CHFClBr current parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_CHFClBr_t):   derived type with the CHFClBr parameters.
!! @param nio            integer          :   file unit to print the parameters.
  SUBROUTINE Write_QML_CHFClBr(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_CHFClBr_t),   intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'CHFClBr default parameters:'
    write(nio,*)
    write(nio,*) ' QFF potential at ??? level:'
    write(nio,*)
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) '  norder: ',QModel%norder
    write(nio,*)
    write(nio,*) 'end CHFClBr current parameters'

  END SUBROUTINE Write_QML_CHFClBr

!> @brief Subroutine wich calculates the CHFClBr potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      Potential with derivatives,.
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_CHFClBr_t):   derived type with the CHFClBr parameters.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_CHFClBr(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_CHFClBr_t), intent(in)     :: QModel
    TYPE(dnS_t),          intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),          intent(in)     :: dnQ(:)
    integer,              intent(in)     :: nderiv


  REAL(Rkind), PARAMETER :: quadratic(9) = (/ &
    0.001049484847345_Rkind, 0.001467338744296_Rkind, 0.001977082805837_Rkind, 0.003126571192830_Rkind, 0.003693818620119_Rkind,  &
    0.005038640790029_Rkind, 0.005581060837652_Rkind, 0.006042199190490_Rkind, 0.014928556478443_Rkind /)

  REAL(Rkind), PARAMETER :: cubic(152) = (/ &
    -0.000099150958196_Rkind, 0.000077709254566_Rkind, -0.000153243405810_Rkind, 0.000053752134934_Rkind, 0.000234921820704_Rkind,  &
    0.000259866526003_Rkind, 0.000181202764764_Rkind, 0.000175386466124_Rkind, -0.000247581780458_Rkind, -0.000176487367848_Rkind,  &
    0.000213337367092_Rkind, 0.000180536947494_Rkind, -0.000132708959656_Rkind, -0.000355954168881_Rkind, -0.000402641888259_Rkind,  &
    0.000330320955763_Rkind, -0.000346084782217_Rkind, -0.000099948863625_Rkind, -0.000131097520567_Rkind, -0.000135930607623_Rkind,  &
    -0.000150603465164_Rkind, 0.000415057992950_Rkind, 0.000252404068999_Rkind, 0.000070257869457_Rkind, -0.000696071154311_Rkind,  &
    -0.000084160888594_Rkind, -0.000812261257198_Rkind, -0.000201264491136_Rkind, 0.000052536823632_Rkind, 0.000860183561176_Rkind,  &
    -0.001250344143587_Rkind, 0.000804450286111_Rkind, -0.000047764700342_Rkind, -0.000254975072292_Rkind, 0.000265298178821_Rkind,  &
    -0.000177943663722_Rkind, -0.000172095516298_Rkind, 0.000283040274915_Rkind, -0.000186537094262_Rkind, -0.000677783801577_Rkind,  &
    -0.000039200248094_Rkind, 0.000239151284466_Rkind, 0.000222798779497_Rkind, -0.000202415558111_Rkind, -0.000960300280191_Rkind,  &
    0.000078378079018_Rkind, -0.000091968898063_Rkind, 0.000370778205626_Rkind, 0.000064736137891_Rkind, -0.000010971154092_Rkind,  &
    0.002081122483942_Rkind, 0.000069650692221_Rkind, -0.000051545000574_Rkind, -0.000221037664796_Rkind, -0.000078041821476_Rkind,  &
    0.000247961505438_Rkind, 0.000306919481215_Rkind, -0.000065091805421_Rkind, -0.000248875005092_Rkind, 0.000448924048236_Rkind,  &
    0.000121099736375_Rkind, -0.000408123979709_Rkind, -0.000205797543515_Rkind, 0.000064420520548_Rkind, 0.000335020496633_Rkind,  &
    0.000108134137653_Rkind, 0.000198042296408_Rkind, 0.000336245512929_Rkind, 0.000211561763244_Rkind, 0.000197329047688_Rkind,  &
    -0.000052266997458_Rkind, 0.000122315184367_Rkind, -0.001240199144203_Rkind, 0.000311740949616_Rkind, -0.000058799962070_Rkind,  &
    0.000068790456126_Rkind, -0.000088786616833_Rkind, 0.000185948506874_Rkind, 0.000031702525055_Rkind, 0.000144096198283_Rkind,  &
    0.000096606791717_Rkind, -0.000372495716200_Rkind, 0.000233426704855_Rkind, 0.000253364043274_Rkind, -0.000023851868288_Rkind,  &
    -0.000465964377575_Rkind, -0.000443010562965_Rkind, -0.000180499494418_Rkind, -0.000124852561906_Rkind, 0.000061622338379_Rkind,  &
    -0.000090942629111_Rkind, 0.000496072549798_Rkind, 0.000256325843441_Rkind, -0.000080660256219_Rkind, 0.000209189097224_Rkind,  &
    0.000123726645902_Rkind, -0.000453504258686_Rkind, -0.000178553301378_Rkind, 0.000439020853564_Rkind, 0.001508596223320_Rkind,  &
    -0.000270363365595_Rkind, 0.000131122944918_Rkind, -0.000132133767893_Rkind, 0.000050231819191_Rkind, -0.000314417887703_Rkind,  &
    -0.000502578176399_Rkind, 0.000388206324658_Rkind, 0.001217411408013_Rkind, 0.000292296105474_Rkind, -0.000436206496405_Rkind,  &
    0.000035111438401_Rkind, 0.000114198528748_Rkind, 0.000172577075371_Rkind, -0.000093620660719_Rkind, 0.000130741579657_Rkind,  &
    -0.000153227640890_Rkind, -0.000230180862747_Rkind, 0.000222499154891_Rkind, 0.000146123175147_Rkind, 0.000059820490040_Rkind,  &
    0.000194903345926_Rkind, 0.000121099189615_Rkind, 0.000189717325141_Rkind, 0.000190349106588_Rkind, -0.000178409776818_Rkind,  &
    -0.000201143520435_Rkind, 0.000064948280860_Rkind, 0.000334284830743_Rkind, -0.000576759870649_Rkind, -0.002164627579203_Rkind,  &
    -0.000438611603532_Rkind, 0.001084430571833_Rkind, -0.001272773068375_Rkind, -0.000986453598979_Rkind, 0.007411174311322_Rkind,  &
    0.000699219946477_Rkind, 0.000443046375760_Rkind, -0.001546390477481_Rkind, -0.002034681535681_Rkind, 0.000429683008971_Rkind,  &
    0.002790005916293_Rkind, 0.000229624807593_Rkind, 0.005772339117688_Rkind, 0.000008321143945_Rkind, 0.000034329981342_Rkind,  &
    -0.000131842618071_Rkind, -0.000058291748436_Rkind, 0.000023126044082_Rkind, -0.000052580154380_Rkind, -0.000102682391397_Rkind,  &
    0.000180924873877_Rkind, -0.008647136929444_Rkind /)

  REAL(Rkind), PARAMETER :: quartic(300) = (/ &
    0.000019085121473_Rkind, 0.000027348536651_Rkind, -0.000041068436674_Rkind, 0.000065562565979_Rkind, 0.000042278508190_Rkind,  &
    0.000066845812240_Rkind, 0.000019870269164_Rkind, 0.000060225639371_Rkind, 0.000045126855610_Rkind, 0.000160815305986_Rkind,  &
    0.000127016502208_Rkind, -0.000145383545245_Rkind, -0.000110701450315_Rkind, 0.000034239219144_Rkind, 0.000093959834315_Rkind,  &
    0.000078432208281_Rkind, -0.000067687640741_Rkind, 0.000055664337713_Rkind, 0.000060927132747_Rkind, -0.000020406914330_Rkind,  &
    -0.000067230184682_Rkind, -0.000040333955431_Rkind, -0.000058991601531_Rkind, 0.000034342010067_Rkind, 0.000170179668451_Rkind,  &
    -0.000055277459286_Rkind, 0.000268418812836_Rkind, -0.000255113949391_Rkind, 0.000317631060887_Rkind, -0.000208734739473_Rkind,  &
    0.000029638778566_Rkind, 0.000084576699749_Rkind, -0.000179659944085_Rkind, 0.000152459169386_Rkind, 0.000124611850715_Rkind,  &
    0.000293187598030_Rkind, 0.000086585496835_Rkind, -0.000454962468220_Rkind, 0.000164419002664_Rkind, 0.000064440249480_Rkind,  &
    0.000356647870923_Rkind, -0.000277557545579_Rkind, 0.000218247729594_Rkind, 0.000058300496600_Rkind, -0.000046157863151_Rkind,  &
    0.000051656813041_Rkind, 0.000041263447823_Rkind, 0.000037643895098_Rkind, 0.000058535238992_Rkind, -0.000078894220675_Rkind,  &
    0.000108775214023_Rkind, 0.000138296439140_Rkind, -0.000092707571135_Rkind, 0.000035386139854_Rkind, 0.000039018997077_Rkind,  &
    0.000115731644434_Rkind, 0.000172002567059_Rkind, 0.000155990146953_Rkind, -0.000063891666715_Rkind, -0.000232565557492_Rkind,  &
    0.000248972738484_Rkind, -0.000198824801425_Rkind, 0.000079159399387_Rkind, -0.000511717273649_Rkind, 0.000099131457081_Rkind,  &
    -0.000265425847335_Rkind, -0.000309543110180_Rkind, 0.000066014189930_Rkind, 0.000667853405584_Rkind, -0.000418515431266_Rkind,  &
    -0.000651144412945_Rkind, -0.000159852461220_Rkind, 0.000868306094467_Rkind, 0.000559248233992_Rkind, -0.000123442057202_Rkind,  &
    -0.000569147054582_Rkind, 0.000034698315484_Rkind, -0.000031929886184_Rkind, 0.000225216370983_Rkind, 0.000068527646708_Rkind,  &
    -0.000018834249654_Rkind, 0.000761752822894_Rkind, -0.000044902683916_Rkind, -0.000028445200983_Rkind, 0.000045568637876_Rkind,  &
    -0.000123077368128_Rkind, -0.000195865370551_Rkind, 0.000101386296271_Rkind, 0.000075825437756_Rkind, 0.000048873803469_Rkind,  &
    -0.000037806829647_Rkind, -0.000028620711017_Rkind, 0.000038721559512_Rkind, -0.000039654879225_Rkind, -0.000148088914865_Rkind,  &
    -0.000078089936377_Rkind, -0.000091383682363_Rkind, -0.000060839651110_Rkind, -0.000059256780243_Rkind, -0.000173316978652_Rkind,  &
    0.000019675713649_Rkind, 0.000015749428435_Rkind, -0.000146745524979_Rkind, 0.000008834551802_Rkind, 0.000157512327435_Rkind,  &
    0.000171093851556_Rkind, 0.000099745651073_Rkind, -0.000044492340363_Rkind, 0.000114707562522_Rkind, 0.000172843757673_Rkind,  &
    -0.000041175966186_Rkind, 0.000062472277157_Rkind, -0.000007402313378_Rkind, -0.000143538229468_Rkind, -0.000238806278761_Rkind,  &
    0.000015707874657_Rkind, 0.000372747772666_Rkind, 0.000163326028964_Rkind, 0.000107152429660_Rkind, 0.000032222220654_Rkind,  &
    -0.000297799520573_Rkind, 0.000006811265569_Rkind, 0.000136907668154_Rkind, 0.000249258147324_Rkind, -0.000040257955759_Rkind,  &
    -0.000300658529817_Rkind, 0.000037767462910_Rkind, 0.000040857204972_Rkind, -0.000230780203086_Rkind, -0.000869654405195_Rkind,  &
    -0.000149553867775_Rkind, 0.000630574381813_Rkind, -0.000318634548163_Rkind, -0.000445441367822_Rkind, 0.001559376260768_Rkind,  &
    -0.000045229099773_Rkind, 0.000039392981075_Rkind, -0.000025792320345_Rkind, -0.000013963891776_Rkind, -0.000065652781417_Rkind,  &
    -0.000056444792378_Rkind, -0.000020022906395_Rkind, 0.000154779619803_Rkind, -0.000142112278787_Rkind, -0.000064830089524_Rkind,  &
    0.000037412068760_Rkind, 0.000030456731870_Rkind, 0.000106795395229_Rkind, -0.000276912186254_Rkind, 0.000037460730421_Rkind,  &
    -0.000080406377219_Rkind, 0.000091512900031_Rkind, -0.000020421494603_Rkind, 0.000079259456509_Rkind, 0.000070981506622_Rkind,  &
    0.000020496947515_Rkind, 0.000031201966064_Rkind, 0.000022463644063_Rkind, -0.000071473590829_Rkind, 0.000018415978078_Rkind,  &
    0.000221536674632_Rkind, -0.000064034371135_Rkind, -0.000026161930261_Rkind, 0.000066181498560_Rkind, -0.000280952744356_Rkind,  &
    0.000154781806844_Rkind, -0.000120965233358_Rkind, 0.000294404686303_Rkind, 0.000146471051343_Rkind, -0.000775133230377_Rkind,  &
    -0.001323452820885_Rkind, 0.000429544131872_Rkind, 0.001042947326398_Rkind, 0.000112334258616_Rkind, 0.000061702438753_Rkind,  &
    -0.000109836475630_Rkind, -0.000063280388778_Rkind, -0.000125023834548_Rkind, 0.000012677456081_Rkind, -0.000015817226703_Rkind,  &
    0.000023013684855_Rkind, 0.000262279242209_Rkind, 0.000020057078909_Rkind, -0.000050231409121_Rkind, 0.000309067155399_Rkind,  &
    0.000008320597185_Rkind, -0.000147677751172_Rkind, -0.000048680250347_Rkind, 0.000174122903232_Rkind, 0.000130937593199_Rkind,  &
    -0.000491873704609_Rkind, -0.000840911766647_Rkind, 0.000348651502555_Rkind, 0.000341619072472_Rkind, -0.000201329692293_Rkind,  &
    -0.000522854779541_Rkind, -0.000149976331180_Rkind, 0.000238250223607_Rkind, -0.000328144713356_Rkind, -0.000183924127121_Rkind,  &
    0.002362594969460_Rkind, 0.000200594117530_Rkind, 0.000168096876480_Rkind, -0.000472563044541_Rkind, -0.000659374976945_Rkind,  &
    0.000144505083808_Rkind, 0.000625761980519_Rkind, 0.000005466508783_Rkind, 0.001023752397245_Rkind, -0.000020122963517_Rkind,  &
    -0.000027985922390_Rkind, -0.000019611195942_Rkind, -0.000026033441607_Rkind, 0.000036863668249_Rkind, 0.000024624987254_Rkind,  &
    -0.000023921306837_Rkind, -0.000021073779558_Rkind, 0.000095790205313_Rkind, 0.000022066149376_Rkind, 0.000009682394666_Rkind,  &
    0.000022085832744_Rkind, -0.000028941659272_Rkind, -0.000090013136719_Rkind, -0.000009562289668_Rkind, -0.000017274342717_Rkind,  &
    0.000097025336673_Rkind, 0.000163147238368_Rkind, -0.000052496636755_Rkind, 0.000081170383514_Rkind, 0.000031441447045_Rkind,  &
    -0.000024542973219_Rkind, 0.000052966303793_Rkind, 0.000327171480146_Rkind, -0.000209016047611_Rkind, -0.000033906788924_Rkind,  &
    -0.000076516360434_Rkind, 0.000134815763513_Rkind, -0.000268434122122_Rkind, 0.000868670783541_Rkind, -0.000208537632410_Rkind,  &
    -0.000199609949115_Rkind, 0.000055277459286_Rkind, -0.000143448560790_Rkind, -0.000032174652514_Rkind, -0.000250645278028_Rkind,  &
    -0.000178751228582_Rkind, 0.000269721195704_Rkind, 0.000265623774538_Rkind, 0.000158264669511_Rkind, 0.000032629010265_Rkind,  &
    0.000249985338430_Rkind, -0.000289764878988_Rkind, -0.001118373993692_Rkind, -0.000349508822596_Rkind, 0.000646905198626_Rkind,  &
    -0.000034935518297_Rkind, -0.000165699515123_Rkind, -0.000216035082068_Rkind, 0.000107222414969_Rkind, -0.000017452039792_Rkind,  &
    -0.000121260210503_Rkind, 0.000271977675175_Rkind, 0.000387996915490_Rkind, -0.000340264200622_Rkind, -0.000462521428404_Rkind,  &
    -0.000101806754888_Rkind, -0.000222498608131_Rkind, 0.000062496334608_Rkind, 0.000333866012406_Rkind, -0.000081110513269_Rkind,  &
    -0.000290220877020_Rkind, -0.000365268639475_Rkind, 0.000396252994968_Rkind, 0.000870844702217_Rkind, -0.000272704319521_Rkind,  &
    -0.000376715064897_Rkind, 0.001010939982514_Rkind, 0.003797357511416_Rkind, 0.000826279916123_Rkind, -0.002393752921326_Rkind,  &
    0.001264514255096_Rkind, 0.001902606407823_Rkind, -0.015675889183712_Rkind, -0.001348496079715_Rkind, -0.000730195007046_Rkind,  &
    0.002644875521035_Rkind, 0.004092073304289_Rkind, -0.001101153233478_Rkind, -0.004148630182515_Rkind, -0.000405638681082_Rkind,  &
    -0.013618671284857_Rkind, -0.000029163279419_Rkind, -0.000041101606795_Rkind, 0.000137352184222_Rkind, 0.000217757103414_Rkind,  &
    -0.000086101067271_Rkind, -0.000038091873980_Rkind, 0.000105679457599_Rkind, -0.000177108396343_Rkind, 0.004714184065553_Rkind /)

  INTEGER, PARAMETER :: idx3(3, 152) = RESHAPE((/ &
    1 , 1 , 1 , 2 , 2 , 1 , 2 , 2 , 2 , 3, & 
    1 , 1 , 3 , 2 , 1 , 3 , 2 , 2 , 3 , 3, & 
    2 , 3 , 3 , 3 , 4 , 2 , 1 , 4 , 2 , 2, & 
    4 , 3 , 1 , 4 , 3 , 2 , 4 , 3 , 3 , 4, & 
    4 , 1 , 4 , 4 , 2 , 4 , 4 , 3 , 4 , 4, & 
    4 , 5 , 1 , 1 , 5 , 2 , 1 , 5 , 2 , 2, & 
    5 , 3 , 1 , 5 , 3 , 2 , 5 , 3 , 3 , 5, & 
    4 , 1 , 5 , 4 , 2 , 5 , 4 , 3 , 5 , 4, & 
    4 , 5 , 5 , 1 , 5 , 5 , 2 , 5 , 5 , 3, & 
    5 , 5 , 4 , 5 , 5 , 5 , 6 , 1 , 1 , 6, & 
    2 , 2 , 6 , 3 , 1 , 6 , 3 , 2 , 6 , 4, & 
    1 , 6 , 4 , 2 , 6 , 4 , 3 , 6 , 4 , 4, & 
    6 , 5 , 1 , 6 , 5 , 2 , 6 , 5 , 3 , 6, & 
    5 , 4 , 6 , 5 , 5 , 6 , 6 , 1 , 6 , 6, & 
    2 , 6 , 6 , 3 , 6 , 6 , 4 , 6 , 6 , 5, & 
    6 , 6 , 6 , 7 , 2 , 2 , 7 , 3 , 1 , 7, & 
    3 , 2 , 7 , 3 , 3 , 7 , 4 , 1 , 7 , 4, & 
    2 , 7 , 4 , 3 , 7 , 4 , 4 , 7 , 5 , 1, & 
    7 , 5 , 2 , 7 , 5 , 3 , 7 , 5 , 4 , 7, & 
    5 , 5 , 7 , 6 , 2 , 7 , 6 , 3 , 7 , 6, & 
    4 , 7 , 6 , 5 , 7 , 6 , 6 , 7 , 7 , 1, & 
    7 , 7 , 2 , 7 , 7 , 3 , 7 , 7 , 4 , 7, & 
    7 , 5 , 7 , 7 , 6 , 7 , 7 , 7 , 8 , 1, & 
    1 , 8 , 2 , 2 , 8 , 3 , 1 , 8 , 3 , 2, & 
    8 , 4 , 1 , 8 , 4 , 2 , 8 , 4 , 3 , 8, & 
    4 , 4 , 8 , 5 , 1 , 8 , 5 , 2 , 8 , 5, & 
    3 , 8 , 5 , 4 , 8 , 5 , 5 , 8 , 6 , 1, & 
    8 , 6 , 2 , 8 , 6 , 3 , 8 , 6 , 4 , 8, & 
    6 , 5 , 8 , 6 , 6 , 8 , 7 , 1 , 8 , 7, & 
    2 , 8 , 7 , 3 , 8 , 7 , 4 , 8 , 7 , 5, & 
    8 , 7 , 6 , 8 , 7 , 7 , 8 , 8 , 1 , 8, & 
    8 , 2 , 8 , 8 , 3 , 8 , 8 , 4 , 8 , 8, & 
    5 , 8 , 8 , 6 , 8 , 8 , 7 , 8 , 8 , 8, & 
    9 , 1 , 1 , 9 , 2 , 1 , 9 , 2 , 2 , 9, & 
    3 , 1 , 9 , 3 , 3 , 9 , 4 , 1 , 9 , 4, & 
    2 , 9 , 4 , 3 , 9 , 4 , 4 , 9 , 5 , 1, & 
    9 , 5 , 2 , 9 , 5 , 5 , 9 , 6 , 1 , 9, & 
    6 , 2 , 9 , 6 , 3 , 9 , 6 , 4 , 9 , 6, & 
    5 , 9 , 6 , 6 , 9 , 7 , 1 , 9 , 7 , 2, & 
    9 , 7 , 3 , 9 , 7 , 4 , 9 , 7 , 5 , 9, & 
    7 , 6 , 9 , 7 , 7 , 9 , 8 , 1 , 9 , 8, & 
    2 , 9 , 8 , 3 , 9 , 8 , 4 , 9 , 8 , 5, & 
    9 , 8 , 6 , 9 , 8 , 7 , 9 , 8 , 8 , 9, & 
    9 , 1 , 9 , 9 , 2 , 9 , 9 , 3 , 9 , 9, & 
    4 , 9 , 9 , 5 , 9 , 9 , 6 , 9 , 9 , 7, & 
    9 , 9 , 8 , 9 , 9 , 9 /), (/3, 152/))

  INTEGER, PARAMETER :: idx4(4, 300) = RESHAPE((/ &
    1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 3 , 2, & 
    2 , 2 , 3 , 3 , 2 , 1 , 3 , 3 , 2 , 2, & 
    3 , 3 , 3 , 2 , 3 , 3 , 3 , 3 , 4 , 2, & 
    2 , 1 , 4 , 3 , 3 , 1 , 4 , 4 , 2 , 1, & 
    4 , 4 , 2 , 2 , 4 , 4 , 3 , 1 , 4 , 4, & 
    3 , 2 , 4 , 4 , 3 , 3 , 4 , 4 , 4 , 1, & 
    4 , 4 , 4 , 2 , 4 , 4 , 4 , 3 , 4 , 4, & 
    4 , 4 , 5 , 2 , 2 , 1 , 5 , 2 , 2 , 2, & 
    5 , 3 , 1 , 1 , 5 , 3 , 2 , 2 , 5 , 3, & 
    3 , 1 , 5 , 3 , 3 , 2 , 5 , 4 , 2 , 2, & 
    5 , 4 , 3 , 3 , 5 , 4 , 4 , 2 , 5 , 4, & 
    4 , 3 , 5 , 4 , 4 , 4 , 5 , 5 , 1 , 1, & 
    5 , 5 , 2 , 1 , 5 , 5 , 2 , 2 , 5 , 5, & 
    3 , 1 , 5 , 5 , 3 , 2 , 5 , 5 , 3 , 3, & 
    5 , 5 , 4 , 1 , 5 , 5 , 4 , 2 , 5 , 5, & 
    4 , 3 , 5 , 5 , 4 , 4 , 5 , 5 , 5 , 2, & 
    5 , 5 , 5 , 3 , 5 , 5 , 5 , 4 , 5 , 5, & 
    5 , 5 , 6 , 2 , 2 , 1 , 6 , 2 , 2 , 2, & 
    6 , 3 , 2 , 2 , 6 , 3 , 3 , 1 , 6 , 3, & 
    3 , 2 , 6 , 3 , 3 , 3 , 6 , 4 , 3 , 3, & 
    6 , 4 , 4 , 1 , 6 , 4 , 4 , 2 , 6 , 4, & 
    4 , 3 , 6 , 4 , 4 , 4 , 6 , 5 , 1 , 1, & 
    6 , 5 , 3 , 3 , 6 , 5 , 4 , 4 , 6 , 5, & 
    5 , 1 , 6 , 5 , 5 , 2 , 6 , 5 , 5 , 3, & 
    6 , 5 , 5 , 4 , 6 , 5 , 5 , 5 , 6 , 6, & 
    2 , 1 , 6 , 6 , 2 , 2 , 6 , 6 , 3 , 1, & 
    6 , 6 , 3 , 2 , 6 , 6 , 3 , 3 , 6 , 6, & 
    4 , 1 , 6 , 6 , 4 , 2 , 6 , 6 , 4 , 3, & 
    6 , 6 , 4 , 4 , 6 , 6 , 5 , 1 , 6 , 6, & 
    5 , 2 , 6 , 6 , 5 , 3 , 6 , 6 , 5 , 4, & 
    6 , 6 , 5 , 5 , 6 , 6 , 6 , 1 , 6 , 6, & 
    6 , 2 , 6 , 6 , 6 , 3 , 6 , 6 , 6 , 4, & 
    6 , 6 , 6 , 5 , 6 , 6 , 6 , 6 , 7 , 2, & 
    2 , 2 , 7 , 3 , 2 , 2 , 7 , 4 , 1 , 1, & 
    7 , 4 , 4 , 1 , 7 , 4 , 4 , 2 , 7 , 4, & 
    4 , 3 , 7 , 4 , 4 , 4 , 7 , 5 , 1 , 1, & 
    7 , 5 , 3 , 3 , 7 , 5 , 4 , 4 , 7 , 5, & 
    5 , 1 , 7 , 5 , 5 , 2 , 7 , 5 , 5 , 3, & 
    7 , 5 , 5 , 4 , 7 , 5 , 5 , 5 , 7 , 6, & 
    2 , 2 , 7 , 6 , 3 , 3 , 7 , 6 , 4 , 4, & 
    7 , 6 , 5 , 5 , 7 , 6 , 6 , 1 , 7 , 6, & 
    6 , 2 , 7 , 6 , 6 , 3 , 7 , 6 , 6 , 4, & 
    7 , 6 , 6 , 5 , 7 , 6 , 6 , 6 , 7 , 7, & 
    1 , 1 , 7 , 7 , 2 , 1 , 7 , 7 , 2 , 2, & 
    7 , 7 , 3 , 1 , 7 , 7 , 3 , 2 , 7 , 7, & 
    3 , 3 , 7 , 7 , 4 , 1 , 7 , 7 , 4 , 2, & 
    7 , 7 , 4 , 3 , 7 , 7 , 4 , 4 , 7 , 7, & 
    5 , 1 , 7 , 7 , 5 , 2 , 7 , 7 , 5 , 3, & 
    7 , 7 , 5 , 4 , 7 , 7 , 5 , 5 , 7 , 7, & 
    6 , 1 , 7 , 7 , 6 , 2 , 7 , 7 , 6 , 3, & 
    7 , 7 , 6 , 4 , 7 , 7 , 6 , 5 , 7 , 7, & 
    6 , 6 , 7 , 7 , 7 , 1 , 7 , 7 , 7 , 2, & 
    7 , 7 , 7 , 3 , 7 , 7 , 7 , 4 , 7 , 7, & 
    7 , 5 , 7 , 7 , 7 , 6 , 7 , 7 , 7 , 7, & 
    8 , 3 , 2 , 2 , 8 , 3 , 3 , 1 , 8 , 3, & 
    3 , 2 , 8 , 3 , 3 , 3 , 8 , 4 , 2 , 2, & 
    8 , 4 , 3 , 3 , 8 , 4 , 4 , 1 , 8 , 4, & 
    4 , 2 , 8 , 4 , 4 , 3 , 8 , 4 , 4 , 4, & 
    8 , 5 , 2 , 2 , 8 , 5 , 3 , 3 , 8 , 5, & 
    4 , 4 , 8 , 5 , 5 , 3 , 8 , 5 , 5 , 4, & 
    8 , 5 , 5 , 5 , 8 , 6 , 2 , 2 , 8 , 6, & 
    3 , 3 , 8 , 6 , 4 , 4 , 8 , 6 , 5 , 5, & 
    8 , 6 , 6 , 1 , 8 , 6 , 6 , 2 , 8 , 6, & 
    6 , 3 , 8 , 6 , 6 , 4 , 8 , 6 , 6 , 5, & 
    8 , 6 , 6 , 6 , 8 , 7 , 1 , 1 , 8 , 7, & 
    2 , 2 , 8 , 7 , 3 , 3 , 8 , 7 , 4 , 4, & 
    8 , 7 , 5 , 5 , 8 , 7 , 6 , 6 , 8 , 7, & 
    7 , 1 , 8 , 7 , 7 , 2 , 8 , 7 , 7 , 3, & 
    8 , 7 , 7 , 4 , 8 , 7 , 7 , 5 , 8 , 7, & 
    7 , 6 , 8 , 7 , 7 , 7 , 8 , 8 , 2 , 1, & 
    8 , 8 , 2 , 2 , 8 , 8 , 3 , 1 , 8 , 8, & 
    3 , 2 , 8 , 8 , 3 , 3 , 8 , 8 , 4 , 1, & 
    8 , 8 , 4 , 2 , 8 , 8 , 4 , 3 , 8 , 8, & 
    4 , 4 , 8 , 8 , 5 , 1 , 8 , 8 , 5 , 2, & 
    8 , 8 , 5 , 3 , 8 , 8 , 5 , 4 , 8 , 8, & 
    5 , 5 , 8 , 8 , 6 , 1 , 8 , 8 , 6 , 2, & 
    8 , 8 , 6 , 3 , 8 , 8 , 6 , 4 , 8 , 8, & 
    6 , 5 , 8 , 8 , 6 , 6 , 8 , 8 , 7 , 1, & 
    8 , 8 , 7 , 2 , 8 , 8 , 7 , 3 , 8 , 8, & 
    7 , 4 , 8 , 8 , 7 , 5 , 8 , 8 , 7 , 6, & 
    8 , 8 , 7 , 7 , 8 , 8 , 8 , 1 , 8 , 8, & 
    8 , 2 , 8 , 8 , 8 , 3 , 8 , 8 , 8 , 4, & 
    8 , 8 , 8 , 5 , 8 , 8 , 8 , 6 , 8 , 8, & 
    8 , 7 , 8 , 8 , 8 , 8 , 9 , 4 , 2 , 2, & 
    9 , 4 , 3 , 3 , 9 , 4 , 4 , 2 , 9 , 4, & 
    4 , 3 , 9 , 4 , 4 , 4 , 9 , 5 , 3 , 3, & 
    9 , 5 , 4 , 4 , 9 , 5 , 5 , 1 , 9 , 5, & 
    5 , 3 , 9 , 5 , 5 , 4 , 9 , 5 , 5 , 5, & 
    9 , 6 , 3 , 3 , 9 , 6 , 4 , 4 , 9 , 6, & 
    5 , 5 , 9 , 6 , 6 , 1 , 9 , 6 , 6 , 2, & 
    9 , 6 , 6 , 3 , 9 , 6 , 6 , 4 , 9 , 6, & 
    6 , 5 , 9 , 6 , 6 , 6 , 9 , 7 , 1 , 1, & 
    9 , 7 , 2 , 2 , 9 , 7 , 3 , 3 , 9 , 7, & 
    4 , 4 , 9 , 7 , 5 , 5 , 9 , 7 , 6 , 6, & 
    9 , 7 , 7 , 1 , 9 , 7 , 7 , 2 , 9 , 7, & 
    7 , 3 , 9 , 7 , 7 , 4 , 9 , 7 , 7 , 5, & 
    9 , 7 , 7 , 7 , 9 , 8 , 1 , 1 , 9 , 8, & 
    2 , 2 , 9 , 8 , 3 , 3 , 9 , 8 , 4 , 4, & 
    9 , 8 , 5 , 5 , 9 , 8 , 6 , 6 , 9 , 8, & 
    7 , 7 , 9 , 8 , 8 , 1 , 9 , 8 , 8 , 2, & 
    9 , 8 , 8 , 4 , 9 , 8 , 8 , 5 , 9 , 8, & 
    8 , 6 , 9 , 8 , 8 , 7 , 9 , 8 , 8 , 8, & 
    9 , 9 , 1 , 1 , 9 , 9 , 2 , 1 , 9 , 9, & 
    2 , 2 , 9 , 9 , 3 , 1 , 9 , 9 , 3 , 2, & 
    9 , 9 , 3 , 3 , 9 , 9 , 4 , 1 , 9 , 9, & 
    4 , 2 , 9 , 9 , 4 , 3 , 9 , 9 , 4 , 4, & 
    9 , 9 , 5 , 1 , 9 , 9 , 5 , 2 , 9 , 9, & 
    5 , 3 , 9 , 9 , 5 , 4 , 9 , 9 , 5 , 5, & 
    9 , 9 , 6 , 1 , 9 , 9 , 6 , 2 , 9 , 9, & 
    6 , 3 , 9 , 9 , 6 , 4 , 9 , 9 , 6 , 5, & 
    9 , 9 , 6 , 6 , 9 , 9 , 7 , 1 , 9 , 9, & 
    7 , 2 , 9 , 9 , 7 , 3 , 9 , 9 , 7 , 4, & 
    9 , 9 , 7 , 5 , 9 , 9 , 7 , 6 , 9 , 9, & 
    7 , 7 , 9 , 9 , 8 , 1 , 9 , 9 , 8 , 2, & 
    9 , 9 , 8 , 3 , 9 , 9 , 8 , 4 , 9 , 9, & 
    8 , 5 , 9 , 9 , 8 , 6 , 9 , 9 , 8 , 7, & 
    9 , 9 , 8 , 8 , 9 , 9 , 9 , 1 , 9 , 9, & 
    9 , 2 , 9 , 9 , 9 , 3 , 9 , 9 , 9 , 4, & 
    9 , 9 , 9 , 5 , 9 , 9 , 9 , 6 , 9 , 9, & 
    9 , 7 , 9 , 9 , 9 , 8 , 9 , 9 , 9 , 9 /), (/4, 300/))



    integer          :: i
    TYPE (dnS_t)     :: dnDR
    TYPE (dnS_t)     :: E2,E3,E4

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_CHFClBr'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) ' nderiv:',nderiv
      write(out_unit,*) ' Q(:):',(get_d0(dnQ(i)),i=1,size(dnQ))
    END IF


    E2 = ZERO
    DO i = 1, size(quadratic)
      E2 = E2 + HALF * quadratic(i) * dnQ(i)**2
    END DO

    E3 = ZERO
    IF (QModel%norder >= 3) THEN
      DO i = 1,size(cubic)
        E3 = E3 + cubic(i)/SIX * dnQ(idx3(1,i)) * dnQ(idx3(2,i)) * dnQ(idx3(3,i))
      END DO
    END IF

    E4 = ZERO
    IF (QModel%norder >= 4) THEN
      DO i = 1, size(quartic)
        E4 = E4 + quartic(i)/24._Rkind * dnQ(idx4(1,i)) * dnQ(idx4(2,i)) * dnQ(idx4(3,i)) * dnQ(idx4(4,i))
      END DO
    END IF

    Mat_OF_PotDia(1,1) = E2 + E3 + E4

    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),out_unit)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
  END SUBROUTINE EvalPot_QML_CHFClBr

   SUBROUTINE RefValues_QML_CHFClBr(QModel,err,nderiv,Q0,dnMatV,d0GGdef,option)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_CHFClBr_t), intent(in)              :: QModel
    integer,              intent(inout)           :: err
    integer,              intent(in)              :: nderiv
    real (kind=Rkind),    intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),       intent(inout), optional :: dnMatV
    real (kind=Rkind),    intent(inout), optional :: d0GGdef(:,:)
    integer,              intent(in),    optional :: option

    real (kind=Rkind), allocatable :: d0(:,:),d1(:,:,:),d2(:,:,:,:),d3(:,:,:,:,:),V(:)
    real (kind=Rkind), allocatable :: masses(:)
    integer        :: i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_CHFClBr'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      flush(out_unit)
    END IF

    IF (.NOT. QModel%Init) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'The model is not initialized!'
      err = -1
      RETURN
    ELSE
      err = 0
    END IF

    IF (present(Q0)) THEN
      IF (size(Q0) /= QModel%ndim) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'incompatible Q0 size:'
        write(out_unit,*) 'size(Q0), ndimQ:',size(Q0),QModel%ndim
        err = 1
        Q0(:) = HUGE(ONE)
        RETURN
      END IF
      Q0(:) = [0.1_Rkind,-0.1_Rkind,0.2_Rkind,-0.2_Rkind,0.3_Rkind,-0.3_Rkind,0.4_Rkind,-0.4_Rkind,0.5_Rkind]
    END IF

    IF (present(dnMatV)) THEN
      err = 0

      IF (nderiv >= 0) THEN ! no derivative
        V  = [3.2766587723962216E-003_Rkind]
        d0 = reshape(V,shape=[QModel%nsurf,QModel%nsurf])
      END IF

      SELECT CASE (nderiv)
      CASE(0)
        CALL set_dnMat(dnMatV,d0=d0)
      CASE Default
        STOP 'ERROR in RefValues_QML_CHFClBr: nderiv = 0'
      END SELECT

    END IF

    IF (present(d0GGdef)) THEN 
      d0GGdef = Identity_Mat(QModel%ndim)
      masses  = [952.84844038464450_Rkind,681.50589213793319_Rkind,505.79570923770643_Rkind, &
                 319.83919070617901_Rkind,270.72255106228903_Rkind,198.46622168004251_Rkind, &
                 179.17740535161565_Rkind,165.50265366523007_Rkind,66.985713015455389_Rkind]
      DO i=1,QModel%ndim
        d0GGdef(i,i) = ONE/masses(i)
      END DO
    END IF


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_CHFClBr
END MODULE QML_CHFClBr_m
