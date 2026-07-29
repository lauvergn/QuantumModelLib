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
     integer           :: norder = 4
     logical           :: bQFF   = .TRUE.

     real (kind=Rkind) ::  quadratic(9) = [0.001049576065177_Rkind, 0.001461717593495_Rkind,  &
                                           0.001969821921105_Rkind, 0.003125330502741_Rkind,  &
                                           0.003688471806263_Rkind,  0.004989751221640_Rkind, &
                                           0.005689101160362_Rkind, 0.006101642415972_Rkind,  &
                                           0.014536840636557_Rkind ]
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
    QModel%masses = ONE/QModel%quadratic

    QModel%d0GGdef = Identity_Mat(QModel%ndim)
    DO i=1,QModel%ndim
      QModel%d0GGdef(i,i) = QModel%quadratic(i)
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
    logical                       :: bQFF

    namelist /CHFClBr/ norder,bQFF

    write(out_unit,*) 'read CHFClBr namelist ...'
    flush(out_unit)

    norder = 4
    bQFF   = .TRUE.

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
    Qmodel%bQFF        = bQFF


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
    write(nio,*) ' QFF potential at MP2/aug-cc-pVTZ level of theory:'
    write(nio,*)
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) '  norder: ',QModel%norder
    write(nio,*) '  bQFF:   ',QModel%bQFF
    write(nio,*)
    write(nio,*) ' at Q={0.1,-0.1,0.2,-0.2,0.3,-0.3,0.4,-0.4,0.5}, ' 
    write(nio,*) ' ... the bQFF potential value is: 3.2310845043165954E-003 Hartree'
    write(nio,*)
    write(nio,*) 'The ZPE (bQFF) is 4627.47 cm-1 (optained with ElVibRot and MCTDH).'
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

    !QFF parameters
    ! quartic(9) is defined in the derived type
    real(kind=Rkind), parameter :: cubic(151) = [ &
    -0.000100764037565_Rkind, 0.000076395253043_Rkind, -0.000146512012797_Rkind, &
    0.000055135575007_Rkind, 0.000232670808836_Rkind,  &
    0.000255456221295_Rkind, 0.000181748021404_Rkind, 0.000174069411856_Rkind, &
    -0.000240855536104_Rkind, -0.000176277548610_Rkind,  &
    0.000211792496061_Rkind, 0.000185594206244_Rkind, -0.000136822920319_Rkind, &
    -0.000353110468923_Rkind, -0.000406044923933_Rkind,  &
    0.000336308527045_Rkind, -0.000346294464765_Rkind, -0.000096799251318_Rkind, &
    -0.000119766188173_Rkind, -0.000140616889557_Rkind,  &
    -0.000147459593840_Rkind, 0.000411424497839_Rkind, 0.000251851567786_Rkind, &
    0.000068597905398_Rkind, -0.000689947166351_Rkind,  &
    -0.000082862606427_Rkind, -0.000806671453983_Rkind, -0.000199583476808_Rkind, &
    0.000048586754348_Rkind, 0.000854253673098_Rkind,  &
    -0.001272208675128_Rkind, 0.000817619461892_Rkind, -0.000047107084475_Rkind, &
    -0.000246088988339_Rkind, 0.000260892384885_Rkind,  &
    -0.000180889607843_Rkind, -0.000170018100803_Rkind, 0.000277294098274_Rkind, &
    -0.000180459034161_Rkind, -0.000676318894230_Rkind,  &
    -0.000029910245020_Rkind, 0.000243838113161_Rkind, 0.000220018230346_Rkind, &
    -0.000204307621888_Rkind, -0.000947182408745_Rkind,  &
    0.000070094798219_Rkind, -0.000084118377986_Rkind, 0.000362168372279_Rkind, &
    0.000060375041604_Rkind, -0.000022887930002_Rkind,  &
    0.002091500676564_Rkind, 0.000073655300838_Rkind, -0.000055445861437_Rkind, &
    -0.000219832605248_Rkind, -0.000082445428371_Rkind,  &
    0.000232076207086_Rkind, 0.000283834170770_Rkind, -0.000062826851167_Rkind, &
    -0.000214426240090_Rkind, 0.000450298056695_Rkind,  &
    0.000114537793471_Rkind, -0.000402092667608_Rkind, -0.000209832907395_Rkind, &
    0.000065364092016_Rkind, 0.000339815037092_Rkind,  &
    0.000095149949083_Rkind, 0.000142180623816_Rkind, 0.000330731982766_Rkind, &
    0.000140205999244_Rkind, 0.000175075632933_Rkind,  &
    -0.000060530321509_Rkind, 0.000126944876618_Rkind, -0.001197720658458_Rkind, &
    0.000308096109230_Rkind, -0.000063930486692_Rkind,  &
    0.000087465507426_Rkind, -0.000086470540497_Rkind, 0.000184893943079_Rkind, &
    0.000034972424613_Rkind, 0.000139912115620_Rkind,  &
    0.000105176984947_Rkind, -0.000352164437035_Rkind, 0.000225536954731_Rkind, &
    0.000235074366809_Rkind, -0.000454977504126_Rkind,  &
    -0.000436117830121_Rkind, -0.000189082809894_Rkind, -0.000121848661201_Rkind, &
    0.000065751745019_Rkind, -0.000107231163133_Rkind,  &
    0.000521311002030_Rkind, 0.000296584892721_Rkind, -0.000074223794788_Rkind, &
    0.000227178055562_Rkind, 0.000126329771358_Rkind,  &
    -0.000431168009756_Rkind, -0.000168475963573_Rkind, 0.000383513208221_Rkind, &
    0.001498991559733_Rkind, -0.000301765992664_Rkind,  &
    0.000141697834532_Rkind, -0.000132247630711_Rkind, 0.000038874105616_Rkind, &
    -0.000288635545732_Rkind, -0.000482750600096_Rkind,  &
    0.000395063381397_Rkind, 0.001243473827913_Rkind, 0.000286576446705_Rkind, &
    -0.000352454265520_Rkind, 0.000033062044368_Rkind,  &
    0.000107000977076_Rkind, 0.000163169928918_Rkind, -0.000087084415546_Rkind, &
    0.000122376968273_Rkind, -0.000138880652446_Rkind,  &
    -0.000216978516846_Rkind, 0.000199340851955_Rkind, 0.000128540596350_Rkind, &
    0.000051645604457_Rkind, 0.000185213934504_Rkind,  &
    0.000020678471911_Rkind, 0.000121917689680_Rkind, 0.000166151685834_Rkind, &
    0.000152346263398_Rkind, -0.000155943945714_Rkind,  &
    -0.000138982623229_Rkind, 0.000031342210063_Rkind, 0.000286572346003_Rkind, &
    -0.000567522083177_Rkind, -0.002114661075414_Rkind,  &
    -0.000408379043356_Rkind, 0.001104340298861_Rkind, -0.001285018037114_Rkind, &
    -0.000832330091452_Rkind, 0.007495931988708_Rkind,  &
    0.000666470375581_Rkind, 0.000435761069071_Rkind, -0.001511182672589_Rkind, &
    -0.001971072544005_Rkind, 0.000358073821604_Rkind,  &
    0.002690976429996_Rkind, 0.000173636833387_Rkind, 0.005967756418230_Rkind, &
    0.000027069278863_Rkind, -0.000123704092042_Rkind,  &
    -0.000070803536167_Rkind, 0.000027150882828_Rkind, -0.000062100753581_Rkind, &
    -0.000090505767687_Rkind, 0.000182754470298_Rkind,  &
    -0.009055543812009_Rkind ]
    real(kind=Rkind), parameter :: quartic(301) = [ &
    0.000017765196714_Rkind, 0.000027769952099_Rkind, -0.000042812419555_Rkind,&
     0.000066614532662_Rkind, 0.000040568242190_Rkind,  &
    0.000068268846866_Rkind, 0.000019124670463_Rkind, 0.000068159130313_Rkind, &
    -0.000044108241301_Rkind, 0.000052617470766_Rkind,  &
    -0.000011253965821_Rkind, 0.000159967827629_Rkind, 0.000115961010350_Rkind,&
     -0.000144166456972_Rkind, -0.000114169550456_Rkind,  &
    0.000023476790770_Rkind, 0.000100608074209_Rkind, 0.000087741165709_Rkind, &
    -0.000079979904239_Rkind, 0.000052277978226_Rkind,  &
    0.000053749811203_Rkind, -0.000020818624784_Rkind, -0.000062430723380_Rkind,&
     -0.000037121739078_Rkind, -0.000066446130511_Rkind,  &
    0.000029733914846_Rkind, 0.000172792088831_Rkind, -0.000052186076944_Rkind, &
    0.000278058195696_Rkind, -0.000255380221623_Rkind,  &
    0.000327765626263_Rkind, -0.000215986420408_Rkind, 0.000080599019074_Rkind, &
    -0.000178576812068_Rkind, 0.000145196553246_Rkind,  &
    0.000116233570325_Rkind, 0.000299453470269_Rkind, 0.000099123255677_Rkind, &
    -0.000472313904130_Rkind, 0.000172107545023_Rkind,  &
    0.000062094648092_Rkind, 0.000357142688932_Rkind, -0.000284514522749_Rkind,&
     0.000222360505609_Rkind, 0.000055470465648_Rkind,  &
    -0.000047471727985_Rkind, 0.000055164826679_Rkind, 0.000037129393721_Rkind,&
     0.000033999738163_Rkind, 0.000059822859335_Rkind,  &
    -0.000076824186443_Rkind, 0.000108578380340_Rkind, 0.000130725996990_Rkind,&
     -0.000087706355308_Rkind, 0.000046997869119_Rkind,  &
    0.000110051352401_Rkind, 0.000175205488488_Rkind, 0.000159530419445_Rkind, &
    -0.000064558714196_Rkind, -0.000225118136395_Rkind,  &
    0.000253467107577_Rkind, -0.000195555357501_Rkind, 0.000083418661581_Rkind,&
     -0.000518583762001_Rkind, 0.000096181138878_Rkind,  &
    -0.000269361427473_Rkind, -0.000314502772229_Rkind, 0.000063775753547_Rkind,&
     0.000666922819672_Rkind, -0.000422757743893_Rkind,  &
    -0.000654939202323_Rkind, -0.000162161976433_Rkind, 0.000871744122796_Rkind,&
     0.000559210507536_Rkind, -0.000119464376526_Rkind,  &
    -0.000575029647900_Rkind, 0.000031159136513_Rkind, -0.000033899498787_Rkind,&
     0.000235035820212_Rkind, 0.000072852884637_Rkind,  &
    -0.000029860580966_Rkind, 0.000783781838137_Rkind, -0.000041888394766_Rkind,&
     -0.000029588476625_Rkind, -0.000120058158136_Rkind,  &
    -0.000194769663050_Rkind, 0.000101581489673_Rkind, 0.000051260047367_Rkind, &
    -0.000049050407023_Rkind, -0.000036453051316_Rkind,  &
    0.000037404414117_Rkind, -0.000035710550924_Rkind, -0.000148588653715_Rkind,&
     -0.000091695517948_Rkind, -0.000094648569952_Rkind,  &
    -0.000026710877532_Rkind, -0.000039982935363_Rkind, -0.000120614213291_Rkind,&
     0.000045242768779_Rkind, 0.000017304414530_Rkind,  &
    -0.000143113943529_Rkind, -0.000010296041897_Rkind, 0.000140886989111_Rkind, &
    0.000176145916084_Rkind, 0.000048557776056_Rkind,  &
    -0.000077241364499_Rkind, 0.000099017184193_Rkind, 0.000146671712348_Rkind, &
    -0.000039145298691_Rkind, 0.000053277957124_Rkind,  &
    -0.000034118931893_Rkind, -0.000112480608104_Rkind, -0.000207110041448_Rkind,&
     -0.000013063195423_Rkind, 0.000299426679018_Rkind,  &
    0.000154669174237_Rkind, 0.000083495754774_Rkind, 0.000042732045802_Rkind, &
    -0.000280103078958_Rkind, 0.000006882891160_Rkind,  &
    0.000104970309580_Rkind, 0.000196398461776_Rkind, -0.000020048057366_Rkind,&
     -0.000212015027475_Rkind, 0.000012368263170_Rkind,  &
    -0.000212166480058_Rkind, -0.000812351928269_Rkind, -0.000145266174049_Rkind,&
     0.000554581817680_Rkind, -0.000293500162628_Rkind,  &
    -0.000325717462440_Rkind, 0.001513137158159_Rkind, 0.000011794347182_Rkind, &
    -0.000042066638601_Rkind, 0.000036305972814_Rkind,  &
    -0.000008897064721_Rkind, 0.000052988720962_Rkind, -0.000053402071696_Rkind,&
     -0.000040149150473_Rkind, -0.000023904904030_Rkind,  &
    0.000152346536778_Rkind, -0.000133857292829_Rkind, -0.000050430247591_Rkind,&
     0.000028392165241_Rkind, 0.000023749624125_Rkind,  &
    0.000104952813253_Rkind, -0.000273656775842_Rkind, 0.000045163488546_Rkind, &
    -0.000083851331177_Rkind, 0.000085392466013_Rkind,  &
    -0.000029827410845_Rkind, 0.000041415447167_Rkind, 0.000049927410433_Rkind, &
    0.000016525281202_Rkind, 0.000022227990404_Rkind,  &
    0.000044116989465_Rkind, -0.000051439202470_Rkind, 0.000015324595736_Rkind, &
    0.000273019071160_Rkind, -0.000032628463505_Rkind,  &
    0.000068909285349_Rkind, -0.000239291255085_Rkind, 0.000163034059001_Rkind, &
    -0.000121172455486_Rkind, 0.000265729299263_Rkind,  &
    0.000145648177197_Rkind, -0.000718300784995_Rkind, -0.001213329842894_Rkind,&
     0.000407231393633_Rkind, 0.000925066914259_Rkind,  &
    0.000121336665808_Rkind, -0.000022508478402_Rkind, 0.000058493502961_Rkind, &
    -0.000120681464799_Rkind, -0.000052068523495_Rkind,  &
    -0.000120209610720_Rkind, 0.000029002349658_Rkind, 0.000227044372686_Rkind, &
    -0.000008692940902_Rkind, -0.000062810174980_Rkind,  &
    0.000294917000638_Rkind, 0.000017201623606_Rkind, -0.000122018840322_Rkind, &
    -0.000054934367242_Rkind, 0.000140809895918_Rkind,  &
    0.000105918938580_Rkind, -0.000458594596430_Rkind, -0.000773696891252_Rkind,&
     0.000323781019958_Rkind, 0.000314155032723_Rkind,  &
    -0.000196324102385_Rkind, -0.000509081342578_Rkind, -0.000150323523926_Rkind,&
     0.000221162326128_Rkind, -0.000299211802248_Rkind,  &
    -0.000092884721449_Rkind, 0.002349999254034_Rkind, 0.000180474434574_Rkind, &
    0.000150043218181_Rkind, -0.000410454180684_Rkind,  &
    -0.000592838266489_Rkind, 0.000117200971426_Rkind, 0.000505001599993_Rkind, &
    -0.000055649620750_Rkind, 0.001096192887978_Rkind,  &
    -0.000061590899666_Rkind, -0.000042384306295_Rkind, -0.000040709579710_Rkind,&
     0.000014266067930_Rkind, -0.000012876203424_Rkind,  &
    -0.000026527166095_Rkind, 0.000028834676520_Rkind, 0.000018902047923_Rkind, &
    -0.000030824701505_Rkind, -0.000017620441943_Rkind,  &
    0.000094324341136_Rkind, 0.000008289431852_Rkind, 0.000010494880368_Rkind, &
    0.000018844091338_Rkind, 0.000011571086755_Rkind,  &
    0.000026428202493_Rkind, -0.000016979638953_Rkind, -0.000089727181119_Rkind,&
     -0.000008031361023_Rkind, -0.000019228463780_Rkind,  &
    0.000100486328931_Rkind, 0.000164750339364_Rkind, -0.000058410395406_Rkind, &
    0.000084357995657_Rkind, -0.000102481457012_Rkind,  &
    -0.000099766245708_Rkind, 0.000272613557323_Rkind, -0.000233886530208_Rkind,&
     -0.000059285758535_Rkind, -0.000080925434931_Rkind,  &
    0.000119621843473_Rkind, -0.000259919971815_Rkind, 0.000809199126528_Rkind, &
    -0.000224145996705_Rkind, -0.000026293152716_Rkind,  &
    -0.000216919466741_Rkind, 0.000048969486509_Rkind, -0.000159571973222_Rkind,&
     -0.000041733661621_Rkind, -0.000246525986453_Rkind,  &
    -0.000190819867146_Rkind, 0.000269791181014_Rkind, 0.000246301267998_Rkind, &
    0.000149939698244_Rkind, 0.000021287016047_Rkind,  &
    0.000019804749063_Rkind, 0.000259216291399_Rkind, -0.000310253625099_Rkind, &
    -0.001189421111542_Rkind, -0.000359752375511_Rkind,  &
    0.000575103004898_Rkind, -0.000030374717836_Rkind, -0.000149686001498_Rkind,&
     -0.000202738693407_Rkind, 0.000094790727612_Rkind,  &
    -0.000012878390465_Rkind, -0.000107403939366_Rkind, 0.000241795416940_Rkind,&
     0.000355186927588_Rkind, -0.000300437638684_Rkind,  &
    -0.000441658151535_Rkind, -0.000087751736407_Rkind, -0.000201519964853_Rkind,&
     0.000046681294945_Rkind, 0.000294721807236_Rkind,  &
    -0.000073190964713_Rkind, -0.000232284522733_Rkind, -0.000273779796894_Rkind,&
     0.000336524360646_Rkind, 0.000713744631996_Rkind,  &
    -0.000207363191435_Rkind, -0.000335896679902_Rkind, 0.000980657120397_Rkind, &
    0.003659988924387_Rkind, 0.000754547160945_Rkind,  &
    -0.002398055924338_Rkind, 0.001205394711671_Rkind, 0.001432180466768_Rkind, &
    -0.015978888940835_Rkind, -0.001280514099716_Rkind,  &
    -0.000703525136531_Rkind, 0.002546269500536_Rkind, 0.003924823177153_Rkind, &
    -0.000973165958481_Rkind, -0.003711558255792_Rkind,  &
    -0.000287805837083_Rkind, -0.014046294921420_Rkind, -0.000016712273200_Rkind, &
    -0.000033863959372_Rkind, 0.000132216647632_Rkind,  &
    0.000235519703016_Rkind, -0.000089973952236_Rkind, -0.000021428262440_Rkind, &
    0.000101478698749_Rkind, -0.000187191566257_Rkind,  &
    0.005071762385699_Rkind ]


    integer, parameter :: idx3(3, 151) = RESHAPE([ &
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
      4 , 4 , 8 , 5 , 2 , 8 , 5 , 3 , 8 , 5, &
      4 , 8 , 5 , 5 , 8 , 6 , 1 , 8 , 6 , 2, &
      8 , 6 , 3 , 8 , 6 , 4 , 8 , 6 , 5 , 8, &
      6 , 6 , 8 , 7 , 1 , 8 , 7 , 2 , 8 , 7, &
      3 , 8 , 7 , 4 , 8 , 7 , 5 , 8 , 7 , 6, &
      8 , 7 , 7 , 8 , 8 , 1 , 8 , 8 , 2 , 8, &
      8 , 3 , 8 , 8 , 4 , 8 , 8 , 5 , 8 , 8, &
      6 , 8 , 8 , 7 , 8 , 8 , 8 , 9 , 1 , 1, &
      9 , 2 , 1 , 9 , 2 , 2 , 9 , 3 , 1 , 9, &
      3 , 3 , 9 , 4 , 1 , 9 , 4 , 2 , 9 , 4, &
      3 , 9 , 4 , 4 , 9 , 5 , 1 , 9 , 5 , 2, &
      9 , 5 , 4 , 9 , 5 , 5 , 9 , 6 , 1 , 9, &
      6 , 2 , 9 , 6 , 3 , 9 , 6 , 4 , 9 , 6, &
      5 , 9 , 6 , 6 , 9 , 7 , 1 , 9 , 7 , 2, &
      9 , 7 , 3 , 9 , 7 , 4 , 9 , 7 , 5 , 9, &
      7 , 6 , 9 , 7 , 7 , 9 , 8 , 1 , 9 , 8, &
      2 , 9 , 8 , 3 , 9 , 8 , 4 , 9 , 8 , 5, &
      9 , 8 , 6 , 9 , 8 , 7 , 9 , 8 , 8 , 9, &
      9 , 2 , 9 , 9 , 3 , 9 , 9 , 4 , 9 , 9, &
      5 , 9 , 9 , 6 , 9 , 9 , 7 , 9 , 9 , 8, &
      9 , 9 , 9], shape=[3, 151])

    integer, parameter :: idx4(4, 301) = reshape([ &
      1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 3 , 2, &
      2 , 2 , 3 , 3 , 2 , 1 , 3 , 3 , 2 , 2, &
      3 , 3 , 3 , 2 , 3 , 3 , 3 , 3 , 4 , 2, &
      2 , 1 , 4 , 3 , 2 , 2 , 4 , 3 , 3 , 1, &
      4 , 3 , 3 , 3 , 4 , 4 , 2 , 1 , 4 , 4, &
      2 , 2 , 4 , 4 , 3 , 1 , 4 , 4 , 3 , 2, &
      4 , 4 , 3 , 3 , 4 , 4 , 4 , 1 , 4 , 4, &
      4 , 2 , 4 , 4 , 4 , 3 , 4 , 4 , 4 , 4, &
      5 , 2 , 2 , 1 , 5 , 2 , 2 , 2 , 5 , 3, &
      1 , 1 , 5 , 3 , 2 , 2 , 5 , 3 , 3 , 1, &
      5 , 3 , 3 , 2 , 5 , 4 , 2 , 2 , 5 , 4, &
      3 , 3 , 5 , 4 , 4 , 2 , 5 , 4 , 4 , 3, &
      5 , 4 , 4 , 4 , 5 , 5 , 1 , 1 , 5 , 5, &
      2 , 2 , 5 , 5 , 3 , 1 , 5 , 5 , 3 , 2, &
      5 , 5 , 3 , 3 , 5 , 5 , 4 , 1 , 5 , 5, &
      4 , 2 , 5 , 5 , 4 , 3 , 5 , 5 , 4 , 4, &
      5 , 5 , 5 , 2 , 5 , 5 , 5 , 3 , 5 , 5, &
      5 , 4 , 5 , 5 , 5 , 5 , 6 , 2 , 2 , 1, &
      6 , 2 , 2 , 2 , 6 , 3 , 2 , 2 , 6 , 3, &
      3 , 1 , 6 , 3 , 3 , 2 , 6 , 3 , 3 , 3, &
      6 , 4 , 3 , 3 , 6 , 4 , 4 , 1 , 6 , 4, &
      4 , 2 , 6 , 4 , 4 , 3 , 6 , 4 , 4 , 4, &
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
      2 , 2 , 7 , 3 , 2 , 2 , 7 , 4 , 4 , 1, &
      7 , 4 , 4 , 2 , 7 , 4 , 4 , 3 , 7 , 4, &
      4 , 4 , 7 , 5 , 3 , 3 , 7 , 5 , 4 , 4, &
      7 , 5 , 5 , 1 , 7 , 5 , 5 , 2 , 7 , 5, &
      5 , 3 , 7 , 5 , 5 , 4 , 7 , 5 , 5 , 5, &
      7 , 6 , 2 , 2 , 7 , 6 , 3 , 3 , 7 , 6, &
      4 , 4 , 7 , 6 , 5 , 5 , 7 , 6 , 6 , 1, &
      7 , 6 , 6 , 2 , 7 , 6 , 6 , 3 , 7 , 6, &
      6 , 4 , 7 , 6 , 6 , 5 , 7 , 6 , 6 , 6, &
      7 , 7 , 1 , 1 , 7 , 7 , 2 , 1 , 7 , 7, &
      2 , 2 , 7 , 7 , 3 , 1 , 7 , 7 , 3 , 2, &
      7 , 7 , 3 , 3 , 7 , 7 , 4 , 1 , 7 , 7, &
      4 , 2 , 7 , 7 , 4 , 3 , 7 , 7 , 4 , 4, &
      7 , 7 , 5 , 1 , 7 , 7 , 5 , 2 , 7 , 7, &
      5 , 3 , 7 , 7 , 5 , 4 , 7 , 7 , 5 , 5, &
      7 , 7 , 6 , 1 , 7 , 7 , 6 , 2 , 7 , 7, &
      6 , 3 , 7 , 7 , 6 , 4 , 7 , 7 , 6 , 6, &
      7 , 7 , 7 , 1 , 7 , 7 , 7 , 2 , 7 , 7, &
      7 , 3 , 7 , 7 , 7 , 4 , 7 , 7 , 7 , 5, &
      7 , 7 , 7 , 6 , 7 , 7 , 7 , 7 , 8 , 2, &
      2 , 2 , 8 , 3 , 2 , 2 , 8 , 3 , 3 , 1, &
      8 , 3 , 3 , 3 , 8 , 4 , 1 , 1 , 8 , 4, &
      2 , 2 , 8 , 4 , 3 , 3 , 8 , 4 , 4 , 1, &
      8 , 4 , 4 , 2 , 8 , 4 , 4 , 3 , 8 , 4, &
      4 , 4 , 8 , 5 , 2 , 2 , 8 , 5 , 3 , 3, &
      8 , 5 , 4 , 4 , 8 , 5 , 5 , 3 , 8 , 5, &
      5 , 4 , 8 , 5 , 5 , 5 , 8 , 6 , 2 , 2, &
      8 , 6 , 3 , 3 , 8 , 6 , 4 , 4 , 8 , 6, &
      5 , 5 , 8 , 6 , 6 , 1 , 8 , 6 , 6 , 2, &
      8 , 6 , 6 , 3 , 8 , 6 , 6 , 4 , 8 , 6, &
      6 , 5 , 8 , 6 , 6 , 6 , 8 , 7 , 1 , 1, &
      8 , 7 , 3 , 3 , 8 , 7 , 4 , 4 , 8 , 7, &
      5 , 5 , 8 , 7 , 6 , 6 , 8 , 7 , 7 , 1, &
      8 , 7 , 7 , 2 , 8 , 7 , 7 , 3 , 8 , 7, &
      7 , 4 , 8 , 7 , 7 , 5 , 8 , 7 , 7 , 6, &
      8 , 7 , 7 , 7 , 8 , 8 , 1 , 1 , 8 , 8, &
      2 , 1 , 8 , 8 , 2 , 2 , 8 , 8 , 3 , 1, &
      8 , 8 , 3 , 2 , 8 , 8 , 4 , 2 , 8 , 8, &
      4 , 3 , 8 , 8 , 4 , 4 , 8 , 8 , 5 , 1, &
      8 , 8 , 5 , 2 , 8 , 8 , 5 , 3 , 8 , 8, &
      5 , 4 , 8 , 8 , 5 , 5 , 8 , 8 , 6 , 1, &
      8 , 8 , 6 , 2 , 8 , 8 , 6 , 3 , 8 , 8, &
      6 , 4 , 8 , 8 , 6 , 5 , 8 , 8 , 6 , 6, &
      8 , 8 , 7 , 1 , 8 , 8 , 7 , 2 , 8 , 8, &
      7 , 3 , 8 , 8 , 7 , 4 , 8 , 8 , 7 , 5, &
      8 , 8 , 7 , 6 , 8 , 8 , 7 , 7 , 8 , 8, &
      8 , 1 , 8 , 8 , 8 , 2 , 8 , 8 , 8 , 3, &
      8 , 8 , 8 , 4 , 8 , 8 , 8 , 5 , 8 , 8, &
      8 , 6 , 8 , 8 , 8 , 7 , 8 , 8 , 8 , 8, &
      9 , 4 , 1 , 1 , 9 , 4 , 2 , 2 , 9 , 4, &
      3 , 3 , 9 , 4 , 4 , 1 , 9 , 4 , 4 , 2, &
      9 , 4 , 4 , 3 , 9 , 4 , 4 , 4 , 9 , 5, &
      3 , 3 , 9 , 5 , 4 , 4 , 9 , 5 , 5 , 1, &
      9 , 5 , 5 , 3 , 9 , 5 , 5 , 4 , 9 , 5, &
      5 , 5 , 9 , 6 , 1 , 1 , 9 , 6 , 2 , 2, &
      9 , 6 , 3 , 3 , 9 , 6 , 4 , 4 , 9 , 6, &
      5 , 5 , 9 , 6 , 6 , 1 , 9 , 6 , 6 , 2, &
      9 , 6 , 6 , 3 , 9 , 6 , 6 , 4 , 9 , 6, &
      6 , 5 , 9 , 6 , 6 , 6 , 9 , 7 , 1 , 1, &
      9 , 7 , 2 , 2 , 9 , 7 , 4 , 4 , 9 , 7, &
      5 , 5 , 9 , 7 , 6 , 6 , 9 , 7 , 7 , 1, &
      9 , 7 , 7 , 2 , 9 , 7 , 7 , 3 , 9 , 7, &
      7 , 4 , 9 , 7 , 7 , 5 , 9 , 7 , 7 , 6, &
      9 , 7 , 7 , 7 , 9 , 8 , 1 , 1 , 9 , 8, &
      2 , 2 , 9 , 8 , 3 , 3 , 9 , 8 , 4 , 4, &
      9 , 8 , 5 , 5 , 9 , 8 , 6 , 6 , 9 , 8, &
      7 , 7 , 9 , 8 , 8 , 1 , 9 , 8 , 8 , 2, &
      9 , 8 , 8 , 3 , 9 , 8 , 8 , 4 , 9 , 8, &
      8 , 5 , 9 , 8 , 8 , 6 , 9 , 8 , 8 , 7, &
      9 , 8 , 8 , 8 , 9 , 9 , 1 , 1 , 9 , 9, &
      2 , 1 , 9 , 9 , 2 , 2 , 9 , 9 , 3 , 1, &
      9 , 9 , 3 , 2 , 9 , 9 , 3 , 3 , 9 , 9, &
      4 , 1 , 9 , 9 , 4 , 2 , 9 , 9 , 4 , 3, &
      9 , 9 , 4 , 4 , 9 , 9 , 5 , 1 , 9 , 9, &
      5 , 2 , 9 , 9 , 5 , 3 , 9 , 9 , 5 , 4, &
      9 , 9 , 5 , 5 , 9 , 9 , 6 , 1 , 9 , 9, &
      6 , 2 , 9 , 9 , 6 , 3 , 9 , 9 , 6 , 4, &
      9 , 9 , 6 , 5 , 9 , 9 , 6 , 6 , 9 , 9, &
      7 , 1 , 9 , 9 , 7 , 2 , 9 , 9 , 7 , 3, &
      9 , 9 , 7 , 4 , 9 , 9 , 7 , 5 , 9 , 9, &
      7 , 6 , 9 , 9 , 7 , 7 , 9 , 9 , 8 , 1, &
      9 , 9 , 8 , 2 , 9 , 9 , 8 , 3 , 9 , 9, &
      8 , 4 , 9 , 9 , 8 , 5 , 9 , 9 , 8 , 6, &
      9 , 9 , 8 , 7 , 9 , 9 , 8 , 8 , 9 , 9, &
      9 , 1 , 9 , 9 , 9 , 2 , 9 , 9 , 9 , 3, &
      9 , 9 , 9 , 4 , 9 , 9 , 9 , 5 , 9 , 9, &
      9 , 6 , 9 , 9 , 9 , 7 , 9 , 9 , 9 , 8, &
      9 , 9 , 9 , 9], shape=[4, 301])

    real(kind=Rkind), parameter :: Q_extremes(9) = [ &
      15.00000_Rkind, 15.00000_Rkind, 15.00000_Rkind,    &
      15.00000_Rkind, 15.00000_Rkind, 15.00000_Rkind,    &
      15.00000_Rkind, 15.00000_Rkind, 15.00000_Rkind]  ! for 1D part

    real(kind=Rkind), parameter :: Qc = 4.25_Rkind
    integer,          parameter :: n1 = 24
    integer,          parameter :: n3 = 8


    integer          :: i,mode
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
    DO i = 1, size(QModel%quadratic)
      E2 = E2 + HALF * QModel%quadratic(i) * dnQ(i)**2
    END DO

    E3 = ZERO
    IF (QModel%norder >= 3) THEN
      IF (QModel%bQFF) THEN
        DO i = 1, size(cubic)
            IF (idx3(1,i) == idx3(2,i) .AND. idx3(2,i) == idx3(3,i)) THEN
                mode = idx3(1,i)
                IF (Q_extremes(mode) == 15.0_Rkind) THEN
                    E3 = E3 + (ONE/SIX) * cubic(i) * dnQ(mode)**3
                ELSE
                    E3 = E3 + (ONE/SIX) * cubic(i) * dnQ(mode)**3 * &
                          EXP(-(dnQ(mode)/Q_extremes(mode))**n1)
                END IF
            ELSE
                dnDR = dnQ(idx3(1,i))**n3 + dnQ(idx3(2,i))**n3 + dnQ(idx3(3,i))**n3
                E3 = E3 + (ONE/SIX) * cubic(i) * dnQ(idx3(1,i)) * &
                         dnQ(idx3(2,i)) * dnQ(idx3(3,i)) * EXP(-dnDR/Qc**n3)
            END IF
        END DO
      ELSE
        DO i = 1, size(cubic)
          E3 = E3 + (ONE/SIX) * cubic(i) * dnQ(idx3(1,i))*dnQ(idx3(2,i))*dnQ(idx3(3,i))
        END DO
      END IF
    END IF

    E4 = ZERO
    IF (QModel%norder >= 4) THEN
      IF (QModel%bQFF) THEN
        DO i = 1,size(quartic)
          IF (idx4(1,i) == idx4(2,i) .AND. idx4(2,i) == idx4(3,i) .AND. idx4(3,i) == idx4(4,i)) THEN
                mode = idx4(1,i)
                IF (Q_extremes(mode) == 15.0_Rkind) THEN
                    E4 = E4 + (ONE/24.0_Rkind) * quartic(i) * dnQ(mode)**4
                ELSE
                    E4 = E4 + (ONE/24.0_Rkind) * quartic(i) * dnQ(mode)**4 * &
                          EXP(-(dnQ(mode)/Q_extremes(mode))**n1)
                END IF
          ELSE
              dnDR = dnQ(idx4(1,i))**n3 + dnQ(idx4(2,i))**n3 + dnQ(idx4(3,i))**n3 + dnQ(idx4(4,i))**n3
              E4 = E4 + (ONE/24.0_Rkind) * quartic(i) * dnQ(idx4(1,i)) * &
                         dnQ(idx4(2,i)) * dnQ(idx4(3,i)) * dnQ(idx4(4,i)) * EXP(-dnDR/Qc**n3)
          END IF
        END DO
      ELSE
        DO i = 1,size(quartic)
          E4 = E4 + (ONE/24.0_Rkind) * quartic(i) * dnQ(idx4(1,i))*dnQ(idx4(2,i))*dnQ(idx4(3,i))*dnQ(idx4(4,i))
        END DO
      END IF
    END IF

    Mat_OF_PotDia(1,1) =  E2 + E3 + E4


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
    real (kind=Rkind), allocatable :: wi(:)
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
      Q0(:) = [0.1_Rkind,-0.1_Rkind,0.2_Rkind,-0.2_Rkind,0.3_Rkind,-0.3_Rkind,0.4_Rkind, &
              -0.4_Rkind,0.5_Rkind]
    END IF

    IF (present(dnMatV)) THEN
      err = 0

      IF (nderiv >= 0) THEN ! no derivative
        V  = [3.2310845043165954E-003_Rkind]
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
      wi = [0.001049576065177_Rkind, 0.001461717593495_Rkind,  &
            0.001969821921105_Rkind, 0.003125330502741_Rkind,  &
            0.003688471806263_Rkind,  0.004989751221640_Rkind, &
            0.005689101160362_Rkind, 0.006101642415972_Rkind,  &
            0.014536840636557_Rkind ]
      DO i=1,QModel%ndim
        d0GGdef(i,i) = wi(i)
      END DO
    END IF


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_CHFClBr
END MODULE QML_CHFClBr_m
