
MODULE QML_dmabn_m
    USE QDUtil_NumParameters_m
    USE QML_Empty_m
    IMPLICIT NONE
  
    PRIVATE
  
  !> @brief Derived type in which the XXX parameters are set-up.
    TYPE, EXTENDS (QML_Empty_t) ::  QML_dmabn_t
  
     PRIVATE
     integer :: i
     real(kind=Rkind), allocatable  :: omega_(:)
     real(kind=Rkind), allocatable  :: kappa1_(:),kappa2_(:),kappa3_(:)
     real(kind=Rkind), allocatable  :: lambda1_2_(:),lambda1_3_(:),lambda2_3_(:)

     real(kind=Rkind)   :: E1,E2,E3

     

     CONTAINS
      PROCEDURE :: EvalPot_QModel  => EvalPot_QML_dmabn
      PROCEDURE :: Write_QModel    => Write_QML_dmabn
    END TYPE QML_dmabn_t
  
    PUBLIC :: QML_dmabn_t,Init_QML_dmabn
  
    CONTAINS
  !> @brief Function which makes the initialization of the XXX parameters.
  !!
  !! @param QModel             TYPE(QML_XXX_t):   result derived type in which the parameters are set-up.
  !! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
  !! @param nio_param_file     integer:             file unit to read the parameters.
  !! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
    FUNCTION Init_QML_dmabn(QModel_in,read_param,nio_param_file) RESULT(QModel)
      USE QDUtil_m,         ONLY : Identity_Mat
      IMPLICIT NONE
  
      TYPE (QML_dmabn_t)                           :: QModel ! RESULT
  
      TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
      integer,                     intent(in)      :: nio_param_file
      logical,                     intent(in)      :: read_param
  
  
      !----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Init_QML_dmabn'
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        flush(out_unit)
      END IF
  
      QModel%QML_Empty_t = QModel_in
      QModel%nsurf    = 3
      QModel%ndim     = 57
      QModel%pot_name = 'dmabn'
      allocate(QModel%omega_(QModel%ndim))
      allocate(QModel%kappa1_(QModel%ndim),QModel%kappa2_(QModel%ndim),QModel%kappa3_(QModel%ndim))
      allocate(QModel%lambda1_2_(QModel%ndim),QModel%lambda2_3_(QModel%ndim),QModel%lambda1_3_(QModel%ndim))
  
      IF (debug) write(out_unit,*) 'Paramters of dmabn'
      
      QModel%omega_      =0d0
      QModel%kappa1_    =0d0
      QModel%kappa2_    =0d0
      QModel%kappa3_    =0d0
      QModel%lambda1_2_ =0d0
      QModel%lambda2_3_ =0d0
      QModel%lambda1_3_ =0d0  

      !> Parameters are from Phys. Chem. Chem. Phys., 2024,26, 1829-1844, doi.org/10.1039/D3CP03964A

      QModel%omega_(1)                 =   0.00622_Rkind! , ev
      QModel%omega_(2)                 =   0.00909_Rkind! , ev
      QModel%omega_(3)                 =   0.01033_Rkind! , ev
      QModel%omega_(4)                 =   0.01673_Rkind! , ev
      QModel%omega_(5)                 =   0.02188_Rkind! , ev
      QModel%omega_(6)                 =   0.02317_Rkind! , ev
      QModel%omega_(7)                 =   0.03186_Rkind! , ev
      QModel%omega_(8)                 =   0.03575_Rkind! , ev
      QModel%omega_(9)                 =   0.04213_Rkind! , ev
      QModel%omega_(10)                =   0.05324_Rkind! , ev
      QModel%omega_(11)                =   0.05919_Rkind! , ev
      QModel%omega_(12)                =   0.06109_Rkind! , ev
      QModel%omega_(13)                =   0.06234_Rkind! , ev
      QModel%omega_(14)                =   0.07081_Rkind! , ev
      QModel%omega_(15)                =   0.07103_Rkind! , ev
      QModel%omega_(16)                =   0.08237_Rkind! , ev
      QModel%omega_(17)                =   0.08386_Rkind! , ev
      QModel%omega_(18)                =   0.09314_Rkind! , ev
      QModel%omega_(19)                =   0.10163_Rkind! , ev
      QModel%omega_(20)                =   0.10377_Rkind! , ev
      QModel%omega_(21)                =   0.10575_Rkind! , ev
      QModel%omega_(22)                =   0.12272_Rkind! , ev
      QModel%omega_(23)                =   0.12294_Rkind! , ev
      QModel%omega_(24)                =   0.12408_Rkind! , ev
      QModel%omega_(25)                =   0.12707_Rkind! , ev
      QModel%omega_(26)                =   0.13479_Rkind! , ev
      QModel%omega_(27)                =   0.14102_Rkind! , ev
      QModel%omega_(28)                =   0.14110_Rkind! , ev
      QModel%omega_(29)                =   0.14298_Rkind! , ev
      QModel%omega_(30)                =   0.14786_Rkind! , ev
      QModel%omega_(31)                =   0.14935_Rkind! , ev
      QModel%omega_(32)                =   0.15535_Rkind! , ev
      QModel%omega_(33)                =   0.16070_Rkind! , ev
      QModel%omega_(34)                =   0.16529_Rkind! , ev
      QModel%omega_(35)                =   0.17258_Rkind! , ev
      QModel%omega_(36)                =   0.17568_Rkind! , ev
      QModel%omega_(37)                =   0.17782_Rkind! , ev
      QModel%omega_(38)                =   0.18109_Rkind! , ev
      QModel%omega_(39)                =   0.18219_Rkind! , ev
      QModel%omega_(40)                =   0.18324_Rkind! , ev
      QModel%omega_(41)                =   0.18354_Rkind! , ev
      QModel%omega_(42)                =   0.18661_Rkind! , ev
      QModel%omega_(43)                =   0.18749_Rkind! , ev
      QModel%omega_(44)                =   0.19708_Rkind! , ev
      QModel%omega_(45)                =   0.20350_Rkind! , ev
      QModel%omega_(46)                =   0.21190_Rkind! , ev
      QModel%omega_(47)                =   0.29621_Rkind! , ev
      QModel%omega_(48)                =   0.37412_Rkind! , ev
      QModel%omega_(49)                =   0.37509_Rkind! , ev
      QModel%omega_(50)                =   0.38359_Rkind! , ev
      QModel%omega_(51)                =   0.38373_Rkind! , ev
      QModel%omega_(52)                =   0.39291_Rkind! , ev
      QModel%omega_(53)                =   0.39407_Rkind! , ev
      QModel%omega_(54)                =   0.39931_Rkind! , ev
      QModel%omega_(55)                =   0.39943_Rkind! , ev
      QModel%omega_(56)                =   0.40299_Rkind! , ev
      QModel%omega_(57)                =   0.40309_Rkind! , ev
      QModel%E1                      =   0.00000_Rkind! , ev
      QModel%E2                      =   4.93397_Rkind! , ev
      QModel%E3                      =   5.31142_Rkind! , ev
      QModel%kappa1_(3)                =  -0.00036_Rkind! , ev
      QModel%kappa1_(4)                =  -0.00082_Rkind! , ev
      QModel%kappa1_(5)                =  -0.00045_Rkind! , ev
      QModel%kappa1_(7)                =   0.00032_Rkind! , ev
      QModel%kappa1_(8)                =   0.00037_Rkind! , ev
      QModel%kappa1_(9)                =  -0.00043_Rkind! , ev
      QModel%kappa1_(11)               =  -0.00049_Rkind! , ev
      QModel%kappa1_(12)               =  -0.00071_Rkind! , ev
      QModel%kappa1_(13)               =  -0.00038_Rkind! , ev
      QModel%kappa1_(15)               =   0.00018_Rkind! , ev
      QModel%kappa1_(17)               =  -0.00019_Rkind! , ev
      QModel%kappa1_(19)               =  -0.00044_Rkind! , ev
      QModel%kappa1_(22)               =   0.00033_Rkind! , ev
      QModel%kappa1_(23)               =  -0.00038_Rkind! , ev
      QModel%kappa1_(25)               =  -0.00026_Rkind! , ev
      QModel%kappa1_(26)               =  -0.00015_Rkind! , ev
      QModel%kappa1_(29)               =   0.00014_Rkind! , ev
      QModel%kappa1_(31)               =  -0.00033_Rkind! , ev
      QModel%kappa1_(32)               =   0.00046_Rkind! , ev
      QModel%kappa1_(36)               =  -0.00066_Rkind! , ev
      QModel%kappa1_(37)               =  -0.00025_Rkind! , ev
      QModel%kappa1_(41)               =   0.00036_Rkind! , ev
      QModel%kappa1_(43)               =   0.00014_Rkind! , ev
      QModel%kappa1_(46)               =   0.00012_Rkind! , ev
      QModel%kappa1_(48)               =  -0.00014_Rkind! , ev
      QModel%kappa1_(49)               =   0.00011_Rkind! , ev
      QModel%kappa1_(53)               =  -0.00038_Rkind! , ev
      QModel%kappa1_(54)               =  -0.00018_Rkind! , ev
      QModel%kappa1_(56)               =  -0.00028_Rkind! , ev
      QModel%kappa1_(57)               =  -0.00027_Rkind! , ev
      QModel%kappa2_(1)                =   0.01981_Rkind! , ev
      QModel%kappa2_(2)                =   0.00114_Rkind! , ev
      QModel%kappa2_(3)                =   0.00679_Rkind! , ev
      QModel%kappa2_(4)                =  -0.00129_Rkind! , ev
      QModel%kappa2_(5)                =  -0.01080_Rkind! , ev
      QModel%kappa2_(6)                =  -0.00051_Rkind! , ev
      QModel%kappa2_(7)                =   0.00074_Rkind! , ev
      QModel%kappa2_(8)                =   0.01364_Rkind! , ev
      QModel%kappa2_(9)                =   0.01326_Rkind! , ev
      QModel%kappa2_(10)               =   0.00025_Rkind! , ev
      QModel%kappa2_(11)               =  -0.00067_Rkind! , ev
      QModel%kappa2_(12)               =  -0.02516_Rkind! , ev
      QModel%kappa2_(13)               =  -0.01572_Rkind! , ev
      QModel%kappa2_(14)               =   0.00029_Rkind! , ev
      QModel%kappa2_(15)               =   0.00120_Rkind! , ev
      QModel%kappa2_(16)               =  -0.00017_Rkind! , ev
      QModel%kappa2_(17)               =   0.00630_Rkind! , ev
      QModel%kappa2_(18)               =   0.00095_Rkind! , ev
      QModel%kappa2_(19)               =  -0.09430_Rkind! , ev
      QModel%kappa2_(21)               =  -0.00292_Rkind! , ev
      QModel%kappa2_(22)               =  -0.03331_Rkind! , ev
      QModel%kappa2_(23)               =   0.04050_Rkind! , ev
      QModel%kappa2_(25)               =  -0.00984_Rkind! , ev
      QModel%kappa2_(26)               =  -0.00036_Rkind! , ev
      QModel%kappa2_(27)               =  -0.01563_Rkind! , ev
      QModel%kappa2_(28)               =   0.00299_Rkind! , ev
      QModel%kappa2_(29)               =  -0.00022_Rkind! , ev
      QModel%kappa2_(30)               =  -0.01587_Rkind! , ev
      QModel%kappa2_(31)               =   0.00569_Rkind! , ev
      QModel%kappa2_(32)               =   0.07393_Rkind! , ev
      QModel%kappa2_(33)               =   0.00135_Rkind! , ev
      QModel%kappa2_(34)               =  -0.00099_Rkind! , ev
      QModel%kappa2_(35)               =   0.00021_Rkind! , ev
      QModel%kappa2_(36)               =   0.07246_Rkind! , ev
      QModel%kappa2_(37)               =   0.00079_Rkind! , ev
      QModel%kappa2_(38)               =  -0.00031_Rkind! , ev
      QModel%kappa2_(39)               =   0.00674_Rkind! , ev
      QModel%kappa2_(40)               =   0.01043_Rkind! , ev
      QModel%kappa2_(41)               =  -0.00021_Rkind! , ev
      QModel%kappa2_(42)               =  -0.00014_Rkind! , ev
      QModel%kappa2_(43)               =   0.02299_Rkind! , ev
      QModel%kappa2_(44)               =   0.00436_Rkind! , ev
      QModel%kappa2_(45)               =  -0.00038_Rkind! , ev
      QModel%kappa2_(46)               =   0.10678_Rkind! , ev
      QModel%kappa2_(47)               =  -0.02123_Rkind! , ev
      QModel%kappa2_(48)               =  -0.00070_Rkind! , ev
      QModel%kappa2_(49)               =  -0.01697_Rkind! , ev
      QModel%kappa2_(50)               =   0.00086_Rkind! , ev
      QModel%kappa2_(51)               =   0.00237_Rkind! , ev
      QModel%kappa2_(52)               =  -0.00024_Rkind! , ev
      QModel%kappa2_(53)               =  -0.00153_Rkind! , ev
      QModel%kappa2_(54)               =  -0.02192_Rkind! , ev
      QModel%kappa2_(55)               =   0.00756_Rkind! , ev
      QModel%kappa2_(56)               =  -0.00137_Rkind! , ev
      QModel%kappa2_(57)               =   0.03034_Rkind! , ev
      QModel%kappa3_(1)                =   0.00797_Rkind! , ev
      QModel%kappa3_(2)                =  -0.00096_Rkind! , ev
      QModel%kappa3_(3)                =  -0.00024_Rkind! , ev
      QModel%kappa3_(4)                =  -0.00094_Rkind! , ev
      QModel%kappa3_(5)                =  -0.00628_Rkind! , ev
      QModel%kappa3_(6)                =  -0.00062_Rkind! , ev
      QModel%kappa3_(7)                =   0.00110_Rkind! , ev
      QModel%kappa3_(8)                =   0.00920_Rkind! , ev
      QModel%kappa3_(9)                =  -0.00695_Rkind! , ev
      QModel%kappa3_(10)               =  -0.00019_Rkind! , ev
      QModel%kappa3_(11)               =  -0.00010_Rkind! , ev
      QModel%kappa3_(12)               =  -0.01937_Rkind! , ev
      QModel%kappa3_(13)               =  -0.01093_Rkind! , ev
      QModel%kappa3_(14)               =  -0.00014_Rkind! , ev
      QModel%kappa3_(15)               =  -0.00166_Rkind! , ev
      QModel%kappa3_(17)               =  -0.00475_Rkind! , ev
      QModel%kappa3_(18)               =  -0.00549_Rkind! , ev
      QModel%kappa3_(19)               =  -0.07216_Rkind! , ev
      QModel%kappa3_(20)               =   0.00024_Rkind! , ev
      QModel%kappa3_(21)               =   0.00189_Rkind! , ev
      QModel%kappa3_(22)               =  -0.03457_Rkind! , ev
      QModel%kappa3_(23)               =   0.04590_Rkind! , ev
      QModel%kappa3_(25)               =  -0.01089_Rkind! , ev
      QModel%kappa3_(26)               =  -0.00013_Rkind! , ev
      QModel%kappa3_(27)               =  -0.01695_Rkind! , ev
      QModel%kappa3_(28)               =   0.00394_Rkind! , ev
      QModel%kappa3_(29)               =   0.00141_Rkind! , ev
      QModel%kappa3_(30)               =  -0.02115_Rkind! , ev
      QModel%kappa3_(31)               =  -0.07638_Rkind! , ev
      QModel%kappa3_(32)               =   0.04662_Rkind! , ev
      QModel%kappa3_(33)               =   0.00065_Rkind! , ev
      QModel%kappa3_(34)               =   0.00011_Rkind! , ev
      QModel%kappa3_(36)               =  -0.00036_Rkind! , ev
      QModel%kappa3_(37)               =  -0.00018_Rkind! , ev
      QModel%kappa3_(38)               =  -0.00029_Rkind! , ev
      QModel%kappa3_(39)               =   0.00487_Rkind! , ev
      QModel%kappa3_(40)               =  -0.00584_Rkind! , ev
      QModel%kappa3_(41)               =   0.00064_Rkind! , ev
      QModel%kappa3_(43)               =  -0.01030_Rkind! , ev
      QModel%kappa3_(44)               =   0.06696_Rkind! , ev
      QModel%kappa3_(45)               =   0.00062_Rkind! , ev
      QModel%kappa3_(46)               =  -0.12549_Rkind! , ev
      QModel%kappa3_(47)               =  -0.10953_Rkind! , ev
      QModel%kappa3_(48)               =   0.00058_Rkind! , ev
      QModel%kappa3_(49)               =   0.02129_Rkind! , ev
      QModel%kappa3_(50)               =   0.00104_Rkind! , ev
      QModel%kappa3_(51)               =   0.00253_Rkind! , ev
      QModel%kappa3_(52)               =  -0.00031_Rkind! , ev
      QModel%kappa3_(53)               =   0.00905_Rkind! , ev
      QModel%kappa3_(54)               =  -0.01492_Rkind! , ev
      QModel%kappa3_(55)               =   0.00526_Rkind! , ev
      QModel%kappa3_(56)               =  -0.00034_Rkind! , ev
      QModel%kappa3_(57)               =   0.00257_Rkind! , ev
      QModel%lambda1_2_(1)             =   0.00023_Rkind! , ev
      QModel%lambda1_2_(2)             =   0.00716_Rkind! , ev
      QModel%lambda1_2_(3)             =  -0.00237_Rkind! , ev
      QModel%lambda1_2_(4)             =   0.00655_Rkind! , ev
      QModel%lambda1_2_(5)             =  -0.00235_Rkind! , ev
      QModel%lambda1_2_(6)             =   0.02356_Rkind! , ev
      QModel%lambda1_2_(7)             =  -0.01311_Rkind! , ev
      QModel%lambda1_2_(8)             =  -0.00098_Rkind! , ev
      QModel%lambda1_2_(9)             =   0.00015_Rkind! , ev
      QModel%lambda1_2_(10)            =  -0.02099_Rkind! , ev
      QModel%lambda1_2_(11)            =   0.04703_Rkind! , ev
      QModel%lambda1_2_(12)            =  -0.00047_Rkind! , ev
      QModel%lambda1_2_(13)            =  -0.00013_Rkind! , ev
      QModel%lambda1_2_(14)            =  -0.05253_Rkind! , ev
      QModel%lambda1_2_(15)            =  -0.00131_Rkind! , ev
      QModel%lambda1_2_(16)            =   0.01382_Rkind! , ev
      QModel%lambda1_2_(17)            =   0.00073_Rkind! , ev
      QModel%lambda1_2_(18)            =  -0.00052_Rkind! , ev
      QModel%lambda1_2_(20)            =  -0.00319_Rkind! , ev
      QModel%lambda1_2_(23)            =  -0.00066_Rkind! , ev
      QModel%lambda1_2_(24)            =   0.00046_Rkind! , ev
      QModel%lambda1_2_(25)            =  -0.00050_Rkind! , ev
      QModel%lambda1_2_(26)            =  -0.03919_Rkind! , ev
      QModel%lambda1_2_(27)            =   0.00331_Rkind! , ev
      QModel%lambda1_2_(28)            =   0.01165_Rkind! , ev
      QModel%lambda1_2_(29)            =  -0.07418_Rkind! , ev
      QModel%lambda1_2_(30)            =  -0.00155_Rkind! , ev
      QModel%lambda1_2_(31)            =   0.00043_Rkind! , ev
      QModel%lambda1_2_(32)            =   0.00221_Rkind! , ev
      QModel%lambda1_2_(33)            =  -0.23453_Rkind! , ev
      QModel%lambda1_2_(34)            =  -0.26477_Rkind! , ev
      QModel%lambda1_2_(35)            =   0.33483_Rkind! , ev
      QModel%lambda1_2_(36)            =  -0.00165_Rkind! , ev
      QModel%lambda1_2_(37)            =   0.05187_Rkind! , ev
      QModel%lambda1_2_(38)            =   0.02173_Rkind! , ev
      QModel%lambda1_2_(39)            =  -0.00409_Rkind! , ev
      QModel%lambda1_2_(40)            =   0.01359_Rkind! , ev
      QModel%lambda1_2_(41)            =   0.13492_Rkind! , ev
      QModel%lambda1_2_(42)            =   0.06252_Rkind! , ev
      QModel%lambda1_2_(43)            =  -0.00181_Rkind! , ev
      QModel%lambda1_2_(44)            =   0.00021_Rkind! , ev
      QModel%lambda1_2_(45)            =   0.03665_Rkind! , ev
      QModel%lambda1_2_(46)            =  -0.00047_Rkind! , ev
      QModel%lambda1_2_(47)            =   0.00018_Rkind! , ev
      QModel%lambda1_2_(48)            =  -0.01867_Rkind! , ev
      QModel%lambda1_2_(49)            =   0.00063_Rkind! , ev
      QModel%lambda1_2_(50)            =   0.00470_Rkind! , ev
      QModel%lambda1_2_(51)            =  -0.00202_Rkind! , ev
      QModel%lambda1_2_(52)            =  -0.00113_Rkind! , ev
      QModel%lambda1_2_(54)            =  -0.00052_Rkind! , ev
      QModel%lambda1_2_(55)            =  -0.00103_Rkind! , ev
      QModel%lambda1_2_(56)            =   0.00272_Rkind! , ev
      QModel%lambda1_2_(57)            =   0.00041_Rkind! , ev
      QModel%lambda1_3_(1)             =  -0.06569_Rkind! , ev
      QModel%lambda1_3_(3)             =  -0.02458_Rkind! , ev
      QModel%lambda1_3_(5)             =   0.02188_Rkind! , ev
      QModel%lambda1_3_(6)             =   0.00162_Rkind! , ev
      QModel%lambda1_3_(7)             =  -0.00029_Rkind! , ev
      QModel%lambda1_3_(8)             =  -0.02782_Rkind! , ev
      QModel%lambda1_3_(9)             =  -0.02332_Rkind! , ev
      QModel%lambda1_3_(10)            =   0.00032_Rkind! , ev
      QModel%lambda1_3_(11)            =  -0.00013_Rkind! , ev
      QModel%lambda1_3_(12)            =  -0.00192_Rkind! , ev
      QModel%lambda1_3_(13)            =   0.01311_Rkind! , ev
      QModel%lambda1_3_(14)            =   0.00024_Rkind! , ev
      QModel%lambda1_3_(15)            =  -0.00651_Rkind! , ev
      QModel%lambda1_3_(16)            =   0.00185_Rkind! , ev
      QModel%lambda1_3_(17)            =  -0.05622_Rkind! , ev
      QModel%lambda1_3_(18)            =  -0.00342_Rkind! , ev
      QModel%lambda1_3_(19)            =   0.01769_Rkind! , ev
      QModel%lambda1_3_(20)            =  -0.00020_Rkind! , ev
      QModel%lambda1_3_(21)            =   0.00219_Rkind! , ev
      QModel%lambda1_3_(22)            =  -0.02000_Rkind! , ev
      QModel%lambda1_3_(23)            =   0.02829_Rkind! , ev
      QModel%lambda1_3_(25)            =  -0.08859_Rkind! , ev
      QModel%lambda1_3_(26)            =  -0.00015_Rkind! , ev
      QModel%lambda1_3_(27)            =   0.03954_Rkind! , ev
      QModel%lambda1_3_(28)            =  -0.00879_Rkind! , ev
      QModel%lambda1_3_(29)            =  -0.00048_Rkind! , ev
      QModel%lambda1_3_(30)            =   0.02818_Rkind! , ev
      QModel%lambda1_3_(31)            =   0.08258_Rkind! , ev
      QModel%lambda1_3_(32)            =  -0.04090_Rkind! , ev
      QModel%lambda1_3_(33)            =   0.00041_Rkind! , ev
      QModel%lambda1_3_(34)            =   0.00041_Rkind! , ev
      QModel%lambda1_3_(35)            =  -0.00069_Rkind! , ev
      QModel%lambda1_3_(36)            =  -0.13269_Rkind! , ev
      QModel%lambda1_3_(37)            =  -0.00241_Rkind! , ev
      QModel%lambda1_3_(38)            =   0.00037_Rkind! , ev
      QModel%lambda1_3_(39)            =  -0.01557_Rkind! , ev
      QModel%lambda1_3_(40)            =   0.04191_Rkind! , ev
      QModel%lambda1_3_(41)            =  -0.00465_Rkind! , ev
      QModel%lambda1_3_(42)            =   0.00085_Rkind! , ev
      QModel%lambda1_3_(43)            =  -0.04983_Rkind! , ev
      QModel%lambda1_3_(44)            =   0.09946_Rkind! , ev
      QModel%lambda1_3_(46)            =   0.16533_Rkind! , ev
      QModel%lambda1_3_(47)            =   0.05411_Rkind! , ev
      QModel%lambda1_3_(48)            =   0.00277_Rkind! , ev
      QModel%lambda1_3_(49)            =   0.07859_Rkind! , ev
      QModel%lambda1_3_(50)            =  -0.00343_Rkind! , ev
      QModel%lambda1_3_(51)            =  -0.00810_Rkind! , ev
      QModel%lambda1_3_(52)            =  -0.00015_Rkind! , ev
      QModel%lambda1_3_(53)            =   0.01106_Rkind! , ev
      QModel%lambda1_3_(54)            =  -0.00814_Rkind! , ev
      QModel%lambda1_3_(55)            =   0.00304_Rkind! , ev
      QModel%lambda1_3_(57)            =  -0.00324_Rkind! , ev
      QModel%lambda2_3_(2)             =   0.00195_Rkind! , ev
      QModel%lambda2_3_(4)             =   0.00617_Rkind! , ev
      QModel%lambda2_3_(5)             =  -0.00123_Rkind! , ev
      QModel%lambda2_3_(6)             =   0.01004_Rkind! , ev
      QModel%lambda2_3_(7)             =  -0.00077_Rkind! , ev
      QModel%lambda2_3_(9)             =   0.00014_Rkind! , ev
      QModel%lambda2_3_(10)            =  -0.00338_Rkind! , ev
      QModel%lambda2_3_(11)            =   0.01277_Rkind! , ev
      QModel%lambda2_3_(12)            =   0.00013_Rkind! , ev
      QModel%lambda2_3_(14)            =   0.02641_Rkind! , ev
      QModel%lambda2_3_(15)            =   0.00099_Rkind! , ev
      QModel%lambda2_3_(16)            =   0.00530_Rkind! , ev
      QModel%lambda2_3_(17)            =   0.00031_Rkind! , ev
      QModel%lambda2_3_(20)            =  -0.00208_Rkind! , ev
      QModel%lambda2_3_(24)            =  -0.00278_Rkind! , ev
      QModel%lambda2_3_(25)            =  -0.00030_Rkind! , ev
      QModel%lambda2_3_(26)            =  -0.01234_Rkind! , ev
      QModel%lambda2_3_(27)            =   0.00099_Rkind! , ev
      QModel%lambda2_3_(28)            =   0.00252_Rkind! , ev
      QModel%lambda2_3_(29)            =  -0.01659_Rkind! , ev
      QModel%lambda2_3_(30)            =  -0.00026_Rkind! , ev
      QModel%lambda2_3_(31)            =  -0.00020_Rkind! , ev
      QModel%lambda2_3_(32)            =   0.00057_Rkind! , ev
      QModel%lambda2_3_(33)            =  -0.05087_Rkind! , ev
      QModel%lambda2_3_(34)            =  -0.00962_Rkind! , ev
      QModel%lambda2_3_(35)            =   0.03118_Rkind! , ev
      QModel%lambda2_3_(36)            =  -0.00021_Rkind! , ev
      QModel%lambda2_3_(37)            =   0.01499_Rkind! , ev
      QModel%lambda2_3_(38)            =   0.00439_Rkind! , ev
      QModel%lambda2_3_(39)            =  -0.00069_Rkind! , ev
      QModel%lambda2_3_(40)            =   0.00299_Rkind! , ev
      QModel%lambda2_3_(41)            =   0.02879_Rkind! , ev
      QModel%lambda2_3_(42)            =   0.02474_Rkind! , ev
      QModel%lambda2_3_(43)            =  -0.00032_Rkind! , ev
      QModel%lambda2_3_(44)            =   0.00058_Rkind! , ev
      QModel%lambda2_3_(45)            =  -0.09530_Rkind! , ev
      QModel%lambda2_3_(46)            =  -0.00035_Rkind! , ev
      QModel%lambda2_3_(47)            =  -0.00014_Rkind! , ev
      QModel%lambda2_3_(48)            =  -0.00565_Rkind! , ev
      QModel%lambda2_3_(49)            =   0.00019_Rkind! , ev
      QModel%lambda2_3_(50)            =   0.00185_Rkind! , ev
      QModel%lambda2_3_(51)            =  -0.00073_Rkind! , ev
      QModel%lambda2_3_(52)            =  -0.00199_Rkind! , ev
      QModel%lambda2_3_(54)            =   0.00635_Rkind! , ev
      QModel%lambda2_3_(55)            =   0.01709_Rkind! , ev
      QModel%lambda2_3_(56)            =  -0.00350_Rkind! , ev
  
      IF (debug) write(out_unit,*) 'init d0GGdef of dmabn'
      QModel%d0GGdef = Identity_Mat(QModel%ndim)
  
      IF (debug) THEN
        write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
  
    END FUNCTION Init_QML_dmabn
  !> @brief Subroutine wich prints the current QML_XXX parameters.
  !!
  !! @param QModel            CLASS(QML_XXX_t):   derived type in which the parameters are set-up.
  !! @param nio               integer:              file unit to print the parameters.
    SUBROUTINE Write_QML_dmabn(QModel,nio)
      IMPLICIT NONE
  
      CLASS(QML_dmabn_t),   intent(in) :: QModel
      integer,                intent(in) :: nio
      integer  :: i

      write(nio,*)'X------------------------------------X'
      write(nio,*)'--------------------------------------'
      write(nio,*) 'model name: ',QModel%pot_name
      write(nio,*) "option 1 is for eV, option 2 for au. you chose option", QModel%option
      write(nio,*) "note that that the position input vector has no units."
      write(nio,*) "Parameters are from Phys. Chem. Chem. Phys., 2024,26, 1829-1844, doi.org/10.1039/D3CP03964A"
      write(nio,*)'--------------------------------------'
      write(nio,*)'X------------------------------------X'
      write(nio,*)'--------------------------------------'

      do i=1,QModel%ndim
      write(nio,*) 'omega_ (ev)',i,QModel%omega_(i)
      end do
      write(nio,*)'--------------------------------------'
      write(nio,*) 'E1-3 (ev)',QModel%E1,QModel%E2,QModel%E3
      write(nio,*)'--------------------------------------'  
      do i=1,QModel%ndim
        if (QModel%kappa1_(i).gt.0d000001) then
        write(nio,*) 'kappa1_ (ev)',i,QModel%kappa1_(i)
        endif
      end do
      write(nio,*)'--------------------------------------' 
      do i=1,QModel%ndim
        if (QModel%kappa2_(i).gt.0d000001) then
        write(nio,*) 'kappa2_ (ev)',i,QModel%kappa2_(i)
        endif
      end do
      write(nio,*)'--------------------------------------' 
      do i=1,QModel%ndim
        if (QModel%kappa3_(i).gt.0d000001) then
        write(nio,*) 'kappa3_ (ev)',i,QModel%kappa3_(i)
        endif
      end do
      write(nio,*)'--------------------------------------' 
      do i=1,QModel%ndim
        if (QModel%lambda1_2_(i).gt.0d000001) then
        write(nio,*) 'lambda1_2_ (ev)',i,QModel%lambda1_2_(i)
        endif
      end do
      write(nio,*)'--------------------------------------' 
      do i=1,QModel%ndim
        if (QModel%lambda1_3_(i).gt.0d000001) then
        write(nio,*) 'lambda1_3_ (ev)',i,QModel%lambda1_3_(i)
        endif
      end do
      write(nio,*)'--------------------------------------' 
      do i=1,QModel%ndim
        if (QModel%lambda2_3_(i).gt.0d000001) then
        write(nio,*) 'lambda2_3_ (ev)',i,QModel%lambda2_3_(i)
        endif
      end do

      write(out_unit,*)'X------------------------------------X'
  
    END SUBROUTINE Write_QML_dmabn
  !> @brief Subroutine wich calculates the XXX potential with derivatives up to the 2d order.
  !!
  !! @param QModel             CLASS(QML_XXX_t):   derived type in which the parameters are set-up.
  !! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
  !! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
  !! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
  !!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
    SUBROUTINE EvalPot_QML_dmabn(QModel,Mat_OF_PotDia,dnQ,nderiv)
      USE ADdnSVM_m
      IMPLICIT NONE
  
      CLASS(QML_dmabn_t),     intent(in)    :: QModel
      TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
      TYPE (dnS_t),         intent(in)    :: dnQ(:)    !>  DIMENSIONLESS (no units)
      integer,              intent(in)    :: nderiv
      TYPE (dnS_t) :: dnH0,  dnW0_1,dnW0_2,dnW0_3,  dnW1
      real (kind=Rkind),parameter :: au_to_eV      = 27.2114_Rkind
      !< Conversion factor from Hartree to electronvolts

      !H0 Hamiltionian
      dnH0=0.5_Rkind*sum(QModel%omega_(:)*dnQ(:)**2 )

      !W1 diagonal
      dnW0_1=sum(QModel%kappa1_(:)*dnQ(:))
      dnW0_2=sum(QModel%kappa2_(:)*dnQ(:))
      dnW0_3=sum(QModel%kappa3_(:)*dnQ(:))

      !Diagonal elements total:
      Mat_OF_PotDia(1,1) = dnH0+dnW0_1+QModel%E1
      Mat_OF_PotDia(2,2) = dnH0+dnW0_2+QModel%E2
      Mat_OF_PotDia(3,3) = dnH0+dnW0_3+QModel%E3

      !Upper off-diagonal
      Mat_OF_PotDia(1,2) =sum(QModel%lambda1_2_(:)*dnQ(:))
      Mat_OF_PotDia(1,3) =sum(QModel%lambda1_3_(:)*dnQ(:))
      Mat_OF_PotDia(2,3) =sum(QModel%lambda2_3_(:)*dnQ(:))
                 

      !Lower off-diagonal
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2) 
      Mat_OF_PotDia(3,1) = Mat_OF_PotDia(1,3)
      Mat_OF_PotDia(3,2) = Mat_OF_PotDia(2,3) 
                    

      Select case (QModel%option)
      case(2)
        Mat_OF_PotDia=Mat_OF_PotDia/au_to_eV
      case default
      end select


    END SUBROUTINE EvalPot_QML_dmabn

  END MODULE QML_dmabn_m
  