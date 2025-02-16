
MODULE QML_fulvene_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the XXX parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_fulvene_t

   PRIVATE
   integer :: i
   real(kind=Rkind), allocatable  :: omega(:)
   real(kind=Rkind)   :: E1,E2

   REAL(kind=Rkind) :: kappa1_6
   REAL(kind=Rkind) :: kappa1_13
   REAL(kind=Rkind) :: kappa1_15
   REAL(kind=Rkind) :: kappa1_17
   REAL(kind=Rkind) :: kappa1_20
   REAL(kind=Rkind) :: kappa1_21
   REAL(kind=Rkind) :: kappa1_22
   REAL(kind=Rkind) :: kappa1_24
   REAL(kind=Rkind) :: kappa1_25
   REAL(kind=Rkind) :: kappa1_28
   REAL(kind=Rkind) :: kappa1_30
   REAL(kind=Rkind) :: kappa2_6
   REAL(kind=Rkind) :: kappa2_13
   REAL(kind=Rkind) :: kappa2_15
   REAL(kind=Rkind) :: kappa2_17
   REAL(kind=Rkind) :: kappa2_20
   REAL(kind=Rkind) :: kappa2_21
   REAL(kind=Rkind) :: kappa2_22
   REAL(kind=Rkind) :: kappa2_24
   REAL(kind=Rkind) :: kappa2_25
   REAL(kind=Rkind) :: kappa2_28
   REAL(kind=Rkind) :: kappa2_30
   REAL(kind=Rkind) :: lambda1_2_2
   REAL(kind=Rkind) :: lambda1_2_9
   REAL(kind=Rkind) :: lambda1_2_14
   REAL(kind=Rkind) :: lambda1_2_16
   REAL(kind=Rkind) :: lambda1_2_18
   REAL(kind=Rkind) :: lambda1_2_19
   REAL(kind=Rkind) :: lambda1_2_23
   REAL(kind=Rkind) :: lambda1_2_26
   REAL(kind=Rkind) :: lambda1_2_27
   REAL(kind=Rkind) :: lambda1_2_29

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_fulvene
    PROCEDURE :: Write_QModel    => Write_QML_fulvene
  END TYPE QML_fulvene_t

  PUBLIC :: QML_fulvene_t,Init_QML_fulvene

  CONTAINS
!> @brief Function which makes the initialization of the XXX parameters.
!!
!! @param QModel             TYPE(QML_XXX_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_fulvene(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_fulvene_t)                         :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_fulvene'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in
    QModel%nsurf    = 2
    QModel%ndim     = 30
    QModel%pot_name = 'fulvene'
    allocate(QModel%omega(QModel%ndim))


    IF (debug) write(out_unit,*) 'Paramters of fulvene'

    
    !> Parameters are from Phys. Chem. Chem. Phys., 2024,26, 1829-1844, doi.org/10.1039/D3CP03964A
    QModel%omega(1)                 =   0.02610_Rkind! , ev
    QModel%omega(2)                 =   0.04530_Rkind! , ev
    QModel%omega(3)                 =   0.06174_Rkind! , ev
    QModel%omega(4)                 =   0.07760_Rkind! , ev
    QModel%omega(5)                 =   0.08560_Rkind! , ev
    QModel%omega(6)                 =   0.08750_Rkind! , ev
    QModel%omega(7)                 =   0.09489_Rkind! , ev
    QModel%omega(8)                 =   0.09766_Rkind! , ev
    QModel%omega(9)                 =   0.10603_Rkind! , ev
    QModel%omega(10)                =   0.11074_Rkind! , ev
    QModel%omega(11)                =   0.11154_Rkind! , ev
    QModel%omega(12)                =   0.11242_Rkind! , ev
    QModel%omega(13)                =   0.11739_Rkind! , ev
    QModel%omega(14)                =   0.12811_Rkind! , ev
    QModel%omega(15)                =   0.12909_Rkind! , ev
    QModel%omega(16)                =   0.14578_Rkind! , ev
    QModel%omega(17)                =   0.14629_Rkind! , ev
    QModel%omega(18)                =   0.16851_Rkind! , ev
    QModel%omega(19)                =   0.18011_Rkind! , ev
    QModel%omega(20)                =   0.18285_Rkind! , ev
    QModel%omega(21)                =   0.19361_Rkind! , ev
    QModel%omega(22)                =   0.20310_Rkind! , ev
    QModel%omega(23)                =   0.20964_Rkind! , ev
    QModel%omega(24)                =   0.22047_Rkind! , ev
    QModel%omega(25)                =   0.41409_Rkind! , ev
    QModel%omega(26)                =   0.42437_Rkind! , ev
    QModel%omega(27)                =   0.42759_Rkind! , ev
    QModel%omega(28)                =   0.42790_Rkind! , ev
    QModel%omega(29)                =   0.43081_Rkind! , ev
    QModel%omega(30)                =   0.43248_Rkind! , ev
    QModel%E1                      =   0.00000_Rkind! , ev
    QModel%E2                      =   4.16142_Rkind! , ev
    QModel%kappa1_6                =  -0.00799_Rkind! , ev
    QModel%kappa1_13               =   0.00608_Rkind! , ev
    QModel%kappa1_15               =   0.00487_Rkind! , ev
    QModel%kappa1_17               =  -0.01463_Rkind! , ev
    QModel%kappa1_20               =  -0.01265_Rkind! , ev
    QModel%kappa1_21               =  -0.00096_Rkind! , ev
    QModel%kappa1_22               =  -0.00421_Rkind! , ev
    QModel%kappa1_24               =  -0.00090_Rkind! , ev
    QModel%kappa1_25               =   0.00464_Rkind! , ev
    QModel%kappa1_28               =   0.01287_Rkind! , ev
    QModel%kappa1_30               =   0.05083_Rkind! , ev
    QModel%kappa2_6                =  -0.16463_Rkind! , ev
    QModel%kappa2_13               =   0.16902_Rkind! , ev
    QModel%kappa2_15               =   0.20870_Rkind! , ev
    QModel%kappa2_17               =  -0.19129_Rkind! , ev
    QModel%kappa2_20               =   0.12386_Rkind! , ev
    QModel%kappa2_21               =  -0.06968_Rkind! , ev
    QModel%kappa2_22               =  -0.38302_Rkind! , ev
    QModel%kappa2_24               =   0.49293_Rkind! , ev
    QModel%kappa2_25               =  -0.05766_Rkind! , ev
    QModel%kappa2_28               =   0.01640_Rkind! , ev
    QModel%kappa2_30               =   0.03264_Rkind! , ev
    QModel%lambda1_2_2             =  -0.02235_Rkind! , ev
    QModel%lambda1_2_9             =  -0.14659_Rkind! , ev
    QModel%lambda1_2_14            =   0.04999_Rkind! , ev
    QModel%lambda1_2_16            =   0.02678_Rkind! , ev
    QModel%lambda1_2_18            =  -0.07546_Rkind! , ev
    QModel%lambda1_2_19            =  -0.14017_Rkind! , ev
    QModel%lambda1_2_23            =   0.04228_Rkind! , ev
    QModel%lambda1_2_26            =   0.00122_Rkind! , ev
    QModel%lambda1_2_27            =   0.00416_Rkind! , ev
    QModel%lambda1_2_29            =  -0.00788_Rkind! , ev

    IF (debug) write(out_unit,*) 'init d0GGdef of fulvene'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_fulvene
!> @brief Subroutine wich prints the current QML_XXX parameters.
!!
!! @param QModel            CLASS(QML_XXX_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_fulvene(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_fulvene_t),   intent(in) :: QModel
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
    write(nio,*) 'omega (ev)',i,QModel%omega(i)
    end do

    write(nio,*)'--------------------------------------'
    write(nio,*) 'E1-2 (ev)',QModel%E1,QModel%E2
    write(nio,*)'--------------------------------------'            
    ! Printing the variables
    WRITE(nio,*) 'kappa1_6 = ', QModel%kappa1_6
    WRITE(nio,*) 'kappa1_13 = ', QModel%kappa1_13
    WRITE(nio,*) 'kappa1_15 = ', QModel%kappa1_15
    WRITE(nio,*) 'kappa1_17 = ', QModel%kappa1_17
    WRITE(nio,*) 'kappa1_20 = ', QModel%kappa1_20
    WRITE(nio,*) 'kappa1_21 = ', QModel%kappa1_21
    WRITE(nio,*) 'kappa1_22 = ', QModel%kappa1_22
    WRITE(nio,*) 'kappa1_24 = ', QModel%kappa1_24
    WRITE(nio,*) 'kappa1_25 = ', QModel%kappa1_25
    WRITE(nio,*) 'kappa1_28 = ', QModel%kappa1_28
    WRITE(nio,*) 'kappa1_30 = ', QModel%kappa1_30
    write(nio,*)'--------------------------------------'
    WRITE(nio,*) 'kappa2_6 = ', QModel%kappa2_6
    WRITE(nio,*) 'kappa2_13 = ', QModel%kappa2_13
    WRITE(nio,*) 'kappa2_15 = ', QModel%kappa2_15
    WRITE(nio,*) 'kappa2_17 = ', QModel%kappa2_17
    WRITE(nio,*) 'kappa2_20 = ', QModel%kappa2_20
    WRITE(nio,*) 'kappa2_21 = ', QModel%kappa2_21
    WRITE(nio,*) 'kappa2_22 = ', QModel%kappa2_22
    WRITE(nio,*) 'kappa2_24 = ', QModel%kappa2_24
    WRITE(nio,*) 'kappa2_25 = ', QModel%kappa2_25
    WRITE(nio,*) 'kappa2_28 = ', QModel%kappa2_28
    WRITE(nio,*) 'kappa2_30 = ', QModel%kappa2_30
    write(nio,*)'--------------------------------------'
    WRITE(nio,*) 'lambda1_2_2 = ', QModel%lambda1_2_2
    WRITE(nio,*) 'lambda1_2_9 = ', QModel%lambda1_2_9
    WRITE(nio,*) 'lambda1_2_14 = ', QModel%lambda1_2_14
    WRITE(nio,*) 'lambda1_2_16 = ', QModel%lambda1_2_16
    WRITE(nio,*) 'lambda1_2_18 = ', QModel%lambda1_2_18
    WRITE(nio,*) 'lambda1_2_19 = ', QModel%lambda1_2_19
    WRITE(nio,*) 'lambda1_2_23 = ', QModel%lambda1_2_23
    WRITE(nio,*) 'lambda1_2_26 = ', QModel%lambda1_2_26
    WRITE(nio,*) 'lambda1_2_27 = ', QModel%lambda1_2_27
    WRITE(nio,*) 'lambda1_2_29 = ', QModel%lambda1_2_29
    write(nio,*)'X------------------------------------X'

  END SUBROUTINE Write_QML_fulvene
!> @brief Subroutine wich calculates the XXX potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_XXX_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_fulvene(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_fulvene_t), intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)  !>  DIMENSIONLESS (no units)
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t) :: dnH0,  dnW0_1,dnW0_2,dnW0_3,  dnW1
    real (kind=Rkind),parameter :: au_to_eV      = 27.2114_Rkind
    !< Conversion factor from Hartree to electronvolts

    !H0 Hamiltionian
    dnH0=0.5_Rkind*sum(QModel%omega(:)*dnQ(:)**2 )

    !W1 diagonal
    dnW0_1=QModel%kappa1_6*dnQ(6)+QModel%kappa1_13*dnQ(13)+QModel%kappa1_15*dnQ(15)&
           & +QModel%kappa1_17*dnQ(17)+QModel%kappa1_20*dnQ(20)&
           & +QModel%kappa1_21*dnQ(21)+QModel%kappa1_22*dnQ(22)&
           & +QModel%kappa1_24*dnQ(24)+QModel%kappa1_25*dnQ(25)&
           & +QModel%kappa1_28*dnQ(28)+QModel%kappa1_30*dnQ(30)      
    dnW0_2=QModel%kappa2_6*dnQ(6)+QModel%kappa2_13*dnQ(13)+QModel%kappa2_15*dnQ(15)&
           & +QModel%kappa2_17*dnQ(17)+QModel%kappa2_20*dnQ(20)&
           & +QModel%kappa2_21*dnQ(21)+QModel%kappa2_22*dnQ(22)&
           & +QModel%kappa2_24*dnQ(24)+QModel%kappa2_25*dnQ(25)&
           & +QModel%kappa2_28*dnQ(28)+QModel%kappa2_30*dnQ(30)      

    !Diagonal elements total:
    Mat_OF_PotDia(1,1) = dnH0+dnW0_1+QModel%E1
    Mat_OF_PotDia(2,2) = dnH0+dnW0_2+QModel%E2

    !Upper off-diagonal
    Mat_OF_PotDia(1,2) =  QModel%lambda1_2_2*   dnQ(2) &
                        &+QModel%lambda1_2_9*   dnQ(9) &
                        &+QModel%lambda1_2_14*  dnQ(14)&
                        &+QModel%lambda1_2_16*  dnQ(16)&
                        &+QModel%lambda1_2_18*  dnQ(18)&
                        &+QModel%lambda1_2_19*  dnQ(19)&
                        &+QModel%lambda1_2_23*  dnQ(23)&
                        &+QModel%lambda1_2_26*  dnQ(26)&
                        &+QModel%lambda1_2_27*  dnQ(27)&
                        &+QModel%lambda1_2_29*  dnQ(29)

    !Lower off-diagonal
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)

    Select case (QModel%option)
    case(2)
      Mat_OF_PotDia=Mat_OF_PotDia/au_to_eV
    end select


  END SUBROUTINE EvalPot_QML_fulvene

END MODULE QML_fulvene_m

