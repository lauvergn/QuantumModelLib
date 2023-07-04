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
!> @brief Module which makes the initialization, calculation of the ClH2p potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_ClH2p_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the ClH2p parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE, EXTENDS (QML_Empty_t) ::  QML_ClH2p_t

     PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

     integer                       :: nb_funcModel
     real(kind=Rkind), allocatable :: F(:)
     integer, allocatable          :: tab_func(:,:)
     real(kind=Rkind)              :: E0 = ZERO

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_ClH2p
    PROCEDURE :: Write_QModel    => Write_QML_ClH2p
  END TYPE QML_ClH2p_t

  PUBLIC :: QML_ClH2p_t,Init_QML_ClH2p

  CONTAINS
!> @brief Subroutine which makes the initialization of the ClH2p parameters.
!!
!! @param ClH2pPot           TYPE(QML_ClH2p_t):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_ClH2p(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : file_open2
    USE QMLLib_UtilLib_m, ONLY : make_QMLInternalFileName
    IMPLICIT NONE

    TYPE (QML_ClH2p_t)                           :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !local variable
    integer                         :: err_read,nio_fit,i,j,k,idum
    character (len=50)              :: name_Qref
    character (len=:), allocatable  :: FileName


    integer                  :: ifit,nFit,iFunc

    !-----------------------------------------------------------------
    ! for the namelist
    integer                  :: MinCoupling,MaxCoupling
    real (kind=Rkind)        :: MinNorm,MaxNorm,max_valB
    integer                  :: ndim,nb_B,nb_WB,MR_order
    logical                  :: svd
    real (kind=Rkind)        :: epsi,epsi_inter
    integer                  :: ind_val,nb_val,max_nb
    integer                  :: Col_FOR_WeightOFFit
    real (kind=Rkind)        :: Scal_FOR_WeightOFFit
    integer                  :: nb_Fit


    namelist /nDFitW/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,        &
                      ind_val,nb_val,svd,                             &
                      epsi,epsi_inter,max_nb,ndim,MR_order,           &
                      Col_FOR_WeightOFFit,Scal_FOR_WeightOFFit,       &
                      nb_B,nb_WB,nb_Fit
    !-----------------------------------------------------------------

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_ClH2p'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf    = 1
    QModel%pot_name = 'ClH2p'


    IF (QModel%option < 1 .OR. QModel%option > 6) QModel%option = 6

    SELECT CASE (QModel%option)
    CASE (1)

      QModel%d0GGdef = reshape(                                                 &
        [0.0001817892_Rkind,-0.0000062919_Rkind,0.0000000000_Rkind,             &
        -0.0000062919_Rkind, 0.0002793633_Rkind,0.0000000000_Rkind,             &
        -0.0000000000_Rkind,-0.0000000000_Rkind,0.0002806450_Rkind],shape=[3,3])

      FileName = make_QMLInternalFileName('InternalData/ClH2p/ClH2p_B3LYP_cc-pVTZ.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      QModel%E0 = ZERO

    CASE (2)

      QModel%d0GGdef = reshape(                                                 &
        [0.0005600083_Rkind,-0.0000012817_Rkind,-0.0000062919_Rkind,            &
        -0.0000012817_Rkind, 0.0005600083_Rkind,-0.0000062919_Rkind,            &
        -0.0000062919_Rkind,-0.0000062919_Rkind,0.0001817892_Rkind],shape=[3,3])

        FileName = make_QMLInternalFileName('InternalData/ClH2p/ClH2p_B3LYP_cc-pVTZ.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      QModel%E0 = ZERO

    CASE (3)

      QModel%d0GGdef = reshape(                                                 &
        [0.0001817892_Rkind,-0.0000062919_Rkind,0.0000000000_Rkind,             &
        -0.0000062919_Rkind, 0.0002793633_Rkind,0.0000000000_Rkind,             &
        -0.0000000000_Rkind,-0.0000000000_Rkind,0.0002806450_Rkind],shape=[3,3])

      FileName = make_QMLInternalFileName('InternalData/ClH2p/ClH2p_B3LYP_cc-pVTZ_v2.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      QModel%E0 = -461.06327752734671_Rkind

    CASE (4)

      QModel%d0GGdef = reshape(                                                 &
        [0.0005600083_Rkind,-0.0000012817_Rkind,-0.0000062919_Rkind,            &
        -0.0000012817_Rkind, 0.0005600083_Rkind,-0.0000062919_Rkind,            &
        -0.0000062919_Rkind,-0.0000062919_Rkind,0.0001817892_Rkind],shape=[3,3])

        FileName = make_QMLInternalFileName('InternalData/ClH2p/ClH2p_B3LYP_cc-pVTZ_v2.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      QModel%E0 = -461.06327752734671_Rkind

    CASE (5)

      QModel%d0GGdef = reshape(                                                 &
        [0.0001843648_Rkind,-0.0000063402_Rkind,0.0000000000_Rkind,             &
        -0.0000063402_Rkind, 0.0002794155_Rkind,0.0000000000_Rkind,             &
        -0.0000000000_Rkind,-0.0000000000_Rkind,0.0002805928_Rkind],shape=[3,3])

      FileName = make_QMLInternalFileName('InternalData/ClH2p/ClH2p_CCSDTF12_cc-pVTZ.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      QModel%E0 = -460.59052688_Rkind

    CASE (6)

      QModel%d0GGdef = reshape(                                                 &
        [0.0005600083_Rkind,-0.0000011773_Rkind,-0.0000063402_Rkind,            &
        -0.0000011773_Rkind, 0.0005600083_Rkind,-0.0000063402_Rkind,            &
        -0.0000063402_Rkind,-0.0000063402_Rkind, 0.0001843648_Rkind],shape=[3,3])

        FileName = make_QMLInternalFileName('InternalData/ClH2p/ClH2p_CCSDTF12_cc-pVTZ.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      QModel%E0 = -460.59052688_Rkind

    CASE Default

      write(out_unit,*) ' ERROR in Init_QML_ClH2p '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1,2,3,4,5,6'
      STOP 'ERROR in Init_QML_ClH2p: wrong option'

    END SELECT
    IF (debug) write(out_unit,*) 'FileName',FileName
    flush(out_unit)


    read(nio_fit,*) name_Qref,nFit

    IF (debug) write(out_unit,*) 'nFit',nFit
    flush(out_unit)

    ! first the number of points
    QModel%nb_funcModel = 0
    DO ifit=1,nFit

      read(nio_fit,nDFitW,IOSTAT=err_read)
      IF (debug) write(out_unit,nDFitW)
      flush(out_unit)

      IF (err_read < 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' End-of-file or End-of-record'
        write(out_unit,*) ' The namelist "nDFitW" is probably absent'
        write(out_unit,*) ' from the file: ',trim(adjustl(FileName))
        write(out_unit,*) ' ERROR in ',name_sub
        STOP
      ELSE IF (err_read > 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Some parameter name of the namelist "nDFitW" are probaly wrong'
        write(out_unit,*) ' in the file: ',trim(adjustl(FileName))
        write(out_unit,*) ' It should never append !!'
        write(out_unit,*) ' Check the fortran'
        write(out_unit,*) ' ERROR in ',name_sub
        STOP
      END IF

      QModel%ndim         = ndim
      ! nb_WB must be used instead of nb_B because some basis functions might be skip
      QModel%nb_funcModel = QModel%nb_funcModel + nb_WB
    END DO

    IF (debug) write(out_unit,*) 'nb_funcModel',QModel%nb_funcModel

    close(nio_fit)
    CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)

    allocate(QModel%Qref(QModel%ndim))
    allocate(QModel%F(QModel%nb_funcModel))
    allocate(QModel%tab_func(QModel%ndim,QModel%nb_funcModel))

    iFunc = 0
    DO ifit=1,nFit

      read(nio_fit,nDFitW)
      IF (debug) write(out_unit,nDFitW)
      IF (err_read < 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' End-of-file or End-of-record'
        write(out_unit,*) ' The namelist "nDFitW" is probably absent'
        write(out_unit,*) ' from the file: ',trim(adjustl(FileName))
        write(out_unit,*) ' ERROR in ',name_sub
        STOP
      ELSE IF (err_read > 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Some parameter name of the namelist "nDFitW" are probaly wrong'
        write(out_unit,*) ' in the file: ',trim(adjustl(FileName))
        write(out_unit,*) ' It should never append !!'
        write(out_unit,*) ' Check the fortran'
        write(out_unit,*) ' ERROR in ',name_sub
        STOP
      END IF

      read(nio_fit,*) name_Qref,QModel%Qref(:)
      read(nio_fit,*)
      read(nio_fit,*)
      read(nio_fit,*)

      DO i=1,nb_WB
        iFunc = iFunc + 1
        read(nio_fit,*) idum,QModel%tab_func(:,iFunc),QModel%F(iFunc)
        QModel%tab_func(3,iFunc) = 2*QModel%tab_func(3,iFunc) ! for even momomials
        IF (debug) write(out_unit,*) 'iFunc,tab_func,F',                       &
                                  iFunc,QModel%tab_func(:,iFunc),QModel%F(iFunc)
        flush(out_unit)
      END DO

    END DO
    flush(out_unit)

    close(nio_fit)
    deallocate(FileName)


    IF (debug) write(out_unit,*) 'init Q0 of ClH2p'
    allocate(QModel%Q0(QModel%ndim))
    CALL get_Q0_QML_ClH2p(QModel%Q0,QModel,option=0)
    IF (debug) write(out_unit,*) 'QModel%Q0',QModel%Q0

    IF (debug) write(out_unit,*) 'init d0GGdef of ClH2p'

    SELECT CASE (QModel%option)
    CASE (1,3)
      IF (QModel%PubliUnit) THEN
        write(out_unit,*) 'PubliUnit=.TRUE.,  Q:[Rad,Bohr,Bohr], Energy: [Hartree]'
      ELSE
        write(out_unit,*) 'PubliUnit=.FALSE., Q:[Rad,Bohr,Bohr], Energy: [Hartree]'
      END IF
    CASE (2)
      IF (QModel%PubliUnit) THEN
        write(out_unit,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad], Energy: [Hartree]'
      ELSE
        write(out_unit,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad]:, Energy: [Hartree]'
      END IF
    END SELECT
    flush(out_unit)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_ClH2p
!> @brief Subroutine wich prints the QML_ClH2p parameters.
!!
!! @param QModel            CLASS(QML_ClH2p_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_QML_ClH2p(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_ClH2p_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'ClH2p current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '    R2  \ a                            '
    write(nio,*) '         Cl------------H               '
    write(nio,*) '              R1                       '
    write(nio,*) '                                       '
    write(nio,*) '  Coordinates (option 1,3,5):          '
    write(nio,*) '  Q(1) = a          (Radian)           '
    write(nio,*) '  Q(2) = 1/2(R1+R2) (Bohr)             '
    write(nio,*) '  Q(3) = 1/2(R1-R2) (Bohr)             '
    write(nio,*) '                                       '
    write(nio,*) '  Coordinates (option 2,4,6):          '
    write(nio,*) '  Q(1) = R1         (Bohr)             '
    write(nio,*) '  Q(2) = R2         (Bohr)             '
    write(nio,*) '  Q(3) = a          (Radian)           '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) ' Options (default: 6)                  '
    write(nio,*) ' 1,2: B3LYP/cc-pVTZ 1st version (do not use)'
    write(nio,*) ' 3,4: B3LYP/cc-pVTZ 2d  version        '
    write(nio,*) ' 5,6: CCSD(T)-F12b/cc-pVTZ-F12         '
    write(nio,*) '                                       '
    write(nio,*) 'Ref: unpublished'
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)
    write(nio,*) '---------------------------------------'

    SELECT CASE (QModel%option)

    CASE (1)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: B3LYP/cc-pVTZ                 '
      write(nio,*) '                                       '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = a    = 1.6525859462 (Rad)     '
      write(nio,*) '  Q(2) = r+   = 2.4849898659 (Bohr)    '
      write(nio,*) '  Q(3) = r-   = 0.           (Bohr)    '
      write(nio,*) '                                       '
      write(nio,*) '  V = 0.0004486447 Hartree             '
      write(nio,*) ' grad(:) =[ 0.0000558107,-0.0004452755,'
      write(nio,*) '            0.0000000000]              '
      write(nio,*) ' hess    =[0.158044529,0.009093148,0.000000000'
      write(nio,*) '           0.009093148,0.520862288,0.000000000'
      write(nio,*) '           0.000000000,0.000000000,0.511935490]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '

    CASE (2)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: B3LYP/cc-pVTZ                 '
      write(nio,*) '                                       '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = R1   = 2.4849898659 (Bohr)    '
      write(nio,*) '  Q(2) = R2   = 2.4849898659 (Bohr)    '
      write(nio,*) '  Q(3) = a    = 1.6525859462 (Rad)     '
      write(nio,*) '                                       '
      write(nio,*) '  V = 0.0004486447 Hartree             '
      write(nio,*) ' grad(:) =[-0.0002226378,-0.0002226378,'
      write(nio,*) '            0.0000558107]              '
      write(nio,*) ' hess    =[0.258199445,0.002231699,0.004546574'
      write(nio,*) '           0.002231699,0.258199445,0.004546574'
      write(nio,*) '           0.004546574,0.004546574,0.158044529]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '
    CASE (3)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: B3LYP/cc-pVTZ_v2              '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = a    = 1.6525859462 (Rad)     '
      write(nio,*) '  Q(2) = r+   = 2.4849898659 (Bohr)    '
      write(nio,*) '  Q(3) = r-   = 0.           (Bohr)    '
      write(nio,*) '                                       '
      write(nio,*) '  V = -461.0632775273 Hartree          '
      write(nio,*) ' grad(:) =[ 0.0000000000, 0.0000726996,'
      write(nio,*) '            0.0000000000]              '
      write(nio,*) ' hess    =[0.165834373,0.008873536,0.000000000'
      write(nio,*) '           0.008873536,0.544253313,0.000000000'
      write(nio,*) '           0.000000000,0.000000000,0.534634932]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '
    CASE (4)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: B3LYP/cc-pVTZ_v2              '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = R1   = 2.4849898659 (Bohr)    '
      write(nio,*) '  Q(2) = R2   = 2.4849898659 (Bohr)    '
      write(nio,*) '  Q(3) = a    = 1.6525859462 (Rad)     '
      write(nio,*) '                                       '
      write(nio,*) '  V = -461.0632775273 Hartree          '
      write(nio,*) ' grad(:) =[ 0.0000363498, 0.0000363498,'
      write(nio,*) '            0.0000000000]              '
      write(nio,*) ' hess    =[0.269721061,0.002405595,0.004436768'
      write(nio,*) '           0.002405595,0.269721061,0.004436768'
      write(nio,*) '           0.004436768,0.004436768,0.165834373]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '
    CASE (5)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: CCSD(T)-F12b/cc-pVTZ-F12      '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = a    = 1.6459159559 (Rad)     '
      write(nio,*) '  Q(2) = r+   = 2.4673412575 (Bohr)    '
      write(nio,*) '  Q(3) = r-   = 0.           (Bohr)    '
      write(nio,*) '                                       '
      write(nio,*) '  V = -460.5905268800 Hartree          '
      write(nio,*) ' grad(:) =[ 0.0000060177, 0.0000283675,'
      write(nio,*) '            0.0000000000]              '
      write(nio,*) ' hess    =[0.167709143,0.010804767,0.000000000'
      write(nio,*) '           0.010804767,0.567302406,0.000000000'
      write(nio,*) '           0.000000000,0.000000000,0.560708185]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '
    CASE (6)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: CCSD(T)-F12b/cc-pVTZ-F12      '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = R1   = 2.4673412575 (Bohr)    '
      write(nio,*) '  Q(2) = R2   = 2.4673412575 (Bohr)    '
      write(nio,*) '  Q(3) = a    = 1.6459159559 (Rad)     '
      write(nio,*) '                                       '
      write(nio,*) '  V = -460.5905268800 Hartree          '
      write(nio,*) ' grad(:) =[ 0.0000141837, 0.0000141837,'
      write(nio,*) '            0.0000060177]              '
      write(nio,*) ' hess    =[0.282002648,0.001648555,0.005402383'
      write(nio,*) '           0.001648555,0.282002648,0.005402383'
      write(nio,*) '           0.005402383,0.005402383,0.167709143]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '
    CASE Default
        write(out_unit,*) ' ERROR in write_QModel '
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1,2,3,4,5,6'
        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end ClH2p current parameters'

  END SUBROUTINE Write_QML_ClH2p

  SUBROUTINE get_Q0_QML_ClH2p(Q0,QModel,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (QML_ClH2p_t),          intent(in)    :: QModel
    integer,                     intent(in)    :: option

    IF (size(Q0) /= 3) THEN
      write(out_unit,*) ' ERROR in get_Q0_QML_ClH2p '
      write(out_unit,*) ' The size of Q0 is not ndim=3: '
      write(out_unit,*) ' size(Q0)',size(Q0)
      STOP 'ERROR in get_Q0_QML_ClH2p: wrong Q0 size'
    END IF


    SELECT CASE (QModel%option)
    CASE (1,3,5) ! a,r+,r-
      Q0(:) = QModel%Qref(:)

    CASE (2,4,6) ! R1,R2,a
      Q0(:) = [QModel%Qref(2),QModel%Qref(2),QModel%Qref(1)]

    CASE Default
      write(out_unit,*) ' ERROR in get_Q0_QML_ClH2p '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1,2,3,4,5,6'
      STOP 'ERROR in get_Q0_QML_ClH2p: wrong option'
    END SELECT

  END SUBROUTINE get_Q0_QML_ClH2p
!> @brief Subroutine wich calculates the ClH2p potential (unpublished model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_ClH2p_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_ClH2p(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_ClH2p_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable :: dnQsym(:)


    SELECT CASE (QModel%option)
    CASE (1,3,5) ! a,r+,r-
      CALL EvalPot1_QML_ClH2p(QModel,Mat_OF_PotDia,dnQ,nderiv)

    CASE (2,4,6) ! R1,R2,a
      dnQsym = [dnQ(3), HALF*(dnQ(1)+dnQ(2)), HALF*(dnQ(1)-dnQ(2))]
      CALL EvalPot1_QML_ClH2p(QModel,Mat_OF_PotDia,dnQsym,nderiv)

    CASE Default
      write(out_unit,*) ' ERROR in EvalPot_QML_ClH2p '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1,2,3,4,5,6'
      STOP 'ERROR in EvalPot_QML_ClH2p: wrong option'
    END SELECT


  END SUBROUTINE EvalPot_QML_ClH2p

  SUBROUTINE EvalPot1_QML_ClH2p(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_ClH2p_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable :: DQ(:,:)
    TYPE (dnS_t)              :: Vtemp
    integer                   :: i,j,max_exp

    !write(out_unit,*) ' sub EvalPot1_QML_ClH2p' ; flush(6)
    max_exp = maxval(QModel%tab_func)
    !write(out_unit,*) ' max_exp',max_exp ; flush(6)

    allocate(DQ(QModel%ndim,max_exp))
    DQ(:,1) = dnQ(:) - QModel%Qref(:)
    DO j=2,size(DQ,dim=2)
      DQ(:,j) = DQ(:,j-1) * DQ(:,1)
    END DO
    !write(out_unit,*) ' DQ done' ; flush(6)


    Mat_OF_PotDia(1,1) = QModel%E0
    Vtemp = Mat_OF_PotDia(1,1) ! to have a correct initialization
    !write(out_unit,*) ' Vtemp init done' ; flush(6)



      DO j=1,QModel%nb_funcModel
        Vtemp = QModel%F(j)
        DO i=1,QModel%ndim
          !Vtemp = Vtemp * DQ(i,1)**QModel%tab_func(i,j)
          IF (QModel%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,QModel%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO
      !CALL Write_dnS(Mat_OF_PotDia(1,1),nio=out_unit)


!-----------------------------------------------------------------------!

   CALL dealloc_dnS(Vtemp)
   CALL dealloc_dnS(DQ)

   !write(out_unit,*) ' end EvalPot1_QML_ClH2p' ; flush(6)

 END SUBROUTINE EvalPot1_QML_ClH2p
END MODULE QML_ClH2p_m
