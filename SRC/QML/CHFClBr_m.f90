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

    STOP 'CHFClBr QFF potential is removed!'
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
      Q0(:) = [0.1_Rkind,-0.1_Rkind,0.2_Rkind,-0.2_Rkind,0.3_Rkind,-0.3_Rkind,0.4_Rkind, &
              -0.4_Rkind,0.5_Rkind]
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
