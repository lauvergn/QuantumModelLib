!===========================================================================
!===========================================================================
!This file is part of QuantumModelLib (QML).
!==============================================================================
! MIT License
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in al
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
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, Fran
![2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France
![3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France
![4]: Durham University, Durham, UK
!* Originally, it has been developed during the Quantum-Dynamics E-CAM project 
!     https://www.e-cam2020.eu/quantum-dynamics
!
!===========================================================================
!===========================================================================

!> @brief Module which makes the initialization, calculation of the DoubleWell pote
!!
!> @author Rabiou ISSA
!! @date 18/11/2025
!!
MODULE QML_DoubleWell_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

  ! default parameter values (Rkind precision)
  real(kind=Rkind), parameter :: DEFAULT_a = ONE
  real(kind=Rkind), parameter :: DEFAULT_b = ONE

  !> @brief Derived type in which the DoubleWell parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) :: QML_DoubleWell_t
    PRIVATE
    real(kind=Rkind), PUBLIC :: a = DEFAULT_a
    real(kind=Rkind), PUBLIC :: b = DEFAULT_b
  CONTAINS
    PROCEDURE :: EvalPot_QModel      => EvalPot_QML_DoubleWell
    PROCEDURE :: Write_QModel        => Write_QML_DoubleWell
    PROCEDURE :: Ref_FOR_Test_QModel => Ref_FOR_Test_QML_DoubleWell
  END TYPE QML_DoubleWell_t

  PUBLIC :: QML_DoubleWell_t, Init_QML_DoubleWell, Read_QML_DoubleWell, Write_QML_DoubleWell

CONTAINS

  !======================================================================
  ! Init_QML_DoubleWell
  !======================================================================
  FUNCTION Init_QML_DoubleWell(QModel_in, read_param, nio_param_file) RESULT(QModel)
    IMPLICIT NONE

    TYPE(QML_DoubleWell_t) :: QModel

    TYPE(QML_Empty_t), INTENT(IN) :: QModel_in
    integer, intent(in)           :: nio_param_file
    logical, intent(in)           :: read_param

    ! debugging
    character(len=*), parameter :: name_sub = 'Init_QML_DoubleWell'
    logical, parameter :: debug = .FALSE.

    ! copy base empty type contents
    QModel%QML_Empty_t = QModel_in

    ! model identification
    QModel%pot_name = 'DoubleWell'
    QModel%ndim     = 1
    QModel%nsurf    = 1

    ! default parameters
    IF (debug) write(out_unit,*) 'init default DoubleWell parameters'
    QModel%a = DEFAULT_a
    QModel%b = DEFAULT_b

    ! default Q0 (test coordinate)
    QModel%Q0 = [ZERO] 

    ! default metric / effective mass (identity by default)
    QModel%d0GGdef = reshape([ ONE ], shape=[1,1])

    ! optionally read parameters from namelist
    IF (read_param) THEN
      CALL Read_QML_DoubleWell(QModel, nio_param_file)
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of DoubleWell'
    QModel%Q0 = [ZERO] 

    IF (debug) THEN
      write(out_unit,*) 'Init_QML_DoubleWell: a=',QModel%a,' b=',QModel%b
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_DoubleWell

  !======================================================================
  ! Read_QML_DoubleWell
  !  Read namelist: /DoubleWell/ a, b
  !======================================================================
  SUBROUTINE Read_QML_DoubleWell(QModel, nio)
    IMPLICIT NONE

    TYPE(QML_DoubleWell_t), INTENT(INOUT) :: QModel
    INTEGER, INTENT(IN) :: nio

    real(kind=Rkind) :: a, b
    integer :: err_read

    namelist /DoubleWell/ a, b

    ! defaults from current model
    a = QModel%a
    b = QModel%b

    write(out_unit,*) 'read DoubleWell namelist ...'
    flush(out_unit)

    read(nio, nml=DoubleWell, IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_DoubleWell'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "DoubleWell" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_DoubleWell'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_DoubleWell'
      write(out_unit,*) ' Some parameter names of the namelist "DoubleWell" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=DoubleWell)
      STOP ' ERROR in Read_QML_DoubleWell'
    END IF

    ! commit read values
    QModel%a = a
    QModel%b = b

  END SUBROUTINE Read_QML_DoubleWell

  !======================================================================
  ! Write_QML_DoubleWell
  !======================================================================
  SUBROUTINE Write_QML_DoubleWell(QModel, nio)
    IMPLICIT NONE

    CLASS(QML_DoubleWell_t), INTENT(IN) :: QModel
    INTEGER, INTENT(IN) :: nio

    write(nio,*) 'DoubleWell potential expression:'
    write(nio,*)
    write(nio,*) ' V(s) = -1/2 * a * s^2 + 1/4 * b * s^4'
    write(nio,*)
    write(nio,*) ' The minima are in s_m=+/- sqrt(a/b), V(s_m)=-a^2/(4b)'
    write(nio,*) ' the curvature at the minima is V(s_m)"=2a'

    write(nio,*) 'end DoubleWell potential expression'
    write(nio,*)
    write(nio,*) 'DoubleWell current parameters:'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) '  pot_name : ', TRIM(QModel%pot_name)
    write(nio,*) '  ndim     : ', QModel%ndim
    write(nio,*) '  nsurf    : ', QModel%nsurf
    write(nio,*) '  a        : ', QModel%a
    write(nio,*) '  b        : ', QModel%b
    write(nio,*)
    write(nio,*) 'end DoubleWell current parameters'
    write(nio,*)
  END SUBROUTINE Write_QML_DoubleWell

  !======================================================================
  ! EvalPot_QML_DoubleWell
  ! Compute potential (and its AD derivatives when dnS_t supports it)
  ! Mat_OF_PotDia is TYPE(dnS_t)(:,:), dnQ is TYPE(dnS_t)(:)
  ! nderiv controls derivative order (kept for API compatibility)
  !======================================================================
  SUBROUTINE EvalPot_QML_DoubleWell(QModel, Mat_OF_PotDia, dnQ, nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_DoubleWell_t), INTENT(IN)    :: QModel
    TYPE(dnS_t),         INTENT(INOUT)     :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         INTENT(IN)        :: dnQ(:)
    INTEGER,             INTENT(IN)        :: nderiv

    TYPE(dnS_t) :: s
    integer :: ns

    ! debugging
    character(len=*), parameter :: name_sub='EvalPot_QML_DoubleWell'
    logical, parameter :: debug = .FALSE.

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ', name_sub
      write(out_unit,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) ' nderiv: ', nderiv
      flush(out_unit)
    END IF

    ns = SIZE(Mat_OF_PotDia,1)
    IF (ns < 1) RETURN

    ! reaction coordinate (assume first coordinate)
    s = dnQ(1)

    ! V(s) = -1/2 * a * s^2 + 1/4 * b * s^4 
    Mat_OF_PotDia(1,1) = - HALF * QModel%a * s**2 + FOURTH * QModel%b * s **4

    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia(1,1):'
      CALL Write_dnS( Mat_OF_PotDia(1,1), out_unit )
      write(out_unit,*)
      write(out_unit,*) 'END ', name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE EvalPot_QML_DoubleWell

    SUBROUTINE Ref_FOR_Test_QML_DoubleWell(QModel,err,Q0,dnMatV,d0GGdef,nderiv)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_DoubleWell_t),intent(in)              :: QModel

    integer,                intent(inout)           :: err

    integer,                intent(in)              :: nderiv
    real (kind=Rkind),      intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),         intent(inout), optional :: dnMatV

    real (kind=Rkind),      intent(inout), optional :: d0GGdef(:,:)

    real (kind=Rkind), allocatable :: Mat(:,:)
    integer        :: i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Ref_FOR_Test_QML_DoubleWell'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
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
      Q0(:) = sqrt(QModel%a/QModel%b)
    END IF

    IF (present(dnMatV)) THEN
      err = 0
      CALL alloc_dnMat(dnMatV,nsurf=QModel%nsurf,nVar=QModel%ndim,nderiv=nderiv)

      allocate(Mat(QModel%nsurf,QModel%nsurf))

      IF (nderiv >= 0) THEN ! no derivative
        Mat = -QModel%a**2/(FOUR*QModel%b)
        CALL Mat_wADDTO_dnMat2_ider(Mat,ONE,dnMatV)
      END IF

      IF (nderiv >= 1) THEN ! 1st order derivatives
        Mat = ZERO
        CALL Mat_wADDTO_dnMat2_ider(Mat,ONE,dnMatV,ider=[1])
      END IF

      IF (nderiv >= 2) THEN ! 2d order derivatives
        Mat = TWO*QModel%a
        CALL Mat_wADDTO_dnMat2_ider(Mat,ONE,dnMatV,ider=[1,1])
      END IF

      IF (nderiv >= 3) THEN ! 3d derivative
        Mat = SIX* sqrt(QModel%a*QModel%b)
        CALL Mat_wADDTO_dnMat2_ider(Mat,ONE,dnMatV,ider=[0,0,0]) ! all 3d order derivatives
      END IF

      IF (allocated(Mat)) deallocate(Mat)
    END IF

    IF (present(d0GGdef)) d0GGdef = ONE


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Ref_FOR_Test_QML_DoubleWell
END MODULE QML_DoubleWell_m