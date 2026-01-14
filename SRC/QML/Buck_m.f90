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
!> @brief Module which makes the initialization, calculation of the Buckingham potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
MODULE QML_Buck_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Buckingham parameters are set-up.
!> @brief V(R) = A.Exp(-B*r)-C/r^6
!> @brief Default parameters for Ar-Ar
!> @brief Reference: R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283. doi:10.1098/rspa.1938.0173.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param A,B,C              real: Buckingham parameters
  TYPE, EXTENDS (QML_Empty_t) :: QML_Buck_t
     private
     real (kind=Rkind) :: A   = 387.63744459726228783977_Rkind  ! for Ar2 (eq 27)
     real (kind=Rkind) :: B   =   1.93837805257707347985_Rkind  ! for Ar2 (eq 27)
     real (kind=Rkind) :: C   = 106.54483566475760255666_Rkind  ! for Ar2 (eq 27)
     ! A= 1.69 10^-8 erg      = 387.63744459726228783977 Hartree
     ! B= 1/0.273 A^-1        = 1.93837805257707347985 bohr^-1
     ! C= 102 10^-12 erg A^-6 = 106.54483566475760255666 Hartree bohr^-6
     real (kind=Rkind) :: mu  = 36423.484024390622_Rkind        ! au
  CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_Buck
    PROCEDURE :: Write_QModel     => Write_QML_Buck
    PROCEDURE :: RefValues_QModel => RefValues_QML_Buck
  END TYPE QML_Buck_t

  PUBLIC :: QML_Buck_t,Init_QML_Buck,Init0_QML_Buck,Write_QML_Buck,QML_dnBuck

CONTAINS
!> @brief Subroutine which makes the initialization of the Buckingham parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel          TYPE(QML_Buck_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param A,B,C              real (optional):    parameters
  FUNCTION Init_QML_Buck(QModel_in,read_param,nio_param_file,A,B,C) RESULT(QModel)
    IMPLICIT NONE

    TYPE (QML_Buck_t)                            :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: A,B,C

    !local variables
    real (kind=Rkind) :: R0

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Buck'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    !Default for Ar-Ar
    CALL Init0_QML_Buck(QModel,A=387.63744459726228783977_Rkind,      &
                               B=1.93837805257707347985_Rkind,        &
                               C=106.54483566475760255666_Rkind,      &
                               model_name='buck')

    IF (read_param) THEN
      CALL Read_QML_Buck(QModel,nio_param_file)
    ELSE
      IF (present(A))   QModel%A = A
      IF (present(B))   QModel%B = B
      IF (present(C))   QModel%C = C
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of Buck'
    CALL get_Q0_QML_Buck(R0,QModel) ! no analytical value
    QModel%Q0 = [R0]

    IF (debug) write(out_unit,*) 'init d0GGdef of Buck'
    QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])

    IF (debug) THEN
      write(out_unit,*) 'A,B,C: ',QModel%A,QModel%B,QModel%C
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Buck
  SUBROUTINE Init0_QML_Buck(QModel,A,B,C,model_name)
    IMPLICIT NONE

    TYPE (QML_Buck_t),            intent(inout)   :: QModel
    real (kind=Rkind), optional,   intent(in)      :: A,B,C
    character (len=*), optional,   intent(in)      :: model_name

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init0_QML_Buck'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%In_a_Model = .TRUE.

    QModel%ndim       = 1
    QModel%nsurf      = 1
    QModel%pot_name   = 'Buck'
    IF (present(model_name)) QModel%pot_name = model_name

    IF (debug) write(out_unit,*) 'init Buck parameters (A,B,C), if present'

    IF (present(A))   QModel%A   = A
    IF (present(B))   QModel%B   = B
    IF (present(C))   QModel%C   = C

    IF (debug) THEN
      write(out_unit,*) 'A,B,C: ',QModel%A,QModel%B,QModel%C
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Init0_QML_Buck
!> @brief Subroutine wich reads the Buckingham parameters with a namelist.
!!   This can be called only from the "Init_QML_Buck" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel          TYPE(QML_Buck_t):    derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to read the parameters.
  SUBROUTINE Read_QML_Buck(QModel,nio)
    IMPLICIT NONE

    TYPE (QML_Buck_t), intent(inout) :: QModel
    integer,            intent(in)    :: nio

    real (kind=Rkind) :: A,B,C
    integer           :: err_read
    namelist /Buck/ A,B,C

    A = QModel%A
    B = QModel%B
    C = QModel%C

    read(nio,nml=Buck,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Buck'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "Buck" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Buck'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Buck'
      write(out_unit,*) ' Some parameter names of the namelist "Buck" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=Buck)
      STOP ' ERROR in Read_QML_Buck'
    END IF

    !write(out_unit,nml=Buck)
    QModel%A = A
    QModel%B = B
    QModel%C = C

  END SUBROUTINE Read_QML_Buck
!> @brief Subroutine wich prints the Buckingham parameters (from his publication) + the ones used
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel          TYPE(QML_Buck_t):    derived type with the Buckingham parameters.
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write_QML_Buck(QModel,nio)
    IMPLICIT NONE

    CLASS (QML_Buck_t), intent(in) :: QModel
    integer,             intent(in) :: nio

    write(nio,*) 'Buckingham default parameters:'
    write(nio,*) ' For Ar-Ar, eq 27 of reference:'
    write(nio,*) '  R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283.'
    write(nio,*)
    write(nio,*) 'Buckingham(r) = A.Exp(-B*r)-C/r^6'
    write(nio,*)
    write(nio,*) '   A= 1.69 10^-8 erg      = 387.63744459726228783977 Hartree'
    write(nio,*) '   B= 1/0.273 A^-1        = 1.93837805257707347985   bohr^-1'
    write(nio,*) '   C= 102 10^-12 erg A^-6 = 106.54483566475760255666 Hartree bohr^-6'
    write(nio,*)
    write(nio,*) 'Value at: R=7.0 Bohr'
    write(nio,*) 'V        = -0.000409 Hartree'
    write(nio,*) 'gradient = -0.000186'
    write(nio,*) 'hessian  =  0.001088'
    write(nio,*)
    write(nio,*) 'end Buckingham parameters'
    write(nio,*) 'Buckingham current parameters:'
    write(nio,*)
    write(nio,*) '    V(R) = A.Exp(-B.r) - C/r^6'
    write(nio,*) '  A:   ',QModel%A
    write(nio,*) '  B:   ',QModel%B
    write(nio,*) '  C:   ',QModel%C
    write(nio,*)
    write(nio,*) 'end Buckingham current parameters'

  END SUBROUTINE Write_QML_Buck

  SUBROUTINE get_Q0_QML_Buck(R0,QModel)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: R0
    TYPE (QML_Buck_t),           intent(in)    :: QModel

    real (kind=Rkind) :: Rt1,Rt2
    integer           :: i

    Rt1 = TEN/QModel%B
    DO i=1,100
      !Rt2 = (QModel%B*QModel%A/(SIX*QModel%C)*exp(-QModel%B*Rt1))**(-ONE/SEVEN)
      Rt2 = -ONE/QModel%B* log(SIX*QModel%C/(QModel%B*QModel%A)*Rt1**(-7))
      !write(out_unit,*) i,RT2
      IF (abs(Rt1-Rt2) < ONETENTH**10) EXIT
      Rt1 = Rt2
    END DO
    R0 = Rt1

  END SUBROUTINE get_Q0_QML_Buck

!> @brief Subroutine wich calculates the Buckingham potential with derivatives.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnS_t):         Matrix of dnS with the potential (pot),.
!! @param dnR                TYPE (dnS_t):         derived type wich contain the value for which the potential is calculated: dnR%d0
!! @param QModel          TYPE(QML_Buck_t):    derived type with the Buckingham parameters.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Buck(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_Buck_t),   intent(in)     :: QModel
    TYPE(dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         intent(in)     :: dnQ(:)
    integer,             intent(in)     :: nderiv

    !write(out_unit,*) 'BEGINNING in EvalPot_QML_Buck'

    Mat_OF_PotDia(1,1) = QML_dnBuck(dnQ(1),QModel)

    !write(out_unit,*) 'END in EvalPot_QML_Buck'
    !flush(out_unit)

  END SUBROUTINE EvalPot_QML_Buck

!> @brief Function wich calculates the Buckingham potential with derivatives.
!> @brief V(R) = A.Exp(-B*r)-C/r^6
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return QML_dnBuck            TYPE (dnS_t):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE (dnS_t):           derived type with the value of "r" and,if required, its derivatives.
!! @param BuckPot          TYPE(QML_Buck_t):    derived type with the Buckingham parameters.
  FUNCTION QML_dnBuck(dnR,BuckPot)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t)                          :: QML_dnBuck
    TYPE (dnS_t),          intent(in)     :: dnR

    TYPE (QML_Buck_t),    intent(in)     :: BuckPot


    !write(out_unit,*) 'BEGINNING in QML_dnBuck'
    !write(out_unit,*) 'dnR'
    !CALL Write_dnS(dnR)

    QML_dnBuck = BuckPot%A * exp(-BuckPot%B*dnR) - BuckPot%C * dnR**(-6)

    !write(out_unit,*) 'Buckingham at',get_d0_FROM_dnS(dnR)
    !CALL Write_dnS(QML_dnBuck)
    !write(out_unit,*) 'END in QML_dnBuck'
    !flush(out_unit)

  END FUNCTION QML_dnBuck

  SUBROUTINE RefValues_QML_Buck(QModel,err,nderiv,Q0,dnMatV,d0GGdef,option)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_Buck_t), intent(in)              :: QModel
    integer,           intent(inout)           :: err
    integer,           intent(in)              :: nderiv

    real (kind=Rkind), intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),    intent(inout), optional :: dnMatV
    real (kind=Rkind), intent(inout), optional :: d0GGdef(:,:)
    integer,           intent(in),    optional :: option

    real (kind=Rkind), allocatable :: d0(:,:),d1(:,:,:),d2(:,:,:,:),d3(:,:,:,:,:),V(:)

    integer        :: i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_Buck'
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
      Q0(:) = 7.000000000000000_Rkind
    END IF

    IF (present(dnMatV)) THEN
      err = 0
      CALL alloc_dnMat(dnMatV,nsurf=QModel%nsurf,nVar=QModel%ndim,nderiv=nderiv)

      IF (nderiv >= 0) THEN ! no derivative
        V  = [-4.0943819112303800E-004_Rkind]
        d0 = reshape(V,shape=[QModel%nsurf,QModel%nsurf])
      END IF

      IF (nderiv >= 1) THEN ! 1st order derivatives
        V  = [-1.8553806270932315E-004_Rkind]
        d1 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim])
      END IF

      IF (nderiv >= 2) THEN ! 2d order derivatives
        V  = [1.0880517621509414E-003_Rkind]
        d2 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim,QModel%ndim])
      END IF
      IF (allocated(V)) deallocate(V)

      SELECT CASE (nderiv)
      CASE(0)
        CALL set_dnMat(dnMatV,d0=d0)
        deallocate(d0)
      CASE(1)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1)
        deallocate(d0)
        deallocate(d1)
      CASE(2)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1,d2=d2)
        deallocate(d0)
        deallocate(d1)
        deallocate(d2)
      CASE Default
        STOP 'ERROR in RefValues_QML_HONO0: nderiv MUST < 3'
      END SELECT
    END IF

    IF (present(d0GGdef)) d0GGdef = 2.7454814573212161E-005_Rkind


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_Buck
END MODULE QML_Buck_m
