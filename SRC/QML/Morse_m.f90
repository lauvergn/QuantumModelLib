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

!> @brief Module which makes the initialization, calculation of the Morse potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!> @author Félix MOUHAT
!! @date 10/07/2019
!!
MODULE QML_Morse_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Morse parameters are set-up.
!> @brief morse(R) = D*(1-exp(-a*(r-Req))**2
!> @brief Default parameters for H-F
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param D              real: Dissociation energy (in Hartree)
!! @param a              real: Scaling parameter (in bohr^-1)
!! @param req            real: Equilibrium distance (in bohr)
!! @param mu             real: Reduced mass of HF (in au)
  TYPE, EXTENDS (QML_Empty_t) :: QML_Morse_t ! V(R) = D*(1-exp(-a*(r-Req))**2
     PRIVATE
     real (kind=Rkind) :: D   = 0.225_Rkind  !< Dissociation energy for HF (in Hartree)
     real (kind=Rkind) :: a   = 1.1741_Rkind !< Scaling parameter for HF (in bohr^-1)
     real (kind=Rkind) :: req = 1.7329_Rkind !< Equilibrium HF distance (in bohr)
     real (kind=Rkind), PUBLIC :: mu  = 1744.60504565084306291455_Rkind !< Reduced mass of HF (in au)
  CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_Morse
    PROCEDURE :: Write_QModel    => Write_QML_Morse
    PROCEDURE :: Write0_QModel   => Write0_QML_Morse
  END TYPE QML_Morse_t

  PUBLIC :: QML_Morse_t,Init_QML_Morse,Init0_QML_Morse,Write_QML_Morse,QML_dnMorse

CONTAINS
!> @brief Subroutine which makes the initialization of the Morse parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param MorsePot         TYPE(QML_Morse_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param D                  real (optional):    Dissociation energy (in Hartree)
!! @param a                  real (optional):    Scaling parameter (in bohr^-1)
!! @param req                real (optional):    Equilibrium distance (in bohr)
  FUNCTION Init_QML_Morse(QModel_in,read_param,nio_param_file,D,a,req) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_Morse_t)                          :: QModel

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: D,a,req

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Morse'
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

    QModel%pot_name = 'morse'
    QModel%ndim     = 1
    QModel%nsurf    = 1


    IF (debug) write(out_unit,*) 'init default morse parameters (HF)'
    ! the initialization with QML_Morse_t does not work. Why ???
    !QModel = QML_Morse_t(D=0.225_Rkind,a=1.1741_Rkind,Req=1.7329_Rkind) ! default values (HF)
    CALL Init0_QML_Morse(QModel,D=0.225_Rkind,a=1.1741_Rkind,req=1.7329_Rkind)

    IF (read_param) THEN
      CALL Read_QML_Morse(QModel%D,QModel%a,QModel%req,nio_param_file)
    ELSE
      IF (debug) write(out_unit,*) 'init morse parameters (D,a,req), if present'
      IF (present(D))   QModel%D   = D
      IF (present(a))   QModel%a   = a
      IF (present(req)) QModel%req = req
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of morse'
    QModel%Q0 = [QModel%req]

    IF (debug) write(out_unit,*) 'init d0GGdef of morse'
    QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Morse

  SUBROUTINE Init0_QML_Morse(QModel,D,a,req,model_name)
  IMPLICIT NONE

    TYPE (QML_Morse_t),           intent(inout)   :: QModel
    real (kind=Rkind), optional,   intent(in)      :: D,a,req
    character (len=*), optional,   intent(in)      :: model_name

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init0_QML_Morse'
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
    QModel%pot_name   = 'morse'

    IF (present(model_name)) QModel%pot_name = model_name

    IF (debug) write(out_unit,*) 'init morse parameters (D,a,req), if present'

    IF (present(D))   QModel%D   = D
    IF (present(a))   QModel%a   = a
    IF (present(req)) QModel%req = req

    IF (debug) THEN
      write(out_unit,*) 'D,a,req: ',QModel%D,QModel%a,QModel%req
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Init0_QML_Morse

!> @brief Subroutine wich reads the Morse parameters with a namelist.
!!   This can be called only from the "Init_QML_Morse" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_Morse_t):   derived type in which the parameters are set-up.
!! @param nio                integer          :   file unit to read the parameters.
  SUBROUTINE Read_QML_Morse(D_inout,a_inout,req_inout,nio)
    IMPLICIT NONE

    integer,           intent(in)    :: nio
    real (kind=Rkind), intent(inout) :: D_inout,a_inout,req_inout

    !local variable
    real (kind=Rkind) :: D,a,req
    integer           :: err_read

    namelist /Morse/ D,a,req

    D   = D_inout
    a   = a_inout
    req = req_inout
    write(out_unit,*) 'read Morse namelist ...'

    read(nio,nml=Morse,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Morse'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "Morse" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Morse'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Morse'
      write(out_unit,*) ' Some parameter names of the namelist "Morse" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=Morse)
      STOP ' ERROR in Read_QML_Morse'
    END IF
    !write(out_unit,nml=Morse)

    D_inout    = D
    a_inout    = a
    req_inout  = req

  END SUBROUTINE Read_QML_Morse
!> @brief Subroutine wich prints the Morse parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer          :   file unit to print the parameters.
  SUBROUTINE Write0_QML_Morse(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_Morse_t),    intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'Morse parameters:'
    write(nio,*)
    write(nio,*) ' For H-F (Default values):'
    write(nio,*) '    V(R) = D.( 1 - exp(-a.(r-req)) )^2'
    write(nio,*) '  D   = 0.225  Hartree'
    write(nio,*) '  a   = 1.1741 bohr^-1'
    write(nio,*) '  req = 1.7329 bohr'
    write(nio,*) '  mu  = 1744.60504565084306291455 au'
    write(nio,*)
    write(nio,*) 'Value at: R=2.0 Bohr'
    write(nio,*) 'V        = 0.016304 Hartree'
    write(nio,*) 'gradient = 0.103940'
    write(nio,*) 'hessian  = 0.209272'
    write(nio,*)
    write(nio,*) 'end Morse parameters'

  END SUBROUTINE Write0_QML_Morse
!> @brief Subroutine wich prints the Morse current parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_Morse_t):   derived type with the Morse parameters.
!! @param nio                integer          :   file unit to print the parameters.
  SUBROUTINE Write_QML_Morse(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_Morse_t),    intent(in) :: QModel
    integer,               intent(in) :: nio


    write(nio,*) 'Morse current parameters:'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) '    V(R) = D.( 1 - exp(-a.(r-req)) )^2'
    write(nio,*) '  D:   ',QModel%D
    write(nio,*) '  a:   ',QModel%a
    write(nio,*) '  req: ',QModel%req
    write(nio,*)
    write(nio,*) 'end Morse current parameters'

  END SUBROUTINE Write_QML_Morse

!> @brief Subroutine wich calculates the Morse potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel         TYPE(QML_Morse_t):   derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Morse(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m, ONLY :  dnS_t,Write_dnS,get_d0
    IMPLICIT NONE

    CLASS(QML_Morse_t), intent(in)      :: QModel
    TYPE(dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         intent(in)     :: dnQ(:)
    integer,             intent(in)     :: nderiv

    integer :: i
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_Morse'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) ' nderiv:',nderiv
      write(out_unit,*) ' Q(:):',(get_d0(dnQ(i)),i=1,size(dnQ))
    END IF

    Mat_OF_PotDia(1,1) = QML_dnMorse(dnQ(1),QModel)

    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
  END SUBROUTINE EvalPot_QML_Morse

!> @brief Function wich calculates the Morse potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return QML_dnMorse           TYPE (dnS_t):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE (dnS_t):           derived type with the value of "r" and,if required, its derivatives.
!! @param QModel         TYPE(QML_Morse_t):   derived type with the Morse parameters.
  FUNCTION QML_dnMorse(dnR,QModel)
    USE ADdnSVM_m, ONLY :  dnS_t,exp,Write_dnS,get_d0,dealloc_dnS
    IMPLICIT NONE

    TYPE (dnS_t)                        :: QML_dnMorse

    TYPE (dnS_t),         intent(in)    :: dnR
    TYPE (QML_Morse_t),   intent(in)    :: QModel

    !local variable
    TYPE (dnS_t)     :: dnbeta

    !write(out_unit,*) 'BEGINNING in QML_dnMorse'
    !write(out_unit,*) 'dnR'
    !CALL Write_dnS(dnR)

    dnbeta  = exp(-QModel%a*(dnR-QModel%req))
    !write(out_unit,*) 'dnbeta'
    !CALL Write_dnS(dnbeta)

    QML_dnMorse = QModel%D * (ONE-dnbeta)**2

     CALL dealloc_dnS(dnbeta)

    !write(out_unit,*) 'Morse at',get_d0(dnR)
    !CALL Write_dnS(QML_dnMorse)
    !write(out_unit,*) 'END in QML_dnMorse'
    !flush(out_unit)

  END FUNCTION QML_dnMorse

END MODULE QML_Morse_m
