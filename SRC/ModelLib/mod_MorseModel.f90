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
!    Copyright 2016 David Lauvergnat [1]
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
MODULE mod_MorseModel
  USE mod_EmptyModel
  USE mod_QML_NumParameters
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
  TYPE, EXTENDS (EmptyModel_t) :: MorseModel_t ! V(R) = D*(1-exp(-a*(r-Req))**2
     PRIVATE
     real (kind=Rkind) :: D   = 0.225_Rkind  !< Dissociation energy for HF (in Hartree)
     real (kind=Rkind) :: a   = 1.1741_Rkind !< Scaling parameter for HF (in bohr^-1)
     real (kind=Rkind) :: req = 1.7329_Rkind !< Equilibrium HF distance (in bohr)
     real (kind=Rkind), PUBLIC :: mu  = 1744.60504565084306291455_Rkind !< Reduced mass of HF (in au)
  CONTAINS
    PROCEDURE :: Eval_QModel_Pot => eval_MorsePot
    PROCEDURE :: Write_QModel    => Write_MorseModel
    PROCEDURE :: Write0_QModel   => Write0_MorseModel
  END TYPE MorseModel_t

  PUBLIC :: MorseModel_t,Init_MorseModel,Init0_MorseModel,Write_MorseModel,dnMorse

CONTAINS
!> @brief Subroutine which makes the initialization of the Morse parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param MorsePot         TYPE(MorseModel_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param D                  real (optional):    Dissociation energy (in Hartree)
!! @param a                  real (optional):    Scaling parameter (in bohr^-1)
!! @param req                real (optional):    Equilibrium distance (in bohr)
  FUNCTION Init_MorseModel(QModel_in,read_param,nio_param_file,D,a,req) RESULT(QModel)
  IMPLICIT NONE

    TYPE (MorseModel_t)                          :: QModel

    TYPE(EmptyModel_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: D,a,req

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_MorseModel'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) '  read_param:     ',read_param
      write(out_unitp,*) '  nio_param_file: ',nio_param_file
      flush(out_unitp)
    END IF


    !QModel_loc%EmptyModel_t = Init_EmptyModel(QModel_in) ! it does not work with nagfor
    CALL Init0_EmptyModel(QModel%EmptyModel_t,QModel_in)

    QModel%pot_name = 'morse'
    QModel%ndim     = 1
    QModel%nsurf    = 1


    IF (debug) write(out_unitp,*) 'init default morse parameters (HF)'
    ! the initialization with MorseModel_t does not work. Why ???
    !QModel = MorseModel_t(D=0.225_Rkind,a=1.1741_Rkind,Req=1.7329_Rkind) ! default values (HF)
    CALL Init0_MorseModel(QModel,D=0.225_Rkind,a=1.1741_Rkind,req=1.7329_Rkind)

    IF (read_param) THEN
      CALL Read_MorseModel(QModel%D,QModel%a,QModel%req,nio_param_file)
    ELSE
      IF (debug) write(out_unitp,*) 'init morse parameters (D,a,req), if present'
      IF (present(D))   QModel%D   = D
      IF (present(a))   QModel%a   = a
      IF (present(req)) QModel%req = req
    END IF

    IF (debug) write(out_unitp,*) 'init Q0 of morse'
    QModel%Q0 = [QModel%req]

    IF (debug) write(out_unitp,*) 'init d0GGdef of morse'
    QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_MorseModel

  SUBROUTINE Init0_MorseModel(QModel,D,a,req,model_name)
  IMPLICIT NONE

    TYPE (MorseModel_t),           intent(inout)   :: QModel
    real (kind=Rkind), optional,   intent(in)      :: D,a,req
    character (len=*), optional,   intent(in)      :: model_name

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init0_MorseModel'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    QModel%ndim     = 1
    QModel%nsurf    = 1
    QModel%pot_name = 'morse'
    IF (present(model_name)) QModel%pot_name = model_name

    IF (debug) write(out_unitp,*) 'init morse parameters (D,a,req), if present'

    IF (present(D))   QModel%D   = D
    IF (present(a))   QModel%a   = a
    IF (present(req)) QModel%req = req

    IF (debug) THEN
      write(out_unitp,*) 'D,a,req: ',QModel%D,QModel%a,QModel%req
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Init0_MorseModel

!> @brief Subroutine wich reads the Morse parameters with a namelist.
!!   This can be called only from the "Init_MorseModel" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(MorseModel_t):   derived type in which the parameters are set-up.
!! @param nio                integer          :   file unit to read the parameters.
  SUBROUTINE Read_MorseModel(D_inout,a_inout,req_inout,nio)
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
    write(out_unitp,*) 'read Morse namelist ...'

    read(nio,nml=Morse,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_MorseModel'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Morse" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_MorseModel'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_MorseModel'
      write(out_unitp,*) ' Some parameter names of the namelist "Morse" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Morse)
      STOP ' ERROR in Read_MorseModel'
    END IF
    !write(out_unitp,nml=Morse)

    D_inout    = D
    a_inout    = a
    req_inout  = req

  END SUBROUTINE Read_MorseModel
!> @brief Subroutine wich prints the Morse parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer          :   file unit to print the parameters.
  SUBROUTINE Write0_MorseModel(QModel,nio)
  IMPLICIT NONE

    CLASS(MorseModel_t),    intent(in) :: QModel
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

  END SUBROUTINE Write0_MorseModel
!> @brief Subroutine wich prints the Morse current parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(MorseModel_t):   derived type with the Morse parameters.
!! @param nio                integer          :   file unit to print the parameters.
  SUBROUTINE Write_MorseModel(QModel,nio)
  IMPLICIT NONE

    CLASS(MorseModel_t),    intent(in) :: QModel
    integer,                intent(in) :: nio


    write(nio,*) 'Morse current parameters:'
    CALL Write_EmptyModel(QModel%EmptyModel_t,nio)
    write(nio,*)
    write(nio,*) '    V(R) = D.( 1 - exp(-a.(r-req)) )^2'
    write(nio,*) '  D:   ',QModel%D
    write(nio,*) '  a:   ',QModel%a
    write(nio,*) '  req: ',QModel%req
    write(nio,*)
    write(nio,*) 'end Morse current parameters'

  END SUBROUTINE Write_MorseModel

!> @brief Subroutine wich calculates the Morse potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel         TYPE(MorseModel_t):   derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE Eval_MorsePot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(MorseModel_t), intent(in)     :: QModel
    TYPE(dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         intent(in)     :: dnQ(:)
    integer,             intent(in)     :: nderiv

    integer :: i
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_MorsePot'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) ' nderiv:',nderiv
      write(out_unitp,*) ' Q(:):',(QML_get_d0_FROM_dnS(dnQ(i)),i=1,size(dnQ))
    END IF

    Mat_OF_PotDia(1,1) = dnMorse(dnQ(1),QModel)

    IF (debug) THEN
      write(out_unitp,*) 'Mat_OF_PotDia'
      CALL QML_Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unitp,*)
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF
  END SUBROUTINE Eval_MorsePot

!> @brief Function wich calculates the Morse potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return dnMorse           TYPE (dnS_t):           derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnR                TYPE (dnS_t):           derived type with the value of "r" and,if required, its derivatives.
!! @param QModel         TYPE(MorseModel_t):   derived type with the Morse parameters.
  FUNCTION dnMorse(dnR,QModel)
  USE mod_dnS
  IMPLICIT NONE

    TYPE (dnS_t)                         :: dnMorse

    TYPE (dnS_t),          intent(in)    :: dnR
    TYPE (MorseModel_t),   intent(in)    :: QModel

    !local variable
    TYPE (dnS_t)     :: dnbeta

    !write(out_unitp,*) 'BEGINNING in dnMorse'
    !write(out_unitp,*) 'dnR'
    !CALL QML_Write_dnS(dnR)

    dnbeta  = exp(-QModel%a*(dnR-QModel%req))
    !write(out_unitp,*) 'dnbeta'
    !CALL QML_Write_dnS(dnbeta)

    dnMorse = QModel%D * (ONE-dnbeta)**2

     CALL QML_dealloc_dnS(dnbeta)

    !write(out_unitp,*) 'Morse at',QML_get_d0_FROM_dnS(dnR)
    !CALL QML_Write_dnS(dnMorse)
    !write(out_unitp,*) 'END in dnMorse'
    !flush(out_unitp)

  END FUNCTION dnMorse

END MODULE mod_MorseModel
