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

!> @brief Module which makes the initialization, calculation of the Poly1D potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!> @author Félix MOUHAT
!! @date 10/07/2019
!!
MODULE QML_Poly1D_m
  USE QML_Empty_m
  USE QMLLib_NumParameters_m
  IMPLICIT NONE

  PRIVATE

  integer, parameter :: max_order = 100

!> @brief Derived type in which the Poly1D parameters are set-up.
!> @brief Poly1D(R) = Sum_i C_i * r-Req)**i
!> @brief Default parameters for H-F
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param norder      integer: order of the expansion (terms from 0 to norder)
!! @param coef(:)        real: table of coefficient
!! @param req            real: Equilibrium distance (in bohr)
!! @param mu             real: Reduced mass of HF (in au)
  TYPE, EXTENDS (QML_Empty_t) :: QML_Poly1D_t ! V(R) = Sum_i C_i * r-Req)**i
     PRIVATE
     integer           :: norder = 2
     real (kind=Rkind), allocatable :: coef(:)
     real (kind=Rkind) :: req = 1.7329_Rkind !< Equilibrium HF distance (in bohr)
     real (kind=Rkind), PUBLIC :: mu  = 1744.60504565084306291455_Rkind !< Reduced mass of HF (in au)
  CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_Poly1D
    PROCEDURE :: Write_QModel    => Write_QML_Poly1D
    PROCEDURE :: Write0_QModel   => Write0_QML_Poly1D
  END TYPE QML_Poly1D_t

  PUBLIC :: QML_Poly1D_t,Init_QML_Poly1D,Write_QML_Poly1D

CONTAINS
!> @brief Subroutine which makes the initialization of the Poly1D parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Poly1DPot         TYPE(QML_Poly1D_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param D                  real (optional):    Dissociation energy (in Hartree)
!! @param a                  real (optional):    Scaling parameter (in bohr^-1)
!! @param req                real (optional):    Equilibrium distance (in bohr)
  FUNCTION Init_QML_Poly1D(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_Poly1D_t)                          :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    real (kind=Rkind), parameter :: D   = 0.225_Rkind  !< Dissociation energy for HF (in Hartree)
    real (kind=Rkind), parameter :: a   = 1.1741_Rkind !< Scaling parameter for HF (in bohr^-1)
    real (kind=Rkind), parameter :: k = TWO * a**2 * D
    !real (kind=Rkind), parameter :: k = 1.16533_Rkind

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Poly1D'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) '  read_param:     ',read_param
      write(out_unitp,*) '  nio_param_file: ',nio_param_file
      flush(out_unitp)
    END IF


    !QModel_loc%QML_Empty_t = Init_QML_Empty(QModel_in) ! it does not work with nagfor
    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%pot_name = 'Poly1D'
    QModel%ndim     = 1
    QModel%nsurf    = 1


    IF (debug) write(out_unitp,*) 'init default Poly1D parameters'
    QModel%norder = 2
    allocate(QModel%coef(0:QModel%norder))
    QModel%coef   = [ZERO,ZERO,k/two]
    QModel%req    = 1.7329_Rkind

    IF (read_param) THEN
      CALL Read_QML_Poly1D(QModel,nio_param_file)
    END IF

    IF (debug) write(out_unitp,*) 'init Q0 of Poly1D'
    QModel%Q0 = [QModel%req]

    IF (debug) write(out_unitp,*) 'init d0GGdef of Poly1D'
    QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_Poly1D

!> @brief Subroutine wich reads the Poly1D parameters with a namelist.
!!   This can be called only from the "Init_QML_Poly1D" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_Poly1D_t):   derived type in which the parameters are set-up.
!! @param nio            integer           :   file unit to read the parameters.
  SUBROUTINE Read_QML_Poly1D(QModel,nio)
  IMPLICIT NONE

    TYPE (QML_Poly1D_t), intent(inout) :: QModel
    integer,             intent(in)    :: nio

    !local variable
    real (kind=Rkind),allocatable :: coef(:)
    integer                       :: err_read
    real (kind=Rkind)             :: req
    integer                       :: n

    namelist /Poly1D/ coef,req

    allocate(coef(0:max_order))
    coef = huge(one)
    coef(0:Qmodel%norder)    = Qmodel%coef
    req     = QModel%req
    write(out_unitp,*) 'read Poly1D namelist ...'
    flush(out_unitp)

    read(nio,nml=Poly1D,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_QML_Poly1D'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Poly1D" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_QML_Poly1D'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_QML_Poly1D'
      write(out_unitp,*) ' Some parameter names of the namelist "Poly1D" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Poly1D)
      STOP ' ERROR in Read_QML_Poly1D'
    END IF
    !write(out_unitp,nml=Poly1D)

    n = count(coef /= huge(one))

    deallocate(Qmodel%coef)
    allocate(Qmodel%coef(0:n-1))

    Qmodel%coef(0:n-1) = coef(0:n-1)
    Qmodel%norder      = n-1
    QModel%req         = req


  END SUBROUTINE Read_QML_Poly1D
!> @brief Subroutine wich prints the Poly1D parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer          :   file unit to print the parameters.
  SUBROUTINE Write0_QML_Poly1D(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Poly1D_t),    intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'Poly1D parameters:'
    write(nio,*)
    write(nio,*) ' For H-F (Default values, from Morse parameters):'
    write(nio,*) '    V(R) = 0.620329864 * (R-req)^2'
    write(nio,*) '  req = 1.7329 bohr'
    write(nio,*) '  mu  = 1744.60504565084306291455 au'
    write(nio,*)
    write(nio,*) 'Value at: R=req Bohr'
    write(nio,*) 'V        = 0.0 Hartree'
    write(nio,*) 'gradient = 0.0'
    write(nio,*) 'hessian  = 0.620329864'
    write(nio,*)
    write(nio,*) 'end Poly1D parameters'

  END SUBROUTINE Write0_QML_Poly1D
!> @brief Subroutine wich prints the Poly1D current parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_Poly1D_t):   derived type with the Poly1D parameters.
!! @param nio                integer          :   file unit to print the parameters.
  SUBROUTINE Write_QML_Poly1D(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Poly1D_t),    intent(in) :: QModel
    integer,                intent(in) :: nio


    write(nio,*) 'Poly1D current parameters:'
    CALL Write_QML_Empty(QModel%QML_Empty_t,nio)
    write(nio,*)
    write(nio,*) '    V(R) = Sum_i C_i * (R-Req)^i'
    write(nio,*) '  norder: ',QModel%norder
    write(nio,*) '  coef:   ',QModel%coef(:)
    write(nio,*) '  req:    ',QModel%req
    write(nio,*)
    write(nio,*) 'end Poly1D current parameters'

  END SUBROUTINE Write_QML_Poly1D

!> @brief Subroutine wich calculates the Poly1D potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel         TYPE(QML_Poly1D_t):   derived type with the Poly1D parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Poly1D(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS(QML_Poly1D_t), intent(in)     :: QModel
    TYPE(dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         intent(in)     :: dnQ(:)
    integer,             intent(in)     :: nderiv

    integer          :: i
    TYPE (dnS_t)     :: dnDR

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_Poly1D'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) ' nderiv:',nderiv
      write(out_unitp,*) ' Q(:):',(QML_get_d0_FROM_dnS(dnQ(i)),i=1,size(dnQ))
    END IF

    dnDR  = dnQ(1)-QModel%req


    Mat_OF_PotDia(1,1) = QModel%coef(QModel%norder)
    DO i=QModel%norder-1,0,-1
      Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1)*dnDR + QModel%coef(i)
    END DO

    IF (debug) THEN
      write(out_unitp,*) 'Mat_OF_PotDia'
      CALL QML_Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unitp,*)
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF
  END SUBROUTINE EvalPot_QML_Poly1D

END MODULE QML_Poly1D_m
