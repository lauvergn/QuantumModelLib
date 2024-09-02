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

!> @brief Module which makes the initialization, calculation of the Poly1D potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 30/11/2021
!!
MODULE QML_Poly1D_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
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
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) '  read_param:     ',read_param
      write(out_unit,*) '  nio_param_file: ',nio_param_file
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%pot_name = 'Poly1D'
    QModel%ndim     = 1
    QModel%nsurf    = 1


    IF (debug) write(out_unit,*) 'init default Poly1D parameters'
    QModel%norder = 2
    allocate(QModel%coef(0:QModel%norder))
    QModel%coef   = [ZERO,ZERO,k/two]
    QModel%req    = 1.7329_Rkind

    IF (read_param) THEN
      CALL Read_QML_Poly1D(QModel,nio_param_file)
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of Poly1D'
    QModel%Q0 = [QModel%req]

    IF (debug) write(out_unit,*) 'init d0GGdef of Poly1D'
    QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
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
    write(out_unit,*) 'read Poly1D namelist ...'
    flush(out_unit)

    read(nio,nml=Poly1D,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Poly1D'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "Poly1D" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Poly1D'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Poly1D'
      write(out_unit,*) ' Some parameter names of the namelist "Poly1D" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=Poly1D)
      STOP ' ERROR in Read_QML_Poly1D'
    END IF
    !write(out_unit,nml=Poly1D)

    n = count(coef /= huge(one))

    deallocate(Qmodel%coef)
    allocate(Qmodel%coef(0:n-1))

    Qmodel%coef(0:n-1) = coef(0:n-1)
    Qmodel%norder      = n-1
    QModel%req         = req


  END SUBROUTINE Read_QML_Poly1D
!! === README ==
!! Polynomial potential: $V(R) = \sum_i coef(i) \cdot (r-Req)^i$
!! pot_name  = 'Poly1D'
!! ndim      = 1
!! nsurf     = 1
!! reduced mass      = 1744.60504565084306291455 au
!! remark: Default parameters for H-F
!! === END README ==
!> @brief Subroutine wich prints the Poly1D current parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_Poly1D_t):   derived type with the Poly1D parameters.
!! @param nio            integer          :   file unit to print the parameters.
  SUBROUTINE Write_QML_Poly1D(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_Poly1D_t),    intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'Poly1D default parameters:'
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
    write(nio,*) 'end Poly1D default parameters'
    write(nio,*)
    write(nio,*) 'Poly1D current parameters:'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) '    V(R) = Sum_i Coef_i * (R-Req)^i'
    write(nio,*) '  norder: ',QModel%norder
    write(nio,*) '  Coef:   ',QModel%coef(:)
    write(nio,*) '  Req:    ',QModel%req
    write(nio,*)
    write(nio,*) 'end Poly1D current parameters'

  END SUBROUTINE Write_QML_Poly1D

!> @brief Subroutine wich calculates the Poly1D potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      Potential with derivatives,.
!! @param r                  real:                value for which the potential is calculated
!! @param QModel         TYPE(QML_Poly1D_t):   derived type with the Poly1D parameters.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Poly1D(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
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
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) ' nderiv:',nderiv
      write(out_unit,*) ' Q(:):',(get_d0(dnQ(i)),i=1,size(dnQ))
    END IF

    dnDR  = dnQ(1)-QModel%req


    Mat_OF_PotDia(1,1) = QModel%coef(QModel%norder)
    DO i=QModel%norder-1,0,-1
      Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1)*dnDR + QModel%coef(i)
    END DO

    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
  END SUBROUTINE EvalPot_QML_Poly1D

END MODULE QML_Poly1D_m
