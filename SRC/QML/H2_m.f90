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

!> @brief Module which makes the initialization, calculation of the H2 potential (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 18/10/2021
!!
MODULE QML_H2_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the H2 parameters are set-up.
!> @brief Default parameters for H-H
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param mu             real: Reduced mass of HH (in au)
  TYPE, EXTENDS (QML_Empty_t) :: QML_H2_t
     PRIVATE
     real (kind=Rkind), PUBLIC      :: mu  = 1837.1526464003414_Rkind/TWO !< Reduced mass of H2 (in au)
     real (kind=Rkind)              :: Req = 1.4_Rkind
     real (kind=Rkind)              :: E0  = -1.17404460_Rkind ! CCSD(T)-F12B/VTZ-F12
     real (kind=Rkind), allocatable :: TaylorPot(:)

  CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_H2
    PROCEDURE :: Write_QModel    => Write_QML_H2
    PROCEDURE :: Write0_QModel   => Write0_QML_H2
  END TYPE QML_H2_t

  PUBLIC :: QML_H2_t,Init_QML_H2,Write_QML_H2,QML_dnH2

CONTAINS
!> @brief Subroutine which makes the initialization of the H2 parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param H2Pot              TYPE(QML_H2_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_H2(QModel_in,read_param,nio_param_file) RESULT(QModel)
    IMPLICIT NONE

    TYPE (QML_H2_t)                              :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H2'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) '  read_param:     ',read_param
      write(out_unit,*) '  nio_param_file: ',nio_param_file
      flush(out_unit)
    END IF


    !QModel_loc%QML_Empty_t = Init_QML_Empty(QModel_in) ! it does not work with nagfor
    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%pot_name = 'h2'
    QModel%ndim     = 1
    QModel%nsurf    = 1


    IF (debug) write(out_unit,*) 'init default H2 parameters (HF)'
    QModel%E0  = -1.17404460_Rkind
    QModel%Req = 1.4_Rkind

    IF (read_param) THEN
      CALL Read_QML_H2(QModel%Req,QModel%E0,nio_param_file)
    ELSE
      IF (debug) write(out_unit,*) 'init H2 parameters, if present'
    END IF

    SELECT CASE (QModel%option)
    CASE (1)
      write(out_unit,*) 'Fourth order expansion in (R-Req)'

      QModel%TaylorPot = [-0.000498716666667_Rkind,0.185668083333_Rkind,-0.219403333333_Rkind,0.182041666667_Rkind]

      IF (debug) write(out_unit,*) 'init Q0 of H2'
      QModel%Q0 = [QModel%req]

      IF (debug) write(out_unit,*) 'init d0GGdef of H2'
      QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])
    CASE(2) ! The coordinate x is defined as x = Req/R
      write(out_unit,*) 'Fourth order expansion in (x-1) with x = Req/R'

      ! with DQ=+/-0.1 and +/- 0.2
      !QModel%TaylorPot = [-0.000753325_Rkind,0.363249875_Rkind,-0.1412875_Rkind,0.0165625_Rkind]

      ! with DQ=+/-0.05 and +/- 0.1
      QModel%TaylorPot = [0.00077465_Rkind,0.363216166667_Rkind,-0.14342_Rkind,0.019933333338_Rkind]

      QModel%Q0 = [ONE]
      ! 1/mu d2/dR^2  =>   1/mu [x" d/dx + x'^2 d2/dx^2]
      !   => G(x) = 1/mu * x'^2 = 1/mu * (-Req/R^2)^2 = 1/mu * (x^2/Req)^2
      !  and xeq = 1.
      QModel%d0GGdef = reshape([ONE/QModel%req**2 * ONE/QModel%mu],shape=[1,1])
    CASE default
      write(out_unit,*) 'ERROR in Init_QML_H2'
      write(out_unit,*) 'Possible option: 1,2'
      write(out_unit,*) 'Actual value, option=',QModel%option
      STOP 'ERROR in Init_QML_H2: wrong option value'
    END SELECT

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_H2

!> @brief Subroutine wich reads the H2 parameters with a namelist.
!!   This can be called only from the "Init_QML_H2" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_H2_t):   derived type in which the parameters are set-up.
!! @param nio                integer          :   file unit to read the parameters.
  SUBROUTINE Read_QML_H2(Req_inout,E0_inout,nio)
    IMPLICIT NONE

    integer,           intent(in)    :: nio
    real (kind=Rkind), intent(inout) :: Req_inout,E0_inout

    !local variable
    real (kind=Rkind) :: Req,E0
    integer           :: err_read

    namelist /H2/ Req,E0

    E0  = E0_inout
    Req = Req_inout
    write(out_unit,*) 'read H2 namelist ...'

    read(nio,nml=H2,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_H2'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "H2" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_H2'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_H2'
      write(out_unit,*) ' Some parameter names of the namelist "H2" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=H2)
      STOP ' ERROR in Read_QML_H2'
    END IF
    !write(out_unit,nml=H2)

    E0_inout   = E0
    Req_inout  = Req

  END SUBROUTINE Read_QML_H2
!> @brief Subroutine wich prints the H2 parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer          :   file unit to print the parameters.
  SUBROUTINE Write0_QML_H2(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_H2_t),    intent(in) :: QModel
    integer,            intent(in) :: nio

    write(nio,*) 'H2 parameters:'
    write(nio,*)
    write(nio,*) ' For H-H (Default values):'
    write(nio,*) '    V(R) = E0 + Sum_i a_i(R-Req)^i'
    write(nio,*) '  Req = 1.4 bohr'
    write(nio,*) '  E0  = -1.17404460 Hartree'

    write(nio,*) '  mu  = 1744.60504565084306291455 au'
    write(nio,*)
    write(nio,*) 'Value at: R=2.0 Bohr'
    write(nio,*) 'V        = 0.016304 Hartree'
    write(nio,*) 'gradient = 0.103940'
    write(nio,*) 'hessian  = 0.209272'
    write(nio,*)
    write(nio,*) 'end H2 parameters'

  END SUBROUTINE Write0_QML_H2
!> @brief Subroutine wich prints the H2 current parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_H2_t):   derived type with the H2 parameters.
!! @param nio            integer          :   file unit to print the parameters.
  SUBROUTINE Write_QML_H2(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_H2_t),    intent(in) :: QModel
    integer,            intent(in) :: nio


    write(nio,*) 'H2 current parameters:'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) '    V(R) = E0 + Sum_i a_i(R-Req)^i'
    write(nio,*) '  E0:            ',QModel%E0
    write(nio,*) '  Req:           ',QModel%Req
    write(nio,*) '  TaylorPot(:)   ',QModel%TaylorPot(:)
    write(nio,*)
    write(nio,*) 'end H2 current parameters'

  END SUBROUTINE Write_QML_H2

!> @brief Subroutine wich calculates the H2 potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_H2_t):   derived type with the H2 parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H2(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2_t), intent(in)     :: QModel
    TYPE(dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         intent(in)     :: dnQ(:)
    integer,             intent(in)     :: nderiv

    integer :: i
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_H2'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) ' nderiv:',nderiv
      write(out_unit,*) ' Q(:):',(get_d0(dnQ(i)),i=1,size(dnQ))
    END IF

    Mat_OF_PotDia(1,1) = QML_dnH2(dnQ(1),QModel)

    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF
  END SUBROUTINE EvalPot_QML_H2

!> @brief Function wich calculates the H2 potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @return QML_dnH2      TYPE (dnS_t):     derived type with a value (pot),,if required, its derivatives (gradient (grad) and hessian (hess)).
!! @param dnQ            TYPE (dnS_t):     derived type with the value of "r" and,if required, its derivatives.
!! @param QModel         TYPE(QML_H2_t):   derived type with the H2 parameters.
  FUNCTION QML_dnH2(dnQ,QModel)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t)                     :: QML_dnH2

    TYPE (dnS_t),      intent(in)    :: dnQ
    TYPE (QML_H2_t),   intent(in)    :: QModel

    !local variable
    TYPE (dnS_t)        :: dnDeltaQ
    integer             :: i
    !logical, parameter  :: debug = .TRUE.
    logical, parameter  :: debug = .FALSE.

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING in QML_dnH2'
      CALL Write_dnS(dnQ,info='dnQ')
    END IF

    SELECT CASE (QModel%option)
    CASE (1)
      dnDeltaQ  = dnQ-QModel%Req
    CASE(2)
      dnDeltaQ  = dnQ-ONE
    END SELECT
    IF (debug) CALL Write_dnS(dnDeltaQ,info='dnDeltaQ')

    QML_dnH2 = QModel%E0
    DO i=1,size(QModel%TaylorPot)
      QML_dnH2 = QML_dnH2 + QModel%TaylorPot(i)*dnDeltaQ**i
    END DO

    CALL dealloc_dnS(dnDeltaQ)

    IF (debug) THEN
      write(out_unit,*) 'H2 at',get_d0(dnQ)
      CALL Write_dnS(QML_dnH2)
      write(out_unit,*) 'END in QML_dnH2'
      flush(out_unit)
    END IF

  END FUNCTION QML_dnH2

END MODULE QML_H2_m
