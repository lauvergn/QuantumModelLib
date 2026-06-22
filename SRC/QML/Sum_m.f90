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

!> @brief Module which makes the initialization, calculation of the Sum of potentials (the model has to be read).
!!
!> @author David Lauvergnat
!! @date 15/06/2026
!!
MODULE QML_Sum_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Sum parameters are set-up.
!> @brief Sum(R) = Sum_i C_i * r-Req)**i
!> @brief Default parameters for H-F
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param norder      integer: order of the expansion (terms from 0 to norder)
!! @param coef(:)        real: table of coefficient
!! @param req            real: Equilibrium distance (in bohr)
!! @param mu             real: Reduced mass of HF (in au)
  TYPE, EXTENDS (QML_Empty_t) :: QML_Sum_t ! V(R) = Sum_i C_i * r-Req)**i
    PRIVATE
    integer                   :: nb_model = -1
    TYPE (QML_t), allocatable :: tab_Op(:) ! Table containing the models
  CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_Sum
    PROCEDURE :: Write_QModel     => Write_QML_Sum
    PROCEDURE :: RefValues_QModel => RefValues_QML_Sum
  END TYPE QML_Sum_t

  PUBLIC :: QML_Sum_t,Init_QML_Sum,Write_QML_Sum

CONTAINS
!> @brief Subroutine which makes the initialization of the Sum parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param SumPot         TYPE(QML_Sum_t):   derived type in which the parameters are set-up.
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param D                  real (optional):    Dissociation energy (in Hartree)
!! @param a                  real (optional):    Scaling parameter (in bohr^-1)
!! @param req                real (optional):    Equilibrium distance (in bohr)
  FUNCTION Init_QML_Sum(QModel_in,read_param,nio_param_file) RESULT(QModel)
    IMPLICIT NONE

    TYPE (QML_Sum_t)                          :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    real (kind=Rkind), parameter :: D   = 0.225_Rkind  !< Dissociation energy for HF (in Hartree)
    real (kind=Rkind), parameter :: a   = 1.1741_Rkind !< Scaling parameter for HF (in bohr^-1)
    real (kind=Rkind), parameter :: k = TWO * a**2 * D
    !real (kind=Rkind), parameter :: k = 1.16533_Rkind

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Sum'
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

    QModel%pot_name = 'Sum'
    QModel%ndim     = 1
    QModel%nsurf    = 1


    IF (debug) write(out_unit,*) 'init default Sum parameters'

    CALL Read_QML_Sum(QModel,nio_param_file)

    IF (debug) write(out_unit,*) 'init Q0 of Sum'
    !QModel%Q0 = [QModel%req]

    IF (debug) write(out_unit,*) 'init d0GGdef of Sum'
    !QModel%d0GGdef = reshape([ONE/QModel%mu],shape=[1,1])

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Sum

!> @brief Subroutine wich reads the Sum parameters with a namelist.
!!   This can be called only from the "Init_QML_Sum" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_Sum_t):   derived type in which the parameters are set-up.
!! @param nio            integer           :   file unit to read the parameters.
  SUBROUTINE Read_QML_Sum(QModel,nio)
    IMPLICIT NONE

    TYPE (QML_Sum_t), intent(inout) :: QModel
    integer,          intent(in)    :: nio

    !local variable
    integer                       :: err_read
    integer                       :: nb_model

    namelist /Sum/ nb_model

    nb_model     = 0
    write(out_unit,*) 'read Sum namelist ...'
    flush(out_unit)

    read(nio,nml=Sum,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Sum'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "Sum" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Sum'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Sum'
      write(out_unit,*) ' Some parameter names of the namelist "Sum" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=Sum)
      STOP ' ERROR in Read_QML_Sum'
    END IF
    !write(out_unit,nml=Sum)

    Qmodel%nb_model      = nb_model


  END SUBROUTINE Read_QML_Sum
!! === README ==
!! Polynomial potential: $V(R) = \sum_i coef(i) \cdot (r-Req)^i$
!! pot_name  = 'Sum'
!! ndim      = 1
!! nsurf     = 1
!! reduced mass      = 1744.60504565084306291455 au
!! remark: Default parameters for H-F
!! === END README ==
!> @brief Subroutine wich prints the Sum current parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel         TYPE(QML_Sum_t):   derived type with the Sum parameters.
!! @param nio            integer          :   file unit to print the parameters.
  SUBROUTINE Write_QML_Sum(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_Sum_t),    intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'Sum default parameters:'
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
    write(nio,*) 'end Sum default parameters'
    write(nio,*)
    write(nio,*) 'Sum current parameters:'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) '    V(R) = Sum_i Coef_i * (R-Req)^i'
    write(nio,*) '  norder: ',QModel%norder
    write(nio,*) '  Coef:   ',QModel%coef(:)
    write(nio,*) '  Req:    ',QModel%req
    write(nio,*)
    write(nio,*) 'end Sum current parameters'

  END SUBROUTINE Write_QML_Sum

!> @brief Subroutine wich calculates the Sum potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      Potential with derivatives,.
!! @param r                  real:                value for which the potential is calculated
!! @param QModel         TYPE(QML_Sum_t):   derived type with the Sum parameters.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Sum(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_Sum_t), intent(in)     :: QModel
    TYPE(dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnS_t),         intent(in)     :: dnQ(:)
    integer,             intent(in)     :: nderiv

    integer          :: i
    TYPE (dnS_t)     :: dnDR

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_Sum'
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
  END SUBROUTINE EvalPot_QML_Sum



   SUBROUTINE RefValues_QML_Sum(QModel,err,nderiv,Q0,dnMatV,d0GGdef,option)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_Sum_t), intent(in)              :: QModel
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
    character (len=*), parameter :: name_sub='RefValues_QML_Sum'
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
      Q0(:) = [1.7329000000000001_Rkind]
    END IF

    IF (present(dnMatV)) THEN
      err = 0

      IF (nderiv >= 0) THEN ! no derivative
        V  = [0.000000000000000_Rkind]
        d0 = reshape(V,shape=[QModel%nsurf,QModel%nsurf])
      END IF
 
      IF (nderiv >= 1) THEN ! 1st order derivatives
        V  = [0.000000000000000_Rkind]
        d1 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim])
      END IF

      IF (nderiv >= 2) THEN ! 2d order derivatives
        V  = [0.62032986449999994_Rkind]
        d2 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim,QModel%ndim])
      END IF
      SELECT CASE (nderiv)
      CASE(0)
        CALL set_dnMat(dnMatV,d0=d0)
      CASE(1)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1)
      CASE(2)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1,d2=d2)
      CASE Default
        STOP 'ERROR in RefValues_QML_Sum: nderiv MUST < 3'
      END SELECT

    END IF

    IF (present(d0GGdef)) THEN 
      d0GGdef = 5.7319563673905317E-004_Rkind * Identity_Mat(QModel%ndim)
    END IF


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_Sum
END MODULE QML_Sum_m
