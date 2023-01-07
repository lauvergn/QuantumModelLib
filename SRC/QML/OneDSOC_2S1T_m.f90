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
!> @brief Module which makes the initialization, calculation of the OneDSOC_2S1T potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_OneDSOC_2S1T_m
  USE QDUtil_NumParameters_m, out_unitp => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the OneDSOC_2S1T parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_OneDSOC_2S1T_t

   PRIVATE

      real (kind=Rkind) :: a1=0.034520_Rkind
      real (kind=Rkind) :: a2=0.50_Rkind

      real (kind=Rkind) :: alpha1=0.35_Rkind
      real (kind=Rkind) :: alpha2=0.25_Rkind
      real (kind=Rkind) :: gamma_1=ONETENTH**3 * sqrt(TWO)
      real (kind=Rkind) :: gamma_2=ONETENTH**3 * sqrt(TWO)

      real (kind=Rkind) :: DE=0.04_Rkind
      real (kind=Rkind) :: DEs=0.035_Rkind

      real (kind=Rkind) :: Phi_12=pi/FOUR ! only this paramter can be changed

      real (kind=Rkind) :: mu  = 20000._Rkind !< Reduced mass from Granucci et al. paper (in au)


   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_OneDSOC_2S1T
    PROCEDURE :: Write_QModel    => Write_QML_OneDSOC_2S1T
    PROCEDURE :: Write0_QModel   => Write0_QML_OneDSOC_2S1T
  END TYPE QML_OneDSOC_2S1T_t

  PUBLIC :: QML_OneDSOC_2S1T_t,Init_QML_OneDSOC_2S1T

  CONTAINS
!> @brief Function which makes the initialization of the OneDSOC_2S1T parameters.
!!
!! @param QModel             TYPE(QML_OneDSOC_2S1T_t): result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):         type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:                    file unit to read the parameters.
!! @param read_param         logical:                    when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_OneDSOC_2S1T(QModel_in,read_param,nio_param_file,&
                                   Phi_12_in) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_OneDSOC_2S1T_t)                  :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: Phi_12_in

    real (kind=Rkind)      :: Phi_12
    integer                :: err_read

    namelist /OneD_SOC_Model/ Phi_12


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_OneDSOC_2S1T'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 4
    QModel%ndim     = 1
    QModel%pot_name = '1dsoc_2s1t'


    IF (read_param) THEN
      Phi_12   = pi/FOUR
      read(nio_param_file,OneD_SOC_Model,IOSTAT=err_read)

      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "OneD_SOC_Model" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*)
        STOP ' ERROR in Init_QML_OneDSOC_2S1T'
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Some parameter names of the namelist "OneD_SOC_Model" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,nml=OneD_SOC_Model)
        STOP ' ERROR in Init_QML_OneDSOC_2S1T'
      END IF
      !write(out_unitp,OneD_SOC_Model)

      QModel%Phi_12     = Phi_12
    ELSE
      IF (present(Phi_12_in)) THEN
        QModel%Phi_12     = Phi_12_in
      ELSE
        QModel%Phi_12     = pi/FOUR
      END IF
    END IF



    IF (debug) write(out_unitp,*) 'init Q0 of OneDSOC_2S1T'
    QModel%Q0 = [8.5_Rkind]

    IF (debug) write(out_unitp,*) 'init d0GGdef of OneDSOC_2S1T'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,1) = ONE / QModel%mu

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_OneDSOC_2S1T
!> @brief Subroutine wich prints the current QML_OneDSOC_2S1T parameters.
!!
!! @param QModel            CLASS(QML_OneDSOC_2S1T_t): derived type in which the parameters are set-up.
!! @param nio               integer:                    file unit to print the parameters.
  SUBROUTINE Write_QML_OneDSOC_2S1T(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_OneDSOC_2S1T_t),   intent(in) :: QModel
    integer,                       intent(in) :: nio

    write(nio,*) 'QML_OneDSOC_2S1T current parameters'
    write(nio,*)
    write(nio,*) '  a1,a2:          ',QModel%a1,QModel%a2
    write(nio,*) '  alpha1,alpha2:  ',QModel%alpha1,QModel%alpha2
    write(nio,*) '  gamma_1,gamma_2:',QModel%gamma_1,QModel%gamma_2
    write(nio,*) '  DE,DEs:         ',QModel%DE,QModel%DEs
    write(nio,*) '  Phi_12:         ',QModel%Phi_12
    !write(nio,*)
    !write(nio,*) ' option:          ',QModel%option
    write(nio,*)
    write(nio,*) 'end QML_OneDSOC_2S1T current parameters'

  END SUBROUTINE Write_QML_OneDSOC_2S1T
!> @brief Subroutine wich prints the default QML_OneDSOC_2S1T parameters.
!!
!! @param QModel            CLASS(QML_OneDSOC_2S1T_t): derived type in which the parameters are set-up.
!! @param nio               integer:                     file unit to print the parameters.
  SUBROUTINE Write0_QML_OneDSOC_2S1T(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_OneDSOC_2S1T_t),   intent(in) :: QModel
    integer,                       intent(in) :: nio

    write(nio,*) 'OneDSOC_2S1T default parameters, from reference:'
    write(nio,*) 'Granucci et al., J. Chem. Phys. V137, p22A501 (2012)'
    write(nio,*)
    write(nio,*) 'Remark: to avoid complex diabatic potential, ...'
    write(nio,*) '        the complex triplet state components are transformed according to Eq 12:'
    write(nio,*) '        |RT1>=|T,+> = 1/sqrt(2) * [ |T,1> + |T,-1> ]'
    write(nio,*) '        |RT2>=|T,-> = i/sqrt(2) * [ |T,1> - |T,-1> ]'
    write(nio,*) '        |RT3>=|T,z> = i *           |T,0> '
    write(nio,*) 'Furthermore, the triplet components are transformed, so that 4 states are used (see Eq 31).'
    write(nio,*)
    write(nio,*) 'Diabatic Potential, QML_OneDSOC_2S1T (with Phi_12=pi/4=0.7853981....)'
    write(nio,*) 'Value at: R=8.5 Bohr'
    write(nio,*) 'Vdia (Hartree)   = [ 0.006762157       0.001414214       0.000000000       0.000000000]'
    write(nio,*) '                   [ 0.001414214       0.059716484       0.000000000       0.001000000]'
    write(nio,*) '                   [ 0.000000000       0.000000000       0.059716484       0.001000000]'
    write(nio,*) '                   [ 0.000000000       0.001000000       0.001000000       0.041762157]'
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, QML_OneDSOC_2S1T (with Phi_12=pi/4=0.7853981....)'
    write(nio,*) 'Value at: R=8.5 Bohr'
    write(nio,*) 'Vadia (Hartree)   = [ 0.006724396      0.041651620       0.059732248       0.059849020]'
    write(nio,*)
    write(nio,*) 'Non Adiabatic Coupling (with Phi_12=pi/4=0.7853981....)'
    write(nio,*) 'Value at: R=8.5 Bohr'
    write(nio,*) 'NAC              = [ 0.000000000      -0.000602048      -0.004194444       0.005831219]'
    write(nio,*) '                   [ 0.000602048       0.000000000       0.010136691       0.060254264]'
    write(nio,*) '                   [ 0.004194444      -0.010136691       0.000000000       0.079570422]'
    write(nio,*) '                   [-0.005831219      -0.060254264      -0.079570422       0.000000000]'
    write(nio,*)
    write(nio,*) 'end OneDSOC_2S1T default parameters'


  END SUBROUTINE Write0_QML_OneDSOC_2S1T

!> @brief Subroutine wich calculates the OneDSOC_2S1T potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_OneDSOC_2S1T_t): derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):                derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)                 value for which the potential is calculated
!! @param nderiv             integer:                    it enables to specify up to which derivatives the potential is calculated:
!!                                                       the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_OneDSOC_2S1T(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_OneDSOC_2S1T_t), intent(in)    :: QModel
    TYPE (dnS_t),                intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                intent(in)    :: dnQ(:)
    integer,                     intent(in)    :: nderiv


    integer :: i,j

    ! mat_V(4,4)=a1*DEXP(-alpha1*R)+DE                    !Eq 17 without DEs
    Mat_OF_PotDia(4,4)   = QModel%a1 * exp(-QModel%alpha1*dnQ(1)) + &
                           QModel%DE

    ! mat_V(1,1)=a1*DEXP(-alpha1*R)+DE-DEs                !Eq 17
    Mat_OF_PotDia(1,1)   = Mat_OF_PotDia(4,4) - QModel%DEs

    ! mat_V(2,2)=a2*DEXP(-alpha2*R)                       !Eq 18
    ! mat_V(3,3)=a2*DEXP(-alpha2*R)                       !Eq 18
    Mat_OF_PotDia(2,2)   = QModel%a2 * exp(-QModel%alpha2*dnQ(1))
    Mat_OF_PotDia(3,3)   = Mat_OF_PotDia(2,2)


    !Generate non diagonal constant terms
    !mat_V(1,2)=gamma_1
    !mat_V(2,1)=gamma_1
    Mat_OF_PotDia(1,2) = QModel%gamma_1
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


    !mat_V(2,4)=gamma_2*dcos(Phi_12)
    !mat_V(4,2)=gamma_2*dcos(Phi_12)
    Mat_OF_PotDia(2,4) = QModel%gamma_2*cos(QModel%Phi_12)
    Mat_OF_PotDia(4,2) = Mat_OF_PotDia(2,4)


    !mat_V(3,4)=gamma_2*dsin(Phi_12)
    !mat_V(4,3)=gamma_2*dsin(Phi_12)
    Mat_OF_PotDia(3,4) = QModel%gamma_2*sin(QModel%Phi_12)
    Mat_OF_PotDia(4,3) = Mat_OF_PotDia(3,4)

    ! ZERO couplings
    i=1 ; j=3
    Mat_OF_PotDia(i,j) = ZERO
    Mat_OF_PotDia(j,i) = ZERO
    i=1 ; j=4
    Mat_OF_PotDia(i,j) = ZERO
    Mat_OF_PotDia(j,i) = ZERO
    i=2 ; j=3
    Mat_OF_PotDia(i,j) = ZERO
    Mat_OF_PotDia(j,i) = ZERO


  END SUBROUTINE EvalPot_QML_OneDSOC_2S1T

END MODULE QML_OneDSOC_2S1T_m
