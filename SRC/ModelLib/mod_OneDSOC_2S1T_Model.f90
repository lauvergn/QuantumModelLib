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
!    Copyright 2016  David LAUVERGNAT, Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================

!> @brief Module which makes the initialization, calculation of the OneDSOC_2S1T potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE mod_OneDSOC_2S1T_Model
  USE mod_NumParameters
  USE mod_EmptyModel
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the OneDSOC_2S1T parameters are set-up.
  TYPE, EXTENDS (EmptyModel_t) ::  OneDSOC_2S1T_Model_t

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
    PROCEDURE :: Eval_QModel_Pot => eval_OneDSOC_2S1T_Pot
    PROCEDURE :: Write_QModel    => Write_OneDSOC_2S1T_Model
    PROCEDURE :: Write0_QModel   => Write0_OneDSOC_2S1T_Model
  END TYPE OneDSOC_2S1T_Model_t
 
  PUBLIC :: OneDSOC_2S1T_Model_t,Init_OneDSOC_2S1T_Model
 
  CONTAINS
!> @brief Function which makes the initialization of the OneDSOC_2S1T parameters.
!!
!! @param QModel             TYPE(OneDSOC_2S1T_Model_t): result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(EmptyModel_t):         type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:                    file unit to read the parameters.
!! @param read_param         logical:                    when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_OneDSOC_2S1T_Model(QModel_in,read_param,nio_param_file,&
                                   Phi_12_in) RESULT(QModel)
  IMPLICIT NONE

    TYPE (OneDSOC_2S1T_Model_t)                  :: QModel ! RESULT

    TYPE(EmptyModel_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: Phi_12_in

    real (kind=Rkind)      :: Phi_12
    integer                :: err_read

    namelist /OneD_SOC_Model/ Phi_12


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_OneDSOC_2S1T_Model'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_EmptyModel(QModel%EmptyModel_t,QModel_in)

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
        STOP ' ERROR in Init_OneDSOC_2S1T_Model'
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Some parameter names of the namelist "OneD_SOC_Model" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,nml=OneD_SOC_Model)
        STOP ' ERROR in Init_OneDSOC_2S1T_Model'
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
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
    QModel%d0GGdef(1,1) = ONE / QModel%mu

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_OneDSOC_2S1T_Model
!> @brief Subroutine wich prints the current OneDSOC_2S1T_Model parameters.
!!
!! @param QModel            CLASS(OneDSOC_2S1T_Model_t): derived type in which the parameters are set-up.
!! @param nio               integer:                    file unit to print the parameters.
  SUBROUTINE Write_OneDSOC_2S1T_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(OneDSOC_2S1T_Model_t),   intent(in) :: QModel
    integer,                       intent(in) :: nio

    write(nio,*) 'OneDSOC_2S1T_Model current parameters'
    write(nio,*)
    write(nio,*) '  a1,a2:          ',QModel%a1,QModel%a2
    write(nio,*) '  alpha1,alpha2:  ',QModel%alpha1,QModel%alpha2
    write(nio,*) '  gamma_1,gamma_2:',QModel%gamma_1,QModel%gamma_2
    write(nio,*) '  DE,DEs:         ',QModel%DE,QModel%DEs
    write(nio,*) '  Phi_12:         ',QModel%Phi_12
    !write(nio,*)
    !write(nio,*) ' option:          ',QModel%option
    write(nio,*)
    write(nio,*) 'end OneDSOC_2S1T_Model current parameters'

  END SUBROUTINE Write_OneDSOC_2S1T_Model
!> @brief Subroutine wich prints the default OneDSOC_2S1T_Model parameters.
!!
!! @param QModel            CLASS(OneDSOC_2S1T_Model_t): derived type in which the parameters are set-up.
!! @param nio               integer:                     file unit to print the parameters.
  SUBROUTINE Write0_OneDSOC_2S1T_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(OneDSOC_2S1T_Model_t),   intent(in) :: QModel
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
    write(nio,*) 'Diabatic Potential, OneDSOC_2S1T_Model (with Phi_12=pi/4=0.7853981....)'
    write(nio,*) 'Value at: R=8.5 Bohr'
    write(nio,*) 'Vdia (Hartree)   = [ 0.006762157       0.001414214       0.000000000       0.000000000]'
    write(nio,*) '                   [ 0.001414214       0.059716484       0.000000000       0.001000000]'
    write(nio,*) '                   [ 0.000000000       0.000000000       0.059716484       0.001000000]'
    write(nio,*) '                   [ 0.000000000       0.001000000       0.001000000       0.041762157]'
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, OneDSOC_2S1T_Model (with Phi_12=pi/4=0.7853981....)'
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


  END SUBROUTINE Write0_OneDSOC_2S1T_Model

!> @brief Subroutine wich calculates the OneDSOC_2S1T potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(OneDSOC_2S1T_Model_t): derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):                derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)                 value for which the potential is calculated
!! @param nderiv             integer:                    it enables to specify up to which derivatives the potential is calculated:
!!                                                       the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_OneDSOC_2S1T_Pot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(OneDSOC_2S1T_Model_t), intent(in)    :: QModel
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


  END SUBROUTINE eval_OneDSOC_2S1T_Pot

END MODULE mod_OneDSOC_2S1T_Model
