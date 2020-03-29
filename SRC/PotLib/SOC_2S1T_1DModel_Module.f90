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

!> @brief Module which makes the initialization, calculation of the OneDSOC_Model potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 27/06/2019
!!
MODULE mod_OneDSOC_2S1T_Model
  USE mod_NumParameters
  IMPLICIT NONE

!> @brief Derived type in which the 1D-Spin Orbit Coupling model parameters are set-up.
!> @brief Reference: Giovanni Granucci, Maurizio Persico, and Gloria Spighi, J. Chem. Phys. V137, p22A501 (2012)
!!
!> @author David Lauvergnat
!! @date 27/06/2019
!!
!! @param a1,a2,alpha1,alpha2,DE  real:    Publication parameters (true fortran parameters)
!! @param gamma_1,gamma_2         real:    Publication parameters (true fortran parameters)
!! @param DE,DEs                  real:    Publication parameters (true fortran parameters)
!! @param Phi_12                  real:    Publication parameters (this parameter can be changed)
  TYPE OneDSOC_2S1T_Pot_t
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

     integer            :: option = 1

     real (kind=Rkind), PUBLIC :: mu  = 20000._Rkind !< Reduced mass from Granucci et al. paper (in au)


  END TYPE OneDSOC_2S1T_Pot_t

CONTAINS
  SUBROUTINE Init_OneDSOC_2S1T(OneDSOC_2S1T_Pot,option,nio,read_param,Phi_12_in)
    TYPE (OneDSOC_2S1T_Pot_t),     intent(inout)   :: OneDSOC_2S1T_Pot
    integer,                     intent(in)      :: option
    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: Phi_12_in

    real (kind=Rkind)      :: Phi_12
    logical                :: read_param_loc
    integer                :: err_read

    namelist /OneD_SOC_Model/ Phi_12


    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_OneDSOC_2S1T'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'STOP in Init_OneDSOC_2S1T: impossible to read the input file. The file unit (nio) is not present'
    END IF


    IF (read_param_loc) THEN
      Phi_12   = pi/FOUR
      read(nio,OneD_SOC_Model,IOSTAT=err_read)

      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in Init_OneDSOC_2S1T'
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "OneD_SOC_Model" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*)
        STOP ' ERROR in Init_OneDSOC_2S1T'
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in Init_OneDSOC_2S1T'
        write(out_unitp,*) ' Some parameter names of the namelist "OneD_SOC_Model" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,nml=OneD_SOC_Model)
        STOP ' ERROR in Init_OneDSOC_2S1T'
      END IF
      write(out_unitp,OneD_SOC_Model)

      OneDSOC_2S1T_Pot%Phi_12     = Phi_12
    ELSE
      IF (present(Phi_12_in)) THEN
        OneDSOC_2S1T_Pot%Phi_12     = Phi_12_in
      ELSE
        OneDSOC_2S1T_Pot%Phi_12     = pi/FOUR
      END IF
    END IF

    OneDSOC_2S1T_Pot%option   = option
    !IF (OneDSOC_2S1T_Pot%option < 1 .OR. OneDSOC_2S1T_Pot%option > 2) OneDSOC_2S1T_Pot%option = 1

    CALL Write_OneDSOC_2S1T(OneDSOC_2S1T_Pot,nio=out_unitp)


  END SUBROUTINE Init_OneDSOC_2S1T

!> @brief Subroutine wich prints the 1D-SOC Model parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer:                 file unit to print the parameters.
  SUBROUTINE Write0_OneDSOC_2S1T(nio)
    integer,            intent(in) :: nio

    write(nio,*) 'OneDSOC_2S1T_Model parameters, from reference:'
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
    write(nio,*) 'end OneDSOC_2S1T_Model parameters'

  END SUBROUTINE Write0_OneDSOC_2S1T
!> @brief Subroutine wich prints the 1D-SOC Model parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param OneDSOC_2S1T_Pot    TYPE(OneDSOC_2S1T_Pot_t):  derived type in which the parameters are set-up.
!! @param nio                integer:                 file unit to print the parameters.
  SUBROUTINE Write_OneDSOC_2S1T(OneDSOC_2S1T_Pot,nio)
    TYPE (OneDSOC_2S1T_Pot_t), intent(in) :: OneDSOC_2S1T_Pot
    integer,            intent(in) :: nio

    write(nio,*) 'OneDSOC_2S1T_Model current parameters'
    write(nio,*)
    write(nio,*) '  a1,a2:          ',OneDSOC_2S1T_Pot%a1,OneDSOC_2S1T_Pot%a2
    write(nio,*) '  alpha1,alpha2:  ',OneDSOC_2S1T_Pot%alpha1,OneDSOC_2S1T_Pot%alpha2
    write(nio,*) '  gamma_1,gamma_2:',OneDSOC_2S1T_Pot%gamma_1,OneDSOC_2S1T_Pot%gamma_2
    write(nio,*) '  DE,DEs:         ',OneDSOC_2S1T_Pot%DE,OneDSOC_2S1T_Pot%DEs
    write(nio,*) '  Phi_12:         ',OneDSOC_2S1T_Pot%Phi_12
    !write(nio,*)
    !write(nio,*) ' option:          ',OneDSOC_2S1T_Pot%option
    write(nio,*)
    write(nio,*) 'end OneDSOC_2S1T_Model current parameters'

  END SUBROUTINE Write_OneDSOC_2S1T

  SUBROUTINE get_Q0_OneDSOC_2S1T(R0,OneDSOC_2S1T_Pot)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: R0
    TYPE (OneDSOC_2S1T_Pot_t),     intent(in)    :: OneDSOC_2S1T_Pot

    R0 = ZERO

  END SUBROUTINE get_Q0_OneDSOC_2S1T

!> @brief Subroutine wich calculates the 1D-SOC Model potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 27/06/2019
!!
!! @param PotVal            TYPE (dnMat_t):     derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                 real:               value for which the potential is calculated
!! @param OneDSOC_2S1T_Pot   TYPE(OneDSOC_2S1T_Pot_t):  derived type in which the parameters are set-up.
!! @param nderiv            integer:            it enables to specify up to which derivatives the potential is calculated:
!!                                              the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_OneDSOC_2S1T(Mat_OF_PotDia,dnR,OneDSOC_2S1T_Pot,nderiv)
    USE mod_dnS

    TYPE (OneDSOC_2S1T_Pot_t),  intent(in)     :: OneDSOC_2S1T_Pot
    TYPE (dnS_t),              intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),              intent(in)     :: dnR
    integer,                  intent(in)     :: nderiv


    integer :: i,j

    ! mat_V(4,4)=a1*DEXP(-alpha1*R)+DE                    !Eq 17 without DEs
    Mat_OF_PotDia(4,4)   = OneDSOC_2S1T_Pot%a1 * exp(-OneDSOC_2S1T_Pot%alpha1*dnR) + &
                           OneDSOC_2S1T_Pot%DE

    ! mat_V(1,1)=a1*DEXP(-alpha1*R)+DE-DEs                !Eq 17
    Mat_OF_PotDia(1,1)   = Mat_OF_PotDia(4,4) - OneDSOC_2S1T_Pot%DEs

    ! mat_V(2,2)=a2*DEXP(-alpha2*R)                       !Eq 18
    ! mat_V(3,3)=a2*DEXP(-alpha2*R)                       !Eq 18
    Mat_OF_PotDia(2,2)   = OneDSOC_2S1T_Pot%a2 * exp(-OneDSOC_2S1T_Pot%alpha2*dnR)
    Mat_OF_PotDia(3,3)   = Mat_OF_PotDia(2,2)




    !Generate non diagonal constant terms
    !mat_V(1,2)=gamma_1
    !mat_V(2,1)=gamma_1
    Mat_OF_PotDia(1,2) = OneDSOC_2S1T_Pot%gamma_1
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


    !mat_V(2,4)=gamma_2*dcos(Phi_12)
    !mat_V(4,2)=gamma_2*dcos(Phi_12)
    Mat_OF_PotDia(2,4) = OneDSOC_2S1T_Pot%gamma_2*cos(OneDSOC_2S1T_Pot%Phi_12)
    Mat_OF_PotDia(4,2) = Mat_OF_PotDia(2,4)


    !mat_V(3,4)=gamma_2*dsin(Phi_12)
    !mat_V(4,3)=gamma_2*dsin(Phi_12)
    Mat_OF_PotDia(3,4) = OneDSOC_2S1T_Pot%gamma_2*sin(OneDSOC_2S1T_Pot%Phi_12)
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

  END SUBROUTINE eval_OneDSOC_2S1T

END MODULE mod_OneDSOC_2S1T_Model
