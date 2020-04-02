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
!! @date 15/04/2019
!!
MODULE mod_OneDSOC_Model
  USE mod_NumParameters
  IMPLICIT NONE

!> @brief Derived type in which the 1D-Spin Orbit Coupling model parameters are set-up.
!> @brief Reference: Giovanni Granucci, Maurizio Persico, and Gloria Spighi, J. Chem. Phys. V137, p22A501 (2012)
!!
!> @author David Lauvergnat
!! @date 15/04/2019
!!
!! @param a1,a2,alpha1,alpha2,DE  real:    Publication parameters (true fortran parameters)
!! @param C0,RC1,IC1              real:    Publication parameters (true fortran parameters)
!! @param DRsig                   real:    Publication parameters (true fortran parameters)
!! @param Rsig                    real:    Publication parameters (this parameter can be changed)
  TYPE OneDSOCPot_t
     PRIVATE
      real (kind=Rkind) :: a1     = 0.03452_Rkind
      real (kind=Rkind) :: a2     = 0.5_Rkind
      real (kind=Rkind) :: alpha1 = 0.35_Rkind
      real (kind=Rkind) :: alpha2 = 0.25_Rkind
      real (kind=Rkind) :: DE     = 0.04_Rkind
      real (kind=Rkind) :: C0     = 0.001_Rkind
      real (kind=Rkind) :: RC1    = 0.0005_Rkind
      real (kind=Rkind) :: IC1    = 0.0005_Rkind

      ! for Sigma(R)
      real (kind=Rkind) :: DRsig  = 2._Rkind
      real (kind=Rkind) :: Rsig   = 8._Rkind ! only this paramter can be changed

     integer           :: option = 1


     real (kind=Rkind), PUBLIC :: mu  = 20000._Rkind !< Reduced mass from Granucci et al. paper (in au)


  END TYPE OneDSOCPot_t

!  PRIVATE Read_OneDSOC_ModelPot,Init0_OneDSOC_ModelPot,eval_OneDSOC_ModelPot1,eval_OneDSOC_ModelPot2,eval_OneDSOC_ModelPot3
!  PRIVATE eval_OneDSOC_ModelPot1_old,eval_OneDSOC_ModelPot2_old,eval_OneDSOC_ModelPot3_old

CONTAINS
!> @brief Subroutine which makes the initialization of the 1D-SOC Model parameters.
!!
!> @author David Lauvergnat
!! @date 15/04/2019
!!
!! @param OneDSOCPot         TYPE(OneDSOCPot_t):  derived type in which the parameters are set-up.
!! @param option1            integer:            to be able to chose between the models (default 1, as in the pubication).
!! @param nsurf              integer:            4 or 2 (see the publication).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param Rsig               real (optional):    OneDSOC_Model parameter (the only one, which can be changed).
  SUBROUTINE Init_OneDSOC(OneDSOCPot,option,nsurf,nio,read_param,Rsig_in)
    TYPE (OneDSOCPot_t),          intent(inout)   :: OneDSOCPot
    integer,                     intent(in)      :: option
    integer,                     intent(in)      :: nsurf
    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: Rsig_in

    real (kind=Rkind)      :: Rsig
    logical                :: read_param_loc
    integer                :: err_read

    namelist /OneD_SOC_Model/ Rsig

    IF (nsurf /= 2 .AND. nsurf /= 4) THEN
       write(out_unitp,*) ' ERROR in Init_OneDSOC'
       write(out_unitp,*) ' nsurf MUST equal to 4 or 2. nusrf: ',nsurf
       STOP 'ERROR in Init_OneDSOC: nsurf MUST equal to 4 or 2'
    END IF


    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_OneDSOC'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'ERROR in Init_OneDSOC: impossible to read the input file. The file unit (nio) is not present'
    END IF


    IF (read_param_loc) THEN
      Rsig   = 8._Rkind
      read(nio,OneD_SOC_Model,IOSTAT=err_read)

      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in Init_OneDSOC'
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "OneD_SOC_Model" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*)
        STOP ' ERROR in Read_OneDSOC_ModelPot'
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in Init_OneDSOC'
        write(out_unitp,*) ' Some parameter names of the namelist "OneD_SOC_Model" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,nml=OneD_SOC_Model)
        STOP ' ERROR in Init_OneDSOC'
      END IF

      OneDSOCPot%Rsig     = Rsig

    ELSE
      IF (present(Rsig_in)) THEN
        OneDSOCPot%Rsig     = Rsig_in
      ELSE
        OneDSOCPot%Rsig     = 8._Rkind
      END IF
    END IF

    OneDSOCPot%option   = option
    IF (OneDSOCPot%option < 1 .OR. OneDSOCPot%option > 2) OneDSOCPot%option = 1

    !CALL Write_OneDSOC(OneDSOCPot,nio=out_unitp)


  END SUBROUTINE Init_OneDSOC
!> @brief Subroutine wich prints the 1D-SOC Model parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer:            file unit to print the parameters.
  SUBROUTINE Write0_OneDSOC(nio)
    integer,            intent(in) :: nio

    write(nio,*) 'OneDSOC_Model parameters, from reference:'
    write(nio,*)
    write(nio,*) 'Granucci et al., J. Chem. Phys. V137, p22A501 (2012)'
    write(nio,*)
    write(nio,*) 'Remark: to avoid complex diabatic potential, ...'
    write(nio,*) '        the complex triplet state components are transformed according to Eq 12:'
    write(nio,*) '        |RT1>=|T,+> = 1/sqrt(2) * [ |T,1> + |T,-1> ]'
    write(nio,*) '        |RT2>=|T,-> = i/sqrt(2) * [ |T,1> - |T,-1> ]'
    write(nio,*) '        |RT3>=|T,z> = i *           |T,0> '
    write(nio,*)
    write(nio,*) 'Two options (1 and 2):'
    write(nio,*) '     -The first one, as in the publication (the default).'
    write(nio,*) '     -The second one, the function sigma(R) is defined as tanh(-4*(R-Rsig)/DRsig).'
    write(nio,*)
    write(nio,*) 'Values for the first option:'
    write(nio,*) 'Diabatic Potential, OneDSOC_Model (with Rsig=8.0)'
    write(nio,*) 'Value at: R=10. Bohr'
    write(nio,*) 'Vdia (Hartree)   = [ 0.041042414      -0.000707107       0.000707107       0.001000000]'
    write(nio,*) '                   [-0.000707107       0.041042499       0.000000000       0.000000000]'
    write(nio,*) '                   [ 0.000707107       0.000000000       0.041042499       0.000000000]'
    write(nio,*) '                   [ 0.001000000       0.000000000       0.000000000       0.041042499]'
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, OneDSOC_Model (with Rsig=8.0)'
    write(nio,*) 'Value at: R=10. Bohr'
    write(nio,*) 'Vadia (Hartree)   = [ 0.039628243      0.041042499       0.041042499       0.042456670]'

    write(nio,*)
    write(nio,*) 'Non Adiabatic Coupling (with Rsig=8.0)'
    write(nio,*) 'Value at: R=10. Bohr'
    write(nio,*) 'NAC              = [ 0.000000000       0.000000000       0.000000000       1.749343292]'
    write(nio,*) '                   [-0.000000000       0.000000000      -0.125000000       0.000000000]'
    write(nio,*) '                   [-0.000000000       0.125000000       0.000000000       0.000000000]'
    write(nio,*) '                   [-1.749343292      -0.000000000      -0.000000000       0.000000000]'
    write(nio,*) 'WARNING: The NAC, associated to the 2 degenerated vectors, are numerically not well defined !!!!'
    write(nio,*)
    write(nio,*) 'end OneDSOC_Model parameters'

  END SUBROUTINE Write0_OneDSOC

!> @brief Subroutine wich prints the 1D-SOC Model parameters.
!!
!> @author David Lauvergnat
!! @date 15/04/2019
!!
!! @param OneDSOCPot         TYPE(OneDSOCPot_t):  derived type in which the parameters are set-up.
!! @param nio                integer:            file unit to print the parameters.
  SUBROUTINE Write_OneDSOC(OneDSOCPot,nio)
    TYPE (OneDSOCPot_t), intent(in) :: OneDSOCPot
    integer,            intent(in) :: nio

    write(nio,*) 'OneDSOC_Model current parameters'
    write(nio,*)
    write(nio,*) '  a1,a2:          ',OneDSOCPot%a1,OneDSOCPot%a2
    write(nio,*) '  alpha1,alpha2:  ',OneDSOCPot%alpha1,OneDSOCPot%alpha2
    write(nio,*) '  DE:             ',OneDSOCPot%DE

    write(nio,*) '  C0:             ',OneDSOCPot%C0
    write(nio,*) '  C1=RC1 + i IC1: ',OneDSOCPot%RC1,OneDSOCPot%IC1

    write(nio,*) '   Rsig:          ',OneDSOCPot%Rsig
    write(nio,*) '  DRsig:          ',OneDSOCPot%DRsig
    write(nio,*) ' option:          ',OneDSOCPot%option
    write(nio,*)
    write(nio,*) 'end OneDSOC_Model current parameters'

  END SUBROUTINE Write_OneDSOC

  SUBROUTINE get_Q0_OneDSOC(R0,OneDSOCPot)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: R0
    TYPE (OneDSOCPot_t),          intent(in)    :: OneDSOCPot

    R0 = OneDSOCPot%Rsig

  END SUBROUTINE get_Q0_OneDSOC

!> @brief Subroutine wich calculates the 1D-SOC Model potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 15/04/2019
!!
!! @param PotVal       TYPE (dnMat_t):     derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r            real:               value for which the potential is calculated
!! @param OneDSOCPot   TYPE(OneDSOCPot_t):  derived type in which the parameters are set-up.
!! @param nderiv       integer:            it enables to specify up to which derivatives the potential is calculated:
!!                                         the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_OneDSOC(Mat_OF_PotDia,dnR,OneDSOCPot,nderiv)
    USE mod_dnS

    TYPE (OneDSOCPot_t),       intent(in)     :: OneDSOCPot
    TYPE (dnS_t),              intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),              intent(in)     :: dnR
    integer,                  intent(in)     :: nderiv

    !local variables (derived type). They have to be deallocated
    TYPE (dnS_t)      :: dnSig,dnx
    integer          :: i,j,nsurf
    real(kind=Rkind) :: RC01



    ! get nsurf from the matrix: Mat_OF_PotDia
    nsurf = size(Mat_OF_PotDia,dim=1)

    ! Sig(R) calculation
    IF (OneDSOCPot%option == 1) THEN ! as in the publication
      dnSig = dnR ! to have the correct initialization for dnSig = +/- ONE
      IF (dnR <= OneDSOCPot%Rsig-HALF*OneDSOCPot%DRsig) THEN
        dnSig = ONE
      ELSE IF (dnR >= OneDSOCPot%Rsig+HALF*OneDSOCPot%DRsig) THEN
        dnSig = -ONE
      ELSE
        dnx   = (dnR - OneDSOCPot%Rsig) / OneDSOCPot%DRsig
        dnSig = dnx*(FOUR*dnx**2 -THREE)
      END IF
    ELSE
      dnx   = (dnR - OneDSOCPot%Rsig) / OneDSOCPot%DRsig
      dnSig = tanh(-FOUR*dnx)
    END IF

    IF (nsurf == 4) THEN
      !singlet Energy
      Mat_OF_PotDia(1,1) = OneDSOCPot%a1 * exp(-OneDSOCPot%alpha1*dnR) + OneDSOCPot%DE

      !Triplet Energy
      Mat_OF_PotDia(2,2) = OneDSOCPot%a2 * exp(-OneDSOCPot%alpha2*dnR)
      Mat_OF_PotDia(3,3) = Mat_OF_PotDia(2,2)
      Mat_OF_PotDia(4,4) = Mat_OF_PotDia(2,2)


      !Singley-triplet coupling
      Mat_OF_PotDia(1,2) =  sqrt(TWO) * OneDSOCPot%RC1 * dnSig
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)

      Mat_OF_PotDia(1,3) =  -sqrt(TWO) * OneDSOCPot%IC1 * dnSig
      Mat_OF_PotDia(3,1) = Mat_OF_PotDia(1,3)

      Mat_OF_PotDia(1,4) =  -OneDSOCPot%C0 * dnSig
      Mat_OF_PotDia(4,1) = Mat_OF_PotDia(1,4)


     ! triplet-triplet component couplings (ZERO)
     DO i=2,4
     DO j=2,4
       IF (j == i) CYCLE
       Mat_OF_PotDia(j,i) = ZERO
     END DO
     END DO

   ELSE
      !singlet Energy
      Mat_OF_PotDia(1,1) = OneDSOCPot%a1 * exp(-OneDSOCPot%alpha1*dnR) + OneDSOCPot%DE

      !Triplet Energy
      Mat_OF_PotDia(2,2) = OneDSOCPot%a2 * exp(-OneDSOCPot%alpha2*dnR)

      ! S-T coupling
      ! V(S,T) = gamma.Exp[ i.theta]
      ! V(T,S) = gamma.Exp[-i.theta]
      ! gamma(R) = |sigma(R)|.sqrt(2abs(C1)^2+C0^2)
      ! theta(R) = pi.h(R-Rsig) =>
      !    if (R<Rsig) then
      !      theta(R)=0    =>   Exp[ i.theta]=Exp[-i.theta]=  1
      !      sigma(R) > 0  =>   |sigma(R)| =   sigma(R)
      !      => V(S,T)=V(T,S) = sigma(R) . sqrt(2abs(C1)^2+C0^2)
      !
      !    else
      !      theta(R)=pi   =>   Exp[ i.theta]=Exp[-i.theta]= -1
      !      sigma(R) < 0  =>   |sigma(R)| =  - sigma(R)
      !      => V(S,T)=V(T,S) = sigma(R) . sqrt(2abs(C1)^2+C0^2)

      RC01 = sqrt(TWO*(OneDSOCPot%RC1**2+OneDSOCPot%IC1**2)+OneDSOCPot%C0**2)
      Mat_OF_PotDia(1,2) = dnSig * RC01
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)
   END IF


    CALL QML_dealloc_dnS(dnx)
    CALL QML_dealloc_dnS(dnSig)


  END SUBROUTINE eval_OneDSOC

END MODULE mod_OneDSOC_Model
