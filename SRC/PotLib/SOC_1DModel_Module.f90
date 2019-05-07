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

!> @brief Module which makes the initialization, calculation of the 1DSOC_Model potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 15/04/2019
!!
MODULE mod_1DSOC_Model
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
  TYPE Param_1DSOC
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


  END TYPE Param_1DSOC

!  PRIVATE Read_1DSOC_ModelPot,Init0_1DSOC_ModelPot,eval_1DSOC_ModelPot1,eval_1DSOC_ModelPot2,eval_1DSOC_ModelPot3
!  PRIVATE eval_1DSOC_ModelPot1_old,eval_1DSOC_ModelPot2_old,eval_1DSOC_ModelPot3_old

CONTAINS
!> @brief Subroutine which makes the initialization of the 1D-SOC Model parameters.
!!
!> @author David Lauvergnat
!! @date 15/04/2019
!!
!! @param Para_1DSOC         TYPE(Param_1DSOC):  derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the models (default 1, as in the pubication).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param Rsig               real (optional):    1DSOC_Model parameter (the only one, which can be changed).
  SUBROUTINE Init_1DSOC(Para_1DSOC,option,nio,read_param,Rsig_in)
    TYPE (Param_1DSOC),          intent(inout)   :: Para_1DSOC
    integer,                     intent(in)      :: option
    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: Rsig_in

    real (kind=Rkind)      :: Rsig
    logical                :: read_param_loc
    integer                :: err_read

    namelist /OneD_SOC_Model/ Rsig


    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_1DSOC_ModelPot'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'STOP in Init_1DSOC_ModelPot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    Para_1DSOC%option = option

    IF (Para_1DSOC%option < 1 .OR. Para_1DSOC%option > 2) Para_1DSOC%option = 1

    !write(out_unitp,*) 'read_param',read_param,nio


    IF (read_param_loc) THEN
      Rsig = 8._Rkind
      read(nio,OneD_SOC_Model,IOSTAT=err_read)

      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in Init_1DSOC_Model'
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "OneD_SOC_Model" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*)
        STOP ' ERROR in Read_1DSOC_ModelPot'
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in Init_1DSOC_Model'
        write(out_unitp,*) ' Some parameter names of the namelist "OneD_SOC_Model" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,nml=OneD_SOC_Model)
        STOP ' ERROR in Init_1DSOC_Model'
      END IF

      Para_1DSOC%Rsig     = Rsig

    ELSE
      IF (present(Rsig_in)) THEN
        Para_1DSOC%Rsig     = Rsig_in
      ELSE
        Para_1DSOC%Rsig     = 8._Rkind
      END IF
    END IF


  END SUBROUTINE Init_1DSOC


!> @brief Subroutine wich prints the 1D-SOC Model parameters.
!!
!> @author David Lauvergnat
!! @date 15/04/2019
!!
!! @param Para_1DSOC         TYPE(Param_1DSOC):  derived type in which the parameters are set-up.
!! @param nio                integer:            file unit to print the parameters.
  SUBROUTINE Write_1DSOC(Para_1DSOC,nio)
    TYPE (Param_1DSOC), intent(in) :: Para_1DSOC
    integer,            intent(in) :: nio

    write(nio,*) '1DSOC_Model parameters, from reference:'
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
    write(nio,*) 'Diabatic Potential, 1DSOC_Model (with Rsig=8.0)'
    write(nio,*) 'Value at: R=10. Bohr'
    write(nio,*) 'Vdia (Hartree)   = [ 0.041042414      -0.000707107       0.000707107       0.001000000]'
    write(nio,*) '                   [-0.000707107       0.041042499       0.000000000       0.000000000]'
    write(nio,*) '                   [ 0.000707107       0.000000000       0.041042499       0.000000000]'
    write(nio,*) '                   [ 0.001000000       0.000000000       0.000000000       0.041042499]'

    write(nio,*)
    write(nio,*) 'Adiabatic Potential, 1DSOC_Model (with Rsig=8.0)'
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


    write(nio,*) 'Current parameters:'
    write(nio,*) '  a1,a2:          ',Para_1DSOC%a1,Para_1DSOC%a2
    write(nio,*) '  alpha1,alpha2:  ',Para_1DSOC%alpha1,Para_1DSOC%alpha2
    write(nio,*) '  DE:             ',Para_1DSOC%DE

    write(nio,*) '  C0:             ',Para_1DSOC%C0
    write(nio,*) '  C1=RC1 + i IC1: ',Para_1DSOC%RC1,Para_1DSOC%IC1

    write(nio,*) '   Rsig:          ',Para_1DSOC%Rsig
    write(nio,*) '  DRsig:          ',Para_1DSOC%DRsig

    write(nio,*) 'end 1DSOC_Model parameters'

  END SUBROUTINE Write_1DSOC

!> @brief Subroutine wich calculates the 1D-SOC Model potential with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 15/04/2019
!!
!! @param PotVal       TYPE(dnMatPot):     derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r            real:               value for which the potential is calculated
!! @param Para_1DSOC   TYPE(Param_1DSOC):  derived type in which the parameters are set-up.
!! @param nderiv       integer:            it enables to specify up to which derivatives the potential is calculated:
!!                                         the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_1DSOC(PotVal,r,Para_1DSOC,nderiv)
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_1DSOC),       intent(in)     :: Para_1DSOC
    real (kind=Rkind),        intent(in)     :: r
    TYPE(dnMatPot),           intent(inout)  :: PotVal
    integer,                  intent(in)     :: nderiv

    !local variables (derived type). They have to be deallocated
    TYPE(dnSca)     :: dnPot,dnR,dnSig,dnx

    PotVal = ZERO

    dnR     = init_dnSca(r,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives

    ! Sig(R) calculation
    IF (Para_1DSOC%option == 1) THEN ! as in the publication
      dnSig = dnR ! to have the correct initialization for dnSig = +/- ONE
      IF (dnR <= Para_1DSOC%Rsig-HALF*Para_1DSOC%DRsig) THEN
        dnSig = ONE
      ELSE IF (dnR >= Para_1DSOC%Rsig+HALF*Para_1DSOC%DRsig) THEN
        dnSig = -ONE
      ELSE
        dnx   = (dnR - Para_1DSOC%Rsig) / Para_1DSOC%DRsig
        dnSig = dnx*(FOUR*dnx**2 -THREE)
      END IF
    ELSE
      dnx   = (dnR - Para_1DSOC%Rsig) / Para_1DSOC%DRsig
      dnSig = tanh(-FOUR*dnx)
    END IF

    !singlet Energy
    dnPot =  Para_1DSOC%a1 * exp(-Para_1DSOC%alpha1*dnR) + Para_1DSOC%DE
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=1)

    !Triplet Energy
    dnPot =  Para_1DSOC%a2 * exp(-Para_1DSOC%alpha2*dnR)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=2,j=2)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=3,j=3)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=4,j=4)

    !Singley-triplet coupling
    dnPot =  sqrt(TWO) * Para_1DSOC%RC1 * dnSig
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=2)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=2,j=1)

    dnPot =  -sqrt(TWO) * Para_1DSOC%IC1 * dnSig
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=3)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=3,j=1)

    dnPot =  -Para_1DSOC%C0 * dnSig
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=4)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=4,j=1)



    CALL dealloc_dnSca(dnR)
    CALL dealloc_dnSca(dnx)
    CALL dealloc_dnSca(dnSig)
    CALL dealloc_dnSca(dnPot)


  END SUBROUTINE eval_1DSOC

END MODULE mod_1DSOC_Model
