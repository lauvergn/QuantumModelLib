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

!> @brief Module which makes the initialization, calculation of the Tully potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!> @author Félix MOUHAT
!! @date 03/08/2017
!!
MODULE mod_TullyPot
  USE mod_NumParameters
  IMPLICIT NONE

!> @brief Derived type in which the Tully parameters are set-up.
!> @brief Reference: Tully, J. Chem. Phys. V93, pp15, 1990
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param A,B,C,D,e0              real:    Tully parameters
!! @param mu                      real:   mass used in Tully paper

!! @param option                  integer: it enables to chose between the 3 models (default 1, Simple avoided crossing)
  TYPE Param_Tully
     PRIVATE
     real (kind=Rkind) :: A      = 0.01_Rkind
     real (kind=Rkind) :: B      = 1.6_Rkind
     real (kind=Rkind) :: C      = 0.005_Rkind
     real (kind=Rkind) :: D      = 1.0_Rkind
     real (kind=Rkind) :: E0     = 0.0_Rkind

     real (kind=Rkind), PUBLIC :: mu     = 2000._Rkind


     integer           :: option = 1
     
  END TYPE Param_Tully

  PRIVATE Read_TullyPot,Init0_TullyPot,eval_TullyPot1,eval_TullyPot2,eval_TullyPot3
  PRIVATE eval_TullyPot1_old,eval_TullyPot2_old,eval_TullyPot3_old

CONTAINS
!> @brief Subroutine which makes the initialization of the Tully parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Tully         TYPE(Param_Tully):  derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param A,B,C,D,E0         real (optional):    Tully parameters.
  SUBROUTINE Init_TullyPot(Para_Tully,option,nio,read_param,A,B,C,D,E0)
    TYPE (Param_Tully),          intent(inout)   :: Para_Tully
    integer,                     intent(in)      :: option
    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: A,B,C,D,E0

    logical :: read_param_loc

    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_TullyPot'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'STOP in Init_TullyPot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    Para_Tully%option = option

    IF (Para_Tully%option < 1 .OR. Para_Tully%option > 3) Para_Tully%option = 1

    !write(out_unitp,*) 'option',option
    !write(out_unitp,*) 'read_param',read_param,nio


    IF (read_param_loc) THEN
      CALL Read_TullyPot(Para_Tully,nio)
      CALL Init0_TullyPot(Para_Tully)
    ELSE

      Para_Tully%A      = Huge(ONE)
      Para_Tully%B      = Huge(ONE)
      Para_Tully%C      = Huge(ONE)
      Para_Tully%D      = Huge(ONE)
      Para_Tully%E0     = Huge(ONE)

      CALL Init0_TullyPot(Para_Tully)

      SELECT CASE (Para_Tully%option)
      CASE (1) ! A. Simple avoided crossing

        IF (present(A)) Para_Tully%A = A
        IF (present(B)) Para_Tully%B = B
        IF (present(C)) Para_Tully%C = C
        IF (present(D)) Para_Tully%D = D

      CASE (2) ! B. Dual avoided crossing

        IF (present(A))  Para_Tully%A  = A
        IF (present(B))  Para_Tully%B  = B
        IF (present(E0)) Para_Tully%E0 = E0
        IF (present(C))  Para_Tully%C  = C
        IF (present(D))  Para_Tully%D  = D

      CASE (3) !C. Extended coupling with reflection

        IF (present(A)) Para_Tully%A = A
        IF (present(B)) Para_Tully%B = B
        IF (present(C)) Para_Tully%C = C

      CASE Default
          write(out_unitp,*) 'ERROR in Init_TullyPot'
          write(out_unitp,*) ' This option is not possible. option:',Para_Tully%option
          write(out_unitp,*) ' Its value MUST be 1 or 3 or 3'
          STOP
      END SELECT
    END IF


  END SUBROUTINE Init_TullyPot
!> @brief Subroutine which makes the initialization of the Tully parameters from his paper.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Tully         TYPE(Param_Tully):  derived type in which the parameters are set-up.
  SUBROUTINE Init0_TullyPot(Para_Tully)
    TYPE (Param_Tully),      intent(inout)   :: Para_Tully

    character (len=*), parameter :: name_sub='Init0_TullyPot'

    !write(out_unitp,*) 'BEGINNING ',name_sub
    !write(out_unitp,*) 'Para_Tully%option',Para_Tully%option

    IF (Para_Tully%option < 1 .OR. Para_Tully%option > 3) Para_Tully%option = 1

    SELECT CASE (Para_Tully%option)
    CASE (1) ! A. Simple avoided crossing
      IF (Para_Tully%A == huge(ONE)) Para_Tully%A = 0.01_Rkind
      IF (Para_Tully%B == huge(ONE)) Para_Tully%B = 1.6_Rkind
      IF (Para_Tully%C == huge(ONE)) Para_Tully%C = 0.005_Rkind
      IF (Para_Tully%D == huge(ONE)) Para_Tully%D = 1.0_Rkind

    CASE (2) ! B. Dual avoided crossing
      IF (Para_Tully%A  == huge(ONE)) Para_Tully%A  = 0.1_Rkind
      IF (Para_Tully%B  == huge(ONE)) Para_Tully%B  = 0.28_Rkind
      IF (Para_Tully%E0 == huge(ONE)) Para_Tully%E0 = 0.05_Rkind
      IF (Para_Tully%C  == huge(ONE)) Para_Tully%C  = 0.015_Rkind
      IF (Para_Tully%D  == huge(ONE)) Para_Tully%D  = 0.06_Rkind

    CASE (3) !C. Extended coupling with reflection
      IF (Para_Tully%A == huge(ONE)) Para_Tully%A = 0.0006_Rkind
      IF (Para_Tully%B == huge(ONE)) Para_Tully%B = 0.1_Rkind
      IF (Para_Tully%C == huge(ONE)) Para_Tully%C = 0.90_Rkind

    CASE Default
        write(out_unitp,*) 'ERROR in Init0_TullyPot'
        write(out_unitp,*) ' This option is not possible. option:',Para_Tully%option
        write(out_unitp,*) ' Its value MUST be 1 or 2 or 3'
        STOP
    END SELECT

    !CALL Write_TullyPot(Para_Tully,out_unitp)
    !write(out_unitp,*) 'END ',name_sub


  END SUBROUTINE Init0_TullyPot

!> @brief Subroutine wich reads the Tully parameters with a namelist.
!!   This can be called only from the "Init_TullyPot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Tully         TYPE(Param_Tully):   derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to read the parameters.
  SUBROUTINE Read_TullyPot(Para_Tully,nio)
    TYPE (Param_Tully),      intent(inout)   :: Para_Tully
    integer,                 intent(in)      :: nio

    real (kind=Rkind)      :: A,B,C,D,E0
    integer                :: err_read

    namelist /Tully/ A,B,C,D,E0


    A      = Huge(ONE)
    B      = Huge(ONE)
    C      = Huge(ONE)
    D      = Huge(ONE)
    E0     = Huge(ONE)

    read(nio,Tully,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_TullyPot'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "Tully" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_TullyPot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_TullyPot'
      write(out_unitp,*) ' Some parameter names of the namelist "Tully" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=Tully)
      STOP ' ERROR in Read_TullyPot'
    END IF

    Para_Tully%A       = A
    Para_Tully%B       = B
    Para_Tully%E0      = E0
    Para_Tully%C       = C
    Para_Tully%D       = D

  END SUBROUTINE Read_TullyPot

!> @brief Subroutine wich prints the Tully parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Tully         TYPE(Param_Tully):   derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write_TullyPot(Para_Tully,nio)
    TYPE (Param_Tully), intent(in) :: Para_Tully
    integer, intent(in) :: nio

    write(nio,*) 'Tully parameters, from reference:'
    write(nio,*) '  Reference: Tully, J. Chem. Phys. V93, pp15, 1990'
    write(nio,*)

    write(nio,*) 'Adiabatic Potential, Tully 1'
    write(nio,*) 'Value at: R=-1.1 Bohr'
    write(nio,*) 'V        = [-0.008413,  0.008413] Hartree'
    write(nio,*) 'gradient = [ 0.002128, -0.002128] '
    write(nio,*) 'hessian  = [ 0.001943, -0.001943] '
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, Tully 2'
    write(nio,*) 'Value at: R=-1.1 Bohr'
    write(nio,*) 'V        = [-0.028170,  0.006908] Hartree'
    write(nio,*) 'gradient = [-0.036718, -0.007180] '
    write(nio,*) 'hessian  = [-0.003754,  0.016620] '
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, Tully 3'
    write(nio,*) 'Value at: R=-1.1 Bohr'
    write(nio,*) 'V        = [-0.037163,  0.037163] Hartree'
    write(nio,*) 'gradient = [-0.033438,  0.033438] '
    write(nio,*) 'hessian  = [-0.030102,  0.030102] '
    write(nio,*)

    write(nio,*) 'Current parameters:'
    write(nio,*) '  A:      ',Para_Tully%A
    write(nio,*) '  B:      ',Para_Tully%B
    write(nio,*) '  C:      ',Para_Tully%C
    write(nio,*) '  D:      ',Para_Tully%D
    write(nio,*) '  E0:     ',Para_Tully%E0
    write(nio,*) '  option: ',Para_Tully%option

    write(nio,*) 'end Tully parameters'

  END SUBROUTINE Write_TullyPot

!> @brief Subroutine wich calculates the Tully potential (for the 3 models) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_Tully         TYPE(Param_Tully):   derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TullyPot(PotVal,r,Para_Tully,nderiv)
    USE mod_dnMatPot

    TYPE (Param_Tully), intent(in)     :: Para_Tully
    real (kind=Rkind),  intent(in)     :: r
    TYPE(dnMatPot),     intent(inout)  :: PotVal
    integer, intent(in)                :: nderiv

    SELECT CASE (Para_Tully%option)
    CASE (1)
      CALL eval_TullyPot1(PotVal,r,Para_Tully,nderiv)
    CASE (2)
      CALL eval_TullyPot2(PotVal,r,Para_Tully,nderiv)
    CASE (3)
      CALL eval_TullyPot3(PotVal,r,Para_Tully,nderiv)
    CASE Default
        write(out_unitp,*) 'ERROR in eval_TullyPot'
        write(out_unitp,*) ' This option is not possible. option:',Para_Tully%option
        write(out_unitp,*) ' Its value MUST be 1 or 3 or 3'
        STOP
    END SELECT

  END SUBROUTINE eval_TullyPot
!> @brief Subroutine wich calculates the Tully potential (for the 1st model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_Tully         TYPE(Param_Tully):   derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TullyPot1(PotVal,r,Para_Tully,nderiv)
  ! A. Simple avoided crossing
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_Tully), intent(in)     :: Para_Tully
    real (kind=Rkind),  intent(in)     :: r
    TYPE(dnMatPot),     intent(inout)  :: PotVal
    integer, intent(in)                :: nderiv


    !local variables (derived type). They have to be deallocated
    TYPE(dnSca)     :: dnPot,dnR

    dnR     = init_dnSca(r,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives

    ! Potential calculation
    IF (dnR >= ZERO) THEN
      dnPot =  Para_Tully%A*(ONE-exp(-Para_Tully%B*dnR))
    ELSE
      dnPot = -Para_Tully%A*(ONE-exp( Para_Tully%B*dnR))
    END IF
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=1)

    CALL sub_dnSca_TO_dnMatPot(-dnPot,PotVal,i=2,j=2)


    dnPot = Para_Tully%C*exp(-Para_Tully%D*dnR**2)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=2)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=2,j=1)


    CALL dealloc_dnSca(dnR)
    CALL dealloc_dnSca(dnPot)

  END SUBROUTINE eval_TullyPot1
!> @brief Subroutine wich calculates the Tully potential (for the 2d model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_Tully         TYPE(Param_Tully):   derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TullyPot2(PotVal,r,Para_Tully,nderiv) !2d Tully's potential
  !  B. Dual avoided crossing
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_Tully), intent(in)     :: Para_Tully
    real (kind=Rkind),  intent(in)     :: r
    TYPE(dnMatPot),     intent(inout)  :: PotVal
    integer, intent(in)                :: nderiv


    !local variables (derived type). They have to be deallocated
    TYPE(dnSca)     :: dnPot,dnR

    dnR     = init_dnSca(r,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives


! Potential calculation
    PotVal         = ZERO ! to initialized the PotVal%d0(1,1) and its derivatives

    dnPot = -Para_Tully%A*exp(-Para_Tully%B*dnR**2)+Para_Tully%E0
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=2,j=2)

    dnPot =  Para_Tully%C*exp(-Para_Tully%D*dnR**2)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=2)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=2,j=1)


    CALL dealloc_dnSca(dnR)
    CALL dealloc_dnSca(dnPot)

  END SUBROUTINE eval_TullyPot2
!> @brief Subroutine wich calculates the Tully potential (for the 3d model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param Para_Tully         TYPE(Param_Tully):   derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TullyPot3(PotVal,r,Para_Tully,nderiv) !3d Tully's potential
    !  C. Extended coupling with reflection
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_Tully), intent(in)     :: Para_Tully
    real (kind=Rkind),  intent(in)     :: r
    TYPE(dnMatPot),     intent(inout)  :: PotVal
    integer, intent(in)                :: nderiv


    !local variables (derived type). They have to be deallocated
    TYPE(dnSca)     :: dnPot,dnR

    dnR     = init_dnSca(r,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives


! Potential calculation
    PotVal         = ZERO ! to initialized the PotVal%d0(1,1) and PotVal%d0(2,2) and their derivatives
    PotVal%d0(1,1) =  Para_Tully%A
    PotVal%d0(2,2) = -Para_Tully%A

    IF (dnR >= ZERO) THEN
      dnPot = Para_Tully%B*(TWO-exp(-Para_Tully%C*dnR))
    ELSE
      dnPot = Para_Tully%B*exp(Para_Tully%C*dnr)
    END IF
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=1,j=2)
    CALL sub_dnSca_TO_dnMatPot(dnPot,PotVal,i=2,j=1)


    CALL dealloc_dnSca(dnR)
    CALL dealloc_dnSca(dnPot)

  END SUBROUTINE eval_TullyPot3


  SUBROUTINE eval_TullyPot1_old(PotVal,r,Para_Tully,nderiv)
  ! A. Simple avoided crossing
    USE mod_dnMatPot

    TYPE (Param_Tully), intent(in)     :: Para_Tully
    real (kind=Rkind),  intent(in)     :: r
    TYPE(dnMatPot),     intent(inout)  :: PotVal
    integer, intent(in)                :: nderiv

!write(out_unitp,*) 'coucou Tully 1',r
! Potential calculation
    if (r >= ZERO) then
      PotVal%d0(1,1) =  Para_Tully%A*(ONE-exp(-Para_Tully%B*r))
    else
      PotVal%d0(1,1) = -Para_Tully%A*(ONE-exp(Para_Tully%B*r))
    endif
    PotVal%d0(2,2) = -PotVal%d0(1,1)

    PotVal%d0(1,2) = Para_Tully%C*exp(-Para_Tully%D*r**2)
    PotVal%d0(2,1) = PotVal%d0(1,2)

! gradient calculation
    if (nderiv >= 1) then
         if (r >= ZERO) then
           PotVal%d1(1,1,1) = Para_Tully%A*Para_Tully%B*exp(-Para_Tully%B*r)
         else
           PotVal%d1(1,1,1) = Para_Tully%A*Para_Tully%B*exp(Para_Tully%B*r)
         endif
         PotVal%d1(2,2,1) = - PotVal%d1(1,1,1)

         PotVal%d1(1,2,1) = -TWO*Para_Tully%D*r*PotVal%d0(1,2)
         PotVal%d1(2,1,1) = PotVal%d1(1,2,1)
    endif

! Hessian calculation
    if (nderiv >= 2) then
        if (r > ZERO) then
          PotVal%d2(1,1,1,1) = -Para_Tully%B*PotVal%d1(1,1,1)
        else
          PotVal%d2(1,1,1,1) = Para_Tully%B*PotVal%d1(1,1,1)
        endif
        PotVal%d2(2,2,1,1) = -PotVal%d2(1,1,1,1)

        PotVal%d2(1,2,1,1) = (-TWO*Para_Tully%D+FOUR*Para_Tully%D**2*r**2)*PotVal%d0(1,2)
        PotVal%d2(2,1,1,1) = PotVal%d2(1,2,1,1)
    endif

  END SUBROUTINE eval_TullyPot1_old
  SUBROUTINE eval_TullyPot2_old(PotVal,r,Para_Tully,nderiv) !2nd Tully's potential
  !  B. Dual avoided crossing
    USE mod_dnMatPot

    TYPE (Param_Tully), intent(in)     :: Para_Tully
    real (kind=Rkind),  intent(in)     :: r
    TYPE(dnMatPot),     intent(inout)  :: PotVal
    integer, intent(in)                :: nderiv

! Potential calculation
    PotVal%d0(1,1) =  ZERO
    PotVal%d0(2,2) = -Para_Tully%A*exp(-Para_Tully%B*r**2)+Para_Tully%E0
    PotVal%d0(1,2) =  Para_Tully%C*exp(-Para_Tully%D*r**2)
    PotVal%d0(2,1) =  PotVal%d0(1,2)

! Gradient calculation
    if (nderiv .ge. 1) then
      PotVal%d1(1,1,1) =  ZERO
      PotVal%d1(2,2,1) =  TWO*Para_Tully%A*Para_Tully%B*r*exp(-Para_Tully%B*r**2)
      PotVal%d1(1,2,1) = -TWO*Para_Tully%C*Para_Tully%D*r*exp(-Para_Tully%D*r**2)
      PotVal%d1(2,1,1) =  PotVal%d1(1,2,1)
    endif

! Hessian calculation
    if (nderiv .eq. 2) then
      PotVal%d2(1,1,1,1) =  ZERO
      PotVal%d2(2,2,1,1) =  TWO*Para_Tully%A*Para_Tully%B*(ONE-TWO*Para_Tully%B*r**2)*exp(-Para_Tully%B*r**2)
      PotVal%d2(1,2,1,1) = -TWO*Para_Tully%C*Para_Tully%D*(ONE-TWO*Para_Tully%D*r**2)*exp(-Para_Tully%D*r**2)
      PotVal%d2(2,1,1,1) =  PotVal%d2(1,2,1,1)
    endif

  END SUBROUTINE eval_TullyPot2_old

  SUBROUTINE eval_TullyPot3_old(PotVal,r,Para_Tully,nderiv)   ! 3rd Tully's potential
    !  C. Extended coupling with reflection
    USE mod_dnMatPot

    TYPE (Param_Tully), intent(in)     :: Para_Tully
    real (kind=Rkind),  intent(in)     :: r
    TYPE(dnMatPot),     intent(inout)  :: PotVal
    integer, intent(in)                :: nderiv

! Potential calculation
    PotVal%d0(1,1) =  Para_Tully%A
    PotVal%d0(2,2) = -Para_Tully%A
    if (r >= ZERO) then
      PotVal%d0(1,2) = Para_Tully%B*(TWO-exp(-Para_Tully%C*r))
    else 
      PotVal%d0(1,2) = Para_Tully%B*exp(Para_Tully%C*r)
    endif
    PotVal%d0(2,1) = PotVal%d0(1,2)

! gradient calculation
    if (nderiv .ge. 1) then
      PotVal%d1(1,1,1) = ZERO
      PotVal%d1(2,2,1) = ZERO
      PotVal%d1(1,2,1) = Para_Tully%B*Para_Tully%C*exp(-Para_Tully%C*abs(r))
      PotVal%d1(2,1,1) = PotVal%d1(1,2,1)
    endif

! Hessian calculation
    if (nderiv .eq. 2) then
        PotVal%d2(1,1,1,1) = ZERO
        PotVal%d2(2,2,1,1) = ZERO
        if (r >= ZERO) then
          PotVal%d2(1,2,1,1) = -Para_Tully%C*PotVal%d1(1,2,1)
        else
          PotVal%d2(1,2,1,1) = Para_Tully%C*PotVal%d1(1,2,1)
        endif 
        PotVal%d2(2,1,1,1) = PotVal%d2(1,2,1,1)
    endif

  END SUBROUTINE eval_TullyPot3_old

END MODULE mod_TullyPot
