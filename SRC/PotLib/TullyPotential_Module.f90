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
  TYPE TullyPot_t
     PRIVATE
     real (kind=Rkind) :: A      = 0.01_Rkind
     real (kind=Rkind) :: B      = 1.6_Rkind
     real (kind=Rkind) :: C      = 0.005_Rkind
     real (kind=Rkind) :: D      = 1.0_Rkind
     real (kind=Rkind) :: E0     = 0.0_Rkind

     real (kind=Rkind), PUBLIC :: mu     = 2000._Rkind


     integer           :: option = 1
     
  END TYPE TullyPot_t

  PRIVATE Read_TullyPot,Init0_TullyPot,eval_TullyPot1,eval_TullyPot2,eval_TullyPot3

CONTAINS
!> @brief Subroutine which makes the initialization of the Tully parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param TullyPot         TYPE(TullyPot_t):  derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
!! @param A,B,C,D,E0         real (optional):    Tully parameters.
  SUBROUTINE Init_TullyPot(TullyPot,option,nio,read_param,A,B,C,D,E0)
    TYPE (TullyPot_t),          intent(inout)   :: TullyPot
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

    TullyPot%option = option

    IF (TullyPot%option < 1 .OR. TullyPot%option > 3) TullyPot%option = 1

    !write(out_unitp,*) 'option',option
    !write(out_unitp,*) 'read_param',read_param,nio


    IF (read_param_loc) THEN
      CALL Read_TullyPot(TullyPot,nio)
      CALL Init0_TullyPot(TullyPot)
    ELSE

      TullyPot%A      = Huge(ONE)
      TullyPot%B      = Huge(ONE)
      TullyPot%C      = Huge(ONE)
      TullyPot%D      = Huge(ONE)
      TullyPot%E0     = Huge(ONE)

      CALL Init0_TullyPot(TullyPot)

      SELECT CASE (TullyPot%option)
      CASE (1) ! A. Simple avoided crossing

        IF (present(A)) TullyPot%A = A
        IF (present(B)) TullyPot%B = B
        IF (present(C)) TullyPot%C = C
        IF (present(D)) TullyPot%D = D

      CASE (2) ! B. Dual avoided crossing

        IF (present(A))  TullyPot%A  = A
        IF (present(B))  TullyPot%B  = B
        IF (present(E0)) TullyPot%E0 = E0
        IF (present(C))  TullyPot%C  = C
        IF (present(D))  TullyPot%D  = D

      CASE (3) !C. Extended coupling with reflection

        IF (present(A)) TullyPot%A = A
        IF (present(B)) TullyPot%B = B
        IF (present(C)) TullyPot%C = C

      CASE Default
          write(out_unitp,*) 'ERROR in Init_TullyPot'
          write(out_unitp,*) ' This option is not possible. option:',TullyPot%option
          write(out_unitp,*) ' Its value MUST be 1 or 2 or 3'
          STOP
      END SELECT
    END IF


  END SUBROUTINE Init_TullyPot
!> @brief Subroutine which makes the initialization of the Tully parameters from his paper.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param TullyPot         TYPE(TullyPot_t):  derived type in which the parameters are set-up.
  SUBROUTINE Init0_TullyPot(TullyPot)
    TYPE (TullyPot_t),      intent(inout)   :: TullyPot

    character (len=*), parameter :: name_sub='Init0_TullyPot'

    !write(out_unitp,*) 'BEGINNING ',name_sub
    !write(out_unitp,*) 'TullyPot%option',TullyPot%option

    IF (TullyPot%option < 1 .OR. TullyPot%option > 3) TullyPot%option = 1

    SELECT CASE (TullyPot%option)
    CASE (1) ! A. Simple avoided crossing
      IF (TullyPot%A == huge(ONE)) TullyPot%A = 0.01_Rkind
      IF (TullyPot%B == huge(ONE)) TullyPot%B = 1.6_Rkind
      IF (TullyPot%C == huge(ONE)) TullyPot%C = 0.005_Rkind
      IF (TullyPot%D == huge(ONE)) TullyPot%D = 1.0_Rkind

    CASE (2) ! B. Dual avoided crossing
      IF (TullyPot%A  == huge(ONE)) TullyPot%A  = 0.1_Rkind
      IF (TullyPot%B  == huge(ONE)) TullyPot%B  = 0.28_Rkind
      IF (TullyPot%E0 == huge(ONE)) TullyPot%E0 = 0.05_Rkind
      IF (TullyPot%C  == huge(ONE)) TullyPot%C  = 0.015_Rkind
      IF (TullyPot%D  == huge(ONE)) TullyPot%D  = 0.06_Rkind

    CASE (3) !C. Extended coupling with reflection
      IF (TullyPot%A == huge(ONE)) TullyPot%A = 0.0006_Rkind
      IF (TullyPot%B == huge(ONE)) TullyPot%B = 0.1_Rkind
      IF (TullyPot%C == huge(ONE)) TullyPot%C = 0.90_Rkind

    CASE Default
        write(out_unitp,*) 'ERROR in Init0_TullyPot'
        write(out_unitp,*) ' This option is not possible. option:',TullyPot%option
        write(out_unitp,*) ' Its value MUST be 1 or 2 or 3'
        STOP
    END SELECT

    !CALL Write_TullyPot(TullyPot,out_unitp)
    !write(out_unitp,*) 'END ',name_sub


  END SUBROUTINE Init0_TullyPot

!> @brief Subroutine wich reads the Tully parameters with a namelist.
!!   This can be called only from the "Init_TullyPot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param TullyPot         TYPE(TullyPot_t):   derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to read the parameters.
  SUBROUTINE Read_TullyPot(TullyPot,nio)
    TYPE (TullyPot_t),      intent(inout)   :: TullyPot
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

    TullyPot%A       = A
    TullyPot%B       = B
    TullyPot%E0      = E0
    TullyPot%C       = C
    TullyPot%D       = D

  END SUBROUTINE Read_TullyPot

!> @brief Subroutine wich prints the Tully parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write0_TullyPot(nio)
    integer, intent(in) :: nio

    write(nio,*) 'Tully parameters, from reference:'
    write(nio,*) '  Reference: Tully, J. Chem. Phys. V93, pp15, 1990'
    write(nio,*)
    write(nio,*) '  mu     = 2000 au'
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, Tully A (option=1)'
    write(nio,*)
    write(nio,*) 'V11(R) =  A.[1-Exp(-B.R)]    R>0'
    write(nio,*) 'V11(R) = -A.[1-Exp( B.R)]    R<0'
    write(nio,*) 'V22(R) = -V11(R)'
    write(nio,*) 'V12(R) = C.Exp(-D.R^2)'
    write(nio,*) 'A = 0.01, B = 1.6, C = 0.005, D = 1.0'
    write(nio,*)
    write(nio,*) 'Value at: R=-1.1 Bohr'
    write(nio,*) 'V        = [-0.008413,  0.008413] Hartree'
    write(nio,*) 'gradient = [ 0.002128, -0.002128] '
    write(nio,*) 'hessian  = [ 0.001943, -0.001943] '
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, Tully B (option=2)'
    write(nio,*)
    write(nio,*) 'V11(R) =  0'
    write(nio,*) 'V22(R) = -A.Exp(-B.R^2) + E0'
    write(nio,*) 'V12(R) = C.Exp(-D.R^2)'
    write(nio,*) 'A = 0.10, B = 0.28, E0=0.05, C = 0.015, D = 0.06'
    write(nio,*)
    write(nio,*) 'Value at: R=-1.1 Bohr'
    write(nio,*) 'V        = [-0.028170,  0.006908] Hartree'
    write(nio,*) 'gradient = [-0.036718, -0.007180] '
    write(nio,*) 'hessian  = [-0.003754,  0.016620] '
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, Tully C (option=3)'
    write(nio,*)
    write(nio,*) 'V11(R) =  A'
    write(nio,*) 'V22(R) = -A'
    write(nio,*) 'V12(R) = B.   Exp( C.R)    R<0'
    write(nio,*) 'V12(R) = B.[2-Exp(-C.R)]   R>0'
    write(nio,*) 'A = 0.0006, B = 0.1, C = 0.90'
    write(nio,*)
    write(nio,*) 'Value at: R=-1.1 Bohr'
    write(nio,*) 'V        = [-0.037163,  0.037163] Hartree'
    write(nio,*) 'gradient = [-0.033438,  0.033438] '
    write(nio,*) 'hessian  = [-0.030102,  0.030102] '
    write(nio,*)
    write(nio,*) 'end Tully parameters'

  END SUBROUTINE Write0_TullyPot
!> @brief Subroutine wich prints the Tully parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param TullyPot         TYPE(TullyPot_t):   derived type in which the parameters are set-up.
!! @param nio                integer:             file unit to print the parameters.
  SUBROUTINE Write_TullyPot(TullyPot,nio)
    TYPE (TullyPot_t), intent(in) :: TullyPot
    integer, intent(in) :: nio

    write(nio,*) 'Tully current parameters:'
    write(nio,*)
    write(nio,*) '  A:      ',TullyPot%A
    write(nio,*) '  B:      ',TullyPot%B
    write(nio,*) '  C:      ',TullyPot%C
    write(nio,*) '  D:      ',TullyPot%D
    write(nio,*) '  E0:     ',TullyPot%E0
    write(nio,*) '  option: ',TullyPot%option
    write(nio,*)
    write(nio,*) 'end Tully parameters'

  END SUBROUTINE Write_TullyPot
  SUBROUTINE get_Q0_Tully(R0,TullyPot)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: R0
    TYPE (TullyPot_t),          intent(in)    :: TullyPot

    R0 = ZERO

  END SUBROUTINE get_Q0_Tully
!> @brief Subroutine wich calculates the Tully potential (for the 3 models) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param TullyPot         TYPE(TullyPot_t):   derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TullyPot(Mat_OF_PotDia,dnR,TullyPot,nderiv)
    USE mod_dnS

    TYPE (TullyPot_t), intent(in)     :: TullyPot
    TYPE (dnS_t),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)     :: dnR
    integer, intent(in)                :: nderiv

    SELECT CASE (TullyPot%option)
    CASE (1)
      CALL eval_TullyPot1(Mat_OF_PotDia,dnR,TullyPot,nderiv)
    CASE (2)
      CALL eval_TullyPot2(Mat_OF_PotDia,dnR,TullyPot,nderiv)
    CASE (3)
      CALL eval_TullyPot3(Mat_OF_PotDia,dnR,TullyPot,nderiv)
    CASE Default
        write(out_unitp,*) 'ERROR in eval_TullyPot'
        write(out_unitp,*) ' This option is not possible. option:',TullyPot%option
        write(out_unitp,*) ' Its value MUST be 1 or 3 or 3'
        STOP
    END SELECT

  END SUBROUTINE eval_TullyPot
!> @brief Subroutine wich calculates the Tully potential (for the 1st model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param TullyPot         TYPE(TullyPot_t):   derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TullyPot1(Mat_OF_PotDia,dnR,TullyPot,nderiv)
  ! A. Simple avoided crossing
    USE mod_dnS

    TYPE (TullyPot_t), intent(in)     :: TullyPot
    TYPE (dnS_t),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)     :: dnR
    integer, intent(in)                :: nderiv


    ! Potential calculation
    IF (dnR >= ZERO) THEN
      Mat_OF_PotDia(1,1) =  TullyPot%A*(ONE-exp(-TullyPot%B*dnR))
    ELSE
      Mat_OF_PotDia(1,1) = -TullyPot%A*(ONE-exp( TullyPot%B*dnR))
    END IF
    Mat_OF_PotDia(2,2)   = -Mat_OF_PotDia(1,1)


    Mat_OF_PotDia(1,2) = TullyPot%C*exp(-TullyPot%D*dnR**2)
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


  END SUBROUTINE eval_TullyPot1
!> @brief Subroutine wich calculates the Tully potential (for the 2d model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param TullyPot         TYPE(TullyPot_t):   derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TullyPot2(Mat_OF_PotDia,dnR,TullyPot,nderiv) !2d Tully's potential
  !  B. Dual avoided crossing
    USE mod_dnS

    TYPE (TullyPot_t), intent(in)     :: TullyPot
    TYPE (dnS_t),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)     :: dnR
    integer, intent(in)                :: nderiv

! Potential calculation
    Mat_OF_PotDia(1,1)  = ZERO ! to initialized the PotVal%d0(1,1) and its derivatives

    Mat_OF_PotDia(2,2)  = -TullyPot%A*exp(-TullyPot%B*dnR**2)+TullyPot%E0

    Mat_OF_PotDia(1,2) = TullyPot%C*exp(-TullyPot%D*dnR**2)
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


  END SUBROUTINE eval_TullyPot2
!> @brief Subroutine wich calculates the Tully potential (for the 3d model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param TullyPot         TYPE(TullyPot_t):   derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TullyPot3(Mat_OF_PotDia,dnR,TullyPot,nderiv) !3d Tully's potential
    !  C. Extended coupling with reflection
    USE mod_dnS

    TYPE (TullyPot_t), intent(in)     :: TullyPot
    TYPE (dnS_t),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)     :: dnR
    integer, intent(in)                :: nderiv


! Potential calculation
    Mat_OF_PotDia(1,1) =  TullyPot%A
    Mat_OF_PotDia(2,2) = -TullyPot%A

    IF (dnR >= ZERO) THEN
      Mat_OF_PotDia(1,2) = TullyPot%B*(TWO-exp(-TullyPot%C*dnR))
    ELSE
      Mat_OF_PotDia(1,2) = TullyPot%B*exp(TullyPot%C*dnR)
    END IF
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


  END SUBROUTINE eval_TullyPot3

END MODULE mod_TullyPot
