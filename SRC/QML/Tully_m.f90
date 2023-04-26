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

!> @brief Module which makes the initialization, calculation of the Tully potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_Tully_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Tully parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_Tully_t
   PRIVATE

     real (kind=Rkind) :: A      = 0.01_Rkind
     real (kind=Rkind) :: B      = 1.6_Rkind
     real (kind=Rkind) :: C      = 0.005_Rkind
     real (kind=Rkind) :: D      = 1.0_Rkind
     real (kind=Rkind) :: E0     = 0.0_Rkind

     real (kind=Rkind) :: mu     = 2000._Rkind


   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_Tully
    PROCEDURE :: Write_QModel    => Write_QML_Tully
  END TYPE QML_Tully_t

  PUBLIC :: QML_Tully_t,Init_QML_Tully

  CONTAINS
!> @brief Function which makes the initialization of the Tully parameters.
!!
!! @param QModel             TYPE(QML_Tully_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_Tully(QModel_in,read_param,nio_param_file,A,B,C,D,E0) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_Tully_t)                           :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: A,B,C,D,E0


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Tully'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf    = 2
    QModel%ndim     = 1
    QModel%pot_name = 'tully'

    IF (QModel%option < 1 .OR. QModel%option > 3) QModel%option = 1

    ! Default values as function of option
    SELECT CASE (QModel%option)
    CASE (1) ! A. Simple avoided crossing

      QModel%A = 0.01_Rkind
      QModel%B = 1.6_Rkind
      QModel%C = 0.005_Rkind
      QModel%D = 1.0_Rkind

    CASE (2) ! B. Dual avoided crossing

      QModel%A  = 0.1_Rkind
      QModel%B  = 0.28_Rkind
      QModel%E0 = 0.05_Rkind
      QModel%C  = 0.015_Rkind
      QModel%D  = 0.06_Rkind

    CASE (3) !C. Extended coupling with reflection

      QModel%A = 0.0006_Rkind
      QModel%B = 0.1_Rkind
      QModel%C = 0.90_Rkind

    CASE Default
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) ' This option is not possible. option:',QModel%option
      write(out_unit,*) ' Its value MUST be 1 or 2 or 3'
      STOP
    END SELECT


    IF (read_param) THEN
      CALL Read_QML_Tully(Asub=QModel%A,Bsub=QModel%B,Csub=QModel%C,   &
                       Dsub=QModel%D,E0sub=QModel%E0,nio=nio_param_file)
    ELSE

      SELECT CASE (QModel%option)
      CASE (1) ! A. Simple avoided crossing

        IF (present(A)) QModel%A = A
        IF (present(B)) QModel%B = B
        IF (present(C)) QModel%C = C
        IF (present(D)) QModel%D = D

      CASE (2) ! B. Dual avoided crossing

        IF (present(A))  QModel%A  = A
        IF (present(B))  QModel%B  = B
        IF (present(E0)) QModel%E0 = E0
        IF (present(C))  QModel%C  = C
        IF (present(D))  QModel%D  = D

      CASE (3) !C. Extended coupling with reflection

        IF (present(A)) QModel%A = A
        IF (present(B)) QModel%B = B
        IF (present(C)) QModel%C = C

      CASE Default
          write(out_unit,*) 'ERROR in ',name_sub
          write(out_unit,*) ' This option is not possible. option:',QModel%option
          write(out_unit,*) ' Its value MUST be 1 or 2 or 3'
          STOP
      END SELECT
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of Tully'
    QModel%Q0 = [ZERO]

    IF (debug) write(out_unit,*) 'init d0GGdef of Tully'
    QModel%d0GGdef      = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%mu

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Tully
!> @brief Subroutine wich reads the Tully parameters with a namelist.
!!   This can be called only from the "Init_QML_Tully" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param nio                       integer:   file unit to read the parameters.
!! @param Asub,Bsub,Csub,Dsub,E0sub real:      model parameters.
  SUBROUTINE Read_QML_Tully(Asub,Bsub,Csub,Dsub,E0sub,nio)
    IMPLICIT NONE

    real (kind=Rkind),    intent(inout) :: Asub,Bsub,Csub,Dsub,E0sub
    integer,              intent(in)    :: nio


    real (kind=Rkind)      :: A,B,C,D,E0
    integer                :: err_read

    namelist /Tully/ A,B,C,D,E0

    A      = Asub
    B      = Bsub
    C      = Csub
    D      = Dsub
    E0     = E0sub

    read(nio,Tully,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Tully'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "Tully" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_QML_Tully'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Tully'
      write(out_unit,*) ' Some parameter names of the namelist "Tully" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=Tully)
      STOP ' ERROR in Read_QML_Tully'
    END IF

    Asub       = A
    Bsub       = B
    E0sub      = E0
    Csub       = C
    Dsub       = D

  END SUBROUTINE Read_QML_Tully
!> @brief Subroutine wich prints the current QML_Tully parameters.
!!
!! @param QModel            CLASS(QML_Tully_t):  derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_Tully(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_Tully_t),  intent(in) :: QModel
    integer,              intent(in) :: nio
    write(nio,*) 'Tully default parameters, from reference:'
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
    write(nio,*) 'end Tully default parameters'
    write(nio,*) 'Tully current parameters:'
    write(nio,*)
    write(nio,*) '  A:      ',QModel%A
    write(nio,*) '  B:      ',QModel%B
    write(nio,*) '  C:      ',QModel%C
    write(nio,*) '  D:      ',QModel%D
    write(nio,*) '  E0:     ',QModel%E0
    write(nio,*) '  option: ',QModel%option
    write(nio,*)
    write(nio,*) 'end Tully parameters'

  END SUBROUTINE Write_QML_Tully
!> @brief Subroutine wich calculates the Tully potential (for the 3 models) with derivatives up to the 2d order is required.
!!
!! @param QModel             CLASS(QML_Tully_t):  derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Tully(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_Tully_t),  intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    SELECT CASE (QModel%option)
    CASE (1)
      CALL EvalPot1_QML_Tully(Mat_OF_PotDia,dnQ(1),QModel,nderiv)
    CASE (2)
      CALL EvalPot2_QML_Tully(Mat_OF_PotDia,dnQ(1),QModel,nderiv)
    CASE (3)
      CALL EvalPot3_QML_Tully(Mat_OF_PotDia,dnQ(1),QModel,nderiv)
    CASE Default
        write(out_unit,*) 'ERROR in EvalPot_QML_Tully'
        write(out_unit,*) ' This option is not possible. option:',QModel%option
        write(out_unit,*) ' Its value MUST be 1 or 3 or 3'
        STOP
    END SELECT

  END SUBROUTINE EvalPot_QML_Tully
!> @brief Subroutine wich calculates the Tully potential (for the 1st model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_Tully_t):  derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot1_QML_Tully(Mat_OF_PotDia,dnR,QModel,nderiv)
  !A. Simple avoided crossing
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (QML_Tully_t), intent(in)     :: QModel
    TYPE (dnS_t),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)     :: dnR
    integer,             intent(in)     :: nderiv


    ! Potential calculation
    IF (dnR >= ZERO) THEN
      Mat_OF_PotDia(1,1) =  QModel%A*(ONE-exp(-QModel%B*dnR))
    ELSE
      Mat_OF_PotDia(1,1) = -QModel%A*(ONE-exp( QModel%B*dnR))
    END IF
    Mat_OF_PotDia(2,2)   = -Mat_OF_PotDia(1,1)


    Mat_OF_PotDia(1,2) = QModel%C*exp(-QModel%D*dnR**2)
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


  END SUBROUTINE EvalPot1_QML_Tully
!> @brief Subroutine wich calculates the Tully potential (for the 2d model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_Tully_t):  derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot2_QML_Tully(Mat_OF_PotDia,dnR,QModel,nderiv) !2d Tully's potential
  !B. Dual avoided crossing
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (QML_Tully_t), intent(in)     :: QModel
    TYPE (dnS_t),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)     :: dnR
    integer,             intent(in)     :: nderiv

! Potential calculation
    Mat_OF_PotDia(1,1)  = ZERO ! to initialized the PotVal%d0(1,1) and its derivatives

    Mat_OF_PotDia(2,2)  = -QModel%A*exp(-QModel%B*dnR**2)+QModel%E0

    Mat_OF_PotDia(1,2) = QModel%C*exp(-QModel%D*dnR**2)
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


  END SUBROUTINE EvalPot2_QML_Tully
!> @brief Subroutine wich calculates the Tully potential (for the 3d model) with derivatives up to the 2d order is required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_Tully_t):  derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot3_QML_Tully(Mat_OF_PotDia,dnR,QModel,nderiv) !3d Tully's potential
  !C. Extended coupling with reflection
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (QML_Tully_t), intent(in)     :: QModel
    TYPE (dnS_t),        intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)     :: dnR
    integer,             intent(in)     :: nderiv


! Potential calculation
    Mat_OF_PotDia(1,1) =  QModel%A
    Mat_OF_PotDia(2,2) = -QModel%A

    IF (dnR >= ZERO) THEN
      Mat_OF_PotDia(1,2) = QModel%B*(TWO-exp(-QModel%C*dnR))
    ELSE
      Mat_OF_PotDia(1,2) = QModel%B*exp(QModel%C*dnR)
    END IF
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


  END SUBROUTINE EvalPot3_QML_Tully

END MODULE QML_Tully_m
