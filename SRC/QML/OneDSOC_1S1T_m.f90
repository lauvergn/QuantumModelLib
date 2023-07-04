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
!> @brief Module which makes the initialization, calculation of the OneDSOC_1S1T potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_OneDSOC_1S1T_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the OneDSOC_1S1T parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_OneDSOC_1S1T_t
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

     real (kind=Rkind), PUBLIC :: mu  = 20000._Rkind !< Reduced mass from Granucci et al. paper (in au)

   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_OneDSOC_1S1T
    PROCEDURE :: Write_QModel    => Write_QML_OneDSOC_1S1T
  END TYPE QML_OneDSOC_1S1T_t

  PUBLIC :: QML_OneDSOC_1S1T_t,Init_QML_OneDSOC_1S1T

  CONTAINS
!> @brief Function which makes the initialization of the OneDSOC_1S1T parameters.
!!
!! @param QModel             TYPE(QML_OneDSOC_1S1T_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_OneDSOC_1S1T(QModel_in,read_param,nio_param_file,&
                                   Rsig_in) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_OneDSOC_1S1T_t)                  :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: Rsig_in

    !local variable
    real (kind=Rkind)      :: Rsig
    integer                :: err_read

    namelist /OneD_SOC_Model/ Rsig


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_OneDSOC_1S1T'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      CALL QModel_in%Write_QModel(out_unit)
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    IF (QModel%nsurf == 0) THEN
      ! it means that nsurf was not present in the CALL Init_Model().
      ! => nsurf is set to the default value (4)
      QModel%nsurf  = 4
    END IF
    QModel%ndim     = 1
    QModel%pot_name = '1dsoc_1s1t'


    !The value of QModel%nsurf must be 2 or 4
    IF (QModel%nsurf /= 2 .AND. QModel%nsurf /= 4) THEN
       write(out_unit,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unit)
       write(out_unit,*) ' ERROR in ',name_sub
       write(out_unit,*) ' nsurf MUST equal to 4 or 2. nusrf: ',QModel%nsurf
       STOP 'ERROR in Init_QML_OneDSOC_1S1T: nsurf MUST equal to 4 or 2'
    END IF

    ! to be able to change the sigmoid function (default 1, the original one)
    IF (QModel%option < 1 .OR. QModel%option > 2) QModel%option = 1


    IF (read_param) THEN
      Rsig   = 8._Rkind
      read(nio_param_file,OneD_SOC_Model,IOSTAT=err_read)

      IF (err_read < 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' End-of-file or End-of-record'
        write(out_unit,*) ' The namelist "OneD_SOC_Model" is probably absent'
        write(out_unit,*) ' check your data!'
        write(out_unit,*)
        STOP ' ERROR in Init_QML_OneDSOC_1S1T'
      ELSE IF (err_read > 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Some parameter names of the namelist "OneD_SOC_Model" are probaly wrong'
        write(out_unit,*) ' check your data!'
        write(out_unit,nml=OneD_SOC_Model)
        STOP ' ERROR in Init_QML_OneDSOC_1S1T'
      END IF

      QModel%Rsig     = Rsig

    ELSE
      IF (present(Rsig_in)) THEN
        QModel%Rsig     = Rsig_in
      ELSE
        QModel%Rsig     = 8._Rkind
      END IF
    END IF


    IF (debug) write(out_unit,*) 'init Q0 of OneDSOC_1S1T'
    QModel%Q0 = [9.5_Rkind]

    IF (debug) write(out_unit,*) 'init d0GGdef of OneDSOC_1S1T'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%mu

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_OneDSOC_1S1T
!> @brief Subroutine wich prints the current QML_OneDSOC_1S1T parameters.
!!
!! @param QModel            CLASS(QML_OneDSOC_1S1T_t): derived type in which the parameters are set-up.
!! @param nio               integer:                     file unit to print the parameters.
  SUBROUTINE Write_QML_OneDSOC_1S1T(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_OneDSOC_1S1T_t), intent(in) :: QModel
    integer,                     intent(in) :: nio


    write(nio,*) 'OneDSOC_1S1T default parameters, from reference:'
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
    write(nio,*) 'Diabatic Potential, QML_OneDSOC_1S1T (with Rsig=8.0)'
    write(nio,*) 'Value at: R=10. Bohr'
    write(nio,*) 'Vdia (Hartree)   = [ 0.041042414      -0.000707107       0.000707107       0.001000000]'
    write(nio,*) '                   [-0.000707107       0.041042499       0.000000000       0.000000000]'
    write(nio,*) '                   [ 0.000707107       0.000000000       0.041042499       0.000000000]'
    write(nio,*) '                   [ 0.001000000       0.000000000       0.000000000       0.041042499]'
    write(nio,*)
    write(nio,*) 'Adiabatic Potential, QML_OneDSOC_1S1T (with Rsig=8.0)'
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
    write(nio,*) 'end OneDSOC_1S1T default parameters'
    write(nio,*) 'OneDSOC_1S1T current parameters:'
    write(nio,*)
    write(nio,*) '  a1,a2:          ',QModel%a1,QModel%a2
    write(nio,*) '  alpha1,alpha2:  ',QModel%alpha1,QModel%alpha2
    write(nio,*) '  DE:             ',QModel%DE

    write(nio,*) '  C0:             ',QModel%C0
    write(nio,*) '  C1=RC1 + i IC1: ',QModel%RC1,QModel%IC1

    write(nio,*) '   Rsig:          ',QModel%Rsig
    write(nio,*) '  DRsig:          ',QModel%DRsig
    write(nio,*)
    write(nio,*) ' option:          ',QModel%option
    IF (QModel%option == 1) THEN
      write(nio,*) ' Original sigmoid function'
    ELSE
      write(nio,*) ' Modified sigmoid function:'
      write(nio,*) ' sigmoid(R) = tanh(-4 * (R-Rsig)/DRsig )'
    END IF
    write(nio,*) '  nsurf:          ',QModel%nsurf
    write(nio,*)
    write(nio,*) 'end OneDSOC_1S1T current parameters'

  END SUBROUTINE Write_QML_OneDSOC_1S1T
!> @brief Subroutine wich calculates the OneDSOC_1S1T potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_OneDSOC_1S1T_t): derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):                derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)                 value for which the potential is calculated
!! @param nderiv             integer:                     it enables to specify up to which derivatives the potential is calculated:
!!                                                        the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_OneDSOC_1S1T(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_OneDSOC_1S1T_t), intent(in)    :: QModel
    TYPE (dnS_t),                intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                intent(in)    :: dnQ(:)
    integer,                     intent(in)    :: nderiv

    !local variables (derived type). They have to be deallocated
    TYPE (dnS_t)      :: dnSig,dnx
    integer          :: i,j,nsurf
    real(kind=Rkind) :: RC01



    ! get nsurf from the matrix: Mat_OF_PotDia
    nsurf = size(Mat_OF_PotDia,dim=1)

    ! Sig(R) calculation
    IF (QModel%option == 1) THEN ! as in the publication
      !dnSig = dnQ(1) ! to have the correct initialization for dnSig = +/- ONE
      IF (dnQ(1) <= QModel%Rsig-HALF*QModel%DRsig) THEN
        dnSig = ONE
      ELSE IF (dnQ(1) >= QModel%Rsig+HALF*QModel%DRsig) THEN
        dnSig = -ONE
      ELSE
        dnx   = (dnQ(1) - QModel%Rsig) / QModel%DRsig
        dnSig = dnx*(FOUR*dnx**2 -THREE)
      END IF
    ELSE
      dnx   = (dnQ(1) - QModel%Rsig) / QModel%DRsig
      dnSig = tanh(-FOUR*dnx)
    END IF

    IF (nsurf == 4) THEN
      !singlet Energy
      Mat_OF_PotDia(1,1) = QModel%a1 * exp(-QModel%alpha1*dnQ(1)) + QModel%DE

      !Triplet Energy
      Mat_OF_PotDia(2,2) = QModel%a2 * exp(-QModel%alpha2*dnQ(1))
      Mat_OF_PotDia(3,3) = Mat_OF_PotDia(2,2)
      Mat_OF_PotDia(4,4) = Mat_OF_PotDia(2,2)


      !Singley-triplet coupling
      Mat_OF_PotDia(1,2) =  sqrt(TWO) * QModel%RC1 * dnSig
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)

      Mat_OF_PotDia(1,3) =  -sqrt(TWO) * QModel%IC1 * dnSig
      Mat_OF_PotDia(3,1) = Mat_OF_PotDia(1,3)

      Mat_OF_PotDia(1,4) =  -QModel%C0 * dnSig
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
      Mat_OF_PotDia(1,1) = QModel%a1 * exp(-QModel%alpha1*dnQ(1)) + QModel%DE

      !Triplet Energy
      Mat_OF_PotDia(2,2) = QModel%a2 * exp(-QModel%alpha2*dnQ(1))

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

      RC01 = sqrt(TWO*(QModel%RC1**2+QModel%IC1**2)+QModel%C0**2)
      Mat_OF_PotDia(1,2) = dnSig * RC01
      Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)
   END IF


    CALL dealloc_dnS(dnx)
    CALL dealloc_dnS(dnSig)


  END SUBROUTINE EvalPot_QML_OneDSOC_1S1T

END MODULE QML_OneDSOC_1S1T_m
