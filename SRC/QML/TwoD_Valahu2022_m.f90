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
!> @brief Module which makes the initialization, calculation of the TwoD_Valahu2022 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 21/12/2022
!!
MODULE QML_TwoD_Valahu2022_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the TwoD_Valahu2022 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_TwoD_Valahu2022_t
   PRIVATE

   real(kind=Rkind)     :: k     = TWO*PI * 1000._Rkind ! 1000 Hz
   real(kind=Rkind)     :: w     = TWO*PI * 1000._Rkind/1.5_Rkind ! 1000/1.5~667 is the frequency in Hz

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_TwoD_Valahu2022
    PROCEDURE :: Write_QModel    => Write_QML_TwoD_Valahu2022
  END TYPE QML_TwoD_Valahu2022_t

  PUBLIC :: QML_TwoD_Valahu2022_t,Init_QML_TwoD_Valahu2022

  CONTAINS
!> @brief Subroutine which makes the initialization of the TwoD_Valahu2022 parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 21/12/2022
!!
!! @param QModel             TYPE(QML_TwoD_Valahu2022_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_TwoD_Valahu2022(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_TwoD_Valahu2022_t)                 :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    real (kind=Rkind),   parameter     :: EnergyConv = 6.62606896e-34_Rkind/(TWO*pi)/4.3597439408138934e-18_Rkind ! Hz => Hartree

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_TwoD_Valahu2022'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) 'EnergyConv ',EnergyConv
      flush(out_unit)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 2
    QModel%ndim     = 2
    QModel%pot_name = 'TwoD_Valahu2022'

    IF(.NOT. QModel%PubliUnit) THEN 
      QModel%w = QModel%w*EnergyConv
      QModel%k = QModel%k*EnergyConv
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of TwoD_Valahu2022'

    SELECT CASE (QModel%option)
    CASE (1) ! minimum of V(1,1)
      QModel%Q0 = [-QModel%k/QModel%w,ZERO]
    CASE (2)  ! minimum of V(2,2)
      QModel%Q0 = [QModel%k/QModel%w,ZERO]
    CASE Default ! ! minimum of V(1,1)
      QModel%Q0 = [-QModel%k/QModel%w,ZERO]
    END SELECT

    IF (debug) write(out_unit,*) 'init d0GGdef of TwoD_Valahu2022'
    QModel%d0GGdef      = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,1) = QModel%w
    QModel%d0GGdef(2,2) = QModel%w
    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_TwoD_Valahu2022
!> @brief Subroutine wich prints the current QML_TwoD_Valahu2022 parameters.
!!
!! @param QModel            CLASS(QML_TwoD_Valahu2022_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_TwoD_Valahu2022(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_TwoD_Valahu2022_t),  intent(in) :: QModel
    integer,                       intent(in) :: nio

    write(nio,*) 'TwoD_Valahu2022 current parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Published model from: '
    write(nio,*) ' C. H. Valahu, V. C. Olaya-Agudelo, R. J. MacDonell, ...'
    write(nio,*) '... T. Navickas, A. D. Rao, M. J. Millican, J. B. Pérez-Sánchez, ...'
    write(nio,*) '... J. Yuen-Zhou, M. J. Biercuk, C. Hempel, T. R. Tan, and I. Kassal, ...'
    write(nio,*) '  .... XXXXX; '
    write(nio,*) 'with the X=Q(1), Y=Q(2) in bohr.'
    write(nio,*) '     and the energy in Hartree.'
    write(nio,*) 'PubliUnit: ',QModel%PubliUnit
    write(nio,*)

    SELECT CASE (QModel%option)
    CASE (1,2)
      write(nio,*) 'Current parameters (from the publication):'
    CASE DEFAULT
      write(nio,*) 'Current parameters (from the publication):'
    END SELECT

    write(nio,*) '-----------------------------------------'
    IF (QModel%PubliUnit) THEN
      write(nio,*) ' Unit in (2*pi) Hz'
    ELSE
      write(nio,*) ' Unit in Hartree'
    END IF
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'w (X and Y) =',QModel%w
    write(nio,*) 'kappa       =',QModel%k
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Q0     =',QModel%Q0
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'end TwoD_Valahu2022 parameters'

  END SUBROUTINE Write_QML_TwoD_Valahu2022

!> @brief Subroutine wich calculates the TwoD_Valahu2022 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_TwoD_Valahu2022_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_TwoD_Valahu2022(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_TwoD_Valahu2022_t), intent(in)    :: QModel
    TYPE (dnS_t),                 intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                 intent(in)    :: dnQ(:)
    integer,                      intent(in)    :: nderiv

    !Hel(1,1,x,y)=0.5*w * (x**2+y**2) + k*x
    !Hel(2,2,x,y)=0.5*w * (x**2+y**2) - k*x
    Mat_OF_PotDia(1,1) = HALF*QModel%w* (dnQ(1)**2 + dnQ(2)**2) + QModel%k * dnQ(1)
    Mat_OF_PotDia(2,2) = HALF*QModel%w* (dnQ(1)**2 + dnQ(2)**2) - QModel%k * dnQ(1)

    !Hel(1,2,x,y) = k * Y
    Mat_OF_PotDia(1,2) = QModel%k * dnQ(2)
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)


  END SUBROUTINE EvalPot_QML_TwoD_Valahu2022

END MODULE QML_TwoD_Valahu2022_m
