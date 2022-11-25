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
!> @brief Module which makes the initialization, calculation of the TwoD_RJDI2014 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_TwoD_RJDI2014_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the TwoD_RJDI2014 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_TwoD_RJDI2014_t
   PRIVATE

   real(kind=Rkind)     :: w1    = 9.557e-3_Rkind
   real(kind=Rkind)     :: w2    = 3.3515e-3_Rkind
   real(kind=Rkind)     :: Delta = 0.01984_Rkind
   real(kind=Rkind)     :: a     = 20.07_Rkind
   real(kind=Rkind)     :: c     = 6.127e-4_Rkind

   real (kind=Rkind)    :: muX  = 1._Rkind
   real (kind=Rkind)    :: muY  = 1._Rkind

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_TwoD_RJDI2014
    PROCEDURE :: Write_QModel    => Write_QML_TwoD_RJDI2014
    PROCEDURE :: Write0_QModel   => Write0_QML_TwoD_RJDI2014
  END TYPE QML_TwoD_RJDI2014_t

  PUBLIC :: QML_TwoD_RJDI2014_t,Init_QML_TwoD_RJDI2014

  CONTAINS
!> @brief Subroutine which makes the initialization of the TwoD_RJDI2014 parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 12/07/2019
!!
!! @param QModel             TYPE(QML_TwoD_RJDI2014_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_TwoD_RJDI2014(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_TwoD_RJDI2014_t)                          :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_TwoD_RJDI2014'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 2
    QModel%ndim     = 2
    QModel%pot_name = 'TwoD_RJDI2014'


    IF (debug) write(out_unitp,*) 'init Q0 of TwoD_RJDI2014'
    SELECT CASE (QModel%option)
    CASE (1) ! minimum of V(1,1)
      QModel%Q0 = [-QModel%a/TWO,ZERO]
    CASE (2)  ! minimum of V(2,2)
      QModel%Q0 = [QModel%a/TWO,ZERO]
    CASE (0)  ! minimum of V(1,1), other parameters
      QModel%w1    = 0.007743_Rkind
      QModel%w2    = 0.00668_Rkind
      QModel%Delta = ZERO
      QModel%a     = 51.05_Rkind
      QModel%c     = 0.00008092_Rkind

      QModel%Q0    = [-QModel%a/TWO,ZERO]
    CASE Default ! ! minimum of V(1,1), other parameters
      QModel%w1    = 0.007743_Rkind
      QModel%w2    = 0.00668_Rkind
      QModel%Delta = ZERO
      QModel%a     = 51.05_Rkind
      QModel%c     = 0.00008092_Rkind

      QModel%Q0    = [-QModel%a/TWO,ZERO]

    END SELECT

    IF (debug) write(out_unitp,*) 'init d0GGdef of TwoD_RJDI2014'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%muX
    QModel%d0GGdef(2,2) = ONE/QModel%muY

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_TwoD_RJDI2014
!> @brief Subroutine wich prints the current QML_TwoD_RJDI2014 parameters.
!!
!! @param QModel            CLASS(QML_TwoD_RJDI2014_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_TwoD_RJDI2014(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_TwoD_RJDI2014_t),  intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'TwoD_RJDI2014 current parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Published model from: '
    write(nio,*) ' Ilya G. Ryabinkin, Loïc Joubert-Doriol, and Artur F. Izmaylov, ...'
    write(nio,*) '  .... J. Chem. Phys. 140, 214116 (2014); https://doi.org/10.1063/1.4881147'
    write(nio,*) 'with the X=Q(1), Y=Q(2) in bohr.'
    write(nio,*) '     and the energy in Hartree.'
    write(nio,*) 'PubliUnit: ',QModel%PubliUnit
    write(nio,*)

    SELECT CASE (QModel%option)
    CASE (1,2)
      write(nio,*) 'Current parameters (from the publication):'
    CASE (0)
      write(nio,*) 'Current parameters (unpublished):'
    CASE DEFAULT
      write(nio,*) 'Current parameters (unpublished):'
    END SELECT

    write(nio,*) '-----------------------------------------'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'w1 (X) =',QModel%w1
    write(nio,*) 'w2 (Y) =',QModel%w2
    write(nio,*) 'DELTA  =',QModel%DELTA
    write(nio,*) 'a      =',QModel%a
    write(nio,*) 'c      =',QModel%c
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'MuX    =',QModel%MuX
    write(nio,*) 'MuY    =',QModel%MuY
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Q0     =',QModel%Q0
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'end TwoD_RJDI2014 parameters'

  END SUBROUTINE Write_QML_TwoD_RJDI2014
!> @brief Subroutine wich prints the default QML_TwoD_RJDI2014 parameters.
!!
!! @param QModel            CLASS(QML_TwoD_RJDI2014_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_TwoD_RJDI2014(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_TwoD_RJDI2014_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'TwoD_RJDI2014 default parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Published model from: '
    write(nio,*) ' Ilya G. Ryabinkin, Loïc Joubert-Doriol, and Artur F. Izmaylov, ...'
    write(nio,*) '  .... J. Chem. Phys. 140, 214116 (2014); https://doi.org/10.1063/1.4881147'
    write(nio,*) 'with the X=Q(1), Y=Q(2) in bohr.'
    write(nio,*) '     and the energy in Hartree.'
    write(nio,*) 'end TwoD_RJDI2014 default parameters'

  END SUBROUTINE Write0_QML_TwoD_RJDI2014

!> @brief Subroutine wich calculates the TwoD_RJDI2014 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_TwoD_RJDI2014_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_TwoD_RJDI2014(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_TwoD_RJDI2014_t),  intent(in)    :: QModel
    TYPE (dnS_t),                intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                intent(in)    :: dnQ(:)
    integer,                     intent(in)    :: nderiv

   !Hel(1,1,x,y)=0.5*w1**2 * (X+a/2)**2 + 0.5d0*w2**2 * Y**2 + Delta/2
   Mat_OF_PotDia(1,1) = HALF*( (QModel%w1 * (dnQ(1)+QModel%a/2))**2 + (QModel%w2 * dnQ(2))**2) + QModel%Delta/2

   !Hel(2,2,x,y)=0.5*w1**2 * (X-a/2)**2 + 0.5d0*w2**2 * Y**2 - Delta/2
   Mat_OF_PotDia(2,2) = HALF*( (QModel%w1 * (dnQ(1)-QModel%a/2))**2 + (QModel%w2 * dnQ(2))**2) - QModel%Delta/2

   !Hel(1,2,x,y) = c * Y
   Mat_OF_PotDia(1,2) = QModel%c * dnQ(2)
   Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)

  END SUBROUTINE EvalPot_QML_TwoD_RJDI2014

END MODULE QML_TwoD_RJDI2014_m
