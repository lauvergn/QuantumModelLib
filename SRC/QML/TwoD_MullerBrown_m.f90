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

!> @brief Module which makes the initialization, calculation of the TwoD_MullerBrown potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_TwoD_MullerBrown_m
  USE QDUtil_NumParameters_m, out_unitp => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the TwoD_MullerBrown parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_TwoD_MullerBrown_t
   PRIVATE

   real(kind=Rkind)     :: AA(4) = [-200._Rkind,-100._Rkind,-170._Rkind,15._Rkind]
   real(kind=Rkind)     :: a(4)  = [-ONE,-ONE,-6.5_Rkind,0.7_Rkind]
   real(kind=Rkind)     :: b(4)  = [ZERO,ZERO,11._Rkind,0.6_Rkind]
   real(kind=Rkind)     :: c(4)  = [-TEN,-TEN,-6.5_Rkind,0.7_Rkind]
   real(kind=Rkind)     :: x0(4) = [ONE,ZERO,-HALF,-ONE]
   real(kind=Rkind)     :: y0(4) = [ZERO,HALF,1.5_Rkind,ONE]

   ! the three minima (A,B,C) and the two TSs (from the gradient c) from table1
   real(kind=Rkind)     :: tab_Q0(2,5) = reshape(                               &
                                                [-0.558_Rkind,1.442_Rkind,      &
                                                  0.623_Rkind,0.028_Rkind,      &
                                                 -0.050_Rkind,0.467_Rkind,      &
                                                 -0.822_Rkind,0.624_Rkind,      &
                                                 -0.212_Rkind,0.293_Rkind],     &
                                                 shape=[2,5])
   real(kind=Rkind)     :: tab_Ene(5) = [-146.700_Rkind,                        &
                                         -108.167_Rkind,                        &
                                          -80.768_Rkind,                        &
                                          -40.665_Rkind,                        &
                                          -72.249_Rkind]


   real (kind=Rkind)    :: muX  = 1000._Rkind
   real (kind=Rkind)    :: muY  = 1000._Rkind

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_TwoD_MullerBrown
    PROCEDURE :: Write_QModel    => Write_QML_TwoD_MullerBrown
    PROCEDURE :: Write0_QModel   => Write0_QML_TwoD_MullerBrown
  END TYPE QML_TwoD_MullerBrown_t

  PUBLIC :: QML_TwoD_MullerBrown_t,Init_QML_TwoD_MullerBrown

  CONTAINS
!> @brief Subroutine which makes the initialization of the TwoD_MullerBrown parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 12/07/2019
!!
!! @param QModel             TYPE(QML_TwoD_MullerBrown_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_TwoD_MullerBrown(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_TwoD_MullerBrown_t)                :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_TwoD_MullerBrown'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 1
    QModel%ndim     = 2
    QModel%pot_name = 'twod_mullerbrown'

    IF (QModel%option < 1 .OR. QModel%option > 5) QModel%option = 1

    IF (debug) write(out_unitp,*) 'init Q0 of TwoD_MullerBrown'
    QModel%Q0 = QModel%tab_Q0(:,QModel%option)

    IF (debug) write(out_unitp,*) 'init d0GGdef of TwoD_MullerBrown'
    QModel%d0GGdef      = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%muX
    QModel%d0GGdef(2,2) = ONE/QModel%muY

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_TwoD_MullerBrown
!> @brief Subroutine wich prints the current QML_TwoD_MullerBrown parameters.
!!
!! @param QModel            CLASS(QML_TwoD_MullerBrown_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_TwoD_MullerBrown(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_TwoD_MullerBrown_t),  intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'TwoD_MullerBrown current parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'The parameters are different from the published ones: '
    write(nio,*) ' Klaus Müller and Leo D. Brown, ...'
    write(nio,*) '  ....Theoret. Chim. Acta (Berl.) 53, 75-93 (1979; https://doi.org/10.1007/BF00547608'
    write(nio,*) 'no unit, we assume atomic units'


    write(nio,*) 'PubliUnit: ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) ' Potential with published parameters.'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'MuX =',QModel%MuX
    write(nio,*) 'MuY =',QModel%MuY
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Q0     =',QModel%Q0
    write(nio,*) 'E(Q0)  =',QModel%tab_Ene(QModel%option)
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'end TwoD_MullerBrown parameters'

  END SUBROUTINE Write_QML_TwoD_MullerBrown
!> @brief Subroutine wich prints the default QML_TwoD_MullerBrown parameters.
!!
!! @param QModel            CLASS(QML_TwoD_MullerBrown_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_TwoD_MullerBrown(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_TwoD_MullerBrown_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'TwoD_MullerBrown default parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'The parameters are different from the published ones: '
    write(nio,*) ' Klaus Müller and Leo D. Brown, ...'
    write(nio,*) '  ....Theoret. Chim. Acta (Berl.) 53, 75-93 (1979; https://doi.org/10.1007/BF00547608'
    write(nio,*) 'no unit, we assume atomic units'


    write(nio,*) 'PubliUnit: ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) ' Potential with published parameters.'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'MuX =',QModel%MuX
    write(nio,*) 'MuY =',QModel%MuY
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'Q0     =',QModel%Q0
    write(nio,*) 'E(Q0)  =',QModel%tab_Ene(QModel%option)
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'end TwoD_MullerBrown default parameters'

  END SUBROUTINE Write0_QML_TwoD_MullerBrown

!> @brief Subroutine wich calculates the TwoD_MullerBrown potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_TwoD_MullerBrown_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_TwoD_MullerBrown(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_TwoD_MullerBrown_t),  intent(in)    :: QModel
    TYPE (dnS_t),                   intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                   intent(in)    :: dnQ(:)
    integer,                        intent(in)    :: nderiv

    !local variables
    TYPE (dnS_t)     :: dnDX,dnDY
    integer :: i

    real (kind=Rkind),   parameter     :: EnergyConv = 627.509_Rkind


    i = 1
    dnDX = dnQ(1)-QModel%x0(i)
    dnDY = dnQ(2)-QModel%y0(i)

    Mat_OF_PotDia(1,1) = QModel%AA(i) *                                         &
          exp(QModel%a(i)*dnDX**2 + QModel%b(i)*dnDX*dnDY + QModel%c(i)*dnDY**2)


    DO i=2,size(QModel%AA)
      dnDX = dnQ(1)-QModel%x0(i)
      dnDY = dnQ(2)-QModel%y0(i)
      Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) +   QModel%AA(i) *                &
          exp(QModel%a(i)*dnDX**2 + QModel%b(i)*dnDX*dnDY + QModel%c(i)*dnDY**2)
    END DO

   IF(.NOT. QModel%PubliUnit) Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1)/EnergyConv


    CALL dealloc_dnS(dnDX)
    CALL dealloc_dnS(dnDY)

  END SUBROUTINE EvalPot_QML_TwoD_MullerBrown

END MODULE QML_TwoD_MullerBrown_m
