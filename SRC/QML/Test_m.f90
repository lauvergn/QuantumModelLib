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

!> @brief Module which makes the initialization, calculation of the test potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_Test_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE
  ! add potential paramters (as Fortran parameters) here
  real (kind =Rkind), parameter :: a = 0.5_Rkind
  real (kind =Rkind), parameter :: b = -0.5_Rkind

!> @brief Derived type in which the test parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_Test_t

    PRIVATE

    ! or add potential paramters here
    real (kind =Rkind) :: c

  CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_Test
    PROCEDURE :: Write_QModel     => Write_QML_Test
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_QML_Test ! optional
  END TYPE QML_Test_t

  PUBLIC :: QML_Test_t,Init_QML_Test

  CONTAINS
!> @brief Function which makes the initialization of the test parameters.
!!
!! @param QModel             TYPE(QML_Test_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_Test(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m, ONLY : Identity_Mat
  IMPLICIT NONE

    TYPE (QML_Test_t)                            :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Test'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF


    QModel%QML_Empty_t = QModel_in

    QModel%nsurf    = 2
    QModel%ndimCart = 9 ! 3 atoms
    QModel%ndimQ    = 3 ! 3 internal coordinates (coordinates of the potential)
    QModel%pot_name = 'test'

    IF (QModel%Cart_TO_Q) THEN
      QModel%ndim       = QModel%ndimCart
    ELSE
      QModel%ndim       = QModel%ndimQ
    END IF

    ! add the value of c in Qmodel 
    QModel%c = 3._Rkind

    IF (debug) write(out_unit,*) 'init Q0 of test'
    QModel%Q0 = [ZERO,ZERO,ZERO] ! change the values here

    IF (debug) write(out_unit,*) 'init d0GGdef of test'
    QModel%d0GGdef = Identity_Mat(QModel%ndim) ! change the values here


    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Test
!> @brief Subroutine wich prints the current QML_Test parameters.
!!
!! @param QModel            CLASS(QML_Test_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_Test(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Test_t),   intent(in) :: QModel
    integer,             intent(in) :: nio

    ! add somthing you want to print about your model (reference, system ...)
    write(out_unit,*) 'model name: ',QModel%pot_name
    write(out_unit,*) 'a,b',a,b
    write(out_unit,*) 'c',QModel%c

  END SUBROUTINE Write_QML_Test
!> @brief Subroutine wich calculates the test potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_Test_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Test(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  USE QDUtil_m, ONLY : TO_string
  IMPLICIT NONE

    CLASS(QML_Test_t),    intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    TYPE (dnS_t)   :: dnHarmo
    integer        :: i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_Test'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) ' nderiv    ',nderiv
      !write(out_unit,*) ' dnQ(:)%d0 ',get_d0(dnQ)
      DO i=1,size(dnQ)
        CALL Write_dnS(dnQ(i),info='dnQ(' // TO_string(i) // ')')
      END DO
      flush(out_unit)
    END IF

    ! You can also add the potential parameters here
    ! the Type dnS_t deal with automatic differentiation 
    ! If you want to use intermediate variables, they must be of dnS_t type

    dnHarmo = a * (dnQ(1)**2 + dnQ(2)**2 + dnQ(3)**2)
    IF (debug) CALL Write_dnS(dnHarmo,info='dnHarmo') ; flush(out_unit)

    Mat_OF_PotDia(1,1) =  HALF*dnHarmo
    IF (debug)  CALL Write_dnS(Mat_OF_PotDia(1,1),info='Mat_OF_PotDia(1,1)') ; flush(out_unit)

    Mat_OF_PotDia(2,2) =  HALF*dnHarmo + b
    IF (debug) CALL Write_dnS(Mat_OF_PotDia(2,2),info='Mat_OF_PotDia(2,2)') ; flush(out_unit)

    Mat_OF_PotDia(1,2) =  dnQ(1) * QModel%c
    IF (debug) CALL Write_dnS(Mat_OF_PotDia(1,2),info='Mat_OF_PotDia(1,2)') ; flush(out_unit)

    Mat_OF_PotDia(2,1) =  Mat_OF_PotDia(1,2)
    IF (debug) CALL Write_dnS(Mat_OF_PotDia(2,1),info='Mat_OF_PotDia(2,1)') ; flush(out_unit)

    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE EvalPot_QML_Test



  SUBROUTINE Cart_TO_Q_QML_Test(QModel,dnX,dnQ,nderiv)
    USE QDUtil_m, ONLY : TO_String
    USE ADdnSVM_m
    IMPLICIT NONE
  
      CLASS(QML_Test_t),       intent(in)    :: QModel
      TYPE (dnS_t),            intent(in)    :: dnX(:,:)
      TYPE (dnS_t),            intent(inout) :: dnQ(:)
      integer,                 intent(in)    :: nderiv
  
      ! local vectors
      integer                   :: i,j
      real (kind=Rkind)         :: sm1,sm2,sm3
      TYPE (dnS_t)              :: Vec12(3),Vec13(3),Vec23(3)
  
      !----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Cart_TO_Q_QML_Test'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unit,*) 'BEGINNING ',name_sub
        write(out_unit,*) 'size(dnQ)',size(dnQ)
        write(out_unit,*) 'dnQ:'
        DO j=1,size(dnQ,dim=1)
          CALL Write_dnS(dnQ(j),out_unit,info='dnQ('// TO_String(j) // ')')
        END DO
        write(out_unit,*) 'shape dnX',shape(dnX)
        write(out_unit,*) 'dnX'
        DO i=1,size(dnX,dim=2)
        DO j=1,size(dnX,dim=1)
          CALL Write_dnS(dnX(j,i),out_unit)
        END DO
        END DO
        flush(out_unit)
      END IF

      Vec23(:) = dnX(:,3)-dnX(:,2)
      Vec12(:) = dnX(:,2)-dnX(:,1)
      Vec13(:) = dnX(:,3)-dnX(:,1)
  
      IF (debug) THEN
        write(out_unit,*) 'in ',name_sub,' vect done'
        flush(out_unit)
        DO j=1,size(Vec23,dim=1)
          CALL Write_dnS(Vec23(j),out_unit,info='Vec23')
        END DO
        DO j=1,size(Vec12,dim=1)
          CALL Write_dnS(Vec23(j),out_unit,info='Vec12')
        END DO
        DO j=1,size(Vec23,dim=1)
          CALL Write_dnS(Vec13(j),out_unit,info='Vec13')
        END DO
        flush(out_unit)
      END IF
  
      dnQ(1) = sqrt(dot_product(Vec23,Vec23))
      dnQ(2) = sqrt(dot_product(Vec12,Vec12))
      dnQ(3) = dot_product(Vec12,Vec23)/(dnQ(1)*dnQ(2)) ! cos(th) between Vec12 and Vec23
  
      IF (debug) THEN
        CALL Write_dnS(dnQ(1),out_unit,info='dnQ(1)')
        CALL Write_dnS(dnQ(2),out_unit,info='dnQ(2)')
        CALL Write_dnS(dnQ(3),out_unit,info='dnQ(3)')
        write(out_unit,*) 'END ',name_sub
        flush(out_unit)
      END IF
    END SUBROUTINE Cart_TO_Q_QML_Test
END MODULE QML_Test_m
