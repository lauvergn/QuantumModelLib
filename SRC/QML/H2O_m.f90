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

!> @brief Module which makes the initialization, calculation of the H2O potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 24/11/2022
!!
MODULE QML_H2O_m
  USE QDUtil_NumParameters_m, out_unitp => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the H2O parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE, EXTENDS (QML_Empty_t) ::  QML_H2O_t

     PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_H2O
    PROCEDURE :: Write_QModel    => Write_QML_H2O
    PROCEDURE :: Write0_QModel   => Write_QML_H2O
  END TYPE QML_H2O_t

  PUBLIC :: QML_H2O_t,Init_QML_H2O

  CONTAINS
!> @brief Subroutine which makes the initialization of the H2O parameters.
!!
!! @param H2OPot           TYPE(QML_H2O_t):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_H2O(QModel_in,read_param,nio_param_file) RESULT(QModel)
    IMPLICIT NONE

    TYPE (QML_H2O_t)                             :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !local variable
    integer                         :: err_read,nio_fit,i,j,k,idum

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H2O'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 1
    QModel%pot_name = 'H2O'
    QModel%ndim     = 3


    IF (QModel%option /= 1) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1)

      QModel%d0GGdef = reshape(                                                  &
        [0.0005786177_Rkind,-0.0000085613_Rkind,-0.0000183408_Rkind,             &
        -0.0000085613_Rkind, 0.0005786177_Rkind,-0.0000183408_Rkind,             &
        -0.0000183408_Rkind,-0.0000183408_Rkind, 0.0003581481_Rkind],shape=[3,3])

        QModel%Qref = [1.8107934895553828_Rkind,1.8107934895553828_Rkind,1.8230843491285216_Rkind]

        IF (QModel%PubliUnit) THEN
          write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad], Energy: [Hartree]'
        ELSE
          write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad]:, Energy: [Hartree]'
        END IF

      CASE Default

      write(out_unitp,*) ' ERROR in Init_QML_H2O '
      write(out_unitp,*) ' This option is not possible. option: ',QModel%option
      write(out_unitp,*) ' Its value MUST be 1'
      STOP 'ERROR in Init_QML_H2O: wrong option'

    END SELECT


    IF (debug) write(out_unitp,*) 'init Q0 of H2O'
    allocate(QModel%Q0(QModel%ndim))
    CALL get_Q0_QML_H2O(QModel%Q0,QModel,option=0)
    IF (debug) write(out_unitp,*) 'QModel%Q0',QModel%Q0

    IF (debug) write(out_unitp,*) 'init d0GGdef of H2O'
    flush(out_unitp)

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_H2O
!> @brief Subroutine wich prints the QML_H2O parameters.
!!
!! @param QModel            CLASS(QML_H2O_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_QML_H2O(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_H2O_t),     intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'H2O current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '    R2  \ a                            '
    write(nio,*) '         O------------H               '
    write(nio,*) '              R1                       '
    write(nio,*) '                                       '
    write(nio,*) '  Coordinates (option 1):          '
    write(nio,*) '  Q(1) = R1         (Bohr)             '
    write(nio,*) '  Q(2) = R2         (Bohr)             '
    write(nio,*) '  Q(3) = a          (Radian)           '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) ' Quadratic model potential for H2O;    '
    write(nio,*) ' TIPS force constants taken from:      '  
    write(nio,*) '  Dang and Pettitt, J. Chem. Phys. 91 (1987)'
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)
    write(nio,*) '---------------------------------------'

    SELECT CASE (QModel%option)

    CASE (1)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: B3LYP/cc-pVTZ                 '
      write(nio,*) '                                       '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = R1   = 1.8107934895553828 (Bohr)'
      write(nio,*) '  Q(2) = R2   = 1.8107934895553828 (Bohr)'
      write(nio,*) '  Q(3) = a    = 1.8230843491285216 (Rad) '
      write(nio,*) '                                       '
      write(nio,*) '  V = 0.0          Hartree             '
      write(nio,*) ' grad(:) =[0.0,0.0,0.0]                '
      write(nio,*) ' hess    =[0.4726717042854987, 0.00 ,0.00'
      write(nio,*) '           0.00, 0.4726717042854987, 0.00'
      write(nio,*) '           0.00, 0.00, 0.1085242579142319]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '
    CASE Default
        write(out_unitp,*) ' ERROR in write_QModel '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1'
        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end H2O current parameters'

  END SUBROUTINE Write_QML_H2O

  SUBROUTINE get_Q0_QML_H2O(Q0,QModel,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (QML_H2O_t),            intent(in)    :: QModel
    integer,                     intent(in)    :: option

    IF (size(Q0) /= 3) THEN
      write(out_unitp,*) ' ERROR in get_Q0_QML_H2O '
      write(out_unitp,*) ' The size of Q0 is not ndim=3: '
      write(out_unitp,*) ' size(Q0)',size(Q0)
      STOP 'ERROR in get_Q0_QML_H2O: wrong Q0 size'
    END IF


    SELECT CASE (QModel%option)
    CASE (1,3,5) ! a,r+,r-
      Q0(:) = QModel%Qref(:)

    CASE (2,4,6) ! R1,R2,a
      Q0(:) = [QModel%Qref(2),QModel%Qref(2),QModel%Qref(1)]

    CASE Default
      write(out_unitp,*) ' ERROR in get_Q0_QML_H2O '
      write(out_unitp,*) ' This option is not possible. option: ',QModel%option
      write(out_unitp,*) ' Its value MUST be 1,2,3,4,5,6'
      STOP 'ERROR in get_Q0_QML_H2O: wrong option'
    END SELECT

  END SUBROUTINE get_Q0_QML_H2O
!> @brief Subroutine wich calculates the H2O potential (unpublished model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_H2O_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H2O(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2O_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable :: dnQsym(:)


    SELECT CASE (QModel%option)
    CASE (1) ! R1,R2,a
      CALL EvalPot1_QML_H2O(QModel,Mat_OF_PotDia,dnQ,nderiv)
    CASE Default
      write(out_unitp,*) ' ERROR in EvalPot_QML_H2O '
      write(out_unitp,*) ' This option is not possible. option: ',QModel%option
      write(out_unitp,*) ' Its value MUST be 1'
      STOP 'ERROR in EvalPot_QML_H2O: wrong option'
    END SELECT


  END SUBROUTINE EvalPot_QML_H2O

  SUBROUTINE EvalPot1_QML_H2O(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2O_t),     intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable :: DQ(:)
    integer                   :: i
    real (kind=Rkind), parameter :: a = 0.4726717042854987_Rkind
    real (kind=Rkind), parameter :: b = 0.1085242579142319_Rkind

    allocate(DQ(QModel%ndim))
    DQ(:) = dnQ(:) - QModel%Qref(:)


    Mat_OF_PotDia(1,1) = HALF*a * (DQ(1)**2 + DQ(2)**2) + HALF*b * DQ(3)**2
    !CALL Write_dnS(Mat_OF_PotDia(1,1),nio=out_unitp)

   CALL dealloc_dnS(DQ)

   !write(out_unitp,*) ' end EvalPot1_QML_H2O' ; flush(6)

 END SUBROUTINE EvalPot1_QML_H2O
END MODULE QML_H2O_m
