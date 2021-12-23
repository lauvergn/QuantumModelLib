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
!    Copyright 2016 David Lauvergnat [1]
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

!> @brief Module which makes the initialization, calculation of the HONO potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_HONO_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the HONO parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_HONO_t
   PRIVATE

      real(kind=Rkind) :: Qref(6)=[2.696732586_Rkind,1.822912197_Rkind,           &
                                              1.777642018_Rkind,2.213326419_Rkind,&
                                              1.9315017_Rkind,ZERO]

      real(kind=Rkind) :: Qtrans(6)=[2.696732586_Rkind,1.822912197_Rkind,           &
                                                1.777642018_Rkind,2.213326419_Rkind,&
                                                1.9315017_Rkind,ZERO]

      real(kind=Rkind) :: Qcis(6)=[2.696732586_Rkind,1.822912197_Rkind,           &
                                              1.777642018_Rkind,2.213326419_Rkind,&
                                              1.9315017_Rkind,ZERO]

   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_HONO
    PROCEDURE :: Write_QModel    => Write_QML_HONO
    PROCEDURE :: Write0_QModel   => Write0_QML_HONO
  END TYPE QML_HONO_t

  PUBLIC :: QML_HONO_t,Init_QML_HONO

  CONTAINS
!> @brief Function which makes the initialization of the HONO parameters.
!!
!! @param QModel             TYPE(QML_HONO_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_HONO(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_HONO_t)                            :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_HONO'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 1
    QModel%ndim     = 6
    QModel%pot_name = 'hono'

    IF (QModel%option < 0 .OR. QModel%option > 2) QModel%option = 1

    IF (debug) write(out_unitp,*) 'init Q0 of HONO'
    SELECT CASE (QModel%option)
    CASE (0) ! ref
      QModel%Q0 = QModel%Qref
    CASE (1) ! trans
      QModel%Q0 = QModel%Qtrans
    CASE (2) ! cis
      QModel%Q0 = QModel%Qcis
    CASE Default

      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' This option is not possible. option: ',QModel%option
      write(out_unitp,*) ' Its value MUST be 0,1,2'

      STOP
    END SELECT


    IF (debug) write(out_unitp,*) 'init d0GGdef of HONO'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    IF (QModel%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_HONO
!> @brief Subroutine wich prints the QML_HONO parameters.
!!
!! @param QModel            CLASS(QML_HONO_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_HONO(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HONO_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'HONO current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '                            +Q(6)      '
    write(nio,*) '      H                  O             '
    write(nio,*) '       \                /              '
    write(nio,*) '   Q(2) \ Q(3)    Q(5) /  Q(4)         '
    write(nio,*) '         O------------N                '
    write(nio,*) '             Q(1)                      '
    write(nio,*) '                                       '
    write(nio,*) '  Q(1) = dON (Bohr)                    '
    write(nio,*) '  Q(2) = dOH (Bohr)                    '
    write(nio,*) '  Q(3) = aHON (Radian)                 '
    write(nio,*) '  Q(4) = dNO (Bohr)                    '
    write(nio,*) '  Q(5) = aONO (Radian)                 '
    write(nio,*) '  Q(6) = Torsion (HONO) (Radian)       '
    write(nio,*) '                                       '
    write(nio,*) '  Trans minimum:                       '
    write(nio,*) '  Q(:)=[2.696732586,1.822912197,       '
    write(nio,*) '        1.777642018,2.213326419,       '
    write(nio,*) '        1.9315017,pi]                  '
    write(nio,*) '  V = 0.0000000000 Hartree             '
    write(nio,*) ' grad(:) = [0.0000000000,-0.0000000005,'
    write(nio,*) '            0.0000000000, 0.0000000000,'
    write(nio,*) '            0.0000000000,-0.0000000000]'
    write(nio,*) '                                       '
    write(nio,*) '  Cis minimum:                         '
    write(nio,*) '  Q(:)=[2.63122,1.84164,1.822274,      '
    write(nio,*) '        2.23738,1.975200,0.]           '
    write(nio,*) '  V = 0.0004240788 Hartree             '
    write(nio,*) ' grad(:) =[-0.0000003952,-0.0000018766,'
    write(nio,*) '           -0.0000000574, 0.0000030809,'
    write(nio,*) '            0.0000003894, 0.0000000000]'
    write(nio,*) 'Remark:'
    write(nio,*) '  V=0.000424079 a.u. (from Richter)'
    write(nio,*) 'Refs:'
    write(nio,*) 'F. Richter, M. Hochlaf, P. Rosmus, F. Gatti, and H.-D. Meyer, JCP. 120, 1306 (2004)'
    write(nio,*) '     https://doi.org/10.1063/1.1632471'
    write(nio,*) 'F. Richter, F. Gatti, C. Léonard, F. Le Quéré, and H.-D. Meyer, JCP 127, 164315 (2007)'
    write(nio,*) '     https://doi.org/10.1063/1.2784553'
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)

    SELECT CASE (QModel%option)

    CASE (0,1,2)

    CONTINUE

    CASE Default
        write(out_unitp,*) ' ERROR in Write_QML_HONO '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 0,1,2'

        STOP
    END SELECT

    write(nio,*)
    write(nio,*) 'end HONO current parameters'

  END SUBROUTINE Write_QML_HONO
  SUBROUTINE Write0_QML_HONO(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HONO_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'HONO default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end HONO default parameters'


  END SUBROUTINE Write0_QML_HONO

!> @brief Subroutine wich calculates the HONO potential with derivatives up to the 2d order.
!!
!! @param QModel             TYPE(QML_HONO_t):    derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_HONO(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS(QML_HONO_t),  intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    SELECT CASE (QModel%option)

    CASE (0,1,2)
      CALL EvalPot1_QML_HONO(Mat_OF_PotDia,dnQ)

    CASE Default
      write(out_unitp,*) ' ERROR in EvalPot_QML_HONO '
      write(out_unitp,*) ' This option is not possible. option: ',QModel%option
      write(out_unitp,*) ' Its value MUST be 1'

      STOP
    END SELECT

  END SUBROUTINE EvalPot_QML_HONO

!> @brief Subroutine wich calculates the HONO potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated

  SUBROUTINE EvalPot1_QML_HONO(Mat_OF_PotDia,dnQ)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    TYPE (dnS_t),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)    :: dnQ(:)


    TYPE (dnS_t)        :: Qw(6)
    TYPE (dnS_t)        :: d6(4),q1(4),q2(4),q3(6),t1(4),t2(4)
    TYPE (dnS_t)        :: Vtemp
    real(kind=Rkind), parameter :: Q0(6)=[2.696732586_Rkind,1.822912197_Rkind,&
                                          1.777642018_Rkind,2.213326419_Rkind,&
                                          1.9315017_Rkind,ZERO]

    Qw(:) = dnQ(:)-Q0(:)
    !Qw(1)=RN=r3,Qw(2)=RH=r1,Qw(3)=AH=ad6(1)(u1),
    !Qw(4)=RO=r2,Qw(5)=AO=ad6(1)(u2),Qw(6)=DO=t
    q1(1) = ONE-exp(-0.70_Rkind*Qw(4))
    q1(2) = q1(1) * q1(1)
    q1(3) = q1(2) * q1(1)
    q1(4) = q1(3) * q1(1)

    q2(1) = ONE-exp(-0.70_Rkind*Qw(1))
    q2(2) = q2(1) * q2(1)
    q2(3) = q2(2) * q2(1)
    q2(4) = q2(3) * q2(1)

    q3(1) = ONE-exp(-0.70_Rkind*Qw(2))
    q3(2) = q3(1) * q3(1)
    q3(3) = q3(2) * q3(1)
    q3(4) = q3(3) * q3(1)
    q3(5) = q3(4) * q3(1)
    q3(6) = q3(5) * q3(1)

    t1(1) = Qw(3)
    t1(2) = t1(1) * t1(1)
    t1(3) = t1(2) * t1(1)
    t1(4) = t1(3) * t1(1)

    t2(1) = Qw(5)
    t2(2) = t2(1) * t2(1)
    t2(3) = t2(2) * t2(1)
    t2(4) = t2(3) * t2(1)

    d6(1) = Cos(Qw(6))
    d6(2) = Cos(Qw(6)*TWO)
    d6(3) = Cos(Qw(6)*THREE)
    d6(4) = Cos(Qw(6)*FOUR)


Vtemp = -0.189952695_Rkind + 0.00109424731_Rkind*d6(1) - &
           0.00967354192_Rkind*d6(2) - &
           0.000358767622_Rkind*d6(3) + &
           0.000122463662_Rkind*d6(4) - &
           0.00980743104_Rkind*q1(1) - &
           0.0114986795_Rkind*d6(1)*q1(1) - &
           0.00182412832_Rkind*d6(2)*q1(1) + &
           0.000113292697_Rkind*d6(3)*q1(1) + &
           0.000246172625_Rkind*d6(4)*q1(1) + &
           0.456806247_Rkind*q1(2) + &
           0.000137509436_Rkind*d6(1)*q1(2) + &
           0.0161360489_Rkind*d6(2)*q1(2) + &
           0.00310754905_Rkind*d6(3)*q1(2) + &
           0.00254032729_Rkind*d6(4)*q1(2) - &
           1.24629979_Rkind*q1(3) - &
           0.0420624218_Rkind*d6(1)*q1(3) - &
           0.0596051696_Rkind*d6(2)*q1(3) + &
           0.0172363933_Rkind*d6(3)*q1(3) + &
           0.991842775_Rkind*q1(4) - &
           0.0209354874_Rkind*d6(1)*q1(4) - &
           0.0141578541_Rkind*q2(1) + &
           0.00219268323_Rkind*d6(1)*q2(1) + &
           0.0181054877_Rkind*d6(2)*q2(1) + &
           0.00156869562_Rkind*d6(3)*q2(1) - &
           0.000186254704_Rkind*d6(4)*q2(1) + &
           0.258630192_Rkind*q1(1)*q2(1) + &
           0.00557018607_Rkind*d6(1)*q1(1)*q2(1) + &
           0.0222968815_Rkind*d6(2)*q1(1)*q2(1) + &
           0.00155489001_Rkind*d6(3)*q1(1)*q2(1) - &
           0.00174933896_Rkind*d6(4)*q1(1)*q2(1) - &
           0.510328693_Rkind*q1(2)*q2(1) - &
           0.0328232234_Rkind*d6(1)*q1(2)*q2(1) + &
           0.281091678_Rkind*d6(2)*q1(2)*q2(1) + &
           0.0632351065_Rkind*d6(3)*q1(2)*q2(1) + &
           0.354308631_Rkind*q1(3)*q2(1) - &
           0.0607342525_Rkind*d6(1)*q1(3)*q2(1) + &
           0.154885992_Rkind*q2(2)
Vtemp = Vtemp - &
           0.00311089401_Rkind*d6(1)*q2(2) - &
           0.0122791104_Rkind*d6(2)*q2(2) - &
           0.00137671372_Rkind*d6(3)*q2(2) + &
           0.00147097007_Rkind*d6(4)*q2(2) - &
           0.0192894747_Rkind*q1(1)*q2(2) - &
           0.00930879651_Rkind*d6(1)*q1(1)*q2(2) + &
           0.0379335753_Rkind*d6(2)*q1(1)*q2(2) - &
           0.000739825929_Rkind*d6(3)*q1(1)*q2(2) + &
           0.104269782_Rkind*q1(2)*q2(2) - &
           0.0566142298_Rkind*d6(1)*q1(2)*q2(2) - &
           0.219359332_Rkind*q2(3) + &
           0.00342798063_Rkind*d6(1)*q2(3) + &
           0.0185860784_Rkind*d6(2)*q2(3) + &
           0.00224338813_Rkind*d6(3)*q2(3) - &
           0.0205495263_Rkind*q1(1)*q2(3) - &
           0.0268389994_Rkind*d6(1)*q1(1)*q2(3) + &
           0.0991511624_Rkind*q2(4) - &
           0.0026913524_Rkind*d6(1)*q2(4) - &
           0.00242284791_Rkind*q3(1) - &
           0.00459921985_Rkind*d6(1)*q3(1) - &
           0.00279190874_Rkind*d6(2)*q3(1) - &
           0.000949245129_Rkind*d6(3)*q3(1) - &
           0.00033370905_Rkind*d6(4)*q3(1) + &
           0.00894355549_Rkind*q1(1)*q3(1) + &
           0.0111935833_Rkind*d6(1)*q1(1)*q3(1) - &
           0.0238707666_Rkind*d6(2)*q1(1)*q3(1) - &
           0.0149416543_Rkind*d6(3)*q1(1)*q3(1) - &
           0.0023285987_Rkind*d6(4)*q1(1)*q3(1) + &
           0.0767590076_Rkind*q1(2)*q3(1) + &
           0.0926272027_Rkind*d6(1)*q1(2)*q3(1) - &
           0.0804465429_Rkind*d6(2)*q1(2)*q3(1) - &
           0.0487228345_Rkind*d6(3)*q1(2)*q3(1) + &
           0.513425067_Rkind*q1(3)*q3(1) + &
           0.187565606_Rkind*d6(1)*q1(3)*q3(1) + &
           0.0255587875_Rkind*q2(1)*q3(1) + &
           0.00915843655_Rkind*d6(1)*q2(1)*q3(1) - &
           0.00774282075_Rkind*d6(2)*q2(1)*q3(1) + &
           0.00252316335_Rkind*d6(3)*q2(1)*q3(1) + &
           0.00131147609_Rkind*d6(4)*q2(1)*q3(1) + &
           0.15324151_Rkind*q1(1)*q2(1)*q3(1) + &
           0.146670106_Rkind*d6(1)*q1(1)*q2(1)*q3(1) - &
           0.19887133_Rkind*d6(2)*q1(1)*q2(1)*q3(1) - &
           0.11810427_Rkind*d6(3)*q1(1)*q2(1)*q3(1) + &
           0.160855899_Rkind*q1(2)*q2(1)*q3(1) + &
           0.706327354_Rkind*d6(1)*q1(2)*q2(1)*q3(1) + &
           0.00159657756_Rkind*q2(2)*q3(1) - &
           0.018933282_Rkind*d6(1)*q2(2)*q3(1) - &
           0.0334347598_Rkind*d6(2)*q2(2)*q3(1) - &
           0.00260068331_Rkind*d6(3)*q2(2)*q3(1) + &
           0.0659378351_Rkind*q1(1)*q2(2)*q3(1) + &
           0.247527495_Rkind*d6(1)*q1(1)*q2(2)*q3(1) + &
           0.0184202212_Rkind*q2(3)*q3(1) + &
           0.0382656334_Rkind*d6(1)*q2(3)*q3(1) + &
           0.419798417_Rkind*q3(2)
Vtemp = Vtemp - &
           0.00448341589_Rkind*d6(1)*q3(2) + &
           0.00130845091_Rkind*d6(2)*q3(2) + &
           0.000116324158_Rkind*d6(3)*q3(2) + &
           0.000670620368_Rkind*d6(4)*q3(2) - &
           0.0505717655_Rkind*q1(1)*q3(2) - &
           0.0487270231_Rkind*d6(1)*q1(1)*q3(2) + &
           0.0211145205_Rkind*d6(2)*q1(1)*q3(2) + &
           0.0203382911_Rkind*d6(3)*q1(1)*q3(2) + &
           0.838195571_Rkind*q1(2)*q3(2) - &
           0.162187822_Rkind*d6(1)*q1(2)*q3(2) - &
           0.0793448881_Rkind*q2(1)*q3(2) + &
           0.0126418131_Rkind*d6(1)*q2(1)*q3(2) + &
           0.113110811_Rkind*d6(2)*q2(1)*q3(2) + &
           0.0193718767_Rkind*d6(3)*q2(1)*q3(2) + &
           0.0556884966_Rkind*q1(1)*q2(1)*q3(2) - &
           0.283143742_Rkind*d6(1)*q1(1)*q2(1)*q3(2) - &
           0.0379021543_Rkind*q2(2)*q3(2) - &
           0.0515344482_Rkind*d6(1)*q2(2)*q3(2) - &
           0.469248343_Rkind*q3(3) + &
           0.00683450913_Rkind*d6(1)*q3(3) - &
           0.0302440097_Rkind*d6(2)*q3(3) - &
           0.00785685084_Rkind*d6(3)*q3(3) + &
           0.0317331596_Rkind*q1(1)*q3(3) + &
           0.0224311461_Rkind*d6(1)*q1(1)*q3(3) - &
           0.014109039_Rkind*q2(1)*q3(3) + &
           0.0784615695_Rkind*d6(1)*q2(1)*q3(3) + &
           0.264118373_Rkind*q3(4) - &
           0.0284378299_Rkind*d6(1)*q3(4) - &
           0.411868833_Rkind*q3(5) - &
           0.0465284723_Rkind*d6(1)*q3(5) + &
           0.684631535_Rkind*q3(6) + &
           0.136390997_Rkind*d6(1)*q3(6) - &
           0.00455237659_Rkind*t1(1) - &
           0.00342497933_Rkind*d6(1)*t1(1) + &
           0.00176367547_Rkind*d6(2)*t1(1) + &
           0.00104093751_Rkind*d6(3)*t1(1) + &
           0.00040465928_Rkind*d6(4)*t1(1) + &
           0.0228697532_Rkind*q1(1)*t1(1) - &
           0.00583129063_Rkind*d6(1)*q1(1)*t1(1) - &
           0.00525965687_Rkind*d6(2)*q1(1)*t1(1) - &
           0.00726039823_Rkind*d6(3)*q1(1)*t1(1) - &
           0.000676636312_Rkind*d6(4)*q1(1)*t1(1) + &
           0.00327215125_Rkind*q1(2)*t1(1) + &
           0.0438207551_Rkind*d6(1)*q1(2)*t1(1) - &
           0.0149841179_Rkind*d6(2)*q1(2)*t1(1) - &
           0.00203628351_Rkind*d6(3)*q1(2)*t1(1) + &
           0.191709174_Rkind*q1(3)*t1(1) + &
           0.093491373_Rkind*d6(1)*q1(3)*t1(1) + &
           0.076350072_Rkind*q2(1)*t1(1) + &
           0.00820471508_Rkind*d6(1)*q2(1)*t1(1) - &
           0.017158215_Rkind*d6(2)*q2(1)*t1(1) - &
           0.00807705701_Rkind*d6(3)*q2(1)*t1(1) - &
           0.000438687742_Rkind*d6(4)*q2(1)*t1(1) + &
           0.0383427126_Rkind*q1(1)*q2(1)*t1(1) - &
           0.00998234682_Rkind*d6(1)*q1(1)*q2(1)*t1(1) - &
           0.0576687154_Rkind*d6(2)*q1(1)*q2(1)*t1(1) - &
           0.0249928303_Rkind*d6(3)*q1(1)*q2(1)*t1(1) + &
           0.0448059036_Rkind*q1(2)*q2(1)*t1(1) + &
           0.149487953_Rkind*d6(1)*q1(2)*q2(1)*t1(1) - &
           0.0337328346_Rkind*q2(2)*t1(1) + &
           0.0155544417_Rkind*d6(1)*q2(2)*t1(1) + &
           0.00134378626_Rkind*d6(2)*q2(2)*t1(1) + &
           0.00316193689_Rkind*d6(3)*q2(2)*t1(1) + &
           0.00767483439_Rkind*q1(1)*q2(2)*t1(1) + &
           0.0968434164_Rkind*d6(1)*q1(1)*q2(2)*t1(1) - &
           0.0305733582_Rkind*q2(3)*t1(1) - &
           0.00197350506_Rkind*d6(1)*q2(3)*t1(1) + &
           0.0087370579_Rkind*q3(1)*t1(1)
Vtemp = Vtemp - &
           0.00815737367_Rkind*d6(1)*q3(1)*t1(1) + &
           0.00712637704_Rkind*d6(2)*q3(1)*t1(1) + &
           0.00869544019_Rkind*d6(3)*q3(1)*t1(1) + &
           0.00152060166_Rkind*d6(4)*q3(1)*t1(1) - &
           0.00265855701_Rkind*q1(1)*q3(1)*t1(1) + &
           0.0455076119_Rkind*d6(1)*q1(1)*q3(1)*t1(1) + &
           0.053347839_Rkind*d6(2)*q1(1)*q3(1)*t1(1) + &
           0.0183506506_Rkind*d6(3)*q1(1)*q3(1)*t1(1) + &
           0.00880111659_Rkind*q1(2)*q3(1)*t1(1) - &
           0.110881402_Rkind*d6(1)*q1(2)*q3(1)*t1(1) - &
           0.0887534042_Rkind*q2(1)*q3(1)*t1(1) - &
           0.0471604585_Rkind*d6(1)*q2(1)*q3(1)*t1(1) - &
           0.0108536348_Rkind*d6(2)*q2(1)*q3(1)*t1(1) - &
           0.00249012995_Rkind*d6(3)*q2(1)*q3(1)*t1(1) + &
           0.384527757_Rkind*q1(1)*q2(1)*q3(1)*t1(1) + &
           0.375548669_Rkind*d6(1)*q1(1)*q2(1)*q3(1)*t1(1) + &
           0.095004823_Rkind*q2(2)*q3(1)*t1(1) + &
           0.0258547162_Rkind*d6(1)*q2(2)*q3(1)*t1(1) + &
           0.0358528992_Rkind*q3(2)*t1(1) + &
           0.0243534749_Rkind*d6(1)*q3(2)*t1(1) + &
           0.000643935255_Rkind*d6(2)*q3(2)*t1(1) + &
           0.00309505006_Rkind*d6(3)*q3(2)*t1(1) + &
           0.122479677_Rkind*q1(1)*q3(2)*t1(1) + &
           0.0992704283_Rkind*d6(1)*q1(1)*q3(2)*t1(1) - &
           0.0688250676_Rkind*q2(1)*q3(2)*t1(1) - &
           0.0236041344_Rkind*d6(1)*q2(1)*q3(2)*t1(1) - &
           0.00198364584_Rkind*q3(3)*t1(1) + &
           0.0116008785_Rkind*d6(1)*q3(3)*t1(1) + &
           0.0430971662_Rkind*t1(2) + &
           0.0042138202_Rkind*d6(1)*t1(2) + &
           0.00872268723_Rkind*d6(2)*t1(2) - &
           0.000354403712_Rkind*d6(3)*t1(2) - &
           0.000828993629_Rkind*d6(4)*t1(2) - &
           0.0108297159_Rkind*q1(1)*t1(2) - &
           0.0254890658_Rkind*d6(1)*q1(1)*t1(2) + &
           0.00899383548_Rkind*d6(2)*q1(1)*t1(2) + &
           0.0185913527_Rkind*d6(3)*q1(1)*t1(2) + &
           0.142859188_Rkind*q1(2)*t1(2) - &
           0.105381071_Rkind*d6(1)*q1(2)*t1(2) - &
           0.144786448_Rkind*q2(1)*t1(2) - &
           0.014820077_Rkind*d6(1)*q2(1)*t1(2) + &
           0.023385895_Rkind*d6(2)*q2(1)*t1(2) + &
           0.0160959682_Rkind*d6(3)*q2(1)*t1(2) + &
           0.174206943_Rkind*q1(1)*q2(1)*t1(2) + &
           0.0721759177_Rkind*d6(1)*q1(1)*q2(1)*t1(2) - &
           0.00784015518_Rkind*q2(2)*t1(2) - &
           0.019216646_Rkind*d6(1)*q2(2)*t1(2) + &
           0.00183572502_Rkind*q3(1)*t1(2) + &
           0.0183112454_Rkind*d6(1)*q3(1)*t1(2) - &
           0.0303663705_Rkind*d6(2)*q3(1)*t1(2) - &
           0.0186101314_Rkind*d6(3)*q3(1)*t1(2) - &
           0.000295514736_Rkind*q1(1)*q3(1)*t1(2) - &
           0.0440716885_Rkind*d6(1)*q1(1)*q3(1)*t1(2) + &
           0.154074326_Rkind*q2(1)*q3(1)*t1(2) + &
           0.0953703164_Rkind*d6(1)*q2(1)*q3(1)*t1(2) + &
           0.0134607961_Rkind*q3(2)*t1(2) - &
           0.0716332099_Rkind*d6(1)*q3(2)*t1(2) - &
           0.0520195143_Rkind*t1(3)
Vtemp = Vtemp - &
           0.00793244153_Rkind*d6(1)*t1(3) - &
           0.000902958062_Rkind*d6(2)*t1(3) - &
           0.000295822999_Rkind*d6(3)*t1(3) + &
           0.0296745986_Rkind*q1(1)*t1(3) + &
           0.0140004163_Rkind*d6(1)*q1(1)*t1(3) + &
           0.0319784133_Rkind*q2(1)*t1(3) + &
           0.000471738266_Rkind*d6(1)*q2(1)*t1(3) + &
           0.040323914_Rkind*q3(1)*t1(3) + &
           0.0107477207_Rkind*d6(1)*q3(1)*t1(3) + &
           0.0191579153_Rkind*t1(4) + &
           0.00714854415_Rkind*d6(1)*t1(4) - &
           0.0101898277_Rkind*t2(1) - &
           0.0119240248_Rkind*d6(1)*t2(1) - &
           0.00141049749_Rkind*d6(2)*t2(1) + &
           0.000683539214_Rkind*d6(3)*t2(1) + &
           0.000359839618_Rkind*d6(4)*t2(1) + &
           0.140439429_Rkind*q1(1)*t2(1) + &
           0.0272860896_Rkind*d6(1)*q1(1)*t2(1) - &
           0.00534097195_Rkind*d6(2)*q1(1)*t2(1) - &
           0.00972014792_Rkind*d6(3)*q1(1)*t2(1) + &
           0.0000887204399_Rkind*d6(4)*q1(1)*t2(1) + &
           0.0787252772_Rkind*q1(2)*t2(1) + &
           0.0383644305_Rkind*d6(1)*q1(2)*t2(1) - &
           0.0898012116_Rkind*d6(2)*q1(2)*t2(1) - &
           0.0187434964_Rkind*d6(3)*q1(2)*t2(1) + &
           0.332543937_Rkind*q1(3)*t2(1) + &
           0.0620785849_Rkind*d6(1)*q1(3)*t2(1) + &
           0.10702459_Rkind*q2(1)*t2(1) + &
           0.00898341305_Rkind*d6(1)*q2(1)*t2(1) - &
           0.018752043_Rkind*d6(2)*q2(1)*t2(1) - &
           0.00841260501_Rkind*d6(3)*q2(1)*t2(1) - &
           0.00218651601_Rkind*d6(4)*q2(1)*t2(1) - &
           0.146551654_Rkind*q1(1)*q2(1)*t2(1) - &
           0.0541718939_Rkind*d6(1)*q1(1)*q2(1)*t2(1) - &
           0.160969317_Rkind*d6(2)*q1(1)*q2(1)*t2(1) - &
           0.0237237205_Rkind*d6(3)*q1(1)*q2(1)*t2(1) + &
           0.29026056_Rkind*q1(2)*q2(1)*t2(1) - &
           0.156599513_Rkind*d6(1)*q1(2)*q2(1)*t2(1) - &
           0.108076274_Rkind*q2(2)*t2(1) - &
           0.020691359_Rkind*d6(1)*q2(2)*t2(1) - &
           0.00610632304_Rkind*d6(2)*q2(2)*t2(1) - &
           0.000629650888_Rkind*d6(3)*q2(2)*t2(1) + &
           0.247912781_Rkind*q1(1)*q2(2)*t2(1) + &
           0.103390714_Rkind*d6(1)*q1(1)*q2(2)*t2(1) - &
           0.0493915552_Rkind*q2(3)*t2(1) - &
           0.00468066434_Rkind*d6(1)*q2(3)*t2(1) - &
           0.0243885454_Rkind*q3(1)*t2(1) - &
           0.0372179913_Rkind*d6(1)*q3(1)*t2(1) + &
           0.0174036663_Rkind*d6(2)*q3(1)*t2(1) + &
           0.00985983973_Rkind*d6(3)*q3(1)*t2(1) + &
           0.00139311559_Rkind*d6(4)*q3(1)*t2(1) - &
           0.0558244531_Rkind*q1(1)*q3(1)*t2(1) + &
           0.0563401662_Rkind*d6(1)*q1(1)*q3(1)*t2(1) + &
           0.10247224_Rkind*d6(2)*q1(1)*q3(1)*t2(1) + &
           0.0374984626_Rkind*d6(3)*q1(1)*q3(1)*t2(1) - &
           0.029018897_Rkind*q1(2)*q3(1)*t2(1) - &
           0.194648012_Rkind*d6(1)*q1(2)*q3(1)*t2(1) - &
           0.0237689328_Rkind*q2(1)*q3(1)*t2(1) + &
           0.0463792516_Rkind*d6(1)*q2(1)*q3(1)*t2(1) + &
           0.015788297_Rkind*d6(2)*q2(1)*q3(1)*t2(1) - &
           0.00555425428_Rkind*d6(3)*q2(1)*q3(1)*t2(1) + &
           0.616253746_Rkind*q1(1)*q2(1)*q3(1)*t2(1) + &
           0.357637604_Rkind*d6(1)*q1(1)*q2(1)*q3(1)*t2(1) - &
           0.022142501_Rkind*q2(2)*q3(1)*t2(1) - &
           0.040878152_Rkind*d6(1)*q2(2)*q3(1)*t2(1) + &
           0.0473943519_Rkind*q3(2)*t2(1) - &
           0.0079409029_Rkind*d6(1)*q3(2)*t2(1) - &
           0.0256743505_Rkind*d6(2)*q3(2)*t2(1) - &
           0.00110769813_Rkind*d6(3)*q3(2)*t2(1) + &
           0.161255254_Rkind*q1(1)*q3(2)*t2(1) + &
           0.107292346_Rkind*d6(1)*q1(1)*q3(2)*t2(1) + &
           0.0102204362_Rkind*q2(1)*q3(2)*t2(1) - &
           0.0188518043_Rkind*d6(1)*q2(1)*q3(2)*t2(1) + &
           0.03073191_Rkind*q3(3)*t2(1) - &
           0.0372863889_Rkind*d6(1)*q3(3)*t2(1) + &
           0.00605146357_Rkind*t1(1)*t2(1)
Vtemp = Vtemp - &
           0.0625016337_Rkind*d6(1)*t1(1)*t2(1) + &
           0.00545661919_Rkind*d6(2)*t1(1)*t2(1) + &
           0.00213458974_Rkind*d6(3)*t1(1)*t2(1) - &
           0.00123142747_Rkind*d6(4)*t1(1)*t2(1) - &
           0.0386098149_Rkind*q1(1)*t1(1)*t2(1) + &
           0.0322323164_Rkind*d6(1)*q1(1)*t1(1)*t2(1) + &
           0.0172990392_Rkind*d6(2)*q1(1)*t1(1)*t2(1) + &
           0.0194093374_Rkind*d6(3)*q1(1)*t1(1)*t2(1) - &
           0.0165358157_Rkind*q1(2)*t1(1)*t2(1) - &
           0.251732966_Rkind*d6(1)*q1(2)*t1(1)*t2(1) - &
           0.0291813159_Rkind*q2(1)*t1(1)*t2(1) - &
           0.0285841202_Rkind*d6(1)*q2(1)*t1(1)*t2(1) + &
           0.035623705_Rkind*d6(2)*q2(1)*t1(1)*t2(1) + &
           0.0106192912_Rkind*d6(3)*q2(1)*t1(1)*t2(1) + &
           0.332547773_Rkind*q1(1)*q2(1)*t1(1)*t2(1) + &
           0.244875852_Rkind*d6(1)*q1(1)*q2(1)*t1(1)*t2(1) - &
           0.0384616389_Rkind*q2(2)*t1(1)*t2(1) + &
           0.0629614664_Rkind*d6(1)*q2(2)*t1(1)*t2(1) - &
           0.00564053761_Rkind*q3(1)*t1(1)*t2(1) + &
           0.0129224721_Rkind*d6(1)*q3(1)*t1(1)*t2(1) - &
           0.0227842179_Rkind*d6(2)*q3(1)*t1(1)*t2(1) - &
           0.0120424733_Rkind*d6(3)*q3(1)*t1(1)*t2(1) - &
           0.211198019_Rkind*q1(1)*q3(1)*t1(1)*t2(1) - &
           0.240333121_Rkind*d6(1)*q1(1)*q3(1)*t1(1)*t2(1) - &
           0.00656922172_Rkind*q2(1)*q3(1)*t1(1)*t2(1) - &
           0.0487814318_Rkind*d6(1)*q2(1)*q3(1)*t1(1)*t2(1) - &
           0.0154154286_Rkind*q3(2)*t1(1)*t2(1) - &
           0.0857145723_Rkind*d6(1)*q3(2)*t1(1)*t2(1) - &
           0.0272882706_Rkind*t1(2)*t2(1) - &
           0.00850651755_Rkind*d6(1)*t1(2)*t2(1) - &
           0.00299544666_Rkind*d6(2)*t1(2)*t2(1) - &
           0.000147364731_Rkind*d6(3)*t1(2)*t2(1) + &
           0.0615294787_Rkind*q1(1)*t1(2)*t2(1) + &
           0.027017626_Rkind*d6(1)*q1(1)*t1(2)*t2(1) - &
           0.0288241546_Rkind*q2(1)*t1(2)*t2(1) - &
           0.0486327764_Rkind*d6(1)*q2(1)*t1(2)*t2(1) + &
           0.0488017381_Rkind*q3(1)*t1(2)*t2(1) + &
           0.0102779306_Rkind*d6(1)*q3(1)*t1(2)*t2(1) + &
           0.00673112764_Rkind*t1(3)*t2(1) + &
           0.0273162997_Rkind*d6(1)*t1(3)*t2(1) + &
           0.233236597_Rkind*t2(2) + &
           0.0185509055_Rkind*d6(1)*t2(2) + &
           0.00637435171_Rkind*d6(2)*t2(2) + &
           0.00246218687_Rkind*d6(3)*t2(2) + &
           0.000386653144_Rkind*d6(4)*t2(2) - &
           0.326720561_Rkind*q1(1)*t2(2) - &
           0.0381322227_Rkind*d6(1)*q1(1)*t2(2) + &
           0.00674174458_Rkind*d6(2)*q1(1)*t2(2) + &
           0.0326043209_Rkind*d6(3)*q1(1)*t2(2) + &
           0.259344792_Rkind*q1(2)*t2(2) - &
           0.205327174_Rkind*d6(1)*q1(2)*t2(2) - &
           0.413490225_Rkind*q2(1)*t2(2) - &
           0.053435954_Rkind*d6(1)*q2(1)*t2(2) + &
           0.0668849012_Rkind*d6(2)*q2(1)*t2(2) + &
           0.0225664395_Rkind*d6(3)*q2(1)*t2(2) + &
           0.530653278_Rkind*q1(1)*q2(1)*t2(2) + &
           0.171724944_Rkind*d6(1)*q1(1)*q2(1)*t2(2) - &
           0.0520286772_Rkind*q2(2)*t2(2) - &
           0.00679308378_Rkind*d6(1)*q2(2)*t2(2) + &
           0.0420264848_Rkind*q3(1)*t2(2) + &
           0.0523998353_Rkind*d6(1)*q3(1)*t2(2) - &
           0.0398707139_Rkind*d6(2)*q3(1)*t2(2) - &
           0.0203991009_Rkind*d6(3)*q3(1)*t2(2) + &
           0.0667289964_Rkind*q1(1)*q3(1)*t2(2) - &
           0.0342822152_Rkind*d6(1)*q1(1)*q3(1)*t2(2) + &
           0.0343237583_Rkind*q2(1)*q3(1)*t2(2) + &
           0.0412379456_Rkind*d6(1)*q2(1)*q3(1)*t2(2) + &
           0.0871897999_Rkind*q3(2)*t2(2) - &
           0.0574355743_Rkind*d6(1)*q3(2)*t2(2) - &
           0.0159867411_Rkind*t1(1)*t2(2) - &
           0.00265186155_Rkind*d6(1)*t1(1)*t2(2) - &
           0.0102316103_Rkind*d6(2)*t1(1)*t2(2) + &
           0.0000685730434_Rkind*d6(3)*t1(1)*t2(2) + &
           0.0784264654_Rkind*q1(1)*t1(1)*t2(2) + &
           0.0278765022_Rkind*d6(1)*q1(1)*t1(1)*t2(2) + &
           0.0182766244_Rkind*q2(1)*t1(1)*t2(2) + &
           0.0707377216_Rkind*d6(1)*q2(1)*t1(1)*t2(2) - &
           0.00745614437_Rkind*q3(1)*t1(1)*t2(2) - &
           0.0179292705_Rkind*d6(1)*q3(1)*t1(1)*t2(2) + &
           0.0820457732_Rkind*t1(2)*t2(2) - &
           0.00379042877_Rkind*d6(1)*t1(2)*t2(2) - &
           0.204903499_Rkind*t2(3) - &
           0.0397918602_Rkind*d6(1)*t2(3) - &
           0.0112635345_Rkind*d6(2)*t2(3) - &
           0.0031871158_Rkind*d6(3)*t2(3) + &
           0.165897228_Rkind*q1(1)*t2(3) + &
           0.00715467491_Rkind*d6(1)*q1(1)*t2(3) + &
           0.184620167_Rkind*q2(1)*t2(3) + &
           0.0150775061_Rkind*d6(1)*q2(1)*t2(3) + &
           0.00919759591_Rkind*q3(1)*t2(3) - &
           0.0263665552_Rkind*d6(1)*q3(1)*t2(3) + &
           0.0221488545_Rkind*t1(1)*t2(3) - &
           0.000663315562_Rkind*d6(1)*t1(1)*t2(3) + &
           0.203571926_Rkind*t2(4) + &
           0.0344417073_Rkind*d6(1)*t2(4)


    Mat_OF_PotDia(1,1) =  0.2002292529406_Rkind + &
        0.00001_Rkind*exp(FOUR*Qw(2)**2 + FOUR*Qw(4)**2 + &
              Qw(1)**2 + FOUR*Qw(3)**2 + FOUR*Qw(5)**2) + &
        exp(-0.25_Rkind*Qw(2)**2 - Qw(4)**2 - Qw(1)**2/9._Rkind - &
            0.25_Rkind*Qw(3)**2 - 0.25_Rkind*Qw(5)**2)* Vtemp
!-----------------------------------------------------------------------!

   CALL QML_dealloc_dnS(Vtemp)

   CALL QML_dealloc_dnS(d6)
   CALL QML_dealloc_dnS(q1)
   CALL QML_dealloc_dnS(q2)
   CALL QML_dealloc_dnS(q3)
   CALL QML_dealloc_dnS(t1)
   CALL QML_dealloc_dnS(t2)
   CALL QML_dealloc_dnS(Qw)

  END SUBROUTINE EvalPot1_QML_HONO

END MODULE QML_HONO_m
