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
module QML_Template_m
  USE QML_Empty_m
  USE QML_Morse_m
  IMPLICIT NONE

  PRIVATE

  TYPE, EXTENDS (QML_Empty_t) :: QML_Template_t
     PRIVATE

     TYPE (QML_Morse_t)  :: morseXpY
     TYPE (QML_Morse_t)  :: morseZ

     real(kind=Rkind)     :: kXmY=0.1_Rkind

  CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_Template
    PROCEDURE :: Write_QModel    => Write_QML_Template
    PROCEDURE :: Write0_QModel   => Write0_QML_Template
  END TYPE QML_Template_t

  PUBLIC :: QML_Template_t,Init_QML_Template


contains

  FUNCTION Init_QML_Template(QModel_in,read_param,nio_param_file) RESULT(QModel)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    TYPE (QML_Template_t)             :: QModel

    TYPE(QML_Empty_t),   intent(in)   :: QModel_in ! variable to transfer info to the init
    logical,              intent(in)   :: read_param
    integer,              intent(in)   :: nio_param_file

    real(kind=Rkind)      :: XpYeq,Zeq
    integer               :: i

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Template'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) '  read_param:     ',read_param
      write(out_unitp,*) '  nio_param_file: ',nio_param_file
      flush(out_unitp)
    END IF

    !QModel%QML_Empty_t = Init_QML_Empty(QModel_in) ! it does not work with nagfor
    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)
    QModel%pot_name = 'template'
    QModel%nsurf    = 1
    QModel%ndim     = 3

    IF (debug) write(out_unitp,*) 'init QModel%morseXpY'
    ! V(1,1) term, all parameters in atomic unit (Hartree, bohr)
    !QModel%morseXpY = Init_QML_Morse(D=0.1_Rkind, a=1._Rkind,req=2._Rkind) ! does not work !!
    XpYeq = 2._Rkind
    CALL Init0_QML_Morse(QModel%morseXpY,D=0.1_Rkind,a=1._Rkind,req=XpYeq,model_name='morseXpY')

    IF (debug) write(out_unitp,*) 'init QModel%morseZ'
    !QModel%morseZ   = Init_QML_Morse(D=0.08_Rkind,a=1._Rkind,req=2._Rkind) ! does not work !!
    Zeq  = 2._Rkind
    CALL Init0_QML_Morse(QModel%morseZ,D=0.08_Rkind,a=1._Rkind,req=Zeq,model_name='morseZ')


    IF (debug) write(out_unitp,*) 'init Q0 of Template'
    QModel%Q0 = [XpYeq/TWO,XpYeq/TWO,Zeq]


    IF (debug) write(out_unitp,*) 'init d0GGdef of Template'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)
    QModel%d0GGdef(:,:) = QModel%d0GGdef / 2000._Rkind

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_Template


  SUBROUTINE EvalPot_QML_Template(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS (QML_Template_t),  intent(in)     :: QModel
    TYPE (dnS_t),             intent(in)     :: dnQ(:)
    TYPE (dnS_t),             intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,                  intent(in)     :: nderiv

    TYPE (dnS_t)  :: mXpY,vXmY,mZ
    integer       :: i

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_Template'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) ' nderiv:',nderiv
      write(out_unitp,*) ' Q(:):',(QML_get_d0_FROM_dnS(dnQ(i)),i=1,size(dnQ))
      flush(out_unitp)
    END IF

    Mat_OF_PotDia(1,1) = QModel%kXmY*HALF


    Mat_OF_PotDia(1,1) = QModel%kXmY*HALF * (dnQ(1)-dnQ(2))**2

    mXpY = QML_dnMorse(dnQ(1)+dnQ(2),QModel%morseXpY)

    mZ   = QML_dnMorse(dnQ(3),QModel%morseZ)

    vXmY = (QModel%kXmY*HALF) * (dnQ(1)-dnQ(2))**2

    Mat_OF_PotDia(1,1) = mXpY+mZ+vXmY

    CALL QML_dealloc_dnS(mXpY)
    CALL QML_dealloc_dnS(mZ)
    CALL QML_dealloc_dnS(vXmY)

    IF (debug) THEN
      write(out_unitp,*) 'Mat_OF_PotDia'
      CALL QML_Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unitp,*)
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE EvalPot_QML_Template

  SUBROUTINE Write_QML_Template(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Template_t), intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) '========================================'
    write(nio,*) 'QML_Template current parameters'
    CALL Write_QML_Empty(QModel%QML_Empty_t,nio)
    write(nio,*)
    CALL Write_QML_Morse(QModel%morseXpY,nio)
    write(nio,*)
    CALL Write_QML_Morse(QModel%morseZ,nio)
    write(nio,*)
    write(nio,*) ' kXmY:',QModel%kXmY
    write(nio,*)
    write(nio,*) 'end QML_Template current parameters'
    write(nio,*) '========================================'
    flush(nio)

  END SUBROUTINE Write_QML_Template
  SUBROUTINE Write0_QML_Template(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_Template_t), intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'QML_Template default parameters'
    write(nio,*)
    write(nio,*) 'V(X,Y,Z) = morez(Z) + morsexpy(X+Y) + 1/2.kXmY.(X-Y)^2'
    write(nio,*)
    write(nio,*) '  Morsez (Dz*(1-exp(-az*(z-zeq)))**2) parameters:'
    write(nio,*) '    Dz   = 0.08 Hartree'
    write(nio,*) '    az   = 1.   Bohr^-2'
    write(nio,*) '    zeq  = 2.   Bohr'
    write(nio,*)
    write(nio,*) '  morsexpy (Dxy*(1-exp(-axy*(xy-xyeq)))**2) parameters:'
    write(nio,*) '    Dxy   = 0.10 Hartree'
    write(nio,*) '    axy   = 1.   Bohr^-2'
    write(nio,*) '    xyeq  = 2.   Bohr'
    write(nio,*)
    write(nio,*) '  kXmY = 0.1 Hartree.bohr^2'
    write(nio,*)
    write(nio,*) '  Potential Value at: Q:'
    write(nio,*) '      2.0000000000       2.0000000000       2.0000000000'
    write(nio,*) '    V = 0.0747645072'
    write(nio,*) '    gradient = [0.0234039289       0.0234039289       0.0000000000]'
    write(nio,*) '    hessian'
    write(nio,*) '      1        0.080259199      -0.119740801       0.000000000'
    write(nio,*) '      2       -0.119740801       0.080259199       0.000000000'
    write(nio,*) '      3        0.000000000       0.000000000       0.160000000'
    write(nio,*)
    write(nio,*) 'end QML_Template default parameters'
    flush(nio)

  END SUBROUTINE Write0_QML_Template

end module QML_Template_m
