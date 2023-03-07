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
module QML_Template_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
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
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_Template
    PROCEDURE :: Write_QModel    => Write_QML_Template
  END TYPE QML_Template_t

  PUBLIC :: QML_Template_t,Init_QML_Template


contains

  FUNCTION Init_QML_Template(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m, ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_Template_t)              :: QModel

    TYPE(QML_Empty_t),    intent(in)   :: QModel_in ! variable to transfer info to the init
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
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) '  read_param:     ',read_param
      write(out_unit,*) '  nio_param_file: ',nio_param_file
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in
    QModel%pot_name = 'template'
    QModel%nsurf    = 1
    QModel%ndim     = 3

    IF (debug) write(out_unit,*) 'init QModel%morseXpY'
    ! V(1,1) term, all parameters in atomic unit (Hartree, bohr)
    !QModel%morseXpY = Init_QML_Morse(D=0.1_Rkind, a=1._Rkind,req=2._Rkind) ! does not work !!
    XpYeq = 2._Rkind
    CALL Init0_QML_Morse(QModel%morseXpY,D=0.1_Rkind,a=1._Rkind,req=XpYeq,model_name='morseXpY')

    IF (debug) write(out_unit,*) 'init QModel%morseZ'
    !QModel%morseZ   = Init_QML_Morse(D=0.08_Rkind,a=1._Rkind,req=2._Rkind) ! does not work !!
    Zeq  = 2._Rkind
    CALL Init0_QML_Morse(QModel%morseZ,D=0.08_Rkind,a=1._Rkind,req=Zeq,model_name='morseZ')


    IF (debug) write(out_unit,*) 'init Q0 of Template'
    QModel%Q0 = [XpYeq/TWO,XpYeq/TWO,Zeq]


    IF (debug) write(out_unit,*) 'init d0GGdef of Template'
    QModel%d0GGdef = Identity_Mat(QModel%ndim) / 2000._Rkind

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Template


  SUBROUTINE EvalPot_QML_Template(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS (QML_Template_t),   intent(in)     :: QModel
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
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) ' nderiv:',nderiv
      write(out_unit,*) ' Q(:):',(get_d0(dnQ(i)),i=1,size(dnQ))
      flush(out_unit)
    END IF

    Mat_OF_PotDia(1,1) = QModel%kXmY*HALF


    Mat_OF_PotDia(1,1) = QModel%kXmY*HALF * (dnQ(1)-dnQ(2))**2

    mXpY = QML_dnMorse(dnQ(1)+dnQ(2),QModel%morseXpY)

    mZ   = QML_dnMorse(dnQ(3),QModel%morseZ)

    vXmY = (QModel%kXmY*HALF) * (dnQ(1)-dnQ(2))**2

    Mat_OF_PotDia(1,1) = mXpY+mZ+vXmY

    CALL dealloc_dnS(mXpY)
    CALL dealloc_dnS(mZ)
    CALL dealloc_dnS(vXmY)

    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE EvalPot_QML_Template

  SUBROUTINE Write_QML_Template(QModel,nio)
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
    write(nio,*) '========================================'
    write(nio,*) 'QML_Template current parameters'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
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

end module QML_Template_m
