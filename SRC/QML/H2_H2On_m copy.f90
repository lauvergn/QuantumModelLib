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

!> @brief Module which makes the initialization, calculation of the H2_H2On potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 13/11/2025
!!
MODULE QML_H2_H2On_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the H2_H2On parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  real(kind=Rkind), parameter :: ang2au = 0.529177_Rkind
  real(kind=Rkind), parameter :: au2cm  = 219474.63_Rkind
  real(kind=Rkind), parameter :: calkcm = 83.592_Rkind

  TYPE, EXTENDS (QML_Empty_t) ::  QML_H2_H2On_t

    PRIVATE

    real(kind=Rkind) :: mH  = 1526464003414_Rkind ! NIST2012
    real(kind=Rkind) :: dHH = 0.750902458653675_Rkind / ang2au
    real(kind=Rkind), allocatable :: Qref(:)

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_H2_H2On
    PROCEDURE :: Write_QModel    => Write_QML_H2_H2On
  END TYPE QML_H2_H2On_t

  PUBLIC :: QML_H2_H2On_t,Init_QML_H2_H2On

  CONTAINS
!> @brief Subroutine which makes the initialization of the H2-H2On parameters.
!!
!! @param H2OPot             TYPE(QML_H2_H2On_t):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_H2_H2On(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE (QML_H2_H2On_t)                         :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !local variable
    integer                         :: err_read,nio_fit,i,j,k,idum

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H2_H2On'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf    = 1
    QModel%pot_name = 'H2_H2On'
    QModel%ndim     = 6 ! H2 in cartesian coordinates


    IF (QModel%option /= 1) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1)

      QModel%d0GGdef = Identity_Mat(QModel%ndim)/QModel%mH

      QModel%Qref = [ZERO,ZERO,-QModel%dHH/TWO,  ZERO,ZERO,+QModel%dHH/TWO]

      IF (QModel%PubliUnit) THEN
        write(out_unit,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Bohr,Bohr,Bohr,Bohr], Energy: [Hartree]'
      ELSE
        write(out_unit,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Bohr,Bohr,Bohr,Bohr]:, Energy: [Hartree]'
      END IF

    CASE Default

      write(out_unit,*) ' ERROR in Init_QML_H2_H2On '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1'
      STOP 'ERROR in Init_QML_H2_H2On: wrong option'

    END SELECT


    IF (debug) write(out_unit,*) 'init Q0 of H2_H2On'
    allocate(QModel%Q0(QModel%ndim))
    CALL get_Q0_QML_H2_H2On(QModel%Q0,QModel,option=0)
    IF (debug) write(out_unit,*) 'QModel%Q0',QModel%Q0

    IF (debug) write(out_unit,*) 'init d0GGdef of H2_H2On'
    flush(out_unit)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_H2_H2On
!> @brief Subroutine wich prints the QML_H2_H2On parameters.
!!
!! @param QModel            CLASS(QML_H2_H2On_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_QML_H2_H2On(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_H2_H2On_t),     intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'H2_H2On current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) ' Clathrate: H2@Wn                      '
    write(nio,*) ' H2 in Carteian coordinates            '
    write(nio,*) '  encapsulated by n-water molecules (W)'
    write(nio,*) '  The water cage is rigid.             '
    write(nio,*) '                                       '
    write(nio,*) '  Coordinates:                         '
    write(nio,*) '  H1 = [x1,y1,z1]         (Bohr)       '
    write(nio,*) '  H2 = [x2,y2,z2]         (Bohr)       '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) ' The potential is the SPC/E            '
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
      write(nio,*) '  Level: SPC/E                         '
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
        write(out_unit,*) ' ERROR in write_QModel '
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1'
        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end H2_H2On current parameters'

  END SUBROUTINE Write_QML_H2_H2On

  SUBROUTINE get_Q0_QML_H2_H2On(Q0,QModel,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (QML_H2_H2On_t),        intent(in)    :: QModel
    integer,                     intent(in)    :: option

    IF (size(Q0) /= 6) THEN
      write(out_unit,*) ' ERROR in get_Q0_QML_H2_H2On '
      write(out_unit,*) ' The size of Q0 is not ndim=6: '
      write(out_unit,*) ' size(Q0)',size(Q0)
      STOP 'ERROR in get_Q0_QML_H2_H2On: wrong Q0 size'
    END IF

    Q0(:) = QModel%Qref(:)

  END SUBROUTINE get_Q0_QML_H2_H2On
!> @brief Subroutine wich calculates the H2_H2On potential with derivatives.
!!
!! @param PotVal             TYPE (dnMat_t):      Potential with derivatives,.
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_H2_H2On_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H2_H2On(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2_H2On_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable :: dnQsym(:)


    SELECT CASE (QModel%option)
    CASE (1)
      CALL EvalPot1_QML_H2_H2On(QModel,Mat_OF_PotDia,dnQ,nderiv)
    CASE Default
      write(out_unit,*) ' ERROR in EvalPot_QML_H2_H2On '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1'
      STOP 'ERROR in EvalPot_QML_H2_H2On: wrong option'
    END SELECT


  END SUBROUTINE EvalPot_QML_H2_H2On

  SUBROUTINE EvalPot1_QML_H2_H2On(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2_H2On_t), intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv


    real(kind=Rkind), allocatable :: d0Q(:)
    real(kind=Rkind) :: v


    d0Q = get_d0(dnQ)

    CALL sc_sp_QML(d0Q,v)

    CALL set_d0S(Mat_OF_PotDia(1,1),v)
  
   !write(out_unit,*) ' end EvalPot1_QML_H2_H2On' ; flush(6)

  END SUBROUTINE EvalPot1_QML_H2_H2On

  !***********************************************************************
  SUBROUTINE sc_sp_QML(h2coor,vpot)
    !***********************************************************************
    !     routine to evaluate the interaction potential for one H2 molecule
    !     inside the small cage cavity (20 water molecules)
    !     of the sII structure of clathrate hydrate
    !
    !     input :  array h2coor(3,2)  [bohr] (H atoms cartesian coordinates
    !                                          of the H2 molecule)
    !     output:  vpot H2-large_cage [cm-1]
    !-------------------
    !     bigcag(3,ntot): coordinates
    !                     for the water molecules in the large cage [bohr]
    !***********************************************************************
    real(kind=Rkind), intent(out) :: vpot
    real(kind=Rkind), intent(in)  :: h2coor(3,2)

    !nwater = number of water molecules for the cage
    !ntot   = total number of atoms of the cage
    integer, parameter :: nwater=20
    integer, parameter :: ntot=nwater*3

    real(kind=Rkind) :: vfit,coord(3,5),smallc(3,ntot)
    integer :: i

    !cccccc Small cage
    smallc(1, 1)=   -0.31649d0
    smallc(2, 1)=    0.38003d0
    smallc(3, 1)=   -7.07291d0
    smallc(1, 2)=    3.70808d0
    smallc(2, 2)=   -2.65146d0
    smallc(3, 2)=   -5.67704d0
    smallc(1, 3)=    6.94833d0
    smallc(2, 3)=    0.28393d0
    smallc(3, 3)=   -2.77947d0
    smallc(1, 4)=    4.85546d0
    smallc(2, 4)=    5.15550d0
    smallc(3, 4)=   -2.41820d0
    smallc(1, 5)=    0.34577d0
    smallc(2, 5)=    5.17496d0
    smallc(3, 5)=   -5.09665d0
    smallc(1, 6)=    6.94765d0
    smallc(2, 6)=   -1.90366d0
    smallc(3, 6)=    2.02974d0
    smallc(1, 7)=    4.78465d0
    smallc(2, 7)=    1.66045d0
    smallc(3, 7)=    5.21254d0
    smallc(1, 8)=    3.53305d0
    smallc(2, 8)=    6.04447d0
    smallc(3, 8)=    2.61917d0
    smallc(1, 9)=    0.34023d0
    smallc(2, 9)=   -0.38558d0
    smallc(3, 9)=    7.05546d0
    smallc(1,10)=   -3.68434d0
    smallc(2,10)=    2.64591d0
    smallc(3,10)=    5.65960d0
    smallc(1,11)=   -6.92460d0
    smallc(2,11)=   -0.28948d0
    smallc(3,11)=    2.76202d0
    smallc(1,12)=   -4.83172d0
    smallc(2,12)=   -5.16106d0
    smallc(3,12)=    2.40075d0
    smallc(1,13)=   -0.32204d0
    smallc(2,13)=   -5.18051d0
    smallc(3,13)=    5.07920d0
    smallc(1,14)=   -4.76091d0
    smallc(2,14)=   -1.66601d0
    smallc(3,14)=   -5.22999d0
    smallc(1,15)=   -6.92391d0
    smallc(2,15)=    1.89811d0
    smallc(3,15)=   -2.04719d0
    smallc(1,16)=   -3.50932d0
    smallc(2,16)=   -6.05003d0
    smallc(3,16)=   -2.63661d0
    smallc(1,17)=    1.76223d0
    smallc(2,17)=   -6.66343d0
    smallc(3,17)=   -2.91488d0
    smallc(1,18)=   -1.73849d0
    smallc(2,18)=    6.65787d0
    smallc(3,18)=    2.89744d0
    smallc(1,19)=    3.76898d0
    smallc(2,19)=   -6.16184d0
    smallc(3,19)=    1.94674d0
    smallc(1,20)=   -3.74524d0
    smallc(2,20)=    6.15629d0
    smallc(3,20)=   -1.96419d0
    smallc(1,21)=   -0.08391d0
    smallc(2,21)=    2.06393d0
    smallc(3,21)=   -6.37888d0
    smallc(1,22)=   -0.38327d0
    smallc(2,22)=    0.63072d0
    smallc(3,22)=   -8.89059d0
    smallc(1,23)=    2.26809d0
    smallc(2,23)=   -1.61376d0
    smallc(3,23)=   -6.14702d0
    smallc(1,24)=    2.98751d0
    smallc(2,24)=   -4.03968d0
    smallc(3,24)=   -4.71534d0
    smallc(1,25)=    5.81406d0
    smallc(2,25)=   -0.74362d0
    smallc(3,25)=   -3.79378d0
    smallc(1,26)=    8.55826d0
    smallc(2,26)=    0.20693d0
    smallc(3,26)=   -3.65894d0
    smallc(1,27)=    5.57003d0
    smallc(2,27)=    3.46539d0
    smallc(3,27)=   -2.48266d0
    smallc(1,28)=    4.41648d0
    smallc(2,28)=    5.40842d0
    smallc(3,28)=   -0.65338d0
    smallc(1,29)=    1.86723d0
    smallc(2,29)=    5.18989d0
    smallc(3,29)=   -4.06892d0
    smallc(1,30)=   -1.01662d0
    smallc(2,30)=    5.52546d0
    smallc(3,30)=   -3.91669d0
    smallc(1,31)=    6.94789d0
    smallc(2,31)=   -1.14343d0
    smallc(3,31)=    0.35843d0
    smallc(1,32)=    8.71272d0
    smallc(2,32)=   -2.27583d0
    smallc(3,32)=    2.37222d0
    smallc(1,33)=    5.50892d0
    smallc(2,33)=    0.38643d0
    smallc(3,33)=    4.10640d0
    smallc(1,34)=    3.23140d0
    smallc(2,34)=    0.90630d0
    smallc(3,34)=    5.83701d0
    smallc(1,35)=    3.97118d0
    smallc(2,35)=    4.50982d0
    smallc(3,35)=    3.52700d0
    smallc(1,36)=    4.50123d0
    smallc(2,36)=    7.36725d0
    smallc(3,36)=    3.44631d0
    smallc(1,37)=   -1.07313d0
    smallc(2,37)=    0.67903d0
    smallc(3,37)=    6.56527d0
    smallc(1,38)=    0.30262d0
    smallc(2,38)=   -0.39329d0
    smallc(3,38)=    8.89117d0
    smallc(1,39)=   -4.79969d0
    smallc(2,39)=    1.68568d0
    smallc(3,39)=    4.56178d0
    smallc(1,40)=   -3.06076d0
    smallc(2,40)=    4.01515d0
    smallc(3,40)=    4.60718d0
    smallc(1,41)=   -6.20152d0
    smallc(2,41)=   -1.97259d0
    smallc(3,41)=    2.63721d0
    smallc(1,42)=   -8.59075d0
    smallc(2,42)=   -0.55116d0
    smallc(3,42)=    3.48779d0
    smallc(1,43)=   -3.24481d0
    smallc(2,43)=   -5.17079d0
    smallc(3,43)=    3.32426d0
    smallc(1,44)=   -4.35616d0
    smallc(2,44)=   -5.47084d0
    smallc(3,44)=    0.65457d0
    smallc(1,45)=   -0.08946d0
    smallc(2,45)=   -3.49661d0
    smallc(3,45)=    5.77323d0
    smallc(1,46)=   -0.37321d0
    smallc(2,46)=   -6.28524d0
    smallc(3,46)=    6.54489d0
    smallc(1,47)=   -3.20767d0
    smallc(2,47)=   -0.91185d0
    smallc(3,47)=   -5.85446d0
    smallc(1,48)=   -5.48518d0
    smallc(2,48)=   -0.39198d0
    smallc(3,48)=   -4.12384d0
    smallc(1,49)=   -6.92415d0
    smallc(2,49)=    1.13788d0
    smallc(3,49)=   -0.37587d0
    smallc(1,50)=   -8.68898d0
    smallc(2,50)=    2.27028d0
    smallc(3,50)=   -2.38967d0
    smallc(1,51)=   -3.94744d0
    smallc(2,51)=   -4.51537d0
    smallc(3,51)=   -3.54444d0
    smallc(1,52)=   -4.47749d0
    smallc(2,52)=   -7.37280d0
    smallc(3,52)=   -3.46376d0
    smallc(1,53)=   -0.05385d0
    smallc(2,53)=   -6.44310d0
    smallc(3,53)=   -2.75819d0
    smallc(1,54)=    2.40286d0
    smallc(2,54)=   -6.47979d0
    smallc(3,54)=   -1.20400d0
    smallc(1,55)=    0.07759d0
    smallc(2,55)=    6.43755d0
    smallc(3,55)=    2.74075d0
    smallc(1,56)=   -2.37912d0
    smallc(2,56)=    6.47424d0
    smallc(3,56)=    1.18655d0
    smallc(1,57)=    4.83460d0
    smallc(2,57)=   -4.66786d0
    smallc(3,57)=    2.00789d0
    smallc(1,58)=    2.35681d0
    smallc(2,58)=   -5.77211d0
    smallc(3,58)=    3.05359d0
    smallc(1,59)=   -4.84346d0
    smallc(2,59)=    4.68511d0
    smallc(3,59)=   -1.99286d0
    smallc(1,60)=   -4.81635d0
    smallc(2,60)=    7.54743d0
    smallc(3,60)=   -2.50152d0
    !cccccc
    !     calculate total potential inside the cage + h2
    vpot=0.0d0
    DO i=1,nwater
      ! oxygen
      coord(1,1)=smallc(1,i)
      coord(2,1)=smallc(2,i)
      coord(3,1)=smallc(3,i)
      ! first hydrogen of water molecules
      coord(1,2)=smallc(1,nwater+2*i-1)
      coord(2,2)=smallc(2,nwater+2*i-1)
      coord(3,2)=smallc(3,nwater+2*i-1)
      ! second hydrogen of water molecules
      coord(1,3)=smallc(1,nwater+2*i)
      coord(2,3)=smallc(2,nwater+2*i)
      coord(3,3)=smallc(3,nwater+2*i)
      ! H2 molecule at the center of the cage
      coord(1,4)=h2coor(1,1)
      coord(2,4)=h2coor(2,1)
      coord(3,4)=h2coor(3,1)
      coord(1,5)=h2coor(1,2)
      coord(2,5)=h2coor(2,2)
      coord(3,5)=h2coor(3,2)
      ! now call SPC/E subroutine for WATER-H2:
      ! subroutine oh2_h2(zcoord,Eint) 
      CALL oh2_h2_QML(coord,vfit)
      vpot=vpot+vfit
    END DO
  END SUBROUTINE sc_sp_QML
  !***********************************************************************
  SUBROUTINE oh2_h2_QML(zcoord,Eint)
    !***********************************************************************
    !     routine to evaluate to interaction potential between
    !     one H2 molecule and one H2O molecule within the SPC/E model,
    !     see Alavi et al. J. Chem. Phys. 123, 024507 (2005)
    !     and Berendsen et al. J. Phys. Chem. 91, 6269 (1987)
    !
    !     distances in [a.u.]
    !     angles    in [radians]
    !     Eint      in [cm-1]
    !***********************************************************************
    real(kind=Rkind), intent(in)  :: zcoord(3,5)
    real(kind=Rkind), intent(out) :: Eint

    ! zcoord(i,1) = Oxygen coords
    ! zcoord(i,2) = H atom coords
    ! zcoord(i,3) = H atom coords
    ! zcoord(i,4) = H2 first atom coords
    ! zcoord(i,5) = H2 second atom coords

    ! SPC/E parameters
    real(kind=Rkind), parameter :: sigma=3.10134d0
    real(kind=Rkind), parameter :: epsi=0.430624d0
    real(kind=Rkind), parameter :: hwater=0.4238d0
    real(kind=Rkind), parameter :: oxygen=-0.8476d0
    real(kind=Rkind), parameter :: hcm=-0.9864d0
    real(kind=Rkind), parameter :: h2=0.4932d0
 
    real(kind=Rkind) :: R,x(3,6),d(9)
    real(kind=Rkind) :: v_coul,v_lj

    integer :: j

    d(:) = ZERO

    ! x(i,1) are the coordinates of the Oxygen atom
    x(1,1)=zcoord(1,1)
    x(2,1)=zcoord(2,1)
    x(3,1)=zcoord(3,1)
    ! H atoms in WATER molecule
    x(1,2)=zcoord(1,2)
    x(2,2)=zcoord(2,2)
    x(3,2)=zcoord(3,2)
    x(1,3)=zcoord(1,3)
    x(2,3)=zcoord(2,3)
    x(3,3)=zcoord(3,3)

    x(1,4)=zcoord(1,4)
    x(2,4)=zcoord(2,4)
    x(3,4)=zcoord(3,4)
    x(1,5)=zcoord(1,5)
    x(2,5)=zcoord(2,5)
    x(3,5)=zcoord(3,5)

    ! x(i,6) are the coordinates of the H_2 center of mass
    x(1,6)=(zcoord(1,4)+zcoord(1,5))/2.0d0
    x(2,6)=(zcoord(2,4)+zcoord(2,5))/2.0d0
    x(3,6)=(zcoord(3,4)+zcoord(3,5))/2.0d0

    R=0.0d0
    DO j=1,3
      R=R+(x(j,1)-x(j,6))**2
    END DO
    R=dsqrt(R)
    DO j=1,3
      d(1)=d(1)+(x(j,1)-x(j,4))**2
      d(2)=d(2)+(x(j,1)-x(j,5))**2
      d(3)=d(3)+(x(j,1)-x(j,6))**2
      d(4)=d(4)+(x(j,2)-x(j,4))**2
      d(5)=d(5)+(x(j,2)-x(j,5))**2
      d(6)=d(6)+(x(j,2)-x(j,6))**2
      d(7)=d(7)+(x(j,3)-x(j,4))**2
      d(8)=d(8)+(x(j,3)-x(j,5))**2
      d(9)=d(9)+(x(j,3)-x(j,6))**2
    END DO
    d = sqrt(d)

    v_lj=(sigma/(R*ang2au))**12-(sigma/(R*ang2au))**6
    !v_lj now is kcalmol-1: below
    v_lj=4.0d0*epsi*v_lj
    !v_lj now is cm-1: below
    v_lj=v_lj*calkcm
    !v_coul now is in a.u.
    v_coul = oxygen*h2/d(1)+h2*oxygen/d(2)+hcm*oxygen/d(3)+  &
             hwater*h2/d(4)+hwater*h2/d(5)+hwater*hcm/d(6)+  &
             hwater*h2/d(7)+hwater*h2/d(8)+hwater*hcm/d(9)

    !v_coul now is in cm-1: below
    v_coul=v_coul*au2cm
    !Vpot in cm-1: below
    Eint= v_lj + v_coul

  END SUBROUTINE oh2_h2_QML

END MODULE QML_H2_H2On_m
