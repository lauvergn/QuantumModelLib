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

!> @brief Module which makes the initialization, calculation of the Phenol potential (value, gradient and hessian).
!> @brief Reference: Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218
!> @brief The potential is a 2D (R,theta) with 3 electronic diabatic surfaces (S_0, Pi-Sigma* and Pi-Pi*)
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
MODULE QML_Phenol_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  USE QML_Morse_m
  USE QML_Sigmoid_m
  IMPLICIT NONE

  PRIVATE


!> @brief Derived type in which the parameters of the Phenol potential are set-up.
!> @brief Reference: Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218
!> @brief the parameter names are from the previous reference and are taken from Eqs 4-23 and tables I-IV.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
  TYPE, EXTENDS (QML_Empty_t) :: QML_Phenol_t
     PRIVATE

     ! V(1,1) term
     TYPE (QML_Morse_t)   :: v10
     TYPE (QML_Sigmoid_t) :: v11

     ! V(3,3) term
     TYPE (QML_Morse_t)   :: v30
     TYPE (QML_Sigmoid_t) :: v31
     real(kind=Rkind)      :: a30=4.85842_Rkind ! eV

     ! V(2,2) term
     TYPE (QML_Morse_t)   :: v201
     real(kind=Rkind)      :: b204=5.50696_Rkind ! eV

     real(kind=Rkind)      :: b205=4.70601_Rkind ! eV
     real(kind=Rkind)      :: b206=2.49826_Rkind ! Angs^-1
     real(kind=Rkind)      :: b207=0.988188_Rkind ! Angs
     real(kind=Rkind)      :: b208=3.3257_Rkind ! eV

     TYPE (QML_Sigmoid_t) :: v211,v212,v221,v222
     real(kind=Rkind)      :: b217=-0.00055_Rkind ! eV

     real(kind=Rkind)      :: X20=0.326432_Rkind ! eV^2
     real(kind=Rkind)      :: X21=0.021105_Rkind ! eV^2
     real(kind=Rkind)      :: X22=0._Rkind ! eV^2


     ! V(1,3), and V(1,2) terms
     TYPE (QML_Sigmoid_t) :: lambda12,lambda23


      ! The metric tensor of Tnum with rigid_type=100 from B3LYP/6-31G** of the ground state (in au)
     real (kind=Rkind), PUBLIC :: G_RR    = 0.0005786177_Rkind
     real (kind=Rkind), PUBLIC :: G_ThTh  = 0.0002550307_Rkind
  CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_Phenol
    PROCEDURE :: Write_QModel    => Write_QML_Phenol
  END TYPE QML_Phenol_t

  PUBLIC :: QML_Phenol_t,Init_QML_Phenol


CONTAINS
!> @brief Subroutine which makes the initialization of the phenol parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PhenolPot          TYPE(PhenolPot_t):   derived type in which the parameters are set-up.
  FUNCTION Init_QML_Phenol(QModel_in,read_param,nio_param_file) RESULT(QModel)
  USE QDUtil_m,         ONLY : Identity_Mat
  IMPLICIT NONE

    TYPE (QML_Phenol_t)               :: QModel

    TYPE(QML_Empty_t),    intent(in)   :: QModel_in ! variable to transfer info to the init
    logical,              intent(in)   :: read_param
    integer,              intent(in)   :: nio_param_file

    real (kind=Rkind) :: Req
    real (kind=Rkind) :: a0      = 0.52917720835354106_Rkind
    real (kind=Rkind) :: auTOeV  = 27.211384_Rkind
!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Phenol'
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
    QModel%pot_name = 'phenol'
    QModel%ndim     = 2
    QModel%nsurf    = 3
    flush(out_unit)

    ! Warning the parameters are given as in the publication.
    !   Therefore, the OH distance R(=Q(1)) is in Angstrom and the energy is in eV.

    ! V(1,1) term
    !De1=4.26302 eV r1=0.96994 Å a1=2.66021 Å−1
    Req = 0.96994_Rkind
    CALL Init0_QML_Morse(QModel%v10,D=4.26302_Rkind,a=2.66021_Rkind,    &
                          req=0.96994_Rkind,model_name='morse-v10')
    !A1=0.27037 eV A2=1.96606 Å A3=0.685264 Å
    CALL Init0_QML_Sigmoid(QModel%v11,A=0.27037_Rkind,B=1.96606_Rkind,  &
                            C=0.685264_Rkind,e=-ONE,model_name='sigmoid-v11')

    ! V(2,2) term
    !B201=0.192205 eV B202=5.67356 Å−1 B203=1.03171 Å
    CALL Init0_QML_Morse(QModel%v201,D=0.192205_Rkind,a=5.67356_Rkind,  &
                          req=1.03171_Rkind,model_name='morse-v201')
    !the exp term is given with the coef's: a205,a206,a207,a208

    !B211=−0.2902 eV B212=2.05715 Å B213=1.01574 Å
    CALL Init0_QML_Sigmoid(QModel%v211,A=-0.2902_Rkind,B=2.05715_Rkind, &
                            C=1.01574_Rkind,e=-ONE,model_name='sigmoid-v211')
    !B214=−73.329 eV B215=1.48285 Å B216=−0.1111 Å
    CALL Init0_QML_Sigmoid(QModel%v212,A=-73.329_Rkind,B=1.48285_Rkind, &
                            C=-0.1111_Rkind,e=-ONE,model_name='sigmoid-v212')
    !B221=27.3756 eV B222=1.66881 Å B223=0.20557 Å
    CALL Init0_QML_Sigmoid(QModel%v221,A=27.3756_Rkind,B=1.66881_Rkind, &
                            C=0.20557_Rkind,e=ONE,model_name='sigmoid-v221')
    !B224=0.35567 Å B225=1.43492 eV B226=0.56968 Å (unit problem between B224 and B225)
    CALL Init0_QML_Sigmoid(QModel%v222,A=0.35567_Rkind,B=1.43492_Rkind, &
                            C=0.56968_Rkind,e=-ONE,model_name='sigmoid-v222')

    ! V(3,3) term
    !De3=4.47382 eV r3=0.96304 Å a3=2.38671 Å−1 a30=4.85842 eV
    CALL Init0_QML_Morse(QModel%v30,D=4.47382_Rkind,a=2.38671_Rkind,    &
                          req=0.96304_Rkind,model_name='morse-v30')
    !C1=0.110336 eV C2=1.21724 Å C3=0.06778 Å̊
    CALL Init0_QML_Sigmoid(QModel%v31,A=0.110336_Rkind,B=1.21724_Rkind, &
                            C=0.06778_Rkind,e=-ONE,model_name='sigmoid-v31')

    ! V(2,3), and V(1,2) terms
    !l12,max=1.47613 eV d12=1.96984 Å l12=0.494373 Å
    CALL Init0_QML_Sigmoid(QModel%lambda12,A=1.47613_Rkind,B=1.96984_Rkind,&
                            C=0.494373_Rkind,e=-ONE,model_name='sigmoid-lambda12')
    !l23,max=0.327204 eV d23=1.22594 Å l23=0.0700604 Å
    CALL Init0_QML_Sigmoid(QModel%lambda23,A=0.327204_Rkind,B=1.22594_Rkind,&
                            C=0.0700604_Rkind,e=-ONE,model_name='sigmoid-lambda23')

    IF (debug) THEN
      IF (QModel%PubliUnit) THEN
        write(out_unit,*) 'init Q0 of Phenol [Angs,Rad]'
      ELSE
        write(out_unit,*) 'init Q0 of Phenol [Bohr,Rad]'
      END IF
    END IF
    QModel%Q0 = [Req,ZERO]
    IF (.NOT. QModel%PubliUnit) QModel%Q0(1) = QModel%Q0(1) /a0


    IF (debug) write(out_unit,*) 'init d0GGdef of Phenol [au]'
    QModel%d0GGdef      = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,1) = QModel%G_RR
    QModel%d0GGdef(2,2) = QModel%G_ThTh


    IF (QModel%PubliUnit) THEN
      write(out_unit,*) 'PubliUnit=.TRUE.,  Q:[Angs,Rad], Energy: [eV]'
    ELSE
      write(out_unit,*) 'PubliUnit=.FALSE., Q:[Bohr,Rad], Energy: [Hartree]'
    END IF

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Phenol
!> @brief Subroutine wich prints the Phenol potential parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel        CLASS(PhenolPot_t):   derived type with the Phenol potential parameters.
!! @param nio           integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_Phenol(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_Phenol_t),    intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'Phenol parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) '---- WARNNING ---------------------------'
    write(nio,*) 'The parameters are given as in the publication: '
    write(nio,*) '  Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, ...'
    write(nio,*) '  .... J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218'
    write(nio,*) 'with the OH distance R=Q(1)     in Angstrom ...'
    write(nio,*) '     the OH angle    theta=Q(2) in Radian   ...'
    write(nio,*) '     and the energy in eV.'
    write(nio,*)
    write(nio,*) 'Diabatic Potential Values (in eV) at: R=1.2 Angs and theta=0.2 Radian'
    write(nio,*) '1        0.912491       0.280793       0.000000'
    write(nio,*) '2        0.280793       5.437178       0.044016'
    write(nio,*) '3        0.000000       0.044016       5.698608'
    write(nio,*) 'Adiabatic Potential Values (in eV) at: R=1.2 Angs and theta=0.2 Radian'
    write(nio,*) '1        0.895130       0.000000       0.000000'
    write(nio,*) '2        0.000000       5.446870       0.000000'
    write(nio,*) '3        0.000000       0.000000       5.706276'
    write(nio,*)
    write(nio,*) 'end Phenol parameters'
    write(nio,*) 'Phenol current parameters'

    write(nio,*) 'PubliUnit: ',QModel%PubliUnit
    write(nio,*)
    write(nio,*)
    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(1,1):'
    CALL Write_QML_Morse(QModel%v10,nio)
    CALL Write_QML_Sigmoid(QModel%v11,nio)

    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(2,2):'
    CALL Write_QML_Morse(QModel%v201,nio)
    write(nio,*) ' b204:',QModel%b204
    write(nio,*) ' v202=B205*exp(-B206*(R-B207)) + B208'
    write(nio,*) ' b205...b208:',QModel%b205,QModel%b206,QModel%b207,QModel%b208
    CALL Write_QML_Sigmoid(QModel%v211,nio)
    write(nio,*) ' b217:',QModel%b217
    CALL Write_QML_Sigmoid(QModel%v212,nio)
    CALL Write_QML_Sigmoid(QModel%v221,nio)
    CALL Write_QML_Sigmoid(QModel%v222,nio)

    write(nio,*) ' X20,X21,X22:',QModel%X20,QModel%X21,QModel%X22


    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(3,3):'
    write(nio,*) ' a30:',QModel%a30

    CALL Write_QML_Morse(QModel%v30,nio)
    CALL Write_QML_Sigmoid(QModel%v31,nio)

    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(1,2):'
    CALL Write_QML_Sigmoid(QModel%lambda12,nio)
    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(2,3):'
    CALL Write_QML_Sigmoid(QModel%lambda23,nio)

    write(nio,*) 'end Phenol current parameters'

  END SUBROUTINE Write_QML_Phenol

!> @brief Subroutine wich calculates the Phenol potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                table of two values for which the potential is calculated (R,theta)
!! @param QModel        TYPE(PhenolPot_t):  derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Phenol(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_Phenol_t),    intent(in) :: QModel

    TYPE (dnS_t),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)     :: dnQ(:) ! R and th
    integer,             intent(in)     :: nderiv

    TYPE (dnS_t)  :: dnR,dnth
    TYPE (dnS_t)  :: v10R,v11R,v11th
    TYPE (dnS_t)  :: v30R,v31R,v31th
    TYPE (dnS_t)  :: lambda12R,lambda23R

    TYPE (dnS_t)  :: v20pR,v20mR,v201R,v202R
    TYPE (dnS_t)  :: v21pR,v21mR,v211R,v212R
    TYPE (dnS_t)  :: v22pR,v22mR,v221R,v222R
    TYPE (dnS_t)  :: v20R,v21R,v22R,v21th,v22th

    integer      :: i,j

    real (kind=Rkind) :: a0      = 0.52917720835354106_Rkind
    real (kind=Rkind) :: auTOeV  = 27.211384_Rkind


   !write(out_unit,*) 'phenol pot in:'

   dnR     = dnQ(1)
   dnth    = dnQ(2)

   IF (.NOT. QModel%PubliUnit) THEN
      dnR = a0*dnR ! to convert the bhor into Angstrom
   END IF

   !--------------------------------------------------------------------
   ! for V(1,1): first diabatic state
   !write(out_unit,*) 'morse:'
   v10R = QML_dnMorse(dnR,QModel%v10)
   !CALL Write_dnS(v10R,6)
   !write(out_unit,*) 'sigmoid:'
   v11R = QML_dnSigmoid(dnR,QModel%v11)
   !CALL Write_dnS(v11R,6)

   !write(out_unit,*) 'f(th):'
   v11th = ONE-cos(dnth+dnth)
   !CALL Write_dnS(v11th,6)

   Mat_OF_PotDia(1,1) = v10R+v11R*v11th

   CALL dealloc_dnS(v10R)
   CALL dealloc_dnS(v11R)
   CALL dealloc_dnS(v11th)
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! for V(2,2): 2d diabatic state
   !v201R = QML_dnMorse(dnR,QModel%v201) + QModel%B204 ! pb with ifort19
   v201R = QML_dnMorse(dnR,QModel%v201)
   v201R = v201R + QModel%B204

   v202R = QModel%B205*exp(-QModel%B206*(dnR-QModel%B207)) + QModel%B208

   v20pR = v201R + v202R

   v20mR = v201R - v202R

   v20R  = HALF*(v20pR - (v20mR**TWO + QModel%X20)**HALF)


   !write(out_unit,*) 'sigmoid:'
   v211R = QML_dnSigmoid(dnR,QModel%v211)
   
   !v212R = QML_dnSigmoid(dnR,QModel%v212) + QModel%B217 ! pb with ifort19
   v212R = QML_dnSigmoid(dnR,QModel%v212)
   v212R = v212R + QModel%B217

   v21pR = v211R + v212R
   v21mR = v211R - v212R
   v21R  = HALF * (v21pR + (v21mR**TWO + QModel%X21)**HALF)

   v221R = QML_dnSigmoid(dnR,QModel%v221)
   v222R = QML_dnSigmoid(dnR,QModel%v222)
   v22pR = v221R + v222R
   v22mR = v221R - v222R
   v22R  = HALF * (v22pR - sqrt(v22mR**TWO + QModel%X22) )

   v21th = ONE-cos(dnth+dnth)
   v22th = v21th*v21th

   Mat_OF_PotDia(2,2) = v20R+v21R*v21th+v22R*v22th


   CALL dealloc_dnS(v20R)
   CALL dealloc_dnS(v20pR)
   CALL dealloc_dnS(v20mR)
   CALL dealloc_dnS(v201R)
   CALL dealloc_dnS(v202R)

   CALL dealloc_dnS(v21R)
   CALL dealloc_dnS(v21th)
   CALL dealloc_dnS(v21pR)
   CALL dealloc_dnS(v21mR)
   CALL dealloc_dnS(v211R)
   CALL dealloc_dnS(v212R)

   CALL dealloc_dnS(v22R)
   CALL dealloc_dnS(v22th)
   CALL dealloc_dnS(v22pR)
   CALL dealloc_dnS(v22mR)
   CALL dealloc_dnS(v221R)
   CALL dealloc_dnS(v222R)
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! for V(3,3): 3d diabatic state
   !write(out_unit,*) 'morse:'
   !v30R = QML_dnMorse(dnR,QModel%v30) + QModel%a30 ! with ifort19 this doesn't work, it has to be split
   v30R = QML_dnMorse(dnR,QModel%v30)
   v30R = v30R + QModel%a30

   !CALL Write_dnMat(v30R,6)
   !write(out_unit,*) 'sigmoid:'
   v31R = QML_dnSigmoid(dnR,QModel%v31)
   !CALL Write_dnMat(v31R,6)

   !write(out_unit,*) 'f(th):'
   v31th = ONE-cos(dnth+dnth)
   !CALL Write_dnS(v11th,6)

   Mat_OF_PotDia(3,3) = v30R+v31R*v31th

   !write(out_unit,*) 'phenol pot diabatic:',nderiv
   !CALL Write_dnMat(PotVal,6)

   CALL dealloc_dnS(v30R)
   CALL dealloc_dnS(v31R)
   CALL dealloc_dnS(v31th)
   !--------------------------------------------------------------------


   !--------------------------------------------------------------------
   lambda12R = QML_dnSigmoid(dnR,QModel%lambda12) * sin(dnth)
   !CALL Write_dnS(lambda12R,6)

   Mat_OF_PotDia(1,2) = lambda12R
   Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)



   lambda23R = QML_dnSigmoid(dnR,QModel%lambda23) * sin(dnth)

   Mat_OF_PotDia(2,3) = lambda23R
   Mat_OF_PotDia(3,2) = Mat_OF_PotDia(2,3)

   Mat_OF_PotDia(1,3) = ZERO
   Mat_OF_PotDia(3,1) = ZERO


   CALL dealloc_dnS(lambda12R)
   CALL dealloc_dnS(lambda23R)
   !--------------------------------------------------------------------


   CALL dealloc_dnS(dnth)
   CALL dealloc_dnS(dnR)

   IF (.NOT. QModel%PubliUnit) THEN ! to convert the eV into Hartree
     DO i=1,3
     DO j=1,3
       Mat_OF_PotDia(j,i) = Mat_OF_PotDia(j,i) * (ONE/auTOev)
     END DO
     END DO
   END IF


   !write(out_unit,*) 'phenol pot diabatic:',nderiv
   !CALL Write_dnMat(PotVal,6)
   !write(out_unit,*)
   !flush(out_unit)

  END SUBROUTINE EvalPot_QML_Phenol

END MODULE QML_Phenol_m
