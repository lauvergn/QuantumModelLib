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
!    Copyright 2016  David LAUVERGNAT, Félix MOUHAT and Liang LIANG
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
MODULE mod_PhenolPot
  USE mod_NumParameters
  USE mod_MorsePot
  USE mod_SigmoidPot

  IMPLICIT NONE

!> @brief Derived type in which the parameters of the Phenol potential are set-up.
!> @brief Reference: Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218
!> @brief the parameter names are from the previous reference and are taken from Eqs 4-23 and tables I-IV.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
  TYPE Param_Phenol
     PRIVATE
     logical              :: PubliUnit = .FALSE. ! when PubliUnit=.TRUE., the units (Angstrom and Ev) are used. Default (atomic unit)

     ! Warning the parameters are given as in the publication.
     !   Therefore, the OH distance R(=Q(1)) is in Angstrom and the energy is in eV.

     ! V(1,1) term
     TYPE (Param_Morse)   :: v10
     TYPE (Param_Sigmoid) :: v11

     ! V(3,3) term
     TYPE (Param_Morse)   :: v30
     TYPE (Param_Sigmoid) :: v31
     real(kind=Rkind)     :: a30=4.85842_Rkind ! eV

     ! V(2,2) term
     TYPE (Param_Morse)   :: v201
     real(kind=Rkind)     :: b204=5.50696_Rkind ! eV

     real(kind=Rkind)     :: b205=4.70601_Rkind ! eV
     real(kind=Rkind)     :: b206=2.49826_Rkind ! A^-1
     real(kind=Rkind)     :: b207=0.988188_Rkind ! A
     real(kind=Rkind)     :: b208=3.3257_Rkind ! eV

     TYPE (Param_Sigmoid) :: v211,v212,v221,v222
     real(kind=Rkind)     :: b217=-0.00055_Rkind ! eV

     real(kind=Rkind)     :: X20=0.326432_Rkind ! eV^2
     real(kind=Rkind)     :: X21=0.021105_Rkind ! eV^2
     real(kind=Rkind)     :: X22=0._Rkind ! eV^2


     ! V(1,3), and V(1,2) terms
     TYPE (Param_Sigmoid) :: lambda12,lambda13


      ! The metric tensor of Tnum with rigid_type=100 from B3LYP/6-31G** of the ground state (in au)
     real (kind=Rkind), PUBLIC :: G_RR    = 0.0005786177_Rkind
     real (kind=Rkind), PUBLIC :: G_ThTh  = 0.0002550307_Rkind

  END TYPE Param_Phenol


CONTAINS
!> @brief Subroutine which makes the initialization of the phenol parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Phenol        TYPE(Param_Phenol):   derived type in which the parameters are set-up.
!! @param PubliUnit          logical (optional):   when PubliUnit=.TRUE., the units (Angstrom and eV) are used. Default (atomic unit).
  SUBROUTINE Init_PhenolPot(Para_Phenol,PubliUnit)
    TYPE (Param_Phenol),      intent(inout)   :: Para_Phenol
     logical, optional,       intent(in)      :: PubliUnit

     IF (present(PubliUnit)) Para_Phenol%PubliUnit = PubliUnit

    ! V(1,1) term
    !De1=4.26302 eV r1=0.96994 Å a1=2.66021 Å−1
    CALL Init_MorsePot(Para_Phenol%v10,D=4.26302_Rkind,a=2.66021_Rkind,req=0.96994_Rkind)
    !A1=0.27037 eV A2=1.96606 Å A3=0.685264 Å
    CALL Init_SigmoidPot(Para_Phenol%v11,nio=5,read_param=.FALSE.,      &
                         A=0.27037_Rkind,B=1.96606_Rkind,C=0.685264_Rkind,e=-ONE)



    ! V(2,2) term
    !B201=0.192205 eV B202=5.67356 Å−1 B203=1.03171 Å
    CALL Init_MorsePot(Para_Phenol%v201,D=0.192205_Rkind,a=5.67356_Rkind,req=1.03171_Rkind)
    !the exp term is given with the coef's: a205,a206,a207,a208

    !B211=−0.2902 eV B212=2.05715 Å B213=1.01574 Å
    CALL Init_SigmoidPot(Para_Phenol%v211,nio=5,read_param=.FALSE.,     &
                         A=-0.2902_Rkind,B=2.05715_Rkind,C=1.01574_Rkind,e=-ONE)
    !B214=−73.329 eV B215=1.48285 Å B216=−0.1111 Å
    CALL Init_SigmoidPot(Para_Phenol%v212,nio=5,read_param=.FALSE.,     &
                         A=-73.329_Rkind,B=1.48285_Rkind,C=-0.1111_Rkind,e=-ONE)
    !B221=27.3756 eV B222=1.66881 Å B223=0.20557 Å
    CALL Init_SigmoidPot(Para_Phenol%v221,nio=5,read_param=.FALSE.,     &
                         A=27.3756_Rkind,B=1.66881_Rkind,C=0.20557_Rkind,e=ONE)
    !B224=0.35567 Å B225=1.43492 eV B226=0.56968 Å (unit problem between B224 and B225)
    CALL Init_SigmoidPot(Para_Phenol%v222,nio=5,read_param=.FALSE.,     &
                         A=0.35567_Rkind,B=1.43492_Rkind,C=0.56968_Rkind,e=-ONE)

    ! V(3,3) term
    !De3=4.47382 eV r3=0.96304 Å a3=2.38671 Å−1 a30=4.85842 eV
    CALL Init_MorsePot(Para_Phenol%v30,D=4.47382_Rkind,a=2.38671_Rkind,req=0.96304_Rkind)
    !C1=0.110336 eV C2=1.21724 Å C3=0.06778 Å̊
    CALL Init_SigmoidPot(Para_Phenol%v31,nio=5,read_param=.FALSE.,      &
                         A=0.110336_Rkind,B=1.21724_Rkind,C=0.06778_Rkind,e=-ONE)

    ! V(1,3), and V(1,2) terms
    !l12,max=1.47613 eV d12=1.96984 Å l12=0.494373 Å
    CALL Init_SigmoidPot(Para_Phenol%lambda12,nio=5,read_param=.FALSE., &
                         A=1.47613_Rkind,B=1.96984_Rkind,C=0.494373_Rkind,e=-ONE)
    !l23,max=0.327204 eV d23=1.22594 Å l23=0.0700604 Å
    CALL Init_SigmoidPot(Para_Phenol%lambda13,nio=5,read_param=.FALSE., &
                         A=0.327204_Rkind,B=1.22594_Rkind,C=0.0700604_Rkind,e=-ONE)

    IF (Para_Phenol%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Angs,Rad], Energy: [eV]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Rad], Energy: [Hartree]'
    END IF

  END SUBROUTINE Init_PhenolPot
!> @brief Subroutine wich prints the Phenol potential parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer:              file unit to print the parameters.
  SUBROUTINE Write0_PhenolPot(nio)
    integer, intent(in) :: nio

    write(nio,*) 'Phenol parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) '---- WARNNING ---------------------------'
    write(nio,*) 'The parameters are given as in the publication: '
    write(nio,*) '  Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, ...'
    write(nio,*) '  .... J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218'
    write(nio,*) 'with the OH distance R=Q(1)     in Angstrom ...'
    write(nio,*) '     the OH angle    thet(=Q(2) in Radian   ...'
    write(nio,*) '     and the energy in eV.'
    write(nio,*)
    write(nio,*) 'Diabatic Potential Values (in eV) at: R=1.2 Angs and theta=0.2 Radian'
    write(nio,*) '1        0.912491       0.280793       0.044016'
    write(nio,*) '2        0.280793       5.437178       0.000000'
    write(nio,*) '3        0.044016       0.000000       5.698608'
    write(nio,*) 'Adiabatic Potential Values (in eV) at: R=1.2 Angs and theta=0.2 Radian'
    write(nio,*) '1        0.894730       0.000000       0.000000'
    write(nio,*) '2        0.000000       5.454506       0.000000'
    write(nio,*) '3        0.000000       0.000000       5.699040'
    write(nio,*)
    write(nio,*) 'end Phenol parameters'

  END SUBROUTINE Write0_PhenolPot
!> @brief Subroutine wich prints the Phenol potential parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Phenol        TYPE(Param_Phenol):   derived type with the Phenol potential parameters.
!! @param nio                integer:              file unit to print the parameters.
  SUBROUTINE Write_PhenolPot(Para_Phenol,nio)
    TYPE (Param_Phenol), intent(in) :: Para_Phenol
    integer, intent(in) :: nio

    write(nio,*) 'Phenol current parameters'


    write(nio,*) 'PubliUnit: ',Para_Phenol%PubliUnit
    write(nio,*)
    write(nio,*)
    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(1,1):'
    CALL Write_MorsePot(Para_Phenol%v10,nio)
    CALL Write_SigmoidPot(Para_Phenol%v11,nio)

    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(2,2):'
    CALL Write_MorsePot(Para_Phenol%v201,nio)
    write(nio,*) ' b204:',Para_Phenol%b204
    write(nio,*) ' v202=B205*exp(-B206*(R-B207)) + B208'
    write(nio,*) ' b205...b208:',Para_Phenol%b205,Para_Phenol%b206,Para_Phenol%b207,Para_Phenol%b208
    CALL Write_SigmoidPot(Para_Phenol%v211,nio)
    write(nio,*) ' b217:',Para_Phenol%b217
    CALL Write_SigmoidPot(Para_Phenol%v212,nio)
    CALL Write_SigmoidPot(Para_Phenol%v221,nio)
    CALL Write_SigmoidPot(Para_Phenol%v222,nio)

    write(nio,*) ' X20,X21,X22:',Para_Phenol%X20,Para_Phenol%X21,Para_Phenol%X22


    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(3,3):'
    write(nio,*) ' a30:',Para_Phenol%a30

    CALL Write_MorsePot(Para_Phenol%v30,nio)
    CALL Write_SigmoidPot(Para_Phenol%v31,nio)

    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(1,2):'
    CALL Write_SigmoidPot(Para_Phenol%lambda12,nio)
    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(1,3):'
    CALL Write_SigmoidPot(Para_Phenol%lambda13,nio)

    write(nio,*) 'end Phenol current parameters'

  END SUBROUTINE Write_PhenolPot

!> @brief Subroutine wich calculates the Phenol potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                table of two values for which the potential is calculated (R,theta)
!! @param Para_Phenol        TYPE(Param_Phenol):  derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_PhenolPot(Mat_OF_PotDia,dnQ,Para_Phenol,nderiv)
    USE mod_dnSca

    TYPE (Param_Phenol), intent(in)     :: Para_Phenol
    TYPE(dnSca),         intent(inout)  :: Mat_OF_PotDia(:,:)
    TYPE(dnSca),         intent(in)     :: dnQ(:) ! R and th
    integer,             intent(in)     :: nderiv

    TYPE(dnSca)  :: dnR,dnth
    TYPE(dnSca)  :: v10R,v11R,v11th
    TYPE(dnSca)  :: v30R,v31R,v31th
    TYPE(dnSca)  :: lambda12R,lambda13R

    TYPE(dnSca)  :: v20pR,v20mR,v201R,v202R
    TYPE(dnSca)  :: v21pR,v21mR,v211R,v212R
    TYPE(dnSca)  :: v22pR,v22mR,v221R,v222R
    TYPE(dnSca)  :: v20R,v21R,v22R,v21th,v22th

    integer      :: i,j

    real (kind=Rkind) :: a0      = 0.52917720835354106_Rkind
    real (kind=Rkind) :: auTOeV  = 27.211384_Rkind


   !write(out_unitp,*) 'phenol pot in:'

   dnR     = dnQ(1)
   dnth    = dnQ(2)

   IF (.NOT. Para_Phenol%PubliUnit) THEN
      dnR = a0*dnR ! to convert the bhor into Angstrom
   END IF

   !--------------------------------------------------------------------
   ! for V(1,1): first diabatic state
   !write(out_unitp,*) 'morse:'
   v10R = dnMorse(dnR,Para_Phenol%v10)
   !CALL Write_dnSca(v10R,6)
   !write(out_unitp,*) 'sigmoid:'
   v11R = dnScaigmoid(dnR,Para_Phenol%v11)
   !CALL Write_dnSca(v11R,6)

   !write(out_unitp,*) 'f(th):'
   v11th = ONE-cos(dnth+dnth)
   !CALL Write_dnSca(v11th,6)

   Mat_OF_PotDia(1,1) = v10R+v11R*v11th

   CALL dealloc_dnSca(v10R)
   CALL dealloc_dnSca(v11R)
   CALL dealloc_dnSca(v11th)
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! for V(2,2): 2d diabatic state
   v201R = dnMorse(dnR,Para_Phenol%v201) + Para_Phenol%B204
   v202R = Para_Phenol%B205*exp(-Para_Phenol%B206*(dnR-Para_Phenol%B207)) + Para_Phenol%B208
   v20pR = v201R + v202R
   v20mR = v201R - v202R
   v20R  = HALF*(v20pR - (v20mR**TWO + Para_Phenol%X20)**HALF)


   !write(out_unitp,*) 'sigmoid:'
   v211R = dnScaigmoid(dnR,Para_Phenol%v211)
   v212R = dnScaigmoid(dnR,Para_Phenol%v212) + Para_Phenol%B217
   v21pR = v211R + v212R
   v21mR = v211R - v212R
   v21R  = HALF * (v21pR + (v21mR**TWO + Para_Phenol%X21)**HALF)

   v221R = dnScaigmoid(dnR,Para_Phenol%v221)
   v222R = dnScaigmoid(dnR,Para_Phenol%v222)
   v22pR = v221R + v222R
   v22mR = v221R - v222R
   v22R  = HALF * (v22pR - sqrt(v22mR**TWO + Para_Phenol%X22) )

   v21th = ONE-cos(dnth+dnth)
   v22th = v21th*v21th

   Mat_OF_PotDia(2,2) = v20R+v21R*v21th+v22R*v22th


   CALL dealloc_dnSca(v20R)
   CALL dealloc_dnSca(v20pR)
   CALL dealloc_dnSca(v20mR)
   CALL dealloc_dnSca(v201R)
   CALL dealloc_dnSca(v202R)

   CALL dealloc_dnSca(v21R)
   CALL dealloc_dnSca(v21th)
   CALL dealloc_dnSca(v21pR)
   CALL dealloc_dnSca(v21mR)
   CALL dealloc_dnSca(v211R)
   CALL dealloc_dnSca(v212R)

   CALL dealloc_dnSca(v22R)
   CALL dealloc_dnSca(v22th)
   CALL dealloc_dnSca(v22pR)
   CALL dealloc_dnSca(v22mR)
   CALL dealloc_dnSca(v221R)
   CALL dealloc_dnSca(v222R)
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! for V(3,3): 3d diabatic state
   !write(out_unitp,*) 'morse:'
   v30R = dnMorse(dnR,Para_Phenol%v30) + Para_Phenol%a30
   !CALL Write_dnMatPot(v30R,6)
   !write(out_unitp,*) 'sigmoid:'
   v31R = dnScaigmoid(dnR,Para_Phenol%v31)
   !CALL Write_dnMatPot(v31R,6)

   !write(out_unitp,*) 'f(th):'
   v31th = ONE-cos(dnth+dnth)
   !CALL Write_dnSca(v11th,6)

   Mat_OF_PotDia(3,3) = v30R+v31R*v31th

   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL Write_dnMatPot(PotVal,6)

   CALL dealloc_dnSca(v30R)
   CALL dealloc_dnSca(v31R)
   CALL dealloc_dnSca(v31th)
   !--------------------------------------------------------------------


   !--------------------------------------------------------------------
   lambda12R = dnScaigmoid(dnR,Para_Phenol%lambda12) * sin(dnth)
   !CALL Write_dnSca(lambda12R,6)

   Mat_OF_PotDia(1,2) = lambda12R
   Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)



   lambda13R = dnScaigmoid(dnR,Para_Phenol%lambda13) * sin(dnth)

   Mat_OF_PotDia(1,3) = lambda13R
   Mat_OF_PotDia(3,1) = Mat_OF_PotDia(1,3)

   Mat_OF_PotDia(2,3) = ZERO
   Mat_OF_PotDia(3,2) = ZERO


   CALL dealloc_dnSca(lambda12R)
   CALL dealloc_dnSca(lambda13R)
   !--------------------------------------------------------------------


   CALL dealloc_dnSca(dnth)
   CALL dealloc_dnSca(dnR)

   IF (.NOT. Para_Phenol%PubliUnit) THEN ! to convert the eV into Hartree
     DO i=1,3
     DO j=1,3
       Mat_OF_PotDia(j,i) = Mat_OF_PotDia(j,i) * (ONE/auTOev)
     END DO
     END DO
   END IF


   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL Write_dnMatPot(PotVal,6)
   !write(out_unitp,*)
   !flush(out_unitp)

  END SUBROUTINE eval_PhenolPot

END MODULE mod_PhenolPot
