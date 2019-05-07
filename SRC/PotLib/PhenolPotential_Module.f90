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

  END SUBROUTINE Init_PhenolPot

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

    write(nio,*) 'PubliUnit: ',Para_Phenol%PubliUnit

    write(nio,*)
    write(nio,*) 'Current parameters:'
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

    write(nio,*) 'end Phenol parameters'

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
  SUBROUTINE eval_PhenolPot(PotVal,Q,Para_Phenol,nderiv)
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_Phenol), intent(in)     :: Para_Phenol
    real (kind=Rkind),   intent(in)     :: Q(2)
    TYPE(dnMatPot),      intent(inout)  :: PotVal
    integer,             intent(in)     :: nderiv

    TYPE(dnSca)  :: dnR,dnth
    TYPE(dnSca)  :: v10R,v11R,v11th
    TYPE(dnSca)  :: v30R,v31R,v31th
    TYPE(dnSca)  :: lambda12R,lambda13R

    TYPE(dnSca)  :: v20pR,v20mR,v201R,v202R
    TYPE(dnSca)  :: v21pR,v21mR,v211R,v212R
    TYPE(dnSca)  :: v22pR,v22mR,v221R,v222R
    TYPE(dnSca)  :: v20R,v21R,v22R,v21th,v22th

    real (kind=Rkind) :: a0      = 0.52917720835354106_Rkind
    real (kind=Rkind) :: auTOeV  = 27.211384_Rkind

    IF ( Check_NotAlloc_dnMatPot(PotVal,nderiv) ) THEN
      CALL alloc_dnMatPot(PotVal,nsurf=3,ndim=2,nderiv=nderiv)
    END IF

   !write(out_unitp,*) 'phenol pot in:'
   PotVal = ZERO
   !CALL Write_dnMatPot(PotVal)

   dnR     = init_dnSca(Q(1),ndim=2,nderiv=nderiv,iQ=1) ! to set up the derivatives
   dnth    = init_dnSca(Q(2),ndim=2,nderiv=nderiv,iQ=2) ! to set up the derivatives

   IF (.NOT. Para_Phenol%PubliUnit) dnR = a0*dnR ! to convert the bhor into Angstrom

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

   CALL sub_dnSca_TO_dnMatPot(v10R+v11R*v11th,PotVal,i=1,j=1)

   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL Write_dnMatPot(PotVal,6)

   CALL dealloc_dnSca(v10R)
   CALL dealloc_dnSca(v11R)
   CALL dealloc_dnSca(v11th)

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

   CALL sub_dnSca_TO_dnMatPot(v20R+v21R*v21th+v22R*v22th,PotVal,i=2,j=2)


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

   CALL sub_dnSca_TO_dnMatPot(v30R+v31R*v31th,PotVal,i=3,j=3)

   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL Write_dnMatPot(PotVal,6)

   CALL dealloc_dnSca(v30R)
   CALL dealloc_dnSca(v31R)
   CALL dealloc_dnSca(v31th)


   lambda12R = dnScaigmoid(dnR,Para_Phenol%lambda12) * sin(dnth)
   !CALL Write_dnSca(lambda12R,6)

   CALL sub_dnSca_TO_dnMatPot(lambda12R,PotVal,i=1,j=2)
   CALL sub_dnSca_TO_dnMatPot(lambda12R,PotVal,i=2,j=1)


   lambda13R = dnScaigmoid(dnR,Para_Phenol%lambda13) * sin(dnth)

   CALL sub_dnSca_TO_dnMatPot(lambda13R,PotVal,i=1,j=3)
   CALL sub_dnSca_TO_dnMatPot(lambda13R,PotVal,i=3,j=1)


   CALL dealloc_dnSca(lambda12R)
   CALL dealloc_dnSca(lambda13R)


   IF (.NOT. Para_Phenol%PubliUnit) PotVal = PotVal*(ONE/auTOev) ! to convert the eV into Hartree


   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL Write_dnMatPot(PotVal,6)
   !write(out_unitp,*)
   !flush(out_unitp)

  END SUBROUTINE eval_PhenolPot

!> @brief Subroutine wich calculates the Phenol potential with derivatives up to the 2d order if required.
!> @brief This old subroutine is the same as the "eval_PhenolPot" subroutine. However, it does not use the "dnSca module".
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                table of two values for which the potential is calculated (R,theta)
!! @param Para_Phenol        TYPE(Param_Phenol):  derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_PhenolPot_old(PotVal,Q,Para_Phenol,nderiv)
    USE mod_dnMatPot

    TYPE (Param_Phenol), intent(in)     :: Para_Phenol
    real (kind=Rkind),       intent(in)     :: Q(2)
    TYPE(dnMatPot),      intent(inout)  :: PotVal
    integer,             intent(in)     :: nderiv

    TYPE(dnMatPot)  :: v10R,v11R
    TYPE(dnMatPot)  :: v30R,v31R
    TYPE(dnMatPot)  :: lambda12R,lambda13R

    TYPE(dnMatPot)  :: v20R,v20pR,v20mR,vsq20mR,v201R,v202R
    TYPE(dnMatPot)  :: v21R,v21pR,v21mR,vsq21mR,v211R,v212R
    TYPE(dnMatPot)  :: v22R,v22pR,v22mR,vsq22mR,v221R,v222R

    real (kind=Rkind)   :: R,th
    real (kind=Rkind)   :: d0vth,d1vth,d2vth
    real (kind=Rkind)   :: d0v2th,d1v2th,d2v2th

    IF ( Check_NotAlloc_dnMatPot(PotVal,nderiv) ) THEN
      CALL alloc_dnMatPot(PotVal,nsurf=1,ndim=2,nderiv=nderiv)
    END IF

   !write(out_unitp,*) 'phenol pot in:'
   PotVal = ZERO
   !CALL Write_dnMatPot(PotVal)

   R  = Q(1)
   th = Q(2)

   ! for V(1,1): first diabatic state
   !write(out_unitp,*) 'morse:'
   CALL Eval_MorsePot(v10R,R,Para_Phenol%v10,nderiv)
   !CALL Write_dnMatPot(v10R,6)

   !write(out_unitp,*) 'sigmoid:'
   CALL Eval_SigmoidPot(v11R,R,Para_Phenol%v11,nderiv)
   !CALL Write_dnMatPot(v11R,6)

   d0vth = ONE-cos(th+th)
   d1vth = TWO*sin(th+th)
   d2vth = FOUR*cos(th+th)
   !write(out_unitp,*) 'f(th):',d0vth,d1vth,d2vth


   PotVal%d0(1,1) = v10R%d0(1,1) + v11R%d0(1,1) * d0vth

   IF (nderiv > 0) THEN
     ! derivative with repspect to R:
     PotVal%d1(1,1,1) = v10R%d1(1,1,1) + v11R%d1(1,1,1) * d0vth
     ! derivative with repspect to th:
     PotVal%d1(1,1,2) =                  v11R%d0(1,1)   * d1vth
   END IF

   IF (nderiv > 1) THEN
     ! derivative with repspect to R,R:
     PotVal%d2(1,1,1,1) = v10R%d2(1,1,1,1) + v11R%d2(1,1,1,1) * d0vth
     ! derivative with repspect to th,th:
     PotVal%d2(1,1,2,2) =                    v11R%d0(1,1)     * d2vth
     ! derivative with repspect to R,th or th,R:
     PotVal%d2(1,1,1,2) =                    v11R%d1(1,1,1)   * d1vth
     PotVal%d2(1,1,2,1) = PotVal%d2(1,1,1,2)
   END IF

   CALL dealloc_dnMatPot(v10R)
   CALL dealloc_dnMatPot(v11R)

   ! for V(2,2): 2d diabatic state
   !write(out_unitp,*) 'morse:'
   CALL Eval_MorsePot(v201R,R,Para_Phenol%v201,nderiv)
   v201R%d0 = v201R%d0 + Para_Phenol%B204
   !CALL Write_dnMatPot(v201R,6)

   CALL alloc_dnMatPot(v202R,nsurf=1,ndim=1,nderiv=nderiv,name_var='v202R',name_sub='eval_PhenolPot')
   ! write(nio,*) ' v202=B205*exp(-B206*(r-B207)) + B208'
   v202R%d0 =  Para_Phenol%B205*exp(-Para_Phenol%B206*(R-Para_Phenol%B207))
   IF (nderiv > 0) v202R%d1 = -Para_Phenol%B206*v202R%d0(1,1)
   IF (nderiv > 1) v202R%d2 = -Para_Phenol%B206*v202R%d1(1,1,1)
   v202R%d0 =  v202R%d0 + Para_Phenol%B208

   v20pR   = v201R + v202R
   v20mR   = v201R - v202R
   v20mR   = v20mR**TWO + Para_Phenol%X20

   vsq20mR = v20mR**HALF
   v20R    = HALF*(v20pR - vsq20mR)

   !write(out_unitp,*) 'sigmoid:'
   CALL Eval_SigmoidPot(v211R,R,Para_Phenol%v211,nderiv)
   CALL Eval_SigmoidPot(v212R,R,Para_Phenol%v212,nderiv)
   v212R%d0 = v212R%d0 + Para_Phenol%B217
   v21pR   = v211R + v212R
   v21mR   = v211R - v212R
   v21mR   = v21mR**TWO + Para_Phenol%X21
   vsq21mR = v21mR**HALF
   v21R    = HALF*(v21pR + vsq21mR)

   CALL Eval_SigmoidPot(v221R,R,Para_Phenol%v221,nderiv)
   CALL Eval_SigmoidPot(v222R,R,Para_Phenol%v222,nderiv)
   v22pR   = v221R + v222R
   v22mR   = v221R - v222R
   v22mR   = v22mR**TWO + Para_Phenol%X22
   vsq22mR = v22mR**HALF
   v22R    = HALF*(v22pR - vsq22mR)

   d0vth = ONE-cos(th+th)
   d1vth = TWO*sin(th+th)
   d2vth = FOUR*cos(th+th)

   d0v2th = d0vth**2
   d1v2th = TWO*d1vth*d0vth
   d2v2th = TWO*(d2vth*d0vth+d1vth**2)

   PotVal%d0(2,2)    = v20R%d0(1,1)    + v21R%d0(1,1)   * d0vth + v22R%d0(1,1)   * d0v2th

   IF (nderiv > 0) THEN
     ! derivative with repspect to R:
     PotVal%d1(2,2,1) = v20R%d1(1,1,1) + v21R%d1(1,1,1) * d0vth + v22R%d1(1,1,1) * d0v2th
     ! derivative with repspect to th:
     PotVal%d1(2,2,2) =                  v21R%d0(1,1)   * d1vth + v22R%d0(1,1)   * d1v2th
   END IF

   IF (nderiv > 1) THEN
     ! derivative with repspect to R,R:
     PotVal%d2(2,2,1,1) = v20R%d2(1,1,1,1) + v21R%d2(1,1,1,1) * d0vth + v22R%d2(1,1,1,1) * d0v2th
     ! derivative with repspect to th,th:
     PotVal%d2(2,2,2,2) =                    v21R%d0(1,1)     * d2vth + v22R%d0(1,1)     * d2v2th
     ! derivative with repspect to R,th or th,R:
     PotVal%d2(2,2,1,2) =                    v21R%d1(1,1,1)   * d1vth + v22R%d1(1,1,1)   * d1v2th
     PotVal%d2(2,2,2,1) = PotVal%d2(2,2,1,2)
   END IF

   CALL dealloc_dnMatPot(v20R)
   CALL dealloc_dnMatPot(v20pR)
   CALL dealloc_dnMatPot(v20mR)
   CALL dealloc_dnMatPot(vsq20mR)
   CALL dealloc_dnMatPot(v201R)
   CALL dealloc_dnMatPot(v202R)

   CALL dealloc_dnMatPot(v21R)
   CALL dealloc_dnMatPot(v21pR)
   CALL dealloc_dnMatPot(v21mR)
   CALL dealloc_dnMatPot(vsq21mR)
   CALL dealloc_dnMatPot(v211R)
   CALL dealloc_dnMatPot(v212R)

   CALL dealloc_dnMatPot(v22R)
   CALL dealloc_dnMatPot(v22pR)
   CALL dealloc_dnMatPot(v22mR)
   CALL dealloc_dnMatPot(vsq22mR)
   CALL dealloc_dnMatPot(v221R)
   CALL dealloc_dnMatPot(v222R)


   ! for V(3,3): 3d diabatic state
   !write(out_unitp,*) 'morse:'
   CALL Eval_MorsePot(v30R,R,Para_Phenol%v30,nderiv)
   !CALL Write_dnMatPot(v30R,6)
   !write(out_unitp,*) 'sigmoid:'
   CALL Eval_SigmoidPot(v31R,R,Para_Phenol%v31,nderiv)
   !CALL Write_dnMatPot(v31R,6)

   PotVal%d0(3,3) = v30R%d0(1,1) + v31R%d0(1,1) * d0vth + Para_Phenol%a30

   IF (nderiv > 0) THEN
     ! derivative with repspect to R:
     PotVal%d1(3,3,1) = v30R%d1(1,1,1) + v31R%d1(1,1,1) * d0vth
     ! derivative with repspect to th:
     PotVal%d1(3,3,2) =                  v31R%d0(1,1)   * d1vth
   END IF

   IF (nderiv > 1) THEN
     ! derivative with repspect to R,R:
     PotVal%d2(3,3,1,1) = v30R%d2(1,1,1,1) + v31R%d2(1,1,1,1) * d0vth
     ! derivative with repspect to th,th:
     PotVal%d2(3,3,2,2) =                    v31R%d0(1,1)     * d2vth
     ! derivative with repspect to R,th or th,R:
     PotVal%d2(3,3,1,2) =                    v31R%d1(1,1,1)   * d1vth
     PotVal%d2(3,3,2,1) = PotVal%d2(3,3,1,2)
   END IF

   CALL dealloc_dnMatPot(v30R)
   CALL dealloc_dnMatPot(v31R)


   CALL Eval_SigmoidPot(lambda12R,R,Para_Phenol%lambda12,nderiv)
   CALL Eval_SigmoidPot(lambda13R,R,Para_Phenol%lambda13,nderiv)
   !CALL Write_dnMatPot(lambda12R,6)
   d0vth =  sin(th)
   d1vth =  cos(th)
   d2vth = -sin(th)

   PotVal%d0(1,2) = lambda12R%d0(1,1) * d0vth
   PotVal%d0(2,1) = PotVal%d0(1,2)

   PotVal%d0(1,3) = lambda13R%d0(1,1) * d0vth
   PotVal%d0(3,1) = PotVal%d0(1,3)

   IF (nderiv > 0) THEN
     ! derivative with repspect to R:
     PotVal%d1(1,2,1) = lambda12R%d1(1,1,1) * d0vth
     PotVal%d1(2,1,1) = PotVal%d1(1,2,1)

     PotVal%d1(1,3,1) = lambda13R%d1(1,1,1) * d0vth
     PotVal%d1(3,1,1) = PotVal%d1(1,3,1)

     ! derivative with repspect to th:
     PotVal%d1(1,2,2) = lambda12R%d0(1,1)   * d1vth
     PotVal%d1(2,1,2) = PotVal%d1(1,2,2)

     PotVal%d1(1,3,2) = lambda13R%d0(1,1)   * d1vth
     PotVal%d1(3,1,2) = PotVal%d1(1,3,2)
   END IF

   IF (nderiv > 1) THEN
     ! derivative with repspect to R,R:
     PotVal%d2(1,2,1,1) = lambda12R%d2(1,1,1,1) * d0vth
     PotVal%d2(2,1,1,1) = PotVal%d2(1,2,1,1)

     PotVal%d2(1,3,1,1) = lambda13R%d2(1,1,1,1) * d0vth
     PotVal%d2(3,1,1,1) = PotVal%d2(1,3,1,1)

     ! derivative with repspect to th,th:
     PotVal%d2(1,2,2,2) = lambda12R%d0(1,1)     * d2vth
     PotVal%d2(2,1,2,2) = PotVal%d2(1,2,2,2)

     PotVal%d2(1,3,2,2) = lambda13R%d0(1,1)     * d2vth
     PotVal%d2(3,1,2,2) = PotVal%d2(1,3,2,2)

     ! derivative with repspect to R,th or th,R:
     PotVal%d2(1,2,1,2) = lambda12R%d1(1,1,1)   * d1vth
     PotVal%d2(1,2,2,1) = PotVal%d2(1,2,1,2)
     PotVal%d2(2,1,1,2) = PotVal%d2(1,2,1,2)
     PotVal%d2(2,1,2,1) = PotVal%d2(1,2,1,2)

     PotVal%d2(1,3,1,2) = lambda13R%d1(1,1,1)   * d1vth
     PotVal%d2(1,3,2,1) = PotVal%d2(1,3,1,2)
     PotVal%d2(3,1,1,2) = PotVal%d2(1,3,1,2)
     PotVal%d2(3,1,2,1) = PotVal%d2(1,3,1,2)

   END IF

   CALL dealloc_dnMatPot(lambda12R)
   CALL dealloc_dnMatPot(lambda13R)


   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL Write_dnMatPot(PotVal,6)
   !write(out_unitp,*)
   !flush(out_unitp)

  END SUBROUTINE eval_PhenolPot_old

END MODULE mod_PhenolPot
