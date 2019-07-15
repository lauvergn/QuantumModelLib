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

!> @brief Module which makes the initialization, calculation of the TwoD potential (value, gradient and hessian).
!> @brief Reference:A. Ferretti, G. Granucci, A. Lami, M. Persico, G. Villani, J. Chem. Phys. 104, 5517 (1996); https://doi.org/10.1063/1.471791
!> @brief The potential is a 2D (X,Y) with 2 electronic diabatic surfaces (S1,S2)
!!
!> @author David Lauvergnat
!! @date 12/07/2019
!!
MODULE mod_TwoDPot
  USE mod_NumParameters

  IMPLICIT NONE

!> @brief Derived type in which the parameters of the 2D model are set-up.
!> @brief Reference: A. Ferretti, G. Granucci, A. Lami, M. Persico, G. Villani, J. Chem. Phys. 104, 5517 (1996); https://doi.org/10.1063/1.471791
!> @brief the parameters have been modified.
!!
!> @author David Lauvergnat
!! @date 12/07/2019
  TYPE Param_TwoD
     PRIVATE
     logical              :: PubliUnit = .FALSE. ! when PubliUnit=.TRUE., the units (Angstrom and Ev) are used. Default (atomic unit)

    real(kind=Rkind)     :: KX    = 0.02_Rkind
    real(kind=Rkind)     :: KY    = 0.1_Rkind
    real(kind=Rkind)     :: DELTA = 0.01_Rkind
    real(kind=Rkind)     :: X1    = 6._Rkind
    real(kind=Rkind)     :: X2    = 2._Rkind
    real(kind=Rkind)     :: X3    = 31._Rkind/8.0_Rkind
    real(kind=Rkind)     :: GAMMA = 0.01_Rkind
    real(kind=Rkind)     :: ALPHA = 3._Rkind
    real(kind=Rkind)     :: BETA  = 1.5_Rkind

     real (kind=Rkind), PUBLIC :: muX  = 20000._Rkind
     real (kind=Rkind), PUBLIC :: muY  = 6667._Rkind


  END TYPE Param_TwoD


CONTAINS
!> @brief Subroutine which makes the initialization of the TwoD parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 12/07/2019
!!
!! @param Para_TwoD          TYPE(Param_TwoD):     derived type in which the parameters are set-up.
!! @param PubliUnit          logical (optional):   when PubliUnit=.TRUE., the units (Angstrom and eV) are used. Default (atomic unit).
  SUBROUTINE Init_TwoDPot(Para_TwoD,PubliUnit)
    TYPE (Param_TwoD),        intent(inout)   :: Para_TwoD
     logical, optional,       intent(in)      :: PubliUnit

     IF (present(PubliUnit)) Para_TwoD%PubliUnit = PubliUnit

     CALL Write_TwoDPot(Para_TwoD,6)

  END SUBROUTINE Init_TwoDPot

!> @brief Subroutine wich prints the TwoD potential parameters.
!!
!> @author David Lauvergnat
!! @date 12/07/2019
!!
!! @param Para_TwoD        TYPE(Param_TwoD):     derived type with the TwoD potential parameters.
!! @param nio              integer:              file unit to print the parameters.
  SUBROUTINE Write_TwoDPot(Para_TwoD,nio)
    TYPE (Param_TwoD), intent(in) :: Para_TwoD
    integer, intent(in) :: nio

    write(nio,*) 'TwoD parameters'
    write(nio,*) '-----------------------------------------'
    write(nio,*) '---- WARNNING ---------------------------'
    write(nio,*) 'The parameters are different from the published ones: '
    write(nio,*) ' A. Ferretti, G. Granucci, A. Lami, M. Persico, G. Villani, ...'
    write(nio,*) '  .... J. Chem. Phys. 104, 5517 (1996); https://doi.org/10.1063/1.471791'
    write(nio,*) 'with the X=Q(1), Y=Q(2) in bohr.'
    write(nio,*) '     and the energy in Hartree.'
    write(nio,*)
    write(nio,*) 'Diabatic Potential Values (in Hartree) at: R=3.875 bohr and theta=0.5 bohr'
    write(nio,*) '1        0.05765625  0.00343645'
    write(nio,*) '2        0.00343645  0.05765625'


    write(nio,*) 'PubliUnit: ',Para_TwoD%PubliUnit

    write(nio,*)
    write(nio,*) 'Current parameters:'
    write(nio,*) '-----------------------------------------'
    write(nio,*) 'KX    =',Para_TwoD%KX
    write(nio,*) 'KY    =',Para_TwoD%KY
    write(nio,*) 'DELTA =',Para_TwoD%DELTA
    write(nio,*) 'X1    =',Para_TwoD%X1
    write(nio,*) 'X2    =',Para_TwoD%X2
    write(nio,*) 'X3    =',Para_TwoD%X3
    write(nio,*) 'GAMMA =',Para_TwoD%GAMMA
    write(nio,*) 'ALPHA =',Para_TwoD%ALPHA
    write(nio,*) 'BETA  =',Para_TwoD%BETA
    write(nio,*) '-----------------------------------------'



    write(nio,*) 'end TwoD parameters'

  END SUBROUTINE Write_TwoDPot

!> @brief Subroutine wich calculates the TwoD potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 12/07/2019
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                table of two values for which the potential is calculated (R,theta)
!! @param Para_TwoD          TYPE(Param_TwoD):    derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TwoDPot(PotVal,Q,Para_TwoD,nderiv)
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_TwoD),   intent(in)     :: Para_TwoD
    real (kind=Rkind),   intent(in)     :: Q(2)
    TYPE(dnMatPot),      intent(inout)  :: PotVal
    integer,             intent(in)     :: nderiv

    TYPE(dnSca)  :: dnX,dnY
    TYPE(dnSca)  :: vTemp



    IF ( Check_NotAlloc_dnMatPot(PotVal,nderiv) ) THEN
      CALL alloc_dnMatPot(PotVal,nsurf=2,ndim=2,nderiv=nderiv)
    END IF

   !write(out_unitp,*) 'TwoD pot in:'
   PotVal = ZERO
   !CALL Write_dnMatPot(PotVal)

   dnX     = init_dnSca(Q(1),ndim=2,nderiv=nderiv,iQ=1) ! to set up the derivatives
   dnY     = init_dnSca(Q(2),ndim=2,nderiv=nderiv,iQ=2) ! to set up the derivatives


   !Hel(1,1,x,y)=0.5d0*KX*(R(1,x)-X1)**2 + 0.5d0*KY*(R(2,y))**2
   vTemp = HALF * Para_TwoD%KX*(dnX-Para_TwoD%X1)**2 + HALF * Para_TwoD%KY*(dnY)**2
   CALL sub_dnSca_TO_dnMatPot(vTemp,PotVal,i=1,j=1)

   !Hel(2,2,x,y)=0.5d0*KX*(R(1,x)-X2)**2 + 0.5d0*KY*(R(2,y))**2 +DELTA
   vTemp = HALF * Para_TwoD%KX*(dnX-Para_TwoD%X2)**2 + HALF * Para_TwoD%KY*(dnY)**2 + Para_TwoD%DELTA
   CALL sub_dnSca_TO_dnMatPot(vTemp,PotVal,i=2,j=2)

   !Hel(1,2,x,y)=GAMMA*R(2,y)*dexp(-ALPHA*(R(1,x)-X3)**2)*dexp(-BETA*(R(2,y))**2)
   vTemp = Para_TwoD%GAMMA * dnY * exp(-Para_TwoD%ALPHA*(dnX-Para_TwoD%X3)**2) * exp(-Para_TwoD%BETA*(dnY)**2)

   !write(out_unitp,*) 'vTemp for V(1,2):'
   !CALL Write_dnSca(vTemp,6)


   CALL sub_dnSca_TO_dnMatPot(vTemp,PotVal,i=1,j=2)
   CALL sub_dnSca_TO_dnMatPot(vTemp,PotVal,i=2,j=1)


   CALL dealloc_dnSca(dnX)
   CALL dealloc_dnSca(dnY)
   CALL dealloc_dnSca(vTemp)


   !write(out_unitp,*) 'TwoD pot diabatic:',nderiv
   !CALL Write_dnMatPot(PotVal,6)
   !write(out_unitp,*)
   !flush(out_unitp)

  END SUBROUTINE eval_TwoDPot


END MODULE mod_TwoDPot
