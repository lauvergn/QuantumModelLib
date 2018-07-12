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

!> @brief Module to show how to implement a 3D-potential with ONE electronic surface.
!!
!> @author David Lauvergnat
!! @date 30/05/2018
!!
MODULE mod_TemplatePot
  USE mod_NumParameters
  USE mod_MorsePot

  IMPLICIT NONE

!> @brief Derived type in which the parameters of potential are set-up (3D-potential)
!> @brief V(1,1) = morseXpY(x+y) + k11/2(x-y)^2 + morseZ(z)

!!
!> @author David Lauvergnat
!! @date 03/08/2017
  TYPE Param_Template
     PRIVATE

     ! V(1,1) term
     TYPE (Param_Morse)   :: morseXpY,morseZ
     real(kind=Rkind)     :: k11=0.1_Rkind

  END TYPE Param_Template


CONTAINS
!> @brief Subroutine which makes the initialization of the Morse parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Pot        TYPE(Param_Template):   derived type in which the parameters are set-up.
  SUBROUTINE Init_TemplatePot(Para_Pot)
    TYPE (Param_Template),      intent(inout)   :: Para_Pot


    ! V(1,1) term, all parameters in atomic unit (Hartree, bohr)

    CALL Init_MorsePot(Para_Pot%morseXpY,D=0.1_Rkind,a=1._Rkind,req=2._Rkind)

    CALL Init_MorsePot(Para_Pot%morseZ,D=0.08_Rkind,a=1._Rkind,req=2._Rkind)

  END SUBROUTINE Init_TemplatePot

!> @brief Subroutine wich prints the potential parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Pot           TYPE(Param_Template):   derived type with the potential parameters
!! @param nio                integer:                file unit to print the parameters.
  SUBROUTINE Write_TemplatePot(Para_Pot,nio)
    TYPE (Param_Template), intent(in) :: Para_Pot
    integer, intent(in) :: nio

    write(nio,*) 'TemplatePot parameters'

    write(nio,*)
    write(nio,*) 'Current parameters:'
    write(nio,*) '-----------------------------------------'
    write(nio,*) ' V(1,1):'
    CALL Write_MorsePot(Para_Pot%morseXpY,nio)
    CALL Write_MorsePot(Para_Pot%morseZ,nio)

    write(nio,*) ' k11:',Para_Pot%k11


    write(nio,*) 'end TemplatePot parameters'

  END SUBROUTINE Write_TemplatePot

!> @brief Subroutine wich calculates the Phenol potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE(dnMatPot):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                table of two values for which the potential is calculated (R,theta)
!! @param Para_Pot           TYPE(Param_Pot):     Potential parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TemplatePot(PotVal,Q,Para_Pot,nderiv)
    USE mod_dnMatPot
    USE mod_dnSca

    TYPE (Param_Template), intent(in)     :: Para_Pot
    real (kind=Rkind),     intent(in)     :: Q(3)
    TYPE(dnMatPot),        intent(inout)  :: PotVal
    integer,               intent(in)     :: nderiv

    TYPE(dnSca)  :: dnX,dnY,dnZ
    !TYPE(dnSca)  :: dnXpY ! x+y
    !TYPE(dnSca)  :: dnXmY ! x-y

    TYPE(dnSca)  :: mXpY,vXmY,mZ



    IF ( Check_NotAlloc_dnMatPot(PotVal,nderiv) ) THEN
      CALL alloc_dnMatPot(PotVal,nsurf=1,ndim=2,nderiv=nderiv)
    END IF

   !write(out_unitp,*) 'TemplatePot in:'
   PotVal = ZERO
   !CALL Write_dnMatPot(PotVal)

   dnX     = init_dnSca(Q(1),ndim=3,nderiv=nderiv,iQ=1) ! to set up the derivatives
   dnY     = init_dnSca(Q(2),ndim=3,nderiv=nderiv,iQ=2) ! to set up the derivatives
   dnZ     = init_dnSca(Q(3),ndim=3,nderiv=nderiv,iQ=3) ! to set up the derivatives

   !dnXpY  = dnX+dnY ! x+y including the derivatives
   !dnXmY  = dnX-dnY ! x-y including the derivatives

   ! for V(1,1): first diabatic state
   !write(out_unitp,*) 'morse:'
   mXpY = dnMorse(dnX+dnY,Para_Pot%morseXpY)
   mZ   = dnMorse(dnZ,Para_Pot%morseZ)

   vXmY = (Para_Pot%k11*HALF) * (dnX-dnY)**2


   CALL sub_dnSca_TO_dnMatPot(mXpY+mZ+vXmY,PotVal,i=1,j=1)

   CALL dealloc_dnSca(dnX)
   CALL dealloc_dnSca(dnY)
   CALL dealloc_dnSca(dnZ)
   CALL dealloc_dnSca(mXpY)
   CALL dealloc_dnSca(mZ)
   CALL dealloc_dnSca(vXmY)


   write(out_unitp,*) 'TemplatePot, nderiv:',nderiv
   CALL Write_dnMatPot(PotVal,6)
   write(out_unitp,*)
   flush(out_unitp)

  END SUBROUTINE eval_TemplatePot

END MODULE mod_TemplatePot
