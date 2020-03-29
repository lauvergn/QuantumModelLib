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
  TYPE TemplatePot_t
     PRIVATE

     ! V(1,1) term
     TYPE (MorsePot_t)   :: morseXpY,morseZ
     real(kind=Rkind)     :: kXmY=0.1_Rkind

  END TYPE TemplatePot_t


CONTAINS
!> @brief Subroutine which makes the initialization of the Morse parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Pot        TYPE(TemplatePot_t):   derived type in which the parameters are set-up.
  SUBROUTINE Init_TemplatePot(Para_Pot)
    TYPE (TemplatePot_t),      intent(inout)   :: Para_Pot


    ! V(1,1) term, all parameters in atomic unit (Hartree, bohr)

    CALL Init_MorsePot(Para_Pot%morseXpY,D=0.1_Rkind,a=1._Rkind,req=2._Rkind)

    CALL Init_MorsePot(Para_Pot%morseZ,D=0.08_Rkind,a=1._Rkind,req=2._Rkind)

  END SUBROUTINE Init_TemplatePot

!> @brief Subroutine wich prints the potential parameters.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer:                file unit to print the parameters.
  SUBROUTINE Write0_TemplatePot(nio)
    integer, intent(in) :: nio

    write(nio,*) 'TemplatePot default parameters'
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
    write(nio,*) 'end TemplatePot default parameters'

  END SUBROUTINE Write0_TemplatePot

!> @brief Subroutine wich prints the potential parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Pot           TYPE(TemplatePot_t):   derived type with the potential parameters
!! @param nio                integer:                file unit to print the parameters.
  SUBROUTINE Write_TemplatePot(Para_Pot,nio)
    TYPE (TemplatePot_t), intent(in) :: Para_Pot
    integer, intent(in) :: nio

    write(nio,*) 'TemplatePot current parameters'
    write(nio,*)
    CALL Write_MorsePot(Para_Pot%morseXpY,nio)
    write(nio,*)
    CALL Write_MorsePot(Para_Pot%morseZ,nio)
    write(nio,*)
    write(nio,*) ' kXmY:',Para_Pot%kXmY
    write(nio,*)
    write(nio,*) 'end TemplatePot current parameters'

  END SUBROUTINE Write_TemplatePot

!> @brief Subroutine wich calculates the Phenol potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                table of two values for which the potential is calculated (R,theta)
!! @param Para_Pot           TYPE(Param_Pot):     Potential parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_TemplatePot(Mat_OF_PotDia,dnQ,Para_Pot,nderiv)
    USE mod_dnS

    TYPE (TemplatePot_t), intent(in)     :: Para_Pot
    TYPE (dnS_t),           intent(in)     :: dnQ(3)
    TYPE (dnS_t),           intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,               intent(in)     :: nderiv

    TYPE (dnS_t)  :: mXpY,vXmY,mZ
    integer      :: i

   !write(out_unitp,*) 'TemplatePot in:'

   mXpY = dnMorse(dnQ(1)+dnQ(2),Para_Pot%morseXpY)

   mZ   = dnMorse(dnQ(3),Para_Pot%morseZ)

   vXmY = (Para_Pot%kXmY*HALF) * (dnQ(1)-dnQ(2))**2

   Mat_OF_PotDia(1,1) = mXpY+mZ+vXmY

   CALL QML_dealloc_dnS(mXpY)
   CALL QML_dealloc_dnS(mZ)
   CALL QML_dealloc_dnS(vXmY)


   !write(out_unitp,*) 'TemplatePot, nderiv:',nderiv
   !write(out_unitp,*) 'Q(:):',(get_d0_FROM_dnS(dnQ(i)),i=1,size(dnQ))
   !CALL QML_Write_dnS( Mat_OF_PotDia(1,1),6)
   !write(out_unitp,*)
   !flush(out_unitp)

  END SUBROUTINE eval_TemplatePot

END MODULE mod_TemplatePot
