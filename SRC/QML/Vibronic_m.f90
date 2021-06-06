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
!!
!> @author David Lauvergnat
!! @date 02/12/2017
!!
MODULE QML_Vibronic_m
  USE QML_NumParameters_m
  USE QML_dnMat_m
  USE QML_Empty_m
  IMPLICIT NONE

!> @brief Derived type in which the parameters of the Phenol potential are set-up.
!> @brief Reference: Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218
!> @brief the parameter names are from the previous reference and are taken from Eqs 4-23 and tables I-IV.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
  TYPE, EXTENDS (QML_Empty_t) ::  VibronicPot_t
     PRIVATE
     ! Vij, the diabatic potentials (i=j) or couplings (i/=j), are
     !   expanded as a taylor expension at second order
     !   around a reference geometry QijO (each component is QijkO):
     !Vij(Q) = Vij0 + Sum_k g_ijk * (Qk-QijkO) + 1/2*Sum_kl h_ijkj * (Qk-QijkO)(Ql-QijlO)

     ! V, where:
     ! V%dO(i,j) is Vij0
     ! V%d1(i,j,k) is g_ijk
     ! V%d2(i,j,k,l) is h_ijkl
     TYPE (dnMat_t)   :: V

     !Q0(:,i,j) is the reference geometry for the diabatic (j=i) or compling (i/=j) terms
     real (kind=Rkind), allocatable     :: Q0(:,:,:)
  CONTAINS
    PROCEDURE :: Eval_QModel_Pot => eval_VibronicPot
    PROCEDURE :: Write_QModel    => Write_VibronicModel
    PROCEDURE :: Write0_QModel   => Write0_VibronicModel
  END TYPE VibronicPot_t


CONTAINS
!> @brief Subroutine which makes the initialization of the Morse parameters.
!!      Those parameters cannot be modified (PRIVATE).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_Phenol        TYPE(Param_Phenol):   derived type in which the parameters are set-up.
!! @param PubliUnit          logical (optional):   when PubliUnit=.TRUE., the units (Angstrom and eV) are used. Default (atomic unit).
  SUBROUTINE Init_Vibronic(VibronicPot,ndim,nsurf,nio,PubliUnit)
    TYPE (VibronicPot_t),      intent(inout)   :: VibronicPot
     logical, optional,         intent(in)      :: PubliUnit

     integer :: i,j,IOerr,order,nbcol
     real (kind=Rkind)     :: V0,V0coupling

     namelist /VibronicSurf/ order,nbcol,V0

     IF (present(PubliUnit)) VibronicPot%PubliUnit = PubliUnit

     CALL QML_alloc_dnMat(VibronicPot%V,nsurf,ndim,nderiv, &
             name_var='VibronicPot%V',name_sub='Init_Vibronic',IOerr)
     IF (IOerr /= 0 ) THEN
       write(out_unitp,*) ' ERROR in Init_Vibronic'
       write(out_unitp,*) '  Problem with allocate of VibronicPot%V'
       STOP
     END IF
     VibronicPot%V = ZERO


     allocate(VibronicPot%Q0(ndim,nsurf,nsurf),stat=IOerr)
     IF (IOerr /= 0 .OR. nsurf < 1 .OR. ndim < 1) THEN
       write(out_unitp,*) ' ERROR in Init_Vibronic'
       write(out_unitp,*) '  Problem with allocate of VibronicPot%Q0'
       write(out_unitp,*) '  nsurf > 0?',nsurf
       write(out_unitp,*) '  ndim > 0?',ndim
       STOP
     END IF
     VibronicPot%Q0 = ZERO


     DO i=1,nsurf
       V0         = ZERO
       V0coupling = ZERO
       order      = -1
       nbcol      = 5
       read(VibronicSurf,nio)
       VibronicPot%VOrder(i,i) = order

       IF (order < 0) CYCLE

       CALL Read_RVec(VibronicPot%Q0(:,i,i),nio,nbcol,IOerr)

       VibronicPot%V%d0(i,i) = V0coupling
       IF (order > 0) THEN
         CALL Read_RVec(VibronicPot%V%d1(:,i,i),nio,nbcol,IOerr)
       END IF
       IF (order > 1) THEN
         CALL Read_RMat(VibronicPot%V%d2(:,:,i,i),nio,nbcol,IOerr)
       END IF

     END DO

     DO i=1,nsurf
     DO j=i+1,nsurf
       V0         = ZERO
       V0coupling = ZERO
       order      = -1
       nbcol      = 5
       read(VibronicSurf,nio)
       VibronicPot%VOrder(i,j) = order
       IF (order < 0) CYCLE

       CALL Read_RVec(VibronicPot%Q0(:,i,j),nio,nbcol,IOerr)

       VibronicPot%V%d0(i,i) = V0
       IF (order > 0) THEN
         CALL Read_RVec(VibronicPot%V%d1(:,i,j),nio,nbcol,IOerr)
       END IF
       IF (order > 1) THEN
         CALL Read_RMat(VibronicPot%V%d2(:,:,i,j),nio,nbcol,IOerr)
       END IF

     END DO


  END SUBROUTINE Init_Vibronic

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
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                table of two values for which the potential is calculated (R,theta)
!! @param Para_Phenol        TYPE(Param_Phenol):  derived type with the Morse parameters.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_Vibronic(PotVal,Q,Para_Phenol,nderiv)
    USE QML_dnMat_m
    USE QML_dnS_m

    TYPE (Param_Phenol), intent(in)     :: Para_Phenol
    real (kind=Rkind),   intent(in)     :: Q(2)
    TYPE (dnMat_t),      intent(inout)  :: PotVal
    integer,             intent(in)     :: nderiv

    TYPE (dnS_t)  :: dnR,dnth
    TYPE (dnS_t)  :: v10R,v11R,v11th
    TYPE (dnS_t)  :: v30R,v31R,v31th
    TYPE (dnS_t)  :: lambda12R,lambda13R

    TYPE (dnS_t)  :: v20pR,v20mR,v201R,v202R
    TYPE (dnS_t)  :: v21pR,v21mR,v211R,v212R
    TYPE (dnS_t)  :: v22pR,v22mR,v221R,v222R
    TYPE (dnS_t)  :: v20R,v21R,v22R,v21th,v22th

    real (kind=Rkind) :: a0      = 0.52917720835354106_Rkind
    real (kind=Rkind) :: auTOeV  = 27.211384_Rkind

    IF ( Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL alloc_dnMat(PotVal,nsurf=1,ndim=2,nderiv=nderiv)
    END IF

   !write(out_unitp,*) 'phenol pot in:'
   PotVal = ZERO
   !CALL QML_Write_dnMat(PotVal)

   dnR     = init_dnS(Q(1),ndim=2,nderiv=nderiv,iQ=1) ! to set up the derivatives
   dnth    = init_dnS(Q(2),ndim=2,nderiv=nderiv,iQ=2) ! to set up the derivatives

   IF (.NOT. Para_Phenol%PubliUnit) dnR = a0*dnR ! to convert the bhor into Angstrom

   ! for V(1,1): first diabatic state
   !write(out_unitp,*) 'morse:'
   v10R = QML_dnMorse(dnR,Para_Phenol%v10)
   !CALL QML_Write_dnS(v10R,6)
   !write(out_unitp,*) 'sigmoid:'
   v11R = dnSigmoid(dnR,Para_Phenol%v11)
   !CALL QML_Write_dnS(v11R,6)

   !write(out_unitp,*) 'f(th):'
   v11th = ONE-cos(dnth+dnth)
   !CALL QML_Write_dnS(v11th,6)

   CALL sub_dnS_TO_dnMat(v10R+v11R*v11th,PotVal,i=1,j=1)

   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL QML_Write_dnMat(PotVal,6)

   CALL dealloc_dnS(v10R)
   CALL dealloc_dnS(v11R)
   CALL dealloc_dnS(v11th)

   ! for V(2,2): 2d diabatic state
   v201R = QML_dnMorse(dnR,Para_Phenol%v201) + Para_Phenol%B204
   v202R = Para_Phenol%B205*exp(-Para_Phenol%B206*(dnR-Para_Phenol%B207)) + Para_Phenol%B208
   v20pR = v201R + v202R
   v20mR = v201R - v202R
   v20R  = HALF*(v20pR - (v20mR**TWO + Para_Phenol%X20)**HALF)


   !write(out_unitp,*) 'sigmoid:'
   v211R = dnSigmoid(dnR,Para_Phenol%v211)
   v212R = dnSigmoid(dnR,Para_Phenol%v212) + Para_Phenol%B217
   v21pR = v211R + v212R
   v21mR = v211R - v212R
   v21R  = HALF * (v21pR + (v21mR**TWO + Para_Phenol%X21)**HALF)

   v221R = dnSigmoid(dnR,Para_Phenol%v221)
   v222R = dnSigmoid(dnR,Para_Phenol%v222)
   v22pR = v221R + v222R
   v22mR = v221R - v222R
   v22R  = HALF * (v22pR - sqrt(v22mR**TWO + Para_Phenol%X22) )

   v21th = ONE-cos(dnth+dnth)
   v22th = v21th*v21th

   CALL sub_dnS_TO_dnMat(v20R+v21R*v21th+v22R*v22th,PotVal,i=2,j=2)


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


   ! for V(3,3): 3d diabatic state
   !write(out_unitp,*) 'morse:'
   v30R = QML_dnMorse(dnR,Para_Phenol%v30) + Para_Phenol%a30
   !CALL QML_Write_dnMat(v30R,6)
   !write(out_unitp,*) 'sigmoid:'
   v31R = dnSigmoid(dnR,Para_Phenol%v31)
   !CALL QML_Write_dnMat(v31R,6)

   !write(out_unitp,*) 'f(th):'
   v31th = ONE-cos(dnth+dnth)
   !CALL QML_Write_dnS(v11th,6)

   CALL sub_dnS_TO_dnMat(v30R+v31R*v31th,PotVal,i=3,j=3)

   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL QML_Write_dnMat(PotVal,6)

   CALL dealloc_dnS(v30R)
   CALL dealloc_dnS(v31R)
   CALL dealloc_dnS(v31th)


   lambda12R = dnSigmoid(dnR,Para_Phenol%lambda12) * sin(dnth)
   !CALL QML_Write_dnS(lambda12R,6)

   CALL sub_dnS_TO_dnMat(lambda12R,PotVal,i=1,j=2)
   CALL sub_dnS_TO_dnMat(lambda12R,PotVal,i=2,j=1)


   lambda13R = dnSigmoid(dnR,Para_Phenol%lambda13) * sin(dnth)

   CALL sub_dnS_TO_dnMat(lambda13R,PotVal,i=1,j=3)
   CALL sub_dnS_TO_dnMat(lambda13R,PotVal,i=3,j=1)


   CALL dealloc_dnS(lambda12R)
   CALL dealloc_dnS(lambda13R)


   IF (.NOT. Para_Phenol%PubliUnit) PotVal = PotVal*(ONE/auTOev) ! to convert the eV into Hartree


   !write(out_unitp,*) 'phenol pot diabatic:',nderiv
   !CALL QML_Write_dnMat(PotVal,6)
   !write(out_unitp,*)
   !flush(out_unitp)

  END SUBROUTINE eval_Vibronic


END MODULE QML_Vibronic_m
