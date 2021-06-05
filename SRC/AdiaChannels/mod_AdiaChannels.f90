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

MODULE mod_AdiaChannels
  USE QML_NumParameters_m

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Make_Hinact(Qact,QModel)
  USE QML_NumParameters_m
  USE QML_UtilLib_m
  USE QML_diago_m
  USE QML_dnS_m
  USE QML_dnMat_m
  USE mod_Model
  USE QML_Basis_m
  IMPLICIT NONE

  real (kind=Rkind), allocatable, intent(in)    :: Qact(:)
  TYPE (Model_t),                 intent(inout) :: QModel

  real(kind=Rkind), parameter :: auTOcm_inv = 219474.631443_Rkind

  integer                        :: i,iq,ib,jb,nb,nq
  !real (kind=Rkind), allocatable :: V(:),H(:,:),HB(:),Vec(:,:)
  real (kind=Rkind), allocatable :: Ene(:)

  TYPE (dnMat_t)                 :: PotVal
  real (kind=Rkind), allocatable :: Q(:),d0GGdef(:,:)


  integer :: nderiv
  TYPE (dnMat_t)              :: dnH ! derivative of the Hamiltonian
  TYPE (dnMat_t)              :: dnVec,dnVecProj ! Eigenvectors
  TYPE (dnMat_t)              :: dnHdiag ! derivative of the Hamiltonian
  integer                     :: ndim_act

  TYPE (dnS_t), allocatable :: dnV(:),dnHB(:)
  TYPE (dnS_t) :: dnVfull,dnHij


  ndim_act = size(QModel%QM%list_act)
  nderiv = 1

  nb = QModel%Basis%nb
  nq = QModel%Basis%nq


  CALL Eval_dnHVib_ana(QModel,Qact,dnH,nderiv)

  allocate(Ene(nb))


  CALL QML_alloc_dnMat(dnVec,    nsurf=nb,ndim=ndim_act,nderiv=nderiv,name_var='dnVec')
  CALL QML_alloc_dnMat(dnHdiag,  nsurf=nb,ndim=ndim_act,nderiv=nderiv,name_var='dnHdiag')
  CALL QML_alloc_dnMat(dnVecProj,nsurf=nb,ndim=ndim_act,nderiv=nderiv,name_var='dnVecProj')

  IF (.NOT. allocated(QModel%QM%Vec0)) allocate(QModel%QM%Vec0)
  IF (QML_Check_NotAlloc_dnMat(QModel%QM%Vec0,nderiv=0)) THEN
!$OMP CRITICAL (CRIT_Make_Hinact)
    CALL QML_alloc_dnMat(QModel%QM%Vec0,nsurf=nb,ndim=ndim_act,nderiv=0)
    CALL diagonalization(dnH%d0,Ene,QModel%QM%Vec0%d0,nb,sort=1,phase=.TRUE.)
    write(out_unitp,*) 'init Vec0 done'
!$OMP END CRITICAL (CRIT_Make_Hinact)
  END IF

  CALL QML_DIAG_dnMat(dnH,dnHdiag,dnVec,dnVecProj,dnVec0=QModel%QM%Vec0)
  !CALL Write_RMat(dnVec%d0(:,1:QModel%QM%nb_Channels),6,6,name_info='Vec')
  write(out_unitp,*) 'NAC'
  !CALL QML_Write_dnMat(dnVecProj,nio=out_unitp)
  CALL Write_RMat(dnVecProj%d1(1:QModel%QM%nb_Channels,1:QModel%QM%nb_Channels,1),out_unitp,6,name_info='NAC')

  DO ib=1,nb
    Ene(ib) = dnHdiag%d0(ib,ib)
  END DO
  write(out_unitp,*) 'Ene'
  CALL Write_RVec(Ene(1:QModel%QM%nb_Channels)*auTOcm_inv,out_unitp,6,name_info='Ene')

  END SUBROUTINE Make_Hinact

END MODULE mod_AdiaChannels
