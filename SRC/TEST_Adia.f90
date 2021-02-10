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
PROGRAM TEST_adia
  USE mod_NumParameters
  USE mod_Lib
  USE mod_diago
  USE mod_dnS
  USE mod_dnMat
  USE mod_Model
  USE mod_AdiaChannels
  IMPLICIT NONE

  real(kind=Rkind), parameter :: auTOcm_inv = 219474.631443_Rkind
  TYPE(BasisPrim_t)  :: BasisPrim

  real (kind=Rkind)              :: A,B
  integer                        :: iq,ib,jb,nb,nq
  real (kind=Rkind), allocatable :: V(:),H(:,:),HB(:)
  real (kind=Rkind), allocatable :: Ene(:),Vec(:,:)

  TYPE (Model_t)            :: QModel
  TYPE (dnMat_t)            :: PotVal
  real (kind=Rkind), allocatable :: Q(:),d0GGdef(:,:)


  integer :: nderiv
  TYPE (dnMat_t)              :: dnH ! derivative of the Hamiltonian
  TYPE (dnMat_t)              :: dnVec,dnVecProj ! Eigenvectors
  TYPE (dnMat_t)              :: dnHdiag ! derivative of the Hamiltonian
  integer, allocatable        :: list_act(:)
  integer                     :: ndim_act

  TYPE (dnS_t), allocatable :: dnV(:),dnHB(:)
  TYPE (dnS_t) :: dnVfull,dnHij


  !initialization of the model
  CALL Init_Model(QModel,pot_name='HBond',Print_init=.FALSE.)
  allocate(Q(QModel%ndim))

  list_act = [1]
  ndim_act = size(list_act)
  nderiv = 1
  !Basis set in Q(2)
  nb = 64
  nq = 64
  A=-2.5_Rkind
  B=2.5_Rkind

  allocate(V(nq))
  allocate(HB(nq))
  allocate(H(nb,nb))
  allocate(Vec(nb,nb))
  allocate(Ene(nb))

  CALL Make_BasisPrim(BasisPrim,nb=nb,nq=nq,A=A,B=B)

  CALL QML_alloc_dnMat(dnH,      nsurf=nb,ndim=ndim_act,nderiv=nderiv,name_var='dnH')
  CALL QML_alloc_dnMat(dnVec,    nsurf=nb,ndim=ndim_act,nderiv=nderiv,name_var='dnVec')
  CALL QML_alloc_dnMat(dnHdiag,  nsurf=nb,ndim=ndim_act,nderiv=nderiv,name_var='dnHdiag')
  CALL QML_alloc_dnMat(dnVecProj,nsurf=nb,ndim=ndim_act,nderiv=nderiv,name_var='dnVecProj')

  allocate(dnV(nq))
  allocate(dnHB(nq))

  Q(1) = 6._Rkind
  !Buid H
  DO iq=1,BasisPrim%nq
    Q(2) =  BasisPrim%x(iq)
    CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
    CALL QML_sub_dnMat_TO_dnS(PotVal,dnVfull,i=1,j=1)
    CALL QML_ReduceDerivatives_dnS2_TO_dnS1(dnV(iq),dnVfull,list_act)
    !V(iq) = PotVal%d0(1,1)
  END DO
  !CALL QML_Write_dnMat(PotVal,6,info='PotVal')
  !CALL QML_Write_dnS(dnV(nq),6,info='dnV',all_type=.TRUE.)

  !CALL Write_RVec(V,6,5,name_info='V')
  !V(:) = QML_get_d0_FROM_dnS(dnV)
  !CALL Write_RVec(V,6,5,name_info='dnV%d0')
  !write(6,*) 'coucou pot: done' ; stop

  d0GGdef = QModel%QM%get_d0GGdef_QModel()
  DO ib=1,nb
    ! H B(:,ib)>
    !HB(:) = -HALF*d0GGdef(2,2)*BasisPrim%d2gb(:,ib,1,1) + V(:)*BasisPrim%d0gb(:,ib)
    !HB(:) = HB(:) * BasisPrim%w(:)
    DO iq=1,nq
      dnHB(iq) = -HALF*d0GGdef(2,2)*BasisPrim%d2gb(iq,ib,1,1) + dnV(iq)*BasisPrim%d0gb(iq,ib)
      dnHB(iq) = dnHB(iq) * BasisPrim%w(iq)
    END DO
    !CALL QML_Write_dnS(dnHB(1),6,info='dnHB',all_type=.TRUE.)
    !write(6,*) 'coucou dnHB: done',ib ; flush(6)
    DO jb=1,nb
      !H(jb,ib) = dot_product(BasisPrim%d0gb(:,jb),HB(:))
      dnHij = dot_product(BasisPrim%d0gb(:,jb),dnHB(:))
      CALL QML_sub_dnS_TO_dnMat(dnHij,dnH,jb,ib)
    END DO
  END DO
  !CALL Write_RMat(H,6,5,name_info='H')


  CALL QML_DIAG_dnMat(dnH,dnHdiag,dnVec,dnVecProj)
  !CALL Write_RMat(dnHdiag%d0*auTOcm_inv,6,5,name_info='Ene')
  CALL Write_RMat(dnVec%d0(:,1:6),6,5,name_info='Vec')
  write(out_unitp,*) 'NAC'
  !CALL QML_Write_dnMat(dnVecProj,nio=out_unitp)
  CALL Write_RMat(dnVecProj%d1(1:6,1:6,1),6,5,name_info='NAC')


  !CALL diagonalization(H,Ene,Vec,nb)
  DO ib=1,nb
    Ene(ib) = dnHdiag%d0(ib,ib)
  END DO
  write(out_unitp,*) 'Ene'
  CALL Write_RVec(Ene(1:6)*auTOcm_inv,6,5,name_info='Ene')



END PROGRAM TEST_adia
