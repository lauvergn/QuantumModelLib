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
MODULE AdiaChannels_MakeHinact_m
  USE QDUtil_NumParameters_m
  IMPLICIT NONE

CONTAINS

  SUBROUTINE QML_MakeHinact(Qact,QModel)
    USE QDUtil_m,         ONLY : diagonalization, Write_Mat, Write_Vec
    USE ADdnSVM_m
    USE Model_m
    USE AdiaChannels_Basis_m
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


  CALL alloc_dnMat(dnVec,    nsurf=nb,nVar=ndim_act,nderiv=nderiv,name_var='dnVec')
  CALL alloc_dnMat(dnHdiag,  nsurf=nb,nVar=ndim_act,nderiv=nderiv,name_var='dnHdiag')
  CALL alloc_dnMat(dnVecProj,nsurf=nb,nVar=ndim_act,nderiv=nderiv,name_var='dnVecProj')

  IF (.NOT. allocated(QModel%QM%Vec0)) allocate(QModel%QM%Vec0)
  IF (Check_NotAlloc_dnMat(QModel%QM%Vec0,nderiv=0)) THEN
!$OMP CRITICAL (CRIT_QML_MakeHinact)
    CALL alloc_dnMat(QModel%QM%Vec0,nsurf=nb,nVar=ndim_act,nderiv=0)
    CALL diagonalization(dnH%d0,Ene,QModel%QM%Vec0%d0,nb,sort=1,phase=.TRUE.)
    write(out_unit,*) 'init Vec0 done'
!$OMP END CRITICAL (CRIT_QML_MakeHinact)
  END IF

  CALL DIAG_dnMat(dnH,dnHdiag,dnVec,dnVecProj,dnVec0=QModel%QM%Vec0)
  !CALL Write_Mat(dnVec%d0(:,1:QModel%QM%nb_Channels),out_unit,6,info='Vec')
  write(out_unit,*) 'NAC'
  !CALL Write_dnMat(dnVecProj,nio=out_unit)
  CALL Write_Mat(dnVecProj%d1(1:QModel%QM%nb_Channels,1:QModel%QM%nb_Channels,1),out_unit,6,info='NAC')

  DO ib=1,nb
    Ene(ib) = dnHdiag%d0(ib,ib)
  END DO
  write(out_unit,*) 'Ene'
  CALL Write_Vec(Ene(1:QModel%QM%nb_Channels)*auTOcm_inv,out_unit,6,info='Ene')

  END SUBROUTINE QML_MakeHinact

END MODULE AdiaChannels_MakeHinact_m
