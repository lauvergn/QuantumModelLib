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
PROGRAM TEST_gridNonAdia
  !USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
  USE QDUtil_m
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: Q0(:)
  real (kind=Rkind),      allocatable     :: GGdef(:,:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: Vec0(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,iQ,ndim,nsurf,option,nb_eval,maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  nb_eval = maxth
  write(out_unit,*) '============================================================'
  write(out_unit,*) 'NTEST_driver. number of threads:',maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  pot_name = 'phenol'
  write(out_unit,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL QML_time_perso('Test ' // pot_name)

  ndim     = 2
  nsurf    = 3
  option   = 0
  adiabatic = .TRUE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
  CALL set_Qmodel_Phase_Following(.TRUE.)
  CALL set_Qmodel_Phase_Checking(.TRUE.)
  write(out_unit,*) 'ndim,nsurf',ndim,nsurf

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,pot_name,option) &
!$OMP   PRIVATE(i,iQ,Q,Q0,V,g,NAC,Vec0) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(Q0(ndim))
  allocate(V(nsurf,nsurf))
  allocate(Vec0(nsurf,nsurf))
  allocate(g(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  CALL  random_number(Q0)
  Q0 = ONETENTH*Q0 + [1.0_Rkind,-0.5_Rkind]
  Q0(1) = 1.5_Rkind
  CALL sub_Qmodel_VVec(V,Vec0,Q0)

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    DO iQ=0,1000
      Q = Q0 + [ONETENTH**2*iq,ZERO]
      CALL sub_Qmodel_VG_NAC_Vec0(V,G,NAC,Vec0,Q)
      write(660+i,*) Q,V(1,1),V(2,2),V(3,3) ; flush(660+i)
      write(770+i,*) Q,NAC(1,2,:),NAC(1,3,:),NAC(2,3,:) ; flush(770+i)
      write(880+i,*) Q,Vec0 ; flush(660+i)
    END DO
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(Q0)
  deallocate(V)
  deallocate(Vec0)
  deallocate(g)
  deallocate(NAC)

!$OMP   END PARALLEL

  CALL QML_time_perso('Test ' // pot_name)
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'


END PROGRAM TEST_gridNonAdia
