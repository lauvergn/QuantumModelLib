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
PROGRAM main_pot
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=8),      allocatable     :: Q(:)
  real (kind=8),      allocatable     :: Q0(:)
  real (kind=8),      allocatable     :: V(:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,ndim,nsurf,option,nb_eval,maxth


  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unitp,*) 'NTEST_driver. number of threads:',maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

  ndim      = 6
  nsurf     = 1
  option    = 0
  adiabatic = .TRUE.
  nb_eval   = 10**6
  pot_name  = 'HONO'
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)  ! a new initialization

  allocate(Q0(ndim))
  CALL get_Qmodel_Q0(Q0,0) ! the option of get_Qmodel_Q0 enables to select different reference geometries.

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,Q0) &
!$OMP   PRIVATE(i,Q,V) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  write(out_unitp,*) 'alloc Q and V'
!$OMP BARRIER

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval

    CALL  random_number(Q)
    Q = Q0 + (Q-0.5d0)/10.d0

    CALL sub_Qmodel_V(V,Q)

    !write(out_unitp,*) Q,(V(k,k),k=1,nsurf)
  END DO
!$OMP   END DO

!$OMP BARRIER
  deallocate(Q)
  deallocate(V)
  write(out_unitp,*) 'dealloc Q and V'

!$OMP   END PARALLEL

  deallocate(Q0)


  CALL time_perso('Test ' // pot_name)
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

END PROGRAM main_pot
SUBROUTINE QML_time_perso(name_sub)
  USE QDUtil_m, ONLY : time_perso
  IMPLICIT NONE

  character (len=*) :: name_sub


  !$OMP    CRITICAL (QML_time_perso_CRIT)
  CALL time_perso(name_sub)
  !$OMP   END CRITICAL (QML_time_perso_CRIT)

END SUBROUTINE QML_time_perso