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
!    Copyright 2017 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================
PROGRAM main_pot
!$ USE omp_lib
USE mod_dnS
  IMPLICIT NONE

  real (kind=8),      allocatable     :: Q(:)
  real (kind=8),      allocatable     :: V(:,:)
  character (len=16)                  :: pot_name

  integer                             :: i,k,ndim,nsurf,option,nb_eval,maxth


  maxth = 1
  !$ maxth           = omp_get_max_threads()

  write(6,*) 'TEST_driver. number of threads:',maxth

  ndim     = 2
  nsurf    = 3
  option   = 0
  nb_eval  = 10**6
  pot_name = 'phenol'

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,pot_name,option) &
!$OMP   PRIVATE(i,k,Q,V) &
!$OMP   NUM_THREADS(maxth)

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    CALL  random_number(Q)
    Q = Q + (/1._8,-0.5_8 /)
    CALL sub_model1_V(V,Q,ndim,nsurf,pot_name,option)

    !write(6,*) Q,(V(k,k),k=1,nsurf)
  END DO
!$OMP   END DO

!$OMP   END PARALLEL

  deallocate(Q)
  deallocate(V)

END PROGRAM main_pot
