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
PROGRAM main_pot
  IMPLICIT NONE

  CALL test2_Vib_adia(100)

stop
  CALL test_HBond()
  CALL test_1DSOC_1S1T()
  CALL test_PSB3()
  CALL test_Phenol_Dia(10**7)
  CALL test_Phenol_ADia()
  CALL test_henonheiles(10**7)
  CALL test_Vib_adia(1000)

END PROGRAM main_pot
SUBROUTINE test2_Vib_adia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'
  maxth = 1
  write(6,*) '============================================================'
  write(6,*) '  Vibrational adiabatic separation'
  write(6,*) 'TEST_driver. number of threads:',maxth

  pot_name  = 'read_model'
  ndim      = 0 ! it would be initialized
  nsurf     = 0 ! it would be initialized
  option    = -1
  adiabatic = .FALSE.
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
  write(6,*) 'ndim,nsurf',ndim,nsurf

  write(6,*) ' Test V, G and NAC'
  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  read(5,*) Q

  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)
  write(6,*) 'Q',Q
  write(6,*) 'Vadia',(V(i,i),i=1,nsurf)
  write(6,*) 'Gadia',(G(i,i,1:ndim),i=1,nsurf)
  write(6,*) 'NAC'
  DO i=1,nsurf
    write(6,*) i,NAC(:,i,1:ndim)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(G)
  deallocate(NAC)

  write(6,*) ' Test of evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

  Q = 2.5d0
  DO i=1,nb_eval
    Q = Q + 0.05d0

    CALL sub_Qmodel_V(V,Q)

    write(6,*) 'pot:',Q,(V(k,k),k=1,nsurf)
  END DO

  deallocate(Q)
  deallocate(V)

  CALL time_perso('Test ' // pot_name)
  write(6,*) '============================================================'
  write(6,*) '============================================================'

END SUBROUTINE test2_Vib_adia
SUBROUTINE test_1DSOC_1S1T
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'
  nb_eval = 1

  pot_name = '1DSOC_1S1T'
  write(6,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name

  ndim      = 1
  nsurf     = 4
  option    = 0
  adiabatic = .TRUE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(Vec(nsurf,nsurf))


  Q(:) = (/8.5_real64 /)

  CALL sub_Qmodel_VVec(V,Vec,Q)
  write(6,*) ' Diabatic potential as a 4x4 matrix:'
  write(6,'(4f12.8)') V
  write(6,*) ' Adiabatic vectors (in column) as a 4x4 matrix:'
  write(6,'(4f12.8)') transpose(Vec)

  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(Q)
  deallocate(Vec)


  write(6,*) '============================================================'
  write(6,*) '============================================================'


END SUBROUTINE test_1DSOC_1S1T
SUBROUTINE test_PSB3
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(6,*) '============================================================'
  write(6,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(6,*) '============================================================'
  write(6,*) 'NTEST_driver. number of threads:',maxth
  write(6,*) '============================================================'
  write(6,*) '============================================================'
  pot_name = 'PSB3'
  write(6,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  ndim     = 3
  nsurf    = 2
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))


  Q(:) = (/0.1_real64,-3.14_real64,0.0_real64/)

  CALL sub_Qmodel_V(V,Q)
  write(6,*) ' Diabatic potential as a 2x2 matrix:'
  write(6,'(2f12.8)') V

  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(Q)


  write(6,*) '============================================================'
  write(6,*) '============================================================'

END SUBROUTINE test_PSB3

SUBROUTINE test_HBond
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(6,*) '============================================================'
  write(6,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(6,*) '============================================================'
  write(6,*) 'NTEST_driver. number of threads:',maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'
  pot_name = 'HBond'
  write(6,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  ndim     = 2
  nsurf    = 1
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))


  Q(:) = (/2._real64,0._real64/)

  CALL sub_Qmodel_V(V,Q)
  write(6,*) V


  allocate(G(nsurf,nsurf,ndim))
  allocate(h(nsurf,nsurf,ndim,ndim))

  CALL sub_Qmodel_VGH(V,G,H,Q)
  write(6,*) 'gradiant',G
  write(6,*) 'hessian',H


  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(Q)


  write(6,*) '============================================================'
  write(6,*) '============================================================'



END SUBROUTINE test_HBond

SUBROUTINE test_Phenol_Dia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(6,*) '============================================================'
  write(6,*) 'NTEST_driver. number of threads:',maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'
  pot_name = 'phenol'
  write(6,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  ndim     = 2
  nsurf    = 3
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

  write(6,*) 'ndim,nsurf',ndim,nsurf

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(GGdef(ndim,ndim))

  ! write GGdef
  write(6,*) 'Get and Write the metric Tensor'
  CALL get_Qmodel_GGdef(GGdef)
  DO i=1,ndim
    write(6,*) i,GGdef(:,i)
  END DO


  Q = (/1.0_real64,-0.5_real64 /)
  CALL sub_Qmodel_V(V,Q)
  write(6,*) ' Diabatic potential as a 3x3 matrix:'
  write(6,'(3f12.8)') V

  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(GGdef)
  deallocate(Q)


  CALL time_perso('Test ' // pot_name)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,pot_name,option) &
!$OMP   PRIVATE(i,Q,V) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    CALL  random_number(Q)
    Q = Q + (/1.0_real64,-0.5_real64 /)
    CALL sub_Qmodel_V(V,Q)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL time_perso('Test ' // pot_name)
  write(6,*) '============================================================'
  write(6,*) '============================================================'


END SUBROUTINE test_Phenol_Dia
SUBROUTINE test_Phenol_ADia
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(6,*) '============================================================'
  write(6,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(6,*) '============================================================'
  write(6,*) 'NTEST_driver. number of threads:',maxth



  write(6,*) '============================================================'
  write(6,*) '============================================================'
  pot_name = 'phenol'
  write(6,*) ' Test of the adiatic potential ',pot_name


  ndim     = 2
  nsurf    = 3
  option   = 0
  adiabatic = .TRUE.

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(g(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)  ! a new initialization


  Q = (/1._real64,-0.5_real64 /)
  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)

  write(6,*) ' Adiabatic potential as a 3x3 matrix:'
  write(6,'(3f12.8)') V
  write(6,*) ' Adiabatic gradient. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(6,*) ' Component:',i
    write(6,'(3f12.8)') g(:,:,i)
  END DO

  write(6,*) ' NAC. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(6,*) ' Component:',i
    write(6,'(3f12.8)') NAC(:,:,i)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(g)
  deallocate(NAC)

  write(6,*) '============================================================'
  write(6,*) '============================================================'


END SUBROUTINE test_Phenol_ADia
SUBROUTINE test_henonheiles(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(6,*) '============================================================'
  write(6,*) 'NTEST_driver. number of threads:',maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'

  ndim     = 6
  nsurf    = 1
  option   = 0
  adiabatic = .TRUE.
  pot_name = 'henonheiles'
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)  ! a new initialization

  write(6,*) '============================================================'
  write(6,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,pot_name,option) &
!$OMP   PRIVATE(i,k,Q,V) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    CALL  random_number(Q)

    CALL sub_Qmodel_V(V,Q)

    !write(6,*) Q,(V(k,k),k=1,nsurf)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL time_perso('Test ' // pot_name)
  write(6,*) '============================================================'
  write(6,*) '============================================================'

END SUBROUTINE test_henonheiles

SUBROUTINE test_Vib_adia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(6,*) '============================================================'
  write(6,*) '  Vibrational adiabatic separation'
  write(6,*) 'TEST_driver. number of threads:',maxth

  pot_name  = 'read_model'
  ndim      = 0 ! it would be initialized
  nsurf     = 0 ! it would be initialized
  option    = -1
  adiabatic = .FALSE.
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
  write(6,*) 'ndim,nsurf',ndim,nsurf

  write(6,*) ' Test V, G and NAC'
  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  CALL  random_number(Q)
  Q = Q + 4.d0
  !Q=7.d0

  ! calculation of the adiabatic potential (as a matrix)
  ! calculation of the gradient and the NAC

  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)
  write(6,*) 'Q',Q
  write(6,*) 'Vadia',(V(i,i),i=1,nsurf)
  write(6,*) 'Gadia',(G(i,i,1:ndim),i=1,nsurf)
  write(6,*) 'NAC'
  DO i=1,nsurf
    write(6,*) i,NAC(:,i,1:ndim)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(G)
  deallocate(NAC)

  write(6,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,pot_name,option) &
!$OMP   PRIVATE(i,k,Q,V) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    CALL  random_number(Q)
    Q = Q + 4.d0

    CALL sub_Qmodel_V(V,Q)

    !write(6,*) Q,(V(k,k),k=1,nsurf)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL time_perso('Test ' // pot_name)
  write(6,*) '============================================================'
  write(6,*) '============================================================'

END SUBROUTINE test_Vib_adia

  SUBROUTINE time_perso(name)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
  IMPLICIT NONE

    character (len=*) :: name


    integer             :: tab_time(8) = 0
    real (kind=real64)  :: t_real
    integer             :: count,count_work,freq
    real                :: t_cpu
    integer             :: seconds,minutes,hours,days

    integer, save       :: count_old,count_ini
    real,    save       :: t_cpu_old,t_cpu_ini
    logical, save       :: begin = .TRUE.


!$OMP    CRITICAL (time_perso_CRIT)

    CALL date_and_time(values=tab_time)
    write(6,21) name,tab_time(5:8),tab_time(3:1:-1)
 21 format('     Time and date in ',a,' : ',i2,'h:',                &
            i2,'m:',i2,'.',i3,'s, the ',i2,'/',i2,'/',i4)

     CALL system_clock(count=count,count_rate=freq)
     call cpu_time(t_cpu)

     IF (begin) THEN
       begin = .FALSE.
       count_old = count
       count_ini = count
       t_cpu_old = t_cpu
       t_cpu_ini = t_cpu
     END IF


     !============================================
     !cpu time in the subroutine: "name"

     count_work = count-count_old
     seconds = count_work/freq

     minutes = seconds/60
     seconds = mod(seconds,60)
     hours   = minutes/60
     minutes = mod(minutes,60)
     days    = hours/24
     hours   = mod(hours,24)


     t_real = real(count_work,kind=real64)/real(freq,kind=real64)
     write(6,31) t_real,name
 31  format('        real (s): ',f18.3,' in ',a)
     write(6,32) days,hours,minutes,seconds,name
 32  format('        real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s in ',a)

     write(6,33) t_cpu-t_cpu_old,name
 33  format('        cpu (s): ',f18.3,' in ',a)


     !============================================
     !Total cpu time

     count_work = count-count_ini
     seconds = count_work/freq

     minutes = seconds/60
     seconds = mod(seconds,60)
     hours   = minutes/60
     minutes = mod(minutes,60)
     days    = hours/24
     hours   = mod(hours,24)

     t_real = real(count_work,kind=real64)/real(freq,kind=real64)
     write(6,41) t_real
 41  format('  Total real (s): ',f18.3)
     write(6,42) days,hours,minutes,seconds
 42  format('  Total real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s')
     write(6,43) t_cpu-t_cpu_ini
 43  format('  Total cpu (s): ',f18.3)

 51  format(a,i10,a,a)


     flush(6)
     !============================================

     count_old = count
     t_cpu_old = t_cpu

!$OMP   END CRITICAL (time_perso_CRIT)


  END SUBROUTINE time_perso
