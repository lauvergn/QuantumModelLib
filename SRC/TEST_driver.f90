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
!      last modification, 11/07/2019 DML
!
!===========================================================================
!===========================================================================
PROGRAM main_pot
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=8),      allocatable     :: Q(:)
  real (kind=8),      allocatable     :: GGdef(:,:)
  real (kind=8),      allocatable     :: V(:,:)
  real (kind=8),      allocatable     :: g(:,:,:)
  real (kind=8),      allocatable     :: h(:,:,:,:)
  real (kind=8),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth


  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(6,*) 'NTEST_driver. number of threads:',maxth

  write(6,*) '============================================================'
  write(6,*) '============================================================'



  nb_eval  = 10**7

  write(6,*) '============================================================'
  write(6,*) '============================================================'
  pot_name = 'phenol'
  write(6,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  ndim     = 2
  nsurf    = 3
  option   = 0
  adiabatic = .FALSE.
  pot_name = 'phenol'

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(GGdef(ndim,ndim))

  ! write GGdef
  write(6,*) 'Get and Write the metric Tensor'
  CALL get_Qmodel_GGdef(GGdef)
  DO i=1,ndim
    write(6,*) i,GGdef(:,i)
  END DO


  Q = (/1._8,-0.5_8 /)
  CALL sub_Qmodel_V(V,Q)
  write(6,*) ' Diabatic potential as a 3x3 matrix:'
  write(6,'(3f12.8)') V

  deallocate(Q)
  deallocate(V)
  deallocate(GGdef)

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
    Q = Q + (/1._8,-0.5_8 /)
    CALL sub_Qmodel_V(V,Q)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL time_perso('Test ' // pot_name)
  write(6,*) '============================================================'
  write(6,*) '============================================================'


  write(6,*) '============================================================'
  write(6,*) '============================================================'
  pot_name = 'phenol'
  write(6,*) ' Test of the adiatic potential ',pot_name


  ndim     = 2
  nsurf    = 3
  option   = 0
  adiabatic = .TRUE.
  pot_name = 'phenol'

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(g(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)  ! a new initialization


  Q = (/1._8,-0.5_8 /)
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


  write(6,*) '============================================================'
  write(6,*) '============================================================'

  ndim     = 6
  nsurf    = 1
  option   = 0
  adiabatic = .TRUE.
  nb_eval  = 10**8
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

END PROGRAM main_pot
  SUBROUTINE time_perso(name)
  IMPLICIT NONE

    character (len=*) :: name


    integer           :: tab_time(8) = 0
    real (kind=8)     :: t_real
    integer           :: count,count_work,freq
    real              :: t_cpu
    integer           :: seconds,minutes,hours,days

    integer, save     :: count_old,count_ini
    real,    save     :: t_cpu_old,t_cpu_ini
    logical, save     :: begin = .TRUE.


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


     t_real = real(count_work,kind=8)/real(freq,kind=8)
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

     t_real = real(count_work,kind=8)/real(freq,kind=8)
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
