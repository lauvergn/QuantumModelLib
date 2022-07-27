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
    write(out_unitp,21) name,tab_time(5:8),tab_time(3:1:-1)
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
     write(out_unitp,31) t_real,name
 31  format('        real (s): ',f18.3,' in ',a)
     write(out_unitp,32) days,hours,minutes,seconds,name
 32  format('        real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s in ',a)

     write(out_unitp,33) t_cpu-t_cpu_old,name
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
     write(out_unitp,41) t_real
 41  format('  Total real (s): ',f18.3)
     write(out_unitp,42) days,hours,minutes,seconds
 42  format('  Total real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s')
     write(out_unitp,43) t_cpu-t_cpu_ini
 43  format('  Total cpu (s): ',f18.3)

 51  format(a,i10,a,a)


     flush(6)
     !============================================

     count_old = count
     t_cpu_old = t_cpu

!$OMP   END CRITICAL (time_perso_CRIT)


  END SUBROUTINE time_perso
