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
MODULE QMLLib_UtilLib_m
  USE QDUtil_NumParameters_m, out_unitp => out_unit, in_unitp => in_unit
  !$ USE omp_lib
  IMPLICIT NONE

  PRIVATE


  PUBLIC :: time_perso
  PUBLIC :: Find_Label,NFind_Label,FFind_Label


INTERFACE time_perso
  MODULE PROCEDURE QML_time_perso
END INTERFACE
INTERFACE Find_Label
  MODULE PROCEDURE QML_Find_Label
END INTERFACE
INTERFACE NFind_Label
  MODULE PROCEDURE QML_NFind_Label
END INTERFACE
INTERFACE FFind_Label
  MODULE PROCEDURE QML_FFind_Label
END INTERFACE
CONTAINS

  SUBROUTINE QML_time_perso(name)
  IMPLICIT NONE

    character (len=*) :: name


    integer           :: tab_time(8) = 0
    real (kind=Rkind) :: t_real
    integer           :: count,count_work,freq
    real              :: t_cpu

    integer, save     :: count_old,count_ini
    real, save        :: t_cpu_old,t_cpu_ini
    integer           :: seconds,minutes,hours,days
    logical, save     :: begin = .TRUE.


!$OMP    CRITICAL (QML_time_perso_CRIT)

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


     t_real = real(count_work,kind=Rkind)/real(freq,kind=Rkind)
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

     t_real = real(count_work,kind=Rkind)/real(freq,kind=Rkind)
     write(out_unitp,41) t_real
 41  format('  Total real (s): ',f18.3)
     write(out_unitp,42) days,hours,minutes,seconds
 42  format('  Total real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s')
     write(out_unitp,43) t_cpu-t_cpu_ini
 43  format('  Total cpu (s): ',f18.3)

     flush(out_unitp)
     !============================================

     count_old = count
     t_cpu_old = t_cpu

!$OMP   END CRITICAL (QML_time_perso_CRIT)


  END SUBROUTINE QML_time_perso

  SUBROUTINE QML_Find_Label(nio,label,located)
    USE QDUtil_m,         ONLY : TO_string
    IMPLICIT NONE

    character (len=*) :: label
        logical :: located
        integer :: nio

        character (len=len(label)) :: labelR
        character (len=name_len) :: format_label

        integer, save :: i_line = 0

        located = .FALSE.
        format_label='(A' // TO_string(len(label)) // ')'

        !write(out_unitp,*) 'format_label',format_label
  !     write(out_unitp,*) 'label to find:',label,located
        DO
          i_line = i_line + 1

          read(nio,format_label,end=999,eor=888,err=888,advance='no') labelR

          located = (labelR .EQ. label)
          !located = verify(label,labelR) == 0
          !write(out_unitp,*) i_line,located,labelR
          IF (located) EXIT
          read(nio,*,end=999,err=888)

   888    CONTINUE
        END DO


   999  CONTINUE
        i_line = 0

        !write(out_unitp,*) 'Find_Label: ',label,located
  END SUBROUTINE QML_Find_Label

  SUBROUTINE QML_NFind_Label(nio,label,located,iformat)
    USE QDUtil_m,         ONLY : TO_string
    IMPLICIT NONE

    character (len=*) :: label
        integer :: iformat
        character (len=Name_len) :: fformat
        logical :: located
        integer :: nio

        character (len=iformat) :: labelR
        integer, save :: i_line = 0


        located = .FALSE.

        IF (iformat < 1) THEN
          write(out_unitp,*) ' ERROR in NFind_Label'
          write(out_unitp,*) ' iformat < 1',iformat
          STOP
        END IF
        fformat='(A' // TO_string(iformat) // ')'
        !write(out_unitp,*) iformat,fformat
        !write(out_unitp,*) 'label to find:',label,located
        DO
          i_line = i_line + 1
          read(nio,fformat,end=999,eor=888,err=888,advance='no') labelR

          located = verify(label,labelR) == 0
          !write(out_unitp,*) i_line,located,labelR
          IF (located) EXIT
          read(nio,*,end=999,err=888)

   888    CONTINUE
        END DO


   999  CONTINUE
        i_line = 0
  end subroutine QML_NFind_Label

  SUBROUTINE QML_FFind_Label(nio,label,located,fformat)
        IMPLICIT NONE
        character (len=*)   :: label
        character (len=*)   :: fformat
        logical             :: located
        integer             :: nio

        character (len=256) :: labelR
        integer, save       :: i_line = 0


        located = .FALSE.

        !write(out_unitp,*) 'label to find:',label,located
        DO
          i_line = i_line + 1
          read(nio,fformat,end=999,eor=888,err=888,advance='no') labelR


          located = verify(label,labelR) == 0
          !write(out_unitp,*) i_line,located,labelR
          IF (located) EXIT
          read(nio,*,end=999,err=888)

   888    CONTINUE
        END DO


   999  CONTINUE
        i_line = 0
        IF (located) read(nio,*)

  END SUBROUTINE QML_FFind_Label

END MODULE QMLLib_UtilLib_m
