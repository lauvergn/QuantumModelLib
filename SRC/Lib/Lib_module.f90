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
!    Copyright 2016 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!    This particular module has been modified from ElVibRot-Tnum:
!       <http://pagesperso.lcp.u-psud.fr/lauvergnat/ElVibRot/ElVibRot.html>
!===========================================================================
!===========================================================================
MODULE mod_Lib
USE mod_NumParameters
!$ USE omp_lib
IMPLICIT NONE

CONTAINS

  !!@description: Change the case of a string (default lowercase)
  !!@param: name_string the string
  !!@param: lower If the variable is present and its value is F,
  !!              the string will be converted into a uppercase string, otherwise,
  !!              it will be convert into a lowercase string.
  SUBROUTINE string_uppercase_TO_lowercase(name_string,lower)
  IMPLICIT NONE

   character (len=*), intent(inout)  :: name_string
   logical, optional  :: lower

   logical  :: lower_loc
   integer  :: i,ascii_char

   IF (present(lower)) THEN
     lower_loc = lower
   ELSE
     lower_loc = .TRUE.
   END IF

   !write(out_unitp,*) 'name_string: ',name_string
   IF (lower_loc) THEN ! uppercase => lowercase
     DO i=1,len_trim(name_string)
       ascii_char = iachar(name_string(i:i))
       IF (ascii_char >= 65 .AND. ascii_char <= 90)                 &
                           name_string(i:i) = achar(ascii_char+32)
     END DO
   ELSE ! lowercase => uppercase
     DO i=1,len_trim(name_string)
       ascii_char = iachar(name_string(i:i))
       IF (ascii_char >= 97 .AND. ascii_char <= 122)                 &
                            name_string(i:i) = achar(ascii_char-32)

     END DO
   END IF
   !write(out_unitp,*) 'name_string: ',name_string


  END SUBROUTINE string_uppercase_TO_lowercase

  FUNCTION strdup(string)
  IMPLICIT NONE

   character (len=*), intent(in)   :: string
   character (len=:), allocatable  :: strdup

   allocate(character(len=len_trim(string)) :: strdup)
   strdup = trim(string)

  END FUNCTION strdup
  FUNCTION int_TO_char(i)
  USE mod_NumParameters

    character (len=:), allocatable  :: int_TO_char
    integer, intent(in)             :: i

    character (len=:), allocatable  :: name_int
    integer :: clen

    ! first approximated size of name_int
    IF (i == 0) THEN
      clen = 1
    ELSE IF (i < 0) THEN
      clen = int(log10(abs(real(i,kind=Rkind))))+2
    ELSE
      clen = int(log10(real(i,kind=Rkind)))+1
    END IF

    ! allocate name_int
    allocate(character(len=clen) :: name_int)

    ! write i in name_int
    write(name_int,'(i0)') i

    ! transfert name_int in int_TO_char
    int_TO_char = strdup(name_int)

    ! deallocate name_int
    deallocate(name_int)

  END FUNCTION int_TO_char
  SUBROUTINE sub_Format_OF_Line(wformat,nb_line,max_col,cplx,       &
                                Rformat,name_info)
  USE mod_NumParameters
  IMPLICIT NONE

   character (len=*), optional  :: Rformat
   character (len=*), optional  :: name_info

   character (len=Name_longlen) :: wformat,name_info_loc
   integer                      :: nb_line,max_col
   logical                      :: cplx


   character (len=Name_longlen) :: iformat,cformat
   character (len=Name_longlen) :: NMatformat


   IF (present(name_info)) THEN
     name_info_loc = '(2x,"' // trim(adjustl(name_info)) // ' ",'
   ELSE
     name_info_loc = '('
   END IF

   IF (present(Rformat)) THEN
     IF (len_trim(Rformat) > 10) THEN
       write(out_unitp,*) ' ERROR in sub_Format_OF_Line'
       write(out_unitp,*) ' The format (len_trim) in "Rformat" is too long',len_trim(Rformat)
       write(out_unitp,*) ' Rformat: ',Rformat
       STOP
     END IF
       IF (cplx) THEN
         NMatformat = "'('," // trim(adjustl(Rformat)) //           &
                      ",' +i'," // trim(adjustl(Rformat)) // ",')'"
       ELSE
         NMatformat = trim(adjustl(Rformat))
       END IF
   ELSE
       IF (cplx) THEN
         NMatformat = trim(adjustl(CMatIO_format))
       ELSE
         NMatformat = trim(adjustl(RMatIO_format))
       END IF
   END IF

   IF (nb_line > 0) THEN

       write(iformat,*) int(log10(real(nb_line,kind=Rkind)))+1
       write(cformat,*) max_col

       wformat = trim(adjustl(name_info_loc)) // '1x,' //           &
                 'i' // trim(adjustl(iformat)) // ',2x,' //         &
                 trim(adjustl(cformat)) // '(' //                   &
                 trim(adjustl(NMatformat)) // ',1x))'
   ELSE
       write(cformat,*) max_col

         wformat = trim(adjustl(name_info_loc)) //                  &
                   trim(adjustl(cformat)) // '(' //                 &
                   trim(adjustl(NMatformat)) // ',1x))'


   END IF

  !write(out_unitp,*) 'format?: ',trim(wformat)
  END SUBROUTINE sub_Format_OF_Line

  SUBROUTINE Write_RMat(f,nio,nbcol1,Rformat,name_info)
  USE mod_NumParameters
  IMPLICIT NONE

     character (len=*), optional :: Rformat
     character (len=*), optional :: name_info

     integer, intent(in)         :: nio,nbcol1
     real(kind=8), intent(in) :: f(:,:)

     integer         :: nl,nc
     integer i,j,nb,nbblocs,nfin,nbcol
     character (len=:), allocatable :: wformat

     nl = size(f,dim=1)
     nc = size(f,dim=2)
     !write(out_unitp,*) 'nl,nc,nbcol',nl,nc,nbcol
     nbcol = nbcol1
     IF (nbcol > 10) nbcol=10
     nbblocs=int(nc/nbcol)
     IF (nbblocs*nbcol == nc) nbblocs=nbblocs-1


     IF (present(name_info)) THEN
       wformat = strdup( '(" ' // trim(name_info) // ': ",' )
     ELSE
       wformat = strdup( '(' )
     END IF
     IF (present(Rformat)) THEN
       wformat = strdup( wformat // 'i0,x,10' // trim(Rformat) // ')' )
     ELSE
       wformat = strdup( wformat // 'i0,x,10f15.6)' )
     END IF

     DO nb=0,nbblocs-1
       DO j=1,nl
         write(nio,wformat) j,(f(j,i+nb*nbcol),i=1,nbcol)
       END DO
       IF (nl > 1 ) write(nio,*)
     END DO
     DO j=1,nl
       nfin=nc-nbcol*nbblocs
       write(nio,wformat) j,(f(j,i+nbcol*nbblocs),i=1,nfin)
     END DO

  END SUBROUTINE Write_RMat
  SUBROUTINE Write_RVec(l,nio,nbcol1,Rformat,name_info)
  USE mod_NumParameters
  IMPLICIT NONE

    character (len=*), optional  :: Rformat
    character (len=*), optional  :: name_info

    integer, intent(in)          :: nio,nbcol1
    real(kind=Rkind), intent(in) :: l(:)

    integer           :: n,i,nb,nbblocs,nfin,nbcol
    character (len=Name_longlen) :: wformat

    n = size(l)
    !write(out_unitp,*) 'n,nbcol',n,nbcol
    nbcol = nbcol1
    IF (nbcol > 10) nbcol=10
    nbblocs=int(n/nbcol)
    IF (nbblocs*nbcol == n) nbblocs=nbblocs-1


    IF (present(Rformat)) THEN
      IF (present(name_info)) THEN
         CALL sub_Format_OF_Line(wformat,0,nbcol,.FALSE.,Rformat,name_info)
       ELSE
         CALL sub_Format_OF_Line(wformat,0,nbcol,.FALSE.,Rformat=Rformat)
       END IF
     ELSE
       IF (present(name_info)) THEN
         CALL sub_Format_OF_Line(wformat,0,nbcol,.FALSE.,name_info=name_info)
       ELSE
         CALL sub_Format_OF_Line(wformat,0,nbcol,.FALSE.)
       END IF
     END IF


     DO nb=0,nbblocs-1
       write(nio,wformat) (l(i+nb*nbcol),i=1,nbcol)
     END DO
     nfin=n-nbcol*nbblocs
     write(nio,wformat) (l(i+nbcol*nbblocs),i=1,nfin)

  END SUBROUTINE Write_RVec

  SUBROUTINE Read_RMat(f,nio,nbcol,err)

  integer, intent(in)             :: nio,nbcol
  integer, intent(inout)          :: err

  real(kind=Rkind), intent(inout) :: f(:,:)

  integer i,j,jj,nb,nbblocs,nfin,nl,nc

!$OMP    CRITICAL (Read_RMat_CRIT)

  nl = size(f,dim=1)
  nc = size(f,dim=2)
! write(out_unitp,*) 'nl,nc,nbcol',nl,nc,nbcol


  nbblocs=int(nc/nbcol)

  IF (nbblocs*nbcol == nc) nbblocs=nbblocs-1
  err = 0

  DO nb=0,nbblocs-1

      DO j=1,nl
        read(nio,*,IOSTAT=err) jj,(f(j,i+nb*nbcol),i=1,nbcol)
        IF (err /= 0) EXIT
      END DO

      IF (err /= 0) EXIT

      IF (nl > 1) read(nio,*,IOSTAT=err)
      IF (err /= 0) EXIT

  END DO

  IF (err == 0) THEN
    DO j=1,nl
      nfin=nc-nbcol*nbblocs
      read(nio,*,IOSTAT=err) jj,(f(j,i+nbcol*nbblocs),i=1,nfin)
      IF (err /= 0) EXIT
    END DO
  END IF

  IF (err /= 0) THEN
    CALL Write_RMat(f,out_unitp,nbcol)
    write(out_unitp,*) ' ERROR in Read_RMat'
    write(out_unitp,*) '  while reading a matrix'
    write(out_unitp,*) '  end of file or end of record'
    write(out_unitp,*) '  The matrix paramters: nl,nc,nbcol',nl,nc,nbcol
    write(out_unitp,*) ' Check your data !!'
  END IF

!$OMP    END CRITICAL (Read_RMat_CRIT)

  END SUBROUTINE Read_RMat

  SUBROUTINE Read_RVec(l,nio,nbcol,err)

  integer, intent(in)                :: nio,nbcol
  real(kind=Rkind), intent(inout)    :: l(:)
  integer, intent(inout)             :: err

  integer :: n,i,nb,nbblocs,nfin

!$OMP    CRITICAL (Read_RVec_CRIT)

  n = size(l,dim=1)
  nbblocs=int(n/nbcol)
  err = 0


  IF (nbblocs*nbcol == n) nbblocs=nbblocs-1

  DO nb=0,nbblocs-1
    read(nio,*,IOSTAT=err) (l(i+nb*nbcol),i=1,nbcol)
    IF (err /= 0) EXIT
  END DO

  nfin=n-nbcol*nbblocs
  read(nio,*,IOSTAT=err) (l(i+nbcol*nbblocs),i=1,nfin)

  IF (err /= 0) THEN
    write(out_unitp,*) ' ERROR in Read_RVec'
    write(out_unitp,*) '  while reading a vector'
    write(out_unitp,*) '  end of file or end of record'
    write(out_unitp,*) '  The vector paramters: n,nbcol',n,nbcol
    write(out_unitp,*) ' Check your data !!'
  END IF

!$OMP    END CRITICAL (Read_RVec_CRIT)

  END SUBROUTINE Read_RVec

  SUBROUTINE time_perso(name)
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


     flush(out_unitp)
     !============================================

     count_old = count
     t_cpu_old = t_cpu

!$OMP   END CRITICAL (time_perso_CRIT)


  END SUBROUTINE time_perso

END MODULE mod_Lib
