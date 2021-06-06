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
MODULE QMLLib_UtilLib_m
USE QMLLib_NumParameters_m
!$ USE omp_lib
IMPLICIT NONE

  !PRIVATE

  character (len=Line_len), public :: File_path = ''

  PUBLIC :: string_uppercase_TO_lowercase,strdup,string_IS_empty,               &
            int_TO_char,sub_Format_OF_Line
  PUBLIC :: Write_RMat,Write_RVec,Read_RMat,Read_RVec,Init_IdMat
  PUBLIC :: time_perso
  PUBLIC :: err_file_name,make_FileName,file_open2

INTERFACE string_uppercase_TO_lowercase
  MODULE PROCEDURE QML_string_uppercase_TO_lowercase
END INTERFACE
INTERFACE strdup
  MODULE PROCEDURE QML_strdup
END INTERFACE
INTERFACE int_TO_char
  MODULE PROCEDURE QML_int_TO_char
END INTERFACE
INTERFACE string_IS_empty
  MODULE PROCEDURE QML_string_IS_empty
END INTERFACE
INTERFACE sub_Format_OF_Line
  MODULE PROCEDURE QML_sub_Format_OF_Line
END INTERFACE

INTERFACE Write_RMat
  MODULE PROCEDURE QML_Write_RMat
END INTERFACE
INTERFACE Write_RVec
  MODULE PROCEDURE QML_Write_RVec
END INTERFACE
INTERFACE Read_RMat
  MODULE PROCEDURE QML_Read_RMat
END INTERFACE
INTERFACE Read_RVec
  MODULE PROCEDURE QML_Read_RVec
END INTERFACE
INTERFACE Init_IdMat
  MODULE PROCEDURE QML_Init_IdMat
END INTERFACE

INTERFACE time_perso
  MODULE PROCEDURE QML_time_perso
END INTERFACE
INTERFACE err_file_name
  MODULE PROCEDURE QML_err_file_name
END INTERFACE
INTERFACE make_FileName
  MODULE PROCEDURE QML_make_FileName
END INTERFACE
INTERFACE file_open2
  MODULE PROCEDURE QML_file_open2
END INTERFACE

CONTAINS

  !!@description: Change the case of a string (default lowercase)
  !!@param: name_string the string
  !!@param: lower If the variable is present and its value is F,
  !!              the string will be converted into a uppercase string, otherwise,
  !!              it will be convert into a lowercase string.
  SUBROUTINE QML_string_uppercase_TO_lowercase(name_string,lower)
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


  END SUBROUTINE QML_string_uppercase_TO_lowercase

  PURE FUNCTION QML_strdup(string)
  IMPLICIT NONE

   character (len=*), intent(in)   :: string
   character (len=:), allocatable  :: QML_strdup

   allocate(character(len=len_trim(string)) :: QML_strdup)
   QML_strdup = trim(string)

  END FUNCTION QML_strdup
  PURE FUNCTION QML_int_TO_char(i)
  USE QMLLib_NumParameters_m

    character (len=:), allocatable  :: QML_int_TO_char
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

    ! transfert name_int in QML_int_TO_char
    QML_int_TO_char = QML_strdup(name_int)

    ! deallocate name_int
    deallocate(name_int)

  END FUNCTION QML_int_TO_char
  FUNCTION QML_string_IS_empty(String)
    logical                          :: QML_string_IS_empty
    character(len=*), intent(in)     :: String

    QML_string_IS_empty = (len_trim(String) == 0)

  END FUNCTION QML_string_IS_empty

  SUBROUTINE QML_sub_Format_OF_Line(wformat,nb_line,max_col,cplx,       &
                                Rformat,name_info)
  USE QMLLib_NumParameters_m
  IMPLICIT NONE

   character (len=*), optional  :: Rformat
   character (len=*), optional  :: name_info

   character (len=Name_longlen) :: wformat,name_info_loc
   integer                      :: nb_line,max_col
   logical                      :: cplx


   character (len=Name_longlen) :: iformat,cformat
   character (len=Name_longlen) :: NMatformat


   IF (present(name_info)) THEN
     name_info_loc = '(2x,"' // trim(adjustl(name_info)) // ': ",'
   ELSE
     name_info_loc = '('
   END IF

   IF (present(Rformat)) THEN
     IF (len_trim(Rformat) > 10) THEN
       write(out_unitp,*) ' ERROR in QML_sub_Format_OF_Line'
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
  END SUBROUTINE QML_sub_Format_OF_Line

  SUBROUTINE QML_Write_RMat(f,nio,nbcol1,Rformat,name_info)
  USE QMLLib_NumParameters_m
  IMPLICIT NONE

     character (len=*), optional :: Rformat
     character (len=*), optional :: name_info

     integer, intent(in)         :: nio,nbcol1
     real(kind=Rkind), intent(in) :: f(:,:)

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
       wformat = QML_strdup( '(" ' // trim(name_info) // ': ",' )
     ELSE
       wformat = QML_strdup( '(' )
     END IF
     IF (present(Rformat)) THEN
       wformat = QML_strdup( wformat // 'i0,1x,10' // trim(Rformat) // ')' )
     ELSE
       wformat = QML_strdup( wformat // 'i0,1x,10f18.9)' )
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
     flush(nio)

  END SUBROUTINE QML_Write_RMat
  SUBROUTINE QML_Write_RVec(l,nio,nbcol1,Rformat,name_info)
  USE QMLLib_NumParameters_m
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
     flush(nio)

  END SUBROUTINE QML_Write_RVec

  SUBROUTINE QML_Read_RMat(f,nio,nbcol,err)

  integer, intent(in)             :: nio,nbcol
  integer, intent(inout)          :: err

  real(kind=Rkind), intent(inout) :: f(:,:)

  integer i,j,jj,nb,nbblocs,nfin,nl,nc

!$OMP    CRITICAL (QML_Read_RMat_CRIT)

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
    write(out_unitp,*) ' ERROR in QML_Read_RMat'
    write(out_unitp,*) '  while reading a matrix'
    write(out_unitp,*) '  end of file or end of record'
    write(out_unitp,*) '  The matrix paramters: nl,nc,nbcol',nl,nc,nbcol
    write(out_unitp,*) ' Check your data !!'
  END IF

!$OMP    END CRITICAL (QML_Read_RMat_CRIT)

  END SUBROUTINE QML_Read_RMat

  SUBROUTINE QML_Read_RVec(l,nio,nbcol,err)

  integer, intent(in)                :: nio,nbcol
  real(kind=Rkind), intent(inout)    :: l(:)
  integer, intent(inout)             :: err

  integer :: n,i,nb,nbblocs,nfin

!$OMP    CRITICAL (QML_Read_RVec_CRIT)

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
    write(out_unitp,*) ' ERROR in QML_Read_RVec'
    write(out_unitp,*) '  while reading a vector'
    write(out_unitp,*) '  end of file or end of record'
    write(out_unitp,*) '  The vector paramters: n,nbcol',n,nbcol
    write(out_unitp,*) ' Check your data !!'
  END IF

!$OMP    END CRITICAL (QML_Read_RVec_CRIT)

  END SUBROUTINE QML_Read_RVec

  SUBROUTINE QML_Init_IdMat(Mat,ndim)
  IMPLICIT NONE

  integer,                        intent(in)    :: ndim
  real (kind=Rkind), allocatable, intent(inout) :: mat(:,:)

  integer :: i

    IF (allocated(mat)) deallocate(mat)

    allocate(mat(ndim,ndim))
    mat(:,:) = ZERO
    DO i=1,ndim
      mat(i,i) = ONE
    END DO

  END SUBROUTINE QML_Init_IdMat


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

 51  format(a,i10,a,a)


     flush(out_unitp)
     !============================================

     count_old = count
     t_cpu_old = t_cpu

!$OMP   END CRITICAL (QML_time_perso_CRIT)


  END SUBROUTINE QML_time_perso

  FUNCTION QML_err_file_name(file_name,name_sub)
  integer                                 :: QML_err_file_name
  character (len=*), intent(in)           :: file_name
  character (len=*), intent(in), optional :: name_sub

    IF (string_IS_empty(file_name) ) THEN
      IF (present(name_sub)) THEN
        write(out_unitp,*) ' ERROR in QML_err_file_name, called from: ',name_sub
      ELSE
        write(out_unitp,*) ' ERROR in QML_err_file_name'
      END IF
      write(out_unitp,*) '   The file name is empty'
      QML_err_file_name = 1
    ELSE
      QML_err_file_name = 0
    END IF

  END FUNCTION QML_err_file_name

  FUNCTION QML_make_FileName(FileName)
    character(len=*), intent(in)    :: FileName

    character (len=:), allocatable  :: QML_make_FileName
    integer :: ilast_char

    ilast_char = len_trim(File_path)

    IF (FileName(1:1) == "/" .OR. FileName(1:1) == "" .OR. ilast_char == 0) THEN
      QML_make_FileName = trim(adjustl(FileName))
    ELSE
      IF (File_path(ilast_char:ilast_char) == "/") THEN
        QML_make_FileName = trim(adjustl(File_path)) // trim(adjustl(FileName))
      ELSE
        QML_make_FileName = trim(adjustl(File_path)) // '/' // trim(adjustl(FileName))
      END IF
    END IF
    !write(666,*) 'QML_make_FileName: ',QML_make_FileName
    !stop
  END FUNCTION QML_make_FileName
  SUBROUTINE QML_file_open2(name_file,iunit,lformatted,append,old,err_file)

  character (len=*),   intent(in)              :: name_file
  integer,             intent(inout)           :: iunit
  logical,             intent(in),    optional :: lformatted,append,old
  integer,             intent(out),   optional :: err_file

  character (len=20)   :: fform,fstatus,fposition
  logical              :: unit_opened
  integer              :: err_file_loc

! - default for the open ---------------------------

! - test if optional arguments are present ---------
  IF (present(lformatted)) THEN
    IF (.NOT. lformatted) THEN
      fform = 'unformatted'
    ELSE
      fform = 'formatted'
    END IF
  ELSE
    fform = 'formatted'
  END IF

  IF (present(append)) THEN
    IF (append) THEN
      fposition = 'append'
    ELSE
      fposition = 'asis'
    END IF
  ELSE
    fposition = 'asis'
  END IF

  IF (present(old)) THEN
    IF (old) THEN
      fstatus = 'old'
    ELSE
      fstatus = 'unknown'
    END IF
  ELSE
      fstatus = 'unknown'
  END IF

  err_file_loc = err_file_name(name_file,'QML_file_open2')
  IF (.NOT. present(err_file) .AND. err_file_loc /= 0) STOP ' ERROR, the file name is empty!'
  IF (present(err_file)) err_file = err_file_loc

! - check if the file is already open ------------------
! write(out_unitp,*) 'name,unit,unit_opened ',name_file,unit,unit_opened

  inquire(FILE=name_file,NUMBER=iunit,OPENED=unit_opened)
! write(out_unitp,*) 'name,unit,unit_opened ',name_file,unit,unit_opened


! - the file is not open, find an unused UNIT ---------
  IF (unit_opened) RETURN ! the file is already open

  iunit = 9
  DO
    iunit = iunit + 1
    inquire(UNIT=iunit,OPENED=unit_opened)
!   write(out_unitp,*) 'name,iunit,unit_opened ',name_file,iunit,unit_opened
    IF (.NOT. unit_opened) exit
  END DO


! -- open the file
  IF (present(err_file)) THEN
    open(UNIT=iunit,FILE=name_file,FORM=fform,STATUS=fstatus,       &
         POSITION=fposition,ACCESS='SEQUENTIAL',IOSTAT=err_file)
  ELSE
    open(UNIT=iunit,FILE=name_file,FORM=fform,STATUS=fstatus,       &
         POSITION=fposition,ACCESS='SEQUENTIAL')
  END IF

! write(out_unitp,*) 'open ',name_file,iunit

  END SUBROUTINE QML_file_open2
END MODULE QMLLib_UtilLib_m
