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
PROGRAM Compare_txt
  IMPLICIT NONE

  character(len=:), allocatable :: arg,arg2
  integer :: long , i,n

  character (len=:), allocatable  :: file_name_old,file_name_new
  character (len=20)               :: name_old,name_new

  integer :: itest,nb_error,nio_old=10,nio_new=11,err_file_old,err_file_new,err_sub



  real (kind=8), allocatable :: tab_old(:),tab_new(:)
  real (kind=8)              :: max_diff,max_diff_test



  DO i=1,COMMAND_ARGUMENT_COUNT(),2
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, LENGTH=long )
    allocate( character(len=long) :: arg )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i, VALUE=arg )

    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, LENGTH=long )
    allocate( character(len=long) :: arg2 )
    CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, VALUE=arg2 )

    SELECT CASE(arg)
    CASE("-o","-old")
      file_name_old = arg2
    CASE("-n","-new")
      file_name_new = arg2
    CASE Default
      STOP 'no default argument'
    END SELECT

    deallocate(arg)
    deallocate(arg2)

  END DO

  write(6,*) '---------------------------------------------'
  write(6,*) '   Differences between:'
  write(6,*) 'file_name_old: "',file_name_old,'"'
  write(6,*) 'file_name_new: "',file_name_new,'"'
  write(6,*) '---------------------------------------------'

  open(unit=nio_old,file=file_name_old,STATUS='old',IOSTAT=err_file_old)
  open(unit=nio_new,file=file_name_new,STATUS='old',IOSTAT=err_file_new)

  IF (err_file_old /= 0) THEN
    write(6,*) ' WARNING while openning the old file. err_file_old: ',err_file_old
    write(6,*) 'file_name_old: ',file_name_old
  END IF
  IF (err_file_new /= 0) THEN
    write(6,*) ' WARNING while openning the new file. err_file_new: ',err_file_new
    write(6,*) 'file_name_new: ',file_name_new
  END IF

  IF (err_file_old /= 0 .AND. err_file_new == 0 .OR. err_file_old == 0 .AND. err_file_new /= 0) THEN
    write(6,*) ' Only ONE file can be open!!'
    write(6,*) ' ERROR'
    STOP
  END IF
  IF (err_file_old /= 0 .AND. err_file_new /= 0) THEN
    write(6,*) ' No file can be open!!'
    write(6,*) ' NO ERROR'
    STOP
  END IF

  itest    = 0
  nb_error = 0
  DO
    itest = itest + 1
    max_diff_test = 0.d0

    read(nio_old,*,IOSTAT=err_file_old) name_old
    read(nio_new,*,IOSTAT=err_file_new) name_new
    CALL compa_name(name_old,name_new,err_file_old,err_file_new,err_sub)
    IF (err_sub /= 0 .OR. err_file_old /= 0) EXIT ! 'END OF FILE or ERROR'

    IF (name_old == 'TEST') THEN
      DO
        read(nio_old,*,IOSTAT=err_file_old) name_old
        read(nio_new,*,IOSTAT=err_file_new) name_new
        !write(6,*) 'name_old ',name_old
        CALL compa_name(name_old,name_new,err_file_old,err_file_new,err_sub)
        IF (err_sub /= 0 .OR. name_old == 'END_TEST') EXIT

        read(nio_old,*) n
        read(nio_new,*) n

        allocate(tab_old(n))
        allocate(tab_new(n))
        read(nio_old,*) tab_old
        read(nio_new,*) tab_new
        max_diff = maxval(abs(tab_old-tab_new))
        IF (max_diff > 1.d-6) THEN
          write(6,*) name_old,n,'max_diff',max_diff,' large difference'
          write(6,*) 'OLD values',n
          write(6,*) tab_old
          write(6,*) 'NEW values',n
          write(6,*) tab_new
        ELSE
          write(6,*) name_old,n,'max_diff',max_diff
        END IF
        deallocate(tab_old)
        deallocate(tab_new)

        max_diff_test = max(max_diff_test,max_diff)

      END DO
      IF (max_diff_test > 1.d-6) THEN
        nb_error = nb_error + 1
        write(6,'(a,i0,a,a,a,e9.3,a)') 'TEST: ',itest,', file_name_new: "',file_name_new,'", max_diff: ',max_diff_test,' ERROR'
      ELSE
        write(6,'(a,i0,a,a,a,e9.3)') 'TEST: ',itest,', file_name_new: "',file_name_new,'", max_diff: ',max_diff_test
      END IF
    END IF
  END DO

  IF (nb_error > 0) THEN
    write(6,'(a,a,a,i0,a)') 'In file_name_new: "',file_name_new,'", ',nb_error,' PROBLEM(S)'
  ELSE
    write(6,'(a,a,a)') 'In file_name_new: "',file_name_new,'", NO PROBLEM'
  END IF
  write(6,*) '---------------------------------------------'

  close(nio_old)
  close(nio_new)

END PROGRAM Compare_txt
SUBROUTINE compa_name(name_old,name_new,err_file_old,err_file_new,err_sub)
  IMPLICIT NONE

  character (len=*)               :: name_old,name_new
  integer                         :: err_file_old,err_file_new,err_sub

  err_sub = 0
  IF (name_old /= name_new .OR. err_file_old /= err_file_new) THEN
    write(6,*) ' The structure of both files are different!'
    write(6,*) 'name_old,err_file_old: ',name_old,err_file_old
    write(6,*) 'name_new,err_file_new: ',name_new,err_file_new
    write(6,*) ' ERROR'
    err_sub = 1
  END IF

END SUBROUTINE compa_name
