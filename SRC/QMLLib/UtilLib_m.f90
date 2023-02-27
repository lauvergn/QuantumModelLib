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
  USE QDUtil_NumParameters_m
  !$ USE omp_lib
  IMPLICIT NONE

  PRIVATE

  character (len=*), parameter :: QML_path   =                         &
#if defined(__QMLPATH)
    __QMLPATH
#else
    '~/QuantumModelLib'
#endif

  PUBLIC :: Find_Label,NFind_Label,FFind_Label,Pot_Name_Analysis,      &
            Read_alloc_Vect,make_QMLInternalFileName,QML_path

  
  INTERFACE Pot_Name_Analysis
    MODULE PROCEDURE QML_Pot_Name_Analysis
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
  INTERFACE Read_alloc_Vect
    MODULE PROCEDURE QML_Read_alloc_Vect
  END INTERFACE

  CONTAINS

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

        !write(out_unit,*) 'format_label',format_label
  !     write(out_unit,*) 'label to find:',label,located
        DO
          i_line = i_line + 1

          read(nio,format_label,end=999,eor=888,err=888,advance='no') labelR

          located = (labelR .EQ. label)
          !located = verify(label,labelR) == 0
          !write(out_unit,*) i_line,located,labelR
          IF (located) EXIT
          read(nio,*,end=999,err=888)

   888    CONTINUE
        END DO


   999  CONTINUE
        i_line = 0

        !write(out_unit,*) 'Find_Label: ',label,located
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
          write(out_unit,*) ' ERROR in NFind_Label'
          write(out_unit,*) ' iformat < 1',iformat
          STOP
        END IF
        fformat='(A' // TO_string(iformat) // ')'
        !write(out_unit,*) iformat,fformat
        !write(out_unit,*) 'label to find:',label,located
        DO
          i_line = i_line + 1
          read(nio,fformat,end=999,eor=888,err=888,advance='no') labelR

          located = verify(label,labelR) == 0
          !write(out_unit,*) i_line,located,labelR
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

        !write(out_unit,*) 'label to find:',label,located
        DO
          i_line = i_line + 1
          read(nio,fformat,end=999,eor=888,err=888,advance='no') labelR


          located = verify(label,labelR) == 0
          !write(out_unit,*) i_line,located,labelR
          IF (located) EXIT
          read(nio,*,end=999,err=888)

   888    CONTINUE
        END DO


   999  CONTINUE
        i_line = 0
        IF (located) read(nio,*)

  END SUBROUTINE QML_FFind_Label

  SUBROUTINE QML_Pot_Name_Analysis(pot_name,tab_pot_name)
    USE QDUtil_m
    IMPLICIT NONE

    character (len=*),              intent(in)    :: pot_name
    character (len=:), allocatable, intent(inout) :: tab_pot_name(:)

    character (len=:), allocatable :: name
    integer :: nb_word,len_word,ioerr,i,iblank

    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING QML_Pot_Name_Analysis'
      write(out_unit,*) 'pot_name: ',pot_name
      flush(out_unit)
    END IF
    !first the number of word
    name     = pot_name ! for the allocation
    nb_word  = 0
    len_word = 0
    DO
      iblank = index(name,' ')
      IF (iblank == 0) THEN
        nb_word = nb_word + 1
        len_word = max(len_word,len_trim(name))
        EXIT
      END IF

      name = name(iblank+1:len(name))
      name = trim(adjustl(name))

      nb_word = nb_word + 1
      len_word = max(len_word,iblank-1)
      IF (debug) write(out_unit,*) 'nb_word,len_word',nb_word,len_word,' ',name
      flush(out_unit)
    END DO
    IF (debug) write(out_unit,*) 'nb_word,len_word',nb_word,len_word

    !now create the table of string
    allocate(character(len=len_word) :: tab_pot_name(nb_word))

    name     = pot_name ! for the allocation
    nb_word  = 0
    len_word = 0
    DO
      iblank = index(name,' ')
      IF (iblank == 0) THEN
        nb_word = nb_word + 1
        len_word = max(len_word,len_trim(name))
        tab_pot_name(nb_word) = trim(adjustl(name))
        IF (debug) write(out_unit,*) 'nb_word',nb_word,' ',tab_pot_name(nb_word)
        EXIT
      END IF

      nb_word = nb_word + 1
      len_word = max(len_word,iblank-1)

      tab_pot_name(nb_word) = trim(adjustl(name(1:iblank)))
      IF (debug)write(out_unit,*) 'nb_word',nb_word,' ',tab_pot_name(nb_word)

      name = name(iblank+1:len(name))
      name = trim(adjustl(name))
      flush(out_unit)
    END DO
    IF (debug) THEN
      write(out_unit,*) 'END QML_Pot_Name_Analysis'
      flush(out_unit)
    END IF
  END SUBROUTINE QML_Pot_Name_Analysis
  FUNCTION QML_Read_alloc_Vect(nio,err_io) RESULT (Vec)
    IMPLICIT NONE

    integer, intent(in)           :: nio
    integer, intent(inout)        :: err_io
    real(kind=Rkind), allocatable :: Vec(:) ! result

    integer :: nsize
    character (len=Line_len) :: name_dum

    read(nio,*,IOSTAT=err_io)

    read(nio,*,IOSTAT=err_io) nsize
    IF (err_io /= 0) RETURN
    allocate(Vec(nsize))
    read(nio,*,IOSTAT=err_io) Vec
    !write(out_unit,*) 'read ',nsize,' real'
    !flush(out_unit)

  END FUNCTION QML_Read_alloc_Vect

  FUNCTION make_QMLInternalFileName(FileName,FPath) RESULT(make_FileName)
    USE QDUtil_m, ONLY : err_FileName
    IMPLICIT NONE

    character (len=:), allocatable          :: make_FileName

    character(len=*), intent(in)            :: FileName
    character(len=*), intent(in), optional  :: FPath

    character (len=:), allocatable          :: FPath_loc


    integer :: ilast_char,err

    IF (present(FPath)) THEN
      FPath_loc = FPath
    ELSE
      FPath_loc = QML_path
    END IF

    err = err_FileName(FileName,name_sub='make_QMLInternalFileName')
    IF (err /= 0) STOP 'ERROR in make_QMLInternalFileName: problem with the FileName'

    ilast_char = len_trim(FPath_loc)

    IF (FileName(1:1) == "/" .OR. FileName(1:1) == "~" .OR. ilast_char == 0) THEN
      make_FileName = trim(adjustl(FileName))
    ELSE
      IF (FPath_loc(ilast_char:ilast_char) == "/") THEN
        make_FileName = trim(adjustl(FPath_loc)) // trim(adjustl(FileName))
      ELSE
        make_FileName = trim(adjustl(FPath_loc)) // '/' // trim(adjustl(FileName))
      END IF
    END IF

    IF (allocated(FPath_loc)) deallocate(FPath_loc)

  END FUNCTION make_QMLInternalFileName
END MODULE QMLLib_UtilLib_m
