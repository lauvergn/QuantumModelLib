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

  PUBLIC :: Find_Label,NFind_Label,FFind_Label,Pot_Name_Analysis

  
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
      write(out_unitp,*) 'BEGINNING QML_Pot_Name_Analysis'
      write(out_unitp,*) 'pot_name: ',pot_name
      flush(out_unitp)
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
      IF (debug) write(out_unitp,*) 'nb_word,len_word',nb_word,len_word,' ',name
      flush(out_unitp)
    END DO
    IF (debug) write(out_unitp,*) 'nb_word,len_word',nb_word,len_word

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
        IF (debug) write(out_unitp,*) 'nb_word',nb_word,' ',tab_pot_name(nb_word)
        EXIT
      END IF

      nb_word = nb_word + 1
      len_word = max(len_word,iblank-1)

      tab_pot_name(nb_word) = trim(adjustl(name(1:iblank)))
      IF (debug)write(out_unitp,*) 'nb_word',nb_word,' ',tab_pot_name(nb_word)

      name = name(iblank+1:len(name))
      name = trim(adjustl(name))
      flush(out_unitp)
    END DO
    IF (debug) THEN
      write(out_unitp,*) 'END QML_Pot_Name_Analysis'
      flush(out_unitp)
    END IF
  END SUBROUTINE QML_Pot_Name_Analysis
END MODULE QMLLib_UtilLib_m
