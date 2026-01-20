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
MODULE QMLValues_m
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: QMLValues_t, dealloc_QMLValues, Write_QMLValues

  TYPE :: QMLValues_t
    logical                        :: adiabatic = .FALSE.
    integer                        :: nderiv    = -1

    real (kind=Rkind), allocatable :: Q(:)

    TYPE (dnMat_t)                 :: PotAdia
    TYPE (dnMat_t)                 :: ImagPotAdia
    TYPE (dnMat_t)                 :: PotDia
    TYPE (dnMat_t)                 :: ImagPotDia
 
    TYPE (dnMat_t)                 :: NAC
    TYPE (dnMat_t)                 :: Vec
    TYPE (dnMat_t)                 :: Vec0

  END TYPE QMLValues_t

CONTAINS

  SUBROUTINE dealloc_QMLValues(QMLValues)
    IMPLICIT NONE

    TYPE(QMLValues_t),      intent(inout)           :: QMLValues

    IF (allocated(QMLValues%Q)) deallocate(QMLValues%Q)
    QMLValues%nderiv       = -1
    QMLValues%adiabatic    = .FALSE.

    CALL dealloc_dnMat(QMLValues%PotAdia)
    CALL dealloc_dnMat(QMLValues%ImagPotAdia)
    CALL dealloc_dnMat(QMLValues%PotDia)
    CALL dealloc_dnMat(QMLValues%ImagPotDia) 

    CALL dealloc_dnMat(QMLValues%Vec)
    CALL dealloc_dnMat(QMLValues%Vec0)
    CALL dealloc_dnMat(QMLValues%NAC)

  END SUBROUTINE dealloc_QMLValues

  SUBROUTINE Write_QMLValues(QMLValues,nio)
    IMPLICIT NONE

    TYPE(QMLValues_t),  intent(in)           :: QMLValues
    integer,            intent(in), optional :: nio

    integer :: nio_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unit
    END IF

    write(nio_loc,*) '-----------------------------------------------'
    write(nio_loc,*) 'QML Values'

    write(nio_loc,*) 'nderiv   : ',QMLValues%nderiv
    write(nio_loc,*) 'adiabatic: ',QMLValues%adiabatic

    IF (allocated(QMLValues%Q)) THEN
      write(nio_loc,*) 'Q (coordinates): ',QMLValues%Q
    ELSE
      write(nio_loc,*) 'Q (coordinates): not defined!'
    END IF
    
    CALL Write_dnMat(QMLValues%PotAdia,     nio_loc,info='PotAdia')
    CALL Write_dnMat(QMLValues%ImagPotAdia, nio_loc,info='ImagPotAdia')
    CALL Write_dnMat(QMLValues%PotDia,      nio_loc,info='PotDia')
    CALL Write_dnMat(QMLValues%ImagPotDia,  nio_loc,info='ImagPotDia') 

    CALL Write_dnMat(QMLValues%Vec,         nio_loc,info='Vec')
    CALL Write_dnMat(QMLValues%Vec0,        nio_loc,info='Vec0')
    CALL Write_dnMat(QMLValues%NAC,         nio_loc,info='NAC') 
  
    write(nio_loc,*) '-----------------------------------------------'
    flush(nio_loc)

  END SUBROUTINE Write_QMLValues

END MODULE QMLValues_m
