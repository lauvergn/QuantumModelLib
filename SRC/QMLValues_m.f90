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
  PUBLIC :: QMLValues_t, alloc_QMLValues, dealloc_QMLValues, Write_QMLValues
  PUBLIC :: WxQMLValuesd0_ADDTO_QMLValues2_ider

  TYPE :: QMLValues_t
    logical                        :: alloc     = .FALSE.
    logical                        :: adiabatic = .FALSE.
    integer                        :: nderiv    = -1
    integer                        :: ndim      = -1
    integer                        :: nsurf     = -1
    integer                        :: nb_ScalOp = -1

    real (kind=Rkind), allocatable :: Q(:)

    TYPE (dnMat_t)                 :: PotAdia
    TYPE (dnMat_t)                 :: PotDia
    TYPE (dnMat_t)                 :: NAC
    TYPE (dnMat_t)                 :: Vec
    TYPE (dnMat_t)                 :: Vec0

    TYPE (dnMat_t)                 :: ImagPotAdia
    TYPE (dnMat_t)                 :: ImagPotDia
    TYPE (dnMat_t)                 :: ImagNAC
    TYPE (dnMat_t)                 :: ImagVec
    TYPE (dnMat_t)                 :: ImagVec0

    TYPE (dnMat_t), allocatable    :: ScalOpAdia(:)
    TYPE (dnMat_t), allocatable    :: ScalOpDia(:)
    TYPE (dnMat_t), allocatable    :: ImagScalOpAdia(:)
    TYPE (dnMat_t), allocatable    :: ImagScalOpDia(:)

  END TYPE QMLValues_t

CONTAINS

  SUBROUTINE alloc_QMLValues(QMLValues,adiabatic,cplx,ndim,nsurf,nb_ScalOp,nderiv)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(QMLValues_t),      intent(inout)           :: QMLValues
    logical,                intent(in)              :: adiabatic,cplx
    integer,                intent(in)              :: ndim,nsurf,nb_ScalOp,nderiv

    integer :: i

    CALL dealloc_QMLValues(QMLValues)

    QMLValues%nderiv       = nderiv
    QMLValues%nsurf        = nsurf
    QMLValues%ndim         = ndim
    QMLValues%nb_ScalOp    = nb_ScalOp

    QMLValues%adiabatic    = adiabatic

    allocate(QMLValues%Q(ndim))

    CALL alloc_dnMat(QMLValues%PotDia,nsurf=nsurf,nVar=ndim,nderiv=nderiv)
    QMLValues%PotDia = ZERO

    IF (cplx) THEN 
      CALL alloc_dnMat(QMLValues%ImagPotDia,nsurf=nsurf,nVar=ndim,nderiv=nderiv)
      QMLValues%ImagPotDia = ZERO
    END IF

    IF (nb_ScalOp >0) THEN
      allocate(QMLValues%ScalOpDia(nb_ScalOp))
      DO i=1,nb_ScalOp
        CALL alloc_dnMat(QMLValues%ScalOpDia(i),nsurf=nsurf,nVar=ndim,nderiv=nderiv)
        QMLValues%ScalOpDia(i) = ZERO
      END DO
      IF (cplx) THEN
        allocate(QMLValues%ImagScalOpDia(nb_ScalOp))
        DO i=1,nb_ScalOp
          CALL alloc_dnMat(QMLValues%ImagScalOpDia(i),nsurf=nsurf,nVar=ndim,nderiv=nderiv)
          QMLValues%ImagScalOpDia(i) = ZERO
        END DO
      END IF
    END IF

    IF (adiabatic) THEN
      CALL alloc_dnMat(QMLValues%PotAdia,nsurf=nsurf,nVar=ndim,nderiv=nderiv)
      QMLValues%PotAdia = ZERO
      CALL alloc_dnMat(QMLValues%Vec,nsurf=nsurf,nVar=ndim,nderiv=nderiv)
      CALL alloc_dnMat(QMLValues%NAC,nsurf=nsurf,nVar=ndim,nderiv=nderiv)
      CALL alloc_dnMat(QMLValues%Vec0,nsurf=nsurf,nVar=ndim,nderiv=0)
      QMLValues%Vec  = ZERO
      QMLValues%Vec0 = ZERO
      QMLValues%NAC  = ZERO

      IF (cplx) THEN 
        CALL alloc_dnMat(QMLValues%ImagPotAdia,nsurf=nsurf,nVar=ndim,nderiv=nderiv)
        QMLValues%ImagPotAdia = ZERO
        CALL alloc_dnMat(QMLValues%ImagVec,nsurf=nsurf,nVar=ndim,nderiv=nderiv)
        CALL alloc_dnMat(QMLValues%ImagNAC,nsurf=nsurf,nVar=ndim,nderiv=nderiv)
        CALL alloc_dnMat(QMLValues%ImagVec0,nsurf=nsurf,nVar=ndim,nderiv=0)
        QMLValues%ImagVec  = ZERO
        QMLValues%ImagVec0 = ZERO
        QMLValues%ImagNAC  = ZERO
      END IF

      IF (nb_ScalOp >0) THEN
        allocate(QMLValues%ScalOpAdia(nb_ScalOp))
        DO i=1,nb_ScalOp
          CALL alloc_dnMat(QMLValues%ScalOpAdia(i),nsurf=nsurf,nVar=ndim,nderiv=nderiv)
          QMLValues%ScalOpAdia(i) = ZERO
        END DO
        IF (cplx) THEN
          allocate(QMLValues%ImagScalOpAdia(nb_ScalOp))
          DO i=1,nb_ScalOp
            CALL alloc_dnMat(QMLValues%ImagScalOpAdia(i),nsurf=nsurf,nVar=ndim,nderiv=nderiv)
            QMLValues%ImagScalOpAdia(i) = ZERO
          END DO
        END IF
      END IF
    END IF

    QMLValues%alloc    = .TRUE.

  END SUBROUTINE alloc_QMLValues
  SUBROUTINE dealloc_QMLValues(QMLValues)
    IMPLICIT NONE

    TYPE(QMLValues_t),      intent(inout)           :: QMLValues

    integer :: i

    QMLValues%alloc        = .FALSE.
    IF (allocated(QMLValues%Q)) deallocate(QMLValues%Q)
    QMLValues%nderiv       = -1
    QMLValues%nsurf        = -1
    QMLValues%ndim         = -1
    QMLValues%nb_ScalOp    = -1
    QMLValues%adiabatic    = .FALSE.

    CALL dealloc_dnMat(QMLValues%PotAdia)
    CALL dealloc_dnMat(QMLValues%PotDia)

    CALL dealloc_dnMat(QMLValues%Vec)
    CALL dealloc_dnMat(QMLValues%Vec0)
    CALL dealloc_dnMat(QMLValues%NAC)

    CALL dealloc_dnMat(QMLValues%ImagPotAdia)
    CALL dealloc_dnMat(QMLValues%ImagPotDia) 
    CALL dealloc_dnMat(QMLValues%ImagVec)
    CALL dealloc_dnMat(QMLValues%ImagVec0)
    CALL dealloc_dnMat(QMLValues%ImagNAC)

    IF (allocated(QMLValues%ScalOpDia)) THEN
      DO i=1,size(QMLValues%ScalOpDia)
        CALL dealloc_dnMat(QMLValues%ScalOpDia(i))
      END DO
      deallocate(QMLValues%ScalOpDia)
    END IF
    IF (allocated(QMLValues%ImagScalOpDia)) THEN
      DO i=1,size(QMLValues%ImagScalOpDia)
        CALL dealloc_dnMat(QMLValues%ImagScalOpDia(i))
      END DO
      deallocate(QMLValues%ImagScalOpDia)
    END IF

    IF (allocated(QMLValues%ScalOpAdia)) THEN
      DO i=1,size(QMLValues%ScalOpAdia)
        CALL dealloc_dnMat(QMLValues%ScalOpDia(i))
      END DO
      deallocate(QMLValues%ScalOpAdia)
    END IF
    IF (allocated(QMLValues%ImagScalOpAdia)) THEN
      DO i=1,size(QMLValues%ImagScalOpDia)
        CALL dealloc_dnMat(QMLValues%ImagScalOpDia(i))
      END DO
      deallocate(QMLValues%ImagScalOpAdia)
    END IF
  END SUBROUTINE dealloc_QMLValues

  SUBROUTINE Write_QMLValues(QMLValues,nio,Flatten)
    USE QDUtil_m
    IMPLICIT NONE

    TYPE(QMLValues_t),  intent(in)           :: QMLValues
    integer,            intent(in), optional :: nio
    logical,            intent(in), optional :: Flatten


    integer :: i,nio_loc
    logical :: Flatten_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unit
    END IF
    IF (present(Flatten)) THEN
      Flatten_loc = Flatten
    ELSE
      Flatten_loc = .FALSE.
    END IF

    write(nio_loc,*) '==============================================='
    write(nio_loc,*) 'QML Values'
    write(nio_loc,*) 'alloc     : ',QMLValues%alloc
    write(nio_loc,*) 'nderiv    : ',QMLValues%nderiv
    write(nio_loc,*) 'ndim      : ',QMLValues%ndim
    write(nio_loc,*) 'nsurf     : ',QMLValues%nsurf
    write(nio_loc,*) 'nb_ScalOp : ',QMLValues%nb_ScalOp
    write(nio_loc,*) 'adiabatic : ',QMLValues%adiabatic

    IF (allocated(QMLValues%Q)) THEN
      write(nio_loc,*) 'Q (coordinates): ',QMLValues%Q
    ELSE
      write(nio_loc,*) 'Q (coordinates): not defined!'
    END IF
    write(nio_loc,*)

    IF (Flatten_loc) THEN
      write(nio_loc,*) 'PotAdia:    ',get_Flatten(QMLValues%PotAdia)
      write(nio_loc,*) 'PotDia:     ',get_Flatten(QMLValues%PotDia)
      write(nio_loc,*) 'Vec:        ',get_Flatten(QMLValues%Vec)
      write(nio_loc,*) 'Vec0:       ',get_Flatten(QMLValues%Vec0)
      write(nio_loc,*) 'NAC:        ',get_Flatten(QMLValues%NAC)
      write(nio_loc,*) 'ImagPotAdia:',get_Flatten(QMLValues%ImagPotAdia)
      write(nio_loc,*) 'ImagPotDia: ',get_Flatten(QMLValues%ImagPotDia)
      write(nio_loc,*) 'ImagVec:    ',get_Flatten(QMLValues%ImagVec)
      write(nio_loc,*) 'ImagVec0:   ',get_Flatten(QMLValues%ImagVec0)
      write(nio_loc,*) 'ImagNAC:    ',get_Flatten(QMLValues%ImagNAC)

      IF (allocated(QMLValues%ScalOpAdia)) THEN
        DO i=1,size(QMLValues%ScalOpAdia)
          write(nio_loc,*) 'ScalOpAdia:     ',i,get_Flatten(QMLValues%ScalOpAdia(i))
        END DO
      END IF
      IF (allocated(QMLValues%ImagScalOpAdia)) THEN
        DO i=1,size(QMLValues%ImagScalOpAdia)
          write(nio_loc,*) 'ImagScalOpAdia: ',i,get_Flatten(QMLValues%ImagScalOpAdia(i))
        END DO
      END IF
      IF (allocated(QMLValues%ScalOpDia)) THEN
        DO i=1,size(QMLValues%ScalOpDia)
          write(nio_loc,*) 'ScalOpDia:      ',i,get_Flatten(QMLValues%ScalOpDia(i))
        END DO
      END IF
      IF (allocated(QMLValues%ImagScalOpDia)) THEN
       DO i=1,size(QMLValues%ImagScalOpDia)
          write(nio_loc,*) 'ImagScalOpDia:  ',i,get_Flatten(QMLValues%ImagScalOpDia(i))
        END DO
      END IF

    ELSE
      write(nio_loc,*) '-----------------------------------------------'
      CALL Write_dnMat(QMLValues%PotAdia,     nio_loc,info='PotAdia')
      CALL Write_dnMat(QMLValues%ImagPotAdia, nio_loc,info='ImagPotAdia')
      write(nio_loc,*)
      CALL Write_dnMat(QMLValues%PotDia,      nio_loc,info='PotDia')
      CALL Write_dnMat(QMLValues%ImagPotDia,  nio_loc,info='ImagPotDia') 
      write(nio_loc,*)

      write(nio_loc,*) '-----------------------------------------------'
      CALL Write_dnMat(QMLValues%Vec,         nio_loc,info='Vec')
      CALL Write_dnMat(QMLValues%ImagVec,     nio_loc,info='ImagVec')
      write(nio_loc,*)
      CALL Write_dnMat(QMLValues%Vec0,        nio_loc,info='Vec0')
      CALL Write_dnMat(QMLValues%ImagVec0,    nio_loc,info='ImagVec0')
      write(nio_loc,*)
      CALL Write_dnMat(QMLValues%NAC,         nio_loc,info='NAC')
      CALL Write_dnMat(QMLValues%ImagNAC,     nio_loc,info='ImagNAC')  
      write(nio_loc,*)

      IF (allocated(QMLValues%ScalOpAdia)) THEN
        write(nio_loc,*) '-----------------------------------------------'
        DO i=1,size(QMLValues%ScalOpAdia)
          CALL Write_dnMat(QMLValues%ScalOpAdia(i),nio_loc,info='ScalOpAdia('//TO_string(i)//')')
        END DO
      END IF
      IF (allocated(QMLValues%ImagScalOpAdia)) THEN
        write(nio_loc,*) '-----------------------------------------------'
        DO i=1,size(QMLValues%ImagScalOpDia)
          CALL Write_dnMat(QMLValues%ImagScalOpDia(i),nio_loc,info='ImagScalOpDia('//TO_string(i)//')')
        END DO
      END IF
      IF (allocated(QMLValues%ScalOpDia)) THEN
        write(nio_loc,*) '-----------------------------------------------'
        DO i=1,size(QMLValues%ScalOpDia)
          CALL Write_dnMat(QMLValues%ScalOpDia(i),nio_loc,info='ScalOpDia('//TO_string(i)//')')
        END DO
      END IF
      IF (allocated(QMLValues%ImagScalOpDia)) THEN
        write(nio_loc,*) '-----------------------------------------------'
        DO i=1,size(QMLValues%ImagScalOpDia)
          CALL Write_dnMat(QMLValues%ImagScalOpDia(i),nio_loc,info='ImagScalOpDia('//TO_string(i)//')')
        END DO
      END IF
      write(nio_loc,*) '==============================================='
    END IF
    flush(nio_loc)

  END SUBROUTINE Write_QMLValues

  SUBROUTINE WxQMLValuesd0_ADDTO_QMLValues2_ider(QMLValues,W,QMLValues2,ider)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE(QMLValues_t),  intent(in)            :: QMLValues
    TYPE(QMLValues_t),  intent(inout)         :: QMLValues2
    integer,            intent(in),  optional :: ider(:)
    real (kind=Rkind),  intent(in)            :: W

    integer :: i
    character (len=*), parameter :: name_sub='WxQMLValuesd0_ADDTO_QMLValues2_ider'

    IF (present(ider)) THEN
      CALL Mat_wADDTO_dnMat2_ider(QMLValues%PotDia%d0,W,QMLValues2%PotDia,ider)

      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagPotDia,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagPotDia%d0,W,QMLValues2%ImagPotDia,ider)

      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%PotAdia,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%PotAdia%d0,W,QMLValues2%PotAdia,ider)
      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagPotAdia,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagPotAdia%d0,W,QMLValues2%ImagPotAdia,ider)

      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%Vec,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%Vec%d0,W,QMLValues2%Vec,ider)
      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagVec,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagVec%d0,W,QMLValues2%ImagVec,ider)

      IF (allocated(QMLValues%ScalOpDia)) THEN
        DO i=1,size(QMLValues%ScalOpDia)
          IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ScalOpDia(i),QMLValues%nderiv))  &
            CALL Mat_wADDTO_dnMat2_ider(QMLValues%ScalOpDia(i)%d0,W,QMLValues2%ScalOpDia(i),ider)
        END DO
      END IF
      IF (allocated(QMLValues%ImagScalOpDia)) THEN
        DO i=1,size(QMLValues%ImagScalOpDia)
          IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagScalOpDia(i),QMLValues%nderiv))  &
            CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagScalOpDia(i)%d0,W,QMLValues2%ImagScalOpDia(i),ider)
        END DO
      END IF
      IF (allocated(QMLValues%ScalOpAdia)) THEN
        DO i=1,size(QMLValues%ScalOpAdia)
          IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ScalOpAdia(i),QMLValues%nderiv))  &
            CALL Mat_wADDTO_dnMat2_ider(QMLValues%ScalOpAdia(i)%d0,W,QMLValues2%ScalOpAdia(i),ider)
        END DO
      END IF
      IF (allocated(QMLValues%ImagScalOpAdia)) THEN
        DO i=1,size(QMLValues%ImagScalOpAdia)
          IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagScalOpAdia(i),QMLValues%nderiv))  &
            CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagScalOpAdia(i)%d0,W,QMLValues2%ImagScalOpAdia(i),ider)
        END DO
      END IF
    ELSE
      CALL Mat_wADDTO_dnMat2_ider(QMLValues%PotDia%d0,W,QMLValues2%PotDia)

      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagPotDia,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagPotDia%d0,W,QMLValues2%ImagPotDia)

      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%PotAdia,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%PotAdia%d0,W,QMLValues2%PotAdia)
      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagPotAdia,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagPotAdia%d0,W,QMLValues2%ImagPotAdia)

      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%Vec,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%Vec%d0,W,QMLValues2%Vec)
      IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagVec,QMLValues%nderiv))  &
        CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagVec%d0,W,QMLValues2%ImagVec) 


      IF (allocated(QMLValues%ScalOpDia)) THEN
        DO i=1,size(QMLValues%ScalOpDia)
          IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ScalOpDia(i),QMLValues%nderiv))  &
            CALL Mat_wADDTO_dnMat2_ider(QMLValues%ScalOpDia(i)%d0,W,QMLValues2%ScalOpDia(i))
        END DO
      END IF
      IF (allocated(QMLValues%ImagScalOpDia)) THEN
        DO i=1,size(QMLValues%ImagScalOpDia)
          IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagScalOpDia(i),QMLValues%nderiv))  &
            CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagScalOpDia(i)%d0,W,QMLValues2%ImagScalOpDia(i))
        END DO
      END IF
      IF (allocated(QMLValues%ScalOpAdia)) THEN
        DO i=1,size(QMLValues%ScalOpAdia)
          IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ScalOpAdia(i),QMLValues%nderiv))  &
            CALL Mat_wADDTO_dnMat2_ider(QMLValues%ScalOpAdia(i)%d0,W,QMLValues2%ScalOpAdia(i))
        END DO
      END IF
      IF (allocated(QMLValues%ImagScalOpAdia)) THEN
        DO i=1,size(QMLValues%ImagScalOpAdia)
          IF (.NOT. Check_NotAlloc_dnMat(QMLValues%ImagScalOpAdia(i),QMLValues%nderiv))  &
            CALL Mat_wADDTO_dnMat2_ider(QMLValues%ImagScalOpAdia(i)%d0,W,QMLValues2%ImagScalOpAdia(i))
        END DO
      END IF
    END IF

  END SUBROUTINE WxQMLValuesd0_ADDTO_QMLValues2_ider

END MODULE QMLValues_m
