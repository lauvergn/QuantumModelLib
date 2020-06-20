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
!> @brief Module which deals with derivatives of a matrix (as function of coordinates).
!!
!! This module deals with operations or functions of a matrix function and its derivatives, dnMat.
!!
!! There is a mapping between the matrix function M, its derivatives and the dnM derived type components:
!!
!! @li M(:,:)                   => M%d0(:,:)
!! @li dM(:,:)/dQ_i             => M%d1(:,:,i)
!! @li d^2M(:,:)/dQ_idQ_j       => M%d2(:,:,i,j)
!! @li d^3M(:,:)/dQ_idQ_jdQ_k   => M%d3(:,:,i,j,k)
!!
!! with M defined as:
!!  TYPE (dnMat_t) :: M
!!
!!
!! Some standard fortran operators (= + - * **) are overloaded (!!! not /):
!!
!! For instance the sum (+) of two dnMat variables, M1 and M2 correspond to:
!! @li (M1+M2)                 => M1%d0    + M2%d0
!! @li d(M1+M2)/dQ_i           => M1%d1(:,:,i) + M2%d1(:,:,i)
!! @li ....
!!
!!
!! @author David Lauvergnat
!! @date 09/08/2017
!!
MODULE mod_dnMat
  USE mod_NumParameters
  IMPLICIT NONE

  TYPE dnMat_t
     integer                        :: nderiv = -1

     real (kind=Rkind), allocatable :: d0(:,:)
     real (kind=Rkind), allocatable :: d1(:,:,:)
     real (kind=Rkind), allocatable :: d2(:,:,:,:)
     real (kind=Rkind), allocatable :: d3(:,:,:,:,:)

  CONTAINS
    PROCEDURE, PRIVATE :: QML_sub_dnMat2_TO_dnMat1
    PROCEDURE, PRIVATE :: QML_set_dnMat_TO_R
    PROCEDURE, PRIVATE :: QML_set_dnMat_FROM_MatOFdnS
    GENERIC,   PUBLIC  :: assignment(=) => QML_sub_dnMat2_TO_dnMat1,            &
                                           QML_set_dnMat_TO_R,                  &
                                           QML_set_dnMat_FROM_MatOFdnS
  END TYPE dnMat_t

  PRIVATE :: QML_sub_dnMat2_TO_dnMat1,QML_set_dnMat_TO_R,QML_set_dnMat_FROM_MatOFdnS
  PRIVATE :: QML_TRANSPOSE_dnMat,QML_MATMUL_dnMat1_dnMat2,              &
             QML_MATMUL_dnMat1_Mat2,QML_MATMUL_Mat1_dnMat2

   INTERFACE transpose
      MODULE PROCEDURE QML_TRANSPOSE_dnMat
   END INTERFACE
   INTERFACE matmul
      MODULE PROCEDURE QML_MATMUL_dnMat1_dnMat2,QML_MATMUL_dnMat1_Mat2, &
                       QML_MATMUL_Mat1_dnMat2
   END INTERFACE

   INTERFACE operator (*)
      MODULE PROCEDURE QML_sub_dnMat_TIME_R,QML_sub_R_TIME_dnMat
   END INTERFACE

   INTERFACE operator (**)
      MODULE PROCEDURE QML_sub_dnMat_EXP_R
   END INTERFACE

   INTERFACE operator (+)
      MODULE PROCEDURE QML_dnMat2_PLUS_dnMat1,QML_sub_dnMat_PLUS_R,QML_sub_R_PLUS_dnMat
   END INTERFACE

   INTERFACE operator (-)
      MODULE PROCEDURE QML_dnMat2_MINUS_dnMat1,QML_sub_dnMat_MINUS_R,QML_sub_R_MINUS_dnMat
   END INTERFACE

CONTAINS
!> @brief Public subroutine which allocates a derived type dnMat.
!!
!> @author David Lauvergnat
!! @date 21/06/2018
!!
!! @param Mat                TYPE (dnMat_t):        derived type which deals with the derivatives of a matrix (as function of coordinates).
!! @param nsurf              integer (optional):    number of electronic surfaces.
!! @param ndim               integer (optional):    number of coordinates (for the derivatives).
!! @param nderiv             integer (optional):    it enables to chose the derivative order (from 0 to 2).
!! @param err_dnMat          integer (optional):    to handle the errors errors (0: no error).
!! @param name_var           character (optional):  Name of the variable from the calling subroutine (debuging purpose).
!! @param name_sub           character (optional):  Name of the calling subroutine (debuging purpose).
  SUBROUTINE QML_alloc_dnMat(Mat,nsurf,ndim,nderiv,name_var,name_sub,err_dnMat)
  IMPLICIT NONE

    TYPE (dnMat_t),    intent(inout)         :: Mat   !< derived type, which contains, matrix potential, its derivatives
    integer,           intent(in),  optional :: nsurf !< number of electronic surfaces
    integer,           intent(in),  optional :: ndim  !< number of coordinates (for the derivatives)
    integer,           intent(in),  optional :: nderiv  !< order of the derivatives [0,1,2]
    character (len=*), intent(in),  optional :: name_var,name_sub
    integer,           intent(out), optional :: err_dnMat  !< to handle the errors

    ! local variables
    integer :: nsurf_loc,ndim_loc,err_dnMat_loc,nderiv_loc



    err_dnMat_loc = 0 ! no error

    CALL QML_dealloc_dnMat(Mat,err_dnMat_loc)
    IF (err_dnMat_loc /= 0) THEN
      write(out_unitp,*) ' ERROR in QML_alloc_dnMat'
      write(out_unitp,*) ' Problem in QML_dealloc_dnMat CALL in QML_alloc_dnMat'
      IF (present(name_var)) write(out_unitp,*) '  for the variable: ',name_var
      IF (present(name_sub)) write(out_unitp,*) '  call from the subroutine: ',name_sub
      IF (present(err_dnMat)) THEN
        err_dnMat = err_dnMat_loc
        RETURN
      ELSE
        STOP 'Problem in QML_dealloc_dnMat CALL in QML_alloc_dnMat'
      END IF
    END IF

    ! test nsurf
    IF (present(nsurf)) THEN
      nsurf_loc = nsurf
    ELSE
      nsurf_loc = 1
    END IF

    ! test ndim
    IF (present(ndim)) THEN
      ndim_loc = ndim
    ELSE
      ndim_loc = 1
    END IF

    ! test nderiv
    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(3,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF
    Mat%nderiv = nderiv_loc

    !write(out_unitp,*) 'Mat%nderiv in alloc_dnMat',Mat%nderiv

    allocate(Mat%d0(nsurf_loc,nsurf_loc),stat=err_dnMat_loc)
    IF (err_dnMat_loc /= 0 .OR. nsurf_loc < 1) THEN
      write(out_unitp,*) ' ERROR in QML_alloc_dnMat'
      write(out_unitp,*) '  Problem with allocate of Mat%d0'
      write(out_unitp,*) '  nsurf > 0?',nsurf_loc
      IF (present(name_var)) write(out_unitp,*) '  for the variable: ',name_var
      IF (present(name_sub)) write(out_unitp,*) '  call from the subroutine: ',name_sub
      IF (present(err_dnMat)) THEN
        err_dnMat = err_dnMat_loc
        RETURN
      ELSE
        STOP 'Problem with allocate in QML_alloc_dnMat'
      END IF
    END IF

    IF (nderiv_loc >= 1) THEN
      allocate(Mat%d1(nsurf_loc,nsurf_loc,ndim_loc),stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0 .OR. nsurf_loc < 1 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in QML_alloc_dnMat'
        write(out_unitp,*) '  Problem with allocate of Mat%d1'
        write(out_unitp,*) '  nsurf > 0?',nsurf_loc
        write(out_unitp,*) '  ndim > 0?',ndim_loc
        IF (present(name_var)) write(out_unitp,*) '  for the variable: ',name_var
        IF (present(name_sub)) write(out_unitp,*) '  call from the subroutine: ',name_sub
        IF (present(err_dnMat)) THEN
          err_dnMat = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in QML_alloc_dnMat'
        END IF
      END IF
    END IF

    IF (nderiv_loc >= 2) THEN
      allocate(Mat%d2(nsurf_loc,nsurf_loc,ndim_loc,ndim_loc),stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0 .OR. nsurf_loc < 1 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in QML_alloc_dnMat'
        write(out_unitp,*) '  Problem with allocate of Mat%d2'
        write(out_unitp,*) '  nsurf > 0',nsurf_loc
        write(out_unitp,*) '  ndim > 0',ndim_loc
        IF (present(name_var)) write(out_unitp,*) '  for the variable: ',name_var
        IF (present(name_sub)) write(out_unitp,*) '  call from the subroutine: ',name_sub
        IF (present(err_dnMat)) THEN
          err_dnMat = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in QML_alloc_dnMat'
        END IF
      END IF
    END IF

    IF (nderiv_loc >= 3) THEN
      allocate(Mat%d3(nsurf_loc,nsurf_loc,ndim_loc,ndim_loc,ndim_loc),stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0 .OR. nsurf_loc < 1 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in QML_alloc_dnMat'
        write(out_unitp,*) '  Problem with allocate of Mat%d2'
        write(out_unitp,*) '  nsurf > 0',nsurf_loc
        write(out_unitp,*) '  ndim > 0',ndim_loc
        IF (present(name_var)) write(out_unitp,*) '  for the variable: ',name_var
        IF (present(name_sub)) write(out_unitp,*) '  call from the subroutine: ',name_sub
        IF (present(err_dnMat)) THEN
          err_dnMat = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in QML_alloc_dnMat'
        END IF
      END IF
    END IF

  END SUBROUTINE QML_alloc_dnMat
!> @brief Public subroutine which deallocates a derived type dnMat.
!!
!> @author David Lauvergnat
!! @date 21/06/2018
!!
!! @param Mat                TYPE (dnMat_t):        derived type which deals with the derivatives of a matrix (as function of coordinates).
!! @param err_dnMat       integer (optional):    to handle the errors errors (0: no error).
  SUBROUTINE QML_dealloc_dnMat(Mat,err_dnMat)
  IMPLICIT NONE

    TYPE (dnMat_t), intent(inout)         :: Mat        !< derived type, which contains, matrix potential, its derivatives
    integer,        intent(out), optional :: err_dnMat  !< to handle the errors

    ! local variables
    integer :: err_dnMat_loc

    err_dnMat_loc = 0
    IF (present(err_dnMat)) err_dnMat = 0

    IF (allocated(Mat%d0)) THEN
      deallocate(Mat%d0,stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0) THEN
        write(out_unitp,*) ' ERROR in dealloc_dnMat'
        write(out_unitp,*) '  Problem with deallocate of Mat%d0'
        IF (present(err_dnMat)) THEN
          err_dnMat = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnMat'
        END IF
      END IF
    END IF

    IF (allocated(Mat%d1)) THEN
      deallocate(Mat%d1,stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0) THEN
        write(out_unitp,*) ' ERROR in dealloc_dnMat'
        write(out_unitp,*) '  Problem with deallocate of Mat%d1'
        IF (present(err_dnMat)) THEN
          err_dnMat = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnMat'
        END IF
      END IF
    END IF

    IF (allocated(Mat%d2)) THEN
      deallocate(Mat%d2,stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0) THEN
        write(out_unitp,*) ' ERROR in dealloc_dnMat'
        write(out_unitp,*) '  Problem with deallocate of Mat%d2'
        IF (present(err_dnMat)) THEN
          err_dnMat = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnMat'
        END IF
      END IF
    END IF

    IF (allocated(Mat%d3)) THEN
      deallocate(Mat%d3,stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0) THEN
        write(out_unitp,*) ' ERROR in dealloc_dnMat'
        write(out_unitp,*) '  Problem with deallocate of Mat%d3'
        IF (present(err_dnMat)) THEN
          err_dnMat = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnMat'
        END IF
      END IF
    END IF

    Mat%nderiv = -1

  END SUBROUTINE QML_dealloc_dnMat
!> @brief Public subroutine which copies two "dnMat" derived types.
!!
!> @author David Lauvergnat
!! @date 21/06/2018
!!
!! @param dnMat1                TYPE (dnMat_t):     derived type which deals with the derivatives of a matrix (as function of coordinates).
!! @param dnMat2                TYPE (dnMat_t):     derived type which deals with the derivatives of a matrix (as function of coordinates).
  SUBROUTINE QML_sub_dnMat2_TO_dnMat1(dnMat1,dnMat2)
    CLASS (dnMat_t), intent(inout) :: dnMat1
    CLASS (dnMat_t), intent(in)    :: dnMat2

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_dnMat2_TO_dnMat1'

    nderiv_loc = QML_get_nderiv_FROM_dnMat(dnMat2)
    nsurf_loc  = QML_get_nsurf_FROM_dnMat(dnMat2)
    ndim_loc   = QML_get_ndim_FROM_dnMat(dnMat2)

    !write(out_unitp,*) 'in ',name_sub,' ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc


    IF (nderiv_loc < 0 .OR. nsurf_loc < 1 .OR. (nderiv_loc > 0  .AND. ndim_loc < 1)) RETURN


    CALL QML_alloc_dnMat(dnMat1,nsurf_loc,ndim_loc,nderiv_loc,name_var='dnMat1',name_sub=name_sub)


    IF (nderiv_loc == 0) THEN
       dnMat1%d0 = dnMat2%d0
    ELSE IF (nderiv_loc == 1) THEN
       dnMat1%d0 = dnMat2%d0
       dnMat1%d1 = dnMat2%d1
    ELSE IF (nderiv_loc == 2) THEN
       dnMat1%d0 = dnMat2%d0
       dnMat1%d1 = dnMat2%d1
       dnMat1%d2 = dnMat2%d2
    ELSE IF (nderiv_loc == 3) THEN
       dnMat1%d0 = dnMat2%d0
       dnMat1%d1 = dnMat2%d1
       dnMat1%d2 = dnMat2%d2
       dnMat1%d3 = dnMat2%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END SUBROUTINE QML_sub_dnMat2_TO_dnMat1
!> @brief Public subroutine which copies a dnS derived type to one element of dnMat derived type.
!!
!> @author David Lauvergnat
!! @date 25/06/2018
!!
!! @param Mat                   TYPE (dnMat_t):    derived type which deals with the derivatives of a matrix.
!! @param S                     TYPE(dnS):       derived type which deals with the derivatives of a scalar.
!! @param i,j                   integer (optional) indices of the matrix element. If not present i=j=1
  SUBROUTINE QML_sub_dnS_TO_dnMat(S,Mat,i,j)
    USE mod_dnS
    TYPE (dnMat_t),     intent(inout) :: Mat
    TYPE (dnS_t),       intent(in)    :: S
    integer, optional,  intent(in)    :: i,j

    integer :: nderiv_dnMat,nsurf_dnMat,ndim_dnMat,nderiv_dnS,ndim_dnS
    integer :: i_loc,j_loc

    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_dnS_TO_dnMat'


    nderiv_dnS = QML_get_nderiv_FROM_dnS(S)
    ndim_dnS   = QML_get_ndim_FROM_dnS(S)

    nderiv_dnMat = QML_get_nderiv_FROM_dnMat(Mat)
    nsurf_dnMat  = QML_get_nsurf_FROM_dnMat(Mat)
    ndim_dnMat   = QML_get_ndim_FROM_dnMat(Mat)

    i_loc = 1
    j_loc = 1
    IF (present(i)) i_loc = i
    IF (present(j)) j_loc = j


    IF (i_loc < 1 .OR. i_loc > nsurf_dnMat .OR. j_loc < 1 .OR. j_loc > nsurf_dnMat) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The matrix indexes, (',i_loc,j_loc,') are out of range [1...',nsurf_dnMat,']'
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF

    IF (nderiv_dnS == -1) THEN
      IF (nderiv_dnMat == -1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' dnMat is not allocated.'
        write(out_unitp,*) 'It should never append! Check the source.'
        STOP 'dnMat is not allocated.'
      END IF
      ! S (dnS) is a constant
      ! value
      Mat%d0(i_loc,j_loc) = QML_get_d0_FROM_dnS(S)

      ! 1st order derivatives
      IF (nderiv_dnMat >= 1) Mat%d1(i_loc,j_loc,:) = ZERO

      ! 2d order derivatives
      IF (nderiv_dnMat >= 2) Mat%d2(i_loc,j_loc,:,:) = ZERO
    ELSE

      IF ( QML_check_notalloc_dnmat(Mat,nderiv_dnS) .OR.                  &
           nderiv_dnS /= nderiv_dnMat  .OR.  ndim_dnS /= ndim_dnMat .OR.  &
           nsurf_dnMat < 1 ) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' dnMat is not allocated or ...'
        write(out_unitp,*) '  ... nderiv from dnMat or dnS are different or ...'
        write(out_unitp,*) '  ... ndim from dnMat or dnS are different or ...'
        write(out_unitp,*) '  ... nsurf from dnMat is < 1'

        write(out_unitp,*) 'nderiv from dnMat and dnS:',nderiv_dnMat,nderiv_dnS
        write(out_unitp,*) 'ndim   from dnMat and dnS:',ndim_dnMat,ndim_dnS
        write(out_unitp,*) 'nsurf  from dnMat        :',nsurf_dnMat

        write(out_unitp,*) 'It should never append! Check the source'
        STOP 'dnMat is not allocated or inconsistent ndim,nderiv parameters.'
      END IF

      ! value
      Mat%d0(i_loc,j_loc) = QML_get_d0_FROM_dnS(S)

      ! 1st order derivatives
      IF (nderiv_dnS >= 1) THEN
        CALL QML_sub_get_dn_FROM_dnS(S,d1=Mat%d1(i_loc,j_loc,:))
      END IF

      ! 2d order derivatives
      IF (nderiv_dnS >= 2) then
        CALL QML_sub_get_dn_FROM_dnS(S,d2=Mat%d2(i_loc,j_loc,:,:))
      END IF

      ! 3d order derivatives
      IF (nderiv_dnS >= 3) then
        CALL QML_sub_get_dn_FROM_dnS(S,d3=Mat%d3(i_loc,j_loc,:,:,:))
      END IF

    END IF

  END SUBROUTINE QML_sub_dnS_TO_dnMat
!> @brief Public subroutine which copies a dnS derived type to one element of dnMat derived type.
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param Mat                   TYPE (dnMat_t):    derived type which deals with the derivatives of a matrix.
!! @param MatOFS                TYPE(dnS):       matrix of derived type which deals with the derivatives of a scalar.
  SUBROUTINE QML_set_dnMat_FROM_MatOFdnS(Mat,MatOFS)
    USE mod_dnS
    CLASS (dnMat_t),   intent(inout) :: Mat
    TYPE (dnS_t),      intent(in)    :: MatOFS(:,:)

    integer :: nderiv_dnMat,nsurf_dnMat,ndim_dnMat,nderiv_dnS,ndim_dnS
    integer :: i,j

    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_set_dnMat_FROM_MatOFdnS'

    IF (lbound(Mat%d0,dim=1) /= lbound(MatOFS,dim=1) .OR. ubound(Mat%d0,dim=1) /= ubound(MatOFS,dim=1) .OR. &
        lbound(Mat%d0,dim=2) /= lbound(MatOFS,dim=2) .OR. ubound(Mat%d0,dim=2) /= ubound(MatOFS,dim=2) ) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  the matrices have not the same dimensions'
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF

    DO i=lbound(MatOFS,dim=2),ubound(MatOFS,dim=2)
    DO j=lbound(MatOFS,dim=1),ubound(MatOFS,dim=1)
      CALL QML_sub_dnS_TO_dnMat(MatOFS(i,j),Mat,i,j)
    END DO
    END DO

  END SUBROUTINE QML_set_dnMat_FROM_MatOFdnS
  SUBROUTINE QML_Mat_wADDTO_dnMat2_ider(Mat1,w1,dnMat2,ider)
    real (kind=Rkind),  intent(in)            :: Mat1(:,:)
    TYPE (dnMat_t),     intent(inout)         :: dnMat2
    integer,            intent(in),  optional :: ider(:)
    real (kind=Rkind),  intent(in)            :: w1

    integer :: nderiv,nsurf,ndim
    integer :: i1,i1i,i1f
    integer :: i2,i2i,i2f
    integer :: i3,i3i,i3f

    character (len=*), parameter :: name_sub='QML_Mat_wADDTO_dnMat2_ider'

    nderiv = QML_get_nderiv_FROM_dnMat(dnMat2)
    nsurf  = QML_get_nsurf_FROM_dnMat(dnMat2)
    ndim   = QML_get_ndim_FROM_dnMat(dnMat2)

    IF (.NOT. allocated(dnMat2%d0)) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  dnMat2%d0 is not allocated.'
      write(out_unitp,*) ' CHECK the fortran source!!'
      STOP
    END IF

    IF (.NOT. all(shape(Mat1) == shape(dnMat2%d0))) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  The shape of Mat1 dnMat2%d0 must be equal.'
      write(out_unitp,*) '  shape(Mat1):      ',shape(Mat1)
      write(out_unitp,*) '  shape(dnMat2%d0): ',shape(dnMat2%d0)
      write(out_unitp,*) ' CHECK the fortran source!!'
      STOP
    END IF
    IF (present(ider)) THEN
      IF (size(ider) > nderiv) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' size(ider) cannot be > and nderiv.'
        write(out_unitp,*) ' size(ider)',size(ider)
        write(out_unitp,*) ' nderiv    ',nderiv
        write(out_unitp,*) ' CHECK the fortran source!!'
        STOP
      END IF
      IF (any(ider < 0) .OR. any(ider > ndim)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Some ider(:) values are out-of-range.'
        write(out_unitp,*) ' ider(:)',ider
        write(out_unitp,'(a,i0,a)') ' derivative range [0:',ndim,']'
        write(out_unitp,*) ' CHECK the fortran source!!'
        STOP
      END IF
    END IF



    IF (present(ider)) THEN

      IF (size(ider) > 0) THEN
        IF (ider(1) == 0) THEN
          i1i = 1
          i1f = ndim
        ELSE
          i1i = ider(1)
          i1f = ider(1)
        END IF
      END IF
      IF (size(ider) > 1) THEN
        IF (ider(2) == 0) THEN
          i2i = 1
          i2f = ndim
        ELSE
          i2i = ider(2)
          i2f = ider(2)
        END IF
      END IF
      IF (size(ider) > 2) THEN
        IF (ider(3) == 0) THEN
          i3i = 1
          i3f = ndim
        ELSE
          i3i = ider(3)
          i3f = ider(3)
        END IF
      END IF


      SELECT CASE (size(ider))
      CASE (0)
        dnMat2%d0(:,:) = w1*Mat1 + dnMat2%d0

      CASE (1)
        DO i1=i1i,i1f
          dnMat2%d1(:,:,i1) = w1*Mat1 + dnMat2%d1(:,:,i1)
        END DO

      CASE (2)
        DO i2=i2i,i2f
        DO i1=i1i,i1f
          dnMat2%d2(:,:,i1,i2) = w1*Mat1 + dnMat2%d2(:,:,i1,i2)
        END DO
        END DO

      CASE (3)

        !IF (present(ider)) write(6,*) 'ider',ider

        DO i3=i3i,i3f
        DO i2=i2i,i2f
        DO i1=i1i,i1f
          dnMat2%d3(:,:,i1,i2,i3) = w1*Mat1 + dnMat2%d3(:,:,i1,i2,i3)
        END DO
        END DO
        END DO

      CASE Default
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' size(ider) > 3 is NOT possible.'
        write(out_unitp,*) '   ider',ider
        write(out_unitp,*) 'It should never append! Check the source'
        STOP
      END SELECT
    ELSE
      dnMat2%d0(:,:) = w1*Mat1 + dnMat2%d0
    END IF

  END SUBROUTINE QML_Mat_wADDTO_dnMat2_ider

!> @brief Public function which calculate set dnMat to zero (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                   TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param set_dnMat_TO_zero  TYPE (dnMat_t) (result):  dnMat derived type
  SUBROUTINE QML_set_dnMat_TO_zero(dnMat)
    TYPE (dnMat_t), intent(inout) :: dnMat

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_set_dnMat_TO_zero'


    CALL QML_set_dnMat_TO_R(dnMat,ZERO)

  END SUBROUTINE QML_set_dnMat_TO_zero
!> @brief Public function which calculate set dnMat to R (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param R                  real:                     some real number
!! @param set_dnMat_TO_R  TYPE (dnMat_t) (result):  dnMat derived type
  SUBROUTINE QML_set_dnMat_TO_R(dnMat,R)

    CLASS (dnMat_t), intent(inout) :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_set_dnMat_TO_R'

    nderiv_loc = QML_get_nderiv_FROM_dnMat(dnMat)

    !write(out_unitp,*) 'nderiv',nderiv_loc


    IF (nderiv_loc == 0) THEN
       dnMat%d0 = R
    ELSE IF (nderiv_loc == 1) THEN
       dnMat%d0 = R
       dnMat%d1 = ZERO
    ELSE IF (nderiv_loc == 2) THEN
       dnMat%d0 = R
       dnMat%d1 = ZERO
       dnMat%d2 = ZERO
    ELSE IF (nderiv_loc == 3) THEN
       dnMat%d0 = R
       dnMat%d1 = ZERO
       dnMat%d2 = ZERO
       dnMat%d3 = ZERO
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 or nderiv < 0 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END SUBROUTINE QML_set_dnMat_TO_R
!> @brief Public function which calculate dnMat*R (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param R                  real:                     some real number
!! @param sub_dnMat_TIME_R TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_sub_dnMat_TIME_R(dnMat,R) RESULT (sub_dnMat_TIME_R)

    TYPE (dnMat_t)                 :: sub_dnMat_TIME_R
    TYPE (dnMat_t),    intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_dnMat_TIME_R'

    nderiv_loc = QML_get_nderiv_FROM_dnMat(dnMat)
    nsurf_loc  = QML_get_nsurf_FROM_dnMat(dnMat)
    ndim_loc   = QML_get_ndim_FROM_dnMat(dnMat)

    !write(out_unitp,*) 'ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc

    CALL QML_alloc_dnMat(sub_dnMat_TIME_R,nsurf_loc,ndim_loc,nderiv_loc,&
                         name_var='sub_dnMat_TIME_R',name_sub=name_sub)

    !write(out_unitp,*) 'nderiv',nderiv_loc


    IF (nderiv_loc == 0) THEN
       sub_dnMat_TIME_R%d0 = dnMat%d0 * R

    ELSE IF (nderiv_loc == 1) THEN
       sub_dnMat_TIME_R%d0 = dnMat%d0 * R
       sub_dnMat_TIME_R%d1 = dnMat%d1 * R

    ELSE IF (nderiv_loc == 2) THEN
       sub_dnMat_TIME_R%d0 = dnMat%d0 * R
       sub_dnMat_TIME_R%d1 = dnMat%d1 * R
       sub_dnMat_TIME_R%d2 = dnMat%d2 * R
    ELSE IF (nderiv_loc == 3) THEN
       sub_dnMat_TIME_R%d0 = dnMat%d0 * R
       sub_dnMat_TIME_R%d1 = dnMat%d1 * R
       sub_dnMat_TIME_R%d2 = dnMat%d2 * R
       sub_dnMat_TIME_R%d3 = dnMat%d3 * R
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION QML_sub_dnMat_TIME_R
!> @brief Public function which calculate R*dnMat (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param R                  real:                     some real number
!! @param sub_R_TIME_dnMat TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_sub_R_TIME_dnMat(R,dnMat)  RESULT(sub_R_TIME_dnMat)

    TYPE (dnMat_t)                :: sub_R_TIME_dnMat
    TYPE (dnMat_t),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_R_TIME_dnMat'

    nderiv_loc = QML_get_nderiv_FROM_dnMat(dnMat)
    nsurf_loc  = QML_get_nsurf_FROM_dnMat(dnMat)
    ndim_loc   = QML_get_ndim_FROM_dnMat(dnMat)

    !write(out_unitp,*) 'ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc

    CALL QML_alloc_dnMat(sub_R_TIME_dnMat,nsurf_loc,ndim_loc,nderiv_loc,&
                         name_var='sub_R_TIME_dnMat',name_sub=name_sub)

    !write(out_unitp,*) 'nderiv',nderiv_loc


    IF (nderiv_loc == 0) THEN
       sub_R_TIME_dnMat%d0 = dnMat%d0 * R

    ELSE IF (nderiv_loc == 1) THEN
       sub_R_TIME_dnMat%d0 = dnMat%d0 * R
       sub_R_TIME_dnMat%d1 = dnMat%d1 * R

    ELSE IF (nderiv_loc == 2) THEN
       sub_R_TIME_dnMat%d0 = dnMat%d0 * R
       sub_R_TIME_dnMat%d1 = dnMat%d1 * R
       sub_R_TIME_dnMat%d2 = dnMat%d2 * R
    ELSE IF (nderiv_loc == 3) THEN
       sub_R_TIME_dnMat%d0 = dnMat%d0 * R
       sub_R_TIME_dnMat%d1 = dnMat%d1 * R
       sub_R_TIME_dnMat%d2 = dnMat%d2 * R
       sub_R_TIME_dnMat%d3 = dnMat%d3 * R
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION QML_sub_R_TIME_dnMat
!> @brief Public function which calculate dnMat1+dnMat2 (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param dnMat1                    TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param dnMat2                    TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param dnMat2_PLUS_dnMat1 TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_dnMat2_PLUS_dnMat1(dnMat1,dnMat2)  RESULT(dnMat2_PLUS_dnMat1)
    TYPE (dnMat_t)                :: dnMat2_PLUS_dnMat1
    TYPE (dnMat_t), intent(in)    :: dnMat1,dnMat2

    integer :: nderiv,nsurf,ndim
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_dnMat2_PLUS_dnMat1'

    nderiv = min(QML_get_nderiv_FROM_dnMat(dnMat1),QML_get_nderiv_FROM_dnMat(dnMat2))
    nsurf  = min(QML_get_nsurf_FROM_dnMat(dnMat1), QML_get_nsurf_FROM_dnMat(dnMat2))
    ndim   = min(QML_get_ndim_FROM_dnMat(dnMat1),  QML_get_ndim_FROM_dnMat(dnMat2))

    !write(out_unitp,*) 'in ',name_sub,' nsurf,ndim,nderiv',nsurf,ndim,nderiv

    CALL QML_dealloc_dnMat(dnMat2_PLUS_dnMat1)

    IF (nderiv < 0 .OR. nsurf < 1 .OR. (nderiv > 0  .AND. ndim < 1)) RETURN

    CALL QML_alloc_dnMat(dnMat2_PLUS_dnMat1,nsurf,ndim,nderiv,          &
                         name_var='dnMat2_PLUS_dnMat1',name_sub=name_sub)

    IF (nderiv == 0) THEN
       dnMat2_PLUS_dnMat1%d0 = dnMat1%d0 + dnMat2%d0
    ELSE IF (nderiv == 1) THEN
       dnMat2_PLUS_dnMat1%d0 = dnMat1%d0 + dnMat2%d0
       dnMat2_PLUS_dnMat1%d1 = dnMat1%d1 + dnMat2%d1
    ELSE IF (nderiv == 2) THEN
       dnMat2_PLUS_dnMat1%d0 = dnMat1%d0 + dnMat2%d0
       dnMat2_PLUS_dnMat1%d1 = dnMat1%d1 + dnMat2%d1
       dnMat2_PLUS_dnMat1%d2 = dnMat1%d2 + dnMat2%d2
    ELSE IF (nderiv == 3) THEN
       dnMat2_PLUS_dnMat1%d0 = dnMat1%d0 + dnMat2%d0
       dnMat2_PLUS_dnMat1%d1 = dnMat1%d1 + dnMat2%d1
       dnMat2_PLUS_dnMat1%d2 = dnMat1%d2 + dnMat2%d2
       dnMat2_PLUS_dnMat1%d3 = dnMat1%d3 + dnMat2%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION QML_dnMat2_PLUS_dnMat1
!> @brief Public function which calculate dnMat+R (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param R                  real:                     some real number
!! @param sub_dnMat_EXP_R TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_sub_dnMat_PLUS_R(dnMat,R)  RESULT (sub_dnMat_PLUS_R)

    TYPE (dnMat_t)                :: sub_dnMat_PLUS_R
    TYPE (dnMat_t),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_dnMat_PLUS_R'


    sub_dnMat_PLUS_R    = dnMat

    sub_dnMat_PLUS_R%d0 = sub_dnMat_PLUS_R%d0 + R

    ! the derivatives of R are zero => nothing to be add to %d1 and %d2

  END FUNCTION QML_sub_dnMat_PLUS_R
!> @brief Public function which calculate R+dnMat (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param R                  real:                     some real number
!! @param sub_R_PLUS_dnMat TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_sub_R_PLUS_dnMat(R,dnMat) RESULT (sub_R_PLUS_dnMat)

    TYPE (dnMat_t)                :: sub_R_PLUS_dnMat
    TYPE (dnMat_t),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_R_PLUS_dnMat'


    sub_R_PLUS_dnMat    = dnMat

    sub_R_PLUS_dnMat%d0 = sub_R_PLUS_dnMat%d0 + R

    ! the derivatives of R are zero

  END FUNCTION QML_sub_R_PLUS_dnMat
!> @brief Public function which calculate dnMat1-dnMat2 (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param dnMat1                    TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param dnMat2                    TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param dnMat2_MINUS_dnMat1 TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_dnMat2_MINUS_dnMat1(dnMat1,dnMat2) RESULT (dnMat2_MINUS_dnMat1)
    TYPE (dnMat_t)                :: dnMat2_MINUS_dnMat1
    TYPE (dnMat_t), intent(in)    :: dnMat1,dnMat2

    integer :: nderiv,nsurf,ndim
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_dnMat2_MINUS_dnMat1'

    nderiv = min(QML_get_nderiv_FROM_dnMat(dnMat1),QML_get_nderiv_FROM_dnMat(dnMat2))
    nsurf  = min(QML_get_nsurf_FROM_dnMat(dnMat1), QML_get_nsurf_FROM_dnMat(dnMat2))
    ndim   = min(QML_get_ndim_FROM_dnMat(dnMat1),  QML_get_ndim_FROM_dnMat(dnMat2))


    CALL QML_dealloc_dnMat(dnMat2_MINUS_dnMat1)

    IF (nderiv < 0 .OR. nsurf < 1 .OR. (nderiv > 0  .AND. ndim < 1)) RETURN

    CALL QML_alloc_dnMat(dnMat2_MINUS_dnMat1,nsurf,ndim,nderiv,         &
                         name_var='dnMat2_MINUS_dnMat1',name_sub=name_sub)

    IF (nderiv == 0) THEN
       dnMat2_MINUS_dnMat1%d0 = dnMat1%d0 - dnMat2%d0
    ELSE IF (nderiv == 1) THEN
       dnMat2_MINUS_dnMat1%d0 = dnMat1%d0 - dnMat2%d0
       dnMat2_MINUS_dnMat1%d1 = dnMat1%d1 - dnMat2%d1
    ELSE IF (nderiv == 2) THEN
       dnMat2_MINUS_dnMat1%d0 = dnMat1%d0 - dnMat2%d0
       dnMat2_MINUS_dnMat1%d1 = dnMat1%d1 - dnMat2%d1
       dnMat2_MINUS_dnMat1%d2 = dnMat1%d2 - dnMat2%d2
    ELSE IF (nderiv == 3) THEN
       dnMat2_MINUS_dnMat1%d0 = dnMat1%d0 - dnMat2%d0
       dnMat2_MINUS_dnMat1%d1 = dnMat1%d1 - dnMat2%d1
       dnMat2_MINUS_dnMat1%d2 = dnMat1%d2 - dnMat2%d2
       dnMat2_MINUS_dnMat1%d3 = dnMat1%d3 - dnMat2%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION QML_dnMat2_MINUS_dnMat1
!> @brief Public function which calculate dnMat-R (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                  TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param R                    real:                     some real number
!! @param sub_dnMat_MINUS_R TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_sub_dnMat_MINUS_R(dnMat,R) RESULT (sub_dnMat_MINUS_R)

    TYPE (dnMat_t)                :: sub_dnMat_MINUS_R
    TYPE (dnMat_t),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_dnMat_MINUS_R'


    sub_dnMat_MINUS_R = dnMat

    sub_dnMat_MINUS_R%d0 = dnMat%d0 - R

    ! the derivatives of R are zero

  END FUNCTION QML_sub_dnMat_MINUS_R
!> @brief Public function which calculate R-dnMat (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                  TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param R                    real:                     some real number
!! @param sub_R_MINUS_dnMat TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_sub_R_MINUS_dnMat(R,dnMat) RESULT (sub_R_MINUS_dnMat)

    TYPE (dnMat_t)                :: sub_R_MINUS_dnMat
    TYPE (dnMat_t),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_R_MINUS_dnMat'

    nderiv_loc = QML_get_nderiv_FROM_dnMat(dnMat)
    nsurf_loc  = QML_get_nsurf_FROM_dnMat(dnMat)
    ndim_loc   = QML_get_ndim_FROM_dnMat(dnMat)

    !write(out_unitp,*) 'ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc

    CALL QML_alloc_dnMat(sub_R_MINUS_dnMat,nsurf_loc,ndim_loc,nderiv_loc,&
                         name_var='sub_R_MINUS_dnMat',name_sub=name_sub)

    !write(out_unitp,*) 'nderiv',nderiv_loc
    IF (nderiv_loc == 0) THEN
       sub_R_MINUS_dnMat%d0 = R - dnMat%d0

    ELSE IF (nderiv_loc == 1) THEN
       sub_R_MINUS_dnMat%d0 = R - dnMat%d0
       sub_R_MINUS_dnMat%d1 =   - dnMat%d1


    ELSE IF (nderiv_loc == 2) THEN
       sub_R_MINUS_dnMat%d0 = R - dnMat%d0
       sub_R_MINUS_dnMat%d1 =   - dnMat%d1
       sub_R_MINUS_dnMat%d2 =   - dnMat%d2

    ELSE IF (nderiv_loc == 3) THEN
       sub_R_MINUS_dnMat%d0 = R - dnMat%d0
       sub_R_MINUS_dnMat%d1 =   - dnMat%d1
       sub_R_MINUS_dnMat%d2 =   - dnMat%d2
       sub_R_MINUS_dnMat%d3 =   - dnMat%d3

    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF

  END FUNCTION QML_sub_R_MINUS_dnMat
!> @brief Public function which calculate dnMat**R (and derivatives).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                TYPE (dnMat_t):           derived type which deals with the derivatives of a matrix.
!! @param R                  real:                     exponent
!! @param sub_dnMat_EXP_R TYPE (dnMat_t) (result):  dnMat derived type
  FUNCTION QML_sub_dnMat_EXP_R(dnMat,R) RESULT (sub_dnMat_EXP_R)

    TYPE (dnMat_t)                 :: sub_dnMat_EXP_R
    TYPE (dnMat_t),    intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc,id,jd,kd
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_sub_dnMat_EXP_R'

    nderiv_loc = QML_get_nderiv_FROM_dnMat(dnMat)
    nsurf_loc  = QML_get_nsurf_FROM_dnMat(dnMat)
    ndim_loc   = QML_get_ndim_FROM_dnMat(dnMat)

    !write(out_unitp,*) 'ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc

    CALL QML_alloc_dnMat(sub_dnMat_EXP_R,nsurf_loc,ndim_loc,            &
                nderiv_loc,name_var='sub_dnMat_EXP_R',name_sub=name_sub)

    !write(out_unitp,*) 'nderiv',nderiv_loc


    IF (nderiv_loc == 0) THEN
       sub_dnMat_EXP_R%d0 = dnMat%d0 ** R

    ELSE IF (nderiv_loc == 1) THEN
       sub_dnMat_EXP_R%d0 = dnMat%d0 ** R

       DO id=1,ndim_loc
         sub_dnMat_EXP_R%d1(:,:,id) = R * dnMat%d0 ** (R-ONE) * dnMat%d1(:,:,id)
       END DO

    ELSE IF (nderiv_loc == 2) THEN
       sub_dnMat_EXP_R%d0 = dnMat%d0 ** R

       DO id=1,ndim_loc
         sub_dnMat_EXP_R%d1(:,:,id) = R * dnMat%d0 ** (R-ONE) * dnMat%d1(:,:,id)
       END DO

       DO id=1,ndim_loc
       DO jd=1,ndim_loc
         sub_dnMat_EXP_R%d2(:,:,jd,id) = R*(R-ONE) * dnMat%d0 ** (R-TWO) * dnMat%d1(:,:,id) * dnMat%d1(:,:,jd) + &
                                            R * dnMat%d0 ** (R-ONE) * dnMat%d2(:,:,jd,id)
       END DO
       END DO

    ELSE IF (nderiv_loc == 3) THEN
       sub_dnMat_EXP_R%d0 = dnMat%d0 ** R

       DO id=1,ndim_loc
         sub_dnMat_EXP_R%d1(:,:,id) = R * dnMat%d0 ** (R-ONE) * dnMat%d1(:,:,id)
       END DO

       DO id=1,ndim_loc
       DO jd=1,ndim_loc
         sub_dnMat_EXP_R%d2(:,:,jd,id) = R*(R-ONE) * dnMat%d0 ** (R-TWO) * dnMat%d1(:,:,id) * dnMat%d1(:,:,jd) + &
                                            R * dnMat%d0 ** (R-ONE) * dnMat%d2(:,:,jd,id)
       END DO
       END DO

       DO id=1,ndim_loc
       DO jd=1,ndim_loc
       DO kd=1,ndim_loc
         sub_dnMat_EXP_R%d3(:,:,kd,jd,id) =                             &
                            R*(R-ONE)*(R-TWO) * dnMat%d0**(R-THREE) *   &
                   dnMat%d1(:,:,id)*dnMat%d1(:,:,jd)*dnMat%d1(:,:,kd) + &
                            R*(R-ONE) * dnMat%d0**(R-TWO) * (           &
                               dnMat%d2(:,:,jd,id)*dnMat%d1(:,:,kd) +   &
                               dnMat%d2(:,:,kd,id)*dnMat%d1(:,:,jd) +   &
                               dnMat%d2(:,:,kd,jd)*dnMat%d1(:,:,id) ) + &
                             R * dnMat%d0**(R-ONE) * dnMat%d3(:,:,kd,jd,id)

       END DO
       END DO
       END DO

    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION QML_sub_dnMat_EXP_R
  FUNCTION QML_TRANSPOSE_dnMat(dnMat)  RESULT(TransdnMat) ! check with t(t(dnmat))-dnMat
    TYPE (dnMat_t)                :: TransdnMat
    TYPE (dnMat_t), intent(in)    :: dnMat

    integer :: nderiv,nsurf,ndim,id,jd,kd
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_TRANSPOSE_dnMat'

    nderiv = QML_get_nderiv_FROM_dnMat(dnMat)
    nsurf  = QML_get_nsurf_FROM_dnMat(dnMat)
    ndim   = QML_get_ndim_FROM_dnMat(dnMat)

    !write(out_unitp,*) 'in ',name_sub,' nsurf,ndim,nderiv',nsurf,ndim,nderiv

    CALL QML_dealloc_dnMat(TransdnMat)

    IF (nderiv < 0 .OR. nsurf < 1 .OR. (nderiv > 0  .AND. ndim < 1)) RETURN

    CALL QML_alloc_dnMat(TransdnMat,nsurf,ndim,nderiv,                  &
                         name_var='TransdnMat',name_sub=name_sub)

    IF (nderiv >= 0) THEN
      TransdnMat%d0(:,:) = transpose(dnMat%d0)
    END IF

    IF (nderiv >= 1) THEN
      DO id=1,ndim
        TransdnMat%d1(:,:,id) = transpose(dnMat%d1(:,:,id))
      END DO
    END IF

    IF (nderiv >= 2) THEN
      DO id=1,ndim
      DO jd=1,ndim
        TransdnMat%d2(:,:,jd,id) = transpose(dnMat%d2(:,:,jd,id))
      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN
      DO id=1,ndim
      DO jd=1,ndim
      DO kd=1,ndim
        TransdnMat%d3(:,:,kd,jd,id) = transpose(dnMat%d3(:,:,kd,jd,id))
      END DO
      END DO
      END DO
    END IF

  END FUNCTION QML_TRANSPOSE_dnMat

  FUNCTION QML_MATMUL_dnMat1_dnMat2(dnMat1,dnMat2)  RESULT(MatmuldnMat)
    TYPE (dnMat_t)                :: MatmuldnMat
    TYPE (dnMat_t), intent(in)    :: dnMat1,dnMat2

    integer :: nderiv,nsurf,ndim,id,jd,kd
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_MATMUL_dnMat1_dnMat2'


    nderiv = min(QML_get_nderiv_FROM_dnMat(dnMat1),QML_get_nderiv_FROM_dnMat(dnMat2))
    nsurf  = min(QML_get_nsurf_FROM_dnMat(dnMat1), QML_get_nsurf_FROM_dnMat(dnMat2))
    ndim   = min(QML_get_ndim_FROM_dnMat(dnMat1),  QML_get_ndim_FROM_dnMat(dnMat2))


    !write(out_unitp,*) 'in ',name_sub,' nsurf,ndim,nderiv',nsurf,ndim,nderiv

    CALL QML_dealloc_dnMat(MatmuldnMat)

    IF (nderiv < 0 .OR. nsurf < 1 .OR. (nderiv > 0  .AND. ndim < 1)) RETURN

    CALL QML_alloc_dnMat(MatmuldnMat,nsurf,ndim,nderiv,                 &
                         name_var='MatmuldnMat',name_sub=name_sub)

    IF (nderiv >= 0) THEN
      MatmuldnMat%d0(:,:) = matmul(dnMat1%d0,dnMat2%d0)
    END IF

    IF (nderiv >= 1) THEN
      DO id=1,ndim
        MatmuldnMat%d1(:,:,id) = matmul(dnMat1%d0,dnMat2%d1(:,:,id)) +  &
                                 matmul(dnMat1%d1(:,:,id),dnMat2%d0)
      END DO
    END IF

    IF (nderiv >= 2) THEN
      DO id=1,ndim
      DO jd=1,ndim
        MatmuldnMat%d2(:,:,jd,id) =                                     &
                    matmul(dnMat1%d0,           dnMat2%d2(:,:,jd,id)) + &
                    matmul(dnMat1%d1(:,:,jd),   dnMat2%d1(:,:,id))    + &
                    matmul(dnMat1%d1(:,:,id),   dnMat2%d1(:,:,jd))    + &
                    matmul(dnMat1%d2(:,:,jd,id),dnMat2%d0)

      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN
      DO id=1,ndim
      DO jd=1,ndim
      DO kd=1,ndim
        MatmuldnMat%d3(:,:,kd,jd,id) =                                  &
             matmul(dnMat1%d0,              dnMat2%d3(:,:,kd,jd,id))  + &
             matmul(dnMat1%d1(:,:,kd),      dnMat2%d2(:,:,jd,id))     + &
             matmul(dnMat1%d1(:,:,id),      dnMat2%d2(:,:,kd,jd))     + &
             matmul(dnMat1%d1(:,:,jd),      dnMat2%d2(:,:,id,kd))     + &
             matmul(dnMat1%d2(:,:,jd,id),   dnMat2%d1(:,:,kd))        + &
             matmul(dnMat1%d2(:,:,kd,jd),   dnMat2%d1(:,:,id))        + &
             matmul(dnMat1%d2(:,:,id,kd),   dnMat2%d1(:,:,jd))        + &
             matmul(dnMat1%d3(:,:,kd,jd,id),dnMat2%d0)

      END DO
      END DO
      END DO
    END IF

  END FUNCTION QML_MATMUL_dnMat1_dnMat2

  FUNCTION QML_MATMUL_dnMat1_Mat2(dnMat1,Mat2)  RESULT(MatmuldnMat)
    TYPE (dnMat_t)                  :: MatmuldnMat
    TYPE (dnMat_t),   intent(in)    :: dnMat1
    real(kind=Rkind), intent(in)    :: Mat2(:,:)

    integer :: nderiv,nsurf,ndim,id,jd,kd
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_MATMUL_dnMat1_Mat2'


    nderiv = QML_get_nderiv_FROM_dnMat(dnMat1)
    nsurf  = QML_get_nsurf_FROM_dnMat(dnMat1)
    ndim   = QML_get_ndim_FROM_dnMat(dnMat1)


    !write(out_unitp,*) 'in ',name_sub,' nsurf,ndim,nderiv',nsurf,ndim,nderiv

    CALL QML_dealloc_dnMat(MatmuldnMat)

    IF (nderiv < 0 .OR. nsurf < 1 .OR. (nderiv > 0  .AND. ndim < 1)) RETURN

    CALL QML_alloc_dnMat(MatmuldnMat,nsurf,ndim,nderiv,                 &
                         name_var='MatmuldnMat',name_sub=name_sub)

    IF (nderiv >= 0) THEN
      MatmuldnMat%d0(:,:) = matmul(dnMat1%d0,Mat2)
    END IF

    IF (nderiv >= 1) THEN
      DO id=1,ndim
        MatmuldnMat%d1(:,:,id) = matmul(dnMat1%d1(:,:,id),Mat2)
      END DO
    END IF

    IF (nderiv >= 2) THEN
      DO id=1,ndim
      DO jd=1,ndim
        MatmuldnMat%d2(:,:,jd,id) = matmul(dnMat1%d2(:,:,jd,id),Mat2)
      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN
      DO id=1,ndim
      DO jd=1,ndim
      DO kd=1,ndim
        MatmuldnMat%d3(:,:,kd,jd,id) = matmul(dnMat1%d3(:,:,kd,jd,id),Mat2)
      END DO
      END DO
      END DO
    END IF

  END FUNCTION QML_MATMUL_dnMat1_Mat2
  FUNCTION QML_MATMUL_Mat1_dnMat2(Mat1,dnMat2)  RESULT(MatmuldnMat)
    TYPE (dnMat_t)                  :: MatmuldnMat
    real(kind=Rkind), intent(in)    :: Mat1(:,:)
    TYPE (dnMat_t),   intent(in)    :: dnMat2

    integer :: nderiv,nsurf,ndim,id,jd,kd
    integer :: err_dnMat_loc
    character (len=*), parameter :: name_sub='QML_MATMUL_Mat1_dnMat2'


    nderiv = QML_get_nderiv_FROM_dnMat(dnMat2)
    nsurf  = QML_get_nsurf_FROM_dnMat(dnMat2)
    ndim   = QML_get_ndim_FROM_dnMat(dnMat2)


    !write(out_unitp,*) 'in ',name_sub,' nsurf,ndim,nderiv',nsurf,ndim,nderiv

    CALL QML_dealloc_dnMat(MatmuldnMat)

    IF (nderiv < 0 .OR. nsurf < 1 .OR. (nderiv > 0  .AND. ndim < 1)) RETURN

    CALL QML_alloc_dnMat(MatmuldnMat,nsurf,ndim,nderiv,                 &
                         name_var='MatmuldnMat',name_sub=name_sub)

    IF (nderiv >= 0) THEN
      MatmuldnMat%d0(:,:) = matmul(Mat1,dnMat2%d0)
    END IF

    IF (nderiv >= 1) THEN
      DO id=1,ndim
        MatmuldnMat%d1(:,:,id) = matmul(Mat1,dnMat2%d1(:,:,id))
      END DO
    END IF

    IF (nderiv >= 2) THEN
      DO id=1,ndim
      DO jd=1,ndim
        MatmuldnMat%d2(:,:,jd,id) = matmul(Mat1,dnMat2%d2(:,:,jd,id))
      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN
      DO id=1,ndim
      DO jd=1,ndim
      DO kd=1,ndim
        MatmuldnMat%d3(:,:,kd,jd,id) = matmul(Mat1,dnMat2%d3(:,:,kd,jd,id))
      END DO
      END DO
      END DO
    END IF

  END FUNCTION QML_MATMUL_Mat1_dnMat2
!> @brief Public subroutine which prints a derived type dnMat.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                TYPE (dnMat_t):      derived type which deals with the derivatives of a matrix.
!! @param nio                integer (optional):  when present unit to print S, otherwise it is the default unit:out_unitp
  SUBROUTINE QML_Write_dnMat(Mat,nio,info)
    USE mod_Lib

    TYPE (dnMat_t),   intent(in)           :: Mat
    integer,          intent(in), optional :: nio
    character(len=*), intent(in), optional :: info

    integer :: i,j,k,nio_loc,nsurf_loc,ndim_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    nsurf_loc  = QML_get_nsurf_FROM_dnMat(Mat)
    ndim_loc   = QML_get_ndim_FROM_dnMat(Mat)

    IF (nsurf_loc == 1 .AND. ndim_loc > 1) THEN
      IF (allocated(Mat%d0)) THEN
        write(nio_loc,'(a,' // RMatIO_format // ')') ' no derivative',Mat%d0
      END IF

      IF (allocated(Mat%d1)) THEN
        write(nio_loc,*) ' 1st derivative'
        CALL Write_RVec(Mat%d1(1,1,:),nio_loc,5)
      END IF

      IF (allocated(Mat%d2)) THEN
        write(nio_loc,*) ' 2d derivative'
        CALL Write_RMat(Mat%d2(1,1,:,:),nio_loc,5)
      END IF
      IF (allocated(Mat%d3)) THEN
        write(nio_loc,*) ' 3d derivative'
        DO i=1,ubound(Mat%d3,dim=5)
        DO j=1,ubound(Mat%d3,dim=4)
        DO k=1,ubound(Mat%d3,dim=3)
          write(nio_loc,'(3(1x,i0)," : ",' // RMatIO_format // ')') k,j,i,Mat%d3(1,1,k,j,i)
        END DO
        END DO
        END DO
      END IF
    ELSE
      IF (allocated(Mat%d0)) THEN
         IF (present(info)) THEN
           write(nio_loc,*) ' no derivative of ',info
         ELSE
           write(nio_loc,*) ' no derivative'
         END IF
        CALL Write_RMat(Mat%d0,nio_loc,5)
      END IF

      IF (allocated(Mat%d1)) THEN
        DO i=1,ubound(Mat%d1,dim=3)
          IF (present(info)) THEN
            write(nio_loc,*) ' 1st derivative of ',info,i
          ELSE
            write(nio_loc,*) ' 1st derivative',i
          END IF
          CALL Write_RMat(Mat%d1(:,:,i),nio_loc,5)
        END DO
      END IF

      IF (allocated(Mat%d2)) THEN
        DO i=1,ubound(Mat%d2,dim=4)
        DO j=1,ubound(Mat%d2,dim=3)
          IF (present(info)) THEN
            write(nio_loc,*) ' 2d derivative of ',info,j,i
          ELSE
            write(nio_loc,*) ' 2d derivative',j,i
          END IF
          CALL Write_RMat(Mat%d2(:,:,j,i),nio_loc,5)
        END DO
        END DO
      END IF

      IF (allocated(Mat%d3)) THEN
        DO i=1,ubound(Mat%d3,dim=5)
        DO j=1,ubound(Mat%d3,dim=4)
        DO k=1,ubound(Mat%d3,dim=3)
          IF (present(info)) THEN
            write(nio_loc,*) ' 3d derivative of ',info,k,j,i
          ELSE
            write(nio_loc,*) ' 3d derivative',k,j,i
          END IF
          CALL Write_RMat(Mat%d3(:,:,k,j,i),nio_loc,5)
        END DO
        END DO
        END DO
      END IF
    END IF

  END SUBROUTINE QML_Write_dnMat
!> @brief Public function to get nderiv from a derived type dnMat.
!!
!> @author David Lauvergnat
!! @date 25/06/2018
!!
!! @param Mat                         TYPE (dnMat_t):     derived type which deals with the derivatives of a matrix (as function of coordinates).
!! @param get_nderiv_FROM_dnMat    integer  (result):  nderiv value, check against Mat%nederiv
  FUNCTION QML_get_nderiv_FROM_dnMat(Mat) RESULT(nderiv)

    integer                       :: nderiv
    TYPE (dnMat_t), intent(in)    :: Mat

    nderiv = Mat%nderiv

    IF (.NOT. allocated(Mat%d0)) THEN
      nderiv = -1
    ELSE IF (.NOT. allocated(Mat%d1)) THEN
      nderiv = 0
    ELSE IF (.NOT. allocated(Mat%d2)) THEN
      nderiv = 1
    ELSE IF (.NOT. allocated(Mat%d3)) THEN
      nderiv = 2
    ELSE
      nderiv = 3
    END IF

    IF (Mat%nderiv /= nderiv) THEN
      write(out_unitp,*) ' ERROR in QML_get_nderiv_FROM_dnMat'
      write(out_unitp,*) '  Problem with nderiv in Mat'
      CALL QML_Write_dnMat(Mat)
      STOP 'ERROR in QML_get_nderiv_FROM_dnMat'
    END IF

    END FUNCTION QML_get_nderiv_FROM_dnMat

!> @brief Public function to get nsurf (the number of electronic surfaces, dimension of the matrix) from a derived type dnMat.
!!
!> @author David Lauvergnat
!! @date 25/06/2018
!!
!! @param Mat                         TYPE (dnMat_t):     derived type which deals with the derivatives of a matrix (as function of coordinates).
!! @param get_nderiv_FROM_dnMat    integer  (result):  nderiv value
  FUNCTION QML_get_nsurf_FROM_dnMat(Mat) RESULT(nsurf)

    integer                       :: nsurf
    TYPE (dnMat_t), intent(in)    :: Mat

    IF (.NOT. allocated(Mat%d0)) THEN
      nsurf = 0
    ELSE
      nsurf = size(Mat%d0(:,1))
    END IF

    END FUNCTION QML_get_nsurf_FROM_dnMat

!> @brief Public function to get ndim (number of coordinates) from a derived type dnMat.
!!
!> @author David Lauvergnat
!! @date 25/06/2018
!!
!! @param Mat                       TYPE (dnMat_t):      derived type which deals with the derivatives of a scalar functions.
!! @param get_ndim_FROM_dnMat    integer  (result):   ndim value from the size of Mat%d1.
  FUNCTION QML_get_ndim_FROM_dnMat(Mat) RESULT(ndim)

    integer                       :: ndim
    TYPE (dnMat_t), intent(in)    :: Mat

    IF (.NOT. allocated(Mat%d1)) THEN
      ndim = 0
    ELSE
      ndim = size(Mat%d1,dim=3)
    END IF

    END FUNCTION QML_get_ndim_FROM_dnMat
!> @brief Public function which ckecks a derived type dnMat is zero (all components).
!!
!> @author David Lauvergnat
!! @date 25/06/2018
!!
!! @param Check_dnMat_IS_ZERO   logical  (result):   result of the comparison
!! @param Mat                      TYPE (dnMat_t):      derived type which deals with the derivatives of a matrix.
!! @param epsi                     real (optional):     when present zero limit, otherwise 10^-10
  FUNCTION QML_Check_dnMat_IS_ZERO(Mat,epsi) RESULT(Check_dnMat_IS_ZERO)
    USE mod_NumParameters

    logical                                  :: Check_dnMat_IS_ZERO
    TYPE (dnMat_t),     intent(in)           :: Mat
    real(kind=Rkind),   intent(in), optional :: epsi


    IF (present(epsi)) THEN
      Check_dnMat_IS_ZERO = QML_get_maxval_OF_dnMat(Mat) <= epsi
    ELSE
      Check_dnMat_IS_ZERO = QML_get_maxval_OF_dnMat(Mat) <= ONETENTH**10
    END IF


    END FUNCTION QML_Check_dnMat_IS_ZERO
!> @brief Public function which gets the largest value of a derived type get_maxval_OF_dnMat (all components).
!!
!> @author David Lauvergnat
!! @date 25/06/2018
!!
!! @param get_maxval_OF_dnMat   real  (result):      largest value (all components)
!! @param Mat                      TYPE (dnMat_t):      derived type which deals with the derivatives of a matrix.
  FUNCTION QML_get_maxval_OF_dnMat(Mat,nderiv) RESULT(get_maxval_OF_dnMat)
    USE mod_NumParameters

    real(kind=Rkind)                     :: get_maxval_OF_dnMat
    TYPE (dnMat_t), intent(in)           :: Mat
    integer,        intent(in), optional :: nderiv

    real(kind=Rkind) :: e0,e1,e2,e3
    integer          :: nderiv_loc

    nderiv_loc = QML_get_nderiv_FROM_dnMat(Mat)
    IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

    IF (allocated(Mat%d0) .AND. nderiv_loc >= 0) THEN
      e0 = maxval(abs(Mat%d0))
    ELSE
      e0 = ZERO
    END IF

    IF (allocated(Mat%d1) .AND. nderiv_loc >= 1) THEN
      e1 = maxval(abs(Mat%d1))
    ELSE
      e1 = ZERO
    END IF

    IF (allocated(Mat%d2) .AND. nderiv_loc >= 2) THEN
      e2 = maxval(abs(Mat%d2))
    ELSE
      e2 = ZERO
    END IF

    IF (allocated(Mat%d3) .AND. nderiv_loc >= 3) THEN
      e3 = maxval(abs(Mat%d3))
    ELSE
      e3 = ZERO
    END IF

    get_maxval_OF_dnMat = max(e0,e1,e2,e3)

    END FUNCTION QML_get_maxval_OF_dnMat
!! @brief Public subroutine which checks if the derived type dnMat is (correctly) allocated.
!!
!> @author David Lauvergnat
!! @date 25/06/2018
!!
!! @param Mat                TYPE (dnMat_t):  derived type which deals with the derivatives of a matrix.
!! @param nderiv             integer:         the derivative order.
  FUNCTION QML_Check_NotAlloc_dnMat(Mat,nderiv) RESULT(NotAlloc)

    logical                       :: NotAlloc
    TYPE (dnMat_t), intent(in)    :: Mat
    integer,        intent(in)    :: nderiv



    NotAlloc = (nderiv >= 0 .AND. .NOT. allocated(Mat%d0))
    NotAlloc = NotAlloc .OR. (nderiv >= 1 .AND. .NOT. allocated(Mat%d1))
    NotAlloc = NotAlloc .OR. (nderiv >= 2 .AND. .NOT. allocated(Mat%d2))
    NotAlloc = NotAlloc .OR. (nderiv >= 3 .AND. .NOT. allocated(Mat%d3))

  END FUNCTION QML_Check_NotAlloc_dnMat

  SUBROUTINE QML_DIAG_dnMat(dnMat,dnMatDiag,dnVec,dnVecProj,dnVec0)
    USE mod_Lib
    USE mod_diago
    IMPLICIT NONE

    TYPE (dnMat_t),     intent(in)              :: dnMat
    TYPE (dnMat_t),     intent(inout)           :: dnMatDiag ! we keep it as a matrix
    TYPE (dnMat_t),     intent(inout)           :: dnVec
    TYPE (dnMat_t),     intent(inout), optional :: dnVecProj
    TYPE (dnMat_t),     intent(inout), optional :: dnVec0

    integer                       :: ndim,nderiv,nsurf
    real(kind=Rkind), allocatable :: Vec(:,:),tVec(:,:),Eig(:),Mtemp(:,:)
    TYPE (dnMat_t)                :: dnMat_OnVec
    integer                       :: i,j,k,id,jd,kd
    real (kind=Rkind)             :: ai,aj,aii,aij,aji,ajj,th,cc,ss

  real (kind=Rkind)               :: epsi = ONETENTH**10


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_DIAG_dnMat'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    nderiv = QML_get_nderiv_FROM_dnMat(dnMat)
    ndim   = QML_get_ndim_FROM_dnMat(dnMat)
    nsurf  = QML_get_nsurf_FROM_dnMat(dnMat)

    CALL QML_dealloc_dnMat(dnMatDiag)
    CALL QML_dealloc_dnMat(dnVec)
    IF (present(dnVecProj)) CALL QML_dealloc_dnMat(dnVecProj)

    IF (nderiv < 0) RETURN

    CALL QML_alloc_dnMat(dnMatDiag,nsurf,ndim,nderiv)
    dnMatDiag = ZERO
    CALL QML_alloc_dnMat(dnVec,nsurf,ndim,nderiv)
    dnVec     = ZERO

    ! the zero order: normal diagonalization
    allocate(Eig(nsurf))
    allocate(Vec(nsurf,nsurf))
    allocate(tVec(nsurf,nsurf))
    allocate(Mtemp(nsurf,nsurf))

    CALL diagonalization(dnMat%d0,Eig,Vec,nsurf,sort=1,phase=.TRUE.)

    IF (present(dnVec0)) THEN
       IF (debug) write(out_unitp,*) 'Change phase?'
       flush(out_unitp)

       DO i=1,nsurf
         IF (dot_product(dnVec0%d0(:,i),Vec(:,i)) < ZERO) Vec(:,i) = -Vec(:,i)
       END DO

       IF (debug) THEN
         write(out_unitp,*) 'Vec before rotation'
         CALL Write_RMat(Vec,nio=out_unitp,nbcol1=5)
       END IF
       !For degenerated eigenvectors (works only with 2 vectors)
       DO i=1,nsurf-1
         IF ( abs(Eig(i)-Eig(i+1)) < epsi) THEN
           j = i+1
           IF (debug) write(out_unitp,*) 'degenerated vectors',i,j

           aii = dot_product(dnVec0%d0(:,i),Vec(:,i))
           aji = dot_product(dnVec0%d0(:,j),Vec(:,i))
           aij = dot_product(dnVec0%d0(:,i),Vec(:,j))
           ajj = dot_product(dnVec0%d0(:,j),Vec(:,j))

           th = ( atan2(aij,ajj) -atan2(aji,aii) ) * HALF
           IF (debug) write(out_unitp,*) 'theta',th

           cc = cos(th)
           ss = sin(th)

           DO k=1,nsurf
             ai = Vec(k,i)
             aj = Vec(k,j)
             Vec(k,i) =  cc * ai + ss * aj
             Vec(k,j) = -ss * ai + cc * aj
           END DO
         END IF
       END DO
    END IF


    tVec         = transpose(Vec)
    dnVec%d0     = matmul(tVec,Vec)                  ! Identity matrix
    dnMatDiag%d0 = matmul(tVec,matmul(dnMat%d0,Vec)) ! diagonal matrix

    dnMat_OnVec = matmul(tVec,matmul(dnMat,Vec)) ! dnMat on the Eigenvector basis
    !  for dnMat_OnVec%d0: Eigenvalues on the diagonal

    IF (nderiv > 0) THEN

      DO id=1,ndim
        Mtemp = dnMat_OnVec%d1(:,:,id)

        DO i=1,nsurf
          ! d1Eig
          dnMatDiag%d1(i,i,id) = Mtemp(i,i)

          ! d1Vec: projection on <i|
          dnVec%d1(i,i,id) = ZERO

          ! d1Vec: projection on <j|
          DO j=1,nsurf
            IF (j == i) CYCLE ! already done
            IF (abs(Eig(i)-Eig(j)) < epsi) CYCLE ! for degenerated eigenvalues

            dnVec%d1(j,i,id) = Mtemp(j,i)/(Eig(i)-Eig(j))

          END DO

        END DO

      END DO


    END IF

    IF (nderiv > 1) THEN

      DO id=1,ndim
      DO jd=1,ndim
        Mtemp = dnMat_OnVec%d2(:,:,jd,id) +                             &
                      matmul(dnMat_OnVec%d1(:,:,id),dnVec%d1(:,:,jd)) + &
                      matmul(dnMat_OnVec%d1(:,:,jd),dnVec%d1(:,:,id))
        DO i=1,nsurf
          Mtemp(:,i) = Mtemp(:,i) -                                     &
                                dnMatDiag%d1(i,i,id)*dnVec%d1(:,i,jd) - &
                                dnMatDiag%d1(i,i,jd)*dnVec%d1(:,i,id)
        END DO

        DO i=1,nsurf
          ! d1Eig
          dnMatDiag%d2(i,i,jd,id) = Mtemp(i,i)

          ! d1Vec: projection on <i|
          dnVec%d2(i,i,jd,id) = -dot_product(dnVec%d1(:,i,id),dnVec%d1(:,i,jd))

          ! d1Vec: projection on <j|
          DO j=1,nsurf
            IF (j == i) CYCLE ! already done
            IF (abs(Eig(i)-Eig(j)) < epsi) CYCLE ! for degenerated eigenvalues

            dnVec%d2(j,i,jd,id) = Mtemp(j,i)/(Eig(i)-Eig(j))

          END DO

        END DO

      END DO
      END DO

    END IF

    IF (nderiv > 2) THEN

      DO id=1,ndim
      DO jd=1,ndim
      DO kd=1,ndim

        Mtemp = dnMat_OnVec%d3(:,:,kd,jd,id) +                          &
                   matmul(dnMat_OnVec%d2(:,:,kd,id),dnVec%d1(:,:,jd)) + &
                   matmul(dnMat_OnVec%d2(:,:,jd,kd),dnVec%d1(:,:,id)) + &
                   matmul(dnMat_OnVec%d2(:,:,id,jd),dnVec%d1(:,:,kd)) + &
                   matmul(dnMat_OnVec%d1(:,:,id),dnVec%d2(:,:,jd,kd)) + &
                   matmul(dnMat_OnVec%d1(:,:,kd),dnVec%d2(:,:,id,jd)) + &
                   matmul(dnMat_OnVec%d1(:,:,jd),dnVec%d2(:,:,kd,id))

        DO i=1,nsurf
          Mtemp(:,i) = Mtemp(:,i) -                                     &
                             dnMatDiag%d2(i,i,kd,id)*dnVec%d1(:,i,jd) - &
                             dnMatDiag%d2(i,i,jd,kd)*dnVec%d1(:,i,id) - &
                             dnMatDiag%d2(i,i,id,jd)*dnVec%d1(:,i,kd) - &
                             dnMatDiag%d1(i,i,id)*dnVec%d2(:,i,jd,kd) - &
                             dnMatDiag%d1(i,i,kd)*dnVec%d2(:,i,id,jd) - &
                             dnMatDiag%d1(i,i,jd)*dnVec%d2(:,i,kd,id)
        END DO

        DO i=1,nsurf
          ! d1Eig
          dnMatDiag%d3(i,i,kd,jd,id) = Mtemp(i,i)

          ! d1Vec: projection on <i|
          dnVec%d3(i,i,kd,jd,id) = - &
               dot_product(dnVec%d1(:,i,kd),dnVec%d2(:,i,jd,id)) - &
               dot_product(dnVec%d1(:,i,jd),dnVec%d2(:,i,id,kd)) - &
               dot_product(dnVec%d1(:,i,id),dnVec%d2(:,i,kd,jd))

          ! d1Vec: projection on <j|
          DO j=1,nsurf
            IF (j == i) CYCLE ! already done
            IF (abs(Eig(i)-Eig(j)) < epsi) CYCLE ! for degenerated eigenvalues

            dnVec%d3(j,i,kd,jd,id) = Mtemp(j,i)/(Eig(i)-Eig(j))

          END DO

        END DO

      END DO
      END DO
      END DO

    END IF


    IF (present(dnVecProj)) dnVecProj = dnVec ! since here dnVec are the projected vectors

    ! unproject the dnVec: correct ???
    dnVec%d0(:,:) = Vec
    IF (nderiv > 0) THEN
      DO id=1,ndim
        dnVec%d1(:,:,id) = matmul(Vec,dnVec%d1(:,:,id))
      END DO
    END IF
    IF (nderiv > 1) THEN
      DO id=1,ndim
      DO jd=1,ndim
        dnVec%d2(:,:,jd,id) = matmul(Vec,dnVec%d2(:,:,jd,id))
      END DO
      END DO
    END IF
    IF (nderiv > 2) THEN
      DO id=1,ndim
      DO jd=1,ndim
      DO kd=1,ndim
        dnVec%d3(:,:,kd,jd,id) = matmul(Vec,dnVec%d3(:,:,kd,jd,id))
      END DO
      END DO
      END DO
    END IF

    IF (allocated(Eig))   deallocate(Eig)
    IF (allocated(Vec))   deallocate(Vec)
    IF (allocated(tVec))  deallocate(tVec)
    IF (allocated(Mtemp)) deallocate(Mtemp)

    CALL QML_dealloc_dnMat(dnMat_OnVec)

    IF (debug) THEN
      IF (present(dnVecProj)) CALL QML_Write_dnMat(dnVecProj,info='dnVecProj')
      CALL QML_Write_dnMat(dnVec,info='dnVec')
      CALL QML_Write_dnMat(dnMatDiag,info='dnMatDiag')
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE QML_DIAG_dnMat

END MODULE mod_dnMat
