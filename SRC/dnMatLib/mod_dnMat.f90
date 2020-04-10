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
!! This module deals with operations or functions of a scalar function and its derivatives, dnS.
!!
!! There is a mapping between the saclar function S, its derivatives and the dnS derived type components:
!!
!! @li M(:,:)                 => M%d0(:,:)
!! @li dM(:,:)/dQ_i           => M%d1(:,:,i)
!! @li d^2M(:,:)/dQ_idQ_j     => M%d2(:,:,i,j)
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

  CONTAINS
    PROCEDURE, PRIVATE :: QML_sub_dnMat2_TO_dnMat1
    PROCEDURE, PRIVATE :: QML_set_dnMat_TO_R
    PROCEDURE, PRIVATE :: QML_set_dnMat_FROM_MatOFdnS
    GENERIC,   PUBLIC  :: assignment(=) => QML_sub_dnMat2_TO_dnMat1,    &
                          QML_set_dnMat_TO_R,QML_set_dnMat_FROM_MatOFdnS
  END TYPE dnMat_t

  PRIVATE :: QML_sub_dnMat2_TO_dnMat1,QML_set_dnMat_TO_R,QML_set_dnMat_FROM_MatOFdnS

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

    TYPE (dnMat_t), intent(inout)  :: Mat   !< derived type, which contains, matrix potential, its derivatives
    integer, intent(in), optional  :: nsurf !< number of electronic surfaces
    integer, intent(in), optional  :: ndim  !< number of coordinates (for the derivatives)
    integer, intent(in), optional  :: nderiv  !< order of the derivatives [0,1,2]
    character (len=*), intent(in), optional  :: name_var,name_sub

    integer, intent(out), optional :: err_dnMat  !< to handle the errors

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
      nderiv_loc = min(2,nderiv_loc)
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

    TYPE (dnMat_t), intent(inout)     :: Mat !< derived type, which contains, matrix potential, its derivatives
    integer, intent(out), optional :: err_dnMat  !< to handle the errors

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
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 is NOT possible',nderiv_loc
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
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 or nderiv < 0 is NOT possible',nderiv_loc
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

    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 is NOT possible',nderiv_loc
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

    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 is NOT possible',nderiv_loc
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
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 is NOT possible',nderiv
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
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 is NOT possible',nderiv
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

    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 is NOT possible',nderiv_loc
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

    integer :: nderiv_loc,nsurf_loc,ndim_loc,id,jd
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
       DO jd=1,ndim_loc
       DO id=1,ndim_loc
         sub_dnMat_EXP_R%d2(:,:,jd,id) = R*(R-ONE) * dnMat%d0 ** (R-TWO) * dnMat%d1(:,:,id) * dnMat%d1(:,:,jd) + &
                                            R * dnMat%d0 ** (R-ONE) * dnMat%d2(:,:,jd,id)
       END DO
       END DO

    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION QML_sub_dnMat_EXP_R
!> @brief Public subroutine which prints a derived type dnMat.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Mat                TYPE (dnMat_t):      derived type which deals with the derivatives of a matrix.
!! @param nio                integer (optional):  when present unit to print S, otherwise it is the default unit:out_unitp
  SUBROUTINE QML_Write_dnMat(Mat,nio)
    USE mod_Lib

    TYPE (dnMat_t), intent(in)    :: Mat
    integer, intent(in), optional :: nio

    integer :: i,j,nio_loc,nsurf_loc,ndim_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = 6
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
    ELSE
      IF (allocated(Mat%d0)) THEN
         write(nio_loc,*) ' no derivative'
        CALL Write_RMat(Mat%d0,nio_loc,5)
      END IF

      IF (allocated(Mat%d1)) THEN
        DO i=1,ubound(Mat%d1,dim=3)
          write(nio_loc,*) ' 1st derivative',i
          CALL Write_RMat(Mat%d1(:,:,i),nio_loc,5)
        END DO
      END IF

      IF (allocated(Mat%d2)) THEN
        DO i=1,ubound(Mat%d2,dim=4)
        DO j=1,ubound(Mat%d2,dim=3)
          write(nio_loc,*) ' 2d derivative',i,j
          CALL Write_RMat(Mat%d2(:,:,j,i),nio_loc,5)
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
    ELSE
      nderiv = 2
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
      ndim = size(Mat%d1(1,1,:))
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

    real(kind=Rkind) :: epsi_loc,e0,e1,e2


    IF (present(epsi)) THEN
      epsi_loc = epsi
    ELSE
      epsi_loc = ONETENTH**10
    END IF


   IF (.NOT. allocated(Mat%d0)) THEN
      Check_dnMat_IS_ZERO = .TRUE.
    ELSE IF (.NOT. allocated(Mat%d1)) THEN
      e0 = maxval(abs(Mat%d0))
      Check_dnMat_IS_ZERO = e0 <= epsi_loc
    ELSE IF (.NOT. allocated(Mat%d2)) THEN
      e0 = maxval(abs(Mat%d0))
      e1 = maxval(abs(Mat%d1))
      Check_dnMat_IS_ZERO = max(e0,e1) <= epsi_loc
    ELSE
      e0 = maxval(abs(Mat%d0))
      e1 = maxval(abs(Mat%d1))
      e2 = maxval(abs(Mat%d2))
      Check_dnMat_IS_ZERO = max(e0,e1,e2) <= epsi_loc
    END IF

    END FUNCTION QML_Check_dnMat_IS_ZERO
!> @brief Public function which gets the largest value of a derived type get_maxval_OF_dnMat (all components).
!!
!> @author David Lauvergnat
!! @date 25/06/2018
!!
!! @param get_maxval_OF_dnMat   real  (result):      largest value (all components)
!! @param Mat                      TYPE (dnMat_t):      derived type which deals with the derivatives of a matrix.
  FUNCTION QML_get_maxval_OF_dnMat(Mat) RESULT(get_maxval_OF_dnMat)
    USE mod_NumParameters

    real(kind=Rkind)              :: get_maxval_OF_dnMat
    TYPE (dnMat_t), intent(in)    :: Mat

    real(kind=Rkind) :: e0,e1,e2

   IF (.NOT. allocated(Mat%d0)) THEN
      get_maxval_OF_dnMat = ZERO
    ELSE IF (.NOT. allocated(Mat%d1)) THEN
      e0 = maxval(abs(Mat%d0))
      get_maxval_OF_dnMat = e0
    ELSE IF (.NOT. allocated(Mat%d2)) THEN
      e0 = maxval(abs(Mat%d0))
      e1 = maxval(abs(Mat%d1))
      get_maxval_OF_dnMat = max(e0,e1)
    ELSE
      e0 = maxval(abs(Mat%d0))
      e1 = maxval(abs(Mat%d1))
      e2 = maxval(abs(Mat%d2))
      get_maxval_OF_dnMat = max(e0,e1,e2)
    END IF

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

  END FUNCTION QML_Check_NotAlloc_dnMat

END MODULE mod_dnMat
