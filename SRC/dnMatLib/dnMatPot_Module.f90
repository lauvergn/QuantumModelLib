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
MODULE mod_dnMatPot
  USE mod_NumParameters
  IMPLICIT NONE

  TYPE dnMatPot
     integer                    :: nderiv = -1

     real (kind=Rkind), allocatable :: d0(:,:)
     real (kind=Rkind), allocatable :: d1(:,:,:)
     real (kind=Rkind), allocatable :: d2(:,:,:,:)
  END TYPE dnMatPot


  INTERFACE assignment (=)
     MODULE PROCEDURE sub_dnMatPot2_TO_dnMatPot1,set_dnMatPot_TO_R
  END INTERFACE

   INTERFACE operator (*)
      MODULE PROCEDURE sub_dnMatPot_TIME_R,sub_R_TIME_dnMatPot
   END INTERFACE

   INTERFACE operator (**)
      MODULE PROCEDURE sub_dnMatPot_EXP_R
   END INTERFACE

   INTERFACE operator (+)
      MODULE PROCEDURE dnMatPot2_PLUS_dnMatPot1,sub_dnMatPot_PLUS_R,sub_R_PLUS_dnMatPot
   END INTERFACE

   INTERFACE operator (-)
      MODULE PROCEDURE dnMatPot2_MINUS_dnMatPot1,sub_dnMatPot_MINUS_R,sub_R_MINUS_dnMatPot
   END INTERFACE

CONTAINS

  SUBROUTINE alloc_dnMatPot(Mat,nsurf,ndim,nderiv,name_var,name_sub,err_dnMatPot)
  IMPLICIT NONE

    TYPE(dnMatPot), intent(inout)  :: Mat !< derived type, which contains, matrix potential, its derivatives
    integer, intent(in), optional  :: nsurf !< number of electronic surfaces
    integer, intent(in), optional  :: ndim  !< number of coordinates (for the derivatives)
    integer, intent(in), optional  :: nderiv  !< order of the derivatives [0,1,2]
    character (len=*), intent(in), optional  :: name_var,name_sub

    integer, intent(out), optional :: err_dnMatPot  !< to handle the errors

    ! local variables
    integer :: nsurf_loc,ndim_loc,err_dnMat_loc,nderiv_loc



    err_dnMat_loc = 0 ! no error

    CALL dealloc_dnMatPot(Mat,err_dnMat_loc)
    IF (err_dnMat_loc /= 0) THEN
      write(6,*) ' ERROR in alloc_dnMatPot'
      write(6,*) ' Problem in dealloc_dnMatPot CALL in alloc_dnMatPot'
      IF (present(name_var)) write(6,*) '  for the variable: ',name_var
      IF (present(name_sub)) write(6,*) '  call from the subroutine: ',name_sub
      IF (present(err_dnMatPot)) THEN
        RETURN
      ELSE
        STOP 'Problem in dealloc_dnMatPot CALL in alloc_dnMatPot'
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

    !write(6,*) 'Mat%nderiv in alloc_dnMatPot',Mat%nderiv

    allocate(Mat%d0(nsurf_loc,nsurf_loc),stat=err_dnMat_loc)
    IF (err_dnMat_loc /= 0 .OR. nsurf_loc < 1) THEN
      write(6,*) ' ERROR in alloc_dnMatPot'
      write(6,*) '  Problem with allocate of Mat%d0'
      write(6,*) '  nsurf > 0?',nsurf_loc
      IF (present(name_var)) write(6,*) '  for the variable: ',name_var
      IF (present(name_sub)) write(6,*) '  call from the subroutine: ',name_sub
      IF (present(err_dnMatPot)) THEN
        err_dnMatPot = err_dnMat_loc
        RETURN
      ELSE
        STOP 'Problem with allocate in alloc_dnMatPot'
      END IF
    END IF

    IF (nderiv_loc >= 1) THEN
      allocate(Mat%d1(nsurf_loc,nsurf_loc,ndim_loc),stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0 .OR. nsurf_loc < 1 .OR. ndim_loc < 1) THEN
        write(6,*) ' ERROR in alloc_dnMatPot'
        write(6,*) '  Problem with allocate of Mat%d1'
        write(6,*) '  nsurf > 0?',nsurf_loc
        write(6,*) '  ndim > 0?',ndim_loc
        IF (present(name_var)) write(6,*) '  for the variable: ',name_var
        IF (present(name_sub)) write(6,*) '  call from the subroutine: ',name_sub
        IF (present(err_dnMatPot)) THEN
          err_dnMatPot = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in alloc_dnMatPot'
        END IF
      END IF
    END IF

    IF (nderiv_loc >= 2) THEN
      allocate(Mat%d2(nsurf_loc,nsurf_loc,ndim_loc,ndim_loc),stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0 .OR. nsurf_loc < 1 .OR. ndim_loc < 1) THEN
        write(6,*) ' ERROR in alloc_dnMatPot'
        write(6,*) '  Problem with allocate of Mat%d2'
        write(6,*) '  nsurf > 0',nsurf_loc
        write(6,*) '  ndim > 0',ndim_loc
        IF (present(name_var)) write(6,*) '  for the variable: ',name_var
        IF (present(name_sub)) write(6,*) '  call from the subroutine: ',name_sub
        IF (present(err_dnMatPot)) THEN
          err_dnMatPot = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in alloc_dnMatPot'
        END IF
      END IF
    END IF

  END SUBROUTINE alloc_dnMatPot

  SUBROUTINE dealloc_dnMatPot(Mat,err_dnMatPot)
  IMPLICIT NONE

    TYPE(dnMatPot), intent(inout)     :: Mat !< derived type, which contains, matrix potential, its derivatives
    integer, intent(out), optional :: err_dnMatPot  !< to handle the errors

    ! local variables
    integer :: err_dnMat_loc

    IF (allocated(Mat%d0)) THEN
      deallocate(Mat%d0,stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0) THEN
        write(6,*) ' ERROR in dealloc_dnMatPot'
        write(6,*) '  Problem with deallocate of Mat%d0'
        IF (present(err_dnMatPot)) THEN
          err_dnMatPot = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnMatPot'
        END IF
      END IF
    END IF

    IF (allocated(Mat%d1)) THEN
      deallocate(Mat%d1,stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0) THEN
        write(6,*) ' ERROR in dealloc_dnMatPot'
        write(6,*) '  Problem with deallocate of Mat%d1'
        IF (present(err_dnMatPot)) THEN
          err_dnMatPot = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnMatPot'
        END IF
      END IF
    END IF

    IF (allocated(Mat%d2)) THEN
      deallocate(Mat%d2,stat=err_dnMat_loc)
      IF (err_dnMat_loc /= 0) THEN
        write(6,*) ' ERROR in dealloc_dnMatPot'
        write(6,*) '  Problem with deallocate of Mat%d2'
        IF (present(err_dnMatPot)) THEN
          err_dnMatPot = err_dnMat_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnMatPot'
        END IF
      END IF
    END IF
    Mat%nderiv = -1

  END SUBROUTINE dealloc_dnMatPot

  SUBROUTINE sub_dnMatPot2_TO_dnMatPot1(dnMat1,dnMat2)
    TYPE (dnMatPot), intent(inout) :: dnMat1
    TYPE (dnMatPot), intent(in)    :: dnMat2

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_dnMatPot2_TO_dnMatPot1'

    nderiv_loc = get_nderiv_FROM_dnMatPot(dnMat2)
    nsurf_loc  = get_nsurf_FROM_dnMatPot(dnMat2)
    ndim_loc   = get_ndim_FROM_dnMatPot(dnMat2)

    !write(6,*) 'in ',name_sub,' ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc


    IF (nderiv_loc < 0 .OR. nsurf_loc < 1 .OR. (nderiv_loc > 0  .AND. ndim_loc < 1)) RETURN


    CALL alloc_dnMatPot(dnMat1,nsurf_loc,ndim_loc,nderiv_loc,name_var='dnMat1',name_sub=name_sub)


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
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' nderiv > 2 is NOT possible',nderiv_loc
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF
  END SUBROUTINE sub_dnMatPot2_TO_dnMatPot1

  SUBROUTINE sub_dnSca_TO_dnMatPot(S,Mat,i,j)
    USE mod_dnSca
    TYPE (dnMatPot),    intent(inout) :: Mat
    TYPE (dnSca),         intent(in)    :: S
    integer, optional,  intent(in)    :: i,j

    integer :: nderiv_dnMat,nsurf_dnMat,ndim_dnMat,nderiv_dnSca,ndim_dnSca
    integer :: i_loc,j_loc

    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_dnSca_TO_dnMatPot'


    nderiv_dnSca = get_nderiv_FROM_dnSca(S)
    ndim_dnSca   = get_ndim_FROM_dnSca(S)

    nderiv_dnMat = get_nderiv_FROM_dnMatPot(Mat)
    nsurf_dnMat  = get_nsurf_FROM_dnMatPot(Mat)
    ndim_dnMat   = get_ndim_FROM_dnMatPot(Mat)

    IF ( Check_NotAlloc_dnMatPot(Mat,nderiv_dnSca) .OR.                   &
         nderiv_dnSca /= nderiv_dnMat  .OR.  ndim_dnSca /= ndim_dnMat .OR.  &
         nsurf_dnMat < 1 ) THEN
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' dnMat is not allocated or ...'
      write(6,*) '  ... nderiv from dnMat or dnSca are different or ...'
      write(6,*) '  ... ndim from dnMat or dnSca are different or ...'
      write(6,*) '  ... nsurf from dnMat is < 1'

      write(6,*) 'nderiv from dnMat and dnSca:',nderiv_dnMat,nderiv_dnSca
      write(6,*) 'ndim   from dnMat and dnSca:',ndim_dnMat,ndim_dnSca
      write(6,*) 'nsurf  from dnMat        :',nsurf_dnMat

      write(6,*) 'It should never append! Check the source'
      STOP
    END IF

    i_loc = 1
    j_loc = 1
    IF (present(i)) i_loc = i
    IF (present(j)) j_loc = j


    IF (i_loc < 1 .OR. i_loc > nsurf_dnMat .OR. j_loc < 1 .OR. j_loc > nsurf_dnMat) THEN
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' The matrix indexes, (',i_loc,j_loc,') are out of range [1...',nsurf_dnMat,']'
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF


    ! Potential
    Mat%d0(i_loc,j_loc) = get_d0_FROM_dnSca(S)

    ! gradient
    IF (nderiv_dnSca >= 1) THEN
      Mat%d1(i_loc,j_loc,:) = get_d1_FROM_dnSca(S)
    END IF

    ! Hessian
    IF (nderiv_dnSca >= 2) then
      Mat%d2(i_loc,j_loc,:,:) = get_d2_FROM_dnSca(S)
    END IF


  END SUBROUTINE sub_dnSca_TO_dnMatPot

  SUBROUTINE set_dnMatPot_TO_zero(dnMat)
    TYPE (dnMatPot), intent(inout) :: dnMat

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='set_dnMatPot_TO_zero'


    CALL set_dnMatPot_TO_R(dnMat,ZERO)

  END SUBROUTINE set_dnMatPot_TO_zero
  SUBROUTINE set_dnMatPot_TO_R(dnMat,R)

    TYPE (dnMatPot), intent(inout) :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='set_dnMatPot_TO_R'

    nderiv_loc = get_nderiv_FROM_dnMatPot(dnMat)

    !write(6,*) 'nderiv',nderiv_loc


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
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' nderiv > 2 is NOT possible',nderiv_loc
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF
  END SUBROUTINE set_dnMatPot_TO_R
  FUNCTION sub_dnMatPot_TIME_R(dnMat,R)

    TYPE (dnMatPot)                :: sub_dnMatPot_TIME_R
    TYPE (dnMatPot),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_dnMatPot_TIME_R'

    nderiv_loc = get_nderiv_FROM_dnMatPot(dnMat)
    nsurf_loc  = get_nsurf_FROM_dnMatPot(dnMat)
    ndim_loc   = get_ndim_FROM_dnMatPot(dnMat)

    !write(6,*) 'ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc

    CALL alloc_dnMatPot(sub_dnMatPot_TIME_R,nsurf_loc,ndim_loc,nderiv_loc,name_var='sub_dnMatPot_TIME_R',name_sub=name_sub)

    !write(6,*) 'nderiv',nderiv_loc


    IF (nderiv_loc == 0) THEN
       sub_dnMatPot_TIME_R%d0 = dnMat%d0 * R

    ELSE IF (nderiv_loc == 1) THEN
       sub_dnMatPot_TIME_R%d0 = dnMat%d0 * R
       sub_dnMatPot_TIME_R%d1 = dnMat%d1 * R

    ELSE IF (nderiv_loc == 2) THEN
       sub_dnMatPot_TIME_R%d0 = dnMat%d0 * R
       sub_dnMatPot_TIME_R%d1 = dnMat%d1 * R
       sub_dnMatPot_TIME_R%d2 = dnMat%d2 * R

    ELSE
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' nderiv > 2 is NOT possible',nderiv_loc
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION sub_dnMatPot_TIME_R
  FUNCTION sub_R_TIME_dnMatPot(R,dnMat)

    TYPE (dnMatPot)                :: sub_R_TIME_dnMatPot
    TYPE (dnMatPot),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_R_TIME_dnMatPot'

    nderiv_loc = get_nderiv_FROM_dnMatPot(dnMat)
    nsurf_loc  = get_nsurf_FROM_dnMatPot(dnMat)
    ndim_loc   = get_ndim_FROM_dnMatPot(dnMat)

    !write(6,*) 'ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc

    CALL alloc_dnMatPot(sub_R_TIME_dnMatPot,nsurf_loc,ndim_loc,nderiv_loc,name_var='sub_R_TIME_dnMatPot',name_sub=name_sub)

    !write(6,*) 'nderiv',nderiv_loc


    IF (nderiv_loc == 0) THEN
       sub_R_TIME_dnMatPot%d0 = dnMat%d0 * R

    ELSE IF (nderiv_loc == 1) THEN
       sub_R_TIME_dnMatPot%d0 = dnMat%d0 * R
       sub_R_TIME_dnMatPot%d1 = dnMat%d1 * R

    ELSE IF (nderiv_loc == 2) THEN
       sub_R_TIME_dnMatPot%d0 = dnMat%d0 * R
       sub_R_TIME_dnMatPot%d1 = dnMat%d1 * R
       sub_R_TIME_dnMatPot%d2 = dnMat%d2 * R

    ELSE
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' nderiv > 2 is NOT possible',nderiv_loc
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION sub_R_TIME_dnMatPot
 FUNCTION dnMatPot2_PLUS_dnMatPot1(dnMat1,dnMat2)
    TYPE (dnMatPot)                :: dnMatPot2_PLUS_dnMatPot1
    TYPE (dnMatPot), intent(in)    :: dnMat1,dnMat2

    integer :: nderiv,nsurf,ndim
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='dnMatPot2_PLUS_dnMatPot1'

    nderiv = min(get_nderiv_FROM_dnMatPot(dnMat1),get_nderiv_FROM_dnMatPot(dnMat2))
    nsurf  = min(get_nsurf_FROM_dnMatPot(dnMat1), get_nsurf_FROM_dnMatPot(dnMat2))
    ndim   = min(get_ndim_FROM_dnMatPot(dnMat1),  get_ndim_FROM_dnMatPot(dnMat2))

    !write(6,*) 'in ',name_sub,' nsurf,ndim,nderiv',nsurf,ndim,nderiv

    CALL dealloc_dnMatPot(dnMatPot2_PLUS_dnMatPot1)

    IF (nderiv < 0 .OR. nsurf < 1 .OR. (nderiv > 0  .AND. ndim < 1)) RETURN

    CALL alloc_dnMatPot(dnMatPot2_PLUS_dnMatPot1,nsurf,ndim,nderiv,name_var='dnMatPot2_PLUS_dnMatPot1',name_sub=name_sub)

    IF (nderiv == 0) THEN
       dnMatPot2_PLUS_dnMatPot1%d0 = dnMat1%d0 + dnMat2%d0
    ELSE IF (nderiv == 1) THEN
       dnMatPot2_PLUS_dnMatPot1%d0 = dnMat1%d0 + dnMat2%d0
       dnMatPot2_PLUS_dnMatPot1%d1 = dnMat1%d1 + dnMat2%d1
    ELSE IF (nderiv == 2) THEN
       dnMatPot2_PLUS_dnMatPot1%d0 = dnMat1%d0 + dnMat2%d0
       dnMatPot2_PLUS_dnMatPot1%d1 = dnMat1%d1 + dnMat2%d1
       dnMatPot2_PLUS_dnMatPot1%d2 = dnMat1%d2 + dnMat2%d2
    ELSE
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' nderiv > 2 is NOT possible',nderiv
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION dnMatPot2_PLUS_dnMatPot1
  FUNCTION sub_dnMatPot_PLUS_R(dnMat,R)

    TYPE (dnMatPot)                :: sub_dnMatPot_PLUS_R
    TYPE (dnMatPot),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_dnMatPot_PLUS_R'


    sub_dnMatPot_PLUS_R    = dnMat

    sub_dnMatPot_PLUS_R%d0 = sub_dnMatPot_PLUS_R%d0 + R

    ! the derivatives of R are zero => nothing to be add to %d1 and %d2

  END FUNCTION sub_dnMatPot_PLUS_R
  FUNCTION sub_R_PLUS_dnMatPot(R,dnMat)

    TYPE (dnMatPot)                :: sub_R_PLUS_dnMatPot
    TYPE (dnMatPot),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_R_PLUS_dnMatPot'


    sub_R_PLUS_dnMatPot    = dnMat

    sub_R_PLUS_dnMatPot%d0 = sub_R_PLUS_dnMatPot%d0 + R

    ! the derivatives of R are zero

  END FUNCTION sub_R_PLUS_dnMatPot
 FUNCTION dnMatPot2_MINUS_dnMatPot1(dnMat1,dnMat2)
    TYPE (dnMatPot)                :: dnMatPot2_MINUS_dnMatPot1
    TYPE (dnMatPot), intent(in)    :: dnMat1,dnMat2

    integer :: nderiv,nsurf,ndim
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='dnMatPot2_MINUS_dnMatPot1'

    nderiv = min(get_nderiv_FROM_dnMatPot(dnMat1),get_nderiv_FROM_dnMatPot(dnMat2))
    nsurf  = min(get_nsurf_FROM_dnMatPot(dnMat1), get_nsurf_FROM_dnMatPot(dnMat2))
    ndim   = min(get_ndim_FROM_dnMatPot(dnMat1),  get_ndim_FROM_dnMatPot(dnMat2))


    CALL dealloc_dnMatPot(dnMatPot2_MINUS_dnMatPot1)

    IF (nderiv < 0 .OR. nsurf < 1 .OR. (nderiv > 0  .AND. ndim < 1)) RETURN

    CALL alloc_dnMatPot(dnMatPot2_MINUS_dnMatPot1,nsurf,ndim,nderiv,name_var='dnMatPot2_MINUS_dnMatPot1',name_sub=name_sub)

    IF (nderiv == 0) THEN
       dnMatPot2_MINUS_dnMatPot1%d0 = dnMat1%d0 - dnMat2%d0
    ELSE IF (nderiv == 1) THEN
       dnMatPot2_MINUS_dnMatPot1%d0 = dnMat1%d0 - dnMat2%d0
       dnMatPot2_MINUS_dnMatPot1%d1 = dnMat1%d1 - dnMat2%d1
    ELSE IF (nderiv == 2) THEN
       dnMatPot2_MINUS_dnMatPot1%d0 = dnMat1%d0 - dnMat2%d0
       dnMatPot2_MINUS_dnMatPot1%d1 = dnMat1%d1 - dnMat2%d1
       dnMatPot2_MINUS_dnMatPot1%d2 = dnMat1%d2 - dnMat2%d2
    ELSE
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' nderiv > 2 is NOT possible',nderiv
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION dnMatPot2_MINUS_dnMatPot1
  FUNCTION sub_dnMatPot_MINUS_R(dnMat,R)

    TYPE (dnMatPot)                :: sub_dnMatPot_MINUS_R
    TYPE (dnMatPot),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_dnMatPot_MINUS_R'


    sub_dnMatPot_MINUS_R = dnMat

    sub_dnMatPot_MINUS_R%d0 = dnMat%d0 - R

    ! the derivatives of R are zero

  END FUNCTION sub_dnMatPot_MINUS_R
  FUNCTION sub_R_MINUS_dnMatPot(R,dnMat)

    TYPE (dnMatPot)                :: sub_R_MINUS_dnMatPot
    TYPE (dnMatPot),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_R_MINUS_dnMatPot'

    nderiv_loc = get_nderiv_FROM_dnMatPot(dnMat)
    nsurf_loc  = get_nsurf_FROM_dnMatPot(dnMat)
    ndim_loc   = get_ndim_FROM_dnMatPot(dnMat)

    !write(6,*) 'ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc

    CALL alloc_dnMatPot(sub_R_MINUS_dnMatPot,nsurf_loc,ndim_loc,nderiv_loc,name_var='sub_R_MINUS_dnMatPot',name_sub=name_sub)

    !write(6,*) 'nderiv',nderiv_loc
    IF (nderiv_loc == 0) THEN
       sub_R_MINUS_dnMatPot%d0 = R - dnMat%d0

    ELSE IF (nderiv_loc == 1) THEN
       sub_R_MINUS_dnMatPot%d0 = R - dnMat%d0
       sub_R_MINUS_dnMatPot%d1 =   - dnMat%d1


    ELSE IF (nderiv_loc == 2) THEN
       sub_R_MINUS_dnMatPot%d0 = R - dnMat%d0
       sub_R_MINUS_dnMatPot%d1 =   - dnMat%d1
       sub_R_MINUS_dnMatPot%d2 =   - dnMat%d2

    ELSE
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' nderiv > 2 is NOT possible',nderiv_loc
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF

  END FUNCTION sub_R_MINUS_dnMatPot

  FUNCTION sub_dnMatPot_EXP_R(dnMat,R)

    TYPE (dnMatPot)                :: sub_dnMatPot_EXP_R
    TYPE (dnMatPot),   intent(in)  :: dnMat
    real (kind=Rkind), intent(in)  :: R

    integer :: nderiv_loc,nsurf_loc,ndim_loc,id,jd
    integer :: err_dnMatPot_loc
    character (len=*), parameter :: name_sub='sub_dnMatPot_EXP_R'

    nderiv_loc = get_nderiv_FROM_dnMatPot(dnMat)
    nsurf_loc  = get_nsurf_FROM_dnMatPot(dnMat)
    ndim_loc   = get_ndim_FROM_dnMatPot(dnMat)

    !write(6,*) 'ndim,nsurf,nderiv',ndim_loc,nsurf_loc,nderiv_loc

    CALL alloc_dnMatPot(sub_dnMatPot_EXP_R,nsurf_loc,ndim_loc,nderiv_loc,name_var='sub_dnMatPot_EXP_R',name_sub=name_sub)

    !write(6,*) 'nderiv',nderiv_loc


    IF (nderiv_loc == 0) THEN
       sub_dnMatPot_EXP_R%d0 = dnMat%d0 ** R

    ELSE IF (nderiv_loc == 1) THEN
       sub_dnMatPot_EXP_R%d0 = dnMat%d0 ** R
       DO id=1,ndim_loc
         sub_dnMatPot_EXP_R%d1(:,:,id) = R * dnMat%d0 ** (R-ONE) * dnMat%d1(:,:,id)
       END DO

    ELSE IF (nderiv_loc == 2) THEN
       sub_dnMatPot_EXP_R%d0 = dnMat%d0 ** R
       DO id=1,ndim_loc
         sub_dnMatPot_EXP_R%d1(:,:,id) = R * dnMat%d0 ** (R-ONE) * dnMat%d1(:,:,id)
       END DO
       DO jd=1,ndim_loc
       DO id=1,ndim_loc
         sub_dnMatPot_EXP_R%d2(:,:,jd,id) = R*(R-ONE) * dnMat%d0 ** (R-TWO) * dnMat%d1(:,:,id) * dnMat%d1(:,:,jd) + &
                                            R * dnMat%d0 ** (R-ONE) * dnMat%d2(:,:,jd,id)
       END DO
       END DO

    ELSE
      write(6,*) ' ERROR in ',name_sub
      write(6,*) ' nderiv > 2 is NOT possible',nderiv_loc
      write(6,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION sub_dnMatPot_EXP_R

  SUBROUTINE Write_dnMatPot(Mat,nio)
    USE mod_Lib

    TYPE(dnMatPot), intent(in)    :: Mat
    integer, intent(in), optional :: nio

    integer :: i,j,nio_loc,nsurf_loc,ndim_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = 6
    END IF

    nsurf_loc  = get_nsurf_FROM_dnMatPot(Mat)
    ndim_loc   = get_ndim_FROM_dnMatPot(Mat)

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

  END SUBROUTINE Write_dnMatPot

  FUNCTION get_nderiv_FROM_dnMatPot(Mat)
    integer :: get_nderiv_FROM_dnMatPot

    TYPE(dnMatPot), intent(in)    :: Mat
    integer :: nderiv

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
      write(6,*) ' ERROR in get_nderiv_FROM_dnMatPot'
      write(6,*) '  Problem with nderiv in Mat'
      CALL Write_dnMatPot(Mat)
      STOP 'ERROR in get_nderiv_FROM_dnMatPot'
    END IF

    get_nderiv_FROM_dnMatPot = nderiv

    END FUNCTION get_nderiv_FROM_dnMatPot

  FUNCTION get_nsurf_FROM_dnMatPot(Mat)
    integer :: get_nsurf_FROM_dnMatPot

    TYPE(dnMatPot), intent(in)    :: Mat
    integer :: nsurf

    IF (.NOT. allocated(Mat%d0)) THEN
      get_nsurf_FROM_dnMatPot = 0
    ELSE
      get_nsurf_FROM_dnMatPot = size(Mat%d0(:,1))
    END IF

    END FUNCTION get_nsurf_FROM_dnMatPot

  FUNCTION get_ndim_FROM_dnMatPot(Mat)
    integer :: get_ndim_FROM_dnMatPot

    TYPE(dnMatPot), intent(in)    :: Mat
    integer :: nsurf

    IF (.NOT. allocated(Mat%d1)) THEN
      get_ndim_FROM_dnMatPot = 0
    ELSE
      get_ndim_FROM_dnMatPot = size(Mat%d1(1,1,:))
    END IF

    END FUNCTION get_ndim_FROM_dnMatPot

  FUNCTION Check_dnMatPot_IS_ZERO(Mat,epsi)
    USE mod_NumParameters
    logical :: Check_dnMatPot_IS_ZERO

    TYPE(dnMatPot),     intent(in)           :: Mat
    real(kind=Rkind),   intent(in), optional :: epsi

    real(kind=Rkind) :: epsi_loc,e0,e1,e2


    IF (present(epsi)) THEN
      epsi_loc = epsi
    ELSE
      epsi_loc = ONETENTH**10
    END IF


   IF (.NOT. allocated(Mat%d0)) THEN
      Check_dnMatPot_IS_ZERO = .TRUE.
    ELSE IF (.NOT. allocated(Mat%d1)) THEN
      e0 = maxval(abs(Mat%d0))
      Check_dnMatPot_IS_ZERO = e0 <= epsi_loc
    ELSE IF (.NOT. allocated(Mat%d2)) THEN
      e0 = maxval(abs(Mat%d0))
      e1 = maxval(abs(Mat%d1))
      Check_dnMatPot_IS_ZERO = max(e0,e1) <= epsi_loc
    ELSE
      e0 = maxval(abs(Mat%d0))
      e1 = maxval(abs(Mat%d1))
      e2 = maxval(abs(Mat%d2))
      Check_dnMatPot_IS_ZERO = max(e0,e1,e2) <= epsi_loc
    END IF

    END FUNCTION Check_dnMatPot_IS_ZERO

  FUNCTION get_maxval_OF_dnMatPot(Mat)
    USE mod_NumParameters

    real(kind=Rkind) :: get_maxval_OF_dnMatPot

    TYPE(dnMatPot), intent(in)    :: Mat

    real(kind=Rkind) :: e0,e1,e2



   IF (.NOT. allocated(Mat%d0)) THEN
      get_maxval_OF_dnMatPot = ZERO
    ELSE IF (.NOT. allocated(Mat%d1)) THEN
      e0 = maxval(abs(Mat%d0))
      get_maxval_OF_dnMatPot = e0
    ELSE IF (.NOT. allocated(Mat%d2)) THEN
      e0 = maxval(abs(Mat%d0))
      e1 = maxval(abs(Mat%d1))
      get_maxval_OF_dnMatPot = max(e0,e1)
    ELSE
      e0 = maxval(abs(Mat%d0))
      e1 = maxval(abs(Mat%d1))
      e2 = maxval(abs(Mat%d2))
      get_maxval_OF_dnMatPot = max(e0,e1,e2)
    END IF

    END FUNCTION get_maxval_OF_dnMatPot

  FUNCTION Check_NotAlloc_dnMatPot(Mat,nderiv)

    TYPE(dnMatPot), intent(in)    :: Mat
    integer, intent(in) :: nderiv
    logical :: Check_NotAlloc_dnMatPot


    logical :: NotAlloc

    NotAlloc = (nderiv >= 0 .AND. .NOT. allocated(Mat%d0))
    NotAlloc = NotAlloc .OR. (nderiv >= 1 .AND. .NOT. allocated(Mat%d1))
    NotAlloc = NotAlloc .OR. (nderiv >= 2 .AND. .NOT. allocated(Mat%d2))

    Check_NotAlloc_dnMatPot = NotAlloc

    END FUNCTION Check_NotAlloc_dnMatPot

END MODULE mod_dnMatPot
