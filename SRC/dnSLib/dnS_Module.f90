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
!===========================================================================
!===========================================================================
!> @brief Module which deals with derivatives of a scalar functions.
!!
!! This module deals with operations or functions of a scalar function and its derivatives, dnSca.
!!
!! There is a mapping between the saclar function S, its derivatives and the dnSca derived type components:
!!
!! @li S                 => S%d0
!! @li dS/dQ_i           => S%d1(i)
!! @li d^2S/dQ_idQ_j     => S%d2(i,j)
!! @li d^3S/dQ_idQ_jdQ_k => S%d3(i,j,k)
!!
!! with S defined as:
!!  TYPE(dnSca) :: S
!!
!!
!! All standard fortran operators (= + - * / **) are overloaded:
!!
!! For instance the sum (+) of two dnSca variables, S1 and S2 corresponds to:
!! @li (S1+S2)                 => S1%d0    + S2%d0
!! @li d(S1+S2)/dQ_i           => S1%d1(i) + S2%d1(i)
!! @li ....
!!
!! The product (*) of two dnSca variables, S1 and S2 correspond to:
!! @li (S1*S2)                 => S1%d0 * S1%d0
!! @li d(S1*S2)/dQ_i           => S1%d0 * S2%d1(i) + S1%d1(i) * S2%d0    (derivative of a product)
!! @li ....
!!
!! All standard fortran functions (exp, sqrt, log ... sin ... sinh) are overloaded
!!
!! For instance the function, f, of dnSca variables, S, corresponds to:
!! @li f(S)                    =>           f(S%d0)
!! @li d(f(S))/dQ_i            => S%d1(i) * f'(S%d0)
!! @li ....
!!
!! All fortran comparison operators (== /= > >= < <=) and (.EQ. .LT. ....) are overloaded, as well.
!! The comparison are done on the zero-order component:
!! S1 == S2                   =>   S1%d0 == S2%d0
!! S1 > S2                    =>   S1%d0  > S2%d0
!!
!!
!! @author David Lauvergnat
!! @date 09/08/2017
!!
MODULE mod_dnSca
  USE mod_NumParameters
  IMPLICIT NONE

!> @brief Derived type which deals with the derivatives of a scalar functions. It is a semi-private type.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param nderiv                  integer: it enables to chose the derivative order (from 0 to 3)
!! @param d0                      real:    0 order derivative (no derivative)
!! @param d1                      real:    1st order derivative (gradient: table of ndim derivatives)
!! @param d2                      real:    2d  order derivative (hessian: matrix of ndim*ndim derivatives)
!! @param d3                      real:    3d  order derivative (ndim*ndim*ndim derivatives)
  TYPE, PUBLIC :: dnSca
     PRIVATE
     integer                        :: nderiv = -1

     real (kind=Rkind)              :: d0
     real (kind=Rkind), allocatable :: d1(:)
     real (kind=Rkind), allocatable :: d2(:,:)
     real (kind=Rkind), allocatable :: d3(:,:,:)
  END TYPE dnSca

  PRIVATE sub_dnSca2_TO_dnSca1,set_dnSca_TO_R
  INTERFACE assignment (=)
     MODULE PROCEDURE sub_dnSca2_TO_dnSca1,set_dnSca_TO_R
  END INTERFACE
  PRIVATE dnSca2_PLUS_dnSca1,dnSca_PLUS_R,R_PLUS_dnSca,PLUS_dnSca
  INTERFACE operator (+)
     MODULE PROCEDURE dnSca2_PLUS_dnSca1,dnSca_PLUS_R,R_PLUS_dnSca,PLUS_dnSca
  END INTERFACE
  PRIVATE dnSca2_MINUS_dnSca1,dnSca_MINUS_R,R_MINUS_dnSca,MINUS_dnSca
  INTERFACE operator (-)
     MODULE PROCEDURE dnSca2_MINUS_dnSca1,dnSca_MINUS_R,R_MINUS_dnSca,MINUS_dnSca
  END INTERFACE
  PRIVATE dnSca2_TIME_dnSca1,dnSca_TIME_R,R_TIME_dnSca
  INTERFACE operator (*)
     MODULE PROCEDURE dnSca2_TIME_dnSca1,dnSca_TIME_R,R_TIME_dnSca
  END INTERFACE
  PRIVATE dnSca2_OVER_dnSca1,dnSca_OVER_R,R_OVER_dnSca
  INTERFACE operator (/)
     MODULE PROCEDURE dnSca2_OVER_dnSca1,dnSca_OVER_R,R_OVER_dnSca
  END INTERFACE
  PRIVATE dnSca_EXP_R,dnSca_EXP_I
  INTERFACE operator (**)
     MODULE PROCEDURE dnSca_EXP_R,dnSca_EXP_I
  END INTERFACE
  PRIVATE dnSca_EQ_dnSca,dnSca_EQ_R,R_EQ_dnSca
  INTERFACE operator (==)
     MODULE PROCEDURE dnSca_EQ_dnSca,dnSca_EQ_R,R_EQ_dnSca
  END INTERFACE
  PRIVATE dnSca_NEQ_dnSca,dnSca_NEQ_R,R_NEQ_dnSca
  INTERFACE operator (/=)
     MODULE PROCEDURE dnSca_NEQ_dnSca,dnSca_NEQ_R,R_NEQ_dnSca
  END INTERFACE
  PRIVATE dnSca_LE_dnSca,dnSca_LE_R,R_LE_dnSca
  INTERFACE operator (<=)
     MODULE PROCEDURE dnSca_LE_dnSca,dnSca_LE_R,R_LE_dnSca
  END INTERFACE
  PRIVATE dnSca_LT_dnSca,dnSca_LT_R,R_LT_dnSca
  INTERFACE operator (<)
     MODULE PROCEDURE dnSca_LT_dnSca,dnSca_LT_R,R_LT_dnSca
  END INTERFACE
  PRIVATE dnSca_GE_dnSca,dnSca_GE_R,R_GE_dnSca
  INTERFACE operator (>=)
     MODULE PROCEDURE dnSca_GE_dnSca,dnSca_GE_R,R_GE_dnSca
  END INTERFACE
  PRIVATE dnSca_GT_dnSca,dnSca_GT_R,R_GT_dnSca
  INTERFACE operator (>)
     MODULE PROCEDURE dnSca_GT_dnSca,dnSca_GT_R,R_GT_dnSca
  END INTERFACE

  PRIVATE get_SQRT_dnSca,get_ABS_dnSca,get_EXP_dnSca,get_LOG_dnSca,get_LOG10_dnSca
  INTERFACE sqrt
     MODULE PROCEDURE get_SQRT_dnSca
  END INTERFACE
  INTERFACE abs
     MODULE PROCEDURE get_ABS_dnSca
  END INTERFACE
  INTERFACE exp
     MODULE PROCEDURE get_EXP_dnSca
  END INTERFACE
  INTERFACE log
     MODULE PROCEDURE get_LOG_dnSca
  END INTERFACE
  INTERFACE log10
     MODULE PROCEDURE get_LOG10_dnSca
  END INTERFACE
  PRIVATE get_COS_dnSca,get_ACOS_dnSca,get_SIN_dnSca,get_ASIN_dnSca,get_TAN_dnSca,get_ATAN_dnSca
  INTERFACE cos
     MODULE PROCEDURE get_COS_dnSca
  END INTERFACE
  INTERFACE acos
     MODULE PROCEDURE get_ACOS_dnSca
  END INTERFACE
  INTERFACE sin
     MODULE PROCEDURE get_SIN_dnSca
  END INTERFACE
  INTERFACE asin
     MODULE PROCEDURE get_ASIN_dnSca
  END INTERFACE
  INTERFACE tan
     MODULE PROCEDURE get_TAN_dnSca
  END INTERFACE
  INTERFACE atan
     MODULE PROCEDURE get_ATAN_dnSca
  END INTERFACE
  PRIVATE get_COSH_dnSca,get_ACOSH_dnSca,get_SINH_dnSca,get_ASINH_dnSca,get_TANH_dnSca,get_ATANH_dnSca
  INTERFACE cosh
     MODULE PROCEDURE get_COSH_dnSca
  END INTERFACE
  INTERFACE acosh
     MODULE PROCEDURE get_ACOSH_dnSca
  END INTERFACE
  INTERFACE sinh
     MODULE PROCEDURE get_SINH_dnSca
  END INTERFACE
  INTERFACE asinh
     MODULE PROCEDURE get_ASINH_dnSca
  END INTERFACE
  INTERFACE tanh
     MODULE PROCEDURE get_TANH_dnSca
  END INTERFACE
  INTERFACE atanh
     MODULE PROCEDURE get_ATANH_dnSca
  END INTERFACE
CONTAINS
!> @brief Public subroutine which allocates a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
!! @param ndim               integer (optional):    number of variables (coordiantes) for the derivatives.
!! @param nderiv             integer (optional):    it enables to chose the derivative order (from 0 to 3).
!! @param err_dnSca            integer (optional):  to handle the errors errors (0: no error).
  SUBROUTINE alloc_dnSca(S,ndim,nderiv,err_dnSca)
    TYPE(dnSca),   intent(inout)          :: S       !< derived type, which contains, matrix potential, its derivatives
    integer,     intent(in),  optional  :: ndim      !< number of variables (coordiantes)
    integer,     intent(in),  optional  :: nderiv    !< order of the derivatives [0,1,3]
    integer,     intent(out), optional  :: err_dnSca !< to handle the errors

    ! local variables
    integer :: ndim_loc,err_dnSca_loc,nderiv_loc

    err_dnSca_loc = 0 ! no error
    IF (present(err_dnSca)) err_dnSca = 0

    CALL dealloc_dnSca(S,err_dnSca_loc)
    IF (err_dnSca_loc /= 0) THEN
      write(out_unitp,*) ' ERROR in alloc_dnSca'
      write(out_unitp,*) ' Problem in dealloc_dnSca CALL in alloc_dnSca'
      IF (present(err_dnSca)) THEN
        err_dnSca = err_dnSca_loc
        RETURN
      ELSE
        STOP 'Problem in dealloc_dnSca CALL in alloc_dnSca'
      END IF
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
    S%nderiv = nderiv_loc

    !write(out_unitp,*) 'S%nderiv in alloc_dnSca',S%nderiv


    IF (nderiv_loc >= 1) THEN
      allocate(S%d1(ndim_loc),stat=err_dnSca_loc)
      IF (err_dnSca_loc /= 0 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in alloc_dnSca'
        write(out_unitp,*) '  Problem with allocate of S%d1'
        write(out_unitp,*) '  ndim > 0?',ndim_loc
        IF (present(err_dnSca)) THEN
          err_dnSca = err_dnSca_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in alloc_dnSca'
        END IF
      END IF
    END IF

    IF (nderiv_loc >= 2) THEN
      allocate(S%d2(ndim_loc,ndim_loc),stat=err_dnSca_loc)
      IF (err_dnSca_loc /= 0 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in alloc_dnSca'
        write(out_unitp,*) '  Problem with allocate of S%d2'
        write(out_unitp,*) '  ndim > 0',ndim_loc
        IF (present(err_dnSca)) THEN
          err_dnSca = err_dnSca_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in alloc_dnSca'
        END IF
      END IF
    END IF

    IF (nderiv_loc >= 3) THEN
      allocate(S%d3(ndim_loc,ndim_loc,ndim_loc),stat=err_dnSca_loc)
      IF (err_dnSca_loc /= 0 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in alloc_dnSca'
        write(out_unitp,*) '  Problem with allocate of S%d3'
        write(out_unitp,*) '  ndim > 0',ndim_loc
        IF (present(err_dnSca)) THEN
          err_dnSca = err_dnSca_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in alloc_dnSca'
        END IF
      END IF
    END IF

    !write(out_unitp,*) 'err_dnSca_loc',err_dnSca_loc
    !IF (present(err_dnSca)) write(out_unitp,*) 'err_dnSca',err_dnSca

  END SUBROUTINE alloc_dnSca

!> @brief Public subroutine which deallocates a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnSca):             derived type which deals with the derivatives of a scalar functions.
!! @param err_dnSca            integer (optional):    to handle the errors errors (0: no error).
  SUBROUTINE dealloc_dnSca(S,err_dnSca)
    TYPE(dnSca), intent(inout)       :: S !< derived type, which contains, matrix potential, its derivatives
    integer, intent(out), optional :: err_dnSca  !< to handle the errors

    ! local variables
    integer :: err_dnSca_loc

    IF (allocated(S%d1)) THEN
      deallocate(S%d1,stat=err_dnSca_loc)
      IF (err_dnSca_loc /= 0) THEN
        write(out_unitp,*) ' ERROR in dealloc_dnSca'
        write(out_unitp,*) '  Problem with deallocate of S%d1'
        IF (present(err_dnSca)) THEN
          err_dnSca = err_dnSca_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnSca'
        END IF
      END IF
    END IF

    IF (allocated(S%d2)) THEN
      deallocate(S%d2,stat=err_dnSca_loc)
      IF (err_dnSca_loc /= 0) THEN
        write(out_unitp,*) ' ERROR in dealloc_dnSca'
        write(out_unitp,*) '  Problem with deallocate of S%d2'
        IF (present(err_dnSca)) THEN
          err_dnSca = err_dnSca_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnSca'
        END IF
      END IF
    END IF

    IF (allocated(S%d3)) THEN
      deallocate(S%d3,stat=err_dnSca_loc)
      IF (err_dnSca_loc /= 0) THEN
        write(out_unitp,*) ' ERROR in dealloc_dnSca'
        write(out_unitp,*) '  Problem with deallocate of S%d3'
        IF (present(err_dnSca)) THEN
          err_dnSca = err_dnSca_loc
          RETURN
        ELSE
          STOP 'Problem with deallocate in dealloc_dnSca'
        END IF
      END IF
    END IF

    S%nderiv = -1

  END SUBROUTINE dealloc_dnSca
!> @brief Public subroutine which checks if the derived type dnSca is (correctly) allocated.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnSca):  derived type which deals with the derivatives of a scalar functions.
!! @param nderiv             integer:    the derivative order.
  FUNCTION Check_NotAlloc_dnSca(S,nderiv)

    logical :: Check_NotAlloc_dnSca
    TYPE(dnSca), intent(in)    :: S
    integer,   intent(in)    :: nderiv


    logical :: NotAlloc

    NotAlloc = (nderiv < 0)
    NotAlloc = NotAlloc .OR. (nderiv >= 1 .AND. .NOT. allocated(S%d1))
    NotAlloc = NotAlloc .OR. (nderiv >= 2 .AND. .NOT. allocated(S%d2))
    NotAlloc = NotAlloc .OR. (nderiv >= 3 .AND. .NOT. allocated(S%d3))

    Check_NotAlloc_dnSca = NotAlloc

  END FUNCTION Check_NotAlloc_dnSca
!> @brief Public function which initializes a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Sres               TYPE(dnSca) (result):    derived type which deals with the derivatives of a scalar functions.
!! @param R                  real:                  Value of Sres%d0.
!! @param ndim               integer (optional):    number of variables (coordiantes) for the derivatives.
!! @param nderiv             integer (optional):    it enables to chose the derivative order (from 0 to 3).
!! @param iQ                 integer (optional):    when ndim > 1, dSres/dQ_iQ = Sres%d1(iQ)= 1, the other derivatives are zero
  FUNCTION init_dnSca(R,ndim,nderiv,iQ) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    real (kind=Rkind), intent(in)    :: R
    integer, optional, intent(in)    :: nderiv,ndim,iQ

    integer :: err_dnSca_loc,nderiv_loc,ndim_loc,iQ_loc
    character (len=*), parameter :: name_sub='init_dnSca'

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

    ! test nderiv
    IF (present(iQ)) THEN
      iQ_loc = iQ
    ELSE
      iQ_loc = 1
    END IF

    IF (iQ_loc < 1 .OR. iQ_loc > ndim_loc) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' the iQ value (',iQ_loc,') is out of range: [1:',ndim_loc,']'
      STOP 'Problem with the iQ value in init_dnSca'
    END IF

    !write(out_unitp,*) 'iQ_loc',iQ_loc ; flush(out_unitp)
    !write(out_unitp,*) 'ndim_loc,nderiv_loc',ndim_loc,nderiv_loc ; flush(out_unitp)

    CALL alloc_dnSca(Sres,ndim_loc,nderiv_loc,err_dnSca_loc)
    IF (err_dnSca_loc /= 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Problem in alloc_dnSca CALL'
      STOP 'Problem Problem in alloc_dnSca CALL in init_dnSca'
    END IF

    Sres = ZERO

    Sres%d0 = R
    IF (nderiv_loc > 0) Sres%d1(iQ_loc) = ONE

  END FUNCTION init_dnSca
!> @brief Public subroutine which initializes a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
!! @param d0,d1,d2,d3        real  (optional):    real value (d0) or table to initialize S.
  SUBROUTINE set_dnSca(S,d0,d1,d2,d3)
    real (kind=Rkind), optional,   intent(in)     :: d0
    real (kind=Rkind), optional,   intent(in)     :: d1(:)
    real (kind=Rkind), optional,   intent(in)     :: d2(:,:)
    real (kind=Rkind), optional,   intent(in)     :: d3(:,:,:)

    TYPE (dnSca),                    intent(inout)  :: S


    character (len=*), parameter :: name_sub='set_dnSca'

    IF (present(d0)) THEN
      S%d0 = d0
      S%nderiv = 0
    END IF

    IF (present(d1)) THEN
      S%d1     = d1
      S%nderiv = 1
      IF (.NOT. present(d0)) THEN
        write(out_unitp,*) ' ERROR in set_dnSca'
        write(out_unitp,*) ' d1 is present but not d0'
        write(out_unitp,*) ' CHECK the fortran!!'
        STOP 'ERROR in set_dnSca'
      END IF
    END IF

    IF (present(d2)) THEN
      S%d2     = d2
      S%nderiv = 2
      IF (.NOT. present(d1)) THEN
        write(out_unitp,*) ' ERROR in set_dnSca'
        write(out_unitp,*) ' d2 is present but not d1'
        write(out_unitp,*) ' CHECK the fortran!!'
        STOP 'ERROR in set_dnSca'
      END IF
    END IF

    IF (present(d3)) THEN
      S%d3     = d3
      S%nderiv = 3
      IF (.NOT. present(d2)) THEN
        write(out_unitp,*) ' ERROR in set_dnSca'
        write(out_unitp,*) ' d3 is present but not d2'
        write(out_unitp,*) ' CHECK the fortran!!'
        STOP 'ERROR in set_dnSca'
      END IF
    END IF

  END SUBROUTINE set_dnSca
!> @brief Public function to get d0 from a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
!! @param d0                 real  (result):      d0=S%d0
  FUNCTION get_d0_FROM_dnSca(S) RESULT(d0)
    real (kind=Rkind)                :: d0
    TYPE (dnSca),        intent(in)    :: S

    character (len=*), parameter :: name_sub='get_d0_FROM_dnSca'

    d0 = S%d0

  END FUNCTION get_d0_FROM_dnSca
!> @brief Public function to get d1 from a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
!! @param d1                 real  (result):        d1=S%d1, d1 is allocatable
  FUNCTION get_d1_FROM_dnSca(S) RESULT(d1)
    real (kind=Rkind), allocatable   :: d1(:)
    TYPE (dnSca),        intent(in)    :: S

    character (len=*), parameter :: name_sub='get_d1_FROM_dnSca'

    !IF (allocated(d1)) deallocate(d1)
    IF (allocated(S%d1)) d1 = S%d1

  END FUNCTION get_d1_FROM_dnSca
  !> @brief Public function to get d1 from a derived type dnSca. Similar to get_d1_FROM_dnSca, but d1 is NOT allocatable.
  !!
  !> @author David Lauvergnat
  !! @date 25/05/2018
  !!
  !! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
  !! @param d1                 real       :           d1=S%d1, d1 is NOT allocatable
    SUBROUTINE sub_get_d1_FROM_dnSca(d1,S)
      real (kind=Rkind),   intent(inout) :: d1(:)
      TYPE (dnSca),        intent(in)    :: S

      character (len=*), parameter :: name_sub='sub_get_d1_FROM_dnSca'

      IF (allocated(S%d1)) d1(:) = S%d1

    END SUBROUTINE sub_get_d1_FROM_dnSca
  !> @brief Public function to get d2 from a derived type dnSca.
  !!
  !> @author David Lauvergnat
  !! @date 03/08/2017
  !!
  !! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
  !! @param d2                 real  (result):        d2=S%d2, d2 is allocatable
    FUNCTION get_d2_FROM_dnSca(S) RESULT(d2)
      real (kind=Rkind), allocatable   :: d2(:,:)
      TYPE (dnSca),        intent(in)    :: S

      character (len=*), parameter :: name_sub='get_d2_FROM_dnSca'

      !IF (allocated(d2)) deallocate(d2)
      IF (allocated(S%d2)) d2 = S%d2

    END FUNCTION get_d2_FROM_dnSca
  !> @brief Public function to get d2 from a derived type dnSca. Similar to get_d2_FROM_dnSca, but d2 is NOT allocatable.
  !!
  !> @author David Lauvergnat
  !! @date 25/05/2018
  !!
  !! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
  !! @param d2                 real:                  d2=S%d2, d2 is NOT allocatable
    SUBROUTINE sub_get_d2_FROM_dnSca(d2,S)
      real (kind=Rkind),   intent(inout) :: d2(:,:)
      TYPE (dnSca),        intent(in)    :: S

      character (len=*), parameter :: name_sub='sub_get_d2_FROM_dnSca'

      IF (allocated(S%d2)) d2(:,:) = S%d2

    END SUBROUTINE sub_get_d2_FROM_dnSca
  !> @brief Public function to get d3 from a derived type dnSca.
  !!
  !> @author David Lauvergnat
  !! @date 03/08/2017
  !!
  !! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
  !! @param d3                 real  (result):        d3=S%d3, d3 is allocatable
    FUNCTION get_d3_FROM_dnSca(S) RESULT(d3)
      real (kind=Rkind), allocatable   :: d3(:,:,:)
      TYPE (dnSca),        intent(in)    :: S

      character (len=*), parameter :: name_sub='get_d3_FROM_dnSca'

      !IF (allocated(d3)) deallocate(d3)
      IF (allocated(S%d3)) d3 = S%d3

    END FUNCTION get_d3_FROM_dnSca
  !> @brief Public function to get d3 from a derived type dnSca. Similar to get_d3_FROM_dnSca, but d3 is NOT allocatable.
  !!
  !> @author David Lauvergnat
  !! @date 25/05/2018
  !!
  !! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
  !! @param d3                 real:                  d3=S%d3, d3 is NOT allocatable
    SUBROUTINE sub_get_d3_FROM_dnSca(d3,S)
      real (kind=Rkind),   intent(inout) :: d3(:,:,:)
      TYPE (dnSca),        intent(in)    :: S

      character (len=*), parameter :: name_sub='sub_get_d3_FROM_dnSca'

      IF (allocated(S%d3)) d3(:,:,:) = S%d3

    END SUBROUTINE sub_get_d3_FROM_dnSca
!> @brief Public subroutine which prints a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
!! @param nio                integer (optional):  when present unit to print S, otherwise it is the default unit:out_unitp
  SUBROUTINE Write_dnSca(S,nio,info)
    USE mod_Lib

    TYPE(dnSca),        intent(in)           :: S
    integer,          intent(in), optional :: nio
    character(len=*), intent(in), optional :: info

    integer :: i,j,k,nio_loc,ndim
    character (len=50) :: fformat

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    IF (present(info)) write(nio_loc,*) info

    ndim = 0
    IF (allocated(S%d1)) ndim = size(S%d1)

    write(nio_loc,'(a,3(3x),x,sp,e12.3)') ' 0   derivative',S%d0

    IF (allocated(S%d1)) THEN
      IF (ndim > 99) THEN
        fformat = '(a,(x,i0),x,sp,e12.3)'
      ELSE
        fformat = '(a,(x,i2),2(3x),x,sp,e12.3)'
      END IF
      DO i=1,ubound(S%d1,dim=1)
        write(nio_loc,fformat) ' 1st derivative',i,S%d1(i)
      END DO
    END IF

    IF (allocated(S%d2)) THEN
      IF (ndim > 99) THEN
        fformat = '(a,2(x,i0),x,sp,e12.3)'
      ELSE
        fformat = '(a,2(x,i2),1(3x),x,sp,e12.3)'
      END IF
      DO i=1,ubound(S%d2,dim=2)
      DO j=1,ubound(S%d2,dim=1)
        write(nio_loc,fformat) ' 2d  derivative',i,j,S%d2(j,i)
      END DO
      END DO
    END IF

    IF (allocated(S%d3)) THEN
      IF (ndim > 99) THEN
        fformat = '(a,3(x,i0),x,sp,e12.3)'
      ELSE
        fformat = '(a,3(x,i2),x,sp,e12.3)'
      END IF
      DO i=1,ubound(S%d3,dim=3)
      DO j=1,ubound(S%d3,dim=2)
      DO k=1,ubound(S%d3,dim=1)
        write(nio_loc,fformat) ' 3d  derivative',i,j,k,S%d3(k,j,i)
      END DO
      END DO
      END DO
    END IF

  END SUBROUTINE Write_dnSca
!> @brief Public function to get nderiv from a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                     TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
!! @param get_nderiv_FROM_dnSca   integer  (result):   nderiv value, check against S%nederiv and the allocated d1,d2 or d3
  FUNCTION get_nderiv_FROM_dnSca(S)
    integer :: get_nderiv_FROM_dnSca

    TYPE(dnSca), intent(in)    :: S
    integer :: nderiv

    nderiv = S%nderiv

    IF (.NOT. allocated(S%d1)) THEN
      nderiv = 0
    ELSE IF (.NOT. allocated(S%d2)) THEN
      nderiv = 1
    ELSE IF (.NOT. allocated(S%d3)) THEN
      nderiv = 2
    ELSE
      nderiv = 3
    END IF

    IF (S%nderiv /= nderiv) THEN
      write(out_unitp,*) ' ERROR in get_nderiv_FROM_dnSca'
      write(out_unitp,*) '  Problem with nderiv in S'
      CALL Write_dnSca(S)
      STOP 'ERROR in get_nderiv_FROM_dnSca'
    END IF

    get_nderiv_FROM_dnSca = nderiv

    END FUNCTION get_nderiv_FROM_dnSca
!> @brief Public function to get ndim from a derived type dnSca.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                     TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
!! @param get_ndim_FROM_dnSca     integer  (result):   ndim value from the size of S%d1.
  FUNCTION get_ndim_FROM_dnSca(S)
    integer :: get_ndim_FROM_dnSca

    TYPE(dnSca), intent(in)    :: S

    IF (.NOT. allocated(S%d1)) THEN
      get_ndim_FROM_dnSca = 0
    ELSE
      get_ndim_FROM_dnSca = size(S%d1(:))
    END IF

  END FUNCTION get_ndim_FROM_dnSca
!> @brief Public function which ckecks a derived type dnSca is zero (all components).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Check_dnSca_IS_ZERO     logical  (result):   result of the comparison
!! @param S                     TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
!! @param epsi                  real (optional):     when present zero limit, otherwise 10^-10
  FUNCTION Check_dnSca_IS_ZERO(S,epsi)
    USE mod_NumParameters
    logical :: Check_dnSca_IS_ZERO

    TYPE(dnSca),          intent(in)           :: S
    real(kind=Rkind),   intent(in), optional :: epsi

    real(kind=Rkind) :: epsi_loc


    IF (present(epsi)) THEN
      epsi_loc = epsi
    ELSE
      epsi_loc = ONETENTH**10
    END IF


    Check_dnSca_IS_ZERO = (get_maxval_OF_dnSca(S) <= epsi_loc)

  END FUNCTION Check_dnSca_IS_ZERO
!> @brief Public function which gets the largest value of a derived type dnSca (all components).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param get_maxval_OF_dnSca     real  (result):      largest value
!! @param S                     TYPE(dnSca):           derived type which deals with the derivatives of a scalar functions.
  FUNCTION get_maxval_OF_dnSca(S)
    USE mod_NumParameters

    real(kind=Rkind) :: get_maxval_OF_dnSca

    TYPE(dnSca), intent(in)    :: S

    real(kind=Rkind) :: e0,e1,e2,e3

    e1 = ZERO
    e2 = ZERO
    e3 = ZERO
    e0 = abs(S%d0)
    IF (allocated(S%d1)) e1 = maxval(abs(S%d1))
    IF (allocated(S%d2)) e2 = maxval(abs(S%d2))
    IF (allocated(S%d3)) e3 = maxval(abs(S%d3))

    get_maxval_OF_dnSca = max(e0,e1,e2,e3)

  END FUNCTION get_maxval_OF_dnSca
!> @brief Public function which calculates numerical derivative of a function
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!

!! @param Snum                     TYPE(dnSca):       function value, f(x), and derivatives f' f" f'".
!! @param f                        real:            function (intrinsic or external)
!! @param x                        real:            abciss
!! @param nderiv                   integer:         order of the derivative
  FUNCTION get_Num_dnSca_FROM_f_x(x,f,nderiv) RESULT(Snum)
    TYPE (dnSca)                                 :: Snum
    real (kind=Rkind)                          :: f ! an intrinsic function: sin exp ....
    real (kind=Rkind),           intent(in)    :: x
    integer,           optional, intent(in)    :: nderiv

    !local variables:
    TYPE (dnSca)            :: Sloc
    real (kind=Rkind)     :: xloc,step=ONETENTH**4
    integer               :: nderiv_loc,i,j,k
    real (kind=Rkind)     :: f0,fp,fm,fpp,fmm,fppp,fmmm


    character (len=*), parameter :: name_sub='get_Num_dnSca_FROM_f_x'

    ! test nderiv
    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(3,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF

    CALL alloc_dnSca(Snum,ndim=1,nderiv=nderiv_loc)
    Snum = ZERO

    SELECT CASE (nderiv_loc)
    CASE(0)
      Snum%d0 = f(x)
    CASE(1)
      f0  = f(x)
      fp  = f(x+step)
      fm  = f(x-step)
      !fpp = f(x+step+step)
      !fmm = f(x-step-step)

      Snum%d0 = f0
      Snum%d1 = (fp-fm)/(step+step)
    CASE(2)
      f0  = f(x)
      fp  = f(x+step)
      fm  = f(x-step)
      fpp = f(x+step+step)
      fmm = f(x-step-step)

      Snum%d0 = f0
      Snum%d1 = (fp-fm)/(step+step)
      Snum%d2 = (fp+fm-TWO*f0)/step**2
    CASE(3)
      f0   = f(x)
      fp   = f(x+step)
      fm   = f(x-step)
      fpp  = f(x+step+step)
      fmm  = f(x-step-step)
      fppp = f(x+step+step+step)
      fmmm = f(x-step-step-step)

      Snum%d0 = f0
      Snum%d1 = ( THREE/FOUR*(fp-fm) - &
                 THREE/20_Rkind*(fpp-fmm) + &
                 ONE/60_Rkind*(fppp-fmmm) )/step

      Snum%d2 = (-30_Rkind*f0+16_Rkind*(fp+fm)-(fpp+fmm))/(TWELVE*step**2)

      Snum%d3 = (-13_Rkind*(fp-fm)+EIGHT*(fpp-fmm)-(fppp-fmmm))/(EIGHT*step**3)

    END SELECT

  END FUNCTION get_Num_dnSca_FROM_f_x

!=========================================================
! operators ==,/=,>=,>,<=,<
!=========================================================
  FUNCTION dnSca_EQ_dnSca(S1,S2) RESULT(lres)
    TYPE (dnSca), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_EQ_dnSca'

    lres = (S1%d0 == S2%d0)

    IF (allocated(S1%d1) .AND. allocated(S2%d1)) lres = lres .AND. all(S1%d1 == S2%d1)
    IF (allocated(S1%d2) .AND. allocated(S2%d2)) lres = lres .AND. all(S1%d2 == S2%d2)
    IF (allocated(S1%d3) .AND. allocated(S2%d3)) lres = lres .AND. all(S1%d3 == S2%d3)

  END FUNCTION dnSca_EQ_dnSca
  FUNCTION dnSca_EQ_R(S1,R) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_EQ_R'

    lres = (S1%d0 == R)

  END FUNCTION dnSca_EQ_R
  FUNCTION R_EQ_dnSca(R,S1) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_EQ_dnSca'

    lres = (R == S1%d0)

  END FUNCTION R_EQ_dnSca

  FUNCTION dnSca_NEQ_dnSca(S1,S2) RESULT(lres)
    TYPE (dnSca), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_NEQ_dnSca'

    lres = (S1%d0 /= S2%d0)

    IF (allocated(S1%d1) .AND. allocated(S2%d1)) lres = lres .OR. all(S1%d1 /= S2%d1)
    IF (allocated(S1%d2) .AND. allocated(S2%d2)) lres = lres .OR. all(S1%d2 /= S2%d2)
    IF (allocated(S1%d3) .AND. allocated(S2%d3)) lres = lres .OR. all(S1%d3 /= S2%d3)

  END FUNCTION dnSca_NEQ_dnSca
  FUNCTION dnSca_NEQ_R(S1,R) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_NEQ_R'

    lres = (S1%d0 /= R)

  END FUNCTION dnSca_NEQ_R
  FUNCTION R_NEQ_dnSca(R,S1) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_NEQ_dnSca'

    lres = (R /= S1%d0)

  END FUNCTION R_NEQ_dnSca

  FUNCTION dnSca_LE_dnSca(S1,S2) RESULT(lres)
    TYPE (dnSca), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_LE_dnSca'

    lres = (S1%d0 <= S2%d0)

  END FUNCTION dnSca_LE_dnSca
  FUNCTION dnSca_LE_R(S1,R) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_LE_R'

    lres = (S1%d0 <= R)

  END FUNCTION dnSca_LE_R
  FUNCTION R_LE_dnSca(R,S1) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_LE_dnSca'

    lres = (R <= S1%d0)

  END FUNCTION R_LE_dnSca

  FUNCTION dnSca_LT_dnSca(S1,S2) RESULT(lres)
    TYPE (dnSca), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_LT_dnSca'

    lres = (S1%d0 < S2%d0)

  END FUNCTION dnSca_LT_dnSca
  FUNCTION dnSca_LT_R(S1,R) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_LT_R'

    lres = (S1%d0 < R)

  END FUNCTION dnSca_LT_R
  FUNCTION R_LT_dnSca(R,S1) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_LT_dnSca'

    lres = (R < S1%d0)

  END FUNCTION R_LT_dnSca

  FUNCTION dnSca_GE_dnSca(S1,S2) RESULT(lres)
    TYPE (dnSca), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_GE_dnSca'

    lres = (S1%d0 >= S2%d0)

  END FUNCTION dnSca_GE_dnSca
  FUNCTION dnSca_GE_R(S1,R) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_GE_R'

    lres = (S1%d0 >= R)

  END FUNCTION dnSca_GE_R
  FUNCTION R_GE_dnSca(R,S1) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_GE_dnSca'

    lres = (R >= S1%d0)

  END FUNCTION R_GE_dnSca

  FUNCTION dnSca_GT_dnSca(S1,S2) RESULT(lres)
    TYPE (dnSca), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_GT_dnSca'

    lres = (S1%d0 > S2%d0)

  END FUNCTION dnSca_GT_dnSca
  FUNCTION dnSca_GT_R(S1,R) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_GT_R'

    lres = (S1%d0 > R)

  END FUNCTION dnSca_GT_R
  FUNCTION R_GT_dnSca(R,S1) RESULT(lres)
    TYPE (dnSca),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_GT_dnSca'

    lres = (R > S1%d0)

  END FUNCTION R_GT_dnSca
!=========================================================
! operators =,+,-,*,/,**
!=========================================================
  SUBROUTINE sub_dnSca2_TO_dnSca1(S1,S2)
    TYPE (dnSca), intent(inout) :: S1
    TYPE (dnSca), intent(in)    :: S2

    integer :: nderiv_loc,ndim_loc
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='sub_dnSca2_TO_dnSca1'

    nderiv_loc = get_nderiv_FROM_dnSca(S2)
    ndim_loc   = get_ndim_FROM_dnSca(S2)

    !write(out_unitp,*) 'ndim,nsurf,nderiv',ndim_loc,nderiv_loc

    CALL alloc_dnSca(S1,ndim_loc,nderiv_loc)

    IF (nderiv_loc < 0 .OR. (nderiv_loc > 0 .AND. ndim_loc < 1)) RETURN


    IF (nderiv_loc == 0) THEN
       S1%d0 = S2%d0
    ELSE IF (nderiv_loc == 1) THEN
       S1%d0 = S2%d0
       S1%d1 = S2%d1
    ELSE IF (nderiv_loc == 2) THEN
       S1%d0 = S2%d0
       S1%d1 = S2%d1
       S1%d2 = S2%d2
    ELSE IF (nderiv_loc == 3) THEN
       S1%d0 = S2%d0
       S1%d1 = S2%d1
       S1%d2 = S2%d2
       S1%d3 = S2%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END SUBROUTINE sub_dnSca2_TO_dnSca1

  SUBROUTINE set_dnSca_TO_R(S,R)

    TYPE (dnSca),      intent(inout) :: S
    real (kind=Rkind), intent(in)    :: R

    integer :: nderiv_loc,ndim_loc
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='set_dnSca_TO_R'

    nderiv_loc = get_nderiv_FROM_dnSca(S)

    !write(out_unitp,*) 'nderiv',nderiv_loc

    IF (nderiv_loc == 0) THEN
       S%d0 = R
    ELSE IF (nderiv_loc == 1) THEN
       S%d0 = R
       S%d1 = ZERO
    ELSE IF (nderiv_loc == 2) THEN
       S%d0 = R
       S%d1 = ZERO
       S%d2 = ZERO
    ELSE IF (nderiv_loc == 3) THEN
       S%d0 = R
       S%d1 = ZERO
       S%d2 = ZERO
       S%d3 = ZERO
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv_loc
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END SUBROUTINE set_dnSca_TO_R
!=========================================================


 FUNCTION dnSca2_PLUS_dnSca1(S1,S2) RESULT(Sres)
    TYPE (dnSca)                :: Sres
    TYPE (dnSca), intent(in)    :: S1,S2

    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca2_PLUS_dnSca1'

    nderiv = min(get_nderiv_FROM_dnSca(S1),get_nderiv_FROM_dnSca(S2))
    ndim   = min(get_ndim_FROM_dnSca(S1),  get_ndim_FROM_dnSca(S2))


    CALL dealloc_dnSca(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    CALL alloc_dnSca(Sres,ndim,nderiv)

    IF (nderiv == 0) THEN
       Sres%d0 = S1%d0 + S2%d0
    ELSE IF (nderiv == 1) THEN
       Sres%d0 = S1%d0 + S2%d0
       Sres%d1 = S1%d1 + S2%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = S1%d0 + S2%d0
       Sres%d1 = S1%d1 + S2%d1
       Sres%d2 = S1%d2 + S2%d2
   ELSE IF (nderiv == 3) THEN
       Sres%d0 = S1%d0 + S2%d0
       Sres%d1 = S1%d1 + S2%d1
       Sres%d2 = S1%d2 + S2%d2
       Sres%d3 = S1%d3 + S2%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 2 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION dnSca2_PLUS_dnSca1
  FUNCTION dnSca_PLUS_R(S,R) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_PLUS_R'

    Sres = S
    Sres%d0 = Sres%d0 + R

  END FUNCTION dnSca_PLUS_R
  FUNCTION R_PLUS_dnSca(R,S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_PLUS_dnSca'

    Sres = S
    Sres%d0 = Sres%d0 + R

  END FUNCTION R_PLUS_dnSca
  FUNCTION PLUS_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='PLUS_dnSca'

    Sres = S

  END FUNCTION PLUS_dnSca
  FUNCTION dnSca2_MINUS_dnSca1(S1,S2) RESULT(Sres)
    TYPE (dnSca)                :: Sres
    TYPE (dnSca), intent(in)    :: S1,S2

    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca2_MINUS_dnSca1'

    nderiv = min(get_nderiv_FROM_dnSca(S1),get_nderiv_FROM_dnSca(S2))
    ndim   = min(get_ndim_FROM_dnSca(S1),  get_ndim_FROM_dnSca(S2))

    CALL dealloc_dnSca(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    CALL alloc_dnSca(Sres,ndim,nderiv)

    IF (nderiv == 0) THEN
       Sres%d0 = S1%d0 - S2%d0
    ELSE IF (nderiv == 1) THEN
       Sres%d0 = S1%d0 - S2%d0
       Sres%d1 = S1%d1 - S2%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = S1%d0 - S2%d0
       Sres%d1 = S1%d1 - S2%d1
       Sres%d2 = S1%d2 - S2%d2
    ELSE IF (nderiv == 3) THEN
       Sres%d0 = S1%d0 - S2%d0
       Sres%d1 = S1%d1 - S2%d1
       Sres%d2 = S1%d2 - S2%d2
       Sres%d3 = S1%d3 - S2%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION dnSca2_MINUS_dnSca1
  FUNCTION dnSca_MINUS_R(S,R) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_MINUS_R'

    Sres = S
    Sres%d0 = Sres%d0 - R

  END FUNCTION dnSca_MINUS_R
  FUNCTION R_MINUS_dnSca(R,S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_MINUS_dnSca'

    nderiv = get_nderiv_FROM_dnSca(S)
    ndim   = get_ndim_FROM_dnSca(S)


    CALL dealloc_dnSca(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    CALL alloc_dnSca(Sres,ndim,nderiv)

    IF (nderiv == 0) THEN
       Sres%d0 = R -S%d0
    ELSE IF (nderiv == 1) THEN
       Sres%d0 = R -S%d0
       Sres%d1 =   -S%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = R - S%d0
       Sres%d1 =   -S%d1
       Sres%d2 =   -S%d2
    ELSE IF (nderiv == 3) THEN
       Sres%d0 = R - S%d0
       Sres%d1 =   -S%d1
       Sres%d2 =   -S%d2
       Sres%d3 =   -S%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF

  END FUNCTION R_MINUS_dnSca
  FUNCTION MINUS_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='MINUS_dnSca'

    nderiv = get_nderiv_FROM_dnSca(S)
    ndim   = get_ndim_FROM_dnSca(S)

    CALL dealloc_dnSca(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    CALL alloc_dnSca(Sres,ndim,nderiv)

    IF (nderiv == 0) THEN
       Sres%d0 = -S%d0
    ELSE IF (nderiv == 1) THEN
       Sres%d0 = -S%d0
       Sres%d1 = -S%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = -S%d0
       Sres%d1 = -S%d1
       Sres%d2 = -S%d2
    ELSE IF (nderiv == 3) THEN
       Sres%d0 = -S%d0
       Sres%d1 = -S%d1
       Sres%d2 = -S%d2
       Sres%d3 = -S%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF

  END FUNCTION MINUS_dnSca

 FUNCTION dnSca2_TIME_dnSca1(S1,S2) RESULT(Sres)
    TYPE (dnSca)                :: Sres
    TYPE (dnSca), intent(in)    :: S1,S2

    integer :: nderiv,ndim,id,jd,kd
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca2_TIME_dnSca1'

    nderiv = min(get_nderiv_FROM_dnSca(S1),get_nderiv_FROM_dnSca(S2))
    ndim   = min(get_ndim_FROM_dnSca(S1),  get_ndim_FROM_dnSca(S2))


    CALL dealloc_dnSca(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    CALL alloc_dnSca(Sres,ndim,nderiv)

    IF (nderiv == 0) THEN
       Sres%d0 = S1%d0 * S2%d0
    ELSE IF (nderiv == 1) THEN
       Sres%d0 = S1%d0 * S2%d0
       Sres%d1 = S1%d1 * S2%d0 + S1%d0 * S2%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = S1%d0 * S2%d0
       Sres%d1 = S1%d1 * S2%d0 + S1%d0 * S2%d1
       Sres%d2 = S1%d2 * S2%d0 + S1%d0 * S2%d2
       DO id=1,ndim
       DO jd=1,ndim
         Sres%d2(jd,id) = Sres%d2(jd,id) + S1%d1(id) * S2%d1(jd) + S1%d1(jd) * S2%d1(id)
       END DO
       END DO
    ELSE IF (nderiv == 3) THEN
       Sres%d0 = S1%d0 * S2%d0
       Sres%d1 = S1%d1 * S2%d0 + S1%d0 * S2%d1

       Sres%d2 = S1%d2 * S2%d0 + S1%d0 * S2%d2
       DO id=1,ndim
       DO jd=1,ndim
         Sres%d2(jd,id) = Sres%d2(jd,id) + S1%d1(id) * S2%d1(jd) + S1%d1(jd) * S2%d1(id)
       END DO
       END DO

       Sres%d3 = S1%d3 * S2%d0 + S1%d0 * S2%d3
       DO id=1,ndim
       DO jd=1,ndim
       DO kd=1,ndim
         Sres%d3(kd,jd,id) = Sres%d3(kd,jd,id) + S1%d1(id) * S2%d2(kd,jd) + &
                                                 S1%d1(jd) * S2%d2(kd,id) + &
                                                 S1%d1(kd) * S2%d2(jd,id) + &
                                                 S1%d2(kd,jd) * S2%d1(id) + &
                                                 S1%d2(kd,id) * S2%d1(jd) + &
                                                 S1%d2(jd,id) * S2%d1(kd)

       END DO
       END DO
       END DO
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF
  END FUNCTION dnSca2_TIME_dnSca1
  FUNCTION dnSca_TIME_R(S,R) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_TIME_R'

    nderiv = get_nderiv_FROM_dnSca(S)
    ndim   = get_ndim_FROM_dnSca(S)


    CALL dealloc_dnSca(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    CALL alloc_dnSca(Sres,ndim,nderiv)

    IF (nderiv == 0) THEN
       Sres%d0 = R * S%d0
    ELSE IF (nderiv == 1) THEN
       Sres%d0 = R * S%d0
       Sres%d1 = R * S%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = R * S%d0
       Sres%d1 = R * S%d1
       Sres%d2 = R * S%d2
    ELSE IF (nderiv == 3) THEN
       Sres%d0 = R * S%d0
       Sres%d1 = R * S%d1
       Sres%d2 = R * S%d2
       Sres%d3 = R * S%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF

  END FUNCTION dnSca_TIME_R

  FUNCTION R_TIME_dnSca(R,S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_TIME_dnSca'

    nderiv = get_nderiv_FROM_dnSca(S)
    ndim   = get_ndim_FROM_dnSca(S)


    CALL dealloc_dnSca(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    CALL alloc_dnSca(Sres,ndim,nderiv)

    IF (nderiv == 0) THEN
       Sres%d0 = R * S%d0
    ELSE IF (nderiv == 1) THEN
       Sres%d0 = R * S%d0
       Sres%d1 = R * S%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = R * S%d0
       Sres%d1 = R * S%d1
       Sres%d2 = R * S%d2
    ELSE IF (nderiv == 3) THEN
       Sres%d0 = R * S%d0
       Sres%d1 = R * S%d1
       Sres%d2 = R * S%d2
       Sres%d3 = R * S%d3
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF

  END FUNCTION R_TIME_dnSca

 FUNCTION dnSca2_OVER_dnSca1(S1,S2) RESULT(Sres)
    TYPE (dnSca)                :: Sres
    TYPE (dnSca), intent(in)    :: S1,S2

    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca2_OVER_dnSca1'


    Sres = S1 * S2**(-ONE)


  END FUNCTION dnSca2_OVER_dnSca1
  FUNCTION dnSca_OVER_R(S,R) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_OVER_R'


    Sres = S * R**(-ONE)

  END FUNCTION dnSca_OVER_R

  FUNCTION R_OVER_dnSca(R,S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='R_OVER_dnSca'


    Sres = R * S**(-ONE)


  END FUNCTION R_OVER_dnSca

!=========================================================
! mathematical functions: cos, sin exp, log, cosh ....
! All functions in the fortran norm except atan2 betause it has two arguments
!=========================================================

  FUNCTION get_F_dnSca(S,d0f,d1f,d2f,d3f) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: d0f,d1f,d2f,d3f

    integer :: nderiv,ndim,id,jd,kd
    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='get_F_dnSca'

    nderiv = get_nderiv_FROM_dnSca(S)
    ndim   = get_ndim_FROM_dnSca(S)

    CALL dealloc_dnSca(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    CALL alloc_dnSca(Sres,ndim,nderiv)

    IF (nderiv == 0) THEN
       Sres%d0 = d0f
    ELSE IF (nderiv == 1) THEN
       Sres%d0 =  d0f
       Sres%d1 =  d1f * S%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = d0f
       Sres%d1 = d1f * S%d1
       DO id=1,ndim
       DO jd=1,ndim
         Sres%d2(jd,id) = d1f * S%d2(jd,id) + d2f * S%d1(id)*S%d1(jd)
       END DO
       END DO
    ELSE IF (nderiv == 3) THEN
       Sres%d0 = d0f
       Sres%d1 = d1f * S%d1
       DO id=1,ndim
       DO jd=1,ndim
         Sres%d2(jd,id) = d1f * S%d2(jd,id) + d2f * S%d1(id)*S%d1(jd)
       END DO
       END DO

       DO id=1,ndim
       DO jd=1,ndim
       DO kd=1,ndim
         Sres%d3(kd,jd,id) = d1f * S%d3(kd,jd,id) + &
                             d2f * (S%d1(id)*S%d2(kd,jd) + S%d1(jd)*S%d2(kd,id) + S%d1(kd)*S%d2(jd,id)) + &
                             d3f * S%d1(id)*S%d1(jd)*S%d1(kd)
       END DO
       END DO
       END DO

    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nderiv > 3 is NOT possible',nderiv
      write(out_unitp,*) 'It should never append! Check the source'
      STOP
    END IF

  END FUNCTION get_F_dnSca

  FUNCTION dnSca_EXP_R(S,R) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='dnSca_EXP_R'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f = S%d0**R
    IF (nderiv >= 1) d1f = S%d0**(R-ONE)   * R
    IF (nderiv >= 2) d2f = S%d0**(R-TWO)   * R*(R-ONE)
    IF (nderiv >= 3) d3f = S%d0**(R-THREE) * R*(R-ONE)*(R-TWO)

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION dnSca_EXP_R
  FUNCTION dnSca_EXP_I(S,I) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S
    integer,           intent(in)    :: I


    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='dnSca_EXP_I'

    Sres = S**real(I,kind=Rkind)

  END FUNCTION dnSca_EXP_I
  FUNCTION get_SQRT_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S


    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='get_SQRT_dnSca'

    Sres = S**HALF

  END FUNCTION get_SQRT_dnSca

  FUNCTION get_ABS_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S


    integer :: err_dnSca_loc
    character (len=*), parameter :: name_sub='get_ABS_dnSca'

    Sres = (S*S)**HALF

  END FUNCTION get_ABS_dnSca

  FUNCTION get_EXP_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f

    character (len=*), parameter :: name_sub='get_EXP_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  exp(S%d0)
    IF (nderiv >= 1) d1f =  exp(S%d0)
    IF (nderiv >= 2) d2f =  exp(S%d0)
    IF (nderiv >= 3) d3f =  exp(S%d0)

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_EXP_dnSca
  FUNCTION get_LOG_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f

    character (len=*), parameter :: name_sub='get_LOG_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  log(S%d0)
    IF (nderiv >= 1) d1f =  ONE/S%d0
    IF (nderiv >= 2) d2f = -ONE/S%d0**2
    IF (nderiv >= 3) d3f =  TWO/S%d0**3

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_LOG_dnSca
  FUNCTION get_LOG10_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f

    character (len=*), parameter :: name_sub='get_LOG10_dnSca'

    Sres = log(S)/log(TEN)

  END FUNCTION get_LOG10_dnSca
  FUNCTION get_COS_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S


    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_COS_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  cos(S%d0)
    IF (nderiv >= 1) d1f = -sin(S%d0)
    IF (nderiv >= 2) d2f = -cos(S%d0)
    IF (nderiv >= 3) d3f =  sin(S%d0)

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_COS_dnSca
  FUNCTION get_ACOS_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S


    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ACOS_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  acos(S%d0)
    IF (nderiv >= 1) d1f = -ONE/sqrt(ONE-S%d0**2)
    IF (nderiv >= 2) d2f = S%d0*d1f**3
    IF (nderiv >= 3) d3f = (ONE+TWO*S%d0**2)*d1f**5

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ACOS_dnSca
  FUNCTION get_SIN_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_SIN_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  sin(S%d0)
    IF (nderiv >= 1) d1f =  cos(S%d0)
    IF (nderiv >= 2) d2f = -sin(S%d0)
    IF (nderiv >= 3) d3f = -cos(S%d0)

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_SIN_dnSca
  FUNCTION get_ASIN_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ASIN_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  asin(S%d0)
    IF (nderiv >= 1) d1f = ONE/sqrt(ONE-S%d0**2)
    IF (nderiv >= 2) d2f = S%d0*d1f**3
    IF (nderiv >= 3) d3f = (ONE+TWO*S%d0**2)*d1f**5

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ASIN_dnSca
  FUNCTION get_TAN_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_TAN_dnSca'

    Sres = Sin(S)/cos(S)

  END FUNCTION get_TAN_dnSca

  FUNCTION get_ATAN_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ATAN_dnSca'

    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  atan(S%d0)
    IF (nderiv >= 1) d1f = ONE/(ONE+S%d0**2)
    IF (nderiv >= 2) d2f = -TWO*S%d0 * d1f**2
    IF (nderiv >= 3) d3f = (-TWO+SIX*S%d0**2) * d1f**3

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ATAN_dnSca

  FUNCTION get_COSH_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S


    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_COSH_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  cosh(S%d0)
    IF (nderiv >= 1) d1f =  sinh(S%d0)
    IF (nderiv >= 2) d2f =  cosh(S%d0)
    IF (nderiv >= 3) d3f =  sinh(S%d0)

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_COSH_dnSca
  FUNCTION get_ACOSH_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S


    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ACOSH_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  acosh(S%d0)
    IF (nderiv >= 1) d1f = ONE/sqrt(-ONE+S%d0**2)
    IF (nderiv >= 2) d2f = -S%d0*d1f**3
    IF (nderiv >= 3) d3f = (ONE+TWO*S%d0**2)*d1f**5

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ACOSH_dnSca
  FUNCTION get_SINH_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_SINH_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  sinh(S%d0)
    IF (nderiv >= 1) d1f =  cosh(S%d0)
    IF (nderiv >= 2) d2f =  sinh(S%d0)
    IF (nderiv >= 3) d3f =  cosh(S%d0)

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_SINH_dnSca
  FUNCTION get_ASINH_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ASINH_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  asinh(S%d0)
    IF (nderiv >= 1) d1f = ONE/sqrt(ONE+S%d0**2)
    IF (nderiv >= 2) d2f = -S%d0*d1f**3
    IF (nderiv >= 3) d3f = (-ONE+TWO*S%d0**2)*d1f**5

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ASINH_dnSca

  FUNCTION get_TANH_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_TANH_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  tanh(S%d0)
    IF (nderiv >= 1) d1f =  ONE/cosh(S%d0)**2
    IF (nderiv >= 2) d2f = -TWO*tanh(S%d0) * d1f
    IF (nderiv >= 3) d3f = (-FOUR+TWO*cosh(TWO*S%d0)) * d1f**2

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_TANH_dnSca

  FUNCTION get_ATANH_dnSca(S) RESULT(Sres)
    TYPE (dnSca)                       :: Sres
    TYPE (dnSca),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnSca_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ATANH_dnSca'


    nderiv = get_nderiv_FROM_dnSca(S)

    IF (nderiv >= 0) d0f =  atanh(S%d0)
    IF (nderiv >= 1) d1f =  ONE/(ONE-S%d0**2)
    IF (nderiv >= 2) d2f =  TWO*S%d0 * d1f**2
    IF (nderiv >= 3) d3f =  (TWO+SIX*S%d0**2) * d1f**3

    Sres = get_F_dnSca(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ATANH_dnSca

END MODULE mod_dnSca
