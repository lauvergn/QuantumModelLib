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
!! This module deals with operations or functions of a scalar function and its derivatives, dnS.
!!
!! There is a mapping between the saclar function S, its derivatives and the dnS derived type components:
!!
!! @li S                 => S%d0
!! @li dS/dQ_i           => S%d1(i)
!! @li d^2S/dQ_idQ_j     => S%d2(i,j)
!! @li d^3S/dQ_idQ_jdQ_k => S%d3(i,j,k)
!!
!! with S defined as:
!!  TYPE(dnS) :: S
!!
!!
!! All standard fortran operators (= + - * / **) are overloaded:
!!
!! For instance the sum (+) of two dnS variables, S1 and S2 corresponds to:
!! @li (S1+S2)                 => S1%d0    + S2%d0
!! @li d(S1+S2)/dQ_i           => S1%d1(i) + S2%d1(i)
!! @li ....
!!
!! The product (*) of two dnS variables, S1 and S2 correspond to:
!! @li (S1*S2)                 => S1%d0 * S1%d0
!! @li d(S1*S2)/dQ_i           => S1%d0 * S2%d1(i) + S1%d1(i) * S2%d0    (derivative of a product)
!! @li ....
!!
!! All standard fortran functions (exp, sqrt, log ... sin ... sinh) are overloaded
!!
!! For instance the function, f, of dnS variables, S, corresponds to:
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
MODULE mod_dnS
  USE mod_NumParameters
  IMPLICIT NONE
  PRIVATE

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
  TYPE, PUBLIC :: dnS
     PRIVATE
     integer                        :: nderiv = -1

     real (kind=Rkind)              :: d0
     real (kind=Rkind), allocatable :: d1(:)
     real (kind=Rkind), allocatable :: d2(:,:)
     real (kind=Rkind), allocatable :: d3(:,:,:)
  CONTAINS
    PROCEDURE, PRIVATE :: set_dnS_TO_R
    PROCEDURE, PRIVATE :: sub_dnS2_TO_dnS1
    GENERIC,   PUBLIC  :: assignment(=) => sub_dnS2_TO_dnS1,set_dnS_TO_R
  END TYPE dnS

  ! overloded operators, functions
  ! Rk: the assignment(=) cannot be here, since it is bounded to the dnS type
  PUBLIC :: operator (+),operator (-),operator (*),operator (/),operator (**)
  PUBLIC :: operator (==),operator (/=),operator (<=),operator (>=),operator (<),operator (>)
  PUBLIC :: dot_product
  PUBLIC :: sqrt,exp,abs,log,log10
  PUBLIC :: cos,sin,tan,acos,asin,atan,cosh,sinh,tanh,acosh,asinh,atanh

  PUBLIC :: alloc_dnS,dealloc_dnS,init_dnS,set_dnS
  PUBLIC :: Write_dnS,WriteAll_dnS,Write_dnS_FOR_test
  PUBLIC :: get_nderiv_FROM_dnS,get_ndim_FROM_dnS
  PUBLIC :: sub_get_d1_FROM_dnS,sub_get_d2_FROM_dnS,sub_get_d3_FROM_dnS
  PUBLIC :: get_d0_FROM_dnS,get_d1_FROM_dnS,get_d2_FROM_dnS,get_d3_FROM_dnS

  PUBLIC :: get_Num_dnS_FROM_f_x,Check_dnS_IS_ZERO,d0S_TIME_R


  INTERFACE operator (+)
     MODULE PROCEDURE dnS2_PLUS_dnS1,dnS_PLUS_R,R_PLUS_dnS,PLUS_dnS
  END INTERFACE
  INTERFACE operator (-)
     MODULE PROCEDURE dnS2_MINUS_dnS1,dnS_MINUS_R,R_MINUS_dnS,MINUS_dnS
  END INTERFACE
  INTERFACE operator (*)
     MODULE PROCEDURE dnS2_TIME_dnS1,dnS_TIME_R,R_TIME_dnS
  END INTERFACE
  INTERFACE operator (/)
     MODULE PROCEDURE dnS2_OVER_dnS1,dnS_OVER_R,R_OVER_dnS
  END INTERFACE
  INTERFACE operator (**)
     MODULE PROCEDURE dnS_EXP_R,dnS_EXP_I
  END INTERFACE

  INTERFACE operator (==)
     MODULE PROCEDURE dnS_EQ_dnS,dnS_EQ_R,R_EQ_dnS
  END INTERFACE
  INTERFACE operator (/=)
     MODULE PROCEDURE dnS_NEQ_dnS,dnS_NEQ_R,R_NEQ_dnS
  END INTERFACE
  INTERFACE operator (<=)
     MODULE PROCEDURE dnS_LE_dnS,dnS_LE_R,R_LE_dnS
  END INTERFACE
  INTERFACE operator (<)
     MODULE PROCEDURE dnS_LT_dnS,dnS_LT_R,R_LT_dnS
  END INTERFACE
  INTERFACE operator (>=)
     MODULE PROCEDURE dnS_GE_dnS,dnS_GE_R,R_GE_dnS
  END INTERFACE
  INTERFACE operator (>)
     MODULE PROCEDURE dnS_GT_dnS,dnS_GT_R,R_GT_dnS
  END INTERFACE

  INTERFACE sqrt
     MODULE PROCEDURE get_SQRT_dnS
  END INTERFACE
  INTERFACE abs
     MODULE PROCEDURE get_ABS_dnS
  END INTERFACE
  INTERFACE exp
     MODULE PROCEDURE get_EXP_dnS
  END INTERFACE
  INTERFACE log
     MODULE PROCEDURE get_LOG_dnS
  END INTERFACE
  INTERFACE log10
     MODULE PROCEDURE get_LOG10_dnS
  END INTERFACE
  INTERFACE cos
     MODULE PROCEDURE get_COS_dnS
  END INTERFACE
  INTERFACE acos
     MODULE PROCEDURE get_ACOS_dnS
  END INTERFACE
  INTERFACE sin
     MODULE PROCEDURE get_SIN_dnS
  END INTERFACE
  INTERFACE asin
     MODULE PROCEDURE get_ASIN_dnS
  END INTERFACE
  INTERFACE tan
     MODULE PROCEDURE get_TAN_dnS
  END INTERFACE
  INTERFACE atan
     MODULE PROCEDURE get_ATAN_dnS
  END INTERFACE
  INTERFACE cosh
     MODULE PROCEDURE get_COSH_dnS
  END INTERFACE
  INTERFACE acosh
     MODULE PROCEDURE get_ACOSH_dnS
  END INTERFACE
  INTERFACE sinh
     MODULE PROCEDURE get_SINH_dnS
  END INTERFACE
  INTERFACE asinh
     MODULE PROCEDURE get_ASINH_dnS
  END INTERFACE
  INTERFACE tanh
     MODULE PROCEDURE get_TANH_dnS
  END INTERFACE
  INTERFACE atanh
     MODULE PROCEDURE get_ATANH_dnS
  END INTERFACE


  INTERFACE dot_product
     MODULE PROCEDURE dot_product_VecOFdnS
  END INTERFACE

CONTAINS
!> @brief Public subroutine which allocates a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
!! @param ndim               integer (optional):    number of variables (coordiantes) for the derivatives.
!! @param nderiv             integer (optional):    it enables to chose the derivative order (from 0 to 3).
!! @param err_dnS            integer (optional):  to handle the errors errors (0: no error).
  SUBROUTINE alloc_dnS(S,ndim,nderiv,err_dnS)
    USE mod_NumParameters
    TYPE(dnS), intent(inout)          :: S         !< derived type, which contains, matrix potential, its derivatives
    integer,     intent(in),  optional  :: ndim      !< number of variables (coordiantes)
    integer,     intent(in),  optional  :: nderiv    !< order of the derivatives [0,1,3]
    integer,     intent(out), optional  :: err_dnS !< to handle the errors

    ! local variables
    integer :: ndim_loc,err_dnS_loc,nderiv_loc

    err_dnS_loc = 0 ! no error
    IF (present(err_dnS)) err_dnS = 0

    CALL dealloc_dnS(S)

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

    !write(out_unitp,*) 'S%nderiv in alloc_dnS',S%nderiv


    IF (nderiv_loc >= 1) THEN
      allocate(S%d1(ndim_loc),stat=err_dnS_loc)
      IF (err_dnS_loc /= 0 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in alloc_dnS'
        write(out_unitp,*) '  Problem with allocate of S%d1'
        write(out_unitp,*) '  ndim > 0?',ndim_loc
        IF (present(err_dnS)) THEN
          err_dnS = err_dnS_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in alloc_dnS'
        END IF
      END IF
    END IF

    IF (nderiv_loc >= 2) THEN
      allocate(S%d2(ndim_loc,ndim_loc),stat=err_dnS_loc)
      IF (err_dnS_loc /= 0 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in alloc_dnS'
        write(out_unitp,*) '  Problem with allocate of S%d2'
        write(out_unitp,*) '  ndim > 0',ndim_loc
        IF (present(err_dnS)) THEN
          err_dnS = err_dnS_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in alloc_dnS'
        END IF
      END IF
    END IF

    IF (nderiv_loc >= 3) THEN
      allocate(S%d3(ndim_loc,ndim_loc,ndim_loc),stat=err_dnS_loc)
      IF (err_dnS_loc /= 0 .OR. ndim_loc < 1) THEN
        write(out_unitp,*) ' ERROR in alloc_dnS'
        write(out_unitp,*) '  Problem with allocate of S%d3'
        write(out_unitp,*) '  ndim > 0',ndim_loc
        IF (present(err_dnS)) THEN
          err_dnS = err_dnS_loc
          RETURN
        ELSE
          STOP 'Problem with allocate in alloc_dnS'
        END IF
      END IF
    END IF

    !write(out_unitp,*) 'err_dnS_loc',err_dnS_loc
    !IF (present(err_dnS)) write(out_unitp,*) 'err_dnS',err_dnS

  END SUBROUTINE alloc_dnS

!> @brief Public subroutine which deallocates a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnS):             derived type which deals with the derivatives of a scalar functions.
  ELEMENTAL SUBROUTINE dealloc_dnS(S)
    USE mod_NumParameters

    TYPE(dnS), intent(inout)       :: S !< derived type, which contains, matrix potential, its derivatives

    IF (allocated(S%d1)) deallocate(S%d1)
    IF (allocated(S%d2)) deallocate(S%d2)
    IF (allocated(S%d3)) deallocate(S%d3)

    S%nderiv = -1

  END SUBROUTINE dealloc_dnS
!> @brief Public subroutine which checks if the derived type dnS is (correctly) allocated.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnS):  derived type which deals with the derivatives of a scalar functions.
!! @param nderiv             integer:    the derivative order.
  ELEMENTAL FUNCTION Check_NotAlloc_dnS(S,nderiv)
    USE mod_NumParameters

    logical :: Check_NotAlloc_dnS
    TYPE(dnS), intent(in)    :: S
    integer,     intent(in)    :: nderiv


    logical :: NotAlloc

    NotAlloc = (nderiv < 0)
    NotAlloc = NotAlloc .OR. (nderiv >= 1 .AND. .NOT. allocated(S%d1))
    NotAlloc = NotAlloc .OR. (nderiv >= 2 .AND. .NOT. allocated(S%d2))
    NotAlloc = NotAlloc .OR. (nderiv >= 3 .AND. .NOT. allocated(S%d3))

    Check_NotAlloc_dnS = NotAlloc

  END FUNCTION Check_NotAlloc_dnS
!> @brief Public function which initializes a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Sres               TYPE(dnS) (result):    derived type which deals with the derivatives of a scalar functions.
!! @param R                  real:                  Value of Sres%d0.
!! @param ndim               integer (optional):    number of variables (coordiantes) for the derivatives.
!! @param nderiv             integer (optional):    it enables to chose the derivative order (from 0 to 3).
!! @param iQ                 integer (optional):    when ndim > 1, dSres/dQ_iQ = Sres%d1(iQ)= 1, the other derivatives are zero
  FUNCTION init_dnS(R,ndim,nderiv,iQ) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                     :: Sres
    real (kind=Rkind), intent(in)    :: R
    integer, optional, intent(in)    :: nderiv,ndim,iQ

    integer :: err_dnS_loc,nderiv_loc,ndim_loc,iQ_loc
    character (len=*), parameter :: name_sub='init_dnS'

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
      STOP 'Problem with the iQ value in init_dnS'
    END IF

    !write(out_unitp,*) 'iQ_loc',iQ_loc ; flush(out_unitp)
    !write(out_unitp,*) 'ndim_loc,nderiv_loc',ndim_loc,nderiv_loc ; flush(out_unitp)

    CALL alloc_dnS(Sres,ndim_loc,nderiv_loc,err_dnS_loc)
    IF (err_dnS_loc /= 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Problem in alloc_dnS CALL'
      STOP 'Problem Problem in alloc_dnS CALL in init_dnS'
    END IF

    Sres = ZERO

    Sres%d0 = R
    IF (nderiv_loc > 0) Sres%d1(iQ_loc) = ONE

  END FUNCTION init_dnS
!> @brief Public subroutine which initializes a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
!! @param d0,d1,d2,d3        real  (optional):    real value (d0) or table to initialize S.
  SUBROUTINE set_dnS(S,d0,d1,d2,d3)
    USE mod_NumParameters

    real (kind=Rkind), optional,   intent(in)     :: d0
    real (kind=Rkind), optional,   intent(in)     :: d1(:)
    real (kind=Rkind), optional,   intent(in)     :: d2(:,:)
    real (kind=Rkind), optional,   intent(in)     :: d3(:,:,:)

    TYPE (dnS),                    intent(inout)  :: S


    character (len=*), parameter :: name_sub='set_dnS'

    IF (present(d0)) THEN
      S%d0 = d0
      S%nderiv = 0
    END IF

    IF (present(d1)) THEN
      S%d1     = d1
      S%nderiv = 1
      IF (.NOT. present(d0)) THEN
        write(out_unitp,*) ' ERROR in set_dnS'
        write(out_unitp,*) ' d1 is present but not d0'
        write(out_unitp,*) ' CHECK the fortran!!'
        STOP 'ERROR in set_dnS'
      END IF
    END IF

    IF (present(d2)) THEN
      S%d2     = d2
      S%nderiv = 2
      IF (.NOT. present(d1)) THEN
        write(out_unitp,*) ' ERROR in set_dnS'
        write(out_unitp,*) ' d2 is present but not d1'
        write(out_unitp,*) ' CHECK the fortran!!'
        STOP 'ERROR in set_dnS'
      END IF
    END IF

    IF (present(d3)) THEN
      S%d3     = d3
      S%nderiv = 3
      IF (.NOT. present(d2)) THEN
        write(out_unitp,*) ' ERROR in set_dnS'
        write(out_unitp,*) ' d3 is present but not d2'
        write(out_unitp,*) ' CHECK the fortran!!'
        STOP 'ERROR in set_dnS'
      END IF
    END IF

  END SUBROUTINE set_dnS
!> @brief Public function to get d0 from a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
!! @param d0                 real  (result):      d0=S%d0
  ELEMENTAL FUNCTION get_d0_FROM_dnS(S) RESULT(d0)
    USE mod_NumParameters

    real (kind=Rkind)                :: d0
    TYPE (dnS),        intent(in)  :: S

    character (len=*), parameter :: name_sub='get_d0_FROM_dnS'

    d0 = S%d0

  END FUNCTION get_d0_FROM_dnS
!> @brief Public function to get d1 from a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
!! @param d1                 real  (result):        d1=S%d1, d1 is allocatable
  FUNCTION get_d1_FROM_dnS(S) RESULT(d1)
    USE mod_NumParameters

    real (kind=Rkind), allocatable   :: d1(:)
    TYPE (dnS),        intent(in)    :: S

    character (len=*), parameter :: name_sub='get_d1_FROM_dnS'

    !IF (allocated(d1)) deallocate(d1)
    IF (allocated(S%d1)) d1 = S%d1

  END FUNCTION get_d1_FROM_dnS
  !> @brief Public function to get d1 from a derived type dnS. Similar to get_d1_FROM_dnS, but d1 is NOT allocatable.
  !!
  !> @author David Lauvergnat
  !! @date 25/05/2018
  !!
  !! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
  !! @param d1                 real       :           d1=S%d1, d1 is NOT allocatable
  SUBROUTINE sub_get_d1_FROM_dnS(d1,S)
    USE mod_NumParameters

      real (kind=Rkind),   intent(inout) :: d1(:)
      TYPE (dnS),        intent(in)    :: S

      character (len=*), parameter :: name_sub='sub_get_d1_FROM_dnS'

      IF (allocated(S%d1)) d1(:) = S%d1

    END SUBROUTINE sub_get_d1_FROM_dnS
  !> @brief Public function to get d2 from a derived type dnS.
  !!
  !> @author David Lauvergnat
  !! @date 03/08/2017
  !!
  !! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
  !! @param d2                 real  (result):        d2=S%d2, d2 is allocatable
  FUNCTION get_d2_FROM_dnS(S) RESULT(d2)
    USE mod_NumParameters

      real (kind=Rkind), allocatable   :: d2(:,:)
      TYPE (dnS),        intent(in)    :: S

      character (len=*), parameter :: name_sub='get_d2_FROM_dnS'

      !IF (allocated(d2)) deallocate(d2)
      IF (allocated(S%d2)) d2 = S%d2

    END FUNCTION get_d2_FROM_dnS
  !> @brief Public function to get d2 from a derived type dnS. Similar to get_d2_FROM_dnS, but d2 is NOT allocatable.
  !!
  !> @author David Lauvergnat
  !! @date 25/05/2018
  !!
  !! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
  !! @param d2                 real:                  d2=S%d2, d2 is NOT allocatable
  SUBROUTINE sub_get_d2_FROM_dnS(d2,S)
    USE mod_NumParameters

      real (kind=Rkind),   intent(inout) :: d2(:,:)
      TYPE (dnS),        intent(in)    :: S

      character (len=*), parameter :: name_sub='sub_get_d2_FROM_dnS'

      IF (allocated(S%d2)) d2(:,:) = S%d2

    END SUBROUTINE sub_get_d2_FROM_dnS
  !> @brief Public function to get d3 from a derived type dnS.
  !!
  !> @author David Lauvergnat
  !! @date 03/08/2017
  !!
  !! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
  !! @param d3                 real  (result):        d3=S%d3, d3 is allocatable
  FUNCTION get_d3_FROM_dnS(S) RESULT(d3)
    USE mod_NumParameters

      real (kind=Rkind), allocatable   :: d3(:,:,:)
      TYPE (dnS),        intent(in)    :: S

      character (len=*), parameter :: name_sub='get_d3_FROM_dnS'

      !IF (allocated(d3)) deallocate(d3)
      IF (allocated(S%d3)) d3 = S%d3

    END FUNCTION get_d3_FROM_dnS
  !> @brief Public function to get d3 from a derived type dnS. Similar to get_d3_FROM_dnS, but d3 is NOT allocatable.
  !!
  !> @author David Lauvergnat
  !! @date 25/05/2018
  !!
  !! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
  !! @param d3                 real:                  d3=S%d3, d3 is NOT allocatable
  SUBROUTINE sub_get_d3_FROM_dnS(d3,S)
    USE mod_NumParameters

      real (kind=Rkind),   intent(inout) :: d3(:,:,:)
      TYPE (dnS),        intent(in)    :: S

      character (len=*), parameter :: name_sub='sub_get_d3_FROM_dnS'

      IF (allocated(S%d3)) d3(:,:,:) = S%d3

    END SUBROUTINE sub_get_d3_FROM_dnS
!> @brief Public subroutine which prints a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                  TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
!! @param nio                integer (optional):  when present unit to print S, otherwise it is the default unit:out_unitp
  SUBROUTINE Write_dnS(S,nio,info)
    USE mod_Lib

    TYPE(dnS),        intent(in)           :: S
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

  END SUBROUTINE Write_dnS
  SUBROUTINE WriteAll_dnS(S,nio)
    USE mod_Lib

    TYPE(dnS),      intent(in)           :: S
    integer,          intent(in), optional :: nio

    integer :: i,j,k,nio_loc,ndim

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    write(nio_loc,*) '-------------------------------------------'
    write(nio_loc,*) 'WriteAll_dnS'

    write(nio_loc,*) 'nderiv',S%nderiv

    write(nio_loc,*) 'S%d0',S%d0

    IF (allocated(S%d1)) THEN
      write(nio_loc,*) 'S%d1:',S%d1
    ELSE
      write(nio_loc,*) 'S%d1: not allocated'
    END IF

    IF (allocated(S%d2)) THEN
      write(nio_loc,*) 'S%d2:',S%d2
    ELSE
      write(nio_loc,*) 'S%d2: not allocated'
    END IF

    IF (allocated(S%d3)) THEN
      write(nio_loc,*) 'S%d3:',S%d3
    ELSE
      write(nio_loc,*) 'S%d3: not allocated'
    END IF

    write(nio_loc,*) 'END WriteAll_dnS'
    write(nio_loc,*) '-------------------------------------------'

  END SUBROUTINE WriteAll_dnS
  SUBROUTINE Write_dnS_FOR_test(S,nio,info)
    USE mod_Lib

    TYPE(dnS),        intent(in)           :: S
    integer,            intent(in), optional :: nio
    character(len=*),   intent(in), optional :: info

    integer :: nio_loc,ndim

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    IF (present(info)) THEN
      write(nio_loc,*) 'TEST dnSca: ',trim(adjustl(info))
    ELSE
      write(nio_loc,*) 'TEST dnSca: '
    END IF

    write(nio_loc,*) 'S%d0'
    write(nio_loc,*) 1
    write(nio_loc,*) S%d0

    IF (allocated(S%d1)) THEN
      write(nio_loc,*) 'S%d1'
      write(nio_loc,*) size(S%d1)
      write(nio_loc,*) S%d1
    END IF

    IF (allocated(S%d2)) THEN
      write(nio_loc,*) 'S%d2'
      write(nio_loc,*) size(S%d2)
      write(nio_loc,*) S%d2
    END IF

    IF (allocated(S%d3)) THEN
      write(nio_loc,*) 'S%d3'
      write(nio_loc,*) size(S%d3)
      write(nio_loc,*) S%d3
    END IF

    write(nio_loc,*) 'END_TEST dnSca: '

  END SUBROUTINE Write_dnS_FOR_test
!> @brief Public function to get nderiv from a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                       TYPE(dnS):         derived type which deals with the derivatives of a scalar functions.
!! @param get_nderiv_FROM_dnS   integer  (result):   nderiv value, check against S%nederiv and the allocated d1,d2 or d3
  ELEMENTAL FUNCTION get_nderiv_FROM_dnS(S)
    USE mod_NumParameters

    integer :: get_nderiv_FROM_dnS

    TYPE(dnS), intent(in)    :: S

    IF (.NOT. allocated(S%d1)) THEN
      get_nderiv_FROM_dnS = 0
    ELSE IF (.NOT. allocated(S%d2)) THEN
      get_nderiv_FROM_dnS = 1
    ELSE IF (.NOT. allocated(S%d3)) THEN
      get_nderiv_FROM_dnS = 2
    ELSE
      get_nderiv_FROM_dnS = 3
    END IF

    END FUNCTION get_nderiv_FROM_dnS
!> @brief Public function to get ndim from a derived type dnS.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param S                     TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
!! @param get_ndim_FROM_dnS     integer  (result):   ndim value from the size of S%d1.
  ELEMENTAL FUNCTION get_ndim_FROM_dnS(S)
    USE mod_NumParameters

    integer :: get_ndim_FROM_dnS

    TYPE(dnS), intent(in)    :: S

    IF (.NOT. allocated(S%d1)) THEN
      get_ndim_FROM_dnS = 0
    ELSE
      get_ndim_FROM_dnS = size(S%d1(:))
    END IF

  END FUNCTION get_ndim_FROM_dnS
!> @brief Public function which ckecks a derived type dnS is zero (all components).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Check_dnS_IS_ZERO   logical  (result):   result of the comparison
!! @param S                     TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
!! @param epsi                  real (optional):     when present zero limit, otherwise 10^-10
  ELEMENTAL FUNCTION Check_dnS_IS_ZERO(S,epsi)
    USE mod_NumParameters
    logical :: Check_dnS_IS_ZERO

    TYPE(dnS),        intent(in)           :: S
    real(kind=Rkind),   intent(in), optional :: epsi

    real(kind=Rkind) :: epsi_loc


    IF (present(epsi)) THEN
      epsi_loc = epsi
    ELSE
      epsi_loc = ONETENTH**10
    END IF


    Check_dnS_IS_ZERO = (get_maxval_OF_dnS(S) <= epsi_loc)

  END FUNCTION Check_dnS_IS_ZERO
!> @brief Public function which gets the largest value of a derived type dnS (all components).
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param get_maxval_OF_dnS     real  (result):      largest value
!! @param S                     TYPE(dnS):           derived type which deals with the derivatives of a scalar functions.
  ELEMENTAL FUNCTION get_maxval_OF_dnS(S)
    USE mod_NumParameters

    real(kind=Rkind) :: get_maxval_OF_dnS

    TYPE(dnS), intent(in)    :: S

    real(kind=Rkind) :: e0,e1,e2,e3

    e1 = ZERO
    e2 = ZERO
    e3 = ZERO
    e0 = abs(S%d0)
    IF (allocated(S%d1)) e1 = maxval(abs(S%d1))
    IF (allocated(S%d2)) e2 = maxval(abs(S%d2))
    IF (allocated(S%d3)) e3 = maxval(abs(S%d3))

    get_maxval_OF_dnS = max(e0,e1,e2,e3)

  END FUNCTION get_maxval_OF_dnS
!> @brief Public function which calculates numerical derivative of a function
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!

!! @param Snum                     TYPE(dnS):       function value, f(x), and derivatives f' f" f'".
!! @param f                        real:            function (intrinsic or external)
!! @param x                        real:            abciss
!! @param nderiv                   integer:         order of the derivative
  FUNCTION get_Num_dnS_FROM_f_x(x,f,nderiv) RESULT(Snum)
    USE mod_NumParameters

    TYPE (dnS)                                   :: Snum
    real (kind=Rkind), external                  :: f ! an intrinsic function: sin exp ....
    real (kind=Rkind),           intent(in)      :: x
    integer,           optional, intent(in)      :: nderiv

    !local variables:
    TYPE (dnS)            :: Sloc
    real (kind=Rkind)     :: xloc,step=ONETENTH**4
    integer               :: nderiv_loc,i,j,k
    real (kind=Rkind)     :: f0,fp,fm,fpp,fmm,fppp,fmmm


    character (len=*), parameter :: name_sub='get_Num_dnS_FROM_f_x'

    ! test nderiv
    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(3,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF

    CALL alloc_dnS(Snum,ndim=1,nderiv=nderiv_loc)
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

  END FUNCTION get_Num_dnS_FROM_f_x

!=========================================================
! operators ==,/=,>=,>,<=,<
!=========================================================
  ELEMENTAL FUNCTION dnS_EQ_dnS(S1,S2) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_EQ_dnS'

    lres = (S1%d0 == S2%d0)

    IF (allocated(S1%d1) .AND. allocated(S2%d1)) lres = lres .AND. all(S1%d1 == S2%d1)
    IF (allocated(S1%d2) .AND. allocated(S2%d2)) lres = lres .AND. all(S1%d2 == S2%d2)
    IF (allocated(S1%d3) .AND. allocated(S2%d3)) lres = lres .AND. all(S1%d3 == S2%d3)

  END FUNCTION dnS_EQ_dnS
  ELEMENTAL FUNCTION dnS_EQ_R(S1,R) RESULT(lres)
    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_EQ_R'

    lres = (S1%d0 == R)

  END FUNCTION dnS_EQ_R
  ELEMENTAL FUNCTION R_EQ_dnS(R,S1) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='R_EQ_dnS'

    lres = (R == S1%d0)

  END FUNCTION R_EQ_dnS

  ELEMENTAL FUNCTION dnS_NEQ_dnS(S1,S2) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_NEQ_dnS'

    lres = (S1%d0 /= S2%d0)

    IF (allocated(S1%d1) .AND. allocated(S2%d1)) lres = lres .OR. all(S1%d1 /= S2%d1)
    IF (allocated(S1%d2) .AND. allocated(S2%d2)) lres = lres .OR. all(S1%d2 /= S2%d2)
    IF (allocated(S1%d3) .AND. allocated(S2%d3)) lres = lres .OR. all(S1%d3 /= S2%d3)

  END FUNCTION dnS_NEQ_dnS
  ELEMENTAL FUNCTION dnS_NEQ_R(S1,R) RESULT(lres)
    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_NEQ_R'

    lres = (S1%d0 /= R)

  END FUNCTION dnS_NEQ_R
  ELEMENTAL FUNCTION R_NEQ_dnS(R,S1) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='R_NEQ_dnS'

    lres = (R /= S1%d0)

  END FUNCTION R_NEQ_dnS

  ELEMENTAL FUNCTION dnS_LE_dnS(S1,S2) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_LE_dnS'

    lres = (S1%d0 <= S2%d0)

  END FUNCTION dnS_LE_dnS
  ELEMENTAL FUNCTION dnS_LE_R(S1,R) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_LE_R'

    lres = (S1%d0 <= R)

  END FUNCTION dnS_LE_R
  ELEMENTAL FUNCTION R_LE_dnS(R,S1) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='R_LE_dnS'

    lres = (R <= S1%d0)

  END FUNCTION R_LE_dnS

  ELEMENTAL FUNCTION dnS_LT_dnS(S1,S2) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_LT_dnS'

    lres = (S1%d0 < S2%d0)

  END FUNCTION dnS_LT_dnS
  ELEMENTAL FUNCTION dnS_LT_R(S1,R) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_LT_R'

    lres = (S1%d0 < R)

  END FUNCTION dnS_LT_R
  ELEMENTAL FUNCTION R_LT_dnS(R,S1) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='R_LT_dnS'

    lres = (R < S1%d0)

  END FUNCTION R_LT_dnS

  ELEMENTAL FUNCTION dnS_GE_dnS(S1,S2) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_GE_dnS'

    lres = (S1%d0 >= S2%d0)

  END FUNCTION dnS_GE_dnS
  ELEMENTAL FUNCTION dnS_GE_R(S1,R) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_GE_R'

    lres = (S1%d0 >= R)

  END FUNCTION dnS_GE_R
  ELEMENTAL FUNCTION R_GE_dnS(R,S1) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='R_GE_dnS'

    lres = (R >= S1%d0)

  END FUNCTION R_GE_dnS

  ELEMENTAL FUNCTION dnS_GT_dnS(S1,S2) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS), intent(in)    :: S1,S2
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_GT_dnS'

    lres = (S1%d0 > S2%d0)

  END FUNCTION dnS_GT_dnS
  ELEMENTAL FUNCTION dnS_GT_R(S1,R) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind), intent(in)    :: R
    logical                   :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS_GT_R'

    lres = (S1%d0 > R)

  END FUNCTION dnS_GT_R
  ELEMENTAL FUNCTION R_GT_dnS(R,S1) RESULT(lres)
    USE mod_NumParameters

    TYPE (dnS),        intent(in)    :: S1
    real (kind=Rkind),   intent(in)    :: R
    logical                            :: lres
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='R_GT_dnS'

    lres = (R > S1%d0)

  END FUNCTION R_GT_dnS
!=========================================================
! operators =,+,-,*,/,**
!=========================================================
  ELEMENTAL SUBROUTINE sub_dnS2_TO_dnS1(S1,S2)
    USE mod_NumParameters
    CLASS (dnS), intent(inout) :: S1
    CLASS (dnS), intent(in)    :: S2

    integer :: nderiv_loc,ndim_loc
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='sub_dnS2_TO_dnS1'

    nderiv_loc = get_nderiv_FROM_dnS(S2)
    ndim_loc   = get_ndim_FROM_dnS(S2)

    !write(out_unitp,*) 'ndim,nsurf,nderiv',ndim_loc,nderiv_loc

    CALL dealloc_dnS(S1)
    IF (nderiv_loc < 0 .OR. (nderiv_loc > 0 .AND. ndim_loc < 1)) RETURN

    S1%d0 = S2%d0
    IF (allocated(S2%d1)) S1%d1 = S2%d1
    IF (allocated(S2%d2)) S1%d2 = S2%d2
    IF (allocated(S2%d3)) S1%d3 = S2%d3

  END SUBROUTINE sub_dnS2_TO_dnS1

  ELEMENTAL SUBROUTINE set_dnS_TO_R(S,R)
    USE mod_NumParameters

    CLASS (dnS),       intent(inout) :: S
    real (kind=Rkind), intent(in)    :: R

    integer :: nderiv_loc,ndim_loc
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='set_dnS_TO_R'

    nderiv_loc = get_nderiv_FROM_dnS(S)

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
    END IF
  END SUBROUTINE set_dnS_TO_R
!=========================================================


  ELEMENTAL FUNCTION dnS2_PLUS_dnS1(S1,S2) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                :: Sres
    TYPE (dnS), intent(in)    :: S1,S2

    integer :: nderiv,ndim
    character (len=*), parameter :: name_sub='dnS2_PLUS_dnS1'

    nderiv = min(get_nderiv_FROM_dnS(S1),get_nderiv_FROM_dnS(S2))
    ndim   = min(get_ndim_FROM_dnS(S1),  get_ndim_FROM_dnS(S2))

    CALL dealloc_dnS(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

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
    END IF
  END FUNCTION dnS2_PLUS_dnS1
  ELEMENTAL FUNCTION dnS_PLUS_R(S,R) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R


    integer :: nderiv,ndim
    character (len=*), parameter :: name_sub='dnS_PLUS_R'

    Sres    = S
    Sres%d0 = Sres%d0 + R

  END FUNCTION dnS_PLUS_R
  ELEMENTAL FUNCTION R_PLUS_dnS(R,S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R


    character (len=*), parameter :: name_sub='R_PLUS_dnS'

    Sres    = S
    Sres%d0 = Sres%d0 + R

  END FUNCTION R_PLUS_dnS
  ELEMENTAL FUNCTION PLUS_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S


    character (len=*), parameter :: name_sub='PLUS_dnS'

    Sres = S

  END FUNCTION PLUS_dnS
  ELEMENTAL FUNCTION dnS2_MINUS_dnS1(S1,S2) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                :: Sres
    TYPE (dnS), intent(in)    :: S1,S2

    integer :: nderiv,ndim
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='dnS2_MINUS_dnS1'

    nderiv = min(get_nderiv_FROM_dnS(S1),get_nderiv_FROM_dnS(S2))
    ndim   = min(get_ndim_FROM_dnS(S1),  get_ndim_FROM_dnS(S2))

    CALL dealloc_dnS(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

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
    END IF
  END FUNCTION dnS2_MINUS_dnS1
  ELEMENTAL FUNCTION dnS_MINUS_R(S,R) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R


    character (len=*), parameter :: name_sub='dnS_MINUS_R'

    Sres    = S
    Sres%d0 = Sres%d0 - R

  END FUNCTION dnS_MINUS_R
  ELEMENTAL FUNCTION R_MINUS_dnS(R,S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R


    integer :: nderiv,ndim
    character (len=*), parameter :: name_sub='R_MINUS_dnS'

    nderiv = get_nderiv_FROM_dnS(S)
    ndim   = get_ndim_FROM_dnS(S)

    CALL dealloc_dnS(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

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
    END IF

  END FUNCTION R_MINUS_dnS
  ELEMENTAL FUNCTION MINUS_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S


    integer :: nderiv,ndim
    character (len=*), parameter :: name_sub='MINUS_dnS'

    nderiv = get_nderiv_FROM_dnS(S)
    ndim   = get_ndim_FROM_dnS(S)

    CALL dealloc_dnS(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

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
    END IF

  END FUNCTION MINUS_dnS

  ELEMENTAL FUNCTION dnS2_TIME_dnS1(S1,S2) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                :: Sres
    TYPE (dnS), intent(in)    :: S1,S2

    integer :: nderiv,ndim,id,jd,kd
    character (len=*), parameter :: name_sub='dnS2_TIME_dnS1'

    nderiv = min(get_nderiv_FROM_dnS(S1),get_nderiv_FROM_dnS(S2))
    ndim   = min(get_ndim_FROM_dnS(S1),  get_ndim_FROM_dnS(S2))


    CALL dealloc_dnS(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

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
    END IF
  END FUNCTION dnS2_TIME_dnS1
  ELEMENTAL FUNCTION d0S_TIME_R(S,R) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                         :: Sres
    TYPE (dnS),          intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='d0S_TIME_R'

    nderiv = get_nderiv_FROM_dnS(S)
    ndim   = get_ndim_FROM_dnS(S)

    Sres = S

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN
    Sres%d0 = R * S%d0

  END FUNCTION d0S_TIME_R
  ELEMENTAL FUNCTION dnS_TIME_R(S,R) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind), intent(in)    :: R


    integer :: nderiv,ndim
    character (len=*), parameter :: name_sub='dnS_TIME_R'

    nderiv = get_nderiv_FROM_dnS(S)
    ndim   = get_ndim_FROM_dnS(S)


    CALL dealloc_dnS(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

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
    END IF

  END FUNCTION dnS_TIME_R

  ELEMENTAL FUNCTION R_TIME_dnS(R,S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R


    integer :: nderiv,ndim
    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='R_TIME_dnS'

    nderiv = get_nderiv_FROM_dnS(S)
    ndim   = get_ndim_FROM_dnS(S)


    CALL dealloc_dnS(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

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
    END IF

  END FUNCTION R_TIME_dnS

  ELEMENTAL FUNCTION dnS2_OVER_dnS1(S1,S2) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                :: Sres
    TYPE (dnS), intent(in)    :: S1,S2

    character (len=*), parameter :: name_sub='dnS2_OVER_dnS1'


    Sres = S1 * S2**(-ONE)


  END FUNCTION dnS2_OVER_dnS1
  ELEMENTAL FUNCTION dnS_OVER_R(S,R) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R


    character (len=*), parameter :: name_sub='dnS_OVER_R'


    Sres = S * R**(-ONE)

  END FUNCTION dnS_OVER_R

  ELEMENTAL FUNCTION R_OVER_dnS(R,S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R


    character (len=*), parameter :: name_sub='R_OVER_dnS'


    Sres = R * S**(-ONE)


  END FUNCTION R_OVER_dnS

!=========================================================
! mathematical functions: cos, sin exp, log, cosh ....
! All functions in the fortran norm except atan2 betause it has two arguments
!=========================================================

  ELEMENTAL FUNCTION get_F_dnS(S,d0f,d1f,d2f,d3f) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: d0f,d1f,d2f,d3f

    integer :: nderiv,ndim,id,jd,kd
    character (len=*), parameter :: name_sub='get_F_dnS'

    nderiv = get_nderiv_FROM_dnS(S)
    ndim   = get_ndim_FROM_dnS(S)

    CALL dealloc_dnS(Sres)

    IF (nderiv < 0 .OR. (nderiv > 0 .AND. ndim < 1)) RETURN

    IF (nderiv == 0) THEN
       Sres%d0 = d0f
    ELSE IF (nderiv == 1) THEN
       Sres%d0 =  d0f
       Sres%d1 =  d1f * S%d1
    ELSE IF (nderiv == 2) THEN
       Sres%d0 = d0f
       Sres%d1 = d1f * S%d1
       Sres%d2 = d1f * S%d2

       DO id=1,ndim
       DO jd=1,ndim
         Sres%d2(jd,id) = Sres%d2(jd,id) + d2f * S%d1(id)*S%d1(jd)
       END DO
       END DO
    ELSE IF (nderiv == 3) THEN
       Sres%d0 = d0f
       Sres%d1 = d1f * S%d1
       Sres%d2 = d1f * S%d2
       Sres%d3 = d1f * S%d3

       DO id=1,ndim
       DO jd=1,ndim
         Sres%d2(jd,id) = Sres%d2(jd,id) + d2f * S%d1(id)*S%d1(jd)
       END DO
       END DO

       DO id=1,ndim
       DO jd=1,ndim
       DO kd=1,ndim
         Sres%d3(kd,jd,id) = Sres%d3(kd,jd,id) + &
                             d2f * (S%d1(id)*S%d2(kd,jd) + S%d1(jd)*S%d2(kd,id) + S%d1(kd)*S%d2(jd,id)) + &
                             d3f * S%d1(id)*S%d1(jd)*S%d1(kd)
       END DO
       END DO
       END DO
    END IF

  END FUNCTION get_F_dnS

  ELEMENTAL FUNCTION dnS_EXP_R(S,R) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    real (kind=Rkind),   intent(in)    :: R

    integer :: nderiv
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='dnS_EXP_R'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (R == ZERO) THEN
      Sres = S ! to have the right initialization
      Sres = ONE
    ELSE IF (R == ONE) THEN
      Sres = S
    ELSE IF (R == TWO) THEN
      Sres = S*S
    ELSE
      IF (nderiv >= 0) d0f = S%d0**R
      IF (nderiv >= 1) d1f = S%d0**(R-ONE)   * R
      IF (nderiv >= 2) d2f = S%d0**(R-TWO)   * R*(R-ONE)
      IF (nderiv >= 3) d3f = S%d0**(R-THREE) * R*(R-ONE)*(R-TWO)
      Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)
    END IF


  END FUNCTION dnS_EXP_R
  ELEMENTAL FUNCTION dnS_EXP_I(S,I) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S
    integer,           intent(in)    :: I


    character (len=*), parameter :: name_sub='dnS_EXP_I'

    IF (I == 0) THEN
      Sres = S ! to have the right initialization
      Sres = ONE
    ELSE IF (I == 1) THEN
      Sres = S
    ELSE IF (I == 2) THEN
      Sres = S*S
    ELSE
      Sres = S**real(I,kind=Rkind)
    END IF

  END FUNCTION dnS_EXP_I
  ELEMENTAL FUNCTION get_SQRT_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S


    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='get_SQRT_dnS'

    Sres = S**HALF

  END FUNCTION get_SQRT_dnS

  ELEMENTAL FUNCTION get_ABS_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S


    integer :: err_dnS_loc
    character (len=*), parameter :: name_sub='get_ABS_dnS'

    Sres = (S*S)**HALF

  END FUNCTION get_ABS_dnS

  ELEMENTAL FUNCTION get_EXP_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f

    character (len=*), parameter :: name_sub='get_EXP_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  exp(S%d0)
    IF (nderiv >= 1) d1f =  exp(S%d0)
    IF (nderiv >= 2) d2f =  exp(S%d0)
    IF (nderiv >= 3) d3f =  exp(S%d0)

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_EXP_dnS
  ELEMENTAL FUNCTION get_LOG_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f

    character (len=*), parameter :: name_sub='get_LOG_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  log(S%d0)
    IF (nderiv >= 1) d1f =  ONE/S%d0
    IF (nderiv >= 2) d2f = -ONE/S%d0**2
    IF (nderiv >= 3) d3f =  TWO/S%d0**3

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_LOG_dnS
  ELEMENTAL FUNCTION get_LOG10_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f

    character (len=*), parameter :: name_sub='get_LOG10_dnS'

    Sres = log(S)/log(TEN)

  END FUNCTION get_LOG10_dnS
  ELEMENTAL FUNCTION get_COS_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S


    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_COS_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  cos(S%d0)
    IF (nderiv >= 1) d1f = -sin(S%d0)
    IF (nderiv >= 2) d2f = -cos(S%d0)
    IF (nderiv >= 3) d3f =  sin(S%d0)

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_COS_dnS
  ELEMENTAL FUNCTION get_ACOS_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S


    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ACOS_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  acos(S%d0)
    IF (nderiv >= 1) d1f = -ONE/sqrt(ONE-S%d0**2)
    IF (nderiv >= 2) d2f = S%d0*d1f**3
    IF (nderiv >= 3) d3f = (ONE+TWO*S%d0**2)*d1f**5

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ACOS_dnS
  ELEMENTAL FUNCTION get_SIN_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_SIN_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  sin(S%d0)
    IF (nderiv >= 1) d1f =  cos(S%d0)
    IF (nderiv >= 2) d2f = -sin(S%d0)
    IF (nderiv >= 3) d3f = -cos(S%d0)

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_SIN_dnS
  ELEMENTAL FUNCTION get_ASIN_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ASIN_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  asin(S%d0)
    IF (nderiv >= 1) d1f = ONE/sqrt(ONE-S%d0**2)
    IF (nderiv >= 2) d2f = S%d0*d1f**3
    IF (nderiv >= 3) d3f = (ONE+TWO*S%d0**2)*d1f**5

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ASIN_dnS
  ELEMENTAL FUNCTION get_TAN_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_TAN_dnS'

    Sres = Sin(S)/cos(S)

  END FUNCTION get_TAN_dnS

  ELEMENTAL FUNCTION get_ATAN_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ATAN_dnS'

    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  atan(S%d0)
    IF (nderiv >= 1) d1f = ONE/(ONE+S%d0**2)
    IF (nderiv >= 2) d2f = -TWO*S%d0 * d1f**2
    IF (nderiv >= 3) d3f = (-TWO+SIX*S%d0**2) * d1f**3

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ATAN_dnS

  ELEMENTAL FUNCTION get_COSH_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S


    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_COSH_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  cosh(S%d0)
    IF (nderiv >= 1) d1f =  sinh(S%d0)
    IF (nderiv >= 2) d2f =  cosh(S%d0)
    IF (nderiv >= 3) d3f =  sinh(S%d0)

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_COSH_dnS
  ELEMENTAL FUNCTION get_ACOSH_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S


    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ACOSH_dnS'


    nderiv = get_nderiv_FROM_dnS(S)
#if __INVHYP == 1
    IF (nderiv >= 0) d0f =   acosh(S%d0)
#else
    IF (nderiv >= 0) d0f =  acosh_perso(S%d0)
#endif
    IF (nderiv >= 1) d1f = ONE/sqrt(-ONE+S%d0**2)
    IF (nderiv >= 2) d2f = -S%d0*d1f**3
    IF (nderiv >= 3) d3f = (ONE+TWO*S%d0**2)*d1f**5

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ACOSH_dnS
  ELEMENTAL FUNCTION get_SINH_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_SINH_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  sinh(S%d0)
    IF (nderiv >= 1) d1f =  cosh(S%d0)
    IF (nderiv >= 2) d2f =  sinh(S%d0)
    IF (nderiv >= 3) d3f =  cosh(S%d0)

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_SINH_dnS
  ELEMENTAL FUNCTION get_ASINH_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ASINH_dnS'


    nderiv = get_nderiv_FROM_dnS(S)
#if __INVHYP == 1
    IF (nderiv >= 0) d0f =  asinh(S%d0)
#else
    IF (nderiv >= 0) d0f =  asinh_perso(S%d0)
#endif
    IF (nderiv >= 1) d1f = ONE/sqrt(ONE+S%d0**2)
    IF (nderiv >= 2) d2f = -S%d0*d1f**3
    IF (nderiv >= 3) d3f = (-ONE+TWO*S%d0**2)*d1f**5

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ASINH_dnS

  ELEMENTAL FUNCTION get_TANH_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_TANH_dnS'


    nderiv = get_nderiv_FROM_dnS(S)

    IF (nderiv >= 0) d0f =  tanh(S%d0)
    IF (nderiv >= 1) d1f =  ONE/cosh(S%d0)**2
    IF (nderiv >= 2) d2f = -TWO*tanh(S%d0) * d1f
    IF (nderiv >= 3) d3f = (-FOUR+TWO*cosh(TWO*S%d0)) * d1f**2

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_TANH_dnS

  ELEMENTAL FUNCTION get_ATANH_dnS(S) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: S

    integer :: nderiv
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='get_ATANH_dnS'


    nderiv = get_nderiv_FROM_dnS(S)
#if __INVHYP == 1
    IF (nderiv >= 0) d0f =  atanh(S%d0)
#else
    IF (nderiv >= 0) d0f =  atanh_perso(S%d0)
#endif
    IF (nderiv >= 1) d1f =  ONE/(ONE-S%d0**2)
    IF (nderiv >= 2) d2f =  TWO*S%d0 * d1f**2
    IF (nderiv >= 3) d3f =  (TWO+SIX*S%d0**2) * d1f**3

    Sres = get_F_dnS(S,d0f,d1f,d2f,d3f)

  END FUNCTION get_ATANH_dnS

  ELEMENTAL FUNCTION atanh_perso(x)
    real (kind=Rkind)             :: atanh_perso
    real (kind=Rkind), intent(in) :: x

#if __INVHYP == 1
    atanh_perso = atanh(x)
#else
    atanh_perso = HALF*log((ONE+x)/(ONE-x))
#endif

  END FUNCTION atanh_perso
  ELEMENTAL FUNCTION asinh_perso(x)
    real (kind=Rkind)             :: asinh_perso
    real (kind=Rkind), intent(in) :: x

#if __INVHYP == 1
    asinh_perso = asinh(x)
#else
    asinh_perso = log(x+sqrt(x*x+ONE))
#endif

  END FUNCTION asinh_perso
  ELEMENTAL FUNCTION acosh_perso(x)
    real (kind=Rkind)             :: acosh_perso
    real (kind=Rkind), intent(in) :: x

#if __INVHYP == 1
    acosh_perso = acosh(x)
#else
    acosh_perso = log(x+sqrt(x*x-ONE))
#endif

  END FUNCTION acosh_perso

  FUNCTION dot_product_VecOFdnS(VecA,VecB) RESULT(Sres)
    USE mod_NumParameters

    TYPE (dnS)                       :: Sres
    TYPE (dnS),        intent(in)    :: VecA(:),VecB(:)

    integer :: i
    integer :: err_dnS_loc
    real(kind=Rkind) :: d0f,d1f,d2f,d3f
    character (len=*), parameter :: name_sub='dot_product_VecOFdnS'


    IF (size(VecA) /= size(VecB)) THEN
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) '  size of both vectors are different'
       write(out_unitp,*) '  size(VecA),size(VecB)',size(VecA),size(VecB)
       STOP 'Problem in dot_product_VecOFdnS'
    END IF

    Sres = VecA(lbound(VecA,dim=1)) ! for the initialization
    Sres = ZERO
    DO i=lbound(VecA,dim=1),ubound(VecA,dim=1)
      Sres = Sres + VecA(i) * VecB(lbound(VecB,dim=1)+i-1)
    END DO

  END FUNCTION dot_product_VecOFdnS

END MODULE mod_dnS
