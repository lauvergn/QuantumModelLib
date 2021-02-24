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
!!  TYPE (dnS_t) :: S
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
!! @date 26/04/2020
!!
MODULE mod_dnPoly
  USE mod_QML_NumParameters
  USE mod_dnS
  IMPLICIT NONE


CONTAINS
  ELEMENTAL FUNCTION QML_dnMonomial(x,i) RESULT(Sres)
    USE mod_QML_NumParameters

    TYPE (dnS_t)                       :: Sres
    TYPE (dnS_t),        intent(in)    :: x
    integer,             intent(in)    :: i


    character (len=*), parameter :: name_sub='QML_dnMonomial'

    Sres = x**i

  END FUNCTION QML_dnMonomial
  ELEMENTAL FUNCTION QML_dnBox(x,i) RESULT(Sres)
    USE mod_QML_NumParameters

    TYPE (dnS_t)                       :: Sres
    TYPE (dnS_t),        intent(in)    :: x
    integer,             intent(in)    :: i

    character (len=*), parameter :: name_sub='QML_dnBox'

    Sres = sin(x*real(i,kind=Rkind)) / sqrt(pi*HALF)

  END FUNCTION QML_dnBox
  ELEMENTAL FUNCTION QML_dnFourier(x,i) RESULT(Sres)
    USE mod_QML_NumParameters

    TYPE (dnS_t)                       :: Sres
    TYPE (dnS_t),        intent(in)    :: x
    integer,             intent(in)    :: i

    real(kind=Rkind), parameter :: sqpi  = ONE/sqrt(pi)
    real(kind=Rkind), parameter :: sq2pi = ONE/sqrt(pi+pi)

    integer          :: ii
    TYPE (dnS_t)     :: xx
    character (len=*), parameter :: name_sub='QML_dnFourier'


    ii = int(i/2)
    xx = x*real(ii,kind=Rkind)

    IF (ii == 0) THEN
      Sres = sq2pi
    ELSE IF (mod(i,2) == 0) THEN
      Sres = sin(xx) * sqpi
    ELSE
      Sres = cos(xx) * sqpi
    END IF


  END FUNCTION QML_dnFourier
  ELEMENTAL FUNCTION QML_dnLegendre0(x,i) RESULT(Sres)
    USE mod_QML_NumParameters

    TYPE (dnS_t)                       :: Sres
    TYPE (dnS_t),        intent(in)    :: x
    integer,             intent(in)    :: i

    TYPE (dnS_t)     :: P2,P1,P0
    integer          :: j
    real(kind=Rkind) :: Pnorm2

    character (len=*), parameter :: name_sub='QML_dnLegendre0'

    IF ( i<= 0) THEN
      Sres = ONE
    ELSE IF ( i== 1) THEN
      Sres = x
    ELSE
      P0 = ONE
      P1 = x
      DO j=2,i
        P2 = (real(2*j-1,kind=Rkind)*x*P1 -real(j-1,kind=Rkind)*P0 )/real(j,kind=Rkind)
        P0 = P1
        P1 = P2
      END DO
      Sres = P2
    END IF
    Pnorm2 = TWO/real(2*i+1,kind=Rkind)
    Sres = Sres/sqrt(Pnorm2)

  END FUNCTION QML_dnLegendre0
!===================================================
!
!   Normalized Hermite polynomial Hm(x,l)
!   with x in ( -inf =< x =< inf )
!
!===================================================
  ELEMENTAL FUNCTION QML_dnHermite(x,l)  RESULT(Sres)
    USE mod_QML_NumParameters

    TYPE (dnS_t)                       :: Sres
    TYPE (dnS_t),        intent(in)    :: x
    integer,             intent(in)    :: l


    ! Polynomial for  l, l-1 et l-2
    TYPE (dnS_t)      :: pl0,pl1,pl2
    real (kind=Rkind) :: norm
    integer           :: i

    norm=sqrt(pi)

    IF (l == 0) THEN
       Sres = ONE/sqrt(norm)
    ELSE IF (l == 1) THEN
       Sres = TWO*x/sqrt(TWO*norm)
    ELSE

     pl2  = ONE
     pl1  = TWO*x
     norm = norm*TWO

     DO i=2,l
       norm = norm*TWO*real(i,kind=Rkind)
       pl0  = TWO*( x*pl1 - real(i-1,kind=Rkind)*pl2 )
       pl2  = pl1
       pl1  = pl0
     END DO
     Sres = pl0/sqrt(norm)
   END IF

  END FUNCTION QML_dnHermite
!===================================================
!
!   Normalized Hermite polynomial with gaussian part:Hm(x,l).exp(-x**2/2)
!   with x in ( -inf =< x =< inf )
!
!===================================================
  ELEMENTAL FUNCTION QML_dnExpHermite(x,l)  RESULT(Sres)
    USE mod_QML_NumParameters

    TYPE (dnS_t)                       :: Sres
    TYPE (dnS_t),        intent(in)    :: x
    integer,             intent(in)    :: l

    Sres = QML_dnHermite(x,l) * exp(-x*x*HALF)

 END FUNCTION QML_dnExpHermite
END MODULE mod_dnPoly
