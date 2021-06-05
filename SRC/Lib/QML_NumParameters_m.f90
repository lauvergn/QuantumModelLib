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
  MODULE QML_NumParameters_m
!$ USE omp_lib
      USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64,real128,int32,int64
      IMPLICIT NONE

      PUBLIC
      PRIVATE :: INPUT_UNIT,OUTPUT_UNIT,real64,real128,int32,int64

      integer, parameter :: Rkind        = real64 ! 8
      integer, parameter :: Ikind        = int32  ! 4
      integer, parameter :: ILkind       = int64  ! 8
      integer, parameter :: Name_len     = 20
      integer, parameter :: Name_longlen = 50
      integer, parameter :: Line_len     = 255
      integer, parameter :: error_l      = 80


      real (kind=Rkind), parameter :: ZERO    = 0._Rkind
      real (kind=Rkind), parameter :: ONE     = 1._Rkind
      real (kind=Rkind), parameter :: TWO     = 2._Rkind
      real (kind=Rkind), parameter :: THREE   = 3._Rkind
      real (kind=Rkind), parameter :: FOUR    = 4._Rkind
      real (kind=Rkind), parameter :: FIVE    = 5._Rkind
      real (kind=Rkind), parameter :: SIX     = 6._Rkind
      real (kind=Rkind), parameter :: SEVEN   = 7._Rkind
      real (kind=Rkind), parameter :: EIGHT   = 8._Rkind
      real (kind=Rkind), parameter :: NINE    = 9._Rkind
      real (kind=Rkind), parameter :: TEN     = 10._Rkind
      real (kind=Rkind), parameter :: ELEVEN  = 11._Rkind
      real (kind=Rkind), parameter :: TWELVE  = 12._Rkind
      real (kind=Rkind), parameter :: HUNDRED = 100._Rkind

      real (kind=Rkind), parameter :: HALF      = ONE/TWO
      real (kind=Rkind), parameter :: THIRD     = ONE/THREE
      real (kind=Rkind), parameter :: FOURTH    = ONE/FOUR
      real (kind=Rkind), parameter :: QUARTER   = ONE/FOUR
      real (kind=Rkind), parameter :: FIFTH     = ONE/FIVE
      real (kind=Rkind), parameter :: SIXTH     = ONE/SIX
      real (kind=Rkind), parameter :: ONETENTH  = ONE/TEN
      real (kind=Rkind), parameter :: TWOTENTHS = TWO/TEN

      real (kind=Rkind), parameter ::                                   &
       pi = 3.14159265358979323846264338327950288419716939937511_Rkind

      complex (kind=Rkind), parameter :: EYE      = (0._Rkind,1._Rkind)
      complex (kind=Rkind), parameter :: CZERO    = (0._Rkind,0._Rkind)
      complex (kind=Rkind), parameter :: CONE     = (1._Rkind,0._Rkind)
      complex (kind=Rkind), parameter :: CHALF    = (0.5_Rkind,0._Rkind)

      integer :: print_level = 0        ! 0 minimal, 1 default, 2 large, -1 nothing

      character (len=Name_longlen) :: EneIO_format = "f20.5"
      character (len=Name_longlen) :: RMatIO_format = "f18.10"
      character (len=Name_longlen) :: CMatIO_format = "'(',f15.7,' +i',f15.7,')'"

      !integer :: in_unitp  = 5 ! Unit for input and the ouptput files
      !integer :: out_unitp = 6 ! Unit for input and the ouptput files
      integer :: in_unitp  = INPUT_UNIT  ! Unit for the ouptput files, with the ISO_FORTRAN_ENV
      integer :: out_unitp = OUTPUT_UNIT ! Unit for the input, with the ISO_FORTRAN_ENV

  END MODULE QML_NumParameters_m
