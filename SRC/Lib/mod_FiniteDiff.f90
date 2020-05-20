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
MODULE mod_FiniteDiff
!$ USE omp_lib
  IMPLICIT NONE

  TYPE FiniteDiff_t
    integer             :: nderiv           ! order of the derivatives (0,1,2,3)
    integer             :: ndim             ! dimension of the points
    integer             :: nb_pts           ! number of points
    integer             :: list_pts(:,:)    ! list_pts(ndim,nb_pts)
    real (kind=Rkind)   :: CoefD1(:)        ! coefficient: Coef(nb_pts) for d./dx

    real (kind=Rkind)   :: CoefD2(:)        ! coefficient: Coef(nb_pts) for d2./dx2
    real (kind=Rkind)   :: CoefDD(:)        ! coefficient: Coef(nb_pts) for d2./dxdy

    real (kind=Rkind)   :: CoefD3(:)        ! coefficient: Coef(nb_pts) for d3./dx3
    real (kind=Rkind)   :: CoefD2D(:)       ! coefficient: Coef(nb_pts) for d3./dx2dy
    real (kind=Rkind)   :: CoefDDD(:)       ! coefficient: Coef(nb_pts) for d3./dxdydz

  END TYPE FiniteDiff_t

  TYPE PrimList_t
    integer             :: nderiv           ! order of the derivatives (0,1,2,3)
    integer             :: ndim             ! dimension of the points
    integer             :: nb_pts           ! number of points
    integer             :: list_pts(:,:)    ! list_pts(ndim,nb_pts)
    real (kind=Rkind)   :: CoefD1(:)        ! coefficient: Coef(nb_pts) for d./dx
    real (kind=Rkind)   :: CoefD2(:)        ! coefficient: Coef(nb_pts) for d2./dx2
    real (kind=Rkind)   :: CoefD3(:)        ! coefficient: Coef(nb_pts) for d3./dx3
  END TYPE PrimList_t

  TYPE (PrimList_t), PRIVATE :: PrimList_D   ! along 1 coordinate
  TYPE (PrimList_t), PRIVATE :: PrimList_DD  ! along 2 coordinates
  TYPE (PrimList_t), PRIVATE :: PrimList_DDD ! along 3 coordinates

  CONTAINS

  SUBROUTINE Init_FiniteDiff_D(nderiv,step,option)
    USE mod_NumParameters
    IMPLICIT NONE

    integer,             intent(in)           :: nderiv
    real (kind=Rkind),   intent(in)           :: step
    integer,             intent(in)           :: option


    integer           :: ndim

    ndim = 1

    ! option = 2,3,5 (number of points)
    IF (nderiv == 3 .AND. option /= 5) THEN
      write(out_unitp,*) ' ERROR in Init_FiniteDiff_D'
      write(out_unitp,*) ' nderiv,option',nderiv,option
      write(out_unitp,*) '  => For nderiv=3, you MUST use option=5'
      STOP ' ERROR in Init_FiniteDiff_D: incompatible nderiv and option'
    END IF
    IF (nderiv == 2 .AND. option == 2) THEN
      write(out_unitp,*) ' ERROR in Init_FiniteDiff_D'
      write(out_unitp,*) ' nderiv,option',nderiv,option
      write(out_unitp,*) '  => For nderiv=2, you MUST use option=3 or 5'
      STOP ' ERROR in Init_FiniteDiff_D: incompatible nderiv and option'
    END IF
    PrimList_D%ndim   = ndim
    PrimList_D%nderiv = nderiv
      SELECT CASE (option)
      CASE (2)
         PrimList_D%nb_pts = 3
         allocate(PrimList_D%list_pts(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD1(ndim,PrimList_D%nb_pts))

         PrimList_D%list_pts(1,1) =  0
         PrimList_D%CoefD1(1)     = ZERO

         PrimList_D%list_pts(1,2) =  1
         PrimList_D%CoefD1(1)     = HALF

        PrimList_D%CoefD1(:) = PrimList_D%CoefD1/step
      CASE (3)
         PrimList_D%nb_pts = 3
         allocate(PrimList_D%list_pts(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD1(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD2(ndim,PrimList_D%nb_pts))

         PrimList_D%list_pts(1,1) = -1
         PrimList_D%CoefD1(1)     = HALF
         PrimList_D%CoefD2(1)     = ONE

         PrimList_D%list_pts(1,2) =  0
         PrimList_D%CoefD1(2)     = ZERO
         PrimList_D%CoefD2(2)     = -TWO

         PrimList_D%list_pts(1,3) =  1
         PrimList_D%CoefD1(3)     = HALF
         PrimList_D%CoefD2(3)     = ONE

        PrimList_D%CoefD1(:) = PrimList_D%CoefD1/step
        PrimList_D%CoefD2(:) = PrimList_D%CoefD2/step**2

      CASE (5)
         PrimList_D%nb_pts = 5
         allocate(PrimList_D%list_pts(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD1(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD2(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD3(ndim,PrimList_D%nb_pts))

! d/dx:   FullSimplify[(8/12 (Dfp00 - Dfm00) - 1/12 (Dfp200 - Dfm200))/h]
! d2/dx2: FullSimplify[(12/9 (Dfp00 + Dfm00) - 1/12 (Dfp200 + Dfm200))/h^2]
! d3/dx3: FullSimplify[(-(Dfp00 - Dfm00) + 1/2 (Dfp200 - Dfm200))/h^3]

         i = 1
         PrimList_D%list_pts(1,i) = -2
         PrimList_D%CoefD1(i)     =  ONE/TWELVE
         PrimList_D%CoefD2(i)     = -ONE/TWELVE
         PrimList_D%CoefD3(i)     = -HALF

         i = 2
         PrimList_D%list_pts(1,i) = -1
         PrimList_D%CoefD1(i)     = -EIGHT/TWELVE
         PrimList_D%CoefD2(i)     =  TWELVE/NINE
         PrimList_D%CoefD3(i)     =  ONE

         i = 3
         PrimList_D%list_pts(1,i) =  0
         PrimList_D%CoefD1(i)     =  ZERO
         PrimList_D%CoefD2(i)     =  FIVE/TWO
         PrimList_D%CoefD3(i)     =  ZERO

         i = 4
         PrimList_D%list_pts(1,i) =  1
         PrimList_D%CoefD1(i)     =  EIGHT/TWELVE
         PrimList_D%CoefD2(i)     =  TWELVE/NINE
         PrimList_D%CoefD3(i)     = -ONE

         i = 5
         PrimList_D%list_pts(1,i) =  2
         PrimList_D%CoefD1(i)     = -ONE/TWELVE
         PrimList_D%CoefD2(i)     = -ONE/TWELVE
         PrimList_D%CoefD3(i)     =  HALF

        PrimList_D%CoefD1(:) = PrimList_D%CoefD1/step
        PrimList_D%CoefD2(:) = PrimList_D%CoefD2/step**2
        PrimList_D%CoefD3(:) = PrimList_D%CoefD3/step**3
      CASE Default
      END SELECT

  END SUBROUTINE Init_FiniteDiff_D

  SUBROUTINE Init_FiniteDiff_DD(nderiv,step,option)
    USE mod_NumParameters
    IMPLICIT NONE

    integer,             intent(in)           :: nderiv
    real (kind=Rkind),   intent(in)           :: step
    integer,             intent(in)           :: option


    integer           :: ndim

    ndim = 1

    ! option = 2,3,5 (number of points)
    IF (nderiv == 3 .AND. option /= 5) THEN
      write(out_unitp,*) ' ERROR in Init_FiniteDiff_DD'
      write(out_unitp,*) ' nderiv,option',nderiv,option
      write(out_unitp,*) '  => For nderiv=3, you MUST use option=5'
      STOP ' ERROR in Init_FiniteDiff_D: incompatible nderiv and option'
    END IF

    PrimList_D%ndim   = ndim
    PrimList_D%nderiv = nderiv
    IF (ndim == 1) THEN
      SELECT CASE (nderiv)
      CASE (1)
         PrimList_D%nb_pts = 3
         allocate(PrimList_D%list_pts(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD1(ndim,PrimList_D%nb_pts))

         PrimList_D%list_pts(1,1) = -1
         PrimList_D%CoefD1(1)     = HALF

         PrimList_D%list_pts(1,2) =  0
         PrimList_D%CoefD1(2)     = ZERO

         PrimList_D%list_pts(1,3) =  1
         PrimList_D%CoefD1(3)     = HALF

        PrimList_D%CoefD1(:) = PrimList_D%CoefD1/step
      CASE (2)
         PrimList_D%nb_pts = 3
         allocate(PrimList_D%list_pts(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD1(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD2(ndim,PrimList_D%nb_pts))

         PrimList_D%list_pts(1,1) = -1
         PrimList_D%CoefD1(1)     = HALF
         PrimList_D%CoefD2(1)     = ONE

         PrimList_D%list_pts(1,2) =  0
         PrimList_D%CoefD1(2)     = ZERO
         PrimList_D%CoefD2(2)     = -TWO

         PrimList_D%list_pts(1,3) =  1
         PrimList_D%CoefD1(3)     = HALF
         PrimList_D%CoefD2(3)     = ONE

        PrimList_D%CoefD1(:) = PrimList_D%CoefD1/step
        PrimList_D%CoefD2(:) = PrimList_D%CoefD2/step**2

      CASE (3)
         PrimList_D%nb_pts = 5
         allocate(PrimList_D%list_pts(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD1(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD2(ndim,PrimList_D%nb_pts))
         allocate(PrimList_D%CoefD3(ndim,PrimList_D%nb_pts))

! d/dx:   FullSimplify[(8/12 (Dfp00 - Dfm00) - 1/12 (Dfp200 - Dfm200))/h]
! d2/dx2: FullSimplify[(12/9 (Dfp00 + Dfm00) - 1/12 (Dfp200 + Dfm200))/h^2]
! d3/dx3: FullSimplify[(-(Dfp00 - Dfm00) + 1/2 (Dfp200 - Dfm200))/h^3]

         i = 1
         PrimList_D%list_pts(1,i) = -2
         PrimList_D%CoefD1(i)     =  ONE/TWELVE
         PrimList_D%CoefD2(i)     = -ONE/TWELVE
         PrimList_D%CoefD3(i)     = -HALF

         i = 2
         PrimList_D%list_pts(1,i) = -1
         PrimList_D%CoefD1(i)     = -EIGHT/TWELVE
         PrimList_D%CoefD2(i)     =  TWELVE/NINE
         PrimList_D%CoefD3(i)     =  ONE

         i = 3
         PrimList_D%list_pts(1,i) =  0
         PrimList_D%CoefD1(i)     =  ZERO
         PrimList_D%CoefD2(i)     =  FIVE/TWO
         PrimList_D%CoefD3(i)     =  ZERO

         i = 4
         PrimList_D%list_pts(1,i) =  1
         PrimList_D%CoefD1(i)     =  EIGHT/TWELVE
         PrimList_D%CoefD2(i)     =  TWELVE/NINE
         PrimList_D%CoefD3(i)     = -ONE

         i = 5
         PrimList_D%list_pts(1,i) =  2
         PrimList_D%CoefD1(i)     = -ONE/TWELVE
         PrimList_D%CoefD2(i)     = -ONE/TWELVE
         PrimList_D%CoefD3(i)     =  HALF

        PrimList_D%CoefD1(:) = PrimList_D%CoefD1/step
        PrimList_D%CoefD2(:) = PrimList_D%CoefD2/step**2
        PrimList_D%CoefD3(:) = PrimList_D%CoefD3/step**3
      CASE Default
      END SELECT
    END IF

  END SUBROUTINE Init_FiniteDiff_DD

END MODULE mod_FiniteDiff
