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
!    Copyright 2017 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================
PROGRAM TEST_grid
  USE mod_NumParameters
  USE mod_Lib
  USE mod_Model
  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind)              :: a
  real (kind=Rkind), allocatable :: Qmin(:),Qmax(:)
  character (len=:), allocatable :: grid_file
  integer                        :: i,j

  a=SIX

  CALL Init_Model(Para_Model,pot_name='HenonHeiles',ndim=6)

  allocate(Qmin(Para_Model%ndim))
  allocate(Qmax(Para_Model%ndim))

  DO i=1,Para_Model%ndim
    Qmax = ZERO
    Qmax(i) = a
    Qmin = -Qmax
    grid_file = strdup('grid1D_' // int_TO_char(i) )
    CALL Eval_pot_ON_Grid(Para_Model,Qmin=Qmin,Qmax=Qmax,       &
                          nb_points=101,nderiv=0,grid_file=grid_file)
  END DO

  IF (Para_Model%ndim == 1) STOP '1D'

  DO i=1,Para_Model%ndim
  DO j=i+1,Para_Model%ndim
    Qmax = ZERO
    Qmax(i) = a
    Qmax(j) = a

    Qmin = -Qmax

    grid_file = strdup('grid2D_' // int_TO_char(i) // '-' // int_TO_char(j))

    CALL Eval_pot_ON_Grid(Para_Model,Qmin=Qmin,Qmax=Qmax,       &
                          nb_points=101,nderiv=0,grid_file=grid_file)
  END DO
  END DO

  deallocate(Qmin)
  deallocate(Qmax)
  deallocate(grid_file)


END PROGRAM TEST_grid
