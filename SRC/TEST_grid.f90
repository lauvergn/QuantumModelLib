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
!    Copyright 2016 David Lauvergnat [1]
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
PROGRAM TEST_grid
  USE QML_NumParameters_m
  USE QML_UtilLib_m
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind)              :: a
  real (kind=Rkind), allocatable :: Qmin(:),Qmax(:)
  character (len=:), allocatable :: grid_file
  integer                        :: i,j

  a=SIX

  CALL Init_Model(QModel,pot_name='HenonHeiles',ndim=6)

  allocate(Qmin(QModel%QM%ndim))
  allocate(Qmax(QModel%QM%ndim))

  DO i=1,QModel%QM%ndim
    Qmax = ZERO
    Qmax(i) = a
    Qmin = -Qmax
    grid_file = strdup('grid1D_' // int_TO_char(i) )
    CALL Eval_pot_ON_Grid(QModel,Qmin=Qmin,Qmax=Qmax,       &
                          nb_points=101,nderiv=0,grid_file=grid_file)
  END DO

  IF (QModel%QM%ndim == 1) STOP '1D'

  DO i=1,QModel%QM%ndim
  DO j=i+1,QModel%QM%ndim
    Qmax = ZERO
    Qmax(i) = a
    Qmax(j) = a

    Qmin = -Qmax

    grid_file = strdup('grid2D_' // int_TO_char(i) // '-' // int_TO_char(j))

    CALL Eval_pot_ON_Grid(QModel,Qmin=Qmin,Qmax=Qmax,       &
                          nb_points=101,nderiv=0,grid_file=grid_file)
  END DO
  END DO

  deallocate(Qmin)
  deallocate(Qmax)
  deallocate(grid_file)


END PROGRAM TEST_grid
