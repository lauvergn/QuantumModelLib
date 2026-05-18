!===========================================================================
!===========================================================================
!This file is part of QuantumModelLib (QML).
!===============================================================================
! MIT License
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!    Copyright (c) 2022 David Lauvergnat [1]
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
  USE QDUtil_NumParameters_m
  USE QDUtil_m, ONLY : TO_string
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind)              :: a
  real (kind=Rkind), allocatable :: Qmin(:),Qmax(:)
  character (len=:), allocatable :: grid_file
  integer                        :: i,j

  ! CALL Init_Model(QModel,pot_name='H3',ndim=2)
  ! Qmin = [0.7_Rkind,0.7_Rkind]
  ! Qmax = [FIVE,FIVE]

  ! grid_file = 'grid2D_H3'

  ! CALL Eval_pot_ON_Grid(QModel,Qmin=Qmin,Qmax=Qmax,       &
  !                       nb_points=41,nderiv=0,grid_file=grid_file)

  ! deallocate(Qmin)
  ! deallocate(Qmax)
  ! deallocate(grid_file)

  ! STOP

  a=SIX

  CALL Init_Model(QModel,pot_name='HenonHeiles',ndim=6)

  allocate(Qmin(QModel%QM%ndim))
  allocate(Qmax(QModel%QM%ndim))

  DO i=1,QModel%QM%ndim
    Qmax = ZERO
    Qmax(i) = a
    Qmin = -Qmax
    grid_file = 'grid1D_' // TO_string(i)
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

    grid_file = 'grid2D_' // TO_string(i) // '-' // TO_string(j)

    CALL Eval_pot_ON_Grid(QModel,Qmin=Qmin,Qmax=Qmax,       &
                          nb_points=101,nderiv=0,grid_file=grid_file)
  END DO
  END DO

  deallocate(Qmin)
  deallocate(Qmax)
  deallocate(grid_file)




END PROGRAM TEST_grid
