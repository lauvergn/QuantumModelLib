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
PROGRAM TEST_HONO
  USE QDUtil_NumParameters_m
  USE QDUtil_Test_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: Model
  real (kind=Rkind), allocatable :: Q(:),Qref(:),xyz(:)
  !real (kind=Rkind), allocatable :: G(:,:),Gref(:,:) ! metric tensor
  integer                        :: ndim,nsurf,nderiv,i,option,err
  logical                        :: Lerr
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotValref
  TYPE (dnMat_t)                 :: dnErr
  TYPE (test_t)                  :: test_var
  real (kind=Rkind), parameter   :: epsi = 1.e-10_Rkind
  real (kind=Rkind), parameter   :: a0 = 0.52917720835354106_Rkind ! from Tnum


  CALL Initialize_Test(test_var,test_name='QModel_HONO')

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' HONO potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(Model,pot_name='HONO',read_param=.FALSE.)
  Q = [2.696732586_Rkind,1.822912197_Rkind,1.777642018_Rkind,2.213326419_Rkind,1.9315017_Rkind,pi]

  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,'(a,6f12.6)') 'Q (Bohr)',Q(:)
  CALL Check_analytical_numerical_derivatives(Model,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives (option=0)'
  Q = [2.696732586_Rkind,1.822912197_Rkind,1.777642018_Rkind,2.213326419_Rkind,1.9315017_Rkind,pi]

  CALL Test_QVG_FOR_Model(Model,Q,test_var,nderiv,option=0)

 write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  Q = [2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,2.23738_Rkind,1.975200_Rkind,ZERO]

  CALL Test_QVG_FOR_Model(Model,Q,test_var,nderiv,option=2)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of q(6), the torsion)'
  write(out_unit,*) '   file name: "grid_HONO"'

  CALL Eval_pot_ON_Grid(Model,                                      &
                        Qmin=[2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,&
                              2.23738_Rkind,1.975200_Rkind,ZERO],        &
                        Qmax=[2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,&
                              2.23738_Rkind,1.975200_Rkind,pi],          &
                        nb_points=1001, grid_file='grid_HONO')

  CALL dealloc_Model(Model)
 
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ Cartesian coordinates ----------'
  write(out_unit,*) '---------------------------------------------'
      !1          8           0        0.000000    0.000000    0.000000
      !2          7           0        0.000000    0.000000    1.427049
      !3          1           0        0.944081    0.000000   -0.198113
      !4          8           0       -1.095870    0.000000    1.840421

  CALL Init_Model(Model,pot_name='HONO',Cart_TO_Q=.TRUE.)

  allocate(xyz(Model%ndim))
  xyz = [ZERO,           ZERO,  ZERO,           &
         ZERO,           ZERO,  1.427049_Rkind, &
         0.944081_Rkind, ZERO, -0.198113_Rkind, &
        -1.095870_Rkind, 0.1_Rkind,  1.840421_Rkind]/a0


  CALL Eval_Pot(Model,xyz,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  deallocate(xyz)

  CALL dealloc_dnMat(PotVal)
  CALL dealloc_dnMat(PotValref)
  CALL dealloc_dnMat(dnErr)

  deallocate(Q)

  CALL dealloc_Model(Model)


  CALL Finalize_Test(test_var)

END PROGRAM TEST_HONO