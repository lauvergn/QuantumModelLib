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
PROGRAM TEST_CHFClBr
  USE QDUtil_NumParameters_m
  USE QDUtil_Test_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: Model
  real (kind=Rkind), allocatable :: Q(:),Qref(:),xyz(:)
  integer                        :: ndim,nsurf,nderiv,i,option,err
  logical                        :: Lerr
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotValref
  TYPE (dnMat_t)                 :: dnErr
  TYPE (test_t)                  :: test_var
  real (kind=Rkind), parameter   :: epsi = 1.e-10_Rkind
  real (kind=Rkind), parameter   :: a0 = 0.52917720835354106_Rkind ! from Tnum


  CALL Initialize_Test(test_var,test_name='QModel_CHFClBr')

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' CHFClBr potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(Model,pot_name='CHFClBr',read_param=.FALSE.)
  Q = [0.1_Rkind,-0.1_Rkind,0.2_Rkind,-0.2_Rkind,0.3_Rkind,-0.3_Rkind,0.4_Rkind,-0.4_Rkind,0.5_Rkind]

  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,'(a,9f12.6)') 'Q (Bohr)',Q(:)
  CALL Check_analytical_numerical_derivatives(Model,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential'
  Q = [0.1_Rkind,-0.1_Rkind, 0.2_Rkind,-0.2_Rkind, 0.3_Rkind,-0.3_Rkind,0.4_Rkind,-0.4_Rkind,0.5_Rkind]

  CALL Test_QVG_FOR_Model(Model,Q,test_var,option=0,nderiv=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------' 
 

  deallocate(Q)
  CALL dealloc_Model(Model)


  CALL Finalize_Test(test_var)

END PROGRAM TEST_CHFClBr