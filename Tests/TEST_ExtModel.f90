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
PROGRAM TEST_model
  USE QDUtil_Test_m
  USE QDUtil_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE


  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: nderiv
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotValref
  TYPE (dnMat_t)                 :: dnErr

  TYPE (test_t)                  :: test_var
  logical                        :: Lerr
  integer                        :: err

  CALL Initialize_Test(test_var,test_name='QModel_ExtModel')

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' ExtModel potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)

  CALL Init_Model(QModel,pot_name='ExtModel',adiabatic=.FALSE.)
  nderiv=2

  allocate(q(QModel%QM%ndim))
  CALL QModel%QM%RefValues_QModel(err,Q0=Q,nderiv=nderiv)

  CALL Logical_Test(test_var,test1=(err == 0),info='Q0 initialization')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  flush(out_unit)
  nderiv = 3
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL QModel%QM%RefValues_QModel(err,dnMatV=PotValref,nderiv=nderiv)
  dnErr = PotValref-PotVal
  Lerr  = Check_dnMat_IS_ZERO(dnErr)
  
  CALL Logical_Test(test_var,test1=Lerr,info='dnVMat')


  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'



  CALL Finalize_Test(test_var)

END PROGRAM TEST_model