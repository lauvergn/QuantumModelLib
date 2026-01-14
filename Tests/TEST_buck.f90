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
PROGRAM TEST_buck
  USE QDUtil_NumParameters_m
  USE QDUtil_Test_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: Model
  real (kind=Rkind), allocatable :: Q(:),Qref(:)
  real (kind=Rkind), allocatable :: G(:,:),Gref(:,:) ! metric tensor
  integer                        :: ndim,nsurf,nderiv,i,option,err
  logical                        :: Lerr
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotValref
  TYPE (dnMat_t)                 :: dnErr
  TYPE (test_t)                  :: test_var
  real (kind=Rkind), parameter   :: epsi = 1.e-10_Rkind


  CALL Initialize_Test(test_var,test_name='QModel_Buck')

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Buckingham potential (Ar-Ar parameters)'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(Model,pot_name='Buck',read_param=.FALSE.)

  Q = [7._Rkind]

  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,'(a,f12.6)') 'R (Bohr)',Q(:)
  CALL Check_analytical_numerical_derivatives(Model,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'

  CALL Eval_Pot(Model,Q,PotVal,nderiv=nderiv)
  write(out_unit,'(a,f12.6)') 'R (Bohr)',q(:)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)
  flush(out_unit)

  ! For testing the model
  allocate(Qref(Model%ndim))
  allocate(Gref(Model%ndim,Model%ndim))
  CALL Model%QM%RefValues_QModel(err,Q0=Qref,d0GGdef=Gref,dnMatV=PotValref,nderiv=nderiv)
  write(out_unit,*) 'Reference Energy (Hartree)'
  CALL Write_dnMat(PotValref,nio=out_unit)
  flush(out_unit)
  dnErr = PotValref-PotVal
  Lerr  = Check_dnMat_IS_ZERO(dnErr)

  CALL Logical_Test(test_var,test1=Lerr,info='dnMatV')

  Lerr = all(abs(Q-Qref) < epsi)
  CALL Logical_Test(test_var,test1=Lerr,info='Q(:)')

  G = get_d0GGdef_Model(Model=Model)
  Lerr = all(abs(G-Gref) < epsi)
  CALL Logical_Test(test_var,test1=Lerr,info='G (metrix tensor)')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)
  CALL dealloc_dnMat(PotValref)
  CALL dealloc_dnMat(dnErr)

  deallocate(Q)
  deallocate(Qref)
  deallocate(G)
  deallocate(Gref)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of R)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '   file name: "grid_Buck"'

  CALL Eval_pot_ON_Grid(Model,Qmin=[6._Rkind],Qmax=[20._Rkind],nb_points=1001,grid_file='grid_Buck')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_Model(Model)

  CALL Finalize_Test(test_var)

END PROGRAM TEST_buck