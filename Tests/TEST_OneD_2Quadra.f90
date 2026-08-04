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
PROGRAM TEST_OneD_Photons
  USE QDUtil_NumParameters_m
  USE QDUtil_m
  USE QDUtil_Test_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:),Qref(:),xyz(:)
  real (kind=Rkind), allocatable :: G(:,:),Gref(:,:) ! metric tensor
  integer                        :: ndim,nsurf,nderiv,i,option,err
  logical                        :: Lerr
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: ScalOp(2)
  TYPE (dnMat_t)                 :: PotValref
  TYPE (dnMat_t)                 :: dnErr
  TYPE (test_t)                  :: test_var
  real (kind=Rkind), parameter   :: epsi = 1.e-10_Rkind


  CALL Initialize_Test(test_var,test_name='QModel_OneD_2Quadra')


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' 1D-diabatic model + dipole moment transition '
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='OneD_2Quadra',adiabatic=.FALSE.,Print_init=.TRUE.,nb_ScalOp=2)

  Q = [THREE]

  nderiv=2

 write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Diabatic Potential and derivatives',nderiv

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  write(out_unit,*) 'Energy (Hartree)',get_Flatten(PotVal)
  CALL Write_dnMat(PotVal,nio=out_unit)


  ! For testing the model
  allocate(Qref(QModel%ndim))
  allocate(Gref(QModel%ndim,QModel%ndim))
  CALL QModel%QM%RefValues_QModel(err,nderiv=nderiv,Q0=Qref,d0GGdef=Gref,dnMatV=PotValref,option=1)

  write(out_unit,*) 'Reference Energy (Hartree)' ; flush(out_unit)
  CALL Write_dnMat(PotValref,nio=out_unit)
  flush(out_unit)
  dnErr = PotValref-PotVal
  Lerr  = Check_dnMat_IS_ZERO(dnErr)

  CALL Logical_Test(test_var,test1=Lerr,info='dnMatV')
  IF (.NOT. Lerr) CALL Write_dnMat(dnErr,nio=out_unit,info='dnErr')

  Lerr = all(abs(Q-Qref) < epsi)
  CALL Logical_Test(test_var,test1=Lerr,info='Q(:)')
  IF (.NOT. Lerr) Write(out_unit,*) 'Q-Qref',Q-Qref

  G = get_d0GGdef_Model(Model=QModel)
  Lerr = all(abs(G-Gref) < epsi)
  CALL Logical_Test(test_var,test1=Lerr,info='G (metrix tensor)')
  IF (.NOT. Lerr) Write(out_unit,*) 'G-Gref',G-Gref

  deallocate(Qref)
  deallocate(Gref)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- Scalar Operators --------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Eval_ScalOp(QModel,Q,ScalOp,nderiv=nderiv)

  CALL Write_dnMat(ScalOp(1),nio=out_unit,info='pot (ScalOp)')
  dnErr = PotValref-PotVal
  Lerr  = Check_dnMat_IS_ZERO(dnErr)
  CALL Logical_Test(test_var,test1=Lerr,info='dnMatV (ScalOp)')
  IF (.NOT. Lerr) CALL Write_dnMat(dnErr,nio=out_unit,info='dnErr')

  CALL Write_dnMat(ScalOp(2),nio=out_unit,info='mu')
  flush(out_unit)

  CALL dealloc_dnMat(PotValref)
  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_dnMat(dnErr)
  deallocate(G)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' 1D-adiabatic model + dipole moment transition '
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='OneD_2Quadra',adiabatic=.TRUE.,Print_init=.TRUE.,nb_ScalOp=2)

  Q = [THREE]

  nderiv=2

 write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Adiabatic Potential and derivatives',nderiv

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)


  ! For testing the model
  allocate(Qref(QModel%ndim))
  allocate(Gref(QModel%ndim,QModel%ndim))
  CALL QModel%QM%RefValues_QModel(err,nderiv=nderiv,Q0=Qref,d0GGdef=Gref,dnMatV=PotValref,option=2)

  write(out_unit,*) 'Reference Energy (Hartree)' ; flush(out_unit)
  CALL Write_dnMat(PotValref,nio=out_unit)
  flush(out_unit)
  dnErr = PotValref-PotVal
  Lerr  = Check_dnMat_IS_ZERO(dnErr)

  CALL Logical_Test(test_var,test1=Lerr,info='dnMatV')
  IF (.NOT. Lerr) CALL Write_dnMat(dnErr,nio=out_unit,info='dnErr')

  Lerr = all(abs(Q-Qref) < epsi)
  CALL Logical_Test(test_var,test1=Lerr,info='Q(:)')
  IF (.NOT. Lerr) Write(out_unit,*) 'Q-Qref',Q-Qref

  G = get_d0GGdef_Model(Model=QModel)
  Lerr = all(abs(G-Gref) < epsi)
  CALL Logical_Test(test_var,test1=Lerr,info='G (metrix tensor)')
  IF (.NOT. Lerr) Write(out_unit,*) 'G-Gref',G-Gref

  deallocate(Qref)
  deallocate(Gref)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- Scalar Operators --------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Eval_ScalOp(QModel,Q,ScalOp,nderiv=nderiv)
 
  CALL Write_dnMat(ScalOp(1),nio=out_unit,info='pot (ScalOp)')
  dnErr = PotValref-PotVal
  Lerr  = Check_dnMat_IS_ZERO(dnErr)
  CALL Logical_Test(test_var,test1=Lerr,info='dnMatV (ScalOp)')
  IF (.NOT. Lerr) CALL Write_dnMat(dnErr,nio=out_unit,info='dnErr')

  CALL Write_dnMat(ScalOp(2),nio=out_unit,info='mu')
  flush(out_unit)

  CALL dealloc_dnMat(PotValref)
  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_dnMat(dnErr)
  deallocate(G)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_Model(QModel)


  CALL Finalize_Test(test_var)

END PROGRAM TEST_OneD_Photons