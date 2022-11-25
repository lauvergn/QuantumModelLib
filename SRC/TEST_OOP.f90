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
program prog
  USE QML_Empty_m
  USE QML_Template_m
  USE QML_Morse_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE


  real(kind=Rkind),  allocatable    :: Q(:)
  TYPE (dnS_t),      allocatable    :: dnQ(:)
  TYPE (dnS_t),      allocatable    :: Mat_OF_PotDia(:,:)
    TYPE (dnMat_t)                  :: PotVal

  TYPE(Model_t)                     :: QModel

  integer :: i,nderiv

  nderiv     = 2

  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  CALL Init_Model(QModel,pot_name='read_model')
  !CALL Write_Model(QModel,nio=out_unitp)

  write(out_unitp,*) 'Gdef',QModel%QM%get_d0GGdef_QModel()

  allocate(Q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)
  IF (allocated(Q)) write(out_unitp,*) 'Q0',Q
  flush(out_unitp)

  Q(:) = [ 1.2_Rkind,0.2_Rkind ]
  CALL Eval_Pot(QModel,Q,PotVal,nderiv)
  write(out_unitp,*) 'PotVal'
  CALL Write_dnMat(PotVal,6)

  CALL dealloc_dnMat(PotVal)

  write(out_unitp,*) "---------------------------------------------------"
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  deallocate(QModel%QM)
  deallocate(Q)
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  write(out_unitp,*)
  flush(out_unitp)
stop
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  CALL Init_Model(QModel,pot_name='morse')
  !CALL Write_Model(QModel,nio=out_unitp)

  write(out_unitp,*) 'Gdef',QModel%QM%get_d0GGdef_QModel()

  allocate(Q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)
  IF (allocated(Q)) write(out_unitp,*) 'Q0',Q

  CALL Eval_Pot(QModel,Q,PotVal,nderiv)
  write(out_unitp,*) 'PotVal'
  CALL Write_dnMat(PotVal,6)

  CALL dealloc_dnMat(PotVal)
  deallocate(QModel%QM)
  deallocate(Q)
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  write(out_unitp,*)
  flush(out_unitp)

  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  CALL Init_Model(QModel,pot_name='template')
  !CALL Write_Model(QModel,nio=out_unitp)

  write(out_unitp,*) 'Gdef',QModel%QM%get_d0GGdef_QModel()

  allocate(Q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)
  IF (allocated(Q)) write(out_unitp,*) 'Q0',Q

  CALL Eval_Pot(QModel,Q,PotVal,nderiv)
  write(out_unitp,*) 'PotVal'
  CALL Write_dnMat(PotVal,6)

  CALL dealloc_dnMat(PotVal)

  write(out_unitp,*) "---------------------------------------------------"
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  deallocate(QModel%QM)
  deallocate(Q)
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  write(out_unitp,*)
  flush(out_unitp)


end program prog
