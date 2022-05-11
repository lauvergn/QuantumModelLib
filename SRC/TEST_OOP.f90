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

  CALL QML_dealloc_dnMat(PotVal)

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

  CALL QML_dealloc_dnMat(PotVal)
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

  CALL QML_dealloc_dnMat(PotVal)

  write(out_unitp,*) "---------------------------------------------------"
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  deallocate(QModel%QM)
  deallocate(Q)
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  write(out_unitp,*)
  flush(out_unitp)


end program prog
