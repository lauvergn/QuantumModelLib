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
PROGRAM TEST_model
  IMPLICIT NONE

  !CALL test_Retinal_JPCB2000() ; stop
  !CALL test_Tully_test() ; stop
  !CALL test_PSB3_test() ; stop
  !CALL test_H3 ; stop
  !CALL test_PSB3_Retinal2000_test ; stop
  !CALL test_Test() ; stop
  !CALL test_CH5() ; stop


  ! One electronic surface
  CALL test_LinearHBond()
  CALL test_HenonHeiles()
  CALL test_Buckingham()
  CALL test_Morse()

  ! Several electronic surfaces
  CALL test_Tully()
  CALL test_OneDSOC_1S1T()
  CALL test_OneDSOC_2S1T()

  CALL test_Phenol()
  CALL test_PSB3()
  CALL test_TwoD()
  CALL test_Retinal_JPCB2000()

  ! 6D (full-D), One electronic surface
  CALL test_HONO()
  CALL test_HNNHp()

  CALL test_H2SiN()
  CALL test_H2NSi()

  ! 3D (full-D), One electronic surface for collision
  CALL test_HOO_DMBE()

  ! 12D (full-D), One electronic surface for collision H+CH4 -> H2+CH3
  CALL test_CH5()

  ! A template with one electronic surface
  CALL test_template()

END PROGRAM TEST_model

SUBROUTINE test_Tully
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Tully potential'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='Tully',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  DO option=1,3
    CALL Init_Model(QModel,pot_name='Tully',option=option)
    allocate(q(QModel%QM%ndim))
    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,*) '----- CHECK POT -----------------------------'
    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,'(a,i1)') ' ----- TULLY',option

    write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

    q(:) = 1.1_Rkind
    write(out_unitp,'(a,f12.6)') ' R (Bohr)',q(:)
    CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

    q(:) = -1.1_Rkind
    write(out_unitp,'(a,f12.6)') ' R (Bohr)',q(:)
    CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)


    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,*) ' Potential and derivatives'

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
    write(out_unitp,'(a,f12.6)') 'R (Bohr)',q(:)
    write(out_unitp,*) 'Energy (Hartree)'
    CALL QML_Write_dnMat(PotVal,nio=out_unitp)
    ! For testing the model
    CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='TULLY ' // int_TO_char(option))

    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,*) '- END CHECK POT -----------------------------'
    write(out_unitp,*) '---------------------------------------------'

    CALL QML_dealloc_dnMat(PotVal)
    deallocate(q)

  END DO



  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '-----------TULLY 1 on a grid ----------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Tully1"'
  write(out_unitp,*) ' You should get the curves (columns: 2,3 5) of the panel (a) of figure 3 of ... '
  write(out_unitp,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unitp,*) ' Be carrefull, the non-adiabatic coupling is scaled by -1/50.'
  CALL Init_Model(QModel,pot_name='Tully',option=1)
  CALL Eval_pot_ON_Grid(QModel,Qmin=(/-TEN/),Qmax=(/TEN/),nb_points=1001, &
                        grid_file='grid_Tully1')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '-----------TULLY 2 on a grid ----------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Tully2"'
  write(out_unitp,*) ' You should get the curves (columns: 2,3 5) of the panel (b) of figure 3 of ... '
  write(out_unitp,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unitp,*) ' Be carrefull, the non-adiabatic coupling is scaled by -1/12.'

  CALL Init_Model(QModel,pot_name='Tully',option=2)
  CALL Eval_pot_ON_Grid(QModel,Qmin=(/-TEN/),Qmax=(/TEN/),nb_points=1001, &
                        grid_file='grid_Tully2')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '-----------TULLY 3 on a grid ----------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Tully3"'
  write(out_unitp,*) ' You should get the curves (columns: 2,3 5) of the panel (c) of figure 3 of ... '
  write(out_unitp,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unitp,*) ' Be carrefull, the sign of non-adiabatic coupling is changed.'

  CALL Init_Model(QModel,pot_name='Tully',option=3)
  CALL Eval_pot_ON_Grid(QModel,Qmin=(/-TEN/),Qmax=(/TEN/),nb_points=1001, &
                        grid_file='grid_Tully3')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Tully

SUBROUTINE test_OneDSOC_1S1T
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' 1D-SOC_1S1T potential'
  write(out_unitp,*) ' With units: Bohr Hartree (au)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='1DSOC_1S1T',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,nsurf=4,pot_name='1DSOC_1S1T')

  allocate(q(QModel%QM%ndim))
  q(:) = (/ 9.5_Rkind /) ! one evaluation to get vec0%d0(:,:)
  write(out_unitp,*) ' ref Eigenvectors at Q:',Q
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=0)
  ! ref Eigenvectors
  CALL QML_Write_dnMat(QModel%QM%vec0,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,Vec=QModel%QM%vec0,info='1DSOC_1S1T')



  q(:) = (/ 8.5_Rkind /)
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Q (Bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  write(out_unitp,*) 'Non adiatic couplings:'
  CALL QML_Write_dnMat(NAC,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='1DSOC_1S1T')

  write(out_unitp,*) 'DIABATIC potential'
  QModel%QM%adiabatic = .FALSE.
  write(out_unitp,*) 'Q (Bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)
  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='1DSOC_1S1T')

  CALL QML_dealloc_dnMat(PotVal)
  CALL QML_dealloc_dnMat(NAC)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  flush(out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----------- 1D-SOC_1S1T on a grid -----------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of R in Bohr)'
  write(out_unitp,*) '   file name: "grid_1DSOC_1S1T"'
  write(out_unitp,*) ' You should get the 3 curves (left panels) of figure 1 of ... '
  write(out_unitp,*) '  ... Granucci et al. J. Chem. Phys. V137, p22A501 (2012)'

  CALL Init_Model(QModel,nsurf=4,pot_name='1DSOC_1S1T',option=1)
  QModel%QM%adiabatic = .TRUE.
  flush(out_unitp)

  q(:) = (/ 9.5_Rkind /) ! one evaluation to get vec0%d0(:,:)
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=1)

  ! ref Eigenvectors
  write(out_unitp,*) ' ref Eigenvectors at Q:',Q
  CALL QML_Write_dnMat(QModel%QM%vec0,nio=out_unitp)


  CALL Eval_pot_ON_Grid(QModel,Qmin=(/3._Rkind/),Qmax=(/20._Rkind/), &
                        nb_points=1001,nderiv=0,grid_file='grid_1DSOC_1S1T')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  deallocate(q)
  CALL QML_dealloc_dnMat(PotVal)
  CALL QML_dealloc_dnMat(NAC)
END SUBROUTINE test_OneDSOC_1S1T
SUBROUTINE test_OneDSOC_2S1T
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' 1D-SOC_2S1T potential'
  write(out_unitp,*) ' With units: Bohr Hartree (au)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='1DSOC_2S1T')
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='1DSOC_2S1T',Print_init=.FALSE.)

  allocate(q(QModel%QM%ndim))
  q(:) = (/ 9.5_Rkind /) ! one evaluation to get vec0%d0(:,:)
  write(out_unitp,*) ' ref Eigenvectors at Q:',Q
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=0)
  ! ref Eigenvectors
  CALL QML_Write_dnMat(QModel%QM%vec0,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,Vec=QModel%QM%vec0,info='1DSOC_2S1T')



  q(:) = (/ 8.5_Rkind /)
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Q (Bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  write(out_unitp,*) 'Non adiatic couplings:'
  CALL QML_Write_dnMat(NAC,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='1DSOC_2S1T')

  write(out_unitp,*) 'DIABATIC potential'
  QModel%QM%adiabatic = .FALSE.
  write(out_unitp,*) 'Q (Bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)
  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='1DSOC_2S1T')

  CALL QML_dealloc_dnMat(PotVal)
  CALL QML_dealloc_dnMat(NAC)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  flush(out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----------- 1D-SOC_2S1T on a grid -----------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of R in Bohr)'
  write(out_unitp,*) '   file name: "grid_1DSOC_2S1T"'

  CALL Init_Model(QModel,pot_name='1DSOC_2S1T')
  QModel%QM%adiabatic = .TRUE.
  flush(out_unitp)

  q(:) = (/ 9.5_Rkind /) ! one evaluation to get vec0%d0(:,:)
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=1)

  ! ref Eigenvectors
  write(out_unitp,*) ' ref Eigenvectors at Q:',Q
  CALL QML_Write_dnMat(QModel%QM%vec0,nio=out_unitp)


  CALL Eval_pot_ON_Grid(QModel,Qmin=(/3._Rkind/),Qmax=(/20._Rkind/), &
                        nb_points=1001,nderiv=0,grid_file='grid_1DSOC_2S1T')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  deallocate(q)
  CALL QML_dealloc_dnMat(PotVal)
  CALL QML_dealloc_dnMat(NAC)
END SUBROUTINE test_OneDSOC_2S1T
SUBROUTINE test_Morse
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Morse potential (H-F parameters)'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='Morse',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Morse',read_param=.FALSE.)

  allocate(q(QModel%QM%ndim))
  q(:) = 2._Rkind

  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,'(a,f12.6)') 'R (Bohr)',q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unitp,'(a,f12.6)') 'R (Bohr)',q(:)
  write(out_unitp,*) 'Energy (Hartree)'
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='Morse_HF')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of R)'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Morse"'

  CALL Eval_pot_ON_Grid(QModel,Qmin=(/1._Rkind/),Qmax=(/5._Rkind/),nb_points=1001,&
                        grid_file='grid_Morse')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Morse

SUBROUTINE test_Buckingham
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Buckingham potential (Ar-Ar parameters)'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='Buck',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Buck',read_param=.FALSE.)

  allocate(q(QModel%QM%ndim))
  q(:) = 7._Rkind
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,'(a,f12.6)') 'R (Bohr)',q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unitp,'(a,f12.6)') 'R (Bohr)',q(:)
  write(out_unitp,*) 'Energy (Hartree)'
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='Buckingham_Ar-Ar')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of R)'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Buck"'

  CALL Eval_pot_ON_Grid(QModel,Qmin=(/6._Rkind/),Qmax=(/20._Rkind/),nb_points=1001,&
                        grid_file='grid_Buck')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Buckingham

SUBROUTINE test_Phenol
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Phenol potential'
  write(out_unitp,*) ' With units: Angs and eV'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='Phenol')
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Phenol',PubliUnit=.TRUE.)

  allocate(q(QModel%QM%ndim))

  q(:) = (/ 1.2_Rkind,0.2_Rkind /)
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q (Angs, radian):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Q (Angs, radian):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  write(out_unitp,*) 'Non adiatic couplings:'
  CALL QML_Write_dnMat(NAC,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='Phenol')

  write(out_unitp,*) 'DIABATIC potential'
  QModel%QM%adiabatic = .FALSE.
  write(out_unitp,*) 'Q (Angs, radian):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='Phenol')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with PubliUnit=.FALSE.'
  CALL Init_Model(QModel,pot_name='Phenol',PubliUnit=.FALSE.,Print_init=.FALSE.)
  CALL get_Q0_Model(q,QModel,option=0)
  write(out_unitp,'(a,2f12.6)') 'q (au)',q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL QML_dealloc_dnMat(PotVal)
  CALL QML_dealloc_dnMat(NAC)
  deallocate(q)


  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----------- Phenol on a grid ----------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of R in Angstrom)'
  write(out_unitp,*) '   file name: "grid_Phenol"'
  write(out_unitp,*) ' You should get the 3 curves (columns: 2,6, 10) of figure 2 of ... '
  write(out_unitp,*) '   Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, ...'
  write(out_unitp,*) '  .... J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218)'

  CALL Init_Model(QModel,pot_name='phenol',PubliUnit=.TRUE.)
  QModel%QM%adiabatic = .FALSE.

  CALL Eval_pot_ON_Grid(QModel,Qmin=(/0.5_Rkind,ZERO/),Qmax=(/5._Rkind,ZERO/), &
                        nb_points=1001,grid_file='grid_Phenol')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Phenol
SUBROUTINE test_HenonHeiles
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 4D-HenonHeiles -----------------'
  CALL Init_Model(QModel,pot_name='HenonHeiles',ndim=4,Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HenonHeiles',ndim=4)

  allocate(q(QModel%QM%ndim))
  q(:) = (/ (0.1_Rkind*real(i,kind=Rkind),i=1,QModel%QM%ndim) /)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='4D-HenonHeiles')

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_HenonHeiles
SUBROUTINE test_LinearHBond
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Linear H-Bond (symmetric one and the default parameters)'
  write(out_unitp,*) ' With units: Angs and kcal.mol^-1'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='HBond',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HBond',PubliUnit=.TRUE.)

  allocate(q(QModel%QM%ndim))

  q(:) = (/ 2.75_Rkind,0._Rkind /)
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,'(a,2f12.6)') 'QQ,q (Angs)',q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unitp,'(a,2f12.6)') 'QQ,q (Angs)',q(:)
  write(out_unitp,*) 'Energy (kcal.mol^-1)'
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='Linear H-Bond (symmetric one)')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with PubliUnit=.FALSE.'
  CALL Init_Model(QModel,pot_name='HBond',PubliUnit=.FALSE.,Print_init=.FALSE.)
  CALL get_Q0_Model(q,QModel,option=0)
  write(out_unitp,'(a,2f12.6)') 'QQ,q (Bohr)',q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'


  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  CALL QML_dealloc_dnMat(PotVal)

  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unitp,*) '   file name: "grid_Hbond-sym"'
  write(out_unitp,*) ' You should get the green curve (QQ=2.75 Angs) of ... '
  write(out_unitp,*) '...the bottom left panel of figure 3.8 (Julien Beutier thesis)'

  CALL Eval_pot_ON_Grid(QModel, &
                        Qmin=(/2.75_Rkind,-0.6_Rkind/),Qmax=(/2.75_Rkind,0.6_Rkind/),nb_points=1001,&
                        grid_file='grid_Hbond-sym')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Linear H-Bond (asymmetric one)'
  write(out_unitp,*) ' With units: Angs and kcal.mol^-1'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HBond',PubliUnit=.TRUE.,         &
                  read_param=.TRUE.,param_file_name='Hbond.dat')
  CALL Write_Model(QModel)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  allocate(q(QModel%QM%ndim))

  q(:) = (/ 3._Rkind,0._Rkind /)
  nderiv=2
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unitp,'(a,2f12.6)') 'QQ,q (Angs)',q(:)
  write(out_unitp,*) 'Energy (kcal.mol^-1)'
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='Linear H-Bond (asymmetric one)')

  CALL QML_dealloc_dnMat(PotVal)

  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unitp,*) '   file name: "grid_Hbond-asym"'
  write(out_unitp,*) ' You should get the green curve (QQ=3. Angs) of ... '
  write(out_unitp,*) '...the bottom right panel of figure 3.8 (Julien Beutier thesis)'

  CALL Eval_pot_ON_Grid(QModel, &
                        Qmin=(/3._Rkind,-0.8_Rkind/),Qmax=(/3._Rkind,0.7_Rkind/),nb_points=1001,&
                        grid_file='grid_Hbond-asym')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_LinearHBond
SUBROUTINE test_PSB3

  USE mod_Lib
  USE mod_dnMat
  USE mod_Model

  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' PSB3 potential'
  write(out_unitp,*) ' With units: Atomic Units (Angstrom, Rad, Rad, kcal.mol^-1)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='PSB3',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='PSB3',PubliUnit=.FALSE.)

  allocate(q(QModel%QM%ndim))

  q(:) = (/0.172459_Rkind,-3.14_Rkind,0._Rkind/)
  nderiv=3

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Evaluated in', q
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  write(out_unitp,*) '---------------------------------------------'

  QModel%QM%adiabatic = .FALSE.
  write(out_unitp,*) 'DIABATIC potential'
  write(out_unitp,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='PSB3')

  QModel%QM%adiabatic = .TRUE.
  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  write(out_unitp,*) 'Non adiatic couplings:'
  CALL QML_Write_dnMat(NAC,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='PSB3')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with PubliUnit=.FALSE.'
  DO option=1,2
    CALL Init_Model(QModel,pot_name='PSB3',PubliUnit=.FALSE.,Print_init=.FALSE.,option=option)
    CALL get_Q0_Model(q,QModel,option=0)
    write(out_unitp,*) 'q (au)',q(:)
    CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)
    CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
    CALL QML_Write_dnMat(PotVal,nio=out_unitp)
  END DO
  write(out_unitp,*) '---------- END CHECK POT --------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL QML_dealloc_dnMat(PotVal)

  deallocate(q)

END SUBROUTINE test_PSB3
SUBROUTINE test_Retinal_JPCB2000

  USE mod_Lib
  USE mod_dnMat
  USE mod_Model

  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings
  real (kind=Rkind)              :: DQ2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Retinal_JPCB2000 potential'
  write(out_unitp,*) ' With units: Atomic Units'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='Retinal_JPCB2000',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Retinal_JPCB2000',PubliUnit=.FALSE.)

  allocate(q(QModel%QM%ndim))

  q(:) = [ZERO,ONE]
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Evaluated in', q
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  write(out_unitp,*) '---------------------------------------------'

  QModel%QM%adiabatic = .FALSE.
  write(out_unitp,*) 'DIABATIC potential'
  write(out_unitp,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='Retinal_JPCB2000')

  QModel%QM%adiabatic = .TRUE.
  q(:) = [ZERO,ZERO]
  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  write(out_unitp,*) 'Non adiatic couplings:'
  CALL QML_Write_dnMat(NAC,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='Retinal_JPCB2000')

  write(out_unitp,*) '---------- END CHECK POT --------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '-----------Retinal_JPCB2000 on a 2D-grid ----'
  write(out_unitp,*) '---------------------------------------------'

  QModel%QM%adiabatic = .TRUE.
  DQ2 = TWO
  CALL Eval_pot_ON_Grid(QModel,Qmin=[-pi/TWO,-DQ2],Qmax=[THREE*pi/TWO,DQ2],         &
                        nb_points=101,grid_file='grid_Retinal_JPCB2000')


  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Retinal_JPCB2000 potential+Bath'
  write(out_unitp,*) ' With units: Atomic Units'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,ndim=4,pot_name='Retinal_JPCB2000',PubliUnit=.FALSE.)

  allocate(q(QModel%QM%ndim))

  q(:) = ZERO
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Evaluated in',q

  QModel%QM%adiabatic = .TRUE.
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  write(out_unitp,*) 'Non adiatic couplings:'
  CALL QML_Write_dnMat(NAC,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='Retinal_JPCB2000')

  write(out_unitp,*) '---------- END CHECK POT --------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

END SUBROUTINE test_Retinal_JPCB2000
SUBROUTINE test_HONO
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 6D-HONO ------------------------'
  CALL Init_Model(QModel,pot_name='HONO',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HONO')

  allocate(q(QModel%QM%ndim))
  q(:) = [2.696732586_Rkind,1.822912197_Rkind,1.777642018_Rkind,2.213326419_Rkind,1.9315017_Rkind,pi]
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives at the trans minimum'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='HONO')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives at the cis minimum'
  q(:) = [2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,2.23738_Rkind,1.975200_Rkind,ZERO]

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='HONO')

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)



  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of q(6), the torsion)'
  write(out_unitp,*) '   file name: "grid_HONO"'

  CALL Eval_pot_ON_Grid(QModel,                                      &
                        Qmin=[2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,&
                              2.23738_Rkind,1.975200_Rkind,ZERO],        &
                        Qmax=[2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,&
                              2.23738_Rkind,1.975200_Rkind,pi],          &
                        nb_points=1001, grid_file='grid_HONO')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_HONO
SUBROUTINE test_HNNHp
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 6D-HNNH+ -----------------------'
  CALL Init_Model(QModel,pot_name='HNNHp',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HNNHp')

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives at the trans minimum'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='HNNHp')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_HNNHp
SUBROUTINE test_H2SiN
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 6D-H2SiN -----------------------'
  CALL Init_Model(QModel,pot_name='H2SiN',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
DO option=1,3

  CALL Init_Model(QModel,pot_name='H2SiN',option=option)

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives at:'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='H2SiN')

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
END DO

END SUBROUTINE test_H2SiN
SUBROUTINE test_H2NSi
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 6D-H2NSi -----------------------'
  CALL Init_Model(QModel,pot_name='H2NSi',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
DO option=1,2

  CALL Init_Model(QModel,pot_name='H2NSi',option=option)

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives at:'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='H2NSi')

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
END DO

END SUBROUTINE test_H2NSi
SUBROUTINE test_template
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Template potential'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='template',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  flush(out_unitp)

  CALL Init_Model(QModel,pot_name='template')

  allocate(q(QModel%QM%ndim))
  Q(:) = 2._Rkind

  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q(:) (bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unitp,*) 'Q(:) (bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  write(out_unitp,*) 'Energy (Hartree)'
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='Template')

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'



END SUBROUTINE test_template
SUBROUTINE test_TwoD
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' TwoD potential'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='TwoD',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  flush(out_unitp)

  CALL Init_Model(QModel,pot_name='TwoD',adiabatic=.FALSE.)
  !CALL Init_Model(QModel,pot_name='TwoD',adiabatic=.TRUE.)

  allocate(q(QModel%ndim))
  Q(:) = [3.875_Rkind,0.5_Rkind]

  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q(:) (bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives',nderiv

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unitp,*) 'Q(:) (bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%ndim)
  write(out_unitp,*) 'Energy (Hartree)'
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='twod')

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'



END SUBROUTINE test_TwoD
SUBROUTINE test_HNO3
  USE mod_Lib
  USE mod_dnS
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnS_t),      allocatable :: Func(:)

  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 9D-HNO3 ------------------------'
  CALL Init_Model(QModel,pot_name='HNO3',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,ndim=9,pot_name='HNO3')

  allocate(q(QModel%ndim))
  Q(:) = [2.65450_Rkind,2.28_Rkind,2.0_Rkind,Pi/TWO,                    &
          ZERO,ZERO,1.83492_Rkind,1.77157_Rkind,Pi/TWO]
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='HNO3')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Functions'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%ndim)

  CALL Eval_Func(QModel,Q,Func,nderiv=nderiv)

  DO i=1,size(Func)
    CALL QML_Write_dnS(Func(i))
  END DO
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_HNO3
SUBROUTINE test_CH5
  USE mod_Lib
  USE mod_dnS
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,i1,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnS_t),      allocatable :: Func(:)

  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 12D-CH5: H+CH4 -> H2+CH3 -------'
  write(out_unitp,*) ' Qmodel, init: ',check_Init_QModel(QModel)
  CALL Init_Model(QModel,pot_name='CH5',Print_init=.FALSE.)
  write(out_unitp,*) ' Qmodel, init: ',check_Init_QModel(QModel)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,ndim=12,pot_name='CH5',option=5)

  allocate(Q(QModel%ndim))
  CALL get_Q0_Model(Q,QModel,option=5)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  ! For testing the model
  CALL Write_QdnV_FOR_Model(Q,PotVal,QModel,info='CH5')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Functions'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%ndim)

  DO i=0,200
    Q(1)= -FIVE + real(i,kind=Rkind)*TEN/real(200,kind=Rkind)
    CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
    write(6,*) Q,QML_get_d0_FROM_dnS(Func)
  END DO

  Q(1)= 0.5_Rkind
  CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
  write(6,*) Q(1),'Energy',QML_get_d0_FROM_dnS(Func(1))
  write(6,*) Q(1),'Qopt',QML_get_d0_FROM_dnS(Func(2:12))
  write(6,*) Q(1),'hessian'
  DO i=1,121
    write(6,*) '                       ',QML_get_d0_FROM_dnS(Func(12+i)),',     &'
  END DO

  deallocate(Q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential in the asymptotic region'

  allocate(Q(QModel%ndim))
  Q(1)= 5._Rkind
  CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
  Q(2:12) = QML_get_d0_FROM_dnS(Func(2:12))

  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,QModel%ndim)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_CH5
SUBROUTINE test_HOO_DMBE
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:,:),x(:,:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 3D-HOO_DMBE --------------------'
  CALL Init_Model(QModel,pot_name='HOO_DMBE',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HOO_DMBE')

  allocate(Q(QModel%ndim,3))

  ! OH...O in TableVII of JCP paper, E=-0.1738 Hartree
  Q(:,1) = [5.663_Rkind,3.821_Rkind,1.842_Rkind]

  ! H...O-O in TableVII of JCP paper, E=-0.1916 Hartree
  Q(:,2) = [2.282_Rkind,7.547_Rkind,9.829_Rkind]

  ! HO2 in TableVII of JCP paper, E=-0.2141 Hartree
  Q(:,3) = [2.806_Rkind,2.271_Rkind,2.271_Rkind]

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  DO i=1,size(Q,dim=2)
    write(out_unitp,*) 'Q:'
    CALL Write_RVec(Q(:,i),out_unitp,QModel%ndim)

    CALL Eval_Pot(QModel,Q(:,i),PotVal,nderiv=nderiv)
    CALL QML_Write_dnMat(PotVal,nio=out_unitp)

    ! For testing the model
    CALL Write_QdnV_FOR_Model(Q(:,i),PotVal,QModel,info='HOO_DMBE')
  END DO
  CALL QML_dealloc_dnMat(PotVal)
  deallocate(Q)
  write(out_unitp,*) ' END Potential and derivatives'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HOO_DMBE',Cart_TO_Q=.TRUE.)
  nderiv = 1

  allocate(x(3,3))
  x(:,1) = [ZERO,ZERO,-7.547_Rkind] ! H1
  x(:,2) = [ZERO,ZERO,ZERO]         ! O2
  x(:,3) = [ZERO,ZERO,2.282_Rkind]  ! O3

  CALL Eval_Pot(QModel,reshape(x,shape=[9]),PotVal,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)

  deallocate(x)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_HOO_DMBE
SUBROUTINE test_H3
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:,:),x(:,:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ 3D-H+H2 ------------------------'
  CALL Init_Model(QModel,pot_name='H3_LSTH',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='H3_LSTH')

  allocate(Q(QModel%ndim,3))

  ! TableII of JCP paper, E=9.802 kcal/mol
  Q(:,1) = [1.757_Rkind,1.757_Rkind,1.757_Rkind]

  Q(:,2) = [1.4_Rkind,1000._Rkind,1001.4_Rkind]

  Q(:,3) = [1000._Rkind,1000._Rkind,2000._Rkind]

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  DO i=1,size(Q,dim=2)
    write(out_unitp,*) 'Q:'
    CALL Write_RVec(Q(:,i),out_unitp,QModel%ndim)

    CALL Eval_Pot(QModel,Q(:,i),PotVal,nderiv=nderiv)
    CALL QML_Write_dnMat(PotVal,nio=out_unitp)

    ! For testing the model
    CALL Write_QdnV_FOR_Model(Q(:,i),PotVal,QModel,info='H3_SLTH')
  END DO
  CALL QML_dealloc_dnMat(PotVal)
  deallocate(Q)
  write(out_unitp,*) ' END Potential and derivatives'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_H3
SUBROUTINE test_Test
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Test potential'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  flush(out_unitp)

  CALL Init_Model(QModel,pot_name='Test',adiabatic=.TRUE.)

  allocate(q(QModel%ndim))
  Q(:) = ZERO

  nderiv=3

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unitp,*) 'Q(:) (bohr):'
  CALL Write_RVec(Q,out_unitp,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Test
SUBROUTINE test_Tully_test
  USE mod_Lib
  USE mod_dnMat
  USE mod_Model
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 3
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Tully potential'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'

    CALL Init_Model(QModel,pot_name='Tully',option=option,adiabatic=.TRUE.)
    allocate(q(QModel%QM%ndim))
    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,*) '----- CHECK POT -----------------------------'
    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,'(a,i1)') ' ----- TULLY',option

    write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

    q(:) = 1.1_Rkind
    write(out_unitp,'(a,f12.6)') ' R (Bohr)',q(:)
    CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)



  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Tully_test
SUBROUTINE test_PSB3_test

  USE mod_Lib
  USE mod_dnMat
  USE mod_Model

  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' PSB3 potential'
  write(out_unitp,*) ' With units: Atomic Units (Angstrom, Rad, Rad, kcal.mol^-1)'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='PSB3',Print_init=.FALSE.)
  CALL Write0_Model(QModel)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='PSB3',PubliUnit=.FALSE.)

  allocate(q(QModel%QM%ndim))

  q(:) = [0.172459_Rkind,PI,ZERO]
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  write(out_unitp,*) '---------------------------------------------'

  QModel%QM%adiabatic = .TRUE.
  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal,nio=out_unitp)


  CALL QML_dealloc_dnMat(PotVal)
  deallocate(q)

END SUBROUTINE test_PSB3_test
SUBROUTINE test_PSB3_Retinal2000_test

  USE mod_Lib
  USE mod_dnMat
  USE mod_Model

  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel_PSB3
  TYPE (Model_t)                 :: QModel_Retinal2000

  real (kind=Rkind), allocatable :: q_PSB3(:)
  real (kind=Rkind), allocatable :: q_Retinal2000(:)
  integer                        :: ndim,nsurf,nderiv,i,option

  TYPE (dnMat_t)                 :: PotVal_PSB3,PotVal_Retinal2000
  TYPE (dnMat_t)                 :: NAC_PSB3,NAC_Retinal2000

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' PSB3 potential'
  write(out_unitp,*) ' With units: Atomic Units'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel_PSB3,pot_name='PSB3',PubliUnit=.FALSE.)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Retinal_JPCB2000 potential'
  write(out_unitp,*) ' With units: Atomic Units'
  write(out_unitp,*) '---------------------------------------------'
  CALL Init_Model(QModel_Retinal2000,pot_name='Retinal_JPCB2000',PubliUnit=.FALSE.)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'


  allocate(q_Retinal2000(QModel_Retinal2000%QM%ndim))
  q_Retinal2000(:) = [ZERO,PI]

  allocate(q_PSB3(QModel_PSB3%QM%ndim))
  q_PSB3(:) = [0.172459_Rkind,PI,ZERO]
  nderiv=1

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'
  write(out_unitp,*) '    ADIABATIC potential'
  write(out_unitp,*) '---------------------------------------------'

  QModel_Retinal2000%QM%adiabatic = .TRUE.
  QModel_PSB3%QM%adiabatic = .TRUE.

  write(out_unitp,*) 'Evaluated in', q_PSB3
  CALL Eval_Pot(QModel_PSB3,q_PSB3,PotVal_PSB3,NAC=NAC_PSB3,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal_PSB3,nio=out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) 'Evaluated in', q_Retinal2000
  CALL Eval_Pot(QModel_Retinal2000,q_Retinal2000,PotVal_Retinal2000,            &
                    NAC=NAC_Retinal2000,nderiv=nderiv)
  CALL QML_Write_dnMat(PotVal_Retinal2000,nio=out_unitp)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  QModel_PSB3%QM%adiabatic = .TRUE.
  CALL Eval_pot_ON_Grid(QModel_PSB3,Qmin=[0.172459_Rkind,-PI,ZERO],             &
                                    Qmax=[0.172459_Rkind,PI,ZERO],nb_points=1001, &
                        grid_file='grid_PSB3')

  QModel_Retinal2000%QM%adiabatic = .TRUE.
  CALL Eval_pot_ON_Grid(QModel_Retinal2000,Qmin=[-PI,ZERO],             &
                                           Qmax=[PI,ZERO],nb_points=1001, &
                        grid_file='grid_Retinal2000')

  CALL QML_dealloc_dnMat(PotVal_Retinal2000)
  CALL QML_dealloc_dnMat(NAC_Retinal2000)
  deallocate(q_Retinal2000)

  CALL QML_dealloc_dnMat(PotVal_PSB3)
  CALL QML_dealloc_dnMat(NAC_PSB3)
  deallocate(q_PSB3)

END SUBROUTINE test_PSB3_Retinal2000_test
