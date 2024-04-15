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
  IMPLICIT NONE

  TYPE (test_t)                  :: test_var

  CALL Initialize_Test(test_var,test_name='QModel')
  !CALL test_TwoD_RJDI2014() ; CALL Finalize_Test(test_var) ; stop

  !CALL test_PSB3() ; CALL Finalize_Test(test_var) ; stop

  !CALL test_Bottleneck ; CALL Finalize_Test(test_var) ; stop
  !CALL test_PH4Jo ; stop
  !CALL test_H3() ; stop

  ! One electronic surface
  CALL test_LinearHBond()
  CALL test_Opt_LinearHBond()
  CALL test_HenonHeiles()
  CALL test_Buckingham()
  CALL test_Morse()
  CALL test_Poly1D()

  ! One electronic surface + optimization + IRC
  CALL test_Opt_MullerBrown()
  CALL test_IRC_MullerBrown()

  ! One electronic surface + optimization
  !CALL test_H3()

  ! Several electronic surfaces
  CALL test_Tully()
  !CALL test_OneDSOC_1S1T() ! this test is removed, because it gives always 2 errors (degenerate eigenvectors)
  CALL test_OneDSOC_2S1T()

  CALL test_Phenol()
  CALL test_PSB3()
  CALL test_TwoD()
  CALL test_TwoD_RJDI2014()
  CALL test_Vibronic()
  CALL test_TwoD_Valahu2022()
  CALL test_OneD_Photons() 
  CALL test_Retinal_JPCB2000()
  CALL test_NO3()

  ! 6D (full-D), One electronic surface (spectro)
  CALL test_HONO()
  CALL test_HNNHp()

  CALL test_H2SiN()
  CALL test_H2NSi()

  CALL test_H2O() ! for testing (PES: quadratic expansion)

  ! 3D (full-D), One electronic surface (spectro); ClH2+ (unpublished)
  CALL test_ClH2p_op12()
  CALL test_ClH2p_op34()
  CALL test_ClH2p_op56()
  ! 3D (full-D), One electronic surface (spectro); ClH2+ (From Botschwina 1988)
  CALL test_ClH2p_Botschwina()


  ! 3D (full-D), One electronic surface (spectro): CNH Murrell
  CALL test_CNH_Murrell

  ! 9D (full-D), One electronic surface (spectro, torsional levels): HO-NO2 
  CALL test_HNO3()

  ! 3D (full-D), One electronic surface for collision
  CALL test_HOO_DMBE()

  ! 12D (full-D), One electronic surface for collision H+CH4 -> H2+CH3
  CALL test_CH5()
  ! 9D (full-D), One electronic surface for collision H+PH3 -> H2+PH2
  CALL test_PH4Jo()

  ! A template with one electronic surface
  CALL test_template()

  ! vibrational adiabatic separation (on HBond potential)
  CALL test_Vib_adia()
  CALL test2_Vib_adia()

  CALL Finalize_Test(test_var)

CONTAINS

SUBROUTINE test_Tully
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : TO_string
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Tully potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'

  DO option=1,3
    CALL Init_Model(QModel,pot_name='Tully',Print_init=.FALSE.,option=option)
    allocate(q(QModel%QM%ndim))
    write(out_unit,*) '---------------------------------------------'
    write(out_unit,*) '----- CHECK POT -----------------------------'
    write(out_unit,*) '---------------------------------------------'
    write(out_unit,'(a,i1)') ' ----- TULLY',option

    write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

    q(:) = 1.1_Rkind
    write(out_unit,'(a,f12.6)') ' R (Bohr)',q(:)
    CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

    q(:) = -1.1_Rkind
    write(out_unit,'(a,f12.6)') ' R (Bohr)',q(:)
    CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)


    write(out_unit,*) '---------------------------------------------'
    write(out_unit,*) ' Potential and derivatives'

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
    write(out_unit,'(a,f12.6)') 'R (Bohr)',q(:)
    write(out_unit,*) 'Energy (Hartree)'
    CALL Write_dnMat(PotVal,nio=out_unit)
    ! For testing the model
    CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Tully ' // TO_string(option), &
                             test_var=test_var,last_test=(option == 3))

    write(out_unit,*) '---------------------------------------------'
    write(out_unit,*) '- END CHECK POT -----------------------------'
    write(out_unit,*) '---------------------------------------------'

    CALL dealloc_dnMat(PotVal)
    deallocate(q)
    CALL dealloc_Model(QModel)
  END DO

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '-----------TULLY 1 on a grid ----------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '   file name: "grid_Tully1"'
  write(out_unit,*) ' You should get the curves (columns: 2,3 5) of the panel (a) of figure 3 of ... '
  write(out_unit,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unit,*) ' Be carrefull, the non-adiabatic coupling is scaled by -1/50.'
  CALL Init_Model(QModel,pot_name='Tully',Print_init=.FALSE.,option=1)
  CALL Eval_pot_ON_Grid(QModel,Qmin=[-TEN],Qmax=[TEN],nb_points=1001,grid_file='grid_Tully1')

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '-----------TULLY 2 on a grid ----------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '   file name: "grid_Tully2"'
  write(out_unit,*) ' You should get the curves (columns: 2,3 5) of the panel (b) of figure 3 of ... '
  write(out_unit,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unit,*) ' Be carrefull, the non-adiabatic coupling is scaled by -1/12.'

  CALL Init_Model(QModel,pot_name='Tully',Print_init=.FALSE.,option=2)
  CALL Eval_pot_ON_Grid(QModel,Qmin=[-TEN],Qmax=[TEN],nb_points=1001,grid_file='grid_Tully2')

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '-----------TULLY 3 on a grid ----------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '   file name: "grid_Tully3"'
  write(out_unit,*) ' You should get the curves (columns: 2,3 5) of the panel (c) of figure 3 of ... '
  write(out_unit,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unit,*) ' Be carrefull, the sign of non-adiabatic coupling is changed.'

  CALL Init_Model(QModel,pot_name='Tully',Print_init=.FALSE.,option=3)
  CALL Eval_pot_ON_Grid(QModel,Qmin=[-TEN],Qmax=[TEN],nb_points=1001,grid_file='grid_Tully3')

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_Tully

SUBROUTINE test_OneDSOC_1S1T
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' 1D-SOC_1S1T potential'
  write(out_unit,*) ' With units: Bohr Hartree (au)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,nsurf=4,pot_name='1DSOC_1S1T')

  allocate(q(QModel%QM%ndim))
  q(:) = [ 9.5_Rkind ] ! one evaluation to get vec0%d0(:,:)
  write(out_unit,*) ' ref Eigenvectors at Q:',Q
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=0)
  ! ref Eigenvectors
  CALL Write_dnMat(QModel%QM%vec0,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,Vec=QModel%QM%vec0,info='1DSOC_1S1T', &
                           test_var=test_var,last_test=.FALSE.)


  q(:) = [8.5_Rkind]
  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Q (Bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) 'Non adiatic couplings:'
  CALL Write_dnMat(NAC,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='1DSOC_1S1T', &
                           test_var=test_var,last_test=.FALSE.)

  write(out_unit,*) 'DIABATIC potential'
  QModel%QM%adiabatic = .FALSE.
  write(out_unit,*) 'Q (Bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)
  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='1DSOC_1S1T', &
                           test_var=test_var,last_test=.TRUE.)

  CALL dealloc_dnMat(PotVal)
  CALL dealloc_dnMat(NAC)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----------- 1D-SOC_1S1T on a grid -----------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of R in Bohr)'
  write(out_unit,*) '   file name: "grid_1DSOC_1S1T"'
  write(out_unit,*) ' You should get the 3 curves (left panels) of figure 1 of ... '
  write(out_unit,*) '  ... Granucci et al. J. Chem. Phys. V137, p22A501 (2012)'

  CALL Init_Model(QModel,nsurf=4,pot_name='1DSOC_1S1T',option=1)
  QModel%QM%adiabatic = .TRUE.
  flush(out_unit)

  q(:) = [ 9.5_Rkind ] ! one evaluation to get vec0%d0(:,:)
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=1)

  ! ref Eigenvectors
  write(out_unit,*) ' ref Eigenvectors at Q:',Q
  CALL Write_dnMat(QModel%QM%vec0,nio=out_unit)


  CALL Eval_pot_ON_Grid(QModel,Qmin=[3._Rkind],Qmax=[20._Rkind], &
                        nb_points=1001,nderiv=0,grid_file='grid_1DSOC_1S1T')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  deallocate(q)
  CALL dealloc_dnMat(PotVal)
  CALL dealloc_dnMat(NAC)
  CALL dealloc_Model(QModel)

END SUBROUTINE test_OneDSOC_1S1T
SUBROUTINE test_OneDSOC_2S1T
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' 1D-SOC_2S1T potential'
  write(out_unit,*) ' With units: Bohr Hartree (au)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='1DSOC_2S1T',Print_init=.FALSE.)

  allocate(q(QModel%QM%ndim))
  q(:) = [ 9.5_Rkind ] ! one evaluation to get vec0%d0(:,:)
  write(out_unit,*) ' ref Eigenvectors at Q:',Q
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=0)
  ! ref Eigenvectors
  CALL Write_dnMat(QModel%QM%vec0,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,Vec=QModel%QM%vec0,info='1DSOC_2S1T', &
                           test_var=test_var,last_test=.FALSE.)


  q(:) = [8.5_Rkind]
  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Q (Bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) 'Non adiatic couplings:'
  CALL Write_dnMat(NAC,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='1DSOC_2S1T', &
                           test_var=test_var,last_test=.FALSE.)
  write(out_unit,*) 'DIABATIC potential'
  QModel%QM%adiabatic = .FALSE.
  write(out_unit,*) 'Q (Bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)
  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='1DSOC_2S1T', &
                           test_var=test_var,last_test=.TRUE.)

  CALL dealloc_dnMat(PotVal)
  CALL dealloc_dnMat(NAC)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----------- 1D-SOC_2S1T on a grid -----------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of R in Bohr)'
  write(out_unit,*) '   file name: "grid_1DSOC_2S1T"'

  CALL Init_Model(QModel,pot_name='1DSOC_2S1T')
  QModel%QM%adiabatic = .TRUE.
  flush(out_unit)

  q(:) = [ 9.5_Rkind ] ! one evaluation to get vec0%d0(:,:)
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=1)

  ! ref Eigenvectors
  write(out_unit,*) ' ref Eigenvectors at Q:',Q
  CALL Write_dnMat(QModel%QM%vec0,nio=out_unit)


  CALL Eval_pot_ON_Grid(QModel,Qmin=[3._Rkind],Qmax=[20._Rkind], &
                        nb_points=1001,nderiv=0,grid_file='grid_1DSOC_2S1T')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  deallocate(q)
  CALL dealloc_dnMat(PotVal)
  CALL dealloc_dnMat(NAC)
  CALL dealloc_Model(QModel)

END SUBROUTINE test_OneDSOC_2S1T
SUBROUTINE test_Morse
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Morse potential (H-F parameters)'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Morse',read_param=.FALSE.)

  allocate(q(QModel%QM%ndim))
  q(:) = 2._Rkind

  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,'(a,f12.6)') 'R (Bohr)',q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unit,'(a,f12.6)') 'R (Bohr)',q(:)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Morse_HF', &
                           test_var=test_var,last_test=.TRUE.)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)
  deallocate(q)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of R)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '   file name: "grid_Morse"'

  CALL Eval_pot_ON_Grid(QModel,Qmin=[1._Rkind],Qmax=[5._Rkind],nb_points=1001,&
                        grid_file='grid_Morse')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL dealloc_Model(QModel)

END SUBROUTINE test_Morse
SUBROUTINE test_Poly1D
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: nderiv
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Polynomial potential (H-F parameters)'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='Poly1D',Print_init=.TRUE.)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(Q(QModel%QM%ndim))
  Q(:) = QModel%QM%Q0

  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,'(a,f12.6)') 'R (Bohr)',Q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unit,'(a,f12.6)') 'R (Bohr)',Q(:)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Poly1D_HF',test_var=test_var,last_test=.TRUE.)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

END SUBROUTINE test_Poly1D
SUBROUTINE test_H2
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: nderiv
  TYPE (dnMat_t)                 :: PotVal

  TYPE (QML_Opt_t)               :: Opt_p

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' H2 potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='H2',Print_init=.TRUE.,option=3)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(Q(QModel%QM%ndim))
  Q(:) = QModel%QM%Q0

  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,'(a,f12.6)') 'R (Bohr)',Q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'

  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=9)
  CALL QML_Opt(Q,QModel,Opt_p,Q0=[1.4_Rkind])

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unit,'(a,f12.6)') 'R (Bohr)',Q(:)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  !CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Poly1D_HF',test_var=test_var,last_test=.TRUE.)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of R)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '   file name: "grid_H2"'

  CALL Eval_pot_ON_Grid(QModel,Qmin=[0.7_Rkind],Qmax=[20._Rkind],nb_points=1001,grid_file='grid_H2')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

END SUBROUTINE test_H2
SUBROUTINE test_Buckingham
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Buckingham potential (Ar-Ar parameters)'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='Buck')
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Buck',read_param=.FALSE.)

  allocate(q(QModel%QM%ndim))
  q(:) = 7._Rkind
  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,'(a,f12.6)') 'R (Bohr)',q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unit,'(a,f12.6)') 'R (Bohr)',q(:)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Buckingham_Ar-Ar', &
                           test_var=test_var,last_test=.TRUE.)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)
  deallocate(q)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of R)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '   file name: "grid_Buck"'

  CALL Eval_pot_ON_Grid(QModel,Qmin=[6._Rkind],Qmax=[20._Rkind],nb_points=1001,&
                        grid_file='grid_Buck')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL dealloc_Model(QModel)

END SUBROUTINE test_Buckingham

SUBROUTINE test_Phenol
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Phenol potential'
  write(out_unit,*) ' With units: Angs and eV'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Phenol',PubliUnit=.TRUE.,Print_init=.FALSE.)

  allocate(q(QModel%QM%ndim))

  q(:) = [1.2_Rkind,0.2_Rkind]
  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q (Angs, radian):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Q (Angs, radian):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) 'Non adiatic couplings:'
  CALL Write_dnMat(NAC,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='Phenol', &
      test_var=test_var,last_test=.FALSE.)


  write(out_unit,*) 'DIABATIC potential'
  QModel%QM%adiabatic = .FALSE.
  write(out_unit,*) 'Q (Angs, radian):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Phenol', &
      test_var=test_var,last_test=.TRUE.)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with PubliUnit=.FALSE.'
  CALL Init_Model(QModel,pot_name='Phenol',PubliUnit=.FALSE.,Print_init=.FALSE.)
  CALL get_Q0_Model(q,QModel,option=0)
  write(out_unit,'(a,2f12.6)') 'q (au)',q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)
  CALL dealloc_dnMat(NAC)
  deallocate(q)
  CALL dealloc_Model(QModel)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----------- Phenol on a grid ----------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of R in Angstrom)'
  write(out_unit,*) '   file name: "grid_Phenol"'
  write(out_unit,*) ' You should get the 3 curves (columns: 2,6, 10) of figure 2 of ... '
  write(out_unit,*) '   Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, ...'
  write(out_unit,*) '  .... J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218)'

  CALL Init_Model(QModel,pot_name='phenol',PubliUnit=.TRUE.,Print_init=.FALSE.,adiabatic = .FALSE.)

  CALL Eval_pot_ON_Grid(QModel,Qmin=[0.5_Rkind,ZERO],Qmax=[5._Rkind,ZERO], &
                        nb_points=1001,grid_file='grid_Phenol')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL dealloc_Model(QModel)

END SUBROUTINE test_Phenol
SUBROUTINE test_HenonHeiles
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 4D-HenonHeiles -----------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HenonHeiles',ndim=4)

  allocate(q(QModel%QM%ndim))
  q(:) = [ (0.1_Rkind*real(i,kind=Rkind),i=1,QModel%QM%ndim) ]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='4D-HenonHeiles', &
      test_var=test_var,last_test=.TRUE.)

  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  nderiv = 3
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 2D-HenonHeiles -----------------'
  write(out_unit,*) '-------------option 1, 2 and 3 --------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HenonHeiles',ndim=2,option=1)

  Q = [ ZERO,ZERO ]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='HenonHeiles',ndim=2,option=2)
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='HenonHeiles',ndim=2,option=3)
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
END SUBROUTINE test_HenonHeiles
SUBROUTINE test_Bottleneck
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 2D-test_Bottleneck -------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Bottleneck',ndim=2, option=3)

  Q = [ (0.1_Rkind*real(i,kind=Rkind),i=1,QModel%QM%ndim) ]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  !CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='4D-HenonHeiles', &
   !   test_var=test_var,last_test=.TRUE.)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' 1D-grid'
  write(out_unit,*) '---------------------------------------------'
  CALL Eval_pot_ON_Grid(QModel, &
  Qmin=[-5._Rkind,ZERO],Qmax=[5._Rkind,ZERO],nb_points=101,grid_file='grid_bottleneck')

  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_Bottleneck
SUBROUTINE test_LinearHBond
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Linear H-Bond (symmetric one and the default parameters)'
  write(out_unit,*) ' With units: Angs and kcal.mol^-1'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HBond',PubliUnit=.TRUE.,Print_init=.FALSE.)

  allocate(Q(QModel%QM%ndim))

  Q(:) = [ 2.75_Rkind,0._Rkind ]
  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,'(a,2f12.6)') 'QQ,q (Angs)',Q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unit,'(a,2f12.6)') 'QQ,q (Angs)',q(:)
  write(out_unit,*) 'Energy (kcal.mol^-1)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Linear H-Bond (symmetric one)',test_var=test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with PubliUnit=.FALSE.'
  CALL Init_Model(QModel,pot_name='HBond',PubliUnit=.FALSE.,Print_init=.FALSE.)
  CALL get_Q0_Model(Q,QModel,option=0)
  write(out_unit,'(a,2f12.6)') 'QQ,q (Bohr)',Q(:)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL dealloc_dnMat(PotVal)

  deallocate(Q)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unit,*) '   file name: "grid_Hbond-sym"'
  write(out_unit,*) ' You should get the green curve (QQ=2.75 Angs) of ... '
  write(out_unit,*) '...the bottom left panel of figure 3.8 (Julien Beutier thesis)'

  CALL Eval_pot_ON_Grid(QModel, &
                        Qmin=[2.75_Rkind,-0.6_Rkind],Qmax=[2.75_Rkind,0.6_Rkind],nb_points=1001,&
                        grid_file='grid_Hbond-sym')

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Linear H-Bond (asymmetric one)'
  write(out_unit,*) ' With units: Angs and kcal.mol^-1'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HBond',PubliUnit=.TRUE.,Print_init=.FALSE.,  &
                  read_param=.TRUE.,param_file_name='DAT_files/Hbond.dat')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  allocate(Q(QModel%QM%ndim))

  Q(:) = [3._Rkind,0._Rkind]
  nderiv=2
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  write(out_unit,'(a,2f12.6)') 'QQ,q (Angs)',q(:)
  write(out_unit,*) 'Energy (kcal.mol^-1)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Linear H-Bond (asymmetric one)', &
      test_var=test_var,last_test=.TRUE.)

  CALL dealloc_dnMat(PotVal)

  deallocate(Q)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unit,*) '   file name: "grid_Hbond-asym"'
  write(out_unit,*) ' You should get the green curve (QQ=3. Angs) of ... '
  write(out_unit,*) '...the bottom right panel of figure 3.8 (Julien Beutier thesis)'

  CALL Eval_pot_ON_Grid(QModel, &
                        Qmin=[3._Rkind,-0.8_Rkind],Qmax=[3._Rkind,0.7_Rkind],nb_points=1001,&
                        grid_file='grid_Hbond-asym')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL dealloc_Model(QModel)

END SUBROUTINE test_LinearHBond
SUBROUTINE test_Opt_LinearHBond
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  TYPE (QML_Opt_t)               :: Opt_param

  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Optimization with HBond potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)
  CALL Init_Model(QModel,pot_name='HBond',PubliUnit=.FALSE.)


  CALL Init_QML_Opt(Opt_param,QModel,read_param=.TRUE.,param_file_name='DAT_files/opt_hbond.dat')

  allocate(Q(QModel%ndim))


  CALL QML_Opt(Q,QModel,Opt_param,Q0=[ 5._Rkind,0.1_Rkind ])

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  write(out_unit,*) 'Potential+derivatives:'
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_Opt_LinearHBond
SUBROUTINE test_Vib_adia
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Mat
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Qact(:)
  TYPE (dnMat_t)                 :: PotVal,NAC
  integer                        :: nderiv

  real(kind=Rkind) :: dQ
  integer :: i,iq
  integer, parameter :: nq=100
  real(kind=Rkind), parameter    :: auTOcm_inv = 219474.631443_Rkind

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Vibrational adiabatic separation with ...'
  write(out_unit,*) ' Linear H-Bond (symmetric one and the default parameters)'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,Print_init=.TRUE.,Vib_adia=.TRUE.,PubliUnit=.TRUE.,   &
                  param_file_name='DAT_files/Vibadia_HBond.dat')
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  Qact = [4._Rkind]
  CALL Eval_Pot(QModel,Qact,PotVal,NAC=NAC,nderiv=1)
  write(out_unit,*) 'NAC'
  CALL Write_Mat(NAC%d1(:,:,1),out_unit,6,info='NAC')
  write(out_unit,*) Qact,'Ene',(PotVal%d0(i,i)*auTOcm_inv,i=1,get_nsurf(PotVal))

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Qact,PotVal,QModel,info='Pot',name_file='Vib_adia.txt', &
      test_var=test_var,last_test=.FALSE.)

  CALL Test_QdnV_FOR_Model(Qact,NAC,QModel,info='NAC',name_file='Vib_adia.txt', &
      test_var=test_var,last_test=.TRUE.)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Adiatic potential on a 1D grid (as a function of QQ)'
  write(out_unit,*) '---------------------------------------------'

  dQ = TWO / real(nq,kind=Rkind)
  DO iq=0,nq
    Qact = [4.0_Rkind+iq*dQ]
    CALL Eval_Pot(QModel,Qact,PotVal,nderiv=0)
    write(out_unit,*) Qact,'Ene',(PotVal%d0(i,i)*auTOcm_inv,i=1,get_nsurf(PotVal))
  END DO
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)
  deallocate(Qact)
  CALL dealloc_Model(QModel)

END SUBROUTINE test_Vib_adia
SUBROUTINE test2_Vib_adia
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Mat
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Qact(:)
  real (kind=Rkind), allocatable :: tab_MatH(:,:,:)

  TYPE (dnMat_t)                 :: PotVal,NAC
  integer                        :: nderiv

  integer            :: i,iq,iterm
  integer, parameter :: nq=100
  real(kind=Rkind)   :: dQ
  real(kind=Rkind), parameter    :: auTOcm_inv = 219474.631443_Rkind

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Vibrational adiabatic separation with ...'
  write(out_unit,*) ' Linear H-Bond (symmetric one and the default parameters)'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,Print_init=.TRUE.,Vib_adia=.TRUE.,PubliUnit=.TRUE.,   &
                  param_file_name='DAT_files/Vibadia_HBond.dat')
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  Qact = [4.0_Rkind]
  CALL Eval_tab_HMatVibAdia(QModel,Qact,tab_MatH)

  DO iterm=1,size(tab_MatH,dim=3)
    CALL Write_Mat(tab_MatH(:,:,iterm),nio=out_unit,nbcol=5)
  END DO

  dQ = TWO / real(nq,kind=Rkind)
  DO iq=0,nq
    Qact = [4.0_Rkind+iq*dQ]
    CALL Eval_Pot(QModel,Qact,PotVal,NAC=NAC,nderiv=1)
    write(out_unit,*) Qact,'Ene',(PotVal%d0(i,i),i=1,get_nsurf(PotVal))
    write(out_unit,*) Qact,'NAC1',NAC%d1
  END DO

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  deallocate(Qact)
  CALL dealloc_Model(QModel)

END SUBROUTINE test2_Vib_adia
SUBROUTINE test_PSB3
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m

  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' PSB3 potential'
  write(out_unit,*) ' With units: Atomic Units (Bohr, Rad, Rad, Hartree)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='PSB3',PubliUnit=.FALSE.)

  allocate(q(QModel%QM%ndim))

  q(:) = [0.172459_Rkind,-3.14_Rkind,0._Rkind]
  nderiv=3

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Evaluated in', q
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'

  QModel%QM%adiabatic = .FALSE.
  write(out_unit,*) 'DIABATIC potential'
  write(out_unit,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='PSB3', &
      test_var=test_var,last_test=.FALSE.)

  QModel%QM%adiabatic = .TRUE.
  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) 'Non adiatic couplings:'
  CALL Write_dnMat(NAC,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='PSB3', &
      test_var=test_var,last_test=.TRUE.)

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with PubliUnit=.FALSE.'
  DO option=1,2
    CALL Init_Model(QModel,pot_name='PSB3',PubliUnit=.FALSE.,Print_init=.FALSE.,option=option)
    CALL get_Q0_Model(q,QModel,option=0)
    write(out_unit,*) 'q (au)',q(:)
    CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)
    CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
    CALL Write_dnMat(PotVal,nio=out_unit)
    CALL dealloc_Model(QModel)
  END DO
  write(out_unit,*) '---------- END CHECK POT --------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)

  deallocate(q)

END SUBROUTINE test_PSB3
SUBROUTINE test_Retinal_CP2000
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m

  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Retinal_CB2000 potential (2+nD)'
  write(out_unit,*) ' With units: Atomic Units'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,PubliUnit=.FALSE.,read_param=.TRUE.,  &
                  param_file_name='DAT_files/retinal_cp2000.dat')

  allocate(q(QModel%QM%ndim))

  q(:) = [HALF,ONE,-ONE]
  nderiv=1

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Evaluated in',q
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'

  QModel%QM%adiabatic = .FALSE.
  write(out_unit,*) 'DIABATIC potential'
  write(out_unit,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,q,PotVal,nderiv=nderiv)

  CALL Write_dnMat(PotVal,nio=out_unit)


  QModel%QM%adiabatic = .TRUE.
  q(:) = [ZERO,ZERO,ZERO]
  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) 'Non adiatic couplings:'
  CALL Write_dnMat(NAC,nio=out_unit)

  q(:) = [THREE,ZERO,ZERO]
  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Evaluated in', q
  CALL Eval_Pot(QModel,q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL dealloc_dnMat(PotVal)
  CALL dealloc_dnMat(NAC)
  deallocate(q)
  CALL dealloc_Model(QModel)

END SUBROUTINE test_Retinal_CP2000
SUBROUTINE test_Retinal_JPCB2000
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m

  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings
  real (kind=Rkind)              :: DQ2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Retinal_JPCB2000 potential'
  write(out_unit,*) ' With units: Atomic Units'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='Retinal_JPCB2000')
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='Retinal_JPCB2000',PubliUnit=.FALSE.,Print_init=.FALSE.)

  allocate(q(QModel%QM%ndim))

  q(:) = [ZERO,ONE]
  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Evaluated in', q
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'

  QModel%QM%adiabatic = .FALSE.
  write(out_unit,*) 'DIABATIC potential'
  write(out_unit,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Retinal_JPCB2000',test_var=test_var)

  QModel%QM%adiabatic = .TRUE.
  q(:) = [ZERO,ZERO]
  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) 'Non adiatic couplings:'
  CALL Write_dnMat(NAC,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='Retinal_JPCB2000',test_var=test_var)

  write(out_unit,*) '---------- END CHECK POT --------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '-----------Retinal_JPCB2000 on a 2D-grid ----'
  write(out_unit,*) '---------------------------------------------'

  QModel%QM%adiabatic = .TRUE.
  DQ2 = TWO
  CALL Eval_pot_ON_Grid(QModel,Qmin=[-pi/TWO,-DQ2],Qmax=[THREE*pi/TWO,DQ2],         &
                        nb_points=101,grid_file='grid_Retinal_JPCB2000')


  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Retinal_JPCB2000 potential+Bath'
  write(out_unit,*) ' With units: Atomic Units'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,ndim=4,pot_name='Retinal_JPCB2000',PubliUnit=.FALSE.,Print_init=.FALSE.)

  allocate(q(QModel%QM%ndim))

  q(:) = ZERO
  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Evaluated in',q

  QModel%QM%adiabatic = .TRUE.
  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) 'Non adiatic couplings:'
  CALL Write_dnMat(NAC,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,NAC=NAC,info='Retinal_JPCB2000', &
      test_var=test_var,last_test=.TRUE.)

  write(out_unit,*) '---------- END CHECK POT --------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

END SUBROUTINE test_Retinal_JPCB2000
SUBROUTINE test_HONO
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 6D-HONO ------------------------'
  CALL Init_Model(QModel,pot_name='HONO',Print_init=.TRUE.)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HONO')

  allocate(q(QModel%QM%ndim))
  q(:) = [2.696732586_Rkind,1.822912197_Rkind,1.777642018_Rkind,2.213326419_Rkind,1.9315017_Rkind,pi]
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the trans minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='HONO',test_var=test_var)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the cis minimum'
  q(:) = [2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,2.23738_Rkind,1.975200_Rkind,ZERO]

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='HONO', &
      test_var=test_var,last_test=.TRUE.)


  CALL dealloc_dnMat(PotVal)
  deallocate(q)



  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of q)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential on a 1D grid (as a function of q(6), the torsion)'
  write(out_unit,*) '   file name: "grid_HONO"'

  CALL Eval_pot_ON_Grid(QModel,                                      &
                        Qmin=[2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,&
                              2.23738_Rkind,1.975200_Rkind,ZERO],        &
                        Qmax=[2.63122_Rkind,1.84164_Rkind,1.822274_Rkind,&
                              2.23738_Rkind,1.975200_Rkind,pi],          &
                        nb_points=1001, grid_file='grid_HONO')

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_HONO
SUBROUTINE test_HNNHp
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:),xyz(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  real (kind=Rkind) :: th


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 6D-HNNH+ -----------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ Z-matrix coordinates -----------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HNNHp',Print_init=.FALSE.)

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the trans minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='HNNHp', &
      test_var=test_var,last_test=.FALSE.)

  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ Cartesian coordinates ----------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HNNHp',Cart_TO_Q=.TRUE.)

  allocate(xyz(QModel%QM%ndim))
  xyz = [ZERO,              ZERO,  ZERO,              &
         ZERO,              ZERO,  2.200000000_Rkind, &
         1.640097797_Rkind, ZERO, -0.959207599_Rkind, &
         0.869189803_Rkind, 1.353682913_Rkind,  2.749592264_Rkind]

  CALL Eval_Pot(QModel,xyz,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

    ! For testing the model
  CALL Test_QdnV_FOR_Model(xyz,PotVal,QModel,info='HNNHp', &
      test_var=test_var,last_test=.TRUE.)

  CALL dealloc_dnMat(PotVal)
  deallocate(xyz)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_HNNHp
SUBROUTINE test_H2SiN
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  DO option=1,3
  write(out_unit,*) '------------ 6D-H2SiN -----------------------'
  write(out_unit,'(a,i0,a)') ' ------------ option: ',option,' ----------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='H2SiN',option=option,Print_init=.FALSE.)

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at:'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='H2SiN', &
      test_var=test_var,last_test=(option == 3))


  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
END DO

END SUBROUTINE test_H2SiN
SUBROUTINE test_H2NSi
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
DO option=1,2
  write(out_unit,*) '------------ 6D-H2NSi -----------------------'
  write(out_unit,'(a,i0,a)') ' ------------ option: ',option,' ----------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL Init_Model(QModel,pot_name='H2NSi',option=option,Print_init=.FALSE.)

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at:'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='H2NSi', &
      test_var=test_var,last_test=(option == 2))


  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
END DO

END SUBROUTINE test_H2NSi

SUBROUTINE test_H2O
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotVal_gaussian
  real (kind=Rkind), allocatable :: qtest(:,:),EAbInitio(:)


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-H2O (option 1: R1,R2,a)------'
  CALL Init_Model(QModel,pot_name='H2O',Print_init=.TRUE.,option=1)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='H2O_op1', &
      test_var=test_var,last_test=.TRUE.)

  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_H2O
SUBROUTINE test_NO3
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotVal_gaussian
  real (kind=Rkind), allocatable :: qtest(:,:),EAbInitio(:)


  nderiv = 1
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 6D-NO3 -------------------------'
  CALL Init_Model(QModel,pot_name='NO3',Print_init=.TRUE.,adiabatic=.FALSE.)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(Q(QModel%QM%ndim))
  Q( 1: 3) = [-0.28276917777022903_Rkind,    -0.41734715646392218_Rkind,-1.0997467054058241e-002_Rkind]
  Q( 4: 6) = [-9.1098922338089527e-002_Rkind,-0.13445551772839726_Rkind, 2.2271508290107942_Rkind]
  Q( 7: 9) = [ 2.2747337906099050_Rkind,     -0.13445551772839726_Rkind,-0.81179512186018443_Rkind]
  Q(10:12) = [-1.9360788283195443_Rkind,      0.63428610976387545_Rkind,-1.4057277504727415_Rkind]
  !CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='NO3', &
      test_var=test_var,last_test=.TRUE.)

  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_NO3
SUBROUTINE test_ClH2p_Botschwina
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotVal_gaussian
  real (kind=Rkind), allocatable :: qtest(:,:),EAbInitio(:)


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-ClH2+ Botschwina ------------'
  write(out_unit,*) '------------ option 1: R1, R2, A ------------'
  CALL Init_Model(QModel,pot_name='ClH2p_Botschwina',Print_init=.TRUE.,option=1)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='ClH2p_Botschwina', &
      test_var=test_var,last_test=.FALSE.)

  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-ClH2+ Botschwina corrected --'
  write(out_unit,*) '------------ option 2: R1, R2, A ------------'
  CALL Init_Model(QModel,pot_name='ClH2p_Botschwina',Print_init=.TRUE.,option=2)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='ClH2p_Botschwina', &
      test_var=test_var,last_test=.TRUE.)

  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'


END SUBROUTINE test_ClH2p_Botschwina
SUBROUTINE test_ClH2p_op12
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotVal_gaussian
  real (kind=Rkind), allocatable :: qtest(:,:),EAbInitio(:)


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-ClH2+ (option 1: a,R+,R-)----'
  CALL Init_Model(QModel,pot_name='ClH2p',Print_init=.TRUE.,option=1)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='ClH2p_op12',test_var=test_var)

  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-ClH2+ (option 2: R1,R2,a)----'
  CALL Init_Model(QModel,pot_name='ClH2p',Print_init=.TRUE.,option=2)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='ClH2p_op12',test_var=test_var)

  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'



END SUBROUTINE test_ClH2p_op12
SUBROUTINE test_ClH2p_op34
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotVal_gaussian
  real (kind=Rkind), allocatable :: qtest(:,:),EAbInitio(:)

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-ClH2+ (option 3: a,R+,R-)----'
  CALL Init_Model(QModel,pot_name='ClH2p',Print_init=.TRUE.,option=3)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='ClH2p_op34',test_var=test_var)

  CALL get_Q0_Model(Q,QModel,option=0)

  CALL Eval_pot_ON_Grid(QModel,Qmin=[0.8_Rkind,Q(2),Q(3)], &
                               Qmax=[2.6_Rkind,Q(2),Q(3)], &
                        nb_points=101,grid_file='grid_ClH2+_op34_Q1.tab')
  CALL Eval_pot_ON_Grid(QModel,Qmin=[Q(1),2.1_Rkind,Q(3)], &
                               Qmax=[Q(1),3.1_Rkind,Q(3)], &
                        nb_points=101,grid_file='grid_ClH2+_op34_Q2.tab')

  CALL Eval_pot_ON_Grid(QModel,Qmin=[Q(1),Q(2),-0.5_Rkind], &
                               Qmax=[Q(1),Q(2), 0.5_Rkind], &
                        nb_points=101,grid_file='grid_ClH2+_op34_Q3.tab')

  CALL Eval_pot_ON_Grid(QModel,Qmin=[0.8_Rkind,2.1_Rkind,Q(3)], &
                               Qmax=[2.6_Rkind,3.1_Rkind,Q(3)], &
                        nb_points=101,grid_file='grid_ClH2+_op34_Q12.tab')
  CALL Eval_pot_ON_Grid(QModel,Qmin=[Q(1),2.1_Rkind,-0.5_Rkind], &
                               Qmax=[Q(1),3.1_Rkind, 0.5_Rkind], &
                        nb_points=101,grid_file='grid_ClH2+_op34_Q23.tab')

  CALL Eval_pot_ON_Grid(QModel,Qmin=[0.8_Rkind,Q(2),-0.5_Rkind], &
                               Qmax=[2.6_Rkind,Q(2), 0.5_Rkind], &
                        nb_points=101,grid_file='grid_ClH2+_op34_Q13.tab')
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-ClH2+ (option 4: R1,R2,a)----'
  CALL Init_Model(QModel,pot_name='ClH2p',Print_init=.TRUE.,option=4)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='ClH2p_op34',test_var=test_var)


  CALL dealloc_dnMat(PotVal)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential test against the ab initio values '
  EAbInitio = [-461.057937359_Rkind,-461.027304042_Rkind,-460.988071490_Rkind]
  Qtest     = reshape([2.7_Rkind,2.5_Rkind,1.7_Rkind, &
                       2.3_Rkind,2.7_Rkind,1.1_Rkind, &
                       3.0_Rkind,2.5_Rkind,2.6_Rkind] &
              ,shape=[3,size(EAbInitio)])
  DO i=1,size(EAbInitio)
    Q = Qtest(:,i)
    write(out_unit,*) 'Q:'
    CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=0)
    CALL Write_dnMat(PotVal,nio=out_unit)
    write(out_unit,*) 'Difference with the ab initio value (cm-1)'

    PotVal = (PotVal - EAbInitio(i))*219475._Rkind
    CALL Write_dnMat(PotVal,nio=out_unit)
  END DO

  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_ClH2p_op34

SUBROUTINE test_ClH2p_op56
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  USE QML_ClH2p_m,      ONLY : QML_ClH2p_CCSDTF12, QML_ClH2p_Qsym_CCSDTF12, QML_ClH2p_Qsym_CCSDTF12_bis
  USE Opt_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: PotVal_gaussian
  real (kind=Rkind), allocatable :: qtest(:,:),EAbInitio(:)
  TYPE (QML_Opt_t)               :: Opt_p
  real (kind=Rkind)              :: V


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-ClH2+ (option 5: a,R+,R-)----'
  CALL Init_Model(QModel,pot_name='ClH2p',Print_init=.TRUE.,option=5)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='ClH2p_op56',test_var=test_var)


  CALL get_Q0_Model(Q,QModel,option=0)

  CALL Eval_pot_ON_Grid(QModel,Qmin=[0.8_Rkind,Q(2),Q(3)], &
                               Qmax=[2.6_Rkind,Q(2),Q(3)], &
                        nb_points=101,grid_file='grid_ClH2+_op56_Q1.tab')
  CALL Eval_pot_ON_Grid(QModel,Qmin=[Q(1),2.1_Rkind,Q(3)], &
                               Qmax=[Q(1),3.1_Rkind,Q(3)], &
                        nb_points=101,grid_file='grid_ClH2+_op56_Q2.tab')

  CALL Eval_pot_ON_Grid(QModel,Qmin=[Q(1),Q(2),-0.5_Rkind], &
                               Qmax=[Q(1),Q(2), 0.5_Rkind], &
                        nb_points=101,grid_file='grid_ClH2+_op56_Q3.tab')

  CALL Eval_pot_ON_Grid(QModel,Qmin=[0.8_Rkind,2.1_Rkind,Q(3)], &
                               Qmax=[2.6_Rkind,3.1_Rkind,Q(3)], &
                        nb_points=101,grid_file='grid_ClH2+_op56_Q12.tab')
  CALL Eval_pot_ON_Grid(QModel,Qmin=[Q(1),2.1_Rkind,-0.5_Rkind], &
                               Qmax=[Q(1),3.1_Rkind, 0.5_Rkind], &
                        nb_points=101,grid_file='grid_ClH2+_op56_Q23.tab')

  CALL Eval_pot_ON_Grid(QModel,Qmin=[0.8_Rkind,Q(2),-0.5_Rkind], &
                               Qmax=[2.6_Rkind,Q(2), 0.5_Rkind], &
                        nb_points=101,grid_file='grid_ClH2+_op56_Q13.tab')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK SUBROUTINE ----------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL dealloc_dnMat(PotVal)
  Q = [ 1.3_Rkind,1.2_Rkind, -0.5_Rkind]
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=0)
  !CALL Write_dnMat(PotVal,nio=out_unit)
  write(out_unit,*) 'pot (QML)',get_d0(PotVal)

  CALL QML_ClH2p_Qsym_CCSDTF12(V,Q)
  write(out_unit,*) 'pot (sub)',V

  Q = [1.7_Rkind, 0.7_Rkind, 1.3_Rkind]
  CALL QML_ClH2p_CCSDTF12(V,Q)
  write(out_unit,*) 'pot (sub)',V
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-ClH2+ (option 6: R1,R2,a)----'
  CALL Init_Model(QModel,pot_name='ClH2p',Print_init=.TRUE.,option=6)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives at the minimum'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='ClH2p_op56',test_var=test_var,last_test=.TRUE.)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Find the true minimum'
  CALL get_Q0_Model(Q,QModel,option=0)
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)

  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=7)
  CALL QML_Opt(Q,QModel,Opt_p,Q0=Q)

  deallocate(Q)
  CALL dealloc_dnMat(PotVal)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_ClH2p_op56

SUBROUTINE test_template
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: nderiv
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Template potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)

  CALL Init_Model(QModel,pot_name='template')

  allocate(q(QModel%QM%ndim))
  Q(:) = 2._Rkind

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
  write(out_unit,*) ' Potential and derivatives'

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Template', &
      test_var=test_var,last_test=.TRUE.)

  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'



END SUBROUTINE test_template
SUBROUTINE test_TwoD
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' TwoD potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)

  CALL Init_Model(QModel,pot_name='TwoD',adiabatic=.FALSE.,Print_init=.FALSE.)

  allocate(q(QModel%ndim))
  Q(:) = [3.875_Rkind,0.5_Rkind]

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
  write(out_unit,*) ' Potential and derivatives',nderiv

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='twod', &
      test_var=test_var,last_test=.TRUE.)


  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_TwoD
SUBROUTINE test_TwoD_RJDI2014
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' TwoD_RJDI2014 potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)
  CALL Init_Model(QModel,pot_name='TwoD_RJDI2014',adiabatic=.TRUE.,Print_init=.TRUE.,option=3)

  Q = [ZERO,HALF]

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
  write(out_unit,*) ' Potential and derivatives',nderiv

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='TwoD_RJDI2014', &
      test_var=test_var,last_test=.TRUE.)


  Q = [TWO,HALF]
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)


  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_TwoD_RJDI2014
SUBROUTINE test_TwoD_Valahu2022
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' TwoD_Valahu2022 potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)
  CALL Init_Model(QModel,pot_name='TwoD_Valahu2022',adiabatic=.TRUE.,Print_init=.TRUE.,option=1,PubliUnit=.TRUE.)

  allocate(Q(QModel%ndim))

  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  option = 1
  CALL get_Q0_Model(Q,QModel,option)
  write(out_unit,*) 'Q(:) (no unit):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives',nderiv

  option = 1
  CALL get_Q0_Model(Q,QModel,option)
  write(out_unit,*) 'Q(:) (no unit):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv) 
  write(out_unit,*) 'Energy (2*pi Hz)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='TwoD_Valahu2022', &
      test_var=test_var,last_test=.TRUE.)

  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'


END SUBROUTINE test_TwoD_Valahu2022
SUBROUTINE test_OneD_Photons
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' test_OneD_Photons potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)
  CALL Init_Model(QModel,pot_name='OneD_Photons',adiabatic=.FALSE.,Print_init=.TRUE.)

  Q = [THREE,ZERO]

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
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='OneD_Photons', &
      test_var=test_var,last_test=.TRUE.)


      ! from Eduarda:    -0.55055055    0.65065065    0.07056147   -0.00110611   -0.00110611    0.43460551
  Q = [-0.55055055_Rkind, 0.65065065_Rkind]
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  write(out_unit,*) 'Diabatic Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

  QModel%QM%adiabatic = .TRUE.
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  !        -0.5506        0.6507        0.0706
  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  write(out_unit,*) 'adiabatic Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)


  CALL dealloc_dnMat(PotVal)
  deallocate(Q)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 2D-Grid ------------------------'

  CALL Eval_pot_ON_Grid(QModel,Qmin=[-ONE,-2._Rkind],Qmax=[Ten,2._Rkind],         &
  nb_points=101,grid_file='grid_OneD_photons')
  write(out_unit,*) '---------------------------------------------'

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_OneD_Photons
SUBROUTINE test_HNO3
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnS_t),      allocatable :: Func(:)

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 9D-HNO3 ------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,ndim=9,pot_name='HNO3')

  Q = [2.65450_Rkind,2.28_Rkind,2.0_Rkind,Pi/TWO,                    &
          ZERO,ZERO,1.83492_Rkind,1.77157_Rkind,Pi/TWO]
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='HNO3', &
      test_var=test_var,last_test=.TRUE.)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Functions'
  write(out_unit,*) 'Q:'
  Q(1) = 1.57080_Rkind
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  CALL Eval_Func(QModel,Q,Func,nderiv=nderiv)

  DO i=1,size(Func)
    CALL Write_dnS(Func(i))
  END DO
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_HNO3
SUBROUTINE test_PH4Jo
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE
  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,i1,option,nio
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnS_t),      allocatable :: Func(:)
  TYPE (QML_Opt_t)               :: Opt_p

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 1,2,9D-PH4: H+PH3 -> H2+PH2 ----'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='PH4',option=7,Print_init=.FALSE.)
 
  allocate(Q(QModel%ndim))
  Q(:) = [1000._Rkind]
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=0)
  CALL Write_dnMat(PotVal,nio=out_unit)

  Q(:) = [-1000._Rkind]
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=0)
  CALL Write_dnMat(PotVal,nio=out_unit)
  deallocate(Q)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 2D-PH4: H+PH3 -> H2+PH2 --------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='PH4',option=8,Print_init=.FALSE.)

  allocate(Q(QModel%ndim))
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=4,list_act=[1])
  CALL QML_Opt(Q,QModel,Opt_p,Q0=[1.4_Rkind,1000._Rkind])

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)


  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=4,list_act=[1])
  CALL QML_Opt(Q,QModel,Opt_p,Q0=[2.6_Rkind,-1000._Rkind])

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=4,list_act=[1,2])
  Opt_p%nb_neg = 1
  CALL QML_Opt(Q,QModel,Opt_p,Q0=[1._Rkind,0._Rkind])
  
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---- 2D-Grid --------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL Eval_pot_ON_Grid(QModel,nb_points=100, &
                        Qmin=[0.4_Rkind,-20._Rkind],Qmax=[1.7_Rkind,+20._Rkind],  &
                        nderiv=0,grid_file='grid_PH4_2D')

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---- 1D-functions ---------------------------'
  write(out_unit,*) '---------------------------------------------'
  allocate(Func(QModel%QM%nb_Func))
  open(newunit=nio,file='grid_PH4_Func')
  DO i=-200,200
    Q = [i*0.1_Rkind]
    CALL Eval_Func(QModel,Q,Func,nderiv=0)
    write(nio,*) 'func',Q,get_d0(Func([1,2,10,18]))
  END DO
  close(nio)

  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_PH4jo
SUBROUTINE test_PH4
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,i1,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnS_t),      allocatable :: Func(:)

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 9D-PH4: H+PH3 -> H2+PH2 --------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,ndim=9,pot_name='PH4',option=4,Print_init=.FALSE.)

  allocate(Q(QModel%ndim))
  CALL get_Q0_Model(Q,QModel,option=3)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='PH4', &
      test_var=test_var,last_test=.TRUE.)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Functions'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  DO i=0,200
    Q(1)= -FIVE + real(i,kind=Rkind)*TEN/real(200,kind=Rkind)
    CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
    write(out_unit,*) Q(1:1),get_d0(Func)
    flush(6)
  END DO

  Q(1)= 0.5_Rkind
  CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
  write(out_unit,*) Q(1),'Energy',get_d0(Func(1))
  write(out_unit,*) Q(1),'Qopt',get_d0(Func(2:9))
  write(out_unit,*) Q(1),'Grad',get_d0(Func(10:17))
  write(out_unit,*) Q(1),'hessian'
  DO i=1,64
    write(out_unit,*) '                       ',get_d0(Func(17+i)),',     &'
  END DO

  deallocate(Q)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential in the asymptotic region'

  allocate(Q(QModel%ndim))
  Q(1)= 10._Rkind
  CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
  Q(2:9) = get_d0(Func(2:9))

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_PH4
SUBROUTINE test_CH5
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,i1,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnS_t),      allocatable :: Func(:)

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 12D-CH5: H+CH4 -> H2+CH3 -------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,ndim=12,pot_name='CH5',option=5,Print_init=.FALSE.)

  allocate(Q(QModel%ndim))
  CALL get_Q0_Model(Q,QModel,option=5)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv,test_var)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='CH5', &
      test_var=test_var,last_test=.TRUE.)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Functions'
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)

  DO i=0,200
    Q(1)= -FIVE + real(i,kind=Rkind)*TEN/real(200,kind=Rkind)
    CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
    write(out_unit,*) Q,get_d0(Func)
    flush(6)
  END DO

  Q(1)= 0.5_Rkind
  CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
  write(out_unit,*) Q(1),'Energy',get_d0(Func(1))
  write(out_unit,*) Q(1),'Qopt',get_d0(Func(2:12))
  write(out_unit,*) Q(1),'hessian'
  DO i=1,121
    write(out_unit,*) '                       ',get_d0(Func(12+i)),',     &'
  END DO

  deallocate(Q)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential in the asymptotic region'

  allocate(Q(QModel%ndim))
  Q(1)= 5._Rkind
  CALL Eval_Func(QModel,Q(1:1),Func,nderiv=0)
  Q(2:12) = get_d0(Func(2:12))

  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_CH5
SUBROUTINE test_HOO_DMBE
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:,:),x(:,:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-HOO_DMBE --------------------'

  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HOO_DMBE',Print_init=.FALSE.)

  allocate(Q(QModel%ndim,3))

  ! OH...O in TableVII of JCP paper, E=-0.1738 Hartree
  Q(:,1) = [5.663_Rkind,3.821_Rkind,1.842_Rkind]

  ! H...O-O in TableVII of JCP paper, E=-0.1916 Hartree
  Q(:,2) = [2.282_Rkind,7.547_Rkind,9.829_Rkind]

  ! HO2 in TableVII of JCP paper, E=-0.2141 Hartree
  Q(:,3) = [2.806_Rkind,2.271_Rkind,2.271_Rkind]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  DO i=1,size(Q,dim=2)
    write(out_unit,*) 'Q:'
    CALL Write_Vec(Q(:,i),out_unit,QModel%ndim)

    CALL Eval_Pot(QModel,Q(:,i),PotVal,nderiv=nderiv)
    CALL Write_dnMat(PotVal,nio=out_unit)

    ! For testing the model
    CALL Test_QdnV_FOR_Model(Q(:,i),PotVal,QModel,info='HOO_DMBE', &
       test_var=test_var,last_test=(i == size(Q,dim=2)))
  END DO
  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) ' END Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) ' With Cartesian coordinates'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='HOO_DMBE',Cart_TO_Q=.TRUE.,Print_init=.FALSE.)
  nderiv = 1

  allocate(x(3,3))
  x(:,1) = [ZERO,ZERO,-7.547_Rkind] ! H1
  x(:,2) = [ZERO,ZERO,ZERO]         ! O2
  x(:,3) = [ZERO,ZERO,2.282_Rkind]  ! O3

  CALL Eval_Pot(QModel,reshape(x,shape=[9]),PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  deallocate(x)
  write(out_unit,*) ' END Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_HOO_DMBE
SUBROUTINE test_H3
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : diagonalization, Write_Vec, Write_Mat
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:,:),x(:,:),Qcart(:),Q1D(:),Q2D(:)
  real (kind=Rkind), allocatable :: MW_hessian(:,:),EigVal(:),EigVec(:,:)
  integer                        :: ndim,nsurf,nderiv,i,j,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (QML_Opt_t)               :: Opt_p

  TYPE (dnS_t),    allocatable   :: Func(:)
  real(kind=Rkind), parameter    :: auTOcm_inv = 219474.631443_Rkind


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-H+H2 ------------------------'
  CALL Init_Model(QModel,pot_name='H3_LSTH',Print_init=.TRUE.)
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='H3_LSTH')

  allocate(Q(QModel%ndim,3))

  ! TableII of JCP paper, E=9.802 kcal/mol
  Q(:,1) = [1.757_Rkind,1.757_Rkind,1.757_Rkind]
  Q(:,2) = [1.4_Rkind,1000._Rkind,1001.4_Rkind]
  Q(:,3) = [1000._Rkind,1000._Rkind,2000._Rkind]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  DO i=1,size(Q,dim=2)
    write(out_unit,*) 'Q:'
    CALL Write_Vec(Q(:,i),out_unit,QModel%ndim)

    CALL Eval_Pot(QModel,Q(:,i),PotVal,nderiv=nderiv)
    CALL Write_dnMat(PotVal,nio=out_unit)

    ! For testing the model
    CALL Test_QdnV_FOR_Model(Q(:,i),PotVal,QModel,info='H3_LSTH', &
        test_var=test_var,last_test=.FALSE.)
  END DO

  write(out_unit,*) ' Potential and derivatives at the asymptotic chanel'
  Q(:,2) = [1.4010443631833713_Rkind,1000._Rkind,1001.4010443631833713_Rkind]
  write(out_unit,*) 'Q:'
  CALL Write_Vec(Q(:,2),out_unit,QModel%ndim)

  CALL Eval_Pot(QModel,Q(:,2),PotVal,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)

  write(out_unit,*) ' H2 Harmonic frequency. k,m,w', &
         PotVal%d2(1,1,1,1),QModel%QM%masses(1)/TWO, &
         auTOcm_inv*sqrt( PotVal%d2(1,1,1,1)/(QModel%QM%masses(1)/TWO) )

  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  write(out_unit,*) ' END Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Optimization H2 ..... H'

  CALL Init_Model(QModel,Print_init=.TRUE.,Cart_TO_Q=.TRUE.,pot_name='H3_LSTH')
  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=6,list_act=[9])
  CALL Write_QML_Opt(Opt_p)
  allocate(Qcart(QModel%ndim))

  CALL QML_Opt(Qcart,QModel,Opt_p,Q0=[ZERO,ZERO,-1000._Rkind,ZERO,ZERO,ZERO,ZERO,ZERO,1.4_Rkind])

  write(out_unit,*) 'Potential+derivatives:'
  CALL Eval_Pot(QModel,Qcart,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  !CALL Test_QdnV_FOR_Model(Qcart,PotVal,QModel,info='H3_LSTH', &
  !    test_var=test_var,last_test=.FALSE.)

  ! harmonic frequencies
  MW_hessian = reshape(get_d2(PotVal),shape=[QModel%ndim,QModel%ndim])
  allocate(EigVal(QModel%ndim))
  allocate(EigVec(QModel%ndim,QModel%ndim))

  DO i=1,QModel%ndim
  DO j=1,QModel%ndim
    MW_hessian(i,j) = MW_hessian(i,j) / sqrt( QModel%QM%masses((i-1)/3+1)*QModel%QM%masses((j-1)/3+1) )
  END DO
  END DO
  !CALL Write_Mat(MW_hessian, out_unit, 5, info='MWhessian')
  CALL diagonalization(MW_hessian,EigVal,EigVec,QModel%ndim)
  EigVal = sqrt(abs(EigVal))*auTOcm_inv
  CALL Write_Vec(EigVal, out_unit, 5, info='freq (cm-1):')

  deallocate(MW_hessian)
  deallocate(EigVal)
  deallocate(EigVec)
  ! enb harmonic frequencies

  CALL dealloc_dnMat(PotVal)
  deallocate(Qcart)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Optimization (TS) H-H-H'

  CALL Init_Model(QModel,Print_init=.TRUE.,Cart_TO_Q=.TRUE.,pot_name='H3_LSTH')
  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=4,list_act=[3,9])
  CALL Write_QML_Opt(Opt_p)
  allocate(Qcart(QModel%ndim))

  CALL QML_Opt(Qcart,QModel,Opt_p,Q0=[ZERO,ZERO,-2._Rkind,ZERO,ZERO,ZERO,ZERO,ZERO,2._Rkind])

  write(out_unit,*) 'Potential+derivatives:'
  CALL Eval_Pot(QModel,Qcart,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  !CALL Test_QdnV_FOR_Model(Qcart,PotVal,QModel,info='H3_LSTH', &
  !    test_var=test_var,last_test=.TRUE.)


 ! harmonic frequencies
  MW_hessian = reshape(get_d2(PotVal),shape=[QModel%ndim,QModel%ndim])
  allocate(EigVal(QModel%ndim))
  allocate(EigVec(QModel%ndim,QModel%ndim))

  DO i=1,QModel%ndim
  DO j=1,QModel%ndim
    MW_hessian(i,j) = MW_hessian(i,j) / sqrt( QModel%QM%masses((i-1)/3+1)*QModel%QM%masses((j-1)/3+1) )
  END DO
  END DO
  !CALL Write_Mat(MW_hessian, out_unit, 5, info='MWhessian')
  CALL diagonalization(MW_hessian,EigVal,EigVec,QModel%ndim)
  EigVal = sqrt(abs(EigVal))*auTOcm_inv
  CALL Write_Vec(EigVal, out_unit, 5, info='freq (cm-1):')

  deallocate(MW_hessian)
  deallocate(EigVal)
  deallocate(EigVec)
  ! enb harmonic frequencies

  CALL dealloc_dnMat(PotVal)
  deallocate(Qcart)
  CALL dealloc_Model(QModel)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Optimization H2-H (weak minimum in the asymtotic region)'

  CALL Init_Model(QModel,Print_init=.TRUE.,Cart_TO_Q=.TRUE.,pot_name='H3_LSTH')
  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=6,list_act=[3,9])
  CALL Write_QML_Opt(Opt_p)
  allocate(Qcart(QModel%ndim))

  CALL QML_Opt(Qcart,QModel,Opt_p,Q0=[ZERO,ZERO,-4._Rkind,ZERO,ZERO,ZERO,ZERO,ZERO,1.4_Rkind])

  write(out_unit,*) 'Potential+derivatives:'
  CALL Eval_Pot(QModel,Qcart,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  !CALL Test_QdnV_FOR_Model(Qcart,PotVal,QModel,info='H3_LSTH', &
  !    test_var=test_var,last_test=.FALSE.)

  CALL dealloc_dnMat(PotVal)
  deallocate(Qcart)
  CALL dealloc_Model(QModel)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' MEP (op=61) : (1/R1, 1/R2) => (rho, s)'
  CALL Init_Model(QModel,pot_name='H3_LSTH',option=61)

  Q1D = [ZERO]
  write(out_unit,*) 'Potential+derivatives at:',Q1D
  CALL Eval_Pot(QModel,Q1D,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL Eval_pot_ON_Grid(QModel,Qmin=[-FOUR],Qmax=[FOUR],nb_points=801,          &
                        grid_file='grid_H3_V61')

  allocate(Func(QModel%QM%nb_Func))

  DO i=-400,400
    Q1D = real(i,kind=Rkind)/100
    CALL Eval_Func(QModel,Q1D,Func,nderiv=0)
    write(out_unit,*) 'grid_H3_func61',Q1D,get_d0(Func)
    flush(6)
  END DO

  deallocate(Func)
  CALL dealloc_dnMat(PotVal)
  deallocate(Q1D)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' MEP (op=72) : (1/R1, 1/R2) => (rho, s)'
  CALL Init_Model(QModel,pot_name='H3_LSTH',option=72)

  Q2D = [ZERO,1.242_Rkind]
  write(out_unit,*) 'Potential+derivatives at:',Q2D
  CALL Eval_Pot(QModel,Q2D,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL Eval_pot_ON_Grid(QModel,Qmin=[-FOUR,0.8_Rkind],Qmax=[FOUR,2._Rkind],nb_points=100,          &
                        grid_file='grid_H3_V72')

  CALL dealloc_dnMat(PotVal)
  deallocate(Q2D)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' MEP (op=74) : (1/R1, 1/R2) => (rho, s)'
  CALL Init_Model(QModel,pot_name='H3_LSTH',option=74)

  Q2D = [ZERO,1.242_Rkind]
  write(out_unit,*) 'Potential+derivatives at:',Q2D
  CALL Eval_Pot(QModel,Q2D,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL Eval_pot_ON_Grid(QModel,Qmin=[-FOUR,0.8_Rkind],Qmax=[FOUR,2._Rkind],nb_points=100,          &
                        grid_file='grid_H3_V74')

  CALL dealloc_dnMat(PotVal)
  deallocate(Q2D)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_H3
SUBROUTINE test_IRC_H3
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  USE IRC_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  TYPE (QML_IRC_t)               :: IRC_p

  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' IRC with H3 potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)
  CALL Init_Model(QModel,Print_init=.TRUE.,Cart_TO_Q=.TRUE.,                    &
                  read_param=.TRUE.,param_file_name='DAT_files/irc_h3.dat')

  CALL Init_QML_IRC(IRC_p,QModel,read_param=.TRUE.,param_file_name='DAT_files/irc_h3.dat')

  allocate(Q(QModel%QM%ndim))


  CALL QML_IRC(Q,QModel,IRC_p,Q0=[ZERO,ZERO,-TWO,ZERO,ZERO,ZERO,ZERO,ZERO,TWO])


  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_IRC_H3

SUBROUTINE test_IRC_H3_AbInitio
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  USE IRC_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  TYPE (QML_IRC_t)               :: IRC_p

  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' IRC with H3 ab initio'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)
  CALL Init_Model(QModel,Print_init=.TRUE.,Cart_TO_Q=.TRUE.,                    &
                  read_param=.TRUE.,param_file_name='DAT_files/irc_h3_abinitio.dat')

  CALL Init_QML_IRC(IRC_p,QModel,read_param=.TRUE.,param_file_name='DAT_files/irc_h3_abinitio.dat')

  allocate(Q(QModel%QM%ndim))


  CALL QML_IRC(Q,QModel,IRC_p,Q0=[ZERO,ZERO,-TWO,ZERO,ZERO,ZERO,ZERO,ZERO,TWO])


  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_IRC_H3_AbInitio

SUBROUTINE test_CNH_Murrell
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  USE IRC_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:,:),Qopt(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (QML_Opt_t)               :: Opt_p


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-CNH ------------------------'
  CALL Init_Model(QModel,pot_name='CNH_Murrell',option=1) ! jacobi
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(Q(QModel%ndim,3))


  Q(:,1) = [ZERO,         3.18722_Rkind,2.17926_Rkind]
  Q(:,2) = [Pi,           2.89321_Rkind,2.20061_Rkind]
  Q(:,3) = [1.16809_Rkind,2.28376_Rkind,2.15316_Rkind]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  DO i=1,size(Q,dim=2)
    write(out_unit,*) 'Q:'
    CALL Write_Vec(Q(:,i),out_unit,QModel%ndim)

    CALL Eval_Pot(QModel,Q(:,i),PotVal,nderiv=nderiv)
    CALL Write_dnMat(PotVal,nio=out_unit)

    ! For testing the model
    CALL Test_QdnV_FOR_Model(Q(:,i),PotVal,QModel,info='CNH_Murrell',test_var=test_var)
  END DO
  CALL dealloc_dnMat(PotVal)

  write(out_unit,*) ' END Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  allocate(Qopt(QModel%ndim))
  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=10,list_act=[2,3])
  CALL Write_QML_Opt(Opt_p)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Optimization H-CN'
  CALL QML_Opt(Qopt,QModel,Opt_p,Q0=Q(:,1))

  write(out_unit,*) 'Potential+derivatives:'
  CALL Eval_Pot(QModel,Qopt,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL Test_QdnV_FOR_Model(Qopt,PotVal,QModel,info='CNH_Murrell',test_var=test_var)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Optimization CN-H'
  flush(out_unit)

  CALL QML_Opt(Qopt,QModel,Opt_p,Q0=Q(:,2))

  write(out_unit,*) 'Potential+derivatives:'
  CALL Eval_Pot(QModel,Qopt,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)
  flush(out_unit)

  CALL Test_QdnV_FOR_Model(Qopt,PotVal,QModel,info='CNH_Murrell',test_var=test_var)


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Optimization CNH (TS)'
  Opt_p%nb_neg = 1
  CALL Init_QML_Opt(Opt_p,QModel,read_param=.FALSE.,icv=10,list_act=[1,2,3])
  CALL QML_Opt(Qopt,QModel,Opt_p,Q0=Q(:,3))

  write(out_unit,*) 'Potential+derivatives:'
  CALL Eval_Pot(QModel,Qopt,PotVal,nderiv=2)
  CALL Write_dnMat(PotVal,nio=out_unit)

  CALL Test_QdnV_FOR_Model(Qopt,PotVal,QModel,info='CNH_Murrell',test_var=test_var)



  CALL dealloc_dnMat(PotVal)
  deallocate(Qopt)
  CALL dealloc_Model(QModel)
  deallocate(Q)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 3D-CNH ------------------------'
  CALL Init_Model(QModel,pot_name='CNH_Murrell',option=11) ! jacobi
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(Q(QModel%ndim,3))

  Q(:,1) = [ ONE,              3.18722_Rkind,2.17926_Rkind]
  Q(:,2) = [-ONE,              2.89321_Rkind,2.20061_Rkind]
  Q(:,3) = [cos(1.16809_Rkind),2.28376_Rkind,2.15316_Rkind]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  DO i=1,size(Q,dim=2)
    write(out_unit,*) 'Q:'
    CALL Write_Vec(Q(:,i),out_unit,QModel%ndim)

    CALL Eval_Pot(QModel,Q(:,i),PotVal,nderiv=nderiv)
    CALL Write_dnMat(PotVal,nio=out_unit)

    ! For testing the model
    CALL Test_QdnV_FOR_Model(Q(:,i),PotVal,QModel,info='CNH_Murrell',test_var=test_var,last_test=(i==size(Q,dim=2)))
  END DO
  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) ' END Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'


  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 1D-CNH ------------------------'
  IF (allocated(QModel%QM)) deallocate(QModel%QM)
  CALL Init_Model(QModel,pot_name='CNH_Murrell',option=2) ! jacobi
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) 'IndexFunc_Ene:             ',QModel%QM%IndexFunc_Ene
  write(out_unit,*) 'IndexFunc_Qop:             ',QModel%QM%IndexFunc_Qop
  write(out_unit,*) 'IndexFunc_Hess:            ',QModel%QM%IndexFunc_Hess
  write(out_unit,*) '---------------------------------------------'


  allocate(Q(QModel%ndim,3))

  Q(:,1) = [ZERO]
  Q(:,2) = [Pi]
  Q(:,3) = [1.16809_Rkind]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  DO i=1,size(Q,dim=2)
    write(out_unit,*) 'Q:'
    CALL Write_Vec(Q(:,i),out_unit,QModel%ndim)

    CALL Eval_Pot(QModel,Q(:,i),PotVal,nderiv=nderiv)
    CALL Write_dnMat(PotVal,nio=out_unit)

    ! For testing the model
    CALL Test_QdnV_FOR_Model(Q(:,i),PotVal,QModel,info='CNH_Murrell', &
        test_var=test_var,last_test=.FALSE.)

  END DO
  CALL dealloc_dnMat(PotVal)
  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) ' END Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  nderiv = 2
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '------------ 1D-CNH ------------------------'
  CALL Init_Model(QModel,pot_name='CNH_Murrell',option=21) ! jacobi
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  allocate(Q(QModel%ndim,3))

  Q(:,1) = [ONE]
  Q(:,2) = [-ONE]
  Q(:,3) = [cos(1.16809_Rkind)]

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  DO i=1,size(Q,dim=2)
    write(out_unit,*) 'Q:'
    CALL Write_Vec(Q(:,i),out_unit,QModel%ndim)

    CALL Eval_Pot(QModel,Q(:,i),PotVal,nderiv=nderiv)
    CALL Write_dnMat(PotVal,nio=out_unit)

    ! For testing the model
    CALL Test_QdnV_FOR_Model(Q(:,i),PotVal,QModel,info='CNH_Murrell', &
        test_var=test_var,last_test=(i == size(Q,dim=2)))

  END DO
  CALL dealloc_dnMat(PotVal)
  CALL dealloc_Model(QModel)
  deallocate(Q)
  write(out_unit,*) ' END Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_CNH_Murrell

SUBROUTINE test_Opt_MullerBrown
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : TO_string
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  TYPE (QML_Opt_t)               :: Opt_param

  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Optimization with 2D Müller-Brown potential'
  write(out_unit,*) ' Assumed units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)

  DO option=1,5

    CALL Init_Model(QModel,pot_name='2D_MB',PubliUnit=.TRUE.,Print_init=.FALSE.,option=option)

    Q = QModel%QM%Q0

    CALL Init_QML_Opt(Opt_param,QModel,read_param=.FALSE.)
    IF (option > 3) Opt_param%nb_neg = 1
    CALL QML_Opt(Q,QModel,Opt_param,Q0=QModel%QM%Q0)

    write(out_unit,*) 'Potential+derivatives:'
    CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)
    CALL Write_dnMat(PotVal,nio=out_unit)
    flush(out_unit)

    ! For testing the model
    CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='Müller-Brown' // TO_string(option), &
        test_var=test_var,last_test=(option == 5))

    flush(out_unit)
    CALL dealloc_Model(QModel)

  END DO

  deallocate(Q)
  CALL dealloc_dnMat(PotVal)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_Opt_MullerBrown
SUBROUTINE test_IRC_MullerBrown
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m
  USE Opt_m
  USE IRC_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  TYPE (QML_IRC_t)               :: IRC_p

  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' IRC with the 2D-MullerBrown potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)
  CALL Init_Model(QModel,pot_name='2D_MB',option=4) ! the first TS


  CALL Init_QML_IRC(IRC_p,QModel,read_param=.TRUE.,param_file_name='DAT_files/irc_mb.dat')

  allocate(Q(QModel%QM%ndim))


  CALL QML_IRC(Q,QModel,IRC_p,Q0=QModel%QM%Q0)


  deallocate(Q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_IRC_MullerBrown
SUBROUTINE test_Test
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Test potential'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)

  CALL Init_Model(QModel,pot_name='Test',adiabatic=.TRUE.)

  allocate(q(QModel%ndim))
  Q(:) = ZERO

  nderiv=3

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '----- CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Check analytical derivatives with respect to numerical ones'

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%QM%ndim)
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_Test
SUBROUTINE test_Vibronic
  USE QDUtil_NumParameters_m
  USE QDUtil_m,         ONLY : Write_Vec
  USE ADdnSVM_m
  USE Model_m
  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: Q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings


  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Vibronic potential (TwoD_RJDI2014)'
  write(out_unit,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  flush(out_unit)

  CALL Init_Model(QModel,pot_name='Vibronic  TwoD_RJDI2014',nsurf=2,ndim=1,adiabatic=.TRUE.)


  Q = [ZERO,HALF]

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
  write(out_unit,*) ' Potential and derivatives',nderiv

  CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv)

  write(out_unit,*) 'Q(:) (bohr):'
  CALL Write_Vec(Q,out_unit,QModel%ndim)
  write(out_unit,*) 'Energy (Hartree)'
  CALL Write_dnMat(PotVal,nio=out_unit)

    ! For testing the model
  CALL Test_QdnV_FOR_Model(Q,PotVal,QModel,info='vibronic', &
      test_var=test_var,last_test=.TRUE.)

  CALL dealloc_Model(QModel)

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '- END CHECK POT -----------------------------'
  write(out_unit,*) '---------------------------------------------'

END SUBROUTINE test_Vibronic

SUBROUTINE test_NO3_test
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m
  USE Model_m

  IMPLICIT NONE

  TYPE (Model_t)                 :: QModel
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: nderiv
  TYPE (dnMat_t)                 :: PotVal
  TYPE (dnMat_t)                 :: NAC ! for non adiabatic couplings

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' NO3 potential'
  write(out_unit,*) ' With units: Atomic Units  ??? '
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'

  CALL Init_Model(QModel,pot_name='NO3',PubliUnit=.FALSE.)
STOP
  allocate(q(QModel%QM%ndim))

  q(:) = ZERO
  nderiv=2

  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) '---------------------------------------------'
  write(out_unit,*) ' Potential and derivatives'
  write(out_unit,*) '---------------------------------------------'

  QModel%QM%adiabatic = .TRUE.
  write(out_unit,*) 'ADIABATIC potential'
  write(out_unit,*) 'Evaluated in', q

  CALL Eval_Pot(QModel,Q,PotVal,NAC=NAC,nderiv=nderiv)
  CALL Write_dnMat(PotVal,nio=out_unit)


  CALL dealloc_dnMat(PotVal)
  deallocate(q)
  CALL dealloc_Model(QModel)

END SUBROUTINE test_NO3_test
END PROGRAM TEST_model