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
!    Copyright 2017 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================
PROGRAM TEST_model
  IMPLICIT NONE

  CALL test_LinearHBond()
  CALL test_HenonHeiles()
  CALL test_Phenol()
  CALL test_Buckingham()
  CALL test_Morse()
  CALL test_Tully()
  CALL test_1DSOC()
  CALL test_PSB3()

END PROGRAM TEST_model

SUBROUTINE test_Tully
  USE mod_dnMatPot
  USE mod_Model
  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMatPot)                :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Tully potential'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'

  DO option=1,3
    CALL Init_Model(Para_Model,pot_name='Tully',option=option)
    CALL Write_Model(Para_Model)
    allocate(q(Para_Model%ndim))
    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,*) '----- CHECK POT -----------------------------'
    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,'(a,i1)') ' ----- TULLY',option

    write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

    q(:) = 1.1_Rkind
    write(out_unitp,'(a,f12.6)') ' R (Bohr)',q(:)
    CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

    q(:) = -1.1_Rkind
    write(out_unitp,'(a,f12.6)') ' R (Bohr)',q(:)
    CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)


    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,*) ' Potential and derivatives'

    CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)
    write(out_unitp,'(a,f12.6)') 'R (Bohr)',q(:)
    write(out_unitp,*) 'Energy (Hartree)'
    CALL Write_dnMatPot(PotVal,nio=out_unitp)

    write(out_unitp,*) '---------------------------------------------'
    write(out_unitp,*) '- END CHECK POT -----------------------------'
    write(out_unitp,*) '---------------------------------------------'

    CALL dealloc_dnMatPot(PotVal)
    deallocate(q)

  END DO



  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '-----------TULLY 1 on a grid ----------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Tully1"'
  write(out_unitp,*) ' You should get the curves (columns: 2,3 5) of the panel (a) of figure 3 of ... '
  write(out_unitp,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unitp,*) ' Be carrefull, the non-adiabatic coupling is scaled by -1/50.'
  CALL Init_Model(Para_Model,pot_name='Tully',option=1)
  CALL Eval_pot_ON_Grid(Para_Model,Qmin=(/-TEN/),Qmax=(/TEN/),nb_points=1001, &
                        grid_file='grid_Tully1')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '-----------TULLY 2 on a grid ----------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Tully2"'
  write(out_unitp,*) ' You should get the curves (columns: 2,3 5) of the panel (b) of figure 3 of ... '
  write(out_unitp,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unitp,*) ' Be carrefull, the non-adiabatic coupling is scaled by -1/12.'

  CALL Init_Model(Para_Model,pot_name='Tully',option=2)
  CALL Eval_pot_ON_Grid(Para_Model,Qmin=(/-TEN/),Qmax=(/TEN/),nb_points=1001, &
                        grid_file='grid_Tully2')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '-----------TULLY 3 on a grid ----------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Tully3"'
  write(out_unitp,*) ' You should get the curves (columns: 2,3 5) of the panel (c) of figure 3 of ... '
  write(out_unitp,*) '   Tully, J. Chem. Phys. V93, pp15, 1990.'
  write(out_unitp,*) ' Be carrefull, the sign of non-adiabatic coupling is changed.'

  CALL Init_Model(Para_Model,pot_name='Tully',option=3)
  CALL Eval_pot_ON_Grid(Para_Model,Qmin=(/-TEN/),Qmax=(/TEN/),nb_points=1001, &
                        grid_file='grid_Tully3')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Tully


SUBROUTINE test_1DSOC
  USE mod_Lib
  USE mod_dnMatPot
  USE mod_Model
  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMatPot)                :: PotVal
  TYPE (dnMatPot)                :: Vec ! for non adiabatic couplings

  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' 1D-SOC potential'
  write(out_unitp,*) ' With units: Bohr Hartree (au)'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(Para_Model,pot_name='1DSOC')
  CALL Write_Model(Para_Model)

  allocate(q(Para_Model%ndim))
  q(:) = (/ 9.5_Rkind /) ! one evaluation to get vec0%d0(:,:)
  write(out_unitp,*) ' ref Eigenvectors at Q:',Q
  CALL Eval_Pot(Para_Model,Q,PotVal,Vec=Vec,nderiv=0)
  ! ref Eigenvectors
  CALL Write_dnMatPot(Para_Model%vec0,nio=out_unitp)



  q(:) = (/ 10._Rkind /)
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Q (Bohr):'
  CALL Write_RVec(Q,out_unitp,Para_Model%ndim)

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  CALL Eval_Pot(Para_Model,Q,PotVal,Vec=Vec,nderiv=nderiv)
  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  write(out_unitp,*) 'Non adiatic couplings:'
  CALL Write_dnMatPot(Vec,nio=out_unitp)

  write(out_unitp,*) 'DIABATIC potential'
  Para_Model%adiabatic = .FALSE.
  write(out_unitp,*) 'Q (Bohr):'
  CALL Write_RVec(Q,out_unitp,Para_Model%ndim)

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  CALL dealloc_dnMatPot(PotVal)
  CALL dealloc_dnMatPot(Vec)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  flush(out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----------- 1D-SOC on a grid ----------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of R in Bohr)'
  write(out_unitp,*) '   file name: "grid_1DSOC"'
  write(out_unitp,*) ' You should get the 1 curves of figure 2 of ... '
  write(out_unitp,*) '  ... Granucci et al. J. Chem. Phys. V137, p22A501 (2012)'

  CALL Init_Model(Para_Model,pot_name='1DSOC',option=1)
  Para_Model%adiabatic = .TRUE.
  flush(out_unitp)

  q(:) = (/ 9.5_Rkind /) ! one evaluation to get vec0%d0(:,:)
  CALL Eval_Pot(Para_Model,Q,PotVal,Vec=Vec,nderiv=1)

  ! ref Eigenvectors
  write(out_unitp,*) ' ref Eigenvectors at Q:',Q
  CALL Write_dnMatPot(Para_Model%vec0,nio=out_unitp)


  CALL Eval_pot_ON_Grid(Para_Model,Qmin=(/3._Rkind/),Qmax=(/20._Rkind/), &
                        nb_points=1001,nderiv=0,grid_file='grid_1DSOC')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  deallocate(q)
  CALL dealloc_dnMatPot(PotVal)
  CALL dealloc_dnMatPot(Vec)
END SUBROUTINE test_1DSOC

SUBROUTINE test_Morse
  USE mod_dnMatPot
  USE mod_Model
  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMatPot)                :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Morse potential (H-F parameters)'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(Para_Model,pot_name='Morse',read_param=.FALSE.)
  CALL Write_Model(Para_Model)

  allocate(q(Para_Model%ndim))
  q(:) = 2._Rkind

  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)
  write(out_unitp,'(a,f12.6)') 'R (Bohr)',q(:)
  write(out_unitp,*) 'Energy (Hartree)'
  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL dealloc_dnMatPot(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of R)'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Morse"'

  CALL Eval_pot_ON_Grid(Para_Model,Qmin=(/1._Rkind/),Qmax=(/5._Rkind/),nb_points=1001,&
                        grid_file='grid_Morse')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Morse

SUBROUTINE test_Buckingham
  USE mod_dnMatPot
  USE mod_Model
  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMatPot)                :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Buckingham potential (Ar-Ar parameters)'
  write(out_unitp,*) ' With units: Bohr and Hartree (atomic units)'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(Para_Model,pot_name='Buck',read_param=.FALSE.)
  CALL Write_Model(Para_Model)

  allocate(q(Para_Model%ndim))
  q(:) = 7._Rkind

  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)
  write(out_unitp,'(a,f12.6)') 'R (Bohr)',q(:)
  write(out_unitp,*) 'Energy (Hartree)'
  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL dealloc_dnMatPot(PotVal)
  deallocate(q)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential on a 1D grid (as a function of R)'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '   file name: "grid_Buck"'

  CALL Eval_pot_ON_Grid(Para_Model,Qmin=(/6._Rkind/),Qmax=(/20._Rkind/),nb_points=1001,&
                        grid_file='grid_Buck')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Buckingham

SUBROUTINE test_Phenol
  USE mod_Lib
  USE mod_dnMatPot
  USE mod_Model
  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMatPot)                :: PotVal
  TYPE (dnMatPot)                :: Vec ! for non adiabatic couplings


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Phenol potential'
  write(out_unitp,*) ' With units: Angs and eV'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(Para_Model,pot_name='Phenol',PubliUnit=.TRUE.)
  CALL Write_Model(Para_Model)

  allocate(q(Para_Model%ndim))

  q(:) = (/ 1.2_Rkind,0.2_Rkind /)
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) 'ADIABATIC potential'
  write(out_unitp,*) 'Q (Angs, radian):'
  CALL Write_RVec(Q,out_unitp,Para_Model%ndim)

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  CALL Eval_Pot(Para_Model,Q,PotVal,Vec=Vec,nderiv=nderiv)
  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  write(out_unitp,*) 'Non adiatic couplings:'
  CALL Write_dnMatPot(Vec,nio=out_unitp)

  write(out_unitp,*) 'DIABATIC potential'
  Para_Model%adiabatic = .FALSE.
  write(out_unitp,*) 'Q (Angs, radian):'
  CALL Write_RVec(Q,out_unitp,Para_Model%ndim)

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL dealloc_dnMatPot(PotVal)
  CALL dealloc_dnMatPot(Vec)
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

  CALL Init_Model(Para_Model,pot_name='phenol',PubliUnit=.TRUE.)
  Para_Model%adiabatic = .FALSE.

  CALL Eval_pot_ON_Grid(Para_Model,Qmin=(/0.5_Rkind,ZERO/),Qmax=(/5._Rkind,ZERO/), &
                        nb_points=1001,grid_file='grid_Phenol')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_Phenol
SUBROUTINE test_HenonHeiles
  USE mod_Lib
  USE mod_dnMatPot
  USE mod_Model
  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMatPot)                :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '------------ HenonHeiles --------------------'

  CALL Init_Model(Para_Model,pot_name='HenonHeiles',ndim=4)
  CALL Write_Model(Para_Model)

  allocate(q(Para_Model%ndim))
  q(:) = (/ (0.1_Rkind*real(i,kind=Rkind),i=1,Para_Model%ndim) /)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) 'Q:'
  CALL Write_RVec(Q,out_unitp,Para_Model%ndim)

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)
  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  CALL dealloc_dnMatPot(PotVal)
  deallocate(q)
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_HenonHeiles
SUBROUTINE test_LinearHBond
  USE mod_dnMatPot
  USE mod_Model
  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMatPot)                :: PotVal


  nderiv = 2
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) ' Linear H-Bond (symmetric one and the default parameters)'
  write(out_unitp,*) ' With units: Angs and kcal.mol^-1'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(Para_Model,pot_name='HBond',PubliUnit=.TRUE.)
  CALL Write_Model(Para_Model)

  allocate(q(Para_Model%ndim))

  q(:) = (/ 2.75_Rkind,0._Rkind /)
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)
  write(out_unitp,'(a,2f12.6)') 'QQ,q (Angs)',q(:)
  write(out_unitp,*) 'Energy (kcal.mol^-1)'
  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '- END CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  CALL dealloc_dnMatPot(PotVal)

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

  CALL Eval_pot_ON_Grid(Para_Model, &
                        Qmin=(/2.75_Rkind,-0.6_Rkind/),Qmax=(/2.75_Rkind,0.6_Rkind/),nb_points=1001,&
                        grid_file='grid_Hbond-sym')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Linear H-Bond (asymmetric one)'
  write(out_unitp,*) ' With units: Angs and kcal.mol^-1'
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(Para_Model,pot_name='HBond',PubliUnit=.TRUE.,             &
                read_param=.TRUE.,param_file_name='Hbond.dat')
  CALL Write_Model(Para_Model)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '----- CHECK POT -----------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Check analytical derivatives with respect to numerical ones'

  allocate(q(Para_Model%ndim))

  q(:) = (/ 3._Rkind,0._Rkind /)
  nderiv=2
  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' Potential and derivatives'

  CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)
  write(out_unitp,'(a,2f12.6)') 'QQ,q (Angs)',q(:)
  write(out_unitp,*) 'Energy (kcal.mol^-1)'
  CALL Write_dnMatPot(PotVal,nio=out_unitp)
  CALL dealloc_dnMatPot(PotVal)

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

  CALL Eval_pot_ON_Grid(Para_Model, &
                        Qmin=(/3._Rkind,-0.8_Rkind/),Qmax=(/3._Rkind,0.7_Rkind/),nb_points=1001,&
                        grid_file='grid_Hbond-asym')

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------------------------------------------'

END SUBROUTINE test_LinearHBond
SUBROUTINE test_PSB3

  USE mod_Lib
  USE mod_dnMatPot
  USE mod_Model

  IMPLICIT NONE

  TYPE (Param_Model)             :: Para_Model
  real (kind=Rkind), allocatable :: q(:)
  integer                        :: ndim,nsurf,nderiv,i,option
  TYPE (dnMatPot)                :: PotVal

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) ' PSB3 potential'
  IF(Para_Model%PubliUnit) THEN
     write(out_unitp,*) ' With units: Atomic Units (Angstrom, Rad, Rad, Kcal/mol)'
  END IF
  IF(.NOT. Para_Model%PubliUnit) THEN
     write(out_unitp,*) ' With units: Atomic Units (Bhor, Rad, Rad, Hartree)'
  END IF
  write(out_unitp,*) '---------------------------------------------'

  CALL Init_Model(Para_Model,pot_name='psb3',PubliUnit=.FALSE.)

  CALL Write_Model(Para_Model)

  allocate(q(Para_Model%ndim))

  q(:) = (/0.172459_Rkind,-3.14_Rkind,0._Rkind/)
  nderiv=2

  write(out_unitp,*) '---------------------------------------------'
  write(out_unitp,*) '---------- CHECK POT ------------------------'
  write(out_unitp,*) '---------------------------------------------'

  write(out_unitp,*) 'Evaluated in', q

  CALL Eval_Pot(Para_Model,Q,PotVal,nderiv=nderiv)

  CALL Write_dnMatPot(PotVal,nio=out_unitp)

  CALL Check_analytical_numerical_derivatives(Para_Model,Q,nderiv)

  write(out_unitp,*) '---------- END CHECK POT --------------------'
  write(out_unitp,*) '---------------------------------------------'

  CALL dealloc_dnMatPot(PotVal)

  deallocate(q)

END SUBROUTINE test_PSB3
