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
PROGRAM main_pot
  IMPLICIT NONE


  CALL test_PH4()
  CALL test_HBond()
  CALL test_1DSOC_1S1T()
  CALL test_PSB3()
  CALL test_Phenol_Dia(10**7)
  CALL test_Phenol_ADia()
  CALL test_henonheiles(10**7)
  !CALL test_Vib_adia(1000)
  CALL test_Test()
  !CALL test_Read_Model() ; stop

END PROGRAM main_pot
SUBROUTINE test_Read_Model()
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT, Rkind => real64
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option
  real (kind=Rkind)                  :: DQ

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '  Read model'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

  pot_name  = 'read_model'
  ndim      = 0 ! it would be initialized
  nsurf     = 0 ! it would be initialized
  option    = -1
  adiabatic = .TRUE.
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

  write(out_unit,*) 'ndim,nsurf',ndim,nsurf
  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))

  IF (nsurf > 1) THEN
    write(out_unit,*) ' Test V, G and NAC'

    allocate(NAC(nsurf,nsurf,ndim))

    CALL get_Qmodel_Q0(Q,0)
    DQ = 0.1_Rkind
    DO k=1,4
      write(out_unit,*) '------------------------------------------------------'
      CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)
      write(out_unit,*) 'Q',Q
      DO i=1,nsurf
        write(out_unit,*) 'Vadia',i,':',V(i,i)
      END DO
      DO i=1,nsurf
        write(out_unit,*) 'Gadia',i,':',G(i,i,1:ndim)
      END DO
      DO i=1,nsurf
        DO j=1,nsurf
          write(out_unit,*) 'NAC',i,j,':',NAC(j,i,1:ndim)
        END DO
      END DO

      Q(:) = Q(:) + DQ
    END DO

    deallocate(NAC)
  ELSE
    write(out_unit,*) ' Test V, G'

    CALL get_Qmodel_Q0(Q,0)

    CALL sub_Qmodel_VG(V,G,Q)
    write(out_unit,*) 'Q',Q
    write(out_unit,*) 'Potential',V
    write(out_unit,*) 'Gradient',G
  END IF

  deallocate(Q)
  deallocate(V)
  deallocate(G)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

END SUBROUTINE test_Read_Model
SUBROUTINE test2_Vib_adia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  maxth = 1
  write(out_unit,*) '============================================================'
  write(out_unit,*) '  Vibrational adiabatic separation'
  write(out_unit,*) 'TEST_driver. number of threads:',maxth

  pot_name  = 'read_model'
  ndim      = 0 ! it would be initialized
  nsurf     = 0 ! it would be initialized
  option    = -1
  adiabatic = .FALSE.
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

  write(out_unit,*) 'ndim,nsurf',ndim,nsurf

  write(out_unit,*) ' Test V, G and NAC'
  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  read(5,*) Q

  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)
  write(out_unit,*) 'Q',Q
  write(out_unit,*) 'Vadia',(V(i,i),i=1,nsurf)
  write(out_unit,*) 'Gadia',(G(i,i,1:ndim),i=1,nsurf)
  write(out_unit,*) 'NAC'
  DO i=1,nsurf
    write(out_unit,*) i,NAC(:,i,1:ndim)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(G)
  deallocate(NAC)

  write(out_unit,*) ' Test of evaluations of the potential ',pot_name
  CALL QML_time_perso('Test ' // pot_name)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

  Q = 2.5d0
  DO i=1,nb_eval
    Q = Q + 0.05d0

    CALL sub_Qmodel_V(V,Q)

    write(out_unit,*) 'pot:',Q,(V(k,k),k=1,nsurf)
  END DO

  deallocate(Q)
  deallocate(V)

  CALL QML_time_perso('Test ' // pot_name)
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

END SUBROUTINE test2_Vib_adia
SUBROUTINE test_1DSOC_1S1T
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: GGdef(:,:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: Vec(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: h(:,:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  nb_eval = 1

  pot_name = '1DSOC_1S1T'
  write(out_unit,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name

  ndim      = 1
  nsurf     = 4
  option    = 0
  adiabatic = .TRUE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(Vec(nsurf,nsurf))


  Q(:) = (/8.5_Rkind /)

  CALL sub_Qmodel_VVec(V,Vec,Q)
  write(out_unit,*) ' Diabatic potential as a 4x4 matrix:'
  write(out_unit,'(4f12.8)') V
  write(out_unit,*) ' Adiabatic vectors (in column) as a 4x4 matrix:'
  write(out_unit,'(4f12.8)') transpose(Vec)

  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(Q)
  deallocate(Vec)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'


END SUBROUTINE test_1DSOC_1S1T
SUBROUTINE test_PSB3
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: GGdef(:,:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: Vec(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: h(:,:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  pot_name = 'PSB3'
  ndim     = 3
  nsurf    = 2
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))


  Q(:) = (/0.1_Rkind,-3.14_Rkind,0.0_Rkind/)

  CALL sub_Qmodel_V(V,Q)
  write(out_unit,*) ' Diabatic potential as a 2x2 matrix:'
  write(out_unit,'(2f12.8)') V

  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  write(out_unit,*) '============================================================'
  write(out_unit,*) '== Optimisation ==='
  write(out_unit,*) '============================================================'

       !sub_Qmodel_Opt(Q,i_surf,nb_deg,icv,Max_it)
  CALL sub_Qmodel_Opt(Q,1,-1,3,-1)

  deallocate(V)
  deallocate(Q)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

END SUBROUTINE test_PSB3

SUBROUTINE test_Test
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT, out_unit=>OUTPUT_UNIT, Rkind => real64
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  pot_name = 'Test'
  ndim      = 0
  nsurf     = 0
  option    = 0
  adiabatic = .TRUE.

  CALL sub_Init_Qmodel_Cart(ndim,nsurf,pot_name,adiabatic,option)
  CALL set_Qmodel_Phase_Following(.FALSE.)
 
  write(out_unit,*) 'ndim,nsurf',ndim,nsurf
  flush(out_unit)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(g(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))


  Q = [1._Rkind,0._Rkind,0._Rkind,  0._Rkind,1._Rkind,0._Rkind,  0._Rkind,0._Rkind,1._Rkind]

  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)

  write(out_unit,*) ' Q:',Q

  write(out_unit,*) ' Adiabatic potential as a 2x2 matrix:'
  write(out_unit,'(2f12.8)') V
  write(out_unit,*) ' Adiabatic gradient. Each component as a 2x2 matrix:'
  DO i=1,ndim
    write(out_unit,*) ' Component:',i
    write(out_unit,'(2f12.8)') g(:,:,i)
  END DO

  write(out_unit,*) ' NAC. Each component as a 2x2 matrix:'
  DO i=1,ndim
    write(out_unit,*) ' Component:',i
    write(out_unit,'(2f12.8)') NAC(:,:,i)
  END DO

  deallocate(V)
  deallocate(g)
  deallocate(NAC)
  deallocate(Q)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

END SUBROUTINE test_Test
SUBROUTINE test_HBond
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: GGdef(:,:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: Vec(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: h(:,:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

  pot_name = 'HBond'
  ndim     = 2
  nsurf    = 1
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))


  Q(:) = (/2._Rkind,0._Rkind/)

  CALL sub_Qmodel_V(V,Q)
  write(out_unit,*) V


  allocate(G(nsurf,nsurf,ndim))
  allocate(h(nsurf,nsurf,ndim,ndim))

  CALL sub_Qmodel_VGH(V,G,H,Q)
  write(out_unit,*) 'gradiant',G
  write(out_unit,*) 'hessian',H

  CALL set_Qmodel_step(1.e-2_Rkind)
  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(Q)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'



END SUBROUTINE test_HBond

SUBROUTINE test_Phenol_Dia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: GGdef(:,:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: Vec(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: h(:,:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unit,*) '============================================================'
  write(out_unit,*) 'NTEST_driver. number of threads:',maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  pot_name = 'phenol'
  write(out_unit,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL QML_time_perso('Test ' // pot_name)

  ndim     = 2
  nsurf    = 3
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

  write(out_unit,*) 'ndim,nsurf',ndim,nsurf

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(g(nsurf,nsurf,ndim))
  allocate(h(nsurf,nsurf,ndim,ndim))

  allocate(GGdef(ndim,ndim))

  ! write GGdef
  write(out_unit,*) 'Get and Write the metric Tensor'
  CALL get_Qmodel_GGdef(GGdef)
  DO i=1,ndim
    write(out_unit,*) i,GGdef(:,i)
  END DO


  Q = [1.0_Rkind,-0.5_Rkind]
  CALL sub_Qmodel_V(V,Q)
  write(out_unit,*) ' Diabatic potential as a 3x3 matrix:'
  write(out_unit,'(3f12.8)') V

  CALL sub_Qmodel_Check_anaVSnum(Q,2)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '== Optimisation ==='
  write(out_unit,*) '============================================================'

       !sub_Qmodel_Opt(Q,i_surf,nb_deg,icv,Max_it)
  CALL sub_Qmodel_Opt(Q,1,-1,2,-1)
  CALL sub_Qmodel_VGH(V,G,H,Q)
  write(out_unit,*) 'Diabatic potential'
  write(out_unit,'(3f12.8)') V
  write(out_unit,*) 'Diabatic gradient of V(1,1)'
  write(out_unit,'(2f12.8)') g(1,1,:)
  write(out_unit,*) 'Diabatic hessian of V(1,1)'
  write(out_unit,'(2f12.8)') h(1,1,:,:)

  deallocate(V)
  deallocate(GGdef)
  deallocate(Q)


  CALL QML_time_perso('Test ' // pot_name)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,pot_name,option) &
!$OMP   PRIVATE(i,Q,V) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    CALL  random_number(Q)
    Q = Q + (/1.0_Rkind,-0.5_Rkind /)
    CALL sub_Qmodel_V(V,Q)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL QML_time_perso('Test ' // pot_name)
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'


END SUBROUTINE test_Phenol_Dia
SUBROUTINE test_Phenol_ADia
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: GGdef(:,:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: Vec(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: h(:,:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unit,*) '============================================================'
  write(out_unit,*) 'NTEST_driver. number of threads:',maxth



  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  pot_name = 'phenol'
  write(out_unit,*) ' Test of the adiatic potential ',pot_name


  ndim     = 2
  nsurf    = 3
  option   = 0
  adiabatic = .TRUE.

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(g(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)  ! a new initialization
  CALL set_Qmodel_Phase_Checking(.FALSE.)

  Q = [1._Rkind,-0.5_Rkind ]
  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)

  write(out_unit,*) ' Q:',Q
  write(out_unit,*) ' Adiabatic potential as a 3x3 matrix:'
  write(out_unit,'(3f12.8)') V
  write(out_unit,*) ' Adiabatic gradient. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(out_unit,*) ' Component:',i
    write(out_unit,'(3f12.8)') g(:,:,i)
  END DO

  write(out_unit,*) ' NAC. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(out_unit,*) ' Component:',i
    write(out_unit,'(3f12.8)') NAC(:,:,i)
  END DO

  Q = [1.1_Rkind,-0.5_Rkind ]
  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)

  write(out_unit,*) ' Q:',Q
  write(out_unit,*) ' Adiabatic potential as a 3x3 matrix:'
  write(out_unit,'(3f12.8)') V
  write(out_unit,*) ' Adiabatic gradient. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(out_unit,*) ' Component:',i
    write(out_unit,'(3f12.8)') g(:,:,i)
  END DO

  write(out_unit,*) ' NAC. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(out_unit,*) ' Component:',i
    write(out_unit,'(3f12.8)') NAC(:,:,i)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(g)
  deallocate(NAC)

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'


END SUBROUTINE test_Phenol_ADia
SUBROUTINE test_henonheiles(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: GGdef(:,:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: Vec(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: h(:,:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unit,*) '============================================================'
  write(out_unit,*) 'NTEST_driver. number of threads:',maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

  ndim     = 6
  nsurf    = 1
  option   = 0
  adiabatic = .TRUE.
  pot_name = 'henonheiles'
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)  ! a new initialization

  write(out_unit,*) '============================================================'
  write(out_unit,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL QML_time_perso('Test ' // pot_name)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,pot_name,option) &
!$OMP   PRIVATE(i,k,Q,V) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    CALL  random_number(Q)

    CALL sub_Qmodel_V(V,Q)

    !write(out_unit,*) Q,(V(k,k),k=1,nsurf)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL QML_time_perso('Test ' // pot_name)
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

END SUBROUTINE test_henonheiles

SUBROUTINE test_Vib_adia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: NAC(:,:,:)


  integer                             :: ndim,nsurf,nio_QML
  integer                             :: i,j,k,nb_eval,maxth

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unit,*) '============================================================'
  write(out_unit,*) '  Vibrational adiabatic separation'
  write(out_unit,*) 'TEST_driver. number of threads:',maxth

  ndim      = 0 ! it would be initialized
  nsurf     = 0 ! it would be initialized
  nio_QML   = in_unit
  CALL sub_Read_Qmodel(ndim,nsurf,nio_QML)

  write(out_unit,*) 'ndim,nsurf',ndim,nsurf

  write(out_unit,*) ' Test V, G and NAC'
  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  CALL  random_number(Q)
  Q = Q + 4.d0
  !Q=7.d0

  ! calculation of the adiabatic potential (as a matrix)
  ! calculation of the gradient and the NAC

  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)
  write(out_unit,*) 'Q',Q
  write(out_unit,*) 'Vadia',(V(i,i),i=1,nsurf)
  write(out_unit,*) 'Gadia',(G(i,i,1:ndim),i=1,nsurf)
  write(out_unit,*) 'NAC'
  DO i=1,nsurf
    write(out_unit,*) i,NAC(:,i,1:ndim)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(G)
  deallocate(NAC)

  write(out_unit,*) ' Test of ',nb_eval,' evaluations of the potential (read_model)'
  CALL QML_time_perso('Test read_model')

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth) &
!$OMP   PRIVATE(i,k,Q,V) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    CALL  random_number(Q)
    Q = Q + 4.d0

    CALL sub_Qmodel_V(V,Q)

    !write(out_unit,*) Q,(V(k,k),k=1,nsurf)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL QML_time_perso('Test read_model')
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

END SUBROUTINE test_Vib_adia

SUBROUTINE test_PH4
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT,Rkind => real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     :: GGdef(:,:)
  real (kind=Rkind),      allocatable     :: V(:,:)
  real (kind=Rkind),      allocatable     :: Vec(:,:)
  real (kind=Rkind),      allocatable     :: g(:,:,:)
  real (kind=Rkind),      allocatable     :: h(:,:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  pot_name = 'PH4'

  ndim     = 1
  nsurf    = 1
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

  DO i=-100,100

    Q(1) = 0.1_Rkind * i

    CALL sub_Qmodel_V(V,Q)
    write(out_unit,*) Q,V
  END DO


  deallocate(V)
  deallocate(Q)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'



END SUBROUTINE test_PH4