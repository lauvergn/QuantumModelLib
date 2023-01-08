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
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
  IMPLICIT NONE

  CALL test_PSB3() ; stop

  CALL test_Read_Model() ; stop
  !CALL test_Phenol_ADia() ; stop
  !CALL test_HBond() ; stop
  CALL test_Vib_adia(1) ; stop


  CALL test_PH4()
  CALL test_HBond()
  CALL test_1DSOC_1S1T()
  CALL test_PSB3()
  CALL test_Phenol_Dia(10**7)
  CALL test_Phenol_ADia()
  CALL test_henonheiles(10**7)
  CALL test_Vib_adia(1000)

END PROGRAM main_pot
SUBROUTINE test_Read_Model()
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option
  real (kind=real64)                  :: DQ

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '  Read model'
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

  pot_name  = 'read_model'
  ndim      = 0 ! it would be initialized
  nsurf     = 0 ! it would be initialized
  option    = -1
  adiabatic = .TRUE.
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

  write(out_unitp,*) 'ndim,nsurf',ndim,nsurf
  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))

  IF (nsurf > 1) THEN
    write(out_unitp,*) ' Test V, G and NAC'

    allocate(NAC(nsurf,nsurf,ndim))

    CALL get_Qmodel_Q0(Q,0)
    DQ = 0.1_real64
    DO k=1,4
      write(out_unitp,*) '------------------------------------------------------'
      CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)
      write(out_unitp,*) 'Q',Q
      DO i=1,nsurf
        write(out_unitp,*) 'Vadia',i,':',V(i,i)
      END DO
      DO i=1,nsurf
        write(out_unitp,*) 'Gadia',i,':',G(i,i,1:ndim)
      END DO
      DO i=1,nsurf
        DO j=1,nsurf
          write(out_unitp,*) 'NAC',i,j,':',NAC(j,i,1:ndim)
        END DO
      END DO

      Q(:) = Q(:) + DQ
    END DO

    deallocate(NAC)
  ELSE
    write(out_unitp,*) ' Test V, G'

    CALL get_Qmodel_Q0(Q,0)

    CALL sub_Qmodel_VG(V,G,Q)
    write(out_unitp,*) 'Q',Q
    write(out_unitp,*) 'Potential',V
    write(out_unitp,*) 'Gradient',G
  END IF

  deallocate(Q)
  deallocate(V)
  deallocate(G)


  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

END SUBROUTINE test_Read_Model
SUBROUTINE test2_Vib_adia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  maxth = 1
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '  Vibrational adiabatic separation'
  write(out_unitp,*) 'TEST_driver. number of threads:',maxth

  pot_name  = 'read_model'
  ndim      = 0 ! it would be initialized
  nsurf     = 0 ! it would be initialized
  option    = -1
  adiabatic = .FALSE.
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

  write(out_unitp,*) 'ndim,nsurf',ndim,nsurf

  write(out_unitp,*) ' Test V, G and NAC'
  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))
  allocate(NAC(nsurf,nsurf,ndim))

  read(5,*) Q

  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)
  write(out_unitp,*) 'Q',Q
  write(out_unitp,*) 'Vadia',(V(i,i),i=1,nsurf)
  write(out_unitp,*) 'Gadia',(G(i,i,1:ndim),i=1,nsurf)
  write(out_unitp,*) 'NAC'
  DO i=1,nsurf
    write(out_unitp,*) i,NAC(:,i,1:ndim)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(G)
  deallocate(NAC)

  write(out_unitp,*) ' Test of evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

  Q = 2.5d0
  DO i=1,nb_eval
    Q = Q + 0.05d0

    CALL sub_Qmodel_V(V,Q)

    write(out_unitp,*) 'pot:',Q,(V(k,k),k=1,nsurf)
  END DO

  deallocate(Q)
  deallocate(V)

  CALL time_perso('Test ' // pot_name)
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

END SUBROUTINE test2_Vib_adia
SUBROUTINE test_1DSOC_1S1T
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  nb_eval = 1

  pot_name = '1DSOC_1S1T'
  write(out_unitp,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name

  ndim      = 1
  nsurf     = 4
  option    = 0
  adiabatic = .TRUE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(Vec(nsurf,nsurf))


  Q(:) = (/8.5_real64 /)

  CALL sub_Qmodel_VVec(V,Vec,Q)
  write(out_unitp,*) ' Diabatic potential as a 4x4 matrix:'
  write(out_unitp,'(4f12.8)') V
  write(out_unitp,*) ' Adiabatic vectors (in column) as a 4x4 matrix:'
  write(out_unitp,'(4f12.8)') transpose(Vec)

  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(Q)
  deallocate(Vec)


  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'


END SUBROUTINE test_1DSOC_1S1T
SUBROUTINE test_PSB3
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) 'NTEST_driver. number of threads:',maxth
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  pot_name = 'PSB3'
  write(out_unitp,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  ndim     = 3
  nsurf    = 2
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))


  Q(:) = (/0.1_real64,-3.14_real64,0.0_real64/)

  CALL sub_Qmodel_V(V,Q)
  write(out_unitp,*) ' Diabatic potential as a 2x2 matrix:'
  write(out_unitp,'(2f12.8)') V

  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '== Optimisation ==='
  write(out_unitp,*) '============================================================'

       !sub_Qmodel_Opt(Q,i_surf,nb_deg,icv,Max_it)
  CALL sub_Qmodel_Opt(Q,1,-1,3,-1)

  deallocate(V)
  deallocate(Q)


  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

END SUBROUTINE test_PSB3

SUBROUTINE test_HBond
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) 'NTEST_driver. number of threads:',maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  pot_name = 'HBond'
  write(out_unitp,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  ndim     = 2
  nsurf    = 1
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))


  Q(:) = (/2._real64,0._real64/)

  CALL sub_Qmodel_V(V,Q)
  write(out_unitp,*) V


  allocate(G(nsurf,nsurf,ndim))
  allocate(h(nsurf,nsurf,ndim,ndim))

  CALL sub_Qmodel_VGH(V,G,H,Q)
  write(out_unitp,*) 'gradiant',G
  write(out_unitp,*) 'hessian',H

  CALL set_Qmodel_step(1.e-2_real64)
  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(Q)


  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'



END SUBROUTINE test_HBond

SUBROUTINE test_Phenol_Dia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) 'NTEST_driver. number of threads:',maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  pot_name = 'phenol'
  write(out_unitp,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

  ndim     = 2
  nsurf    = 3
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)

  write(out_unitp,*) 'ndim,nsurf',ndim,nsurf

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))
  allocate(GGdef(ndim,ndim))

  ! write GGdef
  write(out_unitp,*) 'Get and Write the metric Tensor'
  CALL get_Qmodel_GGdef(GGdef)
  DO i=1,ndim
    write(out_unitp,*) i,GGdef(:,i)
  END DO


  Q = (/1.0_real64,-0.5_real64 /)
  CALL sub_Qmodel_V(V,Q)
  write(out_unitp,*) ' Diabatic potential as a 3x3 matrix:'
  write(out_unitp,'(3f12.8)') V

  CALL sub_Qmodel_Check_anaVSnum(Q,2)

  deallocate(V)
  deallocate(GGdef)
  deallocate(Q)


  CALL time_perso('Test ' // pot_name)

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(nb_eval,ndim,nsurf,maxth,pot_name,option) &
!$OMP   PRIVATE(i,Q,V) &
!$OMP   NUM_THREADS(maxth)

  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

!$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    CALL  random_number(Q)
    Q = Q + (/1.0_real64,-0.5_real64 /)
    CALL sub_Qmodel_V(V,Q)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL time_perso('Test ' // pot_name)
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'


END SUBROUTINE test_Phenol_Dia
SUBROUTINE test_Phenol_ADia
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) 'NTEST_driver. number of threads:',maxth



  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  pot_name = 'phenol'
  write(out_unitp,*) ' Test of the adiatic potential ',pot_name


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

  Q = [1._real64,-0.5_real64 ]
  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)

  write(out_unitp,*) ' Q:',Q
  write(out_unitp,*) ' Adiabatic potential as a 3x3 matrix:'
  write(out_unitp,'(3f12.8)') V
  write(out_unitp,*) ' Adiabatic gradient. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(out_unitp,*) ' Component:',i
    write(out_unitp,'(3f12.8)') g(:,:,i)
  END DO

  write(out_unitp,*) ' NAC. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(out_unitp,*) ' Component:',i
    write(out_unitp,'(3f12.8)') NAC(:,:,i)
  END DO

  Q = [1.1_real64,-0.5_real64 ]
  CALL sub_Qmodel_VG_NAC(V,G,NAC,Q)

  write(out_unitp,*) ' Q:',Q
  write(out_unitp,*) ' Adiabatic potential as a 3x3 matrix:'
  write(out_unitp,'(3f12.8)') V
  write(out_unitp,*) ' Adiabatic gradient. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(out_unitp,*) ' Component:',i
    write(out_unitp,'(3f12.8)') g(:,:,i)
  END DO

  write(out_unitp,*) ' NAC. Each component as a 3x3 matrix:'
  DO i=1,ndim
    write(out_unitp,*) ' Component:',i
    write(out_unitp,'(3f12.8)') NAC(:,:,i)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(g)
  deallocate(NAC)

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'


END SUBROUTINE test_Phenol_ADia
SUBROUTINE test_henonheiles(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) 'NTEST_driver. number of threads:',maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

  ndim     = 6
  nsurf    = 1
  option   = 0
  adiabatic = .TRUE.
  pot_name = 'henonheiles'
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)  ! a new initialization

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) ' Test of ',nb_eval,' evaluations of the potential ',pot_name
  CALL time_perso('Test ' // pot_name)

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

    !write(out_unitp,*) Q,(V(k,k),k=1,nsurf)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL time_perso('Test ' // pot_name)
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

END SUBROUTINE test_henonheiles

SUBROUTINE test_Vib_adia(nb_eval)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: NAC(:,:,:)


  integer                             :: ndim,nsurf,nio_QML
  integer                             :: i,j,k,nb_eval,maxth

  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  maxth = 1
  !$ maxth           = omp_get_max_threads()
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '  Vibrational adiabatic separation'
  write(out_unitp,*) 'TEST_driver. number of threads:',maxth

  ndim      = 0 ! it would be initialized
  nsurf     = 0 ! it would be initialized
  nio_QML   = in_unitp
  CALL sub_Read_Qmodel(ndim,nsurf,nio_QML)

  write(out_unitp,*) 'ndim,nsurf',ndim,nsurf

  write(out_unitp,*) ' Test V, G and NAC'
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
  write(out_unitp,*) 'Q',Q
  write(out_unitp,*) 'Vadia',(V(i,i),i=1,nsurf)
  write(out_unitp,*) 'Gadia',(G(i,i,1:ndim),i=1,nsurf)
  write(out_unitp,*) 'NAC'
  DO i=1,nsurf
    write(out_unitp,*) i,NAC(:,i,1:ndim)
  END DO

  deallocate(Q)
  deallocate(V)
  deallocate(G)
  deallocate(NAC)

  write(out_unitp,*) ' Test of ',nb_eval,' evaluations of the potential (read_model)'
  CALL time_perso('Test read_model')

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

    !write(out_unitp,*) Q,(V(k,k),k=1,nsurf)
  END DO
!$OMP   END DO

  deallocate(Q)
  deallocate(V)

!$OMP   END PARALLEL

  CALL time_perso('Test read_model')
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'

END SUBROUTINE test_Vib_adia

SUBROUTINE test_PH4
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
!$ USE omp_lib
  IMPLICIT NONE

  real (kind=real64),      allocatable     :: Q(:)
  real (kind=real64),      allocatable     :: GGdef(:,:)
  real (kind=real64),      allocatable     :: V(:,:)
  real (kind=real64),      allocatable     :: Vec(:,:)
  real (kind=real64),      allocatable     :: g(:,:,:)
  real (kind=real64),      allocatable     :: h(:,:,:,:)

  character (len=16)                  :: pot_name
  logical                             :: adiabatic

  integer                             :: i,j,k,ndim,nsurf,option,nb_eval,maxth

  nb_eval = 1
  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'
  pot_name = 'PH4'

  ndim     = 1
  nsurf    = 1
  option   = 0
  adiabatic = .FALSE.

  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)


  allocate(Q(ndim))
  allocate(V(nsurf,nsurf))

  DO i=-100,100

    Q(1) = 0.1_real64 * i

    CALL sub_Qmodel_V(V,Q)
    write(out_unitp,*) Q,V
  END DO


  deallocate(V)
  deallocate(Q)


  write(out_unitp,*) '============================================================'
  write(out_unitp,*) '============================================================'



END SUBROUTINE test_PH4


  SUBROUTINE time_perso(name)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unitp=>INPUT_UNIT,out_unitp=>OUTPUT_UNIT,real64
  IMPLICIT NONE

    character (len=*) :: name


    integer             :: tab_time(8) = 0
    real (kind=real64)  :: t_real
    integer             :: count,count_work,freq
    real                :: t_cpu
    integer             :: seconds,minutes,hours,days

    integer, save       :: count_old,count_ini
    real,    save       :: t_cpu_old,t_cpu_ini
    logical, save       :: begin = .TRUE.


!$OMP    CRITICAL (time_perso_CRIT)

    CALL date_and_time(values=tab_time)
    write(out_unitp,21) name,tab_time(5:8),tab_time(3:1:-1)
 21 format('     Time and date in ',a,' : ',i2,'h:',                &
            i2,'m:',i2,'.',i3,'s, the ',i2,'/',i2,'/',i4)

     CALL system_clock(count=count,count_rate=freq)
     call cpu_time(t_cpu)

     IF (begin) THEN
       begin = .FALSE.
       count_old = count
       count_ini = count
       t_cpu_old = t_cpu
       t_cpu_ini = t_cpu
     END IF


     !============================================
     !cpu time in the subroutine: "name"

     count_work = count-count_old
     seconds = count_work/freq

     minutes = seconds/60
     seconds = mod(seconds,60)
     hours   = minutes/60
     minutes = mod(minutes,60)
     days    = hours/24
     hours   = mod(hours,24)


     t_real = real(count_work,kind=real64)/real(freq,kind=real64)
     write(out_unitp,31) t_real,name
 31  format('        real (s): ',f18.3,' in ',a)
     write(out_unitp,32) days,hours,minutes,seconds,name
 32  format('        real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s in ',a)

     write(out_unitp,33) t_cpu-t_cpu_old,name
 33  format('        cpu (s): ',f18.3,' in ',a)


     !============================================
     !Total cpu time

     count_work = count-count_ini
     seconds = count_work/freq

     minutes = seconds/60
     seconds = mod(seconds,60)
     hours   = minutes/60
     minutes = mod(minutes,60)
     days    = hours/24
     hours   = mod(hours,24)

     t_real = real(count_work,kind=real64)/real(freq,kind=real64)
     write(out_unitp,41) t_real
 41  format('  Total real (s): ',f18.3)
     write(out_unitp,42) days,hours,minutes,seconds
 42  format('  Total real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s')
     write(out_unitp,43) t_cpu-t_cpu_ini
 43  format('  Total cpu (s): ',f18.3)

 51  format(a,i10,a,a)


     flush(6)
     !============================================

     count_old = count
     t_cpu_old = t_cpu

!$OMP   END CRITICAL (time_perso_CRIT)


  END SUBROUTINE time_perso
