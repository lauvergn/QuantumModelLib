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
PROGRAM test_driver
  IMPLICIT NONE

  CALL test_Read_Model()

END PROGRAM test_driver
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