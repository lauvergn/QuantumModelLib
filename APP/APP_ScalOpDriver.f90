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

  !CALL test_Morse()
  CALL test_OneD_2Quadra()
 
END PROGRAM test_driver
SUBROUTINE test_Morse()
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT, Rkind => real64
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     ::   ScalOp(:,:,:)
  real (kind=Rkind),      allocatable     ::  dScalOp(:,:,:,:)
  real (kind=Rkind),      allocatable     :: ddScalOp(:,:,:,:,:)

  character (len=16)                  :: model_name
  logical                             :: adiabatic

  integer                             :: iOp,ndim,nsurf,nb_ScalOp,option

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '  Morse (HF) model'
  write(out_unit,*) '    with Scalar Operators (potential+dipole)'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

  model_name = 'Morse'
  ndim       = 0 ! it would be initialized
  nsurf      = 0 ! it would be initialized
  nb_ScalOp  = 2
  option     = -1
  adiabatic  = .FALSE.
  CALL sub_Init_Qmodel_ScalOp(ndim,nsurf,nb_ScalOp,model_name,adiabatic,option)

  write(out_unit,*) 'ndim,nsurf,nb_ScalOp',ndim,nsurf,nb_ScalOp
  allocate(Q(ndim))
  allocate(ScalOp(nsurf,nsurf,nb_ScalOp))
  allocate(dScalOp(nsurf,nsurf,ndim,nb_ScalOp))
  allocate(ddScalOp(nsurf,nsurf,ndim,ndim,nb_ScalOp))

  CALL get_Qmodel_Q0(Q,0)


  CALL sub_Qmodel_ddScalOp(ScalOp,dScalOp,ddScalOp,Q)
  write(out_unit,*) ' Test Scalar Operators'
  write(out_unit,*) '   at, Q:',Q(:)

  DO iOp=1,nb_ScalOp
    write(out_unit,*) '------------------------------------------------------'
    write(out_unit,*) 'iOp',iOp
    write(out_unit,*) '------------------------------------------------------'
    write(out_unit,*) iOp,'Value:    ',ScalOp(:,:,iOp)
    write(out_unit,*) iOp,'gradient: ',dScalOp(:,:,:,iOp)
    write(out_unit,*) iOp,'hessian:  ',ddScalOp(:,:,:,:,iOp)
  END DO

  deallocate(Q)
  deallocate(ScalOp)
  deallocate(dScalOp)
  deallocate(ddScalOp)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

END SUBROUTINE test_Morse
SUBROUTINE test_OneD_2Quadra()
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : in_unit=>INPUT_UNIT,out_unit=>OUTPUT_UNIT, Rkind => real64
  IMPLICIT NONE

  real (kind=Rkind),      allocatable     :: Q(:)
  real (kind=Rkind),      allocatable     ::   ScalOp(:,:,:)
  real (kind=Rkind),      allocatable     ::  dScalOp(:,:,:,:)
  real (kind=Rkind),      allocatable     :: ddScalOp(:,:,:,:,:)

  character (len=16)                  :: model_name
  logical                             :: adiabatic
  integer                             :: iOp,ndim,nsurf,nb_ScalOp,option

  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '  OneD_2Quadra model (1D and 2PES)'
  write(out_unit,*) '    with Scalar Operators (potential+dipole)'
  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

  model_name = 'OneD_2Quadra'
  ndim       = 0 ! it would be initialized
  nsurf      = 0 ! it would be initialized
  nb_ScalOp  = 2 ! it would be initialized
  option     = -1
  adiabatic  = .TRUE.
  CALL sub_Init_Qmodel_ScalOp(ndim,nsurf,nb_ScalOp,model_name,adiabatic,option)

  write(out_unit,*) 'ndim,nsurf,nb_ScalOp',ndim,nsurf,nb_ScalOp
  allocate(Q(ndim))
  allocate(ScalOp(nsurf,nsurf,nb_ScalOp))
  ScalOp = 0._Rkind
  allocate(dScalOp(nsurf,nsurf,ndim,nb_ScalOp))
  dScalOp = 0._Rkind
  allocate(ddScalOp(nsurf,nsurf,ndim,ndim,nb_ScalOp))
  ddScalOp = 0._Rkind

  Q = 3._Rkind

  !CALL sub_Qmodel_ddScalOp(ScalOp,dScalOp,ddScalOp,Q)
  CALL sub_Qmodel_ScalOp(ScalOp,Q)
  write(out_unit,*) ' Test Scalar Operators'
  write(out_unit,*) '   at, Q:',Q(:)

  DO iOp=1,nb_ScalOp
    write(out_unit,*) '------------------------------------------------------'
    write(out_unit,*) 'iOp',iOp
    write(out_unit,*) '------------------------------------------------------'
    write(out_unit,*) iOp,'Value:    ',ScalOp(:,:,iOp)
    !write(out_unit,*) iOp,'gradient: ',dScalOp(:,:,:,iOp)
    !write(out_unit,*) iOp,'hessian:  ',ddScalOp(:,:,:,:,iOp)
  END DO

  deallocate(Q)
  deallocate(ScalOp)
  deallocate(dScalOp)
  deallocate(ddScalOp)


  write(out_unit,*) '============================================================'
  write(out_unit,*) '============================================================'

END SUBROUTINE test_OneD_2Quadra