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
PROGRAM TEST_adia
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE ADdnSVM_m
  USE Model_m
  !USE AdiaChannels_MakeHinact_m
  IMPLICIT NONE


  TYPE (Model_t)                  :: QModel
  real (kind=Rkind), allocatable  :: Qact(:)
  TYPE (dnMat_t)                  :: PotVal,NAC

  real(kind=Rkind) :: dQ
  integer :: i,iq
  integer, parameter :: nq=100
  real(kind=Rkind), parameter :: auTOcm_inv = 219474.631443_Rkind


  CALL Init_Model(QModel,Print_init=.FALSE.,Vib_adia=.TRUE.)

  dQ = TWO / real(nq,kind=Rkind)
  DO iq=0,nq
    Qact = [4.0_Rkind+iq*dQ]
    CALL Eval_Pot(QModel,Qact,PotVal,NAC=NAC,nderiv=1)
    !write(out_unitp,*) 'NAC'
    !CALL Write_RMat(NAC%d1(:,:,1),out_unitp,6,name_info='NAC')
    write(out_unitp,*) Qact,'Ene',(PotVal%d0(i,i)*auTOcm_inv,i=1,get_nsurf(PotVal))
  END DO



END PROGRAM TEST_adia
