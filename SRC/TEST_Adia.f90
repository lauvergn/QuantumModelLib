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
PROGRAM TEST_adia
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE ADdnSVM_m
  USE Model_m
  USE mod_AdiaChannels
  IMPLICIT NONE


  TYPE (Model_t)                  :: QModel
  real (kind=Rkind), allocatable  :: Qact(:)
  TYPE (dnMat_t)                  :: PotVal,NAC
  real(kind=Rkind), allocatable   :: Ene(:)

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
    write(out_unitp,*) Qact,'Ene',(PotVal%d0(i,i)*auTOcm_inv,i=1,QML_get_nsurf_FROM_dnMat(PotVal))
  END DO



END PROGRAM TEST_adia
