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
!    Copyright 2016 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!    This particular module has been modified from ElVibRot-Tnum:
!       <http://pagesperso.lcp.u-psud.fr/lauvergnat/ElVibRot/ElVibRot.html>
!===========================================================================
!===========================================================================
!> @brief Program which checks the dnS module
!!
!! @licence GNU Lesser General Public License
!!
!! @section install_sec Installation
!!
!! Dependencies: this module needs the fortran modules in the @e Lib/Lib directory.
!!
!! Build the module and the testing unit (with dependencies):
!!
!!     make dnS.x
!!
!! Build the module documentation (with doxygen):
!!
!!     make doxy
!!
!! @section test_sec Tests
!!
!! To test the installation, you can run tests (dnS and ModLib).
!!
!!     cd Tests ; ./run_tests
!!
!! The results will be compared to previous ones in run_tests/RES_old
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
PROGRAM TEST_dnPoly
  USE mod_QML_NumParameters
  USE mod_dnS
  USE mod_dnPoly
  IMPLICIT NONE

    TYPE (dnS_t)                       :: dnX,dn2X,dnY,dnZ,Sana,Snum,dnXZ
    TYPE (dnS_t), allocatable          :: Vec_dnS(:)

    real (kind=Rkind)                :: x,y,z,err,maxdiff,maxdnS
    integer :: i,nderiv

    character (len=*), parameter :: name_sub='TEST_dnPoly'

   nderiv = 1
   write(out_unitp,'(a,i2)') "== TESTING dnPoly module with nderiv=",nderiv

   !nio_test = 10
   !open(unit=nio_test,file='dnSca.txt')


   x       = 0.5_Rkind
   dnX     = QML_init_dnS(x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
   dn2X    = QML_init_dnS(x+x,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives

   Sana  = QML_dnBox(dnX,2)
   CALL QML_write_dnS(Sana,out_unitp,info='dnBox(x,2)')

   allocate(Vec_dnS(6))

   write(out_unitp,*) '== Box(x,i) [0,Pi] =='
   Vec_dnS(:) = QML_dnBox(dnX,[1,2,3,4,5,6])
   DO i=1,size(Vec_dnS)
     CALL QML_write_dnS(Vec_dnS(i),out_unitp,info='dnBox')
   END DO

   write(out_unitp,*) '== x**i =='
   Vec_dnS(:) = QML_dnMonomial(dnX,[0,1,2,3,4,5])
   DO i=1,size(Vec_dnS)
     CALL QML_write_dnS(Vec_dnS(i),out_unitp,info='x**i')
   END DO

   write(out_unitp,*) '== Pl0(x,i) =='
   Vec_dnS(:) = QML_dnLegendre0(dnX,[0,1,2,3,4,5])
   DO i=1,size(Vec_dnS)
     CALL QML_write_dnS(Vec_dnS(i),out_unitp,info='Pl0')
   END DO

   write(out_unitp,*) '== Fourier(x,i) =='
   Vec_dnS(:) = QML_dnFourier(dnX,[1,2,3,4,5,6])
   DO i=1,size(Vec_dnS)
     CALL QML_write_dnS(Vec_dnS(i),out_unitp,info='Fourier')
   END DO

END PROGRAM TEST_dnPoly
