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
PROGRAM TEST_dnS
  USE mod_NumParameters
  USE mod_dnS
  IMPLICIT NONE

    TYPE (dnS)                       :: dnX,dn2X,dnY,dnZ,Sana,Snum,dnXZ
    TYPE (dnS), allocatable          :: Vec_dnS(:)

    real (kind=Rkind)                :: x,y,z,err,maxdiff,maxdnS
    integer                          :: nderiv,nio_test
    real (kind=Rkind)                :: dnSerr_test = FIVE*ONETENTH**4

    !intrinsic :: sqrt, abs, exp, log, log10, sin, asin, cos, acos, tan, atan, sinh, asinh, cosh, acosh, tanh, atanh
    intrinsic :: dsqrt, dabs, dexp, dlog, dlog10
    intrinsic :: dsin, dasin, dcos, dacos, dtan, datan
    intrinsic :: dsinh, dcosh, dtanh
#if __INVHYP == 1
    intrinsic :: dasinh, dacosh, datanh
#endif
    real (kind=Rkind), external  :: faplusx,faminusx,fatimex,faoverx

    character (len=*), parameter :: name_sub='TEST_dnS'

   nderiv = 3
   write(out_unitp,'(a,i2)') "== TESTING dnS module with nderiv=",nderiv

   nio_test = 10
   open(unit=nio_test,file='dnSca.txt')


    x       = 0.5_Rkind
    dnX     = init_dnS(x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
    dn2X    = init_dnS(x+x,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "TEST with dnX=0.5"
   CALL Write_dnS(dnX,out_unitp,info='dnX')
   CALL Write_dnS_FOR_test(dnX,nio_test,info='dnX')
   CALL Write_dnS(dn2X,out_unitp,info='dn2X')
   CALL Write_dnS_FOR_test(dn2X,nio_test,info='dn2X')
   write(out_unitp,'(a)') "============================================"

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "operators: == /= > >= < <="
   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a,l2)') 'dnX == dnX:   T ?',(dnX == dnX)
   write(out_unitp,'(a,l2)') 'dnX == dn2X:  F ?',(dnX == dn2X)

   write(out_unitp,'(a,l2)') 'dnX /= dnX:   F ?',(dnX /= dnX)
   write(out_unitp,'(a,l2)') 'dnX /= dn2X:  T ?',(dnX /= dn2X)

   write(out_unitp,'(a,l2)') 'dnX  > dnX:   F ?',(dnX > dnX)
   write(out_unitp,'(a,l2)') 'dnX  > dn2X:  F ?',(dnX > dn2X)
   write(out_unitp,'(a,l2)') 'dn2X > dnX:   T ?',(dn2X > dnX)

   write(out_unitp,'(a,l2)') 'dnX  >= dnX:  T ?',(dnX >= dnX)
   write(out_unitp,'(a,l2)') 'dnX  >= dn2X: F ?',(dnX >= dn2X)
   write(out_unitp,'(a,l2)') 'dn2X >= dnX:  T ?',(dn2X >= dnX)

   write(out_unitp,'(a,l2)') 'dnX  < dnX:   F ?',(dnX < dnX)
   write(out_unitp,'(a,l2)') 'dnX  < dn2X:  T ?',(dnX < dn2X)
   write(out_unitp,'(a,l2)') 'dn2X < dnX:   F ?',(dn2X < dnX)

   write(out_unitp,'(a,l2)') 'dnX  <= dnX:  T ?',(dnX <= dnX)
   write(out_unitp,'(a,l2)') 'dnX  <= dn2X: T ?',(dnX <= dn2X)
   write(out_unitp,'(a,l2)') 'dn2X <= dnX:  F ?',(dn2X <= dnX)

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "operators: .EQ."
   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a,l2)') 'dnX .EQ. dnX:   T ?',(dnX .EQ. dnX)
   write(out_unitp,'(a,l2)') 'dnX .NE.dn2X:   T ?',(dnX .NE. dn2X)

   write(out_unitp,'(a,l2)') 'dnX .GT. dnX:   F ?',(dnX .GT. dnX)
   write(out_unitp,'(a,l2)') 'dnX .GE. dnX:   T ?',(dnX .GE. dnX)

   write(out_unitp,'(a,l2)') 'dnX .LT. dnX:   F ?',(dnX .LT. dnX)
   write(out_unitp,'(a,l2)') 'dnX .LE. dnX:   T ?',(dnX .LE. dnX)

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "operators: + - * / **"
   write(out_unitp,'(a)') "============================================"

   Sana = 0.5_Rkind + dnX
   Snum = get_Num_dnS_FROM_f_x(x,faplusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'a+dnX:         (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='0.5 + dnX')

   Sana = dnX + 0.5_Rkind
   Snum = get_Num_dnS_FROM_f_x(x,faplusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX+a:         (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX + 0.5')

   Sana = dnX + dnX
   Snum = get_Num_dnS_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX+dnX=2*dnX: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX + dnX')

   Sana = +(0.5_Rkind - dnX)
   Snum = get_Num_dnS_FROM_f_x(x,faminusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') '+(a-dnX):      (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='+(0.5 - dnX)')

   Sana = -(dnX - 0.5_Rkind)
   Snum = get_Num_dnS_FROM_f_x(x,faminusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') '-(dnX-a):      (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='-(dnX - 0.5)')

   Sana = dnX - dnX
   write(out_unitp,'(a,l2)') 'dnX-dnX                   ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX - dnX')

   Sana = 2._Rkind * dnX
   Snum = get_Num_dnS_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'a*dnX:         (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='2. * dnX')

   Sana =  dnX * 2._Rkind
   Snum = get_Num_dnS_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX*a:         (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX * 2.')

   Sana =  dnX * dnX - dnX**(2._Rkind)
   write(out_unitp,'(a,l2)') 'dnX*dnX-dnX**2.           ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX * dnX - dnX**(2.)')

   Sana =  dnX / 0.5_Rkind -dnX*2._Rkind
   write(out_unitp,'(a,l2)') 'dnX/0.5-dnX*2             ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   !CALL Write_dnS(Sana,out_unitp,info='dnX/0.5-dnX*2')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX/0.5 -dnX*2.')

   Sana =  0.5_Rkind / dnX - 0.5_Rkind*dnX**(-1)
   write(out_unitp,'(a,l2)') '0.5/dnX-0.5.dnX*^-1       ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='0.5 / dnX - 0.5*dnX**(-1)')

   Sana =  dnX / dnX - 1._Rkind
   write(out_unitp,'(a,l2)') 'dnX/dnX-1.                ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX / dnX - 1.')

   Sana =  dnX**0.5_Rkind - sqrt(dnX)
   write(out_unitp,'(a,l2)') 'dnX**0.5-sqrt(dnX)        ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX**0.5-sqrt(dnX)')

   Sana =  dnX**3 - dnX*dnX*dnX
   write(out_unitp,'(a,l2)') 'dnX**3-dnX*dnX*dnX        ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX**3 - dnX*dnX*dnX')

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "functions: sqrt, exp, log, ... sin, asin, ... acosh ..."
   write(out_unitp,'(a)') "============================================"

   Sana = sqrt(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dsqrt,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sqrt: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on sqrt')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='sqrt(dnX)')

   Sana = abs(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dabs,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'abs:  (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on abs')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='abs(dnX)')

   Sana = exp(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dexp,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'exp:  (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on exp')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='exp(dnX)')

   Sana = log(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dlog,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'log:  (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on log')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='log(dnX)')

   Sana = log10(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dlog10,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'log10: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on log10')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='log10(dnX)')

   Sana = sin(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dsin,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sin:  (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on sin')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='sin(dnX)')

   Sana = asin(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dasin,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'asin: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on asin')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='asin(dnX)')

   Sana = cos(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dcos,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'cos:  (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on cos')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='cos(dnX)')

   Sana = acos(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dacos,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'acos: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on acos')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='acos(dnX)')

   Sana = tan(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dtan,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'tan:  (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on tan')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='tan(dnX)')

   Sana = atan(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,datan,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'atan: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on atan')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='atan(dnX)')

   Sana = sinh(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dsinh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sinh: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on sinh')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='sinh(dnX)')

   Sana = asinh(dnX)
#if __INVHYP == 1
   Snum = get_Num_dnS_FROM_f_x(x,dasinh,nderiv=nderiv)
#else
   Snum = get_Num_dnS_FROM_f_x(x,asinh_perso,nderiv=nderiv)
#endif
   write(out_unitp,'(a,l2)') 'asinh: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on asinh')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='asinh(dnX)')

   Sana = cosh(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dcosh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'cosh: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on cosh')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='cosh(dnX)')

   dnY = init_dnS(FOUR*x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
   Sana = acosh(dnY)
#if __INVHYP == 1
   Snum = get_Num_dnS_FROM_f_x(FOUR*x,dacosh,nderiv=nderiv)
#else
   Snum = get_Num_dnS_FROM_f_x(FOUR*x,acosh_perso,nderiv=nderiv)
#endif
   write(out_unitp,'(a,l2)') 'acosh: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on acosh')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='acosh(dn4X)')

   Sana = tanh(dnX)
   Snum = get_Num_dnS_FROM_f_x(x,dtanh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'tanh: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on tanh')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='tanh(dnX)')

   Sana = atanh(dnX)
#if __INVHYP == 1
   Snum = get_Num_dnS_FROM_f_x(x,datanh,nderiv=nderiv)
#else
   Snum = get_Num_dnS_FROM_f_x(x,atanh_perso,nderiv=nderiv)
#endif
   write(out_unitp,'(a,l2)') 'atanh: (Sana-Snum)==0?',Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnS(Sana,info='test on atanh')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='atanh(dnX)')

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "new tests: **, composition"
   write(out_unitp,'(a)') "============================================"

   Sana =  dnX**0 - ONE
   write(out_unitp,'(a,l2)') 'dnX**0-ONE                ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX**0 - ONE')

   Sana =  dnX**1 - dnX
   write(out_unitp,'(a,l2)') 'dnX**1-dnX                ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX**1 - dnX')

   Sana =  dnX**2 - dnX*dnX
   write(out_unitp,'(a,l2)') 'dnX**2-dnX*dnX            ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX**2 - dnX*dnX')

   Sana =  dnX**3 - dnX*dnX*dnX
   write(out_unitp,'(a,l2)') 'dnX**3-dnX*dnX*dnX        ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX**3 - dnX*dnX*dnX')

   Sana =  sqrt(dnX**2) - dnX
   write(out_unitp,'(a,l2)') 'sqrt(dnX**2) - dnX        ==0?',Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnS_FOR_test(Sana,nio_test,info='sqrt(dnX**2) - dnX')

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "new tests: 3D, nderiv=2"
   write(out_unitp,'(a)') "============================================"
   nderiv = 2
   x=0.5_Rkind
   z=2.0_Rkind

   dnX     = init_dnS(x  ,ndim=3,nderiv=nderiv,iQ=1) ! to set up the derivatives
   !CALL Write_dnS(dnX,out_unitp,info='dnX')
   dnZ     = init_dnS(z  ,ndim=3,nderiv=nderiv,iQ=3) ! to set up the derivatives
   !CALL Write_dnS(dnZ,out_unitp,info='dnZ')
   Sana    = dnX*dnZ ! It is equivalent to 3D function f(x,y,z) = x*z
   CALL set_dnS(dnXZ,d0=         x*z,                         &
                     d1=        [z,   ZERO,x],                &
                     d2=reshape([ZERO,ZERO,ONE,               &
                                 ZERO,ZERO,ZERO,              &
                                 ONE, ZERO,ZERO],shape=[3,3]))
   write(out_unitp,'(a,l2)') 'dnX*dnZ -dnS_Result       ==0?',Check_dnS_IS_ZERO(Sana-dnXZ,dnSerr_test)
   !CALL write_dnS(Sana,info='test on dnX*dnZ (3D)')
   !CALL Write_dnS_FOR_test(Sana,nio_test,info='dnX*dnZ (3D)')

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "new tests : Vec_OF_dnS(1:3), 2D, nderiv=1"
   write(out_unitp,'(a)') "============================================"
   nderiv = 2
   x=0.5_Rkind
   y=1.0_Rkind
   dnX     = init_dnS(x  ,ndim=2,nderiv=nderiv,iQ=1) ! to set up the derivatives
   dnY     = init_dnS(y  ,ndim=2,nderiv=nderiv,iQ=2) ! to set up the derivatives

   allocate(Vec_dnS(3))
   Vec_dnS(:) = [dnX,dnY,dnX+dnY]
   Sana = dot_product(Vec_dnS,Vec_dnS) ! 2*(dnX**2+dnY**2+dnX*dnY)
   CALL write_dnS(Sana,info='test on dot_product (2D)')
   CALL Write_dnS_FOR_test(Sana,nio_test,info='dot_product (2D)')
   dnXZ = TWO*(dnX**2+dnY**2+dnX*dnY)
   write(out_unitp,'(a,l2)') 'dot_product - dnS_Result  ==0?',Check_dnS_IS_ZERO(Sana-dnXZ,dnSerr_test)


   close(unit=nio_test)


END PROGRAM TEST_dnS
FUNCTION faplusx(x)
  USE mod_NumParameters
  IMPLICIT NONE

  real (kind=Rkind) :: faplusx
  real (kind=Rkind), intent(in) :: x

  faplusx = 0.5_Rkind + x

END FUNCTION faplusx
FUNCTION faminusx(x)
  USE mod_NumParameters
  IMPLICIT NONE

  real (kind=Rkind) :: faminusx
  real (kind=Rkind), intent(in) :: x

  faminusx = 0.5_Rkind - x

END FUNCTION faminusx
FUNCTION fatimex(x)
  USE mod_NumParameters
  IMPLICIT NONE

  real (kind=Rkind) :: fatimex
  real (kind=Rkind), intent(in) :: x

  fatimex = 2._Rkind * x

END FUNCTION fatimex
FUNCTION faoverx(x)
  USE mod_NumParameters
  IMPLICIT NONE

  real (kind=Rkind) :: faoverx
  real (kind=Rkind), intent(in) :: x

  faoverx = 0.5_Rkind / x

END FUNCTION faoverx
