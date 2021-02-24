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
!!     cd Tests ; ./run_test_dnS
!!
!! The results will be compared to previous ones in run_tests/RES_old
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
PROGRAM TEST_dnS
  USE mod_QML_NumParameters
  USE mod_dnS
  IMPLICIT NONE

    TYPE (dnS_t)                       :: dnX,dn2X,dnY,dnZ,Sana,Snum,dnXZ
    TYPE (dnS_t), allocatable          :: Vec_dnS(:)

    real (kind=Rkind)                :: x,y,z,err,maxdiff,maxdnS
    integer                          :: nderiv,nio_test
    real (kind=Rkind)                :: dnSerr_test = FIVE*ONETENTH**4


    real (kind=Rkind), external  :: faplusx,faminusx,fatimex,faoverx
    real (kind=Rkind), external  :: SQRT_perso,ABS_perso,EXP_perso,LOG_perso,LOG10_perso
    real (kind=Rkind), external  :: SIN_perso,ASIN_perso,COS_perso,ACOS_perso,TAN_perso,ATAN_perso
    real (kind=Rkind), external  :: SINH_perso,ASINH_perso,COSH_perso,ACOSH_perso,TANH_perso,ATANH_perso

    character (len=*), parameter :: name_sub='TEST_dnS'

    CALL TEST_EXCEPTION()

    nderiv = 3
    write(out_unitp,'(a,i2)') "== TESTING dnS module with nderiv=",nderiv

   nio_test = 10
   open(unit=nio_test,file='dnSca.txt')


    x       = 0.5_Rkind
    dnX     = QML_init_dnS(x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
    dn2X    = QML_init_dnS(x+x,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "TEST with dnX=0.5"
   CALL QML_write_dnS(dnX,out_unitp,info='dnX')
   CALL QML_write_dnS(dnX,nio_test,info='dnX',FOR_test=.TRUE.)
   CALL QML_write_dnS(dn2X,out_unitp,info='dn2X')
   CALL QML_write_dnS(dn2X,nio_test,info='dn2X',FOR_test=.TRUE.)
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
   Snum = QML_get_Num_dnS_FROM_f_x(x,faplusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'a+dnX:         (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='0.5 + dnX',FOR_test=.TRUE.)

   Sana = dnX + 0.5_Rkind
   Snum = QML_get_Num_dnS_FROM_f_x(x,faplusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX+a:         (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX + 0.5',FOR_test=.TRUE.)

   Sana = dnX + dnX
   Snum = QML_get_Num_dnS_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX+dnX=2*dnX: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX + dnX',FOR_test=.TRUE.)

   Sana = +(0.5_Rkind - dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,faminusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') '+(a-dnX):      (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='+(0.5 - dnX)',FOR_test=.TRUE.)

   Sana = -(dnX - 0.5_Rkind)
   Snum = QML_get_Num_dnS_FROM_f_x(x,faminusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') '-(dnX-a):      (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='-(dnX - 0.5)',FOR_test=.TRUE.)

   Sana = dnX - dnX
   write(out_unitp,'(a,l2)') 'dnX-dnX                   ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX - dnX',FOR_test=.TRUE.)

   Sana = 2._Rkind * dnX
   Snum = QML_get_Num_dnS_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'a*dnX:         (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='2. * dnX',FOR_test=.TRUE.)

   Sana =  dnX * 2._Rkind
   Snum = QML_get_Num_dnS_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX*a:         (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX * 2.',FOR_test=.TRUE.)

   Sana =  dnX * dnX - dnX**(2._Rkind)
   write(out_unitp,'(a,l2)') 'dnX*dnX-dnX**2.           ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX * dnX - dnX**(2.)',FOR_test=.TRUE.)

   Sana =  dnX / 0.5_Rkind -dnX*2._Rkind
   write(out_unitp,'(a,l2)') 'dnX/0.5-dnX*2             ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   !CALL QML_write_dnS(Sana,out_unitp,info='dnX/0.5-dnX*2')
   CALL QML_write_dnS(Sana,nio_test,info='dnX/0.5 -dnX*2.',FOR_test=.TRUE.)

   Sana =  0.5_Rkind / dnX - 0.5_Rkind*dnX**(-1)
   write(out_unitp,'(a,l2)') '0.5/dnX-0.5.dnX*^-1       ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='0.5 / dnX - 0.5*dnX**(-1)',FOR_test=.TRUE.)

   Sana =  dnX / dnX - 1._Rkind
   write(out_unitp,'(a,l2)') 'dnX/dnX-1.                ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX / dnX - 1.',FOR_test=.TRUE.)

   Sana =  dnX**0.5_Rkind - sqrt(dnX)
   write(out_unitp,'(a,l2)') 'dnX**0.5-sqrt(dnX)        ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX**0.5-sqrt(dnX)',FOR_test=.TRUE.)

   Sana =  dnX**3 - dnX*dnX*dnX
   write(out_unitp,'(a,l2)') 'dnX**3-dnX*dnX*dnX        ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX**3 - dnX*dnX*dnX',FOR_test=.TRUE.)

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "functions: sqrt, exp, log, ... sin, asin, ... acosh ..."
   write(out_unitp,'(a)') "============================================"

   Sana = sqrt(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,SQRT_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sqrt: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on sqrt')
   CALL QML_write_dnS(Sana,nio_test,info='sqrt(dnX)',FOR_test=.TRUE.)

   Sana = abs(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,ABS_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'abs:  (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on abs')
   CALL QML_write_dnS(Sana,nio_test,info='abs(dnX)',FOR_test=.TRUE.)

   Sana = exp(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,EXP_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'exp:  (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on exp')
   CALL QML_write_dnS(Sana,nio_test,info='exp(dnX)',FOR_test=.TRUE.)

   Sana = log(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,log_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'log:  (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on log')
   CALL QML_write_dnS(Sana,nio_test,info='log(dnX)',FOR_test=.TRUE.)

   Sana = log10(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,log10_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'log10: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on log10')
   CALL QML_write_dnS(Sana,nio_test,info='log10(dnX)',FOR_test=.TRUE.)

   Sana = sin(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,sin_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sin:  (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on sin')
   CALL QML_write_dnS(Sana,nio_test,info='sin(dnX)',FOR_test=.TRUE.)

   Sana = asin(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,asin_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'asin: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on asin')
   CALL QML_write_dnS(Sana,nio_test,info='asin(dnX)',FOR_test=.TRUE.)

   Sana = cos(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,cos_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'cos:  (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on cos')
   CALL QML_write_dnS(Sana,nio_test,info='cos(dnX)',FOR_test=.TRUE.)

   Sana = acos(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,acos_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'acos: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on acos')
   CALL QML_write_dnS(Sana,nio_test,info='acos(dnX)',FOR_test=.TRUE.)

   Sana = tan(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,tan_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'tan:  (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on tan')
   CALL QML_write_dnS(Sana,nio_test,info='tan(dnX)',FOR_test=.TRUE.)

   Sana = atan(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,atan_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'atan: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on atan')
   CALL QML_write_dnS(Sana,nio_test,info='atan(dnX)',FOR_test=.TRUE.)

   Sana = sinh(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,sinh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sinh: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on sinh')
   CALL QML_write_dnS(Sana,nio_test,info='sinh(dnX)',FOR_test=.TRUE.)

   Sana = asinh(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,asinh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'asinh: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on asinh')
   CALL QML_write_dnS(Sana,nio_test,info='asinh(dnX)',FOR_test=.TRUE.)

   Sana = cosh(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,cosh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'cosh: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on cosh')
   CALL QML_write_dnS(Sana,nio_test,info='cosh(dnX)',FOR_test=.TRUE.)

   dnY = QML_init_dnS(FOUR*x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
   Sana = acosh(dnY)
   Snum = QML_get_Num_dnS_FROM_f_x(FOUR*x,acosh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'acosh: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on acosh')
   CALL QML_write_dnS(Sana,nio_test,info='acosh(dn4X)',FOR_test=.TRUE.)

   Sana = tanh(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,tanh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'tanh: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on tanh')
   CALL QML_write_dnS(Sana,nio_test,info='tanh(dnX)',FOR_test=.TRUE.)

   Sana = atanh(dnX)
   Snum = QML_get_Num_dnS_FROM_f_x(x,atanh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'atanh: (Sana-Snum)==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL QML_write_dnS(Sana,info='test on atanh')
   CALL QML_write_dnS(Sana,nio_test,info='atanh(dnX)',FOR_test=.TRUE.)

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "new tests: **, composition"
   write(out_unitp,'(a)') "============================================"

   Sana =  dnX**0 - ONE
   write(out_unitp,'(a,l2)') 'dnX**0-ONE                ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX**0 - ONE',FOR_test=.TRUE.)

   Sana =  dnX**1 - dnX
   write(out_unitp,'(a,l2)') 'dnX**1-dnX                ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX**1 - dnX',FOR_test=.TRUE.)

   Sana =  dnX**2 - dnX*dnX
   write(out_unitp,'(a,l2)') 'dnX**2-dnX*dnX            ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX**2 - dnX*dnX',FOR_test=.TRUE.)

   Sana =  dnX**3 - dnX*dnX*dnX
   write(out_unitp,'(a,l2)') 'dnX**3-dnX*dnX*dnX        ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='dnX**3 - dnX*dnX*dnX',FOR_test=.TRUE.)

   Sana =  sqrt(dnX**2) - dnX
   write(out_unitp,'(a,l2)') 'sqrt(dnX**2) - dnX        ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)
   CALL QML_write_dnS(Sana,nio_test,info='sqrt(dnX**2) - dnX',FOR_test=.TRUE.)

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "new tests: 3D, nderiv=2"
   write(out_unitp,'(a)') "============================================"
   nderiv = 2
   x=0.5_Rkind
   z=2.0_Rkind

   dnX     = QML_init_dnS(x  ,ndim=3,nderiv=nderiv,iQ=1) ! to set up the derivatives
   !CALL QML_write_dnS(dnX,out_unitp,info='dnX')
   dnZ     = QML_init_dnS(z  ,ndim=3,nderiv=nderiv,iQ=3) ! to set up the derivatives
   !CALL QML_write_dnS(dnZ,out_unitp,info='dnZ')
   Sana    = dnX*dnZ ! It is equivalent to 3D function f(x,y,z) = x*z
   CALL QML_set_dnS(dnXZ,d0=     x*z,                         &
                     d1=        [z,   ZERO,x],                &
                     d2=reshape([ZERO,ZERO,ONE,               &
                                 ZERO,ZERO,ZERO,              &
                                 ONE, ZERO,ZERO],shape=[3,3]))
   write(out_unitp,'(a,l2)') 'dnX*dnZ -dnS_Result       ==0?',QML_Check_dnS_IS_ZERO(Sana-dnXZ,dnSerr_test)
   !CALL QML_write_dnS(Sana,info='test on dnX*dnZ (3D)')
   !CALL QML_write_dnS(Sana,nio_test,info='dnX*dnZ (3D)',FOR_test=.TRUE.)

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "new tests : Vec_OF_dnS(1:3), 2D, nderiv=2"
   write(out_unitp,'(a)') "============================================"
   nderiv = 2
   x=0.5_Rkind
   y=1.0_Rkind
   dnX     = QML_init_dnS(x  ,ndim=2,nderiv=nderiv,iQ=1) ! to set up the derivatives
   dnY     = QML_init_dnS(y  ,ndim=2,nderiv=nderiv,iQ=2) ! to set up the derivatives

   allocate(Vec_dnS(3))
   Vec_dnS(:) = [dnX,dnY,dnX+dnY]
   Sana = dot_product(Vec_dnS,Vec_dnS) ! 2*(dnX**2+dnY**2+dnX*dnY)
   CALL QML_write_dnS(Sana,info='test on dot_product (2D)')
   CALL QML_write_dnS(Sana,nio_test,info='dot_product (2D)',FOR_test=.TRUE.)
   dnXZ = TWO*(dnX**2+dnY**2+dnX*dnY)
   write(out_unitp,'(a,l2)') 'dot_product - dnS_Result  ==0?',QML_Check_dnS_IS_ZERO(Sana-dnXZ,dnSerr_test)

   Vec_dnS(:) = [-dnX,-dnY,dnX+dnY]
   Sana = sum(Vec_dnS) ! ZERO
   !CALL QML_write_dnS(Sana,info='test on sum (2D)')
   !CALL QML_write_dnS(Sana,nio_test,info='sum (2D)',FOR_test=.TRUE.)
   write(out_unitp,'(a,l2)') 'sum  ==0?',QML_Check_dnS_IS_ZERO(Sana,dnSerr_test)

   CALL QML_dealloc_dnS(dnZ)
   dnZ = ONE
   Vec_dnS(:) = [dnX**2,ONE/dnX,ONE/dnX]
   Sana = product(Vec_dnS) ! ONE
   !CALL QML_write_dnS(Sana,info='test on product (2D)')
   !CALL QML_write_dnS(Sana,nio_test,info='product (2D)',FOR_test=.TRUE.)
   write(out_unitp,'(a,l2)') 'sum - 1 ==0?',QML_Check_dnS_IS_ZERO(Sana-dnZ,dnSerr_test)

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "new tests : two-argument fuction, 2D, nderiv=2"
   write(out_unitp,'(a)') "============================================"
   nderiv = 2
   x=-0.5_Rkind
   y=-1.0_Rkind
   dnX     = QML_init_dnS(x  ,ndim=2,nderiv=nderiv,iQ=1) ! to set up the derivatives
   dnY     = QML_init_dnS(y  ,ndim=2,nderiv=nderiv,iQ=2) ! to set up the derivatives
   Sana = atan2(dnY,dnX)
   CALL QML_write_dnS(Sana,info='test on atan2 (2D)')
   Snum = atan(dnY/dnX)-pi
   CALL QML_write_dnS(Snum,info='test on atan (2D)')
   write(out_unitp,'(a,l2)') 'atan2 - atan ==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)


   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "new tests : init const, 2D, nderiv=2"
   write(out_unitp,'(a)') "============================================"
   nderiv = 2
   x=0.5_Rkind
   y=1.0_Rkind
   dnX     = QML_init_dnS(x  ,ndim=2,nderiv=nderiv,iQ=1) ! to set up the derivatives

   CALL QML_dealloc_dnS(Sana)
   CALL QML_dealloc_dnS(Snum)
   Sana    = x ! here, Sana%nderiv = -1
   CALL QML_write_dnS(Sana,info='test1 (x)')
   Sana    = Sana + dnX ! here, Sana%nderiv = dnX%nderiv
   CALL QML_write_dnS(Sana,info='test2 (x+dnX)')
   CALL QML_write_dnS(Sana,nio_test,info='init const: (x+dnX)',FOR_test=.TRUE.)
   Snum    = dnX + x    ! here, Snum%nderiv = dnX%nderiv
   CALL QML_write_dnS(Snum,info='test3 (dnX+x)')
   write(out_unitp,'(a,l2)') '(x+dnX) - (dnX+x)  ==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)

   CALL QML_dealloc_dnS(Sana)
   CALL QML_dealloc_dnS(Snum)
   Sana    = x ! here, Sana%nderiv = -1
   CALL QML_write_dnS(Sana,info='test1 (x)')
   Sana    = Sana - dnX ! here, Sana%nderiv = dnX%nderiv
   CALL QML_write_dnS(Sana,info='test2 (x-dnX)')
   CALL QML_write_dnS(Sana,nio_test,info='init const: (x-dnX)',FOR_test=.TRUE.)
   Snum    = -dnX + x    ! here, Snum%nderiv = dnX%nderiv
   CALL QML_write_dnS(Snum,info='test3 (-dnX+x)')
   write(out_unitp,'(a,l2)') '(x-dnX) - (-dnX+x)  ==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)

   CALL QML_dealloc_dnS(Sana)
   CALL QML_dealloc_dnS(Snum)
   Sana    = x ! here, Sana%nderiv = -1
   CALL QML_write_dnS(Sana,info='test1 (x)')
   Sana    = Sana * dnX ! here, Sana%nderiv = dnX%nderiv
   CALL QML_write_dnS(Sana,info='test2 (x*dnX)')
   CALL QML_write_dnS(Sana,nio_test,info='init const: (x*dnX)',FOR_test=.TRUE.)
   Snum    = dnX * x    ! here, Snum%nderiv = dnX%nderiv
   CALL QML_write_dnS(Snum,info='test3 (dnX*x)')
   write(out_unitp,'(a,l2)') '(x*dnX) - (dnX*x)  ==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)

   CALL QML_dealloc_dnS(Sana)
   CALL QML_dealloc_dnS(Snum)
   Sana    = x ! here, Sana%nderiv = -1
   CALL QML_write_dnS(Sana,info='test1 (x)')
   Sana    = Sana / dnX ! here, Sana%nderiv = dnX%nderiv
   CALL QML_write_dnS(Sana,info='test2 (x/dnX)')
   CALL QML_write_dnS(Sana,nio_test,info='init const: (x/dnX)',FOR_test=.TRUE.)
   Snum    = dnX**(-1) * x    ! here, Snum%nderiv = dnX%nderiv
   CALL QML_write_dnS(Snum,info='test3 (dnX**(-1)*x)')
   write(out_unitp,'(a,l2)') '(x/dnX) - (dnX**(-1)*x)  ==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)

   CALL QML_dealloc_dnS(Sana)
   CALL QML_dealloc_dnS(Snum)
   Sana    = x ! here, Sana%nderiv = -1
   CALL QML_write_dnS(Sana,info='test1 (x)')
   Sana    = sin(Sana) ! here, Sana%nderiv = dnX%nderiv
   CALL QML_write_dnS(Sana,info='test2 (sin(x)')
   CALL QML_write_dnS(Sana,nio_test,info='init const: (sin(x))',FOR_test=.TRUE.)
   Snum    = dnX ; Snum = ZERO
   Snum    = sin(x)    ! here, Snum%nderiv = dnX%nderiv
   CALL QML_write_dnS(Snum,info='test3 (sin(x))')
   write(out_unitp,'(a,l2)') '(sin(x)) - (sin(x))  ==0?',QML_Check_dnS_IS_ZERO(Sana-Snum,dnSerr_test)

   close(unit=nio_test)


END PROGRAM TEST_dnS
SUBROUTINE TEST_EXCEPTION
  USE, INTRINSIC :: ieee_exceptions

  USE mod_QML_NumParameters
  USE mod_dnS
  IMPLICIT NONE

    TYPE (dnS_t)                     :: dnX
    real (kind=Rkind)                :: x
    integer                          :: nderiv
    logical                          :: flag_NAN


    character (len=*), parameter :: name_sub='TEST_EXCEPTION'

    nderiv = 3
    write(out_unitp,'(a,i2)') "== TESTING EXCEPTION dnS module with nderiv=",nderiv

    CALL ieee_set_halting_mode(ieee_invalid, .FALSE.)

    x = -ONE
    dnX     = QML_init_dnS(x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives

    write(out_unitp,*) 'log(-1)',log(x)
    CALL  IEEE_GET_FLAG(ieee_invalid, flag_NAN)
    IF (flag_NAN) THEN
       write(out_unitp,*) 'a NAN occurs!'
    ELSE
      write(out_unitp,*) 'no NAN occurs!'
    END IF

    ! reset the NAN exception
    CALL  IEEE_SET_FLAG(ieee_invalid, .FALSE.)

    CALL QML_write_dnS(cos(dnX),info='cos(dnX)')
    CALL  IEEE_GET_FLAG(ieee_invalid, flag_NAN)
    IF (flag_NAN) THEN
       write(out_unitp,*) 'a NAN occurs!'
    ELSE
      write(out_unitp,*) 'no NAN occurs!'
    END IF

    ! reset the NAN exception
    CALL  IEEE_SET_FLAG(ieee_invalid, .FALSE.)

    CALL QML_write_dnS(log(dnX),info='log(dnX)')
    CALL QML_write_dnS(cos(dnX),info='cos(dnX)')
    CALL  IEEE_GET_FLAG(ieee_invalid, flag_NAN)
    IF (flag_NAN) THEN
       write(out_unitp,*) 'a NAN occurs!'
    ELSE
      write(out_unitp,*) 'no NAN occurs!'
    END IF

    !write(out_unitp,*) 'dnX',dnX


END SUBROUTINE TEST_EXCEPTION

FUNCTION faplusx(x)
  USE mod_QML_NumParameters
  IMPLICIT NONE

  real (kind=Rkind) :: faplusx
  real (kind=Rkind), intent(in) :: x

  faplusx = 0.5_Rkind + x

END FUNCTION faplusx
FUNCTION faminusx(x)
  USE mod_QML_NumParameters
  IMPLICIT NONE

  real (kind=Rkind) :: faminusx
  real (kind=Rkind), intent(in) :: x

  faminusx = 0.5_Rkind - x

END FUNCTION faminusx
FUNCTION fatimex(x)
  USE mod_QML_NumParameters
  IMPLICIT NONE

  real (kind=Rkind) :: fatimex
  real (kind=Rkind), intent(in) :: x

  fatimex = 2._Rkind * x

END FUNCTION fatimex
FUNCTION faoverx(x)
  USE mod_QML_NumParameters
  IMPLICIT NONE

  real (kind=Rkind) :: faoverx
  real (kind=Rkind), intent(in) :: x

  faoverx = 0.5_Rkind / x

END FUNCTION faoverx
FUNCTION SQRT_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = sqrt(x)

END FUNCTION SQRT_perso
FUNCTION ABS_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = ABS(x)

END FUNCTION ABS_perso
FUNCTION EXP_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = EXP(x)

END FUNCTION EXP_perso
FUNCTION LOG_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = LOG(x)

END FUNCTION LOG_perso
FUNCTION LOG10_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = LOG10(x)

END FUNCTION LOG10_perso
FUNCTION SIN_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = SIN(x)

END FUNCTION SIN_perso
FUNCTION ASIN_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = ASIN(x)

END FUNCTION ASIN_perso
FUNCTION COS_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = COS(x)

END FUNCTION COS_perso
FUNCTION ACOS_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = ACOS(x)

END FUNCTION ACOS_perso
FUNCTION TAN_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = TAN(x)

END FUNCTION TAN_perso
FUNCTION ATAN_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = ATAN(x)

END FUNCTION ATAN_perso
FUNCTION SINH_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = SINH(x)

END FUNCTION SINH_perso
FUNCTION ASINH_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

#if __INVHYP == 1
    f = asinh(x)
#else
    f = log(x+sqrt(x*x+ONE))
#endif

END FUNCTION ASINH_perso
FUNCTION COSH_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = COSH(x)

END FUNCTION COSH_perso
FUNCTION ACOSH_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

#if __INVHYP == 1
    f = acosh(x)
#else
    f = log(x+sqrt(x*x-ONE))
#endif

END FUNCTION ACOSH_perso
FUNCTION TANH_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = TANH(x)

END FUNCTION TANH_perso
FUNCTION ATANH_perso(x) RESULT(f)
  USE mod_QML_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

#if __INVHYP == 1
    f = atanh(x)
#else
    f = HALF*log((ONE+x)/(ONE-x))
#endif

END FUNCTION ATANH_perso
