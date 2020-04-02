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
!> @brief Program which checks the dnSca module
!!
!! @licence GNU Lesser General Public License
!!
!! @section install_sec Installation
!!
!! Dependencies: this module needs the fortran modules in the @e Lib/Lib directory.
!!
!! Build the module and the testing unit (with dependencies):
!!
!!     make dnSca.x
!!
!! Build the module documentation (with doxygen):
!!
!!     make doxy
!!
!! @section test_sec Tests
!!
!! To test the installation, you can run tests (dnSca and ModLib).
!!
!!     cd Tests ; ./run_tests
!!
!! The results will be compared to previous ones in run_tests/RES_old
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
PROGRAM TEST_dnSca
  USE mod_NumParameters
  USE mod_dnSca
  IMPLICIT NONE

    TYPE (dnSca)                     :: dnX,dn2X,dnY,Sana,Snum
    real (kind=Rkind)                :: x,err,maxdiff,maxdnSca
    integer                          :: nderiv,nio_test
    real (kind=Rkind)                :: dnSerr_test = FIVE*ONETENTH**4


    real (kind=Rkind), external  :: faplusx,faminusx,fatimex,faoverx
    real (kind=Rkind), external  :: SQRT_perso,ABS_perso,EXP_perso,LOG_perso,LOG10_perso
    real (kind=Rkind), external  :: SIN_perso,ASIN_perso,COS_perso,ACOS_perso,TAN_perso,ATAN_perso
    real (kind=Rkind), external  :: SINH_perso,ASINH_perso,COSH_perso,ACOSH_perso,TANH_perso,ATANH_perso

    character (len=*), parameter :: name_sub='TEST_dnSca'

   nderiv = 3
   write(out_unitp,'(a,i2)') "== TESTING dnSca module with nderiv=",nderiv

   nio_test = 10
   open(unit=nio_test,file='dnSca.txt')


    x       = 0.5_Rkind
    dnX     = init_dnSca(x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
    dn2X    = init_dnSca(x+x,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "TEST with dnX=0.5"
   CALL Write_dnSca(dnX,out_unitp,info='dnX')
   CALL Write_dnSca_FOR_test(dnX,nio_test,info='dnX')
   CALL Write_dnSca(dn2X,out_unitp,info='dn2X')
   CALL Write_dnSca_FOR_test(dn2X,nio_test,info='dn2X')
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
   ! test 2+1=3
   Sana = 0.5_Rkind + dnX
   Snum = get_Num_dnSca_FROM_f_x(x,faplusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'a+dnX:         (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='0.5 + dnX')

   Sana = dnX + 0.5_Rkind
   Snum = get_Num_dnSca_FROM_f_x(x,faplusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX+a:         (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX + 0.5')

   Sana = dnX + dnX
   Snum = get_Num_dnSca_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX+dnX=2*dnX: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX + dnX')

   Sana = +(0.5_Rkind - dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,faminusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') '+(a-dnX):      (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='+(0.5 - dnX)')

   Sana = -(dnX - 0.5_Rkind)
   Snum = get_Num_dnSca_FROM_f_x(x,faminusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') '-(dnX-a):      (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='-(dnX - 0.5)')

   Sana = dnX - dnX
   write(out_unitp,'(a,l2)') 'dnX-dnX                   ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX - dnX')

   Sana = 2._Rkind * dnX
   Snum = get_Num_dnSca_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'a*dnX:         (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='2. * dnX')

   Sana =  dnX * 2._Rkind
   Snum = get_Num_dnSca_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX*a:         (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX * 2.')

   Sana =  dnX * dnX - dnX**(2._Rkind)
   write(out_unitp,'(a,l2)') 'dnX*dnX-dnX**2.           ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX * dnX - dnX**(2.)')

   Sana =  dnX / 0.5_Rkind -dnX*2._Rkind
   write(out_unitp,'(a,l2)') 'dnX/0.5-dnX*2             ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
   !CALL Write_dnSca(Sana,out_unitp,info='dnX/0.5-dnX*2')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX/0.5 -dnX*2.')

   Sana =  0.5_Rkind / dnX - 0.5_Rkind*dnX**(-1)
   write(out_unitp,'(a,l2)') '0.5/dnX-0.5.dnX*^-1       ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='0.5 / dnX - 0.5*dnX**(-1)')

   Sana =  dnX / dnX - 1._Rkind
   write(out_unitp,'(a,l2)') 'dnX/dnX-1.                ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX / dnX - 1.')
   !test 13
   Sana =  dnX**0.5_Rkind - sqrt(dnX)
   write(out_unitp,'(a,l2)') 'dnX**0.5-sqrt(dnX)        ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX**0.5-sqrt(dnX)')

!   Sana =  dnX**0 - ONE
!   write(out_unitp,'(a,l2)') 'dnX**0-ONE                ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
!   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX**0 - ONE')
!
!   Sana =  dnX**1 - dnX
!   write(out_unitp,'(a,l2)') 'dnX**1-dnX                ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
!   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX**1 - dnX')
!
!   Sana =  dnX**2 - dnX*dnX
!   write(out_unitp,'(a,l2)') 'dnX**2-dnX*dnX            ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
!   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX**2 - dnX*dnX')
   !test 14+2
   Sana =  dnX**3 - dnX*dnX*dnX
   write(out_unitp,'(a,l2)') 'dnX**3-dnX*dnX*dnX        ==0?',Check_dnSca_IS_ZERO(Sana,dnSerr_test)
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='dnX**3 - dnX*dnX*dnX')

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "functions: sqrt, exp, log, ... sin, asin, ... acosh ..."
   write(out_unitp,'(a)') "============================================"
   !test 17
   Sana = sqrt(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,SQRT_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sqrt: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on sqrt')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='sqrt(dnX)')

   Sana = abs(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,ABS_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'abs:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on abs')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='abs(dnX)')

   Sana = exp(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,EXP_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'exp:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on exp')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='exp(dnX)')

   Sana = log(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,log_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'log:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on log')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='log(dnX)')

   Sana = log10(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,log10_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'log10: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on log10')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='log10(dnX)')

   Sana = sin(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,sin_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sin:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on sin')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='sin(dnX)')

   Sana = asin(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,asin_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'asin: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on asin')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='asin(dnX)')

   Sana = cos(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,cos_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'cos:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on cos')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='cos(dnX)')

   Sana = acos(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,acos_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'acos: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on acos')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='acos(dnX)')

   Sana = tan(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,tan_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'tan:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on tan')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='tan(dnX)')

   Sana = atan(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,atan_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'atan: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on atan')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='atan(dnX)')

   Sana = sinh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,sinh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sinh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on sinh')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='sinh(dnX)')

   Sana = asinh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,asinh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'asinh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on asinh')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='asinh(dnX)')

   Sana = cosh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,cosh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'cosh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on cosh')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='cosh(dnX)')

   dnY = init_dnSca(FOUR*x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
   Sana = acosh(dnY)
   Snum = get_Num_dnSca_FROM_f_x(FOUR*x,acosh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'acosh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on acosh')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='acosh(dn4X)')

   Sana = tanh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,tanh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'tanh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on tanh')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='tanh(dnX)')
   !test 33
   Sana = atanh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,atanh_perso,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'atanh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,dnSerr_test)
   CALL write_dnSca(Sana,info='test on atanh')
   CALL Write_dnSca_FOR_test(Sana,nio_test,info='atanh(dnX)')

   close(unit=nio_test)


END PROGRAM TEST_dnSca
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
FUNCTION SQRT_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = sqrt(x)

END FUNCTION SQRT_perso
FUNCTION ABS_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = ABS(x)

END FUNCTION ABS_perso
FUNCTION EXP_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = EXP(x)

END FUNCTION EXP_perso
FUNCTION LOG_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = LOG(x)

END FUNCTION LOG_perso
FUNCTION LOG10_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = LOG10(x)

END FUNCTION LOG10_perso
FUNCTION SIN_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = SIN(x)

END FUNCTION SIN_perso
FUNCTION ASIN_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = ASIN(x)

END FUNCTION ASIN_perso
FUNCTION COS_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = COS(x)

END FUNCTION COS_perso
FUNCTION ACOS_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = ACOS(x)

END FUNCTION ACOS_perso
FUNCTION TAN_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = TAN(x)

END FUNCTION TAN_perso
FUNCTION ATAN_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = ATAN(x)

END FUNCTION ATAN_perso
FUNCTION SINH_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = SINH(x)

END FUNCTION SINH_perso
FUNCTION ASINH_perso(x) RESULT(f)
  USE mod_NumParameters
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
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = COSH(x)

END FUNCTION COSH_perso
FUNCTION ACOSH_perso(x) RESULT(f)
  USE mod_NumParameters
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
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

  f = TANH(x)

END FUNCTION TANH_perso
FUNCTION ATANH_perso(x) RESULT(f)
  USE mod_NumParameters
  IMPLICIT NONE
  real (kind=Rkind) :: f
  real (kind=Rkind), intent(in) :: x

#if __INVHYP == 1
    f = atanh(x)
#else
    f = HALF*log((ONE+x)/(ONE-x))
#endif

END FUNCTION ATANH_perso
