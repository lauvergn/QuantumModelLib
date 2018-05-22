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

    TYPE (dnSca)                       :: dnX,dn2X,dnY,Sana,Snum
    real (kind=Rkind)                :: x,err,maxdiff,maxdnSca
    integer                          :: nderiv

    !intrinsic :: sqrt, abs, exp, log, log10, sin, asin, cos, acos, tan, atan, sinh, asinh, cosh, acosh, tanh, atanh
    intrinsic :: dsqrt, dabs, dexp, dlog, dlog10
    intrinsic :: dsin, dasin, dcos, dacos, dtan, datan
    intrinsic :: dsinh, dasinh, dcosh, dacosh, dtanh, datanh
    real (kind=Rkind), external  :: faplusx,faminusx,fatimex,faoverx

    character (len=*), parameter :: name_sub='TEST_dnSca'

   nderiv = 3
   write(out_unitp,'(a,i2)') "== TESTING dnSca module with nderiv=",nderiv


    x       = 0.5_Rkind
    dnX     = init_dnSca(x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
    dn2X    = init_dnSca(x+x,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
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
   Snum = get_Num_dnSca_FROM_f_x(x,faplusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'a+dnX:         (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   Sana = dnX + 0.5_Rkind
   Snum = get_Num_dnSca_FROM_f_x(x,faplusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX+a:         (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   Sana = dnX + dnX
   Snum = get_Num_dnSca_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX+dnX=2*dnX: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)

   Sana = +(0.5_Rkind - dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,faminusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') '+(a-dnX):      (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   Sana = -(dnX - 0.5_Rkind)
   Snum = get_Num_dnSca_FROM_f_x(x,faminusx,nderiv=nderiv)
   write(out_unitp,'(a,l2)') '-(dnX-a):      (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   Sana = dnX - dnX
   write(out_unitp,'(a,l2)') 'dnX-dnX==0                   ?',Check_dnSca_IS_ZERO(Sana,ONETENTH**4)


   Sana = 2._Rkind * dnX
   Snum = get_Num_dnSca_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'a*dnX:         (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   Sana =  dnX * 2._Rkind
   Snum = get_Num_dnSca_FROM_f_x(x,fatimex,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'dnX*a:         (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   Sana =  dnX * dnX - dnX**(2._Rkind)
   write(out_unitp,'(a,l2)') 'dnX*dnX-dnX**2.==0           ?',Check_dnSca_IS_ZERO(Sana,ONETENTH**4)

   Sana =  dnX / 0.5_Rkind -dnX*2._Rkind
   write(out_unitp,'(a,l2)') 'dnX/0.5-dnX*2==0             ?',Check_dnSca_IS_ZERO(Sana,ONETENTH**4)
   Sana =  0.5_Rkind / dnX - 0.5_Rkind*dnX**(-1)
   write(out_unitp,'(a,l2)') '0.5/dnX-0.5.dnX*^-1==0       ?',Check_dnSca_IS_ZERO(Sana,ONETENTH**4)
   Sana =  dnX / dnX - 1._Rkind
   write(out_unitp,'(a,l2)') 'dnX/dnX-1.==0                ?',Check_dnSca_IS_ZERO(Sana,ONETENTH**4)

   Sana =  dnX**0.5_Rkind - sqrt(dnX)
   write(out_unitp,'(a,l2)') 'dnX**0.5-sqrt(dnX)==0        ?',Check_dnSca_IS_ZERO(Sana,ONETENTH**4)
   Sana =  dnX**3 - dnX*dnX*dnX
   write(out_unitp,'(a,l2)') 'dnX**3 - dnX*dnX*dnX==0      ?',Check_dnSca_IS_ZERO(Sana,ONETENTH**4)

   write(out_unitp,'(a)') "============================================"
   write(out_unitp,'(a)') "functions: sqrt, exp, log, ... sin, asin, ... acosh ..."
   write(out_unitp,'(a)') "============================================"

   Sana = sqrt(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dsqrt,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sqrt: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on sqrt')

   Sana = abs(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dabs,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'abs:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on abs')

   Sana = exp(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dexp,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'exp:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on exp')

   Sana = log(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dlog,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'log:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on log')

   Sana = log10(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dlog10,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'log10:(Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on log10')

   Sana = sin(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dsin,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sin:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on sin')

   Sana = asin(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dasin,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'asin: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on asin')

   Sana = cos(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dcos,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'cos:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on cos')

   Sana = acos(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dacos,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'acos: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,TWO*ONETENTH**4)
   CALL write_dnSca(Sana,info='test on acos')

   Sana = tan(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dtan,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'tan:  (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on tan')

   Sana = atan(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,datan,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'atan: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on atan')

   Sana = sinh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dsinh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'sinh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on sinh')

   Sana = asinh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dasinh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'asinh:(Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on asinh')

   Sana = cosh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dcosh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'cosh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on cosh')

   dnY = init_dnSca(FOUR*x  ,ndim=1,nderiv=nderiv,iQ=1) ! to set up the derivatives
   Sana = acosh(dnY)
   Snum = get_Num_dnSca_FROM_f_x(FOUR*x,dacosh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'acosh:(Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,TWO*ONETENTH**4)
   CALL write_dnSca(Sana,info='test on acosh')

   Sana = tanh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,dtanh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'tanh: (Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on tanh')

   Sana = atanh(dnX)
   Snum = get_Num_dnSca_FROM_f_x(x,datanh,nderiv=nderiv)
   write(out_unitp,'(a,l2)') 'atanh:(Sana-Snum)==0?',Check_dnSca_IS_ZERO(Sana-Snum,ONETENTH**4)
   CALL write_dnSca(Sana,info='test on atanh')

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
