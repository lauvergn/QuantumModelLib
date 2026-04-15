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
FUNCTION get_Qmodel_ndim() RESULT(ndim)
  USE Model_m
  IMPLICIT NONE

  integer     :: ndim

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_ndim in Model_driver.f90')

  ndim  = QuantumModel%ndim

END FUNCTION get_Qmodel_ndim
FUNCTION get_Qmodel_nsurf() RESULT(nsurf)
  USE Model_m
  IMPLICIT NONE

  integer     :: nsurf

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_ndim in Model_driver.f90')

  nsurf  = QuantumModel%nsurf

END FUNCTION get_Qmodel_nsurf
FUNCTION get_Qmodel_NB() RESULT(NB)
  USE Model_m
  IMPLICIT NONE

  integer     :: NB

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_ndim in Model_driver.f90')

  NB  = QuantumModel%NB

END FUNCTION get_Qmodel_NB
FUNCTION get_Qmodel_Vib_Adia() RESULT(Vib_Adia)
  USE Model_m
  IMPLICIT NONE

  logical     :: Vib_Adia

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_VibAdia in Model_driver.f90')

  Vib_Adia  = QuantumModel%QM%Vib_adia

END FUNCTION get_Qmodel_Vib_Adia

SUBROUTINE sub_Qmodel_Check_anaVSnum(Q,nderiv)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  integer,                intent(in)        :: nderiv

  CALL Check_analytical_numerical_derivatives(QuantumModel,Q,nderiv=2)

END SUBROUTINE sub_Qmodel_Check_anaVSnum
SUBROUTINE sub_Qmodel_check_alloc_d0GGdef(check)
  USE Model_m
  IMPLICIT NONE

  logical,            intent(inout)    :: check

  check = check_Init_QModel(QuantumModel) ! check if QuantumModel%QM is allocated and initialized
  IF (check) check = check_alloc_d0GGdef(QuantumModel%QM)

END SUBROUTINE sub_Qmodel_check_alloc_d0GGdef

SUBROUTINE set_Qmodel_step(step_in)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: step_in

  CALL set_step_epsi_Model(step_in=step_in)

END SUBROUTINE set_Qmodel_step
SUBROUTINE set_Qmodel_Print_level(printlevel)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  integer,                intent(in)        :: printlevel

  CALL set_print_level(printlevel,force=.TRUE.) ! from them module QDUtil lib

END SUBROUTINE set_Qmodel_Print_level
SUBROUTINE set_Qmodel_in_unit(inunitp)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  integer,                intent(in)        :: inunitp

  in_unit = inunitp  ! from the module QMLLib_NumParameters_m.f90

END SUBROUTINE set_Qmodel_in_unit
SUBROUTINE set_Qmodel_out_unit(outunitp)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  integer,                intent(in)        :: outunitp

  out_unit = outunitp  ! from the module QMLLib_NumParameters_m.f90

END SUBROUTINE set_Qmodel_out_unit
SUBROUTINE set_Qmodel_Phase_Following(Phase_Following)
  USE Model_m
  IMPLICIT NONE

  logical,                intent(in)        :: Phase_Following

  CALL check_alloc_QM(QuantumModel,name_sub_in='set_Qmodel_Phase_Following in Model_driver.f90')

  QuantumModel%QM%Phase_Following = Phase_Following

END SUBROUTINE set_Qmodel_Phase_Following
SUBROUTINE set_Qmodel_Phase_Checking(Phase_Checking)
  USE Model_m
  IMPLICIT NONE

  logical,                intent(in)        :: Phase_Checking

  CALL check_alloc_QM(QuantumModel,name_sub_in='set_Qmodel_Phase_Checking in Model_driver.f90')

  QuantumModel%QM%Phase_Checking = Phase_Checking

END SUBROUTINE set_Qmodel_Phase_Checking
SUBROUTINE set_Qmodel_Print_Vec_Overlap(Print_Vec_Overlap)
  USE Model_m
  USE ADdnSVM_m
  IMPLICIT NONE

  logical,                intent(inout)        :: Print_Vec_Overlap

  IF (Print_Vec_Overlap) print_level_dia_dnMat = 2

END SUBROUTINE set_Qmodel_Print_Vec_Overlap

SUBROUTINE get_Qmodel_nb_Func_ndimFunc(nb_Func,ndimFunc)
  USE Model_m
  IMPLICIT NONE

  integer,      intent(inout)     :: nb_Func,ndimFunc


  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_nb_Func_ndimFunc in Model_driver.f90')

  nb_Func  = QuantumModel%QM%nb_Func
  ndimFunc = QuantumModel%QM%ndimFunc

END SUBROUTINE get_Qmodel_nb_Func_ndimFunc
SUBROUTINE get_Qmodel_IndexesFunc(IndexFunc_Ene,IndexFunc_Qop,IndexFunc_Grad,IndexFunc_Hess)
  USE Model_m
  IMPLICIT NONE

  integer,      intent(inout)     :: IndexFunc_Ene
  integer,      intent(inout)     :: IndexFunc_Qop
  integer,      intent(inout)     :: IndexFunc_Grad,IndexFunc_Hess


  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_IndexesFunc in Model_driver.f90')

  IndexFunc_Ene  = QuantumModel%QM%IndexFunc_Ene
  IndexFunc_Qop  = QuantumModel%QM%IndexFunc_Qop
  IndexFunc_Grad = QuantumModel%QM%IndexFunc_Grad
  IndexFunc_Hess = QuantumModel%QM%IndexFunc_Hess

END SUBROUTINE get_Qmodel_IndexesFunc

SUBROUTINE get_Qmodel_d0Func(d0Func,Q,nb_Func,ndimFunc)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE ADdnSVM_m, ONLY : dnS_t,sub_get_dn,dealloc_dnS
  IMPLICIT NONE

  integer,          intent(in)      :: nb_Func,ndimFunc
  real(kind=Rkind), intent(inout)   :: d0Func(nb_Func)
  real(kind=Rkind), intent(in)      :: Q(ndimFunc)


  TYPE (dnS_t),     allocatable     :: Func(:)
  integer                           :: i


  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_d0Func in Model_driver.f90')

  CALL Eval_Func(QuantumModel,Q,Func,nderiv=0)

  DO i=1,size(Func)
    CALL sub_get_dn(Func(i),d0=d0Func(i))
  END DO

  DO i=1,size(Func)
    CALL dealloc_dnS(Func(i))
  END DO
  deallocate(Func)

END SUBROUTINE get_Qmodel_d0Func
SUBROUTINE get_Qmodel_d0d1Func(d0Func,d1Func,Q,nb_Func,ndimFunc)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE ADdnSVM_m, ONLY : dnS_t,sub_get_dn,dealloc_dnS
  IMPLICIT NONE

  integer,          intent(in)      :: nb_Func,ndimFunc
  real(kind=Rkind), intent(inout)   :: d0Func(nb_Func)
  real(kind=Rkind), intent(inout)   :: d1Func(ndimFunc,nb_Func)
  !real(kind=Rkind), intent(inout)   :: d2Func(ndimFunc,ndimFunc,nb_Func)
  !real(kind=Rkind), intent(inout)   :: d3Func(ndimFunc,ndimFunc,ndimFunc,nb_Func)
  real(kind=Rkind), intent(in)      :: Q(ndimFunc)


  TYPE (dnS_t),     allocatable     :: Func(:)
  integer                           :: i


  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_d0d1Func in Model_driver.f90')

  CALL Eval_Func(QuantumModel,Q,Func,nderiv=1)

  DO i=1,size(Func)
    CALL sub_get_dn(Func(i),d0=d0Func(i),d1=d1Func(:,i))
  END DO

  DO i=1,size(Func)
    CALL dealloc_dnS(Func(i))
  END DO
  deallocate(Func)

END SUBROUTINE get_Qmodel_d0d1Func
SUBROUTINE get_Qmodel_d0d1d2Func(d0Func,d1Func,d2Func,Q,nb_Func,ndimFunc)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE ADdnSVM_m, ONLY : dnS_t,sub_get_dn,dealloc_dnS
  IMPLICIT NONE

  integer,          intent(in)      :: nb_Func,ndimFunc
  real(kind=Rkind), intent(inout)   :: d0Func(nb_Func)
  real(kind=Rkind), intent(inout)   :: d1Func(ndimFunc,nb_Func)
  real(kind=Rkind), intent(inout)   :: d2Func(ndimFunc,ndimFunc,nb_Func)
  real(kind=Rkind), intent(in)      :: Q(ndimFunc)


  TYPE (dnS_t),     allocatable     :: Func(:)
  integer                           :: i


  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_d0d1d2Func in Model_driver.f90')

  CALL Eval_Func(QuantumModel,Q,Func,nderiv=2)

  DO i=1,size(Func)
    CALL sub_get_dn(Func(i),d0=d0Func(i),d1=d1Func(:,i),           &
                                         d2=d2Func(:,:,i))
  END DO

  DO i=1,size(Func)
    CALL dealloc_dnS(Func(i))
  END DO
  deallocate(Func)

END SUBROUTINE get_Qmodel_d0d1d2Func
SUBROUTINE get_Qmodel_d0d1d2d3Func(d0Func,d1Func,d2Func,d3Func,Q,nb_Func,ndimFunc)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE ADdnSVM_m, ONLY : dnS_t,sub_get_dn,dealloc_dnS
  IMPLICIT NONE

  integer,          intent(in)      :: nb_Func,ndimFunc
  real(kind=Rkind), intent(inout)   :: d0Func(nb_Func)
  real(kind=Rkind), intent(inout)   :: d1Func(ndimFunc,nb_Func)
  real(kind=Rkind), intent(inout)   :: d2Func(ndimFunc,ndimFunc,nb_Func)
  real(kind=Rkind), intent(inout)   :: d3Func(ndimFunc,ndimFunc,ndimFunc,nb_Func)
  real(kind=Rkind), intent(in)      :: Q(ndimFunc)


  TYPE (dnS_t),     allocatable     :: Func(:)
  integer                           :: i


  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_d0d1d2d3Func in Model_driver.f90')

  CALL Eval_Func(QuantumModel,Q,Func,nderiv=3)

  DO i=1,size(Func)
    CALL sub_get_dn(Func(i),d0=d0Func(i),d1=d1Func(:,i),           &
                                         d2=d2Func(:,:,i),d3=d3Func(:,:,:,i))
  END DO

  DO i=1,size(Func)
    CALL dealloc_dnS(Func(i))
  END DO
  deallocate(Func)

END SUBROUTINE get_Qmodel_d0d1d2d3Func
SUBROUTINE QML_time_perso(name_sub)
  USE QDUtil_m, ONLY : time_perso
  IMPLICIT NONE

  character (len=*) :: name_sub


  !$OMP    CRITICAL (QML_time_perso_CRIT)
  CALL time_perso(name_sub)
  !$OMP   END CRITICAL (QML_time_perso_CRIT)

END SUBROUTINE QML_time_perso
