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
SUBROUTINE sub_Read_Qmodel(ndim,nsurf,nio)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(inout)     :: ndim,nsurf
  integer,                intent(in)        :: nio

!$OMP CRITICAL (CRIT_sub_Read_Qmodel)
    !$ write(out_unitp,*) 'begining: max threads?',omp_get_max_threads()
    write(out_unitp,*) 'ndim,nsurf     : ',ndim,nsurf
    write(out_unitp,*) 'nio            : ',nio
    flush(out_unitp)

    CALL Init_Model(QuantumModel,ndim=ndim,nsurf=nsurf,                         &
                    read_param=.TRUE.,nio_param_file=nio)

    CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Init_Qmodel in Model_driver.f90')

    ndim  = QuantumModel%ndim
    nsurf = QuantumModel%nsurf

!$OMP END CRITICAL (CRIT_sub_Read_Qmodel)

END SUBROUTINE sub_Read_Qmodel
SUBROUTINE sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(inout)     :: ndim,nsurf
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  logical,                intent(in)        :: adiabatic

!$OMP CRITICAL (CRIT_sub_Init_Qmodel)
    !$ write(out_unitp,*) 'begining: max threads?',omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf     : ',ndim,nsurf
    flush(out_unitp)

    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),              &
                    ndim=ndim,nsurf=nsurf,adiabatic=adiabatic,option=option)


    CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Init_Qmodel in Model_driver.f90')

    ndim  = QuantumModel%ndim
    nsurf = QuantumModel%nsurf

!$OMP END CRITICAL (CRIT_sub_Init_Qmodel)

END SUBROUTINE sub_Init_Qmodel
SUBROUTINE sub_Init_Qmodel_Cart(ndim,nsurf,pot_name,adiabatic,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(inout)     :: ndim,nsurf
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  logical,                intent(in)        :: adiabatic

!$OMP CRITICAL (CRIT_sub_Init_Qmodel_Cart)
    !$ write(out_unitp,*) 'begining: max threads?',omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf     : ',ndim,nsurf
    flush(out_unitp)

    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),              &
                    ndim=ndim,nsurf=nsurf,adiabatic=adiabatic,Cart_TO_Q=.TRUE., &
                    option=option)

    CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Init_Qmodel_Cart in Model_driver.f90')

    ndim  = QuantumModel%ndim
    nsurf = QuantumModel%nsurf

!$OMP END CRITICAL (CRIT_sub_Init_Qmodel_Cart)

END SUBROUTINE sub_Init_Qmodel_Cart

SUBROUTINE sub_check_Init_Qmodel(check)
  USE Model_m
  IMPLICIT NONE

  logical,            intent(inout)    :: check

  check = check_Init_QModel(QuantumModel)

END SUBROUTINE sub_check_Init_Qmodel

FUNCTION get_Qmodel_ndim() RESULT(ndim)
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  integer     :: ndim

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_ndim in Model_driver.f90')

  ndim  = QuantumModel%ndim

END FUNCTION get_Qmodel_ndim
FUNCTION get_Qmodel_nsurf() RESULT(nsurf)
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  integer     :: nsurf

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_ndim in Model_driver.f90')

  nsurf  = QuantumModel%nsurf

END FUNCTION get_Qmodel_nsurf
FUNCTION get_Qmodel_Vib_Adia() RESULT(Vib_Adia)
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  logical     :: Vib_Adia

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_VibAdia in Model_driver.f90')

  Vib_Adia  = QuantumModel%QM%Vib_adia

END FUNCTION get_Qmodel_Vib_Adia
SUBROUTINE sub_Write_Qmodel(nio)
  USE Model_m
  IMPLICIT NONE

  integer,            intent(in)    :: nio

  CALL Write_Model(QuantumModel,nio)

END SUBROUTINE sub_Write_Qmodel
SUBROUTINE sub_Qmodel_Opt(Q,i_surf,nb_deg,icv,Max_it)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(inout)     :: Q(QuantumModel%ndim) ! intial and final geometry
  integer,                intent(in)        :: i_surf,nb_deg,icv,Max_it

  !local variables
  TYPE (QML_Opt_t)               :: Opt_param
  real (kind=Rkind)              :: Q0(QuantumModel%ndim)


  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_Opt in Model_driver.f90')


  Q0(:) = Q

  Opt_param%i_surf = i_surf
  IF (Opt_param%i_surf < 1) Opt_param%i_surf = 1

  Opt_param%nb_neg = nb_deg

  Opt_param%max_it = Max_it

  CALL Init_QML_Opt(Opt_param,QuantumModel,icv=icv,read_param=.FALSE.)

  CALL QML_Opt(Q,QuantumModel,Opt_param,Q0=Q0)


END SUBROUTINE sub_Qmodel_Opt
! get a reference geometry, Q0
! the option parameter is not used
SUBROUTINE get_Qmodel_Q0(Q0,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(inout)     :: Q0(QuantumModel%ndim)
  integer,                intent(in)        :: option

  CALL get_Q0_Model(Q0,QuantumModel,option)

END SUBROUTINE get_Qmodel_Q0
SUBROUTINE sub_Qmodel_V(V,Q)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE


  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)

  CALL calc_pot(V,QuantumModel,Q(1:QuantumModel%ndim))

END SUBROUTINE sub_Qmodel_V
SUBROUTINE sub_Qmodel_VVec(V,Vec,Q)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vec(QuantumModel%nsurf,QuantumModel%nsurf)

  TYPE (dnMat_t)                  :: Vec_loc,PotVal_loc


  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VVec in Model_driver.f90')


  CALL alloc_dnMat(Vec_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=0)
  CALL alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=0)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=0,Vec=Vec_loc)

  V(:,:)   = PotVal_loc%d0
  Vec(:,:) = Vec_loc%d0

  CALL dealloc_dnMat(PotVal_loc)
  CALL dealloc_dnMat(Vec_loc)

END SUBROUTINE sub_Qmodel_VVec
SUBROUTINE sub_Qmodel_VG(V,G,Q)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     ::                          &
                                V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     ::                          &
               g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)

  CALL calc_pot_grad(V,G,QuantumModel,Q)

END SUBROUTINE sub_Qmodel_VG
SUBROUTINE sub_Qmodel_VG_NAC(V,G,NAC,Q)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)       :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    ::                           &
                                V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)    ::                           &
              g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: NAC(QuantumModel%nsurf,   &
                                   QuantumModel%nsurf,QuantumModel%ndim)


  TYPE (dnMat_t)                  :: NAC_loc,PotVal_loc

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VG_NAC in Model_driver.f90')


  CALL alloc_dnMat(NAC_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=1)
  CALL alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=1)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=1,NAC=NAC_loc)

  V(:,:)     = PotVal_loc%d0
  G(:,:,:)   = PotVal_loc%d1
  NAC(:,:,:) = NAC_loc%d1

  CALL dealloc_dnMat(PotVal_loc)
  CALL dealloc_dnMat(NAC_loc)

END SUBROUTINE sub_Qmodel_VG_NAC
SUBROUTINE sub_Qmodel_VGH(V,G,H,Q)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     ::                                  &
                                        V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     ::                                  &
                      g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: h(QuantumModel%nsurf,            &
                         QuantumModel%nsurf,QuantumModel%ndim,QuantumModel%ndim)

  CALL calc_pot_grad_hess(V,G,H,QuantumModel,Q)

END SUBROUTINE sub_Qmodel_VGH
SUBROUTINE sub_Qmodel_tab_HMatVibAdia(tab_MatH,Q,nb_terms)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  integer,                intent(in)        :: nb_terms
  real (kind=Rkind),      intent(inout)     ::                                  &
                        tab_MatH(QuantumModel%nsurf,QuantumModel%nsurf,nb_terms)

  real (kind=Rkind),  allocatable  :: tab_MatH_loc(:,:,:)


  CALL Eval_tab_HMatVibAdia(QuantumModel,Q,tab_MatH_loc)

  IF ( any(shape(tab_MatH_loc) /= shape(tab_MatH) )) THEN
    write(out_unitp,*) ' ERROR in sub_Qmodel_tab_HMatVibAdia'
    write(out_unitp,*) ' The shapes of tab_MatH and the one from Eval_tab_HMatVibAdia() are different'
    write(out_unitp,*) ' shape(tab_MatH)    ',shape(tab_MatH)
    write(out_unitp,*) ' shape(tab_MatH_loc)',shape(tab_MatH_loc)
    write(out_unitp,*) ' check the fortran!'
    write(out_unitp,*)
    STOP ' ERROR in sub_Qmodel_tab_HMatVibAdia'
  END IF
  tab_MatH(:,:,:) = tab_MatH_loc

END SUBROUTINE sub_Qmodel_tab_HMatVibAdia
SUBROUTINE sub_Qmodel_Check_anaVSnum(Q,nderiv)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  integer,                intent(in)        :: nderiv

  CALL Check_analytical_numerical_derivatives(QuantumModel,Q,nderiv=2)

END SUBROUTINE sub_Qmodel_Check_anaVSnum
SUBROUTINE get_Qmodel_GGdef(GGdef)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(inout)     :: GGdef(QuantumModel%ndim,QuantumModel%ndim)

  real (kind=Rkind),  allocatable  :: GGdef_loc(:,:)

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_GGdef in Model_driver.f90')

  GGdef_loc = QuantumModel%QM%get_d0GGdef_QModel()

  IF (any(shape(GGdef) /= shape(GGdef_loc))) THEN
      write(out_unitp,*) ' ERROR in get_Qmodel_GGdef'
      write(out_unitp,*) ' The shapes of GGdef and the one from get_d0GGdef_QModel() are different'
      write(out_unitp,*) ' shape(GGdef)    ',shape(GGdef)
      write(out_unitp,*) ' shape(GGdef_loc)',shape(GGdef_loc)
      write(out_unitp,*) ' check the fortran!'
      write(out_unitp,*)
      STOP ' ERROR in get_Qmodel_GGdef'
  END IF

  GGdef(:,:) = GGdef_loc(:,:)

END SUBROUTINE get_Qmodel_GGdef
SUBROUTINE set_Qmodel_GGdef(GGdef,ndim)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  integer,                intent(in)        :: ndim
  real (kind=Rkind),      intent(in)        :: GGdef(ndim,ndim)

  CALL check_alloc_QM(QuantumModel,name_sub_in='set_Qmodel_GGdef in Model_driver.f90')

  IF (.NOT. allocated(QuantumModel%QM%d0GGdef) .OR. ndim /= QuantumModel%ndim) THEN
      write(out_unitp,*) ' ERROR in set_Qmodel_GGdef'
      write(out_unitp,*) ' QuantumModel%QM%d0GGdef is not allocated or'
      write(out_unitp,*) ' ndim /= QuantumModel%ndim'
      write(out_unitp,*) ' Probably, the initialization is not done!!'
      write(out_unitp,*) '   => CALL sub_Init_Qmodel(...)'
      write(out_unitp,*) ' check the fortran!'
      write(out_unitp,*)
      STOP ' ERROR in set_Qmodel_GGdef'
  END IF

  IF (allocated(QuantumModel%QM%d0GGdef)) THEN
    deallocate(QuantumModel%QM%d0GGdef)
  END IF
  allocate(QuantumModel%QM%d0GGdef(ndim,ndim))

  QuantumModel%QM%d0GGdef(:,:) = GGdef

END SUBROUTINE set_Qmodel_GGdef
SUBROUTINE set_Qmodel_step(step_in)
  USE QMLLib_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: step_in

  CALL set_step_epsi_Model(step_in=step_in)

END SUBROUTINE set_Qmodel_step
SUBROUTINE set_Qmodel_Print_level(printlevel)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  integer,                intent(in)        :: printlevel

  print_level = printlevel  ! from the module QMLLib_NumParameters_m.f90

END SUBROUTINE set_Qmodel_Print_level
SUBROUTINE set_Qmodel_in_unitp(inunitp)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  integer,                intent(in)        :: inunitp

  in_unitp = inunitp  ! from the module QMLLib_NumParameters_m.f90

END SUBROUTINE set_Qmodel_in_unitp
SUBROUTINE set_Qmodel_out_unitp(outunitp)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  integer,                intent(in)        :: outunitp

  out_unitp = outunitp  ! from the module QMLLib_NumParameters_m.f90

END SUBROUTINE set_Qmodel_out_unitp
SUBROUTINE set_Qmodel_Phase_Following(Phase_Following)
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  logical,                intent(in)        :: Phase_Following

  CALL check_alloc_QM(QuantumModel,name_sub_in='set_Qmodel_Phase_Following in Model_driver.f90')

  QuantumModel%QM%Phase_Following = Phase_Following

END SUBROUTINE set_Qmodel_Phase_Following
SUBROUTINE set_Qmodel_Phase_Checking(Phase_Checking)
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  logical,                intent(in)        :: Phase_Checking

  CALL check_alloc_QM(QuantumModel,name_sub_in='set_Qmodel_Phase_Checking in Model_driver.f90')

  QuantumModel%QM%Phase_Checking = Phase_Checking

END SUBROUTINE set_Qmodel_Phase_Checking
SUBROUTINE sub_model_V(V,Q,ndim,nsurf)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)

  logical, SAVE                    :: begin = .TRUE.


!$OMP CRITICAL (CRIT_sub_model_V)
  IF (begin) THEN
    !$ write(out_unitp,*) 'begining: max threads?',begin,omp_get_max_threads()
    CALL Init_Model(QuantumModel,pot_name='HenonHeiles',ndim=ndim,      &
                    read_param=.FALSE.,adiabatic=.TRUE.)
    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unitp,*) ' ERROR in sub_model_V'
      write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unitp,*) '   .... with the intialized ones !!'
      write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unitp,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model_V: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model_V)

  CALL calc_pot(V,QuantumModel,Q)

END SUBROUTINE sub_model_V
SUBROUTINE sub_model1_V(V,Q,ndim,nsurf,pot_name,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option

  logical, SAVE                    :: begin = .TRUE.


  IF (option < 0) THEN
    begin = .TRUE.
    RETURN
  END IF

!$OMP CRITICAL (CRIT_sub_model1_V)
  IF (begin .OR. option < 0) THEN
    !$ write(out_unitp,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unitp)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
             ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.TRUE.)
    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unitp,*) ' ERROR in sub_model1_V'
      write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unitp,*) '   .... with the intialized ones !!'
      write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unitp,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_V: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_V)

  CALL calc_pot(V,QuantumModel,Q)


END SUBROUTINE sub_model1_V

SUBROUTINE sub_model1_VG(V,G,Q,ndim,nsurf,pot_name,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  real (kind=Rkind),      intent(inout)     :: G(nsurf,nsurf,ndim)

  logical, SAVE                    :: begin = .TRUE.

  IF (option < 0) THEN
    begin = .TRUE.
    RETURN
  END IF

!$OMP CRITICAL (CRIT_sub_model1_VG)
  IF (begin) THEN
    !$ write(out_unitp,*) 'begining: max threads?',begin,omp_get_max_threads()
    QuantumModel%QM%adiabatic = .TRUE.
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unitp)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),ndim=ndim,read_param=.FALSE.,option=option)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unitp,*) ' ERROR in sub_model1_VG'
      write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unitp,*) '   .... with the intialized ones !!'
      write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unitp,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_VG: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_VG)

  CALL calc_pot_grad(V,G,QuantumModel,Q)


END SUBROUTINE sub_model1_VG
SUBROUTINE sub_model1_VG_NAC(V,G,NAC,Q,ndim,nsurf,pot_name,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  real (kind=Rkind),      intent(inout)     :: G(nsurf,nsurf,ndim)
  real (kind=Rkind),      intent(inout)     :: NAC(nsurf,nsurf,ndim)


  logical, SAVE                    :: begin = .TRUE.
  TYPE (dnMat_t)                  :: NAC_loc,PotVal_loc

  IF (option < 0) THEN
    begin = .TRUE.
    RETURN
  END IF

!$OMP CRITICAL (CRIT_sub_model1_VG)
  IF (begin) THEN
    !$ write(out_unitp,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unitp)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
             ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.TRUE.)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unitp,*) ' ERROR in sub_model1_VG'
      write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unitp,*) '   .... with the intialized ones !!'
      write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unitp,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_VG: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_VG)

  CALL alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=1)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=1,NAC=NAC_loc)

  V(:,:)     = PotVal_loc%d0
  G(:,:,:)   = PotVal_loc%d1
  NAC(:,:,:) = NAC_loc%d1

  CALL dealloc_dnMat(PotVal_loc)
  CALL dealloc_dnMat(NAC_loc)


END SUBROUTINE sub_model1_VG_NAC
SUBROUTINE sub_model1_VGH(V,G,H,Q,ndim,nsurf,pot_name,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  real (kind=Rkind),      intent(inout)     :: G(nsurf,nsurf,ndim)
  real (kind=Rkind),      intent(inout)     :: H(nsurf,nsurf,ndim,ndim)

  logical, SAVE                    :: begin = .TRUE.

  IF (option < 0) THEN
    begin = .TRUE.
    RETURN
  END IF

!$OMP CRITICAL (CRIT_sub_model1_VGH)
  IF (begin) THEN
    !$ write(out_unitp,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unitp)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
             ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.TRUE.)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unitp,*) ' ERROR in sub_model1_VGH'
      write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unitp,*) '   .... with the intialized ones !!'
      write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unitp,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_VGH: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_VGH)

  CALL calc_pot_grad_hess(V,G,H,QuantumModel,Q)


END SUBROUTINE sub_model1_VGH
SUBROUTINE sub_model1_DiaV(V,Q,ndim,nsurf,pot_name,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option

  logical, SAVE                    :: begin = .TRUE.

  character (len=:), allocatable   :: pot_name_loc


  IF (option < 0) THEN
    begin = .TRUE.
    RETURN
  END IF

!$OMP CRITICAL (CRIT_sub_model1_DiaV)
  IF (begin .OR. option < 0) THEN
    !$ write(out_unitp,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unitp)

    pot_name_loc = trim(adjustl(pot_name))
    CALL string_uppercase_TO_lowercase(pot_name_loc)
    IF (pot_name_loc == 'read_model') THEN
      CALL Init_Model(QuantumModel,read_param=.TRUE.,adiabatic=.FALSE.)
    ELSE
      CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
                    ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.FALSE.)
      IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
        write(out_unitp,*) ' ERROR in sub_model1_DiaV'
        write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
        write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
        write(out_unitp,*) '   .... with the intialized ones !!'
        write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
        write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
        write(out_unitp,*) '   CHECK your data !!'
        STOP 'ERROR in sub_model1_DiaV: wrong ndim or nsurf'
      END IF
    END IF
    deallocate(pot_name_loc)

    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_DiaV)

  CALL calc_pot(V,QuantumModel,Q)


END SUBROUTINE sub_model1_DiaV
SUBROUTINE sub_model1_DiaVG(V,G,Q,ndim,nsurf,pot_name,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  real (kind=Rkind),      intent(inout)     :: G(nsurf,nsurf,ndim)

  logical, SAVE                    :: begin = .TRUE.


  IF (option < 0) THEN
    begin = .TRUE.
    RETURN
  END IF

!$OMP CRITICAL (CRIT_sub_model1_DiaVG)
  IF (begin) THEN
    !$ write(out_unitp,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unitp)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
                    ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.FALSE.)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unitp,*) ' ERROR in sub_model1_DiaVG'
      write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unitp,*) '   .... with the intialized ones !!'
      write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unitp,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_DiaVG: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_DiaVG)

  CALL calc_pot_grad(V,G,QuantumModel,Q)


END SUBROUTINE sub_model1_DiaVG
SUBROUTINE sub_model1_DiaVGH(V,G,H,Q,ndim,nsurf,pot_name,option)
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  real (kind=Rkind),      intent(inout)     :: G(nsurf,nsurf,ndim)
  real (kind=Rkind),      intent(inout)     :: H(nsurf,nsurf,ndim,ndim)

  logical, SAVE                    :: begin = .TRUE.

  IF (option < 0) THEN
    begin = .TRUE.
    RETURN
  END IF

!$OMP CRITICAL (CRIT_sub_model1_DiaVGH)
  IF (begin) THEN
    !$ write(out_unitp,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unitp)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
                    ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.FALSE.)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unitp,*) ' ERROR in sub_model1_DiaVGH'
      write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unitp,*) '   .... with the intialized ones !!'
      write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unitp,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_DiaVGH: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_DiaVGH)

  CALL calc_pot_grad_hess(V,G,H,QuantumModel,Q)


END SUBROUTINE sub_model1_DiaVGH
SUBROUTINE get_Qmodel_nb_Func_ndimFunc(nb_Func,ndimFunc)
  USE QMLLib_UtilLib_m
  USE Model_m
  IMPLICIT NONE

  integer,      intent(inout)     :: nb_Func,ndimFunc


  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_nb_Func_ndimFunc in Model_driver.f90')

  nb_Func  = QuantumModel%QM%nb_Func
  ndimFunc = QuantumModel%QM%ndimFunc

END SUBROUTINE get_Qmodel_nb_Func_ndimFunc
SUBROUTINE get_Qmodel_IndexesFunc(IndexFunc_Ene,IndexFunc_Qop,IndexFunc_Grad,IndexFunc_Hess)
  USE QMLLib_UtilLib_m
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
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
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
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
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
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
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
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
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
