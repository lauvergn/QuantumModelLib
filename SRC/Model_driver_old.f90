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
SUBROUTINE sub_Qmodel_V_old(V,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE


  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)

  CALL calc_pot(V,QuantumModel,Q(1:QuantumModel%ndim))

END SUBROUTINE sub_Qmodel_V_old
SUBROUTINE sub_Qmodel_Vdia_Vadia_old(Vdia,Vadia,Q)
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m, ONLY : dnMat_t,dealloc_dnMat
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: Vdia(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vadia(QuantumModel%nsurf,QuantumModel%nsurf)

  TYPE (dnMat_t)                         :: PotVal_adia,PotVal_dia

  CALL check_alloc_QM(QuantumModel,'sub_Qmodel_Vdia_Vadia')

  CALL Eval_Pot(QuantumModel,Q,PotVal_adia,nderiv=0,PotVal_dia=PotVal_dia)

  Vadia = PotVal_adia%d0
  Vdia  = PotVal_dia%d0

  CALL dealloc_dnMat(PotVal_dia)
  CALL dealloc_dnMat(PotVal_adia)

END SUBROUTINE sub_Qmodel_Vdia_Vadia_old
SUBROUTINE sub_Qmodel_VVec_old(V,Vec,Q)
  USE QDUtil_NumParameters_m
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

END SUBROUTINE sub_Qmodel_VVec_old
SUBROUTINE sub_Qmodel_VVec_Vec0_old(V,Vec,Vec0,Q)
  USE QDUtil_NumParameters_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vec(QuantumModel%NB,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vec0(QuantumModel%NB,QuantumModel%nsurf)

  TYPE (dnMat_t)                  :: Vec_loc,PotVal_loc,Vec0_loc


  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VVec_Vec0 in Model_driver.f90')


  CALL alloc_dnMat(Vec_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=0)
  CALL alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=0)
  CALL alloc_dnMat(Vec0_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=0)

  Vec0_loc%d0 = Vec0(:,:)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=0,Vec=Vec_loc,Vec0=Vec0_loc)

  V(:,:)    = PotVal_loc%d0
  Vec(:,:)  = Vec_loc%d0
  Vec0(:,:) = Vec0_loc%d0

  CALL dealloc_dnMat(PotVal_loc)
  CALL dealloc_dnMat(Vec_loc)
  CALL dealloc_dnMat(Vec0_loc)

END SUBROUTINE sub_Qmodel_VVec_Vec0_old
SUBROUTINE sub_Qmodel_VG_old(V,G,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     ::                          &
                                V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     ::                          &
               g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)

  CALL calc_pot_grad(V,G,QuantumModel,Q)

END SUBROUTINE sub_Qmodel_VG_old
SUBROUTINE sub_Qmodel_VG_NAC_old(V,G,NAC,Q)
  USE QDUtil_NumParameters_m
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

END SUBROUTINE sub_Qmodel_VG_NAC_old
SUBROUTINE sub_Qmodel_VG_NAC_Vec0_old(V,G,NAC,Vec0,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)       :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)    :: G(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: NAC(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: Vec0(QuantumModel%NB,QuantumModel%NB)

  TYPE (dnMat_t)                  :: NAC_loc,PotVal_loc,Vec0_loc

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VG_NAC_Vec0 in Model_driver.f90')


  CALL alloc_dnMat(NAC_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=1)
  CALL alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=1)
  CALL alloc_dnMat(Vec0_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=0)

  Vec0_loc%d0 = Vec0(:,:)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=1,NAC=NAC_loc,Vec0=Vec0_loc)

  V(:,:)     = PotVal_loc%d0
  G(:,:,:)   = PotVal_loc%d1
  NAC(:,:,:) = NAC_loc%d1
  Vec0(:,:)  = Vec0_loc%d0

  CALL dealloc_dnMat(PotVal_loc)
  CALL dealloc_dnMat(NAC_loc)
  CALL dealloc_dnMat(Vec0_loc)


END SUBROUTINE sub_Qmodel_VG_NAC_Vec0_old
SUBROUTINE sub_Qmodel_VG_NAC_VecVec0_old(V,G,NAC,Vec,Vec0,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE ADdnSVM_m, ONLY : dnMat_t,alloc_dnMat,dealloc_dnMat
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)       :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)    :: G(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: NAC(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: Vec(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)    :: Vec0(QuantumModel%nsurf,QuantumModel%nsurf)

  TYPE (dnMat_t)                  :: NAC_loc,PotVal_loc,Vec_loc,Vec0_loc

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VG_NAC_VecVec0 in Model_driver.f90')


  CALL alloc_dnMat(NAC_loc,   nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=1)
  CALL alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=1)
  CALL alloc_dnMat(Vec0_loc,  nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=0)
  CALL alloc_dnMat(Vec_loc,   nsurf=QuantumModel%nsurf,nVar=QuantumModel%ndim,nderiv=0)

  Vec0_loc%d0 = Vec0(:,:)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=1,NAC=NAC_loc,Vec=Vec_loc,Vec0=Vec0_loc)

  V(:,:)     = PotVal_loc%d0
  G(:,:,:)   = PotVal_loc%d1
  NAC(:,:,:) = NAC_loc%d1
  Vec(:,:)   = Vec_loc%d0
  Vec0(:,:)  = Vec0_loc%d0

  CALL dealloc_dnMat(PotVal_loc)
  CALL dealloc_dnMat(NAC_loc)
  CALL dealloc_dnMat(Vec_loc)
  CALL dealloc_dnMat(Vec0_loc)

END SUBROUTINE sub_Qmodel_VG_NAC_VecVec0_old
SUBROUTINE sub_Qmodel_VGH_old(V,G,H,Q)
  USE QDUtil_NumParameters_m
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

END SUBROUTINE sub_Qmodel_VGH_old

SUBROUTINE sub_model_V(V,Q,ndim,nsurf)
  USE QDUtil_NumParameters_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)

  logical, SAVE                    :: begin = .TRUE.


!$OMP CRITICAL (CRIT_sub_model_V)
  IF (begin) THEN
    !$ write(out_unit,*) 'begining: max threads?',begin,omp_get_max_threads()
    CALL Init_Model(QuantumModel,pot_name='HenonHeiles',ndim=ndim,      &
                    read_param=.FALSE.,adiabatic=.TRUE.)
    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unit,*) ' ERROR in sub_model_V'
      write(out_unit,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unit,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unit,*) '   .... with the intialized ones !!'
      write(out_unit,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unit,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unit,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model_V: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model_V)

  CALL calc_pot(V,QuantumModel,Q)

END SUBROUTINE sub_model_V
SUBROUTINE sub_model1_V(V,Q,ndim,nsurf,pot_name,option)
  USE QDUtil_NumParameters_m
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
    !$ write(out_unit,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unit,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unit)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
             ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.TRUE.)
    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unit,*) ' ERROR in sub_model1_V'
      write(out_unit,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unit,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unit,*) '   .... with the intialized ones !!'
      write(out_unit,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unit,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unit,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_V: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_V)

  CALL calc_pot(V,QuantumModel,Q)


END SUBROUTINE sub_model1_V

SUBROUTINE sub_model1_VG(V,G,Q,ndim,nsurf,pot_name,option)
  USE QDUtil_NumParameters_m
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
    !$ write(out_unit,*) 'begining: max threads?',begin,omp_get_max_threads()
    QuantumModel%QM%adiabatic = .TRUE.
    write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unit,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unit)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),ndim=ndim,read_param=.FALSE.,option=option)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unit,*) ' ERROR in sub_model1_VG'
      write(out_unit,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unit,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unit,*) '   .... with the intialized ones !!'
      write(out_unit,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unit,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unit,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_VG: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_VG)

  CALL calc_pot_grad(V,G,QuantumModel,Q)


END SUBROUTINE sub_model1_VG
SUBROUTINE sub_model1_VG_NAC(V,G,NAC,Q,ndim,nsurf,pot_name,option)
  USE QDUtil_NumParameters_m
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
    !$ write(out_unit,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unit,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unit)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
             ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.TRUE.)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unit,*) ' ERROR in sub_model1_VG'
      write(out_unit,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unit,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unit,*) '   .... with the intialized ones !!'
      write(out_unit,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unit,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unit,*) '   CHECK your data !!'
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
  USE QDUtil_NumParameters_m
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
    !$ write(out_unit,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unit,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unit)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
             ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.TRUE.)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unit,*) ' ERROR in sub_model1_VGH'
      write(out_unit,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unit,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unit,*) '   .... with the intialized ones !!'
      write(out_unit,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unit,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unit,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_VGH: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_VGH)

  CALL calc_pot_grad_hess(V,G,H,QuantumModel,Q)


END SUBROUTINE sub_model1_VGH
SUBROUTINE sub_model1_DiaV(V,Q,ndim,nsurf,pot_name,option)
  USE QDUtil_NumParameters_m
  USE QDUtil_m, ONLY : TO_lowercase
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
    !$ write(out_unit,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unit,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unit)

    pot_name_loc = TO_lowercase(trim(adjustl(pot_name)))
    IF (pot_name_loc == 'read_model') THEN
      CALL Init_Model(QuantumModel,read_param=.TRUE.,adiabatic=.FALSE.)
    ELSE
      CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
                    ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.FALSE.)
      IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
        write(out_unit,*) ' ERROR in sub_model1_DiaV'
        write(out_unit,*) ' ndim, nsurf :',ndim,nsurf
        write(out_unit,*) ' The ndim or nsurf values are incompatible ...'
        write(out_unit,*) '   .... with the intialized ones !!'
        write(out_unit,*) ' model name                  : ',QuantumModel%QM%pot_name
        write(out_unit,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
        write(out_unit,*) '   CHECK your data !!'
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
  USE QDUtil_NumParameters_m
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
    !$ write(out_unit,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unit,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unit)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
                    ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.FALSE.)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unit,*) ' ERROR in sub_model1_DiaVG'
      write(out_unit,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unit,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unit,*) '   .... with the intialized ones !!'
      write(out_unit,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unit,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unit,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_DiaVG: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_DiaVG)

  CALL calc_pot_grad(V,G,QuantumModel,Q)


END SUBROUTINE sub_model1_DiaVG
SUBROUTINE sub_model1_DiaVGH(V,G,H,Q,ndim,nsurf,pot_name,option)
  USE QDUtil_NumParameters_m
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
    !$ write(out_unit,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unit,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(out_unit)
    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
                    ndim=ndim,read_param=.FALSE.,option=option,adiabatic=.FALSE.)

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unit,*) ' ERROR in sub_model1_DiaVGH'
      write(out_unit,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unit,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unit,*) '   .... with the intialized ones !!'
      write(out_unit,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unit,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unit,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model1_DiaVGH: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model1_DiaVGH)

  CALL calc_pot_grad_hess(V,G,H,QuantumModel,Q)


END SUBROUTINE sub_model1_DiaVGH
