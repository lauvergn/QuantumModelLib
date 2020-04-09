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
!    Copyright 2017 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================
SUBROUTINE sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
  USE mod_Lib
  USE mod_QModel
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  logical,                intent(in)        :: adiabatic

!$OMP CRITICAL (CRIT_sub_Init_Qmodel)
    !$ write(out_unitp,*) 'begining: max threads?',omp_get_max_threads()
    write(out_unitp,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unitp,*) 'ndim,nsurf     : ',ndim,nsurf
    flush(out_unitp)

    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),      &
                    ndim=ndim,nsurf=nsurf,adiabatic=adiabatic,option=option)

    CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Init_Qmodel in Model_driver.f90')

    IF (ndim /= QuantumModel%ndim .OR. nsurf /= QuantumModel%nsurf) THEN
      write(out_unitp,*) ' ERROR in sub_Init_Qmodel'
      write(out_unitp,*) ' ndim, nsurf :',ndim,nsurf
      write(out_unitp,*) ' The ndim or nsurf values are incompatible ...'
      write(out_unitp,*) '   .... with the intialized ones !!'
      write(out_unitp,*) ' model name                  : ',QuantumModel%QM%pot_name
      write(out_unitp,*) ' ndim, nsurf (from the model): ',QuantumModel%ndim,QuantumModel%nsurf
      write(out_unitp,*) '   CHECK your data !!'
      STOP 'ERROR in sub_Init_Qmodel: wrong ndim or nsurf'
    END IF
!$OMP END CRITICAL (CRIT_sub_Init_Qmodel)

END SUBROUTINE sub_Init_Qmodel
SUBROUTINE sub_Write_Qmodel(nio)
  USE mod_QModel
  IMPLICIT NONE

  integer,            intent(in)    :: nio

  CALL Write_Model(QuantumModel,nio)

END SUBROUTINE sub_Write_Qmodel
SUBROUTINE get_Qmodel_Q0(Q0,option)
  USE mod_Lib
  USE mod_QModel
  IMPLICIT NONE

  real (kind=Rkind),      intent(inout)     :: Q0(QuantumModel%ndim)
  integer,                intent(in)        :: option

  CALL get_Q0_Model(Q0,QuantumModel,option)

END SUBROUTINE get_Qmodel_Q0
SUBROUTINE sub_Qmodel_V(V,Q)
  USE mod_Lib
  USE mod_QModel
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)

  CALL calc_pot(V,QuantumModel,Q)

END SUBROUTINE sub_Qmodel_V
SUBROUTINE sub_Qmodel_VVec(V,Vec,Q)
  USE mod_Lib
  USE mod_dnMat
  USE mod_QModel
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vec(QuantumModel%nsurf,QuantumModel%nsurf)

  TYPE (dnMat_t)                  :: Vec_loc,PotVal_loc



  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VVec in Model_driver.f90')


  CALL QML_alloc_dnMat(Vec_loc,nsurf=QuantumModel%nsurf,ndim=QuantumModel%ndim,nderiv=0)
  CALL QML_alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,ndim=QuantumModel%ndim,nderiv=0)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=0,Vec=Vec_loc)

  V(:,:)   = PotVal_loc%d0
  Vec(:,:) = Vec_loc%d0

  CALL QML_dealloc_dnMat(PotVal_loc)
  CALL QML_dealloc_dnMat(Vec_loc)

END SUBROUTINE sub_Qmodel_VVec
SUBROUTINE sub_Qmodel_VG(V,G,Q)
  USE mod_Lib
  USE mod_QModel
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     ::                          &
                                V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     ::                          &
               g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)

  CALL calc_pot_grad(V,G,QuantumModel,Q)

END SUBROUTINE sub_Qmodel_VG
SUBROUTINE sub_Qmodel_VG_NAC(V,G,NAC,Q)
  USE mod_Lib
  USE mod_QModel
  USE mod_dnMat
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


  CALL QML_alloc_dnMat(NAC_loc,nsurf=QuantumModel%nsurf,ndim=QuantumModel%ndim,nderiv=1)
  CALL QML_alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,ndim=QuantumModel%ndim,nderiv=1)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=1,NAC=NAC_loc)

  V(:,:)     = PotVal_loc%d0
  G(:,:,:)   = PotVal_loc%d1
  NAC(:,:,:) = NAC_loc%d1

  CALL QML_dealloc_dnMat(PotVal_loc)
  CALL QML_dealloc_dnMat(NAC_loc)

END SUBROUTINE sub_Qmodel_VG_NAC
SUBROUTINE sub_Qmodel_VGH(V,G,H,Q)
  USE mod_Lib
  USE mod_QModel
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     ::                          &
                                V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     ::                          &
              g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: h(QuantumModel%nsurf,    &
                 QuantumModel%nsurf,QuantumModel%ndim,QuantumModel%ndim)

  CALL calc_pot_grad_hess(V,G,H,QuantumModel,Q)

END SUBROUTINE sub_Qmodel_VGH
SUBROUTINE sub_Qmodel_Check_anaVSnum(Q,nderiv)
  USE mod_Lib
  USE mod_QModel
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  integer,                intent(in)        :: nderiv

  CALL Check_analytical_numerical_derivatives(QuantumModel,Q,nderiv=2)

END SUBROUTINE sub_Qmodel_Check_anaVSnum
SUBROUTINE get_Qmodel_GGdef(GGdef)
  USE mod_Lib
  USE mod_QModel
  IMPLICIT NONE

  real (kind=Rkind),      intent(inout)     :: GGdef(QuantumModel%ndim,QuantumModel%ndim)


  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_GGdef in Model_driver.f90')

  GGdef(:,:) = QuantumModel%QM%d0GGdef

END SUBROUTINE get_Qmodel_GGdef
SUBROUTINE set_Qmodel_GGdef(GGdef,ndim)
  USE mod_Lib
  USE mod_QModel
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
SUBROUTINE sub_model_V(V,Q,ndim,nsurf)
  USE mod_Lib
  USE mod_QModel
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
  USE mod_Lib
  USE mod_QModel
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
  USE mod_Lib
  USE mod_QModel
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
  USE mod_Lib
  USE mod_QModel
  USE mod_dnMat
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

  CALL QML_alloc_dnMat(PotVal_loc,nsurf=QuantumModel%nsurf,ndim=QuantumModel%ndim,nderiv=1)

  CALL Eval_Pot(QuantumModel,Q,PotVal_loc,nderiv=1,NAC=NAC_loc)

  V(:,:)     = PotVal_loc%d0
  G(:,:,:)   = PotVal_loc%d1
  NAC(:,:,:) = NAC_loc%d1

  CALL QML_dealloc_dnMat(PotVal_loc)
  CALL QML_dealloc_dnMat(NAC_loc)


END SUBROUTINE sub_model1_VG_NAC
SUBROUTINE sub_model1_VGH(V,G,H,Q,ndim,nsurf,pot_name,option)
  USE mod_Lib
  USE mod_QModel
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
  USE mod_Lib
  USE mod_QModel
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
  USE mod_Lib
  USE mod_QModel
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
  USE mod_Lib
  USE mod_QModel
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
