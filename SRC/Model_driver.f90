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
SUBROUTINE sub_model_V(V,Q,ndim,nsurf)
  USE mod_NumParameters
  USE mod_dnMatPot
  USE mod_Model
  USE mod_Lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  !real (kind=Rkind),      intent(inout)     :: g(nsurf,nsurf,ndim)
  !real (kind=Rkind),      intent(inout)     :: h(nsurf,nsurf,ndim,ndim)

  TYPE (Param_Model), SAVE         :: Para_Model
  logical, SAVE                    :: begin = .TRUE.


!$OMP CRITICAL (CRIT_sub_model_V)
  IF (begin) THEN
    !$ write(6,*) 'begining: max threads?',begin,omp_get_max_threads()
    CALL Init_Model(Para_Model,pot_name='HenonHeiles',ndim=ndim,read_param=.FALSE.)

    IF (ndim /= Para_Model%ndim .OR. nsurf /= Para_Model%nsurf) THEN
      write(6,*) ' ERROR in sub_model_V'
      write(6,*) ' ndim, nsurf :',ndim,nsurf
      write(6,*) ' The ndim or nsurf values are incompatible ...'
      write(6,*) '   .... with the intialized ones !!'
      write(6,*) ' model name                  : ',Para_Model%pot_name
      write(6,*) ' ndim, nsurf (from the model): ',Para_Model%ndim,Para_Model%nsurf
      write(6,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model_V: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model_V)

  CALL calc_pot(V,Para_Model,Q)


END SUBROUTINE sub_model_V
SUBROUTINE sub_model1_V(V,Q,ndim,nsurf,pot_name,option)
  USE mod_NumParameters
  USE mod_dnMatPot
  USE mod_Model
  USE mod_Lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=16),     intent(in)        :: pot_name
  integer,                intent(in)        :: option
  !real (kind=Rkind),      intent(inout)     :: g(nsurf,nsurf,ndim)
  !real (kind=Rkind),      intent(inout)     :: h(nsurf,nsurf,ndim,ndim)

  TYPE (Param_Model), SAVE         :: Para_Model
  logical, SAVE                    :: begin = .TRUE.


!$OMP CRITICAL (CRIT_sub_model_V)
  IF (begin .OR. option < 0) THEN
    !$ write(6,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(6,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(6,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(6)
    CALL Init_Model(Para_Model,pot_name=trim(adjustl(pot_name)),ndim=ndim,read_param=.FALSE.,option=option)

    IF (ndim /= Para_Model%ndim .OR. nsurf /= Para_Model%nsurf) THEN
      write(6,*) ' ERROR in sub_model_V'
      write(6,*) ' ndim, nsurf :',ndim,nsurf
      write(6,*) ' The ndim or nsurf values are incompatible ...'
      write(6,*) '   .... with the intialized ones !!'
      write(6,*) ' model name                  : ',Para_Model%pot_name
      write(6,*) ' ndim, nsurf (from the model): ',Para_Model%ndim,Para_Model%nsurf
      write(6,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model_V: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model_V)

  CALL calc_pot(V,Para_Model,Q)


END SUBROUTINE sub_model1_V
SUBROUTINE sub_model1_VG(V,G,Q,ndim,nsurf,pot_name,option)
  USE mod_NumParameters
  USE mod_dnMatPot
  USE mod_Model
  USE mod_Lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=16),     intent(in)        :: pot_name
  integer,                intent(in)        :: option
  real (kind=Rkind),      intent(inout)     :: G(nsurf,nsurf,ndim)
  !real (kind=Rkind),      intent(inout)     :: h(nsurf,nsurf,ndim,ndim)

  TYPE (Param_Model), SAVE         :: Para_Model
  logical, SAVE                    :: begin = .TRUE.


!$OMP CRITICAL (CRIT_sub_model_V)
  IF (begin) THEN
    !$ write(6,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(6,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(6,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(6)
    CALL Init_Model(Para_Model,pot_name=trim(adjustl(pot_name)),ndim=ndim,read_param=.FALSE.,option=option)

    IF (ndim /= Para_Model%ndim .OR. nsurf /= Para_Model%nsurf) THEN
      write(6,*) ' ERROR in sub_model_V'
      write(6,*) ' ndim, nsurf :',ndim,nsurf
      write(6,*) ' The ndim or nsurf values are incompatible ...'
      write(6,*) '   .... with the intialized ones !!'
      write(6,*) ' model name                  : ',Para_Model%pot_name
      write(6,*) ' ndim, nsurf (from the model): ',Para_Model%ndim,Para_Model%nsurf
      write(6,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model_V: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model_V)

  CALL calc_pot_grad(V,G,Para_Model,Q)


END SUBROUTINE sub_model1_VG
SUBROUTINE sub_model1_VGH(V,G,H,Q,ndim,nsurf,pot_name,option)
  USE mod_NumParameters
  USE mod_dnMatPot
  USE mod_Model
  USE mod_Lib
  IMPLICIT NONE

  integer,                intent(in)        :: ndim,nsurf
  real (kind=Rkind),      intent(in)        :: Q(ndim)
  real (kind=Rkind),      intent(inout)     :: V(nsurf,nsurf)
  character (len=16),     intent(in)        :: pot_name
  integer,                intent(in)        :: option
  real (kind=Rkind),      intent(inout)     :: G(nsurf,nsurf,ndim)
  real (kind=Rkind),      intent(inout)     :: H(nsurf,nsurf,ndim,ndim)

  TYPE (Param_Model), SAVE         :: Para_Model
  logical, SAVE                    :: begin = .TRUE.


!$OMP CRITICAL (CRIT_sub_model_V)
  IF (begin) THEN
    !$ write(6,*) 'begining: max threads?',begin,omp_get_max_threads()
    write(6,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(6,*) 'ndim,nsurf,    : ',ndim,nsurf
    flush(6)
    CALL Init_Model(Para_Model,pot_name=trim(adjustl(pot_name)),ndim=ndim,read_param=.FALSE.,option=option)

    IF (ndim /= Para_Model%ndim .OR. nsurf /= Para_Model%nsurf) THEN
      write(6,*) ' ERROR in sub_model_V'
      write(6,*) ' ndim, nsurf :',ndim,nsurf
      write(6,*) ' The ndim or nsurf values are incompatible ...'
      write(6,*) '   .... with the intialized ones !!'
      write(6,*) ' model name                  : ',Para_Model%pot_name
      write(6,*) ' ndim, nsurf (from the model): ',Para_Model%ndim,Para_Model%nsurf
      write(6,*) '   CHECK your data !!'
      STOP 'ERROR in sub_model_V: wrong ndim or nsurf'
    END IF
    begin = .FALSE.
  END IF
!$OMP END CRITICAL (CRIT_sub_model_V)

  CALL calc_pot_grad_hess(V,G,H,Para_Model,Q)


END SUBROUTINE sub_model1_VGH
