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
  !< Subroutine which enables to set the model (potential) in reading the 'potential' namelist.
  USE QDUtil_NumParameters_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(inout)     :: ndim  !< number of degree(s) of freedom
  integer,                intent(inout)     :: nsurf !< number of electronic states or vibrational adiabatic channels
  integer,                intent(in)        :: nio   !< unit file where the namelist is read

!$OMP CRITICAL (CRIT_sub_Read_Qmodel)
    !$ write(out_unit,*) 'begining: max threads?',omp_get_max_threads()
    write(out_unit,*) 'ndim,nsurf     : ',ndim,nsurf
    write(out_unit,*) 'nio            : ',nio
    flush(out_unit)

    CALL Init_Model(QuantumModel,ndim=ndim,nsurf=nsurf,                         &
                    read_param=.TRUE.,nio_param_file=nio)

    CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Read_Qmodel in Model_driver.f90')

    ndim  = QuantumModel%ndim
    nsurf = QuantumModel%nsurf

!$OMP END CRITICAL (CRIT_sub_Read_Qmodel)

END SUBROUTINE sub_Read_Qmodel
SUBROUTINE sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
  !< Subroutine which enables to set the model (potential).
  USE QDUtil_NumParameters_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(inout)     :: ndim      !< integer, number of degree(s) of freedom
  integer,                intent(inout)     :: nsurf     !< integer, number of electronic states or vibrational adiabatic channels
  character (len=*),      intent(in)        :: pot_name  !< character string, potential of model name
  integer,                intent(in)        :: option    !< integer, option of the model (relevant to some models)
  logical,                intent(in)        :: adiabatic !< logical, flag to turn on/off the adiatic calculations

!$OMP CRITICAL (CRIT_sub_Init_Qmodel)
    !$ write(out_unit,*) 'begining: max threads?',omp_get_max_threads()
    write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
    write(out_unit,*) 'ndim,nsurf     : ',ndim,nsurf
    flush(out_unit)

    CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),              &
                    ndim=ndim,nsurf=nsurf,adiabatic=adiabatic,option=option)


    CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Init_Qmodel in Model_driver.f90')

    ndim  = QuantumModel%ndim
    nsurf = QuantumModel%nsurf

!$OMP END CRITICAL (CRIT_sub_Init_Qmodel)

END SUBROUTINE sub_Init_Qmodel
SUBROUTINE sub_Init_Qmodel_Cart(ndim,nsurf,pot_name,adiabatic,option)
  USE QDUtil_NumParameters_m
  USE Model_m
  !$ USE omp_lib
  IMPLICIT NONE

  integer,                intent(inout)     :: ndim,nsurf
  character (len=*),      intent(in)        :: pot_name
  integer,                intent(in)        :: option
  logical,                intent(in)        :: adiabatic

!$OMP CRITICAL (CRIT_sub_Init_Qmodel_Cart)
  write(out_unit,*) 'In sub_Init_Qmodel_Cart'
  !$ write(out_unit,*) 'begining: max threads?',omp_get_max_threads()
  write(out_unit,*) 'in sub_Init_Qmodel_Cart'
  write(out_unit,*) 'pot_name,option: ',trim(adjustl(pot_name)),option
  write(out_unit,*) 'ndim,nsurf     : ',ndim,nsurf
  flush(out_unit)

  CALL Init_Model(QuantumModel,pot_name=trim(adjustl(pot_name)),              &
                    ndim=ndim,nsurf=nsurf,adiabatic=adiabatic,Cart_TO_Q=.TRUE., &
                    option=option)

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Init_Qmodel_Cart in Model_driver.f90')

  ndim  = QuantumModel%ndim
  nsurf = QuantumModel%nsurf

  write(out_unit,*) 'ndim,nsurf     : ',ndim,nsurf
  write(out_unit,*) 'End sub_Init_Qmodel_Cart'
  flush(out_unit)

!$OMP END CRITICAL (CRIT_sub_Init_Qmodel_Cart)

END SUBROUTINE sub_Init_Qmodel_Cart

SUBROUTINE sub_check_Init_Qmodel(check)
  USE Model_m
  IMPLICIT NONE

  logical,            intent(inout)    :: check

  check = check_Init_QModel(QuantumModel)

END SUBROUTINE sub_check_Init_Qmodel

SUBROUTINE sub_Write_Qmodel(nio)
  USE Model_m
  IMPLICIT NONE

  integer,            intent(in)    :: nio

  CALL Write_Model(QuantumModel,nio)

END SUBROUTINE sub_Write_Qmodel
SUBROUTINE sub_Qmodel_Opt(Q,i_surf,nb_neg,icv,Max_it)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE Opt_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(inout)     :: Q(QuantumModel%ndim) ! intial and final geometry
  integer,                intent(in)        :: i_surf,nb_neg,icv,Max_it

  !local variables
  TYPE (QML_Opt_t)               :: Opt_param
  real (kind=Rkind)              :: Q0(QuantumModel%ndim)


  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_Opt in Model_driver.f90')


  Q0(:) = Q

  Opt_param%i_surf = i_surf
  IF (Opt_param%i_surf < 1) Opt_param%i_surf = 1

  Opt_param%nb_neg = nb_neg

  Opt_param%max_it = Max_it

  CALL Init_QML_Opt(Opt_param,QuantumModel,icv=icv,read_param=.FALSE.)

  CALL QML_Opt(Q,QuantumModel,Opt_param,Q0=Q0)


END SUBROUTINE sub_Qmodel_Opt
! get a reference geometry, Q0
! the option parameter is not used
SUBROUTINE get_Qmodel_Q0(Q0,option)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(inout)     :: Q0(QuantumModel%ndim)
  integer,                intent(in)        :: option

  CALL get_Q0_Model(Q0,QuantumModel,option)

END SUBROUTINE get_Qmodel_Q0
SUBROUTINE sub_Qmodel_V(V,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE QMLValues_m
  IMPLICIT NONE


  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,'sub_Qmodel_V in Model_driver.f90')

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=0)
  
  IF (allocated(QMLValues%PotAdia%d0)) THEN
    V = QMLValues%PotAdia%d0
  ELSE
    V = QMLValues%PotDia%d0
  END IF

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_V
SUBROUTINE sub_Qmodel_VG(V,G,Q)
  USE QDUtil_NumParameters_m
  USE QMLValues_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     ::                          &
                                V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     ::                          &
               g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,'sub_Qmodel_VG in Model_driver.f90')

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=1)
  
  IF (allocated(QMLValues%PotAdia%d0)) THEN
    V = QMLValues%PotAdia%d0
    g = QMLValues%PotAdia%d1
  ELSE
    V = QMLValues%PotDia%d0
    g = QMLValues%PotDia%d1
  END IF

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_VG
SUBROUTINE sub_Qmodel_VGH(V,G,H,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE QMLValues_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     ::                                  &
                                        V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     ::                                  &
                      g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: h(QuantumModel%nsurf,            &
                         QuantumModel%nsurf,QuantumModel%ndim,QuantumModel%ndim)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,'sub_Qmodel_VGH in Model_driver.f90')

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=2)
  
  IF (allocated(QMLValues%PotAdia%d0)) THEN
    V = QMLValues%PotAdia%d0
    g = QMLValues%PotAdia%d1
    h = QMLValues%PotAdia%d2
  ELSE
    V = QMLValues%PotDia%d0
    g = QMLValues%PotDia%d1
    h = QMLValues%PotDia%d2
  END IF

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_VGH
SUBROUTINE sub_Qmodel_Vdia_Vadia(Vdia,Vadia,Q)
  USE QDUtil_NumParameters_m
  USE QMLValues_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: Vdia(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vadia(QuantumModel%nsurf,QuantumModel%nsurf)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,'sub_Qmodel_Vdia_Vadia in Model_driver.f90')

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=0)

  IF (allocated(QMLValues%PotAdia%d0)) THEN
    Vadia = QMLValues%PotAdia%d0
  ELSE
    Vadia = ZERO
  END IF
  Vdia  = QMLValues%PotDia%d0

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_Vdia_Vadia
SUBROUTINE sub_Qmodel_VVec(V,Vec,Q)
  USE QDUtil_NumParameters_m
  USE QMLValues_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vec(QuantumModel%nsurf,QuantumModel%nsurf)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VVec in Model_driver.f90')

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=0)

  V = QMLValues%PotAdia%d0

  IF (allocated(QMLValues%Vec%d0)) THEN
    Vec = QMLValues%Vec%d0
  ELSE
    Vec = ZERO
  END IF

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_VVec
SUBROUTINE sub_Qmodel_VVec_Vec0(V,Vec,Vec0,Q)
  USE QDUtil_NumParameters_m
  USE QMLValues_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)     :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vec(QuantumModel%NB,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)     :: Vec0(QuantumModel%NB,QuantumModel%nsurf)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VVec_Vec0 in Model_driver.f90')

  QMLValues%Vec0%d0 = Vec0(:,:)

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=0)

  V    = QMLValues%PotAdia%d0
  Vec  = QMLValues%Vec%d0
  Vec0 = QMLValues%Vec0%d0

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_VVec_Vec0

SUBROUTINE sub_Qmodel_VG_NAC(V,G,NAC,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE QMLValues_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)       :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    ::                           &
                                V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)    ::                           &
              g(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: NAC(QuantumModel%nsurf,   &
                                   QuantumModel%nsurf,QuantumModel%ndim)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VG_NAC in Model_driver.f90')

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=1)
  
  V   = QMLValues%PotAdia%d0
  g   = QMLValues%PotAdia%d1
  NAC = QMLValues%NAC%d1

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_VG_NAC
SUBROUTINE sub_Qmodel_VG_NAC_Vec0(V,G,NAC,Vec0,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE QMLValues_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)       :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)    :: G(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: NAC(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: Vec0(QuantumModel%NB,QuantumModel%NB)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VG_NAC_Vec0 in Model_driver.f90')

  QMLValues%Vec0%d0 = Vec0(:,:)

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=1)
  
  V   = QMLValues%PotAdia%d0
  g   = QMLValues%PotAdia%d1
  NAC = QMLValues%NAC%d1

  Vec0 = QMLValues%Vec0%d0

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_VG_NAC_Vec0
SUBROUTINE sub_Qmodel_VG_NAC_VecVec0(V,G,NAC,Vec,Vec0,Q)
  USE QDUtil_NumParameters_m
  USE Model_m
  USE QMLValues_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)       :: Q(QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: V(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)    :: G(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: NAC(QuantumModel%nsurf,QuantumModel%nsurf,QuantumModel%ndim)
  real (kind=Rkind),      intent(inout)    :: Vec(QuantumModel%nsurf,QuantumModel%nsurf)
  real (kind=Rkind),      intent(inout)    :: Vec0(QuantumModel%nsurf,QuantumModel%nsurf)

  TYPE (QMLValues_t)         :: QMLValues

  CALL check_alloc_QM(QuantumModel,name_sub_in='sub_Qmodel_VG_NAC_VecVec0 in Model_driver.f90')

  QMLValues%Vec0%d0 = Vec0(:,:)

  CALL Eval_Pot(QuantumModel,Q,QMLValues,nderiv=1)
  
  V    = QMLValues%PotAdia%d0
  g    = QMLValues%PotAdia%d1
  NAC  = QMLValues%NAC%d1
  Vec  = QMLValues%Vec%d0

  Vec0 = QMLValues%Vec0%d0

  CALL dealloc_QMLValues(QMLValues)

END SUBROUTINE sub_Qmodel_VG_NAC_VecVec0

SUBROUTINE sub_Qmodel_tab_HMatVibAdia(tab_MatH,Q,nb_terms)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(in)        :: Q(QuantumModel%ndim)
  integer,                intent(in)        :: nb_terms
  real (kind=Rkind),      intent(inout)     ::                                  &
                        tab_MatH(QuantumModel%nsurf,QuantumModel%nsurf,nb_terms)

  real (kind=Rkind),  allocatable  :: tab_MatH_loc(:,:,:)


  CALL Eval_tab_HMatVibAdia(QuantumModel,Q,tab_MatH_loc)

  IF ( any(shape(tab_MatH_loc) /= shape(tab_MatH) )) THEN
    write(out_unit,*) ' ERROR in sub_Qmodel_tab_HMatVibAdia'
    write(out_unit,*) ' The shapes of tab_MatH and the one from Eval_tab_HMatVibAdia() are different'
    write(out_unit,*) ' shape(tab_MatH)    ',shape(tab_MatH)
    write(out_unit,*) ' shape(tab_MatH_loc)',shape(tab_MatH_loc)
    write(out_unit,*) ' check the fortran!'
    write(out_unit,*)
    STOP ' ERROR in sub_Qmodel_tab_HMatVibAdia'
  END IF
  tab_MatH(:,:,:) = tab_MatH_loc

END SUBROUTINE sub_Qmodel_tab_HMatVibAdia

SUBROUTINE get_Qmodel_GGdef(GGdef)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  real (kind=Rkind),      intent(inout)     :: GGdef(QuantumModel%ndim,QuantumModel%ndim)

  real (kind=Rkind),  allocatable  :: GGdef_loc(:,:)

  CALL check_alloc_QM(QuantumModel,name_sub_in='get_Qmodel_GGdef in Model_driver.f90')

  GGdef_loc = QuantumModel%QM%get_d0GGdef_QModel()

  IF (any(shape(GGdef) /= shape(GGdef_loc))) THEN
      write(out_unit,*) ' ERROR in get_Qmodel_GGdef'
      write(out_unit,*) ' The shapes of GGdef and the one from get_d0GGdef_QModel() are different'
      write(out_unit,*) ' shape(GGdef)    ',shape(GGdef)
      write(out_unit,*) ' shape(GGdef_loc)',shape(GGdef_loc)
      write(out_unit,*) ' check the fortran!'
      write(out_unit,*)
      STOP ' ERROR in get_Qmodel_GGdef'
  END IF

  GGdef(:,:) = GGdef_loc(:,:)

END SUBROUTINE get_Qmodel_GGdef
SUBROUTINE set_Qmodel_GGdef(GGdef,ndim)
  USE QDUtil_NumParameters_m
  USE Model_m
  IMPLICIT NONE

  integer,                intent(in)        :: ndim
  real (kind=Rkind),      intent(in)        :: GGdef(ndim,ndim)

  CALL check_alloc_QM(QuantumModel,name_sub_in='set_Qmodel_GGdef in Model_driver.f90')

  IF (.NOT. allocated(QuantumModel%QM%d0GGdef) .OR. ndim /= QuantumModel%ndim) THEN
      write(out_unit,*) ' ERROR in set_Qmodel_GGdef'
      write(out_unit,*) ' QuantumModel%QM%d0GGdef is not allocated or'
      write(out_unit,*) ' ndim /= QuantumModel%ndim'
      write(out_unit,*) ' Probably, the initialization is not done!!'
      write(out_unit,*) '   => CALL sub_Init_Qmodel(...)'
      write(out_unit,*) ' check the fortran!'
      write(out_unit,*)
      STOP ' ERROR in set_Qmodel_GGdef'
  END IF

  IF (allocated(QuantumModel%QM%d0GGdef)) THEN
    deallocate(QuantumModel%QM%d0GGdef)
  END IF
  allocate(QuantumModel%QM%d0GGdef(ndim,ndim))

  QuantumModel%QM%d0GGdef(:,:) = GGdef

END SUBROUTINE set_Qmodel_GGdef