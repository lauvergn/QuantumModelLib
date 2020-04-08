program prog
  USE mod_EmptyModel
  USE mod_TemplateModel
  USE mod_MorseModel
  USE mod_dnS
  USE mod_dnMat
  USE mod_QModel
  IMPLICIT NONE


  real(kind=Rkind),  allocatable    :: Q(:)
  TYPE (dnS_t),      allocatable    :: dnQ(:)
  TYPE (dnS_t),      allocatable    :: Mat_OF_PotDia(:,:)
    TYPE (dnMat_t)                  :: PotVal

  TYPE(QModel_t)                        :: QModel

  integer :: i,nderiv

  nderiv     = 2

  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  CALL Init_Model(QModel,pot_name='read_model')
  !CALL Write_Model(QModel,nio=out_unitp)

  write(out_unitp,*) 'Gdef',QModel%QM%get_d0GGdef_QModel()

  allocate(Q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)
  IF (allocated(Q)) write(out_unitp,*) 'Q0',Q
  flush(out_unitp)

  Q(:) = [ 1.2_Rkind,0.2_Rkind ]
  CALL Eval_Pot(QModel,Q,PotVal,nderiv)
  write(out_unitp,*) 'PotVal'
  CALL QML_Write_dnMat(PotVal,6)

  CALL QML_dealloc_dnMat(PotVal)

  write(out_unitp,*) "---------------------------------------------------"
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  deallocate(QModel%QM)
  deallocate(Q)
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  write(out_unitp,*)
  flush(out_unitp)
stop
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  CALL Init_Model(QModel,pot_name='morse')
  !CALL Write_Model(QModel,nio=out_unitp)

  write(out_unitp,*) 'Gdef',QModel%QM%get_d0GGdef_QModel()

  allocate(Q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)
  IF (allocated(Q)) write(out_unitp,*) 'Q0',Q

  CALL Eval_Pot(QModel,Q,PotVal,nderiv)
  write(out_unitp,*) 'PotVal'
  CALL QML_Write_dnMat(PotVal,6)

  CALL QML_dealloc_dnMat(PotVal)
  deallocate(QModel%QM)
  deallocate(Q)
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  write(out_unitp,*)
  flush(out_unitp)

  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  CALL Init_Model(QModel,pot_name='template')
  !CALL Write_Model(QModel,nio=out_unitp)

  write(out_unitp,*) 'Gdef',QModel%QM%get_d0GGdef_QModel()

  allocate(Q(QModel%QM%ndim))
  CALL get_Q0_Model(Q,QModel,option=0)
  IF (allocated(Q)) write(out_unitp,*) 'Q0',Q

  CALL Eval_Pot(QModel,Q,PotVal,nderiv)
  write(out_unitp,*) 'PotVal'
  CALL QML_Write_dnMat(PotVal,6)

  CALL QML_dealloc_dnMat(PotVal)

  write(out_unitp,*) "---------------------------------------------------"
  CALL Check_analytical_numerical_derivatives(QModel,Q,nderiv)

  deallocate(QModel%QM)
  deallocate(Q)
  write(out_unitp,*) "===================================================="
  write(out_unitp,*) "===================================================="
  write(out_unitp,*)
  flush(out_unitp)


end program prog
