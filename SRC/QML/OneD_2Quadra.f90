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
module QML_OneD_2Quadra_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE


  TYPE, EXTENDS (QML_Empty_t) :: QML_OneD_2Quadra_t
    PRIVATE

    ! 1 vibrational mode (R) + 2 electronic surfaces
    real(kind=Rkind)     :: MR    = 20000._Rkind ! mass associated to R
    real(kind=Rkind)     :: k     = 0.020_Rkind ! force constant both diabatic surfaces
    real(kind=Rkind)     :: R1    = 2._Rkind ! equilibrium value for the 1st diabatic surface
    real(kind=Rkind)     :: R2    = 6._Rkind ! equilibrium value for the 2d  diabatic surface
    real(kind=Rkind)     :: a     = 3._Rkind ! parameter for the coupling term
    real(kind=Rkind)     :: b     = 0.01_Rkind ! parameter for the coupling term (b=0 => no-coupling)
    real(kind=Rkind)     :: R3    = 3.875_Rkind ! "equilibrium" value for the coupling term
    real(kind=Rkind)     :: Delta = 0._Rkind ! Energy shift for the second quadratic potential

    ! Dipole moment Mu_12
    real(kind=Rkind)     :: Mu12  = 1.0_Rkind     ! transition dipole moment between the 2 diabatic states
    real(kind=Rkind)     :: Z     = 1.0_Rkind     ! charge

  CONTAINS
    PROCEDURE :: EvalPot_QModel    => EvalPot_QML_OneD_2Quadra
    PROCEDURE :: EvalScalOp_QModel => EvalScalOp_QML_OneD_2Quadra
    PROCEDURE :: Write_QModel      => Write_QML_OneD_2Quadra
    PROCEDURE :: RefValues_QModel  => RefValues_QML_OneD_2Quadra
  END TYPE QML_OneD_2Quadra_t

  PUBLIC :: QML_OneD_2Quadra_t,Init_QML_OneD_2Quadra


contains

  FUNCTION Init_QML_OneD_2Quadra(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m, ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_OneD_2Quadra_t)          :: QModel

    TYPE(QML_Empty_t),    intent(in)   :: QModel_in ! variable to transfer info to the init
    logical,              intent(in)   :: read_param
    integer,              intent(in)   :: nio_param_file

    ! 1 vibrational mode (R) + 2 electronic surfaces
    real(kind=Rkind)     :: MR    ! mass associated to R
    real(kind=Rkind)     :: k     ! force constant both diabatic surfaces
    real(kind=Rkind)     :: R1    ! equilibrium value for the 1st diabatic surface
    real(kind=Rkind)     :: R2    ! equilibrium value for the 2d  diabatic surface
    real(kind=Rkind)     :: a     ! parameter for the coupling term
    real(kind=Rkind)     :: b     ! parameter for the coupling term (b=0 => no-coupling)
    real(kind=Rkind)     :: R3    ! "equilibrium" value for the coupling term
    real(kind=Rkind)     :: Delta ! Energy shift for the second quadratic potential

    real(kind=Rkind)     :: Mu12  ! transition dipole moment between the 2 diabatic states
    real(kind=Rkind)     :: Z     ! charge

    integer :: err_read
    namelist /OneD_2Quadra/ MR,k,R1,R2,a,b,R3,Delta,  Z,Mu12

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_OneD_2Quadra'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) '  read_param:     ',read_param
      write(out_unit,*) '  nio_param_file: ',nio_param_file
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in
    QModel%pot_name    = 'OneD_2Quadra'
    QModel%nsurf       = 2

    IF (QModel%ndim == 0 .OR. QModel%ndim == 1) THEN ! old initialization
      QModel%ndim     = 1
      QModel%Z        = ONE
    ELSE ! ndim /= 1 but read_param=.FALSE.
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) '  wrong ndim and/or read_param'
      write(out_unit,*) '  ndim      ',QModel%ndim
      write(out_unit,*) '  read_param',read_param
      write(out_unit,*) '  ndim /=1 is not possible'
      STOP 'ERROR in Init_QML_OneD_2Quadra: wrong ndim and/or read_param'
    END IF

    IF (read_param) THEN
      MR    = QModel%MR
      k     = QModel%k
      R1    = QModel%R1
      R2    = QModel%R2
      a     = QModel%a
      b     = QModel%b
      R3    = QModel%R3
      Delta = QModel%Delta

      Mu12  = QModel%Mu12
      Z     = QModel%Z

      read(nio_param_file,OneD_2Quadra,IOSTAT=err_read)

      IF (err_read < 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' End-of-file or End-of-record'
        write(out_unit,*) ' The namelist "OneD_2Quadra" is probably absent'
        write(out_unit,*) ' check your data!'
        write(out_unit,*)
        STOP ' ERROR in Init_QML_OneD_2Quadra: The namelist "OneD_2Quadra" is probably absent'
      ELSE IF (err_read > 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Some parameter names of the namelist "OneD_2Quadra" are probaly wrong'
        write(out_unit,*) ' check your data!'
        write(out_unit,nml=OneD_2Quadra)
        STOP ' ERROR in Init_QML_OneD_2Quadra: Some parameter names of the namelist "OneD_2Quadra" are probaly wrong'
      END IF

      QModel%MR    = MR
      QModel%k     = k 
      QModel%R1    = R1
      QModel%R2    = R2
      QModel%a     = a
      QModel%b     = b
      QModel%R3    = R3
      QModel%Delta = Delta

      QModel%Mu12  = Mu12
      QModel%Z     = Z

    END IF


    IF (debug) write(out_unit,*) 'init Q0 of OneD_2Quadra'
    QModel%Q0 = [QModel%R1]


    IF (debug) write(out_unit,*) 'init d0GGdef of OneD_2Quadra'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%MR

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_OneD_2Quadra


  SUBROUTINE EvalPot_QML_OneD_2Quadra(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS (QML_OneD_2Quadra_t), intent(in)     :: QModel
    TYPE (dnS_t),               intent(in)     :: dnQ(:)
    TYPE (dnS_t),               intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,                    intent(in)     :: nderiv

    TYPE (dnS_t)  :: R
    integer       :: i

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_OneD_2Quadra'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      write(out_unit,*) ' QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) ' nderiv:',nderiv
      write(out_unit,*) ' Q(:):',(get_d0(dnQ(i)),i=1,size(dnQ))
      flush(out_unit)
    END IF

    R = dnQ(1)
    !vibrational + electronic part
    Mat_OF_PotDia(1,1) = HALF*QModel%k*(R-QModel%R1)**2
    Mat_OF_PotDia(2,2) = HALF*QModel%k*(R-QModel%R2)**2 + QModel%Delta
    Mat_OF_PotDia(1,2) = QModel%b*exp(-QModel%a*(R-QModel%R3)**2)
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)
 
    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE EvalPot_QML_OneD_2Quadra

  SUBROUTINE EvalScalOp_QML_OneD_2Quadra(QModel,Mat_OF_ScalOpDia,list_Op,dnQ,nderiv)
    USE QDUtil_m,  ONLY : ZERO
    USE ADdnSVM_m, ONLY : dnS_t
    IMPLICIT NONE
  
      CLASS (QML_OneD_2Quadra_t),    intent(in)     :: QModel
      TYPE (dnS_t),                  intent(in)     :: dnQ(:)
      integer,                       intent(in)     :: list_Op(:)
      TYPE (dnS_t),                  intent(inout)  :: Mat_OF_ScalOpDia(:,:,:)
      integer,                       intent(in)     :: nderiv
  
      IF (QModel%nb_ScalOp /= size(list_Op)) THEN
        write(out_unit,*) 'ERROR in EvalScalOp_QML_OneD_2Quadra'
        write(out_unit,*) '  QModel%nb_ScalOp and size(list_Op) are inconsistent'
        write(out_unit,*) '  QModel%nb_ScalOp (including potential)',QModel%nb_ScalOp
        write(out_unit,*) ' size(list_Op)',size(list_Op)

        STOP 'ERROR in EvalScalOp_QML_OneD_2Quadra: QModel%nb_ScalOp and size(list_Op) are inconsistent'
      END IF
 
      CALL EvalPot_QML_OneD_2Quadra(QModel,Mat_OF_ScalOpDia(:,:,1),dnQ,nderiv)

      IF (QModel%nb_ScalOp > 1) THEN ! dipole moments (diagonals+transition)
        Mat_OF_ScalOpDia(1,1,2) = dnQ(1) * QModel%Z
        Mat_OF_ScalOpDia(2,2,2) = dnQ(1) * QModel%Z
        Mat_OF_ScalOpDia(1,2,2) = QModel%Mu12
        Mat_OF_ScalOpDia(2,1,2) = QModel%Mu12
      END IF
  
  END SUBROUTINE EvalScalOp_QML_OneD_2Quadra
  SUBROUTINE Write_QML_OneD_2Quadra(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_OneD_2Quadra_t), intent(in) :: QModel
    integer,                   intent(in) :: nio

    write(nio,*) '========================================'
    write(nio,*) 'QML_OneD_2Quadra current parameters'
    CALL QModel%QML_Empty_t%Write_QModel(nio)
    write(nio,*)
    write(nio,*) ' Parameters for the vibrational mode'
    write(nio,*) ' MR:   ',QModel%MR
    write(nio,*) ' k:    ',QModel%k
    write(nio,*) ' R1:   ',QModel%R1
    write(nio,*) ' R2:   ',QModel%R2
    write(nio,*) ' a:    ',QModel%a
    write(nio,*) ' b:    ',QModel%b
    write(nio,*) ' R3:   ',QModel%R3
    write(nio,*) ' Delta:',QModel%Delta

    write(nio,*) ' Parameters for the dipole'
    write(nio,*) ' Mu12:   ',QModel%Mu12
    write(nio,*) ' Z:      ',QModel%Z ! diagonal part
    write(nio,*)
    write(nio,*) 'end QML_OneD_2Quadra current parameters'
    write(nio,*) '========================================'
    flush(nio)

  END SUBROUTINE Write_QML_OneD_2Quadra

  SUBROUTINE RefValues_QML_OneD_2Quadra(QModel,err,nderiv,Q0,dnMatV,d0GGdef,option)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_OneD_2Quadra_t), intent(in)              :: QModel
    integer,                   intent(inout)           :: err
    integer,                   intent(in)              :: nderiv
      
    real (kind=Rkind),         intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),            intent(inout), optional :: dnMatV
    real (kind=Rkind),         intent(inout), optional :: d0GGdef(:,:)
    integer,                   intent(in),    optional :: option

    integer :: option_loc
    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_OneD_2Quadra'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      IF (present(option))  write(out_unit,*) 'option',option
      flush(out_unit)
    END IF

    IF (.NOT. QModel%Init) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'The model is not initialized!'
      err = -1
      RETURN
    ELSE
      err = 0
    END IF

    IF (present(option)) THEN
      option_loc = option
    ELSE
      option_loc = 1
    END IF


    SELECT CASE (option_loc)
    CASE (1)
      IF (present(Q0))      CALL RefValues_QML_OneD_2Quadra1(QModel,err,nderiv=nderiv,Q0=Q0)
      IF (present(dnMatV))  CALL RefValues_QML_OneD_2Quadra1(QModel,err,nderiv=nderiv,dnMatV=dnMatV)
      IF (present(d0GGdef)) CALL RefValues_QML_OneD_2Quadra1(QModel,err,nderiv=nderiv,d0GGdef=d0GGdef)
    CASE (2)
      IF (present(Q0))      CALL RefValues_QML_OneD_2Quadra2(QModel,err,nderiv=nderiv,Q0=Q0)
      IF (present(dnMatV))  CALL RefValues_QML_OneD_2Quadra2(QModel,err,nderiv=nderiv,dnMatV=dnMatV)
      IF (present(d0GGdef)) CALL RefValues_QML_OneD_2Quadra2(QModel,err,nderiv=nderiv,d0GGdef=d0GGdef)
    CASE Default
      STOP 'ERROR in RefValues_QML_OneD_2Quadra: wrong option. Possible value(s): 1, 2'
    END SELECT


    IF (debug) THEN
      write(out_unit,*) 'present Q0 dnMatV d0GGdef',present(Q0),present(dnMatV),present(d0GGdef)
      IF (present(Q0))      write(out_unit,*) 'Q0',Q0
      IF (present(dnMatV))  THEN
        write(out_unit,*) 'dnMatV is present'
        write(out_unit,*) 'dnMatV is allocated',(.NOT. Check_NotAlloc_dnMat(dnMatV,nderiv))
        CALL write_dnMat(dnMatV,info='dnMatV')
      END IF
      IF (present(d0GGdef)) write(out_unit,*) 'd0GGdef',d0GGdef
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_OneD_2Quadra
 SUBROUTINE RefValues_QML_OneD_2Quadra1(QModel,err,Q0,dnMatV,d0GGdef,nderiv)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_OneD_2Quadra_t), intent(in)              :: QModel

    integer,           intent(inout)           :: err

    integer,           intent(in)              :: nderiv
    real (kind=Rkind), intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),    intent(inout), optional :: dnMatV

    real (kind=Rkind), intent(inout), optional :: d0GGdef(:,:)

    real (kind=Rkind), allocatable :: d0(:,:),d1(:,:,:),d2(:,:,:,:),d3(:,:,:,:,:),V(:)
    integer        :: i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_OneD_2Quadra1'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) ' nderiv ',nderiv
      flush(out_unit)
    END IF

    IF (.NOT. QModel%Init) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'The model is not initialized!'
      err = -1
      RETURN
    ELSE
      err = 0
    END IF

    IF (present(Q0)) THEN
      IF (size(Q0) /= QModel%ndim) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'incompatible Q0 size:'
        write(out_unit,*) 'size(Q0), ndimQ:',size(Q0),QModel%ndim
        err = 1
        Q0(:) = HUGE(ONE)
        RETURN
      END IF
      Q0(:) = [3.00000000000000_Rkind]
    END IF

    IF (present(dnMatV)) THEN
      err = 0
 !Energy (Hartree)   1.0000000000000000E-002_Rkind,1.0057264266491435E-003_Rkind,1.0057264266491435E-003_Rkind,8.9999999999999997E-002_Rkind,
 !                   2.0000000000000000E-002_Rkind,5.2800637399080036E-003_Rkind,5.2800637399080036E-003_Rkind,-5.9999999999999998E-002_Rkind,
  !                  2.0000000000000000E-002_Rkind,2.1685976074622155E-002_Rkind,2.1685976074622155E-002_Rkind,2.0000000000000000E-002_Rkind

      IF (nderiv >= 0) THEN ! no derivative
        V  = [1.0000000000000000E-002_Rkind,1.0057264266491435E-003_Rkind, &
              1.0057264266491435E-003_Rkind,8.9999999999999997E-002_Rkind]
        d0 = reshape(V,shape=[QModel%nsurf,QModel%nsurf])
      END IF

      IF (nderiv >= 1) THEN ! 1st order derivatives
        V  = [2.0000000000000000E-002_Rkind,5.2800637399080036E-003_Rkind,  &
              5.2800637399080036E-003_Rkind,-5.9999999999999998E-002_Rkind]
        d1 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim])
      END IF

      IF (nderiv >= 2) THEN ! 2d order derivatives
        V  = [2.0000000000000000E-002_Rkind,2.1685976074622155E-002_Rkind,  &
              2.1685976074622155E-002_Rkind,2.0000000000000000E-002_Rkind]
        d2 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim,QModel%ndim])
      END IF
      SELECT CASE (nderiv)
      CASE(0)
        CALL set_dnMat(dnMatV,d0=d0)
      CASE(1)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1)
      CASE(2)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1,d2=d2)
      CASE Default
        STOP 'ERROR in RefValues_QML_OneD_2Quadra1: nderiv MUST < 3'
      END SELECT

    END IF

    IF (present(d0GGdef)) THEN 
      V  = [5.0000000000000002E-005_Rkind]
      d0GGdef = reshape(V,shape=[QModel%ndim,QModel%ndim])
    END IF


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_OneD_2Quadra1
 SUBROUTINE RefValues_QML_OneD_2Quadra2(QModel,err,Q0,dnMatV,d0GGdef,nderiv)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_OneD_2Quadra_t), intent(in)              :: QModel

    integer,           intent(inout)           :: err

    integer,           intent(in)              :: nderiv
    real (kind=Rkind), intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),    intent(inout), optional :: dnMatV

    real (kind=Rkind), intent(inout), optional :: d0GGdef(:,:)

    real (kind=Rkind), allocatable :: d0(:,:),d1(:,:,:),d2(:,:,:,:),d3(:,:,:,:,:),V(:)
    integer        :: i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_OneD_2Quadra2'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) ' nderiv ',nderiv
      flush(out_unit)
    END IF

    IF (.NOT. QModel%Init) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'The model is not initialized!'
      err = -1
      RETURN
    ELSE
      err = 0
    END IF

    IF (present(Q0)) THEN
      IF (size(Q0) /= QModel%ndim) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'incompatible Q0 size:'
        write(out_unit,*) 'size(Q0), ndimQ:',size(Q0),QModel%ndim
        err = 1
        Q0(:) = HUGE(ONE)
        RETURN
      END IF
      Q0(:) = [3.00000000000000_Rkind]
    END IF

    IF (present(dnMatV)) THEN
      err = 0
 !Energy (Hartree)   9.9873584270513266E-003   2.7105054312137611E-020   0.0000000000000000        9.0012641572948665E-002   
 !                   1.9854646873414719E-002   0.0000000000000000        0.0000000000000000       -5.9854646873414720E-002
 !                   1.8468075275483067E-002   0.0000000000000000        0.0000000000000000        2.1531924724516933E-002


      IF (nderiv >= 0) THEN ! no derivative
        V  = [9.9873584270513266E-003_Rkind,ZERO, ZERO,9.0012641572948665E-002_Rkind]
        d0 = reshape(V,shape=[QModel%nsurf,QModel%nsurf])
      END IF

      IF (nderiv >= 1) THEN ! 1st order derivatives
        V  = [1.9854646873414719E-002_Rkind,ZERO, ZERO,-5.9854646873414720E-002_Rkind]
        d1 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim])
      END IF

      IF (nderiv >= 2) THEN ! 2d order derivatives
        V  = [1.8468075275483067E-002_Rkind,ZERO, ZERO,2.1531924724516933E-002_Rkind]
        d2 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim,QModel%ndim])
      END IF
      SELECT CASE (nderiv)
      CASE(0)
        CALL set_dnMat(dnMatV,d0=d0)
      CASE(1)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1)
      CASE(2)
        CALL set_dnMat(dnMatV,d0=d0,d1=d1,d2=d2)
      CASE Default
        STOP 'ERROR in RefValues_QML_OneD_2Quadra1: nderiv MUST < 3'
      END SELECT

    END IF

    IF (present(d0GGdef)) THEN 
      V  = [5.0000000000000002E-005_Rkind]
      d0GGdef = reshape(V,shape=[QModel%ndim,QModel%ndim])
    END IF


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_OneD_2Quadra2
end module QML_OneD_2Quadra_m
