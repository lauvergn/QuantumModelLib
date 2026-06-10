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
module QML_OneD_Photons_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE


  TYPE, EXTENDS (QML_Empty_t) :: QML_OneD_Photons_t
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


    ! photon modes
    integer                        :: nb_cavity_modes = 1 ! default only one cavity mode
    real(kind=Rkind), allocatable  :: MP(:)  ! mass associated to the photo mode (default 1.)
    real(kind=Rkind), allocatable  :: w(:)   ! photo mode angular frequency (default 0.17)

    ! coupling betwee the ibrational mode and the cavity modes
    real(kind=Rkind)             :: Mu12 = 1.0_Rkind     ! transition dipole moment between the 2 diabatic states
    real(kind=Rkind),allocatable :: Z(:)       ! for the diagonal contribution (default value 1.)
    real(kind=Rkind),allocatable :: lambda(:)  ! global parameter (default value 0.01)

  CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_OneD_Photons
    PROCEDURE :: Write_QModel     => Write_QML_OneD_Photons
    PROCEDURE :: RefValues_QModel => RefValues_QML_OneD_Photons
  END TYPE QML_OneD_Photons_t

  PUBLIC :: QML_OneD_Photons_t,Init_QML_OneD_Photons


contains

  FUNCTION Init_QML_OneD_Photons(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m, ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_OneD_Photons_t)          :: QModel

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
    real(kind=Rkind)     :: Mu12    ! transition dipole moment between the 2 diabatic states

    ! photon modes
    real(kind=Rkind), allocatable  :: MP(:)  !  mass associated to the photo mode
    real(kind=Rkind), allocatable  :: w(:)   ! photo mode angular frequency

    ! coupling betwee the vibrational mode and the photon mode
    real(kind=Rkind),allocatable :: Z(:)       ! for the diagonal contribution
    real(kind=Rkind),allocatable :: lambda(:)  ! global parameter

    integer :: i

    integer :: err_read
    namelist /OneD_Photons/ MR,k,R1,R2,a,b,R3,Delta,   MP,w,  Z,Mu12,lambda

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_OneD_Photons'
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
    QModel%pot_name = 'OneD_Photons'
    QModel%nsurf    = 2

    IF (QModel%ndim <= 0 .OR. QModel%ndim == 2) THEN ! old initialization
      QModel%nb_cavity_modes = 1
      QModel%ndim     = 1+QModel%nb_cavity_modes  ! 1+1: 1 vibrational mode + 1 cavity mode
      QModel%MP       = [ONE]
      QModel%w        = [0.17_Rkind]
      QModel%Z        = [ONE]
      QModel%lambda   = [0.01_Rkind]
    ELSE IF (QModel%ndim == 1) THEN
      QModel%nb_cavity_modes = 0
      write(out_unit,*) 'WARNING in ',name_sub
      write(out_unit,*) '  wrong ndim and/or read_param'
      write(out_unit,*) '  ndim      ',QModel%ndim
      write(out_unit,*) '  read_param',read_param
      write(out_unit,*) '  ndim=1 => no cavity mode!'
    ELSE IF (QModel%ndim > 2 .AND. read_param) THEN ! here data need to be read
      QModel%nb_cavity_modes = QModel%ndim - 1
      allocate(QModel%MP(QModel%nb_cavity_modes))
      QModel%MP = ONE
      allocate(QModel%w(QModel%nb_cavity_modes))
      QModel%w = 0.17_Rkind ! it is silly to have allways the same frequencies
     allocate(QModel%Z(QModel%nb_cavity_modes))
      QModel%Z = ONE
      allocate(QModel%lambda(QModel%nb_cavity_modes))
      QModel%lambda = 0.01_Rkind
    ELSE ! ndim > 2 but read_param=.FALSE.
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) '  wrong ndim and/or read_param'
      write(out_unit,*) '  ndim      ',QModel%ndim
      write(out_unit,*) '  read_param',read_param
      write(out_unit,*) '  ndim=1 is not possible => no cavity mode'
      write(out_unit,*) '  ndim>2 and read_param=.FALSE. => cavity parameters cannot be read'
      STOP 'ERROR in Init_QML_OneD_Photons: wrong ndim and/or read_param'
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

      IF (QModel%nb_cavity_modes > 0) THEN
        MP     = QModel%MP
        w      = QModel%w
        Z      = QModel%Z
        lambda = QModel%lambda
      END IF

      read(nio_param_file,OneD_Photons,IOSTAT=err_read)

      IF (err_read < 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' End-of-file or End-of-record'
        write(out_unit,*) ' The namelist "OneD_Photons" is probably absent'
        write(out_unit,*) ' check your data!'
        write(out_unit,*)
        STOP ' ERROR in Init_QML_OneD_Photons: The namelist "OneD_Photons" is probably absent'
      ELSE IF (err_read > 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Some parameter names of the namelist "OneD_Photons" are probaly wrong'
        write(out_unit,*) ' check your data!'
        write(out_unit,nml=OneD_Photons)
        STOP ' ERROR in Init_QML_OneD_Photons: Some parameter names of the namelist "OneD_Photons" are probaly wrong'
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

      IF (QModel%nb_cavity_modes > 0) THEN
        QModel%MP     = MP
        QModel%w      = w
        QModel%Z      = Z
        QModel%lambda = lambda
      END IF

    END IF


    IF (debug) write(out_unit,*) 'init Q0 of OneD_Photons'
    allocate(QModel%Q0(QModel%ndim))
    QModel%Q0 = ZERO
    QModel%Q0(1) = QModel%R1


    IF (debug) write(out_unit,*) 'init d0GGdef of OneD_Photons'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)
    QModel%d0GGdef(1,1) = ONE/QModel%MR
    DO i=1,QModel%nb_cavity_modes
      QModel%d0GGdef(i+1,i+1) = ONE/QModel%MP(i)
    END DO

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_OneD_Photons


  SUBROUTINE EvalPot_QML_OneD_Photons(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS (QML_OneD_Photons_t), intent(in)     :: QModel
    TYPE (dnS_t),               intent(in)     :: dnQ(:)
    TYPE (dnS_t),               intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,                    intent(in)     :: nderiv

    TYPE (dnS_t)  :: mXpY,vXmY,mZ,R,Q
    integer       :: i

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_OneD_Photons'
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


    DO i=1,QModel%nb_cavity_modes
      Q = dnQ(i+1)
      Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + HALF*(QModel%w(i)*Q)**2 + QModel%Z(i)*QModel%w(i)*QModel%lambda(i)*Q*R
      Mat_OF_PotDia(2,2) = Mat_OF_PotDia(2,2) + HALF*(QModel%w(i)*Q)**2 + QModel%Z(i)*QModel%w(i)*QModel%lambda(i)*Q*R
      Mat_OF_PotDia(1,2) = Mat_OF_PotDia(1,2)                           - QModel%w(i)*QModel%lambda(i)*QModel%Mu12*Q
    END DO
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)

    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE EvalPot_QML_OneD_Photons

  SUBROUTINE Write_QML_OneD_Photons(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_OneD_Photons_t), intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) '========================================'
    write(nio,*) 'QML_OneD_Photons current parameters'
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

    write(nio,*) ' Parameters for the cavity mode(s)'
    write(nio,*) ' MP:',QModel%MP
    write(nio,*) ' w: ',QModel%w
    write(nio,*) ' Parameters for the coupling (vib/photon)'
    write(nio,*) ' Z:      ',QModel%Z
    write(nio,*) ' Mu12:   ',QModel%Mu12
    write(nio,*) ' lambda: ',QModel%lambda
    write(nio,*)
    write(nio,*) 'end QML_OneD_Photons current parameters'
    write(nio,*) '========================================'
    flush(nio)

  END SUBROUTINE Write_QML_OneD_Photons


  SUBROUTINE RefValues_QML_OneD_Photons(QModel,err,nderiv,Q0,dnMatV,d0GGdef,option)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_OneD_Photons_t), intent(in)              :: QModel
    integer,                  intent(inout)           :: err
    integer,                  intent(in)              :: nderiv
      
    real (kind=Rkind),        intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),           intent(inout), optional :: dnMatV
    real (kind=Rkind),        intent(inout), optional :: d0GGdef(:,:)
    integer,                  intent(in),    optional :: option

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_OneD_Photons'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
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

    SELECT CASE (option)
    CASE (1)
      IF (present(Q0))      CALL RefValues_QML_OneD_Photons1(QModel,err,nderiv=nderiv,Q0=Q0)
      IF (present(dnMatV))  CALL RefValues_QML_OneD_Photons1(QModel,err,nderiv=nderiv,dnMatV=dnMatV)
      IF (present(d0GGdef)) CALL RefValues_QML_OneD_Photons1(QModel,err,nderiv=nderiv,d0GGdef=d0GGdef)
    CASE Default
      STOP 'ERROR in RefValues_QML_OneD_Photons1: wrong option. Possible value(s): 1'
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

  END SUBROUTINE RefValues_QML_OneD_Photons
 SUBROUTINE RefValues_QML_OneD_Photons1(QModel,err,Q0,dnMatV,d0GGdef,nderiv)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_OneD_Photons_t), intent(in)              :: QModel

    integer,           intent(inout)           :: err

    integer,           intent(in)              :: nderiv
    real (kind=Rkind), intent(inout), optional :: Q0(:)
    TYPE (dnMat_t),    intent(inout), optional :: dnMatV

    real (kind=Rkind), intent(inout), optional :: d0GGdef(:,:)

    real (kind=Rkind), allocatable :: d0(:,:),d1(:,:,:),d2(:,:,:,:),d3(:,:,:,:,:),V(:)
    integer        :: i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='RefValues_QML_OneD_Photons1'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
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

! Q
!            2
!    3.0000000000000000        0.0000000000000000     
!  V
!            4
!    1.0000000000000000E-002   1.0057264266491435E-003   1.0057264266491435E-003   8.9999999999999997E-002
!  Grad
!            8
!    2.0000000000000000E-002   5.2800637399080036E-003   5.2800637399080036E-003  -5.9999999999999998E-002   5.1000000000000004E-003  -1.7000000000000001E-003  -1.7000000000000001E-003   5.1000000000000004E-003
!  Hess
!           16
!    2.0000000000000000E-002   2.1685976074622155E-002   2.1685976074622155E-002   2.0000000000000000E-002   1.7000000000000001E-003   0.0000000000000000        0.0000000000000000        1.7000000000000001E-003   1.7000000000000001E-003   0.0000000000000000        0.0000000000000000        1.7000000000000001E-003   2.8900000000000006E-002   0.0000000000000000        0.0000000000000000        2.8900000000000006E-002
!  d0GGdef
!            4
!    5.0000000000000002E-005   0.0000000000000000        0.0000000000000000        1.0000000000000000 

    IF (present(Q0)) THEN
      IF (size(Q0) /= QModel%ndim) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'incompatible Q0 size:'
        write(out_unit,*) 'size(Q0), ndimQ:',size(Q0),QModel%ndim
        err = 1
        Q0(:) = HUGE(ONE)
        RETURN
      END IF
      Q0(:) = [3.00000000000000_Rkind,0.0000000000000000_Rkind]
    END IF

    IF (present(dnMatV)) THEN
      err = 0

      IF (nderiv >= 0) THEN ! no derivative
        V  = [1.0000000000000000E-002_Rkind,1.0057264266491435E-003_Rkind,1.0057264266491435E-003_Rkind,8.9999999999999997E-002_Rkind]
        d0 = reshape(V,shape=[QModel%nsurf,QModel%nsurf])
      END IF

      IF (nderiv >= 1) THEN ! 1st order derivatives
        V  = [2.0000000000000000E-002_Rkind, 5.2800637399080036E-003_Rkind, 5.2800637399080036E-003_Rkind,-5.9999999999999998E-002_Rkind, &
              5.1000000000000004E-003_Rkind,-1.7000000000000001E-003_Rkind,-1.7000000000000001E-003_Rkind, 5.1000000000000004E-003_Rkind]
        d1 = reshape(V,shape=[QModel%nsurf,QModel%nsurf,QModel%ndim])
      END IF

      IF (nderiv >= 2) THEN ! 2d order derivatives
        V  = [2.0000000000000000E-002_Rkind,2.1685976074622155E-002_Rkind,2.1685976074622155E-002_Rkind,2.0000000000000000E-002_Rkind,&
              1.7000000000000001E-003_Rkind,0.0000000000000000_Rkind,     0.0000000000000000_Rkind,     1.7000000000000001E-003_Rkind,   &
              1.7000000000000001E-003_Rkind,0.0000000000000000_Rkind,     0.0000000000000000_Rkind,     1.7000000000000001E-003_Rkind,   &
              2.8900000000000006E-002_Rkind,0.0000000000000000_Rkind,     0.0000000000000000_Rkind,     2.8900000000000006E-002_Rkind]
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
        STOP 'ERROR in RefValues_QML_OneD_Photons1: nderiv MUST < 3'
      END SELECT

    END IF

    IF (present(d0GGdef)) THEN 
      V  = [5.0000000000000002E-005_Rkind,0.0000000000000000_Rkind,0.0000000000000000_Rkind,1.0000000000000000_Rkind]
      d0GGdef = reshape(V,shape=[QModel%ndim,QModel%ndim])
    END IF


    IF (debug) THEN
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE RefValues_QML_OneD_Photons1
end module QML_OneD_Photons_m
