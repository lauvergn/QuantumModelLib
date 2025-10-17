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
module QML_OneD_Photons2_m
  USE QDUtil_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE


  TYPE, EXTENDS (QML_Empty_t) :: QML_OneD_Photons2_t
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

    ! photon mode
    real(kind=Rkind)     :: MP  = 1._Rkind !  mass associated to the photo mode
    real(kind=Rkind)     :: w   = 0.17_Rkind ! photo mode angular frequency

    ! coupling betwee the ibrational mode and the photon mode
    real(kind=Rkind)     :: Z      = 1._Rkind ! for the diagonal contribution
    real(kind=Rkind)     :: Mu12   = 1.0_Rkind ! transition dipole moment between the 2 diabatic states
    real(kind=Rkind)     :: lambda = 0.01_Rkind ! global parameter

  CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_OneD_Photons2
    PROCEDURE :: Write_QModel    => Write_QML_OneD_Photons2
  END TYPE QML_OneD_Photons2_t

  PUBLIC :: QML_OneD_Photons2_t,Init_QML_OneD_Photons2


contains

  FUNCTION Init_QML_OneD_Photons2(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m, ONLY : Identity_Mat
    IMPLICIT NONE

    TYPE (QML_OneD_Photons2_t)          :: QModel

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

    ! photon mode
    real(kind=Rkind)     :: MP  !  mass associated to the photo mode
    real(kind=Rkind)     :: w   ! photo mode angular frequency

    ! coupling betwee the ibrational mode and the photon mode
    real(kind=Rkind)     :: Z       ! for the diagonal contribution
    real(kind=Rkind)     :: Mu12    ! transition dipole moment between the 2 diabatic states
    real(kind=Rkind)     :: lambda  ! global parameter

    integer :: err_read
    namelist /OneD_Photons2/ MR,k,R1,R2,a,b,R3,Delta,   MP,w,  Z,Mu12,lambda

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_OneD_Photons2'
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
    QModel%pot_name = 'OneD_Photons2'
    QModel%nsurf    = 4 ! two photons and two electron states
    QModel%ndim     = 1 ! 1 vibrational mode 

    IF (read_param) THEN
      MR    = QModel%MR
      k     = QModel%k
      R1    = QModel%R1
      R2    = QModel%R2
      a     = QModel%a
      b     = QModel%b
      R3    = QModel%R3
      Delta = QModel%Delta

      MP = QModel%MP
      w  = QModel%w

      Z      = QModel%Z
      Mu12   = QModel%Mu12
      lambda = QModel%lambda

      read(nio_param_file,OneD_Photons2,IOSTAT=err_read)

      IF (err_read < 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' End-of-file or End-of-record'
        write(out_unit,*) ' The namelist "OneD_Photons2" is probably absent'
        write(out_unit,*) ' check your data!'
        write(out_unit,*)
        STOP ' ERROR in Init_QML_OneD_Photons2: The namelist "OneD_Photons2" is probably absent'
      ELSE IF (err_read > 0) THEN
        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' Some parameter names of the namelist "OneD_Photons2" are probaly wrong'
        write(out_unit,*) ' check your data!'
        write(out_unit,nml=OneD_Photons2)
        STOP ' ERROR in Init_QML_OneD_Photons2: Some parameter names of the namelist "OneD_Photons2" are probaly wrong'
      END IF

      QModel%MR    = MR
      QModel%k     = k 
      QModel%R1    = R1
      QModel%R2    = R2
      QModel%a     = a
      QModel%b     = b
      QModel%R3    = R3
      QModel%Delta = Delta

      QModel%MP = MP
      QModel%w  = w

      QModel%Z      = Z
      QModel%Mu12   = Mu12
      QModel%lambda = lambda

    END IF


    IF (debug) write(out_unit,*) 'init Q0 of OneD_Photons2'
    QModel%Q0 = [QModel%R1]


    IF (debug) write(out_unit,*) 'init d0GGdef of OneD_Photons2'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)
    QModel%d0GGdef = ONE/QModel%MR


    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_OneD_Photons2


  SUBROUTINE EvalPot_QML_OneD_Photons2(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS (QML_OneD_Photons2_t), intent(in)     :: QModel
    TYPE (dnS_t),               intent(in)     :: dnQ(:)
    TYPE (dnS_t),               intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,                    intent(in)     :: nderiv

    TYPE (dnS_t)  :: mXpY,vXmY,mZ,R,Q
    integer       :: i

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPot_QML_OneD_Photons2'
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

    ! DIAGONAL TERMS
    Mat_OF_PotDia(1,1) = HALF*QModel%k*(R-QModel%R1)**2 + 0.5D0*QModel%w
    Mat_OF_PotDia(2,2) = HALF*QModel%k*(R-QModel%R2)**2 + 0.5D0*QModel%w + QModel%Delta
    Mat_OF_PotDia(3,3) = HALF*QModel%k*(R-QModel%R1)**2 + 1.5D0*QModel%w
    Mat_OF_PotDia(4,4) = HALF*QModel%k*(R-QModel%R2)**2 + 1.5D0*QModel%w + QModel%Delta

    ! ELECTRON NUCLEAR INTERACTIONS
    Mat_OF_PotDia(1,2) = QModel%b*exp(-QModel%a*(R-QModel%R3)**2)  
    Mat_OF_PotDia(2,1) = Mat_OF_PotDia(1,2)
    Mat_OF_PotDia(3,4) = Mat_OF_PotDia(1,2)
    Mat_OF_PotDia(4,3) = Mat_OF_PotDia(1,2)

    ! NUCLEAR PHOTON INTERACTIONS
    Mat_OF_PotDia(1,3) = QModel%Z*SQRT(HALF*QModel%w)*QModel%lambda*R
    Mat_OF_PotDia(2,4) = Mat_OF_PotDia(1,3)
    Mat_OF_PotDia(3,1) = Mat_OF_PotDia(1,3)
    Mat_OF_PotDia(4,2) = Mat_OF_PotDia(1,3)

    ! ELECTRON PHOTON INTERACTIONS
    Mat_OF_PotDia(1,4) =  - SQRT(HALF*QModel%w)*QModel%lambda*QModel%Mu12
    Mat_OF_PotDia(2,3) = Mat_OF_PotDia(1,4)
    Mat_OF_PotDia(3,2) = Mat_OF_PotDia(1,4)
    Mat_OF_PotDia(4,1) = Mat_OF_PotDia(1,4)


    IF (debug) THEN
      write(out_unit,*) 'Mat_OF_PotDia'
      CALL Write_dnS( Mat_OF_PotDia(1,1),6)
      write(out_unit,*)
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE EvalPot_QML_OneD_Photons2

  SUBROUTINE Write_QML_OneD_Photons2(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_OneD_Photons2_t), intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) '========================================'
    write(nio,*) 'QML_OneD_Photons2 current parameters'
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

    write(nio,*) ' Parameters for the photon mode'
    write(nio,*) ' MP:',QModel%MP
    write(nio,*) ' w: ',QModel%w
    write(nio,*) ' Parameters for the coupling (vib/photon)'
    write(nio,*) ' Z:      ',QModel%Z
    write(nio,*) ' Mu12:   ',QModel%Mu12
    write(nio,*) ' lambda: ',QModel%lambda
    write(nio,*)
    write(nio,*) 'end QML_OneD_Photons2 current parameters'
    write(nio,*) '========================================'
    flush(nio)

  END SUBROUTINE Write_QML_OneD_Photons2

end module QML_OneD_Photons2_m
