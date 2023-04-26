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
MODULE QML_Empty_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE ADdnSVM_m, ONLY : dnMat_t
  IMPLICIT NONE

  PRIVATE

  TYPE :: QML_Empty_t
    logical :: Init        = .FALSE.

    integer :: nsurf       = 0
    integer :: ndim        = 0
    integer :: ndimQ       = 0
    integer :: ndimCart    = 0

    ! for functions (used in the fit, Qeq(Q(1:ndimFunc)), hess() ....)
    integer :: ndimFunc         = 0
    integer :: nb_Func          = 0
    integer :: IndexFunc_Ene    = 0
    integer :: IndexFunc_Qop    = 0
    integer :: IndexFunc_Grad   = 0
    integer :: IndexFunc_Hess   = 0


    logical :: numeric     = .FALSE.
    logical :: no_ana_der  = .FALSE. ! to force numerical derivatives
                                     ! for potential without analitical derivatives

    logical :: Cart_TO_Q        = .FALSE. ! to perform the Cartesian to model coordinates
    logical :: MassWeighted     = .FALSE. ! Cartesian with mass Weighted
    logical :: AbInitio         = .FALSE. ! To use abitio calculation (experimental)
    integer :: nb_ScalOp        = 1 ! number scalar operators including the potential
                                    ! numbering: [0: potential, 1...: other operators]

    logical :: Phase_Following  = .TRUE.
    logical :: Phase_Checking   = .TRUE.
    logical :: adiabatic        = .TRUE.
    integer :: option           = 0
    logical :: PubliUnit        = .FALSE. ! when PubliUnit=.TRUE., the units of a reference (publi ...) are used. Default (atomic unit)
    logical :: In_a_Model       = .FALSE.


    logical :: Vib_adia         = .FALSE.
    integer :: nb_Channels      = 0
    integer, allocatable :: list_act(:)
    integer, allocatable :: list_inact(:)

    logical :: print_EigenVec_Grid  = .FALSE.
    logical :: print_EigenVec_Basis = .FALSE.

    character (len=:),  allocatable :: pot_name
    real (kind=Rkind),  allocatable :: d0GGdef(:,:)
    real (kind=Rkind),  allocatable :: masses(:) ! atomic masses (in au)
    real (kind=Rkind),  allocatable :: Q0(:)

    !Vec0 must be allocatable, to be able to deallocate with deallocate of the QML_Empty_t variable.
    TYPE (dnMat_t),     allocatable :: Vec0 ! to get the correct phase of the adiatic couplings
    CONTAINS
      PROCEDURE :: EvalPot_QModel         => EvalPot_QML_Empty
      PROCEDURE :: EvalPotAbInitio_QModel => EvalPotAbInitio_QML_Empty
      PROCEDURE :: EvalScalOp_QModel      => EvalScalOp_QML_Empty
      PROCEDURE :: EvalFunc_QModel        => EvalFunc_QML_Empty
      PROCEDURE :: Write_QModel           => Write_QML_Empty
      PROCEDURE :: Write0_QModel          => Write0_QML_Empty
     !PROCEDURE :: get2_Q0_QModel         => get2_Q0_QML_Empty
      PROCEDURE :: get_d0GGdef_QModel     => get_d0GGdef_QML_Empty
      PROCEDURE :: Cart_TO_Q_QModel       => Cart_TO_Q_QML_Empty
      PROCEDURE, PRIVATE :: Empty2_TO_Empty1_QML_Empty
      GENERIC,   PUBLIC  :: assignment(=) => Empty2_TO_Empty1_QML_Empty
      END TYPE QML_Empty_t

  INTERFACE get_Q0_QModel
    MODULE PROCEDURE get_Q0_QML_Empty
  END INTERFACE

  INTERFACE Qact_TO_Q
    MODULE PROCEDURE Qact_TO_Q_QML_Empty
  END INTERFACE

  !INTERFACE Empty2_TO_Empty1
  !  MODULE PROCEDURE Empty2_TO_Empty1_QML_Empty
  !END INTERFACE

  PUBLIC :: QML_Empty_t
  PUBLIC :: get_Q0_QModel, get2_Q0_QML_Empty, Qact_TO_Q

  CONTAINS

  SUBROUTINE Empty2_TO_Empty1_QML_Empty(QModel,QModel_in)
    USE ADdnSVM_m, ONLY : dealloc_dnMat
    IMPLICIT NONE

    CLASS (QML_Empty_t),  intent(inout)   :: QModel
    TYPE (QML_Empty_t),  intent(in)      :: QModel_in ! variable to transfer info to the init

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Empty2_TO_Empty1_QML_Empty'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF
    QModel%Init             = QModel_in%Init

    QModel%nsurf            = QModel_in%nsurf
    QModel%ndim             = QModel_in%ndim
    QModel%ndimQ            = QModel_in%ndimQ
    QModel%ndimCart         = QModel_in%ndimCart

    QModel%numeric          = QModel_in%numeric
    QModel%adiabatic        = QModel_in%adiabatic
    QModel%Phase_Following  = QModel_in%Phase_Following
    QModel%Phase_Checking   = QModel_in%Phase_Checking
    QModel%option           = QModel_in%option
    QModel%PubliUnit        = QModel_in%PubliUnit

    QModel%no_ana_der       = QModel_in%no_ana_der
    QModel%Cart_TO_Q        = QModel_in%Cart_TO_Q
    QModel%MassWeighted     = QModel_in%MassWeighted

    QModel%AbInitio         = QModel_in%AbInitio
    QModel%nb_ScalOp        = QModel_in%nb_ScalOp

    IF (QModel%adiabatic) THEN
      write(out_unit,*) 'Adiabatic potential . . .'
    ELSE
      write(out_unit,*) 'Non-adiabatic potential . . .'
    END IF
    flush(out_unit)

    IF (QModel%numeric) THEN
      write(out_unit,*) 'You have decided to perform a numeric checking of the analytic formulas.'
    END IF

    IF (allocated(QModel%pot_name)) deallocate(QModel%pot_name )
    QModel%pot_name = 'QML_Empty'

    IF (allocated(QModel%Vec0)) THEN
      CALL dealloc_dnMat(QModel%Vec0)
      deallocate(QModel%Vec0)
    END If

    IF (allocated(QModel%d0GGdef)) deallocate(QModel%d0GGdef)
    IF (allocated(QModel%Q0))      deallocate(QModel%Q0)
    IF (allocated(QModel%masses))  deallocate(QModel%masses)

    QModel%ndimFunc     = QModel_in%ndimFunc
    QModel%nb_Func      = QModel_in%nb_Func



    QModel%Vib_adia     = QModel_in%Vib_adia
    QModel%nb_Channels  = QModel_in%nb_Channels

    QModel%print_EigenVec_Basis = QModel_in%print_EigenVec_Basis
    QModel%print_EigenVec_Grid  = QModel_in%print_EigenVec_Grid

    IF (allocated(QModel_in%list_act)) THEN
      QModel%list_act     = QModel_in%list_act
    END IF

    IF (allocated(QModel%list_inact)) deallocate(QModel%list_inact)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE Empty2_TO_Empty1_QML_Empty
  SUBROUTINE get2_Q0_QML_Empty(QModel,Q0)
    IMPLICIT NONE

    real (kind=Rkind),    intent(inout)  :: Q0(:)
    CLASS(QML_Empty_t),   intent(in)     :: QModel

    IF (size(Q0) /= QModel%ndim) THEN
      STOP 'STOP in get_Q0_QML_Empty, wrong ndim size.'
    END IF

    IF (allocated(QModel%Q0)) THEN
      Q0(:) =  QModel%Q0
    END IF


  END SUBROUTINE get2_Q0_QML_Empty

  SUBROUTINE get_Q0_QML_Empty(QModel,Q0,err_Q0)
    USE QDUtil_m, ONLY : Rkind
    IMPLICIT NONE

    CLASS(QML_Empty_t),   intent(in)              :: QModel
    real (kind=Rkind),    intent(inout)           :: Q0(:)
    integer,              intent(inout), optional ::  err_Q0

    IF (size(Q0) /= QModel%ndim) THEN
      STOP 'STOP in get_Q0_QML_Empty, wrong ndim size.'
    END IF

    IF (allocated(QModel%Q0)) THEN
      Q0(:) =  QModel%Q0
      err_Q0 = 0
    ELSE
      err_Q0 = 1
    END IF

  END SUBROUTINE get_Q0_QML_Empty

  FUNCTION get_d0GGdef_QML_Empty(QModel) RESULT(d0GGdef)
    IMPLICIT NONE

    real (kind=Rkind),   allocatable              :: d0GGdef(:,:)
    CLASS(QML_Empty_t),             intent(in)    :: QModel

    integer :: i,nact

    IF (allocated(d0GGdef)) deallocate(d0GGdef)


    IF (allocated(QModel%d0GGdef)) THEN
      IF (allocated(QModel%list_act)) THEN
        d0GGdef =  QModel%d0GGdef(QModel%list_act,QModel%list_act)
      else
        d0GGdef =  QModel%d0GGdef(:,:)
      END IF
    ELSE
      nact = size(QModel%list_act)
      allocate(d0GGdef(nact,nact))
      d0GGdef = ZERO
      DO i=1,nact
        d0GGdef(i,i) = ONE
      END DO
    END IF

    !write(out_unit,*) 'alloc Q0',allocated(Q0)

  END FUNCTION get_d0GGdef_QML_Empty
  SUBROUTINE EvalPot_QML_Empty(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m, ONLY : dnS_t
    IMPLICIT NONE

    CLASS (QML_Empty_t),    intent(in)     :: QModel
    TYPE (dnS_t),           intent(in)     :: dnQ(:)
    TYPE (dnS_t),           intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,                intent(in)     :: nderiv

    write(out_unit,*) 'ERROR in EvalPot_QModel (EvalPot_QML_Empty)'
    write(out_unit,*) '  The intialized model does not have EvalPot_QModel subroutine!'
    STOP 'ERROR in EvalPot_QML_Empty: the intialized model does not have EvalPot_QModel subroutine'
    Mat_OF_PotDia(:,:) = ZERO

  END SUBROUTINE EvalPot_QML_Empty

  SUBROUTINE EvalScalOp_QML_Empty(QModel,Mat_OF_ScalOpDia,list_Op,dnQ,nderiv)
    USE QDUtil_m,  ONLY : ZERO
    USE ADdnSVM_m, ONLY : dnS_t
    IMPLICIT NONE
  
      CLASS (QML_Empty_t),    intent(in)     :: QModel
      TYPE (dnS_t),           intent(in)     :: dnQ(:)
      integer,                intent(in)     :: list_Op(:)
      TYPE (dnS_t),           intent(inout)  :: Mat_OF_ScalOpDia(:,:,:)
      integer,                intent(in)     :: nderiv
  
      write(out_unit,*) 'ERROR in EvalScalOp_QModel (EvalScalOp_QML_Empty)'
      write(out_unit,*) '  The intialized model does not have EvalScalOp_QModel subroutine!'
      STOP 'ERROR in EvalScalOp_QML_Empty: the intialized model does not have EvalScalOp_QModel subroutine'
      Mat_OF_ScalOpDia(:,:,:) = ZERO
  
  END SUBROUTINE EvalScalOp_QML_Empty

  SUBROUTINE EvalPotAbInitio_QML_Empty(QModel,Mat_OF_PotDia,dnX,nderiv)
    USE QMLLib_UtilLib_m, ONLY : Find_Label
    USE ADdnSVM_m, ONLY :  dnS_t, get_d0, set_dnS, write_dnS
    IMPLICIT NONE

    CLASS (QML_Empty_t),    intent(in)     :: QModel
    TYPE (dnS_t),           intent(in)     :: dnX(:,:)
    TYPE (dnS_t),           intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,                intent(in)     :: nderiv

    integer                       :: i,j,Z,nio_otf,err
    logical                       :: located
    real(kind=Rkind)              :: d0E
    real(kind=Rkind), allocatable :: d1E(:)
    real(kind=Rkind), allocatable :: d2E(:,:)

    character (len=Name_longlen)  :: labelR
    character (len=Name_len)      :: name1_i

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalPotAbInitio_QML_Empty'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) ' BEGINNING ',name_sub
      write(out_unit,*) '   nderiv:       ',nderiv
      flush(out_unit)
    END IF

    ! check if Cart_TO_Q=t
    IF (.NOT. QModel%Cart_TO_Q) THEN
      write(out_unit,*) 'ERROR in EvalPotAbInitio_QML_Empty'
      write(out_unit,*) 'Cartessian coordinates must be used and Cart_TO_Q=F'
      write(out_unit,*) 'Check your data'
      STOP 'ERROR in EvalPotAbInitio_QML_Empty'
    END IF

    ! generate the input file
    open(newunit=nio_otf,file='xx.com',form='formatted')
    write(nio_otf,*) '%chk=xx.chk'

    IF (nderiv == 0) write(nio_otf,*) '#n HF STO-3G'
    IF (nderiv == 1) write(nio_otf,*) '#n HF STO-3G force '
    IF (nderiv == 2) write(nio_otf,*) '#n HF STO-3G freq=noraman '

    write(nio_otf,*) '# unit=(au,rad) nosymm FormCheck=All'
    write(nio_otf,*) '# integral=(grid=ultrafine) '
    write(nio_otf,*) '# iop(1/18=40,2/9=111,2/11=1,7/33=1,1/33=1)'
    write(nio_otf,*) '# iop(3/27=30,11/27=30)'

    write(nio_otf,*) ' '
    write(nio_otf,*) ' xx shape dnX',shape(dnX),'ndimCart',QModel%ndimCart
    write(nio_otf,*) ' '
    write(nio_otf,*) '0 2'


    DO i=1,QModel%ndimCart/3
      Z = 1
      write(nio_otf,13) Z,0,get_d0(dnX(:,i))
13    format(i5,1x,i5,3(1x,f20.15))
    END DO
    write(nio_otf,*) ' '
    close(nio_otf)

    !calculation + checking the normal termination
    CALL EXECUTE_COMMAND_LINE('gauss09.run  xx &> gaussexe.log')

    !checking the normal termination
    located = .FALSE.
    open(newunit=nio_otf,file='xx.log',status='old',position='append',form='formatted')
    backspace(nio_otf,err=999)
    read(nio_otf,'(a32)',err=999) labelR
    !write(out_unit,*) 'last line: ',labelR
    located = verify(labelR,' Normal termination of Gaussian') == 0
    close(nio_otf)
999 CONTINUE
    IF (.NOT. located) THEN
      write(out_unit,*) 'ERROR in EvalPotAbInitio_QML_Empty'
      write(out_unit,*) 'no line: "Normal termination of Gaussian"'
      write(out_unit,*) 'log file last line: ',labelR
      STOP
    END IF


    !- read the energy from the file energy
    open(newunit=nio_otf,file='xx.fchk',status='old',form='formatted')
    CALL Find_Label(nio_otf,'Total Energy',located)
    IF (debug) write(out_unit,*) 'located: Total Energy',located
    IF (located) THEN
      read(nio_otf,*,iostat=err) name1_i,d0E
    ELSE
      err = -1
    END IF

    IF (.NOT. located .OR. err /=0) THEN
      write(out_unit,*) 'ERROR in ',name_sub
      write(out_unit,*) 'I cannot find the energy in : xx.fchk'
      write(out_unit,*) 'located,err',located,err
      STOP
    END IF
    close(nio_otf)

    !- read the gradient
    IF (nderiv >= 1) THEN
      allocate(d1E(QModel%ndimCart))
      open(newunit=nio_otf,file='xx.fchk',status='old',form='formatted')

      CALL Find_Label(nio_otf,'Cartesian Gradient',located)
      IF (debug) write(out_unit,*) 'located: Cartesian Gradient',located
      IF (located) THEN
        read(nio_otf,*,iostat=err)
        read(nio_otf,*,iostat=err) d1E(:)
      END IF
      IF (.NOT. located .OR. err /=0) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'I cannot find the Gradient in : xx.fchk'
        write(out_unit,*) 'located,err',located,err
        STOP
      END IF
      close(nio_otf)
    END IF

    IF (nderiv == 2) THEN
      allocate(d2E(QModel%ndimCart,QModel%ndimCart))
      open(newunit=nio_otf,file='xx.fchk',status='old',form='formatted')
      CALL Find_Label(nio_otf,'Cartesian Force Constants',located)
      IF (debug) write(out_unit,*) 'located: Cartesian Force Constants (hessian)',located
      IF (located) THEN
        read(nio_otf,*,iostat=err)
        read(nio_otf,*,iostat=err) ((d2E(i,j),i=1,j),j=1,QModel%ndimCart)
      END IF
      IF (.NOT. located .OR. err /=0) THEN
        write(out_unit,*) 'ERROR in ',name_sub
        write(out_unit,*) 'I cannot find the hessian in : xx.fchk'
        write(out_unit,*) 'located,err',located,err
        STOP
      END IF

      DO j=1,QModel%ndimCart
      DO i=1,j-1
        d2E(j,i) = d2E(i,j)
      END DO
      END DO
      close(nio_otf)
    END IF

    IF (nderiv == 0) THEN
      CALL set_dnS(Mat_OF_PotDia(1,1),d0=d0E)
    ELSE IF (nderiv == 1) THEN
      CALL set_dnS(Mat_OF_PotDia(1,1),d0=d0E,d1=d1E)
    ELSE IF (nderiv == 2) THEN
      CALL set_dnS(Mat_OF_PotDia(1,1),d0=d0E,d1=d1E,d2=d2E)
    END IF

    IF (debug) THEN
      CALL write_dnS(Mat_OF_PotDia(1,1),info='dnE')
      write(out_unit,*) ' END ',name_sub
      flush(out_unit)
    END IF

  END SUBROUTINE EvalPotAbInitio_QML_Empty

  SUBROUTINE EvalFunc_QML_Empty(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m, ONLY : dnS_t
    IMPLICIT NONE

    CLASS (QML_Empty_t),    intent(in)     :: QModel
    TYPE (dnS_t),           intent(in)     :: dnQ(:)
    TYPE (dnS_t),           intent(inout)  :: Func(:)
    integer,                intent(in)     :: nderiv

    integer :: i

    write(out_unit,*) 'ERROR in EvalFunc_QModel (EvalFunc_QML_Empty)'
    write(out_unit,*) '  The intialized model does not have EvalFunc_QModel subroutine!'
    STOP 'ERROR in EvalFunc_QML_Empty: the intialized model does not have EvalFunc_QModel subroutine'
    DO i=1,size(Func)
      Func(i) = ZERO
    END DO

  END SUBROUTINE EvalFunc_QML_Empty

  SUBROUTINE Write_QML_Empty(QModel,nio)
    USE QDUtil_m,  ONLY : Write_RMat => Write_Mat, Write_RVec => Write_Vec
    IMPLICIT NONE

    CLASS (QML_Empty_t), intent(in)    :: QModel
    integer,             intent(in)    :: nio

    IF (QModel%In_a_Model) RETURN

    write(nio,*) 'Potential parameters are written just below'
    write(nio,*) 'Init:                      ',QModel%Init
    write(nio,*) 'In_a_Model:                ',QModel%In_a_Model

    write(nio,*) 'option:                    ',QModel%option
    write(nio,*) 'AbInitio:                  ',QModel%AbInitio

    write(nio,*)
    write(nio,*) 'nsurf:                     ',QModel%nsurf
    write(nio,*) 'ndim:                      ',QModel%ndim
    write(nio,*) 'numeric:                   ',QModel%numeric
    write(nio,*) 'no analytical derivatives: ',QModel%no_ana_der
    write(nio,*) 'Cartesian => model coord.: ',QModel%Cart_TO_Q
    write(nio,*) 'Cart. With Mass Weighted:  ',QModel%MassWeighted
    write(nio,*) 'ndimQ:                     ',QModel%ndimQ
    write(nio,*) 'ndimCart:                  ',QModel%ndimCart
    write(nio,*)
    write(nio,*) 'ndimFunc:                  ',QModel%ndimFunc
    write(nio,*) 'nb_Func:                   ',QModel%nb_Func
    write(nio,*) 'IndexFunc_Ene:             ',QModel%IndexFunc_Ene
    write(nio,*) 'IndexFunc_Qop:             ',QModel%IndexFunc_Qop
    write(nio,*) 'IndexFunc_Grad:            ',QModel%IndexFunc_Grad
    write(nio,*) 'IndexFunc_Hess:            ',QModel%IndexFunc_Hess

    write(nio,*)
    write(nio,*) 'adiabatic:                 ',QModel%adiabatic
    write(nio,*) 'Vib_adia:                  ',QModel%Vib_adia
    write(nio,*) 'Phase_Following:           ',QModel%Phase_Following
    write(nio,*) 'Phase_Checking:            ',QModel%Phase_Checking

    IF (QModel%Vib_adia) THEN
      write(nio,*) 'nb_Channels:                  ',QModel%nb_Channels
      IF (allocated(QModel%list_act)) &
        write(nio,*) 'list_act(:):               ',QModel%list_act(:)
      IF (allocated(QModel%list_inact)) &
        write(nio,*) 'list_inact(:):             ',QModel%list_inact(:)
    END IF
    write(nio,*) 'print_EigenVec Basis/Grid: ',QModel%print_EigenVec_Basis,QModel%print_EigenVec_Grid


    IF (allocated(QModel%pot_name)) write(nio,*) 'pot_name: ',QModel%pot_name
    write(nio,*)

    IF (allocated(QModel%d0GGdef)) THEN
      write(nio,*) 'Deformation metric tensor (~ 1/Mii)'
      CALL Write_RMat(QModel%d0GGdef,nio,nbcol=5)
    END IF

    IF (allocated(QModel%Q0)) THEN
      write(nio,*) 'Reference Coordinate values, Q0(:)'
      CALL Write_RVec(QModel%Q0,nio,nbcol=5)
    END IF

    IF (allocated(QModel%masses)) THEN
      write(nio,*) 'Masses in au'
      CALL Write_RVec(QModel%masses,nio,nbcol=5)
    END IF
    flush(nio)

  END SUBROUTINE Write_QML_Empty
  SUBROUTINE Write0_QML_Empty(QModel,nio)
    USE QDUtil_m,  ONLY : Write_RMat => Write_Mat
    IMPLICIT NONE

    CLASS (QML_Empty_t), intent(in)    :: QModel
    integer,             intent(in)    :: nio


    write(nio,*) 'QUANTUM MODEL default parameters'
    flush(nio)

    write(nio,*)
    write(nio,*) 'Potential parameters are written just below'
    write(nio,*) 'Init:                      ',QModel%Init
    write(nio,*) 'option:                    ',QModel%option
    write(nio,*)
    write(nio,*) 'nsurf:                     ',QModel%nsurf
    write(nio,*) 'ndim:                      ',QModel%ndim
    write(nio,*) 'numeric:                   ',QModel%numeric
    write(nio,*) 'adiabatic:                 ',QModel%adiabatic
    write(nio,*) 'no analytical derivatives: ',QModel%no_ana_der
    write(nio,*) 'Cartesian => model coord.: ',QModel%Cart_TO_Q
    write(nio,*) 'ndimQ:                     ',QModel%ndimQ
    write(nio,*) 'ndimCart:                  ',QModel%ndimCart
    write(nio,*) 'ndimFunc:                  ',QModel%ndimFunc
    write(nio,*) 'nb_Func:                   ',QModel%nb_Func
    write(nio,*)

    IF (allocated(QModel%d0GGdef)) THEN
     write(nio,*) 'Deformation metric tensor (~ 1/Mii)'
     CALL Write_RMat(QModel%d0GGdef,nio,nbcol=5)
    END IF

    write(nio,*) 'END QUANTUM MODEL default parameters'
    flush(nio)


  END SUBROUTINE Write0_QML_Empty

  SUBROUTINE Cart_TO_Q_QML_Empty(QModel,dnX,dnQ,nderiv)
    USE ADdnSVM_m, ONLY : dnS_t
    IMPLICIT NONE

    CLASS(QML_Empty_t),      intent(in)    :: QModel
    TYPE (dnS_t),            intent(in)    :: dnX(:,:)
    TYPE (dnS_t),            intent(inout) :: dnQ(:)
    integer,                 intent(in)    :: nderiv

    write(out_unit,*) 'ERROR in Cart_TO_Q_QModel (Cart_TO_Q_QML_Empty)'
    write(out_unit,*) '  The intialized model does not have Cart_TO_Q transformation!'
    STOP 'ERROR in Cart_TO_Q_QML_Empty: the intialized model does not have Cart_TO_Q transformation'

  END SUBROUTINE Cart_TO_Q_QML_Empty
  SUBROUTINE Qact_TO_Q_QML_Empty(Qact,Q,list_act)
    USE QDUtil_m,  ONLY : Rkind
    IMPLICIT NONE

    real (kind=Rkind),       intent(in)    :: Qact(:)
    integer,                 intent(in)    :: list_act(:)
    real (kind=Rkind),       intent(inout) :: Q(:)


    integer    :: i

    DO i=1,size(Qact)
      Q(list_act(i)) = Qact(i)
    END DO

  END SUBROUTINE Qact_TO_Q_QML_Empty
END MODULE QML_Empty_m
