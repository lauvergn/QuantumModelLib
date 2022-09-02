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
!    Copyright 2016 David Lauvergnat [1]
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
  USE QMLLib_NumParameters_m
  USE QMLLib_UtilLib_m
  USE ADdnSVM_m, ONLY : dnMat_t
  IMPLICIT NONE

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
    real (kind=Rkind),  allocatable :: Q0(:)

    !Vec0 must be allocatable, to be able to deallocate with deallocate of the QML_Empty_t variable.
    TYPE (dnMat_t),     allocatable :: Vec0 ! to get the correct phase of the adiatic couplings
    CONTAINS
      PROCEDURE :: EvalPot_QModel    => EvalPot_QML_Empty
      PROCEDURE :: Eval_QModel_Func   => EvalFunc_QML_Empty
      PROCEDURE :: Write_QModel       => Write_QML_Empty
      PROCEDURE :: Write0_QModel      => Write0_QML_Empty
     !PROCEDURE :: get2_Q0_QModel     => get2_Q0_QML_Empty
      PROCEDURE :: get_d0GGdef_QModel => get_d0GGdef_QML_Empty
      PROCEDURE :: Cart_TO_Q_QModel   => Cart_TO_Q_QML_Empty
  END TYPE QML_Empty_t

  INTERFACE get_Q0_QModel
    MODULE PROCEDURE get_Q0_QML_Empty
  END INTERFACE

  INTERFACE Qact_TO_Q
    MODULE PROCEDURE Qact_TO_Q_QML_Empty
  END INTERFACE

CONTAINS

  FUNCTION Init_QML_Empty(QModel_in) RESULT(QModel)
  USE QMLLib_UtilLib_m
  IMPLICIT NONE

    TYPE (QML_Empty_t)                  :: QModel

    TYPE(QML_Empty_t),  intent(in)      :: QModel_in ! variable to transfer info to the init

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Empty'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF


    CALL Init0_QML_Empty(QModel,QModel_in)

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_Empty
  SUBROUTINE Init0_QML_Empty(QModel,QModel_in)
  USE QMLLib_UtilLib_m
  USE ADdnSVM_m, ONLY : dealloc_dnMat
  IMPLICIT NONE

    TYPE (QML_Empty_t), intent(inout)   :: QModel

    TYPE(QML_Empty_t),  intent(in)      :: QModel_in ! variable to transfer info to the init

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init0_QML_Empty'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    !QModel = QModel_in   ! it does not work always with nagfor
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

    IF (QModel%adiabatic) THEN
      write(out_unitp,*) 'Adiabatic potential . . .'
    ELSE
      write(out_unitp,*) 'Non-adiabatic potential . . .'
    END IF
    flush(out_unitp)

    IF (QModel%numeric) THEN
      write(out_unitp,*) 'You have decided to perform a numeric checking of the analytic formulas.'
    END IF

    IF (allocated(QModel%pot_name)) deallocate(QModel%pot_name )
    QModel%pot_name = 'QML_Empty'

    IF (allocated(QModel%Vec0)) THEN
      CALL dealloc_dnMat(QModel%Vec0)
      deallocate(QModel%Vec0)
    END If

    IF (allocated(QModel%d0GGdef)) deallocate(QModel%d0GGdef)
    IF (allocated(QModel%Q0))      deallocate(QModel%Q0)

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
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Init0_QML_Empty
  SUBROUTINE get2_Q0_QML_Empty(QModel,Q0)
    IMPLICIT NONE

    real (kind=Rkind),    intent(inout)  :: Q0(:)
    CLASS(QML_Empty_t),  intent(in)     :: QModel

    IF (size(Q0) /= QModel%ndim) THEN
      STOP 'STOP in get_Q0_QML_Empty, wrong ndim size.'
    END IF

    IF (allocated(QModel%Q0)) THEN
      Q0(:) =  QModel%Q0
    END IF


  END SUBROUTINE get2_Q0_QML_Empty

  SUBROUTINE get_Q0_QML_Empty(QModel,Q0,err_Q0)
    IMPLICIT NONE

    CLASS(QML_Empty_t),  intent(in)              :: QModel
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

    real (kind=Rkind),   allocatable               :: d0GGdef(:,:)
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

    !write(out_unitp,*) 'alloc Q0',allocated(Q0)

  END FUNCTION get_d0GGdef_QML_Empty
  SUBROUTINE EvalPot_QML_Empty(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m, ONLY :  dnS_t
  IMPLICIT NONE

    CLASS (QML_Empty_t),    intent(in)     :: QModel
    TYPE (dnS_t),           intent(in)     :: dnQ(:)
    TYPE (dnS_t),           intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,                intent(in)     :: nderiv


    Mat_OF_PotDia(:,:) = ZERO

  END SUBROUTINE EvalPot_QML_Empty

  SUBROUTINE EvalFunc_QML_Empty(QModel,Func,dnQ,nderiv)
  USE ADdnSVM_m, ONLY :  dnS_t
  IMPLICIT NONE

    CLASS (QML_Empty_t),    intent(in)     :: QModel
    TYPE (dnS_t),           intent(in)     :: dnQ(:)
    TYPE (dnS_t),           intent(inout)  :: Func(:)
    integer,                intent(in)     :: nderiv

    integer :: i

    DO i=1,size(Func)
      Func(i) = ZERO
    END DO

  END SUBROUTINE EvalFunc_QML_Empty

  SUBROUTINE Write_QML_Empty(QModel,nio)
  IMPLICIT NONE

    CLASS (QML_Empty_t), intent(in)    :: QModel
    integer,             intent(in)    :: nio

    IF (QModel%In_a_Model) RETURN

    write(nio,*) 'Potential parameters are written just below'
    write(nio,*) 'Init:                      ',QModel%Init
    write(nio,*) 'In_a_Model:                ',QModel%In_a_Model

    write(nio,*) 'option:                    ',QModel%option
    write(nio,*)
    write(nio,*) 'nsurf:                     ',QModel%nsurf
    write(nio,*) 'ndim:                      ',QModel%ndim
    write(nio,*) 'numeric:                   ',QModel%numeric
    write(nio,*) 'no analytical derivatives: ',QModel%no_ana_der
    write(nio,*) 'Cartesian => model coord.: ',QModel%Cart_TO_Q
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
      CALL Write_RMat(QModel%d0GGdef,nio,nbcol1=5)
    END IF

    IF (allocated(QModel%Q0)) THEN
      write(nio,*) 'Reference Coordinate values, Q0(:)'
      CALL Write_RVec(QModel%Q0,nio,nbcol1=5)
    END IF
    flush(nio)


  END SUBROUTINE Write_QML_Empty
  SUBROUTINE Write0_QML_Empty(QModel,nio)
  IMPLICIT NONE

    CLASS (QML_Empty_t), intent(in)    :: QModel
    integer,              intent(in)    :: nio


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
     CALL Write_RMat(QModel%d0GGdef,nio,nbcol1=5)
    END IF

    write(nio,*) 'END QUANTUM MODEL default parameters'
    flush(nio)


  END SUBROUTINE Write0_QML_Empty

  SUBROUTINE Cart_TO_Q_QML_Empty(QModel,dnX,dnQ,nderiv)
  USE ADdnSVM_m, ONLY :  dnS_t
  IMPLICIT NONE

    CLASS(QML_Empty_t),      intent(in)    :: QModel
    TYPE (dnS_t),            intent(in)    :: dnX(:,:)
    TYPE (dnS_t),            intent(inout) :: dnQ(:)
    integer,                 intent(in)    :: nderiv


    dnQ(:) = ZERO

  END SUBROUTINE Cart_TO_Q_QML_Empty
  SUBROUTINE Qact_TO_Q_QML_Empty(Qact,Q,list_act)
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
