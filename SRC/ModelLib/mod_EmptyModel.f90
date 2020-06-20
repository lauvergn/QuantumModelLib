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
MODULE mod_EmptyModel
  USE mod_Lib
  USE mod_dnS
  USE mod_dnMat
  IMPLICIT NONE

  TYPE :: EmptyModel_t
    integer :: nsurf       = 0
    integer :: ndim        = 0
    logical :: numeric     = .FALSE.
    logical :: no_ana_der  = .FALSE. ! to force numerical derivatives
                                     ! for potential without analitical derivatives

    logical :: adiabatic   = .TRUE.
    integer :: option      = 0
    logical :: PubliUnit   = .FALSE. ! when PubliUnit=.TRUE., the units of a reference (publi ...) are used. Default (atomic unit)

    character (len=:),  allocatable :: pot_name
    real (kind=Rkind),  allocatable :: d0GGdef(:,:)
    real (kind=Rkind),  allocatable :: Q0(:)

    !Vec0 must be allocatable, to be able to deallocate with deallocate of the EmptyModel_t variable.
    TYPE (dnMat_t),     allocatable :: Vec0 ! to get the correct phase of the adiatic couplings
    CONTAINS
      PROCEDURE :: Eval_QModel_Pot    => Eval_EmptyModel_Pot
      PROCEDURE :: Write_QModel       => Write_EmptyModel
      PROCEDURE :: Write0_QModel      => Write0_EmptyModel
     !PROCEDURE :: get2_Q0_QModel     => get2_Q0_EmptyModel
      PROCEDURE :: get_d0GGdef_QModel => get_d0GGdef_EmptyModel
  END TYPE EmptyModel_t

  INTERFACE get_Q0_QModel
    MODULE PROCEDURE get_Q0_EmptyModel
  END INTERFACE


CONTAINS

  FUNCTION Init_EmptyModel(QModel_in) RESULT(QModel)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (EmptyModel_t)                  :: QModel

    TYPE(EmptyModel_t),  intent(in)      :: QModel_in ! variable to transfer info to the init

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_EmptyModel'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF


    CALL Init0_EmptyModel(QModel,QModel_in)

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_EmptyModel
  SUBROUTINE Init0_EmptyModel(QModel,QModel_in)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (EmptyModel_t), intent(inout)   :: QModel

    TYPE(EmptyModel_t),  intent(in)      :: QModel_in ! variable to transfer info to the init

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init0_EmptyModel'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    !QModel = QModel_in   ! it does not work always with nagfor
    QModel%nsurf     = QModel_in%nsurf
    QModel%ndim      = QModel_in%ndim
    QModel%numeric   = QModel_in%numeric
    QModel%adiabatic = QModel_in%adiabatic
    QModel%option    = QModel_in%option
    QModel%PubliUnit = QModel_in%PubliUnit


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
    QModel%pot_name = 'emptymodel'

    IF (allocated(QModel%Vec0)) THEN
      CALL QML_dealloc_dnMat(QModel%Vec0)
      deallocate(QModel%Vec0)
    END If

    IF (allocated(QModel%d0GGdef)) deallocate(QModel%d0GGdef)
    IF (allocated(QModel%Q0))      deallocate(QModel%Q0)

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Init0_EmptyModel
  SUBROUTINE get2_Q0_EmptyModel(QModel,Q0)
    IMPLICIT NONE

    real (kind=Rkind),    intent(inout)  :: Q0(:)
    CLASS(EmptyModel_t),  intent(in)     :: QModel

    IF (size(Q0) /= QModel%ndim) THEN
      STOP 'STOP in get_Q0_EmptyModel, wrong ndim size.'
    END IF

    IF (allocated(QModel%Q0)) THEN
      Q0(:) =  QModel%Q0
    END IF


  END SUBROUTINE get2_Q0_EmptyModel

  SUBROUTINE get_Q0_EmptyModel(QModel,Q0,err_Q0)
    IMPLICIT NONE

    CLASS(EmptyModel_t),  intent(in)              :: QModel
    real (kind=Rkind),    intent(inout)           :: Q0(:)
    integer,              intent(inout), optional ::  err_Q0

    IF (size(Q0) /= QModel%ndim) THEN
      STOP 'STOP in get_Q0_EmptyModel, wrong ndim size.'
    END IF

    IF (allocated(QModel%Q0)) THEN
      Q0(:) =  QModel%Q0
      err_Q0 = 0
    ELSE
      err_Q0 = 1
    END IF

  END SUBROUTINE get_Q0_EmptyModel

  FUNCTION get_d0GGdef_EmptyModel(QModel) RESULT(d0GGdef)
    IMPLICIT NONE

    real (kind=Rkind),   allocatable               :: d0GGdef(:,:)
    CLASS(EmptyModel_t),             intent(in)    :: QModel

    integer :: i

    IF (allocated(d0GGdef)) deallocate(d0GGdef)

    IF (allocated(QModel%d0GGdef)) THEN
      d0GGdef =  QModel%d0GGdef
    ELSE
      allocate(d0GGdef(QModel%ndim,QModel%ndim))
      d0GGdef = ZERO
      DO i=1,QModel%ndim
        d0GGdef(i,i) = ONE
      END DO
    END IF

    !write(6,*) 'alloc Q0',allocated(Q0)

  END FUNCTION get_d0GGdef_EmptyModel
  SUBROUTINE Eval_EmptyModel_Pot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS (EmptyModel_t),   intent(in)     :: QModel
    TYPE (dnS_t),           intent(in)     :: dnQ(:)
    TYPE (dnS_t),           intent(inout)  :: Mat_OF_PotDia(:,:)
    integer,                intent(in)     :: nderiv


    Mat_OF_PotDia(1,1) = ZERO

  END SUBROUTINE Eval_EmptyModel_Pot

  SUBROUTINE Write_EmptyModel(QModel,nio)
  !USE mod_Lib
  IMPLICIT NONE

    CLASS (EmptyModel_t), intent(in)    :: QModel
    integer,              intent(in)    :: nio

    write(nio,*) 'nsurf:                     ',QModel%nsurf
    write(nio,*) 'ndim:                      ',QModel%ndim
    write(nio,*) 'numeric:                   ',QModel%numeric
    write(nio,*) 'adiabatic:                 ',QModel%adiabatic
    write(nio,*) 'no analitical derivatives: ',QModel%no_ana_der

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


  END SUBROUTINE Write_EmptyModel
  SUBROUTINE Write0_EmptyModel(QModel,nio)
  !USE mod_Lib
  IMPLICIT NONE

    CLASS (EmptyModel_t), intent(in)    :: QModel
    integer,              intent(in)    :: nio


    write(nio,*) 'QUANTUM MODEL default parameters'
    flush(nio)

    write(nio,*)
    write(nio,*) 'Potential parameters are written just below'
    write(nio,*)
    write(nio,*) 'nsurf:                     ',QModel%nsurf
    write(nio,*) 'ndim:                      ',QModel%ndim
    write(nio,*) 'numeric:                   ',QModel%numeric
    write(nio,*) 'adiabatic:                 ',QModel%adiabatic
    write(nio,*) 'no analitical derivatives: ',QModel%no_ana_der
    write(nio,*)

     IF (allocated(QModel%d0GGdef)) THEN
       write(nio,*) 'Deformation metric tensor (~ 1/Mii)'
       CALL Write_RMat(QModel%d0GGdef,nio,nbcol1=5)
     END IF

    write(nio,*) 'END QUANTUM MODEL default parameters'
    flush(nio)


  END SUBROUTINE Write0_EmptyModel
END MODULE mod_EmptyModel
