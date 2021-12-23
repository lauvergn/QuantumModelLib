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

!> @brief Module which makes the initialization, calculation of the CH5 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_CH5_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE


  integer, parameter :: ndim    = 1
  integer, parameter :: max_fit = 12
  integer, parameter :: max_nn  = 20
                                             !   1  2  3  4  5  6  7  8  9 10 11 12
  integer, parameter :: listQop_fit3(0:12) = [0,-1, 2, 3,-1,-1, 6,-1,-1,-1,-1,-1,-1]

  character (len=*), parameter :: base_fit3_fileName='InternalData/CH5/fit3/inter_'
  character (len=*), parameter :: base_fit3_hess_fileName='InternalData/CH5/fit3/inter_'

  character (len=*), parameter :: base_fit4_fileName='InternalData/CH5/fit4/inter_'
  character (len=*), parameter :: base_fit4_hess_fileName='InternalData/CH5/fit4/inter_mp2_'

  character (len=*), parameter :: base_fit5_fileName='InternalData/CH5/fit5/inter_'
  character (len=*), parameter :: base_fit5_hess_fileName='InternalData/CH5/fit5/inter_mp2_'

!> @brief Derived type in which the CH5 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_CH5_t

    PRIVATE

    real (kind=Rkind)  :: F(max_nn,0:max_fit,0:max_fit)  = ZERO
    integer            :: nn(0:ndim,0:max_fit,0:max_fit) = 0
    integer            :: nt(0:max_fit,0:max_fit)        = 0
    INTEGER            :: largest_nn = 0

    real (kind=Rkind)  :: a(0:max_fit,0:max_fit)  = ZERO
    real (kind=Rkind)  :: b(0:max_fit,0:max_fit)  = ZERO

    logical            :: file_exist(0:max_fit,0:max_fit)  = .FALSE.

    integer, allocatable :: ifunc_TO_i1i2(:,:)
    integer, allocatable :: i1i2_TO_ifunc(:,:)


   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_CH5
    PROCEDURE :: Eval_QModel_Func => EvalFunc_QML_CH5
    PROCEDURE :: Write_QModel     => Write_QML_CH5
    PROCEDURE :: Write0_QModel    => Write0_QML_CH5
  END TYPE QML_CH5_t

  PUBLIC :: QML_CH5_t,Init_QML_CH5

  CONTAINS
!> @brief Function which makes the initialization of the CH5 parameters.
!!
!! @param QModel             TYPE(QML_CH5_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_CH5(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_CH5_t)                           :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    integer :: ii,jj,ifunc
    logical :: exist
    character (len=:), allocatable  :: FileName

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_CH5'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    IF (QModel%ndim == 0) THEN
      ! it means that ndim was not present in the CALL Init_Model().
      ! => ndim is set to the default value (1)
      QModel%ndim  = 1
    END IF

    !The value of QModel%ndim must be 1 or 12
    IF (QModel%ndim /= 1 .AND. QModel%ndim /= 12 .AND. QModel%ndim /= 2 ) THEN
       write(out_unitp,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unitp)
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' ndim MUST equal to 1 or 2 or 12. ndim: ',QModel%ndim
       STOP 'ERROR in Init_QML_CH5: ndim MUST equal to 1 or 2 or 12'
    END IF

    QModel%nsurf    = 1
    QModel%pot_name = 'ch5'

    QModel%ndimFunc = 1
                  !  ene   Qopt(i)       hess(i,j)
    QModel%nb_Func  = 1 + (max_fit-1) + (max_fit-1)**2
    allocate(QModel%ifunc_TO_i1i2(2,QModel%nb_Func))
    allocate(QModel%i1i2_TO_ifunc(0:max_fit,0:max_fit))
    QModel%ifunc_TO_i1i2(:,:) = 0
    QModel%i1i2_TO_ifunc(:,:) = 0

    ! read the parameters
    jj = 0
    ii = 0
    ! for the energy
    ifunc = 1

    SELECT CASE (QModel%option)
    CASE (3)
      FileName = make_FileName(base_fit3_fileName // int_TO_char(ii))
    CASE (4)
      FileName = make_FileName(base_fit4_fileName // int_TO_char(ii))
    CASE (5)
      FileName = make_FileName(base_fit5_fileName // int_TO_char(ii))
    CASE Default
      FileName = make_FileName(base_fit4_fileName // int_TO_char(ii))
    END SELECT

    !write(out_unitp,*) ii,'FileName: ',FileName ; flush(out_unitp)
    CALL QML_read_para4d(QModel%a(ii,jj),QModel%b(ii,jj),QModel%F(:,ii,jj),     &
                         QModel%nn(:,ii,jj),ndim,QModel%nt(ii,jj),max_nn,       &
                         FileName,QModel%file_exist(ii,jj),print_info=debug)
    !write(out_unitp,*) ii,'Read done' ; flush(out_unitp)
    IF ( .NOT. QModel%file_exist(ii,jj)) STOP ' ERROR while reading CH5 energy parameters'
    QModel%ifunc_TO_i1i2(:,ifunc) = [0,0]
    QModel%i1i2_TO_ifunc(0,0)     = ifunc

    QModel%largest_nn = QModel%nn(0,ii,jj)

    DO ii=2,max_fit
      ifunc = ifunc + 1
      QModel%ifunc_TO_i1i2(:,ifunc) = [ii,0]
      QModel%i1i2_TO_ifunc(ii,0)    = ifunc

      IF (listQop_fit3(ii) == -1) CYCLE

      SELECT CASE (QModel%option)
      CASE (3)
        FileName = make_FileName(base_fit3_fileName // int_TO_char(ii))
      CASE (4)
        FileName = make_FileName(base_fit4_fileName // int_TO_char(ii))
      CASE (5)
        FileName = make_FileName(base_fit5_fileName // int_TO_char(ii))
      CASE Default
        FileName = make_FileName(base_fit4_fileName // int_TO_char(ii))
      END SELECT

      !write(out_unitp,*) ii,'FileName: ',FileName ; flush(out_unitp)
      CALL QML_read_para4d(QModel%a(ii,jj),QModel%b(ii,jj),QModel%F(:,ii,jj),   &
                           QModel%nn(:,ii,jj),ndim,QModel%nt(ii,jj),max_nn,     &
                           FileName,QModel%file_exist(ii,jj),print_info=debug)
      !write(out_unitp,*) ii,'Read done' ; flush(out_unitp)

      IF ( .NOT. QModel%file_exist(ii,jj)) STOP ' ERROR while reading CH5 Qop parameters'

       QModel%largest_nn = max(QModel%largest_nn,QModel%nn(0,ii,jj))

    END DO

    ! !hess(i,j)
    DO ii=2,max_fit
    DO jj=2,max_fit
      ifunc = ifunc + 1
      QModel%i1i2_TO_ifunc(ii,jj)    = ifunc

      IF (jj < ii) THEN
        QModel%ifunc_TO_i1i2(:,ifunc) = [jj,ii]
        CYCLE
      ELSE
        QModel%ifunc_TO_i1i2(:,ifunc) = [ii,jj]
      END IF

      SELECT CASE (QModel%option)
      CASE (3)
        FileName = make_FileName(base_fit3_hess_fileName //                     &
                             int_TO_char(ii) // '_' // int_TO_char(jj) )
      CASE (4)
        FileName = make_FileName(base_fit4_hess_fileName //                     &
                             int_TO_char(ii) // '_' // int_TO_char(jj) )
      CASE (5)
        FileName = make_FileName(base_fit5_hess_fileName //                     &
                             int_TO_char(ii) // '_' // int_TO_char(jj) )
      CASE Default
        FileName = make_FileName(base_fit4_hess_fileName //                     &
                             int_TO_char(ii) // '_' // int_TO_char(jj) )
      END SELECT

      CALL QML_read_para4d(QModel%a(ii,jj),QModel%b(ii,jj),QModel%F(:,ii,jj),   &
                           QModel%nn(:,ii,jj),ndim,QModel%nt(ii,jj),max_nn,     &
                           FileName,QModel%file_exist(ii,jj),print_info=debug)

      !IF ( .NOT. QModel%file_exist(ii,jj)) STOP ' ERROR while reading CH5 hessian parameters'

       QModel%largest_nn = max(QModel%largest_nn,QModel%nn(0,ii,jj))

    END DO
    END DO

    deallocate(FileName)

    IF (debug) write(out_unitp,*) 'largest_nn',QModel%largest_nn

    IF (debug) write(out_unitp,*) 'init Q0 of CH5' ! for the rigid constraints

    SELECT CASE (QModel%option)
    CASE (3,4)
      QModel%Q0 = [0.5_Rkind,2.182480603843_Rkind,              &
                   2.069879989614_Rkind,ZERO,ZERO,              &
                   1.797743737992_Rkind,ZERO,ZERO,              &
                   TWO*PI/THREE,-TWO*PI/THREE,PI/TWO,PI/TWO]
    CASE (5)
      QModel%Q0 = [0.5_Rkind,2.182480603843_Rkind,              &
                   2.069879989614_Rkind,ZERO,ZERO,              &
                   tan(1.797743737992_Rkind-Pi/TWO),ZERO,ZERO,  &
                   TWO*PI/THREE,-TWO*PI/THREE,ZERO,PI/TWO]
    CASE Default
      QModel%Q0 = [0.5_Rkind,2.182480603843_Rkind,              &
                   2.069879989614_Rkind,ZERO,ZERO,              &
                   1.797743737992_Rkind,ZERO,ZERO,              &
                   TWO*PI/THREE,-TWO*PI/THREE,PI/TWO,PI/TWO]
    END SELECT

    IF (debug) write(out_unitp,*) 'init d0GGdef of CH5'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_CH5
!> @brief Subroutine wich prints the current QML_CH5 parameters.
!!
!! @param QModel            CLASS(QML_CH5_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_CH5(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_CH5_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

  END SUBROUTINE Write_QML_CH5
!> @brief Subroutine wich prints the default QML_CH5 parameters.
!!
!! @param QModel            CLASS(QML_CH5_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_CH5(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_CH5_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'CH5 parameters'
    write(nio,*)
    write(nio,*) ' H + H-CH3 -> H-H + CH3 approximate 12D-potential:'
    write(nio,*)
    write(nio,*) '   Quadratic potential along the reaction path'
    write(nio,*) '   Reaction coordinate: R- = 1/2(RCH1-RHH)'
    write(nio,*) '   Optimal coordinates along the path at CCSD(T)-F12/cc-pVTZ-F12'
    write(nio,*) '   V0 along the path at CCSD(T)-F12/cc-pVTZ-F12'
    write(nio,*) '   Hessian along the path at MP2/cc-pVDZ'
    write(nio,*)
    write(nio,*) 'Two options are availble, which differ from the coordinates'
    write(nio,*) '-option=4 (default)'
    write(nio,*) '-option=5:'
    write(nio,*) '  Coordinate transformations on some valence angles are added.'
    write(nio,*) '  x = tan(A-Pi/2), so that the tranformed angle range is ]-inf,+inf[.'
    write(nio,*)
    write(nio,*) 'Primitive coordinates (z-matrix):'
    write(nio,*) '   C'
    write(nio,*) '   H 1 RCH1'
    write(nio,*) '   H 1 RCH2 2 ACH2'
    write(nio,*) '   H 1 RCH3 2 ACH3 3 DH3'
    write(nio,*) '   H 1 RCH4 2 ACH4 3 DH4'
    write(nio,*) '   X 2 1.   1 Pi/2 3 0.'
    write(nio,*) '   X 2 1.   1 Pi/2 3 Pi/2'
    write(nio,*) '   H 2 RHH  7 A    6 D'
    write(nio,*)
    write(nio,*) 'If option=5 is selected, the angles (ACH2,ACH3,ACH4,A) are transformed:'
    write(nio,*)     ' x = tan(A-Pi/2)'
    write(nio,*)
    write(nio,*) 'Symetrization of some coordinates (linear combinations):'
    write(nio,*)
    write(nio,*) '   R-= 1/2 (RCH1 - RHH)'
    write(nio,*) '   R+= 1/2 (RCH1 + RHH)'
    write(nio,*)
    write(nio,*) '   Ra1= 1/3      (  RCH1 + RCH2 + RCH3)'
    write(nio,*) '   Re1= 1/sqrt(6)(2.RCH1 - RCH2 - RCH3)'
    write(nio,*) '   Re2= 1/sqrt(2)(         RCH2 - RCH3)'
    write(nio,*)
    write(nio,*) '   Aa1= 1/3      (  ACH1 + ACH2 + ACH3)'
    write(nio,*) '   Ae1= 1/sqrt(6)(2.ACH1 - ACH2 - ACH3)'
    write(nio,*) '   Ae2= 1/sqrt(2)(         ACH2 - ACH3)'
    write(nio,*) '   Remark: when option=5 is selected, ...'
    write(nio,*) '      ... the symetrization is performed on the transformed angle, x.'
    write(nio,*)
    write(nio,*) 'end CH5 parameters'

  END SUBROUTINE Write0_QML_CH5

!> @brief Subroutine wich calculates the CH5 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_CH5_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_CH5(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS(QML_CH5_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv



    TYPE (dnS_t)    :: dnDQ(QModel%ndim),vh,Rm
    integer         :: i1,i2

    Rm  = dnQ(1)

    DO i1=2,QModel%ndim
      dnDQ(i1) = dnQ(i1) - QML_dnvfour_fit3(Rm,i1,0,QModel)
    END DO

    vh = ZERO
    DO i1=2,QModel%ndim
      vh = vh + dnDQ(i1)*dnDQ(i1) * QML_dnvfour_fit3(Rm,i1,i1,QModel)
      DO i2=i1+1,QModel%ndim
        IF (.NOT. QModel%file_exist(i1,i2)) CYCLE
        vh = vh + dnDQ(i1)*dnDQ(i2) * TWO*QML_dnvfour_fit3(Rm,i1,i2,QModel)
      END DO
    END DO
    Mat_OF_PotDia(1,1) = QML_dnvfour_fit3(Rm,0,0,QModel) + vh * HALF

  END SUBROUTINE EvalPot_QML_CH5

  SUBROUTINE EvalFunc_QML_CH5(QModel,Func,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    CLASS(QML_CH5_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    TYPE (dnS_t)    :: Rm
    integer         :: i1,i2,ifunc

    CALL EvalFunc_QML_CH5_fit3(QModel,Func,dnQ,nderiv)

  END SUBROUTINE EvalFunc_QML_CH5

  SUBROUTINE EvalFunc_QML_CH5_fit3(QModel,Func,dnQ,nderiv)
  USE QMLdnSVM_dnS_m
  USE QMLdnSVM_dnPoly_m
  IMPLICIT NONE

    CLASS(QML_CH5_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable   :: dnPoly(:)
    TYPE (dnS_t)                :: Rm,tRm
    integer                     :: i,i1,i2,ifunc,ifunc2

    Rm  = dnQ(1)
    tRm = tanh(Rm)


    allocate(dnPoly(QModel%largest_nn))
    DO i=1,QModel%largest_nn
        dnPoly(i) = QML_dnLegendre0(tRm,i-1)
    END DO

    ! energy
    ifunc = QModel%i1i2_TO_ifunc(0,0)
    Func(ifunc) = QML_dnvfour_fit3_WITH_poly(Rm,0,0,QModel,dnPoly)

    ! Qopt: 11 values
    DO i1=2,max_fit
      ifunc = QModel%i1i2_TO_ifunc(i1,0)
      Func(ifunc) = QML_dnvfour_fit3_WITH_poly(Rm,i1,0,QModel,dnPoly)
    END DO


    ! hess_ij: xx values
    DO i1=2,max_fit
      ifunc = QModel%i1i2_TO_ifunc(i1,i1)
      Func(ifunc) = QML_dnvfour_fit3_WITH_poly(Rm,i1,i1,QModel,dnPoly)
    END DO

    DO i1=2,max_fit
    DO i2=i1+1,max_fit
      ifunc = QModel%i1i2_TO_ifunc(i1,i2)
      Func(ifunc) = QML_dnvfour_fit3_WITH_poly(Rm,i1,i2,QModel,dnPoly)
      ifunc2 = QModel%i1i2_TO_ifunc(i2,i1)
      Func(ifunc2) = Func(ifunc)
    END DO
    END DO

    CALL QML_dealloc_dnS(Rm)
    CALL QML_dealloc_dnS(tRm)
    CALL QML_dealloc_dnS(dnPoly)
    deallocate(dnPoly)


  END SUBROUTINE EvalFunc_QML_CH5_fit3
  FUNCTION QML_dnvfour_fit3(Rm,iq,jq,QModel) RESULT(dnvfour)
  USE QMLdnSVM_dnS_m
  USE QMLdnSVM_dnPoly_m
  IMPLICIT NONE

    TYPE (dnS_t)                        :: dnvfour

    CLASS(QML_CH5_t),   intent(in)    :: QModel
    integer,              intent(in)    :: iq,jq
    TYPE (dnS_t),         intent(in)    :: Rm

    integer         :: i,kl,iiq,jjq
    TYPE (dnS_t)    :: tRm ! transformation of Rm

    !write(6,*) 'in QML_dnvfour_fit3',iq,jq

    IF (iq > max_fit .OR. iq < 0 .OR. jq > max_fit .OR. jq < 0) THEN
      write(out_unitp,*) ' ERROR in dnvfour'
      write(out_unitp,*) ' wrong value for iq or jq',iq,jq
      STOP ' ERROR in QML_dnvfour_fit3: wrong value for iq or jq'
    END IF

    IF (listQop_fit3(iq) == -1 .AND. jq == 0) THEN
      dnvfour = QModel%Q0(iq)
    ELSE
      tRm = tanh(Rm)

      dnvfour = ZERO
      DO i=1,QModel%nn(0,iq,jq)
        dnvfour = dnvfour + QModel%F(i,iq,jq)*QML_dnLegendre0(tRm,i-1)
      END DO

      CALL QML_dealloc_dnS(tRm)

    END IF

    IF (jq == 0  .AND. listQop_fit3(iq) /= -1) THEN
      IF (iq == 2) THEN
        dnvfour = dnvfour + QML_sc2_fit3(Rm,QModel%a(iq,0),QModel%b(iq,0))
      ELSE
        dnvfour = dnvfour + QML_sc_fit3(Rm,QModel%a(iq,0),QModel%b(iq,0))
      END IF
    END IF

    IF (jq > 0 .AND. iq == jq) THEN
      iiq = iq
      IF (iq == 4)  iiq = 5
      IF (iq == 7)  iiq = 8
      IF (iq == 9)  iiq = 10
      IF (iq == 11) iiq = 12

      tRm = tanh(Rm)

      dnvfour = ZERO
      DO i=1,QModel%nn(0,iiq,iiq)
        dnvfour = dnvfour + QModel%F(i,iiq,iiq)*QML_dnLegendre0(tRm,i-1)
      END DO
      dnvfour = dnvfour + QML_sc_fit3(Rm,QModel%a(iiq,iiq),QModel%b(iiq,iiq))

      CALL QML_dealloc_dnS(tRm)
    ELSE IF (jq > 0 .AND. iq > 0 .AND. iq < jq) THEN
      tRm = tanh(Rm)
      iiq = iq
      jjq = jq

      dnvfour = ZERO
      DO i=1,QModel%nn(0,iiq,jjq)
        dnvfour = dnvfour + QModel%F(i,iiq,jjq)*QML_dnLegendre0(tRm,i-1)
      END DO
      dnvfour = dnvfour + QML_sc_fit3(Rm,QModel%a(iiq,jjq),QModel%b(iiq,jjq))

      CALL QML_dealloc_dnS(tRm)
    END IF

    IF (iq == 6 .AND. jq == 0 .AND. QModel%option == 5) THEN
      dnvfour = tan(dnvfour - pi/TWO)
    END IF

  END FUNCTION QML_dnvfour_fit3
  FUNCTION QML_dnvfour_fit3_WITH_poly(Rm,iq,jq,QModel,dnPoly) RESULT(dnvfour)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    TYPE (dnS_t)                        :: dnvfour

    CLASS(QML_CH5_t),   intent(in)    :: QModel
    integer,              intent(in)    :: iq,jq
    TYPE (dnS_t),         intent(in)    :: Rm
    TYPE (dnS_t),         intent(in)    :: dnPoly(:)

    integer         :: i,kl,iiq,jjq,np

    !write(6,*) 'in QML_dnvfour_fit3',iq,jq
    dnvfour = ZERO
    IF (iq > max_fit .OR. iq < 0 .OR. jq > max_fit .OR. jq < 0) THEN
      write(out_unitp,*) ' ERROR in dnvfour'
      write(out_unitp,*) ' wrong value for iq or jq',iq,jq
      STOP ' ERROR in QML_dnvfour_fit3: wrong value for iq or jq'
    END IF

    IF (listQop_fit3(iq) == -1 .AND. jq == 0) THEN
      dnvfour = QModel%Q0(iq)
    ELSE IF (jq == 0  .AND. listQop_fit3(iq) /= -1) THEN

      np = QModel%nn(0,iq,jq)
      IF (np > 0) dnvfour = dot_product(QModel%F(1:np,iq,jq),dnPoly(1:np))

    END IF

    IF (jq == 0  .AND. listQop_fit3(iq) /= -1) THEN
      IF (iq == 2) THEN
        dnvfour = dnvfour + QML_sc2_fit3(Rm,QModel%a(iq,0),QModel%b(iq,0))
      ELSE
        dnvfour = dnvfour + QML_sc_fit3(Rm,QModel%a(iq,0),QModel%b(iq,0))
      END IF
    END IF

    IF (jq > 0 .AND. iq == jq) THEN
      iiq = iq
      IF (iq == 4)  iiq = 5
      IF (iq == 7)  iiq = 8
      IF (iq == 9)  iiq = 10
      IF (iq == 11) iiq = 12

      np = QModel%nn(0,iiq,iiq)

      IF (np > 0) dnvfour = dot_product(QModel%F(1:np,iiq,iiq),dnPoly(1:np))
      dnvfour = dnvfour + QML_sc_fit3(Rm,QModel%a(iiq,iiq),QModel%b(iiq,iiq))

    ELSE IF (jq > 0 .AND. iq > 0 .AND. iq < jq) THEN
      iiq = iq
      jjq = jq


      np = QModel%nn(0,iiq,jjq)

      IF (np > 0) dnvfour = dot_product(QModel%F(1:np,iiq,jjq),dnPoly(1:np))
      dnvfour = dnvfour + QML_sc_fit3(Rm,QModel%a(iiq,jjq),QModel%b(iiq,jjq))

    END IF

    IF (iq == 6 .AND. jq == 0 .AND. QModel%option == 5) THEN
      dnvfour = tan(dnvfour - pi/TWO)
    END IF

  END FUNCTION QML_dnvfour_fit3_WITH_poly
  FUNCTION QML_sc2_fit3(x,a,b)
    TYPE (dnS_t)                          :: QML_sc2_fit3
    TYPE (dnS_t),           intent(in)    :: x
    real(kind=Rkind),       intent(in)    :: a,b

    QML_sc2_fit3 = QML_dnSigmoid_CH5(-x,ONE)*(-x+a) + QML_dnSigmoid_CH5(x,ONE)*(x+b)

  END FUNCTION QML_sc2_fit3
  FUNCTION QML_sc_fit3(x,a,b)
    TYPE (dnS_t)                          :: QML_sc_fit3
    TYPE (dnS_t),           intent(in)    :: x
    real(kind=Rkind),       intent(in)    :: a,b

    QML_sc_fit3 = a + b*tanh(x)

  END FUNCTION QML_sc_fit3

  SUBROUTINE QML_read_para4d(a,b,F,n,ndim,nt,max_points,nom1,exist,print_info)
  IMPLICIT NONE

   integer,           intent(in)    :: max_points,ndim
   integer,           intent(inout) :: n(0:ndim),nt
   real (kind=Rkind), intent(inout) :: a,b,F(max_points)
   character (len=*), intent(in)    :: nom1
   logical,           intent(inout) :: exist
   logical,           intent(in)    :: print_info

   integer :: no,ios,kl,i

   IF (print_info) write(out_unitp,*) 'QML_read_para4d: nom1,max_points: ',nom1,max_points


   CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,                   &
                   old=.TRUE.,err_file=ios)
   IF (ios == 0) THEN

     read(no,*) i ! for nb_fit (not used)

     IF (print_info) write(out_unitp,*) 'nom1,nt,ndim: ',nom1,nt,ndim
     read(no,*) n(0:ndim)
     IF (print_info) write(out_unitp,*) 'nom1,n ',nom1,n(0:ndim)
     IF (n(0) > max_points) THEN
         write(out_unitp,*) ' ERROR : The number of coefficients (',n(0),') >'
         write(out_unitp,*) '         than max_points (',max_points,')'
         write(out_unitp,*) '         STOP in QML_read_para4d'
         STOP 'ERROR in QML_read_para4d'
     END IF
     DO kl=1,n(0)
        read(no,*) F(kl)
!       write(out_unitp,*) F(kl)
     END DO

     read(no,*) a,b

     CLOSE(no)
     exist = .TRUE.
   ELSE
     IF (print_info) write(out_unitp,*) 'The file (',nom1,') does not exist !!'
     exist = .FALSE.
   END IF

  END SUBROUTINE QML_read_para4d
  FUNCTION QML_dnSigmoid_CH5(x,a)
  USE QMLdnSVM_dnS_m
  IMPLICIT NONE

    TYPE (dnS_t)                        :: QML_dnSigmoid_CH5

    TYPE (dnS_t),         intent(in)    :: x
    real(kind=Rkind),     INTENT(IN)    :: a

    QML_dnSigmoid_CH5 = (ONE+tanh(a*x))/TWO

  END FUNCTION QML_dnSigmoid_CH5
END MODULE QML_CH5_m
