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
MODULE mod_CH5_Model
  USE mod_NumParameters
  USE mod_EmptyModel
  IMPLICIT NONE

  PRIVATE


  integer, parameter :: ndim    = 1
  integer, parameter :: max_fit = 12
  integer, parameter :: max_nn  = 50
                                             !   1  2  3  4  5  6  7  8  9 10 11 12
  integer, parameter :: listQop_fit1(0:12) = [0,-1,-1, 3,-1,-1, 6,-1,-1,-1,-1,-1,-1]
  integer, parameter :: listQop_fit3(0:12) = [0,-1, 2, 3,-1,-1, 6,-1,-1,-1,-1,-1,-1]

  integer :: fit = 3
  character (len=*), parameter :: base_fit1_fileName='InternalData/CH5/fit1/inter12_'
  character (len=*), parameter :: base_fit3_fileName='InternalData/CH5/fit3/inter_'
  character (len=*), parameter :: base_fit3_hess_fileName='InternalData/CH5/fit3/inter_mp2_'


!> @brief Derived type in which the CH5 parameters are set-up.
  TYPE, EXTENDS (EmptyModel_t) ::  CH5_Model_t

    PRIVATE

    real (kind=Rkind)  :: F(max_nn,0:max_fit,0:max_fit)  = ZERO
    integer            :: nn(0:ndim,0:max_fit,0:max_fit) = 0
    integer            :: nt(0:max_fit,0:max_fit)        = 0

    real (kind=Rkind)  :: a(0:max_fit,0:max_fit)  = ZERO
    real (kind=Rkind)  :: b(0:max_fit,0:max_fit)  = ZERO

    logical            :: file_exist(0:max_fit,0:max_fit)  = .FALSE.


   CONTAINS
    PROCEDURE :: Eval_QModel_Pot  => eval_CH5_Pot
    PROCEDURE :: Eval_QModel_Func => Eval_CH5_Func
    PROCEDURE :: Write_QModel     => Write_CH5_Model
    PROCEDURE :: Write0_QModel    => Write0_CH5_Model
  END TYPE CH5_Model_t

  PUBLIC :: CH5_Model_t,Init_CH5_Model

  CONTAINS
!> @brief Function which makes the initialization of the CH5 parameters.
!!
!! @param QModel             TYPE(CH5_Model_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(EmptyModel_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_CH5_Model(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (CH5_Model_t)                           :: QModel ! RESULT

    TYPE(EmptyModel_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    integer :: ii,jj
    logical :: exist
    character (len=:), allocatable  :: FileName,base_fit_fileName

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_CH5_Model'
    !logical, parameter :: debug = .FALSE.
    logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_EmptyModel(QModel%EmptyModel_t,QModel_in)

    IF (QModel%ndim == 0) THEN
      ! it means that ndim was not present in the CALL Init_Model().
      ! => ndim is set to the default value (1)
      QModel%ndim  = 1
    END IF

    !The value of QModel%ndim must be 1 or 12
    IF (QModel%ndim /= 1 .AND. QModel%ndim /= 12) THEN
       write(out_unitp,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unitp)
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' ndim MUST equal to 1 or 12. ndim: ',QModel%ndim
       STOP 'ERROR in Init_CH5_Model: ndim MUST equal to 1 or 9'
    END IF

    QModel%nsurf    = 1
    QModel%pot_name = 'ch5'

    QModel%ndimFunc = 1

                  !  ene   Qopt(i)       hess(i,j)
    QModel%nb_Func  = 1 + (max_fit-1) + (max_fit-1)**2

    SELECT CASE (fit)
    CASE(1)
      base_fit_fileName = base_fit1_fileName
                   !  ene   Qopt(i)
      QModel%nb_Func  = 1 + (max_fit-1)
    CASE (3)
      base_fit_fileName = base_fit3_fileName
    END SELECT

    ! read the parameters
    jj = 0
    ii = 0
    ! for the energy
    FileName = make_FileName(base_fit_fileName // int_TO_char(ii))
    !write(out_unitp,*) ii,'FileName: ',FileName ; flush(out_unitp)
    CALL read_para4d(QModel%a(ii,jj),QModel%b(ii,jj),QModel%F(:,ii,jj),QModel%nn(:,ii,jj),  &
                     ndim,QModel%nt(ii,jj),max_nn,FileName,QModel%file_exist(ii,jj))
    !write(out_unitp,*) ii,'Read done' ; flush(out_unitp)
    IF ( .NOT. QModel%file_exist(ii,jj)) STOP ' ERROR while reading CH5 parameters'

    DO ii=2,max_fit

      IF (fit == 1 .AND. listQop_fit1(ii) == -1) CYCLE
      IF (fit == 3 .AND. listQop_fit3(ii) == -1) CYCLE

      FileName = make_FileName(base_fit_fileName // int_TO_char(ii))

      !write(out_unitp,*) ii,'FileName: ',FileName ; flush(out_unitp)
      CALL read_para4d(QModel%a(ii,jj),QModel%b(ii,jj),QModel%F(:,ii,jj),QModel%nn(:,ii,jj),  &
                       ndim,QModel%nt(ii,jj),max_nn,FileName,QModel%file_exist(ii,jj))
      !write(out_unitp,*) ii,'Read done' ; flush(out_unitp)

      IF ( .NOT. QModel%file_exist(ii,jj)) STOP ' ERROR while reading CH5 parameters'

    END DO

    ! !hess(i,j)
    DO ii=2,max_fit
    DO jj=ii,max_fit
      FileName = make_FileName(base_fit3_hess_fileName //                       &
                               int_TO_char(ii) // '_' // int_TO_char(jj) )
      CALL read_para4d(QModel%a(ii,jj),QModel%b(ii,jj),QModel%F(:,ii,jj),QModel%nn(:,ii,jj),  &
                      ndim,QModel%nt(ii,jj),max_nn,FileName,QModel%file_exist(ii,jj))
    END DO
    END DO

    deallocate(FileName)
    deallocate(base_fit_fileName)

    IF (debug) write(out_unitp,*) 'init Q0 of CH5' ! for the rigid constraints
    QModel%Q0 = [0.5_Rkind,2.182480603843_Rkind,              &
                 2.069879989614_Rkind,ZERO,ZERO,              &
                 1.797743737992_Rkind,ZERO,ZERO,              &
                 TWO*PI/THREE,-TWO*PI/THREE,PI/TWO,PI/TWO]


    IF (debug) write(out_unitp,*) 'init d0GGdef of CH5'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_CH5_Model
!> @brief Subroutine wich prints the current CH5_Model parameters.
!!
!! @param QModel            CLASS(CH5_Model_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_CH5_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(CH5_Model_t),  intent(in) :: QModel
    integer,              intent(in) :: nio

  END SUBROUTINE Write_CH5_Model
!> @brief Subroutine wich prints the default CH5_Model parameters.
!!
!! @param QModel            CLASS(CH5_Model_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_CH5_Model(QModel,nio)
  IMPLICIT NONE

    CLASS(CH5_Model_t),  intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'CH5 default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end CH5 default parameters'


  END SUBROUTINE Write0_CH5_Model

!> @brief Subroutine wich calculates the CH5 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(CH5_Model_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE eval_CH5_Pot(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(CH5_Model_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv



    TYPE (dnS_t)    :: dnDQ(QModel%ndim),vh,Rm
    integer         :: i1,i2

    Rm  = dnQ(1)

    SELECT CASE (fit)
    CASE (1)
      Mat_OF_PotDia(1,1) = dnvfour_fit1(Rm,0,0,QModel)
    CASE(3)
      DO i1=2,QModel%ndim
        dnDQ(i1) = dnQ(i1) - dnvfour_fit3(Rm,i1,0,QModel)
      END DO

      vh = ZERO
      DO i1=2,QModel%ndim
        vh = vh + dnDQ(i1)*dnDQ(i1) * dnvfour_fit3(Rm,i1,i1,QModel)
        DO i2=i1+1,QModel%ndim
          IF (.NOT. QModel%file_exist(i1,i2)) CYCLE
          vh = vh + dnDQ(i1)*dnDQ(i2) * TWO*dnvfour_fit3(Rm,i1,i2,QModel)
        END DO
      END DO
      Mat_OF_PotDia(1,1) = dnvfour_fit3(Rm,0,0,QModel) + vh * HALF
    END SELECT

  END SUBROUTINE eval_CH5_Pot

  SUBROUTINE eval_CH5_Func(QModel,Func,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(CH5_Model_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    TYPE (dnS_t)    :: Rm
    integer         :: i1,i2,ifunc

    SELECT CASE (fit)
    CASE (1)
      CALL eval_CH5_Func_fit1(QModel,Func,dnQ,nderiv)
    CASE(3)
      CALL eval_CH5_Func_fit3(QModel,Func,dnQ,nderiv)
    END SELECT

  END SUBROUTINE eval_CH5_Func

  SUBROUTINE eval_CH5_Func_fit3(QModel,Func,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(CH5_Model_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv



    TYPE (dnS_t)    :: Rm
    integer         :: i1,i2,ifunc

    Rm    = dnQ(1)
    ifunc = 1

    ! energy
    Func(ifunc) = dnvfour_fit3(Rm,0,0,QModel)

    ! Qopt: 11 values
    DO i1=2,max_fit
      ifunc = ifunc + 1
      Func(ifunc) = dnvfour_fit3(Rm,i1,0,QModel)
    END DO


    ! hess_ij: xx values
    DO i1=2,max_fit
    DO i2=i1,max_fit
       IF ( .NOT. QModel%file_exist(i1,i2)) CYCLE
       ifunc = ifunc + 1
       Func(ifunc) = dnvfour_fit3(Rm,i1,i2,QModel)
    END DO
    END DO

  END SUBROUTINE eval_CH5_Func_fit3
  FUNCTION dnvfour_fit3(Rm,iq,jq,QModel) RESULT(dnvfour)
  USE mod_dnS
  USE mod_dnPoly
  IMPLICIT NONE

    TYPE (dnS_t)                        :: dnvfour

    CLASS(CH5_Model_t),   intent(in)    :: QModel
    integer,              intent(in)    :: iq,jq
    TYPE (dnS_t),         intent(in)    :: Rm

    integer         :: i,kl,iiq,jjq
    TYPE (dnS_t)    :: tRm,Rm0 ! transformation of Rm

    !write(6,*) 'in dnvfour_fit3',iq,jq

    IF (iq > max_fit .OR. iq < 0 .OR. jq > max_fit .OR. jq < 0) THEN
      write(out_unitp,*) ' ERROR in dnvfour'
      write(out_unitp,*) ' wrong value for iq or jq',iq,jq
      STOP ' ERROR in dnvfour_fit3: wrong value for iq or jq'
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
        dnvfour = dnvfour + sc2_fit3(Rm,QModel%a(iq,0),QModel%b(iq,0))
      ELSE
        dnvfour = dnvfour + sc_fit3(Rm,QModel%a(iq,0),QModel%b(iq,0))
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
      dnvfour = dnvfour + sc_fit3(Rm,QModel%a(iiq,iiq),QModel%b(iiq,iiq))

      CALL QML_dealloc_dnS(tRm)
    ELSE IF (jq > 0 .AND. iq > 0 .AND. iq < jq) THEN
      tRm = tanh(Rm)
      iiq = iq
      jjq = jq

      dnvfour = ZERO
      DO i=1,QModel%nn(0,iiq,jjq)
        dnvfour = dnvfour + QModel%F(i,iiq,jjq)*QML_dnLegendre0(tRm,i-1)
      END DO
      dnvfour = dnvfour + sc_fit3(Rm,QModel%a(iiq,jjq),QModel%b(iiq,jjq))

      CALL QML_dealloc_dnS(tRm)
    END IF

  END FUNCTION dnvfour_fit3
  FUNCTION sc2_fit3(x,a,b)
    TYPE (dnS_t)                          :: sc2_fit3
    TYPE (dnS_t),           intent(in)    :: x
    real(kind=Rkind),       intent(in)    :: a,b

    sc2_fit3 = dnSigmoid(-x,ONE)*(-x+a) + dnSigmoid(x,ONE)*(x+b)

  END FUNCTION sc2_fit3
  FUNCTION sc_fit3(x,a,b)
    TYPE (dnS_t)                          :: sc_fit3
    TYPE (dnS_t),           intent(in)    :: x
    real(kind=Rkind),       intent(in)    :: a,b

    sc_fit3 = a + b*tanh(x)

  END FUNCTION sc_fit3
  SUBROUTINE eval_CH5_Func_fit1(QModel,Func,dnQ,nderiv)
  USE mod_dnS
  IMPLICIT NONE

    CLASS(CH5_Model_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv



    TYPE (dnS_t)    :: Rm
    integer         :: i1,i2,ifunc

    Rm    = dnQ(1)
    ifunc = 1

    Func(ifunc) = dnvfour_fit1(Rm,0,0,QModel)

    DO i1=1,12
      ifunc = ifunc + 1
      Func(ifunc) = dnvfour_fit1(Rm,i1,0,QModel)
    END DO

    ! DO i1=1,QModel%ndim-1
    ! DO i2=1,QModel%ndim-1
    !   ifunc = ifunc + 1
    !   Func(ifunc) = dnvfour_fit1(Rm,i1,i2,QModel)
    ! END DO
    ! END DO

  END SUBROUTINE eval_CH5_Func_fit1

  FUNCTION dnvfour_fit1(Rm,iq,jq,QModel) RESULT(dnvfour)
  USE mod_dnS
  USE mod_dnPoly
  IMPLICIT NONE

    TYPE (dnS_t)                        :: dnvfour

    CLASS(CH5_Model_t),   intent(in)    :: QModel
    integer,              intent(in)    :: iq,jq
    TYPE (dnS_t),         intent(in)    :: Rm

    integer         :: i,kl
    TYPE (dnS_t)    :: tRm,Rm0 ! transformation of Rm

    IF (iq > max_fit .OR. iq < 0 .OR. jq > max_fit .OR. jq < 0) THEN
      write(out_unitp,*) ' ERROR in dnvfour'
      write(out_unitp,*) ' wrong value for iq or jq E [0...8]',iq,jq
      STOP ' ERROR in dnvfour: wrong value for iq or jq E [0...8]'
    END IF

    IF (listQop_fit1(iq) == -1 .AND. jq == 0) THEN
      IF (iq == 1) THEN
        Rm0 = 0.3_Rkind
        dnvfour =               fm_fit1(Rm ,1.3_Rkind)+fp_fit1(Rm ,1.3_Rkind) + &
                  (2.15_Rkind - fm_fit1(Rm0,1.3_Rkind)-fp_fit1(Rm0,1.3_Rkind)) *&
                  exp(-0.9_Rkind*(Rm-0.3_Rkind)**2)
      ELSE
        dnvfour = QModel%Q0(iq)
      END IF
    ELSE
      tRm = tanh(Rm)

      dnvfour = ZERO
      DO i=1,QModel%nn(0,iq,jq)
        dnvfour = dnvfour + QModel%F(i,iq,jq)*QML_dnLegendre0(tRm,i-1)
      END DO

      CALL QML_dealloc_dnS(tRm)
    END IF
  END FUNCTION dnvfour_fit1
  SUBROUTINE read_para4d(a,b,F,n,ndim,nt,max_points,nom1,exist)
  IMPLICIT NONE

   integer,           intent(in)    :: max_points,ndim
   integer,           intent(inout) :: n(0:ndim),nt
   real (kind=Rkind), intent(inout) :: a,b,F(max_points)
   character (len=*), intent(in)    :: nom1
   logical,           intent(inout) :: exist

   integer :: no,ios,kl,i

   write(out_unitp,*) 'read_para4d: nom1,max_points: ',nom1,max_points


   CALL file_open2(name_file=nom1,iunit=no,lformatted=.TRUE.,           &
                   old=.TRUE.,err_file=ios)
   IF (ios == 0) THEN

     read(no,*) i ! for nb_fit (not used)

     write(out_unitp,*) 'nom1,nt,ndim: ',nom1,nt,ndim
     read(no,*) n(0:ndim)
     write(out_unitp,*) 'nom1,n ',nom1,n(0:ndim)
     IF (n(0) > max_points) THEN
         write(out_unitp,*) ' ERROR : The number of coefficients (',n(0),') >'
         write(out_unitp,*) '         than max_points (',max_points,')'
         write(out_unitp,*) '         STOP in read_para4d'
         STOP 'ERROR in read_para4d'
     END IF
     DO kl=1,n(0)
        read(no,*) F(kl)
!       write(out_unitp,*) F(kl)
     END DO

     read(no,*) a,b

     CLOSE(no)
     exist = .TRUE.
   ELSE
     write(out_unitp,*) 'The file (',nom1,') does not exist !!'
     exist = .FALSE.
   END IF

  END SUBROUTINE read_para4d
  FUNCTION dnSigmoid(x,a)
  USE mod_dnS
  IMPLICIT NONE

    TYPE (dnS_t)                        :: dnSigmoid

    TYPE (dnS_t),         intent(in)    :: x
    real(kind=Rkind),     INTENT(IN)    :: a

    dnSigmoid = (ONE+tanh(a*x))/TWO

  END FUNCTION dnSigmoid
  FUNCTION fp_fit1(x,a)
  USE mod_dnS
  IMPLICIT NONE

    TYPE (dnS_t)                        :: fp_fit1

    TYPE (dnS_t),         intent(in)    :: x
    real(kind=Rkind),     INTENT(IN)    :: a

    real(kind=Rkind), parameter    :: x0=0.3_Rkind
    real(kind=Rkind), parameter    :: b=6.401578173119_Rkind-FIVE

    fp_fit1 = dnSigmoid(x-x0,a)*(x+b)

  END FUNCTION fp_fit1
  FUNCTION fm_fit1(x,a)
  USE mod_dnS
  IMPLICIT NONE

    TYPE (dnS_t)                        :: fm_fit1

    TYPE (dnS_t),         intent(in)    :: x
    real(kind=Rkind),     INTENT(IN)    :: a

    real(kind=Rkind), parameter    :: x0=0.3_Rkind
    real(kind=Rkind), parameter    :: b=7.055647895899_Rkind-FIVE

    fm_fit1 = dnSigmoid(-x+x0,a)*(-x+b)

  END FUNCTION fm_fit1
END MODULE mod_CH5_Model
