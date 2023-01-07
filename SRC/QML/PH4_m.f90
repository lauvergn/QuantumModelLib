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
!> @brief Module which makes the initialization, calculation of the PH4 potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 08/04/2021
!!
MODULE QML_PH4_m
  USE QDUtil_NumParameters_m, out_unitp => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE


  integer, parameter :: ndim      = 1 ! ????
  integer, parameter :: max_ndim  = 9 ! dimension of the system, here 9

  integer, parameter :: max_fit   = 1 + (max_ndim-1) + (max_ndim-1) + (max_ndim-1)**2 + 5
  integer, parameter :: max_nn    = 20
                                            !Qopt  1  2  3  4  5  6  7  8  9
  integer, parameter :: listQop_fit3(max_ndim) = [-1, 2, 3,-1, 5,-1, 7,-1,-1]

  character (len=*), parameter :: base_fit3_Ene1_fileName='InternalData/PH4/fit3/interEneMP2'
  character (len=*), parameter :: base_fit3_Ene2_fileName='InternalData/PH4/fit3/interEneCCSDT-F12'
  character (len=*), parameter :: base_fit3_Qopt_fileName='InternalData/PH4/fit3/interQ_'
  character (len=*), parameter :: base_fit3_grad_fileName='InternalData/PH4/fit3/interGrad_'
  character (len=*), parameter :: base_fit3_hess_fileName='InternalData/PH4/fit3/interHess_'
  character (len=*), parameter :: base_fit3_AnHar_fileName='InternalData/PH4/fit3/inter_anharCCSDT-F12_'

!> @brief Derived type in which the PH4 parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_PH4_t

    PRIVATE

    real (kind=Rkind)  :: F(max_nn,max_fit)   = ZERO
    integer            :: nn(0:ndim,max_fit)  = 0
    integer            :: nt(max_fit)         = 0
    integer            :: largest_nn          = 0
    integer            :: ifunc_Ene           = 0
    integer            :: ifunc_Qopt          = 0
    integer            :: ifunc_Grad          = 0
    integer            :: ifunc_Hess          = 0
    integer            :: ifunc_AnHar         = 0


    real (kind=Rkind)  :: a(max_fit)          = ZERO
    real (kind=Rkind)  :: b(max_fit)          = ZERO

    logical            :: file_exist(max_fit) = .FALSE.

    integer, allocatable :: iQopt_TO_ifunc(:)
    integer, allocatable :: iQgrad_TO_ifunc(:)
    integer, allocatable :: iQjQHess_TO_ifunc(:,:)
    integer, allocatable :: AnHar_TO_ifunc(:)

   CONTAINS
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_PH4
    PROCEDURE :: EvalFunc_QModel => EvalFunc_QML_PH4
    PROCEDURE :: Write_QModel     => Write_QML_PH4
    PROCEDURE :: Write0_QModel    => Write0_QML_PH4
  END TYPE QML_PH4_t

  PUBLIC :: QML_PH4_t,Init_QML_PH4

  CONTAINS
!> @brief Function which makes the initialization of the PH4 parameters.
!!
!! @param QModel             TYPE(QML_PH4_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_PH4(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat, TO_string
    USE QMLLib_UtilLib_m, ONLY : make_FileName, file_open2
    IMPLICIT NONE

    TYPE (QML_PH4_t)                             :: QModel ! RESULT

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    integer :: i,j,ifunc
    logical :: exist,read_ab
    character (len=:), allocatable  :: FileName

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_PH4'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    IF (debug) write(out_unitp,*) 'option',QModel%option

    IF (QModel%ndim == 0) THEN
      ! it means that ndim was not present in the CALL Init_Model().
      ! => ndim is set to the default value (1)
      QModel%ndim  = 1
    END IF

    !The value of QModel%ndim must be 1, 2 or 9
    IF (QModel%ndim /= 1 .AND. QModel%ndim /= 9 .AND. QModel%ndim /= 2 ) THEN
       write(out_unitp,*) 'Write_QModel'
       CALL QModel%Write_QModel(out_unitp)
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' ndim MUST equal to 1 or 2 or 9. ndim: ',QModel%ndim
       STOP 'ERROR in Init_QML_PH4: ndim MUST equal to 1 or 2 or 9'
    END IF

    QModel%nsurf    = 1
    QModel%pot_name = 'ph4'

    QModel%ndimFunc  = 1

    ! ene
    QModel%nb_Func    = 1
    QModel%ifunc_Ene  = QModel%nb_Func
    ! Qopt(i)
    QModel%ifunc_Qopt = QModel%nb_Func + 1
    QModel%nb_Func    = QModel%nb_Func + (max_ndim-1)
    ! grad(i)
    QModel%ifunc_Grad = QModel%nb_Func + 1
    QModel%nb_Func    = QModel%nb_Func + (max_ndim-1)
    ! hess(i,j)
    QModel%ifunc_Hess = QModel%nb_Func + 1
    QModel%nb_Func    = QModel%nb_Func + (max_ndim-1)**2

    IF (QModel%option == 5) THEN
      ! anharm on Rp
      QModel%ifunc_AnHar= QModel%nb_Func + 1
      QModel%nb_Func    = QModel%nb_Func + 5 ! DeltaE(DeltaRp) = sum_i=0,4 coef(s,i) DeltaRp^i
    END IF

    !write(out_unitp,*) ' ifunc_Ene ',QModel%ifunc_Ene
    !write(out_unitp,*) ' ifunc_Qopt',QModel%ifunc_Qopt
    !write(out_unitp,*) ' ifunc_Grad',QModel%ifunc_Grad
    !write(out_unitp,*) ' ifunc_Hess',QModel%ifunc_Hess

    allocate(QModel%iQopt_TO_ifunc(max_ndim))
    QModel%iQopt_TO_ifunc(:) = -1
    allocate(QModel%iQgrad_TO_ifunc(max_ndim))
    QModel%iQgrad_TO_ifunc(:) = -1
    allocate(QModel%iQjQHess_TO_ifunc(max_ndim,max_ndim))
    QModel%iQjQHess_TO_ifunc(:,:) = -1
    IF (QModel%option == 5) THEN
      allocate(QModel%AnHar_TO_ifunc(5))
      QModel%AnHar_TO_ifunc(:) = -1
    END IF

    ! read the parameters

    ! for the energy
    ifunc   = 1
    read_ab = .FALSE.
    SELECT CASE (QModel%option)
    CASE (3)
      FileName = make_FileName(base_fit3_Ene1_fileName)
    CASE (4,5) ! 4 harmonic along the path, 5 harmonic + anharmonic(Rp) along the path
      FileName = make_FileName(base_fit3_Ene2_fileName)
    CASE Default
      FileName = make_FileName(base_fit3_Ene2_fileName)
    END SELECT

    !write(out_unitp,*) i,'FileName: ',FileName ; flush(out_unitp)
    CALL QML_read_para4d(QModel%a(ifunc),QModel%b(ifunc),QModel%F(:,ifunc),     &
                         QModel%nn(:,ifunc),ndim,QModel%nt(ifunc),max_nn,       &
                         FileName,QModel%file_exist(ifunc),read_ab,print_info=debug)
    !write(out_unitp,*) i,'Read done' ; flush(out_unitp)
    IF ( .NOT. QModel%file_exist(ifunc)) STOP ' ERROR while reading PH4 energy parameters'

    QModel%largest_nn = QModel%nn(0,ifunc)

    ! for Qopt
    DO i=2,max_ndim
      ifunc = ifunc + 1
      QModel%iQopt_TO_ifunc(i) = ifunc

      IF (listQop_fit3(i) == -1) CYCLE

      SELECT CASE (QModel%option)
      CASE (3,4,5)
        FileName = make_FileName(base_fit3_Qopt_fileName // TO_string(i))
      CASE Default
        FileName = make_FileName(base_fit3_Qopt_fileName // TO_string(i))
      END SELECT

      read_ab = (i == 2)

      !write(out_unitp,*) i,'FileName: ',FileName ; flush(out_unitp)
      CALL QML_read_para4d(QModel%a(ifunc),QModel%b(ifunc),QModel%F(:,ifunc),   &
                           QModel%nn(:,ifunc),ndim,QModel%nt(ifunc),max_nn,     &
                           FileName,QModel%file_exist(ifunc),read_ab,print_info=debug)
      !write(out_unitp,*) i,'Read done' ; flush(out_unitp)

      IF ( .NOT. QModel%file_exist(ifunc)) STOP ' ERROR while reading PH4 Qop parameters'

       QModel%largest_nn = max(QModel%largest_nn,QModel%nn(0,ifunc))

    END DO

    ! for the gradient
    DO i=2,max_ndim
      ifunc = ifunc + 1
      QModel%iQgrad_TO_ifunc(i) = ifunc

      read_ab = .FALSE.

      SELECT CASE (QModel%option)
      CASE (3,4,5)
        FileName = make_FileName(base_fit3_grad_fileName // TO_string(i))
      CASE Default
        FileName = make_FileName(base_fit3_grad_fileName // TO_string(i))
      END SELECT

      !write(out_unitp,*) i,'FileName: ',FileName ; flush(out_unitp)
      CALL QML_read_para4d(QModel%a(ifunc),QModel%b(ifunc),QModel%F(:,ifunc),   &
                           QModel%nn(:,ifunc),ndim,QModel%nt(ifunc),max_nn,     &
                           FileName,QModel%file_exist(ifunc),read_ab,print_info=debug)
      !write(out_unitp,*) i,'Read done' ; flush(out_unitp)

       QModel%largest_nn = max(QModel%largest_nn,QModel%nn(0,ifunc))

    END DO

    ! !hess(i,j)
    DO i=2,max_ndim
    DO j=2,max_ndim
      ifunc = ifunc + 1
      QModel%iQjQHess_TO_ifunc(i,j) = ifunc

      read_ab = .FALSE.

      SELECT CASE (QModel%option)
      CASE (3,4,5)
        FileName = make_FileName(base_fit3_hess_fileName //                     &
                                       TO_string(i) // '_' // TO_string(j) )
      CASE Default
        FileName = make_FileName(base_fit3_hess_fileName //                     &
                                       TO_string(i) // '_' // TO_string(j) )
      END SELECT

      !write(out_unitp,*) i,'FileName: ',FileName ; flush(out_unitp)
      CALL QML_read_para4d(QModel%a(ifunc),QModel%b(ifunc),QModel%F(:,ifunc),   &
                           QModel%nn(:,ifunc),ndim,QModel%nt(ifunc),max_nn,     &
                           FileName,QModel%file_exist(ifunc),read_ab,print_info=debug)
      !write(out_unitp,*) i,'Read done' ; flush(out_unitp)

       QModel%largest_nn = max(QModel%largest_nn,QModel%nn(0,ifunc))

    END DO
    END DO

    IF (QModel%option == 5) THEN
      ! for the Anharmonic contribution along Rp
      DO i=1,5
        ifunc = ifunc + 1
        QModel%AnHar_TO_ifunc(i) = ifunc

        read_ab = .FALSE.

        FileName = make_FileName(base_fit3_AnHar_fileName // TO_string(i-1))

        !write(out_unitp,*) i,'FileName: ',FileName ; flush(out_unitp)
        CALL QML_read_para4d(QModel%a(ifunc),QModel%b(ifunc),QModel%F(:,ifunc), &
                             QModel%nn(:,ifunc),ndim,QModel%nt(ifunc),max_nn,   &
                             FileName,QModel%file_exist(ifunc),read_ab,print_info=debug)
        !write(out_unitp,*) i,'Read done' ; flush(out_unitp)

        QModel%largest_nn = max(QModel%largest_nn,QModel%nn(0,ifunc))

      END DO
    END IF

    deallocate(FileName)

    IF (debug) write(out_unitp,*) 'largest_nn',QModel%largest_nn

    IF (debug) write(out_unitp,*) 'init Q0 of PH4' ! for the rigid constraints

    SELECT CASE (QModel%option)
    CASE (3)
      QModel%Q0 = [0.5_Rkind,2.437585_Rkind,2.671068_Rkind,                     &
                   ZERO,0.810373_Rkind,ZERO,1.586462_Rkind,ZERO,PI]

    CASE Default
      QModel%Q0 = [0.5_Rkind,2.437585_Rkind,2.671068_Rkind,                     &
                   ZERO,0.810373_Rkind,ZERO,1.586462_Rkind,ZERO,PI]
    END SELECT

    IF (debug) write(out_unitp,*) 'init d0GGdef of PH4'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_PH4
!> @brief Subroutine wich prints the current QML_PH4 parameters.
!!
!! @param QModel            CLASS(QML_PH4_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_PH4(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_PH4_t),     intent(in) :: QModel
    integer,              intent(in) :: nio

  END SUBROUTINE Write_QML_PH4
!> @brief Subroutine wich prints the default QML_PH4 parameters.
!!
!! @param QModel            CLASS(QML_PH4_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write0_QML_PH4(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_PH4_t),     intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'PH4 parameters'
    write(nio,*)
    write(nio,*) ' H + H-PH2 -> H-H + PH2 approximate 9D-potential:'
    write(nio,*)
    write(nio,*) '   Quadratic potential along the reaction path'
    write(nio,*) '   Reaction coordinate: R- = 1/2(RCH1-RHH)'
    write(nio,*) '   Optimal coordinates along the path at MP2/cc-pVTZ'
    write(nio,*) '   V0 along the path at CCSD(T)-F12/cc-pVTZ-F12 (not yet)'
    write(nio,*) '   Hessian along the path at MP2/cc-pVTZ'
    write(nio,*)
    write(nio,*) 'One option is availble:'
    write(nio,*) '-option=3 (default)'
    write(nio,*) '  Coordinate transformations on some valence angles are added.'
    write(nio,*) '  x = tan(A-Pi/2), so that the tranformed angle range is ]-inf,+inf[.'
    write(nio,*)
    write(nio,*) 'Primitive coordinates (z-matrix):'
    write(nio,*) '   P'
    write(nio,*) '   H        1 RPH1'
    write(nio,*) '   X        1 RX3     2 Pi/2'
    write(nio,*) '   X        1 RX4     3 Pi/2   2 DX4(#7)'
    write(nio,*) '   H        1 RPH2    4 APH2   2 DH2'
    write(nio,*) '   H        1 RPH3    4 APH3   2 DH3'
    write(nio,*) '   X        2 RX7     1 Pi/2   4 ZERO'
    write(nio,*) '   H        2 RHH     7 AHH    1 DHH(#9)'
    write(nio,*)
    write(nio,*) 'If option=5 is selected, the angles (AHH #17) is transformed:'
    write(nio,*)     ' tAHH = tan(AHH-Pi/2) (#8)'
    write(nio,*)
    write(nio,*) 'Symetrization of some coordinates (linear combinations):'
    write(nio,*)
    write(nio,*) '#1      R-= 1/2 (RPH1 - RHH)'
    write(nio,*) '#2      R+= 1/2 (RPH1 + RHH)'
    write(nio,*)
    write(nio,*) '#3      RPH+= 1/2(RPH2 + RPH3)'
    write(nio,*) '#4      RPH-= 1/2(RPH2 - RPH3)'
    write(nio,*)
    write(nio,*) '#5      APH+= 1/2(APH2 + APH3)'
    write(nio,*) '(cte)   APH-= 1/2(APH2 - APH3)'
    write(nio,*)
    write(nio,*) '#6      DPH+= 1/2(DH2 + DH3)'
    write(nio,*) '(cte)   DPH-= 1/2(DH2 - DH3)'
    write(nio,*)
    write(nio,*) 'end PH4 parameters'

  END SUBROUTINE Write0_QML_PH4

!> @brief Subroutine wich calculates the PH4 potential with derivatives up to the 2d order.
!!
!! @param QModel             CLASS(QML_PH4_t):   derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_PH4(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_PH4_t),     intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable   :: dnPoly(:)
    TYPE (dnS_t)    :: dnDQ(QModel%ndim),vh,Rm,tRm
    integer         :: i1,i2,i,ifunc

    !write(out_unitp,*) 'coucou pot PH4 '

    Rm  = dnQ(1)
    tRm = tanh(Rm)

    allocate(dnPoly(QModel%largest_nn))
    DO i=1,QModel%largest_nn
        dnPoly(i) = dnLegendre0(tRm,i-1)
    END DO

    DO i1=2,QModel%ndim
      ifunc = QModel%iQopt_TO_ifunc(i1)
      dnDQ(i1) = dnQ(i1) - QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly)
      !write(out_unitp,*) 'DQ',i1,get_d0(dnDQ(i1))
    END DO

    vh = ZERO

    SELECT CASE(QModel%option)
    CASE(3,4)
      DO i1=2,QModel%ndim
        ifunc = QModel%iQgrad_TO_ifunc(i1)
        IF (QModel%file_exist(ifunc)) THEN
          vh = vh + dnDQ(i1) * QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly)
          !write(out_unitp,*) 'vh grad',i1,get_d0(vh)
        END IF

        ifunc = QModel%iQjQHess_TO_ifunc(i1,i1)
        IF (QModel%file_exist(ifunc)) THEN
          vh = vh + HALF * dnDQ(i1)**2 * QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly)
          !write(out_unitp,*) 'vh hessii',i1,get_d0(vh)
        END IF
      END DO
    CASE(5) ! with anharmonicity along Rp (Q(2))
      !write(out_unitp,*) 'coucou anhar'

      ! function in DeltaQ2 (poly fourth order)
      i1    = 2
      DO i=1,5
        ifunc = QModel%AnHar_TO_ifunc(i) ! DeltaQ2^(i-1)
        IF (QModel%file_exist(ifunc)) THEN
          vh = vh + dnDQ(i1)**(i-1) * QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly)
          !write(out_unitp,*) 'vh anhar',i1,i,get_d0(vh)
        END IF
      END DO

      DO i1=3,QModel%ndim ! i1=2 is treated before
        ifunc = QModel%iQgrad_TO_ifunc(i1)
        IF (QModel%file_exist(ifunc)) THEN
          vh = vh + dnDQ(i1) * QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly)
          !write(out_unitp,*) 'vh grad',i1,get_d0(vh)
        END IF

        ifunc = QModel%iQjQHess_TO_ifunc(i1,i1)
        IF (QModel%file_exist(ifunc)) THEN
          vh = vh + HALF * dnDQ(i1)**2 * QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly)
        END IF
        !write(out_unitp,*) 'vh hessii',i1,i1,get_d0(vh)
      END DO

    END SELECT

    ! off diagonal contribution
    DO i1=2,QModel%ndim
      DO i2=i1+1,QModel%ndim
        ifunc = QModel%iQjQHess_TO_ifunc(i1,i2)
        IF (.NOT. QModel%file_exist(ifunc)) CYCLE
        vh = vh + dnDQ(i1)*dnDQ(i2) * QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly)
        !write(out_unitp,*) 'vh hessij',i1,i2,get_d0(vh)
      END DO
    END DO
    ifunc = QModel%ifunc_Ene
    Mat_OF_PotDia(1,1) = QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly) + vh

    CALL dealloc_dnS(Rm)
    CALL dealloc_dnS(tRm)
    CALL dealloc_dnS(dnPoly)
    CALL dealloc_dnS(dnDQ)
    deallocate(dnPoly)

  END SUBROUTINE EvalPot_QML_PH4

  SUBROUTINE EvalFunc_QML_PH4(QModel,Func,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_PH4_t),     intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    TYPE (dnS_t)    :: Rm
    integer         :: i1,i2,ifunc

    CALL EvalFunc_QML_PH4_fit3(QModel,Func,dnQ,nderiv)

  END SUBROUTINE EvalFunc_QML_PH4

  SUBROUTINE EvalFunc_QML_PH4_fit3(QModel,Func,dnQ,nderiv)
    USE QDUtil_m,         ONLY : out_unitp => out_unit
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_PH4_t),     intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Func(:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable   :: dnPoly(:)
    TYPE (dnS_t)                :: Rm,tRm
    integer                     :: i,i1,i2,ifunc,ifunc2

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='EvalFunc_QML_PH4_fit3'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING EvalFunc_QML_PH4_fit3'
      write(out_unitp,*) 'Func',size(Func)
      DO i=1,size(Func)
        CALL Write_dnS(Func(i),nio=out_unitp,all_type=.TRUE.)
      END DO
      flush(out_unitp)
    END IF

    Rm  = dnQ(1)
    tRm = tanh(Rm)


    allocate(dnPoly(QModel%largest_nn))
    DO i=1,QModel%largest_nn
        dnPoly(i) = dnLegendre0(tRm,i-1)
    END DO

    DO ifunc=1,QModel%nb_Func
      Func(ifunc) = QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly)
    END DO

    CALL dealloc_dnS(Rm)
    CALL dealloc_dnS(tRm)
    CALL dealloc_dnS(dnPoly)
    deallocate(dnPoly)

    IF (debug) THEN
      DO i=1,size(Func)
        CALL Write_dnS(Func(i),nio=out_unitp,all_type=.TRUE.)
      END DO
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE EvalFunc_QML_PH4_fit3

  FUNCTION QML_dnvfour_fit3_WITH_poly(Rm,ifunc,QModel,dnPoly) RESULT(dnvfour)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t)                        :: dnvfour

    CLASS(QML_PH4_t),     intent(in)    :: QModel
    integer,              intent(in)    :: ifunc
    TYPE (dnS_t),         intent(in)    :: Rm
    TYPE (dnS_t),         intent(in)    :: dnPoly(:)

    integer                             :: iq,np

    !write(out_unitp,*) 'in QML_dnvfour_fit3_WITH_poly',ifunc
    !dnvfour = Rm
    dnvfour = ZERO
    IF (ifunc > max_fit) THEN
      write(out_unitp,*) ' ERROR in dnvfour'
      write(out_unitp,*) ' wrong value for ifunc',ifunc
      write(out_unitp,*) ' It must be smaller than max_fit:',max_fit
      STOP ' ERROR in QML_dnvfour_fit3_WITH_poly: wrong value for ifunc'
    END IF

    !CALL Write_dnS(dnvfour,nio=out_unitp,all_type=.TRUE.)


    np = QModel%nn(0,ifunc)

    ! special treatment for Qopt
    IF (ifunc >= QModel%ifunc_Qopt .AND. ifunc < QModel%ifunc_Grad) THEN
      iq = ifunc - QModel%ifunc_Qopt+2 ! index of the optimal coordinate
      IF (listQop_fit3(iq) == -1 .OR. np < 1) THEN ! Qopt is constant
        dnvfour = QModel%Q0(iq)
      ELSE
        dnvfour = dot_product(QModel%F(1:np,ifunc),dnPoly(1:np))
      END IF

      IF (iq == 2) dnvfour = dnvfour + QML_sc2_fit3(Rm,QModel%a(ifunc),QModel%b(ifunc))

    ELSE
      IF (np > 1) dnvfour = dot_product(QModel%F(1:np,ifunc),dnPoly(1:np))
    END IF
    !CALL Write_dnS(dnvfour,nio=out_unitp,all_type=.TRUE.)

  END FUNCTION QML_dnvfour_fit3_WITH_poly
  FUNCTION QML_sc2_fit3(x,a,b)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t)                          :: QML_sc2_fit3
    TYPE (dnS_t),           intent(in)    :: x
    real(kind=Rkind),       intent(in)    :: a,b

    QML_sc2_fit3 =  QML_dnSigmoid_PH4(-x,ONE)*(-x + a) +                        &
                    QML_dnSigmoid_PH4( x,ONE)*( x + b)

  END FUNCTION QML_sc2_fit3

  SUBROUTINE QML_read_para4d(a,b,F,n,ndim,nt,max_points,file_name,exist,read_ab,print_info)
    USE QMLLib_UtilLib_m, ONLY : file_open2
    IMPLICIT NONE

   integer,           intent(in)    :: max_points,ndim
   integer,           intent(inout) :: n(0:ndim),nt
   real (kind=Rkind), intent(inout) :: a,b,F(max_points)
   character (len=*), intent(in)    :: file_name
   logical,           intent(inout) :: exist
   logical,           intent(in)    :: read_ab
   logical,           intent(in)    :: print_info

   integer :: no,ios,kl,i

   IF (print_info) write(out_unitp,*) 'QML_read_para4d: file_name,max_points: ',file_name,max_points


   CALL file_open2(name_file=file_name,iunit=no,lformatted=.TRUE.,              &
                   old=.TRUE.,err_file=ios)
   IF (ios == 0) THEN

     read(no,*) i ! for nb_fit (not used)

     IF (print_info) write(out_unitp,*) 'file_name,nt,ndim: ',file_name,nt,ndim
     read(no,*) n(0:ndim)
     IF (print_info) write(out_unitp,*) 'file_name,n ',file_name,n(0:ndim)
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

     IF (read_ab) read(no,*) a,b

     CLOSE(no)
     exist = .TRUE.
   ELSE
     IF (print_info) write(out_unitp,*) 'The file (',file_name,') does not exist !!'
     exist = .FALSE.
   END IF

  END SUBROUTINE QML_read_para4d
  FUNCTION QML_dnSigmoid_PH4(x,a)
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t)                        :: QML_dnSigmoid_PH4

    TYPE (dnS_t),         intent(in)    :: x
    real(kind=Rkind),     INTENT(IN)    :: a

    QML_dnSigmoid_PH4 = (ONE+tanh(a*x))/TWO

  END FUNCTION QML_dnSigmoid_PH4
END MODULE QML_PH4_m
