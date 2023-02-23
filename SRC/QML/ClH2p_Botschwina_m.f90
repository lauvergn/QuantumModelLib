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

!> @brief Module which makes the initialization, calculation of the ClH2p_Botschwina potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_ClH2p_Botschwina_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the ClH2p_Botschwina parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

    real(kind=Rkind), parameter :: a0 = 0.52917720835354106_Rkind


    real(kind=Rkind), parameter :: Req_CEPA1 = 1.30207_Rkind / a0
    real(kind=Rkind), parameter :: Aeq_CEPA1 = 94.244_Rkind  /180._Rkind*pi

    real(kind=Rkind), parameter :: FitR_CEPA1(10) = [                           &
       ZERO,0.141206_Rkind, -0.141162_Rkind, 0.095262_Rkind, -0.052840_Rkind,   &
       0.022631_Rkind, -0.006592_Rkind, 0.001370_Rkind, -0.000288_Rkind, 0.000044_Rkind]
    real(kind=Rkind), parameter :: FitA_CEPA1(10) = [                           &
       ZERO,0.086473_Rkind, -0.018289_Rkind, 0.002232_Rkind, 0.016749_Rkind,    &
       -0.051448_Rkind, 0.055196_Rkind, -0.035462_Rkind, 0.01236_Rkind, -0.001703_Rkind]

    real(kind=Rkind), parameter :: FitRRA_CEPA1(16) = [&
        0.001432_Rkind,  0.005146_Rkind, -0.000339_Rkind, -0.002293_Rkind,      &
       -0.009195_Rkind, -0.027619_Rkind,  0.000191_Rkind,  0.000449_Rkind,      &
        0.000964_Rkind, -0.003548_Rkind,  0.005228_Rkind,  0.028442_Rkind,      &
        0.013751_Rkind,  0.009600_Rkind,  0.013060_Rkind, -0.006492_Rkind]

        !rr’     ra      r2r'    r2a     rr’a    ra2     r3r’    r3a
        !r2r'2   r2a     r2r’a   rr’a2   ra3     ra4     r2a3    rr‘a3
    integer, parameter :: Fit_iRiRia(3,16) = reshape([                          &
         1,1,0,  1,0,1,  2,1,0,  2,0,1,  1,1,1,  1,0,2,  3,1,0,  3,0,1,         &
         2,2,0,  2,0,1,  2,1,1,  1,1,2,  1,0,3,  1,0,4,  2,0,3,  1,1,3],shape=[3,16])


         real(kind=Rkind), parameter :: Req_CEPA1_corr = 1.3030_Rkind / a0
         real(kind=Rkind), parameter :: Aeq_CEPA1_corr = 94.244_Rkind  /180._Rkind*pi

         real(kind=Rkind), parameter :: FitR_CEPA1_corr(10) = [                 &
          ZERO, 0.140500_Rkind, -0.140494_Rkind, 0.094799_Rkind, -0.052602_Rkind,&
          0.022550_Rkind, -0.006572_Rkind, 0.001366_Rkind, -0.000288_Rkind, 0.000044_Rkind]
         real(kind=Rkind), parameter :: FitA_CEPA1_corr(10) = [                 &
            ZERO,0.083509_Rkind,-0.018241_Rkind,0.002266_Rkind, 0.016749_Rkind, &
           -0.051448_Rkind, 0.055196_Rkind, -0.035462_Rkind,0.012362_Rkind, -0.001703_Rkind]

         real(kind=Rkind), parameter :: FitRRA_CEPA1_corr(16) = [               &
             0.001475_Rkind,  0.005122_Rkind, -0.000335_Rkind, -0.002281_Rkind, &
            -0.009159_Rkind, -0.027582_Rkind,  0.000191_Rkind,  0.000449_Rkind, &
             0.000964_Rkind, -0.003548_Rkind,  0.005228_Rkind,  0.028442_Rkind, &
             0.013786_Rkind,  0.009600_Rkind,  0.013060_Rkind, -0.006492_Rkind]


  TYPE, EXTENDS (QML_Empty_t) ::  QML_ClH2p_Botschwina_t

     PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

     integer                       :: nb_funcModel
     real(kind=Rkind), allocatable :: F(:)
     integer, allocatable          :: tab_func(:,:)
     real(kind=Rkind)              :: E0 = ZERO

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_ClH2p_Botschwina
    PROCEDURE :: Write_QModel    => Write_QML_ClH2p_Botschwina
    PROCEDURE :: Write0_QModel   => Write_QML_ClH2p_Botschwina
  END TYPE QML_ClH2p_Botschwina_t

  PUBLIC :: QML_ClH2p_Botschwina_t,Init_QML_ClH2p_Botschwina

  CONTAINS
!> @brief Subroutine which makes the initialization of the ClH2p_Botschwina parameters.
!!
!! @param ClH2p_BotschwinaPot  TYPE(QML_ClH2p_Botschwina_t):   derived type in which the parameters are set-up.
!! @param option               integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                  integer (optional): file unit to read the parameters.
!! @param read_param           logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_ClH2p_Botschwina(QModel_in,read_param,nio_param_file) RESULT(QModel)
    IMPLICIT NONE

    TYPE (QML_ClH2p_Botschwina_t)                :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !local variable
    integer                         :: i,ifunc,i1,i2,i3


    !-----------------------------------------------------------------

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_ClH2p_Botschwina'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 1
    QModel%ndim     = 3
    QModel%pot_name = 'ClH2p_Botschwina'


    IF (QModel%option < 1 .OR. QModel%option > 2) QModel%option = 1

    QModel%nb_funcModel = 2*size(FitR_CEPA1) + size(FitA_CEPA1)
    DO i=1,size(Fit_iRiRia,dim=2)
      i1 = Fit_iRiRia(1,i)
      i2 = Fit_iRiRia(2,i)
      QModel%nb_funcModel = QModel%nb_funcModel + 1
      IF (i1 /= i2) QModel%nb_funcModel = QModel%nb_funcModel + 1
    END DO
    write(out_unit,*) 'QModel%nb_funcModel',QModel%nb_funcModel

    allocate(QModel%tab_func(3,QModel%nb_funcModel))
    QModel%tab_func(:,:) = 0
    allocate(QModel%F(QModel%nb_funcModel))
    QModel%F(:) = ZERO

    SELECT CASE (QModel%option)
    CASE (1)

      QModel%d0GGdef = reshape(                                                 &
        [0.0005600083_Rkind,-0.0000011610_Rkind,-0.0000063582_Rkind,            &
        -0.0000011610_Rkind, 0.0005600083_Rkind,-0.0000063582_Rkind,            &
        -0.0000063582_Rkind,-0.0000063582_Rkind, 0.0001853777_Rkind],shape=[3,3])

     QModel%Qref = [Req_CEPA1,Req_CEPA1,Aeq_CEPA1]

     ifunc = 0
     ! coef for r^i and r'^i
     DO i=1,size(FitR_CEPA1)
       ifunc = ifunc + 1
       QModel%tab_func(:,ifunc)  = [i,0,0]
       QModel%F(ifunc)           = FitR_CEPA1(i)

       ifunc = ifunc + 1
       QModel%tab_func(:,ifunc)  = [0,i,0]
       QModel%F(ifunc)           = FitR_CEPA1(i)
     END DO

     ! coef for a^i
     DO i=1,size(FitA_CEPA1)
       ifunc = ifunc + 1
       QModel%tab_func(:,ifunc)  = [0,0,i]
       QModel%F(ifunc)           = FitA_CEPA1(i)
     END DO

     ! coef for r^i1 r'^i2 a^i3
     DO i=1,size(Fit_iRiRia,dim=2)
       i1 = Fit_iRiRia(1,i)
       i2 = Fit_iRiRia(2,i)
       i3 = Fit_iRiRia(3,i)

       ifunc = ifunc + 1
       QModel%tab_func(:,ifunc)  = [i1,i2,i3]
       QModel%F(ifunc)           = FitRRA_CEPA1(i)

       IF (i1 /= i2) THEN
         ifunc = ifunc + 1
         QModel%tab_func(:,ifunc)  = [i2,i1,i3]
         QModel%F(ifunc)           = FitRRA_CEPA1(i)
       END IF
     END DO

    CASE (2)

      QModel%d0GGdef = reshape(                                                 &
        [0.0005600083_Rkind,-0.0000011610_Rkind,-0.0000063536_Rkind,            &
        -0.0000011610_Rkind, 0.0005600083_Rkind,-0.0000063536_Rkind,            &
        -0.0000063536_Rkind,-0.0000063536_Rkind, 0.0001851131_Rkind],shape=[3,3])

        QModel%Qref = [Req_CEPA1_corr,Req_CEPA1_corr,Aeq_CEPA1_corr]

        ifunc = 0
        ! coef for r^i and r'^i
        DO i=1,size(FitR_CEPA1_corr)
          ifunc = ifunc + 1
          QModel%tab_func(:,ifunc)  = [i,0,0]
          QModel%F(ifunc)           = FitR_CEPA1_corr(i)

          ifunc = ifunc + 1
          QModel%tab_func(:,ifunc)  = [0,i,0]
          QModel%F(ifunc)           = FitR_CEPA1_corr(i)
        END DO

        ! coef for a^i
        DO i=1,size(FitA_CEPA1_corr)
          ifunc = ifunc + 1
          QModel%tab_func(:,ifunc)  = [0,0,i]
          QModel%F(ifunc)           = FitA_CEPA1_corr(i)
        END DO

        ! coef for r^i1 r'^i2 a^i3
        DO i=1,size(Fit_iRiRia,dim=2)
          i1 = Fit_iRiRia(1,i)
          i2 = Fit_iRiRia(2,i)
          i3 = Fit_iRiRia(3,i)

          ifunc = ifunc + 1
          QModel%tab_func(:,ifunc)  = [i1,i2,i3]
          QModel%F(ifunc)           = FitRRA_CEPA1_corr(i)

          IF (i1 /= i2) THEN
            ifunc = ifunc + 1
            QModel%tab_func(:,ifunc)  = [i2,i1,i3]
            QModel%F(ifunc)           = FitRRA_CEPA1_corr(i)
          END IF
        END DO


    CASE Default

      write(out_unit,*) ' ERROR in Init_QML_ClH2p_Botschwina '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1,2'
      STOP 'ERROR in Init_QML_ClH2p_Botschwina: wrong option'

    END SELECT
    flush(out_unit)



    IF (debug) write(out_unit,*) 'init Q0 of ClH2p_Botschwina'
    allocate(QModel%Q0(QModel%ndim))
    CALL get_Q0_QML_ClH2p_Botschwina(QModel%Q0,QModel,option=0)
    IF (debug) write(out_unit,*) 'QModel%Q0',QModel%Q0

    IF (debug) write(out_unit,*) 'init d0GGdef of ClH2p_Botschwina'

    SELECT CASE (QModel%option)
    CASE (1,3)
      IF (QModel%PubliUnit) THEN
        write(out_unit,*) 'PubliUnit=.TRUE.,  Q:[Rad,Bohr,Bohr], Energy: [Hartree]'
      ELSE
        write(out_unit,*) 'PubliUnit=.FALSE., Q:[Rad,Bohr,Bohr], Energy: [Hartree]'
      END IF
    CASE (2)
      IF (QModel%PubliUnit) THEN
        write(out_unit,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad], Energy: [Hartree]'
      ELSE
        write(out_unit,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad]:, Energy: [Hartree]'
      END IF
    END SELECT
    flush(out_unit)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_ClH2p_Botschwina
!> @brief Subroutine wich prints the QML_ClH2p_Botschwina parameters.
!!
!! @param QModel            CLASS(QML_ClH2p_Botschwina_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_QML_ClH2p_Botschwina(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_ClH2p_Botschwina_t), intent(in) :: QModel
    integer,                       intent(in) :: nio

    write(nio,*) 'ClH2p_Botschwina current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '    R2  \ a                            '
    write(nio,*) '         Cl------------H               '
    write(nio,*) '              R1                       '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Coordinates (option 1,2):          '
    write(nio,*) '  Q(1) = R1         (Bohr)             '
    write(nio,*) '  Q(2) = R2         (Bohr)             '
    write(nio,*) '  Q(3) = a          (Radian)           '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) ' Options (default: 6)                  '
    write(nio,*) ' 1: CEPA-1                             '
    write(nio,*) ' 2: CEPA-1-corrected                   '
    write(nio,*) '                                       '
    write(nio,*) 'Ref: Peter Botschwina, ...             '
    write(nio,*) '...  J. Chem. Soc., Faraday Trans. 2, 1988, 84(9), 1263-1276'
    write(nio,*) 'DOI: 10.1039/F29888401263'
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)
    write(nio,*) '---------------------------------------'

    SELECT CASE (QModel%option)

    CASE (1)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: CEPA-1                        '
      write(nio,*) '                                       '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = R1   = 2.4605557069 (Bohr)    '
      write(nio,*) '  Q(2) = R2   = 2.4605557069 (Bohr)    '
      write(nio,*) '  Q(3) = a    = 1.6448681002 (Rad)     '
      write(nio,*) '                                       '
      write(nio,*) '  V = 0.0 Hartree                      '
      write(nio,*) ' grad(:) =[ 0.0, 0.0, 0.0]             '
      write(nio,*) ' hess    =[0.282412000,0.001432000,0.005146000'
      write(nio,*) '           0.001432000,0.282412000,0.005146000'
      write(nio,*) '           0.005146000,0.005146000,0.172946000]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '

    CASE (2)
      write(nio,*) '                                       '
      write(nio,*) '                                       '
      write(nio,*) '  Level: CEPA-1 corrected              '
      write(nio,*) '                                       '
      write(nio,*) '  Minimum:                             '
      write(nio,*) '  Q(1) = R1   = 2.4623131523 (Bohr)    '
      write(nio,*) '  Q(2) = R2   = 2.4623131523 (Bohr)    '
      write(nio,*) '  Q(3) = a    = 1.6448681002 (Rad)     '
      write(nio,*) '                                       '
      write(nio,*) '  V = 0.0 Hartree                      '
      write(nio,*) ' grad(:) =[ 0.0, 0.0, 0.0]             '
      write(nio,*) ' hess    =[0.281000000,0.001475000,0.005122000'
      write(nio,*) '           0.001475000,0.281000000,0.005122000'
      write(nio,*) '           0.005122000,0.005146000,0.167018000]'
      write(nio,*) '                                       '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) ' with Tnum at the minimun.             '
      write(nio,*) '                                       '
    CASE Default
        write(out_unit,*) ' ERROR in write_QModel '
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1,2'
        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end ClH2p_Botschwina current parameters'

  END SUBROUTINE Write_QML_ClH2p_Botschwina

  SUBROUTINE get_Q0_QML_ClH2p_Botschwina(Q0,QModel,option)
    IMPLICIT NONE

    real (kind=Rkind),             intent(inout) :: Q0(:)
    TYPE (QML_ClH2p_Botschwina_t), intent(in)    :: QModel
    integer,                       intent(in)    :: option

    IF (size(Q0) /= 3) THEN
      write(out_unit,*) ' ERROR in get_Q0_QML_ClH2p_Botschwina '
      write(out_unit,*) ' The size of Q0 is not ndim=3: '
      write(out_unit,*) ' size(Q0)',size(Q0)
      STOP 'ERROR in get_Q0_QML_ClH2p_Botschwina: wrong Q0 size'
    END IF

    Q0(:) = QModel%Qref

  END SUBROUTINE get_Q0_QML_ClH2p_Botschwina
!> @brief Subroutine wich calculates the ClH2p_Botschwina potential (unpublished model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_ClH2p_Botschwina_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_ClH2p_Botschwina(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_ClH2p_Botschwina_t), intent(in)    :: QModel
    TYPE (dnS_t),                  intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                  intent(in)    :: dnQ(:) !
    integer,                       intent(in)    :: nderiv

    CALL EvalPot1_QML_ClH2p_Botschwina(QModel,Mat_OF_PotDia,dnQ,nderiv)

  END SUBROUTINE EvalPot_QML_ClH2p_Botschwina

  SUBROUTINE EvalPot1_QML_ClH2p_Botschwina(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_ClH2p_Botschwina_t),  intent(in)    :: QModel
    TYPE (dnS_t),                   intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),                   intent(in)    :: dnQ(:) !
    integer,                        intent(in)    :: nderiv


    TYPE (dnS_t), allocatable :: DQ(:,:)
    TYPE (dnS_t)              :: Vtemp
    integer                   :: i,j,max_exp

    !write(out_unit,*) ' sub EvalPot1_QML_ClH2p_Botschwina' ; flush(6)
    max_exp = 10
    !write(out_unit,*) ' max_exp',max_exp ; flush(6)

    allocate(DQ(QModel%ndim,max_exp))
    DQ(:,1) = dnQ(:) - QModel%Qref(:)
    DO j=2,size(DQ,dim=2)
      DQ(:,j) = DQ(:,j-1) * DQ(:,1)
    END DO
    !write(out_unit,*) ' DQ done' ; flush(6)


    Mat_OF_PotDia(1,1) = ZERO
    Vtemp = Mat_OF_PotDia(1,1) ! to have a correct initialization
    !write(out_unit,*) ' Vtemp init done' ; flush(6)



      DO j=1,QModel%nb_funcModel
        Vtemp = QModel%F(j)
        DO i=1,QModel%ndim
          !Vtemp = Vtemp * DQ(i,1)**QModel%tab_func(i,j)
          IF (QModel%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,QModel%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO
      !CALL Write_dnS(Mat_OF_PotDia(1,1),nio=out_unit)


!-----------------------------------------------------------------------!

   CALL dealloc_dnS(Vtemp)
   CALL dealloc_dnS(DQ)

   !write(out_unit,*) ' end EvalPot1_QML_ClH2p_Botschwina' ; flush(6)

 END SUBROUTINE EvalPot1_QML_ClH2p_Botschwina
END MODULE QML_ClH2p_Botschwina_m
