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

!> @brief Module which makes the initialization, calculation of the ClH2p potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_ClH2p_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the ClH2p parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE, EXTENDS (QML_Empty_t) ::  QML_ClH2p_t

     PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

     integer                       :: nb_funcModel
     real(kind=Rkind), allocatable :: F(:)
     integer, allocatable          :: tab_func(:,:)

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_ClH2p
    PROCEDURE :: Write_QModel    => Write_QML_ClH2p
    PROCEDURE :: Write0_QModel   => Write_QML_ClH2p
  END TYPE QML_ClH2p_t

  PUBLIC :: QML_ClH2p_t,Init_QML_ClH2p

  CONTAINS
!> @brief Subroutine which makes the initialization of the ClH2p parameters.
!!
!! @param ClH2pPot           TYPE(QML_ClH2p_t):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_ClH2p(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE
    TYPE (QML_ClH2p_t)                           :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !local variable
    integer                         :: err_read,nio_fit,i,j,k,idum
    character (len=50)              :: name_Qref
    character (len=:), allocatable  :: FileName

    !-----------------------------------------------------------------
    ! for the namelist
    integer                  :: MinCoupling,MaxCoupling
    real (kind=Rkind)        :: MinNorm,MaxNorm,max_valB
    integer                  :: ndim,nb_B,nb_WB,MR_order
    logical                  :: svd
    real (kind=Rkind)        :: epsi,epsi_inter
    integer                  :: ind_val,nb_val,max_nb
    integer                  :: Col_FOR_WeightOFFit
    real (kind=Rkind)        :: Scal_FOR_WeightOFFit
    integer                  :: nb_Fit

    namelist /nDFitW/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,        &
                      ind_val,nb_val,svd,                             &
                      epsi,epsi_inter,max_nb,ndim,MR_order,           &
                      Col_FOR_WeightOFFit,Scal_FOR_WeightOFFit,       &
                      nb_B,nb_WB,                                     &
                      nb_Fit

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_ClH2p'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)



    QModel%nsurf    = 1
    QModel%pot_name = 'ClH2p'


    IF (QModel%option < 1 .OR. QModel%option > 2) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1)

      FileName = make_FileName('InternalData/ClH2p/ClH2p_B3LYP_cc-pVTZ.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)

    CASE Default

      write(out_unitp,*) ' ERROR in Init_QML_ClH2p '
      write(out_unitp,*) ' This option is not possible. option: ',QModel%option
      write(out_unitp,*) ' Its value MUST be 1'
      STOP

    END SELECT

    read(nio_fit,nDFitW)
    IF (debug) write(out_unitp,nDFitW)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "nDFitW" is probably absent'
      write(out_unitp,*) ' from the file: ',trim(adjustl(FileName))
      write(out_unitp,*) ' ERROR in ',name_sub
      STOP
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Some parameter name of the namelist "nDFitW" are probaly wrong'
      write(out_unitp,*) ' in the file: ',trim(adjustl(FileName))
      write(out_unitp,*) ' It should never append !!'
      write(out_unitp,*) ' Check the fortran'
      write(out_unitp,*) ' ERROR in ',name_sub
      STOP
    END IF

    QModel%ndim         = ndim
    QModel%nb_funcModel = nb_B

    allocate(QModel%Qref(QModel%ndim))
    read(nio_fit,*) name_Qref,QModel%Qref(:)
    read(nio_fit,*)
    read(nio_fit,*)
    read(nio_fit,*)

    allocate(QModel%F(QModel%nb_funcModel))
    allocate(QModel%tab_func(QModel%ndim,QModel%nb_funcModel))


    DO i=1,QModel%nb_funcModel
      read(nio_fit,*) idum,QModel%tab_func(:,i),QModel%F(i)
      QModel%tab_func(3,i) = 2*QModel%tab_func(3,i) ! for even momomials
      IF (debug) write(out_unitp,*) 'i,tab_func,F',i,QModel%tab_func(:,i),QModel%F(i)
    END DO
    close(nio_fit)
    deallocate(FileName)


    IF (debug) write(out_unitp,*) 'init Q0 of ClH2p'
    QModel%Q0 = QModel%Qref
    IF (debug) write(out_unitp,*) 'QModel%Q0',QModel%Q0

    IF (debug) write(out_unitp,*) 'init d0GGdef of ClH2p'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)

    IF (QModel%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Rad,Bohr,Bohr], Energy: [Hartree]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Rad,Bohr,Bohr], Energy: [Hartree]'
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_ClH2p
!> @brief Subroutine wich prints the QML_ClH2p parameters.
!!
!! @param QModel            CLASS(QML_ClH2p_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_QML_ClH2p(QModel,nio)

    CLASS(QML_ClH2p_t),   intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'ClH2p current parameters'
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
    write(nio,*) '  Coordinates :                        '
    write(nio,*) '  Q(1) = a          (Radian)           '
    write(nio,*) '  Q(2) = 1/2(R1+R2) (Bohr)             '
    write(nio,*) '  Q(3) = 1/2(R1-R2) (Bohr)             '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) 'Ref: unpublished'
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
    write(nio,*) '  Level: B3LYP/cc-pVTZ                 '
    write(nio,*) '                                       '
    write(nio,*) '  Minimum:                             '
    write(nio,*) '  Q(1) = a    = 1.6525859462 (Rad)     '
    write(nio,*) '  Q(2) = r+   = 2.4849898659 (Bohr)    '
    write(nio,*) '  Q(3) = r-   = 0.           (Bohr)    '
    write(nio,*) '                                       '
    write(nio,*) '  V = 0.0004486447 Hartree             '
    write(nio,*) ' grad(:) =[ 0.0000558107,-0.0004452755,'
    write(nio,*) '            0.0000000000]              '
    write(nio,*) ' hess    =[0.158044529,0.009093148,0.000000000'
    write(nio,*) '           0.009093148,0.520862288,0.000000000'
    write(nio,*) '           0.000000000,0.000000000,0.511935490]'
    write(nio,*) '                                       '

    CASE Default
        write(out_unitp,*) ' ERROR in write_QModel '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1 or 2'

        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end ClH2p current parameters'

  END SUBROUTINE Write_QML_ClH2p

  SUBROUTINE get_Q0_QML_ClH2p(Q0,QModel,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (QML_ClH2p_t),          intent(in)    :: QModel
    integer,                     intent(in)    :: option

    IF (size(Q0) /= 3) THEN
      write(out_unitp,*) ' ERROR in get_Q0_QML_ClH2p '
      write(out_unitp,*) ' The size of Q0 is not ndim=3: '
      write(out_unitp,*) ' size(Q0)',size(Q0)
      STOP
    END IF

    Q0(:) = QModel%Qref(:)

  END SUBROUTINE get_Q0_QML_ClH2p
!> @brief Subroutine wich calculates the ClH2p potential (unpublished model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_ClH2p_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_ClH2p(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE QMLdnSVM_dnS_m

    CLASS(QML_ClH2p_t),   intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv


    TYPE (dnS_t), allocatable :: DQ(:,:)
    TYPE (dnS_t)              :: Vtemp
    integer                   :: i,j,max_exp

    !write(6,*) ' sub EvalPot_QML_ClH2p' ; flush(6)
    max_exp = maxval(QModel%tab_func)
    !write(6,*) ' max_exp',max_exp ; flush(6)

    allocate(DQ(QModel%ndim,max_exp))
    DQ(:,1) = dnQ(:) - QModel%Qref(:)
    DO j=2,size(DQ,dim=2)
      DQ(:,j) = DQ(:,j-1) * DQ(:,1)
    END DO
    !write(6,*) ' DQ done' ; flush(6)


    Mat_OF_PotDia(1,1) = ZERO
    Vtemp = Mat_OF_PotDia(1,1) ! to have a correct initialization
    !write(6,*) ' Vtemp init done' ; flush(6)



      DO j=1,QModel%nb_funcModel
        Vtemp = QModel%F(j)
        DO i=1,QModel%ndim
          !Vtemp = Vtemp * DQ(i,1)**QModel%tab_func(i,j)
          IF (QModel%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,QModel%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO
      !CALL QML_Write_dnS(Mat_OF_PotDia(1,1),nio=out_unitp)


!-----------------------------------------------------------------------!

   CALL QML_dealloc_dnS(Vtemp)
   CALL QML_dealloc_dnS(DQ)

   !write(6,*) ' end EvalPot1_QML_ClH2p' ; flush(6)

  END SUBROUTINE EvalPot_QML_ClH2p

END MODULE QML_ClH2p_m
