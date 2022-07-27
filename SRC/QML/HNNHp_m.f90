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

!> @brief Module which makes the initialization, calculation of the HNNHp potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_HNNHp_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the HNNHp parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_HNNHp_t
   PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

     integer                       :: nb_funcModel
     real(kind=Rkind), allocatable :: F(:)
     integer,          allocatable :: tab_func(:,:)

   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_HNNHp
    PROCEDURE :: Write_QModel    => Write_QML_HNNHp
    PROCEDURE :: Write0_QModel   => Write0_QML_HNNHp
  END TYPE QML_HNNHp_t

  PUBLIC :: QML_HNNHp_t,Init_QML_HNNHp

  CONTAINS
!> @brief Function which makes the initialization of the HNNHp parameters.
!!
!! @param QModel             TYPE(QML_HNNHp_t): result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_HNNHp(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_HNNHp_t)                           :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    integer :: nio_fit,nb_columns,i,j,k
    character (len=:), allocatable  :: FileName

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_HNNHp'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 1
    QModel%ndim     = 6
    QModel%pot_name = 'hnnhp'

    IF (QModel%option < 1 .OR. QModel%option > 1) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1)

      FileName = make_FileName('InternalData/HNNHp/n2h2pcc5.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      read(nio_fit,*) QModel%nb_funcModel
      !write(out_unitp,*) 'nb_funcModel',QModel%nb_funcModel

      allocate(QModel%Qref(QModel%ndim))
      allocate(QModel%F(QModel%nb_funcModel))
      allocate(QModel%tab_func(QModel%ndim,QModel%nb_funcModel))

       k = 0
       DO
        nb_columns = min(6,QModel%nb_funcModel-k)
        !write(out_unitp,*) k+1,k+nb_columns,nb_columns
        IF (nb_columns == 0) EXIT
         read(nio_fit,11) (QModel%tab_func(1:QModel%ndim,j),QModel%F(j),j=k+1,k+nb_columns)
 11      format(6i1,f15.8,5(2x,6i1,f15.8))
         k = k + nb_columns
       END DO
       read(nio_fit,*) QModel%Qref(:)
       !write(out_unitp,*) 'Qref',QModel%Qref

       close(nio_fit)
       deallocate(FileName)


    CASE Default

        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1'

        STOP

    END SELECT


    IF (debug) write(out_unitp,*) 'init Q0 of HNNHp'
    QModel%Q0 = QModel%Qref([2,1,4,3,5,6])

    IF (debug) write(out_unitp,*) 'init d0GGdef of HNNHp'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_HNNHp
!> @brief Subroutine wich prints the QML_HNNHp parameters.
!!
!! @param QModel            CLASS(QML_HNNHp_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_HNNHp(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HNNHp_t),   intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'trans-HNNH+ current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '   Q(2) \ Q(3)                         '
    write(nio,*) '         N------------N                '
    write(nio,*) '            Q(1)   Q(5) \  Q(4)        '
    write(nio,*) '                          \            '
    write(nio,*) '                            H    +Q(6) '
    write(nio,*) '                                       '
    write(nio,*) '  Q(1) = dNN (Bohr)                    '
    write(nio,*) '  Q(2) = dNH (Bohr)                    '
    write(nio,*) '  Q(3) = aHNN (Radian)                 '
    write(nio,*) '  Q(4) = dNH (Bohr)                    '
    write(nio,*) '  Q(5) = aHNN (Radian)                 '
    write(nio,*) '  Q(6) = Torsion (HNNH) (Radian)       '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) 'Ref:'
    write(nio,*) '  D. Lauvergnat and M. Hochlaf, J. Chem. Phys. 130, 224312 (2009).'
    write(nio,*) '      https://doi.org/10.1063/1.3154141'
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
    write(nio,*) '  Level: RCCSD(T)/aug-cc-pV5Z*         '
    write(nio,*) '    spdfg for N and spdf for H         '
    write(nio,*) '                                       '
    write(nio,*) '  Reference Trans-geometry:            '
    write(nio,*) '  Q(1) = rNN  = 2.2           (Bohr)   '
    write(nio,*) '  Q(2) = rH1  = 2.0           (Bohr)   '
    write(nio,*) '  Q(3) = aH1  = 2.1816615395  (rad)    '
    write(nio,*) '  Q(4) = rH2  = 2.0           (Bohr)   '
    write(nio,*) '  Q(5) = aH3  = 2.1816615395  (rad)    '
    write(nio,*) '  Q(6) = phi  = 3.1415926169  (rad)    '
    write(nio,*) '                                       '
    write(nio,*) '  V = -110.16613531 Hartree            '
    write(nio,*) ' grad(:) =[-0.0002097800, 0.0133780500,'
    write(nio,*) '            0.0054179000, 0.0133780500,'
    write(nio,*) '            0.0054179000, 0.0000000000]'
    write(nio,*) '                                       '

    CONTINUE

    CASE Default
        write(out_unitp,*) ' ERROR in Write_QML_HNNHp'
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1'

        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end HNNHp current parameters'

  END SUBROUTINE Write_QML_HNNHp
  SUBROUTINE Write0_QML_HNNHp(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HNNHp_t),   intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'HNNHp default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end HNNHp default parameters'


  END SUBROUTINE Write0_QML_HNNHp

!> @brief Subroutine wich calculates the HNNHp potential with derivatives up to the 2d order.
!!
!! @param QModel             TYPE(QML_HNNHp_t):    derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_HNNHp(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HNNHp_t), intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    SELECT CASE (QModel%option)

    CASE (1)
      CALL EvalPot1_QML_HNNHp(Mat_OF_PotDia,dnQ,QModel)

    CASE Default
        write(out_unitp,*) ' ERROR in eval_HNNHpPot '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1'

        STOP
    END SELECT

  END SUBROUTINE EvalPot_QML_HNNHp

!> @brief Subroutine wich calculates the HNNHp potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_HNNHp_t): derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot1_QML_HNNHp(Mat_OF_PotDia,dnQ,QModel)
  USE ADdnSVM_m
  IMPLICIT NONE

    TYPE (dnS_t),      intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),      intent(in)    :: dnQ(:)
    TYPE(QML_HNNHp_t), intent(in)    :: QModel


    TYPE (dnS_t)        :: DQ(6,6)
    TYPE (dnS_t)        :: Vtemp
    integer             :: i,j

    !write(out_unitp,*) ' sub EvalPot1_QML_HNNHp' ; flush(6)

      ! Warning, the coordinate ordering in the potential data (from the file) is different from the z-matrix one.
      DQ(:,1) = dnQ([2,1,4,3,5,6]) - QModel%Qref(:)
      DO j=2,size(DQ,dim=2)
        DQ(1:5,j) = DQ(1:5,j-1) * DQ(1:5,1)
      END DO
      DO i=size(DQ,dim=2),1,-1
        DQ(6,i) = cos( DQ(6,1) * real(i,kind=Rkind) )
      END DO


      Mat_OF_PotDia(1,1) = ZERO
      Vtemp = Mat_OF_PotDia(1,1) ! to have a correct initialization

      DO j=1,QModel%nb_funcModel
        Vtemp = QModel%F(j)
        DO i=1,QModel%ndim
          IF (QModel%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,QModel%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO


!-----------------------------------------------------------------------!

   CALL dealloc_dnS(Vtemp)
   CALL dealloc_dnS(DQ)

   !write(out_unitp,*) ' end EvalPot1_QML_HNNHp' ; flush(6)

  END SUBROUTINE EvalPot1_QML_HNNHp

  ! here we suppose that the atom ordering: N1-N2-H1-H2
  ! the bounds are N1-N2, N1-H1, n2-H2
  SUBROUTINE Cart_TO_Q_QML_HNNHp(dnX,dnQ,QModel,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    TYPE (dnS_t),        intent(in)    :: dnX(:,:)
    TYPE (dnS_t),        intent(inout) :: dnQ(:)
    TYPE(QML_HNNHp_t), intent(in)    :: QModel
    integer,             intent(in)    :: nderiv


    ! local vector
    TYPE (dnS_t)    :: VecNN(3),VecN1H1(3),VecN2H2(3)


    VecNN(:)   = dnX(:,1)-dnX(:,2)
    VecN1H1(:) = dnX(:,1)-dnX(:,3)
    VecN2H2(:) = dnX(:,2)-dnX(:,4)

    dnQ(1) = sqrt(dot_product(VecNN,VecNN))
    dnQ(2) = sqrt(dot_product(VecN1H1,VecN1H1))
    dnQ(4) = sqrt(dot_product(VecN2H2,VecN2H2))

    dnQ(3) = acos(dot_product(VecNN,VecN1H1)/(dnQ(1)*dnQ(2)))
    dnQ(5) = acos(dot_product(VecNN,VecN2H2)/(dnQ(1)*dnQ(4)))



  END SUBROUTINE Cart_TO_Q_QML_HNNHp

END MODULE QML_HNNHp_m
