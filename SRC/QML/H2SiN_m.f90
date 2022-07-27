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
!    Copyright 2016  David LAUVERGNAT, Félix MOUHAT and Liang LIANG
!
!===========================================================================
!===========================================================================

!> @brief Module which makes the initialization, calculation of the H2SiN potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_H2SiN_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the H2SiN parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_H2SiN_t
   PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

     integer                       :: nb_funcModel
     real(kind=Rkind), allocatable :: F(:)
     integer,          allocatable :: tab_func(:,:)

   CONTAINS
    PROCEDURE :: EvalPot_QModel => EvalPot_QML_H2SiN
    PROCEDURE :: Write_QModel    => Write_QML_H2SiN
    PROCEDURE :: Write0_QModel   => Write0_QML_H2SiN
  END TYPE QML_H2SiN_t

  PUBLIC :: QML_H2SiN_t,Init_QML_H2SiN

  CONTAINS
!> @brief Function which makes the initialization of the H2SiN parameters.
!!
!! @param QModel             TYPE(QML_H2SiN_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_H2SiN(QModel_in,read_param,nio_param_file) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_H2SiN_t)                         :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param


    integer :: nio_fit,nb_columns,j,k
    character (len=:), allocatable  :: FileName

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H2SiN'
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
    QModel%pot_name = 'h2sin'

    IF (QModel%option < 1 .OR. QModel%option > 3) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1)

      FileName = make_FileName('InternalData/H2SiN/h2sinf12a.pot')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      read(nio_fit,*) QModel%nb_funcModel
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
       QModel%Qref(6) = pi ! Qref(6) must be changed to be compatible with a z-matrix.
       close(nio_fit)
       deallocate(FileName)

    CASE (2)

      FileName = make_FileName('InternalData/H2SiN/h2sinf12b.pot')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      read(nio_fit,*) QModel%nb_funcModel
      allocate(QModel%Qref(QModel%ndim))
      allocate(QModel%F(QModel%nb_funcModel))
      allocate(QModel%tab_func(QModel%ndim,QModel%nb_funcModel))

       k = 0
       DO
        nb_columns = min(6,QModel%nb_funcModel-k)
        !write(out_unitp,*) k+1,k+nb_columns,nb_columns
        IF (nb_columns == 0) EXIT
         read(nio_fit,11) (QModel%tab_func(1:QModel%ndim,j),QModel%F(j),j=k+1,k+nb_columns)
         k = k + nb_columns
       END DO
       read(nio_fit,*) QModel%Qref(:)
       QModel%Qref(6) = pi ! Qref(6) must be changed to be compatible with a z-matrix.
       close(nio_fit)
       deallocate(FileName)

    CASE (3)

      FileName = make_FileName('InternalData/H2SiN/h2sincc.pot')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      read(nio_fit,*) QModel%nb_funcModel
      allocate(QModel%Qref(QModel%ndim))
      allocate(QModel%F(QModel%nb_funcModel))
      allocate(QModel%tab_func(QModel%ndim,QModel%nb_funcModel))

       k = 0
       DO
        nb_columns = min(6,QModel%nb_funcModel-k)
        !write(out_unitp,*) k+1,k+nb_columns,nb_columns
        IF (nb_columns == 0) EXIT
         read(nio_fit,11) (QModel%tab_func(1:QModel%ndim,j),QModel%F(j),j=k+1,k+nb_columns)
         k = k + nb_columns
       END DO
       read(nio_fit,*) QModel%Qref(:)
       QModel%Qref(6) = pi ! Qref(6) must be changed to be compatible with a z-matrix.
       close(nio_fit)
       deallocate(FileName)

    CASE Default

      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' This option is not possible. option: ',QModel%option
      write(out_unitp,*) ' Its value MUST be 1,2 or 3'
      STOP
    END SELECT

    IF (debug) write(out_unitp,*) 'init Q0 of H2SiN'
    QModel%Q0 = QModel%Qref([3,1,4,2,5,6])

    IF (debug) write(out_unitp,*) 'init d0GGdef of H2SiN'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_H2SiN
!> @brief Subroutine wich prints the QML_H2SiN parameters.
!!
!! @param QModel            CLASS(QML_H2SiN_t):   derived type in which the parameters are set-up.
!! @param nio               integer:              file unit to print the parameters.
  SUBROUTINE Write_QML_H2SiN(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_H2SiN_t),   intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'H2SiN current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '   Q(2) \ Q(3)                         '
    write(nio,*) '         Si------------N               '
    write(nio,*) '   Q(4) / Q(5)    Q(1)                 '
    write(nio,*) '       /                               '
    write(nio,*) '      H   +Q(6)                        '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Minimum:                             '
    write(nio,*) '  Q(1) = rSiN (Bohr)                   '
    write(nio,*) '  Q(2) = rH1  (Bohr)                   '
    write(nio,*) '  Q(3) = aH1  (Radian)                 '
    write(nio,*) '  Q(4) = rH2  (Bohr)                   '
    write(nio,*) '  Q(5) = aH3  (Radian)                 '
    write(nio,*) '  Q(6) = phi  (Radian)                 '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) '                                       '
    write(nio,*) 'Ref:'
    write(nio,*) '  D. Lauvergnat, M.L. Senent, L. Jutier, and M. Hochlaf, JCP 135, 074301 (2011).'
    write(nio,*) '      https://doi.org/10.1063/1.3624563'
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
    write(nio,*) '  Level: RCCSD(T)-F12a/cc-pVTZ-F12     '
    write(nio,*) '                                       '
    write(nio,*) '  Minimum:                             '
    write(nio,*) '  Q(1) = rSiN = 3.1103305087 (Bohr)    '
    write(nio,*) '  Q(2) = rH1  = 2.7870835493 (Bohr)    '
    write(nio,*) '  Q(3) = aH1  = 2.1336938250 (Radian)  '
    write(nio,*) '  Q(4) = rH2  = 2.7870835493 (Bohr)    '
    write(nio,*) '  Q(5) = aH3  = 2.1336938250 (Radian)  '
    write(nio,*) '  Q(6) = phi  =  pi (Radian)           '
    write(nio,*) '                                       '
    write(nio,*) '  V = -344.84785889 Hartree            '
    write(nio,*) ' grad(:) =[-0.0000000100,-0.0000000100,'
    write(nio,*) '            0.0000000000,-0.0000000100,'
    write(nio,*) '            0.0000000000, 0.0000000000]'
    write(nio,*) '                                       '

    CASE (2)
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Level: RCCSD(T)-F12b/cc-pVTZ-F12     '
    write(nio,*) '                                       '
    write(nio,*) '  Minimum:                             '
    write(nio,*) '  Q(1) = rSiN = 3.1096943169 (Bohr)    '
    write(nio,*) '  Q(2) = rH1  = 2.7871268895 (Bohr)    '
    write(nio,*) '  Q(3) = aH1  = 2.1335577958 (Radian)  '
    write(nio,*) '  Q(4) = rH2  = 2.7871268895 (Bohr)    '
    write(nio,*) '  Q(5) = aH3  = 2.1335577958 (Radian)  '
    write(nio,*) '  Q(6) = phi  =  pi (Radian)           '
    write(nio,*) '                                       '
    write(nio,*) '  V = -344.84209978 Hartree            '
    write(nio,*) ' grad(:) =[ 0.0000000000, 0.0000000000,'
    write(nio,*) '            0.0000000000, 0.0000000000,'
    write(nio,*) '            0.0000000000, 0.0000000000]'
    write(nio,*) '                                       '

    CASE (3)
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Level: RCCSD(T)/Aug-cc-pV5Z          '
    write(nio,*) '                                       '
    write(nio,*) '  Reference geometry:                  '
    write(nio,*) '  Q(1) = rSiN = 3.10         (Bohr)    '
    write(nio,*) '  Q(2) = rH1  = 2.80         (Bohr)    '
    write(nio,*) '  Q(3) = aH1  = 2.1816615395 (Radian)  '
    write(nio,*) '  Q(4) = rH2  = 2.80         (Bohr)    '
    write(nio,*) '  Q(5) = aH3  = 2.1816615395 (Radian)  '
    write(nio,*) '  Q(6) = phi  =  pi          (Radian)  '
    write(nio,*) '                                       '
    write(nio,*) '  V = -344.84115204 Hartree            '
    write(nio,*) ' grad(:) =[-0.0036112000, 0.0018618700,'
    write(nio,*) '            0.0112321500, 0.0018618800,'
    write(nio,*) '            0.0112320500,0.0000000000] '
    write(nio,*) '                                       '
    CONTINUE

    CASE Default
        write(out_unitp,*) ' ERROR in Write_QML_H2SiN '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1,2,3'

        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end H2SiN current parameters'

  END SUBROUTINE Write_QML_H2SiN
  SUBROUTINE Write0_QML_H2SiN(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_H2SiN_t),   intent(in) :: QModel
    integer,                intent(in) :: nio

    write(nio,*) 'H2SiN default parameters'
    write(nio,*)
    write(nio,*)
    write(nio,*) 'end H2SiN default parameters'


  END SUBROUTINE Write0_QML_H2SiN

!> @brief Subroutine wich calculates the H2SiN potential with derivatives up to the 2d order.
!!
!! @param QModel             TYPE(QML_H2SiN_t):    derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H2SiN(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_H2SiN_t), intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:)
    integer,              intent(in)    :: nderiv

    SELECT CASE (QModel%option)

    CASE (1,2,3)
      CALL EvalPot1_QML_H2SiN(Mat_OF_PotDia,dnQ,QModel)

    CASE Default
        write(out_unitp,*) ' ERROR in EvalPot_QML_H2SiN '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1, 2 or 3'

        STOP
    END SELECT


  END SUBROUTINE EvalPot_QML_H2SiN

!> @brief Subroutine wich calculates the H2SiN potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):       derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                 value for which the potential is calculated
!! @param QModel             TYPE(QML_H2SiN_t):  derived type in which the parameters are set-up.
!! @param nderiv             integer:              it enables to specify up to which derivatives the potential is calculated:
!!                                                 the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot1_QML_H2SiN(Mat_OF_PotDia,dnQ,QModel)
  USE ADdnSVM_m
  IMPLICIT NONE

    TYPE (dnS_t),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)    :: dnQ(:)
    TYPE(QML_H2SiN_t), intent(in)    :: QModel


    TYPE (dnS_t)        :: DQ(6,4)
    TYPE (dnS_t)        :: Vtemp
    integer             :: i,j

    !write(out_unitp,*) ' sub EvalPot1_QML_H2SiN' ; flush(6)

      ! Warning, the coordinate ordering in the potential data (from the file) is different from the z-matrix one.
      DQ(:,1) = dnQ([2,4,1,3,5,6]) - QModel%Qref(:)
      DO j=2,4
        DQ(:,j) = DQ(:,j-1) * DQ(:,1)
      END DO


      Mat_OF_PotDia(1,1) = ZERO
      Vtemp = Mat_OF_PotDia(1,1) ! to have a correct initialization

      DO j=1,QModel%nb_funcModel
        Vtemp = QModel%F(j)
        DO i=1,QModel%ndim
          !Vtemp = Vtemp * DQ(i,1)**QModel%tab_func(i,j)
          IF (QModel%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,QModel%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO


!-----------------------------------------------------------------------!

   CALL dealloc_dnS(Vtemp)
   CALL dealloc_dnS(DQ)

   !write(out_unitp,*) ' end EvalPot1_QML_H2SiN' ; flush(6)

  END SUBROUTINE EvalPot1_QML_H2SiN


END MODULE QML_H2SiN_m
