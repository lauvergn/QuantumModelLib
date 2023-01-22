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

!> @brief Module which makes the initialization, calculation of the H2NSi potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_H2NSi_m
  USE QDUtil_NumParameters_m, out_unitp => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the H2NSi parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE, EXTENDS (QML_Empty_t) ::  QML_H2NSi_t

     PRIVATE

     real(kind=Rkind), allocatable :: Qref(:)

     integer                       :: nb_funcModel
     real(kind=Rkind), allocatable :: F(:)
     integer, allocatable          :: tab_func(:,:)

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_H2NSi
    PROCEDURE :: Write_QModel    => Write_QML_H2NSi
    PROCEDURE :: Write0_QModel   => Write_QML_H2NSi
  END TYPE QML_H2NSi_t

  PUBLIC :: QML_H2NSi_t,Init_QML_H2NSi

  CONTAINS
!> @brief Subroutine which makes the initialization of the H2NSi parameters.
!!
!! @param H2NSiPot          TYPE(QML_H2NSi_t):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_H2NSi(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m, ONLY : Identity_Mat, make_FileName, file_open2
    IMPLICIT NONE

    TYPE (QML_H2NSi_t)                           :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !local variable
    integer                         :: nio_fit,nb_columns,i,j,k
    character (len=50)              :: name_Qref
    character (len=:), allocatable  :: FileName

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_H2NSi'
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
    QModel%pot_name = 'h2nsi'

    IF (QModel%option < 1 .OR. QModel%option > 2) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1)

      FileName = make_FileName('InternalData/H2NSi/h2nsif12a.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      allocate(QModel%Qref(QModel%ndim))

      DO i=1,QModel%ndim
        read(nio_fit,*) name_Qref,QModel%Qref(i)
      END DO
      QModel%Qref(4:6) = QModel%Qref(4:6) / (180._Rkind/pi)
      QModel%Qref(6)   = pi

      read(nio_fit,*) QModel%nb_funcModel
      allocate(QModel%F(QModel%nb_funcModel))
      allocate(QModel%tab_func(QModel%ndim,QModel%nb_funcModel))


      DO i=1,QModel%nb_funcModel
         read(nio_fit,*) name_Qref,QModel%tab_func(:,i),QModel%F(i)
         QModel%F(i) = QModel%F(i) * (180._Rkind/pi)**sum(QModel%tab_func(4:QModel%ndim,i))
      END DO
      close(nio_fit)
      deallocate(FileName)

    CASE (2)

      FileName = make_FileName('InternalData/H2NSi/h2nsicc.pot')
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


    CASE Default

          write(out_unitp,*) ' ERROR in Init_QML_H2NSi '
          write(out_unitp,*) ' This option is not possible. option: ',QModel%option
          write(out_unitp,*) ' Its value MUST be 1 or 2'

          STOP

      END SELECT


     IF (debug) write(out_unitp,*) 'init Q0 of H2NSi'
     QModel%Q0 = QModel%Qref([3,1,4,2,5,6])

    IF (debug) write(out_unitp,*) 'init d0GGdef of H2NSi'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)


    IF (QModel%PubliUnit) THEN
      write(out_unitp,*) 'PubliUnit=.TRUE.,  Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    ELSE
      write(out_unitp,*) 'PubliUnit=.FALSE., Q:[Bohr,Bohr,Rad,Bohr,Rad,Rad], Energy: [Hartree]'
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_H2NSi
!> @brief Subroutine wich prints the QML_H2NSi parameters.
!!
!! @param QModel            CLASS(QML_H2NSi_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_QML_H2NSi(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_H2NSi_t), intent(in) :: QModel
    integer,              intent(in) :: nio

    write(nio,*) 'H2NSi current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '         Internal coordinates          '
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '      H                                '
    write(nio,*) '       \                               '
    write(nio,*) '   Q(2) \ Q(3)                         '
    write(nio,*) '         N-------------Si              '
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
    write(nio,*) '  Q(1) = rSiN = 3.2262192852 (Bohr)    '
    write(nio,*) '  Q(2) = rH1  = 1.9084802130 (Bohr)    '
    write(nio,*) '  Q(3) = aH1  = 124.2468848  (°)       '
    write(nio,*) '  Q(4) = rH2  = 1.9084802130 (Bohr)    '
    write(nio,*) '  Q(5) = aH3  = 124.2468848  (°)       '
    write(nio,*) '  Q(6) = phi  = 180.         (°)       '
    write(nio,*) '                                       '
    write(nio,*) '  V = -344.92729937 Hartree            '
    write(nio,*) ' grad(:) =[ 0.0000000000, 0.0000000000,'
    write(nio,*) '            0.0000000000, 0.0000000000,'
    write(nio,*) '            0.0000000000, 0.0000000000]'
    write(nio,*) '                                       '
    write(nio,*) '  The angles are convert in radian.    '

    CASE (2)
    write(nio,*) '                                       '
    write(nio,*) '                                       '
    write(nio,*) '  Level: RCCSD(T)/Aug-cc-pV5Z          '
    write(nio,*) '                                       '
    write(nio,*) '  Reference geometry:                  '
    write(nio,*) '  Q(1) = rSiN = 3.20         (Bohr)    '
    write(nio,*) '  Q(2) = rH1  = 1.90         (Bohr)    '
    write(nio,*) '  Q(3) = aH1  = 2.1816615395 (Radian)  '
    write(nio,*) '  Q(4) = rH2  = 1.90         (Bohr)    '
    write(nio,*) '  Q(5) = aH3  = 2.1816615395 (Radian)  '
    write(nio,*) '  Q(6) = phi  = pi           (Radian)  '
    write(nio,*) '                                       '
    write(nio,*) '  V = -344.92196350 Hartree            '
    write(nio,*) ' grad(:) =[-0.0077950700,-0.0035043800,'
    write(nio,*) '            0.0024132300,-0.0035043800,'
    write(nio,*) '            0.0024132300,0.0000000000] '
    write(nio,*) '                                       '
    write(nio,*) ' Unpublished potential                 '

    CASE Default
        write(out_unitp,*) ' ERROR in write_QModel '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1 or 2'

        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end H2NSi current parameters'

  END SUBROUTINE Write_QML_H2NSi

  SUBROUTINE get_Q0_QML_H2NSi(Q0,QModel,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (QML_H2NSi_t),        intent(in)    :: QModel
    integer,                     intent(in)    :: option

    IF (size(Q0) /= 6) THEN
      write(out_unitp,*) ' ERROR in get_Q0_QML_H2NSi '
      write(out_unitp,*) ' The size of Q0 is not ndim=6: '
      write(out_unitp,*) ' size(Q0)',size(Q0)
      STOP
    END IF

    SELECT CASE (option)
    CASE (0) ! ref
      Q0(:) = QModel%Qref([3,1,4,2,5,6])
    CASE Default ! ref
      Q0(:) = QModel%Qref([3,1,4,2,5,6])
    END SELECT

  END SUBROUTINE get_Q0_QML_H2NSi
!> @brief Subroutine wich calculates the H2NSi potential (for the 3 models) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel         TYPE(QML_H2NSi_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_H2NSi(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_H2NSi_t), intent(in)    :: QModel
    TYPE (dnS_t),         intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),         intent(in)    :: dnQ(:) !
    integer,              intent(in)    :: nderiv

    SELECT CASE (QModel%option)

    CASE (1,2)
      CALL EvalPot1_QML_H2NSi(Mat_OF_PotDia,dnQ,QModel)

    CASE Default
        write(out_unitp,*) ' ERROR in EvalPot_QML_H2NSi '
        write(out_unitp,*) ' This option is not possible. option: ',QModel%option
        write(out_unitp,*) ' Its value MUST be 1 or 2'

        STOP
    END SELECT

  END SUBROUTINE EvalPot_QML_H2NSi

!> @brief Subroutine wich calculates the H2NSi potential (Not published model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel         TYPE(QML_H2NSi_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).

  SUBROUTINE EvalPot1_QML_H2NSi(Mat_OF_PotDia,dnQ,QModel)
    !Unpublished model potential
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),        intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),        intent(in)    :: dnQ(:)
    TYPE(QML_H2NSi_t) , intent(in)    :: QModel


    TYPE (dnS_t)        :: DQ(6,6)
    TYPE (dnS_t)        :: Vtemp
    integer            :: i,j

    !write(out_unitp,*) ' sub EvalPot1_QML_H2NSi' ; flush(6)

      ! Warning, the coordinate ordering in the potential data (from the file) is different from the z-matrix one.
      DQ(:,1) = dnQ([2,4,1,3,5,6]) - QModel%Qref(:)
      DO j=2,size(DQ,dim=2)
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

   !write(out_unitp,*) ' end EvalPot1_QML_H2NSi' ; flush(6)

  END SUBROUTINE EvalPot1_QML_H2NSi

END MODULE QML_H2NSi_m
