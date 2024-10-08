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

!> @brief Module which makes the initialization, calculation of the HNNHp potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_HNNHp_m
  USE QDUtil_NumParameters_m
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
    PROCEDURE :: EvalPot_QModel   => EvalPot_QML_HNNHp
    PROCEDURE :: Write_QModel     => Write_QML_HNNHp
    PROCEDURE :: Cart_TO_Q_QModel => Cart_TO_Q_QML_HNNHp

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
    USE QDUtil_m,         ONLY : Identity_Mat, file_open2
    USE QMLLib_UtilLib_m, ONLY : make_QMLInternalFileName

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
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    QModel%QML_Empty_t = QModel_in

    QModel%nsurf    = 1
    QModel%ndimQ    = 6
    QModel%ndimCart = 12
    QModel%pot_name = 'hnnhp'

    IF (QModel%option < 1 .OR. QModel%option > 1) QModel%option = 1

    SELECT CASE (QModel%option)
    CASE (1)

      FileName = make_QMLInternalFileName('InternalData/HNNHp/n2h2pcc5.txt')
      CALL file_open2(name_file=FileName,iunit=nio_fit,old=.TRUE.)
      read(nio_fit,*) QModel%nb_funcModel
      !write(out_unit,*) 'nb_funcModel',QModel%nb_funcModel

      allocate(QModel%Qref(QModel%ndimQ))
      allocate(QModel%F(QModel%nb_funcModel))
      allocate(QModel%tab_func(QModel%ndimQ,QModel%nb_funcModel))

       k = 0
       DO
        nb_columns = min(6,QModel%nb_funcModel-k)
        !write(out_unit,*) k+1,k+nb_columns,nb_columns
        IF (nb_columns == 0) EXIT
         read(nio_fit,11) (QModel%tab_func(1:QModel%ndimQ,j),QModel%F(j),j=k+1,k+nb_columns)
 11      format(6i1,f15.8,5(2x,6i1,f15.8))
         k = k + nb_columns
       END DO
       read(nio_fit,*) QModel%Qref(:)
       !write(out_unit,*) 'Qref',QModel%Qref

       close(nio_fit)
       deallocate(FileName)


    CASE Default

        write(out_unit,*) ' ERROR in ',name_sub
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1'

        STOP

    END SELECT

    IF (QModel%Cart_TO_Q) THEN
      QModel%ndim       = QModel%ndimCart
    ELSE
      QModel%ndim       = QModel%ndimQ
    END IF

    IF (debug) write(out_unit,*) 'init Q0 of HNNHp'
    QModel%Q0 = QModel%Qref([2,1,4,3,5,6])

    IF (debug) write(out_unit,*) 'init d0GGdef of HNNHp'
    QModel%d0GGdef = Identity_Mat(QModel%ndim)


    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
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
        write(out_unit,*) ' ERROR in Write_QML_HNNHp'
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1'

        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end HNNHp current parameters'

  END SUBROUTINE Write_QML_HNNHp

!> @brief Subroutine wich calculates the HNNHp potential with derivatives up to the 2d order.
!!
!! @param QModel             TYPE(QML_HNNHp_t):    derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):         Potential with derivatives,.
!! @param dnQ(:)             TYPE (dnS_t)          value for which the potential is calculated
!! @param nderiv             integer:              it enables to secify the derivative order:
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
        write(out_unit,*) ' ERROR in eval_HNNHpPot '
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1'

        STOP
    END SELECT

  END SUBROUTINE EvalPot_QML_HNNHp

!> @brief Subroutine wich calculates the HNNHp potential (Not published model) with derivatives.
!!
!! @param PotVal             TYPE (dnMat_t):      Potential with derivatives,.
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_HNNHp_t): derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to secify the derivative order:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot1_QML_HNNHp(Mat_OF_PotDia,dnQ,QModel)
    USE QDUtil_m
    USE ADdnSVM_m
    IMPLICIT NONE

    TYPE (dnS_t),      intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),      intent(in)    :: dnQ(:)
    TYPE(QML_HNNHp_t), intent(in)    :: QModel


    TYPE (dnS_t)        :: DQ(6,6)
    TYPE (dnS_t)        :: Vtemp
    integer             :: i,j
    !logical, parameter  :: debug = .TRUE.
    logical, parameter  :: debug = .FALSE.

    IF (debug) THEN
      write(out_unit,*) ' sub EvalPot1_QML_HNNHp'
      write(out_unit,*) ' dnQ(:)',get_d0(dnQ)

      !DO i=1,size(dnQ)
      !  CALL Write_dnS(dnQ(i),info='dnQ' // TO_string(i),nderiv=0)
      !END DO
      flush(out_unit)
    END IF
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
        DO i=1,size(DQ,dim=1)
          IF (QModel%tab_func(i,j) == 0) CYCLE
          Vtemp = Vtemp * DQ(i,QModel%tab_func(i,j))
        END DO
        Mat_OF_PotDia(1,1) = Mat_OF_PotDia(1,1) + Vtemp
      END DO


!-----------------------------------------------------------------------!

   CALL dealloc_dnS(Vtemp)
   CALL dealloc_dnS(DQ)

   IF (debug) THEN
    write(out_unit,*) ' end EvalPot1_QML_HNNHp'
    flush(out_unit)
   END IF

  END SUBROUTINE EvalPot1_QML_HNNHp

  ! here we suppose that the atom ordering: N1-N2-H1-H2
  ! the bounds are N1-N2, N1-H1, n2-H2
  SUBROUTINE Cart_TO_Q_QML_HNNHp(QModel,dnX,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_HNNHp_t),  intent(in)    :: QModel
    TYPE (dnS_t),        intent(in)    :: dnX(:,:)
    TYPE (dnS_t),        intent(inout) :: dnQ(:)
    integer,             intent(in)    :: nderiv


    ! local vector
    TYPE (dnS_t)      :: V1(3),V2(3),v3(3),v4(3),V5(3)
    TYPE (dnS_t)      :: ca,sa
    real (kind=Rkind) :: phi



    ! atoms : N1 N2 H3 H4

    !dihedral first: zmatrix order : nc1=4 (H4),nc2=2 (N2),nc3=1 (N1),nc4=3 (H3)
    !CALL calc_vector2(v1,norm1,nc2,nc1,dnx%d0,ncart0) ! N2->H4
    V1(:) = dnX(:,4)-dnX(:,2)  ! N2->H4
    dnQ(4) = sqrt(dot_product(V1,V1))

    !CALL calc_vector2(v2,norm2,nc2,nc3,dnx%d0,ncart0) ! N2->N1
    V2(:) = dnX(:,1)-dnX(:,2)  ! N2->N1
    dnQ(1) = sqrt(dot_product(V2,V2))
    dnQ(5) = acos(dot_product(V1,V2)/(dnQ(1)*dnQ(4))) ! aNNH1

    !CALL calc_vector2(v3,norm3,nc3,nc4,dnx%d0,ncart0) ! N1->H3
    V3(:) = dnX(:,3)-dnX(:,1)
    dnQ(2) = sqrt(dot_product(V3,V3)) ! N1->H3
    dnQ(3) = acos(-dot_product(V3,V2)/(dnQ(1)*dnQ(2))) ! aNNH1

    !CALL calc_cross_product(v1,norm1,v2,norm2,v4,norm4)
    !CALL calc_cross_product(v3,norm3,v2,norm2,v5,norm5)
    V4(:) = cross_product(V1,V2)
    V5(:) = cross_product(V3,V2)


    !CALL calc_angle_d(angle_d,v4,norm4,v5,norm5,v2,norm2)
    !SUBROUTINE calc_angle_d(angle_d,v1,norm1,v2,norm2,v3,norm3)
    ca = dot_product(v4,v5)
    sa = (v4(1)*(v5(2)*v2(3)-v5(3)*v2(2))                             &
         -v4(2)*(v5(1)*v2(3)-v5(3)*v2(1))                             &
         +v4(3)*(v5(1)*v2(2)-v5(2)*v2(1)))                            &
         /dnQ(1)

    !write(6,*) 'cos, sin: ',get_d0(ca),get_d0(sa)
    dnQ(6) = atan2(sa,ca)

    ! Transformation to have phi (dnQ(6)) between [0,2pi]
    phi = mod(get_d0(dnQ(6)),TWO*pi)
    IF (phi < ZERO) phi = phi + TWO*pi
    CALL set_dnS(dnQ(6),phi,ider=[0])

    !write(6,*) 'dnQ: ',get_d0(dnQ)

  END SUBROUTINE Cart_TO_Q_QML_HNNHp

  END MODULE QML_HNNHp_m
