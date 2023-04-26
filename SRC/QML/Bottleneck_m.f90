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

!> @brief Module which makes the initialization, calculation of the Bottleneck potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 01/02/2023
!!
MODULE QML_Bottleneck_m
  USE QDUtil_NumParameters_m, out_unit => out_unit
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the Bottleneck parameters are set-up.
!!
!! @param option                  integer: it enables to chose between the 1 model(s) (default 1)

  TYPE, EXTENDS (QML_Empty_t) ::  QML_Bottleneck_t

    PRIVATE

    ! parameter of the reference, option=1
    real(kind=Rkind) :: V0 =  0.036_Rkind ! Hartree
    real(kind=Rkind) :: a =  0.4_Rkind ! bohr^-1
    real(kind=Rkind) :: k0 =  0.1909_Rkind ! au
    real(kind=Rkind) :: sig=0.15_Rkind
    real(kind=Rkind) :: lambda=0.25_Rkind !bohr^-2
    real(kind=Rkind) :: mu = 1837.152_Rkind ! au 

    ! parameter of the reference, option=2
    real(kind=Rkind) :: b  = 0.1_Rkind ! instead of b
    real(kind=Rkind) :: w0 = -ONE ! so that k0=mu*w0^2

   CONTAINS
    PROCEDURE :: EvalPot_QModel  => EvalPot_QML_Bottleneck
    PROCEDURE :: Write_QModel    => Write_QML_Bottleneck
  END TYPE QML_Bottleneck_t

  PUBLIC :: QML_Bottleneck_t,Init_QML_Bottleneck

  CONTAINS
!> @brief Subroutine which makes the initialization of the Bottleneck parameters.
!!
!! @param BottleneckPot      TYPE(QML_Bottleneck_t):   derived type in which the parameters are set-up.
!! @param option             integer:            to be able to chose between the 3 models (default 1, Simple avoided crossing).
!! @param nio                integer (optional): file unit to read the parameters.
!! @param read_param         logical (optional): when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_Bottleneck(QModel_in,read_param,nio_param_file) RESULT(QModel)
    USE QDUtil_m,         ONLY : Identity_Mat
    IMPLICIT NONE

    !TYPE (QML_Bottleneck_t), allocatable         :: QModel
    TYPE (QML_Bottleneck_t)                      :: QModel

    TYPE(QML_Empty_t),           intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param

    !local variable
    integer                         :: err_read,nio_fit,i,j,k,idum

    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_Bottleneck'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unit,*) 'BEGINNING ',name_sub
      flush(out_unit)
    END IF

    !allocate(QML_Bottleneck_t :: QModel)

    QModel%QML_Empty_t = QModel_in
    !CALL Empty2_TO_Empty1(QModel%QML_Empty_t,QModel_in)

    QModel%pot_name = 'Bottleneck'
    QModel%nsurf    = 1
    QModel%ndim     = max(QModel%ndim,1)
    IF (QModel%ndim == 1) QModel%pot_name = 'Eckart'


    IF (QModel%option < 1 .OR. QModel%option > 2) QModel%option = 2

    SELECT CASE (QModel%option)
    CASE (1)
      QModel%V0     =  0.036_Rkind ! Hartree
      QModel%a      =  0.4_Rkind ! bohr^-1
      QModel%k0     =  0.1909_Rkind ! au
      QModel%sig    = 0.15_Rkind
      QModel%lambda = 0.25_Rkind !bohr^-2
      QModel%mu     = 1837.152_Rkind ! au 
    CASE (2)
      QModel%V0     = 0.425_Rkind / 27.211384_Rkind ! 0.425 eV
      QModel%a      = 1._Rkind ! bohr^-1

      QModel%b      = 0.1_Rkind
      QModel%w0     = -ONE
      QModel%mu     = 1060._Rkind ! au 

      QModel%k0     = QModel%mu * QModel%w0**2

      IF (read_param) CALL Read_QML_Bottleneck(QModel,nio_param_file)

      IF (QModel%w0 < ZERO) QModel%w0 = QModel%V0

    CASE Default
      write(out_unit,*) ' ERROR in EvalPot_QML_Bottleneck '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1 or 2'
      STOP 'ERROR in EvalPot_QML_Bottleneck: wrong option'
    END SELECT


    QModel%d0GGdef = Identity_Mat(QModel%ndim)/QModel%mu

    write(out_unit,*) 'PubliUnit=.TRUE./.FALSE.,  Q:in bohr, Energy: [Hartree]'

    IF (debug) write(out_unit,*) 'init Q0 of Bottleneck'
    allocate(QModel%Q0(QModel%ndim))
    CALL get_Q0_QML_Bottleneck(QModel%Q0,QModel,option=0)
    IF (debug) write(out_unit,*) 'QModel%Q0',QModel%Q0

    IF (debug) write(out_unit,*) 'init d0GGdef of Bottleneck'
    flush(out_unit)

    IF (debug) THEN
      write(out_unit,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unit,*) 'END ',name_sub
      flush(out_unit)
    END IF

  END FUNCTION Init_QML_Bottleneck
!> @brief Subroutine wich reads the Hénon-Heiles parameter with a namelist.
!!   This can be called only from the "Init_BottleneckPot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel            TYPE(QML_Bottleneck_t):   derived type in which the parameters are set-up.
!! @param nio               integer:                    file unit to read the parameters.
  SUBROUTINE Read_QML_Bottleneck(QModel,nio)
    IMPLICIT NONE

    TYPE (QML_Bottleneck_t), intent(inout)   :: QModel
    integer,                 intent(in)      :: nio

    real (kind=Rkind) :: sig,b,V0,a,mu,w0
    integer           :: err_read


    namelist /Bottleneck/ sig,b,V0,a,w0,mu

    sig = QModel%sig
    b   = QModel%b
    V0  = QModel%V0
    a   = QModel%a
    mu  = QModel%mu

    w0  = -ONE

    read(nio,nml=Bottleneck,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Bottleneck'
      write(out_unit,*) ' End-of-file or End-of-record'
      write(out_unit,*) ' The namelist "Bottleneck" is probably absent'
      write(out_unit,*) ' check your data!'
      write(out_unit,*)
      STOP ' ERROR in Read_BottleneckPot'
    ELSE IF (err_read > 0) THEN
      write(out_unit,*) ' ERROR in Read_QML_Bottleneck'
      write(out_unit,*) ' Some parameter names of the namelist "Bottleneck" are probaly wrong'
      write(out_unit,*) ' check your data!'
      write(out_unit,nml=Bottleneck)
      STOP ' ERROR in Read_BottleneckPot'
    END IF

    !write(out_unit,nml=Bottleneck)
    QModel%sig = sig
    QModel%b   = b
    QModel%V0  = V0
    QModel%a   = a
    QModel%mu  = mu
    QModel%w0  = w0

  END SUBROUTINE Read_QML_Bottleneck
  !> @brief Subroutine wich prints the QML_Bottleneck parameters.
!!
!! @param QModel            CLASS(QML_Bottleneck_t):   derived type in which the parameters are set-up.
!! @param nio               integer:            file unit to print the parameters.
  SUBROUTINE Write_QML_Bottleneck(QModel,nio)
    IMPLICIT NONE

    CLASS(QML_Bottleneck_t),     intent(in) :: QModel
    integer,                     intent(in) :: nio

    write(nio,*) 'Bottleneck current parameters'
    write(nio,*)
    write(nio,*) '---------------------------------------'
    write(nio,*) '  n dimensional potential              '
    write(nio,*) '                                       '
    write(nio,*) '  1st coordinates: s                   '
    write(nio,*) '             an Eckart Barrier         '
    write(nio,*) '                                       '
    write(nio,*) '  along the other coordinates: q2...qn '
    write(nio,*) '     a quadratic potential             '
    write(nio,*) '                                       '
    write(nio,*) '  V           (Hartree)                '
    write(nio,*) ' option 1:                             '
    write(nio,*) ' ref: Trahan, Wyatt and Poirier, J Chem Phys 122, 164104 (2005)'
    write(nio,*) '   Multidimensional quantum trajectories: Applications of the derivative propagation method.'
    write(nio,*) ' option 2 (default):                             '
    write(nio,*) ' ref: Dupuy, Lauvergnat and Scribano, CPL 787, 139241 (2022)'
    write(nio,*) '   Smolyak representations with absorbing boundary conditions ...'
    write(nio,*) '       for reaction path Hamiltonian model of reactive scattering.'
    write(nio,*) '   DOI: 10.1016/j.cplett.2021.139241'
    write(nio,*) '---------------------------------------'

    write(nio,*) '  PubliUnit:      ',QModel%PubliUnit
    write(nio,*)
    write(nio,*) '  Option   :      ',QModel%option
    write(nio,*)
    write(nio,*) '---------------------------------------'

    SELECT CASE (QModel%option)

    CASE (1)
      write(nio,*) '                                       '
      write(nio,*) ' V1(s) = V0 sech(a s)^2                '
      write(nio,*) '  with V0 =  0.036 Hartree             '
      write(nio,*) '  with a =  0.4 bohr^-1                '
      write(nio,*) '                                       '
      write(nio,*) ' V2(q,s) = sum_i 1/2.k(s).q^2          '
      write(nio,*) '   k(s) = k0(1 - sig.exp(- lambda.s^2))'
      write(nio,*) '  with k0 =  0.1909 au                 '
      write(nio,*) '  with sig=0.15                        '
      write(nio,*) '  with lambda=0.25 bohr^-2             '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) '   mu = 1837.152 au                    '
      write(nio,*) '                                       '
    CASE (2)
      write(nio,*) '                                       '
      write(nio,*) ' V1(s) = V0 sech(a s)^2                '
      write(nio,*) '  with V0 =  0.425 eV (converted in au)'
      write(nio,*) '  with a =  1. bohr^-1                 '
      write(nio,*) '                                       '
      write(nio,*) ' V2(q,s) = sum_i 1/2 k(s) q^2          '
      write(nio,*) '   k(s) = k0(1 + b sech(a s)^2)        '
      write(nio,*) '  with k0 =  mu*w0^2                   '
      write(nio,*) '  with b=0.1                           '
      write(nio,*) '  with w0 =  0.425 eV (converted in au)'
      write(nio,*) '  with lambda=0.25 bohr^-2             '
      write(nio,*) 'The constant metric tensor is obtained '
      write(nio,*) '   mu = 1060. au                       '
      write(nio,*) '                                       '
    CASE Default
        write(out_unit,*) ' ERROR in write_QModel '
        write(out_unit,*) ' This option is not possible. option: ',QModel%option
        write(out_unit,*) ' Its value MUST be 1 or 2'
        STOP
    END SELECT
    write(nio,*) '---------------------------------------'
    write(nio,*)
    write(nio,*) 'end Bottleneck current parameters'

  END SUBROUTINE Write_QML_Bottleneck

  SUBROUTINE get_Q0_QML_Bottleneck(Q0,QModel,option)
    IMPLICIT NONE

    real (kind=Rkind),           intent(inout) :: Q0(:)
    TYPE (QML_Bottleneck_t),     intent(in)    :: QModel
    integer,                     intent(in)    :: option

    Q0(:) = ZERO

  END SUBROUTINE get_Q0_QML_Bottleneck
!> @brief Subroutine wich calculates the Bottleneck potential (unpublished model) with derivatives up to the 2d order is required.
!!
!! @param PotVal             TYPE (dnMat_t):      derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param r                  real:                value for which the potential is calculated
!! @param QModel             TYPE(QML_Bottleneck_t):    derived type in which the parameters are set-up.
!! @param nderiv             integer:             it enables to specify up to which derivatives the potential is calculated:
!!                                                the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_Bottleneck(QModel,Mat_OF_PotDia,dnQ,nderiv)
    USE ADdnSVM_m
    IMPLICIT NONE

    CLASS(QML_Bottleneck_t),   intent(in)    :: QModel
    TYPE (dnS_t),              intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),              intent(in)    :: dnQ(:) !
    integer,                   intent(in)    :: nderiv

    TYPE (dnS_t) :: k

    SELECT CASE (QModel%option)
    CASE (1)
      k = QModel%k0*(ONE - QModel%sig * exp(-QModel%lambda * dnQ(1)**2))
    CASE (2)
      k = QModel%k0*(ONE + QModel%b / cosh(QModel%a * dnQ(1))**2 )
    CASE Default
      write(out_unit,*) ' ERROR in EvalPot_QML_Bottleneck '
      write(out_unit,*) ' This option is not possible. option: ',QModel%option
      write(out_unit,*) ' Its value MUST be 1 or 2'
      STOP 'ERROR in EvalPot_QML_Bottleneck: wrong option'
    END SELECT

    Mat_OF_PotDia(1,1) = QModel%V0 / cosh(QModel%a * dnQ(1))**2 + &
        HALF*k * sum(dnQ(2:QModel%ndim)**2)

  END SUBROUTINE EvalPot_QML_Bottleneck


END MODULE QML_Bottleneck_m
