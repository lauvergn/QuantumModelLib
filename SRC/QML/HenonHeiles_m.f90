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

!> @brief Module which makes the initialization, calculation of the HenonHeiles potentials (value, gradient and hessian).
!!
!> @author David Lauvergnat
!! @date 07/01/2020
!!
MODULE QML_HenonHeiles_m
  USE QMLLib_NumParameters_m
  USE QML_Empty_m
  IMPLICIT NONE

  PRIVATE

!> @brief Derived type in which the HenonHeiles parameters are set-up.
  TYPE, EXTENDS (QML_Empty_t) ::  QML_HenonHeiles_t
   PRIVATE

     real (kind=Rkind) :: lambda = 0.111803_Rkind

   CONTAINS
    !PROCEDURE :: EvalPot_QModel => EvalPot_QML_HenonHeiles
    PROCEDURE :: EvalPot_QModel => EvalPotnew_QML_HenonHeiles
    PROCEDURE :: Write_QModel    => Write_QML_HenonHeiles
    PROCEDURE :: Write0_QModel   => Write0_QML_HenonHeiles
  END TYPE QML_HenonHeiles_t

  PUBLIC :: QML_HenonHeiles_t,Init_QML_HenonHeiles

  CONTAINS
!> @brief Function which makes the initialization of the HenonHeiles parameters.
!!
!! @param QModel             TYPE(QML_HenonHeiles_t):   result derived type in which the parameters are set-up.
!! @param QModel_in          TYPE(QML_Empty_t):  type to transfer ndim, nsurf ...
!! @param nio_param_file     integer:             file unit to read the parameters.
!! @param read_param         logical:             when it is .TRUE., the parameters are read. Otherwise, they are initialized.
  FUNCTION Init_QML_HenonHeiles(QModel_in,read_param,nio_param_file,lambda) RESULT(QModel)
  IMPLICIT NONE

    TYPE (QML_HenonHeiles_t)                    :: QModel ! RESULT

    TYPE(QML_Empty_t),          intent(in)      :: QModel_in ! variable to transfer info to the init
    integer,                     intent(in)      :: nio_param_file
    logical,                     intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: lambda


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Init_QML_HenonHeiles'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      CALL QModel_in%Write_QModel(out_unitp)
      flush(out_unitp)
    END IF

    CALL Init0_QML_Empty(QModel%QML_Empty_t,QModel_in)

    QModel%nsurf    = 1
    QModel%pot_name = 'henonheiles'

    ! initalization of the default values
    QModel%lambda     = 0.111803_Rkind

    IF (read_param) THEN
      CALL Read_QML_HenonHeiles(QModel,nio_param_file)
    ELSE
      IF (present(lambda))  QModel%lambda = lambda
    END IF

    IF (QModel%ndim < 1) THEN
      write(out_unitp,*) ' ERROR in Init_QML_HenonHeiles'
      write(out_unitp,*) ' ndim MUST be > 0. ndim: ',QModel%ndim
      write(out_unitp,*) ' Its value MUST be given in with Init_Model or Read_Model subroutines.'
      STOP 'ERROR in Init_QML_HenonHeiles: Wrong ndim value.'
    END IF


    IF (debug) write(out_unitp,*) 'init Q0 of HenonHeiles'
    allocate(QModel%Q0(QModel%ndim))
    QModel%Q0 = ZERO

    IF (debug) write(out_unitp,*) 'init d0GGdef of HenonHeiles'
    CALL Init_IdMat(QModel%d0GGdef,QModel%ndim)


    IF (debug) THEN
      write(out_unitp,*) 'QModel%pot_name: ',QModel%pot_name
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END FUNCTION Init_QML_HenonHeiles
!> @brief Subroutine wich reads the Hénon-Heiles parameter with a namelist.
!!   This can be called only from the "Init_HenonHeilesPot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param QModel            TYPE(QML_HenonHeiles_t):   derived type in which the parameters are set-up.
!! @param nio               integer:                    file unit to read the parameters.
  SUBROUTINE Read_QML_HenonHeiles(QModel,nio)
    TYPE (QML_HenonHeiles_t), intent(inout)   :: QModel
    integer,                   intent(in)      :: nio

    real (kind=Rkind) :: lambda
    integer           :: err_read


    namelist /HenonHeiles/ lambda

    lambda = QModel%lambda

    read(nio,nml=HenonHeiles,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_QML_HenonHeiles'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "HenonHeiles" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_HenonHeilesPot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_QML_HenonHeiles'
      write(out_unitp,*) ' Some parameter names of the namelist "HenonHeiles" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=HenonHeiles)
      STOP ' ERROR in Read_HenonHeilesPot'
    END IF

    !write(out_unitp,nml=HenonHeiles)
    QModel%lambda = lambda

  END SUBROUTINE Read_QML_HenonHeiles
!> @brief Subroutine wich prints the QML_HenonHeiles parameters.
!!
!! @param QModel            CLASS(QML_HenonHeiles_t): derived type in which the parameters are set-up.
!! @param nio               integer:                   file unit to print the parameters.
  SUBROUTINE Write_QML_HenonHeiles(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HenonHeiles_t),   intent(in) :: QModel
    integer,                     intent(in) :: nio


    write(nio,*) 'HenonHeiles current parameters'
    write(nio,*)
    write(nio,*) '    ndim  :  ',QModel%ndim
    write(nio,*) '    lambda:  ',QModel%lambda
    write(nio,*)
    write(nio,*) 'end HenonHeiles current parameters'

  END SUBROUTINE Write_QML_HenonHeiles
  SUBROUTINE Write0_QML_HenonHeiles(QModel,nio)
  IMPLICIT NONE

    CLASS(QML_HenonHeiles_t),   intent(in) :: QModel
    integer,                     intent(in) :: nio

    write(nio,*) 'HenonHeiles default parameters'
    write(nio,*)
    write(nio,*) ' Potential and parameters from M. Nest, H.-D. Meyer, J. Chem. Phys. 117 (2002) 10499'
    write(nio,*) '    lambda = 0.111803'
    write(nio,*)
    write(nio,*) '4D-Potential Value at: Q:'
    write(nio,*) '      0.1000000000  0.2000000000  0.3000000000  0.4000000000'
    write(nio,*) 'V = 0.151900651'
    write(nio,*) 'gradient = [0.1044721200  0.2100622700  0.3212425700  0.3921737900]'
    write(nio,*) 'hessian'
    write(nio,*) '1        1.044721       0.022361       0.000000       0.000000'
    write(nio,*) '2        0.022361       1.022361       0.044721       0.000000'
    write(nio,*) '3        0.000000       0.044721       1.022361       0.067082'
    write(nio,*) '4        0.000000       0.000000       0.067082       0.910558'
    write(nio,*)
    write(nio,*) ' Energy levels obtained with ElVibRot100.1-Tnum30.37-Tana8.7'
    write(nio,*) '  Using direct-product basis set and grid.'
    write(nio,*) '  Primitive basis/grid for each coordinate: '
    write(nio,*) '    HO, Harmonic oscilator, with 10 basis functions'
    write(nio,*) '    Gauss Hermite quadrature with 10 grid points.'

    write(nio,*)  ' The first five energy levels:'
    write(nio,*)  ' 1.995725  2.972602  2.980991  2.987181  2.989254'
    write(nio,*)  '    precision around 10^-6'
    write(nio,*)
    write(nio,*) 'end HenonHeiles default parameters'

  END SUBROUTINE Write0_QML_HenonHeiles

!> @brief Subroutine wich calculates the HenonHeiles potential with derivatives up to the 2d order.
!!
!! @param QModel             TYPE(QML_HenonHeiles_t): derived type in which the parameters are set-up.
!! @param Mat_OF_PotDia(:,:) TYPE (dnS_t):             derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param dnQ(:)             TYPE (dnS_t)              value for which the potential is calculated
!! @param nderiv             integer:                  it enables to specify up to which derivatives the potential is calculated:
!!                                                     the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE EvalPot_QML_HenonHeiles(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HenonHeiles_t), intent(in)    :: QModel
    TYPE (dnS_t),              intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),              intent(in)    :: dnQ(:)
    integer,                   intent(in)    :: nderiv


    integer :: i
    real(kind=Rkind) :: d0,d1(QModel%ndim),d2(QModel%ndim,QModel%ndim)
    real(kind=Rkind) :: Q(QModel%ndim)

    Q(:) = get_d0(dnQ)

    ! Potential calculation
    d0 = ZERO
    DO i=1,QModel%ndim
      d0 = d0 + HALF * Q(i)**2
    END DO
    DO i=1,QModel%ndim-1
      d0 = d0 + QModel%lambda * &
         ( Q(i)**2 * Q(i+1) - Q(i+1)**3/THREE )
    END DO

    ! Forces calculation
    IF (nderiv >= 1) THEN
      d1(:) = ZERO

      ! first the anharmonic part without lambda
      DO i=1,QModel%ndim-1 ! derivative along the coodinate Q(i)
        d1(i)   = d1(i)   +  TWO * Q(i)*Q(i+1)
        d1(i+1) = d1(i+1) +  Q(i)**2 - Q(i+1)**2
      END DO

      d1(:) = d1(:) * QModel%lambda

      ! Then the harmonic part
      DO i=1,QModel%ndim ! derivative along the coodinate Q(k)
        d1(i) = d1(i) + Q(i)
      END DO

    END IF

    ! Hessian calculation
    IF (nderiv >= 2) THEN
      d2(:,:) = ZERO

      ! first the anharmonic part without lambda
      DO i=1,QModel%ndim-1 ! derivative along the coodinate Q(i)
        d2(i  ,i  ) = d2(i  ,i  ) + TWO * Q(i+1)

        d2(i  ,i+1) = d2(i  ,i+1) + TWO * Q(i)
        d2(i+1,i  ) = d2(i+1,i  ) + TWO * Q(i)

        d2(i+1,i+1) = d2(i+1,i+1) - TWO * Q(i+1)

      END DO

      d2(:,:) = d2(:,:) * QModel%lambda

      ! Then the harmonic part
      DO i=1,QModel%ndim ! derivative along the coodinate Q(i)
        d2(i,i) = d2(i,i) + ONE
      END DO
    END IF

    IF (nderiv == 0) CALL set_dnS(Mat_OF_PotDia(1,1),d0)
    IF (nderiv == 1) CALL set_dnS(Mat_OF_PotDia(1,1),d0,d1)
    IF (nderiv == 2) CALL set_dnS(Mat_OF_PotDia(1,1),d0,d1,d2)

  END SUBROUTINE EvalPot_QML_HenonHeiles

  SUBROUTINE EvalPotnew_QML_HenonHeiles(QModel,Mat_OF_PotDia,dnQ,nderiv)
  USE ADdnSVM_m
  IMPLICIT NONE

    CLASS(QML_HenonHeiles_t), intent(in)    :: QModel
    TYPE (dnS_t),              intent(inout) :: Mat_OF_PotDia(:,:)
    TYPE (dnS_t),              intent(in)    :: dnQ(:)
    integer,                   intent(in)    :: nderiv


    integer          :: i
    TYPE (dnS_t)     :: dnV


    ! Potential calculation
    dnV = dnQ(1) ! for the initialization
    dnV = ZERO
    DO i=1,QModel%ndim
      dnV = dnV + HALF * dnQ(i)**2
    END DO
    DO i=1,QModel%ndim-1
      dnV = dnV + QModel%lambda*( dnQ(i)**2 * dnQ(i+1) - dnQ(i+1)**3/THREE )
    END DO

    Mat_OF_PotDia(1,1) = dnV

    CALL dealloc_dnS(dnV)

  END SUBROUTINE EvalPotnew_QML_HenonHeiles

END MODULE QML_HenonHeiles_m
