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

!> @brief Module which makes the initialization, calculation of the Hénon-Heiles potential (value, gradient and hessian).
!!        Potential in nD (ndim) uses a parameter as defined in the following reference.
!> @brief Reference: M. Nest, H.-D. Meyer, J. Chem. Phys. 117 (2002) 10499. doi:10.1063/1.1521129
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
MODULE mod_HenonHeilesPot
  USE mod_NumParameters
  IMPLICIT NONE

!> @brief Derived type in which the parameter of the Hénon-Heiles potential are set-up.
!!        Potential in nD (ndim) uses a parameter as defined in the following reference.
!> @brief Reference: M. Nest, H.-D. Meyer, J. Chem. Phys. 117 (2002) 10499. doi:10.1063/1.1521129
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param lambda             real:    parameter of the Hénon-Heiles potential
!! @param ndim               integer: number of coordinates
  TYPE HenonHeilesPot_t
     PRIVATE
     real (kind=Rkind) :: lambda
     integer           :: ndim
  END TYPE HenonHeilesPot_t

  PRIVATE Read_HenonHeilesPot

CONTAINS
  SUBROUTINE Init_HenonHeilesPot(Para_HenonHeiles,ndim,nio,read_param,lambda)

    TYPE (HenonHeilesPot_t),    intent(inout)   :: Para_HenonHeiles
    integer,                     intent(in)      :: ndim
    integer,           optional, intent(in)      :: nio
    logical,           optional, intent(in)      :: read_param
    real (kind=Rkind), optional, intent(in)      :: lambda

    logical           :: read_param_loc
    real (kind=Rkind) :: lambda_loc



    read_param_loc = .FALSE.
    IF (present(read_param)) read_param_loc = read_param
    IF (read_param_loc .AND. .NOT. present(nio)) THEN
       write(out_unitp,*) ' ERROR in Init_HenonHeilesPot'
       write(out_unitp,*) ' read_param = t and The file unit (nio) is not present'
       write(out_unitp,*) ' => impossible to read the input file'
       STOP 'ERROR in Init_HenonHeilesPot: impossible to read the input file. The file unit (nio) is not present'
    END IF

    ! initalization of the default values
    lambda_loc     = 0.111803_Rkind

    IF (read_param_loc) THEN
      CALL Read_HenonHeilesPot(Para_HenonHeiles,ndim,nio,lambda_loc)
    ELSE

      IF (present(lambda))       lambda_loc     = lambda

      Para_HenonHeiles = HenonHeilesPot_t(lambda=lambda_loc,ndim=ndim)

    END IF

  END SUBROUTINE Init_HenonHeilesPot

!> @brief Subroutine wich reads the Hénon-Heiles parameter with a namelist.
!!   This can be called only from the "Init_HenonHeilesPot" subroutine.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_HenonHeiles            TYPE(HenonHeilesPot_t):    derived type in which the parameters are set-up.
!! @param nio                         integer:                    file unit to read the parameters.
!! @param lambdasub                   real (optional):            lambda parameters (default value for the namelist)
  SUBROUTINE Read_HenonHeilesPot(Para_HenonHeiles,ndim,nio,lambdasub)
    TYPE (HenonHeilesPot_t), intent(inout)  :: Para_HenonHeiles
    integer,                  intent(in)     :: nio,ndim
    real (kind=Rkind),        intent(in)     :: lambdasub

    real (kind=Rkind) :: lambda
    integer           :: err_read


    namelist /HenonHeiles/ lambda

    lambda = lambdasub

    read(nio,nml=HenonHeiles,IOSTAT=err_read)
    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_HenonHeilesPot'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "HenonHeiles" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_HenonHeilesPot'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_HenonHeilesPot'
      write(out_unitp,*) ' Some parameter names of the namelist "HenonHeiles" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=HenonHeiles)
      STOP ' ERROR in Read_HenonHeilesPot'
    END IF

    !write(out_unitp,nml=HenonHeiles)

    Para_HenonHeiles = HenonHeilesPot_t(lambda=lambda,ndim=ndim)


  END SUBROUTINE Read_HenonHeilesPot
!> @brief Subroutine wich prints the Hénon-Heiles potential parameters and 4D-calculations
!!
!> @author David Lauvergnat
!! @date 30/07/2019
!!
!! @param nio                integer:                   file unit to print the parameters.
  SUBROUTINE Write0_HenonHeilesPot(nio)
    integer,                  intent(in) :: nio


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

  END SUBROUTINE Write0_HenonHeilesPot
!> @brief Subroutine wich prints the Hénon-Heiles potential parameters.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param Para_HenonHeiles   TYPE(HenonHeilesPot_t):   derived type with the Phenol potential parameters.
!! @param nio                integer:                   file unit to print the parameters.
  SUBROUTINE Write_HenonHeilesPot(Para_HenonHeiles,nio)
    TYPE (HenonHeilesPot_t), intent(in) :: Para_HenonHeiles
    integer,                  intent(in) :: nio


    write(nio,*) 'HenonHeiles current parameters'
    write(nio,*)
    write(nio,*) '    ndim  :  ',Para_HenonHeiles%ndim
    write(nio,*) '    lambda:  ',Para_HenonHeiles%lambda
    write(nio,*)
    write(nio,*) 'end HenonHeiles current parameters'

  END SUBROUTINE Write_HenonHeilesPot

!> @brief Subroutine wich calculates the Hénon-Heiles potential with derivatives up to the 2d order if required.
!!
!> @author David Lauvergnat
!! @date 03/08/2017
!!
!! @param PotVal             TYPE (dnMat_t):          derived type with the potential (pot),  the gradient (grad) and the hessian (hess).
!! @param Q                  real:                    table of two values for which the potential is calculated (R,theta)
!! @param Para_HenonHeiles   TYPE(HenonHeilesPot_t): derived type with the Hénon-Heiles parameters.
!! @param nderiv             integer:                 it enables to specify up to which derivatives the potential is calculated:
!!                                                    the pot (nderiv=0) or pot+grad (nderiv=1) or pot+grad+hess (nderiv=2).
  SUBROUTINE Eval_HenonHeilesPot(Mat_OF_PotDia,Q,Para_HenonHeiles,nderiv)
    USE mod_dnS

    TYPE (HenonHeilesPot_t), intent(in)     :: Para_HenonHeiles
    TYPE (dnS_t),              intent(inout)  :: Mat_OF_PotDia(:,:)
    real(kind=Rkind),         intent(in)     :: Q(:)
    integer, intent(in)                      :: nderiv

    integer :: i
    real(kind=Rkind) :: d0,d1(Para_HenonHeiles%ndim),d2(Para_HenonHeiles%ndim,Para_HenonHeiles%ndim)

    ! Potential calculation
    d0 = ZERO
    DO i=1,Para_HenonHeiles%ndim
      d0 = d0 + HALF * Q(i)**2
    END DO
    DO i=1,Para_HenonHeiles%ndim-1
      d0 = d0 + Para_HenonHeiles%lambda * &
         ( Q(i)**2 * Q(i+1) - Q(i+1)**3/THREE )
    END DO

    ! Forces calculation
    IF (nderiv >= 1) THEN
      d1(:) = ZERO

      ! first the anharmonic part without lambda
      DO i=1,Para_HenonHeiles%ndim-1 ! derivative along the coodinate Q(i)
        d1(i)   = d1(i)   +  TWO * Q(i)*Q(i+1)
        d1(i+1) = d1(i+1) +  Q(i)**2 - Q(i+1)**2
      END DO

      d1(:) = d1(:) * Para_HenonHeiles%lambda

      ! Then the harmonic part
      DO i=1,Para_HenonHeiles%ndim ! derivative along the coodinate Q(k)
        d1(i) = d1(i) + Q(i)
      END DO

    END IF

    ! Hessian calculation
    IF (nderiv >= 2) THEN
      d2(:,:) = ZERO

      ! first the anharmonic part without lambda
      DO i=1,Para_HenonHeiles%ndim-1 ! derivative along the coodinate Q(i)
        d2(i  ,i  ) = d2(i  ,i  ) + TWO * Q(i+1)

        d2(i  ,i+1) = d2(i  ,i+1) + TWO * Q(i)
        d2(i+1,i  ) = d2(i+1,i  ) + TWO * Q(i)

        d2(i+1,i+1) = d2(i+1,i+1) - TWO * Q(i+1)

      END DO

      d2(:,:) = d2(:,:) * Para_HenonHeiles%lambda

      ! Then the harmonic part
      DO i=1,Para_HenonHeiles%ndim ! derivative along the coodinate Q(i)
        d2(i,i) = d2(i,i) + ONE
      END DO
    END IF

    IF (nderiv == 0) CALL QML_set_dnS(Mat_OF_PotDia(1,1),d0)
    IF (nderiv == 1) CALL QML_set_dnS(Mat_OF_PotDia(1,1),d0,d1)
    IF (nderiv == 2) CALL QML_set_dnS(Mat_OF_PotDia(1,1),d0,d1,d2)

  END SUBROUTINE Eval_HenonHeilesPot
  SUBROUTINE Eval_HenonHeilesPot_old(PotVal,Q,Para_HenonHeiles,nderiv)
    USE mod_dnMat

    TYPE (HenonHeilesPot_t),  intent(in)     :: Para_HenonHeiles
    TYPE (dnMat_t),           intent(inout)  :: PotVal
    real(kind=Rkind),         intent(in)     :: Q(:)
    integer, intent(in)                      :: nderiv

    integer :: i

    IF ( QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=1,ndim=Para_HenonHeiles%ndim,nderiv=nderiv)
    END IF

    ! Potential calculation
    PotVal%d0(1,1) = ZERO
    DO i=1,Para_HenonHeiles%ndim
      PotVal%d0(1,1) = PotVal%d0(1,1) + HALF * Q(i)**2
    END DO
    DO i=1,Para_HenonHeiles%ndim-1
      PotVal%d0(1,1) = PotVal%d0(1,1) + Para_HenonHeiles%lambda * &
         ( Q(i)**2 * Q(i+1) - Q(i+1)**3/THREE )
    END DO

    ! Forces calculation
    IF (nderiv >= 1) THEN
      PotVal%d1(1,1,:) = ZERO

      ! first the anharmonic part without lambda
      DO i=1,Para_HenonHeiles%ndim-1 ! derivative along the coodinate Q(i)
        PotVal%d1(1,1,i)   = PotVal%d1(1,1,i)   +  TWO * Q(i)*Q(i+1)
        PotVal%d1(1,1,i+1) = PotVal%d1(1,1,i+1) +  Q(i)**2 - Q(i+1)**2
      END DO

      PotVal%d1(1,1,:) = PotVal%d1(1,1,:) * Para_HenonHeiles%lambda

      ! Then the harmonic part
      DO i=1,Para_HenonHeiles%ndim ! derivative along the coodinate Q(k)
        PotVal%d1(1,1,i) = PotVal%d1(1,1,i) + Q(i)
      END DO

    END IF

    ! Hessian calculation
    IF (nderiv >= 2) THEN
      PotVal%d2(1,1,:,:) = ZERO

      ! first the anharmonic part without lambda
      DO i=1,Para_HenonHeiles%ndim-1 ! derivative along the coodinate Q(i)
        PotVal%d2(1,1,i  ,i  ) = PotVal%d2(1,1,i  ,i  ) + TWO * Q(i+1)

        PotVal%d2(1,1,i  ,i+1) = PotVal%d2(1,1,i  ,i+1) + TWO * Q(i)
        PotVal%d2(1,1,i+1,i  ) = PotVal%d2(1,1,i+1,i  ) + TWO * Q(i)

        PotVal%d2(1,1,i+1,i+1) = PotVal%d2(1,1,i+1,i+1) - TWO * Q(i+1)

      END DO

      PotVal%d2(1,1,:,:) = PotVal%d2(1,1,:,:) * Para_HenonHeiles%lambda

      ! Then the harmonic part
      DO i=1,Para_HenonHeiles%ndim ! derivative along the coodinate Q(i)
        PotVal%d2(1,1,i,i) = PotVal%d2(1,1,i,i) + ONE
      END DO
    END IF
  
  END SUBROUTINE Eval_HenonHeilesPot_old
END MODULE mod_HenonHeilesPot
