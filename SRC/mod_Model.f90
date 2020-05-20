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
MODULE mod_FiniteDiff_coeff
  USE mod_NumParameters
  IMPLICIT NONE

  PRIVATE

    !--------------------------------------------------------------------------
    integer, parameter :: list_1D_indDQ(1,4) = reshape([-2,-1, 1, 2],shape=[1,4])
    integer, parameter :: list_2D_indDQ(2,8) = reshape([                &
                                                       -1,-1,           &
                                                        1,-1,           &
                                                       -1, 1,           &
                                                        1, 1,           &
                                                       -2,-2,           &
                                                        2,-2,           &
                                                       -2, 2,           &
                                                        2, 2            &
                                                            ],shape=[2,8])
    integer, parameter :: list_3D_indDQ(3,6) = reshape([                &
                                                       -1,-1, 1,        &
                                                       -1, 1,-1,        &
                                                        1,-1,-1,        &
                                                        1, 1,-1,        &
                                                        1,-1, 1,        &
                                                       -1, 1, 1         &
                                                            ],shape=[3,6])
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !for 1st derivatives
    ! for d./dx
    real (kind=Rkind), parameter ::   wD(-2:2) = [ONE/TWELVE,-EIGHT/TWELVE,ZERO,EIGHT/TWELVE,-ONE/TWELVE]
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !for 2d derivatives
    ! for d2./dx2
    real (kind=Rkind), parameter ::  wDD(-2:2) = [-ONE/TWELVE,FOUR/THREE,-FIVE/TWO,FOUR/THREE,-ONE/TWELVE]

    ! for d2./dxdy
    real (kind=Rkind), parameter ::  w11DD(-1:1,-1:1) = reshape([       &
                                            ONE/THREE,ZERO,-ONE/THREE,  &
                                             ZERO,    ZERO, ZERO,       &
                                           -ONE/THREE,ZERO, ONE/THREE   &
                                                              ],shape=[3,3])

    real (kind=Rkind), parameter ::  w22DD(-1:1,-1:1) = reshape([       &
                                     -ONE/48._Rkind,ZERO, ONE/48._Rkind,&
                                          ZERO,     ZERO, ZERO,         &
                                      ONE/48._Rkind,ZERO,-ONE/48._Rkind &
                                                              ],shape=[3,3])
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !for 3d derivatives
    ! for d3./dx3
    real (kind=Rkind), parameter :: wDDD(-2:2) = [-HALF,ONE,ZERO,-ONE,HALF]

    ! for d3./dx2dy
    real (kind=Rkind), parameter :: a = TWO/THREE
    real (kind=Rkind), parameter :: b = ONE/48._Rkind

    real (kind=Rkind), parameter ::  w11DDD(-1:1,-1:1) = reshape([       &
                                             -a,  ZERO,-a,               &
                                             ZERO,ZERO, ZERO,            &
                                              a,  ZERO, a                &
                                                              ],shape=[3,3])

    real (kind=Rkind), parameter ::  w22DDD(-1:1,-1:1) = reshape([       &
                                              b,  ZERO, b,               &
                                             ZERO,ZERO, ZERO,            &
                                             -b,  ZERO,-b                &
                                                              ],shape=[3,3])
    real (kind=Rkind), parameter :: w0aDDD(-2:2) = [-b-b,a+a,ZERO,-a-a,b+b]

    ! for d3./dxdydz
    real (kind=Rkind), parameter :: wxyz = ONE/SIX

    real (kind=Rkind), parameter :: w100DDD(-2:2) = [ZERO,wxyz,ZERO,-wxyz,ZERO]
    real (kind=Rkind), parameter :: w110DDD(-1:1,-1:1) = reshape([      &
                                            -wxyz,ZERO,ZERO,            &
                                             ZERO,ZERO,ZERO,            &
                                             ZERO,ZERO,wxyz             &
                                                              ],shape=[3,3])

    real (kind=Rkind) :: w111DDD(-1:1,-1:1,-1:1) = reshape([            &
                                                    ZERO, & ! -1 -1 -1
                                                    ZERO, & !  0 -1 -1
                                                    wxyz, & !  1 -1 -1
                                                    ZERO, & ! -1  0 -1
                                                    ZERO, & !  0  0 -1
                                                    ZERO, & !  1  0 -1
                                                    wxyz, & ! -1  1 -1
                                                    ZERO, & !  0  1 -1
                                                   -wxyz, & !  1  1 -1
                                                    ZERO, & ! -1 -1  0
                                                    ZERO, & !  0 -1  0
                                                    ZERO, & !  1 -1  0
                                                    ZERO, & ! -1  0  0
                                                    ZERO, & !  0  0  0
                                                    ZERO, & !  1  0  0
                                                    ZERO, & ! -1  1  0
                                                    ZERO, & !  0  1  0
                                                    ZERO, & !  1  1  0
                                                    wxyz, & ! -1 -1  1
                                                    ZERO, & !  0 -1  1
                                                   -wxyz, & !  1 -1  1
                                                    ZERO, & ! -1  0  1
                                                    ZERO, & !  0  0  1
                                                    ZERO, & !  1  0  1
                                                   -wxyz, & ! -1  1  1
                                                    ZERO, & !  0  1  1
                                                    ZERO  & !  1  1  1
                                                           ],shape=[3,3,3])
    !--------------------------------------------------------------------------

PUBLIC :: Get_nb_pts,Get_indDQ,NumericDerivativesQ_AddMat_TO_dnMat

CONTAINS
  FUNCTION Get_nb_pts(ndim) RESULT(nb_pts)
  IMPLICIT NONE
    integer                       :: nb_pts   ! number of points
    integer,           intent(in) :: ndim     ! indDQ(ndim)

    SELECT CASE (ndim)
    CASE (0)
      nb_pts = 1
    CASE (1)
      nb_pts  = size(list_1D_indDQ,dim=2)
    CASE (2)
      nb_pts  = size(list_2D_indDQ,dim=2)
    CASE (3)
      nb_pts  = size(list_3D_indDQ,dim=2)
    CASE Default
      nb_pts  = -1
    END SELECT

  END FUNCTION Get_nb_pts
  SUBROUTINE Get_indDQ(indDQ,i_pt)
  IMPLICIT NONE

    integer,           intent(inout) :: indDQ(:) ! amplitude of DQ

    integer,           intent(in)    :: i_pt     ! when 0, force the initialization (number of points)


    SELECT CASE (size(indDQ))
    CASE (1)
      indDQ(:) = list_1D_indDQ(:,i_pt)
    CASE (2)
      indDQ(:) = list_2D_indDQ(:,i_pt)
    CASE (3)
      indDQ(:) = list_3D_indDQ(:,i_pt)
    END SELECT

  END SUBROUTINE Get_indDQ

  SUBROUTINE NumericDerivativesQ_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ,option)
  USE mod_dnMat
  IMPLICIT NONE

    TYPE (dnMat_t),    intent(inout)         :: dnMat
    real (kind=Rkind), intent(in)            :: Mat(:,:)
    integer,           intent(in),  optional :: indQ(:)  ! indexes of the variables along DQ is made
    integer,           intent(in),  optional :: indDQ(:) ! amplitude of DQ
    integer,           intent(in),  optional :: option   ! version 3 or 4 (default 3)

    integer:: option_loc

    IF (present(option)) THEN
      option_loc = option
    ELSE
      option_loc = 3
    END IF

    IF (option_loc == 4) THEN
      IF (present(indQ) .AND. present(indDQ)) THEN

        CALL NumericDerivatives4Q_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ)

      ELSE IF (.NOT. present(indQ) .AND. .NOT. present(indDQ)) THEN

        CALL NumericDerivatives4Q_AddMat_TO_dnMat(dnMat,Mat)

      ELSE
        write(out_unitp,*) ' ERROR in NumericDerivativesQ_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' Both indQ and indDQ MUST be present'
        write(out_unitp,*) '     or '
        write(out_unitp,*) ' Both indQ and indDQ MUST be absent'
        write(out_unitp,*) ' present(indQ) ',present(indQ)
        write(out_unitp,*) ' present(indDQ)',present(indDQ)
        STOP 'STOP in NumericDerivativesQ_AddMat_TO_dnMat: Inconsitent parameters.'
      END IF
    ELSE ! option_loc == 3
      IF (present(indQ) .AND. present(indDQ)) THEN

        CALL NumericDerivatives3Q_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ)

      ELSE IF (.NOT. present(indQ) .AND. .NOT. present(indDQ)) THEN

        CALL NumericDerivatives3Q_AddMat_TO_dnMat(dnMat,Mat)

      ELSE
        write(out_unitp,*) ' ERROR in NumericDerivativesQ_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' Both indQ and indDQ MUST be present'
        write(out_unitp,*) '     or '
        write(out_unitp,*) ' Both indQ and indDQ MUST be absent'
        write(out_unitp,*) ' present(indQ) ',present(indQ)
        write(out_unitp,*) ' present(indDQ)',present(indDQ)
        STOP 'STOP in NumericDerivativesQ_AddMat_TO_dnMat: Inconsitent parameters.'
      END IF
    END IF

  END SUBROUTINE NumericDerivativesQ_AddMat_TO_dnMat

  SUBROUTINE NumericDerivatives4Q_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ)
  USE mod_dnMat
  IMPLICIT NONE

    TYPE (dnMat_t),    intent(inout)         :: dnMat
    real (kind=Rkind), intent(in)            :: Mat(:,:)
    integer,           intent(in),  optional :: indQ(:)  ! indexes of the variables along DQ is made
    integer,           intent(in),  optional :: indDQ(:) ! amplitude of DQ

    ! local variable
    integer                            :: size_indQ,ndim,nderiv
    integer                            :: i,j,k,ip,jp,kp

    nderiv = QML_get_nderiv_FROM_dnMat(dnMat)
    ndim   = QML_get_ndim_FROM_dnMat(dnMat)

    IF (nderiv < 0) THEN
      write(out_unitp,*) ' ERROR in NumericDerivatives4Q_AddMat_TO_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' nderiv < 0',nderiv
      write(out_unitp,*) '  => dnMat is not allocated'
      STOP 'STOP in NumericDerivatives4Q_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (.NOT. all(shape(Mat) == shape(dnMat%d0))) THEN
      write(out_unitp,*) ' ERROR in NumericDerivatives4Q_AddMat_TO_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' shape(Mat),shape(dnMat%d0)',shape(Mat),shape(dnMat%d0)
      STOP 'STOP in NumericDerivatives4Q_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (present(indQ) .AND. present(indDQ)) THEN
      IF (size(indQ) /= size(indDQ) .OR. minval(indQ) < 1 .OR. maxval(indQ) > ndim) THEN
        write(out_unitp,*) ' ERROR in NumericDerivatives4Q_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' size(indQ),size(indDQ)',size(indQ),size(indDQ)
        write(out_unitp,*) ' ndim',ndim
        write(out_unitp,*) ' indQ(:)',indQ
        write(out_unitp,*) ' indDQ(:)',indDQ
        STOP 'STOP in NumericDerivatives4Q_AddMat_TO_dnMat: Inconsitent parameters.'
      END IF

    ELSE IF ( (.NOT. present(indQ) .AND.       present(indDQ)) .OR.     &
              (      present(indQ) .AND. .NOT. present(indDQ)) ) THEN
        write(out_unitp,*) ' ERROR in NumericDerivatives4Q_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' Both indQ and indDQ MUST be present'
        write(out_unitp,*) '     or '
        write(out_unitp,*) ' Both indQ and indDQ MUST be absent'
        write(out_unitp,*) ' present(indQ) ',present(indQ)
        write(out_unitp,*) ' present(indDQ)',present(indDQ)
        STOP 'STOP in NumericDerivatives4Q_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (present(indQ)) THEN
      size_indQ = size(indQ)
    ELSE
      size_indQ = 0
    END IF

    SELECT CASE (size_indQ)
    CASE (0)
      dnMat%d0 = Mat

      IF (nderiv >= 2) THEN
        DO i=1,ndim
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,wDD(0),dnMat,ider=[i,i])
        END DO
      END IF

    CASE (1)
      i  = indQ(1)
      ip = indDQ(1)

      CALL QML_Mat_wADDTO_dnMat2_ider(Mat,wD(ip)  ,dnMat,ider=[i])

      IF (nderiv >= 2) THEN
        CALL QML_Mat_wADDTO_dnMat2_ider(Mat,wDD(ip) ,dnMat,ider=[i,i])
      END IF

      IF (nderiv >= 3) THEN
        ! d3/dQidQidQi
        CALL QML_Mat_wADDTO_dnMat2_ider(Mat,wDDD(ip),dnMat,ider=[i,i,i])

        ! d3/dQjdQjdQi
        DO j=1,ndim
          IF (i == j) CYCLE
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w0aDDD(ip),dnMat,ider=[j,j,i])
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w0aDDD(ip),dnMat,ider=[j,i,j])
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w0aDDD(ip),dnMat,ider=[i,j,j])
        END DO

        ! d3/dQkdQjdQi
        DO j=1,ndim
        DO k=1,ndim
          IF (i == j .OR. i == k .OR. j == k) CYCLE
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w100DDD(ip),dnMat,ider=[k,j,i])
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w100DDD(ip),dnMat,ider=[k,i,j])
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w100DDD(ip),dnMat,ider=[i,k,j])
        END DO
        END DO

      END IF

    CASE (2)

      i  = indQ(1)
      j  = indQ(2)
      ip = indDQ(1)
      jp = indDQ(2)

      IF (abs(ip) == 1 .AND. abs(jp) == 1) THEN
        ! d2/dQidQj
        CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w11DD(jp,ip),dnMat,ider=[j,i])

        IF (nderiv >= 3) THEN
          ! d3/dQidQidQj
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w11DDD(ip,jp),dnMat,ider=[i,i,j])
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w11DDD(ip,jp),dnMat,ider=[i,j,i])
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w11DDD(ip,jp),dnMat,ider=[j,i,i])
        END IF

        IF (nderiv >= 3 .AND. ndim > 2) THEN
          ! d3/dQidQjdQk
          DO k=1,ndim
            IF (k == i .OR. k == j) CYCLE
            CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w110DDD(ip,jp),dnMat,ider=[k,j,i])
            CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w110DDD(ip,jp),dnMat,ider=[j,k,i])
            CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w110DDD(ip,jp),dnMat,ider=[j,i,k])
          END DO
        END IF

      ELSE IF (abs(ip) == 2 .AND. abs(jp) == 2) THEN

        ! d2/dQidQj
        CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w22DD(jp/2,ip/2),dnMat,ider=[j,i])

        IF (nderiv >= 3) THEN
          ! d3/dQidQidQj
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w22DDD(ip/2,jp/2),dnMat,ider=[i,i,j])
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w22DDD(ip/2,jp/2),dnMat,ider=[i,j,i])
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w22DDD(ip/2,jp/2),dnMat,ider=[j,i,i])
        END IF
      END IF

    CASE (3)

      i  = indQ(1)
      j  = indQ(2)
      k  = indQ(3)

      ip = indDQ(1)
      jp = indDQ(2)
      kp = indDQ(3)

      CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w111DDD(kp,jp,ip) ,dnMat,ider=[k,j,i])

    CASE Default
    END SELECT

  END SUBROUTINE NumericDerivatives4Q_AddMat_TO_dnMat
  SUBROUTINE NumericDerivatives3Q_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ)
  USE mod_dnMat
  IMPLICIT NONE

    TYPE (dnMat_t),    intent(inout)         :: dnMat
    real (kind=Rkind), intent(in)            :: Mat(:,:)
    integer,           intent(in),  optional :: indQ(:)  ! indexes of the variables along DQ is made
    integer,           intent(in),  optional :: indDQ(:) ! amplitude of DQ

    ! local variable
    integer                            :: size_indQ,ndim,nderiv
    integer                            :: i,j,k,ip,jp,kp

    nderiv = QML_get_nderiv_FROM_dnMat(dnMat)
    ndim   = QML_get_ndim_FROM_dnMat(dnMat)

    IF (nderiv < 0) THEN
      write(out_unitp,*) ' ERROR in NumericDerivatives3Q_AddMat_TO_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' nderiv < 0',nderiv
      write(out_unitp,*) '  => dnMat is not allocated'
      STOP 'STOP in NumericDerivatives3Q_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (.NOT. all(shape(Mat) == shape(dnMat%d0))) THEN
      write(out_unitp,*) ' ERROR in NumericDerivatives3Q_AddMat_TO_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' shape(Mat),shape(dnMat%d0)',shape(Mat),shape(dnMat%d0)
      STOP 'STOP in NumericDerivatives3Q_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (present(indQ) .AND. present(indDQ)) THEN
      IF (size(indQ) /= size(indDQ) .OR. minval(indQ) < 1 .OR. maxval(indQ) > ndim) THEN
        write(out_unitp,*) ' ERROR in NumericDerivatives3Q_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' size(indQ),size(indDQ)',size(indQ),size(indDQ)
        write(out_unitp,*) ' ndim',ndim
        write(out_unitp,*) ' indQ(:)',indQ
        write(out_unitp,*) ' indDQ(:)',indDQ
        STOP 'STOP in NumericDerivatives3Q_AddMat_TO_dnMat: Inconsitent parameters.'
      END IF

    ELSE IF ( (.NOT. present(indQ) .AND.       present(indDQ)) .OR.     &
              (      present(indQ) .AND. .NOT. present(indDQ)) ) THEN
        write(out_unitp,*) ' ERROR in NumericDerivatives3Q_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' Both indQ and indDQ MUST be present'
        write(out_unitp,*) '     or '
        write(out_unitp,*) ' Both indQ and indDQ MUST be absent'
        write(out_unitp,*) ' present(indQ) ',present(indQ)
        write(out_unitp,*) ' present(indDQ)',present(indDQ)
        STOP 'STOP in NumericDerivatives3Q_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (present(indQ)) THEN
      size_indQ = size(indQ)
    ELSE
      size_indQ = 0
    END IF

    SELECT CASE (size_indQ)
    CASE (0)
      dnMat%d0 = Mat

      IF (nderiv >= 2) THEN
        DO i=1,ndim
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,wDD(0),dnMat,ider=[i,i])
        END DO
      END IF

    CASE (1)
      i  = indQ(1)
      ip = indDQ(1)

       CALL QML_Mat_wADDTO_dnMat2_ider(Mat,wD(ip)  ,dnMat,ider=[i])

       IF (nderiv >= 2) THEN
         CALL QML_Mat_wADDTO_dnMat2_ider(Mat,wDD(ip) ,dnMat,ider=[i,i])
       END IF

       IF (nderiv >= 3) THEN
        ! d3/dQidQidQi
         CALL QML_Mat_wADDTO_dnMat2_ider(Mat,wDDD(ip),dnMat,ider=[i,i,i])

         ! d3/dQjdQjdQi
         DO j=1,ndim
           IF (i == j) CYCLE
           CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w0aDDD(ip),dnMat,ider=[j,j,i])
         END DO
       END IF

       ! d3/dQkdQjdQi
       IF (nderiv >= 3 .AND. ndim > 2 .AND. abs(ip) == 1) THEN

         DO j=1,ndim
         DO k=1,ndim
           IF (k > j .AND. j > i)  THEN
             CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w100DDD(ip),dnMat,ider=[k,j,i])
           END IF

           !permutation between i and j (order: k,i,j)
           IF (k > i .AND. i > j)  THEN
             CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w100DDD(ip),dnMat,ider=[k,i,j])
           END IF

           !permutation between i and k (order: i,k,j)
           IF (i > k .AND. k > j)  THEN
             CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w100DDD(ip),dnMat,ider=[i,k,j])
           END IF

         END DO
         END DO
       END IF

    CASE (2)

      i  = indQ(1)
      j  = indQ(2)
      ip = indDQ(1)
      jp = indDQ(2)

      IF (abs(ip) == 1 .AND. abs(jp) == 1) THEN
        ! d2/dQidQj
        CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w11DD(jp,ip),dnMat,ider=[j,i])

        IF (nderiv >= 3) THEN
          ! d3/dQidQidQj
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w11DDD(ip,jp),dnMat,ider=[i,i,j])
          ! d3/dQjdQjdQi
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w11DDD(jp,ip),dnMat,ider=[j,j,i])

          ! d3/dQidQjdQk
          DO k=1,ndim
            !IF (k == i .OR. k == j) CYCLE
            IF (k > j .AND. j > i) THEN ! remark: in this version, j>i always
              CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w110DDD(ip,jp),dnMat,ider=[k,j,i])
            END IF
            IF (j > k .AND. k > i) THEN
              CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w110DDD(ip,jp),dnMat,ider=[j,k,i])
            END IF
            IF (j > i .AND. i > k) THEN
              CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w110DDD(ip,jp),dnMat,ider=[j,i,k])
            END IF
          END DO
        END IF

      ELSE IF (abs(ip) == 2 .AND. abs(jp) == 2) THEN

        ! d2/dQidQj
        CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w22DD(jp/2,ip/2),dnMat,ider=[j,i])

        IF (nderiv >= 3) THEN
          ! d3/dQidQidQj
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w22DDD(ip/2,jp/2),dnMat,ider=[i,i,j])
          ! d3/dQjdQjdQi
          CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w22DDD(jp/2,ip/2),dnMat,ider=[j,j,i])
        END IF
      END IF

    CASE (3)

      i  = indQ(1)
      j  = indQ(2)
      k  = indQ(3)

      ip = indDQ(1)
      jp = indDQ(2)
      kp = indDQ(3)

      CALL QML_Mat_wADDTO_dnMat2_ider(Mat,w111DDD(kp,jp,ip) ,dnMat,ider=[k,j,i])

    CASE Default
    END SELECT

  END SUBROUTINE NumericDerivatives3Q_AddMat_TO_dnMat

END MODULE mod_FiniteDiff_coeff
MODULE mod_Model
!$ USE omp_lib
  USE mod_NumParameters
  USE mod_dnMat

  USE mod_EmptyModel

  USE mod_TemplateModel
  USE mod_TestModel
  USE mod_MorseModel

  USE mod_HenonHeilesModel
  USE mod_TullyModel

  USE mod_PSB3_Model
  USE mod_HONO_Model
  USE mod_HNNHp_Model
  USE mod_H2SiN_Model
  USE mod_H2NSi_Model
  USE mod_HNO3_Model

  USE mod_OneDSOC_1S1T_Model
  USE mod_OneDSOC_2S1T_Model

  USE mod_LinearHBondModel
  USE mod_BuckModel
  USE mod_PhenolModel
  USE mod_SigmoidModel
  USE mod_TwoD_Model

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Model_t,Init_Model,Eval_Pot,check_alloc_QM
  PUBLIC :: Write0_Model,Write_Model,Write_QdnV_FOR_Model
  PUBLIC :: calc_pot,calc_grad,calc_hess,calc_pot_grad,calc_pot_grad_hess
  PUBLIC :: Check_analytical_numerical_derivatives
  PUBLIC :: Eval_pot_ON_Grid,get_Q0_Model

  TYPE :: Model_t
    ! add nsurf and ndim to avoid crash when using the driver without initialization
    ! at the intialization, the variables are set-up to the correct values and are
    !   identical to QM%nsurf and QM%ndim ones respectively.
    integer :: nsurf       = 0
    integer :: ndim        = 0
    CLASS (EmptyModel_t), allocatable :: QM
  END TYPE Model_t

  !real (kind=Rkind)                     :: step = ONETENTH**4 ! model TWOD => 0.4e-7 (nderiv=2)
  real (kind=Rkind)                     :: step = ONETENTH**3 ! model TWOD => 0.6e-9 (nderiv=2)
  !real (kind=Rkind)                     :: step = ONETENTH**2 ! model TWOD => 0.1e-7 (nderiv=2)

  real (kind=Rkind)                     :: epsi = ONETENTH**10

#if defined(__QML_VER)
      character (len=Name_len) :: QML_version = __QML_VER
#else
      character (len=Name_len) :: QML_version = "unknown: -D__QML_VER=?"
#endif

#if defined(__QMLPATH)
      character (len=Line_len) :: QML_path   =                         &
       __QMLPATH
#else
      character (len=Line_len) :: QML_path   = '~/QuantumModelLib'
#endif

#if defined(__COMPILE_DATE)
      character (len=Line_len) :: compile_date = __COMPILE_DATE
#else
      character (len=Line_len) :: compile_date = "unknown: -D__COMPILE_DATE=?"
#endif

#if defined(__COMPILE_HOST)
      character (len=Line_len) :: compile_host = __COMPILE_HOST
#else
      character (len=Line_len) :: compile_host = "unknown: -D__COMPILE_HOST=?"
#endif
#if defined(__COMPILER)
      character (len=Line_len) :: compiler = __COMPILER
#else
      character (len=Line_len) :: compiler = "unknown: -D__COMPILER=?"
#endif
#if defined(__COMPILER_VER)
      character (len=Line_len) :: compiler_ver = __COMPILER_VER
#else
      character (len=Line_len) :: compiler_ver = "unknown: -D__COMPILER_VER=?"
#endif
#if defined(__COMPILER_OPT)
      character (len=Line_len) :: compiler_opt = &
      __COMPILER_OPT
#else
      character (len=Line_len) :: compiler_opt = "unknown: -D__COMPILER_OPT=?"
#endif
#if defined(__COMPILER_LIBS)
      character (len=Line_len) :: compiler_libs = __COMPILER_LIBS
#else
      character (len=Line_len) :: compiler_libs = "unknown: -D__COMPILER_LIBS=?"
#endif




  TYPE(Model_t), PUBLIC  :: QuantumModel

CONTAINS


  SUBROUTINE Read_Model(QModel_inout,nio,read_nml1)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (EmptyModel_t), intent(inout) :: QModel_inout ! variable to transfer info to the init
    integer,             intent(in)    :: nio
    logical,             intent(inout) :: read_nml1

    ! local variable
    integer :: ndim,nsurf,nderiv,option
    logical :: adiabatic,numeric,PubliUnit,read_nml
    character (len=20) :: pot_name
    integer :: err_read

    ! Namelists for input file
    namelist /potential/ ndim,nsurf,pot_name,numeric,adiabatic,option,PubliUnit,read_nml

!    ! Default values defined
    ndim      = QModel_inout%ndim
    nsurf     = QModel_inout%nsurf
    adiabatic = QModel_inout%adiabatic
    pot_name  = 'morse'
    numeric   = .FALSE.
    PubliUnit = .FALSE.
    read_nml  = .TRUE. ! if T, read the namelist in PotLib (HenonHeiles ....)
    option    = -1 ! no option


    write(out_unitp,*) 'Reading input file . . .'
    read(nio,nml=potential,IOSTAT=err_read)

    IF (err_read < 0) THEN
      write(out_unitp,*) ' ERROR in Read_Model'
      write(out_unitp,*) ' End-of-file or End-of-record'
      write(out_unitp,*) ' The namelist "potential" is probably absent'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,*)
      STOP ' ERROR in Read_Model'
    ELSE IF (err_read > 0) THEN
      write(out_unitp,*) ' ERROR in Read_Model'
      write(out_unitp,*) ' Some parameter names of the namelist "potential" are probaly wrong'
      write(out_unitp,*) ' check your data!'
      write(out_unitp,nml=potential)
      STOP ' ERROR in Read_Model'
    END IF

    !write(out_unitp,nml=potential)

    read_nml1              = read_nml

    QModel_inout%option    = option
    QModel_inout%ndim      = ndim
    QModel_inout%nsurf     = nsurf
    QModel_inout%adiabatic = adiabatic
    QModel_inout%numeric   = numeric
    QModel_inout%pot_name  = trim(pot_name)
    !QModel_inout%pot_name  = strdup(pot_name) ! panic with nagfor !!!
    QModel_inout%PubliUnit = PubliUnit

  END SUBROUTINE Read_Model

  SUBROUTINE Init_Model(QModel,pot_name,ndim,nsurf,adiabatic,       &
                        read_param,param_file_name,nio_param_file,      &
                        option,PubliUnit,Print_init)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (Model_t),     intent(inout)         :: QModel

    character (len=*),   intent(in), optional :: pot_name
    integer,             intent(in), optional :: ndim,nsurf
    logical,             intent(in), optional :: adiabatic

    logical,             intent(in), optional :: read_param
    integer,             intent(in), optional :: nio_param_file
    character (len=*),   intent(in), optional :: param_file_name

    integer,             intent(in), optional :: option
    logical,             intent(in), optional :: PubliUnit
    logical,             intent(in), optional :: Print_init

    ! local variables
    TYPE(EmptyModel_t)             :: QModel_in ! variable to transfer info to the init
    integer                        :: i,nio_loc
    logical                        :: read_param_loc,read_nml,Print_init_loc
    character (len=:), allocatable :: param_file_name_loc,pot_name_loc
    real (kind=Rkind), allocatable :: Q0(:)

    Print_init_loc = .TRUE.
    IF (present(Print_init)) Print_init_loc = Print_init

    IF (Print_init_loc) THEN
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== QML: Quantum Model Lib (E-CAM) ==============='
      write(out_unitp,*) '== QML version:       ',QML_version
      write(out_unitp,*) '== QML path:          ',trim(adjustl(QML_path))
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) '== Compiled on       "',trim(compile_host), '" the ',trim(compile_date)
      write(out_unitp,*) '== Compiler:         ',trim(compiler)
      write(out_unitp,*) '== Compiler version: ',trim(compiler_ver)
      write(out_unitp,*) '== Compiler options: ',trim(compiler_opt)
      write(out_unitp,*) '== Compiler libs:     ',trim(compiler_libs)
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) 'QML is under GNU LGPL3 license and '
      write(out_unitp,*) '  is written by David Lauvergnat [1]'
      write(out_unitp,*) '  with contributions of'
      write(out_unitp,*) '     Félix MOUHAT [2]'
      write(out_unitp,*) '     Liang LIANG [3]'
      write(out_unitp,*) '     Emanuele MARSILI [1,4]'
      write(out_unitp,*)
      write(out_unitp,*) '[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
      write(out_unitp,*) '[2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France'
      write(out_unitp,*) '[3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France'
      write(out_unitp,*) '[4]: Durham University, Durham, UK'
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== Initialization of the Model =================='
    END IF

    IF (allocated(QModel%QM)) deallocate(QModel%QM)

    ! set the "File_path" in the Lib_module.f90
    File_path = trim(adjustl(QML_path))

    IF (present(ndim)) THEN
      QModel_in%ndim      = ndim
    ELSE
      QModel_in%ndim      = 0
    END IF
    IF (present(nsurf)) THEN
      QModel_in%nsurf     = nsurf
    ELSE
      QModel_in%nsurf     = 0
    END IF

    IF (present(adiabatic)) THEN
      QModel_in%adiabatic = adiabatic
    ELSE
      QModel_in%adiabatic = .TRUE.
    END IF

    IF (present(option)) THEN
      QModel_in%option = option
    ELSE
      QModel_in%option = -1
    END IF

    IF (present(PubliUnit)) THEN
      QModel_in%PubliUnit      = PubliUnit
    ELSE
      QModel_in%PubliUnit      = .FALSE.
    END IF


    IF (present(read_param)) THEN
      read_param_loc = read_param
    ELSE
      read_param_loc = .FALSE.
    END IF

    IF (present(nio_param_file)) THEN
      nio_loc = nio_param_file
    ELSE
      IF (present(param_file_name)) THEN
        IF (len_trim(param_file_name) == 0) THEN
          param_file_name_loc = strdup("input.dat")
          nio_loc = 99
        ELSE
          param_file_name_loc = strdup(param_file_name)
          nio_loc = 99
        END IF
      ELSE
      nio_loc = in_unitp
      END IF
    END IF

    IF (present(pot_name)) THEN
      pot_name_loc  = strdup(pot_name)
      CALL string_uppercase_TO_lowercase(pot_name_loc)
    ELSE
      IF (.NOT. read_param_loc) THEN
        write(out_unitp,*) 'ERROR in Init_Model'
        write(out_unitp,*) ' pot_name is not present and read_param=F'
        STOP ' pot_name is not present and read_param=F'
      END IF
    END IF

    read_nml       = .FALSE.
    read_param_loc = (read_param_loc .OR.  pot_name_loc == 'read_model')
    IF (read_param_loc) THEN
      IF (nio_loc /= in_unitp) THEN
        open(unit=nio_loc,file=param_file_name_loc,status='old',form='formatted')
      END IF
      CALL Read_Model(QModel_in,nio_loc,read_nml)
      pot_name_loc = QModel_in%pot_name

      IF (allocated(param_file_name_loc)) deallocate(param_file_name_loc)
    END IF


    IF (QModel_in%adiabatic) THEN
      IF (Print_init_loc) write(out_unitp,*) 'Adiabatic potential . . .'
    ELSE
      IF (Print_init_loc) write(out_unitp,*) 'Non-adiabatic potential . . .'
    END IF

    IF (QModel_in%numeric .AND. Print_init_loc) THEN
      write(out_unitp,*) 'You have decided to perform a numeric checking of the analytic formulas.'
    END IF

    read_param_loc = (read_param_loc .AND. read_nml) ! this enables to not read the next namelist when read_param_loc=t

    CALL string_uppercase_TO_lowercase(pot_name_loc)
    IF (Print_init_loc) write(out_unitp,*) 'pot_name_loc: ',pot_name_loc

    SELECT CASE (pot_name_loc)
    CASE ('morse')
      !! === README ==
      !! Morse potential: V(R) = D*(1-exp(-a*(r-Req))**2
      !! pot_name  = 'Morse'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 1744.60504565084306291455 au
      !! remark: Default parameters for H-F
      !! === END README ==

      allocate(MorseModel_t :: QModel%QM)
      QModel%QM = Init_MorseModel(QModel_in,read_param=read_param_loc,nio_param_file=nio_loc)

    CASE ('sigmoid')
      !! sigmoid function: A * 1/2(1+e*tanh((x-B)/C))  remark: e=+/-1
      allocate(SigmoidModel_t :: QModel%QM)
      QModel%QM = Init_SigmoidModel(QModel_in,read_param=read_param_loc,nio_param_file=nio_loc)

    CASE ('buck')
      !! === README ==
      !! Buckingham potential: V(R) = A*exp(-B*R)-C/R^6
      !! pot_name  = 'buck'
      !! ndim      = 1
      !! nsurf     = 1
      !! reduced mass      = 36423.484024390622 au
      !! remark: default parameters for Ar2
      !! ref:  R.A. Buckingham, Proc. R. Soc. A Math. Phys. Eng. Sci. 168 (1938) 264–283. doi:10.1098/rspa.1938.0173
      !! === END README ==
      allocate(BuckModel_t :: QModel%QM)
      QModel%QM = Init_BuckModel(QModel_in,read_param=read_param_loc,nio_param_file=nio_loc)

    CASE ('hbond','linearhbond ')
      !! === README ==
      !! LinearHBond potential: Morse1(QQ/2+q,param1)+Morse2(QQ/2-q,param2)+Eref2+Buckingham(QQ)
      !! pot_name  = 'HBond'
      !! ndim      = 2   (QQ,q)
      !! nsurf     = 1
      !! reduced masses      = (/ 29156.946380706224, 1837.1526464003414 /) au
      !! remark:
      !!    A--------------H-----X--------------------B
      !!     <--------------QQ----------------------->
      !!                    <-q->
      !! ref:  Eq 3.79 of J. Beutier, thesis.
      !! === END README ==
      allocate(LinearHBondModel_t :: QModel%QM)
      QModel%QM = Init_LinearHBondModel(QModel_in,read_param=read_param_loc,   &
                                 nio_param_file=nio_loc)

    CASE ('henonheiles')
      !! === README ==
      !! HenonHeiles potential
      !! pot_name  = 'HenonHeiles'
      !! ndim      > 1
      !! nsurf     = 1
      !! reduced masses(:)      = ONE au
      !! ref:  parameters taken from M. Nest, H.-D. Meyer, J. Chem. Phys. 117 (2002) 10499. doi:10.1063/1.1521129
      !! === END README ==
      allocate(HenonHeilesModel_t :: QModel%QM)
      QModel%QM = Init_HenonHeilesModel(QModel_in,                      &
                        read_param=read_param_loc,nio_param_file=nio_loc)

    CASE ('tully')
      !! === README ==
      !! Tully potential: three options
      !! pot_name  = 'Tully'
      !! ndim      = 1
      !! nsurf     = 2
      !! reduced mass      = 2000. au
      !! remark: three options are possible (option = 1,2,3)
      !! ref:  Tully, J. Chem. Phys. V93, pp15, 1990
      !! === END README ==
      allocate(TullyModel_t :: QModel%QM)
      QModel%QM = Init_TullyModel(QModel_in,read_param=read_param_loc,  &
                                  nio_param_file=nio_loc)

    CASE ('1dsoc','1dsoc_1s1t')
      !! === README ==
      !! Spin Orbit coupling model
      !! pot_name  = '1DSOC_1S1T'
      !! ndim      = 1
      !! nsurf     = 4 or 2
      !! reduced mass      = 20000. au
      !! remarks: 1 singlet and 3 triplet components                           => nsurf     = 4
      !!  or      1 singlet and 1 linear combibation of the triplet components => nsurf     = 2
      !! ref: Giovanni Granucci, Maurizio Persico, and Gloria Spighi, J. Chem. Phys. V137, p22A501 (2012)
      !! === END README ==
      allocate(OneDSOC_1S1T_Model_t :: QModel%QM)
      QModel%QM = Init_OneDSOC_1S1T_Model(QModel_in,                    &
                       read_param=read_param_loc,nio_param_file=nio_loc)

    CASE ('1dsoc_2s1t')
      !! === README ==
      !! Spin Orbit coupling model
      !! pot_name  = '1DSOC_2S1T'
      !! ndim      = 1
      !! nsurf     = 4
      !! reduced mass      = 20000. au
      !! remark: 2 singlets and 1 triplet (2 linear combinations of the triplet components are not included) => nsurf     = 4
      !! ref: Giovanni Granucci, Maurizio Persico, and Gloria Spighi, J. Chem. Phys. V137, p22A501 (2012)
      !! === END README ==
      allocate(OneDSOC_2S1T_Model_t :: QModel%QM)
      QModel%QM = Init_OneDSOC_2S1T_Model(QModel_in,                    &
                       read_param=read_param_loc,nio_param_file=nio_loc)

    CASE ('phenol')
      !! === README ==
      !! Phenol model
      !! pot_name  = 'Phenol'
      !! ndim      = 2 (R=rOH, th=OH-torsion)
      !! nsurf     = 3
      !! Diagonal Metric Tensor(:)      = (/ 0.0005786177, 0.0002550307 /) au
      !! remark:
      !! ref: Z. Lan, W. Domcke, V. Vallet, A.L. Sobolewski, S. Mahapatra, J. Chem. Phys. 122 (2005) 224315. doi:10.1063/1.1906218.
      !! === END README ==

      allocate(PhenolModel_t :: QModel%QM)
      QModel%QM = Init_PhenolModel(QModel_in,read_param=read_param_loc,nio_param_file=nio_loc)


    CASE ('twod')
      !! === README ==
      !! 2D model
      !! pot_name  = 'TwoD'
      !! ndim      = 2 (X,Y)
      !! nsurf     = 2
      !! Reduced masses(:)      = (/ 20000., 6667. /) au
      !! remark: The parameter values have been modified
      !! ref: A. Ferretti, G. Granucci, A. Lami, M. Persico, G. Villani, J. Chem. Phys. 104, 5517 (1996); https://doi.org/10.1063/1.471791
      !! === END README ==
      allocate(TwoD_Model_t :: QModel%QM)
      QModel%QM = Init_TwoD_Model(QModel_in,read_param=read_param_loc,  &
                                  nio_param_file=nio_loc)

    CASE ('psb3')
      !! === README ==
      !! Model for the photo-isomerization of the penta-2,4-dieniminium (PSB3) cation.
      !! pot_name  = 'PSB3'
      !! ndim      = 3
      !! nsurf     = 2
      !! remarks: two options are possible (option = 1,2)
      !! The default is option=1 (unpublished).
      !! The parameters for option=2 come from the following reference.
      !! ref: E. Marsili, M. H. Farag, X. Yang, L. De Vico, and M. Olivucci, JPCA, 123, 1710–1719 (2019). https://doi.org/10.1021/acs.jpca.8b10010
      !! === END README ==
      allocate(PSB3_Model_t :: QModel%QM)
      QModel%QM = Init_PSB3_Model(QModel_in,read_param=read_param_loc,  &
                                  nio_param_file=nio_loc)

    CASE ('hono')
      allocate(HONO_Model_t :: QModel%QM)
      QModel%QM = Init_HONO_Model(QModel_in,read_param=read_param_loc,  &
                                  nio_param_file=nio_loc)

    CASE ('hno3')
      allocate(HNO3_Model_t :: QModel%QM)
      QModel%QM = Init_HNO3_Model(QModel_in,read_param=read_param_loc,  &
                                  nio_param_file=nio_loc)

    CASE ('hnnhp')
      allocate(HNNHp_Model_t :: QModel%QM)
      QModel%QM = Init_HNNHp_Model(QModel_in,read_param=read_param_loc, &
                                   nio_param_file=nio_loc)

    CASE ('h2sin')
      allocate(H2SiN_Model_t :: QModel%QM)
      QModel%QM = Init_H2SiN_Model(QModel_in,read_param=read_param_loc, &
                                   nio_param_file=nio_loc)

    CASE ('h2nsi')
      allocate(H2NSi_Model_t :: QModel%QM)
      QModel%QM = Init_H2NSi_Model(QModel_in,read_param=read_param_loc, &
                                   nio_param_file=nio_loc)

    CASE ('template')
      !! 3D-potential with 1 surface
      allocate(TemplateModel_t :: QModel%QM)
      QModel%QM = Init_TemplateModel(QModel_in,                         &
                        read_param=read_param_loc,nio_param_file=nio_loc)

    CASE ('test')
      !! test-potential
      allocate(TemplateModel_t :: QModel%QM)
      QModel%QM = Init_testModel(QModel_in,                             &
                        read_param=read_param_loc,nio_param_file=nio_loc)

    CASE DEFAULT
        write(out_unitp,*) ' ERROR in Init_Model'
        write(out_unitp,*) ' This model/potential is unknown. pot_name: ',pot_name_loc
        STOP 'STOP in Init_Model: Other potentials have to be done'
    END SELECT

    IF (read_param_loc .AND. nio_loc /= in_unitp) THEN
       close(nio_loc)
    END IF

    IF (present(ndim)) THEN
    IF (ndim /= QModel%QM%ndim) THEN
        write(out_unitp,*) ' ERROR in Init_Model'
        write(out_unitp,*) ' ndim is present and ...'
        write(out_unitp,*) ' its value is not equal to QModel%QM%ndim'
        write(out_unitp,*) ' ndim,QModel%QM%ndim',ndim,QModel%QM%ndim
        write(out_unitp,*) ' check your data!'
        STOP 'STOP in Init_Model: wrong ndim'
    END IF
    END IF
    IF (present(nsurf)) THEN
    IF (nsurf /= QModel%QM%nsurf) THEN
        write(out_unitp,*) ' ERROR in Init_Model'
        write(out_unitp,*) ' nsurf is present and ...'
        write(out_unitp,*) ' its value is not equal to QModel%QM%nsurf'
        write(out_unitp,*) ' nsurf,QModel%QM%nsurf',nsurf,QModel%QM%nsurf
        write(out_unitp,*) ' check your data!'
        STOP 'STOP in Init_Model: wrong nsurf'
    END IF
    END IF

    QModel%ndim  = QModel%QM%ndim
    QModel%nsurf = QModel%QM%nsurf
    IF (Print_init_loc) THEN
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Quantum Model'
      CALL Write_Model(QModel,nio=out_unitp)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '================================================='
    END IF


  END SUBROUTINE Init_Model

  SUBROUTINE get_Q0_Model(Q0,QModel,option)
  USE mod_Lib
  IMPLICIT NONE

    real (kind=Rkind),  intent(inout)            :: Q0(:)
    TYPE (Model_t),    intent(in)               :: QModel
    integer,            intent(in)               :: option


    integer :: err_Q0
!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='get_Q0_Model'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      flush(out_unitp)
    END IF

    CALL check_alloc_QM(QModel,name_sub)

    IF (size(Q0) /= QModel%QM%ndim) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The size of Q0 is not QModel%ndim: '
      write(out_unitp,*) ' size(Q0)',size(Q0)
      write(out_unitp,*) ' ndim',QModel%QM%ndim
      STOP 'STOP in get_Q0_Model: Wrong Q0 size'
    END IF

    Q0(:) = ZERO

    CALL get_Q0_QModel(QModel%QM,Q0,err_Q0)
    IF (err_Q0 /= 0) THEN
      CALL Write_Model(QModel,out_unitp)
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' Q0 is not set-up in the model'
      STOP 'STOP Q0 is not set-up in the model'
    END IF

    IF (debug) THEN
      CALL Write_RVec(Q0,out_unitp,nbcol1=5,name_info='Q0: ')
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE get_Q0_Model

  ! check if the QM [CLASS(EmptyModel_t)] is allocated
  SUBROUTINE check_alloc_QM(QModel,name_sub_in)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (Model_t),    intent(in)     :: QModel
    character (len=*),  intent(in)     :: name_sub_in

    IF ( .NOT. allocated(QModel%QM)) THEN
      write(out_unitp,*) ' ERROR in check_alloc_QM'
      write(out_unitp,*) ' QM is not allocated in QModel.'
      write(out_unitp,*) '  check_alloc_QM is called from ',name_sub_in
      write(out_unitp,*) '  You MUST initialize the model with:'
      write(out_unitp,*) '    CALL init_Model(...) in mod_Model.f90'
      write(out_unitp,*) ' or'
      write(out_unitp,*) '    CALL sub_Init_Qmodel(...) in Model_driver.f90'
      STOP 'STOP in check_alloc_QM: QM is not allocated in QModel.'
    END IF

  END SUBROUTINE check_alloc_QM


  SUBROUTINE Eval_Pot(QModel,Q,PotVal,nderiv,NAC,Vec,numeric)
  USE mod_dnS
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)            :: QModel
    TYPE (dnMat_t),     intent(inout)            :: PotVal
    real (kind=Rkind),  intent(in)               :: Q(:)
    integer,            intent(in),    optional  :: nderiv
    TYPE (dnMat_t),     intent(inout), optional  :: NAC,Vec
    logical,            intent(in),    optional  :: numeric

    ! local variables
    integer                    :: nderiv_loc
    TYPE (dnMat_t)             :: Vec_loc,NAC_loc
    logical                    :: numeric_loc,adia_loc

    integer :: numeric_option = 3   ! 0 old (up to 2d derivatives
                                    ! 3 version up to 3d derivatives less points than 4
                                    ! 4 version up to 3d derivatives more points than 3

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
      flush(out_unitp)
    END IF

    CALL check_alloc_QM(QModel,name_sub)

    IF (debug) THEN
      write(out_unitp,*) '  QModel%QM%numeric   ',QModel%QM%numeric
      write(out_unitp,*) '  QModel%QM%adiabatic ',QModel%QM%adiabatic
      flush(out_unitp)
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(3,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF

    IF (present(numeric)) THEN
      numeric_loc = (numeric .AND. nderiv_loc > 0)
    ELSE
      numeric_loc = (QModel%QM%numeric .AND. nderiv_loc > 0)
    END IF

    adia_loc = (QModel%QM%adiabatic .AND. QModel%QM%nsurf > 1)

    IF (numeric_loc) THEN  ! numerical
      IF (.NOT. adia_loc) THEN
         SELECT CASE (numeric_option)
         CASE (0)
           CALL Eval_Pot_Numeric_dia_old(QModel,Q,PotVal,nderiv_loc)
         CASE (3)
           CALL Eval_Pot_Numeric_dia_v3(QModel,Q,PotVal,nderiv_loc)
         CASE (4)
           CALL Eval_Pot_Numeric_dia_v4(QModel,Q,PotVal,nderiv_loc)
         CASE Default
           CALL Eval_Pot_Numeric_dia_old(QModel,Q,PotVal,nderiv_loc)
         END SELECT
      ELSE

        IF (present(Vec)) THEN
          IF (present(NAC)) THEN
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,      &
                                                 Vec,NAC,numeric_option)
          ELSE
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,      &
                                             Vec,NAC_loc,numeric_option)
            CALL QML_dealloc_dnMat(NAC_loc)
          END IF
        ELSE
          IF (present(NAC)) THEN
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,      &
                                             Vec_loc,NAC,numeric_option)
          ELSE
            CALL Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv_loc,      &
                                         Vec_loc,NAC_loc,numeric_option)
            CALL QML_dealloc_dnMat(NAC_loc)
          END IF
          CALL QML_dealloc_dnMat(Vec_loc)
        END IF
      END IF
    ELSE ! analytical calculation
      IF (present(Vec)) THEN
        IF (present(NAC)) THEN
          CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv_loc,Vec=Vec,Nac=NAC)
        ELSE
          CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv_loc,Vec=Vec)
        END IF
      ELSE
        IF (present(NAC)) THEN
          CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv_loc,Nac=NAC)
        ELSE
          CALL Eval_Pot_ana(QModel,Q,PotVal,nderiv_loc)
        END IF
      END IF

    END IF


    IF (debug) THEN
      IF ( QModel%QM%adiabatic) THEN
        write(out_unitp,*) 'PotVal (adia)'
      ELSE
        write(out_unitp,*) 'PotVal (dia)'
      END IF
      CALL QML_Write_dnMat(PotVal,nio=out_unitp)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Eval_Pot

  SUBROUTINE Eval_Pot_ana(QModel,Q,PotVal,nderiv,NAC,Vec)
  USE mod_dnS
  IMPLICIT NONE

    TYPE (Model_t),       intent(inout)            :: QModel

    TYPE (dnMat_t),        intent(inout)            :: PotVal
    real (kind=Rkind),     intent(in)               :: Q(:)
    integer,               intent(in)               :: nderiv
    TYPE (dnMat_t),        intent(inout), optional  :: NAC,Vec

    ! local variables
    integer                     :: i,j,id
    TYPE (dnMat_t)              :: PotVal_dia,Vec_loc,NAC_loc
    TYPE (dnS_t), allocatable   :: dnQ(:)
    TYPE (dnS_t), allocatable   :: Mat_OF_PotDia(:,:)
    logical :: old = .FALSE.
    !logical :: old = .TRUE.

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot_ana'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      write(out_unitp,*) '   nderiv    ',nderiv
      flush(out_unitp)
    END IF

    CALL check_alloc_QM(QModel,name_sub)

    IF (debug) write(out_unitp,*) '   adiabatic ',QModel%QM%adiabatic


    IF ( QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    ! allocate Mat_OF_PotDia
    allocate(Mat_OF_PotDia(QModel%QM%nsurf,QModel%QM%nsurf))
    DO j=1,size(Mat_OF_PotDia(1,:))
    DO i=1,size(Mat_OF_PotDia(:,1))
        CALL QML_alloc_dnS(Mat_OF_PotDia(i,j),QModel%QM%ndim,nderiv)
    END DO
    END DO

    ! intialization of the dnQ(:)
    allocate(dnQ(QModel%QM%ndim))
    DO i=1,QModel%QM%ndim
      dnQ(i) = QML_init_dnS(Q(i),ndim=QModel%QM%ndim,nderiv=nderiv,iQ=i) ! to set up the derivatives
    END DO

    CALL QModel%QM%Eval_QModel_Pot(Mat_OF_PotDia,dnQ,nderiv=nderiv)


    PotVal = Mat_OF_PotDia ! transfert the potential and its derivatives to the matrix form (PotVal)

    ! deallocation
    DO i=1,size(dnQ)
      CALL QML_dealloc_dnS(dnQ(i))
    END DO
    deallocate(dnQ)

    DO j=1,size(Mat_OF_PotDia(1,:))
    DO i=1,size(Mat_OF_PotDia(:,1))
      CALL QML_dealloc_dnS(Mat_OF_PotDia(i,j))
    END DO
    END DO
    deallocate(Mat_OF_PotDia)
    ! end deallocation

    IF ( QModel%QM%adiabatic .AND. QModel%QM%nsurf > 1) THEN
      IF (debug) THEN
        write(out_unitp,*) 'PotVal (dia)'
        CALL QML_Write_dnMat(PotVal,nio=out_unitp)
        flush(out_unitp)
      END IF
      IF (.NOT. allocated(QModel%QM%Vec0)) allocate(QModel%QM%Vec0)

      PotVal_dia = PotVal
IF (old) THEN
      IF (present(Vec)) THEN
        IF (present(NAC)) THEN
          CALL  dia_TO_adia_old(PotVal_dia,PotVal,Vec,QModel%QM%Vec0,NAC,nderiv)
        ELSE
          CALL  dia_TO_adia_old(PotVal_dia,PotVal,Vec,QModel%QM%Vec0,NAC_loc,nderiv)
          CALL QML_dealloc_dnMat(NAC_loc)
        END IF
      ELSE
        IF (present(NAC)) THEN
          CALL  dia_TO_adia_old(PotVal_dia,PotVal,Vec_loc,QModel%QM%Vec0,NAC,nderiv)
        ELSE
          CALL  dia_TO_adia_old(PotVal_dia,PotVal,Vec_loc,QModel%QM%Vec0,NAC_loc,nderiv)
          CALL QML_dealloc_dnMat(NAC_loc)
        END IF
        CALL QML_dealloc_dnMat(Vec_loc)
      END IF
      CALL QML_dealloc_dnMat(PotVal_dia)
ELSE
      IF (present(Vec)) THEN
        IF (present(NAC)) THEN
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec,QModel%QM%Vec0,NAC,nderiv)
        ELSE
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec,QModel%QM%Vec0,NAC_loc,nderiv)
          CALL QML_dealloc_dnMat(NAC_loc)
        END IF
      ELSE
        IF (present(NAC)) THEN
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec_loc,QModel%QM%Vec0,NAC,nderiv)
        ELSE
          CALL dia_TO_adia(PotVal_dia,PotVal,Vec_loc,QModel%QM%Vec0,NAC_loc,nderiv)
          CALL QML_dealloc_dnMat(NAC_loc)
        END IF
        CALL QML_dealloc_dnMat(Vec_loc)
      END IF
      CALL QML_dealloc_dnMat(PotVal_dia)
END IF
!STOP 'STOP in Eval_Pot_ana'
    END IF


    IF (debug) THEN
      IF ( QModel%QM%adiabatic) THEN
        write(out_unitp,*) 'PotVal (adia)'
      ELSE
        write(out_unitp,*) 'PotVal (dia)'
      END IF
      CALL QML_Write_dnMat(PotVal,nio=out_unitp)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF


  END SUBROUTINE Eval_Pot_ana

  SUBROUTINE Set_QplusDQ(Qout,Qin,indQ,indDQ,step_sub)
  IMPLICIT NONE

    real (kind=Rkind), intent(inout)  :: Qout(:)
    real (kind=Rkind), intent(in)     :: Qin(:)
    integer,           intent(in)     :: indQ(:)
    integer,           intent(in)     :: indDQ(:)
    real (kind=Rkind), intent(in)     :: step_sub

    integer :: i

    IF (size(Qin) /= size(Qout) .OR. size(indQ) /= size(indDQ) .OR.     &
        minval(indQ) < 1 .OR. maxval(indQ) > size(Qin)) THEN
      write(out_unitp,*) ' ERROR in Set_QplusDQ'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' size(Qin),size(Qout)  ',size(Qin),size(Qout)
      write(out_unitp,*) ' size(indQ),size(indDQ)',size(indQ),size(indDQ)
      write(out_unitp,*) ' indQ(:)',indQ
      write(out_unitp,*) ' indDQ(:)',indDQ
      STOP 'STOP in Set_QplusDQ: Inconsitent parameters.'
    END IF

    Qout(:) = Qin
    DO i=1,size(indQ)
      Qout(indQ(i)) = Qin(indQ(i)) + step_sub * real(indDQ(i),kind=Rkind)
    END DO

  END SUBROUTINE Set_QplusDQ


  SUBROUTINE Eval_Pot_Numeric_dia_v4(QModel,Q,PotVal,nderiv)
  USE mod_FiniteDiff_coeff
  IMPLICIT NONE

    TYPE (Model_t),   intent(inout)   :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j,k,ip,jp,kp
    integer                            :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_dia_v4')


    IF (QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL QML_alloc_dnMat(PotVal_loc0,QModel%QM%nsurf,QModel%QM%ndim,nderiv=0)
    PotVal = ZERO

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0)
    CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=4)

    IF (nderiv >= 1) THEN ! along ONE coordinates (first derivatives and higher)

      ! Numeric evaluation of forces
      DO i=1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                            indQ=[i],indDQ=ind1DQ,option=4)
        END DO

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=1,QModel%QM%ndim
        IF (i == j) CYCLE

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                           indQ=[i,j],indDQ=ind2DQ,option=4)
        END DO

      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN ! 3d derivatives:d3/dQidQjdQk

      DO i=1,QModel%QM%ndim
      DO j=1,QModel%QM%ndim
      IF (i == j) CYCLE
      DO k=1,QModel%QM%ndim
        IF (i == k .OR. j == k) CYCLE

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                         indQ=[i,j,k],indDQ=ind3DQ,option=4)
        END DO

      END DO
      END DO
      END DO
    END IF

    IF (nderiv >= 1) PotVal%d1 = PotVal%d1/step
    IF (nderiv >= 2) PotVal%d2 = PotVal%d2/step**2
    IF (nderiv >= 3) PotVal%d3 = PotVal%d3/step**3

    deallocate(Q_loc)
    CALL QML_dealloc_dnMat(PotVal_loc0)

  END SUBROUTINE Eval_Pot_Numeric_dia_v4
  SUBROUTINE Eval_Pot_Numeric_dia_v3(QModel,Q,PotVal,nderiv)
  USE mod_FiniteDiff_coeff
  IMPLICIT NONE

    TYPE (Model_t),   intent(inout)   :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j,k,ip,jp,kp

    integer                            :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)



    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_dia_v3')


    IF (QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL QML_alloc_dnMat(PotVal_loc0,QModel%QM%nsurf,QModel%QM%ndim,nderiv=0)
    PotVal = ZERO

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0)
    PotVal%d0 = PotVal_loc0%d0

    CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=3)


    IF (nderiv >= 1) THEN ! 1st derivatives

      DO i=1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                           indQ=[i],indDQ=ind1DQ,option=3)
        END DO

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                          indQ=[i,j],indDQ=ind2DQ,option=3)
        END DO

        PotVal%d2(:,:,i,j) = PotVal%d2(:,:,j,i)

        IF (nderiv >= 3) THEN
          PotVal%d3(:,:,i,j,i) = PotVal%d3(:,:,i,i,j)
          PotVal%d3(:,:,j,i,i) = PotVal%d3(:,:,i,i,j)

          PotVal%d3(:,:,i,j,j) = PotVal%d3(:,:,j,j,i)
          PotVal%d3(:,:,j,i,j) = PotVal%d3(:,:,j,j,i)
        END IF

      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN ! 3d derivatives: d3/dQidQidQj

      ! d3/dQidQjdQk
      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim
      DO k=j+1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                       indQ=[i,j,k],indDQ=ind3DQ,option=3)
        END DO

        PotVal%d3(:,:,k,i,j) = PotVal%d3(:,:,k,j,i)
        PotVal%d3(:,:,j,k,i) = PotVal%d3(:,:,k,j,i)
        PotVal%d3(:,:,i,j,k) = PotVal%d3(:,:,k,j,i)

        PotVal%d3(:,:,i,k,j) = PotVal%d3(:,:,k,j,i)
        PotVal%d3(:,:,j,i,k) = PotVal%d3(:,:,k,j,i)


      END DO
      END DO
      END DO
    END IF

    IF (nderiv >= 1) PotVal%d1 = PotVal%d1/step
    IF (nderiv >= 2) PotVal%d2 = PotVal%d2/step**2
    IF (nderiv >= 3) PotVal%d3 = PotVal%d3/step**3

    deallocate(Q_loc)
    CALL QML_dealloc_dnMat(PotVal_loc0)

  END SUBROUTINE Eval_Pot_Numeric_dia_v3

  SUBROUTINE Eval_Pot_Numeric_dia_old(QModel,Q,PotVal,nderiv)
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0
    integer                            :: i,j

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_dia_old')


    IF (QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    PotVal = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL QML_alloc_dnMat(PotVal_loc0,QModel%QM%nsurf,QModel%QM%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0)
    PotVal%d0 = PotVal_loc0%d0


    IF (nderiv >= 1) THEN ! 1st derivatives

      ! Numeric evaluation of forces
      DO i=1,QModel%QM%ndim

        Q_loc(i) = Q(i) + step        ! q+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0) ! Ep(q+dq)
        PotVal%d1(:,:,i) = PotVal_loc0%d0

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = PotVal_loc0%d0
        END IF

        Q_loc(i) = Q(i) - step        ! q-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0) ! Ep(q-dq)
        PotVal%d1(:,:,i) = (PotVal%d1(:,:,i)-PotVal_loc0%d0)/(TWO*step)


        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = (PotVal%d2(:,:,i,i) + PotVal_loc0%d0 - TWO*PotVal%d0)/ &
                                step**2
        END IF

        Q_loc(i) = Q(i)

      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) + PotVal_loc0%d0

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0

        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i)/(FOUR*step**2)
        PotVal%d2(:,:,i,j) = PotVal%d2(:,:,j,i)

        Q_loc(i) = Q(i)
        Q_loc(j) = Q(j)
      END DO
      END DO
    END IF

    deallocate(Q_loc)
    CALL QML_dealloc_dnMat(PotVal_loc0)

  END SUBROUTINE Eval_Pot_Numeric_dia_old

  SUBROUTINE Eval_Pot_Numeric_adia(QModel,Q,PotVal,nderiv,Vec,NAC,option)
  IMPLICIT NONE

    TYPE (Model_t),   intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv
    TYPE (dnMat_t),    intent(inout)  :: Vec,NAC
    integer,           intent(in)     :: option


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Eval_Pot_Numeric_adia'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      write(out_unitp,*) '   nderiv',nderiv
      write(out_unitp,*) '   option',option
      flush(out_unitp)
    END IF

    SELECT CASE (option)
    CASE (0)
      CALL Eval_Pot_Numeric_adia_old(QModel,Q,PotVal,nderiv,Vec,NAC)
    CASE (3)
      CALL Eval_Pot_Numeric_adia_v3(QModel,Q,PotVal,nderiv,Vec,NAC)
    CASE (4)
    !  CALL Eval_Pot_Numeric_adia_v4(QModel,Q,PotVal,nderiv,Vec,NAC)
    CASE Default
      CALL Eval_Pot_Numeric_adia_old(QModel,Q,PotVal,nderiv,Vec,NAC)
    END SELECT

    IF (debug) THEN
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF
  END SUBROUTINE Eval_Pot_Numeric_adia
  SUBROUTINE Eval_Pot_Numeric_adia_old(QModel,Q,PotVal,nderiv,Vec,NAC)
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv
    TYPE (dnMat_t),    intent(inout)  :: Vec,NAC

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0,Vec_loc0
    integer                            :: i,j
    real (kind=Rkind), allocatable     :: tVec(:,:)

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_adia_old')


    IF (QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                           nderiv=nderiv)
    END IF
    PotVal = ZERO

    IF (QML_Check_NotAlloc_dnMat(Vec,nderiv) ) THEN
      CALL QML_alloc_dnMat(Vec,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    Vec = ZERO

    IF (QML_Check_NotAlloc_dnMat(NAC,nderiv) ) THEN
      CALL QML_alloc_dnMat(NAC,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    NAC = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL QML_alloc_dnMat(PotVal_loc0,QModel%QM%nsurf,QModel%QM%ndim,nderiv=0)
    CALL QML_alloc_dnMat(Vec_loc0,   QModel%QM%nsurf,QModel%QM%ndim,nderiv=0)

    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0,vec=Vec_loc0)

    PotVal%d0 = PotVal_loc0%d0
    Vec%d0    = Vec_loc0%d0
    CALL Init_IdMat(NAC%d0,QModel%QM%nsurf)

    allocate(tVec(QModel%QM%nsurf,QModel%QM%nsurf))
    tVec(:,:)      = transpose(Vec%d0)

    IF (nderiv >= 1) THEN ! 1st derivatives

      ! Numeric evaluation of forces
      DO i=1,QModel%QM%ndim

        Q_loc(i) = Q(i) + step        ! q+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0) ! Ep(q+dq)
        PotVal%d1(:,:,i) = PotVal_loc0%d0
        Vec%d1(:,:,i)    = Vec_loc0%d0

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = PotVal_loc0%d0
          Vec%d2(:,:,i,i)    = Vec_loc0%d0
        END IF

        Q_loc(i) = Q(i) - step        ! q-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0) ! Ep(q-dq)
        PotVal%d1(:,:,i) = (PotVal%d1(:,:,i)-PotVal_loc0%d0)/(TWO*step)
        Vec%d1(:,:,i)    = (Vec%d1(:,:,i)-Vec_loc0%d0)/(TWO*step)

        NAC%d1(:,:,i)   = matmul(tVec,Vec%d1(:,:,i))

        IF (nderiv >= 2) THEN
          PotVal%d2(:,:,i,i) = (PotVal%d2(:,:,i,i) + PotVal_loc0%d0 - TWO*PotVal%d0)/ &
                                step**2
          Vec%d2(:,:,i,i)    = (Vec%d2(:,:,i,i) + Vec_loc0%d0 - TWO*Vec%d0)/ &
                                step**2
        END IF

        Q_loc(i) = Q(i)

      END DO
    END IF


    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        PotVal%d2(:,:,j,i) = PotVal_loc0%d0
        Vec%d2(:,:,j,i)    = Vec_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) + PotVal_loc0%d0
        Vec%d2(:,:,j,i)    = Vec%d2(:,:,j,i)    + Vec_loc0%d0

        Q_loc(i) = Q(i) + step        ! qi+dq
        Q_loc(j) = Q(j) - step        ! qj-dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0
        Vec%d2(:,:,j,i)   = Vec%d2(:,:,j,i)     - Vec_loc0%d0

        Q_loc(i) = Q(i) - step        ! qi-dq
        Q_loc(j) = Q(j) + step        ! qj+dq
        CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i) - PotVal_loc0%d0
        Vec%d2(:,:,j,i)    = Vec%d2(:,:,j,i)    - Vec_loc0%d0

        PotVal%d2(:,:,j,i) = PotVal%d2(:,:,j,i)/(FOUR*step**2)
        PotVal%d2(:,:,i,j) = PotVal%d2(:,:,j,i)
        Vec%d2(:,:,j,i)    = Vec%d2(:,:,j,i)/(FOUR*step**2)
        Vec%d2(:,:,i,j)    = Vec%d2(:,:,j,i)

        Q_loc(i) = Q(i)
        Q_loc(j) = Q(j)
      END DO
      END DO
    END IF


    deallocate(tVec)
    deallocate(Q_loc)
    CALL QML_dealloc_dnMat(PotVal_loc0)
    CALL QML_dealloc_dnMat(Vec_loc0)

  END SUBROUTINE Eval_Pot_Numeric_adia_old
  SUBROUTINE Eval_Pot_Numeric_adia_v3(QModel,Q,PotVal,nderiv,Vec,NAC)
  USE mod_FiniteDiff_coeff
  IMPLICIT NONE

    TYPE (Model_t),    intent(inout)  :: QModel
    TYPE (dnMat_t),    intent(inout)  :: PotVal
    real (kind=Rkind), intent(in)     :: Q(:)
    integer,           intent(in)     :: nderiv
    TYPE (dnMat_t),    intent(inout)  :: Vec,NAC

    ! local variable
    real (kind=Rkind), allocatable     :: Q_loc(:)
    TYPE (dnMat_t)                     :: PotVal_loc0,Vec_loc0

    integer                            :: i,j,k,ip,jp,kp

    integer                            :: i_pt,nb_pts,ind1DQ(1),ind2DQ(2),ind3DQ(3)

     real (kind=Rkind), allocatable     :: tVec(:,:)

    CALL check_alloc_QM(QModel,'Eval_Pot_Numeric_adia_v3')

    !write(6,*) 'coucou Eval_Pot_Numeric_adia_v3' ; flush(6)


    IF (QML_Check_NotAlloc_dnMat(PotVal,nderiv) ) THEN
      CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                           nderiv=nderiv)
    END IF
    PotVal = ZERO

    IF (QML_Check_NotAlloc_dnMat(Vec,nderiv) ) THEN
      CALL QML_alloc_dnMat(Vec,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    Vec = ZERO

    IF (QML_Check_NotAlloc_dnMat(NAC,nderiv) ) THEN
      CALL QML_alloc_dnMat(NAC,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,&
                          nderiv=nderiv)
    END IF
    NAC = ZERO

    allocate(Q_loc(QModel%QM%ndim))
    Q_loc(:) = Q
    CALL QML_alloc_dnMat(PotVal_loc0,QModel%QM%nsurf,QModel%QM%ndim,nderiv=0)
    CALL QML_alloc_dnMat(Vec_loc0,   QModel%QM%nsurf,QModel%QM%ndim,nderiv=0)


    ! no derivative : PotVal%d0
    CALL Eval_Pot_ana(QModel,Q,PotVal_loc0,nderiv=0,vec=Vec_loc0)

    CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0,option=3)
    CALL NumericDerivativesQ_AddMat_TO_dnMat(Vec,   Vec_loc0%d0,   option=3)

    CALL Init_IdMat(NAC%d0,QModel%QM%nsurf)

    allocate(tVec(QModel%QM%nsurf,QModel%QM%nsurf))
    tVec(:,:)      = transpose(Vec%d0)

    IF (nderiv >= 1) THEN ! 1st derivatives

      DO i=1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(1)
          CALL Get_indDQ(ind1DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i],indDQ=ind1DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                           indQ=[i],indDQ=ind1DQ,option=3)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(Vec,Vec_loc0%d0,       &
                                           indQ=[i],indDQ=ind1DQ,option=3)
        END DO

        NAC%d1(:,:,i)                      = matmul(tVec,Vec%d1(:,:,i))
        IF (nderiv >= 2) NAC%d2(:,:,i,i)   = matmul(tVec,Vec%d2(:,:,i,i))
        IF (nderiv >= 3) NAC%d3(:,:,i,i,i) = matmul(tVec,Vec%d3(:,:,i,i,i))


      END DO
    END IF

    IF (nderiv >= 2) THEN ! 2d derivatives

      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(2)
          CALL Get_indDQ(ind2DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j],indDQ=ind2DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                          indQ=[i,j],indDQ=ind2DQ,option=3)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(Vec,Vec_loc0%d0,       &
                                          indQ=[i,j],indDQ=ind2DQ,option=3)
        END DO

        PotVal%d2(:,:,i,j) = PotVal%d2(:,:,j,i)
        Vec%d2(:,:,i,j)    = Vec%d2(:,:,j,i)

        NAC%d2(:,:,j,i)    = matmul(tVec,Vec%d2(:,:,j,i))
        NAC%d2(:,:,i,j)    = NAC%d2(:,:,j,i)

        IF (nderiv >= 3) THEN
          PotVal%d3(:,:,i,j,i) = PotVal%d3(:,:,i,i,j)
          PotVal%d3(:,:,j,i,i) = PotVal%d3(:,:,i,i,j)
          PotVal%d3(:,:,i,j,j) = PotVal%d3(:,:,j,j,i)
          PotVal%d3(:,:,j,i,j) = PotVal%d3(:,:,j,j,i)

          Vec%d3(:,:,i,j,i)    = Vec%d3(:,:,i,i,j)
          Vec%d3(:,:,j,i,i)    = Vec%d3(:,:,i,i,j)
          Vec%d3(:,:,i,j,j)    = Vec%d3(:,:,j,j,i)
          Vec%d3(:,:,j,i,j)    = Vec%d3(:,:,j,j,i)

          NAC%d3(:,:,i,i,j)    = matmul(tVec,Vec%d3(:,:,i,i,j))
          NAC%d3(:,:,i,j,i)    = NAC%d3(:,:,i,i,j)
          NAC%d3(:,:,j,i,i)    = NAC%d3(:,:,i,i,j)
          NAC%d3(:,:,j,j,i)    = matmul(tVec,Vec%d3(:,:,j,j,i))
          NAC%d3(:,:,i,j,j)    = NAC%d3(:,:,j,j,i)
          NAC%d3(:,:,j,i,j)    = NAC%d3(:,:,j,j,i)
        END IF

      END DO
      END DO
    END IF

    IF (nderiv >= 3) THEN ! 3d derivatives: d3/dQidQidQj

      ! d3/dQidQjdQk
      DO i=1,QModel%QM%ndim
      DO j=i+1,QModel%QM%ndim
      DO k=j+1,QModel%QM%ndim

        DO i_pt=1,Get_nb_pts(3)
          CALL Get_indDQ(ind3DQ,i_pt)
          CALL Set_QplusDQ(Q_loc,Q,indQ=[i,j,k],indDQ=ind3DQ,step_sub=step)
          CALL Eval_Pot_ana(QModel,Q_loc,PotVal_loc0,nderiv=0,vec=Vec_loc0)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(PotVal,PotVal_loc0%d0, &
                                       indQ=[i,j,k],indDQ=ind3DQ,option=3)
          CALL NumericDerivativesQ_AddMat_TO_dnMat(Vec,Vec_loc0%d0,       &
                                       indQ=[i,j,k],indDQ=ind3DQ,option=3)
        END DO

        PotVal%d3(:,:,k,i,j) = PotVal%d3(:,:,k,j,i)
        PotVal%d3(:,:,j,k,i) = PotVal%d3(:,:,k,j,i)
        PotVal%d3(:,:,i,j,k) = PotVal%d3(:,:,k,j,i)
        PotVal%d3(:,:,i,k,j) = PotVal%d3(:,:,k,j,i)
        PotVal%d3(:,:,j,i,k) = PotVal%d3(:,:,k,j,i)

        Vec%d3(:,:,k,i,j)    = Vec%d3(:,:,k,j,i)
        Vec%d3(:,:,j,k,i)    = Vec%d3(:,:,k,j,i)
        Vec%d3(:,:,i,j,k)    = Vec%d3(:,:,k,j,i)
        Vec%d3(:,:,i,k,j)    = Vec%d3(:,:,k,j,i)
        Vec%d3(:,:,j,i,k)    = Vec%d3(:,:,k,j,i)

        NAC%d3(:,:,k,j,i)    = matmul(tVec,Vec%d3(:,:,k,j,i))
        NAC%d3(:,:,k,i,j)    = NAC%d3(:,:,k,j,i)
        NAC%d3(:,:,j,k,i)    = NAC%d3(:,:,k,j,i)
        NAC%d3(:,:,i,j,k)    = NAC%d3(:,:,k,j,i)
        NAC%d3(:,:,i,k,j)    = NAC%d3(:,:,k,j,i)
        NAC%d3(:,:,j,i,k)    = NAC%d3(:,:,k,j,i)

      END DO
      END DO
      END DO
    END IF

    IF (nderiv >= 1) PotVal%d1 = PotVal%d1/step
    IF (nderiv >= 2) PotVal%d2 = PotVal%d2/step**2
    IF (nderiv >= 3) PotVal%d3 = PotVal%d3/step**3

    IF (nderiv >= 1) Vec%d1    = Vec%d1/step
    IF (nderiv >= 2) Vec%d2    = Vec%d2/step**2
    IF (nderiv >= 3) Vec%d3    = Vec%d3/step**3

    IF (nderiv >= 1) NAC%d1    = NAC%d1/step
    IF (nderiv >= 2) NAC%d2    = NAC%d2/step**2
    IF (nderiv >= 3) NAC%d3    = NAC%d3/step**3

    deallocate(tVec)
    deallocate(Q_loc)
    CALL QML_dealloc_dnMat(PotVal_loc0)
    CALL QML_dealloc_dnMat(Vec_loc0)

  END SUBROUTINE Eval_Pot_Numeric_adia_v3
  SUBROUTINE dia_TO_adia(PotVal_dia,PotVal_adia,Vec,Vec0,NAC,nderiv)
    USE mod_diago
    IMPLICIT NONE

    TYPE (dnMat_t), intent(in)               :: PotVal_dia
    TYPE (dnMat_t), intent(inout)            :: PotVal_adia,Vec,Vec0,NAC

    integer, intent(in), optional            :: nderiv

    ! local variable
    integer                        :: i,j,k,id,jd,kd,nderiv_loc,ndim,nsurf
    real (kind=Rkind)              :: ai,aj,aii,aij,aji,ajj,th,cc,ss,DEne
    real (kind=Rkind), allocatable :: Eig(:),tVec(:,:),Vdum(:),Vi(:)

    TYPE (dnMat_t)                :: PotVal_dia_onadia



    !test DIAG_dnMat
    TYPE (dnMat_t)              :: dnVec,dnDiag,dnMat

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='dia_TO_adia'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------


    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
      flush(out_unitp)
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(3,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF


    IF ( QML_Check_NotAlloc_dnMat(PotVal_dia,nderiv_loc) ) THEN
      write(out_unitp,*) ' The diabatic potential MUST be allocated!'
      CALL QML_Write_dnMat(PotVal_dia)
      STOP 'PotVal_dia%dn NOT allocated in "dia_TO_adia"'
    END IF
    IF (debug) THEN
      write(out_unitp,*) 'PotVal_dia'
      CALL QML_Write_dnMat(PotVal_dia,nio=out_unitp)
      flush(out_unitp)
    END IF

    nsurf = QML_get_nsurf_FROM_dnMat(PotVal_dia)
    ndim  = QML_get_ndim_FROM_dnMat(PotVal_dia)


    IF (QML_Check_NotAlloc_dnMat(Vec0,nderiv=0)) THEN
       !$OMP CRITICAL (CRIT_dia_TO_adia_old)
       CALL QML_alloc_dnMat(Vec0,nsurf=nsurf,ndim=ndim,nderiv=0)

       allocate(Eig(nsurf))

       CALL diagonalization(PotVal_dia%d0,Eig,Vec0%d0,nsurf,sort=1,phase=.TRUE.)

       deallocate(Eig)

       IF (debug) write(out_unitp,*) 'init Vec0 done'

       !$OMP END CRITICAL (CRIT_dia_TO_adia_old)
    END IF


    CALL QML_DIAG_dnMat(dnMat=PotVal_dia,dnMatDiag=PotVal_adia,         &
                        dnVec=Vec,dnVecProj=NAC,dnVec0=Vec0)

    IF (debug) THEN

      write(out_unitp,*) 'Eig',(PotVal_adia%d0(i,i),i=1,nsurf)

      write(out_unitp,*) 'PotVal_adia'
      CALL QML_Write_dnMat(PotVal_adia,nio=out_unitp)

      write(out_unitp,*) 'Vec'
      CALL QML_Write_dnMat(Vec,nio=out_unitp)

      write(out_unitp,*) 'NAC'
      CALL QML_Write_dnMat(NAC,nio=out_unitp)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE dia_TO_adia
  SUBROUTINE dia_TO_adia_old(PotVal_dia,PotVal_adia,Vec,Vec0,NAC,nderiv)
    USE mod_diago
    IMPLICIT NONE

    TYPE (dnMat_t), intent(in)               :: PotVal_dia
    TYPE (dnMat_t), intent(inout)            :: PotVal_adia,Vec,Vec0,NAC

    integer, intent(in), optional            :: nderiv

    ! local variable
    integer                        :: i,j,k,id,jd,kd,nderiv_loc,ndim,nsurf
    real (kind=Rkind)              :: ai,aj,aii,aij,aji,ajj,th,cc,ss,DEne
    real (kind=Rkind), allocatable :: Eig(:),tVec(:,:),Vdum(:),Vi(:)

    TYPE (dnMat_t)                :: PotVal_dia_onadia


!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='dia_TO_adia_old'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      IF (present(nderiv)) write(out_unitp,*) '   nderiv',nderiv
      flush(out_unitp)
    END IF

    IF (present(nderiv)) THEN
      nderiv_loc = max(0,nderiv)
      nderiv_loc = min(3,nderiv_loc)
    ELSE
      nderiv_loc = 0
    END IF


    IF ( QML_Check_NotAlloc_dnMat(PotVal_dia,nderiv_loc) ) THEN
      write(out_unitp,*) ' The diabatic potential MUST be allocated!'
      CALL QML_Write_dnMat(PotVal_dia)
      STOP 'PotVal_dia%dn NOT allocated in "dia_TO_adia_old"'
    END IF
    IF (debug) THEN
      write(out_unitp,*) 'PotVal_dia'
      CALL QML_Write_dnMat(PotVal_dia,nio=out_unitp)
      flush(out_unitp)
    END IF

    nsurf = QML_get_nsurf_FROM_dnMat(PotVal_dia)
    ndim  = QML_get_ndim_FROM_dnMat(PotVal_dia)

    IF ( QML_Check_NotAlloc_dnMat(PotVal_adia,nderiv_loc) ) THEN
      CALL QML_alloc_dnMat(PotVal_adia,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    END IF
    PotVal_adia = ZERO


    IF ( QML_Check_NotAlloc_dnMat(Vec,nderiv_loc) ) THEN
      CALL QML_alloc_dnMat(Vec,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    END IF
    Vec = ZERO

    IF ( QML_Check_NotAlloc_dnMat(NAC,nderiv_loc) ) THEN
      CALL QML_alloc_dnMat(NAC,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    END IF
    NAC = ZERO


    ! local variables
    CALL QML_alloc_dnMat(PotVal_dia_onadia,nsurf=nsurf,ndim=ndim,nderiv=nderiv_loc)
    PotVal_dia_onadia = ZERO


    allocate(Eig(nsurf))
    allocate(tVec(nsurf,nsurf))
    CALL diagonalization(PotVal_dia%d0,Eig,Vec%d0,nsurf)
    IF (QML_Check_NotAlloc_dnMat(Vec0,nderiv=0)) THEN
       !$OMP CRITICAL (CRIT_dia_TO_adia_old)
       IF (debug) write(out_unitp,*) 'init Vec0'
       CALL QML_alloc_dnMat(Vec0,nsurf=nsurf,ndim=ndim,nderiv=0)
       Vec0%d0 = Vec%d0
       !$OMP END CRITICAL (CRIT_dia_TO_adia_old)
    ELSE ! change the phase if required
       IF (debug) write(out_unitp,*) 'Change phase?'
       flush(out_unitp)

       DO i=1,nsurf
         IF (dot_product(Vec0%d0(:,i),Vec%d0(:,i)) < ZERO) Vec%d0(:,i) = -Vec%d0(:,i)
       END DO

       IF (debug) THEN
         write(out_unitp,*) 'Vec before rotation'
         CALL QML_Write_dnMat(Vec,nio=out_unitp)
       END IF
       !For degenerated eigenvectors (works only with 2 vectors)
       DO i=1,nsurf-1
        IF ( abs(Eig(i)-Eig(i+1)) < epsi) THEN
          j = i+1
          IF (debug) write(out_unitp,*) 'degenerated vectors',i,j

          aii = dot_product(Vec0%d0(:,i),Vec%d0(:,i))
          aji = dot_product(Vec0%d0(:,j),Vec%d0(:,i))
          aij = dot_product(Vec0%d0(:,i),Vec%d0(:,j))
          ajj = dot_product(Vec0%d0(:,j),Vec%d0(:,j))

          th = ( atan2(aij,ajj) -atan2(aji,aii) ) * HALF
          IF (debug) write(out_unitp,*) 'theta',th

          cc = cos(th)
          ss = sin(th)

          DO k=1,nsurf
           ai = Vec%d0(k,i)
           aj = Vec%d0(k,j)

           Vec%d0(k,i) =  cc * ai + ss * aj
           Vec%d0(k,j) = -ss * ai + cc * aj

          END DO

        END IF
       END DO

       IF (debug) THEN
         write(out_unitp,*) 'Vec after rotation'
         CALL QML_Write_dnMat(Vec,nio=out_unitp)
       END IF

    END IF
    tVec = transpose(Vec%d0)

    IF (debug) write(out_unitp,*) 'Eig',Eig

    ! transformation of PotVal_dia on the adiabatic basis (Vec)
    !PotVal_dia_onadia%d0 = matmul(tVec,matmul(PotVal_dia%d0,Vec%d0))
    DO i=1,nsurf
      PotVal_dia_onadia%d0(i,i) = Eig(i)
    END DO
    IF (nderiv_loc > 0) THEN
      DO id=1,ndim
        PotVal_dia_onadia%d1(:,:,id) = matmul(tVec,matmul(PotVal_dia%d1(:,:,id),Vec%d0))
      END DO
    END IF
    IF (nderiv_loc > 1) THEN
      DO jd=1,ndim
      DO id=1,ndim
        PotVal_dia_onadia%d2(:,:,id,jd) =                               &
                    matmul(tVec,matmul(PotVal_dia%d2(:,:,id,jd),Vec%d0) )
      END DO
      END DO
    END IF
    IF (nderiv_loc > 2) THEN
      DO kd=1,ndim
      DO jd=1,ndim
      DO id=1,ndim
        PotVal_dia_onadia%d3(:,:,id,jd,kd) =                            &
                 matmul(tVec,matmul(PotVal_dia%d3(:,:,id,jd,kd),Vec%d0) )
      END DO
      END DO
      END DO
    END IF
    deallocate(Eig)
    deallocate(tVec)
    IF (debug) THEN
      write(out_unitp,*) "< Psi_j I dnHdia I Psi_i>"
      CALL QML_Write_dnMat(PotVal_dia_onadia,nio=out_unitp)
      flush(out_unitp)
    END IF


    ! project the eigenvector derivatives on the eigenvectors => NAC

    ! no derivative
    PotVal_adia%d0 = ZERO
    NAC%d0         = ZERO
    DO i=1,nsurf
      PotVal_adia%d0(i,i) = PotVal_dia_onadia%d0(i,i)
      NAC%d0(i,i)         = ONE
    END DO

    ! 1st order derivatives
    IF (nderiv_loc > 0) THEN

      ! eigenvalue derivatives
      PotVal_adia%d1 = ZERO
      DO id=1,ndim
      DO i=1,nsurf
        PotVal_adia%d1(i,i,id) = PotVal_dia_onadia%d1(i,i,id)
      END DO
      END DO

      ! eigenvector derivatives projected on the eigenvectors
      DO id=1,ndim
      DO i=1,nsurf ! I Psi_i' >
      DO j=1,nsurf ! projection on < Psi_j I
        DEne = PotVal_adia%d0(j,j) - PotVal_adia%d0(i,i)
        IF (j /= i .AND. abs(DEne) > ONETENTH**10) THEN
          NAC%d1(j,i,id) = - PotVal_dia_onadia%d1(j,i,id) / DEne
        ELSE
          NAC%d1(j,i,id) = ZERO
        END IF
      END DO
      END DO
      END DO

    END IF


    ! 2d order derivatives
    IF (nderiv_loc > 1) THEN
      PotVal_adia%d2 = ZERO
      allocate(Vdum(nsurf))
      ! eigenvalue derivatives
      DO id=1,ndim
      DO jd=1,ndim
        DO i=1,nsurf
          Vdum = matmul(PotVal_dia_onadia%d2(:,:,id,jd),NAC%d0(:,i))    + &
                 matmul(PotVal_dia_onadia%d1(:,:,jd),   NAC%d1(:,i,id)) + &
                 matmul(PotVal_dia_onadia%d1(:,:,id),   NAC%d1(:,i,jd))

          !write(out_unitp,*) 'Dum',id,jd,i,':',Vdum
          PotVal_adia%d2(i,i,id,jd) = dot_product(NAC%d0(:,i),Vdum)
        END DO

      END DO
      END DO
      deallocate(Vdum)
    END IF

    ! 3d order derivatives (it is not working, because NAC%d2 are not calculated)
    IF (nderiv_loc > 2) THEN
      PotVal_adia%d3 = ZERO
      allocate(Vdum(nsurf))
      ! eigenvalue derivatives
      DO id=1,ndim
      DO jd=1,ndim
      DO kd=1,ndim
        DO i=1,nsurf
          Vdum = matmul(PotVal_dia_onadia%d3(:,:,id,jd,kd),NAC%d0(:,i))       + &
                 matmul(PotVal_dia_onadia%d2(:,:,jd,kd),   NAC%d1(:,i,id))    + &
                 matmul(PotVal_dia_onadia%d2(:,:,id,jd),   NAC%d1(:,i,kd))    + &
                 matmul(PotVal_dia_onadia%d2(:,:,kd,id),   NAC%d1(:,i,jd))    + &
                 matmul(PotVal_dia_onadia%d1(:,:,id),      NAC%d2(:,i,jd,kd)) + &
                 matmul(PotVal_dia_onadia%d1(:,:,jd),      NAC%d2(:,i,kd,id)) + &
                 matmul(PotVal_dia_onadia%d1(:,:,kd),      NAC%d2(:,i,id,jd))
          !write(out_unitp,*) 'Dum',id,jd,i,':',Vdum
          PotVal_adia%d3(i,i,id,jd,kd) = dot_product(NAC%d0(:,i),Vdum)
        END DO
      END DO
      END DO
      END DO
      deallocate(Vdum)
    END IF

    CALL QML_dealloc_dnMat(PotVal_dia_onadia)

    IF (debug) THEN
      write(out_unitp,*) 'PotVal_adia'
      CALL QML_Write_dnMat(PotVal_adia,nio=out_unitp)

      write(out_unitp,*) 'Vec'
      CALL QML_Write_dnMat(Vec,nio=out_unitp)

      write(out_unitp,*) 'NAC'
      CALL QML_Write_dnMat(NAC,nio=out_unitp)
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE dia_TO_adia_old
  SUBROUTINE Write_Model(QModel,nio)
  USE mod_Lib
  IMPLICIT NONE

    TYPE(Model_t),     intent(in)              :: QModel
    integer,            intent(in), optional    :: nio

    integer :: nio_loc

    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    CALL check_alloc_QM(QModel,'Write_Model')

    IF (nio_loc /= out_unitp) THEN
      open(nio_loc,file=trim(adjustl(QModel%QM%pot_name))//'.out',form='formatted')
    END IF


    write(nio_loc,*) '-----------------------------------------------'
    write(nio_loc,*) 'Output file for potential library'

    CALL QModel%QM%Write_QModel(nio=nio_loc)

    write(nio_loc,*) '-----------------------------------------------'
    flush(nio_loc)


  END SUBROUTINE Write_Model
  SUBROUTINE Write0_Model(QModel,nio)
  USE mod_Lib
  IMPLICIT NONE

    TYPE(Model_t),   intent(in)              :: QModel
    integer,          intent(in), optional    :: nio

    integer :: nio_loc


    IF (present(nio)) THEN
      nio_loc = nio
    ELSE
      nio_loc = out_unitp
    END IF

    CALL check_alloc_QM(QModel,'Write0_Model')

    IF (nio_loc /= out_unitp) THEN
      open(nio_loc,file=trim(adjustl(QModel%QM%pot_name))//'.out',form='formatted')
    END IF


    CALL QModel%QM%Write0_QModel(nio=nio_loc)


     IF (nio_loc /= out_unitp) THEN
      close(nio_loc)
    END IF


  END SUBROUTINE Write0_Model
  SUBROUTINE Write_QdnV_FOR_Model(Q,PotVal,QModel,Vec,NAC,info)
  USE mod_Lib
  IMPLICIT NONE

    TYPE (Model_t),   intent(in)            :: QModel
    TYPE (dnMat_t),    intent(in)           :: PotVal
    real (kind=Rkind), intent(in)           :: Q(:)
    TYPE (dnMat_t),    intent(in), optional :: Vec ! for non adiabatic couplings
    TYPE (dnMat_t),    intent(in), optional :: NAC ! for non adiabatic couplings
    character(len=*),  intent(in), optional :: info

    integer :: nio_loc,err_io

    CALL file_open2(trim(adjustl(QModel%QM%pot_name))//'.txt',nio_loc,lformatted=.TRUE.,append=.TRUE.,err_file=err_io)
    IF (err_io /= 0) THEN
      write(out_unitp,*) 'ERROR in Write_QdnV_FOR_Model'
      write(out_unitp,*) ' Impossible to open the file "',trim(adjustl(QModel%QM%pot_name))//'.txt','"'
      STOP 'Impossible to open the file'
    END IF


    write(nio_loc,'(a)',advance='no') 'TEST output: '
    IF (present(info)) THEN
      write(nio_loc,'(a)',advance='no') info
    END IF
    IF (QModel%QM%adiabatic) THEN
      write(nio_loc,'(a)') ' Adiabatic'
    ELSE
      write(nio_loc,'(a)') ' Diabatic'
    END IF

    write(nio_loc,*) 'Q'
    write(nio_loc,*) size(Q)
    write(nio_loc,*) Q

    IF (allocated(PotVal%d0)) THEN
      write(nio_loc,*) 'V'
      write(nio_loc,*) size(PotVal%d0)
      write(nio_loc,*) PotVal%d0
    END IF
    IF (allocated(PotVal%d1)) THEN
      write(nio_loc,*) 'Grad'
      write(nio_loc,*) size(PotVal%d1)
      write(nio_loc,*) PotVal%d1
    END IF
    IF (allocated(PotVal%d2)) THEN
      write(nio_loc,*) 'Hess'
      write(nio_loc,*) size(PotVal%d2)
      write(nio_loc,*) PotVal%d2
    END IF

    IF (present(Vec)) THEN
    IF (allocated(Vec%d0)) THEN
      write(nio_loc,*) 'Vec'
      write(nio_loc,*) size(Vec%d0)
      write(nio_loc,*) Vec%d0
    END IF
    IF (allocated(Vec%d1)) THEN
      write(nio_loc,*) 'd1Vec'
      write(nio_loc,*) size(Vec%d1)
      write(nio_loc,*) Vec%d1
    END IF
    IF (allocated(Vec%d2)) THEN
      write(nio_loc,*) 'd2Vec'
      write(nio_loc,*) size(Vec%d2)
      write(nio_loc,*) Vec%d2
    END IF
    END IF

    IF (present(NAC)) THEN
    IF (allocated(NAC%d1)) THEN
      write(nio_loc,*) 'NAC'
      write(nio_loc,*) size(NAC%d1)
      write(nio_loc,*) NAC%d1
    END IF
    END IF

    IF (allocated(QModel%QM%d0GGdef)) THEN
      write(nio_loc,*) 'd0GGdef'
      write(nio_loc,*) size(QModel%QM%d0GGdef)
      write(nio_loc,*) QModel%QM%d0GGdef
    END IF

    write(nio_loc,'(a)',advance='no') 'END_TEST output: '
    IF (present(info)) THEN
      write(nio_loc,'(a)',advance='no') info
    END IF
    IF (QModel%QM%adiabatic) THEN
      write(nio_loc,'(a)') ' Adiabatic'
    ELSE
      write(nio_loc,'(a)') ' Diabatic'
    END IF

    close(nio_loc)

  END SUBROUTINE Write_QdnV_FOR_Model
  SUBROUTINE Check_analytical_numerical_derivatives(QModel,Q,nderiv)
  IMPLICIT NONE

    TYPE (Model_t),       intent(inout)   :: QModel
    real (kind=Rkind),    intent(in)      :: Q(:)
    integer,              intent(in)      :: nderiv

    ! local variables
    TYPE (dnMat_t)            :: Mat_diff
    TYPE (dnMat_t)            :: PotVal_ana,PotVal_num
    TYPE (dnMat_t)            :: NAC_ana,NAC_num
    TYPE (dnMat_t)            :: Vec_ana,Vec_num

    real (kind=Rkind)         :: MaxMat,MaxDiffMat

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Check_analytical_numerical_derivatives'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

    IF (debug) THEN
      write(out_unitp,*) ' BEGINNING ',name_sub
      write(out_unitp,*) '   nderiv    ',nderiv
      flush(out_unitp)
    END IF

    CALL check_alloc_QM(QModel,name_sub)


    CALL QML_alloc_dnMat(PotVal_ana,nsurf=QModel%QM%nsurf,              &
                         ndim=QModel%QM%ndim,nderiv=nderiv)

    CALL QML_alloc_dnMat(PotVal_num,nsurf=QModel%QM%nsurf,              &
                         ndim=QModel%QM%ndim,nderiv=nderiv)


    IF (QModel%QM%adiabatic) THEN
      CALL Eval_Pot(QModel,Q,PotVal_ana,nderiv,NAC_ana,Vec_ana,numeric=.FALSE.)
    ELSE
      CALL Eval_Pot(QModel,Q,PotVal_ana,nderiv,numeric=.FALSE.)
    END IF

    IF (debug) THEN
      write(out_unitp,*)   'PotVal_ana'
      CALL QML_Write_dnMat(PotVal_ana,nio=out_unitp)
      flush(out_unitp)
    END IF

    IF (QModel%QM%adiabatic) THEN
      CALL Eval_Pot(QModel,Q,PotVal_num,nderiv,NAC_num,Vec_num,numeric=.TRUE.)
    ELSE
      CALL Eval_Pot(QModel,Q,PotVal_num,nderiv,numeric=.TRUE.)
    END IF
    IF (debug) THEN
      write(out_unitp,*)   'PotVal_num'
      CALL QML_Write_dnMat(PotVal_num,nio=out_unitp)
      flush(out_unitp)
    END IF


    MaxMat      = QML_get_maxval_OF_dnMat(PotVal_ana)
    IF (MaxMat < ONETENTH**6) MaxMat = ONE
    Mat_diff    = PotVal_num - PotVal_ana
    MaxDiffMat  = QML_get_maxval_OF_dnMat(Mat_diff)

    write(out_unitp,'(3a,e9.2)') 'With ',QModel%QM%pot_name,            &
               ': max of the relative Potential diff:',MaxDiffMat/MaxMat
    write(out_unitp,'(3a,l9)')   'With ',QModel%QM%pot_name,            &
     ': Potential diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

    IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
      write(out_unitp,*)   'Potential diff (ana-numer)'
      CALL QML_Write_dnMat(Mat_diff,nio=out_unitp)
    END IF

    IF (QModel%QM%adiabatic) THEN

      MaxMat      = QML_get_maxval_OF_dnMat(NAC_ana)
      IF (MaxMat < ONETENTH**6) MaxMat = ONE
      Mat_diff    = NAC_num - NAC_ana
      MaxDiffMat  = QML_get_maxval_OF_dnMat(Mat_diff)

      write(out_unitp,'(3a,e9.2)') 'With ',QModel%QM%pot_name,          &
                 ': max of the relative NAC diff:',MaxDiffMat/MaxMat
      write(out_unitp,'(3a,l9)')   'With ',QModel%QM%pot_name,          &
       ': NAC diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

      IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
        write(out_unitp,*)   'NAC diff (ana-numer)'
        CALL QML_Write_dnMat(Mat_diff,nio=out_unitp)
        write(out_unitp,*)   'NAC_ana'
        CALL QML_Write_dnMat(NAC_ana,nio=out_unitp)
        write(out_unitp,*)   'NAC_num'
        CALL QML_Write_dnMat(NAC_num,nio=out_unitp)
      END IF

      MaxMat      = QML_get_maxval_OF_dnMat(Vec_ana)
      IF (MaxMat < ONETENTH**6) MaxMat = ONE
      Mat_diff    = Vec_num - Vec_ana
      MaxDiffMat  = QML_get_maxval_OF_dnMat(Mat_diff)

      write(out_unitp,'(3a,e9.2)') 'With ',QModel%QM%pot_name,            &
                 ': max of the relative Vec diff:',MaxDiffMat/MaxMat
      write(out_unitp,'(3a,l9)')   'With ',QModel%QM%pot_name,            &
       ': Vec diff (numer-ana), ZERO?  ',(MaxDiffMat/MaxMat <= step)

      IF (MaxDiffMat/MaxMat > step .OR. debug) THEN
        write(out_unitp,*)   'Vec diff (ana-numer)'
        CALL QML_Write_dnMat(Mat_diff,nio=out_unitp)
        write(out_unitp,*)   'Vec_ana'
        CALL QML_Write_dnMat(Vec_ana,nio=out_unitp)
        write(out_unitp,*)   'Vec_num'
        CALL QML_Write_dnMat(Vec_num,nio=out_unitp)
      END IF

    END IF

    CALL QML_dealloc_dnMat(PotVal_ana)
    CALL QML_dealloc_dnMat(PotVal_num)
    CALL QML_dealloc_dnMat(NAC_ana)
    CALL QML_dealloc_dnMat(NAC_num)
    CALL QML_dealloc_dnMat(Vec_ana)
    CALL QML_dealloc_dnMat(Vec_num)
    CALL QML_dealloc_dnMat(Mat_diff)

    IF (debug) THEN
      write(out_unitp,*) ' END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE Check_analytical_numerical_derivatives

  SUBROUTINE Eval_pot_ON_Grid(QModel,Qmin,Qmax,nb_points,nderiv,grid_file)
  IMPLICIT NONE

    TYPE (Model_t),              intent(inout)   :: QModel
    real (kind=Rkind),            intent(in)      :: Qmin(:),Qmax(:)
    integer, optional,            intent(in)      :: nb_points,nderiv
    character (len=*), optional,  intent(in)      :: grid_file


    ! local variables
    integer                        :: unit_grid_file
    integer                        :: i,iq,jq,i1,i2,nb_points_loc,nderiv_loc,ndim_loc
    integer,           allocatable :: i_Q(:)
    real (kind=Rkind), allocatable :: dQ(:),Q(:)
    TYPE (dnMat_t)                 :: PotVal,NAC


    CALL check_alloc_QM(QModel,'Eval_pot_ON_Grid')


    IF (size(Qmin) /= QModel%QM%ndim .OR. size(Qmax) /= QModel%QM%ndim) THEN
       write(out_unitp,*) ' ERROR in Eval_pot_ON_Grid'
       write(out_unitp,*) ' The size of Qmin or Qmax are different from QModel%QM%ndim'
       write(out_unitp,*) '  size(Qmin)    ',size(Qmin)
       write(out_unitp,*) '  size(Qmax)    ',size(Qmax)
       write(out_unitp,*) '  QModel%QM%ndim',QModel%QM%ndim
       write(out_unitp,*) ' => Check the fortran'
       STOP 'ERROR in Eval_pot_ON_Grid: problem with QModel%QM%ndim'
    END IF

    IF (present(grid_file)) THEN
      IF (len_trim(grid_file) == 0) THEN
        unit_grid_file = out_unitp
      ELSE
        unit_grid_file = 99
        open(unit=unit_grid_file,file=trim(grid_file) )
      END IF
    ELSE
      unit_grid_file = out_unitp
    END IF

    nb_points_loc = 100
    IF (present(nb_points)) nb_points_loc = nb_points
    nb_points_loc = max(nb_points_loc,2)

    IF (present(nderiv)) THEN
      nderiv_loc = nderiv
    ELSE
      nderiv_loc = 0
    END IF

    allocate(dQ(QModel%QM%ndim))
    allocate(Q(QModel%QM%ndim))
    allocate(i_Q(QModel%QM%ndim))

    dQ(:)       = (Qmax-Qmin) / real(nb_points_loc-1,kind=Rkind)
    ndim_loc    = 0
    i_Q(:)      = 0
    DO i=1,QModel%QM%ndim
      IF (dQ(i) /= ZERO) THEN
        ndim_loc = ndim_loc + 1
        i_Q(ndim_loc) = i
      END IF
    END DO
    write(out_unitp,*) 'QModel%QM%ndim',QModel%QM%ndim
    write(out_unitp,*) 'QModel%QM%numeric',QModel%QM%numeric
    write(out_unitp,*) 'ndim for the grid',ndim_loc
    write(out_unitp,*) 'i_Q',i_Q(1:ndim_loc)


    Q(:) = Qmin

    IF (ndim_loc == 1) THEN
      i1 = i_Q(1)
      DO iq=0,nb_points_loc-1
        Q(i1) = Qmin(i1) + dQ(i1)*real(iq,kind=Rkind)
        IF (QModel%QM%nsurf > 1 .AND. QModel%QM%adiabatic) THEN
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=max(1,nderiv_loc),NAC=NAC)

          IF (nderiv_loc == 0) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%QM%nsurf),NAC%d1
          ELSE IF (nderiv_loc == 1) THEN
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%QM%nsurf),(PotVal%d1(i,i,:),i=1,QModel%QM%nsurf)
          ELSE
            write(unit_grid_file,*) Q(i1),(PotVal%d0(i,i),i=1,QModel%QM%nsurf),(PotVal%d1(i,i,:),i=1,QModel%QM%nsurf), &
                            (PotVal%d2(i,i,:,:),i=1,QModel%QM%nsurf)
          END IF
          flush(unit_grid_file)

        ELSE
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv_loc)

          IF (nderiv_loc == 0) THEN
            write(unit_grid_file,*) Q(i1),PotVal%d0
          ELSE IF (nderiv_loc == 1) THEN
            write(unit_grid_file,*) Q(i1),PotVal%d0,PotVal%d1
          ELSE
            write(unit_grid_file,*) Q(i1),PotVal%d0,PotVal%d1,PotVal%d2
          END IF
          flush(unit_grid_file)

        END IF
      END DO
    ELSE IF (ndim_loc == 2) THEN
      i1 = i_Q(1)
      i2 = i_Q(2)
      DO iq=0,nb_points_loc-1
      DO jq=0,nb_points_loc-1
        Q(i1) = Qmin(i1) + dQ(i1)*real(iq,kind=Rkind)
        Q(i2) = Qmin(i2) + dQ(i2)*real(jq,kind=Rkind)

        IF (QModel%QM%nsurf > 1 .AND. QModel%QM%adiabatic) THEN
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=0,NAC=NAC)
          write(unit_grid_file,*) Q(i1),Q(i2),(PotVal%d0(i,i),i=1,QModel%QM%nsurf)
        ELSE
          CALL Eval_Pot(QModel,Q,PotVal,nderiv=nderiv_loc)
          write(unit_grid_file,*) Q(i1),Q(i2),PotVal%d0
        END IF
        flush(unit_grid_file)

      END DO
      write(unit_grid_file,*)
      END DO


    END IF

    CALL QML_dealloc_dnMat(PotVal)
    CALL QML_dealloc_dnMat(NAC)
    deallocate(dQ)
    deallocate(Q)
    deallocate(i_Q)

    IF (unit_grid_file /= out_unitp) THEN
      close(unit_grid_file)
    END IF


  END SUBROUTINE Eval_pot_ON_Grid


  SUBROUTINE calc_pot(V,QModel,Q)
  IMPLICIT NONE

    TYPE (Model_t),        intent(inout)   :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated

    TYPE (dnMat_t)                         :: PotVal

    CALL check_alloc_QM(QModel,'calc_pot')


    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,nderiv=0)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=0)

    V = PotVal%d0

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot
  SUBROUTINE calc_pot_grad(V,g,QModel,Q)
  IMPLICIT NONE

    TYPE (Model_t),        intent(inout)   :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(QModel,'calc_pot_grad')


    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,nderiv=1)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=1)

    V = PotVal%d0
    g = PotVal%d1

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot_grad
  SUBROUTINE calc_grad(g,QModel,Q)
  IMPLICIT NONE

    TYPE (Model_t),        intent(inout)   :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated

    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(QModel,'calc_grad')


    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,nderiv=1)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=1)

    g = PotVal%d1

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_grad
  SUBROUTINE calc_pot_grad_hess(V,g,h,QModel,Q)
  IMPLICIT NONE

    TYPE (Model_t),        intent(inout)   :: QModel
    real (kind=Rkind),      intent(in)      :: Q(:)
    real (kind=Rkind),      intent(inout)   :: V(:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: g(:,:,:) ! it has to be allocated
    real (kind=Rkind),      intent(inout)   :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(QModel,'calc_pot_grad_hess')


    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,nderiv=2)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)

    V = PotVal%d0
    g = PotVal%d1
    h = PotVal%d2

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_pot_grad_hess
  SUBROUTINE calc_hess(h,QModel,Q)
  IMPLICIT NONE

    TYPE (Model_t),     intent(inout)     :: QModel
    real (kind=Rkind),  intent(in)        :: Q(:)
    real (kind=Rkind),  intent(inout)     :: h(:,:,:,:) ! it has to be allocated


    TYPE (dnMat_t)           :: PotVal

    CALL check_alloc_QM(QModel,'calc_hess')


    CALL QML_alloc_dnMat(PotVal,nsurf=QModel%QM%nsurf,ndim=QModel%QM%ndim,nderiv=2)

    CALL Eval_Pot(QModel,Q,PotVal,nderiv=2)

    h = PotVal%d2

    CALL QML_dealloc_dnMat(PotVal)


  END SUBROUTINE calc_hess
END MODULE mod_Model
