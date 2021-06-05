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
!    Copyright 2016 David Lauvergnat
!      with contributions of Félix MOUHAT and Liang LIANG
!
!    This particular module has been modified from ElVibRot-Tnum:
!       <http://pagesperso.lcp.u-psud.fr/lauvergnat/ElVibRot/ElVibRot.html>
!===========================================================================
!===========================================================================
MODULE QML_FiniteDiff_m
  USE QML_NumParameters_m
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

PUBLIC :: Get_nb_pts,Get_indDQ,Set_QplusDQ
PUBLIC :: FiniteDiff_AddMat_TO_dnMat, FiniteDiff3_SymPerm_OF_dnMat, FiniteDiff_Finalize_dnMat

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


  SUBROUTINE FiniteDiff_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ,option)
  USE QML_dnMat_m
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

        CALL FiniteDiff4_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ)

      ELSE IF (.NOT. present(indQ) .AND. .NOT. present(indDQ)) THEN

        CALL FiniteDiff4_AddMat_TO_dnMat(dnMat,Mat)

      ELSE
        write(out_unitp,*) ' ERROR in FiniteDiff_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' Both indQ and indDQ MUST be present'
        write(out_unitp,*) '     or '
        write(out_unitp,*) ' Both indQ and indDQ MUST be absent'
        write(out_unitp,*) ' present(indQ) ',present(indQ)
        write(out_unitp,*) ' present(indDQ)',present(indDQ)
        STOP 'STOP in FiniteDiff_AddMat_TO_dnMat: Inconsitent parameters.'
      END IF
    ELSE ! option_loc == 3
      IF (present(indQ) .AND. present(indDQ)) THEN

        CALL FiniteDiff3_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ)

      ELSE IF (.NOT. present(indQ) .AND. .NOT. present(indDQ)) THEN

        CALL FiniteDiff3_AddMat_TO_dnMat(dnMat,Mat)

      ELSE
        write(out_unitp,*) ' ERROR in FiniteDiff_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' Both indQ and indDQ MUST be present'
        write(out_unitp,*) '     or '
        write(out_unitp,*) ' Both indQ and indDQ MUST be absent'
        write(out_unitp,*) ' present(indQ) ',present(indQ)
        write(out_unitp,*) ' present(indDQ)',present(indDQ)
        STOP 'STOP in FiniteDiff_AddMat_TO_dnMat: Inconsitent parameters.'
      END IF
    END IF

  END SUBROUTINE FiniteDiff_AddMat_TO_dnMat

  SUBROUTINE FiniteDiff4_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ)
  USE QML_dnMat_m
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
      write(out_unitp,*) ' ERROR in FiniteDiff4_AddMat_TO_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' nderiv < 0',nderiv
      write(out_unitp,*) '  => dnMat is not allocated'
      STOP 'STOP in FiniteDiff4_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (.NOT. all(shape(Mat) == shape(dnMat%d0))) THEN
      write(out_unitp,*) ' ERROR in FiniteDiff4_AddMat_TO_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' shape(Mat),shape(dnMat%d0)',shape(Mat),shape(dnMat%d0)
      STOP 'STOP in FiniteDiff4_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (present(indQ) .AND. present(indDQ)) THEN
      IF (size(indQ) /= size(indDQ) .OR. minval(indQ) < 1 .OR. maxval(indQ) > ndim) THEN
        write(out_unitp,*) ' ERROR in FiniteDiff4_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' size(indQ),size(indDQ)',size(indQ),size(indDQ)
        write(out_unitp,*) ' ndim',ndim
        write(out_unitp,*) ' indQ(:)',indQ
        write(out_unitp,*) ' indDQ(:)',indDQ
        STOP 'STOP in FiniteDiff4_AddMat_TO_dnMat: Inconsitent parameters.'
      END IF

    ELSE IF ( (.NOT. present(indQ) .AND.       present(indDQ)) .OR.     &
              (      present(indQ) .AND. .NOT. present(indDQ)) ) THEN
        write(out_unitp,*) ' ERROR in FiniteDiff4_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' Both indQ and indDQ MUST be present'
        write(out_unitp,*) '     or '
        write(out_unitp,*) ' Both indQ and indDQ MUST be absent'
        write(out_unitp,*) ' present(indQ) ',present(indQ)
        write(out_unitp,*) ' present(indDQ)',present(indDQ)
        STOP 'STOP in FiniteDiff4_AddMat_TO_dnMat: Inconsitent parameters.'
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

  END SUBROUTINE FiniteDiff4_AddMat_TO_dnMat
  SUBROUTINE FiniteDiff3_AddMat_TO_dnMat(dnMat,Mat,indQ,indDQ)
  USE QML_dnMat_m
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
      write(out_unitp,*) ' ERROR in FiniteDiff3_AddMat_TO_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' nderiv < 0',nderiv
      write(out_unitp,*) '  => dnMat is not allocated'
      STOP 'STOP in FiniteDiff3_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (.NOT. all(shape(Mat) == shape(dnMat%d0))) THEN
      write(out_unitp,*) ' ERROR in FiniteDiff3_AddMat_TO_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' shape(Mat),shape(dnMat%d0)',shape(Mat),shape(dnMat%d0)
      STOP 'STOP in FiniteDiff3_AddMat_TO_dnMat: Inconsitent parameters.'
    END IF

    IF (present(indQ) .AND. present(indDQ)) THEN
      IF (size(indQ) /= size(indDQ) .OR. minval(indQ) < 1 .OR. maxval(indQ) > ndim) THEN
        write(out_unitp,*) ' ERROR in FiniteDiff3_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' size(indQ),size(indDQ)',size(indQ),size(indDQ)
        write(out_unitp,*) ' ndim',ndim
        write(out_unitp,*) ' indQ(:)',indQ
        write(out_unitp,*) ' indDQ(:)',indDQ
        STOP 'STOP in FiniteDiff3_AddMat_TO_dnMat: Inconsitent parameters.'
      END IF

    ELSE IF ( (.NOT. present(indQ) .AND.       present(indDQ)) .OR.     &
              (      present(indQ) .AND. .NOT. present(indDQ)) ) THEN
        write(out_unitp,*) ' ERROR in FiniteDiff3_AddMat_TO_dnMat'
        write(out_unitp,*) ' Inconsitent parameters.'
        write(out_unitp,*) ' Both indQ and indDQ MUST be present'
        write(out_unitp,*) '     or '
        write(out_unitp,*) ' Both indQ and indDQ MUST be absent'
        write(out_unitp,*) ' present(indQ) ',present(indQ)
        write(out_unitp,*) ' present(indDQ)',present(indDQ)
        STOP 'STOP in FiniteDiff3_AddMat_TO_dnMat: Inconsitent parameters.'
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

  END SUBROUTINE FiniteDiff3_AddMat_TO_dnMat

  SUBROUTINE FiniteDiff3_SymPerm_OF_dnMat(dnMat,indQ)
  USE QML_dnMat_m
  IMPLICIT NONE

    TYPE (dnMat_t),    intent(inout)   :: dnMat
    integer,           intent(in)      :: indQ(:)  ! indexes of the variables along DQ is made

    ! local variables
    integer                            :: ndim,nderiv
    integer                            :: i,j,k

    nderiv = QML_get_nderiv_FROM_dnMat(dnMat)
    ndim   = QML_get_ndim_FROM_dnMat(dnMat)

    IF (nderiv < 0) THEN
      write(out_unitp,*) ' ERROR in FiniteDiff3_SymPerm_OF_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' nderiv < 0',nderiv
      write(out_unitp,*) '  => dnMat is not allocated'
      STOP 'STOP in FiniteDiff3_SymPerm_OF_dnMat: Inconsitent parameters.'
    END IF

    IF (minval(indQ) < 1 .OR. maxval(indQ) > ndim) THEN
      write(out_unitp,*) ' ERROR in FiniteDiff3_SymPerm_OF_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' ndim',ndim
      write(out_unitp,*) ' indQ(:)',indQ
      STOP 'STOP in FiniteDiff3_SymPerm_OF_dnMat: Inconsitent parameters.'
    END IF

    SELECT CASE (size(indQ))

    CASE (2)
      i  = indQ(1)
      j  = indQ(2)

      dnMat%d2(:,:,i,j)    = dnMat%d2(:,:,j,i)

      IF (nderiv >= 3) THEN
        dnMat%d3(:,:,i,j,i)    = dnMat%d3(:,:,i,i,j)
        dnMat%d3(:,:,j,i,i)    = dnMat%d3(:,:,i,i,j)
        dnMat%d3(:,:,i,j,j)    = dnMat%d3(:,:,j,j,i)
        dnMat%d3(:,:,j,i,j)    = dnMat%d3(:,:,j,j,i)
      END IF

    CASE (3)

      i  = indQ(1)
      j  = indQ(2)
      k  = indQ(3)

      dnMat%d3(:,:,k,i,j) = dnMat%d3(:,:,k,j,i)
      dnMat%d3(:,:,j,k,i) = dnMat%d3(:,:,k,j,i)
      dnMat%d3(:,:,i,j,k) = dnMat%d3(:,:,k,j,i)
      dnMat%d3(:,:,i,k,j) = dnMat%d3(:,:,k,j,i)
      dnMat%d3(:,:,j,i,k) = dnMat%d3(:,:,k,j,i)

    END SELECT

  END SUBROUTINE FiniteDiff3_SymPerm_OF_dnMat
  SUBROUTINE FiniteDiff_Finalize_dnMat(dnMat,step)
  USE QML_dnMat_m
  IMPLICIT NONE

    TYPE (dnMat_t),    intent(inout)         :: dnMat
    real (kind=Rkind), intent(in)            :: step

    integer:: nderiv

    nderiv = QML_get_nderiv_FROM_dnMat(dnMat)

    IF (nderiv < 0) THEN
      write(out_unitp,*) ' ERROR in FiniteDiff_Finalize_dnMat'
      write(out_unitp,*) ' Inconsitent parameters.'
      write(out_unitp,*) ' nderiv < 0',nderiv
      write(out_unitp,*) '  => dnMat is not allocated'
      STOP 'STOP in FiniteDiff_Finalize_dnMat: Inconsitent parameters.'
    END IF

    IF (nderiv >= 1) dnMat%d1 = dnMat%d1/step
    IF (nderiv >= 2) dnMat%d2 = dnMat%d2/step**2
    IF (nderiv >= 3) dnMat%d3 = dnMat%d3/step**3

  END SUBROUTINE FiniteDiff_Finalize_dnMat

END MODULE QML_FiniteDiff_m
