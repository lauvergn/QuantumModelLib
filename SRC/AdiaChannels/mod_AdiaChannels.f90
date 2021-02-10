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

MODULE mod_AdiaChannels
  USE mod_NumParameters

  IMPLICIT NONE

  !PRIVATE
  !PUBLIC :: Model_t,Init_Model,Eval_Pot,Eval_Func

  TYPE :: AdiaChannels_t
    integer                 :: nact        = 0
    integer                 :: ninact      = 0
    integer, allocatable    :: list_act(:)
    integer, allocatable    :: list_inact(:)
    integer, allocatable    :: nb_inact(:)
    integer, allocatable    :: nq_inact(:)
  END TYPE AdiaChannels_t
  TYPE :: BasisPrim_t
    integer                 :: nb      = 0
    integer                 :: nq      = 0
    real(kind=Rkind), allocatable :: x(:)
    real(kind=Rkind), allocatable :: w(:)
    real(kind=Rkind), allocatable :: d0gb(:,:)      ! basis functions d0gb(nq,nb)
    real(kind=Rkind), allocatable :: d1gb(:,:,:)    ! basis functions d2gb(nq,nb,1)
    real(kind=Rkind), allocatable :: d2gb(:,:,:,:)  ! basis functions d2gb(nq,nb,1,1)

  END TYPE BasisPrim_t
CONTAINS
  SUBROUTINE Write_BasisPrim(BasisPrim)
  USE mod_Lib

    TYPE(BasisPrim_t),       intent(in)  :: BasisPrim

    write(6,*) '-------------------------------------------------'
    write(6,*) 'Write_BasisPrim'
    write(6,*) 'nb,nq',BasisPrim%nb,BasisPrim%nq
    write(6,*)
    CALL Write_RVec(BasisPrim%x,6,5,name_info='x')
    write(6,*)
    CALL Write_RVec(BasisPrim%w,6,5,name_info='w')
    write(6,*)
    CALL Write_RMat(BasisPrim%d0gb,6,5,name_info='d0gb')
    write(6,*)
    CALL Write_RMat(BasisPrim%d1gb(:,:,1),6,5,name_info='d1gb')
    write(6,*)
    CALL Write_RMat(BasisPrim%d2gb(:,:,1,1),6,5,name_info='d2gb')
    write(6,*) '-------------------------------------------------'

  END SUBROUTINE Write_BasisPrim
  SUBROUTINE Make_BasisPrim(BasisPrim,nb,nq,A,B) ! boxAB
  USE mod_Lib
  USE mod_dnS
  USE mod_dnPoly

    TYPE(BasisPrim_t),       intent(inout)  :: BasisPrim
    integer,                 intent(in)     :: nb
    integer,                 intent(in)     :: nq
    real(kind=Rkind),        intent(in)     :: A,B

    real(kind=Rkind)          :: dx,sx,x0
    TYPE (dnS_t)              :: dnBox
    TYPE (dnS_t)              :: dnx
    integer :: iq,ib,jb
    real(kind=Rkind), ALLOCATABLE          :: S(:,:)

    BasisPrim%nb = nb
    BasisPrim%nq = nq

    x0 = A
    sx = pi/(B-A)

    allocate(BasisPrim%d0gb(nq,nb))
    allocate(BasisPrim%d1gb(nq,nb,1))
    allocate(BasisPrim%d2gb(nq,nb,1,1))

    ! grid
    dx = pi/nq
    allocate(BasisPrim%x(nq))
    BasisPrim%x(:) = [(dx*(iq-HALF),iq=1,nq)]

    ! weight
    allocate(BasisPrim%w(nq))
    BasisPrim%w(:) = dx

    DO iq=1,nq
      dnx = QML_init_dnS(BasisPrim%x(iq),ndim=1,nderiv=2,iQ=1) ! to set up the derivatives
      DO ib=1,nb
        dnBox = QML_dnBox(dnx,ib)
        CALL QML_sub_get_dn_FROM_dnS(dnBox,BasisPrim%d0gb(iq,ib),               &
                              BasisPrim%d1gb(iq,ib,:),BasisPrim%d2gb(iq,ib,:,:))
      END DO
    END DO
    IF (nb == nq) THEN
      BasisPrim%d0gb(:,nb)      = BasisPrim%d0gb(:,nb)      / sqrt(TWO)
      BasisPrim%d1gb(:,nb,:)    = BasisPrim%d1gb(:,nb,:)    / sqrt(TWO)
      BasisPrim%d2gb(:,nb,:,:)  = BasisPrim%d2gb(:,nb,:,:)  / sqrt(TWO)
    END IF

    CALL Scale_BasisPrim(BasisPrim,x0,sx)

    RETURN
    allocate(S(nb,nb))
    DO ib=1,nb
    DO jb=1,nb
      S(jb,ib) = dot_product(BasisPrim%d0gb(:,jb),BasisPrim%w(:)*BasisPrim%d0gb(:,ib))
    END DO
    END DO
    CALL Write_RMat(S,6,5,name_info='S')

    DO ib=1,nb
    DO jb=1,nb
      S(jb,ib) = dot_product(BasisPrim%d0gb(:,jb),BasisPrim%w(:)*BasisPrim%d2gb(:,ib,1,1))
    END DO
    END DO
    CALL Write_RMat(S,6,5,name_info='<d0b|d2b>')


    CALL Write_BasisPrim(BasisPrim)

  END SUBROUTINE Make_BasisPrim
  SUBROUTINE Scale_BasisPrim(BasisPrim,x0,sx)
  USE mod_Lib
  USE mod_dnS
  USE mod_dnPoly

    TYPE(BasisPrim_t),       intent(inout)  :: BasisPrim
    real(kind=Rkind),        intent(in)     :: x0,sx


    BasisPrim%x(:) = x0 + BasisPrim%x(:) / sx
    BasisPrim%w(:) =      BasisPrim%w(:) / sx

    BasisPrim%d0gb(:,:)     = BasisPrim%d0gb(:,:)     * sqrt(sx)
    BasisPrim%d1gb(:,:,:)   = BasisPrim%d1gb(:,:,:)   * sqrt(sx)*sx
    BasisPrim%d2gb(:,:,:,:) = BasisPrim%d2gb(:,:,:,:) * sqrt(sx)*sx*sx


  END SUBROUTINE Scale_BasisPrim
  SUBROUTINE Make_Hinact(Qact,QModel,AdiaChannels)
  USE mod_Model

  real (kind=Rkind), allocatable, intent(in) :: Qact(:)
  TYPE (Model_t),                 intent(in) :: QModel
  TYPE (AdiaChannels_t),          intent(in) :: AdiaChannels

  END SUBROUTINE Make_Hinact
END MODULE mod_AdiaChannels
