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
MODULE AdiaChannels_Basis_m
  USE QDUtil_NumParameters_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: QML_Basis_t,Read_Basis,Write_Basis, &
            GridTOBasis_Basis,BasisTOGrid_Basis

  TYPE :: QML_Basis_t
    integer                         :: nb_basis   = 0

    integer                         :: nb         = 0
    integer                         :: nq         = 0

    integer                         :: symab   = 0
    integer,  allocatable           :: tab_symab(:)


    character (len=:),  allocatable :: name

    real(kind=Rkind),   allocatable :: x(:)
    real(kind=Rkind),   allocatable :: w(:)
    real(kind=Rkind),   allocatable :: d0gb(:,:)      ! basis functions d0gb(nq,nb)
    real(kind=Rkind),   allocatable :: d1gb(:,:,:)    ! basis functions d2gb(nq,nb,1)
    real(kind=Rkind),   allocatable :: d2gb(:,:,:,:)  ! basis functions d2gb(nq,nb,1,1)

    TYPE (QML_Basis_t), pointer     :: tab_basis(:) => null()
    !TYPE (QML_Basis_t), allocatable :: tab_basis(:)

  END TYPE QML_Basis_t

  INTERFACE Basis_IS_allocated
      MODULE PROCEDURE QML_Basis_IS_allocated
  END INTERFACE
  INTERFACE Write_Basis
      MODULE PROCEDURE QML_Write_Basis
  END INTERFACE
  INTERFACE Read_Basis
      MODULE PROCEDURE QML_Read_Basis
  END INTERFACE
  INTERFACE Construct_Basis_Sin
      MODULE PROCEDURE QML_Construct_Basis_Sin
  END INTERFACE
  INTERFACE CheckOrtho_Basis
      MODULE PROCEDURE QML_CheckOrtho_Basis
  END INTERFACE
  INTERFACE Scale_Basis
      MODULE PROCEDURE QML_Scale_Basis
  END INTERFACE

  INTERFACE Set_symab_Basis
      MODULE PROCEDURE QML_Set_symab_Basis
  END INTERFACE

  INTERFACE BasisTOGrid_Basis
      MODULE PROCEDURE QML_BasisTOGrid_Basis
  END INTERFACE

  INTERFACE GridTOBasis_Basis
      MODULE PROCEDURE QML_GridTOBasis_Basis
  END INTERFACE

CONTAINS
  FUNCTION QML_Basis_IS_allocated(Basis) RESULT(alloc)
    IMPLICIT NONE

    TYPE(QML_Basis_t),   intent(in)  :: Basis
    logical                          :: alloc

    alloc =             allocated(Basis%x)
    alloc = alloc .AND. allocated(Basis%w)
    alloc = alloc .AND. allocated(Basis%d0gb)
    alloc = alloc .AND. allocated(Basis%d1gb)
    alloc = alloc .AND. allocated(Basis%d2gb)

  END FUNCTION QML_Basis_IS_allocated
  RECURSIVE SUBROUTINE QML_Write_Basis(Basis)
    USE QDUtil_m,         ONLY : Write_Vec, Write_Mat
    IMPLICIT NONE

    TYPE(QML_Basis_t),       intent(in)  :: Basis

    integer :: ib

    write(out_unit,*) '-------------------------------------------------'
    write(out_unit,*) 'Write_Basis'
    write(out_unit,*) 'nb,nq',Basis%nb,Basis%nq
    write(out_unit,*) 'Numero of the symmetry plan (symab):',Basis%symab
    IF (Basis_IS_allocated(Basis)) THEN
      CALL Write_Vec(Basis%x,out_unit,5,info='x')
      write(out_unit,*)
      CALL Write_Vec(Basis%w,out_unit,5,info='w')
      write(out_unit,*)
      CALL Write_Mat(Basis%d0gb,out_unit,5,info='d0gb')
      write(out_unit,*)
      CALL Write_Mat(Basis%d1gb(:,:,1),out_unit,5,info='d1gb')
      write(out_unit,*)
      CALL Write_Mat(Basis%d2gb(:,:,1,1),out_unit,5,info='d2gb')
    ELSE
      write(out_unit,*) ' Basis tables (x, w, dngb) are not allocated.'
    END IF

    write(out_unit,*)
    write(out_unit,*)
    write(out_unit,*) 'nb_basis',Basis%nb_basis
    !IF (allocated(Basis%tab_basis)) THEN
    IF (associated(Basis%tab_basis)) THEN
      DO ib=1,size(Basis%tab_basis)
        CALL Write_Basis(Basis%tab_basis(ib))
      END DO
    END IF
    write(out_unit,*) '-------------------------------------------------'

  END SUBROUTINE QML_Write_Basis
  RECURSIVE SUBROUTINE QML_Read_Basis(Basis,nio)
    USE QDUtil_m,         ONLY : TO_lowercase
    IMPLICIT NONE

    TYPE(QML_Basis_t),       intent(inout)  :: Basis
    integer,             intent(in)     :: nio



    integer                         :: ib,err_io

    integer                         :: nb,nq,nb_basis,symab
    character (len=Name_len)        :: name
    real(kind=Rkind)                :: A,B,scaleQ,Q0


    NAMELIST /basis_nD/ name,nb_basis,nb,nq,A,B,scaleQ,Q0,symab


    nb_basis  = 0
    nb        = 0
    nq        = 0
    A         = ZERO
    B         = ZERO
    Q0        = ZERO
    scaleQ    = ONE
    name      = '0'
    symab     = -1


    read(nio,nml=basis_nD,IOSTAT=err_io)
    write(out_unit,nml=basis_nD)
    IF (err_io < 0) THEN
      write(out_unit,basis_nD)
      write(out_unit,*) ' ERROR in QML_Read_Basis'
      write(out_unit,*) '  while reading the namelist "basis_nD"'
      write(out_unit,*) ' end of file or end of record'
      write(out_unit,*) ' Probably, you forget a basis set ...'
      write(out_unit,*) ' Check your data !!'
      STOP ' ERROR in Read_Basis: problems with the namelist.'
    END IF
    IF (err_io > 0) THEN
      write(out_unit,basis_nD)
      write(out_unit,*) ' ERROR in QML_Read_Basis'
      write(out_unit,*) '  while reading the namelist "basis_nD"'
      write(out_unit,*) ' Probably, some arguments of namelist are wrong.'
      write(out_unit,*) ' Check your data !!'
      STOP ' ERROR in Read_Basis: problems with the namelist.'
    END IF

    IF (nb_basis > 0) THEN
      Basis%name      = TO_lowercase('direct-product')

      allocate(Basis%tab_basis(nb_basis))
      DO ib=1,nb_basis
        CALL Read_Basis(Basis%tab_basis(ib),nio)
      END DO
      Basis%nb = product(Basis%tab_basis(:)%nb)
      Basis%nq = product(Basis%tab_basis(:)%nq)
    ELSE
      Basis%nb_basis  = nb_basis
      Basis%nb        = nb
      Basis%nq        = nq
      Basis%symab     = symab
      Basis%name      = TO_lowercase(trim(adjustl(name)))

      SELECT CASE (Basis%name)
      CASE ('boxab')
        CALL Construct_Basis_Sin(Basis)
        Q0      = A
        scaleQ  = pi/(B-A)
      CASE default
        STOP 'ERROR in Read_Basis: no default basis.'
      END SELECT

      CALL Scale_Basis(Basis,Q0,scaleQ)
      CALL QML_Set_symab_Basis(Basis)
      CALL CheckOrtho_Basis(Basis,nderiv=-1)

      !CALL Write_Basis(Basis)

    END IF

  END SUBROUTINE QML_Read_Basis
  SUBROUTINE QML_Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
    USE ADdnSVM_m, ONLY : dnS_t,Variable,sub_get_dn,dnBox
    IMPLICIT NONE

    TYPE(QML_Basis_t),       intent(inout)  :: Basis


    real(kind=Rkind)          :: dx,sx,x0
    TYPE (dnS_t)              :: dnBasis
    TYPE (dnS_t)              :: dnx
    integer                   :: iq,ib,jb,nb,nq

    nb = Basis%nb
    nq = Basis%nq
    dx = pi/nq

    ! grid and weight
    Basis%x = [(dx*(iq-HALF),iq=1,nq)]
    Basis%w = [(dx,iq=1,nq)]

    allocate(Basis%d0gb(nq,nb))
    allocate(Basis%d1gb(nq,nb,1))
    allocate(Basis%d2gb(nq,nb,1,1))

    DO iq=1,nq
      dnx = Variable(Basis%x(iq),nVar=1,nderiv=2,iVar=1) ! to set up the derivatives
      DO ib=1,nb
        dnBasis = dnBox(dnx,ib)
        CALL sub_get_dn(dnBasis,Basis%d0gb(iq,ib),               &
                        Basis%d1gb(iq,ib,:),Basis%d2gb(iq,ib,:,:))
      END DO
    END DO
    IF (nb == nq) THEN
      Basis%d0gb(:,nb)      = Basis%d0gb(:,nb)      / sqrt(TWO)
      Basis%d1gb(:,nb,:)    = Basis%d1gb(:,nb,:)    / sqrt(TWO)
      Basis%d2gb(:,nb,:,:)  = Basis%d2gb(:,nb,:,:)  / sqrt(TWO)
    END IF

  END SUBROUTINE QML_Construct_Basis_Sin
  SUBROUTINE QML_Set_symab_Basis(Basis)
    IMPLICIT NONE

    TYPE(QML_Basis_t),       intent(inout)  :: Basis


    real(kind=Rkind)          :: dx,sx,x0
    !TYPE (dnS_t)              :: dnBox
    !TYPE (dnS_t)              :: dnx
    integer                   :: iq,ib,jb,nb,nq

    IF (.NOT. Basis_IS_allocated(Basis)) THEN
      write(out_unit,*) ' ERROR in QML_Set_symab_Basis'
      write(out_unit,*) ' the basis is not allocated.'
      STOP 'ERROR in QML_Set_symab_Basis: the basis is not allocated'
    END IF

    allocate(Basis%tab_symab(Basis%nb))

    SELECT CASE (Basis%symab)
    CASE (-1)
      Basis%tab_symab(:) = -1
    CASE (0,1,2,3,4,5,6,7)
      Basis%tab_symab(:) = 0
      DO ib=2,Basis%nb,2  ! tab_symab = [0 s 0 s 0 s ....] with s=symab
        Basis%tab_symab(ib) = Basis%symab
      END DO
    CASE DEFAULT
      write(out_unit,*) ' ERROR in QML_Set_symab_Basis'
      write(out_unit,*) '  Wrong symab value:',Basis%symab
      write(out_unit,*) ' Its values must be: [-1,0,1...,7]'
      write(out_unit,*) ' CHECK your data!!'
      STOP 'ERROR in QML_Set_symab_Basis: Wrong symab value'
    END SELECT

  END SUBROUTINE QML_Set_symab_Basis
  SUBROUTINE QML_CheckOrtho_Basis(Basis,nderiv)
    USE QDUtil_m,         ONLY : Write_Mat
    IMPLICIT NONE

    TYPE(QML_Basis_t),       intent(in)     :: Basis
    integer,                 intent(in)     :: nderiv

    integer                         :: iq,ib,jb
    real(kind=Rkind), ALLOCATABLE   :: S(:,:)
    real(kind=Rkind), ALLOCATABLE   :: d0bgw(:,:)
    real(kind=Rkind)                :: Sii,Sij


    IF (Basis_IS_allocated(Basis)) THEN
      d0bgw = transpose(Basis%d0gb)
      DO ib=1,Basis%nb
        d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
      END DO

      S = matmul(d0bgw,Basis%d0gb)
      IF (nderiv > -1) CALL Write_Mat(S,out_unit,5,info='S')
      Sii = ZERO
      Sij = ZERO
      DO ib=1,Basis%nb
        IF (abs(S(ib,ib)-ONE) > Sii) Sii = abs(S(ib,ib)-ONE)
        S(ib,ib) = ZERO
      END DO
      Sij = maxval(S)
      write(out_unit,*) 'Sii,Sij',Sii,Sij

      IF (nderiv > 0) THEN
        S = matmul(d0bgw,Basis%d1gb(:,:,1))
        CALL Write_Mat(S,out_unit,5,info='<d0b|d1b>')
      END IF

      IF (nderiv > 1) THEN
        S = matmul(d0bgw,Basis%d2gb(:,:,1,1))
        CALL Write_Mat(S,out_unit,5,info='<d0b|d2b>')
      END IF

    ELSE
      write(out_unit,*) ' WARNNING in QML_CheckOrtho_Basis'
      write(out_unit,*) ' the basis is not allocated.'
    END IF

  END SUBROUTINE QML_CheckOrtho_Basis
  SUBROUTINE QML_Scale_Basis(Basis,x0,sx)
    IMPLICIT NONE

    TYPE(QML_Basis_t),       intent(inout)  :: Basis
    real(kind=Rkind),    intent(in)     :: x0,sx

    IF (Basis%nb_basis == 0 .AND. abs(sx) > ONETENTH**6 .AND. &
        Basis_IS_allocated(Basis)) THEN

      Basis%x(:) = x0 + Basis%x(:) / sx
      Basis%w(:) =      Basis%w(:) / sx

      Basis%d0gb(:,:)     = Basis%d0gb(:,:)     * sqrt(sx)
      Basis%d1gb(:,:,:)   = Basis%d1gb(:,:,:)   * sqrt(sx)*sx
      Basis%d2gb(:,:,:,:) = Basis%d2gb(:,:,:,:) * sqrt(sx)*sx*sx
    ELSE
      write(out_unit,*) ' ERROR in QML_Scale_Basis'
      write(out_unit,*) ' nb_basis > 0     or ...'
      write(out_unit,*) ' sx is too small  or ...'
      write(out_unit,*) ' the basis is not allocated.'
      STOP 'ERROR in Scale_Basis: nb_basis > 0'
    END IF

  END SUBROUTINE QML_Scale_Basis

  SUBROUTINE QML_BasisTOGrid_Basis(G,B,Basis)
    IMPLICIT NONE

    TYPE(QML_Basis_t),    intent(in)    :: Basis
    real(kind=Rkind),     intent(in)    :: B(:)
    real(kind=Rkind),     intent(inout) :: G(:)

    IF (.NOT. Basis_IS_allocated(Basis)) THEN
      write(out_unit,*) ' ERROR in QML_BasisTOGrid_Basis'
      write(out_unit,*) ' the basis is not allocated.'
      STOP 'ERROR in QML_BasisTOGrid_Basis: the basis is not allocated.'
    END IF
    IF (size(B) /= Basis%nb) THEN
      write(out_unit,*) ' ERROR in QML_BasisTOGrid_Basis'
      write(out_unit,*) ' the size of B is different from nb.'
      write(out_unit,*) ' size(B), Basis%nb',size(B),Basis%nb
      STOP 'ERROR in QML_BasisTOGrid_Basis: wrong B size.'
    END IF

    G = matmul(Basis%d0gb,B)

  END SUBROUTINE QML_BasisTOGrid_Basis

  SUBROUTINE QML_GridTOBasis_Basis(B,G,Basis)
    IMPLICIT NONE

    TYPE(QML_Basis_t),   intent(in)    :: Basis
    real(kind=Rkind),    intent(in)    :: G(:)
    real(kind=Rkind),    intent(inout) :: B(:)

    IF (.NOT. Basis_IS_allocated(Basis)) THEN
      write(out_unit,*) ' ERROR in QML_GridTOBasis_Basis'
      write(out_unit,*) ' the basis is not allocated.'
      STOP 'ERROR in QML_GridTOBasis_Basis: the basis is not allocated.'
    END IF
    IF (size(B) /= Basis%nb) THEN
      write(out_unit,*) ' ERROR in QML_GridTOBasis_Basis'
      write(out_unit,*) ' the size of B is different from nb.'
      write(out_unit,*) ' size(B), Basis%nb',size(B),Basis%nb
      STOP 'ERROR in QML_GridTOBasis_Basis: wrong B size.'
    END IF
    IF (size(G) /= Basis%nq) THEN
      write(out_unit,*) ' ERROR in QML_GridTOBasis_Basis'
      write(out_unit,*) ' the size of G is different from nq.'
      write(out_unit,*) ' size(G), Basis%nq',size(G),Basis%nq
      STOP 'ERROR in QML_GridTOBasis_Basis: wrong G size.'
    END IF

    B(:) = matmul(transpose(Basis%d0gb),Basis%w*G)

  END SUBROUTINE QML_GridTOBasis_Basis

END MODULE AdiaChannels_Basis_m
