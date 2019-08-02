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
MODULE mod_diago
!$ USE omp_lib
  IMPLICIT NONE

   PRIVATE
   PUBLIC diagonalization
   CONTAINS
!============================================================
!
!   diagonalisation par jacobi
!   le vecteur i est V(.,i)
!
!============================================================
!
      SUBROUTINE diagonalization(Mat,Eig,Vec,n,type_diag,sort,phase)
      USE mod_NumParameters
      IMPLICIT NONE

      integer, intent(in)             :: n
      real(kind=Rkind), intent(in)    :: Mat(n,n)
      real(kind=Rkind), intent(inout) :: Eig(n),Vec(n,n)

      integer, intent(in), optional   :: type_diag,sort
      logical, intent(in), optional   :: phase


      integer  :: type_diag_loc

      !real(kind=Rkind) :: trav(n),Mat_save(n,n)
      real(kind=Rkind), allocatable :: trav(:),Mat_save(:,:)

      integer          :: ierr
      integer          :: lwork
      real(kind=Rkind), allocatable :: work(:)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='diagonalization'
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       allocate(Mat_save(n,n))
       Mat_save = Mat ! save mat


      IF (present(type_diag)) THEN
        type_diag_loc = type_diag
      ELSE
        type_diag_loc = 2
      END IF

      !when lapack is used and Rkind/= 8 (not a double), it switch to type_diag=2
      IF (Rkind /= 8 .AND. type_diag_loc == 3) type_diag_loc = 2

      SELECT CASE (type_diag_loc)
      CASE(1)
        CALL jacobi2(Mat_save,n,Eig,Vec)
      CASE(2)
        allocate(trav(n))
        CALL tred2(n,n,Mat_save,Eig,trav,Vec)
        CALL tql2(n,n,Eig,trav,Vec,ierr)
        deallocate(trav)
!      CASE(3) ! lapack77
!
!#if __LAPACK == 1
!        lwork = 3*n-1
!        CALL alloc_NParray(work,(/ lwork /),'work',name_sub)
!        Vec(:,:) = Mat_save(:,:)
!        CALL DSYEV('V','U',n,Vec,n,Eig,work,lwork,ierr)
!        IF (debug) write(out_unitp,*)'ierr=',ierr
!        IF (ierr /= 0) THEN
!           write(out_unitp,*) ' ERROR in ',name_sub
!           write(out_unitp,*) ' DSYEV lapack subroutine has FAILED!'
!           STOP
!        END IF
!        CALL dealloc_NParray(work,'work',name_sub)
!#else
!        allocate(trav(n))
!        CALL tred2(n,n,Mat_save,Eig,trav,Vec)
!        CALL tql2(n,n,Eig,trav,Vec,ierr)
!        deallocate(trav)
!#endif
!

      CASE DEFAULT
        CALL jacobi2(Mat_save,n,Eig,Vec)
      END SELECT

  IF (present(sort)) THEN
      SELECT CASE (sort)
      CASE(1)
        CALL trie(n,Eig,Vec,n)
      CASE(-1)
        Eig = -Eig
        CALL trie(n,Eig,Vec,n)
        Eig = -Eig
      CASE(2)
        CALL trie_abs(n,Eig,Vec,n)
      CASE DEFAULT ! no sort
        CONTINUE
      END SELECT
  ELSE
    CALL trie(n,Eig,Vec,n)
  END IF

  IF (present(phase)) THEN
      IF (phase) CALL Unique_phase(n,Vec,n)
  ELSE
    CALL Unique_phase(n,Vec,n)
  END IF

  deallocate(Mat_save)


      END SUBROUTINE diagonalization

      SUBROUTINE JACOBI(A,N,D,V,B,Z,max_N)
      USE mod_NumParameters
      IMPLICIT NONE

      integer    max_it
      parameter (max_it=500)

      integer   N,max_N
      real(kind=Rkind)    A(max_N,max_N),B(max_N),Z(max_N)
      real(kind=Rkind)    V(max_N,max_N),D(max_N)


      real(kind=Rkind) h,t,g,sm,tresh,tau,s,theta,c

      integer i,j,iq,nrot,ip


      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=ZERO
11      CONTINUE
        V(IP,IP)=ONE
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=ZERO
13    CONTINUE
      NROT=0
      DO 24 I=1,max_it
        SM=ZERO
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+abs(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.ZERO)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2_Rkind*SM/N**2
        ELSE
          TRESH=ZERO
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=HUNDRED*abs(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))                &
               .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=ZERO
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=HALF*H/A(IP,IQ)
                T=ONE/(ABS(THETA)+sqrt(ONE+THETA**2))
                IF(THETA.LT.ZERO)T=-T
              ENDIF
              C=ONE/sqrt(ONE+T**2)
              S=T*C
              TAU=S/(ONE+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=ZERO
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=ZERO
23      CONTINUE
24    CONTINUE
      write(out_unitp,*) max_it,' iterations should never happen'
      STOP

      end subroutine JACOBI

      SUBROUTINE JACOBI2(A,N,D,V)
      USE mod_NumParameters
      IMPLICIT NONE

      integer            :: N
      real(kind=Rkind)   :: A(N,N),V(N,N),D(N)


      integer, parameter :: max_it = 500
      real(kind=Rkind)   :: B(N),Z(N)


      real(kind=Rkind)   :: h,t,g,sm,tresh,tau,s,theta,c

      integer            :: i,j,iq,nrot,ip

!     V(:,:) = Id(:,:)
      V(:,:) = ZERO
      DO IP=1,N
        V(IP,IP)=ONE
      END DO

!     initialization
      DO IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=ZERO
      END DO

      NROT=0
      DO I=1,max_it ! main loop

        ! SM value
        SM = ZERO
        DO IP=1,N-1
          DO IQ=IP+1,N
            SM = SM+abs(A(IP,IQ))
          END DO
        END DO
        IF(SM == ZERO)RETURN

        ! TRESH value
        IF(I < 4)THEN
          TRESH = TWOTENTHS*SM/N**2
        ELSE
          TRESH = ZERO
        ENDIF

        DO IP=1,N-1
          DO IQ=IP+1,N
            G = HUNDRED*abs(A(IP,IQ))
            IF ( I > 4 .AND. ABS(D(IP))+G == ABS(D(IP))                 &
               .AND. ABS(D(IQ))+G == ABS(D(IQ)) ) THEN
              A(IP,IQ)=ZERO
            ELSE IF ( ABS(A(IP,IQ)) > TRESH ) THEN
              H=D(IQ)-D(IP)
              IF ( ABS(H)+G == ABS(H) ) THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=HALF*H/A(IP,IQ)
                T=ONE/(ABS(THETA)+sqrt(ONE+THETA**2))
                IF ( THETA < ZERO) T=-T
              ENDIF
              C=ONE/sqrt(ONE+T**2)
              S=T*C
              TAU=S/(ONE+C)

              H=T*A(IP,IQ)

              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=ZERO
              DO J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
              END DO
              DO J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
              END DO
              DO J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
              END DO
              DO J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
              END DO
              NROT=NROT+1
            ENDIF
          END DO
        END DO

        DO IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=ZERO
        END DO

      END DO ! end main loop

      write(out_unitp,*) max_it,' iterations should never happen'
      STOP

      end subroutine JACOBI2
!
!============================================================
!
!   diagonalisation trigonalisation puis diagonalisation
!
!============================================================
!
!       call tred2(200,idim1,hint,val,tablo,vec)
!       call tql2(200,idim1,val,tablo,vec,ierr)
!       write(*,*)'ierr=',ierr

      SUBROUTINE TRED2(NM,N,A,D,E,Z)
      USE mod_NumParameters
      IMPLICIT NONE

      integer   N,NM
      real(kind=Rkind)    A(NM,N),Z(NM,N),D(N),E(N)

!     already in mod_NumParameters
!     real(kind=Rkind)    h,scale,hh,f,g,one,zero
      real(kind=Rkind)    h,scale,hh,f,g
      integer   ii,l,k,jp1,j,i

!     write(*,*)
!     do 2 i=1,n
!     write(*,*)(a(i,j),j=1,n)
!2    continue
!     DATA ZERO/0./,ONE/1./

      DO 100 I=1,N
      DO 100 J=1,I
  100 Z(I,J)=A(I,J)
!     FOR I=N STEP -1 UNTIL 2 DO --
      DO 300 II=2,N
      I=N+2-II
      L=I-1
      H=ZERO
      SCALE=ZERO
      IF(L.LT.2) GOTO 130
!     SCALE ROW (ALGOL TOL THEN NOT NEEDED)
      DO 120 K=1,L
  120 SCALE=SCALE+ABS(Z(I,K))
      IF(SCALE.NE.ZERO) GOTO 140
  130 E(I)=Z(I,L)
      GOTO 290
  140 DO 150 K=1,L
      Z(I,K)=Z(I,K)/SCALE
      H=H+Z(I,K)*Z(I,K)
  150 CONTINUE
      F=Z(I,L)
      G=-SIGN(SQRT(H),F)
      E(I)=SCALE*G
      H=H-F*G
      Z(I,L)=F-G
      F=ZERO
      DO 240 J=1,L
      Z(J,I)=Z(I,J)/(SCALE*H)
      G=ZERO
!     FORM ELEMENT OF A*U
      DO 180 K=1,J
  180 G=G+Z(J,K)*Z(I,K)
      JP1=J+1
      IF(L.LT.JP1) GOTO 220
      DO 200 K=JP1,L
  200 G=G+Z(K,J)*Z(I,K)
!     FORM ELEMENT OF PP
  220 E(J)=G/H
      F=F+E(J)*Z(I,J)
  240 CONTINUE
      HH=F/(H+H)
!     FORM REDUCED A
      DO 260 J=1,L
      F=Z(I,J)
      G=E(J)-HH*F
      E(J)=G
      DO 260 K=1,J
      Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  260 CONTINUE
      DO 280 K=1,L
  280 Z(I,K)=SCALE*Z(I,K)
  290 D(I) = H
  300 CONTINUE
  320 D(1)=ZERO
      E(1)=ZERO
!     ACCUMULATION OF TRANSFORMATION MATRICES
      DO 500 I=1,N
      L=I-1
      IF(D(I).EQ.ZERO) GOTO 380
      DO 360 J=1,L
      G=ZERO
      DO 340 K=1,L
  340 G=G+Z(I,K)*Z(K,J)
      DO 360 K=1,L
      Z(K,J)=Z(K,J)-G*Z(K,I)
  360 CONTINUE
  380 D(I)=Z(I,I)
      Z(I,I)=ONE
      IF(L.LT.1) GOTO 500
      DO 400 J=1,L
      Z(I,J)=ZERO
      Z(J,I)=ZERO
  400 CONTINUE
  500 CONTINUE
      RETURN
      end subroutine TRED2

      SUBROUTINE TQL2(NZ,N,D,E,Z,IERR)
      USE mod_NumParameters
      IMPLICIT NONE

      integer   NZ,N
      real(kind=Rkind)    MACHEP
      real(kind=Rkind)    Z(NZ,N),D(N),E(N)
!     MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!     THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!     MACHEP = 2.**(-47) FOR SINGLE PRECISION ARITHMETIC
!     BOTH ON CDC AND CRAY

!     already in mod_NumParameters
!     real(kind=Rkind) ZERO,ONE,TWO
!     DATA ZERO/0./,ONE/1./,TWO/2./

      real(kind=Rkind)  g,pp,h,r,c,s,f,b
      integer l1,j,m,mml,ii,k,l,ierr,i
!
!RAY  1
!     MACHEP=2.**(-47)
!IBM  real(kind=Rkind)
!     MACHEP=16.**(-13)
!      MACHEP=epsilon(ONE)
      MACHEP=tiny(ONE)
      !write(out_unitp,*) 'MACHEP',epsilon(ONE),tiny(ONE),MACHEP
!IBM  simple precision
!     MACHEP=16.**(-5)
      IERR=0
      DO 100 I=2,N
  100 E(I-1)=E(I)
      F=ZERO
      B=ZERO
      E(N)=ZERO
      DO 240 L=1,N
      J=0
      H=MACHEP*(ABS(D(L))+ABS(E(L)))
      IF(B.LT.H) B=H
!     LOOK FOR SMALL SUB-DIAGONAL ELEMENT
      DO 110 M=L,N
      IF(ABS(E(M)).LE.B) GOTO 120
!     E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!     THROUGH THE BOTTOM OF THE LOOP
  110 CONTINUE
  120 IF(M.EQ.L) GOTO 220
  130 IF(J.EQ.30) GOTO 1000
      J=J+1
!     FORM SHIFT
      L1=L+1
      G=D(L)
      PP=(D(L1)-G)/(TWO*E(L))
      R=SQRT(PP*PP+ONE)
      D(L)=E(L)/(PP+SIGN(R,PP))
      H=G-D(L)
      DO 140 I=L1,N
  140 D(I)=D(I)-H
      F=F+H
!     QL TRANSFORMATION
      PP=D(M)
      C=ONE
      S=ZERO
      MML=M-L
!     FOR I=M-1 STEP -1 UNTIL L DO --
      DO 200 II=1,MML
      I=M-II
      G=C*E(I)
      H=C*PP
      IF(ABS(PP).LT.ABS(E(I))) GOTO 150
      C=E(I)/PP
      R=SQRT(C*C+ONE)
      E(I+1)=S*PP*R
      S=C/R
      C=ONE/R
      GOTO 160
  150 C=PP/E(I)
      R=SQRT(C*C+ONE)
      E(I+1)=S*E(I)*R
      S=ONE/R
      C=C*S
  160 PP=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
!     FORM VECTOR
      DO 180 K=1,N
      H=Z(K,I+1)
      Z(K,I+1)=S*Z(K,I)+C*H
      Z(K,I)=C*Z(K,I)-S*H
  180 CONTINUE
  200 CONTINUE
      E(L)=S*PP
      D(L)=C*PP
      IF(ABS(E(L)).GT.B) GOTO 130
  220 D(L)=D(L)+F
  240 CONTINUE
!     ORDER EIGENVALUES AND EIGENVECTORS
!     DO 300 II=2,N
!     I=II-1
!     K=I
!     PP=D(I)
!     DO 260 J=II,N
!     IF(D(J).GE.PP) GOTO 260
!     K=J
!     PP=D(J)
! 260 CONTINUE
!     IF(K.EQ.I) GOTO 300
!     D(K)=D(I)
!     D(I)=PP
!     DO 280 J=1,N
!     P=Z(J,I)
!     Z(J,I)=Z(J,K)
!     Z(J,K)=P
! 280 CONTINUE
! 300 CONTINUE
      GOTO 310
!     SET ERROR -- NO CONVERGENCE TO AN
!     EIGENVALUE AFTER 30 ITERATIONS
 1000 IERR=L
  310 RETURN
      end subroutine TQL2



!
!============================================================
!
!   trie des vecteur dans l'ordre croissant
!   le vecteur i est psi(.,i)
!
!============================================================
!
      SUBROUTINE trie_tab(nb_niv,ene,max_niv)
      USE mod_NumParameters
      IMPLICIT NONE

      integer nb_niv,max_niv
      real(kind=Rkind) ene(max_niv)
      real(kind=Rkind) a

        integer i,j,k

      DO i=1,nb_niv
      DO j=i+1,nb_niv
       IF (ene(i) .GT. ene(j)) THEN
!             permutation
          a=ene(i)
          ene(i)=ene(j)
          ene(j)=a
        END IF
      END DO
      END DO

      end subroutine trie_tab
!
!============================================================
!
!   trie des vecteur dans l'ordre croissant
!   le vecteur i est psi(.,i)
!
!============================================================
!
      SUBROUTINE trie(nb_niv,ene,psi,max_niv)
      USE mod_NumParameters
      IMPLICIT NONE

      integer nb_niv,max_niv
      real(kind=Rkind) ene(max_niv),psi(max_niv,max_niv)
      real(kind=Rkind) a

        integer i,j,k

      DO i=1,nb_niv
      DO j=i+1,nb_niv
       IF (ene(i) .GT. ene(j)) THEN
!	      permutation
          a=ene(i)
          ene(i)=ene(j)
          ene(j)=a
          DO k=1,nb_niv
            a=psi(k,i)
            psi(k,i)=psi(k,j)
            psi(k,j)=a
          END DO
        END IF
      END DO
      END DO

      END SUBROUTINE trie
!
!============================================================
!
!   Change the phase of Vec(:,i) shuch its lmargest coef is positive
!
!============================================================
!
      SUBROUTINE Unique_phase(n,Vec,max_n)
      USE mod_NumParameters
      IMPLICIT NONE

      integer          :: n,max_n
      real(kind=Rkind) :: Vec(max_n,max_n)

      real(kind=Rkind) :: max_val
      integer          :: ind_max_val(max_n),nb_max_val

      integer          :: i,j,max_j

      DO i=1,n
        max_val        = maxval(abs(Vec(:,i)))
        DO j=1,n
          IF ( abs(abs(Vec(j,i))-max_val) < ONETENTH**6 ) THEN
            IF (Vec(j,i) < ZERO) THEN
              Vec(:,i) = -Vec(:,i)
              EXIT
            END IF
          END IF
        END DO
      END DO

      END SUBROUTINE Unique_phase
      SUBROUTINE Unique_phase2(n,Vec,max_n)
      USE mod_NumParameters
      IMPLICIT NONE

      integer          :: n,max_n
      real(kind=Rkind) :: Vec(max_n,max_n)


      integer          :: i

      DO i=1,n
        IF (Vec(1,i) < ZERO) Vec(:,i) = -Vec(:,i)
      END DO

      END SUBROUTINE Unique_phase2

      SUBROUTINE Unique_phase_old(n,Vec,max_n)
      USE mod_NumParameters
      IMPLICIT NONE

      integer          :: n,max_n
      real(kind=Rkind) :: Vec(max_n,max_n)

      real(kind=Rkind) :: max_val
      integer          :: i,j,max_j

      DO i=1,n
        max_val = abs(Vec(1,i))
        max_j   = 1
        DO j=1,n
          IF ( abs(Vec(j,i)) > max_val ) THEN
            max_val = abs(Vec(j,i))
            max_j   = j
          END IF
        END DO
        IF (Vec(max_j,i) < ZERO) Vec(:,i) = -Vec(:,i)
      END DO

      END SUBROUTINE Unique_phase_old
      SUBROUTINE trie_abs(nb_niv,ene,psi,max_niv)
      USE mod_NumParameters
      IMPLICIT NONE

        integer nb_niv,max_niv
        real(kind=Rkind) ene(max_niv),psi(max_niv,max_niv)
        real(kind=Rkind) a

        integer i,j,k

        DO i=1,nb_niv
          DO j=i+1,nb_niv
            IF (abs(ene(i)) .GT. abs(ene(j))) THEN
!             permutation
              a=ene(i)
              ene(i)=ene(j)
              ene(j)=a
              DO k=1,nb_niv
                a=psi(k,i)
                psi(k,i)=psi(k,j)
                psi(k,j)=a
              END DO
            END IF
          END DO
        END DO

        RETURN
        end subroutine trie_abs
!
!============================================================
!
!      *******       Extension to complex symmetrical matrices of
!      *******       the <Tql2> Eispack routine, implemented by
!      *******       Claude Leforestier.
!
!
!     on input-
!        n is the order of the matrix,
!        d contains the diagonal elements of the input matrix,
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary,
!        z should have been initialized to the identity matrix.
!
!      on output-
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1,
!        e has been destroyed,
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues,
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after nbiter iterations.
!
!============================================================
!
      Subroutine cTql2(nZ,n,D,E,Z,ierr)
      USE mod_NumParameters
      IMPLICIT NONE

      Integer :: i,j,k,l,m,n,ii,Nbiter,nm,mml,nZ,ierr
      Data Nbiter/60/
      complex(kind=Rkind) :: D(n),E(n),Z(nZ,n)
      real (kind=Rkind) :: machep,norm1,norm2,rsign
      complex(kind=Rkind) :: b,c,f,g,p,r,s
!     complex(kind=Rkind) :: b,c,f,g,p,r,s,zero,one,two
!     Data zero/(0.,0.)/,one/(1.,0.)/,two/(2.,0.)/,Nbiter/60/


      machep=10.d-16
      ierr=0
!     initialize z to e one.
      do i=2,n
         e(i-1)=e(i)
      enddo
      e(n)=zero
      do 240 l=1,n
      j=0
!     ********** look for small sub-diagonal element **********
  105 do 110 m=l,n
      if(m.eq.n) goto 120
      norm1=abs(real(e(m),kind=Rkind))+abs(aimag(e(m)))
      norm2=abs(real(d(m),kind=Rkind))+abs(aimag(d(m))) +               &
            abs(real(d(m+1),kind=Rkind))+abs(aimag(d(m+1)))
      if(norm1.le.machep*norm2) goto 120
  110 continue
  120 p=d(l)
      if(m.eq.l) goto 240
      if(j.eq.nbiter) goto 1000
      j=j+1
!     ********** form shift **********
      g=(d(l+1)-p)/(TWO*e(l))
      r=sqrt(g*g+ONE)
      rsign=1.
      if(real(g,kind=Rkind).lt.0.) rsign=-ONE
      g=d(m)-p+e(l)/(g+rsign*r)
      s=ONE
      c=ONE
      p=ZERO
      mml=m-l
!     ********** for i=m-1 step -1 until l do -- **********
      do 200 ii=1,mml
      i=m-ii
      f=s*e(i)
      b=c*e(i)
      norm1=f*conjg(f)
      norm2=g*conjg(g)
      if(norm1.lt.norm2) goto 150
      c=g/f
      r=sqrt(c*c+ONE)
      e(i+1)=f*r
      s=ONE/r
      c=c*s
      go to 160
  150 s=f/g
      r=sqrt(s*s+ONE)
      e(i+1)=g*r
      c=ONE/r
      s=s*c
  160 g=d(i+1)-p
      r=(d(i)-g)*s+TWO*c*b
      p=s*r
      d(i+1)=g+p
      g=c*r-b
!     ********** form vector **********
      do k=1,n
      f=z(k,i+1)
      z(k,i+1)=s*z(k,i)+c*f
      z(k,i)=c*z(k,i)-s*f
      enddo
  200 continue
      d(l)=d(l)-p
      e(l)=g
      e(m)=ZERO
      go to 105
  240 continue
!     ********** order eigenvalues and eigenvectors **********
      do 300 ii=2,n
      i=ii-1
      k=i
      p=d(i)
      do 260 j=ii,n
      if(real(d(j),kind=Rkind).ge.real(p,kind=Rkind)) goto 260
      k=j
      p=d(j)
  260 continue
      if(k.eq.i) goto 300
      d(k)=d(i)
      d(i)=p
      do m=1,n
         p=z(m,i)
         z(m,i)=z(m,k)
         z(m,k)=p
      enddo
  300 continue
      goto 1001
!     ********** set error -- no convergence to an
!                eigenvalue after Nbiter iterations **********
 1000 write (out_unitp,1010) l
 1010 format(//10x,'$$$ <Cmtql2> return code :',i5,' $$$')
      stop
 1001 Return
      end subroutine cTql2

      Subroutine cTred2(nm,n,A,d,e,Z)
      USE mod_NumParameters
      IMPLICIT NONE

      INTEGER I,II,J,JP1,K,L,N,NM
      complex(kind=Rkind) A(nm,n),Z(nm,n),d(n),e(n)                     &
                ,f,g,h,hh
      real (kind=Rkind) :: scale,rsign
!     real (kind=Rkind) :: one,scale,rsign,zero
!     Data zero/0./,one/1./

      do i=1,n
         do j=1,i
            Z(i,j)=A(i,j)
         enddo
      enddo
!
!     FOR I=N STEP -1 UNTIL 2 DO --
      do ii=2,n
         I=N+2-II
         L=I-1
         H=ZERO
         SCALE=ZERO
         IF(L.LT.2) GOTO 130
!     SCALE ROW (ALGOL TOL THEN NOT NEEDED)
         DO 120 K=1,L
  120    SCALE=SCALE+ABS(Z(I,K))
         IF(SCALE.NE.ZERO) GOTO 140
  130    E(I)=Z(I,L)
         GOTO 290
  140    DO 150 K=1,L
         Z(I,K)=Z(I,K)/SCALE
         H=H+Z(I,K)*Z(I,K)
  150 CONTINUE
         F=Z(I,L)
         rsign=ONE
         if( real(F,kind=Rkind).lt.0.) rsign=-ONE
         G=-rsign*SQRT(H)
         E(I)=SCALE*G
         H=H-F*G
         Z(I,L)=F-G
         F=ZERO
         DO 240 J=1,L
         Z(J,I)=Z(I,J)/(SCALE*H)
         G=ZERO
!     FORM ELEMENT OF A*U
         DO 180 K=1,J
  180    G=G+Z(J,K)*Z(I,K)
         JP1=J+1
         IF(L.LT.JP1) GOTO 220
         DO 200 K=JP1,L
  200    G=G+Z(K,J)*Z(I,K)
!     FORM ELEMENT OF PP
  220    E(J)=G/H
         F=F+E(J)*Z(I,J)
  240    CONTINUE
         HH=F/(H+H)
!     FORM REDUCED A
         DO 260 J=1,L
         F=Z(I,J)
         G=E(J)-HH*F
         E(J)=G
         DO 260 K=1,J
         Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  260 CONTINUE
         DO 280 K=1,L
  280    Z(I,K)=SCALE*Z(I,K)
  290    D(I) = H
      enddo
  320 D(1)=ZERO
      E(1)=ZERO
!     ACCUMULATION OF TRANSFORMATION MATRICES
      DO 500 I=1,N
      L=I-1
      IF(abs(D(I)).EQ.ZERO) GOTO 380
      DO 360 J=1,L
      G=ZERO
      DO 340 K=1,L
  340 G=G+Z(I,K)*Z(K,J)
      DO 360 K=1,L
      Z(K,J)=Z(K,J)-G*Z(K,I)
  360 CONTINUE
  380 D(I)=Z(I,I)
      Z(I,I)=ONE
      IF(L.LT.1) GOTO 500
      DO 400 J=1,L
      Z(I,J)=ZERO
      Z(J,I)=ZERO
  400 CONTINUE
  500 CONTINUE
      RETURN
      end subroutine cTred2

END MODULE mod_diago
