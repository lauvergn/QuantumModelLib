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
!   Driver for the diagonalization
!      Default: tred2+tql2 (type_diag=2)
!            Other possibilities: Jacobi (type_diag=1) or Lapack (type_diag=3)
!            Rk: Lapack diago is not possible
!      Sort: the eigenvalues/eigenvectors:
!            sort=1:  ascending (default)
!            sort=-1: descending
!            sort=2:  ascending on the absolute eigenvalues
!     phase:
!============================================================
!
  SUBROUTINE diagonalization(Mat,REig,Vec,n,type_diag,sort,phase,IEig)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64,int32
    USE mod_NumParameters
    IMPLICIT NONE

    integer,          intent(in)              :: n
    real(kind=Rkind), intent(in)              :: Mat(n,n)
    real(kind=Rkind), intent(inout)           :: REig(n),Vec(n,n)
    real(kind=Rkind), intent(inout), optional :: IEig(n)

    integer,          intent(in),    optional :: type_diag,sort
    logical,          intent(in),    optional :: phase


    !local variables
    integer  :: type_diag_loc
    real(kind=Rkind), allocatable :: trav(:),Mat_save(:,:)
    integer              :: type_diag_default = 2 ! tred+tql

    !for lapack
    integer              :: i
    integer              ::    lwork ,lda ,ldvr ,ierr
    integer(kind=int32)  :: n4,lwork4,lda4,ldvr4,ierr4
    real(kind=Rkind), allocatable :: work(:)
    real(kind=Rkind), allocatable :: IEig_loc(:)

    real(kind=Rkind) :: dummy(1,1)


    !                                    Jacobi tred+tql DSYEV  DGEEV
    integer, parameter :: list_type(7) = [1,    2,202,   3,377, 4,477]


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='diagonalization'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    IF (present(type_diag)) THEN
      type_diag_loc = type_diag
    ELSE
      type_diag_loc = type_diag_default
    END IF

    !when lapack is used and Rkind /= real64 (not a double)
    IF (Rkind /= real64 .AND. type_diag_loc == 3) type_diag_loc = type_diag_default

#if __LAPACK != 1
    IF (count([3,377,395] == type_diag_loc) == 1) type_diag_loc = type_diag_default
    IF (count([4,477] == type_diag_loc) == 1) THEN
      !type_diag_loc = 0
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The diagonalization of non-symmetric needs LAPACK.'
      write(out_unitp,*) '  Try to link LAPACK with the code (use LAPACK=1 in the makfile).'
      write(out_unitp,*) '   type_diag:      ',type_diag_loc
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in diagonalization: Problem with non-symmetric matrix.'
    END IF

#endif

    IF (count(list_type == type_diag_loc) == 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' type_diag is out-of-range.'
      write(out_unitp,*) '   type_diag:      ',type_diag_loc
      write(out_unitp,*) '   Possible values:',list_type(:)
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in diagonalization: type_diag is out-of-range.'
    END IF


    SELECT CASE (type_diag_loc)
    CASE(1) ! jacobi
      IF (debug) write(out_unitp,*) 'Jacobi (symmetric)'
      allocate(Mat_save(n,n))
      Mat_save = Mat ! save mat

      CALL jacobi2(Mat_save,n,REig,Vec)

      deallocate(Mat_save)
    CASE (2) ! tred+tql
      IF (debug) write(out_unitp,*) 'tred+tql, new version (symmetric)'
      allocate(trav(n))

      Vec = Mat
      CALL TRED2_EISPACK(Vec,n,n,REig,trav)
      CALL TQLI_EISPACK(REig,trav,n,n,Vec)

      deallocate(trav)
    CASE(202) ! do not use anymore (old tred+tql)
      IF (debug) write(out_unitp,*) 'tred+tql, old version (symmetric)'

      allocate(Mat_save(n,n))
      Mat_save = Mat ! save mat
      allocate(trav(n))

      CALL tred2(n,n,Mat_save,REig,trav,Vec)
      CALL tql2(n,n,REig,trav,Vec,ierr)



      deallocate(Mat_save)
      deallocate(trav)
    CASE(3,377) ! lapack77
      IF (debug) write(out_unitp,*) 'lapack77: DSYEV (symmetric)'

#if __LAPACK == 1
      lwork = 3*n-1
      allocate(work(lwork))
      Vec(:,:) = Mat(:,:)

      ! lapack subroutines need integer (kind=4 or int32), therefore, we add a conversion, otherwise
      ! it fails when integers (kind=8 or int64) are used (at the compilation).
      n4     = int(n,kind=int32)
      lwork4 = int(lwork,kind=int32)
      CALL DSYEV('V','U',n4,Vec,n4,REig,work,lwork4,ierr4)

      IF (debug) write(out_unitp,*) 'ierr=',ierr4
      flush(out_unitp)

      IF (ierr4 /= 0_int32) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' DSYEV lapack subroutine has FAILED!'
         STOP
      END IF


      deallocate(work)
#else
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  LAPACK is not linked (LAPACK=0 in the makfile).'
      write(out_unitp,*) '  The program should not reach the LAPACK case.'
      write(out_unitp,*) '  => Probabely, wrong type_diag_default.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in diagonalization: LAPACK case impossible'
#endif
!      CASE(395) ! lapack95
!        IF (debug) write(out_unitp,*) 'lapack95: LA_SYEVD'
!        flush(out_unitp)
!        Vec(:,:) = Mat
!        CALL LA_SYEVD(Vec,Eig)

    CASE(4,477) ! lapack77 (non-symmetric)
#if __LAPACK == 1
      IF (debug) write(out_unitp,*) 'lapack77: DGEEV (non-symmetric)'
      flush(out_unitp)

      allocate(Mat_save(n,n))
      Mat_save = Mat ! save mat


      lwork = (2+64)*n
      ldvr  = n
      lda   = n
      allocate(work(lwork))


      n4     = int(n,kind=int32)
      lwork4 = int(lwork,kind=int32)
      lda4   = int(lda,kind=int32)
      ldvr4  = int(ldvr,kind=int32)

      IF (present(IEig)) THEN
        CALL DGEEV('N','V',n4,Mat_save,lda4,REig,IEig,dummy,              &
                   int(1,kind=int32),Vec,ldvr4,work,lwork4,ierr4)
        IF (debug) write(out_unitp,*)'ierr=',ierr4
        IF (ierr4 /= 0_int32) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' DGEEV lapack subroutine has FAILED!'
           STOP
        END IF

        IF (debug) THEN
          DO i=1,n
            write(out_unitp,*) 'Eigenvalue(', i, ') = ', REig(i),'+I ',IEig(i)
          END DO
        END IF
      ELSE
        allocate(IEig_loc(n))

        CALL DGEEV('N','V',n4,Mat_save,lda4,REig,IEig_loc,dummy,        &
                   int(1,kind=int32),Vec,ldvr4,work,lwork4,ierr4)
        IF (debug) write(out_unitp,*)'ierr=',ierr4
        IF (ierr4 /= 0_int32) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' DGEEV lapack subroutine has FAILED!'
           STOP
        END IF

        DO i=1,n
          write(out_unitp,*) 'Eigenvalue(', i, ') = ', REig(i),'+I ',IEig_loc(i)
        END DO

        deallocate(IEig_loc)
      END IF

      deallocate(work)
      deallocate(Mat_save)
#else
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  LAPACK is not linked (LAPACK=0 in the makfile).'
      write(out_unitp,*) '  The program should not reach the LAPACK case.'
      write(out_unitp,*) '  => Probabely, wrong type_diag_default.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in diagonalization: LAPACK case impossible'
#endif

    CASE DEFAULT
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The default CASE is not defined.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in diagonalization: default case impossible'
    END SELECT


    IF (present(sort)) THEN
        SELECT CASE (sort)
        CASE(1)
          CALL trie(n,REig,Vec,n)
          CALL rota_denerated(REig,Vec,n)
        CASE(-1)
          REig = -REig
          CALL trie(n,REig,Vec,n)
          REig = -REig
          CALL rota_denerated(REig,Vec,n)
        CASE(2)
          CALL trie_abs(n,REig,Vec,n)
        CASE DEFAULT ! no sort
          CONTINUE
        END SELECT
    ELSE
      CALL trie(n,REig,Vec,n)
      CALL rota_denerated(REig,Vec,n)
    END IF

    IF (present(phase)) THEN
        !IF (phase) CALL Unique_phase_old(n,Vec,n)
        IF (phase) CALL Unique_phase(Vec)

    ELSE
      CALL Unique_phase(Vec)
      !CALL Unique_phase_old(n,Vec,n)
    END IF

  END SUBROUTINE diagonalization

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

      integer          ::  N,NM
      real(kind=Rkind) ::  A(NM,N),Z(NM,N),D(N),E(N)

!     already in mod_NumParameters
!     real(kind=Rkind)    h,scale,hh,f,g,one,zero
      real(kind=Rkind)    h,scale,hh,f,g
      integer             ii,l,k,jp1,j,i

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

      GOTO 310
!     SET ERROR -- NO CONVERGENCE TO AN
!     EIGENVALUE AFTER 30 ITERATIONS
 1000 IERR=L
  310 RETURN
      end subroutine TQL2


      SUBROUTINE TRED2_EISPACK(A,N,NP,D,E)
      USE mod_NumParameters
      IMPLICIT NONE

      integer          :: N,NP
      real(kind=Rkind) :: A(NP,NP),D(NP),E(NP)

      !local variables
      integer          :: I,J,K,L
      real(kind=Rkind) :: F,G,H,HH,SCALE

      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=0.
          SCALE=0.
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE.EQ.0.)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=0.
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=0.
      E(1)=0.
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.)THEN
          DO 21 J=1,L
            G=0.
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.
            A(J,I)=0.
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
      END SUBROUTINE TRED2_EISPACK

      SUBROUTINE TQLI_EISPACK(D,E,N,NP,Z)
      USE mod_NumParameters
      IMPLICIT NONE

      integer          :: N,NP
      real(kind=Rkind) :: D(NP),E(NP),Z(NP,NP)

      !local variables
      integer          :: I,K,L,M,ITER
      real(kind=Rkind) :: G,R,S,C,P,F,B,DD

      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30) STOP 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END SUBROUTINE TQLI_EISPACK

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

        end subroutine trie_abs
!
!============================================================
!
!   Change the phase of Vec(:,i) shuch its largest coeficient is positive
!
!============================================================
!
      SUBROUTINE Unique_phase(Vec)
      USE mod_NumParameters
      IMPLICIT NONE

      real(kind=Rkind), intent(inout) :: Vec(:,:)

      integer          :: i,jloc

      DO i=lbound(Vec,dim=2),ubound(Vec,dim=2)
        jloc           = maxloc(abs(Vec(:,i)),dim=1)
        IF (abs(Vec(jloc,i)) < ONETENTH**6 ) CYCLE
        IF (Vec(jloc,i) < ZERO) Vec(:,i) = -Vec(:,i)
      END DO

      END SUBROUTINE Unique_phase
      SUBROUTINE Unique_phase_old(n,Vec,max_n)
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

      END SUBROUTINE Unique_phase_old
!=====================================================================
!
!   c_new(:,i) =  cos(th) c(:,i) + sin(th) c(:,j)
!   c_new(:,j) = -sin(th) c(:,j) + cos(th) c(:,j)
!
!    the angle is obtained such ...
!
!      it is working only if 2 vectors are degenerated !!!!
!
!=====================================================================
      SUBROUTINE rota_denerated(v,c,n)
      USE mod_NumParameters
      IMPLICIT NONE


      integer       :: n
      real (kind=Rkind), intent(in)    :: v(n)
      real (kind=Rkind), intent(inout) :: c(n,n)



      real (kind=Rkind) :: cd(n)
      integer           :: i,j,k,kloc
      real (kind=Rkind) :: ai,aj,norm,cc,ss

      real (kind=Rkind), parameter :: epsi = ONETENTH**10

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING rota_denerated'
      write(out_unitp,*) 'v',v
      !write(out_unitp,*) 'c',c
      END IF
!---------------------------------------------------------------------
      DO i=1,n-1

        j = i+1
        IF ( abs(v(i)-v(j)) < epsi) THEN
          !write(6,*) 'i,j',i,j
          !write(6,*) 'vec i',c(:,i)
          !write(6,*) 'vec j',c(:,j)
          cd(:) = c(:,i)**2+c(:,j)**2


          kloc = maxloc(cd(:),dim=1)

          cc   =  c(kloc,i)
          ss   = -c(kloc,j)
          !write(6,*) i,j,'cos sin',kloc,cc,ss
          norm = sqrt(cc*cc+ss*ss)
          cc   = cc/norm
          ss   = ss/norm
          !write(6,*) i,j,'cos sin',cc,ss

          DO k=1,n
           ai = c(k,i)
           aj = c(k,j)

           c(k,i) =  cc * ai + ss * aj
           c(k,j) = -ss * ai + cc * aj

          END DO

        END IF
      END DO



!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'new c',c
      write(out_unitp,*) 'END rota_denerated'
      END IF
!---------------------------------------------------------------------

      end subroutine rota_denerated
      SUBROUTINE rota_denerated_old(v,c,n)
      USE mod_NumParameters
      IMPLICIT NONE


      integer       :: n
      real (kind=Rkind), intent(in)    :: v(n)
      real (kind=Rkind), intent(inout) :: c(n,n)




      integer       :: i,j,k,kloc
      real (kind=Rkind) :: ai,aj,norm,cc,ss

      real (kind=Rkind), parameter :: epsi = ONETENTH**10

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING rota_denerated_old'
      write(out_unitp,*) 'v',v
      !write(out_unitp,*) 'c',c
      END IF
!---------------------------------------------------------------------
      DO i=1,n-1

        j = i+1
        IF ( abs(v(i)-v(j)) < epsi) THEN
          !write(6,*) 'i,j',i,j
          !write(6,*) 'vec i',c(:,i)
          !write(6,*) 'vec j',c(:,j)

          kloc = maxloc(abs(c(:,i)),dim=1)

          cc   =  c(kloc,i)
          ss   = -c(kloc,j)
          !write(6,*) i,j,'cos sin',kloc,cc,ss
          norm = sqrt(cc*cc+ss*ss)
          cc   = cc/norm
          ss   = ss/norm
          !write(6,*) i,j,'cos sin',cc,ss

          DO k=1,n
           ai = c(k,i)
           aj = c(k,j)

           c(k,i) =  cc * ai + ss * aj
           c(k,j) = -ss * ai + cc * aj

          END DO

        END IF
      END DO



!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'new c',c
      write(out_unitp,*) 'END rota_denerated_old'
      END IF
!---------------------------------------------------------------------

      end subroutine rota_denerated_old

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
