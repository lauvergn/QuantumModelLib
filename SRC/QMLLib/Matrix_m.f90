MODULE QMLLib_Matrix_m
!$ USE omp_lib
  IMPLICIT NONE

  PRIVATE

  PUBLIC inv_m1_TO_m2
  INTERFACE inv_m1_TO_m2
     MODULE PROCEDURE QML_inv_m1_TO_m2
     MODULE PROCEDURE QML_inv_m1_TO_m2_cplx
  END INTERFACE

  PUBLIC Linear_Sys
  INTERFACE Linear_Sys
     MODULE PROCEDURE QML_Linear_Sys
  END INTERFACE

  CONTAINS
!================================================================
!    inversion de la matrice m1 : m2=m1^-1
!================================================================
      SUBROUTINE QML_inv_m1_TO_m2(m1,m2,n,inv_type,epsi)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer          :: n
       real(kind=Rkind) :: m1(n,n)
       real(kind=Rkind) :: m2(n,n)
       integer          :: inv_type
       real(kind=Rkind) :: epsi

       integer          :: indx(n)
       real(kind=Rkind) :: trav(n),m1w(n,n)
       real(kind=Rkind) :: vv(n,n)
       real(kind=Rkind) :: b(n)

       real(kind=Rkind) :: wmax,wmin
       real(kind=Rkind) :: d
       integer          :: j



       CALL QML_mat_id(m2,n)
       m1w = m1

       SELECT CASE (inv_type)
       CASE (0) ! ludcmp ...
         CALL QML_ludcmp(m1w,n,trav,indx,d)
         DO j=1,n
           CALL QML_lubksb(m1w,n,indx,m2(:,j))
         END DO
       CASE (1) ! svd
         CALL QML_SVDCMP(m1w,n,n,trav,vv,n)
         ! Find maximum singular value
         !write(out_unitp,*) 'SVD : epsi',epsi
         !write(out_unitp,*) 'SVD : trav',trav

         wmax = maxval(trav(:))
         wmin = wmax * epsi
         !write(out_unitp,*) 'SVD : count non zero',count(trav >= wmin)
         ! Zero the "small" singular values
         WHERE (trav < WMIN) trav = ZERO

         DO j=1,n
           b(:) = m2(:,j)
           CALL QML_SVBKSB(m1w,trav,vv,n,n,b,m2(:,j),n)
         END DO
       CASE Default ! ludcmp ...
          CALL QML_ludcmp(m1w,n,trav,indx,d)
          DO j=1,n
            CALL QML_lubksb(m1w,n,indx,m2(:,j))
          END DO
       END SELECT

       END SUBROUTINE QML_inv_m1_TO_m2
!================================================================
!    inversion de la matrice m1 : m2=m1^-1
!================================================================
      SUBROUTINE QML_inv_m1_TO_m2_cplx(m1,m2,n,inv_type,epsi)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer          :: n
       complex(kind=Rkind) :: m1(n,n)
       complex(kind=Rkind) :: m2(n,n)
       integer          :: inv_type
       real(kind=Rkind) :: epsi

       integer          :: indx(n)
       complex(kind=Rkind) :: trav(n),m1w(n,n)
       complex(kind=Rkind) :: vv(n,n)
       complex(kind=Rkind) :: b(n)

       complex(kind=Rkind) :: wmax,wmin
       complex(kind=Rkind) :: d
       integer          :: j



       CALL QML_Cplx_mat_id(m2,n)
       m1w = m1

       SELECT CASE (inv_type)
       CASE (0) ! ludcmp ...
         CALL QML_ludcmp_cplx(m1w,n,trav,indx,d)
         DO j=1,n
           CALL QML_lubksb_cplx(m1w,n,indx,m2(:,j))
         END DO

       CASE (1) ! svd

          STOP 'SVD not yet in complex'

       CASE Default ! ludcmp ...
         CALL QML_ludcmp_cplx(m1w,n,trav,indx,d)
         DO j=1,n
           CALL QML_lubksb_cplx(m1w,n,indx,m2(:,j))
         END DO
       END SELECT

     END SUBROUTINE QML_inv_m1_TO_m2_cplx
!================================================================
!    Dertermniant of m1
!================================================================
      SUBROUTINE QML_Det_OF_m1(m1,det,n)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer          :: n
       real(kind=Rkind) :: m1(n,n)
       real(kind=Rkind) :: det

       integer          :: index(n)
       real(kind=Rkind) :: trav(n),m1w(n,n)

       real(kind=Rkind) :: d
       integer          :: j

       m1w = m1

       CALL QML_ludcmp(m1w,n,trav,index,d)

       det = d
       DO j=1,n
         det = det * m1w(j,j)
       END DO

       END SUBROUTINE QML_Det_OF_m1
!================================================================
!    Solve,x: a.x=b
!================================================================
      SUBROUTINE QML_Linear_Sys(a,b,x,n)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer          :: n
       real(kind=Rkind) :: a(n,n),vv(n,n)
       real(kind=Rkind) :: b(n),x(n)

       integer          :: indx(n)
       real(kind=Rkind) :: trav(n),d
       real(kind=Rkind) :: aa(n,n)
       real(kind=Rkind) :: wmax,wmin,epsi=ONETENTH**10
       integer          :: k

       logical          :: svd = .TRUE.
       !logical          :: svd = .FALSE.

       x  = b
       aa = a

       IF (svd) THEN
           ! une facon.... SVD
           CALL QML_SVDCMP(aa,n,n,trav,vv,n)
           ! Find maximum singular value
           wmax = maxval(trav(:))
           wmin = wmax * epsi
           ! Zero the "small" singular values
           DO k=1,n
              IF (trav(k) < WMIN) trav(k) = ZERO
           END DO

           CALL QML_SVBKSB(aa,trav,vv,n,n,b,x,n)


           !write(out_unitp,*) 'solve?',sum(abs(matmul(a,x)-b))
           !STOP
        ELSE
          ! une autre ...
          CALL QML_ludcmp(aa,n,trav,indx,d)
          CALL QML_lubksb(aa,n,indx,x)
          !IF (mpro) CALL CALL mprove(a,aa,n,indx,b,x)

        END IF

       END SUBROUTINE QML_Linear_Sys
!================================================================
!    ameliore la solution d un systeme d equations
!    par une iteration
!================================================================
      SUBROUTINE QML_mprove(A,ALUD,N,INDX,B,X)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

      integer          :: n
      real(kind=Rkind) :: a(n,n)
      real(kind=Rkind) :: alud(n,n)
      real(kind=Rkind) :: b(n),x(n),r(n)
      integer          :: indx(n)
      integer          :: i,j

      DO I=1,N
        R(I) = -B(I) + dot_product(A(I,:),X(:))
      END DO
      CALL QML_LUBKSB(ALUD,N,INDX,R)

      X(:) = X(:) - R(:)

      END SUBROUTINE QML_mprove

!
!================================================================
!    resolution de a*x=b apres la procedure ludcmp
!
!================================================================

      SUBROUTINE QML_lubksb(a,n,index,b)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer n
       real(kind=Rkind) a(n,n),b(n)
       integer index(n)
       real(kind=Rkind)  sum

       integer i,j,ii,ll

       ii=0
       DO 12 i=1,n
         ll=index(i)
         sum=b(ll)
         b(ll)=b(i)
         IF (II .NE. 0) THEN
            DO 11 j=ii,i-1
              sum=sum-a(i,j)*b(j)
 11         CONTINUE
         ELSE IF (sum .NE. ZERO) THEN
                ii=i
              ENDIF
         b(i)=sum
 12    CONTINUE
       DO 14 i=n,1,-1
         sum=b(i)
         DO 13 j=i+1,n
           sum=sum-a(i,j)*b(j)
 13      CONTINUE
         b(i)=sum/a(i,i)
 14    CONTINUE

       RETURN
       end subroutine QML_lubksb
!================================================================
!    decomposition de a=l*u (pour la resolution d un systeme d equations
!     l matrice triangulaire inferieur
!     u matrice triangulaire superieur
!
!    a l u matrices n*n
!
!================================================================

      SUBROUTINE QML_ludcmp(a,n,vv,index,d)
      USE QMLLib_NumParameters_m
      USE QMLLib_UtilLib_m
      IMPLICIT NONE

       integer n
       real(kind=Rkind)   tiny
       parameter (tiny=ONETENTH**20)
       real(kind=Rkind) a(n,n),vv(n)
       real(kind=Rkind) aamax,sum,dum,d
       integer index(n)

       integer i,j,k,imax

       d=ONE
       DO 12 i=1,n
        aamax=ZERO
        DO 11 j=1,n
          IF (abs(a(i,j)) .GT. aamax) aamax=abs(a(i,j))
 11     CONTINUE
        IF (aamax < tiny) STOP "matrice singuliere"
        vv(i)=ONE/aamax
 12    CONTINUE


       DO 19 j=1,n

        DO 14 i=1,j-1
         sum=a(i,j)
         DO 13 k=1,i-1
          sum=sum-a(i,k)*a(k,j)
 13      CONTINUE
         a(i,j)=sum
 14     CONTINUE

        aamax=ZERO
        imax=0
        DO 16 i=j,n
         sum=a(i,j)
         DO 15 k=1,j-1
          sum=sum-a(i,k)*a(k,j)
 15      CONTINUE
         a(i,j)=sum
         dum=vv(i)*abs(sum)
         IF (dum .GE. aamax) THEN
           imax=i
           aamax=dum
         ENDIF
 16     CONTINUE
        IF (imax ==0) THEN
          write(out_unitp,*) ' ERROR in ludcmp'
          write(out_unitp,*) ' imax = 0 !!!'
          write(out_unitp,*) ' matrix a:'
          CALL Write_RMat(a,out_unitp,4)
          STOP
        END IF

        IF (j .NE. imax) THEN
          DO 17 k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
 17       CONTINUE
          d=-d
          vv(imax)=vv(j)
        ENDIF

        index(j)=imax
        IF (a(j,j) .EQ. ZERO) a(j,j)=tiny
        IF (j .NE. n) THEN
          dum=ONE/a(j,j)
          DO 18 i=j+1,n
            a(i,j)=a(i,j)*dum
 18       CONTINUE
        ENDIF

 19    CONTINUE


       RETURN
     END SUBROUTINE QML_ludcmp

      SUBROUTINE QML_SVDCMP(A,M,N,W,V,max_n)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

      integer max_n,N,M
      real (kind=Rkind) :: A(max_n,max_n),V(max_n,max_n)
      real (kind=Rkind) :: W(max_n),RV1(max_n)
      real (kind=Rkind) :: G,SCALE,ANORM,S,F,H,C,Y,Z,X

      integer I,K,J,NM,JJ,L,ITS


      G=ZERO
      SCALE=ZERO
      ANORM=ZERO
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=ZERO
        S=ZERO
        SCALE=ZERO
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.ZERO) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(sqrt(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=ZERO
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=ZERO
        S=ZERO
        SCALE=ZERO
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.ZERO) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(sqrt(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=ZERO
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.ZERO) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=ZERO
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=ZERO
            V(J,I)=ZERO
31        CONTINUE
        ENDIF
        V(I,I)=ONE
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=ZERO
33        CONTINUE
        ENDIF
        IF (G.NE.ZERO) THEN
          G=ONE/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=ZERO
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=ZERO
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+ONE
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=ZERO
          S=ONE
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=sqrt(F*F+G*G)
              W(I)=H
              H=ONE/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.ZERO) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.50) STOP 'No convergence in 50 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
          G=sqrt(F*F+ONE)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=ONE
          S=ONE
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=sqrt(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 JJ=1,N
              X=V(JJ,J)
              Z=V(JJ,I)
              V(JJ,J)= (X*C)+(Z*S)
              V(JJ,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=sqrt(F*F+H*H)
            W(J)=Z
            IF (Z.NE.ZERO) THEN
              Z=ONE/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 JJ=1,M
              Y=A(JJ,J)
              Z=A(JJ,I)
              A(JJ,J)= (Y*C)+(Z*S)
              A(JJ,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=ZERO
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
    END SUBROUTINE QML_SVDCMP


      SUBROUTINE QML_SVBKSB(U,W,V,M,N,B,X,max_n)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

      integer max_n,M,N
      real (kind=Rkind) :: U(max_n,max_n),V(max_n,max_n)
      real (kind=Rkind) :: W(max_n),B(max_n),X(max_n),TMP(max_n)
      real (kind=Rkind) :: s
      integer I,J,JJ

      DO 12 J=1,N
        S=ZERO
        IF(W(J).NE.ZERO)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=ZERO
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN

    END SUBROUTINE QML_SVBKSB


!================================================================
!    inversion de la matrice a : c=1/a
!================================================================

      SUBROUTINE QML_inversion_cplx(c,a,trav,index,n)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer n
       complex(kind=Rkind) a(n,n),d
       complex(kind=Rkind) c(n,n)

       integer index(n)
       complex(kind=Rkind) trav(n)

       integer i,j

       DO i=1,n
         DO j=1,n
           c(i,j)=CZERO
         END DO
         c(i,i)=CONE
       END DO
       CALL QML_ludcmp_cplx(a,n,trav,index,d)

       DO j=1,n
         CALL QML_lubksb_cplx(a,n,index,c(1,j))
       END DO

       RETURN
       end subroutine QML_inversion_cplx
!================================================================
!    resolution de a*x=b apres la procedure ludcmp
!================================================================
  SUBROUTINE QML_Driver_LU_solve_cplx(a,n,LU_index,b,type_lu)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64,int32
  USE QMLLib_NumParameters_m
  IMPLICIT NONE

  integer,             intent(in)    :: n,type_lu
  complex(kind=Rkind), intent(inout) :: a(n,n),b(n)
  integer,             intent(in)    :: LU_index(n)

  integer               :: err,type_lu_loc
  integer, parameter    :: type_lu_default = 1
  integer(kind=int32)   :: n4,ierr4



    type_lu_loc = type_lu

    !when lapack is used and Rkind /= real64 (not a double)
    IF (Rkind /= real64 .AND. type_lu_loc == 3) type_lu_loc = type_lu_default

#if __LAPACK != 1
    IF ( type_lu_loc == 3) type_lu_loc = type_lu_default
#endif

    SELECT CASE (type_lu)
    CASE(1) ! ori
      CALL QML_lubksb_cplx(a,n,LU_index,b)
    CASE(3) ! lapack
#if __LAPACK == 1
      n4     = int(n,kind=int32)
      CALL ZGETRS('No transpose',n4,1,a,n4,LU_index,b,n4,ierr4)
      err = int(ierr4)
      IF (err /= 0) STOP 'LU Driver_LU_solve_cplx'
#else
      write(out_unitp,*) ' ERROR in Driver_LU_solve_cplx'
      write(out_unitp,*) '  LAPACK is not linked (LAPACK=0 in the makefile).'
      write(out_unitp,*) '  The program should not reach the LAPACK case.'
      write(out_unitp,*) '  => Probabely, wrong type_diag_default.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in Driver_LU_solve_cplx: LAPACK case impossible'
#endif
    CASE Default
      CALL QML_lubksb_cplx(a,n,LU_index,b)
    END SELECT

  END SUBROUTINE QML_Driver_LU_solve_cplx
  SUBROUTINE QML_Driver_LU_decomp_cplx(a,n,LU_index,d,type_lu)
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64,int32
  USE QMLLib_NumParameters_m
  IMPLICIT NONE

  integer,             intent(in)    :: n,type_lu
  complex(kind=Rkind), intent(inout) :: d,a(n,n)
  integer,             intent(in)    :: LU_index(n)

  integer               :: err,type_lu_loc
  integer, parameter    :: type_lu_default = 1
  integer(kind=int32)  :: n4,ierr4
  complex(kind=Rkind), allocatable :: work(:)



    type_lu_loc = type_lu

    !when lapack is used and Rkind /= real64 (not a double)
    IF (Rkind /= real64 .AND. type_lu_loc == 3) type_lu_loc = type_lu_default

#if __LAPACK != 1
    IF ( type_lu_loc == 3) type_lu_loc = type_lu_default
#endif

    SELECT CASE (type_lu)
    CASE(1) ! ori
      allocate(work(n))
      CALL QML_ludcmp_cplx(a,n,work,LU_index,d)
      deallocate(work)
    CASE(3) ! lapack
#if __LAPACK == 1
      n4     = int(n,kind=int32)
      CALL ZGETRF(n4,n4,a,n4,LU_index,ierr4)
      err = int(ierr4)
      IF (err /= 0) STOP 'Driver_LU_decomp_cplx'
#else
      write(out_unitp,*) ' ERROR in Driver_LU_decomp_cplx'
      write(out_unitp,*) '  LAPACK is not linked (LAPACK=0 in the makefile).'
      write(out_unitp,*) '  The program should not reach the LAPACK case.'
      write(out_unitp,*) '  => Probabely, wrong type_diag_default.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in Driver_LU_decomp_cplx: LAPACK case impossible'
#endif
    CASE Default
      allocate(work(n))
      CALL QML_ludcmp_cplx(a,n,work,LU_index,d)
      deallocate(work)
    END SELECT

  END SUBROUTINE QML_Driver_LU_decomp_cplx
  SUBROUTINE QML_lubksb_cplx(a,n,index,b)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer n
       complex(kind=Rkind) a(n,n),b(n)
       integer index(n)
       complex(kind=Rkind)  sum

       integer i,j,ii,ll

       ii=0
       DO 12 i=1,n
         ll=index(i)
         sum=b(ll)
         b(ll)=b(i)
         IF (II .NE. 0) THEN
            DO 11 j=ii,i-1
              sum=sum-a(i,j)*b(j)
 11         CONTINUE
         ELSE IF (abs(sum) .NE. ZERO) THEN
                ii=i
              ENDIF
         b(i)=sum
 12    CONTINUE
       DO 14 i=n,1,-1
         sum=b(i)
         DO 13 j=i+1,n
           sum=sum-a(i,j)*b(j)
 13      CONTINUE
         b(i)=sum/a(i,i)
 14    CONTINUE

       RETURN
     end subroutine QML_lubksb_cplx
!================================================================
!    decomposition de a=l*u (pour la resolution d un systeme d equations
!     l matrice triangulaire inferieur
!     u matrice triangulaire superieur
!
!    a l u matrices n*n
!
!================================================================

      SUBROUTINE QML_ludcmp_cplx(a,n,vv,index,d)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer n
       real(kind=Rkind)   tiny
       parameter (tiny=ONETENTH**20)
       complex(kind=Rkind) a(n,n),vv(n)
       complex(kind=Rkind) aamax,sum,dum,d
       integer index(n)

       integer i,j,k,imax

       d=CONE
       DO 12 i=1,n
        aamax=CZERO
        DO 11 j=1,n
          IF (abs(a(i,j)) .GT. abs(aamax)) aamax=cmplx(abs(a(i,j)),kind=Rkind)
 11     CONTINUE
        IF (abs(aamax) < tiny) STOP "matrice singuliere"
        vv(i)=CONE/aamax
 12    CONTINUE


       DO 19 j=1,n

        DO 14 i=1,j-1
         sum=a(i,j)
         DO 13 k=1,i-1
          sum=sum-a(i,k)*a(k,j)
 13      CONTINUE
         a(i,j)=sum
 14     CONTINUE

        aamax=CZERO
        imax = 0
        DO 16 i=j,n
         sum=a(i,j)
         DO 15 k=1,j-1
          sum=sum-a(i,k)*a(k,j)
 15      CONTINUE
         a(i,j)=sum
         dum=vv(i)*cmplx(abs(sum),kind=Rkind)
         IF (abs(dum) .GE. abs(aamax)) THEN
           imax=i
           aamax=dum
         ENDIF
 16     CONTINUE

        IF (j .NE. imax) THEN
          DO 17 k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
 17       CONTINUE
          d=-d
          vv(imax)=vv(j)
        ENDIF

        index(j)=imax
        IF (abs(a(j,j)) .EQ. ZERO) a(j,j)=cmplx(tiny,kind=Rkind)
        IF (j .NE. n) THEN
          dum=CONE/a(j,j)
          DO 18 i=j+1,n
            a(i,j)=a(i,j)*dum
 18       CONTINUE
        ENDIF

 19    CONTINUE


       RETURN
       end subroutine QML_ludcmp_cplx


!=====================================================================
!
! ++   A = Id =>  A(i,i)=ONE
!      A : square matrix
!
!=====================================================================
!
  SUBROUTINE QML_mat_id(a,n)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer          :: i,n
       real(kind=Rkind) :: a(n,n)

       a(:,:) = ZERO

       DO i=1,n
         a(i,i) = ONE
       END DO

  END SUBROUTINE QML_mat_id
  SUBROUTINE QML_Cplx_mat_id(a,n)
      USE QMLLib_NumParameters_m
      IMPLICIT NONE

       integer          :: i,n
       complex(kind=Rkind) :: a(n,n)

       a(:,:) = CZERO

       DO i=1,n
         a(i,i) = CONE
       END DO

  END SUBROUTINE QML_Cplx_mat_id

END MODULE QMLLib_Matrix_m
