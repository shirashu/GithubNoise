*
      SUBROUTINE GLU2( A, N, N1, EPS, WK1, WK2, IP, IER )
*
*        GLU2
*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1
*
*               SOLVES SIMULTANEOUS LINEAR EQUATIONS
*               BY GAUSSIAN ELIMINATION METHOD.
*
*        INPUT - -
*             A(N1,N)  R *8  : 2-DIM. ARRAY CONTAINING THE COEFFICIENTS.
*             N        I *4  : ORDER OF MATRIX.
*             N1       I *4  : SIZE OF ARRAY A.
*             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE
*                              MATRIX. ( STANDARD VALUE 3.52D-15 )
*        OUTPUT - -
*             A(N1,N)        : RESULT OF GAUSSIAN ELIMINATION.
*             IP(N1)   I *4  : PIVOT NUMBER.
*             IER      I *4  : = 0,  FOR NORMAL EXECUTION.
*                              = 1,  FOR SINGULAR MATRIX.
*                              = 3,  FOR INVALID ARGUEMENT.
*        WORKING  -
*             WK1(N), WK2(N)
*                            : 1-DIM. ARRAY.
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N1,*), IP(*), WK1(*), WK2(*)
*             LEFT HAND SIDE
      IF( EPS.LT.0.0D0 )  EPS = 3.52D-15
      IF( ( 2*((N+1)/2).GT.N1 ).OR.( N.LE.0 ) )  THEN
         IER = 3
         WRITE(*,*) '  (SUBR. GLU2)  INVALID ARGUMENT.  N1, N =', N1, N
         RETURN
      END IF
*
      IER = 0
      IF( MOD(N,2).NE.0 )  THEN
         DO 10 J = 1, N
   10      A(N+1,J) = 0.0D0
         DO 20 I = 1, N
   20      A(I,N+1) = 0.0D0
         A(N+1,N+1) = 1.0D0
      END IF
*
      DO 100 K = 1, N, 2
*             FIND MAXIMUM ELEMENT IN THE K-TH COLUMN.
        K1 = K+1
        IPK = K
        AMAX = ABS(A(K,K))
        DO 110 I = K+1, N
          AIK = ABS(A(I,K))
          IF( AIK.GT.AMAX )  THEN
             IPK = I
             AMAX = AIK
          END IF
  110   CONTINUE
        IP(K) = IPK
*             EXCHANGE K-TH AND K+1'TH COLUMN FOR K-TH ELIMINATION.
        IF( AMAX.GT.EPS )  THEN
           IF( IPK.NE.K )  THEN
              W = A(IPK,K)
              A(IPK,K) = A(K,K)
              A(K,K) = W
              W = A(IPK,K1)
              A(IPK,K1) = A(K,K1)
              A(K,K1) = W
           END IF
*
           T0 = A(K,K)
           T1 = A(K,K1)
           DO 120 I = K+1, N
             A(I,K) = -A(I,K)/T0
  120        A(I,K1) = A(I,K1)+A(I,K)*T1
*             MATRIX IS SINGULAR.
        ELSE
           IER = 1
           IP(K) = K
           DO 140 I = K+1, N
  140        A(I,K) = 0.0D0
           WRITE(*,*) '  (SUBR. GLU2)  MATRIX IS SINGULAR AT K =', K
           RETURN
        END IF
*             FIND MAXIMUM ELEMENT IN THE K+1'TH COLUMN.
        AMAX = ABS(A(K1,K1))
        IPK1 = K1
        DO 210 I = K1+1, N
          AIK = ABS(A(I,K1))
          IF( AIK.GT.AMAX )  THEN
             IPK1 = I
             AMAX = AIK
          END IF
  210   CONTINUE
        IP(K1) = IPK1
*             EXCHANGE K-TH AND K+1'TH COLUMN FOR K+1'TH ELIMINATION.
        IF( AMAX.GT.EPS )  THEN
           IF( IPK1.NE.K1 )  THEN
              W = A(IPK1,K)
              A(IPK1,K) = A(K1,K)
              A(K1,K) = W
              W = A(IPK1,K1)
              A(IPK1,K1) = A(K1,K1)
              A(K1,K1) = W
           END IF
*             COMPUTE ALFA
           DO 220 I = K1+1, N
  220        A(I,K1) = -A(I,K1)/A(K1,K1)
*             MATRIX IS SINGULAR.
        ELSE
           IER = 1
           IP(K1) = K1
           DO 230 I = K1+1, N
  230        A(I,K1) = 0.0D0
           WRITE(*,*) '  (SUBR. GLU2)  MATRIX IS SINGULAR AT K =', K1
           RETURN
        END IF
*
        DO 290 I = K1+1, N
          WK1(I) = A(I,K)
  290     WK2(I) = A(I,K1)
*
        DO 300 J = K1+1, N
*             EXCHANGE FOR K-TH ELIMINATION.
          IF( IPK.NE.K )  THEN
             W = A(IPK,J)
             A(IPK,J) = A(K,J)
             A(K,J) = W
          END IF
*             EXCHANGE FOR K+1'TH ELIMINATION.
          IF( IPK1.NE.K1 )  THEN
             W = A(IPK1,J)
             A(IPK1,J) = A(K1,J)
             A(K1,J) = W
          END IF
*
          A(K1,J) = A(K1,J)+A(K1,K)*A(K,J)
*             GAUSSIAN ELIMINATION FOR K-TH AND K+1'TH PROCESS.
          T = A(K,J)
          T1 = A(K1,J)
          DO 310 I = K1+1, N
  310       A(I,J) = A(I,J)+WK1(I)*T+WK2(I)*T1
  300   CONTINUE
  100 CONTINUE
      RETURN
      END
*
      SUBROUTINE GLU4( A, N, N1, EPS, WK1, WK2, IP, IER )
*
*        GLU4
*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1
*
*               SOLVES SIMULTANEOUS LINEAR EQUATIONS
*               BY GAUSSIAN ELIMINATION METHOD.
*
*        INPUT - -
*             A(N1,N)  R *8  : 2-DIM. ARRAY CONTAINING THE COEFFICIENTS.
*             N        I *4  : ORDER OF MATRIX.
*             N1       I *4  : SIZE OF ARRAY A.
*             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE
*                              MATRIX. ( STANDARD VALUE 3.52D-15 )
*        OUTPUT - -
*             A(N1,N)        : RESULT OF GAUSSIAN ELIMINATION.
*             IP(N1)   I *4  : PIVOT NUMBER.
*             IER      I *4  : = 0,  FOR NORMAL EXECUTION.
*                              = 1,  FOR SINGULAR MATRIX.
*                              = 3,  FOR INVALID ARGUEMENT.
*        WORKING  -
*             WK1(N), WK2(N)
*                      R *8  : 1-DIM. ARRAY.
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N1,*), IP(*), WK1(*), WK2(*)
*             LEFT HAND SIDE
      IF( EPS.LT.0.0D0 )  EPS = 3.52D-15
      IF( ( 2*((N+1)/2).GT.N1 ).OR.( N.LE.0 ) )  THEN
         IER = 3
         WRITE(*,*) '  (SUBR. GLU4)  INVALID ARGUMENT.  N1, N =', N1, N
         RETURN
      END IF
*
      IER = 0
      IF( MOD(N,2).NE.0 )  THEN
         DO 10 J = 1, N
   10      A(N+1,J) = 0.0D0
         DO 20 I = 1, N
   20      A(I,N+1) = 0.0D0
         A(N+1,N+1) = 1.0D0
      END IF
*
      DO 100 K = 1, N, 2
*             FIND MAXIMUM ELEMENT IN THE K-TH COLUMN.
        K1 = K+1
        AMAX = ABS(A(K,K))
        IPK = K
        DO 110 I = K+1, N
          AIK = ABS(A(I,K))
          IF( AIK.GT.AMAX )  THEN
             IPK = I
             AMAX = AIK
          END IF
  110   CONTINUE
        IP(K) = IPK
*             EXCHANGE K-TH AND K+1'TH COLUMN FOR K-TH ELIMINATION.
        IF( AMAX.GT.EPS )  THEN
           IF( IPK.NE.K )  THEN
              W = A(IPK,K)
              A(IPK,K) = A(K,K)
              A(K,K) = W
              W = A(IPK,K1)
              A(IPK,K1) = A(K,K1)
              A(K,K1) = W
           END IF
*
           T0 = A(K,K)
           T1 = A(K,K1)
           DO 120 I = K+1, N
             A(I,K) = -A(I,K)/T0
  120        A(I,K1) = A(I,K1)+A(I,K)*T1
*             MATRIX IS SINGULAR.
        ELSE
           IER = 1
           IP(K) = K
           DO 140 I = K+1, N
  140        A(I,K) = 0.0D0
           WRITE(*,*) '  (SUBR. GLU4)  MATRIX IS SINGULAR AT K =', K
           RETURN
        END IF
*             FIND MAXIMUM ELEMENT IN THE K+1'TH COLUMN.
        AMAX = ABS(A(K1,K1))
        IPK1 = K1
        DO 210 I = K1+1, N
          AIK = ABS(A(I,K1))
          IF( AIK.GT.AMAX )  THEN
             IPK1 = I
             AMAX = AIK
          END IF
  210   CONTINUE
        IP(K1) = IPK1
*             EXCHANGE K-TH AND K+1'TH COLUMN FOR K+1'TH ELIMINATION.
        IF( AMAX.GT.EPS )  THEN
           IF( IPK1.NE.K1 )  THEN
              W = A(IPK1,K)
              A(IPK1,K) = A(K1,K)
              A(K1,K) = W
              W = A(IPK1,K1)
              A(IPK1,K1) = A(K1,K1)
              A(K1,K1) = W
           END IF
*             COMPUTE ALFA
           DO 220 I = K1+1, N
  220        A(I,K1) = -A(I,K1)/A(K1,K1)
*             MATRIX IS SINGULAR.
        ELSE
           IER = 1
           IP(K1) = K1
           DO 230 I = K1+1, N
  230        A(I,K1) = 0.0D0
           WRITE(*,*) '  (SUBR. GLU4)  MATRIX IS SINGULAR AT K =', K1
           RETURN
        END IF
*
        DO 290 I = K1+1, N
          WK1(I) = A(I,K)
  290     WK2(I) = A(I,K1)
*
        DO 300 J = K1+1, N, 2
*             EXCHANGE FOR K-TH ELIMINATION.
          IF( IPK.NE.K )  THEN
             W = A(IPK,J)
             A(IPK,J) = A(K,J)
             A(K,J) = W
             W = A(IPK,J+1)
             A(IPK,J+1) = A(K,J+1)
             A(K,J+1) = W
          END IF
*             EXCHANGE FOR K+1'TH ELIMINATION.
          IF( IPK1.NE.K1 )  THEN
             W = A(IPK1,J)
             A(IPK1,J) = A(K1,J)
             A(K1,J) = W
             W = A(IPK1,J+1)
             A(IPK1,J+1) = A(K1,J+1)
             A(K1,J+1) = W
          END IF
*
          A(K1,J) = A(K1,J)+A(K1,K)*A(K,J)
          A(K1,J+1) = A(K1,J+1)+A(K1,K)*A(K,J+1)
*             GAUSSIAN ELIMINATION 2 COL. FOR K-TH AND K+1'TH PROCESS.
          T = A(K,J)
          T1 = A(K1,J)
          U = A(K,J+1)
          U1 = A(K1,J+1)
          DO 310 I = K1+1, N
            A(I,J) = A(I,J)+WK1(I)*T+WK2(I)*T1
  310       A(I,J+1) = A(I,J+1)+WK1(I)*U+WK2(I)*U1
  300      CONTINUE
  100 CONTINUE
      RETURN
      END
*
      SUBROUTINE GSLV4( A, N, N1, B, IP )
*
*        GSLV4
*               COPYRIGHT : H.HASEGAWA, OCT.  4 1991 V.1
*
*               SOLVES SIMULTANEOUS LINEAR EQUATIONS
*               BY GAUSSIAN ELIMINATION METHOD.
*
*        INPUT - -
*             A(N1,N)  R *8  : RESULT OF GAUSSIAN ELIMINATION.
*             N        I *4  : ORDER OF MATRIX.
*             N1       I *4  : SIZE OF ARRAY A.
*             B(N1)    R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT HAND
*                              SIDE VECTOR.
*             IP(N1)   I *4  : PIVOT NUMBER.
*        OUTPUT - -
*             B(N1)          : SOLUTION.
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N1,*), B(*), IP(*)
*             FORWARD ELIMINATION PROCESS
      DO 100 K = 1, N, 2
*             EXCHANGE FOR K+1'TH ELIMINATION.
        IF( IP(K).NE.K ) THEN
           W = B(IP(K))
           B(IP(K)) = B(K)
           B(K) = W
        END IF
*             EXCHANGE FOR K+1'TH ELIMINATION.
        K1 = K+1
        IF( IP(K1).NE.K1 ) THEN
           W = B(IP(K1))
           B(IP(K1)) = B(K1)
           B(K1) = W
        END IF
*
        B(K1) = B(K1)+A(K1,K)*B(K)
*             GAUSSIAN ELIMINATION FOR K-TH AND K+1'TH PROCESS.
        T = B(K)
        T1 = B(K1)
        DO 110 I = K1+1, N
  110     B(I) = B(I)+A(I,K)*T+A(I,K1)*T1
  100 CONTINUE
*             BACKWARD SUBSTITUTION PROCESS
      IF( MOD(N,2).EQ.0 )  THEN
         NEND = N
         B(N) = B(N)/A(N,N)
      ELSE
         NEND = N+1
         B(NEND) = 0.0D0
      END IF
*
      B(NEND-1) = (B(NEND-1)-A(NEND-1,NEND)*B(NEND))/A(NEND-1,NEND-1)
      DO 200 K = NEND-2, 1, -2
        T1 = B(K+1)
        T2 = B(K+2)
        DO 210 I = 1, K
  210     B(I) = B(I)-A(I,K+1)*T1-A(I,K+2)*T2
        B(K) = B(K)/A(K,K)
        B(K-1) = (B(K-1)-A(K-1,K)*B(K))/A(K-1,K-1)
  200 CONTINUE
      RETURN
      END
