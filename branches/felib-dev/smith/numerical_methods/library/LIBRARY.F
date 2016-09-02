C         Appendix 2
C         ----------
C         Alphabetic listing of subroutine library in FORTRAN 77
C         ------------------------------------------------------
C
C
      SUBROUTINE BILIN(DER,FUN,SAMP,ISAMP,I,J)
C
C      SHAPE FUNCTIONS AND DERIVATIVES
C      FOR QUADRILATERAL TRANSFORMATION
C
      REAL DER(2,4),FUN(4),SAMP(ISAMP,2)
      ETA = SAMP(I,1)
      XI = SAMP(J,1)
      ETAM = .25* (1.-ETA)
      ETAP = .25* (1.+ETA)
      XIM = .25* (1.-XI)
      XIP = .25* (1.+XI)
      FUN(1) = 4.*XIM*ETAM
      FUN(2) = 4.*XIM*ETAP
      FUN(3) = 4.*XIP*ETAP
      FUN(4) = 4.*XIP*ETAM
      DER(1,1) = -ETAM
      DER(1,2) = -ETAP
      DER(1,3) = ETAP
      DER(1,4) = ETAM
      DER(2,1) = -XIM
      DER(2,2) = XIM
      DER(2,3) = XIP
      DER(2,4) = -XIP
      RETURN
      END
C
C
C
      SUBROUTINE CHECK(X1,X0,TOL,ICON)
C
C        SETS ICON TO 0 IF RELATIVE CHANGE IN X
C        IS GREATER THAN TOL
C
      ICON = 1
      IF (ABS(X1-X0)/ABS(X1).GT.TOL) ICON = 0
      X0 = X1
      RETURN
      END
C
C
C
      SUBROUTINE CHECON(LOADS,OLDLDS,N,TOL,ICON)
C
C      SETS ICON TO ZERO IF THE RELATIVE CHANGE
C      IN VECTORS 'LOADS' AND 'OLDLDS' IS GREATER THAN 'TOL'
C
      REAL LOADS(*),OLDLDS(*)
      ICON = 1
      BIG = 0.
      DO 1 I = 1,N
    1 IF (ABS(LOADS(I)).GT.BIG) BIG = ABS(LOADS(I))
      DO 2 I = 1,N
          IF (ABS(LOADS(I)-OLDLDS(I))/BIG.GT.TOL) ICON = 0
    2 OLDLDS(I) = LOADS(I)
      RETURN
      END
C
C
C
      SUBROUTINE CHOFAC(LB,ILB,N,IW)
C
C      CHOLESKI FACTORISATION OF BANDED LOWER TRIANGLE LB
C
      REAL LB(ILB,*)
      DO 1 I = 1,N
          X = 0.0
          DO 2 J = 1,IW
              X = X + LB(I,J)**2
    2     CONTINUE
          LB(I,IW+1) = SQRT(LB(I,IW+1)-X)
          DO 3 K = 1,IW
              X = 0.0
              IF (I+K.LE.N) THEN
                  IF (K.NE.W) THEN
                      DO 4 L = IW - K,1,-1
                          X = X + LB(I+K,L)*LB(I,L+K)
    4                 CONTINUE
                  END IF
                  IA = I + K
                  IB = IW - K + 1
                  LB(IA,IB) = (LB(IA,IB)-X)/LB(I,IW+1)
              END IF
    3     CONTINUE
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE CHOSUB(LB,ILB,B,N,IW)
C
C      CHOLESKI FORWARD AND BACKWARD SUBSTITUTION COMBINED
C
      REAL LB(ILB,*),B(*)
      B(1) = B(1)/LB(1,IW+1)
      DO 1 I = 2,N
          X = 0.0
          K = 1
          IF (I.LE.IW+1) K = IW - I + 2
          DO 2 J = K,IW
              X = X + LB(I,J)*B(I+J-IW-1)
    2     CONTINUE
          B(I) = (B(I)-X)/LB(I,IW+1)
    1 CONTINUE
      B(N) = B(N)/LB(N,IW+1)
      DO 3 I = N - 1,1,-1
          X = 0.0
          L = I + IW
          IF (I.GT.N-IW) L = N
          M = I + 1
          DO 4 J = M,L
              X = X + LB(J,IW+I-J+1)*B(J)
    4     CONTINUE
          B(I) = (B(I)-X)/LB(I,IW+1)
    3 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE GAULAG(SAMP,ISAMP,NGL)
C
C      WEIGHTS AND SAMPLING POINTS
C      FOR GAUSS-LAGUERRE QUADRATURE
C
      REAL SAMP(ISAMP,*)
      GO TO (1,2,3,4,5),NGL
    1 SAMP(1,1) = 1.
      SAMP(1,2) = 1.
      GO TO 100
    2 SAMP(1,1) = 0.585786437627
      SAMP(2,1) = 3.414213562373
      SAMP(1,2) = 0.853553390593
      SAMP(2,2) = 0.146446609407
      GO TO 100
    3 SAMP(1,1) = 0.415774556783
      SAMP(2,1) = 2.294280360279
      SAMP(3,1) = 6.289945082937
      SAMP(1,2) = 0.711093009929
      SAMP(2,2) = 0.278517733569
      SAMP(3,2) = 0.0103892565016
      GO TO 100
    4 SAMP(1,1) = 0.322547689619
      SAMP(2,1) = 1.745761101158
      SAMP(3,1) = 4.536620296921
      SAMP(4,1) = 9.395070912301
      SAMP(1,2) = 0.603154104342
      SAMP(2,2) = 0.357418692438
      SAMP(3,2) = 0.038887908515
      SAMP(4,2) = 0.000539294705561
      GO TO 100
    5 SAMP(1,1) = 0.263560319718
      SAMP(2,1) = 1.413403059107
      SAMP(3,1) = 3.596425771041
      SAMP(4,1) = 7.085810005859
      SAMP(5,1) = 12.640800844276
      SAMP(1,2) = 0.521755610583
      SAMP(2,2) = 0.398666811083
      SAMP(3,2) = 0.0759424496817
      SAMP(4,2) = 0.00361175867992
      SAMP(5,2) = 0.0000233699723858
  100 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE GAULEG(SAMP,ISAMP,NGP)
C
C      WEIGHTS AND SAMPLING POINTS
C      FOR GAUSS-LEGENDRE QUADRATURE
C
      REAL SAMP(ISAMP,*)
      GO TO (1,2,3,4,5,6,7),NGP
    1 SAMP(1,1) = 0.
      SAMP(1,2) = 2.
      GO TO 100
    2 SAMP(1,1) = 1./SQRT(3.)
      SAMP(2,1) = -SAMP(1,1)
      SAMP(1,2) = 1.
      SAMP(2,2) = 1.
      GO TO 100
    3 SAMP(1,1) = .2*SQRT(15.)
      SAMP(2,1) = .0
      SAMP(3,1) = -SAMP(1,1)
      SAMP(1,2) = 5./9.
      SAMP(2,2) = 8./9.
      SAMP(3,2) = SAMP(1,2)
      GO TO 100
    4 SAMP(1,1) = .861136311594053
      SAMP(2,1) = .339981043584856
      SAMP(3,1) = -SAMP(2,1)
      SAMP(4,1) = -SAMP(1,1)
      SAMP(1,2) = .347854845137454
      SAMP(2,2) = .652145154862546
      SAMP(3,2) = SAMP(2,2)
      SAMP(4,2) = SAMP(1,2)
      GO TO 100
    5 SAMP(1,1) = .906179845938664
      SAMP(2,1) = .538469310105683
      SAMP(3,1) = .0
      SAMP(4,1) = -SAMP(2,1)
      SAMP(5,1) = -SAMP(1,1)
      SAMP(1,2) = .236926885056189
      SAMP(2,2) = .478628670499366
      SAMP(3,2) = .568888888888889
      SAMP(4,2) = SAMP(2,2)
      SAMP(5,2) = SAMP(1,2)
      GO TO 100
    6 SAMP(1,1) = .932469514203152
      SAMP(2,1) = .661209386466265
      SAMP(3,1) = .238619186083197
      SAMP(4,1) = -SAMP(3,1)
      SAMP(5,1) = -SAMP(2,1)
      SAMP(6,1) = -SAMP(1,1)
      SAMP(1,2) = .171324492379170
      SAMP(2,2) = .360761573048139
      SAMP(3,2) = .467913934572691
      SAMP(4,2) = SAMP(3,2)
      SAMP(5,2) = SAMP(2,2)
      SAMP(6,2) = SAMP(1,2)
      GO TO 100
    7 SAMP(1,1) = .949107912342759
      SAMP(2,1) = .741531185599394
      SAMP(3,1) = .405845151377397
      SAMP(4,1) = .0
      SAMP(5,1) = -SAMP(3,1)
      SAMP(6,1) = -SAMP(2,1)
      SAMP(7,1) = -SAMP(1,1)
      SAMP(1,2) = .129484966168870
      SAMP(2,2) = .279705391489277
      SAMP(3,2) = .381830050505119
      SAMP(4,2) = .417959183673469
      SAMP(5,2) = SAMP(3,2)
      SAMP(6,2) = SAMP(2,2)
      SAMP(7,2) = SAMP(1,2)
  100 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE LDLFOR(A,IA,B,N)
C
C      FORWARD SUBSTITUTION ON A LOWER TRIANGLE
C      STORED AS AN UPPER TRIANGLE
C
      REAL A(IA,*),B(*)
      DO 1 I = 1,N
          SUM = B(I)
          IF (I.GT.1) THEN
              DO 2 J = 1,I - 1
                  SUM = SUM - A(J,I)*B(J)
    2         CONTINUE
          END IF
          B(I) = SUM/A(I,I)
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE LDLT(A,IA,D,N)
C
C      L*D*LT FACTORISATION OF A SQUARE MATRIX A
C      LT OVERWRITES A
C
      REAL A(IA,*),D(*)
      DO 1 K = 1,N - 1
          D(1) = A(1,1)
          IF (ABS(A(K,K)).GT.1.E-6) THEN
              DO 2 I = K + 1,N
                  X = A(I,K)/A(K,K)
                  DO 3 J = K + 1,N
                      A(I,J) = A(I,J) - A(K,J)*X
    3             CONTINUE
                  D(I) = A(I,I)
    2         CONTINUE
          ELSE
              WRITE (6,1000)
              WRITE (6,*) K
          END IF
    1 CONTINUE
 1000 FORMAT ('ZERO PIVOT FOUND IN THIS LINE:')
      RETURN
      END
C
C
C
      SUBROUTINE LUFAC(A,UPTRI,LOWTRI,IA,N)                               
C
C         L*U FACTORISATION OF A SQUARE MATRIX A                                
C
      REAL A(IA,*),UPTRI(IA,*),LOWTRI(IA,*)                                     
      CALL NULL(UPTRI,IA,N,N)                                                   
      CALL NULL(LOWTRI,IA,N,N)                                                  
      DO 10 I=1,N                                                               
         UPTRI(1,I)=A(1,I)                                                      
   10    LOWTRI(I,I)=1.0                                                        
      DO 1 K=1,N-1                                                              
         IF(ABS(UPTRI(K,K)).GT.1.E-6)THEN                                       
            DO 2 I=K+1,N                                                        
C      LOWER TRIANGULAR COMPONENTS                                              
               DO 3 J=1,I-1                                                     
                  SUM=0.0                                                       
                     DO 4 L=1,J-1                                               
                        SUM=SUM-LOWTRI(I,L)*UPTRI(L,J)                          
    4                CONTINUE                                                   
                  LOWTRI(I,J)=(A(I,J)+SUM)/UPTRI(J,J)                           
    3          CONTINUE                                                         
C       UPPER TRIANGULAR COMPONENTS                                             
               DO 5 J=I,N                                                       
                  SUM=0.0                                                       
                  DO 6 L=1,I-1                                                  
                     SUM=SUM-LOWTRI(I,L)*UPTRI(L,J)                             
    6             CONTINUE                                                      
                  UPTRI(I,J)=A(I,J)+SUM                                         
    5          CONTINUE                                                         
    2       CONTINUE                                                            
         ELSE                                                                   
         WRITE(6,1000)                                                          
         WRITE(6,*)K                                                            
         STOP                                                                   
         ENDIF                                                                  
    1 CONTINUE                                                                  
 1000 FORMAT('ZERO PIVOT FOUND IN THIS ROW')                                    
      RETURN                                                                    
      END                                                                       
C
C
C
      SUBROUTINE LUPFAC(A,IA,N,ROW)
C
C      L*U FACTORISATION OF A SQUARE MATRIX WITH PIVOTS
C
      REAL A(IA,*)
      INTEGER ROW(*)
      DO 1 I = 1,N
    1 ROW(I) = I
      DO 2 I = 1,N - 1
          IP = I
          PVAL = A(ROW(IP),IP)
          DO 3 J = I + 1,N
              IF (ABS(A(ROW(J),I)).GT.ABS(PVAL)) THEN
                  IP = J
                  PVAL = A(ROW(J),I)
              END IF
    3     CONTINUE
          IF (ABS(PVAL).LT.1.E-10) THEN
              WRITE (6,1000)
              STOP
          END IF
          IH = ROW(IP)
          ROW(IP) = ROW(I)
          ROW(I) = IH
          DO 4 J = I + 1,N
              IE = ROW(J)
              PIVOT = A(IE,I)/PVAL
              A(IE,I) = PIVOT
              IROW = ROW(I)
              DO 5 K = I + 1,N
                  A(IE,K) = A(IE,K) - A(IROW,K)*PIVOT
    5         CONTINUE
    4     CONTINUE
    2 CONTINUE
      IF (ABS(A(ROW(N),N)).LT.1.E-10) THEN
          WRITE (6,1000)
          STOP
      END IF
 1000 FORMAT ('SINGULAR EQUATIONS DETECTED')
      RETURN
      END
C
C
C
      SUBROUTINE LUPSOL(A,IA,B,SOL,N,ROW)
C
C      FORWARD AND BACKWARD SUBSTITUTION COMBINED WITH PIVOTS
C
      REAL A(IA,*),B(*),SOL(*)
      INTEGER ROW(*)
      DO 1 I = 1,N
          IROW = ROW(I)
          SUM = B(IROW)
          IF (I.GT.1) THEN
              DO 2 J = 1,I - 1
                  SUM = SUM - A(IROW,J)*B(ROW(J))
    2         CONTINUE
              B(IROW) = SUM
          END IF
    1 CONTINUE
      DO 3 I = N,1,-1
          IROW = ROW(I)
          SUM = B(IROW)
          IF (I.LT.N) THEN
              DO 4 J = I + 1,N
                  SUM = SUM - A(IROW,J)*B(ROW(J))
    4         CONTINUE
          END IF
          B(IROW) = SUM/A(IROW,I)
    3 CONTINUE
      DO 5 I = 1,N
    5 SOL(I) = B(ROW(I))
      RETURN
      END
C
C
C
      SUBROUTINE MATINV(A,IA,N)
C
C      FORMS THE INVERSE OF A MATRIX
C      USING GAUSS-JORDAN TRANSFORMATION
C
      REAL A(IA,*)
      DO 1 K = 1,N
          CON = A(K,K)
          A(K,K) = 1.
          DO 2 J = 1,N
    2     A(K,J) = A(K,J)/CON
          DO 1 I = 1,N
              IF (I.EQ.K) GO TO 1
              CON = A(I,K)
              A(I,K) = 0.
              DO 3 J = 1,N
    3         A(I,J) = A(I,J) - A(K,J)*CON
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE MATMUL(A,IA,B,IB,C,IC,L,M,N)
C
C      PRODUCT OF TWO MATRICES
C
      REAL A(IA,*),B(IB,*),C(IC,*)
      DO 1 I = 1,L
          DO 1 J = 1,N
              X = 0.0
              DO 2 K = 1,M
    2         X = X + A(I,K)*B(K,J)
              C(I,J) = X
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE MATRAN(A,IA,B,IB,M,N)
C
C      FORMS THE TRANSPOSE OF A MATRIX
C
      REAL A(IA,*),B(IB,*)
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 A(J,I) = B(I,J)
      RETURN
      END
C
C
C
      SUBROUTINE MSMULT(A,IA,C,M,N)
C
C      MULTIPLIES A MATRIX BY A SCALAR
C
      REAL A(IA,*)
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 A(I,J) = A(I,J)*C
      RETURN
      END
C
C
C
      SUBROUTINE MVMULT(M,IM,V,K,L,Y)
C
C      MULTIPLIES A MATRIX BY A VECTOR
C
      REAL M(IM,*),V(*),Y(*)
      DO 1 I = 1,K
          X = 0.
          DO 2 J = 1,L
    2     X = X + M(I,J)*V(J)
          Y(I) = X
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE NEWCOT(WEIGHT,NNC)
C
C      WEIGHTING COEFFICIENTS
C      FOR NEWTON-COTES QUADRATURE
C
      REAL WEIGHT(*)
      GO TO (1,2,3,4,5) NNC
    1 WEIGHT(1) = 1.
      GO TO 100
    2 WEIGHT(1) = .5
      WEIGHT(2) = .5
      GO TO 100
    3 WEIGHT(1) = 0.333333333333333
      WEIGHT(2) = 1.333333333333333
      WEIGHT(3) = 0.333333333333333
      GO TO 100
    4 WEIGHT(1) = 0.375
      WEIGHT(2) = 1.125
      WEIGHT(3) = 1.125
      WEIGHT(4) = 0.375
      GO TO 100
    5 WEIGHT(1) = 0.311111111111111
      WEIGHT(2) = 1.422222222222222
      WEIGHT(3) = 0.533333333333333
      WEIGHT(4) = 1.422222222222222
      WEIGHT(5) = 0.311111111111111
  100 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE NULL(A,IA,M,N)
C
C      NULLS A 2-D ARRAY
C
      REAL A(IA,*)
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 A(I,J) = 0.0
      RETURN
      END
C
C
C
      SUBROUTINE NULVEC(VEC,N)
C
C      NULLS A COLUMN VECTOR
C
      REAL VEC(*)
      DO 1 I = 1,N
    1 VEC(I) = 0.
      RETURN
      END
C
C
C
      SUBROUTINE PRINTA(A,IA,M,N,ICH)
C
C      WRITES A 2-D ARRAY TO OUTPUT CHANNEL 'ICH'
C
      REAL A(IA,*)
      DO 1 I = 1,M
    1 WRITE (ICH,2) (A(I,J),J=1,N)
    2 FORMAT (1X,6E12.4)
      RETURN
      END
C
C
C
      SUBROUTINE PRINTV(VEC,N,ICH)
C
C      WRITES A COLUMN VECTOR TO OUTPUT CHANNEL 'ICH'
C
      REAL VEC(*)
      WRITE (ICH,1) (VEC(I),I=1,N)
    1 FORMAT (1X,6E12.4)
      RETURN
      END
C
C
C
      SUBROUTINE SKYFAC(A,N,KDIAG)
C
C      CHOLESKI FACTORISATION OF VARIABLE BANDWIDTH
C      A STORED AS A VECTOR AND OVERWRITTEN
C
      REAL A(*)
      INTEGER KDIAG(*)
      A(1) = SQRT(A(1))
      DO 1 I = 2,N
          KI = KDIAG(I) - I
          L = KDIAG(I-1) - KI + 1
          DO 2 J = L,I
              X = A(KI+J)
              KJ = KDIAG(J) - J
              IF (J.NE.1) THEN
                  LL = KDIAG(J-1) - KJ + 1
                  LL = MAX0(L,LL)
                  IF (LL.NE.J) THEN
                      M = J - 1
                      DO 3 K = LL,M
                          X = X - A(KI+K)*A(KJ+K)
    3                 CONTINUE
                  END IF
              END IF
              A(KI+J) = X/A(KJ+J)
    2     CONTINUE
          A(KI+I) = SQRT(X)
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE SKYSUB(A,B,N,KDIAG)
C
C      CHOLESKI FORWARD AND BACKWARD SUBSTITUTION COMBINED
C      VARIABLE BANDWIDTH FACTORED A STORED AS A VECTOR
C
      REAL A(*),B(*)
      INTEGER KDIAG(*)
      B(1) = B(1)/A(1)
      DO 1 I = 2,N
          KI = KDIAG(I) - I
          L = KDIAG(I-1) - KI + 1
          X = B(I)
          IF (L.NE.I) THEN
              M = I - 1
              DO 2 J = L,M
                  X = X - A(KI+J)*B(J)
    2         CONTINUE
          END IF
          B(I) = X/A(KI+I)
    1 CONTINUE
      DO 3 IT = 2,N
          I = N + 2 - IT
          KI = KDIAG(I) - I
          X = B(I)/A(KI+I)
          B(I) = X
          L = KDIAG(I-1) - KI + 1
          IF (L.NE.I) THEN
              M = I - 1
              DO 4 K = L,M
                  B(K) = B(K) - X*A(KI+K)
    4         CONTINUE
          END IF
    3 CONTINUE
      B(1) = B(1)/A(1)
      RETURN
      END
C
C
C
      SUBROUTINE SOLVE(K,IK,U,F,N)
C
C      PERFORMS GAUSSIAN ELIMINATION WITH
C      PARTIAL PIVOTING ON A FULL N*N MATRIX
C
      REAL K(IK,*),F(*),U(*)
C
C      PIVOTING STAGE
C
      DO 1 I = 1,N - 1
          BIG = ABS(K(I,I))
          IHOLD = I
          DO 10 J = I + 1,N
              IF (ABS(K(J,I)).GT.BIG) THEN
                  BIG = ABS(K(J,I))
                  IHOLD = J
              END IF
   10     CONTINUE
          IF (IHOLD.NE.I) THEN
              DO 12 J = I,N
                  HOLD = K(I,J)
                  K(I,J) = K(IHOLD,J)
                  K(IHOLD,J) = HOLD
   12         CONTINUE
              HOLD = F(I)
              F(I) = F(IHOLD)
              F(IHOLD) = HOLD
          END IF
C
C      ELIMINATION STAGE
C
          DO 3 J = I + 1,N
              FAC = K(J,I)/K(I,I)
              DO 4 L = I,N
    4         K(J,L) = K(J,L) - K(I,L)*FAC
              F(J) = F(J) - F(I)*FAC
    3     CONTINUE
    1 CONTINUE
C
C      BACK-SUBSTITUTION STAGE
C
      DO 9 I = N,1,-1
          SUM = 0.
          DO 6 L = I + 1,N
    6     SUM = SUM + K(I,L)*U(L)
          U(I) = (F(I)-SUM)/K(I,I)
    9 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE SUBBAC(A,IA,B,N)
C
C      BACKWARD SUBSTITUTION ON AN UPPER TRIANGLE
C
      REAL A(IA,*),B(*)
      DO 1 I = N,1,-1
          SUM = B(I)
          IF (I.LT.N) THEN
              DO 2 J = I + 1,N
                  SUM = SUM - A(I,J)*B(J)
    2         CONTINUE
          END IF
          B(I) = SUM/A(I,I)
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE SUBFOR(A,IA,B,N)
C
C      FORWARD SUBSTITUTION ON A LOWER TRIANGLE
C
      REAL A(IA,*),B(*)
      DO 1 I = 1,N
          SUM = B(I)
          IF (I.GT.1) THEN
              DO 2 J = 1,I - 1
                  SUM = SUM - A(I,J)*B(J)
    2         CONTINUE
          END IF
          B(I) = SUM/A(I,I)
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE VDOTV(V1,V2,DOTPR,N)
C
C            DOT PRODUCT V1*V2
C
      REAL V1(*),V2(*)
      DOTPR = 0.0
      DO 1 I = 1,N
    1 DOTPR = DOTPR + V1(I)*V2(I)
      RETURN
      END
C
C
C
      SUBROUTINE VECADD(A,B,C,N)
C
C      ADDS VECTORS  A+B=C
C
      REAL A(*),B(*),C(*)
      DO 1 I = 1,N
    1 C(I) = A(I) + B(I)
      RETURN
      END
C
C
C
      SUBROUTINE VECCOP(A,B,N)
C
C      COPIES VECTOR A INTO VECTOR B
C
      REAL A(*),B(*)
      DO 1 I = 1,N
    1 B(I) = A(I)
      RETURN
      END
C
C
C
      SUBROUTINE VECSUB(A,B,C,N)
C
C      VECTOR SUBTRACT C:=A-B
C
      REAL A(*),B(*),C(*)
      DO 1 I = 1,N
    1 C(I) = A(I) - B(I)
      RETURN
      END
C
C
C
      SUBROUTINE VSMULT(V,SCAL,N)
C
C            MULTIPLY A VECTOR BY A SCALAR
C
      REAL V(*)
      DO 1 I = 1,N
    1 V(I) = V(I)*SCAL
      RETURN
      END
C
C
C
      SUBROUTINE VVMULT(V1,V2,A,IA,M,N)
C
C      FORMS A VECTOR PRODUCT
C
      REAL V1(*),V2(*),A(IA,*)
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 A(I,J) = V1(I)*V2(J)
      RETURN
      END
