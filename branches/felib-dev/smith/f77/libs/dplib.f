      SUBROUTINE AXIKM(KM,CSA,E,ELL)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX FOR AN
C      AXIALLY LOADED LINE ELEMENT
C
      DOUBLE PRECISION CSA
      DOUBLE PRECISION E
      DOUBLE PRECISION ELL
      DOUBLE PRECISION KM(2,2)
 
      KM(1,1) = CSA*E/ELL
      KM(2,2) = KM(1,1)
      KM(1,2) = -KM(1,1)
      KM(2,1) = KM(1,2)
      RETURN
 
      END
      SUBROUTINE BACSUB(BK,LOADS,N,IW)
C
C      THIS SUBROUTINE PERFORMS THE GAUSSIAN BACK-SUBSTITUTION
C
      DOUBLE PRECISION SUM
      DOUBLE PRECISION BK(*),LOADS(*)
 
      LOADS(1) = LOADS(1)/BK(1)
      DO 1 I = 2,N
          SUM = LOADS(I)
          I1 = I - 1
          NKB = I - IW
          IF (NKB) 2,2,3
    2     NKB = 1
    3     DO 4 K = NKB,I1
              JN = (I-K)*N + K
              SUM = SUM - BK(JN)*LOADS(K)
    4     CONTINUE
          LOADS(I) = SUM/BK(I)
    1 CONTINUE
      DO 5 JJ = 2,N
          I = N - JJ + 1
          SUM = 0.D0
          I1 = I + 1
          NKB = I + IW
          IF (NKB-N) 7,7,6
    6     NKB = N
    7     DO 8 K = I1,NKB
              JN = (K-I)*N + I
    8     SUM = SUM + BK(JN)*LOADS(K)
          LOADS(I) = LOADS(I) - SUM/BK(I)
    5 CONTINUE
      RETURN
 
      END
      SUBROUTINE BANDRD(N,IW,A,IA,D,E,E2)
C
C      THIS SUBROUTINE TRANSFORMS A REAL SYMMETRIC BAND MATRIX A,
C      OF ORDER N AND BAND WIDTH IW,
C      TO TRIDIAGONAL FORM BY AN APPROPRIATE
C      SEQUENCE OF JACOBI ROTATIONS. DURING THE TRANSFORMATION THE
C      PROPERTY OF THE BAND MATRIX IS MAINTAINED. THE METHOD YIELDS
C      A TRIDIAGONAL MATRIX, THE DIAGONAL ELEMENTS OF WHICH ARE IN
C      D(N) AND OFF-DIAGONAL ELEMENTS IN E(N).
C
      INTEGER M,IW,N2,N,K,MAXR,IRR,IR,KR,J,JM,IUGL,J2,L,JL,MAXL,I,IA
      DOUBLE PRECISION G,B,S,C,C2,S2,CS,U,U1,A(IA,*),D(*),E(*),E2(*)
 
      N2 = N - 2
      IF (N2.LT.1) GO TO 180
      DO 160 K = 1,N2
          MAXR = IW
          IF (N-K.LT.IW) MAXR = N - K
          DO 140 IRR = 2,MAXR
              IR = 2 + MAXR - IRR
              KR = K + IR
              DO 120 J = KR,N,IW
                  IF (J.EQ.KR) GO TO 20
                  IF (G.EQ.0.0D0) GO TO 140
                  JM = J - IW
                  B = -A(JM-1,IW+1)/G
                  IUGL = J - IW
                  GO TO 40
 
   20             IF (A(K,IR+1).EQ.0.0D0) GO TO 140
                  B = -A(K,IR)/A(K,IR+1)
                  IUGL = K
   40             S = 1.0D0/SQRT(1.0D0+B*B)
                  C = B*S
                  C2 = C*C
                  S2 = S*S
                  CS = C*S
                  U = C2*A(J-1,1) - 2.0D0*CS*A(J-1,2) + S2*A(J,1)
                  U1 = S2*A(J-1,1) + 2.0D0*CS*A(J-1,2) + C2*A(J,1)
                  A(J-1,2) = CS* (A(J-1,1)-A(J,1)) + (C2-S2)*A(J-1,2)
                  A(J-1,1) = U
                  A(J,1) = U1
                  J2 = J - 2
                  DO 60 L = IUGL,J2
                      JL = J - L
                      U = C*A(L,JL) - S*A(L,JL+1)
                      A(L,JL+1) = S*A(L,JL) + C*A(L,JL+1)
                      A(L,JL) = U
   60             CONTINUE
                  JM = J - IW
                  IF (J.NE.KR) A(JM-1,IW+1) = C*A(JM-1,IW+1) - S*G
                  MAXL = IW - 1
                  IF (N-J.LT.IW-1) MAXL = N - J
                  IF (MAXL.LE.0) GO TO 100
                  DO 80 L = 1,MAXL
                      U = C*A(J-1,L+2) - S*A(J,L+1)
                      A(J,L+1) = S*A(J-1,L+2) + C*A(J,L+1)
                      A(J-1,L+2) = U
   80             CONTINUE
  100             IF (J+IW.GT.N) GO TO 120
                  G = -S*A(J,IW+1)
                  A(J,IW+1) = C*A(J,IW+1)
  120         CONTINUE
  140     CONTINUE
  160 CONTINUE
  180 E(1) = 0.0D0
      DO 200 I = 1,N
          D(I) = A(I,1)
  200 CONTINUE
      IF (2.GT.N) GO TO 240
      DO 220 I = 2,N
          E(I) = A(I-1,2)
  220 CONTINUE
  240 DO 260 I = 1,N
          E2(I) = E(I)*E(I)
  260 CONTINUE
      RETURN
 
      END
      SUBROUTINE BANMUL(KB,IKB,LOADS,ANS,N,IW)
C
C      THIS SUBROUTINE MULTIPLIES A MATRIX BY A VECTOR
C      THE MATRIX IS SYMMETRICAL WITH ITS LOWER TRIANGLE
C      STORED AS A RECTANGLE
C
      DOUBLE PRECISION X
      DOUBLE PRECISION KB(IKB,*),LOADS(*),ANS(*)
 
      DO 1 I = 1,N
          X = 0.D0
          J = IW + 1
    2     IF (I+J.LE.IW+1) GO TO 3
          X = X + KB(I,J)*LOADS(I+J-IW-1)
    3     J = J - 1
          IF (J.NE.0) GO TO 2
          J = IW
    6     IF (I-J.GE.N-IW) GO TO 7
          X = X + KB(I-J+IW+1,J)*LOADS(I-J+IW+1)
    7     J = J - 1
          IF (J.NE.0) GO TO 6
          ANS(I) = X
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE BANRED(BK,N,IW)
C
C      THIS SUBROUTINE PERFORMS GAUSSIAN REDUCTION OF
C      THE STIFFNESS MATRIX STORED AS A VECTOR BK(N*(IW+1))
C
      DOUBLE PRECISION SUM
      DOUBLE PRECISION BK(*)
 
      DO 1 I = 2,N
          IL1 = I - 1
          KBL = IL1 + IW + 1
          IF (KBL-N) 3,3,2
    2     KBL = N
    3     DO 1 J = I,KBL
              IJ = (J-I)*N + I
              SUM = BK(IJ)
              NKB = J - IW
              IF (NKB) 4,4,5
    4         NKB = 1
    5         IF (NKB-IL1) 6,6,8
    6         DO 7 M = NKB,IL1
                  NI = (I-M)*N + M
                  NJ = (J-M)*N + M
    7         SUM = SUM - BK(NI)*BK(NJ)/BK(M)
    8         BK(IJ) = SUM
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE BANTML(KB,IKB,LOADS,ANS,N,IW)
C
C      THIS SUBROUTINE MULTIPLIES AN UNSYMMETRIC BANDED MATRIX
C      'PB' BY THE VECTOR 'LOADS'.
C
      DOUBLE PRECISION X
      DOUBLE PRECISION KB(IKB,*),LOADS(*),ANS(*)
 
      DO 1 I = 1,N
          X = 0.D0
          K = IW + 2
          L = IW + IW + 1
          DO 2 J = 1,L
              K = K - 1
              IF (I-K+1.GT.N) GO TO 2
              IF (I-K+1.LT.1) GO TO 2
              X = X + KB(I,J)*LOADS(I-K+1)
    2     CONTINUE
          ANS(I) = X
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE BCMASS(MM,RHO,AREA,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE CONSISTENT MASS MATRIX OF
C      AN INCLINED 2-D BEAM-COLUMN ELEMENT
C
      DOUBLE PRECISION RHO
      DOUBLE PRECISION AREA
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION ELL
      DOUBLE PRECISION C
      DOUBLE PRECISION S
      DOUBLE PRECISION FAC
      DOUBLE PRECISION MM(6,6),COORD(ICOORD,*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      X2 = COORD(IP,3)
      Y2 = COORD(IP,4)
      ELL = SQRT((Y2-Y1)**2+ (X2-X1)**2)
      C = (X2-X1)/ELL
      S = (Y2-Y1)/ELL
      MM(1,1) = 140.D0*C*C + 156.D0*S*S
      MM(4,4) = MM(1,1)
      MM(2,2) = 156.D0*C*C + 140.D0*S*S
      MM(5,5) = MM(2,2)
      MM(3,3) = 4.D0*ELL*ELL
      MM(6,6) = MM(3,3)
      MM(1,2) = -16.D0*C*S
      MM(4,5) = MM(1,2)
      MM(1,5) = -MM(1,2)
      MM(2,4) = -MM(1,2)
      MM(1,3) = -22.D0*ELL*S
      MM(4,6) = -MM(1,3)
      MM(2,3) = 22.D0*ELL*C
      MM(5,6) = -MM(2,3)
      MM(1,4) = 70.D0*C*C + 54.D0*S*S
      MM(1,6) = 13.D0*ELL*S
      MM(3,4) = -MM(1,6)
      MM(2,5) = 54.D0*C*C + 70.D0*S*S
      MM(2,6) = -13.D0*ELL*C
      MM(3,5) = -MM(2,6)
      MM(3,6) = -3.D0*ELL*ELL
      FAC = RHO*AREA*ELL/420.D0
      DO 1 I = 1,6
          DO 1 J = I,6
              MM(I,J) = MM(I,J)*FAC
    1 MM(J,I) = MM(I,J)
      RETURN
 
      END
      SUBROUTINE BEAMKM(KM,EI,ELL)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX OF A
C      HORIZONTAL BEAM ELEMENT(BENDING ONLY)
C
      DOUBLE PRECISION EI
      DOUBLE PRECISION ELL
      DOUBLE PRECISION KM(4,4)
 
      KM(1,1) = 12.D0*EI/ (ELL*ELL*ELL)
      KM(3,3) = KM(1,1)
      KM(1,2) = 6.D0*EI/ (ELL*ELL)
      KM(2,1) = KM(1,2)
      KM(1,4) = KM(1,2)
      KM(4,1) = KM(1,4)
      KM(1,3) = -KM(1,1)
      KM(3,1) = KM(1,3)
      KM(3,4) = -KM(1,2)
      KM(4,3) = KM(3,4)
      KM(2,3) = KM(3,4)
      KM(3,2) = KM(2,3)
      KM(2,2) = 4.D0*EI/ELL
      KM(4,4) = KM(2,2)
      KM(2,4) = 2.D0*EI/ELL
      KM(4,2) = KM(2,4)
      RETURN
 
      END
      SUBROUTINE BEAMKP(KP,ELL)
C
C      THIS SUBROUTINE FORMS THE TERMS OF THE BEAM STIFFNESS
C      MATRIX DUE TO AXIAL LOADING
C
      DOUBLE PRECISION ELL
      DOUBLE PRECISION KP(4,4)
 
      KP(1,1) = 1.2D0/ELL
      KP(1,2) = 0.1D0
      KP(2,1) = 0.1D0
      KP(1,3) = -1.2D0/ELL
      KP(3,1) = -1.2D0/ELL
      KP(1,4) = 0.1D0
      KP(4,1) = 0.1D0
      KP(2,2) = 2.0D0*ELL/15.0D0
      KP(2,3) = -0.1D0
      KP(3,2) = -0.1D0
      KP(2,4) = -ELL/30.0D0
      KP(4,2) = -ELL/30.0D0
      KP(3,3) = 1.2D0/ELL
      KP(3,4) = -0.1D0
      KP(4,3) = -0.1D0
      KP(4,4) = 2.0D0*ELL/15.0D0
      RETURN
 
      END
      SUBROUTINE BISECT(N,ACHEPS,D,E,IFAIL)
C
C      THIS SUBROUTINE FINDS THE EIGENVALUES OF A TRIDIAGONAL
C      MATRIX,
C      T, GIVEN WITH ITS DIAGONAL ELEMENTS IN THE ARRAY D(N) AND
C      ITS SUBDIAGONAL ELEMENTS IN THE LAST N - 1 STORES OF THE
C      ARRAY E(N), USING QL TRANSFORMATIONS. THE EIGENVALUES ARE
C      OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE ARRAY D IN
C      ASCENDING ORDER. THE SUBROUTINE WILL FAIL IF ANY ONE
C      EIGENVALUE TAKES MORE THAN 30 ITERATIONS.
C
      INTEGER P01AAF,ISAVE,IFAIL,N,I,L,J,M,I1,M1,II
      DOUBLE PRECISION B,F,H,ACHEPS,G,P,R,C,S,D(*),E(*)
 
      ISAVE = IFAIL
      IF (N.EQ.1) GO TO 40
      DO 20 I = 2,N
          E(I-1) = E(I)
   20 CONTINUE
   40 E(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      DO 340 L = 1,N
          J = 0
          H = ACHEPS* (ABS(D(L))+ABS(E(L)))
          IF (B.LT.H) B = H
C     LOOK FOR SMALL SUB DIAGONAL ELEMENT
          DO 60 M = L,N
              IF (ABS(E(M)).LE.B) GO TO 80
   60     CONTINUE
   80     IF (M.EQ.L) GO TO 260
  100     IF (J.EQ.30) GO TO 360
          J = J + 1
C     FORM SHIFT
          G = D(L)
          H = D(L+1) - G
          IF (ABS(H).GE.ABS(E(L))) GO TO 120
          P = H*0.5D0/E(L)
          R = SQRT(P*P+1.0D0)
          H = P + R
          IF (P.LT.0.0D0) H = P - R
          D(L) = E(L)/H
          GO TO 140
 
  120     P = 2.0D0*E(L)/H
          R = SQRT(P*P+1.0D0)
          D(L) = E(L)*P/ (1.0D0+R)
  140     H = G - D(L)
          I1 = L + 1
          IF (I1.GT.N) GO TO 180
          DO 160 I = I1,N
              D(I) = D(I) - H
  160     CONTINUE
  180     F = F + H
C     QL TRANSFORMATION
          P = D(M)
          C = 1.0D0
          S = 0.0D0
          M1 = M - 1
          DO 240 II = L,M1
              I = M1 - II + L
              G = C*E(I)
              H = C*P
              IF (ABS(P).LT.ABS(E(I))) GO TO 200
              C = E(I)/P
              R = SQRT(C*C+1.0D0)
              E(I+1) = S*P*R
              S = C/R
              C = 1.0D0/R
              GO TO 220
 
  200         C = P/E(I)
              R = SQRT(C*C+1.0D0)
              E(I+1) = S*E(I)*R
              S = 1.0D0/R
              C = C/R
  220         P = C*D(I) - S*G
              D(I+1) = H + S* (C*G+S*D(I))
  240     CONTINUE
          E(L) = S*P
          D(L) = C*P
          IF (ABS(E(L)).GT.B) GO TO 100
  260     P = D(L) + F
C     ORDER EIGENVALUE
          IF (L.EQ.1) GO TO 300
          DO 280 II = 2,L
              I = L - II + 2
              IF (P.GE.D(I-1)) GO TO 320
              D(I) = D(I-1)
  280     CONTINUE
  300     I = 1
  320     D(I) = P
  340 CONTINUE
      IFAIL = 0
      RETURN
 
  360 IFAIL = 1
      RETURN
 
      END
      SUBROUTINE BMCOL2(KM,EA,EI,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX OF AN
C      INCLINED 2-D BEAM-COLUMN ELEMENT
C
      DOUBLE PRECISION EA
      DOUBLE PRECISION EI
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION ELL
      DOUBLE PRECISION C
      DOUBLE PRECISION S
      DOUBLE PRECISION E1
      DOUBLE PRECISION E2
      DOUBLE PRECISION E3
      DOUBLE PRECISION E4
      DOUBLE PRECISION KM(6,6),COORD(ICOORD,*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      X2 = COORD(IP,3)
      Y2 = COORD(IP,4)
      ELL = SQRT((Y2-Y1)**2+ (X2-X1)**2)
      C = (X2-X1)/ELL
      S = (Y2-Y1)/ELL
      E1 = EA/ELL
      E2 = 12.D0*EI/ (ELL*ELL*ELL)
      E3 = EI/ELL
      E4 = 6.D0*EI/ (ELL*ELL)
      KM(1,1) = C*C*E1 + S*S*E2
      KM(4,4) = KM(1,1)
      KM(1,2) = S*C* (E1-E2)
      KM(2,1) = KM(1,2)
      KM(4,5) = KM(1,2)
      KM(5,4) = KM(4,5)
      KM(1,3) = -S*E4
      KM(3,1) = KM(1,3)
      KM(1,6) = KM(1,3)
      KM(6,1) = KM(1,6)
      KM(3,4) = S*E4
      KM(4,3) = KM(3,4)
      KM(4,6) = KM(3,4)
      KM(6,4) = KM(4,6)
      KM(1,4) = -KM(1,1)
      KM(4,1) = KM(1,4)
      KM(1,5) = S*C* (-E1+E2)
      KM(5,1) = KM(1,5)
      KM(2,4) = KM(1,5)
      KM(4,2) = KM(2,4)
      KM(2,2) = S*S*E1 + C*C*E2
      KM(5,5) = KM(2,2)
      KM(2,5) = -KM(2,2)
      KM(5,2) = KM(2,5)
      KM(2,3) = C*E4
      KM(3,2) = KM(2,3)
      KM(2,6) = KM(2,3)
      KM(6,2) = KM(2,6)
      KM(3,3) = 4.D0*E3
      KM(6,6) = KM(3,3)
      KM(3,5) = -C*E4
      KM(5,3) = KM(3,5)
      KM(5,6) = KM(3,5)
      KM(6,5) = KM(5,6)
      KM(3,6) = 2.D0*E3
      KM(6,3) = KM(3,6)
      RETURN
 
      END
      SUBROUTINE BMCOL3(KM,EA,EIY,EIZ,GJ,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX OF A
C      GENERAL 3-D BEAM-COLUMN ELEMENT
C
      DOUBLE PRECISION EA
      DOUBLE PRECISION EIY
      DOUBLE PRECISION EIZ
      DOUBLE PRECISION GJ
      DOUBLE PRECISION PI
      DOUBLE PRECISION GAMA
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Z1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION Z2
      DOUBLE PRECISION XL
      DOUBLE PRECISION YL
      DOUBLE PRECISION ZL
      DOUBLE PRECISION ELL
      DOUBLE PRECISION CG
      DOUBLE PRECISION SG
      DOUBLE PRECISION DEN
      DOUBLE PRECISION A1
      DOUBLE PRECISION A2
      DOUBLE PRECISION A3
      DOUBLE PRECISION A4
      DOUBLE PRECISION A5
      DOUBLE PRECISION A6
      DOUBLE PRECISION A7
      DOUBLE PRECISION A8
      DOUBLE PRECISION X
      DOUBLE PRECISION SUM
      DOUBLE PRECISION KM(12,12),COORD(ICOORD,*),T(12,12),TT(12,12),
     +                 R0(3,3),C(12,12)
 
      PI = 4.D0*ATAN(1.D0)
      GAMA = COORD(IP,7)*PI/180.D0
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      Z1 = COORD(IP,3)
      X2 = COORD(IP,4)
      Y2 = COORD(IP,5)
      Z2 = COORD(IP,6)
      XL = X2 - X1
      YL = Y2 - Y1
      ZL = Z2 - Z1
      ELL = SQRT(XL*XL+YL*YL+ZL*ZL)
      CG = COS(GAMA)
      SG = SIN(GAMA)
      DEN = ELL*SQRT(XL*XL+ZL*ZL)
      DO 1 I = 1,12
          DO 1 J = 1,12
              KM(I,J) = 0.D0
              T(I,J) = 0.D0
    1 TT(I,J) = 0.D0
      A1 = EA/ELL
      A2 = 12.D0*EIZ/ (ELL*ELL*ELL)
      A3 = 12.D0*EIY/ (ELL*ELL*ELL)
      A4 = 6.D0*EIZ/ (ELL*ELL)
      A5 = 6.D0*EIY/ (ELL*ELL)
      A6 = 4.D0*EIZ/ELL
      A7 = 4.D0*EIY/ELL
      A8 = GJ/ELL
      KM(1,1) = A1
      KM(7,7) = A1
      KM(1,7) = -A1
      KM(7,1) = -A1
      KM(2,2) = A2
      KM(8,8) = A2
      KM(2,8) = -A2
      KM(8,2) = -A2
      KM(3,3) = A3
      KM(9,9) = A3
      KM(3,9) = -A3
      KM(9,3) = -A3
      KM(4,4) = A8
      KM(10,10) = A8
      KM(4,10) = -A8
      KM(10,4) = -A8
      KM(5,5) = A7
      KM(11,11) = A7
      KM(5,11) = .5D0*A7
      KM(11,5) = .5D0*A7
      KM(6,6) = A6
      KM(12,12) = A6
      KM(6,12) = .5D0*A6
      KM(12,6) = .5D0*A6
      KM(2,6) = A4
      KM(6,2) = A4
      KM(2,12) = A4
      KM(12,2) = A4
      KM(6,8) = -A4
      KM(8,6) = -A4
      KM(8,12) = -A4
      KM(12,8) = -A4
      KM(5,9) = A5
      KM(9,5) = A5
      KM(9,11) = A5
      KM(11,9) = A5
      KM(3,5) = -A5
      KM(5,3) = -A5
      KM(3,11) = -A5
      KM(11,3) = -A5
      IF (DEN.EQ.0.D0) GO TO 50
      R0(1,1) = XL/ELL
      R0(1,2) = YL/ELL
      R0(1,3) = ZL/ELL
      R0(2,1) = (-XL*YL*CG-ELL*ZL*SG)/DEN
      R0(2,2) = DEN*CG/ (ELL*ELL)
      R0(2,3) = (-YL*ZL*CG+ELL*XL*SG)/DEN
      R0(3,1) = (XL*YL*SG-ELL*ZL*CG)/DEN
      R0(3,2) = -DEN*SG/ (ELL*ELL)
      R0(3,3) = (YL*ZL*SG+ELL*XL*CG)/DEN
      GO TO 60
 
   50 R0(1,1) = 0.D0
      R0(1,3) = 0.D0
      R0(2,2) = 0.D0
      R0(3,2) = 0.D0
      R0(1,2) = 1.D0
      R0(2,1) = -CG
      R0(3,3) = CG
      R0(2,3) = SG
      R0(3,1) = SG
   60 CONTINUE
      DO 2 I = 1,3
          DO 2 J = 1,3
              X = R0(I,J)
              DO 2 K = 0,9,3
                  T(I+K,J+K) = X
    2 TT(J+K,I+K) = X
      DO 3 I = 1,12
          DO 3 J = 1,12
              SUM = 0.D0
              DO 4 K = 1,12
    4         SUM = SUM + KM(I,K)*T(K,J)
    3 C(I,J) = SUM
      DO 5 I = 1,12
          DO 5 J = 1,12
              SUM = 0.D0
              DO 6 K = 1,12
    6         SUM = SUM + TT(I,K)*C(K,J)
    5 KM(I,J) = SUM
      RETURN
 
      END
      SUBROUTINE BNONAX(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD,
     +                  IFLAG,LTH)
C
C      THIS SUBROUTINE FORMS THE STRAIN-DISPLACEMENT MATRIX FOR
C      AXISYMMETRIC SOLIDS SUBJECTED TO NON-AXISYMMETRIC LOADING
C
      DOUBLE PRECISION SUM
      DOUBLE PRECISION BEE(IBEE,*),DERIV(IDERIV,*),FUN(*),
     +                 COORD(ICOORD,*)
 
      SUM = 0.D0
      DO 1 K = 1,NOD
    1 SUM = SUM + FUN(K)*COORD(K,1)
      DO 2 M = 1,NOD
          N = 3*M
          K = N - 1
          L = K - 1
          BEE(1,L) = DERIV(1,M)
          BEE(2,K) = DERIV(2,M)
          BEE(3,L) = FUN(M)/SUM
          BEE(3,N) = IFLAG*LTH*BEE(3,L)
          BEE(4,L) = DERIV(2,M)
          BEE(4,K) = DERIV(1,M)
          BEE(5,K) = -IFLAG*LTH*FUN(M)/SUM
          BEE(5,N) = DERIV(2,M)
          BEE(6,L) = BEE(5,K)
    2 BEE(6,N) = DERIV(1,M) - FUN(M)/SUM
      RETURN
 
      END
      SUBROUTINE CHECON(LOADS,OLDLDS,N,TOL,ICON)
C
C      THIS SUBROUTINE SETS ICON TO ZERO IF THE RELATIVE CHANGE
C      IN VECTORS 'LOADS' AND 'OLDLDS' IS GREATER THAN 'TOL'
C
      DOUBLE PRECISION TOL
      DOUBLE PRECISION BIG
      DOUBLE PRECISION LOADS(*),OLDLDS(*)
 
      ICON = 1
      BIG = 0.D0
      DO 1 I = 1,N
    1 IF (ABS(LOADS(I)).GT.BIG) BIG = ABS(LOADS(I))
      DO 2 I = 1,N
          IF (ABS(LOADS(I)-OLDLDS(I))/BIG.GT.TOL) ICON = 0
    2 OLDLDS(I) = LOADS(I)
      RETURN
 
      END
      SUBROUTINE CHOBAC(KB,IKB,LOADS,N,IW)
C
C      THIS SUBROUTINE PERFORMS THE CHOLESKI BACK-SUBSTITUTION
C
      DOUBLE PRECISION X
      DOUBLE PRECISION KB(IKB,*),LOADS(*)
 
      LOADS(1) = LOADS(1)/KB(1,IW+1)
      DO 1 I = 2,N
          X = 0.0D0
          K = 1
          IF (I.LE.IW+1) K = IW - I + 2
          DO 2 J = K,IW
    2     X = X + KB(I,J)*LOADS(I+J-IW-1)
    1 LOADS(I) = (LOADS(I)-X)/KB(I,IW+1)
      LOADS(N) = LOADS(N)/KB(N,IW+1)
      I = N - 1
    3 X = 0.0D0
      L = I + IW
      IF (I.GT.N-IW) L = N
      M = I + 1
      DO 4 J = M,L
    4 X = X + KB(J,IW+I-J+1)*LOADS(J)
      LOADS(I) = (LOADS(I)-X)/KB(I,IW+1)
      I = I - 1
      IF (I) 5,5,3
    5 CONTINUE
      RETURN
 
      END
      SUBROUTINE CHOBK1(KB,IKB,LOADS,N,IW)
C
C      THIS SUBROUTINE PERFORMS CHOLESKI FORWARD SUBSTITUTION
C
      DOUBLE PRECISION X
      DOUBLE PRECISION KB(IKB,*),LOADS(*)
 
      LOADS(1) = LOADS(1)/KB(1,IW+1)
      DO 1 I = 2,N
          X = 0.D0
          K = 1
          IF (I.LE.IW+1) K = IW - I + 2
          DO 2 J = K,IW
    2     X = X + KB(I,J)*LOADS(I+J-IW-1)
    1 LOADS(I) = (LOADS(I)-X)/KB(I,IW+1)
      RETURN
 
      END
      SUBROUTINE CHOBK2(KB,IKB,LOADS,N,IW)
C
C      THIS SUBROUTINE PERFORMS BACKWARD CHOLESKI SUBSTITUTION
C      FOR RIGHT HAND SIDE 'LOADS'
C
      DOUBLE PRECISION X
      DOUBLE PRECISION KB(IKB,*),LOADS(*)
 
      LOADS(N) = LOADS(N)/KB(N,IW+1)
      I = N - 1
    1 X = 0.D0
      L = I + IW
      IF (I.GT.N-IW) L = N
      M = I + 1
      DO 2 J = M,L
    2 X = X + KB(J,IW+I-J+1)*LOADS(J)
      LOADS(I) = (LOADS(I)-X)/KB(I,IW+1)
      I = I - 1
      IF (I.NE.0) GO TO 1
      RETURN
 
      END
      SUBROUTINE CHOLIN(KB,IKB,N,IW)
C
C      THIS SUBROUTINE PERFORMS CHOLESKI REDUCTION OF
C      THE STIFFNESS MATRIX STORED AS AN ARRAY BK(N,IW+1)
C
      DOUBLE PRECISION X
      DOUBLE PRECISION KB(IKB,*)
 
      DO 1 I = 1,N
          X = 0.D0
          DO 2 J = 1,IW
    2     X = X + KB(I,J)**2
          KB(I,IW+1) = SQRT(KB(I,IW+1)-X)
          DO 3 K = 1,IW
              X = 0.D0
              IF (I+K.GT.N) GO TO 3
    6         IF (K.EQ.IW) GO TO 4
    7         L = IW - K
    5         X = X + KB(I+K,L)*KB(I,L+K)
              L = L - 1
              IF (L.NE.0) GO TO 5
    4         IA = I + K
              IB = IW - K + 1
              KB(IA,IB) = (KB(IA,IB)-X)/KB(I,IW+1)
    3     CONTINUE
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE COMBAC(BK,R,L,IW)
C
C      THIS SUBROUTINE PERFORMS THE BACK-SUBSTITUTION OF THE
C      COMPLEX STIFFNESS EQUATIONS
C
      COMPLEX*16 BK(*),R(*),SUM
 
      KB = IW + 1
      R(1) = R(1)/BK(1)
      DO 1 I = 2,L
          SUM = R(I)
          I1 = I - 1
          NKB = I - KB + 1
          IF (NKB) 2,2,3
    2     NKB = 1
    3     DO 4 K = NKB,I1
              JN = (I-K)*L + K
              SUM = SUM - BK(JN)*R(K)
    4     CONTINUE
          R(I) = SUM/BK(I)
    1 CONTINUE
      DO 5 JJ = 2,L
          I = L - JJ + 1
          SUM = 0.0D0
          I1 = I + 1
          NKB = I - 1 + KB
          IF (NKB-L) 7,7,6
    6     NKB = L
    7     DO 8 K = I1,NKB
              JN = (K-I)*L + I
    8     SUM = SUM + BK(JN)*R(K)
          R(I) = R(I) - SUM/BK(I)
    5 CONTINUE
      RETURN
 
      END
      SUBROUTINE COMRED(BK,L,IW)
C
C      THIS SUBROUTINE REDUCES THE COMPLEX STIFFNESS MATRIX
C
      COMPLEX*16 BK(*),SUM
 
      KB = IW + 1
      DO 1 I = 2,L
          IL1 = I - 1
          KBL = IL1 + KB
          IF (KBL-L) 3,3,2
    2     KBL = L
    3     DO 1 J = I,KBL
              IJ = (J-I)*L + I
              SUM = BK(IJ)
              NKB = J - KB + 1
              IF (NKB) 4,4,5
    4         NKB = 1
    5         IF (NKB-IL1) 6,6,8
    6         DO 7 N = NKB,IL1
                  NI = (I-N)*L + N
                  NJ = (J-N)*L + N
    7         SUM = SUM - BK(NI)*BK(NJ)/BK(N)
    8         BK(IJ) = SUM
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE ECMAT(ECM,IECM,TN,ITN,NT,INT,FUN,NOD,NODOF)
C
C      THIS SUBROUTINE FORMS THE CONSISTENT MASS MATRIX
C
      DOUBLE PRECISION X
      DOUBLE PRECISION ECM(IECM,*),TN(ITN,*),NT(INT,*),FUN(*)
 
      IDOF = NOD*NODOF
      DO 1 I = 1,IDOF
          DO 1 J = 1,NODOF
              NT(I,J) = 0.D0
    1 TN(J,I) = NT(I,J)
      DO 2 I = 1,NOD
          DO 2 J = 1,NODOF
              NT((I-1)*NODOF+J,J) = FUN(I)
    2 TN(J, (I-1)*NODOF+J) = FUN(I)
      DO 3 I = 1,IDOF
          DO 3 J = 1,IDOF
              X = 0.0D0
              DO 4 K = 1,NODOF
    4         X = X + NT(I,K)*TN(K,J)
              ECM(I,J) = X
    3 CONTINUE
      RETURN
 
      END
      SUBROUTINE EVECTS(N,ACHEPS,D,E,Z,IZ,IFAIL)
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS OF A
C     TRIDIAGONAL MATRIX, T, GIVEN WITH ITS DIAGONAL ELEMENTS IN
C     THE ARRAY D(N) AND ITS SUB-DIAGONAL ELEMENTS IN THE LAST N
C     - 1 STORES OF THE ARRAY E(N), USING QL TRANSFORMATIONS. THE
C     EIGENVALUES ARE OVERWRITTEN ON THE DIAGONAL ELEMENTS IN THE
C     ARRAY D IN ASCENDING ORDER. THE EIGENVECTORS ARE FORMED IN
C     THE ARRAY Z(N,N), OVERWRITING THE ACCUMULATED
C     TRANSFORMATIONS AS SUPPLIED BY THE SUBROUTINE F01AJF. THE
C     SUBROUTINE WILL FAIL IF ANY ONE EIGENVALUE TAKES MORE THAN 30
C     ITERATIONS.
C
      INTEGER P01AAF,ISAVE,IFAIL,N,I,L,J,M,I1,M1,II,K,IZ
      DOUBLE PRECISION B,F,H,ACHEPS,G,P,R,C,S,D(*),E(*),Z(IZ,*)
 
      ISAVE = IFAIL
      IF (N.EQ.1) GO TO 40
      DO 20 I = 2,N
          E(I-1) = E(I)
   20 CONTINUE
   40 E(N) = 0.0D0
      B = 0.0D0
      F = 0.0D0
      DO 300 L = 1,N
          J = 0
          H = ACHEPS* (ABS(D(L))+ABS(E(L)))
          IF (B.LT.H) B = H
C     LOOK FOR SMALL SUB-DIAG ELEMENT
          DO 60 M = L,N
              IF (ABS(E(M)).LE.B) GO TO 80
   60     CONTINUE
   80     IF (M.EQ.L) GO TO 280
  100     IF (J.EQ.30) GO TO 400
          J = J + 1
C     FORM SHIFT
          G = D(L)
          H = D(L+1) - G
          IF (ABS(H).GE.ABS(E(L))) GO TO 120
          P = H*0.5D0/E(L)
          R = SQRT(P*P+1.0D0)
          H = P + R
          IF (P.LT.0.0D0) H = P - R
          D(L) = E(L)/H
          GO TO 140
 
  120     P = 2.0D0*E(L)/H
          R = SQRT(P*P+1.0D0)
          D(L) = E(L)*P/ (1.0D0+R)
  140     H = G - D(L)
          I1 = L + 1
          IF (I1.GT.N) GO TO 180
          DO 160 I = I1,N
              D(I) = D(I) - H
  160     CONTINUE
  180     F = F + H
C     QL TRANSFORMATION
          P = D(M)
          C = 1.0D0
          S = 0.0D0
          M1 = M - 1
          DO 260 II = L,M1
              I = M1 - II + L
              G = C*E(I)
              H = C*P
              IF (ABS(P).LT.ABS(E(I))) GO TO 200
              C = E(I)/P
              R = SQRT(C*C+1.0D0)
              E(I+1) = S*P*R
              S = C/R
              C = 1.0D0/R
              GO TO 220
 
  200         C = P/E(I)
              R = SQRT(C*C+1.0D0)
              E(I+1) = S*E(I)*R
              S = 1.0D0/R
              C = C/R
  220         P = C*D(I) - S*G
              D(I+1) = H + S* (C*G+S*D(I))
C     FORM VECTOR
              DO 240 K = 1,N
                  H = Z(K,I+1)
                  Z(K,I+1) = S*Z(K,I) + C*H
                  Z(K,I) = C*Z(K,I) - S*H
  240         CONTINUE
  260     CONTINUE
          E(L) = S*P
          D(L) = C*P
          IF (ABS(E(L)).GT.B) GO TO 100
  280     D(L) = D(L) + F
  300 CONTINUE
C     ORDER EIGENVALUES AND EIGENVECTORS
      DO 380 I = 1,N
          K = I
          P = D(I)
          I1 = I + 1
          IF (I1.GT.N) GO TO 340
          DO 320 J = I1,N
              IF (D(J).GE.P) GO TO 320
              K = J
              P = D(J)
  320     CONTINUE
  340     IF (K.EQ.I) GO TO 380
          D(K) = D(I)
          D(I) = P
          DO 360 J = 1,N
              P = Z(J,I)
              Z(J,I) = Z(J,K)
              Z(J,K) = P
  360     CONTINUE
  380 CONTINUE
      IFAIL = 1
      RETURN
 
  400 IFAIL = 0
      RETURN
 
      END
      SUBROUTINE FDIAGV(BK,KM,IKM,G,KDIAG,IDOF)
C
C      THIS SUBROUTINE ASSEMBLES A DIAGONAL ELEMENT MATRIX
C      INTO THE GLOBAL SYSTEM
C
      DOUBLE PRECISION BK(*),KM(IKM,*)
      INTEGER G(*),KDIAG(*)
 
      DO 1 I = 1,IDOF
          J = G(I)
          IF (J.EQ.0) GO TO 1
          K = KDIAG(J)
          BK(K) = BK(K) + KM(I,I)
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FKDIAG(KDIAG,G,IDOF)
C
C      THIS SUBROUTINE FINDS THE MAXIMUM BANDWIDTH
C      FOR EACH FREEDOM
C
      INTEGER KDIAG(*),G(*)
 
      DO 1 I = 1,IDOF
          IWP1 = 1
          IF (G(I).EQ.0) GO TO 1
          DO 2 J = 1,IDOF
              IF (G(J).EQ.0) GO TO 2
              IM = G(I) - G(J) + 1
              IF (IM.GT.IWP1) IWP1 = IM
    2     CONTINUE
          K = G(I)
          IF (IWP1.GT.KDIAG(K)) KDIAG(K) = IWP1
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FMBEAM(DER2,FUN,SAMP,ISAMP,ELL,I)
C
C      THIS SUBROUTINE FORMS THE BEAM SHAPE FUNCTIONS
C      AND THEIR 2ND DERIVATIVES IN LOCAL COORDINATES
C
      DOUBLE PRECISION ELL
      DOUBLE PRECISION XI
      DOUBLE PRECISION XI2
      DOUBLE PRECISION XI3
      DOUBLE PRECISION DER2(*),FUN(*),SAMP(ISAMP,*)
 
      XI = SAMP(I,1)
      XI2 = XI*XI
      XI3 = XI2*XI
      FUN(1) = .25D0* (XI3-3.D0*XI+2.D0)
      FUN(2) = .125D0*ELL* (XI3-XI2-XI+1.D0)
      FUN(3) = .25D0* (-XI3+3.D0*XI+2.D0)
      FUN(4) = .125D0*ELL* (XI3+XI2-XI-1.D0)
      DER2(1) = 1.5D0*XI
      DER2(2) = .25D0*ELL* (3.D0*XI-1.D0)
      DER2(3) = -1.5D0*XI
      DER2(4) = .25D0*ELL* (3.D0*XI+1.D0)
      RETURN
 
      END
      SUBROUTINE FMBIGK(BIGK,IBIGK,KM,IKM,G,IDOF)
C
C       THIS SUBROUTINE ASSEMBLES ELEMENT MATRICES INTO A
C       FULL GLOBAL MATRIX
C
      DOUBLE PRECISION BIGK(IBIGK,*),KM(IKM,*)
      INTEGER G(*)
 
      DO 1 I = 1,IDOF
          IF (G(I).EQ.0) GO TO 1
    2     DO 3 J = 1,IDOF
              IF (G(J).EQ.0) GO TO 3
    4         BIGK(G(I),G(J)) = BIGK(G(I),G(J)) + KM(I,J)
    3     CONTINUE
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FMBRAD(BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)
C
C      THIS SUBROUTINE FORMS THE STRAIN/DISPLACEMENT MATRIX
C      FOR AXISYMMETRIC STRAIN
C
      DOUBLE PRECISION SUM
      DOUBLE PRECISION X
      DOUBLE PRECISION Y
      DOUBLE PRECISION BEE(IBEE,*),DERIV(IDERIV,*),FUN(*),
     +                 COORD(ICOORD,*)
 
      SUM = 0.D0
      DO 1 I = 1,NOD
    1 SUM = SUM + FUN(I)*COORD(I,1)
      DO 2 M = 1,NOD
          K = 2*M
          L = K - 1
          X = DERIV(1,M)
          BEE(1,L) = X
          BEE(3,K) = X
          Y = DERIV(2,M)
          BEE(2,K) = Y
          BEE(3,L) = Y
          BEE(4,L) = FUN(M)/SUM
    2 CONTINUE
      RETURN
 
      END
      SUBROUTINE FMDEPS(DEE,IDEE,E,V)
C
C      THIS SUBROUTINE FORMS THE ELASTIC PLANE STRAIN
C      STRESS/STRAIN MATRIX
C
      DOUBLE PRECISION E
      DOUBLE PRECISION V
      DOUBLE PRECISION V1
      DOUBLE PRECISION C
      DOUBLE PRECISION DEE(IDEE,*)
 
      V1 = 1.D0 - V
      C = E/ ((1.D0+V)* (1.D0-2.D0*V))
      DEE(1,1) = V1*C
      DEE(2,2) = V1*C
      DEE(3,3) = .5D0*C* (1.D0-2.D0*V)
      DEE(1,2) = V*C
      DEE(2,1) = V*C
      DEE(1,3) = 0.D0
      DEE(3,1) = 0.D0
      DEE(2,3) = 0.D0
      DEE(3,2) = 0.D0
      RETURN
 
      END
      SUBROUTINE FMDRAD(DEE,IDEE,E,V)
C
C      THIS SUBROUTINE FORMS THE ELASTIC AXISYMMETRIC
C      STRESS/STRAIN MATRIX
C
      DOUBLE PRECISION E
      DOUBLE PRECISION V
      DOUBLE PRECISION V1
      DOUBLE PRECISION C
      DOUBLE PRECISION DEE(IDEE,*)
 
      V1 = 1.D0 - V
      C = E/ ((1.D0+V)* (1.D0-2.D0*V))
      DEE(1,1) = V1*C
      DEE(2,2) = V1*C
      DEE(3,3) = .5D0*C* (1.D0-2.D0*V)
      DEE(4,4) = V1*C
      DEE(1,2) = V*C
      DEE(2,1) = V*C
      DEE(1,3) = 0.D0
      DEE(3,1) = 0.D0
      DEE(1,4) = V*C
      DEE(4,1) = V*C
      DEE(2,3) = 0.D0
      DEE(3,2) = 0.D0
      DEE(2,4) = V*C
      DEE(4,2) = V*C
      DEE(4,3) = 0.D0
      DEE(3,4) = 0.D0
      RETURN
 
      END
      SUBROUTINE FMDSIG(DEE,IDEE,E,V)
C
C      THIS SUBROUTINE FORMS THE ELASTIC PLANE STRESS
C      STRESS/STRAIN MATRIX
C
      DOUBLE PRECISION E
      DOUBLE PRECISION V
      DOUBLE PRECISION C
      DOUBLE PRECISION DEE(IDEE,*)
 
      C = E/ (1.D0-V*V)
      DEE(1,1) = C
      DEE(2,2) = C
      DEE(3,3) = .5D0*C* (1.D0-V)
      DEE(1,2) = V*C
      DEE(2,1) = V*C
      DEE(1,3) = 0.D0
      DEE(3,1) = 0.D0
      DEE(3,2) = 0.D0
      DEE(2,3) = 0.D0
      RETURN
 
      END
      SUBROUTINE FMKDKE(KM,IKM,KP,IKP,C,IC,KE,IKE,KD,IKD,IDOF,NOD,ITOT,
     +                  THETA)
C
C      THIS SUBROUTINE FORMS THE ELEMENT COUPLED STIFFNESS
C      MATRICES KE AND KD FROM THE ELASTIC STIFFNESS KM,
C      THE FLUID 'STIFFNESS' KP AND COUPLING MATRIX C
C
      DOUBLE PRECISION THETA
      DOUBLE PRECISION KM(IKM,*),KP(IKP,*),C(IC,*),KE(IKE,*),KD(IKD,*)
 
      DO 11 I = 1,IDOF
          DO 12 J = 1,IDOF
   12     KE(I,J) = KM(I,J)
          DO 13 K = 1,NOD
              KE(I,IDOF+K) = C(I,K)
   13     KE(IDOF+K,I) = C(I,K)
   11 CONTINUE
      DO 14 I = 1,NOD
          DO 16 K = 1,NOD
   16     KE(IDOF+I,IDOF+K) = KP(I,K)
   14 CONTINUE
      DO 17 I = 1,IDOF
          DO 17 J = 1,ITOT
              KD(I,J) = (THETA-1.D0)*KE(I,J)
   17 KE(I,J) = THETA*KE(I,J)
      M = IDOF + 1
      DO 18 I = M,ITOT
          DO 18 J = 1,IDOF
              KD(I,J) = KE(I,J)*THETA
   18 KE(I,J) = KD(I,J)
      DO 19 I = M,ITOT
          DO 19 J = M,ITOT
              KD(I,J) = KE(I,J)* (1.D0-THETA)*THETA
   19 KE(I,J) = -KE(I,J)*THETA*THETA
      RETURN
 
      END
      SUBROUTINE FMLAG9(DER,IDER,FUN,SAMP,ISAMP,I,J)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 9-NODED QUADRILATERAL ELEMENTS
C
      DOUBLE PRECISION ETA
      DOUBLE PRECISION XI
      DOUBLE PRECISION ETAM
      DOUBLE PRECISION XIM
      DOUBLE PRECISION ETAP
      DOUBLE PRECISION XIP
      DOUBLE PRECISION X2P1
      DOUBLE PRECISION X2M1
      DOUBLE PRECISION E2P1
      DOUBLE PRECISION E2M1
      DOUBLE PRECISION DER(IDER,*),FUN(*),SAMP(ISAMP,*)
 
      ETA = SAMP(I,1)
      XI = SAMP(J,1)
      ETAM = ETA - 1.D0
      XIM = XI - 1.D0
      ETAP = ETA + 1.D0
      XIP = XI + 1.D0
      X2P1 = 2.D0*XI + 1.D0
      X2M1 = 2.D0*XI - 1.D0
      E2P1 = 2.D0*ETA + 1.D0
      E2M1 = 2.D0*ETA - 1.D0
      FUN(1) = .25D0*XI*XIM*ETA*ETAM
      FUN(2) = -.5D0*XI*XIM*ETAP*ETAM
      FUN(3) = .25D0*XI*XIM*ETA*ETAP
      FUN(4) = -.5D0*XIP*XIM*ETA*ETAP
      FUN(5) = .25D0*XI*XIP*ETA*ETAP
      FUN(6) = -.5D0*XI*XIP*ETAP*ETAM
      FUN(7) = .25D0*XI*XIP*ETA*ETAM
      FUN(8) = -.5D0*XIP*XIM*ETA*ETAM
      FUN(9) = XIP*XIM*ETAP*ETAM
      DER(1,1) = .25D0*X2M1*ETA*ETAM
      DER(1,2) = -.5D0*X2M1*ETAP*ETAM
      DER(1,3) = .25D0*X2M1*ETA*ETAP
      DER(1,4) = -XI*ETA*ETAP
      DER(1,5) = .25D0*X2P1*ETA*ETAP
      DER(1,6) = -.5D0*X2P1*ETAP*ETAM
      DER(1,7) = .25D0*X2P1*ETA*ETAM
      DER(1,8) = -XI*ETA*ETAM
      DER(1,9) = 2.D0*XI*ETAP*ETAM
      DER(2,1) = .25D0*XI*XIM*E2M1
      DER(2,2) = -XI*XIM*ETA
      DER(2,3) = .25D0*XI*XIM*E2P1
      DER(2,4) = -.5D0*XIP*XIM*E2P1
      DER(2,5) = .25D0*XI*XIP*E2P1
      DER(2,6) = -XI*XIP*ETA
      DER(2,7) = .25D0*XI*XIP*E2M1
      DER(2,8) = -.5D0*XIP*XIM*E2M1
      DER(2,9) = 2.D0*XIP*XIM*ETA
      RETURN
 
      END
      SUBROUTINE FMLIN3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES FOR 8-NODED BRICK ELEMENTS
C
      DOUBLE PRECISION ETA
      DOUBLE PRECISION XI
      DOUBLE PRECISION ZETA
      DOUBLE PRECISION ETAM
      DOUBLE PRECISION XIM
      DOUBLE PRECISION ZETAM
      DOUBLE PRECISION ETAP
      DOUBLE PRECISION XIP
      DOUBLE PRECISION ZETAP
      DOUBLE PRECISION DER(IDER,*),FUN(*),SAMP(ISAMP,*)
 
      ETA = SAMP(I,1)
      XI = SAMP(J,1)
      ZETA = SAMP(K,1)
      ETAM = 1.D0 - ETA
      XIM = 1.D0 - XI
      ZETAM = 1.D0 - ZETA
      ETAP = ETA + 1.D0
      XIP = XI + 1.D0
      ZETAP = ZETA + 1.D0
      FUN(1) = .125D0*XIM*ETAM*ZETAM
      FUN(2) = .125D0*XIM*ETAM*ZETAP
      FUN(3) = .125D0*XIP*ETAM*ZETAP
      FUN(4) = .125D0*XIP*ETAM*ZETAM
      FUN(5) = .125D0*XIM*ETAP*ZETAM
      FUN(6) = .125D0*XIM*ETAP*ZETAP
      FUN(7) = .125D0*XIP*ETAP*ZETAP
      FUN(8) = .125D0*XIP*ETAP*ZETAM
      DER(1,1) = -.125D0*ETAM*ZETAM
      DER(1,2) = -.125D0*ETAM*ZETAP
      DER(1,3) = .125D0*ETAM*ZETAP
      DER(1,4) = .125D0*ETAM*ZETAM
      DER(1,5) = -.125D0*ETAP*ZETAM
      DER(1,6) = -.125D0*ETAP*ZETAP
      DER(1,7) = .125D0*ETAP*ZETAP
      DER(1,8) = .125D0*ETAP*ZETAM
      DER(2,1) = -.125D0*XIM*ZETAM
      DER(2,2) = -.125D0*XIM*ZETAP
      DER(2,3) = -.125D0*XIP*ZETAP
      DER(2,4) = -.125D0*XIP*ZETAM
      DER(2,5) = .125D0*XIM*ZETAM
      DER(2,6) = .125D0*XIM*ZETAP
      DER(2,7) = .125D0*XIP*ZETAP
      DER(2,8) = .125D0*XIP*ZETAM
      DER(3,1) = -.125D0*XIM*ETAM
      DER(3,2) = .125D0*XIM*ETAM
      DER(3,3) = .125D0*XIP*ETAM
      DER(3,4) = -.125D0*XIP*ETAM
      DER(3,5) = -.125D0*XIM*ETAP
      DER(3,6) = .125D0*XIM*ETAP
      DER(3,7) = .125D0*XIP*ETAP
      DER(3,8) = -.125D0*XIP*ETAP
      RETURN
 
      END
      SUBROUTINE FMLUMP(DIAG,EMM,IEMM,G,IDOF)
C
C      THIS SUBROUTINE FORMS THE GLOBAL MASS MATRIX
C      AS VECTOR DIAG
C
      DOUBLE PRECISION DIAG(*),EMM(IEMM,*)
      INTEGER G(*)
 
      DO 1 I = 1,IDOF
          IF (G(I).EQ.0) GO TO 1
          DIAG(G(I)) = DIAG(G(I)) + EMM(I,I)
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FMPLAT(FUN,D1X,D1Y,D2X,D2Y,D2XY,SAMP,ISAMP,AA,BB,I,J)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND THEIR 1ST
C      AND 2ND DERIVATIVES FOR RECTANGULAR PLATE BENDING ELEMENTS
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION X
      DOUBLE PRECISION E
      DOUBLE PRECISION XP1
      DOUBLE PRECISION XP12
      DOUBLE PRECISION XP13
      DOUBLE PRECISION EP1
      DOUBLE PRECISION EP12
      DOUBLE PRECISION EP13
      DOUBLE PRECISION P1
      DOUBLE PRECISION Q1
      DOUBLE PRECISION P2
      DOUBLE PRECISION Q2
      DOUBLE PRECISION P3
      DOUBLE PRECISION Q3
      DOUBLE PRECISION P4
      DOUBLE PRECISION Q4
      DOUBLE PRECISION DP1
      DOUBLE PRECISION DQ1
      DOUBLE PRECISION DP2
      DOUBLE PRECISION DQ2
      DOUBLE PRECISION DP3
      DOUBLE PRECISION DQ3
      DOUBLE PRECISION DP4
      DOUBLE PRECISION DQ4
      DOUBLE PRECISION D2P1
      DOUBLE PRECISION D2P2
      DOUBLE PRECISION D2P3
      DOUBLE PRECISION D2P4
      DOUBLE PRECISION D2Q1
      DOUBLE PRECISION D2Q2
      DOUBLE PRECISION D2Q3
      DOUBLE PRECISION D2Q4
      DOUBLE PRECISION FUN(*),D1X(*),D1Y(*),D2X(*),D2Y(*),D2XY(*),
     +                 SAMP(ISAMP,*)
 
      X = SAMP(I,1)
      E = SAMP(J,1)
      XP1 = X + 1.D0
      XP12 = XP1*XP1
      XP13 = XP12*XP1
      EP1 = E + 1.D0
      EP12 = EP1*EP1
      EP13 = EP12*EP1
      P1 = 1.D0 - .75D0*XP12 + .25D0*XP13
      Q1 = 1.D0 - .75D0*EP12 + .25D0*EP13
      P2 = .5D0*AA*XP1* (1.D0-XP1+.25D0*XP12)
      Q2 = .5D0*BB*EP1* (1.D0-EP1+.25D0*EP12)
      P3 = .25D0*XP12* (3.D0-XP1)
      Q3 = .25D0*EP12* (3.D0-EP1)
      P4 = .25D0*AA*XP12* (.5D0*XP1-1.D0)
      Q4 = .25D0*BB*EP12* (.5D0*EP1-1.D0)
      FUN(1) = P1*Q1
      FUN(2) = P2*Q1
      FUN(3) = P1*Q2
      FUN(4) = P2*Q2
      FUN(5) = P1*Q3
      FUN(6) = P2*Q3
      FUN(7) = P1*Q4
      FUN(8) = P2*Q4
      FUN(9) = P3*Q3
      FUN(10) = P4*Q3
      FUN(11) = P3*Q4
      FUN(12) = P4*Q4
      FUN(13) = P3*Q1
      FUN(14) = P4*Q1
      FUN(15) = P3*Q2
      FUN(16) = P4*Q2
      DP1 = 1.5D0*XP1* (.5D0*XP1-1.D0)
      DQ1 = 1.5D0*EP1* (.5D0*EP1-1.D0)
      DP2 = AA* (.5D0-XP1+.375D0*XP12)
      DQ2 = BB* (.5D0-EP1+.375D0*EP12)
      DP3 = 1.5D0*XP1* (1.D0-.5D0*XP1)
      DQ3 = 1.5D0*EP1* (1.D0-.5D0*EP1)
      DP4 = .5D0*AA*XP1* (.75D0*XP1-1.D0)
      DQ4 = .5D0*BB*EP1* (.75D0*EP1-1.D0)
      D2P1 = 1.5D0*X
      D2P2 = .25D0*AA* (3.D0*X-1.D0)
      D2P3 = -D2P1
      D2P4 = .25D0*AA* (3.D0*X+1.D0)
      D2Q1 = 1.5D0*E
      D2Q2 = .25D0*BB* (3.D0*E-1.D0)
      D2Q3 = -D2Q1
      D2Q4 = .25D0*BB* (3.D0*E+1.D0)
      D1X(1) = DP1*Q1
      D1X(2) = DP2*Q1
      D1X(3) = DP1*Q2
      D1X(4) = DP2*Q2
      D1X(5) = DP1*Q3
      D1X(6) = DP2*Q3
      D1X(7) = DP1*Q4
      D1X(8) = DP2*Q4
      D1X(9) = DP3*Q3
      D1X(10) = DP4*Q3
      D1X(11) = DP3*Q4
      D1X(12) = DP4*Q4
      D1X(13) = DP3*Q1
      D1X(14) = DP4*Q1
      D1X(15) = DP3*Q2
      D1X(16) = DP4*Q2
      D1Y(1) = P1*DQ1
      D1Y(2) = P2*DQ1
      D1Y(3) = P1*DQ2
      D1Y(4) = P2*DQ2
      D1Y(5) = P1*DQ3
      D1Y(6) = P2*DQ3
      D1Y(7) = P1*DQ4
      D1Y(8) = P2*DQ4
      D1Y(9) = P3*DQ3
      D1Y(10) = P4*DQ3
      D1Y(11) = P3*DQ4
      D1Y(12) = P4*DQ4
      D1Y(13) = P3*DQ1
      D1Y(14) = P4*DQ1
      D1Y(15) = P3*DQ2
      D1Y(16) = P4*DQ2
      D2X(1) = D2P1*Q1
      D2X(2) = D2P2*Q1
      D2X(3) = D2P1*Q2
      D2X(4) = D2P2*Q2
      D2X(5) = D2P1*Q3
      D2X(6) = D2P2*Q3
      D2X(7) = D2P1*Q4
      D2X(8) = D2P2*Q4
      D2X(9) = D2P3*Q3
      D2X(10) = D2P4*Q3
      D2X(11) = D2P3*Q4
      D2X(12) = D2P4*Q4
      D2X(13) = D2P3*Q1
      D2X(14) = D2P4*Q1
      D2X(15) = D2P3*Q2
      D2X(16) = D2P4*Q2
      D2Y(1) = P1*D2Q1
      D2Y(2) = P2*D2Q1
      D2Y(3) = P1*D2Q2
      D2Y(4) = P2*D2Q2
      D2Y(5) = P1*D2Q3
      D2Y(6) = P2*D2Q3
      D2Y(7) = P1*D2Q4
      D2Y(8) = P2*D2Q4
      D2Y(9) = P3*D2Q3
      D2Y(10) = P4*D2Q3
      D2Y(11) = P3*D2Q4
      D2Y(12) = P4*D2Q4
      D2Y(13) = P3*D2Q1
      D2Y(14) = P4*D2Q1
      D2Y(15) = P3*D2Q2
      D2Y(16) = P4*D2Q2
      D2XY(1) = DP1*DQ1
      D2XY(2) = DP2*DQ1
      D2XY(3) = DP1*DQ2
      D2XY(4) = DP2*DQ2
      D2XY(5) = DP1*DQ3
      D2XY(6) = DP2*DQ3
      D2XY(7) = DP1*DQ4
      D2XY(8) = DP2*DQ4
      D2XY(9) = DP3*DQ3
      D2XY(10) = DP4*DQ3
      D2XY(11) = DP3*DQ4
      D2XY(12) = DP4*DQ4
      D2XY(13) = DP3*DQ1
      D2XY(14) = DP4*DQ1
      D2XY(15) = DP3*DQ2
      D2XY(16) = DP4*DQ2
      RETURN
 
      END
      SUBROUTINE FMQUAD(DER,IDER,FUN,SAMP,ISAMP,I,J)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 8-NODED QUADRILATERAL ELEMENTS
C
      DOUBLE PRECISION ETA
      DOUBLE PRECISION XI
      DOUBLE PRECISION ETAM
      DOUBLE PRECISION ETAP
      DOUBLE PRECISION XIM
      DOUBLE PRECISION XIP
      DOUBLE PRECISION DER(IDER,*),FUN(*),SAMP(ISAMP,*)
 
      ETA = SAMP(I,1)
      XI = SAMP(J,1)
      ETAM = .25D0* (1.D0-ETA)
      ETAP = .25D0* (1.D0+ETA)
      XIM = .25D0* (1.D0-XI)
      XIP = .25D0* (1.D0+XI)
      FUN(1) = 4.D0*ETAM*XIM* (-XI-ETA-1.D0)
      FUN(2) = 32.D0*ETAM*XIM*ETAP
      FUN(3) = 4.D0*ETAP*XIM* (-XI+ETA-1.D0)
      FUN(4) = 32.D0*XIM*XIP*ETAP
      FUN(5) = 4.D0*ETAP*XIP* (XI+ETA-1.D0)
      FUN(6) = 32.D0*ETAP*XIP*ETAM
      FUN(7) = 4.D0*XIP*ETAM* (XI-ETA-1.D0)
      FUN(8) = 32.D0*XIM*XIP*ETAM
      DER(1,1) = ETAM* (2.D0*XI+ETA)
      DER(1,2) = -8.D0*ETAM*ETAP
      DER(1,3) = ETAP* (2.D0*XI-ETA)
      DER(1,4) = -4.D0*ETAP*XI
      DER(1,5) = ETAP* (2.D0*XI+ETA)
      DER(1,6) = 8.D0*ETAP*ETAM
      DER(1,7) = ETAM* (2.D0*XI-ETA)
      DER(1,8) = -4.D0*ETAM*XI
      DER(2,1) = XIM* (XI+2.D0*ETA)
      DER(2,2) = -4.D0*XIM*ETA
      DER(2,3) = XIM* (2.D0*ETA-XI)
      DER(2,4) = 8.D0*XIM*XIP
      DER(2,5) = XIP* (XI+2.D0*ETA)
      DER(2,6) = -4.D0*XIP*ETA
      DER(2,7) = XIP* (2.D0*ETA-XI)
      DER(2,8) = -8.D0*XIM*XIP
      RETURN
 
      END
      SUBROUTINE FMQUA3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES FOR 20-NODED BRICK ELEMENTS
C
      DOUBLE PRECISION XI
      DOUBLE PRECISION ETA
      DOUBLE PRECISION ZETA
      DOUBLE PRECISION XIO
      DOUBLE PRECISION ETAO
      DOUBLE PRECISION ZETAO
      DOUBLE PRECISION DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      INTEGER XII(20),ETAI(20),ZETAI(20)
 
      XI = SAMP(I,1)
      ETA = SAMP(J,1)
      ZETA = SAMP(K,1)
      XII(1) = -1
      XII(2) = -1
      XII(3) = -1
      XII(9) = -1
      XII(10) = -1
      XII(13) = -1
      XII(14) = -1
      XII(15) = -1
      XII(4) = 0
      XII(8) = 0
      XII(16) = 0
      XII(20) = 0
      XII(5) = 1
      XII(6) = 1
      XII(7) = 1
      XII(11) = 1
      XII(12) = 1
      XII(17) = 1
      XII(18) = 1
      XII(19) = 1
      DO 1 L = 1,8
    1 ETAI(L) = -1
      DO 2 L = 9,12
    2 ETAI(L) = 0
      DO 3 L = 13,20
    3 ETAI(L) = 1
      ZETAI(1) = -1
      ZETAI(7) = -1
      ZETAI(8) = -1
      ZETAI(9) = -1
      ZETAI(12) = -1
      ZETAI(13) = -1
      ZETAI(19) = -1
      ZETAI(20) = -1
      ZETAI(2) = 0
      ZETAI(6) = 0
      ZETAI(14) = 0
      ZETAI(18) = 0
      ZETAI(3) = 1
      ZETAI(4) = 1
      ZETAI(5) = 1
      ZETAI(10) = 1
      ZETAI(11) = 1
      ZETAI(15) = 1
      ZETAI(16) = 1
      ZETAI(17) = 1
      DO 4 L = 1,20
          XIO = XI*XII(L)
          ETAO = ETA*ETAI(L)
          ZETAO = ZETA*ZETAI(L)
          IF (L.EQ.4 .OR. L.EQ.8 .OR. L.EQ.16 .OR. L.EQ.20) THEN
              FUN(L) = .25D0* (1.D0-XI*XI)* (1.D0+ETAO)* (1.D0+ZETAO)
              DER(1,L) = -.5D0*XI* (1.D0+ETAO)* (1.D0+ZETAO)
              DER(2,L) = .25D0*ETAI(L)* (1.D0-XI*XI)* (1.D0+ZETAO)
              DER(3,L) = .25D0*ZETAI(L)* (1.D0-XI*XI)* (1.D0+ETAO)
 
          ELSE IF (L.GE.9 .AND. L.LE.12) THEN
              FUN(L) = .25D0* (1.D0+XIO)* (1.D0-ETA*ETA)* (1.D0+ZETAO)
              DER(1,L) = .25D0*XII(L)* (1.D0-ETA*ETA)* (1.D0+ZETAO)
              DER(2,L) = -.5D0*ETA* (1.D0+XIO)* (1.D0+ZETAO)
              DER(3,L) = .25D0*ZETAI(L)* (1.D0+XIO)* (1.D0-ETA*ETA)
 
          ELSE IF (L.EQ.2 .OR. L.EQ.6 .OR. L.EQ.14 .OR. L.EQ.18) THEN
              FUN(L) = .25D0* (1.D0+XIO)* (1.D0+ETAO)* (1.D0-ZETA*ZETA)
              DER(1,L) = .25D0*XII(L)* (1.D0+ETAO)* (1.D0-ZETA*ZETA)
              DER(2,L) = .25D0*ETAI(L)* (1.D0+XIO)* (1.D0-ZETA*ZETA)
              DER(3,L) = -.5D0*ZETA* (1.D0+XIO)* (1.D0+ETAO)
 
          ELSE
              FUN(L) = .125D0* (1.D0+XIO)* (1.D0+ETAO)* (1.D0+ZETAO)*
     +                 (XIO+ETAO+ZETAO-2.D0)
              DER(1,L) = .125D0*XII(L)* (1.D0+ETAO)* (1.D0+ZETAO)*
     +                   (2.D0*XIO+ETAO+ZETAO-1.D0)
              DER(2,L) = .125D0*ETAI(L)* (1.D0+XIO)* (1.D0+ZETAO)*
     +                   (XIO+2.D0*ETAO+ZETAO-1.D0)
              DER(3,L) = .125D0*ZETAI(L)* (1.D0+XIO)* (1.D0+ETAO)*
     +                   (XIO+ETAO+2.D0*ZETAO-1.D0)
          END IF
 
    4 CONTINUE
      RETURN
 
      END
      SUBROUTINE FMTET4(DER,IDER,FUN,SAMP,ISAMP,I)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 4-NODE TETRAHEDRON ELEMENTS
C
      DOUBLE PRECISION DER(IDER,*),SAMP(ISAMP,*),FUN(*)
 
      FUN(1) = SAMP(I,1)
      FUN(2) = SAMP(I,2)
      FUN(3) = SAMP(I,3)
      FUN(4) = 1.D0 - FUN(1) - FUN(2) - FUN(3)
      DO 1 M = 1,3
          DO 1 N = 1,4
    1 DER(M,N) = 0.D0
      DER(1,1) = 1.D0
      DER(2,2) = 1.D0
      DER(3,3) = 1.D0
      DER(1,4) = -1.D0
      DER(2,4) = -1.D0
      DER(3,4) = -1.D0
      RETURN
 
      END
      SUBROUTINE FMTRI3(DER,IDER,FUN,SAMP,ISAMP,I)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 3-NODED TRIANGULAR ELEMENTS
C
      DOUBLE PRECISION DER(IDER,*),SAMP(ISAMP,*),FUN(*)
 
      FUN(1) = SAMP(I,1)
      FUN(2) = SAMP(I,2)
      FUN(3) = 1.D0 - FUN(1) - FUN(2)
      DER(1,1) = 1.D0
      DER(1,2) = 0.D0
      DER(1,3) = -1.D0
      DER(2,1) = 0.D0
      DER(2,2) = 1.D0
      DER(2,3) = -1.D0
      RETURN
 
      END
      SUBROUTINE FMTRI6(DER,IDER,FUN,SAMP,ISAMP,I)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 6-NODED TRIANGULAR ELEMENTS
C
      DOUBLE PRECISION C1
      DOUBLE PRECISION C2
      DOUBLE PRECISION C3
      DOUBLE PRECISION DER(IDER,*),SAMP(ISAMP,*),FUN(*)
 
      C1 = SAMP(I,1)
      C2 = SAMP(I,2)
      C3 = 1.D0 - C1 - C2
      FUN(1) = (2.D0*C1-1.D0)*C1
      FUN(2) = 4.D0*C1*C2
      FUN(3) = (2.D0*C2-1.D0)*C2
      FUN(4) = 4.D0*C2*C3
      FUN(5) = (2.D0*C3-1.D0)*C3
      FUN(6) = 4.D0*C3*C1
      DER(1,1) = 4.D0*C1 - 1.D0
      DER(1,2) = 4.D0*C2
      DER(1,3) = 0.D0
      DER(1,4) = -4.D0*C2
      DER(1,5) = - (4.D0*C3-1.D0)
      DER(1,6) = 4.D0* (C3-C1)
      DER(2,1) = 0.D0
      DER(2,2) = 4.D0*C1
      DER(2,3) = 4.D0*C2 - 1.D0
      DER(2,4) = 4.D0* (C3-C2)
      DER(2,5) = - (4.D0*C3-1.D0)
      DER(2,6) = -4.D0*C1
      RETURN
 
      END
      SUBROUTINE FMTR15(DER,IDER,FUN,SAMP,ISAMP,I)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 15-NODED TRIANGULAR ELEMENTS
C
      DOUBLE PRECISION C1
      DOUBLE PRECISION C2
      DOUBLE PRECISION C3
      DOUBLE PRECISION T1
      DOUBLE PRECISION T2
      DOUBLE PRECISION T3
      DOUBLE PRECISION T4
      DOUBLE PRECISION T5
      DOUBLE PRECISION T6
      DOUBLE PRECISION T7
      DOUBLE PRECISION T8
      DOUBLE PRECISION T9
      DOUBLE PRECISION DER(IDER,*),SAMP(ISAMP,*),FUN(*)
 
      C1 = SAMP(I,1)
      C2 = SAMP(I,2)
      C3 = 1.D0 - C1 - C2
      T1 = C1 - .25D0
      T2 = C1 - .5D0
      T3 = C1 - .75D0
      T4 = C2 - .25D0
      T5 = C2 - .5D0
      T6 = C2 - .75D0
      T7 = C3 - .25D0
      T8 = C3 - .5D0
      T9 = C3 - .75D0
      FUN(1) = 32.D0/3.D0*C1*T1*T2*T3
      FUN(2) = 128.D0/3.D0*C1*C2*T1*T2
      FUN(3) = 64.D0*C1*C2*T1*T4
      FUN(4) = 128.D0/3.D0*C1*C2*T4*T5
      FUN(5) = 32.D0/3.D0*C2*T4*T5*T6
      FUN(6) = 128.D0/3.D0*C2*C3*T4*T5
      FUN(7) = 64.D0*C2*C3*T4*T7
      FUN(8) = 128.D0/3.D0*C2*C3*T7*T8
      FUN(9) = 32.D0/3.D0*C3*T7*T8*T9
      FUN(10) = 128.D0/3.D0*C3*C1*T7*T8
      FUN(11) = 64.D0*C3*C1*T1*T7
      FUN(12) = 128.D0/3.D0*C3*C1*T1*T2
      FUN(13) = 128.D0*C1*C2*T1*C3
      FUN(14) = 128.D0*C1*C2*C3*T4
      FUN(15) = 128.D0*C1*C2*C3*T7
      DER(1,1) = 32.D0/3.D0* (T2*T3* (T1+C1)+C1*T1* (T3+T2))
      DER(1,2) = 128.D0/3.D0*C2* (T2* (T1+C1)+C1*T1)
      DER(1,3) = 64.D0*C2*T4* (T1+C1)
      DER(1,4) = 128.D0/3.D0*C2*T4*T5
      DER(1,5) = 0.D0
      DER(1,6) = -128.D0/3.D0*C2*T4*T5
      DER(1,7) = -64.D0*C2*T4* (T7+C3)
      DER(1,8) = -128.D0/3.D0*C2* (T8* (T7+C3)+C3*T7)
      DER(1,9) = -32.D0/3.D0* (T8*T9* (T7+C3)+C3*T7* (T8+T9))
      DER(1,10) = 128.D0/3.D0* (C3*T7*T8-C1* (T8* (T7+C3)+C3*T7))
      DER(1,11) = 64.D0* (C3*T7* (T1+C1)-C1*T1* (T7+C3))
      DER(1,12) = 128.D0/3.D0* (C3* (T2* (T1+C1)+C1*T1)-C1*T1*T2)
      DER(1,13) = 128.D0*C2* (C3* (T1+C1)-C1*T1)
      DER(1,14) = 128.D0*C2*T4* (C3-C1)
      DER(1,15) = 128.D0*C2* (C3*T7-C1* (T7+C3))
      DER(2,1) = 0.0D0
      DER(2,2) = 128.D0/3.D0*C1*T1*T2
      DER(2,3) = 64.D0*C1*T1* (T4+C2)
      DER(2,4) = 128.D0/3.D0*C1* (T5* (T4+C2)+C2*T4)
      DER(2,5) = 32.D0/3.D0* (T5*T6* (T4+C2)+C2*T4* (T6+T5))
      DER(2,6) = 128.D0/3.D0* ((C3* (T5* (T4+C2)+C2*T4))-C2*T4*T5)
      DER(2,7) = 64.D0* (C3*T7* (T4+C2)-C2*T4* (T7+C3))
      DER(2,8) = 128.D0/3.D0* (C3*T7*T8-C2* (T8* (T7+C3)+C3*T7))
      DER(2,9) = -32.D0/3.D0* (T8*T9* (T7+C3)+C3*T7* (T8+T9))
      DER(2,10) = -128.D0/3.D0*C1* (T8* (T7+C3)+C3*T7)
      DER(2,11) = -64.D0*C1*T1* (T7+C3)
      DER(2,12) = -128.D0/3.D0*C1*T1*T2
      DER(2,13) = 128.D0*C1*T1* (C3-C2)
      DER(2,14) = 128.D0*C1* (C3* (T4+C2)-C2*T4)
      DER(2,15) = 128.D0*C1* (C3*T7-C2* (C3+T7))
      RETURN
 
      END
      SUBROUTINE FORMB(BEE,IBEE,DERIV,IDERIV,NOD)
C
C      THIS SUBROUTINE FORMS THE STRAIN/DISPLACEMENT MATRIX
C      FOR PLANE STRAIN
C
      DOUBLE PRECISION X
      DOUBLE PRECISION Y
      DOUBLE PRECISION BEE(IBEE,*),DERIV(IDERIV,*)
 
      DO 1 M = 1,NOD
          K = 2*M
          L = K - 1
          X = DERIV(1,M)
          BEE(1,L) = X
          BEE(3,K) = X
          Y = DERIV(2,M)
          BEE(2,K) = Y
          BEE(3,L) = Y
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FORMB3(BEE,IBEE,DERIV,IDERIV,NOD)
C
C      THIS SUBROUTINE FORMS THE 3-D STRAIN-DISPLACEMENT MATRIX
C
      DOUBLE PRECISION X
      DOUBLE PRECISION Y
      DOUBLE PRECISION Z
      DOUBLE PRECISION BEE(IBEE,*),DERIV(IDERIV,*)
 
      DO 1 M = 1,NOD
          N = 3*M
          K = N - 1
          L = K - 1
          X = DERIV(1,M)
          BEE(1,L) = X
          BEE(4,K) = X
          BEE(6,N) = X
          Y = DERIV(2,M)
          BEE(2,K) = Y
          BEE(4,L) = Y
          BEE(5,N) = Y
          Z = DERIV(3,M)
          BEE(3,N) = Z
          BEE(5,K) = Z
          BEE(6,L) = Z
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FORMD3(DEE,IDEE,E,V)
C
C      THIS SUBROUTINE FORMS THE 3-D STRAIN
C      STRESS/STRAIN MATRIX
C
      DOUBLE PRECISION E
      DOUBLE PRECISION V
      DOUBLE PRECISION V1
      DOUBLE PRECISION VV
      DOUBLE PRECISION DEE(IDEE,*)
 
      V1 = V/ (1.D0-V)
      VV = (1.D0-2.D0*V)/ (1.D0-V)*.5D0
      DO 1 I = 1,6
          DO 1 J = 1,6
    1 DEE(I,J) = 0.D0
      DEE(1,1) = 1.D0
      DEE(2,2) = 1.D0
      DEE(3,3) = 1.D0
      DEE(1,2) = V1
      DEE(2,1) = V1
      DEE(1,3) = V1
      DEE(3,1) = V1
      DEE(2,3) = V1
      DEE(3,2) = V1
      DEE(4,4) = VV
      DEE(5,5) = VV
      DEE(6,6) = VV
      DO 2 I = 1,6
          DO 2 J = 1,6
    2 DEE(I,J) = DEE(I,J)*E/ (2.D0* (1.D0+V)*VV)
      RETURN
 
      END
      SUBROUTINE FORMGP(IP,IQ,NYE,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE 'STEERING' VECTOR FOR A
C      4-NODE RECTANGULAR PLATE BENDING ELEMENT
C
      INTEGER G(*),NF(INF,*),NUM(4)
 
      I1 = (IP-1)* (NYE+1) + IQ
      I2 = I1 + 1
      I3 = IP* (NYE+1) + IQ
      I4 = I3 + 1
      DO 1 I = 1,4
          G(I) = NF(I1,I)
          G(I+4) = NF(I2,I)
          G(I+8) = NF(I4,I)
          G(I+12) = NF(I3,I)
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FORMKB(KB,IKB,KM,IKM,G,IW,IDOF)
C
C      THIS SUBROUTINE FORMS THE GLOBAL STIFFNESS MATRIX
C      STORING THE LOWER TRIANGLE AS AN ARRAY BK(N,IW+1)
C
      DOUBLE PRECISION KB(IKB,*),KM(IKM,*)
      INTEGER G(*)
 
      DO 1 I = 1,IDOF
          IF (G(I)) 1,1,2
    2     DO 3 J = 1,IDOF
              IF (G(J)) 3,3,4
    4         ICD = G(J) - G(I) + IW + 1
              IF (ICD-IW-1) 5,5,3
    5         KB(G(I),ICD) = KB(G(I),ICD) + KM(I,J)
    3     CONTINUE
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FORMKC(BK,KM,IKM,CM,ICM,G,N,IDOF)
C
C      THIS SUBROUTINE ASSEMBLES COMPLEX ELEMENT MATRICES INTO A
C      SYMMETRICAL BAND GLOBAL MATRIX(STORED AS A VECTOR)
C
      COMPLEX*16 BK(*)
      DOUBLE PRECISION KM(IKM,*),CM(ICM,*)
      INTEGER G(*)
 
      DO 1 I = 1,IDOF
          IF (G(I).EQ.0) GO TO 1
    2     DO 5 J = 1,IDOF
              IF (G(J).EQ.0) GO TO 5
    3         ICD = G(J) - G(I) + 1
              IF (ICD-1) 5,4,4
    4         IVAL = N* (ICD-1) + G(I)
              BK(IVAL) = BK(IVAL) + DCMPLX(KM(I,J),CM(I,J))
    5     CONTINUE
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FORMKU(KU,IKU,KM,IKM,G,IW,IDOF)
C
C      THIS SUBROUTINE ASSEMBLES ELEMENT MATRICES INTO SYMMETRICAL
C      GLOBAL MATRIX(STORED AS AN UPPER RECTANGLE)
C
      DOUBLE PRECISION KU(IKU,*),KM(IKM,*)
      INTEGER G(*)
 
      DO 1 I = 1,IDOF
          IF (G(I).EQ.0) GO TO 1
          DO 2 J = 1,IDOF
              IF (G(J).EQ.0) GO TO 2
              ICD = G(J) - G(I) + 1
              IF (ICD.LT.1) GO TO 2
              KU(G(I),ICD) = KU(G(I),ICD) + KM(I,J)
    2     CONTINUE
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FORMKV(BK,KM,IKM,G,N,IDOF)
C
C      THIS SUBROUTINE FORMS THE GLOBAL STIFFNESS MATRIX
C      STORING THE UPPER TRIANGLE AS A VECTOR BK(N*(IW+1))
C
      DOUBLE PRECISION BK(*),KM(IKM,*)
      INTEGER G(*)
 
      DO 1 I = 1,IDOF
          IF (G(I).EQ.0) GO TO 1
          DO 5 J = 1,IDOF
              IF (G(J).EQ.0) GO TO 5
              ICD = G(J) - G(I) + 1
              IF (ICD-1) 5,4,4
    4         IVAL = N* (ICD-1) + G(I)
              BK(IVAL) = BK(IVAL) + KM(I,J)
    5     CONTINUE
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FORMLN(DER,IDER,FUN,SAMP,ISAMP,I,J)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 4-NODED QUADRILATERAL ELEMENTS
C
      DOUBLE PRECISION ETA
      DOUBLE PRECISION XI
      DOUBLE PRECISION ETAM
      DOUBLE PRECISION ETAP
      DOUBLE PRECISION XIM
      DOUBLE PRECISION XIP
      DOUBLE PRECISION DER(IDER,*),FUN(*),SAMP(ISAMP,*)
 
      ETA = SAMP(I,1)
      XI = SAMP(J,1)
      ETAM = .25D0* (1.D0-ETA)
      ETAP = .25D0* (1.D0+ETA)
      XIM = .25D0* (1.D0-XI)
      XIP = .25D0* (1.D0+XI)
      FUN(1) = 4.D0*XIM*ETAM
      FUN(2) = 4.D0*XIM*ETAP
      FUN(3) = 4.D0*XIP*ETAP
      FUN(4) = 4.D0*XIP*ETAM
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
      SUBROUTINE FORMM(STRESS,M1,M2,M3)
C
C      THIS SUBROUTINE FORMS THE DERIVATIVES OF THE INVARIANTS
C      WITH RESPECT TO THE STRESSES
C
      DOUBLE PRECISION SX
      DOUBLE PRECISION SY
      DOUBLE PRECISION TXY
      DOUBLE PRECISION SZ
      DOUBLE PRECISION DX
      DOUBLE PRECISION DY
      DOUBLE PRECISION DZ
      DOUBLE PRECISION SIGM
      DOUBLE PRECISION STRESS(*),M1(4,4),M2(4,4),M3(4,4)
 
      SX = STRESS(1)
      SY = STRESS(2)
      TXY = STRESS(3)
      SZ = STRESS(4)
      DX = (2.D0*SX-SY-SZ)/3.D0
      DY = (2.D0*SY-SZ-SX)/3.D0
      DZ = (2.D0*SZ-SX-SY)/3.D0
      SIGM = (SX+SY+SZ)/3.D0
      DO 1 I = 1,4
          DO 1 J = 1,4
              M1(I,J) = 0.D0
              M2(I,J) = 0.D0
    1 M3(I,J) = 0.D0
      M1(1,1) = 1.D0
      M1(1,2) = 1.D0
      M1(2,1) = 1.D0
      M1(1,4) = 1.D0
      M1(4,1) = 1.D0
      M1(2,2) = 1.D0
      M1(2,4) = 1.D0
      M1(4,2) = 1.D0
      M1(4,4) = 1.D0
      DO 2 I = 1,4
          DO 2 J = 1,4
    2 M1(I,J) = M1(I,J)/ (9.D0*SIGM)
      M2(1,1) = .6666666666666666D0
      M2(2,2) = .6666666666666666D0
      M2(4,4) = .6666666666666666D0
      M2(3,3) = 2.D0
      M2(2,4) = -.3333333333333333D0
      M2(4,2) = -.3333333333333333D0
      M2(1,2) = -.3333333333333333D0
      M2(2,1) = -.3333333333333333D0
      M2(1,4) = -.3333333333333333D0
      M2(4,1) = -.3333333333333333D0
      M3(1,1) = DX/3.D0
      M3(2,4) = DX/3.D0
      M3(4,2) = DX/3.D0
      M3(2,2) = DY/3.D0
      M3(1,4) = DY/3.D0
      M3(4,1) = DY/3.D0
      M3(4,4) = DZ/3.D0
      M3(1,2) = DZ/3.D0
      M3(2,1) = DZ/3.D0
      M3(3,3) = -DZ
      M3(3,4) = -2.D0*TXY/3.D0
      M3(4,3) = -2.D0*TXY/3.D0
      M3(1,3) = TXY/3.D0
      M3(3,1) = TXY/3.D0
      M3(2,3) = TXY/3.D0
      M3(3,2) = TXY/3.D0
      RETURN
 
      END
      SUBROUTINE FORMM3(STRESS,M1,M2,M3)
C
C      THIS SUBROUTINE FORMS THE DERIVATIVES OF THE INVARIANTS
C      WITH RESPECT TO THE STRESSES (3-D)
C
      DOUBLE PRECISION SX
      DOUBLE PRECISION SY
      DOUBLE PRECISION SZ
      DOUBLE PRECISION TXY
      DOUBLE PRECISION TYZ
      DOUBLE PRECISION TZX
      DOUBLE PRECISION SIGM
      DOUBLE PRECISION DX
      DOUBLE PRECISION DY
      DOUBLE PRECISION DZ
      DOUBLE PRECISION STRESS(*),M1(6,6),M2(6,6),M3(6,6)
 
      SX = STRESS(1)
      SY = STRESS(2)
      SZ = STRESS(3)
      TXY = STRESS(4)
      TYZ = STRESS(5)
      TZX = STRESS(6)
      SIGM = (SX+SY+SZ)/3.D0
      DX = SX - SIGM
      DY = SY - SIGM
      DZ = SZ - SIGM
      DO 1 I = 1,6
          DO 1 J = I,6
              M1(I,J) = 0.D0
    1 M2(I,J) = 0.D0
      DO 2 I = 1,3
          DO 2 J = 1,3
    2 M1(I,J) = 1.D0/ (3.D0*SIGM)
      DO 3 I = 1,3
          M2(I,I) = 2.D0
    3 M2(I+3,I+3) = 6.D0
      M2(1,2) = -1.D0
      M2(1,3) = -1.D0
      M2(2,3) = -1.D0
      M3(1,1) = DX
      M3(1,2) = DZ
      M3(1,3) = DY
      M3(1,4) = TXY
      M3(1,5) = -2.D0*TYZ
      M3(1,6) = TZX
      M3(2,2) = DY
      M3(2,3) = DX
      M3(2,4) = TXY
      M3(2,5) = TYZ
      M3(2,6) = -2.D0*TZX
      M3(3,3) = DZ
      M3(3,4) = -2.D0*TXY
      M3(3,5) = TYZ
      M3(3,6) = TZX
      M3(4,4) = -3.D0*DZ
      M3(4,5) = 3.D0*TZX
      M3(4,6) = 3.D0*TYZ
      M3(5,5) = -3.D0*DX
      M3(5,6) = 3.D0*TXY
      M3(6,6) = -3.D0*DY
      DO 4 I = 1,6
          DO 4 J = I,6
              M1(I,J) = M1(I,J)/3.D0
              M1(J,I) = M1(I,J)
              M2(I,J) = M2(I,J)/3.D0
              M2(J,I) = M2(I,J)
              M3(I,J) = M3(I,J)/3.D0
    4 M3(J,I) = M3(I,J)
      RETURN
 
      END
      SUBROUTINE FORMTB(KB,IKB,KM,IKM,G,IW,IDOF)
C
C      THIS SUBROUTINE ASSEMBLES THE ELEMENT MATRICES INTO AN
C      UNSYMMETRICAL BANDED MATRIX 'KB'.
C
      DOUBLE PRECISION KB(IKB,*),KM(IKM,*)
      INTEGER G(*)
 
      DO 1 I = 1,IDOF
          IF (G(I).EQ.0) GO TO 1
          DO 2 J = 1,IDOF
              IF (G(J).EQ.0) GO TO 2
              ICD = G(J) - G(I) + IW + 1
              KB(G(I),ICD) = KB(G(I),ICD) + KM(I,J)
    2     CONTINUE
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE FORMXI(FSOIL,FMAX,RF,RM,R0,XI)
C
C      THIS SUBROUTINE FORMS PART OF THE SPRING STIFFNESS TERM
C      FOR PROGRAM 12.2
C
      DOUBLE PRECISION FSOIL
      DOUBLE PRECISION FMAX
      DOUBLE PRECISION RF
      DOUBLE PRECISION RM
      DOUBLE PRECISION R0
      DOUBLE PRECISION XI
      DOUBLE PRECISION PHI
 
      PHI = FSOIL*R0*RF/FMAX
      XI = LOG((RM-PHI)/ (R0-PHI)) + PHI* (RM-R0)/ ((RM-PHI)* (R0-PHI))
      RETURN
 
      END
      SUBROUTINE FRMUPV(KE,IKE,C11,IC11,C12,IC12,C21,IC21,C23,IC23,C32,
     +                  IC32,NOD,NODF,ITOT)
C
C      THIS SUBROUTINE FORMS THE UNSYMMETRICAL STIFFNESS MATRIX
C      FOR THE U-V-P VERSION OF THE NAVIER STOKES EQUATIONS
C
      DOUBLE PRECISION KE(IKE,*),C11(IC11,*),C21(IC21,*),C23(IC23,*),
     +                 C32(IC32,*),C12(IC12,*)
 
      K = NOD + NODF
      DO 1 I = 1,NOD
          DO 1 J = 1,NOD
    1 KE(I,J) = C11(I,J)
      DO 2 I = 1,NOD
          DO 2 J = NOD + 1,K
    2 KE(I,J) = C12(I,J-NOD)
      DO 3 I = NOD + 1,K
          DO 3 J = 1,NOD
    3 KE(I,J) = C21(I-NOD,J)
      DO 4 I = NOD + 1,K
          DO 4 J = K + 1,ITOT
    4 KE(I,J) = C23(I-NOD,J-K)
      DO 5 I = K + 1,ITOT
          DO 5 J = NOD + 1,K
    5 KE(I,J) = C32(I-K,J-NOD)
      DO 6 I = K + 1,ITOT
          DO 6 J = K + 1,ITOT
    6 KE(I,J) = C11(I-K,J-K)
      RETURN
 
      END
      SUBROUTINE FSPARV(BK,KM,IKM,G,KDIAG,IDOF)
C
C      THIS SUBROUTINE ASSEMBLES THE ELEMENT STIFFNESS MATRIX INTO
C      THE GLOBAL MATRIX STORED AS A VECTOR ACCOUNTING FOR A
C      VARIABLE BANDWIDTH
C
      INTEGER KDIAG(*),G(*)
      DOUBLE PRECISION BK(*),KM(IKM,*)
 
      DO 1 I = 1,IDOF
          K = G(I)
          IF (K.EQ.0) GO TO 1
          DO 2 J = 1,IDOF
              IF (G(J).EQ.0) GO TO 2
              IW = K - G(J)
              IF (IW.LT.0) GO TO 2
              IVAL = KDIAG(K) - IW
              BK(IVAL) = BK(IVAL) + KM(I,J)
    2     CONTINUE
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE GAUSBA(PB,IPB,WORK,IWORK,N,IW)
C
C      THIS SUBROUTINE PERFORMS GAUSSIAN REDUCTION OF AN
C      UNSYMMETRIC BANDED MATRIX 'PB' .  ARRAY 'WORK'
C      USED AS WORKING SPACE.
C
      DOUBLE PRECISION S
      DOUBLE PRECISION PB(IPB,*),WORK(IWORK,*)
 
      IWP1 = IW + 1
      IQ = 2*IWP1 - 1
      IQP = IWP1
      IWP11 = IWP1 - 1
      DO 1 I = 1,IWP11
          DO 1 J = 1,IQ
              IF (J.GE.IWP1+I) GO TO 2
              PB(I,J) = PB(I,J+IWP1-I)
              GO TO 1
 
    2         PB(I,J) = 0.D0
              PB(N-I+1,J) = 0.D0
    1 CONTINUE
      DO 3 K = 1,N
          L = K + IWP1 - 1
          IF (L.GT.N) L = N
          IP = 0
          S = 1.D-10
          DO 4 I = K,L
              IF (ABS(PB(I,1)).LE.S) GO TO 4
              S = ABS(PB(I,1))
              IP = I
    4     CONTINUE
          IF (IP.EQ.0) GO TO 5
          IF (K.EQ.N) GO TO 11
          WORK(IWP1,K) = IP
          IQP = IQP - 1
          J = IWP1 + IP - K
          IF (IQP.LT.J) IQP = J
          IF (J.EQ.IWP1) GO TO 6
          DO 7 J = 1,IQP
              S = PB(K,J)
              PB(K,J) = PB(IP,J)
              PB(IP,J) = S
    7     CONTINUE
    6     K1 = K + 1
          DO 8 I = K1,L
              S = PB(I,1)/PB(K,1)
              DO 9 J = 2,IQ
                  IF (J.GT.IQP) GO TO 10
                  PB(I,J-1) = PB(I,J) - S*PB(K,J)
                  GO TO 9
 
   10             PB(I,J-1) = PB(I,J)
    9         CONTINUE
              PB(I,IQ) = 0.D0
              WORK(I-K,K) = S
    8     CONTINUE
    3 CONTINUE
    5 WRITE (6,FMT=*) ' SINGULAR'
   11 RETURN
 
      END
      SUBROUTINE GAUSS(SAMP,ISAMP,NGP)
C
C      THIS SUBROUTINE PROVIDES THE WEIGHTS AND SAMPLING POINTS
C      FOR GAUSS-LEGENDRE QUADRATURE
C
      DOUBLE PRECISION SAMP(ISAMP,*)
 
      GO TO (1,2,3,4,5,6,7) NGP
 
    1 SAMP(1,1) = 0.D0
      SAMP(1,2) = 2.D0
      GO TO 100
 
    2 SAMP(1,1) = 1.D0/SQRT(3.D0)
      SAMP(2,1) = -SAMP(1,1)
      SAMP(1,2) = 1.D0
      SAMP(2,2) = 1.D0
      GO TO 100
 
    3 SAMP(1,1) = .2D0*SQRT(15.D0)
      SAMP(2,1) = .0D0
      SAMP(3,1) = -SAMP(1,1)
      SAMP(1,2) = 5.D0/9.D0
      SAMP(2,2) = 8.D0/9.D0
      SAMP(3,2) = SAMP(1,2)
      GO TO 100
 
    4 SAMP(1,1) = .861136311594053D0
      SAMP(2,1) = .339981043584856D0
      SAMP(3,1) = -SAMP(2,1)
      SAMP(4,1) = -SAMP(1,1)
      SAMP(1,2) = .347854845137454D0
      SAMP(2,2) = .652145154862546D0
      SAMP(3,2) = SAMP(2,2)
      SAMP(4,2) = SAMP(1,2)
      GO TO 100
 
    5 SAMP(1,1) = .906179845938664D0
      SAMP(2,1) = .538469310105683D0
      SAMP(3,1) = .0D0
      SAMP(4,1) = -SAMP(2,1)
      SAMP(5,1) = -SAMP(1,1)
      SAMP(1,2) = .236926885056189D0
      SAMP(2,2) = .478628670499366D0
      SAMP(3,2) = .568888888888889D0
      SAMP(4,2) = SAMP(2,2)
      SAMP(5,2) = SAMP(1,2)
      GO TO 100
 
    6 SAMP(1,1) = .932469514203152D0
      SAMP(2,1) = .661209386466265D0
      SAMP(3,1) = .238619186083197D0
      SAMP(4,1) = -SAMP(3,1)
      SAMP(5,1) = -SAMP(2,1)
      SAMP(6,1) = -SAMP(1,1)
      SAMP(1,2) = .171324492379170D0
      SAMP(2,2) = .360761573048139D0
      SAMP(3,2) = .467913934572691D0
      SAMP(4,2) = SAMP(3,2)
      SAMP(5,2) = SAMP(2,2)
      SAMP(6,2) = SAMP(1,2)
      GO TO 100
 
    7 SAMP(1,1) = .949107912342759D0
      SAMP(2,1) = .741531185599394D0
      SAMP(3,1) = .405845151377397D0
      SAMP(4,1) = .0D0
      SAMP(5,1) = -SAMP(3,1)
      SAMP(6,1) = -SAMP(2,1)
      SAMP(7,1) = -SAMP(1,1)
      SAMP(1,2) = .129484966168870D0
      SAMP(2,2) = .279705391489277D0
      SAMP(3,2) = .381830050505119D0
      SAMP(4,2) = .417959183673469D0
      SAMP(5,2) = SAMP(3,2)
      SAMP(6,2) = SAMP(2,2)
      SAMP(7,2) = SAMP(1,2)
  100 CONTINUE
      RETURN
 
      END
      SUBROUTINE GCOORD(FUN,COORD,ICOORD,NOD,IT,GC)
C
C      THIS SUBROUTINE OBTAINS THE CARTESIAN COORDINATES OF THE
C      GAUSS-POINTS FROM THE SHAPE FUNCTIONS
C
      DOUBLE PRECISION GC(*),FUN(*),COORD(ICOORD,*)
 
      DO 1 I = 1,IT
          GC(I) = 0.D0
          DO 1 J = 1,NOD
    1 GC(I) = GC(I) + COORD(J,I)*FUN(J)
      RETURN
 
      END
      SUBROUTINE GEOM(IP,NN,NPILE,G)
C
C      THIS SUBROUTINE FORMS THE NODE NUMBERS OF PILE IP
C
      INTEGER G(*)
 
      DO 1 I = 1,NN
    1 G(I) = IP + (I-1)*NPILE
      RETURN
 
      END
      SUBROUTINE GSTRNG(IP,NODOF,NF,INF,G)
C
C      THIS SUBROUTINE SELECTS THE G-VECTOR FROM THE NF-DATA
C
      INTEGER NF(INF,*),G(*),NODE(2)
 
      NODE(1) = IP
      NODE(2) = IP + 1
      L = 0
      DO 1 I = 1,2
          DO 1 J = 1,NODOF
              L = L + 1
    1 G(L) = NF(NODE(I),J)
      RETURN
 
      END
      SUBROUTINE HING2(IP,HOLDR,COORD,IPROP,REACT,ACTION,BMP)
C
C      THIS SUBROUTINE FORMS THE END FORCES AND MOMENTS TO BE
C      APPLIED TO A MEMBER IF A JOINT HAS GONE PLASTIC
C
      DOUBLE PRECISION BMP
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION ELL
      DOUBLE PRECISION CSCH
      DOUBLE PRECISION SNCH
      DOUBLE PRECISION BM1
      DOUBLE PRECISION BM2
      DOUBLE PRECISION S1
      DOUBLE PRECISION S2
      DOUBLE PRECISION HOLDR(IPROP,*),COORD(IPROP,*),REACT(*),ACTION(*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      X2 = COORD(IP,3)
      Y2 = COORD(IP,4)
      ELL = SQRT((Y2-Y1)**2+ (X2-X1)**2)
      CSCH = (X2-X1)/ELL
      SNCH = (Y2-Y1)/ELL
      BM1 = 0.D0
      BM2 = 0.D0
      S1 = HOLDR(IP,3) + REACT(3)
      S2 = HOLDR(IP,6) + REACT(6)
      IF (ABS(S1).LE.BMP) GO TO 2
      IF (S1.GT.0.D0) BM1 = BMP - S1
      IF (S1.LE.0.D0) BM1 = -BMP - S1
    2 CONTINUE
      IF (ABS(S2).LE.BMP) GO TO 1
      IF (S2.GT.0.D0) BM2 = BMP - S2
      IF (S2.LE.0.D0) BM2 = -BMP - S2
    1 CONTINUE
      ACTION(1) = - (BM1+BM2)*SNCH/ELL
      ACTION(2) = (BM1+BM2)*CSCH/ELL
      ACTION(3) = BM1
      ACTION(4) = -ACTION(1)
      ACTION(5) = -ACTION(2)
      ACTION(6) = BM2
      RETURN
 
      END
      SUBROUTINE INVAR(STRESS,SIGM,DSBAR,THETA)
C
C      THIS SUBROUTINE FORMS THE STRESS INVARIANTS (2-D)
C
      DOUBLE PRECISION SIGM
      DOUBLE PRECISION DSBAR
      DOUBLE PRECISION THETA
      DOUBLE PRECISION SX
      DOUBLE PRECISION SY
      DOUBLE PRECISION TXY
      DOUBLE PRECISION SZ
      DOUBLE PRECISION DX
      DOUBLE PRECISION DY
      DOUBLE PRECISION DZ
      DOUBLE PRECISION XJ3
      DOUBLE PRECISION SINE
      DOUBLE PRECISION STRESS(*)
 
      SX = STRESS(1)
      SY = STRESS(2)
      TXY = STRESS(3)
      SZ = STRESS(4)
      SIGM = (SX+SY+SZ)/3.D0
      DSBAR = SQRT((SX-SY)**2+ (SY-SZ)**2+ (SZ-SX)**2+6.D0*TXY**2)/
     +        SQRT(2.D0)
      IF (DSBAR.EQ.0.D0) THEN
          THETA = 0.D0
 
      ELSE
          DX = (2.D0*SX-SY-SZ)/3.D0
          DY = (2.D0*SY-SZ-SX)/3.D0
          DZ = (2.D0*SZ-SX-SY)/3.D0
          XJ3 = DX*DY*DZ - DZ*TXY**2
          SINE = -13.5D0*XJ3/DSBAR**3
          IF (SINE.GT.1.D0) SINE = 1.D0
          IF (SINE.LT.-1.D0) SINE = -1.D0
          THETA = ASIN(SINE)/3.D0
      END IF
 
      RETURN
 
      END
      SUBROUTINE INVAR3(STRESS,SIGM,DSBAR,THETA)
C
C      THIS SUBROUTINE FORMS THE STRESS INVARIANTS FROM THE
C      STRESS COMPONENTS IN 3-D
C
      DOUBLE PRECISION SIGM
      DOUBLE PRECISION DSBAR
      DOUBLE PRECISION THETA
      DOUBLE PRECISION SQ3
      DOUBLE PRECISION S1
      DOUBLE PRECISION S2
      DOUBLE PRECISION S3
      DOUBLE PRECISION S4
      DOUBLE PRECISION S5
      DOUBLE PRECISION S6
      DOUBLE PRECISION D2
      DOUBLE PRECISION DS1
      DOUBLE PRECISION DS2
      DOUBLE PRECISION DS3
      DOUBLE PRECISION D3
      DOUBLE PRECISION SINE
      DOUBLE PRECISION STRESS(*)
 
      SQ3 = SQRT(3.D0)
      S1 = STRESS(1)
      S2 = STRESS(2)
      S3 = STRESS(3)
      S4 = STRESS(4)
      S5 = STRESS(5)
      S6 = STRESS(6)
      SIGM = (S1+S2+S3)/3.D0
      D2 = ((S1-S2)**2+ (S2-S3)**2+ (S3-S1)**2)/6.D0 + S4*S4 + S5*S5 +
     +     S6*S6
      DS1 = S1 - SIGM
      DS2 = S2 - SIGM
      DS3 = S3 - SIGM
      D3 = DS1*DS2*DS3 - DS1*S5*S5 - DS2*S6*S6 - DS3*S4*S4 +
     +     2.D0*S4*S5*S6
      DSBAR = SQ3*SQRT(D2)
      IF (DSBAR.EQ.0.D0) THEN
          THETA = 0.D0
 
      ELSE
          SINE = -3.D0*SQ3*D3/ (2.D0*SQRT(D2)**3)
          IF (SINE.GT.1.D0) SINE = 1.D0
          IF (SINE.LT.-1.D0) SINE = -1.D0
          THETA = ASIN(SINE)/3.D0
      END IF
 
      RETURN
 
      END
      SUBROUTINE KVDET(KV,N,IW,DET,KSC)
C
C      THIS SUBROUTINE FORMS THE DETERMINANT OF THE GLOBAL STIFFNESS
C      MATRIX STORED AS A VECTOR (UPPER TRIANGLE)
C
      DOUBLE PRECISION DET
      DOUBLE PRECISION CONST
      DOUBLE PRECISION EM
      DOUBLE PRECISION KV(*)
 
      IWP2 = IW + 2
      L = N - 1
      DO 900 J = 1,L
          IF (IW-N+J) 810,810,820
  810     MA = IWP2
  820     MA = MA - 1
          CONST = 1.D0/KV(J)
          II = J
          DO 880 K = 2,MA
              II = II + 1
              IF (KV((K-1)*N+J)) 840,880,840
  840         EM = KV((K-1)*N+J)*CONST
              LL = 0
              DO 860 I = K,MA
                  LL = LL + 1
                  KV((LL-1)*N+II) = KV((LL-1)*N+II) - EM*KV((I-1)*N+J)
  860         CONTINUE
  880     CONTINUE
  900 CONTINUE
      DET = 1.D0
      KSC = 0
      DO 950 J = 1,N
          DET = DET*KV(J)
          IF (KV(J).LT.0.D0) KSC = KSC + 1
  950 CONTINUE
      RETURN
 
      END
      SUBROUTINE LANCZ1(N,EL,ER,ACC,LEIG,LX,LALFA,LP,ITAPE,IFLAG,U,V,
     +                  EIG,JEIG,NEIG,X,DEL,NU,ALFA,BETA)
C USE THE LANCZOS ALGORITHM TO FIND THE SPECTRUM OF A SYMMETRIC MATRIX
C     OR THAT PART OF THE SPECTRUM THAT LIES IN A GIVEN INTERVAL. THIS
C     VERSION RETURNS INFORMATION FOR LATER CALCULATION OF EIGENVECTORS
C     BY LANCZ2.
C  STANDARD FORTRAN 66 (A VERIFIED PFORT SUBROUTINE)
      DOUBLE PRECISION ANORM
      DOUBLE PRECISION G
      DOUBLE PRECISION FA01AS
      DOUBLE PRECISION EL,ER,ACC,U(N),V(N),EIG(LEIG),X(LX),DEL(LX),
     +                 ALFA(LALFA),BETA(LALFA)
      INTEGER NU(LX),JEIG(2,LEIG)
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C EL,ER MUST BE SET TO INDICATE THE RANGE WITHIN WHICH THE SPECTRUM
C     IS WANTED. IF EL.GE.ER THEN IT IS ASSUMED THAT THE WHOLE SPECTRUM
C     IS WANTED.THEY ARE NOT ALTERED.
C ACC MUST BE SET TO THE REQUIRED PRECISION, RELATIVE TO THE LARGEST
C     EIGENVALUE OF A. IF IT IS VERY SMALL OR NEGATIVE THEN AS MUCH
C     ACCURACY AS THE PRECISION ALLOWS WILL BE FOUND. IT IS NOT ALTERED.
C LEIG MUST BE SET TO LENGTH OF ARRAY EIG. IT MUST BE AS LARGE AS THE
C     NUMBER OF DISTINCT EIGENVALUES IN THE RANGE (EL,ER). IT IS NOT
C     ALTERED.
C LX MUST BE SET TO THE LENGTH OF ARRAYS X,DEL,NU. A VALUE THREE TIMES
C     NUMBER OF DISTINCT EIGENVALUES IN THE RANGE (EL,ER) USUALLY
C     SUFFICES. IT IS NOT ALTERED.
C LALFA MUST BE SET TO THE LENGTH OF ARRAYS ALFA AND BETA. IT LIMITS
C      THE NUMBER OF LANCZOS STEPS POSSIBLE. IT IS NOT ALTERED.
C LP MUST BE SET TO THE UNIT NUMBER FOR DIAGNOSTIC MESSAGES. IF LP.LE.0
C     THE MESSAGES ARE SUPPRESSED. IT IS NOT ALTERED.
C ITAPE MUST BE SET TO THE UNIT NUMBER OF A SEQUENTIAL FILE IF
C     EIGENVECTORS ARE WANTED FROM LANCZ2 OR NON-POSITIVE IF THEY ARE
C     NOT WANTED. IT IS NOT ALTERED.
C IFLAG MUST BE SET PRIOR TO THE FIRST CALL FOR A PARTICULAR MATRIX TO
C       -1 IF THE USER DOES NOT WANT TO SPECIFY A START VECTOR, OR
C       -2 IF THE USER HAS PLACED A REQUIRED START VECTOR IN V.
C     IT SHOULD NOT OTHERWISE BE CHANGED.
C     ON RETURN IT HAS THE VALUE
C       0 ON SUCCESSFUL COMPLETION
C       1 ON A INTERMEDIATE RETURN
C       2,3,...,7 ON ERROR CONDITIONS
C U,V HOLD THE WORKING VECTORS OF THE LANCZOS RECURRENCE. ON A RETURN
C     WITH IFLAG=1 THE USED MUST ADD TO U THE VECTOR A*V WITHOUT
C     ALTERING V.
C EIG NEED NOT BE SET BY THE USER. ON ANY RETURN IT CONTAINS NEIG
C     COMPUTED EIGENVALUES, STORED IN INCREASING ORDER.
C JEIG NEED NOT BE SET BY THE USER. ON RETURN JEIG(1,I) CONTAINS THE
C     LANCZOS STEP AT WHICH EIG(I) WAS ACCEPTED AND JEIG(2,I) CONTAINS
C     THE MATCHING POINT FOR THE RECURRENCES TO GET THE EIGENVECTOR.
C NEIG NEED NOT BE SET BY THE USER. ON ANY RETURN IT CONTAINS THE NUMBER
C     OF EIGENVALUES IN EIG.
C X IS USED FOR A WORK ARRAY.  X(1).LT.X(2).LT.X(3).LT. ... ARE POINTS
C     APPROXIMATING EIGENVALUES OF A.
C DEL IS USED AS A WORKARRAY.  DEL(I),I=1,2,3,... CONTAIN THE LAST PIVOT
C     IN THE LU FACTORIZATION OF T-LAMDA*I, WHERE T IS THE LANCZOS
C     TRIDIAGONAL MATRIX, I IS THE IDENTITY MATRIX AND LAMDA=X(I).
C NU IS USED AS A WORKARRAY.  IABS(NU(I))-1,I=1,2,3,... CONTAIN THE
C     NUMBER OF NEGATIVE PIVOTS IN THE LU FACTORIZATION OF T-LAMDA*I,
C     LAMDA=X(I). FIXED INTERVALS X(I),X(I+1)  ABOUT CONVERGED
C     EIGENVALUES ARE FLAGGED BY NEGATIVE VALUES OF NU(I).
C ALFA AND BETA ARE USED AS WORKARRAYS HOLDING THE LANCZOS TRIDIAGONAL
C     MATRIX T. BETA(1)=0.  AND A TYPICAL ROW IS
C         BETA(I)     ALFA(I)     BETA(I+1)
C
      DOUBLE PRECISION DRELPR,VV,UU,ALF,BET,VI,TOL,BE2,OLDEL,OLDER,DR,
     +                 XMID,XA,XC,FA,FC,TEM,DEN,PI,SIG,TAU,DISC,XAV,
     +                 POLE,RITZ,XL,XR,DIFF,ERR,BND,T,TOLC,W,XLI,XRI,EN
      COMMON /EA15BD/ANORM,OLDEL,OLDER,TOLC,KX,MAXRZ,MAXRZO,NLAN,MLAN,
     +       NXTBND
C COMMON BLOCK EA15BD CONTAINS THE FOLLOWING QUANTITIES,
C     WHICH ARE REQUIRED TO BE PRESERVED BETWEEN ENTRIES:
C   ANORM IS A ESTIMATE OF THE NORM OF THE MATRIX A, BASED
C      ON A GERSHGORIN BOUND APPLIED TO T.
C   OLDEL,OLDER ARE THE VALUES OF EL,ER ON THE PREVIOUS ENTRY.
C   TOLC IS THE TOLERANCE FOR AGREEMENT BETWEEN SUCESSIVE RITZ
C      VALUES THAT DECIDES WHETHER TO CALL THE ERROR BOUNDING ROUTINE
C      EA15CD (BUT SEE ALSO NXTBND)
C   KX IS THE NUMBER OF POINTS CURRENTLY STORED IN X, WITH ASSOCIATED
C      INFORMATION IN DEL AND NU.
C   MAXRZ IS THE MAXIMAL NUMBER OF RITZ VALUES IN A GAP BETWEEN
C      FIXED INTERVALS.
C   MAXRZO IS THE PREVIOUS VALUE OF MAXRZ
C   NLAN IS THE NUMBER OF LANCZOS STEPS PERFORMED. ALFA(I)<
C      BETA(I),I=1,2,...,NLAN HAVE BEEN SET.
C   MLAN IS THE ORDER OF TRIDIAGONAL MATRIX SO FAR USED FOR THIS
C      INTERVAL OF THE SPECTRUM.
C   NXTBND IS USED TO DELAY ERROR BOUNDING WHEN TOLC IS AT ROUNDOFF
C      ERROR LEVEL. EA15CD IS NOT CALLED UNTIL WE ARE LOOKING
C      AT THE TRIDIAGONAL MATRIX OF ORDER NXTBND.
      DOUBLE PRECISION ZERO,ONE,TWO,HALF,THREE
C DRELPR IS THE RELATIVE PRECISION
      DATA DRELPR/2.2D-16/,ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/
CRAY  DATA DRELPR/1.4E-14/, ZERO/0.0E0/, ONE/1.0E0/, TWO/2.0E0/
CIBM  DATA DRELPR/2.2E-16/, ZERO/0.0E0/, ONE/1.0E0/, TWO/2.0E0/
      DATA HALF/0.5D0/,THREE/3.0D0/
 
      EN = DBLE(N)
      IF (N.LE.0 .OR. LX.LE.5 .OR. LALFA.LE.0) GO TO 950
      IF (IFLAG.EQ.1) GO TO 35
      VV = ZERO
      NXTBND = 0
      IF (IFLAG.EQ.-2) GO TO 5
      IF (IFLAG.EQ.-1) GO TO 15
      IF (IFLAG.EQ.0) GO TO 70
      GO TO 970
C
C FIND THE NORM OF THE USER~S START VECTOR.
    5 DO 10 I = 1,N
          VV = VV + V(I)**2
   10 CONTINUE
      IF (VV.GT.ZERO) GO TO 25
C
C GENERATE PSEUDO-RANDOM START VECTOR
   15 G = 1431655765.D0
      DO 20 I = 1,N
          G = DMOD(G*9228907.D0,4294967296.D0)
          IF (I.GE.0) FA01AS = G/4294967296.D0
          IF (I.LT.0) FA01AS = 2.D0*G/4294967296.D0 - 1.D0
          V(I) = FA01AS
          VV = VV + V(I)**2
   20 CONTINUE
C
C NORMALIZE START VECTOR AND PERFORM INITIALIZATIONS
   25 VV = ONE/SQRT(VV)
      DO 30 I = 1,N
          V(I) = V(I)*VV
          U(I) = ZERO
   30 CONTINUE
      NLAN = 0
      ANORM = ZERO
      BETA(1) = ZERO
      NEIG = 0
      GO TO 912
C
C PERFORM A LANCZOS STEP
   35 NLAN = NLAN + 1
      CALL EA15FD(N,ITAPE,NLAN,1,V)
      IF (NLAN.GE.LALFA) GO TO 940
      ALF = ZERO
      DO 40 I = 1,N
          ALF = ALF + U(I)*V(I)
   40 CONTINUE
      ALFA(NLAN) = ALF
      UU = ZERO
      DO 50 I = 1,N
          U(I) = U(I) - ALF*V(I)
          UU = UU + U(I)**2
   50 CONTINUE
      UU = SQRT(UU)
      BET = BETA(NLAN)
      ANORM = DMAX1(ANORM,BET+ABS(ALF)+UU)
      BET = DMAX1(UU,ANORM*DRELPR)
      UU = ONE/BET
      BETA(NLAN+1) = BET
      DO 60 I = 1,N
          VI = U(I)*UU
          U(I) = -BET*V(I)
          V(I) = VI
   60 CONTINUE
C
C THIS IS THE BEGINNING OF THE LOOP THAT ANALYSES T
C NORMALLY WE ADVANCE THE ANALYSIS OF T BY ONE ROW, BUT THE
C     USER MAY HAVE REQUESTED A RESTART
   70 JLAN1 = NLAN
      IF (OLDEL.EQ.EL .AND. OLDER.EQ.ER .AND. IFLAG.EQ.1) GO TO 71
      JLAN1 = 1
      NEIG = 0
   71 DO 910 JLAN = JLAN1,NLAN
          MLAN = JLAN
          IF (JLAN-2) 73,76,85
C
C JLAN=1. CHECK FOR THE TRIVIAL CASE.
   73     EIG(1) = ALFA(1)
          JEIG(1,1) = JLAN
          JEIG(2,1) = JLAN
          IF (BETA(2).GT.DRELPR*ANORM) GO TO 910
          IF (EL.GE.ER .OR. (EL.LE.EIG(1).AND.EIG(1).LE.ER)) NEIG = 1
          GO TO 960
C
C JLAN=2. SET UP FOUR POINTS THAT SPAN AND SEPARATE  THE
C     THREE RITZ VALUES AT LEVELS 1 AND 2
   76     W = HALF*ABS(ALFA(1)-ALFA(2))
          T = SQRT(W*W+BETA(2)**2)
          X(2) = HALF* (ALFA(1)+ALFA(2)-T-W)
          X(3) = X(2) + T + W
          X(1) = X(2) + W - 1.1D0*T
          X(4) = X(3) - W + 1.1D0*T
          DO 80 L = 1,4
              DEL(L) = ALFA(2) - X(L) - BETA(2)**2/ (ALFA(1)-X(L))
              NU(L) = L - (L-1)/2
   80     CONTINUE
          KX = 4
          MAXRZ = 2
          MAXRZO = 1
          IF (X(1).LT.X(2) .AND. X(3).LT.X(4)) GO TO 910
          EIG(1) = (ALFA(1)+ALFA(2))*HALF
          JEIG(1,1) = JLAN
          JEIG(2,1) = JLAN
          NEIG = 1
          GO TO 960
C
C JLAN.GT.2
C
C ADD OR REMOVE POINTS TO THE RIGHT, IF APPROPRIATE
   85     NUK = JLAN - 1
          IF (EL.GE.ER) GO TO 100
C FIND INTERVAL (X(K-1),X(K)) CONTAINING ER
          DO 92 I = 2,KX
              K = 2 + KX - I
              IF (X(K-1).LE.ER) GO TO 95
   92     CONTINUE
          K = 1
   95     NUK = MIN0(IABS(NU(K)),JLAN-1)
C LOOK FOR POINT WITH NU VALUE AT LEAST 2 GREATER THAN ANY IN INTERVAL
C     CONTAINING ER
          DO 97 I = K,KX
              IF (IABS(NU(I)).GT.NUK+1) GO TO 98
   97     CONTINUE
          GO TO 100
C  REDUCE KX TO GET RID OF UNNECESSARY POINTS
   98     KX = I
          GO TO 110
C IF NECESSARY ADD EXTRA POINTS TO RIGHT
  100     DO 105 IDUMMY = 1,LX
              IF (IABS(NU(KX-1)).GT.NUK) GO TO 110
              IF (KX.GE.LX) GO TO 930
              KX = KX + 1
              X(KX) = X(KX-1)*THREE - X(KX-2)*TWO
              CALL EA15DD(JLAN-1,ALFA,BETA,X(KX),DEL(KX),NU(KX),DR,NUR)
  105     CONTINUE
C
C COPY INFORMATION TO ENDS OF ARRAYS X,DEL,NU. IF ANY
C     GAP BETWEEN FIXED INTERVALS CONTAINS NO RITZ VALUES,
C     REMOVE ALL POINTS IN THAT GAP
  110     LFP = LX
C LFP IS THE LAST LEFT-END OF A FIXED INTERVAL
          M = LX
          DO 140 IDUMMY = 1,LX
              IF (NU(KX-1).GE.0) GO TO 130
              IF (NU(KX).EQ.IABS(NU(LFP))) M = LFP - 1
              LFP = M - 1
  130         X(M) = X(KX)
              DEL(M) = DEL(KX)
              NU(M) = NU(KX)
              KX = KX - 1
              M = M - 1
              IF (KX.EQ.1) GO TO 150
  140     CONTINUE
  150     X(M) = X(1)
          DEL(M) = DEL(1)
          NU(M) = NU(1)
C
C ADD OR REMOVE POINTS TO THE LEFT, IF APPROPRIATE
          NUK = 2
          IF (EL.GE.ER) GO TO 175
C FIND INTERVAL (X(K-1),X(K)) CONTAINING EL
          M1 = M + 1
          DO 155 K = M1,LX
              IF (X(K).GE.EL) GO TO 160
  155     CONTINUE
          K = LX + 1
  160     NUK = MAX0(IABS(NU(K-1)),2)
C LOOK FOR POINT WITH NU VALUE AT LEAST 2 LESS THAN ANY IN INTERVAL
C     CONTAINING EL
          DO 165 J = M1,K
              I = M + K - J
              IF (IABS(NU(I)).LT.NUK-1) GO TO 170
  165     CONTINUE
          GO TO 175
C INCREASE M TO GET RID OF UNNECESSARY POINTS
  170     M = I
          GO TO 190
C IF NECESSARY, ADD EXTRA POINTS TO LEFT
  175     DO 180 IDUMMY = 1,LX
              IF (IABS(NU(M+1)).LT.NUK) GO TO 190
              IF (M.EQ.3) GO TO 930
              M = M - 1
              X(M) = THREE*X(M+1) - TWO*X(M+2)
              CALL EA15DD(JLAN-1,ALFA,BETA,X(M),DEL(M),NU(M),DR,NUR)
  180     CONTINUE
  190     K = 3
          IF (M.LE.5) GO TO 930
          XRI = X(M)
          ALF = ALFA(JLAN)
          BET = BETA(JLAN)
          BE2 = BET**2
C TOL HOLDS THE RADIUS OF FIXED INTERVALS SET UP AROUND CONVERGED
C     EIGENVALUES.
          TOL = DMAX1(ACC,EN*TWO*DRELPR)*ANORM
C TOLC IS THE TOLERANCE USED FOR ACCEPTING EIGENVALUES.
          IF (JLAN.LT.10) TOLC = DMAX1(TOL**2*EN/ANORM,
     +                           ANORM*DRELPR*5.D0)
C MOVE LEADING THREE POINTS TO BEGINNING OF STORE AND UPDATE
C     ASSOCIATED DATA.
          M = M - 1
          DO 210 I = 1,3
              M = M + 1
              X(I) = X(M)
              NU(I) = NU(M)
              DEL(I) = ALF - X(M) - BE2/DEL(M)
              IF (DEL(I)) 195,200,210
  195         NU(I) = NU(I) + ISIGN(1,NU(I))
              GO TO 210
 
  200         DEL(I) = DRELPR*ANORM
  210     CONTINUE
C
C PROCESS THE POINTS FROM LEFT TO RIGHT. THE CURRENT INTERVAL
C     IS (X(K-1),X(K))=(X(M-1),X(M)), K.LT.M-2.
C NU(I),DEL(I),I=1,2,...,K HOLD NEW VALUES OF NU,DELTA.
C NU(I),DEL(I),I=M-2,M-1,...,LX HOLD OLD VALUES.
          NE = 1
C NE NORMALLY HOLDS THE NUMBER OF EMPTY INTERVALS ADJACENT ON THE
C     LEFT OF THE CURRENT INTERVAL.
          LFP = 1
C LFP POINTS TO THE LAST FIXED POINT ENCOUNTERED.
C      DO 530 UNTIL M.EQ.LX
  220     IF (M.GE.LX) GO TO 600
          IF (NU(K-1).GE.0) GO TO 240
C WE HAVE A FIXED INTERVAL
          IF (IABS(NU(K-1)).EQ.IABS(NU(K))) GO TO 230
C ACCEPT FIXED INTERVAL
          LFP = K
          GO TO 500
C FIXED INTERVAL NO LONGER CONTAINS RITZ VALUE. FREE IT.
  230     NU(K-1) = IABS(NU(K-1))
          DO 235 I = 1,NEIG
              IF (EIG(I).LE.X(K)) GO TO 235
              EIG(I-1) = EIG(I)
              JEIG(1,I-1) = JEIG(1,I)
              JEIG(2,I-1) = JEIG(2,I)
  235     CONTINUE
          NEIG = NEIG - 1
  240     IF (NU(K-1).LT.IABS(NU(K))) GO TO 250
C CURRENT INTERVAL CONTAINS NO RITZ VALUE.
          NE = NE + 1
          IF (NE.LE.3) GO TO 502
C THERE ARE FOUR ADJACENT EMPTY INTERVALS. REMOVE MIDDLE POINT.
          X(K-2) = X(K-1)
          NU(K-2) = NU(K-1)
          DEL(K-2) = DEL(K-1)
          X(K-1) = X(K)
          NU(K-1) = NU(K)
          DEL(K-1) = DEL(K)
          NE = 3
          K = K - 1
          GO TO 502
 
  250     IF (IABS(NU(M-1)).EQ.IABS(NU(M))) GO TO 500
C
C INTERVAL IS ~INTERESTING~. IT CONTAINS AT LEAST ONE NEW AND
C     AT LEAST ONE OLD RITZ VALUE.
C JUMP IF INTERVAL IS A GAP BETWEEN FIXED INTERVALS,
C     AND CONTAINS JUST ONE RITZ VALUE, UNLESS THIS IS THE
C     FIRST TIME THAT ALL SUCH GAPS HAVE LESS THAN TWO
C     RITZ VALUES.
          IF (MAXRZ.LE.1 .AND. MAXRZO.GT.1) GO TO 252
          IF (NU(K-1)+NU(K).EQ.-1) GO TO 500
  252     XMID = (X(M-1)+X(M))*HALF
C BISECT IF INTERVAL CONTAINS MORE THAN ONE RITZ VALUE.
          RITZ = XMID
          IF (IABS(NU(K)).GT.NU(K-1)+1) GO TO 275
C CHECK WHETHER PROGRESS IS SO SLOW THAT BISECTION IS NEEDED
C     THE CRITERION IS THAT THREE STEPS HAVE BEEN TAKEN WITHOUT
C     REDUCING THE INTERVAL LENGTH BY AT LEAST THE FACTOR 0.6.
C     (XLI,XRI) IS THE ORIGINAL INTERVAL OF INTEREST OR THE
C     INTERVAL OF INTEREST WHEN IT WAS LAST REDUCED IN LENGTH
C      BY AT LEAST THE FACTOR 0.6
          IF (X(K-1).GE.XRI) GO TO 255
          IF (X(K)-X(K-1).GE.0.6D0* (XRI-XLI)) GO TO 257
  255     NXTR = 0
          XRI = X(K-1)
          XLI = X(K)
C BISECT IF PROGRESS ON THIS INTERVAL IS SLOW
  257     IF (NXTR.GE.2) GO TO 275
          NXTR = NXTR + 1
C ESTIMATE POSITION (POLE) OF OLD RITZ VALUE BY 2-1 RATIONAL
C     INTERPOLATION AT X(L-1), X(L), X(L+1)
C CHOOSE L SO THAT EXTRA POINT IS NEAR
          L = M
          IF (XMID-X(M-2).LT.X(M+1)-XMID) L = M - 1
          XA = X(L-1) - X(L)
          XC = X(L+1) - X(L)
          FA = XA* (XA+DEL(L-1))
          FC = XC* (XC+DEL(L+1))
          TEM = XA/ (XA-XC)
          DEN = DEL(L) - DEL(L-1) - TEM* (DEL(L+1)-DEL(L-1))
          IF (DEN.EQ.ZERO) GO TO 275
          PI = (TEM* (FA-FC)-FA)/DEN
          SIG = HALF* ((FC-FA)-PI* (DEL(L+1)-DEL(L-1)))/ (XC-XA)
          TAU = FA - TWO*XA*SIG - DEL(L-1)*PI
          DISC = SIG**2 + TAU
          IF (DISC.LT.ZERO) GO TO 275
          XAV = SIG + SIGN(SQRT(DISC),SIG)
          POLE = -TAU/XAV
C JUMP IF POLE IS IN REQUIRED INTERVAL
          IF (ABS(POLE+X(L)-XMID).LE.X(M)-XMID) GO TO 260
C TRY OTHER ROOT.
          POLE = XAV
          IF (ABS(POLE+X(L)-XMID).GT.X(M)-XMID) GO TO 275
C ESTIMATE POSITION (RITZ) OF RITZ VALUE BY 2-1 RATIONAL
C     INTERPOLATION AT X(K-1),X(K) USING (X-POLE) FOR
C     DENOMINATOR.
  260     IF (L.NE.M) GO TO 265
          TAU = -POLE*DEL(K)
          SIG = XA + ((XA-POLE)*DEL(K-1)-TAU)/XA
          GO TO 270
 
  265     TAU = -POLE*DEL(K-1)
          SIG = XC + ((XC-POLE)*DEL(K)-TAU)/XC
  270     SIG = SIG*HALF
          DISC = SIG**2 + TAU
          IF (DISC.LT.ZERO) GO TO 275
          XAV = SIG + SIGN(SQRT(DISC),SIG)
          TAU = -TAU/XAV
          IF (ABS(TAU+X(L)-XMID).LE.X(M)-XMID) GO TO 280
          TAU = XAV
          IF (ABS(TAU+X(L)-XMID).LE.X(M)-XMID) GO TO 280
C CALCULATION HAS FAILED
C IF THIRD POINT IS JUST OUTSIDE CURRENT INTERVAL THEN TAKE
C     POINT TWICE AS FAR INSIDE INTERVAL. OTHERWISE BISECT.
  275     RITZ = XMID
          POLE = XMID
          IF (X(M)-X(M-1).LE.TOLC) GO TO 300
          GO TO 490
 
  280     RITZ = X(L) + TAU
          DIFF = ABS(TAU-POLE)
          POLE = X(L) + POLE
          IF (ABS(TAU).GE.50.D0*DMAX1(DIFF,TOLC)) GO TO 480
          IF (DIFF.GT.TOLC) GO TO 500
C WE MAY HAVE A CONVERGED EIGENVALUE
  300     IF (JLAN.LT.NXTBND) GO TO 500
          CALL EA15CD(JLAN,LALFA,ALFA,BETA,RITZ,ANORM*DRELPR*EN,ERR,BND,
     +                MATCH)
          IF (BND.LT.ERR*0.1D0) TOLC = DIFF* (TOL/ERR)**2
          IF (TOLC.GT.ANORM*DRELPR*5.0D0) GO TO 305
          TOLC = ANORM*DRELPR*5.0D0
C TOLC HAS BECOME TOO SMALL TO GIVE GOOD CRITERION FOR CALLING EA15CD.
C    DO NOT CALL EA15CD AGAIN UNTIL 1% MORE STEPS PERFORMED.
          NXTBND = JLAN + JLAN/100
  305     IF (ERR.LE.TOL) GO TO 310
          IF (BND.GT.ERR*0.1D0 .AND. X(M-1).LT.RITZ .AND.
     +        RITZ.LT.X(M)) GO TO 480
          GO TO 500
C
C WE HAVE AN ACCEPTED POINT
C TEST WHETHER (RITZ-ERR,RITZ+ERR) OVERLAPS A FIXED INTERVAL
  310     IF (RITZ-ERR.LT.X(LFP) .AND. LFP.GT.1) GO TO 500
          DO 410 I = M,LX
              IF (X(I).GT.RITZ+ERR) GO TO 420
              IF (NU(I).LT.0) GO TO 500
  410     CONTINUE
C SET UP NEW FIXED INTERVAL
  420     XL = RITZ - TOL
          IF (XL.LE.X(LFP)) GO TO 450
C REMOVE POINTS IN INTERVAL (XL,RITZ)
          DO 430 IDUMMY = 1,LX
              IF (X(K-1).LT.XL) GO TO 440
              K = K - 1
  430     CONTINUE
  440     X(K) = XL
          CALL EA15DD(JLAN,ALFA,BETA,XL,DEL(K),NU(K),DR,NUR)
          GO TO 455
 
  450     K = LFP
  455     NU(K) = -NU(K)
          XR = RITZ + TOL
C REMOVE POINTS IN INTERVAL (RITZ,XR)
          I = M
          DO 460 M = I,LX
              IF (X(M).GT.XR) GO TO 470
              IF (NU(M).LT.0) GO TO 473
  460     CONTINUE
          M = LX + 1
  470     M = M - 1
          X(M) = XR
          CALL EA15DD(JLAN-1,ALFA,BETA,X(M),DEL(M),NU(M),DR,NUR)
  473     NFIX = NFIX + 1
          IF (NEIG.GE.LEIG) GO TO 920
          IF (EL.GE.ER) GO TO 474
          IF (RITZ.LT.EL .OR. RITZ.GT.ER) GO TO 478
  474     I = NEIG
          IF (NEIG.LE.0) GO TO 477
          DO 475 J = 1,NEIG
              IF (EIG(I).LT.RITZ) GO TO 477
              EIG(I+1) = EIG(I)
              JEIG(1,I+1) = JEIG(1,I)
              JEIG(2,I+1) = JEIG(2,I)
              I = I - 1
  475     CONTINUE
  477     EIG(I+1) = RITZ
          JEIG(1,I+1) = JLAN
          JEIG(2,I+1) = MATCH
          NEIG = NEIG + 1
  478     IF (M.GT.K+2) GO TO 505
          GO TO 930
C EXTRAPOLATE TO ESTIMATE POSITION OF EIGENVALUE
  480     RITZ = RITZ + RITZ - POLE
          RITZ = DMAX1(RITZ,X(M-1)+ (X(M)-X(M-1))*0.01D0)
          RITZ = DMIN1(RITZ,X(M)- (X(M)-X(M-1))*0.01D0)
C INSERT NEW POINT
  490     IF (K.GE.M-3) GO TO 930
          M = M - 1
          X(M-2) = X(M-1)
          DEL(M-2) = DEL(M-1)
          NU(M-2) = NU(M-1)
          X(M-1) = X(M)
          DEL(M-1) = DEL(M)
          NU(M-1) = NU(M)
          X(M) = RITZ
          X(K) = RITZ
          CALL EA15DD(JLAN,ALFA,BETA,RITZ,DEL(K),NU(K),DR,NUR)
          DEL(M) = DR
          NU(M) = NUR
          GO TO 530
C INTERVAL IS ACCEPTABLE. ADVANCE BY ONE POINT
  500     NE = 0
  502     M = M + 1
  505     K = K + 1
          X(K) = X(M)
          DEL(K) = ALF - X(K) - BE2/DEL(M)
          NU(K) = NU(M)
          IF (DEL(K)) 510,520,530
  510     NU(K) = NU(K) + ISIGN(1,NU(K))
          GO TO 530
 
  520     DEL(K) = DRELPR*ANORM
  530     GO TO 220
C
C SCAN TO FIND THE MAXIMUM NUMBER OF RITZ VALUES (MAXRZ) IN A GAP
C     BETWEEN FIXED INTERVALS, THE NUMBER OF GAPS CONTAINING CANDIDATES
C     THAT ARE BEING WATCHED (NCAND) AND THE NUMBER OF FIXED INTERVALS
C     (NFIX).
  600     KX = K
          NCAND = 0
          MAXRZO = MAXRZ
          MAXRZ = 0
          NFIX = 0
          LFP = 2
          K = K - 1
          DO 710 I = 1,K
              IF (NU(I).GT.0) GO TO 710
              NFIX = NFIX + 1
              NRITZ = IABS(NU(I)) - IABS(NU(LFP))
              MAXRZ = MAX0(MAXRZ,NRITZ)
              IF (NRITZ.GT.0 .AND. I.GT.LFP+1) NCAND = NCAND + 1
              DO 680 LFP = I,KX
                  IF (NU(LFP).GT.0) GO TO 710
  680         CONTINUE
              LFP = KX
  710     CONTINUE
          NRITZ = IABS(NU(K)) - IABS(NU(LFP))
          MAXRZ = MAX0(MAXRZ,NRITZ)
          IF (NRITZ.GT.0 .AND. K.GT.LFP+1) NCAND = NCAND + 1
          IF (BET.LE.EN*ANORM*DRELPR) GO TO 960
          IF (MAXRZ.GT.1 .OR. NCAND.GT.0) GO TO 910
          IF (NFIX.GT.0) GO TO 915
  910 CONTINUE
C
C NORMAL RETURNS
  912 IFLAG = 1
      GO TO 1000
 
  915 IFLAG = 0
      GO TO 1000
C
C ERROR RETURNS
  920 IFLAG = 2
      IF (LP.GT.0) WRITE (LP,FMT=925) LEIG
 
  925 FORMAT (34H ERROR RETURN 2 FROM LANCZ1. LEIG=,I6)
 
      GO TO 1000
 
  930 IFLAG = 3
      IF (LP.GT.0) WRITE (LP,FMT=935) LX
 
  935 FORMAT (32H ERROR RETURN 3 FROM LANCZ1. LX=,I6)
 
      GO TO 1000
 
  940 IFLAG = 4
      IF (LP.GT.0) WRITE (LP,FMT=945) LALFA
 
  945 FORMAT (35H ERROR RETURN 4 FROM LANCZ1. LALFA=,I6)
 
      GO TO 1000
 
  950 IFLAG = 5
      IF (LP.GT.0) WRITE (LP,FMT=955) N,LX,LALFA
 
  955 FORMAT (40H ERROR RETURN 5 FROM LANCZ1. N,LX,LALFA=,3I6)
 
      GO TO 1000
 
  960 IFLAG = 6
      IF (LP.GT.0) WRITE (LP,FMT=965) NEIG
 
  965 FORMAT (38H ERROR RETURN 6 FROM LANCZ1. PREMATURE,
     +       13H TERMINATION.,I6,18H EIGENVALUES FOUND)
 
      GO TO 1000
 
  970 IF (LP.GT.0) WRITE (LP,FMT=975) IFLAG
 
  975 FORMAT (46H ERROR RETURN 7 FROM LANCZ1. ON ENTRY IFLAG IS,I4)
 
      IFLAG = 7
 1000 OLDEL = EL
      OLDER = ER
      RETURN
 
      END
      SUBROUTINE EA15CD(J,LALFA,ALFA,BETA,E,ENORM,ERR,BND,MATCH)
C GIVEN A TRIDIAGONAL MATRIX T GENERATED BY THE LANCZOS PROCESS APPLIED
C     TO A GIVEN SYMMETRIC MATRIX A AND AN APPROXIMATE EIGENVALUE OF T,
C     FIND A BOUND FOR ITS ERROR AS AN EIGENVALUE OF A.
C J IS THE INDEX OF THE LANCZOS STEP. IT IS NOT ALTERED.
C LALFA MUST BE SET TO THE LENGTHS OF ARRAYS ALFA AND BETA. IT IS NOT
C     ALTERED.
C ALFA(I),I=1,J MUST BE SET TO THE DIAGONAL ELEMENTS OF THE LANCZOS
C     TRIDIAGONAL MATRIX T. THEY ARE NOT ALTERED.
C BETA(I),I=2,J+1 MUST BE SET TO THE OFF-DIAGONAL ELEMENTS OF THE
C     LANCZOS TRIDIAGONAL MATRIX T. THEY ARE NOT ALTERED.  NOTE THAT A
C     TYPICAL ROW OF T IS
C          BETA(I)  ALFA(I)  BETA(I+1)
C E MUST BE SET TO THE APPROXIMATE EIGENVALUE OF T. IT IS NOT ALTERED.
C ENORM MUST BE SET TO THE RELATIVE PRECISION TIMES THE NORM OF A TIMES
C     THE ORDER OF A. IT IS NOT ALTERED.
C ERR IS SET TO A BOUND FOR THE ERROR OF E.
C BND IS SET TO AN ESTIMATE OF THE AMOUNT THAT ERR MIGHT BE REDUCED BY
C     USING A MORE ACCURATE E AND/OR A MORE ACCURATE EIGENVECTOR
C     CALCULATION.
C MATCH IS SET TO THE MATCHING POINT BETWEEN USE OF FORWARD AND BACKWARD
C     RECURRENCES.
      DOUBLE PRECISION ALFA(LALFA),BETA(LALFA),E,S,PSI,BND,ERR,ENORM
      DOUBLE PRECISION W,W0,W1,W2,SW,Z1,Z2,BNDK,ERRK,WMAX10
      DOUBLE PRECISION ZERO,ONE
      DATA ZERO/0.0D0/,ONE/1.0D0/
 
      J1 = J - 1
C
C START THE FORWARD RECURRENCE FOR SOLVING (T-E)*Z=PSI*EK
      W1 = ZERO
      W2 = ONE
      SW = ZERO
      K = 0
      W = ALFA(1) - E
      WMAX10 = ZERO
C
      DO 100 KTRY = 1,10
          IF (K-J1) 5,15,110
    5     K1 = K + 1
C CONTINUE THE FORWARD RECURRENCE UNTIL THE FIRST MAXIMUM
C     OR A MAXIMUM AT LEAST TEN TIMES BIGGER THAN THE PREVIOUS
C     MAXIMUM IS FOUND.
          DO 10 K = K1,J1
              W0 = W1
              W1 = W2
              W2 = -W/BETA(K+1)
              SW = SW + W1**2
              W = (ALFA(K+1)-E)*W2 + BETA(K+1)*W1
              IF (ABS(W1).LT.WMAX10) GO TO 10
              IF (ABS(W1).GE.DMAX1(ABS(W0),ABS(W2))) GO TO 20
   10     CONTINUE
   15     K = J
          W0 = W1
          W1 = W2
          SW = SW + W1**2
          IF (ABS(W1).LT.WMAX10) GO TO 100
          IF (ABS(W1).LT.ABS(W0)) GO TO 100
C
C FORWARD RECURRENCE SHOWS A MAXIMUM AT COMPONENT K. CONTINUE
C     SOLVING (T-E)*Z=PSI*EK BY BACKWARD RECURRENCE
   20     S = ZERO
          WMAX10 = ABS(W1)*10.D0
          Z1 = ONE
          PSI = ALFA(J) - E
          IF (K.GT.J1) GO TO 40
          DO 30 KK = K,J1
              I = K + J1 - KK
              Z2 = Z1
              Z1 = -PSI/BETA(I+1)
              S = S + Z2**2
              PSI = (ALFA(I)-E)*Z1 + BETA(I+1)*Z2
   30     CONTINUE
   40     IF (W1.EQ.ZERO) GO TO 100
C
C FORM PRESENT ERROR BOUND AND STORE THE BEST SO FAR FOUND
          Z1 = Z1/W1
          S = S + SW*Z1**2
          S = SQRT(S)
          IF (K.GT.1) PSI = PSI + BETA(K)*Z1*W0
          BNDK = 1.25D0*ABS(PSI)/S
          ERRK = 1.25D0* (ENORM+BETA(J+1)/S) + BNDK
          IF (KTRY.EQ.1) GO TO 50
          IF (ERR.LT.ERRK) GO TO 90
   50     ERR = ERRK
          BND = BNDK
          MATCH = K
C
C TEST RESULT FOR ACCEPTABILITY
   90     IF (BNDK.LT.ERRK/5.D0) GO TO 110
  100 CONTINUE
  110 RETURN
 
      END
      SUBROUTINE EA15DD(N,ALFA,BETA,X,DEL,NU,DR,NUR)
C COMPUTES THE TRIANGULAR FACTORIZATION OF T-X*I, WHERE T IS A
C     TRIDIAGONAL MATRIX, X IS A SCALAR AND I IS THE IDENTITY MATRIX. IT
C     RECORDS THE NUMBER OF NEGATIVE PIVOTS AND THE LAST TWO PIVOTS.
      DOUBLE PRECISION ALFA(N),BETA(N),X,DEL,DR
C N MUST BE SET TO THE ORDER OF THE MATRIX. IT IS NOT ALTERED.
C ALFA MUST BE SET OT CONTAIN THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL
C     MATRIX. IT IS NOT ALTERED.
C BETA(1) MUST HAVE THE VALUE ZERO AND BETA(2),...,BETA(NLAN) MUST BE
C     SET TO CONTAIN THE OFF-DIAGONAL ELEMENTS. BETA IS NOT ALTERED.
C X   MUST BE SET TO THE SCALAR. IT IS NOT ALTERED.
C DEL IS SET TO THE LAST PIVOT.
C NU IS SET TO 1+(THE NUMBER OF NEGATIVE PIVOTS).
C DR IS SET TO THE LAST-BUT-ONE PIVOT.
C NUR IS SET TO 1+(THE NUMBER OF NEGATIVE PIVOTS WHEN THE LAST IS
C     EXCLUDED).
      DOUBLE PRECISION ZERO,ONE
C DRELPR IS THE RELATIVE PRECISION.
      DOUBLE PRECISION DRELPR
CIBM  DATA DRELPR/2.2E-16/
CRAY  DATA DRELPR/1.4E-14/
      DATA DRELPR/2.2D-16/
      DATA ZERO/0.0D0/,ONE/1.0D0/
 
      NU = 1
      DEL = ONE
      DO 10 K = 1,N
          DR = DEL
          DEL = ALFA(K) - X - BETA(K)**2/DEL
          IF (DEL) 6,3,10
    3     DEL = DRELPR* (ABS(ALFA(K))+ABS(X))
          GO TO 10
 
    6     NU = NU + 1
   10 CONTINUE
   20 NUR = NU
      IF (DEL.LT.ZERO) NUR = NU - 1
      RETURN
 
      END
      SUBROUTINE LANCZ2(N,LALFA,LP,ITAPE,EIG,JEIG,NEIG,ALFA,BETA,LY,LZ,
     +                  JFLAG,Y,W,Z)
C FIND APPROXIMATE EIGENVECTORS CORRESPONDING TO EIGENVALUES FOUND
C     BY LANCZ1.
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C LALFA MUST BE SET TO THE LENGTH OF ARRAYS ALFA AND BETA. IT
C     IS NOT ALTERED.
C LP MUST BE SET TO THE UNIT NUMBER FOR DIAGNOSTIC MESSAGES. IF LP.LE.0
C     THE MESSAGES ARE SUPPRESSED. IT IS NOT ALTERED.
C ITAPE MUST BE SET TO THE UNIT NUMBER OF THE SEQUENTIAL FILE USED
C     BY LANCZ1. IT IS NOT ALTERED.
C EIG MUST CONTAIN EIGENVALUES AS FOUND BY LANCZ1. IT IS NOT ALTERED.
C JEIG MUST HOLD THE INFORMATION PASSED DOWN BY LANCZ1, NAMELY JEIG(1,I)
C     MUST CONTAIN  THE LANCZOS STEP AT WHICH EIG(I) WAS ACCEPTED AND
C     JEIG(2,I) MUST CONTAIN THE MATCHING POINT FOR THE EIGENVECTOR
C     CORRESPONDING TO EIG(I). JEIG IS NOT ALTERED.
C NEIG MUST HOLD THE NUMBER OF EIGENVECTORS WANTED. IT IS NOT ALTERED.
C ALFA,BETA MUST BE AS LEFT BY LANCZ1. THEY ARE NOT ALTERED.
C LY MUST HOLD THE FIRST DIMENSION OF Y. IT MUST BE AT LEAST N. IT IS
C     NOT ALTERED.
C LZ MUST HOLD THE FIRST DIMENSION OF Z. IT MUST BE AT LEAST
C     MAX(JEIG(1,I),I=1,NEIG). IT IS NOT ALTERED.
C IFLAG NEED NOT BE SET. ON RETURN IT HAS ONE OF THE VALUES
C      0    SUCCESSFUL COMPLETION
C      1    N.LE.0 .OR. LALFA.LE.0 .OR. ITAPE.LE.0 .OR. NEIG.LE.0
C      2    LY TOO SMALL
C      3    LZ TOO SMALL
C Y IS SET TO THE REQUIRED EIGENVECTORS.
C W IS USED AS A WORKVECTOR.
C Z IS USED AS A WORKVECTOR.
C LZ MUST BE SET TO THE FIRST DIMENSION OF Z. IT MUST BE AT LEAST
C     MAX(JEIG(1,I),I=1,2,...,NEIG).
      DOUBLE PRECISION EIG(NEIG),ALFA(LALFA),BETA(LALFA),Y(LY,NEIG),
     +                 W(N),Z(LZ,NEIG)
      INTEGER JEIG(2,NEIG)
      DOUBLE PRECISION S
      DOUBLE PRECISION ZERO,ONE
      DATA ZERO/0.0D0/,ONE/1.0D0/
 
      JFLAG = 0
      IF (N.LE.0 .OR. LALFA.LE.0 .OR. ITAPE.LE.0 .OR.
     +    NEIG.LE.0) GO TO 100
      IF (LY.LT.N) GO TO 120
      M = 0
      DO 20 L = 1,NEIG
          M = MAX0(M,JEIG(1,L))
          IF (LZ.LT.M) GO TO 140
          DO 10 I = 1,N
              Y(I,L) = ZERO
   10     CONTINUE
          CALL EA15GD(JEIG(1,L),LALFA,ALFA,BETA,EIG(L),JEIG(2,L),Z(1,L))
   20 CONTINUE
      DO 50 J = 1,M
          CALL EA15FD(N,ITAPE,J,2,W)
          DO 40 L = 1,NEIG
              IF (J.GT.JEIG(1,L)) GO TO 40
              DO 30 I = 1,N
                  Y(I,L) = Y(I,L) + Z(J,L)*W(I)
   30         CONTINUE
   40     CONTINUE
   50 CONTINUE
C
C NORMALIZE THE VECTORS
      DO 90 L = 1,NEIG
          S = ZERO
          DO 70 I = 1,N
              S = S + Y(I,L)**2
   70     CONTINUE
          S = ONE/SQRT(S)
          DO 80 I = 1,N
              Y(I,L) = Y(I,L)*S
   80     CONTINUE
   90 CONTINUE
      GO TO 150
 
  100 JFLAG = 1
      IF (LP.GT.0) WRITE (LP,FMT=110) N,LALFA,ITAPE,NEIG
 
  110 FORMAT (28H ERROR RETURN 1 FROM LANCZ2.,/,3H N=,I6,7H LALFA=,I6,
     +       7H ITAPE=,I6,6H NEIG=,I6)
 
      GO TO 150
 
  120 JFLAG = 2
      IF (LP.GT.0) WRITE (LP,FMT=130) LY,N
 
  130 FORMAT (32H ERROR RETURN 2 FROM LANCZ2. LY=,I6,/,
     +       23H AND SHOULD BE AT LEAST,I6)
 
      GO TO 150
 
  140 JFLAG = 3
      IF (LP.GT.0) WRITE (LP,FMT=145) LZ,M
 
  145 FORMAT (32H ERROR RETURN 3 FROM LANCZ2. LZ=,I6,/,
     +       23H AND SHOULD BE AT LEAST,I6)
 
  150 RETURN
 
      END
      SUBROUTINE EA15FD(N,ITAPE,NLAN,IO,V)
C STORE OR RECOVER A LANCZOS VECTOR.
      DOUBLE PRECISION V(N)
C N IS THE ORDER OF THE MATRIX (ALSO THE ORDER OF THE VECTOR BEING
C     STORED OR RECOVERED). IT MUST NOT BE ALTERED.
C ITAPE IS THE UNIT NUMBER OF THE TAPE UNIT. IT IS NOT ALTERED.
C NLAN IS THE LANCZOS ITERATION NUMBER. THE VECTORS ARE ALWAYS SENT FOR
C     STORAGE IN ORDER (NLAN=1,2,...) AND ARE RECOVERED IN ORDER
C     (NLAN=1,2,...). IF THE FIRST VECTOR (NLAN=1) IS SENT FOR STORAGE
C     THEN IT CAN BE ASSUMED THAT A NEW PROBLEM IS IN HAND AND THE OLD
C     VECTORS MAY BE OVERWRITTEN. NLAN MUST NOT BE ALTERED.
C IO IS 1 IF A VECTOR IS TO BE STORED AND TO 2 IF ONE IS TO BE
C     RECOVERED. IT MUST NOT BE ALTERED.
C V HOLDS THE VECTOR TO BE STORED OR IS SET TO THE VECTOR RECOVERED.
      COMMON /EA15HD/LAST
C LAST IS THE VALUE OF NLAN ON THE LAST CALL.
      IF (ITAPE.LE.0) GO TO 50
      IF (NLAN.NE.1) GO TO 10
      REWIND ITAPE
      LAST = 0
C IF NECESSARY, SKIP OVER THE VECTORS WRITTEN EARLIER
   10 NGAP = NLAN - LAST - 1
      IF (NGAP.LE.0) GO TO 30
      DO 20 I = 1,NGAP
          READ (ITAPE)
   20 CONTINUE
C PERFORM ACTUAL READ OR WRITE
   30 IF (IO.EQ.1) WRITE (ITAPE) V
      IF (IO.EQ.2) READ (ITAPE) V
   40 LAST = NLAN
   50 RETURN
 
      END
      SUBROUTINE EA15GD(J,LALFA,ALFA,BETA,E,MATCH,Z)
C FIND AN APPROXIMATE EIGENVECTOR OF A LANCZOS TRIDIAGONAL MATRIX BY
C     FORWARD AND BACKWARD RECURRENCE.
C J IS THE ORDER OF THE MATRIX AND IS NOT ALTERED.
C ALFA,BETA MUST BE EXACTLY AS LEFT BY LANCZ1. THEY ARE NOT ALTERED.
C E IS THE EIGENVALUE CORRESPONDING TO WHICH AN EIGENVECTOR IS WANTED.
C     IT IS NOT ALTERED.
C MATCH IS THE MATCHING POINT BETWEEN THE RECURRENCES. IT IS NOT
C     ALTERED.
C Z IS SET TO THE REQUIRED EIGENVECTOR.
      DOUBLE PRECISION ALFA(LALFA),BETA(LALFA),E,Z(J)
      DOUBLE PRECISION PSI,W,W0,W1,W2,Z1,Z2
      DOUBLE PRECISION ZERO,ONE
      DATA ZERO/0.0D0/,ONE/1.0D0/
 
      J1 = J - 1
C
C FORWARD RECURRENCE.
      W1 = ZERO
      W2 = ONE
      W = ALFA(1) - E
      IF (J1.LE.0) GO TO 15
      DO 10 K = 1,J1
          W0 = W1
          W1 = W2
          W2 = -W/BETA(K+1)
          W = (ALFA(K+1)-E)*W2 + BETA(K+1)*W1
          Z(K) = W1
          IF (K.EQ.MATCH) GO TO 20
   10 CONTINUE
   15 K = J
      W0 = W1
      W1 = W2
      Z(K) = W1
C BACKWARD RECURRENCE
   20 Z1 = ONE
      PSI = ALFA(J) - E
      IF (K.GT.J1) GO TO 40
      DO 30 KK = K,J1
          I = K + J1 - KK
          Z2 = Z1
          Z1 = -PSI/BETA(I+1)
          Z(I+1) = Z2
          PSI = (ALFA(I)-E)*Z1 + BETA(I+1)*Z2
   30 CONTINUE
C RESCALE THE LAST SET OF COMPONENTS
   40 W1 = W1/Z1
      IF (K.GT.J1) GO TO 60
      DO 50 I = K,J1
          Z(I+1) = Z(I+1)*W1
   50 CONTINUE
   60 RETURN
 
      END
      SUBROUTINE LBBT(L,BT,IBT,IW,N)
C
C      THIS SUBROUTINE SOLVES LA=BT
C
      DOUBLE PRECISION C
      DOUBLE PRECISION L(IBT,*),BT(IBT,*)
 
      DO 1 J = 1,N
          BT(J,1) = BT(1,J)/L(1,IW+1)
          IF (J.EQ.1) GO TO 1
          DO 2 I = 2,J
              C = 0.0D0
              IP = I - 1
              DO 3 M = 1,IP
                  IF (M.GT.IW) GO TO 3
    4             C = C - BT(J,I-M)*L(I,IW-M+1)
    3         CONTINUE
              C = C + BT(I,J)
              C = C/L(I,IW+1)
    2     BT(J,I) = C
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE LBKBAN(L,B,K,IK,IW,N)
C
C      THIS SUBROUTINE SOLVES LB=K
C
      DOUBLE PRECISION X
      DOUBLE PRECISION Y
      DOUBLE PRECISION L(IK,*),B(IK,*),K(IK,*)
 
      IP = IW + 1
      DO 1 J = 1,IP
    1 B(1,J) = K(J,IW+2-J)/L(1,IW+1)
      DO 2 I = 2,IP
          IQ = IW + I
          DO 2 J = 1,IQ
              IF (J.GT.N) GO TO 2
              X = 0.0D0
              IR = I - 1
              DO 3 M = 1,IR
    3         X = X + L(I,M-I+IW+1)*B(M,J)
              IF (J.LE.I) Y = K(I,J-I+IW+1)
              IF (J.GT.I) Y = K(J,IW+1-J+I)
              B(I,J) = (Y-X)/L(I,IW+1)
    2 CONTINUE
      IP = IW + 2
      IF (IP.GT.N) GO TO 5
      DO 4 I = IP,N
          DO 4 J = 1,N
              IF (IW+1-J+I.LT.1) GO TO 4
              X = .0D0
              DO 6 M = 1,IW
    6         X = X + L(I,M)*B(I-IW+M-1,J)
              IF (I-J.LE.IW) GO TO 8
              Y = .0D0
              GO TO 7
 
    8         IF (J.GT.I) GO TO 10
              Y = K(I,J-I+IW+1)
              GO TO 7
 
   10         Y = K(J,IW+1-J+I)
    7         B(I,J) = (Y-X)/L(I,IW+1)
    4 CONTINUE
    5 CONTINUE
      RETURN
 
      END
      SUBROUTINE LINMLS(BP,DISPS,LOADS,N,KDIAG)
C
C      THIS SUBROUTINE FORMS THE PRODUCT OF A MATRIX AND A VECTOR
C      WHERE THE MATRIX IS STORED IN A SKYLINE VECTOR
C
      DOUBLE PRECISION X
      DOUBLE PRECISION BP(*),LOADS(*),DISPS(*)
      INTEGER KDIAG(*)
 
      DO 2 I = 1,N
          X = 0.D0
          LUP = KDIAG(I)
          IF (I.EQ.1) LOW = LUP
          IF (I.NE.1) LOW = KDIAG(I-1) + 1
          DO 3 J = LOW,LUP
    3     X = X + BP(J)*DISPS(I+J-LUP)
          LOADS(I) = X
          IF (I.EQ.1) GO TO 2
          LUP = LUP - 1
          DO 4 J = LOW,LUP
              K = I + J - LUP - 1
    4     LOADS(K) = LOADS(K) + BP(J)*DISPS(I)
    2 CONTINUE
      RETURN
 
      END
      SUBROUTINE LINMUL(BK,DISPS,LOADS,N,IW)
C
C      THIS SUBROUTINE MULTIPLIES A MATRIX BY A VECTOR
C      THE MATRIX IS SYMMETRICAL WITH ITS UPPER TRIANGLE
C      STORED AS A VECTOR
C
      DOUBLE PRECISION X
      DOUBLE PRECISION BK(*),DISPS(*),LOADS(*)
 
      DO 1 I = 1,N
          X = 0.D0
          DO 2 J = 1,IW + 1
              IF (I+J.LE.N+1) X = X + BK(N* (J-1)+I)*DISPS(I+J-1)
    2     CONTINUE
          DO 3 J = 2,IW + 1
              IF (I-J+1.GE.1) X = X + BK((N-1)* (J-1)+I)*DISPS(I-J+1)
    3     CONTINUE
          LOADS(I) = X
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE LOC2F(LOCAL,GLOBAL,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE TRANSFORMS THE END REACTIONS AND MOMENTS
C      INTO THE ELEMENT'S LOCAL COORDINATE SYSTEM (2-D)
C
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION ELL
      DOUBLE PRECISION C
      DOUBLE PRECISION S
      DOUBLE PRECISION COORD(ICOORD,*),LOCAL(*),GLOBAL(*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      X2 = COORD(IP,3)
      Y2 = COORD(IP,4)
      ELL = SQRT((X2-X1)**2+ (Y2-Y1)**2)
      C = (X2-X1)/ELL
      S = (Y2-Y1)/ELL
      LOCAL(1) = C*GLOBAL(1) + S*GLOBAL(2)
      LOCAL(2) = C*GLOBAL(2) - S*GLOBAL(1)
      LOCAL(3) = GLOBAL(3)
      LOCAL(4) = C*GLOBAL(4) + S*GLOBAL(5)
      LOCAL(5) = C*GLOBAL(5) - S*GLOBAL(4)
      LOCAL(6) = GLOBAL(6)
      RETURN
 
      END
      SUBROUTINE LOC2T(AXIAL,GLOBAL,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE RETRIEVES THE AXIAL FORCE IN A
C      2-D PIN-JOINTED ELEMENT FROM END REACTIONS (COMP -VE)
C
      DOUBLE PRECISION AXIAL
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION ELL
      DOUBLE PRECISION C
      DOUBLE PRECISION S
      DOUBLE PRECISION GLOBAL(*),COORD(ICOORD,*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      X2 = COORD(IP,3)
      Y2 = COORD(IP,4)
      ELL = SQRT((X2-X1)**2+ (Y2-Y1)**2)
      C = (X2-X1)/ELL
      S = (Y2-Y1)/ELL
      AXIAL = C*GLOBAL(3) + S*GLOBAL(4)
      RETURN
 
      END
      SUBROUTINE LOC3F(LOCAL,GLOBAL,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE TRANSFORMS THE END REACTION AND MOMENTS
C      INTO THE ELEMENT'S LOCAL COORDINATE SYSTEM (3-D)
C
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Z1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION Z2
      DOUBLE PRECISION PI
      DOUBLE PRECISION GAMA
      DOUBLE PRECISION CG
      DOUBLE PRECISION SG
      DOUBLE PRECISION XL
      DOUBLE PRECISION YL
      DOUBLE PRECISION ZL
      DOUBLE PRECISION ELL
      DOUBLE PRECISION DEN
      DOUBLE PRECISION X
      DOUBLE PRECISION SUM
      DOUBLE PRECISION LOCAL(*),GLOBAL(*),COORD(ICOORD,*),R0(3,3),
     +                 T(12,12)
 
      DO 1 I = 1,12
          DO 1 J = 1,12
    1 T(I,J) = 0.D0
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      Z1 = COORD(IP,3)
      X2 = COORD(IP,4)
      Y2 = COORD(IP,5)
      Z2 = COORD(IP,6)
      PI = 4.D0*ATAN(1.D0)
      GAMA = COORD(IP,7)*PI/180.D0
      CG = COS(GAMA)
      SG = SIN(GAMA)
      XL = X2 - X1
      YL = Y2 - Y1
      ZL = Z2 - Z1
      ELL = SQRT(XL*XL+YL*YL+ZL*ZL)
      DEN = ELL*SQRT(XL*XL+ZL*ZL)
      IF (DEN.EQ.0.D0) GO TO 50
      R0(1,1) = XL/ELL
      R0(1,2) = YL/ELL
      R0(1,3) = ZL/ELL
      R0(2,1) = (-XL*YL*CG-ELL*ZL*SG)/DEN
      R0(2,2) = DEN*CG/ (ELL*ELL)
      R0(2,3) = (-YL*ZL*CG+ELL*XL*SG)/DEN
      R0(3,1) = (XL*YL*SG-ELL*ZL*CG)/DEN
      R0(3,2) = -DEN*SG/ (ELL*ELL)
      R0(3,3) = (YL*ZL*SG+ELL*XL*CG)/DEN
      GO TO 60
 
   50 R0(1,1) = 0.D0
      R0(1,3) = 0.D0
      R0(2,2) = 0.D0
      R0(3,2) = 0.D0
      R0(1,2) = 1.D0
      R0(2,1) = -CG
      R0(3,3) = CG
      R0(2,3) = SG
      R0(3,1) = SG
   60 CONTINUE
      DO 2 I = 1,3
          DO 2 J = 1,3
              X = R0(I,J)
              DO 2 K = 0,9,3
    2 T(I+K,J+K) = X
      DO 3 I = 1,12
          SUM = 0.D0
          DO 4 J = 1,12
    4     SUM = SUM + T(I,J)*GLOBAL(J)
    3 LOCAL(I) = SUM
      RETURN
 
      END
      SUBROUTINE LOC3T(AXIAL,GLOBAL,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE RETRIEVES THE AXIAL FORCE IN A
C      3-D PIN-JOINTED ELEMENT FROM END REACTIONS (COMP -VE)
C
      DOUBLE PRECISION AXIAL
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Z1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION Z2
      DOUBLE PRECISION XL
      DOUBLE PRECISION YL
      DOUBLE PRECISION ZL
      DOUBLE PRECISION ELL
      DOUBLE PRECISION GLOBAL(*),COORD(ICOORD,*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      Z1 = COORD(IP,3)
      X2 = COORD(IP,4)
      Y2 = COORD(IP,5)
      Z2 = COORD(IP,6)
      XL = X2 - X1
      YL = Y2 - Y1
      ZL = Z2 - Z1
      ELL = SQRT(XL*XL+YL*YL+ZL*ZL)
      XL = XL/ELL
      YL = YL/ELL
      ZL = ZL/ELL
      AXIAL = XL*GLOBAL(4) + YL*GLOBAL(5) + ZL*GLOBAL(6)
      RETURN
 
      END
      SUBROUTINE MATADD(A,IA,B,IB,M,N)
C
C      THIS SUBROUTINE ADDS TWO EQUAL SIZED ARRAYS
C
      DOUBLE PRECISION A(IA,*),B(IB,*)
 
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 A(I,J) = A(I,J) + B(I,J)
      RETURN
 
      END
      SUBROUTINE MATCOP(A,IA,B,IB,M,N)
C
C      THIS SUBROUTINE COPIES ARRAY A INTO ARRAY B
C
      DOUBLE PRECISION A(IA,*),B(IB,*)
 
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 B(I,J) = A(I,J)
      RETURN
 
      END
      SUBROUTINE MATINV(A,IA,N)
C
C      THIS SUBROUTINE FORMS THE INVERSE OF A MATRIX
C      USING GAUSS-JORDAN TRANSFORMATION
C
      DOUBLE PRECISION CON
      DOUBLE PRECISION A(IA,*)
 
      DO 1 K = 1,N
          CON = A(K,K)
          A(K,K) = 1.D0
          DO 2 J = 1,N
    2     A(K,J) = A(K,J)/CON
          DO 1 I = 1,N
              IF (I.EQ.K) GO TO 1
              CON = A(I,K)
              A(I,K) = 0.D0
              DO 3 J = 1,N
    3         A(I,J) = A(I,J) - A(K,J)*CON
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE MATMUL(A,IA,B,IB,C,IC,L,M,N)
C
C      THIS SUBROUTINE FORMS THE PRODUCT OF TWO MATRICES
C
      DOUBLE PRECISION X
      DOUBLE PRECISION A(IA,*),B(IB,*),C(IC,*)
 
      DO 1 I = 1,L
          DO 1 J = 1,N
              X = 0.0D0
              DO 2 K = 1,M
    2         X = X + A(I,K)*B(K,J)
              C(I,J) = X
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE MATRAN(A,IA,B,IB,M,N)
C
C      THIS SUBROUTINE FORMS THE TRANSPOSE OF A MATRIX
C
      DOUBLE PRECISION A(IA,*),B(IB,*)
 
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 A(J,I) = B(I,J)
      RETURN
 
      END
      SUBROUTINE MOCOPL(PHI,PSI,E,V,STRESS,PL)
C
C      THIS SUBROUTINE FORMS THE PLASTIC STRESS/STRAIN MATRIX
C      FOR A MOHR-COULOMB MATERIAL  (PHI,PSI IN DEGREES)
C
      DOUBLE PRECISION PHI
      DOUBLE PRECISION PSI
      DOUBLE PRECISION E
      DOUBLE PRECISION V
      DOUBLE PRECISION SX
      DOUBLE PRECISION SY
      DOUBLE PRECISION TXY
      DOUBLE PRECISION SZ
      DOUBLE PRECISION PI
      DOUBLE PRECISION PHIR
      DOUBLE PRECISION PSIR
      DOUBLE PRECISION SNPH
      DOUBLE PRECISION SNPS
      DOUBLE PRECISION SQ3
      DOUBLE PRECISION CC
      DOUBLE PRECISION DX
      DOUBLE PRECISION DY
      DOUBLE PRECISION DZ
      DOUBLE PRECISION D2
      DOUBLE PRECISION D3
      DOUBLE PRECISION TH
      DOUBLE PRECISION SNTH
      DOUBLE PRECISION SIG
      DOUBLE PRECISION RPH
      DOUBLE PRECISION RPS
      DOUBLE PRECISION CPS
      DOUBLE PRECISION CPH
      DOUBLE PRECISION EE
      DOUBLE PRECISION ALP
      DOUBLE PRECISION CA
      DOUBLE PRECISION SA
      DOUBLE PRECISION DD
      DOUBLE PRECISION S1
      DOUBLE PRECISION S2
      DOUBLE PRECISION STRESS(4),ROW(4),COL(4),PL(4,4)
 
      SX = STRESS(1)
      SY = STRESS(2)
      TXY = STRESS(3)
      SZ = STRESS(4)
      PI = 4.D0*ATAN(1.D0)
      PHIR = PHI*PI/180.D0
      PSIR = PSI*PI/180.D0
      SNPH = SIN(PHIR)
      SNPS = SIN(PSIR)
      SQ3 = SQRT(3.D0)
      CC = 1.D0 - 2.D0*V
      DX = (2.D0*SX-SY-SZ)/3.D0
      DY = (2.D0*SY-SZ-SX)/3.D0
      DZ = (2.D0*SZ-SX-SY)/3.D0
      D2 = SQRT(-DX*DY-DY*DZ-DZ*DX+TXY*TXY)
      D3 = DX*DY*DZ - DZ*TXY*TXY
      TH = -3.D0*SQ3*D3/ (2.D0*D2**3)
      IF (TH.GT.1.D0) TH = 1.D0
      IF (TH.LT.-1.D0) TH = -1.D0
      TH = ASIN(TH)/3.D0
      SNTH = SIN(TH)
      IF (ABS(SNTH).GT..49D0) THEN
          SIG = -1.D0
          IF (SNTH.LT.0.D0) SIG = 1.D0
          RPH = SNPH* (1.D0+V)/3.D0
          RPS = SNPS* (1.D0+V)/3.D0
          CPS = .25D0*SQ3/D2* (1.D0+SIG*SNPS/3.D0)
          CPH = .25D0*SQ3/D2* (1.D0+SIG*SNPH/3.D0)
          COL(1) = RPH + CPH* ((1.D0-V)*DX+V* (DY+DZ))
          COL(2) = RPH + CPH* ((1.D0-V)*DY+V* (DZ+DX))
          COL(3) = CPH*CC*TXY
          COL(4) = RPH + CPH* ((1.D0-V)*DZ+V* (DX+DY))
          ROW(1) = RPS + CPS* ((1.D0-V)*DX+V* (DY+DZ))
          ROW(2) = RPS + CPS* ((1.D0-V)*DY+V* (DZ+DX))
          ROW(3) = CPS*CC*TXY
          ROW(4) = RPS + CPS* ((1.D0-V)*DZ+V* (DX+DY))
          EE = E/ ((1.D0+V)*CC* (RPH*SNPS+2.D0*CPH*CPS*D2*D2*CC))
 
      ELSE
          ALP = ATAN(ABS((SX-SY)/ (2.D0*TXY)))
          CA = COS(ALP)
          SA = SIN(ALP)
          DD = CC*SA
          S1 = 1.D0
          S2 = 1.D0
          IF ((SX-SY).LT..0D0) S1 = -1.D0
          IF (TXY.LT..0D0) S2 = -1.D0
          COL(1) = SNPH + S1*DD
          COL(2) = SNPH - S1*DD
          COL(3) = S2*CC*CA
          COL(4) = 2.D0*V*SNPH
          ROW(1) = SNPS + S1*DD
          ROW(2) = SNPS - S1*DD
          ROW(3) = S2*CC*CA
          ROW(4) = 2.D0*V*SNPS
          EE = E/ (2.D0* (1.D0+V)*CC* (SNPH*SNPS+CC))
      END IF
 
      DO 1 I = 1,4
          DO 1 J = 1,4
    1 PL(I,J) = EE*ROW(I)*COL(J)
      RETURN
 
      END
      SUBROUTINE MOCOUF(PHI,C,SIGM,DSBAR,THETA,F)
C
C      THIS SUBROUTINE CALCULATES THE VALUE OF THE YIELD FUNCTION
C      FOR A MOHR-COULOMB MATERIAL (PHI IN DEGREES)
C
      DOUBLE PRECISION PHI
      DOUBLE PRECISION C
      DOUBLE PRECISION SIGM
      DOUBLE PRECISION DSBAR
      DOUBLE PRECISION THETA
      DOUBLE PRECISION F
      DOUBLE PRECISION PHIR
      DOUBLE PRECISION SNPH
      DOUBLE PRECISION CSPH
      DOUBLE PRECISION CSTH
      DOUBLE PRECISION SNTH
 
      PHIR = PHI*4.D0*ATAN(1.D0)/180.D0
      SNPH = SIN(PHIR)
      CSPH = COS(PHIR)
      CSTH = COS(THETA)
      SNTH = SIN(THETA)
      F = SNPH*SIGM + DSBAR* (CSTH/SQRT(3.D0)-SNTH*SNPH/3.D0) - C*CSPH
      RETURN
 
      END
      SUBROUTINE MOCOUQ(PSI,DSBAR,THETA,DQ1,DQ2,DQ3)
C
C      THIS SUBROUTINE FORMS THE DERIVATIVES OF A MOHR-COULOMB
C      POTENTIAL FUNCTION WITH RESPECT TO THE THREE INVARIANTS
C      PSI IN DEGREES
C
      DOUBLE PRECISION PSI
      DOUBLE PRECISION DSBAR
      DOUBLE PRECISION THETA
      DOUBLE PRECISION DQ1
      DOUBLE PRECISION DQ2
      DOUBLE PRECISION DQ3
      DOUBLE PRECISION PSIR
      DOUBLE PRECISION SNTH
      DOUBLE PRECISION SNPS
      DOUBLE PRECISION SQ3
      DOUBLE PRECISION C1
      DOUBLE PRECISION CSTH
      DOUBLE PRECISION CS3TH
      DOUBLE PRECISION TN3TH
      DOUBLE PRECISION TNTH
 
      PSIR = PSI*4.D0*ATAN(1.D0)/180.D0
      SNTH = SIN(THETA)
      SNPS = SIN(PSIR)
      SQ3 = SQRT(3.D0)
      DQ1 = SNPS
      IF (ABS(SNTH).GT..49D0) THEN
          C1 = 1.D0
          IF (SNTH.LT.0.D0) C1 = -1.D0
          DQ2 = (SQ3*.5D0-C1*SNPS*.5D0/SQ3)*SQ3*.5D0/DSBAR
          DQ3 = 0.D0
 
      ELSE
          CSTH = COS(THETA)
          CS3TH = COS(3.D0*THETA)
          TN3TH = TAN(3.D0*THETA)
          TNTH = SNTH/CSTH
          DQ2 = SQ3*CSTH/DSBAR* ((1.D0+TNTH*TN3TH)+
     +          SNPS* (TN3TH-TNTH)/SQ3)*.5D0
          DQ3 = 1.5D0* (SQ3*SNTH+SNPS*CSTH)/ (CS3TH*DSBAR*DSBAR)
      END IF
 
      RETURN
 
      END
      SUBROUTINE MSMULT(A,IA,C,M,N)
C
C      THIS SUBROUTINE MULTIPLIES A MATRIX BY A SCALAR
C
      DOUBLE PRECISION C
      DOUBLE PRECISION A(IA,*)
 
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 A(I,J) = A(I,J)*C
      RETURN
 
      END
      SUBROUTINE MVMULT(M,IM,V,K,L,Y)
C
C      THIS SUBROUTINE MULTIPLIES A MATRIX BY A VECTOR
C
      DOUBLE PRECISION X
      DOUBLE PRECISION M(IM,*),V(*),Y(*)
 
      DO 1 I = 1,K
          X = 0.D0
          DO 2 J = 1,L
    2     X = X + M(I,J)*V(J)
          Y(I) = X
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE NULL(A,IA,M,N)
C
C      THIS SUBROUTINE NULLS A 2-D ARRAY
C
      DOUBLE PRECISION A(IA,*)
 
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 A(I,J) = 0.0D0
      RETURN
 
      END
      SUBROUTINE NULL3(A,IA1,IA2,L,M,N)
C
C      THIS SUBROUTINE NULLS A 3-D MATRIX
C
      DOUBLE PRECISION A(IA1,IA2,*)
 
      DO 1 I = 1,L
          DO 1 J = 1,M
              DO 1 K = 1,N
    1 A(I,J,K) = 0.D0
      RETURN
 
      END
      SUBROUTINE NULVEC(VEC,N)
C
C      THIS SUBROUTINE NULLS A COLUMN VECTOR
C
      DOUBLE PRECISION VEC(*)
 
      DO 1 I = 1,N
    1 VEC(I) = 0.D0
      RETURN
 
      END
      SUBROUTINE NUMINT(S,IS,WT,NGP)
C
C      THIS SUBROUTINE FORMS THE SAMPLING POINTS AND
C      WEIGHTS FOR INTEGRATION OVER A TRIANGULAR AREA
C
      DOUBLE PRECISION S(IS,*),WT(*)
 
      GO TO (1,1,3,4,4,6,7,7,7,7,7,12,12,12,12,16) NGP
 
    1 S(1,1) = 1.D0/3.D0
      S(1,2) = 1.D0/3.D0
      WT(1) = 1.D0
      GO TO 99
 
    3 S(1,1) = .5D0
      S(1,2) = .5D0
      S(2,1) = .5D0
      S(2,2) = 0.D0
      S(3,1) = 0.D0
      S(3,2) = .5D0
      WT(1) = 1.D0/3.D0
      WT(2) = WT(1)
      WT(3) = WT(1)
      GO TO 99
 
    4 S(1,1) = 1.D0/3.D0
      S(1,2) = 1.D0/3.D0
      S(2,1) = .6D0
      S(2,2) = .2D0
      S(3,1) = .2D0
      S(3,2) = .6D0
      S(4,1) = .2D0
      S(4,2) = .2D0
      WT(1) = -9.D0/16.D0
      WT(2) = 25.D0/48.D0
      WT(3) = WT(2)
      WT(4) = WT(2)
      GO TO 99
 
    6 S(1,1) = .816847572980459D0
      S(1,2) = .091576213509771D0
      S(2,1) = S(1,2)
      S(2,2) = S(1,1)
      S(3,1) = S(1,2)
      S(3,2) = S(1,2)
      S(4,1) = .108103018168070D0
      S(4,2) = .445948490915965D0
      S(5,1) = S(4,2)
      S(5,2) = S(4,1)
      S(6,1) = S(4,2)
      S(6,2) = S(4,2)
      WT(1) = .109951743655322D0
      WT(2) = WT(1)
      WT(3) = WT(1)
      WT(4) = .223381589678011D0
      WT(5) = WT(4)
      WT(6) = WT(4)
      GO TO 99
 
    7 S(1,1) = 1.D0/3.D0
      S(1,2) = 1.D0/3.D0
      S(2,1) = .797426985353087D0
      S(2,2) = .101286507323456D0
      S(3,1) = S(2,2)
      S(3,2) = S(2,1)
      S(4,1) = S(2,2)
      S(4,2) = S(2,2)
      S(5,1) = .470142064105115D0
      S(5,2) = .059715871789770D0
      S(6,1) = S(5,2)
      S(6,2) = S(5,1)
      S(7,1) = S(5,1)
      S(7,2) = S(5,1)
      WT(1) = .225D0
      WT(2) = .125939180544827D0
      WT(3) = WT(2)
      WT(4) = WT(2)
      WT(5) = .132394152788506D0
      WT(6) = WT(5)
      WT(7) = WT(5)
      GO TO 99
 
   12 S(1,1) = .873821971016996D0
      S(1,2) = .063089014491502D0
      S(2,1) = S(1,2)
      S(2,2) = S(1,1)
      S(3,1) = S(1,2)
      S(3,2) = S(1,2)
      S(4,1) = .501426509658179D0
      S(4,2) = .249286745170910D0
      S(5,1) = S(4,2)
      S(5,2) = S(4,1)
      S(6,1) = S(4,2)
      S(6,2) = S(4,2)
      S(7,1) = .636502499121399D0
      S(7,2) = .310352451033785D0
      S(8,1) = S(7,1)
      S(8,2) = .053145049844816D0
      S(9,1) = S(7,2)
      S(9,2) = S(7,1)
      S(10,1) = S(7,2)
      S(10,2) = S(8,2)
      S(11,1) = S(8,2)
      S(11,2) = S(7,1)
      S(12,1) = S(8,2)
      S(12,2) = S(7,2)
      WT(1) = .050844906370207D0
      WT(2) = WT(1)
      WT(3) = WT(1)
      WT(4) = .116786275726379D0
      WT(5) = WT(4)
      WT(6) = WT(4)
      WT(7) = .082851075618374D0
      WT(8) = WT(7)
      WT(9) = WT(7)
      WT(10) = WT(7)
      WT(11) = WT(7)
      WT(12) = WT(7)
      GO TO 99
 
   16 S(1,1) = 1.D0/3.D0
      S(1,2) = 1.D0/3.D0
      S(2,1) = .658861384496478D0
      S(2,2) = .170569307751761D0
      S(3,1) = S(2,2)
      S(3,2) = S(2,1)
      S(4,1) = S(2,2)
      S(4,2) = S(2,2)
      S(5,1) = .898905543365938D0
      S(5,2) = .050547228317031D0
      S(6,1) = S(5,2)
      S(6,2) = S(5,1)
      S(7,1) = S(5,2)
      S(7,2) = S(5,2)
      S(8,1) = .081414823414554D0
      S(8,2) = .459292588292723D0
      S(9,1) = S(8,2)
      S(9,2) = S(8,1)
      S(10,1) = S(8,2)
      S(10,2) = S(8,2)
      S(11,1) = .008394777409958D0
      S(11,2) = .263112829634638D0
      S(12,1) = S(11,1)
      S(12,2) = .728492392955404D0
      S(13,1) = S(11,2)
      S(13,2) = S(11,1)
      S(14,1) = S(11,2)
      S(14,2) = S(12,2)
      S(15,1) = S(12,2)
      S(15,2) = S(11,1)
      S(16,1) = S(12,2)
      S(16,2) = S(11,2)
      WT(1) = .144315607677787D0
      WT(2) = .103217370534718D0
      WT(3) = WT(2)
      WT(4) = WT(2)
      WT(5) = .032458497623198D0
      WT(6) = WT(5)
      WT(7) = WT(5)
      WT(8) = .095091634267284D0
      WT(9) = WT(8)
      WT(10) = WT(8)
      WT(11) = .027230314174435D0
      WT(12) = WT(11)
      WT(13) = WT(11)
      WT(14) = WT(11)
      WT(15) = WT(11)
      WT(16) = WT(11)
   99 CONTINUE
      RETURN
 
      END
      SUBROUTINE NUMIN3(SAMP,ISAMP,WT,NGP)
C
C      THIS SUBROUTINE FORMS THE SAMPLING POINTS AND
C      WEIGHTS FOR INTEGRATION OVER A TETRAHEDRON
C
      DOUBLE PRECISION SAMP(ISAMP,*),WT(*)
 
      IF (NGP.EQ.1) GO TO 10
      IF (NGP.EQ.4) GO TO 40
      IF (NGP.EQ.5) GO TO 50
   10 SAMP(1,1) = .25D0
      SAMP(1,2) = .25D0
      SAMP(1,3) = .25D0
      WT(1) = 1.D0
      GO TO 99
 
   40 SAMP(1,1) = .58541020D0
      SAMP(1,2) = .13819660D0
      SAMP(1,3) = SAMP(1,2)
      SAMP(2,2) = SAMP(1,1)
      SAMP(2,3) = SAMP(1,2)
      SAMP(2,1) = SAMP(1,2)
      SAMP(3,3) = SAMP(1,1)
      SAMP(3,1) = SAMP(1,2)
      SAMP(3,2) = SAMP(1,2)
      SAMP(4,1) = SAMP(1,2)
      SAMP(4,2) = SAMP(1,2)
      SAMP(4,3) = SAMP(1,2)
      WT(1) = .25D0
      WT(2) = .25D0
      WT(3) = .25D0
      WT(4) = .25D0
      GO TO 99
 
   50 SAMP(1,1) = .25D0
      SAMP(1,2) = .25D0
      SAMP(1,3) = .25D0
      SAMP(2,1) = .5D0
      SAMP(2,2) = 1.D0/6.D0
      SAMP(2,3) = SAMP(2,2)
      SAMP(3,2) = .5D0
      SAMP(3,3) = 1.D0/6.D0
      SAMP(3,1) = SAMP(3,3)
      SAMP(4,3) = .5D0
      SAMP(4,1) = 1.D0/6.D0
      SAMP(4,2) = SAMP(4,1)
      SAMP(5,1) = 1.D0/6.D0
      SAMP(5,2) = SAMP(5,1)
      SAMP(5,3) = SAMP(5,1)
      WT(1) = -.8D0
      WT(2) = 9.D0/20.D0
      WT(3) = WT(2)
      WT(4) = WT(2)
      WT(5) = WT(2)
   99 CONTINUE
      RETURN
 
      END
      SUBROUTINE PINJ2(KM,EA,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX FOR AN
C      INCLINED 2-D PIN-JOINTED ELEMENT
C
      DOUBLE PRECISION EA
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION ELL
      DOUBLE PRECISION CS
      DOUBLE PRECISION SN
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION C
      DOUBLE PRECISION KM(4,4),COORD(ICOORD,*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      X2 = COORD(IP,3)
      Y2 = COORD(IP,4)
      ELL = SQRT((Y2-Y1)**2+ (X2-X1)**2)
      CS = (X2-X1)/ELL
      SN = (Y2-Y1)/ELL
      A = CS*CS
      B = SN*SN
      C = CS*SN
      KM(1,1) = A
      KM(3,3) = A
      KM(1,3) = -A
      KM(3,1) = -A
      KM(2,2) = B
      KM(4,4) = B
      KM(2,4) = -B
      KM(4,2) = -B
      KM(1,2) = C
      KM(2,1) = C
      KM(3,4) = C
      KM(4,3) = C
      KM(1,4) = -C
      KM(4,1) = -C
      KM(2,3) = -C
      KM(3,2) = -C
      DO 1 I = 1,4
          DO 1 J = 1,4
    1 KM(I,J) = KM(I,J)*EA/ELL
      RETURN
 
      END
      SUBROUTINE PINJ3(KM,EA,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX FOR A
C      GENERAL 3-D PIN-JOINTED ELEMENT
C
      DOUBLE PRECISION EA
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Z1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION Z2
      DOUBLE PRECISION XL
      DOUBLE PRECISION YL
      DOUBLE PRECISION ZL
      DOUBLE PRECISION ELL
      DOUBLE PRECISION A
      DOUBLE PRECISION B
      DOUBLE PRECISION C
      DOUBLE PRECISION D
      DOUBLE PRECISION E
      DOUBLE PRECISION F
      DOUBLE PRECISION KM(6,6),COORD(ICOORD,*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      Z1 = COORD(IP,3)
      X2 = COORD(IP,4)
      Y2 = COORD(IP,5)
      Z2 = COORD(IP,6)
      XL = X2 - X1
      YL = Y2 - Y1
      ZL = Z2 - Z1
      ELL = SQRT(XL*XL+YL*YL+ZL*ZL)
      XL = XL/ELL
      YL = YL/ELL
      ZL = ZL/ELL
      A = XL*XL
      B = YL*YL
      C = ZL*ZL
      D = XL*YL
      E = YL*ZL
      F = ZL*XL
      KM(1,1) = A
      KM(4,4) = A
      KM(2,2) = B
      KM(5,5) = B
      KM(3,3) = C
      KM(6,6) = C
      KM(1,2) = D
      KM(2,1) = D
      KM(4,5) = D
      KM(5,4) = D
      KM(2,3) = E
      KM(3,2) = E
      KM(5,6) = E
      KM(6,5) = E
      KM(1,3) = F
      KM(3,1) = F
      KM(4,6) = F
      KM(6,4) = F
      KM(1,4) = -A
      KM(4,1) = -A
      KM(2,5) = -B
      KM(5,2) = -B
      KM(3,6) = -C
      KM(6,3) = -C
      KM(1,5) = -D
      KM(5,1) = -D
      KM(2,4) = -D
      KM(4,2) = -D
      KM(2,6) = -E
      KM(6,2) = -E
      KM(3,5) = -E
      KM(5,3) = -E
      KM(1,6) = -F
      KM(6,1) = -F
      KM(3,4) = -F
      KM(4,3) = -F
      DO 1 I = 1,6
          DO 1 J = 1,6
    1 KM(I,J) = KM(I,J)*EA/ELL
      RETURN
 
      END
      SUBROUTINE PRINTA(A,IA,M,N)
C
C      THIS SUBROUTINE WRITES A 2-D ARRAY TO OUTPUT CHANNEL 6
C
      DOUBLE PRECISION A(IA,*)
 
      DO 1 I = 1,M
    1 WRITE (6,FMT=2) (A(I,J),J=1,N)
 
    2 FORMAT (1X,10D12.4)
 
      RETURN
 
      END
      SUBROUTINE PRINTV(VEC,N)
C
C      THIS SUBROUTINE WRITES A COLUMN VECTOR TO OUTPUT CHANNEL 6
C
      DOUBLE PRECISION VEC(*)
 
      WRITE (6,FMT=1) (VEC(I),I=1,N)
 
    1 FORMAT (1X,6D12.4)
 
      RETURN
 
      END
      SUBROUTINE READNF(NF,INF,NN,NODOF,NR)
C
C      THIS SUBROUTINE READS THE NODAL FREEDOM DATA
C
      INTEGER NF(INF,*)
 
      DO 1 I = 1,NN
          DO 1 J = 1,NODOF
    1 NF(I,J) = 1
      IF (NR.GT.0) READ (5,FMT=*) (K, (NF(K,J),J=1,NODOF),I=1,NR)
      N = 0
      DO 2 I = 1,NN
          DO 2 J = 1,NODOF
              IF (NF(I,J).NE.0) THEN
                  N = N + 1
                  NF(I,J) = N
              END IF
 
    2 CONTINUE
      RETURN
 
      END
      SUBROUTINE SOLVBA(PB,IPB,COPY,ICOPY,ANS,N,IW)
C
C      THIS SUBROUTINE PERFORMS THE GAUSSIAN BACK-SUBSTITUTION
C      ON THE REDUCED MATRIX 'PB'.
C
      DOUBLE PRECISION S
      DOUBLE PRECISION PB(IPB,*),COPY(ICOPY,*),ANS(*)
 
      IWP1 = IW + 1
      IQ = 2*IWP1 - 1
      N1 = N - 1
      DO 1 IV = 1,N1
          I = INT(COPY(IWP1,IV)+.5D0)
          IF (I.EQ.IV) GO TO 2
          S = ANS(IV)
          ANS(IV) = ANS(I)
          ANS(I) = S
    2     L = IV + IWP1 - 1
          IF (L.GT.N) L = N
          IV1 = IV + 1
          DO 3 I = IV1,L
    3     ANS(I) = ANS(I) - COPY(I-IV,IV)*ANS(IV)
    1 CONTINUE
      ANS(N) = ANS(N)/PB(N,1)
      IV = N - 1
    6 S = ANS(IV)
      L = IQ
      IF (IV+L-1.GT.N) L = N - IV + 1
      DO 4 I = 2,L
          S = S - PB(IV,I)*ANS(IV+I-1)
    4 ANS(IV) = S/PB(IV,1)
      IV = IV - 1
      IF (IV.NE.0) GO TO 6
    5 CONTINUE
      RETURN
 
      END
      SUBROUTINE SOLVE(K,IK,U,F,N)
C
C      THIS SUBROUTINE PERFORMS GAUSSIAN ELIMINATION WITH
C      PARTIAL PIVOTING ON A FULL N*N MATRIX
C
      DOUBLE PRECISION BIG
      DOUBLE PRECISION HOLD
      DOUBLE PRECISION FAC
      DOUBLE PRECISION SUM
      DOUBLE PRECISION K(IK,*),F(*),U(*)
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
          SUM = 0.D0
          DO 6 L = I + 1,N
    6     SUM = SUM + K(I,L)*U(L)
          U(I) = (F(I)-SUM)/K(I,I)
    9 CONTINUE
      RETURN
 
      END
      SUBROUTINE SPABAC(A,B,N,KDIAG)
C
C      THIS SUBROUTINE PERFORMS THE CHOLESKI BACK-SUBSTITUTION
C      ON THE VARIABLE BANDWIDTH STIFFNESS MATRIX
C
      DOUBLE PRECISION X
      DOUBLE PRECISION A(*),B(*)
      INTEGER KDIAG(*)
 
      B(1) = B(1)/A(1)
      DO 1 I = 2,N
          KI = KDIAG(I) - I
          L = KDIAG(I-1) - KI + 1
          X = B(I)
          IF (L.EQ.I) GO TO 1
          M = I - 1
          DO 2 J = L,M
    2     X = X - A(KI+J)*B(J)
    1 B(I) = X/A(KI+I)
      DO 3 IT = 2,N
          I = N + 2 - IT
          KI = KDIAG(I) - I
          X = B(I)/A(KI+I)
          B(I) = X
          L = KDIAG(I-1) - KI + 1
          IF (L.EQ.I) GO TO 3
          M = I - 1
          DO 4 K = L,M
    4     B(K) = B(K) - X*A(KI+K)
    3 CONTINUE
      B(1) = B(1)/A(1)
      RETURN
 
      END
      SUBROUTINE SPARIN(A,N,KDIAG)
C
C      THIS SUBROUTINE PERFORMS CHOLESKI REDUCTION OF THE
C      VARIABLE-BANDWIDTH STIFFNESS MATRIX STORED AS A VECTOR
C
      DOUBLE PRECISION X
      DOUBLE PRECISION A(*)
      INTEGER KDIAG(*)
 
      A(1) = SQRT(A(1))
      DO 1 I = 2,N
          KI = KDIAG(I) - I
          L = KDIAG(I-1) - KI + 1
          DO 2 J = L,I
              X = A(KI+J)
              KJ = KDIAG(J) - J
              IF (J.EQ.1) GO TO 2
              LBAR = KDIAG(J-1) - KJ + 1
              LBAR = MAX0(L,LBAR)
              IF (LBAR.EQ.J) GO TO 2
              M = J - 1
              DO 3 K = LBAR,M
    3         X = X - A(KI+K)*A(KJ+K)
    2     A(KI+J) = X/A(KJ+J)
    1 A(KI+I) = SQRT(X)
      RETURN
 
      END
      SUBROUTINE STAB2D(KM,EA,EI,IP,COORD,ICOORD,PAX)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX OF AN
C      INCLINED 2-D BEAM-COLUMN ELEMENT TAKING ACCOUNT
C      OF THE EFFECTS OF AXIAL FORCES
C
      DOUBLE PRECISION EA
      DOUBLE PRECISION EI
      DOUBLE PRECISION PAX
      DOUBLE PRECISION X1
      DOUBLE PRECISION Y1
      DOUBLE PRECISION X2
      DOUBLE PRECISION Y2
      DOUBLE PRECISION ELL
      DOUBLE PRECISION C
      DOUBLE PRECISION S
      DOUBLE PRECISION ALP
      DOUBLE PRECISION SBAR
      DOUBLE PRECISION CBAR
      DOUBLE PRECISION BET1
      DOUBLE PRECISION BET2
      DOUBLE PRECISION E1
      DOUBLE PRECISION E2
      DOUBLE PRECISION E3
      DOUBLE PRECISION E4
      DOUBLE PRECISION KM(6,6),COORD(ICOORD,*)
 
      X1 = COORD(IP,1)
      Y1 = COORD(IP,2)
      X2 = COORD(IP,3)
      Y2 = COORD(IP,4)
      ELL = SQRT((X2-X1)**2+ (Y2-Y1)**2)
      C = (X2-X1)/ELL
      S = (Y2-Y1)/ELL
      ALP = .5D0*ELL*SQRT(ABS(PAX)/EI)
      IF (PAX.GT. (5.D-5*EI/ELL**2)) THEN
          SBAR = ALP* (1.D0-2.D0*ALP/TANH(2.D0*ALP))/ (TANH(ALP)-ALP)
          CBAR = (2.D0*ALP-SINH(2.D0*ALP))/
     +           (SINH(2.D0*ALP)-2.D0*ALP*COSH(2.D0*ALP))
 
      ELSE IF (PAX.LT. (-5.D-5*EI/ELL**2)) THEN
          SBAR = ALP* (1.D0-2.D0*ALP/TAN(2.D0*ALP))/ (TAN(ALP)-ALP)
          CBAR = (2.D0*ALP-SIN(2.D0*ALP))/
     +           (SIN(2.D0*ALP)-2.D0*ALP*COS(2.D0*ALP))
 
      ELSE
          SBAR = 4.D0
          CBAR = .5D0
      END IF
 
      BET1 = 2.D0*SBAR* (1.D0+CBAR) + PAX*ELL**2/EI
      BET2 = SBAR* (1.D0+CBAR)
      E1 = EA/ELL
      E2 = BET1*EI/ (ELL*ELL*ELL)
      E3 = EI/ELL
      E4 = BET2*EI/ (ELL*ELL)
      KM(1,1) = C*C*E1 + S*S*E2
      KM(4,4) = KM(1,1)
      KM(1,2) = S*C* (E1-E2)
      KM(2,1) = KM(1,2)
      KM(4,5) = KM(1,2)
      KM(5,4) = KM(4,5)
      KM(1,3) = -S*E4
      KM(3,1) = KM(1,3)
      KM(1,6) = KM(1,3)
      KM(6,1) = KM(1,6)
      KM(3,4) = S*E4
      KM(4,3) = KM(3,4)
      KM(4,6) = KM(3,4)
      KM(6,4) = KM(4,6)
      KM(1,4) = -KM(1,1)
      KM(4,1) = KM(1,4)
      KM(1,5) = S*C* (-E1+E2)
      KM(5,1) = KM(1,5)
      KM(2,4) = KM(1,5)
      KM(4,2) = KM(2,4)
      KM(2,2) = S*S*E1 + C*C*E2
      KM(5,5) = KM(2,2)
      KM(2,5) = -KM(2,2)
      KM(5,2) = KM(2,5)
      KM(2,3) = C*E4
      KM(3,2) = KM(2,3)
      KM(2,6) = KM(2,3)
      KM(6,2) = KM(2,6)
      KM(3,3) = SBAR*E3
      KM(6,6) = KM(3,3)
      KM(3,5) = -C*E4
      KM(5,3) = KM(3,5)
      KM(5,6) = KM(3,5)
      KM(6,5) = KM(5,6)
      KM(3,6) = SBAR*CBAR*E3
      KM(6,3) = KM(3,6)
      RETURN
 
      END
      SUBROUTINE TREEX3(JAC,IJAC,JAC1,IJAC1,DET)
C
C      THIS SUBROUTINE FORMS THE INVERSE OF A 3 BY 3 MATRIX
C
      DOUBLE PRECISION DET
      DOUBLE PRECISION JAC(IJAC,*),JAC1(IJAC1,*)
 
      DET = JAC(1,1)* (JAC(2,2)*JAC(3,3)-JAC(3,2)*JAC(2,3))
      DET = DET - JAC(1,2)* (JAC(2,1)*JAC(3,3)-JAC(3,1)*JAC(2,3))
      DET = DET + JAC(1,3)* (JAC(2,1)*JAC(3,2)-JAC(3,1)*JAC(2,2))
      JAC1(1,1) = JAC(2,2)*JAC(3,3) - JAC(3,2)*JAC(2,3)
      JAC1(2,1) = -JAC(2,1)*JAC(3,3) + JAC(3,1)*JAC(2,3)
      JAC1(3,1) = JAC(2,1)*JAC(3,2) - JAC(3,1)*JAC(2,2)
      JAC1(1,2) = -JAC(1,2)*JAC(3,3) + JAC(3,2)*JAC(1,3)
      JAC1(2,2) = JAC(1,1)*JAC(3,3) - JAC(3,1)*JAC(1,3)
      JAC1(3,2) = -JAC(1,1)*JAC(3,2) + JAC(3,1)*JAC(1,2)
      JAC1(1,3) = JAC(1,2)*JAC(2,3) - JAC(2,2)*JAC(1,3)
      JAC1(2,3) = -JAC(1,1)*JAC(2,3) + JAC(2,1)*JAC(1,3)
      JAC1(3,3) = JAC(1,1)*JAC(2,2) - JAC(2,1)*JAC(1,2)
      DO 1 K = 1,3
          DO 1 L = 1,3
              JAC1(K,L) = JAC1(K,L)/DET
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE TRIDIA(N,ATOL,A,D,E,Z,IZ)
C
C     THIS SUBROUTINE REDUCES THE GIVEN LOWER TRIANGLE OF A
C     SYMMETRIC MATRIX, A, STORED IN THE ARRAY A(N,N), TO
C     TRIDIAGONAL FORM USING HOUSEHOLDERS REDUCTION. THE DIAGONAL
C     OF THE RESULT IS STORED IN THE ARRAY D(N) AND THE
C     SUB-DIAGONAL IN THE LAST N - 1 STORES OF THE ARRAY E(N)
C     (WITH THE ADDITIONAL ELEMENT E(1) = 0). THE TRANSFORMATION
C     MATRICES ARE ACCUMULATED IN THE ARRAY Z(N,N). THE ARRAY
C     A IS LEFT UNALTERED UNLESS THE ACTUAL PARAMETERS
C     CORRESPONDING TO A AND Z ARE IDENTICAL.
C
      DOUBLE PRECISION ATOL
      DOUBLE PRECISION F
      DOUBLE PRECISION G
      DOUBLE PRECISION HH
      INTEGER N,I,L,IM1,IM2,J,K,IA,IZ
      DOUBLE PRECISION H,S,A(IZ,*),D(*),E(*),Z(IZ,*)
 
      DO 40 I = 1,N
          DO 20 J = 1,I
              Z(I,J) = A(I,J)
   20     CONTINUE
   40 CONTINUE
      IF (N.EQ.1) GO TO 280
      DO 260 II = 2,N
          I = N - II + 2
          L = I - 2
          F = Z(I,I-1)
          G = 0.0D0
          IF (L.EQ.0) GO TO 80
          DO 60 K = 1,L
              G = G + Z(I,K)*Z(I,K)
   60     CONTINUE
   80     H = G + F*F
C
C     IF G IS TOO SMALL FOR ORTHOGONALITY TO BE
C     GUARANTEED THE TRANSFORMATION IS SKIPPED
C
          IF (G.GT.ATOL) GO TO 100
          E(I) = F
          H = 0.0D0
          GO TO 240
 
  100     L = L + 1
          G = SQRT(H)
          IF (F.GE.0.0D0) G = -G
          E(I) = G
          H = H - F*G
          Z(I,I-1) = F - G
          F = 0.0D0
          DO 180 J = 1,L
              Z(J,I) = Z(I,J)/H
              G = 0.0D0
C
C     FORM ELEMENT OF A*U
C
              DO 120 K = 1,J
                  G = G + Z(J,K)*Z(I,K)
  120         CONTINUE
              J1 = J + 1
              IF (J1.GT.L) GO TO 160
              DO 140 K = J1,L
                  G = G + Z(K,J)*Z(I,K)
  140         CONTINUE
C
C     FORM ELEMENT OF P
C
  160         E(J) = G/H
              F = F + G*Z(J,I)
  180     CONTINUE
C
C     FORM K
C
          HH = F/ (H+H)
C
C     FORM REDUCED A
C
          DO 220 J = 1,L
              F = Z(I,J)
              G = E(J) - HH*F
              E(J) = G
              DO 200 K = 1,J
                  Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
  200         CONTINUE
  220     CONTINUE
  240     D(I) = H
  260 CONTINUE
  280 E(1) = 0.0D0
      D(1) = 0.0D0
C
C     ACCUMULATION OF TRANSFORMATION MATRICES
C
      DO 400 I = 1,N
          L = I - 1
          IF (D(I).EQ.0.0D0) GO TO 360
          DO 340 J = 1,L
              G = 0.0D0
              DO 300 K = 1,L
                  G = G + Z(I,K)*Z(K,J)
  300         CONTINUE
              DO 320 K = 1,L
                  Z(K,J) = Z(K,J) - G*Z(K,I)
  320         CONTINUE
  340     CONTINUE
  360     D(I) = Z(I,I)
          Z(I,I) = 1.0D0
          IF (L.EQ.0) GO TO 400
          DO 380 J = 1,L
              Z(I,J) = 0.0D0
              Z(J,I) = 0.0D0
  380     CONTINUE
  400 CONTINUE
      RETURN
 
      END
      SUBROUTINE TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)
C
C      THIS SUBROUTINE FORMS THE INVERSE OF A 2 BY 2 MATRIX
C
      DOUBLE PRECISION DET
      DOUBLE PRECISION JAC(IJAC,*),JAC1(IJAC1,*)
 
      DET = JAC(1,1)*JAC(2,2) - JAC(1,2)*JAC(2,1)
      JAC1(1,1) = JAC(2,2)
      JAC1(1,2) = -JAC(1,2)
      JAC1(2,1) = -JAC(2,1)
      JAC1(2,2) = JAC(1,1)
      DO 1 K = 1,2
          DO 1 L = 1,2
    1 JAC1(K,L) = JAC1(K,L)/DET
      RETURN
 
      END
      SUBROUTINE VECADD(A,B,C,N)
C
C      THIS SUBROUTINE ADDS VECTORS  A+B=C
C
      DOUBLE PRECISION A(*),B(*),C(*)
 
      DO 1 I = 1,N
    1 C(I) = A(I) + B(I)
      RETURN
 
      END
      SUBROUTINE VECCOP(A,B,N)
C
C      THIS SUBROUTINE COPIES VECTOR A INTO VECTOR B
C
      DOUBLE PRECISION A(*),B(*)
 
      DO 1 I = 1,N
    1 B(I) = A(I)
      RETURN
 
      END
      SUBROUTINE VMPL(E,V,STRESS,PL)
C
C      THIS SUBROUTINE FORMS THE PLASTIC MATRIX FOR A
C      VON-MISES MATERIAL
C
      DOUBLE PRECISION E
      DOUBLE PRECISION V
      DOUBLE PRECISION SX
      DOUBLE PRECISION SY
      DOUBLE PRECISION TXY
      DOUBLE PRECISION SZ
      DOUBLE PRECISION DSBAR
      DOUBLE PRECISION EE
      DOUBLE PRECISION STRESS(*),TERM(4),PL(4,4)
 
      SX = STRESS(1)
      SY = STRESS(2)
      TXY = STRESS(3)
      SZ = STRESS(4)
      DSBAR = SQRT((SX-SY)**2+ (SY-SZ)**2+ (SZ-SX)**2+6.D0*TXY**2)/
     +        SQRT(2.D0)
      EE = 1.5D0*E/ ((1.D0+V)*DSBAR*DSBAR)
      TERM(1) = (2.D0*SX-SY-SZ)/3.D0
      TERM(2) = (2.D0*SY-SZ-SX)/3.D0
      TERM(3) = TXY
      TERM(4) = (2.D0*SZ-SX-SY)/3.D0
      DO 1 I = 1,4
          DO 1 J = I,4
              PL(I,J) = TERM(I)*TERM(J)*EE
    1 PL(J,I) = PL(I,J)
      RETURN
 
      END
      SUBROUTINE VOL2D(BEE,IBEE,VOL,NOD)
C
C      THIS SUBROUTINE FORMS A VECTOR CONTAINING THE
C      DERIVATIVES OF THE SHAPE FUNCTIONS (PLANE 2-D)
C
      DOUBLE PRECISION BEE(IBEE,*),VOL(*)
 
      DO 1 M = 1,NOD
          K = 2*M
          L = K - 1
          VOL(L) = BEE(1,L)
          VOL(K) = BEE(2,K)
    1 CONTINUE
      RETURN
 
      END
      SUBROUTINE VVMULT(V1,V2,PROD,IPROD,M,N)
C
C      THIS SUBROUTINE FORMS A VECTOR PRODUCT
C
      DOUBLE PRECISION V1(*),V2(*),PROD(IPROD,*)
 
      DO 1 I = 1,M
          DO 1 J = 1,N
    1 PROD(I,J) = V1(I)*V2(J)
      RETURN
 
      END
      SUBROUTINE GENA8X(IP,IQ,NRE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE QUADS COUNTING IN X-DIRECTION (3-FREEDOMS/NODE)
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(8)
 
      NUM(1) = IQ* (3*NRE+2) + 2*IP - 1
      NUM(2) = IQ* (3*NRE+2) + IP - NRE - 1
      NUM(3) = (IQ-1)* (3*NRE+2) + 2*IP - 1
      NUM(4) = NUM(3) + 1
      NUM(5) = NUM(4) + 1
      NUM(6) = NUM(2) + 1
      NUM(7) = NUM(1) + 2
      NUM(8) = NUM(1) + 1
      INC = 0
      DO 1 I = 1,8
          DO 1 J = 1,3
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = AA* (IP-1)
      COORD(2,1) = AA* (IP-1)
      COORD(3,1) = AA* (IP-1)
      COORD(5,1) = AA*IP
      COORD(6,1) = AA*IP
      COORD(7,1) = AA*IP
      COORD(4,1) = .5D0* (COORD(3,1)+COORD(5,1))
      COORD(8,1) = .5D0* (COORD(7,1)+COORD(1,1))
      COORD(1,2) = -BB*IQ
      COORD(8,2) = -BB*IQ
      COORD(7,2) = -BB*IQ
      COORD(3,2) = -BB* (IQ-1)
      COORD(4,2) = -BB* (IQ-1)
      COORD(5,2) = -BB* (IQ-1)
      COORD(2,2) = .5D0* (COORD(1,2)+COORD(3,2))
      COORD(6,2) = .5D0* (COORD(5,2)+COORD(7,2))
      RETURN
 
      END
      SUBROUTINE GEOM3X(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 3-NODE TRIANGLES COUNTING IN X-DIRECTION
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(3)
 
      IF (MOD(IQ,2).EQ.0) GO TO 1
      NUM(1) = (NXE+1)* (IQ-1)/2 + IP
      NUM(2) = (NXE+1)* (IQ+1)/2 + IP
      NUM(3) = NUM(1) + 1
      COORD(1,1) = (IP-1)*AA
      COORD(1,2) = - (IQ-1)/2*BB
      COORD(2,1) = (IP-1)*AA
      COORD(2,2) = - (IQ+1)/2*BB
      COORD(3,1) = IP*AA
      COORD(3,2) = COORD(1,2)
      GO TO 2
 
    1 NUM(1) = (NXE+1)*IQ/2 + IP + 1
      NUM(2) = (NXE+1)* (IQ-2)/2 + IP + 1
      NUM(3) = NUM(1) - 1
      COORD(1,1) = IP*AA
      COORD(1,2) = -IQ/2*BB
      COORD(2,1) = IP*AA
      COORD(2,2) = - (IQ-2)/2*BB
      COORD(3,1) = (IP-1)*AA
      COORD(3,2) = COORD(1,2)
    2 CONTINUE
      INC = 0
      DO 3 I = 1,3
          DO 3 J = 1,2
              INC = INC + 1
    3 G(INC) = NF(NUM(I),J)
      RETURN
 
      END
      SUBROUTINE GEOM4Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 4-NODE QUADS COUNTING IN Y-DIRECTION
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(4)
 
      NUM(1) = (IP-1)* (NYE+1) + IQ + 1
      NUM(2) = NUM(1) - 1
      NUM(3) = IP* (NYE+1) + IQ
      NUM(4) = NUM(3) + 1
      INC = 0
      DO 1 I = 1,4
          DO 1 J = 1,2
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = AA* (IP-1)
      COORD(2,1) = AA* (IP-1)
      COORD(3,1) = AA*IP
      COORD(4,1) = AA*IP
      COORD(1,2) = -BB*IQ
      COORD(2,2) = -BB* (IQ-1)
      COORD(3,2) = -BB* (IQ-1)
      COORD(4,2) = -BB*IQ
      RETURN
 
      END
      SUBROUTINE GEOM6X(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 6-NODE TRIANGLES COUNTING IN X-DIRECTION
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(6)
 
      IF (MOD(IQ,2).EQ.0) GO TO 1
      NUM(1) = (IQ-1)* (2*NXE+1) + 2*IP - 1
      NUM(2) = (IQ-1)* (2*NXE+1) + 2*NXE + 2*IP
      NUM(3) = (IQ+1)* (2*NXE+1) + 2*IP - 1
      NUM(4) = NUM(2) + 1
      NUM(5) = NUM(1) + 2
      NUM(6) = NUM(1) + 1
      COORD(1,1) = (IP-1)*AA
      COORD(1,2) = - (IQ-1)/2*BB
      COORD(3,1) = (IP-1)*AA
      COORD(3,2) = - (IQ+1)/2*BB
      COORD(5,1) = IP*AA
      COORD(5,2) = COORD(1,2)
      GO TO 2
 
    1 NUM(1) = IQ* (2*NXE+1) + 2*IP + 1
      NUM(2) = (IQ-2)* (2*NXE+1) + 2*NXE + 2*IP + 2
      NUM(3) = (IQ-2)* (2*NXE+1) + 2*IP + 1
      NUM(4) = NUM(2) - 1
      NUM(5) = NUM(1) - 2
      NUM(6) = NUM(1) - 1
      COORD(1,1) = IP*AA
      COORD(1,2) = -IQ/2*BB
      COORD(3,1) = IP*AA
      COORD(3,2) = - (IQ-2)/2*BB
      COORD(5,1) = (IP-1)*AA
      COORD(5,2) = COORD(1,2)
    2 DO 3 I = 1,2
          COORD(2,I) = .5D0* (COORD(1,I)+COORD(3,I))
          COORD(4,I) = .5D0* (COORD(3,I)+COORD(5,I))
    3 COORD(6,I) = .5D0* (COORD(5,I)+COORD(1,I))
      INC = 0
      DO 4 I = 1,6
          DO 4 J = 1,2
              INC = INC + 1
    4 G(INC) = NF(NUM(I),J)
      RETURN
 
      END
      SUBROUTINE GEOM8X(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE QUADS COUNTING IN X-DIRECTION
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(8)
 
      NUM(1) = IQ* (3*NXE+2) + 2*IP - 1
      NUM(2) = IQ* (3*NXE+2) + IP - NXE - 1
      NUM(3) = (IQ-1)* (3*NXE+2) + 2*IP - 1
      NUM(4) = NUM(3) + 1
      NUM(5) = NUM(4) + 1
      NUM(6) = NUM(2) + 1
      NUM(7) = NUM(1) + 2
      NUM(8) = NUM(1) + 1
      INC = 0
      DO 1 I = 1,8
          DO 1 J = 1,2
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = AA* (IP-1)
      COORD(2,1) = AA* (IP-1)
      COORD(3,1) = AA* (IP-1)
      COORD(5,1) = AA*IP
      COORD(6,1) = AA*IP
      COORD(7,1) = AA*IP
      COORD(4,1) = .5D0* (COORD(3,1)+COORD(5,1))
      COORD(8,1) = .5D0* (COORD(7,1)+COORD(1,1))
      COORD(1,2) = -BB*IQ
      COORD(8,2) = -BB*IQ
      COORD(7,2) = -BB*IQ
      COORD(3,2) = -BB* (IQ-1)
      COORD(4,2) = -BB* (IQ-1)
      COORD(5,2) = -BB* (IQ-1)
      COORD(2,2) = .5D0* (COORD(1,2)+COORD(3,2))
      COORD(6,2) = .5D0* (COORD(5,2)+COORD(7,2))
      RETURN
 
      END
      SUBROUTINE GEOM8Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE QUADS COUNTING IN THE Y-DIRECTION
C
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER NUM(8),G(*),NF(INF,*)
 
      NUM(1) = (IP-1)* (3*NYE+2) + 2*IQ + 1
      NUM(2) = NUM(1) - 1
      NUM(3) = NUM(2) - 1
      NUM(4) = (IP-1)* (3*NYE+2) + 2*NYE + IQ + 1
      NUM(5) = IP* (3*NYE+2) + 2*IQ - 1
      NUM(6) = NUM(5) + 1
      NUM(7) = NUM(6) + 1
      NUM(8) = NUM(4) + 1
      INC = 0
      DO 1 I = 1,8
          DO 1 J = 1,2
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = (IP-1)*AA
      COORD(2,1) = (IP-1)*AA
      COORD(3,1) = (IP-1)*AA
      COORD(5,1) = IP*AA
      COORD(6,1) = IP*AA
      COORD(7,1) = IP*AA
      COORD(4,1) = (COORD(3,1)+COORD(5,1))*.5D0
      COORD(8,1) = COORD(4,1)
      COORD(3,2) = - (IQ-1)*BB
      COORD(4,2) = - (IQ-1)*BB
      COORD(5,2) = - (IQ-1)*BB
      COORD(1,2) = -IQ*BB
      COORD(8,2) = -IQ*BB
      COORD(7,2) = -IQ*BB
      COORD(2,2) = (COORD(1,2)+COORD(3,2))*.5D0
      COORD(6,2) = COORD(2,2)
      RETURN
 
      END
      SUBROUTINE GEOM9X(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 9-NODE QUADS COUNTING IN X-DIRECTION
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(9)
 
      NUM(1) = IQ* (4*NXE+2) + 2*IP - 1
      NUM(2) = IQ* (4*NXE+2) + 2*IP - NXE - 4
      NUM(3) = (IQ-1)* (4*NXE+2) + 2*IP - 1
      NUM(4) = NUM(3) + 1
      NUM(5) = NUM(4) + 1
      NUM(6) = NUM(2) + 2
      NUM(7) = NUM(1) + 2
      NUM(8) = NUM(1) + 1
      NUM(9) = NUM(2) + 1
      INC = 0
      DO 1 I = 1,9
          DO 1 J = 1,2
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = (IP-1)*AA
      COORD(3,1) = (IP-1)*AA
      COORD(5,1) = IP*AA
      COORD(7,1) = IP*AA
      COORD(1,2) = -IQ*BB
      COORD(3,2) = - (IQ-1)*BB
      COORD(5,2) = - (IQ-1)*BB
      COORD(7,2) = -IQ*BB
      COORD(2,1) = .5D0* (COORD(1,1)+COORD(3,1))
      COORD(2,2) = .5D0* (COORD(1,2)+COORD(3,2))
      COORD(4,1) = .5D0* (COORD(3,1)+COORD(5,1))
      COORD(4,2) = .5D0* (COORD(3,2)+COORD(5,2))
      COORD(6,1) = .5D0* (COORD(5,1)+COORD(7,1))
      COORD(6,2) = .5D0* (COORD(5,2)+COORD(7,2))
      COORD(8,1) = .5D0* (COORD(1,1)+COORD(7,1))
      COORD(8,2) = .5D0* (COORD(1,2)+COORD(7,2))
      COORD(9,1) = .5D0* (COORD(2,1)+COORD(6,1))
      COORD(9,2) = .5D0* (COORD(4,2)+COORD(8,2))
      RETURN
 
      END
      SUBROUTINE GEOUPV(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,COORDF,
     +                  ICORDF,G,NF,INF)
C
C       THIS SUBROUTINE FORMS THE NODAL COORDINATES AND
C       STEERING VECTORFOR A RECTANGULAR MESH OF 4-NODE/8-NODE
C       QUADRILATERAL ELEMENTS NUMBERING IN THE X-DIRECTION
C       (U,P,V NAVIER STOKES FLOW)
C
      DOUBLE PRECISION COORD(ICOORD,*),COORDF(ICORDF,*),WIDTH(*),
     +                 DEPTH(*)
      INTEGER NUM(8),G(*),NF(INF,*)
 
      NUM(1) = IQ* (3*NXE+2) + 2*IP - 1
      NUM(2) = IQ* (3*NXE+2) + IP - NXE - 1
      NUM(3) = (IQ-1)* (3*NXE+2) + 2*IP - 1
      NUM(4) = NUM(3) + 1
      NUM(5) = NUM(4) + 1
      NUM(6) = NUM(2) + 1
      NUM(7) = NUM(1) + 2
      NUM(8) = NUM(1) + 1
      INC = 0
      DO 1 I = 1,8
          INC = INC + 1
    1 G(INC) = NF(NUM(I),1)
      DO 2 I = 1,7,2
          INC = INC + 1
    2 G(INC) = NF(NUM(I),2)
      DO 3 I = 1,8
          INC = INC + 1
    3 G(INC) = NF(NUM(I),3)
      COORD(1,1) = WIDTH(IP)
      COORD(2,1) = WIDTH(IP)
      COORD(3,1) = WIDTH(IP)
      COORDF(1,1) = WIDTH(IP)
      COORDF(2,1) = WIDTH(IP)
      COORD(5,1) = WIDTH(IP+1)
      COORD(6,1) = WIDTH(IP+1)
      COORD(7,1) = WIDTH(IP+1)
      COORDF(3,1) = WIDTH(IP+1)
      COORDF(4,1) = WIDTH(IP+1)
      COORD(4,1) = .5D0* (COORD(3,1)+COORD(5,1))
      COORD(8,1) = .5D0* (COORD(7,1)+COORD(1,1))
      COORD(1,2) = DEPTH(IQ+1)
      COORD(8,2) = DEPTH(IQ+1)
      COORD(7,2) = DEPTH(IQ+1)
      COORDF(1,2) = DEPTH(IQ+1)
      COORDF(4,2) = DEPTH(IQ+1)
      COORD(3,2) = DEPTH(IQ)
      COORD(4,2) = DEPTH(IQ)
      COORD(5,2) = DEPTH(IQ)
      COORDF(2,2) = DEPTH(IQ)
      COORDF(3,2) = DEPTH(IQ)
      COORD(2,2) = .5D0* (COORD(1,2)+COORD(3,2))
      COORD(6,2) = .5D0* (COORD(5,2)+COORD(7,2))
      RETURN
 
      END
      SUBROUTINE GEOUVP(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,COORDF,
     +                  ICORDF,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE NODAL COORDINATES AND STEERING
C      VECTOR FOR A VARIABLE MESH OF 4-NODE/8-NODE
C      QUADRILATERAL ELEMENTS NUMBERING IN THE X-DIRECTION
C      (U,V,P  BIOT CONSOLIDATION)
C
      DOUBLE PRECISION COORD(ICOORD,*),COORDF(ICORDF,*),WIDTH(*),
     +                 DEPTH(*)
      INTEGER NUM(8),G(*),NF(INF,*)
 
      NUM(1) = IQ* (3*NXE+2) + 2*IP - 1
      NUM(2) = IQ* (3*NXE+2) + IP - NXE - 1
      NUM(3) = (IQ-1)* (3*NXE+2) + 2*IP - 1
      NUM(4) = NUM(3) + 1
      NUM(5) = NUM(4) + 1
      NUM(6) = NUM(2) + 1
      NUM(7) = NUM(1) + 2
      NUM(8) = NUM(1) + 1
      INC = 0
      DO 1 I = 1,8
          DO 1 J = 1,2
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      DO 2 I = 1,7,2
          INC = INC + 1
    2 G(INC) = NF(NUM(I),3)
      COORD(1,1) = WIDTH(IP)
      COORD(2,1) = WIDTH(IP)
      COORD(3,1) = WIDTH(IP)
      COORDF(1,1) = WIDTH(IP)
      COORDF(2,1) = WIDTH(IP)
      COORD(5,1) = WIDTH(IP+1)
      COORD(6,1) = WIDTH(IP+1)
      COORD(7,1) = WIDTH(IP+1)
      COORDF(3,1) = WIDTH(IP+1)
      COORDF(4,1) = WIDTH(IP+1)
      COORD(4,1) = .5D0* (COORD(3,1)+COORD(5,1))
      COORD(8,1) = .5D0* (COORD(7,1)+COORD(1,1))
      COORD(1,2) = DEPTH(IQ+1)
      COORD(8,2) = DEPTH(IQ+1)
      COORD(7,2) = DEPTH(IQ+1)
      COORDF(1,2) = DEPTH(IQ+1)
      COORDF(4,2) = DEPTH(IQ+1)
      COORD(3,2) = DEPTH(IQ)
      COORD(4,2) = DEPTH(IQ)
      COORD(5,2) = DEPTH(IQ)
      COORDF(2,2) = DEPTH(IQ)
      COORDF(3,2) = DEPTH(IQ)
      COORD(2,2) = .5D0* (COORD(1,2)+COORD(3,2))
      COORD(6,2) = .5D0* (COORD(5,2)+COORD(7,2))
      RETURN
 
      END
      SUBROUTINE GEOV4Y(IP,IQ,NDE,RAD,DEP,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR EACH ELEMENT (NUMBERING IN THE Y-DIRECTION)
C
      DOUBLE PRECISION COORD(ICOORD,*),RAD(*),DEP(*)
      INTEGER G(*),NF(INF,*),NUM(4)
 
      NUM(1) = (IP-1)* (NDE+1) + IQ + 1
      NUM(2) = NUM(1) - 1
      NUM(3) = IP* (NDE+1) + IQ
      NUM(4) = NUM(3) + 1
      INC = 0
      DO 1 I = 1,4
          DO 1 J = 1,2
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = RAD(IP)
      COORD(2,1) = RAD(IP)
      COORD(3,1) = RAD(IP+1)
      COORD(4,1) = RAD(IP+1)
      COORD(1,2) = DEP(IQ+1)
      COORD(2,2) = DEP(IQ)
      COORD(3,2) = DEP(IQ)
      COORD(4,2) = DEP(IQ+1)
      RETURN
 
      END
      SUBROUTINE GEOV8Y(IP,IQ,NYE,WIDTH,DEPTH,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE QUADRILATERALS NUMBERING IN Y-DIRECTION
C
      DOUBLE PRECISION COORD(ICOORD,*),WIDTH(*),DEPTH(*)
      INTEGER G(*),NF(INF,*),NUM(8)
 
      NUM(1) = (IP-1)* (3*NYE+2) + 2*IQ + 1
      NUM(2) = NUM(1) - 1
      NUM(3) = NUM(1) - 2
      NUM(4) = (IP-1)* (3*NYE+2) + 2*NYE + IQ + 1
      NUM(5) = IP* (3*NYE+2) + 2*IQ - 1
      NUM(6) = NUM(5) + 1
      NUM(7) = NUM(5) + 2
      NUM(8) = NUM(4) + 1
      INC = 0
      DO 1 I = 1,8
          DO 1 J = 1,2
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = WIDTH(IP)
      COORD(2,1) = WIDTH(IP)
      COORD(3,1) = WIDTH(IP)
      COORD(5,1) = WIDTH(IP+1)
      COORD(6,1) = WIDTH(IP+1)
      COORD(7,1) = WIDTH(IP+1)
      COORD(4,1) = .5D0* (COORD(3,1)+COORD(5,1))
      COORD(8,1) = .5D0* (COORD(7,1)+COORD(1,1))
      COORD(1,2) = DEPTH(IQ+1)
      COORD(8,2) = DEPTH(IQ+1)
      COORD(7,2) = DEPTH(IQ+1)
      COORD(3,2) = DEPTH(IQ)
      COORD(4,2) = DEPTH(IQ)
      COORD(5,2) = DEPTH(IQ)
      COORD(2,2) = .5D0* (COORD(1,2)+COORD(3,2))
      COORD(6,2) = .5D0* (COORD(5,2)+COORD(7,2))
      RETURN
 
      END
      SUBROUTINE GEO15Y(IP,IQ,NYE,WID,DEP,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 15-NODE TRIANGLES COUNTING IN Y-DIRECTION
C
      DOUBLE PRECISION FAC1
      DOUBLE PRECISION FAC2
      DOUBLE PRECISION COORD(ICOORD,*),WID(*),DEP(*)
      INTEGER G(*),NF(INF,*),NUM(15)
 
      IF (MOD(IQ,2).EQ.0) GO TO 1
      FAC1 = 4* (2*NYE+1)* (IP-1) + 2*IQ - 1
      NUM(1) = FAC1
      NUM(2) = FAC1 + 1
      NUM(3) = FAC1 + 2
      NUM(4) = FAC1 + 3
      NUM(5) = FAC1 + 4
      NUM(6) = FAC1 + 2*NYE + 4
      NUM(7) = FAC1 + 4*NYE + 4
      NUM(8) = FAC1 + 6*NYE + 4
      NUM(9) = FAC1 + 8*NYE + 4
      NUM(10) = FAC1 + 6*NYE + 3
      NUM(11) = FAC1 + 4*NYE + 2
      NUM(12) = FAC1 + 2*NYE + 1
      NUM(13) = FAC1 + 2*NYE + 2
      NUM(14) = FAC1 + 2*NYE + 3
      NUM(15) = FAC1 + 4*NYE + 3
      COORD(1,1) = WID(IP)
      COORD(1,2) = DEP((IQ+1)/2)
      COORD(5,1) = WID(IP)
      COORD(5,2) = DEP((IQ+3)/2)
      COORD(9,1) = WID(IP+1)
      COORD(9,2) = DEP((IQ+1)/2)
      GO TO 2
 
    1 FAC2 = 4* (2*NYE+1)* (IP-1) + 2*IQ + 8*NYE + 5
      NUM(1) = FAC2
      NUM(2) = FAC2 - 1
      NUM(3) = FAC2 - 2
      NUM(4) = FAC2 - 3
      NUM(5) = FAC2 - 4
      NUM(6) = FAC2 - 2*NYE - 4
      NUM(7) = FAC2 - 4*NYE - 4
      NUM(8) = FAC2 - 6*NYE - 4
      NUM(9) = FAC2 - 8*NYE - 4
      NUM(10) = FAC2 - 6*NYE - 3
      NUM(11) = FAC2 - 4*NYE - 2
      NUM(12) = FAC2 - 2*NYE - 1
      NUM(13) = FAC2 - 2*NYE - 2
      NUM(14) = FAC2 - 2*NYE - 3
      NUM(15) = FAC2 - 4*NYE - 3
      COORD(1,1) = WID(IP+1)
      COORD(1,2) = DEP((IQ+2)/2)
      COORD(5,1) = WID(IP+1)
      COORD(5,2) = DEP(IQ/2)
      COORD(9,1) = WID(IP)
      COORD(9,2) = DEP((IQ+2)/2)
    2 DO 3 I = 1,2
          COORD(3,I) = .5D0* (COORD(1,I)+COORD(5,I))
          COORD(7,I) = .5D0* (COORD(5,I)+COORD(9,I))
          COORD(11,I) = .5D0* (COORD(9,I)+COORD(1,I))
          COORD(2,I) = .5D0* (COORD(1,I)+COORD(3,I))
          COORD(4,I) = .5D0* (COORD(3,I)+COORD(5,I))
          COORD(6,I) = .5D0* (COORD(5,I)+COORD(7,I))
          COORD(8,I) = .5D0* (COORD(7,I)+COORD(9,I))
          COORD(10,I) = .5D0* (COORD(9,I)+COORD(11,I))
          COORD(12,I) = .5D0* (COORD(11,I)+COORD(1,I))
          COORD(15,I) = .5D0* (COORD(7,I)+COORD(11,I))
          COORD(14,I) = .5D0* (COORD(3,I)+COORD(7,I))
          COORD(13,I) = .5D0* (COORD(2,I)+COORD(15,I))
    3 CONTINUE
      INC = 0
      DO 4 I = 1,15
          DO 4 J = 1,2
              INC = INC + 1
    4 G(INC) = NF(NUM(I),J)
      RETURN
 
      END
      SUBROUTINE GEO4X1(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 4-NODE QUADS COUNTING IN X-DIRECTION
C      LAPLACE'S EQUATION   1-FREEDOM PER NODE
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER NUM(4),G(*),NF(INF,*)
 
      NUM(1) = IQ* (NXE+1) + IP
      NUM(2) = (IQ-1)* (NXE+1) + IP
      NUM(3) = NUM(2) + 1
      NUM(4) = NUM(1) + 1
      DO 1 I = 1,4
    1 G(I) = NF(NUM(I),1)
      COORD(1,1) = (IP-1)*AA
      COORD(2,1) = (IP-1)*AA
      COORD(3,1) = IP*AA
      COORD(4,1) = IP*AA
      COORD(1,2) = -IQ*BB
      COORD(2,2) = - (IQ-1)*BB
      COORD(3,2) = - (IQ-1)*BB
      COORD(4,2) = -IQ*BB
      RETURN
 
      END
      SUBROUTINE GEO83D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE BRICK ELEMENTS COUNTING X-Z PLANES IN Y-DIRECTION
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION CC
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(8)
 
      NUM(1) = (IQ-1)* (NXE+1)* (NZE+1) + IS* (NXE+1) + IP
      NUM(2) = NUM(1) - NXE - 1
      NUM(3) = NUM(2) + 1
      NUM(4) = NUM(1) + 1
      NUM(5) = NUM(1) + (NXE+1)* (NZE+1)
      NUM(6) = NUM(5) - NXE - 1
      NUM(7) = NUM(6) + 1
      NUM(8) = NUM(5) + 1
      INC = 0
      DO 1 I = 1,8
          DO 1 J = 1,3
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = (IP-1)*AA
      COORD(2,1) = (IP-1)*AA
      COORD(5,1) = (IP-1)*AA
      COORD(6,1) = (IP-1)*AA
      COORD(3,1) = IP*AA
      COORD(4,1) = IP*AA
      COORD(7,1) = IP*AA
      COORD(8,1) = IP*AA
      COORD(1,2) = (IQ-1)*BB
      COORD(2,2) = (IQ-1)*BB
      COORD(3,2) = (IQ-1)*BB
      COORD(4,2) = (IQ-1)*BB
      COORD(5,2) = IQ*BB
      COORD(6,2) = IQ*BB
      COORD(7,2) = IQ*BB
      COORD(8,2) = IQ*BB
      COORD(1,3) = -IS*CC
      COORD(4,3) = -IS*CC
      COORD(5,3) = -IS*CC
      COORD(8,3) = -IS*CC
      COORD(2,3) = - (IS-1)*CC
      COORD(3,3) = - (IS-1)*CC
      COORD(6,3) = - (IS-1)*CC
      COORD(7,3) = - (IS-1)*CC
      RETURN
 
      END
      SUBROUTINE GEV4X3(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE NODAL COORDINATES AND STEERING
C      VECTOR FOR A VARIABLE MESH OF 4-NODE QUADRILATERAL ELEMENTS
C      NUMBERING IN THE X-DIRECTION
C      (U,V,P  BIOT CONSOLIDATION)
C
      DOUBLE PRECISION COORD(ICOORD,*),WIDTH(*),DEPTH(*)
      INTEGER NUM(4),G(*),NF(INF,*)
 
      NUM(1) = IQ* (NXE+1) + IP
      NUM(2) = (IQ-1)* (NXE+1) + IP
      NUM(3) = NUM(2) + 1
      NUM(4) = NUM(1) + 1
      NINC = 0
      DO 1 I = 1,4
          DO 1 J = 1,2
              NINC = NINC + 1
    1 G(NINC) = NF(NUM(I),J)
      DO 2 I = 1,4
          NINC = NINC + 1
    2 G(NINC) = NF(NUM(I),3)
      COORD(1,1) = WIDTH(IP)
      COORD(2,1) = WIDTH(IP)
      COORD(3,1) = WIDTH(IP+1)
      COORD(4,1) = WIDTH(IP+1)
      COORD(2,2) = DEPTH(IQ)
      COORD(3,2) = DEPTH(IQ)
      COORD(4,2) = DEPTH(IQ+1)
      COORD(1,2) = DEPTH(IQ+1)
      RETURN
 
      END
      SUBROUTINE GE203D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE STEERING VECTOR AND COORDINATES
C      FOR 20-NODE BRICK ELEMENTS COUNTING X-Z PLANES IN Y-DIRECTION
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION BB
      DOUBLE PRECISION CC
      DOUBLE PRECISION FAC1
      DOUBLE PRECISION FAC2
      DOUBLE PRECISION COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(20)
 
      FAC1 = ((2*NXE+1)* (NZE+1)+ (2*NZE+1)* (NXE+1))* (IQ-1)
      FAC2 = ((2*NXE+1)* (NZE+1)+ (2*NZE+1)* (NXE+1))*IQ
      NUM(1) = FAC1 + (3*NXE+2)*IS + 2*IP - 1
      NUM(2) = FAC1 + (3*NXE+2)*IS - NXE + IP - 1
      NUM(3) = NUM(1) - 3*NXE - 2
      NUM(4) = NUM(3) + 1
      NUM(5) = NUM(4) + 1
      NUM(6) = NUM(2) + 1
      NUM(7) = NUM(1) + 2
      NUM(8) = NUM(1) + 1
      NUM(9) = FAC2 - (NXE+1)* (NZE+1) + (NXE+1)*IS + IP
      NUM(10) = NUM(9) - NXE - 1
      NUM(11) = NUM(10) + 1
      NUM(12) = NUM(9) + 1
      NUM(13) = FAC2 + (3*NXE+2)*IS + 2*IP - 1
      NUM(14) = FAC2 + (3*NXE+2)*IS - NXE + IP - 1
      NUM(15) = NUM(13) - 3*NXE - 2
      NUM(16) = NUM(15) + 1
      NUM(17) = NUM(16) + 1
      NUM(18) = NUM(14) + 1
      NUM(19) = NUM(13) + 2
      NUM(20) = NUM(13) + 1
      INC = 0
      DO 1 I = 1,20
          DO 1 J = 1,3
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      COORD(1,1) = (IP-1)*AA
      COORD(2,1) = (IP-1)*AA
      COORD(3,1) = (IP-1)*AA
      COORD(9,1) = (IP-1)*AA
      COORD(10,1) = (IP-1)*AA
      COORD(13,1) = (IP-1)*AA
      COORD(14,1) = (IP-1)*AA
      COORD(15,1) = (IP-1)*AA
      COORD(5,1) = IP*AA
      COORD(6,1) = IP*AA
      COORD(7,1) = IP*AA
      COORD(11,1) = IP*AA
      COORD(12,1) = IP*AA
      COORD(17,1) = IP*AA
      COORD(18,1) = IP*AA
      COORD(19,1) = IP*AA
      COORD(4,1) = .5D0* (COORD(3,1)+COORD(5,1))
      COORD(8,1) = .5D0* (COORD(1,1)+COORD(7,1))
      COORD(16,1) = .5D0* (COORD(15,1)+COORD(17,1))
      COORD(20,1) = .5D0* (COORD(13,1)+COORD(19,1))
      COORD(1,2) = (IQ-1)*BB
      COORD(2,2) = (IQ-1)*BB
      COORD(3,2) = (IQ-1)*BB
      COORD(4,2) = (IQ-1)*BB
      COORD(5,2) = (IQ-1)*BB
      COORD(6,2) = (IQ-1)*BB
      COORD(7,2) = (IQ-1)*BB
      COORD(8,2) = (IQ-1)*BB
      COORD(13,2) = IQ*BB
      COORD(14,2) = IQ*BB
      COORD(15,2) = IQ*BB
      COORD(16,2) = IQ*BB
      COORD(17,2) = IQ*BB
      COORD(18,2) = IQ*BB
      COORD(19,2) = IQ*BB
      COORD(20,2) = IQ*BB
      COORD(9,2) = .5D0* (COORD(1,2)+COORD(13,2))
      COORD(10,2) = .5D0* (COORD(3,2)+COORD(15,2))
      COORD(11,2) = .5D0* (COORD(5,2)+COORD(17,2))
      COORD(12,2) = .5D0* (COORD(7,2)+COORD(19,2))
      COORD(1,3) = -IS*CC
      COORD(7,3) = -IS*CC
      COORD(8,3) = -IS*CC
      COORD(9,3) = -IS*CC
      COORD(12,3) = -IS*CC
      COORD(13,3) = -IS*CC
      COORD(19,3) = -IS*CC
      COORD(20,3) = -IS*CC
      COORD(3,3) = - (IS-1)*CC
      COORD(4,3) = - (IS-1)*CC
      COORD(5,3) = - (IS-1)*CC
      COORD(10,3) = - (IS-1)*CC
      COORD(11,3) = - (IS-1)*CC
      COORD(15,3) = - (IS-1)*CC
      COORD(16,3) = - (IS-1)*CC
      COORD(17,3) = - (IS-1)*CC
      COORD(2,3) = .5D0* (COORD(1,3)+COORD(3,3))
      COORD(6,3) = .5D0* (COORD(5,3)+COORD(7,3))
      COORD(14,3) = .5D0* (COORD(13,3)+COORD(15,3))
      COORD(18,3) = .5D0* (COORD(17,3)+COORD(19,3))
      RETURN
 
      END
      SUBROUTINE SLOGEO(IP,IQ,NYE,TOP,BOT,DEPTH,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE QUADRILATERALS IN A 'SLOPE' GEOMETRY
C      (NUMBERING IN THE Y-DIRECTION)
C
      DOUBLE PRECISION FAC1
      DOUBLE PRECISION FAC2
      DOUBLE PRECISION TOP(*),BOT(*),COORD(ICOORD,*),DEPTH(*)
      INTEGER G(*),NF(INF,*),NUM(8)
 
      NUM(1) = (IP-1)* (3*NYE+2) + 2*IQ + 1
      NUM(2) = NUM(1) - 1
      NUM(3) = NUM(1) - 2
      NUM(4) = (IP-1)* (3*NYE+2) + 2*NYE + IQ + 1
      NUM(5) = IP* (3*NYE+2) + 2*IQ - 1
      NUM(6) = NUM(5) + 1
      NUM(7) = NUM(5) + 2
      NUM(8) = NUM(4) + 1
      INC = 0
      DO 1 I = 1,8
          DO 1 J = 1,2
              INC = INC + 1
    1 G(INC) = NF(NUM(I),J)
      FAC1 = (BOT(IP)-TOP(IP))/ (DEPTH(NYE+1)-DEPTH(1))
      FAC2 = (BOT(IP+1)-TOP(IP+1))/ (DEPTH(NYE+1)-DEPTH(1))
      COORD(1,1) = TOP(IP) + (DEPTH(IQ+1)-DEPTH(1))*FAC1
      COORD(3,1) = TOP(IP) + (DEPTH(IQ)-DEPTH(1))*FAC1
      COORD(5,1) = TOP(IP+1) + (DEPTH(IQ)-DEPTH(1))*FAC2
      COORD(7,1) = TOP(IP+1) + (DEPTH(IQ+1)-DEPTH(1))*FAC2
      COORD(2,1) = .5D0* (COORD(1,1)+COORD(3,1))
      COORD(6,1) = .5D0* (COORD(5,1)+COORD(7,1))
      COORD(4,1) = .5D0* (COORD(3,1)+COORD(5,1))
      COORD(8,1) = .5D0* (COORD(7,1)+COORD(1,1))
      COORD(1,2) = DEPTH(IQ+1)
      COORD(8,2) = DEPTH(IQ+1)
      COORD(7,2) = DEPTH(IQ+1)
      COORD(3,2) = DEPTH(IQ)
      COORD(4,2) = DEPTH(IQ)
      COORD(5,2) = DEPTH(IQ)
      COORD(2,2) = .5D0* (COORD(1,2)+COORD(3,2))
      COORD(6,2) = .5D0* (COORD(5,2)+COORD(7,2))
      RETURN
 
      END
      SUBROUTINE WELGEO(IP,IQ,NXE,NYE,WIDTH,SURF,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 4-NODE QUADS NUMBERING IN THE X-DIRECTION
C      LAPLACE'S EQUATION,VARIABLE MESH, 1-FREEDOM PER NODE
C
      DOUBLE PRECISION BB
      DOUBLE PRECISION BVAR
      DOUBLE PRECISION COORD(ICOORD,*),WIDTH(*),SURF(*)
      INTEGER NUM(4),G(*),NF(INF,*)
 
      NUM(1) = IQ* (NXE+1) + IP
      NUM(2) = (IQ-1)* (NXE+1) + IP
      NUM(3) = NUM(2) + 1
      NUM(4) = NUM(1) + 1
      DO 1 I = 1,4
    1 G(I) = NF(NUM(I),1)
      BB = SURF(IP)/NYE
      BVAR = SURF(IP+1)/NYE
      COORD(1,1) = WIDTH(IP)
      COORD(2,1) = WIDTH(IP)
      COORD(3,1) = WIDTH(IP+1)
      COORD(4,1) = WIDTH(IP+1)
      COORD(1,2) = (NYE-IQ)*BB
      COORD(2,2) = (NYE-IQ+1)*BB
      COORD(3,2) = (NYE-IQ+1)*BVAR
      COORD(4,2) = (NYE-IQ)*BVAR
      RETURN
 
      END
