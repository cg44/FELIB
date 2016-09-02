C  ********************************************************************
C  *                                                                  *
C  *                       FUNCTION DBINT2                            *
C  *                                                                  *
C  ********************************************************************
C  DOUBLE PRECISION VERSION 1.0
C  WRITTEN BY G.A. FENTON, TUNS, 1991
C
C  PURPOSE  TO COMPUTE THE QUADRUPLE INTEGRAL OF A FUNCTION INVOLVING THE
C           DIFFERENCE OF THE COORDINATES, IE F2D(X_1-X_2, Y_1-Y_2)
C
C  RETURNS THE QUADRUPLE INTEGRAL OF A USER SUPPLIED FUNCTION VIA GAUSSIAN
C  QUADRATURE. THE INTEGRAL IS ASSUMED TO HAVE THE FORM
C
C           B    D    F    H
C         INT  INT  INT  INT F2D( X_1 - X_2, Y_1 - Y_2 ) DX_1 DX_2 DY_1 DY_2
C           A    C    E    G
C
C  ARGUMENTS ARE DESCRIBED AS FOLLOWS;
C
C    F2D   EXTERNALLY SUPPLIED SINGLE PRECISION FUNCTION NAME WHICH RETURNS
C          THE FUNCTION VALUES (THIS IS THE FUNCTION TO BE INTEGRATED).
C          F2D(X,Y) IS ASSUMED TO NOT BLOW UP AT ANY OF THE GAUSS POINTS.
C
C  B1,B2   REAL VECTORS OF LENGTH AT LEAST 4 CONTAINING RESPECTIVELY, (INPUT)
C
C             B1(1) = LOWER INTEGRATION BOUND ON X_1 (G)
C             B1(2) = UPPER INTEGRATION BOUND ON X_1 (H)
C             B1(3) = LOWER INTEGRATION BOUND ON Y_1 (C)
C             B1(4) = UPPER INTEGRATION BOUND ON Y_1 (D)
C
C             B2(1) = LOWER INTEGRATION BOUND ON X_2 (E)
C             B2(2) = UPPER INTEGRATION BOUND ON X_2 (F)
C             B2(3) = LOWER INTEGRATION BOUND ON Y_2 (A)
C             B2(4) = UPPER INTEGRATION BOUND ON Y_2 (B)
C
C    NG    NUMBER OF GAUSS POINTS DESIRED (SELECT FROM 1,2,3,4,5,6,7,8,9,10,
C          16, AND 20). (INPUT)
C
C  THIS ROUTINE CALLS `DBNT2A'.
C ========================================================================
 
      REAL FUNCTION DBINT2( F2D, B1, B2, NG )
      DIMENSION B1(*), B2(*)
      DIMENSION W(54), Z(54), W16(8), Z16(8), W20(10), Z20(10)
      INTEGER KEY(10)
      EXTERNAL F2D
      COMMON/GDBLCK/ KEY, ZERO, HALF, TWO, W, W16, W20, Z, Z16, Z20
 
C					INITIALIZE INTEGRATION PARAMETERS
      A1 = B1(3) - B2(4)
      B3 = B1(4) - B2(3)
      IF( (B1(4)-B1(3)) .GE. (B2(4)-B2(3)) ) THEN
         A2 = B1(3) - B2(3)
         A3 = B1(4) - B2(4)
         C2 = B2(4) - B2(3)
      ELSE
         A2 = B1(4) - B2(4)
         A3 = B1(3) - B2(3)
         C2 = B1(4) - B1(3)
      ENDIF
 
      R1 = HALF*(A2 - A1)
      R2 = HALF*(A3 - A2)
      R3 = HALF*(B3 - A3)
      S1 = HALF*(A1 + A2)
      S2 = HALF*(A2 + A3)
      S3 = HALF*(A3 + B3)
C					INITIALIZE DBNT2A
      DBINT2 = DBNT2A( F2D, ZERO, B1, B2, 0 )
 
      D1 = ZERO
      D2 = ZERO
      D3 = ZERO
C					COMPUTE THE INTEGRAL
      IF( NG .EQ. 1 ) THEN
         D1 = TWO*(S1 - A1)*DBNT2A( F2D, S1, B1, B2, NG )
         D2 = TWO  *  C2  * DBNT2A( F2D, S2, B1, B2, NG )
         D3 = TWO*(B3 - S3)*DBNT2A( F2D, S3, B1, B2, NG )
      ELSEIF( NG .LT. 11 ) THEN
         J1   = KEY(NG-1)
         J2   = KEY(NG) - 1
         DO 10 I = J1, J2
            T1 = Z(I)*R1 + S1
            T2 = Z(I)*R2 + S2
            T3 = Z(I)*R3 + S3
            D1 = D1 + W(I)*(T1 - A1)*DBNT2A(F2D, T1, B1, B2, NG )
            D2 = D2 + W(I)  *  C2  * DBNT2A(F2D, T2, B1, B2, NG )
            D3 = D3 + W(I)*(B3 - T3)*DBNT2A(F2D, T3, B1, B2, NG )
  10     CONTINUE
      ELSEIF( NG .LT. 17 ) THEN
         DO 20 I = 1, 8
            T1 = Z16(I)*R1
            T2 = Z16(I)*R2
            T3 = Z16(I)*R3
            TA = S1 + T1
            TC = S3 + T3
            T1 = S1 - T1
            T3 = S3 - T3
            D1 = D1 + W16(I)*((TA-A1)*DBNT2A(F2D, TA,    B1, B2, NG)
     >                      + (T1-A1)*DBNT2A(F2D, T1,    B1, B2, NG))
            D2 = D2 + W16(I) * C2 *  (DBNT2A(F2D, S2+T2, B1, B2, NG)
     >                      +         DBNT2A(F2D, S2-T2, B1, B2, NG))
            D3 = D3 + W16(I)*((B3-TC)*DBNT2A(F2D, TC,    B1, B2, NG)
     >                      + (B3-T3)*DBNT2A(F2D, T3,    B1, B2, NG))
  20     CONTINUE
      ELSE
         DO 30 I = 1, 10
            T1 = Z20(I)*R1
            T2 = Z20(I)*R2
            T3 = Z20(I)*R3
            TA = S1 + T1
            TC = S3 + T3
            T1 = S1 - T1
            T3 = S3 - T3
            D1 = D1 + W20(I)*((TA-A1)*DBNT2A(F2D, TA,    B1, B2, NG)
     >                      + (T1-A1)*DBNT2A(F2D, T1,    B1, B2, NG))
            D2 = D2 + W20(I) * C2 *  (DBNT2A(F2D, S2+T2, B1, B2, NG)
     >                      +         DBNT2A(F2D, S2-T2, B1, B2, NG))
            D3 = D3 + W20(I)*((B3-TC)*DBNT2A(F2D, TC,    B1, B2, NG)
     >                      + (B3-T3)*DBNT2A(F2D, T3,    B1, B2, NG))
  30     CONTINUE
      ENDIF
 
      DBINT2 = R1*D1 + R2*D2 + R3*D3
 
      RETURN
      END
C  ********************************************************************
C  *                                                                  *
C  *                       FUNCTION DBNT2A                            *
C  *                                                                  *
C  ********************************************************************
C  DOUBLE PRECISION VERSION 1.0
C  WRITTEN BY G.A. FENTON, TUNS, 1991
C
C  PURPOSE  TO COMPUTE THE DOUBLE INTEGRAL OF A FUNCTION INVOLVING THE
C           DIFFERENCE OF THE COORDINATES, IE F2D(X_1-X_2, Y). CALLED
C           BY DBINT2.
C
C  RETURNS THE DOUBLE INTEGRAL OF A USER SUPPLIED FUNCTION VIA GAUSSIAN
C  QUADRATURE. THE INTEGRAL IS ASSUMED TO HAVE THE FORM
C
C           B    D
C         INT  INT  F2D( X_1 - X_2, Y ) DX_1 DX_2
C           A    C
C
C  FOR A GIVEN FIXED Y. ARGUMENTS ARE DESCRIBED AS FOLLOWS;
C
C    F2D   EXTERNALLY SUPPLIED SINGLE PRECISION FUNCTION NAME WHICH RETURNS
C          THE FUNCTION VALUES (THIS IS THE FUNCTION TO BE INTEGRATED).
C          F2D(X,Y) IS ASSUMED TO NOT BLOW UP AT ANY OF THE GAUSS POINTS.
C
C    Y     FIXED VALUE AT WHICH THE SECOND COORDINATE OF F2D(X,Y) IS TO BE
C          HELD DURING INTEGRATION. (INPUT)
C
C  B1,B2   REAL VECTORS OF LENGTH AT LEAST 2 CONTAINING RESPECTIVELY,
C
C             B1(1) = LOWER INTEGRATION BOUND ON X_1 (C)
C             B1(2) = UPPER INTEGRATION BOUND ON X_1 (D)
C
C             B2(1) = LOWER INTEGRATION BOUND ON X_2 (A)
C             B2(2) = UPPER INTEGRATION BOUND ON X_2 (B)
C
C    NG    NUMBER OF GAUSS POINTS DESIRED (SELECT FROM 1,2,3,4,5,6,7,8,9,10,
C          16, AND 20). IF NG = 0 THEN ONLY THE INTEGRATION PARAMETERS ARE
C          INITIALIZED AND 0 IS RETURNED. TO ACTUALLY PERFORM THE INTEGRATION,
C          THIS ROUTINE MUST BE CALLED ONCE WITH NG = 0 AND THEN AGAIN WITH
C          THE DESIRED NUMBER OF GAUSS POINTS. (INPUT)
C
C ========================================================================
 
      REAL FUNCTION DBNT2A( F2D, Y, B1, B2, NG )
      DIMENSION B1(*), B2(*)
      DIMENSION W(54), Z(54), W16(8), Z16(8), W20(10), Z20(10)
      INTEGER KEY(10)
      SAVE A1, B3, C2, S1, S2, S3, R1, R2, R3
      EXTERNAL F2D
      COMMON/GDBLCK/ KEY, ZERO, HALF, TWO, W, W16, W20, Z, Z16, Z20
 
C					INITIALIZE IF NECESSARY
      IF( NG .LE. 0 ) THEN
         A1 = B1(1) - B2(2)
         B3 = B1(2) - B2(1)
         IF( (B1(2)-B1(1)) .GE. (B2(2)-B2(1)) ) THEN
            A2 = B1(1) - B2(1)
            A3 = B1(2) - B2(2)
            C2 = B2(2) - B2(1)
         ELSE
            A2 = B1(2) - B2(2)
            A3 = B1(1) - B2(1)
            C2 = B1(2) - B1(1)
         ENDIF
 
         R1 = HALF*(A2 - A1)
         R2 = HALF*(A3 - A2)
         R3 = HALF*(B3 - A3)
         S1 = HALF*(A2 + A1)
         S2 = HALF*(A3 + A2)
         S3 = HALF*(B3 + A3)
         DBNT2A = ZERO
         RETURN
      ENDIF
 
      D1 = ZERO
      D2 = ZERO
      D3 = ZERO
      IF( NG .EQ. 1 ) THEN
         D1 = TWO*(S1 - A1)*F2D(S1,Y)
         D2 = TWO  *  C2  * F2D(S2,Y)
         D3 = TWO*(B3 - S3)*F2D(S3,Y)
      ELSEIF( NG .LT. 11 ) THEN
         J1   = KEY(NG-1)
         J2   = KEY(NG) - 1
         DO 10 I = J1, J2
            T1 = S1 + Z(I)*R1
            T2 = S2 + Z(I)*R2
            T3 = S3 + Z(I)*R3
            D1 = D1 + W(I)*(T1 - A1)*F2D(T1,Y)
            D2 = D2 + W(I)  *  C2  * F2D(T2,Y)
            D3 = D3 + W(I)*(B3 - T3)*F2D(T3,Y)
  10     CONTINUE
      ELSEIF( NG .LT. 17 ) THEN
         DO 20 I = 1, 8
            T1 = Z16(I)*R1
            T2 = Z16(I)*R2
            T3 = Z16(I)*R3
            TA = S1 + T1
            TC = S3 + T3
            T1 = S1 - T1
            T3 = S3 - T3
            D1 = D1 + W16(I)*((TA-A1)*F2D(TA,Y) + (T1-A1)*F2D(T1,Y))
            D2 = D2 + W16(I)*C2*(     F2D(S2+T2,Y)   +    F2D(S2-T2,Y))
            D3 = D3 + W16(I)*((B3-TC)*F2D(TC,Y) + (B3-T3)*F2D(T3,Y))
  20     CONTINUE
      ELSE
         DO 30 I = 1, 10
            T1 = Z20(I)*R1
            T2 = Z20(I)*R2
            T3 = Z20(I)*R3
            TA = S1 + T1
            TC = S3 + T3
            T1 = S1 - T1
            T3 = S3 - T3
            D1 = D1 + W20(I)*((TA-A1)*F2D(TA,Y) + (T1-A1)*F2D(T1,Y))
            D2 = D2 + W20(I)*C2*(     F2D(S2+T2,Y)   +    F2D(S2-T2,Y))
            D3 = D3 + W20(I)*((B3-TC)*F2D(TC,Y) + (B3-T3)*F2D(T3,Y))
  30     CONTINUE
      ENDIF
 
      DBNT2A = R1*D1 + R2*D2 + R3*D3
 
      RETURN
      END
C   *****************************************************************
C   *                                                               *
C   *                      FUNCTION DCHOLL                          *
C   *                                                               *
C   *****************************************************************
C   DOUBLE PRECISION VERSION 2.0, MATRIX STORAGE
C   WRITTEN BY GORDON A. FENTON, TUNS, FEB. 1991
C
C   PURPOSE TO COMPUTE THE CHOLESKY DECOMPOSITION L ACCORDING TO L*L' = A
C           AND RETURN THE DETERMINANT OF A
C
C   COMPUTES THE CHOLESKY DECOMPOSITION L OF A SYMMETRIC POSITIVE DEFINITE
C   REAL MATRIX A = L*L' WHERE THE PRIME INDICATES TRANSPOSE. L IS A REAL
C   NON-SINGULAR LOWER-TRIANGULAR MATRIX. THIS METHOD INVOLVES N DIVISIONS,
C   N SQUARE-ROOTS, AND (0.16*N**3 + 0.5*N**2) MULTIPLICATIONS. IF A IS
C   DISCOVERED TO BE NON-POSITIVE DEFINITE OR SINGULAR, THE DETERMINANT IS
C   SET TO 0 AND THE FUNCTION RETURNED.
C   THE VARIABLE DESCRIPTION IS AS FOLLOWS;
C
C    A    INPUT POSITIVE DEFINITE SYMMETRIC MATRIX WITH NECESSARY VALUES
C         STORED IN THE LOWER TRIANGLE (AT LEAST). ON OUTPUT, A WILL CONTAIN
C         THE CHOLESKY DECOMPOSITION IN ITS LOWER TRIANGLE (INCLUDING THE
C         DIAGONAL). THE UPPER TRIANGLE OF A IS NOT TOUCHED.
C
C    IA   INPUT INTEGER CONTAINING THE COLUMN DIMENSION OF A (AND THUS L) AS
C         SPECIFIED IN THE CALLING ROUTINE
C
C    N    INPUT INTEGER CONTAINING THE ORDER OF MATRIX A.
C
C ===========================================================================
      REAL FUNCTION DCHOLL( A, IA, N )
      DIMENSION A(IA,1)
      DATA ZERO/0.0/, ONE/1.0/
 
      DET    = ONE
      DCHOLL = ZERO
      DO 20 K = 1, N
C						CHECK DIAGONAL ELEMENT
         IF( A(K,K) .LE. ZERO ) RETURN
C						FIND L AND DETERMINANT
         DET    = DET*A(K,K)
         A(K,K) = SQRT( A(K,K) )
         DO 10 J = K+1, N
            A(J,K) = A(J,K)/A(K,K)
  10     CONTINUE
 
         DO 20 J = K+1, N
         DO 20 I = J, N
            A(I,J) = A(I,J) - A(I,K)*A(J,K)
  20  CONTINUE
 
      DCHOLL = DET
      RETURN
      END
C  ***********************************************************************
C  *                                                                     *
C  *                          BLOCK DATA DGQBLK                          *
C  *                                                                     *
C  ***********************************************************************
C  DOUBLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  TO SET UP THE GAUSSIAN INTEGRATION WEIGHTS AND POINTS DATA FOR
C           USE IN GAUSSIAN INTEGRATION ROUTINES.
C
C  THESE ARE THE WEIGHTS AND INTEGRATION POINTS FOR GAUSSIAN QUADRATURE FOR
C  1,2,3,4,5,6,7,8,9,10,16, AND 20 GAUSS POINTS.
C------------------------------------------------------------------------------
      BLOCK DATA DGQBLK
      DIMENSION W(54), Z(54), W16(8), Z16(8), W20(10), Z20(10)
      INTEGER KEY(10)
      COMMON/GDBLCK/ KEY, ZERO, HALF, TWO, W, W16, W20, Z, Z16, Z20
 
      DATA KEY/1, 3, 6, 10, 15, 21, 28, 36, 45, 55/
      DATA ZERO/0.0/, HALF/0.50/, TWO/2.0/
 
      DATA W/1.00000000000000,1.000000000000000,.555555555555556,
     >        .888888888888889,.555555555555556,.347854845137454,
     >        .652145154862546,.652145154862546,.347854845137454,
     >        .236926885056189,.478628670499366,.568888888888889,
     >        .478628670499366,.236926885056189,.171324492379170,
     >        .360761573048139,.467913934572691,.467913934572691,
     >        .360761573048139,.171324492379170,.129484966168870,
     >        .279705391489277,.381830050505119,.417959183673469,
     >        .381830050505119,.279705391489277,.129484966168870,
     >        .101228536290376,.222381034453374,.313706645877887,
     >        .362683783378362,.362683783378362,.313706645877887,
     >        .222381034453374,.101228536290376,.081274388361574,
     >        .180648160694857,.260610696402935,.312347077040003,
     >        .330239355001260,.312347077040003,.260610696402935,
     >        .180648160694857,.081274388361574,.066671344308688,
     >        .149451349150581,.219086362515982,.269266719309996,
     >        .295524224714753,.295524224714753,.269266719309996,
     >        .219086362515982,.149451349150581,.066671344308688/
      DATA W16/0.027152459411754094852, 0.062253523938647892863,
     >         0.095158511682492784810, 0.124628971255533872052,
     >         0.149595988816576732081, 0.169156519395002538189,
     >         0.182603415044923588867, 0.189450610455068496285/
      DATA W20/0.017614007139152118312, 0.040601429800386941331,
     >         0.062672048334109063570, 0.083276741576704748725,
     >         0.101930119817240435037, 0.118194531961518417312,
     >         0.131688638449176626898, 0.142096109318382051329,
     >         0.149172986472603746788, 0.152753387130725850698/
 
      DATA Z/-.577350269189626,.577350269189626,-.774596669241483,
     >       .000000000000000, .774596669241483,-.861136311594053,
     >      -.339981043584856, .339981043584856, .861136311594053,
     >      -.906179845938664,-.538469310105683, .000000000000000,
     >       .538469310105683, .906179845938664,-.932469514203152,
     >      -.661209386466265,-.238619186083197, .238619186083197,
     >       .661209386466265, .932469514203152,-.949107912342759,
     >      -.741531185599394,-.405845151377397, .000000000000000,
     >       .405845151377397, .741531185599394, .949107912342759,
     >      -.960289856497536,-.796666477413627,-.525532409916329,
     >      -.183434642495650, .183434642495650, .525532409916329,
     >       .796666477413627, .960289856497536,-.968160239507626,
     >      -.836031107326636,-.613371432700590,-.324253423403809,
     >       .000000000000000, .324253423403809, .613371432700590,
     >       .836031107326636, .968160239507626,-.973906528517172,
     >      -.865063366688985,-.679409568299024,-.433395394129247,
     >      -.148874338981632, .148874338981632, .433395394129247,
     >       .679409568299024, .865063366688985, .973906528517172/
      DATA Z16/0.989400934991649932596, 0.944575023073232576078,
     >         0.865631202387831743880, 0.755404408355003033895,
     >         0.617876244402643748447, 0.458016777657227386342,
     >         0.281603550779258913230, 0.095012509837637440185/
      DATA Z20/0.993128599185094924786, 0.963971927277913791268,
     >         0.912234428251325905868, 0.839116971822218823395,
     >         0.746331906460150792614, 0.636053680726515025453,
     >         0.510867001950827098004, 0.373706088715419560673,
     >         0.227785851141645078080, 0.076526521133497333755/
 
      END
C  *********************************************************************
C  *                                                                   *
C  *                         SUBROUTINE DPRTMT                         *
C  *                                                                   *
C  *********************************************************************
C  DOUBLE PRECISION VERSION 2.0
C  WRITTEN BY GORDON A. FENTON, 1989, LAST UPDATE MARCH 1991
C
C  PURPOSE   TO PRINT A MATRIX
C
C  THIS ROUTINE PRINTS AN N X M MATRIX TO UNIT IOUT.
C  ARGUMENTS ARE AS FOLLOWS;
C
C IOUT   FORTRAN UNIT NUMBER (ASSUMED ALREADY OPENED). (INPUT)
C
C    U   REAL ARRAY CONTAINING THE MATRIX TO BE PRINTED. (INPUT)
C
C   IU   LEADING DIMENSION OF THE ARRAY U AS SPECIFIED IN THE DIMENSION
C        STATEMENT OF THE CALLING ROUTINE. (INPUT)
C
C    N   COLUMN DIMENSION OF U (NUMBER OF ROWS). (INPUT)
C
C    M   ROW DIMENSION OF U (NUMBER OF COLUMNS). (INPUT)
C
C TITLE  CHARACTER STRING GIVING THE TITLE TO BE PRINTED. (INPUT)
C
C THE PARAMETER NCOL IS THE NUMBER OF COLUMNS WITHIN WHICH TO FORMAT THE
C ROWS OF THE MATRIX.
C--------------------------------------------------------------------------
      SUBROUTINE DPRTMT( IOUT, U, IU, N, M, TITLE )
      PARAMETER (NCOL = 80)
      REAL U(IU,1)
      CHARACTER*(*) TITLE
      CHARACTER*16 FMT
      DATA ZERO/0.0/, ONE/1.0/
 
   1  FORMAT(A)
   2  FORMAT()
   3  FORMAT('(',I2,'e',I2,'.',I1,')')
   4  FORMAT('(',I2,'f',I2,'.',I2,')')
 
      WRITE(IOUT,2)
      WRITE(IOUT,1) TITLE
C					FIND FIELD WIDTH FOR 80 COLUMN PRINTER
      NFLD = NCOL/M
      IF( NFLD .LT. 1 ) THEN
         WRITE(IOUT,1) '    Matrix to big to print ...'
         WRITE(IOUT,2)
         RETURN
      ENDIF
C					CHECK DATA TO FIND RANGE
      UMAX = ABS(U(1,1))
      UMIN = ABS(U(1,1))
      ISYN = 0
      DO 10 J = 1, M
      DO 10 I = 1, N
         U0 = ABS(U(I,J))
         IF( U0 .GT. UMAX ) UMAX = U0
         IF( U0 .LT. UMIN ) UMIN = U0
         IF( U(I,J) .LT. 0.0 ) ISYN = 1
  10  CONTINUE
C					CAN WE USE F FORMAT?
      IF( UMAX .GE. ONE ) THEN
         NRT = 2 + ISYN + INT( ALOG10(UMAX) )
      ELSE
         NRT = 1 + ISYN
      ENDIF
      IF( UMIN .LT. ONE ) THEN
         IF( UMIN .EQ. ZERO ) THEN
            NLT = 0
         ELSE
            NLT = 2 - INT( ALOG10(UMIN) )
         ENDIF
      ELSE
         NLT = 0
      ENDIF
      NREQ = NRT + NLT
      IF( NREQ .GT. NFLD ) THEN
C						NO, USE E FORMAT
         NDEC = MIN0(7, NFLD - 8)
         NFLD = MIN0( NFLD, NDEC + 10 )
         IF( NDEC .LT. 1 ) THEN
            WRITE(IOUT,1) '    Matrix to big to print ...'
            WRITE(IOUT,2)
            RETURN
         ENDIF
C						HARDWIRE SOME LIMITS
         IF( NFLD .GT. 20 ) NFLD = 20
         IF( NDEC .GT. 9  ) NDEC =  9
         WRITE(FMT,3) M,NFLD,NDEC
      ELSE
C						YES, USE F FORMAT
         NSIG = MAX0(1, 9 - (NRT + NLT - 2 - ISYN))
         NDEC = NLT + MIN0( NSIG, NFLD - NREQ )
         NFLD = MIN0( NFLD, NRT + NDEC + 2 )
         IF( NFLD .GT. 30 ) NFLD = 30
         IF( NDEC .LT. 0 ) THEN
            NDEC = 0
         ELSE IF( NDEC .GT. NFLD ) THEN
            NDEC = NFLD - 2
         ENDIF
         WRITE(FMT,4) M, NFLD, NDEC
      ENDIF
 
      DO 20 I = 1, N
         WRITE(IOUT,FMT) ( U(I,J), J = 1, M )
  20  CONTINUE
 
      RETURN
      END
C  ***********************************************************************
C  *                                                                     *
C  *                          SUBROUTINE DSLVCH                          *
C  *                                                                     *
C  ***********************************************************************
C  DOUBLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  TO COMPUTE THE SOLUTION TO THE MATRIX EQUATION [A]{X} = {B}
C           GIVEN THE CHOLESKY LL' DECOMPOSITION OF A AS PRODUCED BY DCHOLL.
C
C  THIS ROUTINE ACCEPTS AS INPUT THE CHOLESKY LL' DECOMPOSITION OF [A] STORED
C  IN PLACE (IE IN [A], SEE BELOW) AND THE RIGHT-HAND-SIDE VECTOR {B}, AND
C  COMPUTES THE SOLUTION {X} OF [A]{X} = {B}. THE SOLUTION IS RETURNED IN THE
C  VECTOR {B}. ARGUMENTS TO THE ROUTINE ARE AS FOLLOWS;
C
C     A   REAL ARRAY OF SIZE AT LEAST N X N WHICH ON INPUT CONTAINS THE
C         CHOLESKY LL' DECOMPOSITION OF THE MATRIX [A]. [L] IS
C         ASSUMED TO BE STORED IN THE LOWER TRIANGLE OF A. (INPUT)
C
C    IA   INTEGER GIVING THE LEADING DIMENSION OF A EXACTLY AS SPECIFIED
C         IN THE CALLING ROUTINE. (INPUT)
C
C     N   INTEGER GIVING THE SIZE OF THE MATRIX A. (INPUT)
C
C     B   REAL VECTOR OF LENGTH AT LEAST N CONTAINING THE RIGHT-HAND-SIDE
C         VECTOR OF [A]{X} = {B}. ON OUTPUT, B WILL CONTAIN THE SOLUTION.
C         (INPUT/OUTPUT)
C
C-----------------------------------------------------------------------------
      SUBROUTINE DSLVCH( A, IA, N, B )
      DIMENSION A(IA,1), B(1)
C					FORWARD SUBSTITUTION (SOLVE LY = B)
      DO 10 J = 1, N-1
         B(J) = B(J)/A(J,J)
         DO 10 I = J+1, N
            B(I) = B(I) - A(I,J)*B(J)
  10  CONTINUE
      B(N) = B(N)/(A(N,N)*A(N,N))
C					BACK-SUB (SOLVE L'X = Y)
      DO 30 I = N-1, 1, -1
         DO 20 J = I+1, N
            B(I) = B(I) - A(J,I)*B(J)
  20     CONTINUE
         B(I) = B(I)/A(I,I)
  30  CONTINUE
 
      RETURN
      END
C  ********************************************************************
C  *                                                                  *
C  *                       FUNCTION DVINT2                            *
C  *                                                                  *
C  ********************************************************************
C  DOUBLE PRECISION VERSION 1.0
C  WRITTEN BY G.A. FENTON, TUNS, 1991
C
C  PURPOSE  TO COMPUTE THE QUADRUPLE INTEGRAL OF A FUNCTION INVOLVING THE
C           DIFFERENCE OF THE COORDINATES, IE F2D(X_1-X_2, Y_1-Y_2), USING THE
C           VARIANCE FUNCTION OF F2D, CALLED V2D.
C
C  RETURNS THE QUADRUPLE INTEGRAL OF A FUNCTION ASSUMED TO HAVE THE FORM
C
C           B    D    F    H
C         INT  INT  INT  INT F2D( X_1 - X_2, Y_1 - Y_2 ) DX_1 DX_2 DY_1 DY_2
C           A    C    E    G
C
C  THE INTEGRAL IS EVALUATED USING THE USER DEFINED `VARIANCE' FUNCTION
C
C                 1        S    S    R    R
C  V2D(R,S) =  ------- * INT  INT  INT  INT F2D( X_1-X_2, Y_1-Y_2 ) DX_1 DX_2
C              R*R*S*S     0    0    0    0                           DY_1 DY_2
C
C  ARGUMENTS ARE DESCRIBED AS FOLLOWS;
C
C    V2D   EXTERNALLY SUPPLIED SINGLE PRECISION FUNCTION NAME WHICH RETURNS
C          THE VARIANCE FUNCTION OF F2D, AS DEFINED ABOVE.
C
C  B1,B2   REAL VECTORS OF LENGTH AT LEAST 4 CONTAINING RESPECTIVELY, (INPUT)
C
C             B1(1) = LOWER INTEGRATION BOUND ON X_1 (G)
C             B1(2) = UPPER INTEGRATION BOUND ON X_1 (H)
C             B1(3) = LOWER INTEGRATION BOUND ON Y_1 (C)
C             B1(4) = UPPER INTEGRATION BOUND ON Y_1 (D)
C
C             B2(1) = LOWER INTEGRATION BOUND ON X_2 (E)
C             B2(2) = UPPER INTEGRATION BOUND ON X_2 (F)
C             B2(3) = LOWER INTEGRATION BOUND ON Y_2 (A)
C             B2(4) = UPPER INTEGRATION BOUND ON Y_2 (B)
C
C   INIT   INTEGER FLAG. IF INIT = 0, THEN ONLY THE PARAMETERS OF THE
C          INTEGRATION ARE DETERMINED AND ZERO IS RETURNED. OTHERWISE
C          THE INTEGRATION IS PERFORMED. THUS TO PERFORM THE INTEGRATION,
C          THIS ROUTINE MUST BE CALLED ONCE WITH INIT = 0 AND AGAIN WITH
C          INIT > 0. (INPUT)
C ========================================================================
 
      REAL FUNCTION DVINT2( V2D, B1, B2 )
      DIMENSION B1(*), B2(*)
      EXTERNAL V2D
      DATA QUART/0.250/
 
      T1  = B1(4) - B2(3)
      T11 = T1*T1
      T2  = B1(4) - B2(4)
      T22 = T2*T2
      T3  = B1(3) - B2(3)
      T33 = T3*T3
      T4  = B1(3) - B2(4)
      T44 = T4*T4
 
      R1  = B1(2) - B2(1)
      R11 = R1*R1
      R2  = B1(2) - B2(2)
      R22 = R2*R2
      R3  = B1(1) - B2(1)
      R33 = R3*R3
      R4  = B1(1) - B2(2)
      R44 = R4*R4
 
      C1 = R11*V2D(R1,T1)-R22*V2D(R2,T1)-R33*V2D(R3,T1)+R44*V2D(R4,T1)
      C2 = R11*V2D(R1,T2)-R22*V2D(R2,T2)-R33*V2D(R3,T2)+R44*V2D(R4,T2)
      C3 = R11*V2D(R1,T3)-R22*V2D(R2,T3)-R33*V2D(R3,T3)+R44*V2D(R4,T3)
      C4 = R11*V2D(R1,T4)-R22*V2D(R2,T4)-R33*V2D(R3,T4)+R44*V2D(R4,T4)
 
      DVINT2 = QUART*(T11*C1 - T22*C2 - T33*C3 + T44*C4)
 
      RETURN
      END
C  *********************************************************************
C  *                                                                   *
C  *                          SUBROUTINE ECHO                          *
C  *                                                                   *
C  *********************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  ECHOES RFLOW'S INPUT DATA TO A SPECIFIED FILE WITH DESCRIPTIONS.
C
C---------------------------------------------------------------------------
      SUBROUTINE ECHO( ISTAT, SUB1, DEBUG, DMPFLD, VERBOS, DMPPH,
     >                 DMPRAT, NCMN, NCSD, NW, NX1, NX2, NX3, NYE,
     >                 IC1, IC2, DX, DY, JF, KSEED, THX, THY, KMNX,
     >                 KMNY, KSDX,KSDY,RHOXY)
      REAL KMNX, KMNY, KSDX, KSDY
      LOGICAL DEBUG, DMPFLD, VERBOS, DMPPH, DMPRAT
      CHARACTER*(*) SUB1
      CHARACTER*24 STRUC
      CHARACTER*16 FMT
 
   1  FORMAT(A)
   2  FORMAT(A,T40,I8)
   3  FORMAT(A,T40,E14.6)
   4  FORMAT(A,T40,L2)
   5  FORMAT(A,T40,I8,' (ignored)')
   6  FORMAT(I1)
C						CONSTRUCT SUBTITLE FROM DATA
      FMT = '(I5,A,I1,A,A)'
      I   = 1 + INT( ALOG10( FLOAT( JF )))
      WRITE(FMT(8:8),6) I
      IF( KMNX .EQ. KMNY ) THEN
         STRUC = 'isotropic permeability'
      ELSE
         STRUC = 'anisotropic permeability'
      ENDIF
      IF( NW .EQ. 0 ) THEN
         NE = NX1*NYE
         J  = 1 + INT( ALOG10( FLOAT( NE )))
         WRITE(FMT(3:3),6) J
         WRITE(SUB1,FMT) NE,' element, zero-wall problem, ',JF,
     >                 ' realizations, ',STRUC
      ELSEIF( NW .EQ. 1 ) THEN
         NE = (NX1+NX2)*NYE
         J  = 1 + INT( ALOG10( FLOAT( NE )))
         WRITE(FMT(3:3),6) J
         WRITE(SUB1,FMT) NE,' element, one-wall problem, ',JF,
     >                 ' realizations, ',STRUC
      ELSEIF( NW .EQ. 2 ) THEN
         NE = (NX1+NX2+NX3)*NYE
         J  = 1 + INT( ALOG10( FLOAT( NE )))
         WRITE(FMT(3:3),6) J
         WRITE(SUB1,FMT) NE,' element, two-wall problem, ',JF,
     >                 ' realizations, ',STRUC
      ENDIF
 
      WRITE(ISTAT,1) SUB1
      WRITE(ISTAT,'()')
      WRITE(ISTAT,1)'                                INPUT DATA'
      WRITE(ISTAT,1)'                              =============='
      WRITE(ISTAT,'()')
 
      WRITE(ISTAT,4)'Dump debug data (t/f)? . . . . . . . . ',DEBUG
      WRITE(ISTAT,4)'Dump first random field (t/f)? . . . . ',DMPFLD
      WRITE(ISTAT,4)'Report progress to user (t/f)? . . . . ',VERBOS
      WRITE(ISTAT,4)'Dump pressure head stats fields (t/f)? ',DMPPH
      WRITE(ISTAT,4)'Dump flow rate samples (t/f)?  . . . . ',DMPRAT
      WRITE(ISTAT,2)'Number of mean PH contour lines  . . . ',NCMN
      WRITE(ISTAT,2)'Number of S.D. PH contour lines  . . . ',NCSD
      WRITE(ISTAT,2)'Number of walls in the problem . . . . ',NW
      IF( NW .EQ. 0 ) THEN
      WRITE(ISTAT,2)'No. of elements in X-direction . . . . ',NX1
      WRITE(ISTAT,5)'                         (nx2) . . . . ',NX2
      WRITE(ISTAT,5)'                         (nx3) . . . . ',NX3
      WRITE(ISTAT,2)'No. of elements in Y-direction . . . . ',NYE
      WRITE(ISTAT,5)'                         (ic1) . . . . ',IC1
      WRITE(ISTAT,5)'                         (ic2) . . . . ',IC2
      ELSEIF( NW .EQ. 1 ) THEN
      WRITE(ISTAT,2)'No. of elements left of wall (X-dir) . ',NX1
      WRITE(ISTAT,2)'No. of elements right of wall (X-dir). ',NX2
      WRITE(ISTAT,5)'                         (nx3) . . . . ',NX3
      WRITE(ISTAT,2)'No. of elements in Y-direction . . . . ',NYE
      WRITE(ISTAT,2)'No. of elements against wall (Y-dir) . ',IC1
      WRITE(ISTAT,5)'                         (ic2) . . . . ',IC2
      ELSEIF( NW .EQ. 2 ) THEN
      WRITE(ISTAT,2)'No. of elements left of walls (X-dir). ',NX1
      WRITE(ISTAT,2)'No. of elements between walls (X-dir). ',NX2
      WRITE(ISTAT,2)'No. of elements right of walls(X-dir). ',NX3
      WRITE(ISTAT,2)'No. of elements in Y-direction . . . . ',NYE
      WRITE(ISTAT,2)'No. of elements against left wall (Y). ',IC1
      WRITE(ISTAT,2)'No. of elements against right wall(Y). ',IC2
      ENDIF
      WRITE(ISTAT,3)'Size of each element in X-direction. . ',DX
      WRITE(ISTAT,3)'Size of each element in Y-direction. . ',DY
      WRITE(ISTAT,2)'Number of realizations . . . . . . . . ',JF
      WRITE(ISTAT,2)'Generator seed . . . . . . . . . . . . ',KSEED
      WRITE(ISTAT,3)'Scale of fluctuation in X-direction. . ',THX
      WRITE(ISTAT,3)'Scale of fluctuation in Y-direction. . ',THY
      WRITE(ISTAT,3)'Mean permeability in X-direction . . . ',KMNX
      WRITE(ISTAT,3)'Mean permeability in Y-direction . . . ',KMNY
      WRITE(ISTAT,3)'Standard Deviation in X-direction  . . ',KSDX
      WRITE(ISTAT,3)'Standard Deviation in Y-direction  . . ',KSDY
      WRITE(ISTAT,3)'Correlation coeff. between X and Y . . ',RHOXY
      WRITE(ISTAT,1)'===================================================
     >============'
 
      RETURN
      END
C  ********************************************************************
C  *                                                                  *
C  *                        SUBROUTINE ERROR                          *
C  *                                                                  *
C  ********************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  REPORTS ERROR MESSAGES FOR RFLOW ON UNIT 6
C
C  ARGUMENTS ARE AS FOLLOWS;
C
C   IERR    INTEGER FLAG CONTAINING THE ERROR CODE ENCOUNTERED. (INPUT)
C
C   JR      INTEGER CONTAINING THE REALIZATION NUMBER FOR WHICH THE ERROR
C           OCCURRED. (INPUT)
C
C   K       UNIT NUMBER CONNECTED TO THE USER'S SCREEN. (INPUT)
C---------------------------------------------------------------------------
      SUBROUTINE ERROR( IERR, JR, K )
      LOGICAL DEBUG
      REAL AD(9,3), CD(3,3), RD(9,9), SD(9,3), Q(5,3)
      COMMON/DBGRFL/ IDBG, DEBUG
      COMMON/DBG2D/ AD, CD, RD, SD, Q, I, MDBG
 
      GO TO (1, 2, 3, 4, 5, 6 ), IABS(IERR)
 
      WRITE(K,99) 'Warning: Unknown error encountered by SIM2DK on reali
     >zation number ',JR
      RETURN
 
   1  WRITE(K,99) 'Error: Temporary vectors Z1 and Z2 are too small to c
     >ontain the required random fields!'
      CLOSE(IDBG)
      STOP
 
   2  WRITE(K,99) 'Error: Random field to FE mesh mapping failure!'
      CLOSE(IDBG)
      STOP
 
   3  WRITE(K,99) 'Error: LAS2D cannot produce such a big field!'
      CLOSE(IDBG)
      STOP
 
   4  WRITE(K,99) 'Error: LAS2D could not Cholesky decompose C matrix!'
      WRITE(K,99) '       (at stage ',I,')'
      WRITE(IDBG,99)'Error: LAS2D couldnt Cholesky decompose C matrix!'
      WRITE(IDBG,99) '      (at stage ',I,')'
      CALL DPRTMT( IDBG, CD, 3, 3, 3, 'C matrix:')
      CLOSE(IDBG)
      STOP
 
   5  WRITE(K,99) 'Error: LAS2D could not Cholesky decompose R matrix!'
      WRITE(K,99) '       (at stage ',I,')'
      WRITE(IDBG,99)'Error: LAS2D couldnt Cholesky decompose R matrix!'
      WRITE(IDBG,99) '      (at stage ',I,')'
      CALL DPRTMT( IDBG, RD, 9, 4, 4, 'R matrix:')
      CLOSE(IDBG)
      STOP
 
   6  WRITE(K,99) 'Error: RFLOW can''t handle',JR,' walls!'
      STOP
 
  99  FORMAT(A,I5,A)
      END
C  *******************************************************************
C  *                                                                 *
C  *                    FUNCTION GAUSV                               *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, PRINCETON, MAR. 19, 1988.
C
C  PURPOSE  RETURN A NORMALLY DISTRIBUTED N(0,VAR) RANDOM REALIZATION.
C
C  RETURNS A NORMALLY DISTRIBUTED, ZERO MEAN, RANDOM VARIABLE.
C  THE ARGUMENT "VAR" IS THE VARIANCE OF THE RANDOM NUMBER.
C  ENSURE THAT THE MACHINE SPECIFIC RANDOM NUMBER GENERATOR "RAND" IS
C  INITIATED IN THE CALLING ROUTINE WITH A SUITABLE SEED.
C---------------------------------------------------------------------------
      REAL FUNCTION GAUSV( VAR )
      REAL VAR, A
      DATA ONE/1.0/, TWO/2.0/, TWOPI/6.2831853071795864769/
 
      A     = TWOPI*G05CAF(A)
      R     = G05CAF(R)
      GAUSV = SQRT(-TWO*VAR*ALOG(ONE-R))*COS(A)
 
      RETURN
      END
C  **********************************************************************
C  *                                                                    *
C  *                   INTEGER FUNCTION ISEED                           *
C  *                                                                    *
C  **********************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, PRINCETON, DEC. 8, 1988.
C
C  PURPOSE  INITIALIZES THE SYSTEM PSEUDO-RANDOM NUMBER GENERATED USING
C           PROCESS ID AS DEFAULT SEED.
C
C  INITIALIZES THE F77 RANDOM NUMBER GENERATOR. IF THE ARGUMENT INTEGER
C  SEED (KSEED) IS ZERO, A RANDOM SEED IS CALCULATED BASED ON THE CURRENT
C  CLOCK TIME. THE FUNCTION RETURNS THE ACTUAL SEED USED.
C--------------------------------------------------------------------------
      INTEGER FUNCTION ISEED( KSEED )
      INTEGER KSEED
C
      ISEED = KSEED
      IF( KSEED .EQ. 0 ) THEN
C					TAKE LEAST SIG. 5 DIGITS OF WALL TIME
C				(FOR THE SCALAR AMDAHL)
          CALL CLOCK( IC )
C				(FOR THE VECTOR AMDAHL)
C         CALL TIME( IC )
          ISEED = 1 + IC - 100000*INT(IC/100000)
      ENDIF
C                                INITIALIZE THE GENERATOR
      CALL G05CBF( ISEED )
 
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                       SUBROUTINE LAS2D                          *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 2.0
C  WRITTEN BY GORDON A. FENTON, JAN. 1990. LATEST REVISION MAY, 1991.
C
C  PURPOSE  SIMULATE A 2-D GAUSSIAN RANDOM PROCESS OF LOCAL AVERAGES.
C
C  THIS ROUTINE CREATES A REALIZATION OF A 2-D RANDOM PROCESS GIVEN
C  ITS CONTINUOUS COVARIANCE FUNCTION OR ITS VARIANCE FUNCTION. EACH DISCRETE
C  VALUE GENERATED REPRESENTS THE LOCAL AVERAGE OF A REALIZATION OF
C  THE PROCESS OVER THE VOLUME DX X DY, WHERE (DX,DY) IS THE CELL SIZE
C  OF ELEMENTS IN THE FINAL RANDOM FIELD. THE CONSTRUCTION OF THE REALIZATION
C  PROCEEDS RECURSIVELY AS FOLLOWS;
C    1) GENERATE A LOW RESOLUTION FIELD (FIRST STAGE IS A 1 X 1 FIELD, BUT
C       A 5 X 5 FIELD IS GENERATED HERE TO AVOID EDGE PROBLEMS),
C    2) SUBDIVIDE THE DOMAIN BY DIVIDING EACH CELL INTO 4 EQUAL PARTS,
C    3) GENERATE NEW RANDOM VALUES FOR EACH NEW CELL. THE PARENT
C       CELLS ARE USED TO CONDITION THE NEW FIELD - MOST NOTABLE TO
C       ENSURE THAT SPATIAL CORRELATION IS APPROXIMATED AND UPWARDS
C       AVERAGING IS PRESERVED. IN THIS ALGORITHM, ONLY PARENT CELLS
C       IN A NEIGHBORHOOD OF 3 X 3 ARE CONSIDERED IN THE CONDITIONING
C       PROCESS (THUS THE SPATIAL CORRELATION IS ONLY APPROXIMATED).
C
C  THE CONDITIONING IS ACCOMPLISHED BY USING THE COVARIANCE BETWEEN LOCAL
C  AVERAGES OVER EACH CELL, CONSISTENT WITH THE GOAL OF PRODUCING A LOCAL
C  AVERAGE FIELD. NOTE THAT THE CONDITIONING PROCESS IMPLIES THAT THE
C  CONSTRUCTION OF CELLS NEAR THE EDGE OF THE BOUNDARY WILL REQUIRE THE USE
C  OF VALUES WHICH ARE, STRICTLY SPEAKING, OUTSIDE THE BOUNDARY OF THE FIELD
C  IN QUESTION. THIS IS HANDLED BY SIMULATING INITIALLY A FIELD LARGER
C  THAN REQUIRED AND DISCARDING THE EXCESS IN THE LAST STAGE. FOR THE
C  3 X 3 NEIGHBORHOOD USED IN THIS ALGORITHM, THE INITIAL FIELD MUST
C  BE OF SIZE 5 X 5. A FURTHER APPROXIMATION IS INTRODUCED AT THIS STAGE
C  SINCE ONLY A 3 X 3 CORE IS SIMULATED (HAVING THE CORRECT COVARIANCE
C  STRUCTURE) AND THE OUTER LAYER IS BUILT UP BY LETTING ELEMENTS IN
C  THE OUTER LAYER BE EQUAL TO THE CLOSEST ELEMENTS IN THE INNER CORE.
C
C  AT ANY STAGE IN THE SUBDIVISION, MEMORY MUST BE ALLOCATED FOR BOTH
C  THE NEW SUBDIVIDED FIELD AND THE PARENT FIELD (SEE Z BELOW).
C  NOTE THAT THIS ROUTINE SETS UP A NUMBER OF PARAMETERS ON THE
C  FIRST CALL AND THUS THE TIME REQUIRED TO PRODUCE THE FIRST REALIZATION
C  IS SUBSTANTIALLY GREATER THAN ON SUBSEQUENT CALLS (SEE INIT). FOR
C  MORE INFORMATION ON LOCAL AVERAGE PROCESSES, SEE EH VANMARCKE,
C  "RANDOM FIELDS: ANALYSIS AND SYNTHESIS", MIT PRESS, 1984, AND
C  GA FENTON, "SIMULATION AND ANALYSIS OF RANDOM FIELDS", PH.D. THESIS,
C  DEPT. OF CIVIL ENG. AND OP. RESEARCH, PRINCETON UNIVERSITY, PRINCETON,
C  NEW JERSEY, 1990.
C
C--------------------------------------------------------------------
C NOTES: 1) SIMULATION TIMING IS AVAILABLE THROUGH COMMON BLOCK TYMLAS WHERE
C           TI IS THE TIME REQUIRED TO SET UP THE PARAMETERS OF THE PROCESS
C           AND TS IS THE CUMULATIVE SIMULATION TIME. THESE TIMINGS ARE
C           OBTAINED THROUGH THE FUNCTION SECOND WHICH RETURNS ELAPSED USER
C           TIME IN SECONDS SINCE THE START OF THE PROGRAM.
C        2) THIS ROUTINE EMPLOYS THE FOLLOWING EXTERNAL MODULES:
C           SECOND, ISEED, DBINT2, DVINT2, DCHOLL, DSLVCH, GAUSV, VNORM,
C           DBNT2A, DGQBLK
C           AS WELL AS THE USER DEFINED (CO)VARIANCE FUNCTION.
C--------------------------------------------------------------------
C
C   ARGUMENTS ARE AS FOLLOWS;
C
C     Z   REAL VECTOR OF LENGTH AT LEAST [N*N + (0.5*N + 4)*(0.5*N + 4)]
C         WHERE N IS THE DESIRED NUMBER OF CELLS IN THE REALIZATION IN
C         EACH DIRECTION; N = 2**M. TO ACCOMPLISH THIS, Z MAY BE DIMENSIONED
C         TO BE A 2-D ARRAY IN THE CALLING PROGRAM OF SIZE (N, 1.25*N + 5)
C         FOR ANY N > 8 (M > 3). FOR VALUES OF N <= 8, USE THE ABOVE
C         RELATIONSHIP TO DETERMINE THE REQUIRED SIZE. IF THE LOGICAL FLAG
C         MEAN IS SET TO TRUE THE FIRST STAGE OR GLOBAL MEAN OF Z IS SET TO THE
C         DESIRED GLOBAL MEAN VALUE. ON OUTPUT, Z WILL CONTAIN A REALIZATION
C         OF THE 2-D PROCESS IN ITS FIRST N*N ELEMENTS - THE FIRST ELEMENT
C         CORRESPONDING TO Z(1,1), THE SECOND TO Z(2,1), ETC. WHERE COLUMNS
C         TRAVERSE THE (POSITIVE) X DIRECTION AND ROWS TRAVERSE THE (POSITIVE)
C         Y DIRECTION. (OUTPUT)
C
C     M   DESIRED NUMBER OF SUBDIVISIONS IN THE PROCESS (SUCH THAT N = 2**M).
C         (INPUT)
C
C XL, YL  PHYSICAL DIMENSIONS OF THE PROCESS. (INPUT)
C
C    COV  EXTERNAL REAL*8 FUNCTION WHICH RETURNS EITHER THE COVARIANCE
C         BETWEEN POINTS AT A GIVEN LAG, OR RETURNS THE VARIANCE OF THE
C         PROCESS AVERAGED OVER A GIVEN AREA. IN EITHER CASE,
C         COV IS REFERENCED AS FOLLOWS
C
C                COV( V1, V2 )
C
C         WHERE (V1,V2) ARE EITHER THE COMPONENTS OF THE LAG OR THE SIZE
C         OF THE RECTANGULAR AVERAGING DOMAIN (SEE ITYPE). NOTE THAT THE
C         VARIANCE OF THE PROCESS AVERAGED OVER THE AREA (V1 X V2) IS THE
C         PRODUCT OF THE POINT VARIANCE AND THE VARIANCE FUNCTION, AS
C         DISCUSSED BY VANMARCKE (PG 186). THE PARAMETERS OF THE FUNCTION
C         SHOULD BE PASSED BY A COMMON BLOCK FROM THE CALLING ROUTINE.
C
C  ITYPE  INTEGER FLAG WHICH DENOTES THE TYPE OF RESULT RETURNED BY THE
C         USER SUPPLIED FUNCTION "COV"; (INPUT)
C          = 1 IF COV RETURNS THE POINT COVARIANCE FUNCTION OF THE PROCESS
C              (V1 AND V2 ARE LAGS IN PHYSICAL SPACE)
C          = 2 IF COV RETURNS THE VARIANCE OF A LOCAL AVERAGE OF THE PROCESS
C              (V1 AND V2 ARE THE DIMENSIONS OF THE AVERAGING AREA)
C
C  KSEED  INTEGER SEED TO BE USED FOR THE PSEUDO-RANDOM NUMBER GENERATOR.
C         IF KSEED = 0, THEN A RANDOM SEED WILL BE USED (BASED ON THE
C         CLOCK TIME WHEN THIS ROUTINE IS CALLED FOR THE FIRST TIME).
C         ON OUTPUT, KSEED IS SET TO THE VALUE OF THE ACTUAL SEED USED.
C
C  INIT   INTEGER FLAG WHICH MUST BE 1 WHEN PARAMETERS OF THE PROCESS
C         ARE TO BE CALCULATED OR RECALCULATED. IF MULTIPLE REALIZATIONS
C         OF THE SAME PROCESS ARE DESIRED, THEN SUBSEQUENT CALLS SHOULD
C         USE INIT EQUAL TO 0 AND M LESS THAN OR EQUAL TO THE VALUE USED
C         INITIALLY.
C
C  MEAN   LOGICAL FLAG WHICH IS TRUE IF THE GLOBAL MEAN OF THE PROCESS IS
C         ALREADY SPECIFIED IN THE FIRST ELEMENT OF Z. IF MEAN IS FALSE,
C         THEN A GLOBAL MEAN IS GENERATED RANDOMLY ACCORDING TO A ZERO MEAN
C         GAUSSIAN DISTRIBUTION WITH VARIANCE EQUAL TO THAT OF THE PROCESS
C         AVERAGED OVER THE ENTIRE DOMAIN, XL X YL.
C
C  IERR   INTEGER FLAG WHICH DENOTES THE STATUS OF THE RUN;
C         =  0  FOR SUCCESSFUL SIMULATION
C         = -2  IF ITYPE IS NOT 1 OR 2
C         = -3  IF THE REQUIRED PROCESS IS TOO LARGE (SEE PARAMETERS BELOW)
C         = -4  IF ANY OF THE C COEFFICIENTS ARE FOUND TO BE IMAGINARY
C         NOTE THAT ERROR -4 IS USUALLY A RESULT OF NUMERICAL
C         ERRORS IN THE INTEGRATION OF THE POINT COVARIANCE FUNCTION
C         (ITYPE = 1). IT IS USUALLY BEST TO USE THE LOCAL AVERAGE VARIANCE
C         FORMULATION (ITYPE = 2) WHERE COV RETURNS THE VARIANCE OF A
C         LOCAL AVERAGE. THE LOCAL AVERAGE VARIANCE CAN BE FOUND BY ANALYTICAL
C         INTEGRATION OF THE COVARIANCE OR (EQUIVALENTLY) THE CORRELATION
C         FUNCTION (SEE VANMARCKE, EQ. 5.1.6, OR FENTON, EQS. 1.31-1.33).
C
C---------------------------------------------------------------------------
C  PARAMETERS:
C    MX   REPRESENTS THE MAXIMUM NUMBER OF SUBDIVISIONS THAT THE ROUTINE
C         CAN CARRY OUT. NOTE THAT THIS MEANS THAT THE MAXIMUM PROCESS SIZE
C         IS (2**MX X 2**MX)
C   NGS   REPRESENTS THE MAXIMUM NUMBER OF RANDOM GAUSSIAN VARIATES THAT CAN
C         BE PRODUCED ON EACH CALL TO VNORM. NGS SHOULD BE 3*[2**(MX-1)], BUT
C         NOT LESS THAN 9.
C ==========================================================================
      SUBROUTINE LAS2D( Z, M, XL, YL, COV, ITYPE, KSEED,
     >                  INIT, MEAN, IERR )
      PARAMETER( MX = 12, NGS = 6144 )
      DIMENSION Z(1)
      REAL AD(9,3), CD(3,3), RD(9,9), SD(9,3)
      REAL QR1(4), QR2(4), QS1(4), QS2(4), QS3(4), Q(5,3)
C     REAL DR, DS, T1, T2, COV, DBLE, DBINT2, DVINT2, DCHOLL
      DIMENSION A(9,3,MX), C(6,MX), R(45)
      DIMENSION U(NGS)
      LOGICAL MEAN
      EXTERNAL COV
C                                 SAVE PARAMETERS FOR LATER REALIZATIONS
      SAVE A, C, R
 
      COMMON/TYMLAS/ TI, TS
      COMMON/DBG2D/ AD, CD, RD, SD, Q, I, MDBG
 
      DATA IFIRST/1/, ZERO/0.0/, THREE/3.0/, FOUR/4.0/, PT01/0.01/
C-------------------------------------------------------------------------
 
      IERR = 0
C				INITIALIZE
 
      IF( IFIRST .EQ. 1 .OR. INIT .EQ. 1 ) THEN
C					CHECK DATA
         IF( ITYPE .NE. 1 .AND. ITYPE .NE. 2 ) THEN
            IERR = -2
            RETURN
         ENDIF
         IF( M .GT. MX ) THEN
            IERR = -3
            RETURN
         ENDIF
         IFIRST = 0
C					FOR DEBUG PURPOSES
         MDBG = M
C					START TIMER
         TI = SECOND()
C					START GENERATOR
         KSEED = ISEED( KSEED )
C					DETERMINE COVARIANCES BETWEEN CELLS
         T1   = XL
         T2   = YL
         DR   = 1.0/(T1*T1*T2*T2)
         DS   = 4.0*DR
         NG   = 10
C						FOR STAGES 1, 2, ..., M
         DO 80 I = 1, M
C						SET INTEGRATION BOUNDS
            QR1(1) = 2.0*T1
            QR1(2) = 3.0*T1
            QR1(3) = 0.0
            QR1(4) = T2
 
            QS1(1) = T1
            QS1(2) = 1.50*T1
            QS1(3) = T2
            QS1(4) = 1.50*T2
 
            QS2(1) = 1.50*T1
            QS2(2) = 2.00*T1
            QS2(3) = T2
            QS2(4) = 1.50*T2
 
            QS3(1) = T1
            QS3(2) = 1.50*T1
            QS3(3) = 1.50*T2
            QS3(4) = 2.00*T2
 
C						NOW INTEGRATE TO FIND COV.
            IF( ITYPE .EQ. 1 ) THEN
C							COV RETURNS COVARIANCE
               L = 0
               DO 20 K = 1, 3
                  QR2(3) = FLOAT(K-1)*T2
                  QR2(4) = QR2(3) + T2
                  DO 10 J = 1, 3
                     L = L + 1
                     QR2(1)  = FLOAT(J-1)*T1
                     QR2(2)  = QR2(1) + T1
                     Q(J,K)  = DR*DBINT2( COV, QR1, QR2, NG )
                     SD(L,1) = DS*DBINT2( COV, QS1, QR2, NG )
                     SD(L,2) = DS*DBINT2( COV, QS2, QR2, NG )
                     SD(L,3) = DS*DBINT2( COV, QS3, QR2, NG )
                     AD(L,1) = SD(L,1)
                     AD(L,2) = SD(L,2)
                     AD(L,3) = SD(L,3)
  10              CONTINUE
                  QR2(1) = 3.0*T1
                  QR2(2) = QR2(1) + T1
                  Q(4,K) = DR*DBINT2( COV, QR1, QR2, NG )
                  QR2(1) = QR2(2)
                  QR2(2) = QR2(1) + T1
                  Q(5,K) = DR*DBINT2( COV, QR1, QR2, NG )
  20           CONTINUE
            ELSE
C							COV RETURNS VAR FUNCT'N
               L = 0
               DO 40 K = 1, 3
                  QR2(3) = FLOAT(K-1)*T2
                  QR2(4) = QR2(3) + T2
                  DO 30 J = 1, 3
                     L = L + 1
                     QR2(1)  = FLOAT(J-1)*T1
                     QR2(2)  = QR2(1) + T1
                     Q(J,K)  = DR*DVINT2( COV, QR1, QR2 )
                     SD(L,1) = DS*DVINT2( COV, QS1, QR2 )
                     SD(L,2) = DS*DVINT2( COV, QS2, QR2 )
                     SD(L,3) = DS*DVINT2( COV, QS3, QR2 )
                     AD(L,1) = SD(L,1)
                     AD(L,2) = SD(L,2)
                     AD(L,3) = SD(L,3)
  30              CONTINUE
                  QR2(1) = 3.0*T1
                  QR2(2) = QR2(1) + T1
                  Q(4,K) = DR*DVINT2( COV, QR1, QR2 )
                  QR2(1) = QR2(2)
                  QR2(2) = QR2(1) + T1
                  Q(5,K) = DR*DVINT2( COV, QR1, QR2 )
  40           CONTINUE
            ENDIF
C						FILL IN COVARIANCE MATRIX
            DO 50 L = 1, 9
               LL = INT( (FLOAT(L)/THREE) - PT01 )
               DO 50 K = L, 9
                  KK = INT( (FLOAT(K)/THREE) - PT01 )
                  II = K - L + 3*(LL - KK + 1)
                  JJ = 1 + KK - LL
                  RD(K,L) = Q(II,JJ)
  50        CONTINUE
C						CHOLESKY FACTORIZATION OF C
            IF( I .GT. 1 ) THEN
               CD(1,1) = RD(1,1) - CD(1,1)
               CD(2,1) = RD(2,1) - CD(2,1)
               CD(3,1) = RD(4,1) - CD(3,1)
               CD(2,2) = RD(2,2) - CD(2,2)
               CD(3,2) = RD(4,2) - CD(3,2)
               CD(3,3) = RD(4,4) - CD(3,3)
               IF( DCHOLL( CD, 3, 3 ) .EQ. 0.0 ) THEN
                  IERR = -4
                  RETURN
               ENDIF
C							SAVE IN REAL*4
               J = I - 1
               C(1,J) = CD(1,1)
               C(2,J) = CD(2,1)
               C(3,J) = CD(2,2)
               C(4,J) = CD(3,1)
               C(5,J) = CD(3,2)
               C(6,J) = CD(3,3)
            ENDIF
C						CHOLESKY FACTORIZATION OF R
            IF( DCHOLL( RD, 9, 9 ) .EQ. 0.0 ) THEN
               IERR = -5
               RETURN
            ENDIF
C						KEEP FIRST STAGE R
            IF( I .EQ. 1 ) THEN
               L = 0
               DO 60 J = 1, 9
               DO 60 K = 1, J
                  L = L + 1
                  R(L) = RD(J,K)
  60           CONTINUE
            ENDIF
C						SOLVE [R]{A} = {S} FOR A
            CALL DSLVCH( RD, 9, 9, AD(1,1) )
            CALL DSLVCH( RD, 9, 9, AD(1,2) )
            CALL DSLVCH( RD, 9, 9, AD(1,3) )
C							SAVE A IN REAL*4
            DO 70 K = 1, 3
               A(1,K,I) = AD(1,K)
               A(2,K,I) = AD(2,K)
               A(3,K,I) = AD(3,K)
               A(4,K,I) = AD(4,K)
               A(5,K,I) = AD(5,K)
               A(6,K,I) = AD(6,K)
               A(7,K,I) = AD(7,K)
               A(8,K,I) = AD(8,K)
               A(9,K,I) = AD(9,K)
C						COMPUTE SECOND HALF OF [C]
               DO 70 J = K, 3
                  CD(J,K) = AD(1,J)*SD(1,K)
                  DO 70 L = 2, 9
                     CD(J,K) = CD(J,K) + AD(L,J)*SD(L,K)
  70        CONTINUE
C						SET LENGTH OF NEXT STAGE CELL
            T1 = 0.50*T1
            T2 = 0.50*T2
            DR = 1.0/(T1*T1*T2*T2)
            DS = 4.0*DR
            NG = NG - 1
            IF( NG .LT. 5 ) NG = 5
  80     CONTINUE
C						SET FINAL [C] MATRIX
         QR1(1) = 0.0
         QR1(2) = T1
         QR1(3) = 0.0
         QR1(4) = T2
         QR2(1) = 0.0
         QR2(2) = T1
         QR2(3) = 0.0
         QR2(4) = T2
 
         IF( ITYPE .EQ. 1 ) THEN
            RD(1,1) = DR*DBINT2( COV, QR1, QR2, NG )
            QR2(1)  = T1
            QR2(2)  = 2.0*T1
            RD(2,1) = DR*DBINT2( COV, QR1, QR2, NG )
            QR2(1)  = 0.0
            QR2(2)  = T1
            QR2(3)  = T2
            QR2(4)  = 2.0*T2
            RD(4,1) = DR*DBINT2( COV, QR1, QR2, NG )
            QR1(1)  = T1
            QR1(2)  = 2.0*T1
            RD(4,2) = DR*DBINT2( COV, QR1, QR2, NG )
         ELSE
            RD(1,1) = DR*DVINT2( COV, QR1, QR2 )
            QR2(1)  = T1
            QR2(2)  = 2.0*T1
            RD(2,1) = DR*DVINT2( COV, QR1, QR2 )
            QR2(1)  = 0.0
            QR2(2)  = T1
            QR2(3)  = T2
            QR2(4)  = 2.0*T2
            RD(4,1) = DR*DVINT2( COV, QR1, QR2 )
            QR1(1)  = T1
            QR1(2)  = 2.0*T1
            RD(4,2) = DR*DVINT2( COV, QR1, QR2 )
         ENDIF
 
         CD(1,1) = RD(1,1) - CD(1,1)
         CD(2,1) = RD(2,1) - CD(2,1)
         CD(3,1) = RD(4,1) - CD(3,1)
         CD(2,2) = RD(1,1) - CD(2,2)
         CD(3,2) = RD(4,2) - CD(3,2)
         CD(3,3) = RD(1,1) - CD(3,3)
         IF( DCHOLL( CD, 3, 3 ) .EQ. 0.0 ) THEN
            IERR = -4
            RETURN
         ENDIF
         C(1,M) = CD(1,1)
         C(2,M) = CD(2,1)
         C(3,M) = CD(2,2)
         C(4,M) = CD(3,1)
         C(5,M) = CD(3,2)
         C(6,M) = CD(3,3)
C						CLEAN UP
         TS   = ZERO
         TI   = SECOND() - TI
      ENDIF
 
C------------------------------ CREATE THE REALIZATION OF THE PROCESS
 
      TT = SECOND()
 
      IF( M .LE. 0 ) THEN
         Z(1) = GAUSV( R(1)*R(1) )
         RETURN
      ENDIF
      N  = 2**M
      NN = N*N
      ZM = Z(1)
      IF( MOD(M,2) .EQ. 0 ) THEN
         IN = 0
         IO = NN
      ELSE
         IN = NN
         IO = 0
      ENDIF
C					GENERATE STAGE 0 FIELD
      CALL VNORM( U, 9 )
C						THE INNER 3X3 FIELD
      Z(IN+ 7) = R( 1)*U(1)
      Z(IN+ 8) = R( 2)*U(1)+R( 3)*U(2)
      Z(IN+ 9) = R( 4)*U(1)+R( 5)*U(2)+R( 6)*U(3)
      Z(IN+12) = R( 7)*U(1)+R( 8)*U(2)+R( 9)*U(3)+R(10)*U(4)
      Z(IN+13) = R(11)*U(1)+R(12)*U(2)+R(13)*U(3)+R(14)*U(4)+R(15)*U(5)
      Z(IN+14) = R(16)*U(1)+R(17)*U(2)+R(18)*U(3)+R(19)*U(4)+R(20)*U(5)
     >         + R(21)*U(6)
      Z(IN+17) = R(22)*U(1)+R(23)*U(2)+R(24)*U(3)+R(25)*U(4)+R(26)*U(5)
     >         + R(27)*U(6)+R(28)*U(7)
      Z(IN+18) = R(29)*U(1)+R(30)*U(2)+R(31)*U(3)+R(32)*U(4)+R(33)*U(5)
     >         + R(34)*U(6)+R(35)*U(7)+R(36)*U(8)
      Z(IN+19) = R(37)*U(1)+R(38)*U(2)+R(39)*U(3)+R(40)*U(4)+R(41)*U(5)
     >         + R(42)*U(6)+R(43)*U(7)+R(44)*U(8)+R(45)*U(9)
C						FUDGE THE OUTER RING TO GET 5X5
      Z(IN+ 1) = Z(IN+ 7)
      Z(IN+ 2) = Z(IN+ 7)
      Z(IN+ 3) = Z(IN+ 8)
      Z(IN+ 4) = Z(IN+ 9)
      Z(IN+ 5) = Z(IN+ 9)
      Z(IN+ 6) = Z(IN+ 7)
      Z(IN+10) = Z(IN+ 9)
      Z(IN+11) = Z(IN+12)
      Z(IN+15) = Z(IN+14)
      Z(IN+16) = Z(IN+17)
      Z(IN+20) = Z(IN+19)
      Z(IN+21) = Z(IN+17)
      Z(IN+22) = Z(IN+17)
      Z(IN+23) = Z(IN+18)
      Z(IN+24) = Z(IN+19)
      Z(IN+25) = Z(IN+19)
C					CONDITION ON THE MEAN?
      IF( MEAN ) THEN
         DIF  = ZM - Z(IN+13)
         DO 100 K = 1, 25
            Z(IN+K) = Z(IN+K) + DIF
 100     CONTINUE
      ENDIF
C					GENERATE 1ST, 2ND, ... M FIELDS
      IP = 1
      DO 120 I = 1, M
C						SWAP CURRENT AND PREV FIELDS
         IT = IO
         IO = IN
         IN = IT
 
         IC = IP + 4
         IF( I .LT. M ) THEN
            IQ = IP + 2
            NC = 2*IP + 4
            K0 = IO - 1
            II = 2
         ELSE
            IQ = IP
            NC = 2*IP
            K0 = IO + IC - 2
            II = 4
         ENDIF
C						BUILD NEW FIELD FROM PREVIOUS
         N1 = IN + 1
         DO 110 K = 1, IQ
            K0 = K0 + II
            K1 = K0 + IC
            K2 = K1 + IC
            N0 = N1
            N1 = N0 + NC
            L   = 1
            CALL VNORM( U, 3*IQ )
         DO 110 J = 1, IQ
            Z(N0  )=A(1,1,I)*Z(K0) + A(2,1,I)*Z(K0+1) + A(3,1,I)*Z(K0+2)
     >             +A(4,1,I)*Z(K1) + A(5,1,I)*Z(K1+1) + A(6,1,I)*Z(K1+2)
     >             +A(7,1,I)*Z(K2) + A(8,1,I)*Z(K2+1) + A(9,1,I)*Z(K2+2)
     >             + C(1,I)*U(L)
            Z(N0+1)=A(1,2,I)*Z(K0) + A(2,2,I)*Z(K0+1) + A(3,2,I)*Z(K0+2)
     >             +A(4,2,I)*Z(K1) + A(5,2,I)*Z(K1+1) + A(6,2,I)*Z(K1+2)
     >             +A(7,2,I)*Z(K2) + A(8,2,I)*Z(K2+1) + A(9,2,I)*Z(K2+2)
     >             + C(2,I)*U(L)   + C(3,I)*U(L+1)
            Z(N1  )=A(1,3,I)*Z(K0) + A(2,3,I)*Z(K0+1) + A(3,3,I)*Z(K0+2)
     >             +A(4,3,I)*Z(K1) + A(5,3,I)*Z(K1+1) + A(6,3,I)*Z(K1+2)
     >             +A(7,3,I)*Z(K2) + A(8,3,I)*Z(K2+1) + A(9,3,I)*Z(K2+2)
     >             + C(4,I)*U(L)   + C(5,I)*U(L+1)    + C(6,I)*U(L+2)
            Z(N1+1)=FOUR*Z(K1+1) - Z(N0) - Z(N0+1) - Z(N1)
 
            K0 = K0 + 1
            K1 = K1 + 1
            K2 = K2 + 1
            N0 = N0 + 2
            N1 = N1 + 2
            L = L + 3
 110     CONTINUE
 
         IP = 2*IP
 120  CONTINUE
C						ALL DONE, COMPUTE ELAPSED TIME
      TS = TS + (SECOND() - TT)
 
      RETURN
      END
C  ********************************************************************
C  *                                                                  *
C  *                       SUBROUTINE MAP2FE                          *
C  *                                                                  *
C  ********************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE   TO MAP THE AVERAGE OF ELEMENTS FROM ONE MATRIX INTO ANOTHER
C
C  THIS ROUTINE TAKES TWO MATRICES A, OF SIZE N X N, AND B, OF SIZE NX X NY,
C  AND MAPS ELEMENTS OF A INTO B. MORE SPECIFICALLY EACH GROUP OF (IX X IY)
C  ELEMENTS IN A ARE MAPPED INTO A SINGLE ELEMENT OF B BY COMPUTING THEIR
C  AVERAGE. IT IS ASSUMED THAT N >= MAX( IX*NX, IY*NY ).
C  ARGUMENTS ARE DESCRIBED AS FOLLOWS;
C
C      A     REAL ARRAY OF SIZE N X N CONTAINING THE ELEMENTS TO BE MAPPED
C            INTO B. (INPUT)
C
C      N     ORDER OF THE MATRIX A. (INPUT)
C
C      B     REAL ARRAY OF SIZE AT LEAST NX X NY WHICH ON OUTPUT WILL CONTAIN
C            THE AVERAGE OF CORRESPONDING (IX X IY) ELEMENTS OF A. (OUTPUT)
C
C      IB    LEADING DIMENSION OF B AS SPECIFIED IN THE CALLING ROUTINE.
C            (INPUT)
C
C      NX    COLUMN (FIRST DIMENSION) SIZE OF B. (INPUT)
C
C      NY    ROW (SECOND DIMENSION) SIZE OF B. (INPUT)
C
C      IX    NUMBER OF ELEMENTS OF A CORRESPONDING TO ONE ELEMENT OF B IN
C            THE X DIRECTION. (INPUT)
C
C      IY    NUMBER OF ELEMENTS OF A CORRESPONDING TO ONE ELEMENT OF B IN
C            THE Y DIRECTION. (INPUT)
C----------------------------------------------------------------------------
      SUBROUTINE MAP2FE( A, N, B, IB, NX, NY, IX, IY )
      DIMENSION A(N,1), B(IB,1)
      DATA ZERO/0.0/, ONE/1.0/
 
      DIV = 1./FLOAT(IX*IY)
      JJ  = 1 - IY
      IX1 = IX - 1
      IY1 = IY - 1
 
      DO 20 J = 1, NY
         II = 1
         JJ = JJ + IY
      DO 20 I = 1, NX
         S = ZERO
         DO 10 L = JJ, JJ + IY1
         DO 10 K = II, II + IX1
            S = S + A(K,L)
  10     CONTINUE
         B(I,J) = S*DIV
         II = II + IX
  20  CONTINUE
 
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                      SUBROUTINE NOCUT                           *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY D. VAUGHAN GRIFFITHS, 1991
C
C  PURPOSE  FORMS THE COORDINATES AND STEERING VECTOR FOR 4-NODE
C           QUADRILATERALS COUNTING IN Y-DIRECTION (LAPLACES' EQUATION,
C           0-WALL PROBLEM = 1-D FLOW ALONG A CHANNEL)
C
C--------------------------------------------------------------------------
      SUBROUTINE NOCUT(IP,IQ,NYE,AA,BB,COORD,G,NF,INF)
      REAL COORD(4,2)
      INTEGER G(*),NF(INF,*), A1, A2, A3, A4
 
      A1=(IP-1)*(NYE+1)+IQ+1
      A2=A1-1
      A3=IP*(NYE+1)+IQ
      A4=A3+1
      G(1)=NF(A1,1)
      G(2)=NF(A2,1)
      G(3)=NF(A3,1)
      G(4)=NF(A4,1)
      COORD(1,1) = FLOAT(IP-1)*AA
      COORD(2,1) = FLOAT(IP-1)*AA
      COORD(3,1) = FLOAT(IP)  *AA
      COORD(4,1) = FLOAT(IP)  *AA
      COORD(1,2) =-FLOAT(IQ)  *BB
      COORD(2,2) =-FLOAT(IQ-1)*BB
      COORD(3,2) =-FLOAT(IQ-1)*BB
      COORD(4,2) =-FLOAT(IQ)  *BB
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                      SUBROUTINE ONECUT                          *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY D. VAUGHAN GRIFFITHS, 1991
C
C  PURPOSE  FORMS THE COORDINATES AND STEERING VECTOR FOR 4-NODE
C           QUADRILATERALS COUNTING IN Y-DIRECTION (LAPLACES' EQUATION,
C           1-WALL PROBLEM)
C
C--------------------------------------------------------------------------
      SUBROUTINE ONECUT(IP,IQ,NX1,NYE,IC1,AA,BB,COORD,G,NF,INF)
      REAL COORD(4,2)
      INTEGER G(*),NF(INF,*), A1, A2, A3, A4
 
      A1=(IP-1)*(NYE+1)+IQ+1
      A2=A1-1
      A3=IP*(NYE+1)+IQ
      A4=A3+1
      IF(IP.EQ.NX1+1)THEN
        IF(IQ.LT.IC1)THEN
          A1=A1+NYE+1
          A2=A2+NYE+1
          A3=A3+IC1
          A4=A4+IC1
        ELSE IF(IQ.EQ.IC1)THEN
          A1=A1
          A2=A2+NYE+1
          A3=A3+IC1
          A4=A4+IC1
        ELSE
          A1=A1
          A2=A2
          A3=A3+IC1
          A4=A4+IC1
        END IF
      ELSE IF(IP.GT.NX1+1)THEN
          A1=A1+IC1
          A2=A2+IC1
          A3=A3+IC1
          A4=A4+IC1
      END IF
      G(1)=NF(A1,1)
      G(2)=NF(A2,1)
      G(3)=NF(A3,1)
      G(4)=NF(A4,1)
      COORD(1,1) = FLOAT(IP-1)*AA
      COORD(2,1) = FLOAT(IP-1)*AA
      COORD(3,1) = FLOAT(IP)  *AA
      COORD(4,1) = FLOAT(IP)  *AA
      COORD(1,2) =-FLOAT(IQ)  *BB
      COORD(2,2) =-FLOAT(IQ-1)*BB
      COORD(3,2) =-FLOAT(IQ-1)*BB
      COORD(4,2) =-FLOAT(IQ)  *BB
      RETURN
      END
C  ********************************************************************
C  *                                                                  *
C  *                       SUBROUTINE OPENFL                          *
C  *                                                                  *
C  ********************************************************************
C  INTEGER VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  TO OPEN A SET OF OUTPUT FILES HAVING THE SAME BASENAME AS
C           THE DATA FILE
C
C  ARGUMENTS ARE AS FOLLOWS;
C
C     ITERM   IS SET BY THIS ROUTINE TO THE VALUE APPROPRIATE TO PRINT
C             TO YOUR TERMINAL. (OUTPUT)
C
C     ISTAT   OUTPUT UNIT FOR JOB STATISTICS AND GLOBAL RESULTS. (INPUT)
C
C     IRATE   OUTPUT UNIT FOR INDIVIDUAL FLOW RATE REALIZATIONS. (INPUT)
C
C     IFLD    OUTPUT UNIT FOR RANDOM FIELD DISPLAY. (INPUT)
C
C     IMDL    OUTPUT UNIT FOR FE MESH LOG-PERMEABILITY DISPLAY. (INPUT)
C
C     DMPFLD  LOGICAL FLAG WHICH IS TRUE IF RANDOM FIELD (IFLD) AND MESH
C             FIELDS (IMDL) ARE TO BE OUTPUT. (INPUT)
C
C     NW      THE NUMBER OF WALLS IN THE PROBLEM. (INPUT)
C
C     DMPPH   LOGICAL FLAG WHICH IS TRUE IF PRESSURE HEAD DATA IS TO BE
C             OUTPUT. (INPUT)
C
C     DMPRAT  LOGICAL FLAG WHICH IS TRUE IF FLOW RATE, EXIT GRADIENT
C             SAMPLES ARE TO BE DUMPED TO A FILE. (INPUT)
C
C     IP???   OUTPUT UNIT NUMBERS FOR A SET OF PRESSURE HEAD
C             STATISTICS FIELDS. (INPUT)
C		IPMNL = PRESSURE MEAN (ALL 0-WALL, LEFT 1 AND 2-WALL)
C		IPMNC = PRESSURE MEAN (RIGHT 1-WALL, CENTER 2-WALL)
C		IPMNR = PRESSURE MEAN (RIGHT 2-WALL)
C		IPSDL = PRESS. STAND. DEV. (ALL 0-WALL, LEFT 1 AND 2-WALL)
C		IPSDC = PRESS. STAND. DEV. (RIGHT 1-WALL, CENTER 2-WALL)
C		IPSDR = PRESS. STAND. DEV. (RIGHT 2-WALL)
C
C---------------------------------------------------------------------------
      SUBROUTINE OPENFL( ITERM, ISTAT, IRATE, IFLD, IMDL, DMPFLD, NW,
     >        DMPPH, DMPRAT, IPMNL, IPMNC, IPMNR, IPSDL, IPSDC, IPSDR )
      LOGICAL DMPFLD, DMPPH, DMPRAT
 
C						SET STANDARD OUTPUT UNIT
      ITERM = 7
 
      RETURN
      END
C  ********************************************************************
C  *                                                                  *
C  *                       SUBROUTINE OPENIN                          *
C  *                                                                  *
C  ********************************************************************
C  INTEGER VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  TO OPEN A DATA FILE GIVEN BY AN ARGUMENT ON THE COMMAND LINE
C           (ON THE AMDAHL, THIS DOES NOTHING)
C
C  THE PARAMETER `KARG' IS PROVIDED BECAUSE NOT ALL COMPILERS NUMBER
C  THEIR COMMAND LINE ARGUMENTS THE SAME WAY (IE HP-UX USES KARG = 1).
C  ARGUMENTS ARE AS FOLLOWS;
C
C     IIN     UNIT NUMBER TO WHICH THE INPUT FILE SHOULD BE OPENED. (INPUT)
C
C---------------------------------------------------------------------------
      SUBROUTINE OPENIN( IIN )
 
      RETURN
      END
C  ******************************************************************
C  *                                                                *
C  *                      SUBROUTINE PHSTAT                         *
C  *                                                                *
C  ******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  TO OUTPUT PRESSURE HEAD STATISTICS AT NODAL POINTS IN
C           A DISPLAY FORMAT FILE
C
C  ARGUMENTS TO THE ROUTINE ARE AS FOLLOWS;
C
C    JOB    CHARACTER STRING CONTAINING THE TITLE OF THE RUN. (INPUT)
C
C    SUB1   CHARACTER STRING CONTAINING THE RUN SUBTITLE. (INPUT)
C
C    IP???  PRESSURE STATISTICS FIELD OUTPUT UNIT NUMBERS DEFINED AS
C           FOLLOWS;
C		IPMNL = PRESSURE MEAN (ALL 0-WALL, LEFT 1 AND 2-WALL)
C		IPMNC = PRESSURE MEAN (RIGHT 1-WALL, CENTER 2-WALL)
C		IPMNR = PRESSURE MEAN (RIGHT 2-WALL)
C		IPSDL = PRESS. STAND. DEV. (ALL 0-WALL, LEFT 1 AND 2-WALL)
C		IPSDC = PRESS. STAND. DEV. (RIGHT 1-WALL, CENTER 2-WALL)
C		IPSDR = PRESS. STAND. DEV. (RIGHT 2-WALL)
C
C    ITERM  UNIT NUMBER CONNECTED TO THE USER'S SCREEN. (INPUT)
C
C    PHMN   ON INPUT, PHMN CONTAINS THE PRESSURE HEAD AT EACH NODE SUMMED
C           OVER THE TOTAL NUMBER OF REALIZATIONS. ON OUPUT, PHMN WILL
C           CONTAIN THE MEAN PRESSURE HEAD AT EACH NODE. (INPUT/OUTPUT)
C
C    PHSD   ON INPUT, PHSD CONTAINS THE SQUARED PRESSURE HEAD AT EACH
C           NODE SUMMED OVER THE TOTAL NUMBER OF REALIZATIONS. ON OUTPUT,
C           PHSD WILL CONTAIN THE STANDARD DEVIATION OF THE PRESSURE
C           HEAD AT EACH NODE. (INPUT/OUTPUT)
C
C    N      THE NUMBER OF NODES AT WHICH PHMN AND PHSD ARE DEFINED. (INPUT)
C
C    NW     THE NUMBER OF WALLS. (INPUT)
C
C    NS     THE NUMBER OF REALIZATIONS OVER WHICH THE STATISTICS ARE
C           TO BE AVERAGED (IF NS = 1, THEN OUTPUT ESTIMATED VARIANCES ARE
C           SET TO ZERO RATHER THAN INFINITY). (INPUT)
C
C    NX1    THE NUMBER OF ELEMENTS IN THE X DIRECTION (0-WALL CASE) OR TO
C           THE LEFT OF THE LEFT-MOST WALL (1 OR 2-WALL CASES). (INPUT)
C
C    NX2    THE NUMBER OF ELEMENTS IN THE X DIRECTION TO THE RIGHT OF
C           THE WALL (1-WALL CASE) OR BETWEEN WALLS (2-WALL CASE). (INPUT)
C
C    NX3    THE NUMBER OF ELEMENTS IN THE X DIRECTION TO THE RIGHT OF
C           THE RIGHT-MOST WALL (2-WALL CASE). (INPUT)
C
C    NYE    THE TOTAL NUMBER OF ELEMENTS IN THE Y DIRECTION. (INPUT)
C
C    IC1    THE NUMBER OF ELEMENTS IN THE Y DIRECTION AGAINST THE LEFT-MOST
C           WALL (1 OR 2-WALL CASES). (INPUT)
C
C    IC2    THE NUMBER OF ELEMENTS IN THE Y DIRECTION AGAINST THE RIGHT-MOST
C           WALL (2-WALL CASE). (INPUT)
C
C     DX    PHYSICAL DIMENSION OF EACH FINITE ELEMENT IN THE X DIRECTION.
C           (INPUT)
C
C     DY    PHYSICAL DIMENSION OF EACH FINITE ELEMENT IN THE Y DIRECTION.
C           (INPUT)
C
C    NCMN   THE NUMBER OF CONTOURING LINES TO SPECIFY IN THE DISPLAY FILE
C           FOR THE MEAN NODAL PRESSURE HEADS. (INPUT)
C
C    NCSD   THE NUMBER OF CONTOURING LINES TO SPECIFY IN THE DISPLAY FILE
C           FOR THE STANDARD DEVIATION OF NODAL PRESSURE HEADS. (INPUT)
C---------------------------------------------------------------------------
      SUBROUTINE PHSTAT( JOB, SUB1, IPMNL, IPMNC, IPMNR, IPSDL, IPSDC,
     >                   IPSDR, ITERM, PHMN, PHSD, N, NW, NS, NX1, NX2,
     >                   NX3, NYE, IC1, IC2, DX, DY, NCMN, NCSD )
      DIMENSION PHMN(1), PHSD(1)
      CHARACTER*(*) JOB, SUB1
      DATA ZERO/0.0/, ONE/1.0/, TWO/2.0/
C						PRESSURE HEAD STATISTICS
      IF( NS .GT. 1 ) THEN
         FNS = FLOAT(NS)
         DNS = ONE/FNS
         DN1 = ONE/FLOAT(NS-1)
      ELSE
         FNS = ONE
         DNS = ONE
         DN1 = ZERO
      ENDIF
 
      DO 10 I = 1, N
         A = PHMN(I)*DNS
         V = (PHSD(I) - TWO*A*PHMN(I) + FNS*A*A)*DN1
         PHMN(I) = A
         PHSD(I) = SQRT(V)
  10  CONTINUE
C						CREATE OUTPUT FILES
C							0-WALL CASE
      IF( NW .EQ. 0 ) THEN
         CALL PRPH0( JOB, SUB1, IPMNL, PHMN, NX1, NYE, DX, DY, NCMN,
     >               'Mean Nodal Pressures' )
         IF( NS .GT. 1 ) THEN
            CALL PRPH0( JOB, SUB1, IPSDL, PHSD, NX1, NYE, DX, DY,
     >                  NCSD, 'Nodal Pressure Standard Deviations' )
         ENDIF
C							1-WALL CASE
      ELSEIF( NW .EQ. 1 ) THEN
         CALL PRPH1( JOB, SUB1, IPMNL, IPMNC, PHMN, NX1, NX2, NYE, IC1,
     >               DX, DY, NCMN,
     >               'Mean Nodal Pressures' )
         IF( NS .GT. 1 ) THEN
            CALL PRPH1( JOB, SUB1, IPSDL, IPSDC, PHSD, NX1, NX2, NYE,
     >                  IC1, DX, DY, NCSD,
     >                  'Nodal Pressure Standard Deviations')
         ENDIF
C							2-WALL CASE
      ELSEIF( NW .EQ. 2 ) THEN
         CALL PRPH2( JOB, SUB1, IPMNL, IPMNC, IPMNR, PHMN, NX1, NX2,
     >               NX3, NYE, IC1, IC2, DX, DY, NCMN,
     >               'Mean Nodal Pressures')
         IF( NS .GT. 1 ) THEN
            CALL PRPH2( JOB, SUB1, IPSDL, IPSDC, IPSDR, PHSD, NX1, NX2,
     >                  NX3, NYE, IC1, IC2, DX, DY, NCSD,
     >                  'Nodal Pressure Standard Deviations')
         ENDIF
C							UNDEFINED CASE
      ELSE
         CALL ERROR( 6, NW, ITERM )
      ENDIF
 
      RETURN
      END
C  ********************************************************************
C  *                                                                  *
C  *                       SUBROUTINE PLTFLD                          *
C  *                                                                  *
C  ********************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  TO OUTPUT AN ARRAY OF DATA TO A FILE IN DISPLAY FORMAT (FOR
C           DISPLAY ON A POSTSCRIPT PRINTER).
C
C  THIS ROUTINE TAKES AN ARRAY (2-D) OF DATA AND DUMPS IT TO A FILE
C  HAVING A FORMAT READABLE BY DISPLAY. ARGUMENTS ARE AS FOLLOWS;
C
C    JOB    CHARACTER STRING CONTAINING THE TITLE OF THE RUN. (INPUT)
C
C   SUB1    CHARACTER STRING CONTAINING THE SUBTITLE OF THE RUN. (INPUT)
C
C      Z    REAL ARRAY OF SIZE AT LEAST N1 X N2 CONTAINING THE DATA TO
C           BE DISPLAYED GRAPHICALLY. (INPUT)
C
C     IZ    LEADING DIMENSION OF Z EXACTLY AS SPECIFIED IN THE CALLING
C           ROUTINE. (INPUT)
C
C     N1    COLUMN (1ST INDEX) DIMENSION OF Z. (INPUT)
C
C     N2    ROW (2ND INDEX) DIMENSION OF Z. (INPUT)
C
C     D1    PHYSICAL DIMENSION OF THE DATA IN THE X (1ST INDEX) DIRECTION.
C           (INPUT)
C
C     D2    PHYSICAL DIMENSION OF THE DATA IN THE Y (2ND INDEX) DIRECTION.
C           (INPUT)
C
C  TITLE    CHARACTER STRING CONTAINING THE TITLE OF THE DISPLAY. (INPUT)
C
C   IFLD    UNIT NUMBER TO WHICH OUTPUT IS SENT. (INPUT)
C
C-----------------------------------------------------------------------------
      SUBROUTINE PLTFLD( JOB, SUB1, Z, IZ, N1, N2, D1, D2, TITLE, IFLD )
      DIMENSION Z(IZ,1)
      CHARACTER*(*) TITLE
 
   1  FORMAT(A)
   2  FORMAT(2E13.5)
   3  FORMAT()
   4  FORMAT(2I8)
   5  FORMAT(8E12.4)
 
      WRITE(IFLD,1) JOB
      WRITE(IFLD,1) SUB1
      WRITE(IFLD,1) TITLE
      WRITE(IFLD,3)
      WRITE(IFLD,1) '2'
      WRITE(IFLD,1) 'x'
      WRITE(IFLD,1) 'y'
      WRITE(IFLD,1) 'z'
      WRITE(IFLD,2) 0., D1
      WRITE(IFLD,2) 0., D2
      WRITE(IFLD,4) N1, N2
      WRITE(IFLD,3)
      WRITE(IFLD,5) ((Z(I,J), I = 1, N1), J = 1, N2)
      CLOSE(IFLD)
 
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                      SUBROUTINE PRPH0                           *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY D. VAUGHAN GRIFFITHS AND GORDON A. FENTON, 1991
C
C  PURPOSE  MAPS NODAL LOADS TO PHYSICAL NODE LOCATIONS FOR THE SINGLE
C           WALL PROBLEM AND DUMPS THE DATA TO A FILE IN DISPLAY FORMAT.
C
C  ARGUMENTS TO THE ROUTINE ARE AS FOLLOWS;
C
C    JOB    CHARACTER STRING CONTAINING THE TITLE OF THE RUN. (INPUT)
C
C    SUB1   CHARACTER STRING CONTAINING THE SUBTITLE OF THE RUN. (INPUT)
C
C    IOUT   NODAL LOADS ARE PRINTED TO UNIT NUMBER IOUT. (INPUT)
C
C    LOADS  THE NODAL DATA TO BE MAPPED TO THE PHYSICAL NODE LOCATIONS
C           AND PRINTED. (INPUT)
C
C    NX1    THE NUMBER OF ELEMENTS IN THE X DIRECTION TO THE LEFT OF
C           THE WALL. (INPUT)
C
C    NYE    THE TOTAL NUMBER OF ELEMENTS IN THE Y DIRECTION. (INPUT)
C
C     DX    PHYSICAL DIMENSION OF EACH FINITE ELEMENT IN THE X DIRECTION.
C           (INPUT)
C
C     DY    PHYSICAL DIMENSION OF EACH FINITE ELEMENT IN THE Y DIRECTION.
C           (INPUT)
C
C    NCNT   THE NUMBER OF CONTOURING LINES TO SPECIFY IN THE DISPLAY FILE.
C           NCNT MAY NOT EXCEED 20. (INPUT)
C
C    TITLE  CHARACTER STRING CONTAINING THE MAIN TITLE OF THE DATA. (INPUT)
C--------------------------------------------------------------------------
      SUBROUTINE PRPH0( JOB, SUB1, IOUT, LOADS, NX1, NYE, DX, DY,
     >                  NCNT, TITLE )
      REAL LOADS(*), CV(20), V(6)
      CHARACTER*(*) JOB, SUB1, TITLE
      DATA ZERO/0.0/, PT2/0.2/
 
   1  FORMAT(A)
   2  FORMAT(6E12.4)
   3  FORMAT(4I5)
 
      NY  = NYE + 1
      N1  = NX1 + 1
      N   = NX1*NY
      SXL = DX*FLOAT(NX1)
      SY  = DY*FLOAT(NYE)
 
      WRITE(IOUT,1) JOB
      WRITE(IOUT,1) SUB1
      WRITE(IOUT,1) TITLE
      WRITE(IOUT,1) ' '
      WRITE(IOUT,1) '2'
      WRITE(IOUT,1) 'X'
      WRITE(IOUT,1) 'Y'
      WRITE(IOUT,1) ' '
      WRITE(IOUT,2) ZERO, SXL
      WRITE(IOUT,2) ZERO, SY
      WRITE(IOUT,3) N1, NY
C					FIND RANGE IN DATA
      DMIN = ZERO
      DMAX = LOADS(1)
      DO 10 I = 2, N
         IF( LOADS(I) .LT. DMIN ) DMIN = LOADS(I)
         IF( LOADS(I) .GT. DMAX ) DMAX = LOADS(I)
  10  CONTINUE
C					PICK SOME CONTOURS
      IF( NCNT .GT. 20 ) THEN
         JCNT = 20
      ELSEIF( NCNT .LT. 1 ) THEN
         JCNT = 1
      ELSE
         JCNT = NCNT
      ENDIF
 
      R = (DMAX - DMIN)/FLOAT(JCNT+1)
      CV(1) = DMIN + R
      DO 20 I = 2, JCNT
         CV(I) = CV(I-1) + R
  20  CONTINUE
      WRITE(IOUT,3) JCNT
      WRITE(IOUT,2) (CV(I),I=1,JCNT)
C					PRINT THE DATA
C	NOTE; LOADS IS STORED SEQUENTIALLY IN THE Y DIRECTION (INCREASING
C             DOWNWARDS), DISPLAY WANTS IT STORED IN THE X DIRECTION WITH
C             Y INCREASING UPWARDS.
 
      IV = 0
      DO 40 J = NY, 1, -1
         DO 30 I = J, N, NY
            IV = IV + 1
            V(IV) = LOADS(I)
            IF( IV .EQ. 6 ) THEN
               WRITE(IOUT,2) (V(K), K = 1, 6)
               IV = 0
            ENDIF
  30     CONTINUE
         IV = IV + 1
         V(IV) = ZERO
         IF( IV .EQ. 6 ) THEN
            WRITE(IOUT,2) (V(K), K = 1, 6)
            IV = 0
         ENDIF
  40  CONTINUE
 
      IF( IV .GT. 0 ) WRITE(IOUT,2) (V(K), K = 1, IV)
C						EXTRA DATA FOR DISPLAY
      WRITE(IOUT,2) PT2, ZERO, SY, SXL, SY
      WRITE(IOUT,2) PT2, SXL, SY, SXL, ZERO
 
      CLOSE(IOUT)
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                      SUBROUTINE PRPH1                           *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY D. VAUGHAN GRIFFITHS AND GORDON A. FENTON, 1991
C
C  PURPOSE  MAPS NODAL LOADS TO PHYSICAL NODE LOCATIONS FOR THE SINGLE
C           WALL PROBLEM AND DUMPS THE DATA TO A FILE IN DISPLAY FORMAT.
C
C  ARGUMENTS TO THE ROUTINE ARE AS FOLLOWS;
C
C    JOB    CHARACTER STRING CONTAINING THE TITLE OF THE RUN. (INPUT)
C
C   SUB1    CHARACTER STRING CONTAINING THE SUBTITLE OF THE RUN. (INPUT)
C
C    IPL    UNIT NUMBER TO WHICH THE NODAL DATA TO THE LEFT OF THE WALL
C           IS PRINTED. (INPUT)
C
C    IPR    UNIT NUMBER TO WHICH THE NODAL DATA TO THE RIGHT OF THE WALL
C           IS PRINTED. (INPUT)
C
C  LOADS    THE NODAL DATA TO BE MAPPED TO THE PHYSICAL NODE LOCATIONS
C           AND PRINTED. (INPUT)
C
C    NX1    THE NUMBER OF ELEMENTS IN THE X DIRECTION TO THE LEFT OF
C           THE WALL. (INPUT)
C
C    NX2    THE NUMBER OF ELEMENTS IN THE Y DIRECTION TO THE RIGHT OF
C           THE WALL. (INPUT)
C
C    NYE    THE TOTAL NUMBER OF ELEMENTS IN THE Y DIRECTION. (INPUT)
C
C    IC1    THE NUMBER OF ELEMENTS IN THE Y DIRECTION AGAINST THE WALL.
C           (INPUT)
C
C    SXL    PHYSICAL DIMENSION OF THE FINITE ELEMENT MODEL IN THE X DIRECTION
C           TO THE LEFT OF THE WALL. (INPUT)
C
C     DX    PHYSICAL DIMENSION OF EACH FINITE ELEMENT IN THE X DIRECTION
C           (INPUT)
C
C     DY    PHYSICAL DIMENSION OF EACH FINITE ELEMENT IN THE Y DIRECTION.
C           (INPUT)
C
C   NCNT    THE NUMBER OF CONTOURING LINES TO SPECIFY IN THE DISPLAY FILE.
C           NCNT MAY NOT EXCEED 20. (INPUT)
C
C    TITLE  CHARACTER STRING CONTAINING THE MAIN TITLE OF THE DATA. (INPUT)
C--------------------------------------------------------------------------
      SUBROUTINE PRPH1( JOB, SUB1, IPL, IPR, LOADS, NX1, NX2, NYE,
     >                  IC1, DX, DY, NCNT, TITLE )
      REAL LOADS(*), CV(20), V(6)
      CHARACTER*(*) JOB, SUB1, TITLE
      DATA ZERO/0.0/, PT2/0.2/, TWO/2.0/
 
   1  FORMAT(A,A)
   2  FORMAT(6E12.4)
   3  FORMAT(4I5)
 
      NY  = NYE + 1
      N1  = NX1 + 1
      N2  = NX2 + 1
      NC  = N1*NY
      N   = N1*NY + NX2*NYE + IC1 - 1
      SXL = DX*FLOAT(NX1)
      SXR = SXL + DX*FLOAT(NX2)
      SY  = DY*FLOAT(NYE)
      WL1 = SY - DY*FLOAT(IC1)
 
      WRITE(IPL,1) JOB
      WRITE(IPL,1) SUB1
      JT = LNBLNK(TITLE)
      WRITE(IPL,1) TITLE(1:JT)
      WRITE(IPL,1) ' '
      WRITE(IPL,1) '2'
      WRITE(IPL,1) 'X'
      WRITE(IPL,1) 'Y'
      WRITE(IPL,1) ' '
      WRITE(IPL,2) ZERO, SXL
      WRITE(IPL,2) ZERO, SY
      WRITE(IPL,3) N1, NY
C					FIND RANGE IN DATA
      DMIN = ZERO
      DMAX = LOADS(1)
      DO 10 I = 2, N
         IF( LOADS(I) .LT. DMIN ) DMIN = LOADS(I)
         IF( LOADS(I) .GT. DMAX ) DMAX = LOADS(I)
  10  CONTINUE
C					PICK SOME CONTOURS
      IF( NCNT .GT. 20 ) THEN
         JCNT = 20
      ELSEIF( NCNT .LT. 1 ) THEN
         JCNT = 1
      ELSE
         JCNT = NCNT
      ENDIF
 
      R = (DMAX - DMIN)/FLOAT(JCNT+1)
      CV(1) = DMIN + R
      DO 20 I = 2, JCNT
         CV(I) = CV(I-1) + R
  20  CONTINUE
      WRITE(IPL,3) JCNT
      WRITE(IPL,2) (CV(I),I=1,JCNT)
 
C					PRINT THE DATA - LEFT SIDE FIRST
C	NOTE; LOADS IS STORED SEQUENTIALLY IN THE Y DIRECTION (INCREASING
C             DOWNWARDS), DISPLAY WANTS IT STORED IN THE X DIRECTION WITH
C             Y INCREASING UPWARDS.
 
      IV = 0
      DO 30 J = NY, 1, -1
      DO 30 I = J, NC, NY
         IV = IV + 1
         V(IV) = LOADS(I)
         IF( IV .EQ. 6 ) THEN
            WRITE(IPL,2) ( V(K), K = 1, 6)
            IV = 0
         ENDIF
  30  CONTINUE
      IF( IV .GT. 0 ) WRITE(IPL,2) (V(K), K = 1, IV)
      WRITE(IPL,2) PT2, ZERO, SY, SXL, SY
      WRITE(IPL,2) TWO, SXL, SY, SXL, WL1
      CLOSE(IPL)
 
C					NOW FOR THE RIGHT HAND SIDE ...
      WRITE(IPR,1) JOB
      WRITE(IPR,1) SUB1
      WRITE(IPR,1) TITLE(1:JT),' (right side)'
      WRITE(IPR,1) ' '
      WRITE(IPR,1) '2'
      WRITE(IPR,1) 'X'
      WRITE(IPR,1) 'Y'
      WRITE(IPR,1) ' '
      WRITE(IPR,2) SXL, SXR
      WRITE(IPR,2) ZERO, SY
      WRITE(IPR,3) N2, NY
C					USE THE SAME CONTOURS
      WRITE(IPR,3) JCNT
      WRITE(IPR,2) (CV(I),I=1,JCNT)
C					EXTRACT AND PRINT THE DATA
      IV = 0
      DO 40 J = NY, 2, -1
         IV = IV + 1
         IF( J .LE. IC1 ) THEN
            V(IV) = LOADS(NC+J-1)
         ELSE
            V(IV) = LOADS(NX1*NY+J)
         ENDIF
         IF( IV .EQ. 6 ) THEN
            WRITE(IPR,2) (V(K), K = 1, 6)
            IV = 0
         ENDIF
         NJ = NC + IC1 + J - 2
         DO 40 I = 2, N2
            IV = IV + 1
            V(IV) = LOADS(NJ+(I-2)*NYE)
            IF( IV .EQ. 6 ) THEN
               WRITE(IPR,2) (V(K), K = 1, 6)
               IV = 0
            ENDIF
  40  CONTINUE
      DO 50 I = 1, N2
         IV = IV + 1
         V(IV) = ZERO
         IF( IV .EQ. 6 ) THEN
            WRITE(IPR,2) ( V(K), K = 1, 6 )
            IV = 0
         ENDIF
  50  CONTINUE
      IF( IV .GT. 0 ) WRITE(IPR,2) (V(K), K = 1, IV)
      WRITE(IPR,2) PT2, SXL, SY, SXR, SY
      WRITE(IPR,2) PT2, SXR, SY, SXR, ZERO
      CLOSE(IPR)
 
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                      SUBROUTINE PRPH2                           *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY D. VAUGHAN GRIFFITHS AND GORDON A. FENTON, 1991
C
C  PURPOSE  MAPS NODAL LOADS TO PHYSICAL NODE LOCATIONS FOR THE DOUBLE
C           WALL PROBLEM AND DUMPS THE DATA TO A FILE IN DISPLAY FORMAT.
C
C  ARGUMENTS TO THE ROUTINE ARE AS FOLLOWS;
C
C    JOB    CHARACTER STRING CONTAINING THE TITLE OF THE RUN. (INPUT)
C
C   SUB1    CHARACTER STRING CONTAINING THE SUBTITLE OF THE RUN. (INPUT)
C
C    IPL    UNIT NUMBER TO WHICH THE NODAL DATA TO THE LEFT OF THE LEFT WALL
C           IS PRINTED. (INPUT)
C
C    IPC    UNIT NUMBER TO WHICH THE NODAL DATA BETWEEN WALLS IS PRINTED.
C           (INPUT)
C
C    IPR    UNIT NUMBER TO WHICH THE NODAL DATA TO THE RIGHT OF THE RIGHT WALL
C           IS PRINTED. (INPUT)
C
C  LOADS    THE NODAL DATA TO BE MAPPED TO THE PHYSICAL NODE LOCATIONS
C           AND PRINTED. (INPUT)
C
C    NX1    THE NUMBER OF ELEMENTS IN THE X DIRECTION TO THE LEFT OF
C           THE LEFT WALL. (INPUT)
C
C    NX2    THE NUMBER OF ELEMENTS IN THE Y DIRECTION BETWEEN WALLS. (INPUT)
C
C    NX3    THE NUMBER OF ELEMENTS IN THE X DIRECTION TO THE RIGHT OF
C           THE RIGHT WALL. (INPUT)
C
C    NYE    THE TOTAL NUMBER OF ELEMENTS IN THE Y DIRECTION. (INPUT)
C
C    IC1    THE NUMBER OF ELEMENTS IN THE Y DIRECTION AGAINST THE LEFT WALL.
C           (INPUT)
C
C    IC2    THE NUMBER OF ELEMENTS IN THE Y DIRECTION AGAINST THE RIGHT WALL.
C           (INPUT)
C
C     DX    PHYSICAL DIMENSION OF EACH FINITE ELEMENT IN THE X DIRECTION.
C           (INPUT)
C
C     DY    PHYSICAL DIMENSION OF EACH FINITE ELEMENT IN THE Y DIRECTION.
C           (INPUT)
C
C   NCNT    THE NUMBER OF CONTOURING LINES TO SPECIFY IN THE DISPLAY FILE.
C           NCNT MAY NOT EXCEED 20. (INPUT)
C
C  TITLE    CHARACTER STRING CONTAINING THE MAIN TITLE OF THE DATA. (INPUT)
C--------------------------------------------------------------------------
      SUBROUTINE PRPH2( JOB, SUB1, IPL, IPC, IPR, LOADS, NX1, NX2, NX3,
     >                  NYE, IC1, IC2, DX, DY, NCNT, TITLE )
      REAL LOADS(*), CV(20), V(6)
      CHARACTER*(*) JOB, SUB1, TITLE
      DATA ZERO/0.0/, PT2/0.2/, TWO/2.0/
 
   1  FORMAT(A,A)
   2  FORMAT(6E12.4)
   3  FORMAT(4I5)
 
      NY  = NYE + 1
      N1  = NX1 + 1
      N2  = NX2 + 1
      N3  = NX3 + 1
      NC  = N1*NY
      ND  = NC + IC1 + NX2*NY
      N   = ND + NX3*NYE + IC2 - 1
      SXL = DX*FLOAT(NX1)
      SXC = SXL + DX*FLOAT(NX2)
      SXR = SXC + DX*FLOAT(NX3)
      SY  = DY*FLOAT(NYE)
      WL1 = SY - DY*FLOAT(IC1)
      WL2 = SY - DY*FLOAT(IC2)
 
      WRITE(IPL,1) JOB
      WRITE(IPL,1) SUB1
      JT = LNBLNK(TITLE)
      WRITE(IPL,1) TITLE(1:JT)
      WRITE(IPL,1) ' '
      WRITE(IPL,1) '2'
      WRITE(IPL,1) 'X'
      WRITE(IPL,1) 'Y'
      WRITE(IPL,1) ' '
      WRITE(IPL,2) ZERO, SXL
      WRITE(IPL,2) ZERO, SY
      WRITE(IPL,3) N1, NY
C					FIND RANGE IN DATA
      DMIN = ZERO
      DMAX = LOADS(1)
      DO 10 I = 2, N
         IF( LOADS(I) .LT. DMIN ) DMIN = LOADS(I)
         IF( LOADS(I) .GT. DMAX ) DMAX = LOADS(I)
  10  CONTINUE
C					PICK SOME CONTOURS (TRY 6)
      IF( NCNT .GT. 20 ) THEN
         JCNT = 20
      ELSEIF( NCNT .LT. 1 ) THEN
         JCNT = 1
      ELSE
         JCNT = NCNT
      ENDIF
 
      R = (DMAX - DMIN)/FLOAT(JCNT+1)
      CV(1) = DMIN + R
      DO 20 I = 2, JCNT
         CV(I) = CV(I-1) + R
  20  CONTINUE
      WRITE(IPL,3) JCNT
      WRITE(IPL,2) (CV(I),I=1,JCNT)
C					PRINT THE DATA - LEFT SIDE FIRST
C	NOTE; LOADS IS STORED SEQUENTIALLY IN THE Y DIRECTION (INCREASING
C             DOWNWARDS), DISPLAY WANTS IT STORED IN THE X DIRECTION WITH
C             Y INCREASING UPWARDS.
      IV = 0
      DO 30 J = NY, 1, -1
      DO 30 I = J, NC, NY
         IV = IV + 1
         V(IV) = LOADS(I)
         IF( IV .EQ. 6 ) THEN
            WRITE(IPL,2) ( V(K), K = 1, 6)
            IV = 0
         ENDIF
  30  CONTINUE
      IF( IV .GT. 0 ) WRITE(IPL,2) (V(K), K = 1, IV)
      WRITE(IPL,2) PT2, ZERO, SY, SXL, SY
      WRITE(IPL,2) TWO, SXL, SY, SXL, WL1
      CLOSE(IPL)
 
C					NOW FOR THE CENTER BLOCK ...
      WRITE(IPC,1) JOB
      WRITE(IPC,1) SUB1
      WRITE(IPC,1) TITLE(1:JT),' (center)'
      WRITE(IPC,1) ' '
      WRITE(IPC,1) '2'
      WRITE(IPC,1) 'X'
      WRITE(IPC,1) 'Y'
      WRITE(IPC,1) ' '
      WRITE(IPC,2) SXL, SXC
      WRITE(IPC,2) ZERO, SY
      WRITE(IPC,3) N2, NY
C					USE THE SAME CONTOURS
      WRITE(IPC,3) JCNT
      WRITE(IPC,2) (CV(I),I=1,JCNT)
C					EXTRACT AND PRINT THE DATA
      IV = 0
      DO 40 J = NY, 1, -1
         IV = IV + 1
         IF( J .LE. IC1 ) THEN
            V(IV) = LOADS(NC+J)
         ELSE
            V(IV) = LOADS(NX1*NY+J)
         ENDIF
         IF( IV .EQ. 6 ) THEN
            WRITE(IPC,2) (V(K), K = 1, 6)
            IV = 0
         ENDIF
         NJ = NC + IC1 + J
         DO 40 I = 2, N2
            IV = IV + 1
            V(IV) = LOADS(NJ+(I-2)*NY)
            IF( IV .EQ. 6 ) THEN
               WRITE(IPC,2) (V(K), K = 1, 6)
               IV = 0
            ENDIF
  40  CONTINUE
      IF( IV .GT. 0 ) WRITE(IPC,2) (V(K), K = 1, IV)
      WRITE(IPC,2) PT2, SXL, SY, SXC, SY
      WRITE(IPC,2) TWO, SXC, SY, SXC, WL2
      CLOSE(IPC)
 
C					NOW FOR THE RIGHT HAND SIDE ...
      WRITE(IPR,1) JOB
      WRITE(IPR,1) SUB1
      WRITE(IPR,1) TITLE(1:JT),' (right side)'
      WRITE(IPR,1) ' '
      WRITE(IPR,1) '2'
      WRITE(IPR,1) 'X'
      WRITE(IPR,1) 'Y'
      WRITE(IPR,1) ' '
      WRITE(IPR,2) SXC, SXR
      WRITE(IPR,2) ZERO, SY
      WRITE(IPR,3) N3, NY
C					USE THE SAME CONTOURS
      WRITE(IPR,3) JCNT
      WRITE(IPR,2) (CV(I),I=1,JCNT)
C					EXTRACT AND PRINT THE DATA
      IV = 0
      ICC = NC + IC1 + (NX2-1)*NY
      DO 50 J = NY, 2, -1
         IV = IV + 1
         IF( J .LE. IC2 ) THEN
            V(IV) = LOADS(ND+J-1)
         ELSE
            V(IV) = LOADS(ICC+J)
         ENDIF
         IF( IV .EQ. 6 ) THEN
            WRITE(IPR,2) (V(K), K = 1, 6)
            IV = 0
         ENDIF
         NJ = ND + IC2 + J - 2
         DO 50 I = 2, N3
            IV = IV + 1
            V(IV) = LOADS(NJ+(I-2)*NYE)
            IF( IV .EQ. 6 ) THEN
               WRITE(IPR,2) (V(K), K = 1, 6)
               IV = 0
            ENDIF
  50  CONTINUE
      DO 60 I = 1, N3
         IV = IV + 1
         V(IV) = ZERO
         IF( IV .EQ. 6 ) THEN
            WRITE(IPR,2) ( V(K), K = 1, 6 )
            IV = 0
         ENDIF
  60  CONTINUE
      IF( IV .GT. 0 ) WRITE(IPR,2) (V(K), K = 1, IV)
      WRITE(IPR,2) PT2, SXC, SY, SXR, SY
      WRITE(IPR,2) PT2, SXR, SY, SXR, ZERO
      CLOSE(IPR)
 
      RETURN
      END
C  *********************************************************************
C  *                                                                   *
C  *                          FUNCTION SECOND                          *
C  *                                                                   *
C  *********************************************************************
C  SINGLE PRECISION VERSION 3.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1990
C
C  PURPOSE RETURNS ELAPSED USER EXECUTION TIME IN SECONDS.
C
C  RETURNS ELAPSED USER EXECUTION TIME OF THE CALLING PROCESS IN SECONDS.
C  THIS ROUTINE TENDS TO BE SYSEM SPECIFIC AND YOU MAY NEED TO CUSTOMIZE IT
C  FOR YOUR ENVIRONMENT. THIS ROUTINE WORKS FOR SUN'S AND VAX'S RUNNING
C  ULTRIX F77.
C---------------------------------------------------------------------------
      REAL FUNCTION SECOND()
C				(FOR THE VECTOR AMDAHL)
C     CALL CLOCK( TYME, 0, 2 )
C     SECOND = TYME
C				(FOR THE SCALAR AMDAHL)
      CALL CLOCKX( TYME )
      SECOND = 1.E-6*TYME
 
      RETURN
      END
      SUBROUTINE SEEP4(AA,BB,PX,PY,KP)
C
C  PURPOSE  THIS SUBROUTINE PRODUCES THE 'STIFFNESS' MATRIX FOR LAPLACE'S
C           EQUATION FOR RECTANGULAR ELEMENTS
C
      REAL KP(4,4)
      FAC1=(BB*BB*PX+AA*AA*PY)/(AA*BB)
      FAC2=-(-BB*BB*PX+2.*AA*AA*PY)/(6.*AA*BB)
      FAC3= (-2.*BB*BB*PX+AA*AA*PY)/(6.*AA*BB)
      KP(1,1)=FAC1/3.
      KP(1,2)=FAC2
      KP(1,3)=-FAC1/6.
      KP(1,4)=FAC3
      KP(2,2)=FAC1/3.
      KP(2,3)=FAC3
      KP(2,4)=-FAC1/6.
      KP(3,3)=FAC1/3.
      KP(3,4)=FAC2
      KP(4,4)=FAC1/3.
 
      KP(2,1) = KP(1,2)
      KP(3,1) = KP(1,3)
      KP(4,1) = KP(1,4)
      KP(3,2) = KP(2,3)
      KP(4,2) = KP(2,4)
      KP(4,3) = KP(3,4)
 
      RETURN
      END
C  ********************************************************************
C  *                                                                  *
C  *                        SUBROUTINE SIM2DK                         *
C  *                                                                  *
C  ********************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  TO GENERATE A REALIZATION OF A SOIL PERMEABILITY FIELD
C
C  THIS ROUTINE CREATES A REALIZATION OF A SOIL PERMEABILITY FIELD WHERE
C  THE PERMEABILITY HAS TWO COMPONENTS: X AND Y DIRECTIONS. THE PERMEABILITIES
C  IN THE TWO DIRECTIONS ARE CORRELATED ACCORDING TO RHOXY AND CORRELATION
C  OF THE FIELD IS PERFORMED ON A POINT BY POINT BASIS. THE GENERATION
C  PROCEEDS AS FOLLOWS;
C
C   1) GENERATE TWO INDEPENDANT ZERO-MEAN, UNIT-VARIANCE, GAUSSIAN RANDOM
C      FIELDS,
C   2) MAP THE GENERATED FIELDS INTO THE REQUIRED PERMEABILITY MATRICES
C      (ACCORDING TO THE PHYSICAL POSITION OF EACH ELEMENT)
C   3) LINEARLY COMBINE TO FORM A PAIR OF JOINTLY CORRELATED RANDOM FIELDS
C   4) TRANSFORM THE PERMEABILITIES INTO THEIR APPROPRIATE DISTRIBUTIONS
C      (IE, LOG-NORMAL).
C
C  ARGUMENTS ARE DESCRIBED AS FOLLOWS;
C
C    PERMX   REAL ARRAY OF SIZE AT LEAST NXE X NYE WHICH ON OUTPUT WILL
C            CONTAIN THE DESIRED REALIZATION (X-COMPONENT OF PERMEABILITY).
C            (OUTPUT)
C
C    PERMY   REAL ARRAY OF SIZE AT LEAST NXE X NYE WHICH ON OUTPUT WILL
C            CONTAIN THE DESIRED REALIZATION (Y-COMPONENT OF PERMEABILITY).
C            (OUTPUT)
C
C    IXY     LEADING DIMENSION OF PERMX AND PERMY AS SPECIFIED IN THE
C            CALLING ROUTINE. (INPUT)
C
C    Z1,Z2   TWO TEMPORARY REAL VECTORS USED TO STORE THE INTERMEDIATE
C            RANDOM FIELDS.
C
C    MZ      SIZE OF THE VECTORS Z1 AND Z2 AS ALLOCATED IN THE CALLING
C            ROUTINE. (INPUT)
C
C    NXE     NUMBER OF ELEMENTS IN THE X DIRECTION. (INPUT)
C
C    NYE     NUMBER OF ELEMENTS IN THE Y DIRECTION. (INPUT)
C
C    SX      PHYSICAL SIZE OF THE SITE MODEL IN THE X DIRECTION. (INPUT)
C
C    SY      PHYSICAL SIZE OF THE SITE MODEL IN THE Y DIRECTION. (INPUT)
C
C    THX     SCALE OF FLUCTUATION OF LOG-PERMEABILITY IN THE X DIRECTION.
C            (INPUT)
C
C    THY     SCALE OF FLUCTUATION OF LOG-PERMEABILITY IN THE Y DIRECTION.
C            (INPUT)
C
C    KMNX    (REAL) MEAN PERMEABILITY IN THE X DIRECTION. (INPUT)
C
C    KMNY    (REAL) MEAN PERMEABILITY IN THE Y DIRECTION. (INPUT)
C
C    KSDX    (REAL) STANDARD DEVIATION OF PERMEABILITY IN THE X DIRECTION.
C            KSDX IS RETURNED AS THE VALUE OF STANDARD DEVIATION EXPECTED
C            THEORETICALLY UNDER LOCAL AVERAGING. (INPUT/OUTPUT)
C
C    KSDY    (REAL) STANDARD DEVIATION OF PERMEABILITY IN THE Y DIRECTION.
C            KSDY IS RETURNED AS THE VALUE OF STANDARD DEVIATION EXPECTED
C            THEORETICALLY UNDER LOCAL AVERAGING. (INPUT/OUTPUT)
C
C    RHOXY   CORRELATION COEFFICIENT BETWEEN X AND Y DIRECTION
C            LOG-PERMEABILITIES. (INPUT)
C
C    INIT    INTEGER FLAG WHICH MUST BE 1 WHEN PARAMETERS OF THE RANDOM
C            FIELD ARE TO BE CALCULATED (OR RECALCULATED). IF MULTIPLE
C            REALIZATIONS OF THE SAME VECTOR FIELD ARE DESIRED, THEN
C            SUBSEQUENT CALLS TO THIS ROUTINE SHOULD USE INIT NOT
C            EQUAL TO 1. (INPUT)
C
C    KSEED   INTEGER SEED TO BE USED FOR THE PSEUDO-RANDOM NUMBER GENERATOR.
C            IF KSEED = 0, THEN A RANDOM SEED WILL BE USED (BASED ON THE
C            CLOCK TIME WHEN LAS2D IS CALLED FOR THE FIRST TIME).
C            ON OUTPUT, KSEED IS SET TO THE VALUE OF THE ACTUAL SEED USED.
C
C    IERR    INTEGER FLAG WHICH DENOTES THE STATUS OF THE RUN;
C            =  0 FOR SUCCESSFUL SIMULATION
C            = -1 IF THE TEMPORARY VECTORS Z1 AND Z2 ARE TOO SMALL
C            = -2 IF THE RANDOM FIELD TO MESH MAPPING FAILS
C            = -3 IF THE REQUIRED PROCESS IS TOO LARGE FOR LAS2D TO HANDLE
C            = -4 IF LAS2D RUNS INTO IMAGINARY SQRT'S
C
C    IFLD    UNIT NUMBER TO WHICH THE FIRST GENERATED RANDOM FIELD IS DUMPED
C            (IN A DISPLAY FORMAT). ONLY ONE FIELD IS PRINTED, SO IFLD IS
C            SET TO ZERO AFTER DUMPING. (INPUT/OUTPUT)
C
C    IMDL    UNIT NUMBER TO WHICH THE FIRST SITE MODEL (MAPPED FROM THE
C            FIRST RANDOM FIELD) OF LOG-PERMEABILITIES IS DUMPED. ONLY ONE
C            FIELD IS PRINTED, SO IMDL IS SET TO ZERO AFTER DUMPING.
C            (INPUT/OUTPUT)
C
C    JOB     CHARACTER STRING CONTAINING THE TITLE OF THE RUN. (INPUT)
C
C    SUB1    CHARACTER STRING CONTAINING THE SUBTITLE OF THE RUN. (INPUT)
C
C  NOTE: STRICTLY SPEAKING THERE SHOULD BE X AND Y SCALES OF FLUCTUATION
C        FOR EACH OF THE X AND Y DIRECTION PERMEABILITIES - I ASSUME THEY
C        ARE EQUIVALENT IN THIS ROUTINE (IE THE X DIRECTION PERMEABILITY
C        HAS THE SAME SCALES OF FLUCTUATION - IN THE X AND Y DIRECTIONS -
C        AS DOES THE Y DIRECTION PERMEABILITY). THE X AND Y DIRECTION
C        PERMEABILITIES MAY HAVE DIFFERING MEANS AND VARIANCES, HOWEVER.
C
C==================================
C PARAMETERS:
C
C    TOL     REAL VALUE DENOTING THE MAXIMUM PERCENTAGE OF ANISOTROPY IN
C            THE FIELD GENERATED BY LAS2D.
C---------------------------------------------------------------------------
      SUBROUTINE SIM2DK(PERMX,PERMY,IXY,Z1,Z2,MZ,NXE,NYE,SX,SY,THX,THY,
     >                  KMNX,KMNY,KSDX,KSDY,RHOXY,INIT,KSEED,IERR,
     >                  IFLD,IMDL, JOB, SUB1)
      PARAMETER ( TOL = 0.05 )
      DIMENSION PERMX(IXY,1), PERMY(IXY,1), Z1(1), Z2(1)
C				EXPORT PARAMETERS TO VARIANCE FUNCTION
      REAL PA, PB, PX, PY, PZ, PIBY2, VARFNC, DX, DY
      REAL KMNX, KMNY, KSDX, KSDY
      LOGICAL DEBUG
      CHARACTER*(*) JOB, SUB1
      EXTERNAL VARFNC
      SAVE N, M, IX, IY, XL, YL, C, PMNX, PMNY, PSDX, PSDY
      COMMON/DPARAM/ PA, PB, PX, PY, PZ
      COMMON/DBGRFL/ ISTAT, DEBUG
      COMMON/DBGLAS/ ZMN, ZVR, ZCV, NN
      DATA ICNT/1/, ZERO/0.0/, HALF/0.5/, ONE/1.0/
      DATA PIBY2/1.570796326794896619231322/
 
   1  FORMAT(A)
   2  FORMAT(A,'(',I4,',',I4,')')
   3  FORMAT(A,I10,A,I10)
   4  FORMAT(A,E12.5,A,E12.5)
C					COMPUTE REQUIRED FIELD SIZE (ONCE)
      IF( INIT .EQ. 1 .OR. ICNT .EQ. 1 ) THEN
         ICNT = 0
         ZMN  = ZERO
         ZVR  = ZERO
         ZCV  = ZERO
C						ASSUME FAILURE (IERR = -2)
         IERR = -2
         IF( DEBUG ) WRITE(ISTAT,1)'SIM2DK: setting field parameters...'
 
C						SET VARIANCE FNC PARAMETERS
         PA = 1.D0
         PX = 1.D0
         PY = 1.D0
C						THIS ONE YOU MAY WANT TO CHANGE
C						(ITS SET FOR EXP. VARIANCE FNC)
         PB = PIBY2
 
         A = SY*THX*FLOAT(NXE)
         B = SX*THY*FLOAT(NYE)
         Q = A/B
         DO 10 IX = 1, 256
            D  = Q*FLOAT(IX)
            IY = INT(D + HALF)
            R  = ABS(D - FLOAT(IY))/FLOAT(IY)
            IF( R .LE. TOL ) GO TO 20
   10    CONTINUE
C						COULDN'T FIND A CELL MAPPING
         RETURN
C						FIND NUMBER OF SUBDIVISIONS
   20    NM = MAX0( IY*NYE, IX*NXE )
         IF( DEBUG ) WRITE(ISTAT,2)'   Field -> model cells (ix,iy ) = '
     >               ,IX, IY
         N  = 1
         DO 30 M = 0, 40
            IF( N .GE. NM ) GO TO 40
            N = 2*N
  30     CONTINUE
C						TOO MANY SUBDIVISIONS
         RETURN
C						CHECK MEMORY
  40     NN = N*N
         MEMREQ = NN + (4 + N/2)*(4 + N/2)
         IF( DEBUG ) THEN
            WRITE(ISTAT,2)'   required field size (N,N) = ', N, N
            WRITE(ISTAT,3)'   memreq = ',MEMREQ,',  memavail = ', MZ
         ENDIF
         IF( MEMREQ .GT. MZ ) THEN
            IERR = -1
            RETURN
         ENDIF
C						DETERMINE FIELD SIZE, ADMIT
C						SLIGHT ANISOTROPY (XL ~= YL)
         XL = FLOAT(N)*SX/(FLOAT(IX)*NXE*THX)
         YL = FLOAT(N)*SY/(FLOAT(IY)*NYE*THY)
C					FIND MEAN AND S.D. OF LOG-PERMEABILITY
         DX   = XL/FLOAT(N)
         DY   = YL/FLOAT(N)
         VF   = VARFNC(DX,DY)
         A1   = ONE + KSDX*KSDX/(KMNX*KMNX)
         A2   = ONE + KSDY*KSDY/(KMNY*KMNY)
         PVRX = VF*ALOG(A1)
         PMNX = ALOG(KMNX) - HALF*PVRX
         PSDX = SQRT(PVRX/VF)
         PVRY = VF*ALOG(A2)
         PMNY = ALOG(KMNY) - HALF*PVRY
         PSDY = SQRT(PVRY/VF)
C					SD OF LOCALLY AVERAGED PERMEABILITY
         KSDX = KMNX*SQRT( A1**VF - ONE )
         KSDY = KMNY*SQRT( A2**VF - ONE )
 
         IF( DEBUG ) THEN
            WRITE(ISTAT,4)'   Required field size = ',XL,' x ',YL
            WRITE(ISTAT,4)'             cell size = ',DX,' x ',DY
            WRITE(ISTAT,4)'   Variance reduction  = ',VF
            WRITE(ISTAT,4)'   Field target mean(X)= ',PMNX,
     >                    ', target var (X)= ',PVRX
            WRITE(ISTAT,4)'   Perm. target SD (X) = ',KSDX,
     >                    ', target SD (Y) = ',KSDY
            WRITE(ISTAT,1)'SIM2DK: parameter setup completed'
            WRITE(ISTAT,'()')
         ENDIF
C						CHOL. DECOMP OF CORR. MATRIX
         IF( RHOXY .LT. ONE ) C  = SQRT(ONE - RHOXY*RHOXY)
      ENDIF
C					NOW GENERATE THE ACTUAL FIELD(S)
      IERR = 0
      IF( RHOXY .EQ. ONE ) THEN
 
         CALL LAS2D( Z1, M, XL, YL, VARFNC, 2, KSEED, INIT,.FALSE.,IERR)
         IF( IERR .NE. 0 ) RETURN
         IF( IFLD .GT. 0 ) THEN
            CALL PLTFLD( JOB, SUB1, Z1, N, N, N, XL, YL,
     >                   'Isotropic Random Field', IFLD )
            IFLD = 0
         ENDIF
         IF( DEBUG ) THEN
            DO 45 I = 1, NN
               ZMN = ZMN + Z1(I)
               ZVR = ZVR + Z1(I)*Z1(I)
  45        CONTINUE
            ZCV = ZVR
         ENDIF
 
         CALL MAP2FE( Z1, N, PERMX, IXY, NXE, NYE, IX, IY )
         IF( IMDL .GT. 0 ) THEN
            CALL PLTFLD( JOB, SUB1, PERMX, IXY, NXE, NYE, SX, SY,
     >               'Log-Permeability Field', IMDL )
            IMDL = 0
         ENDIF
C						NORMAL -> LOG-NORMAL DIST'N
         DO 50 J = 1, NYE
         DO 50 I = 1, NXE
            PERMY(I,J) = EXP( PMNY + PSDY*PERMX(I,J) )
            PERMX(I,J) = EXP( PMNX + PSDX*PERMX(I,J) )
  50     CONTINUE
      ELSE
         CALL LAS2D( Z1, M, XL, YL, VARFNC, 2, KSEED, INIT,.FALSE.,IERR)
         IF( IERR .NE. 0 ) RETURN
         IF( IFLD .GT. 0 ) THEN
            CALL PLTFLD( JOB, SUB1, Z1, N, N, N, XL, YL,
     >                   'Isotropic Random Field', IFLD )
            IFLD = 0
         ENDIF
         CALL LAS2D( Z2, M, XL, YL, VARFNC, 2, KSEED, 0,.FALSE.,IERR)
         IF( IERR .NE. 0 ) RETURN
         IF( DEBUG ) THEN
            DO 55 I = 1, NN
               ZMN = ZMN + 0.5*(Z1(I) + Z2(I))
               ZVR = ZVR + 0.5*(Z1(I)*Z1(I) + Z2(I)*Z2(I))
               ZCV = ZCV + Z1(I)*Z2(I)
  55        CONTINUE
         ENDIF
 
         CALL MAP2FE( Z1, N, PERMX, IXY, NXE, NYE, IX, IY )
         IF( IMDL .GT. 0 ) THEN
            CALL PLTFLD( JOB, SUB1, PERMX, IXY, NXE, NYE, SX, SY,
     >               'Log-Permeability Field', IMDL )
            IMDL = 0
         ENDIF
         CALL MAP2FE( Z2, N, PERMY, IXY, NXE, NYE, IX, IY )
         DO 60 J = 1, NYE
         DO 60 I = 1, NXE
            TEMP       = RHOXY*PERMX(I,J) + C*PERMY(I,J)
            PERMX(I,J) = EXP( PMNX + PSDX*PERMX(I,J) )
            PERMY(I,J) = EXP( PMNY + PSDY*TEMP       )
  60     CONTINUE
      ENDIF
 
      RETURN
      END
C  *********************************************************************
C  *                                                                   *
C  *                        SUBROUTINE STATS                           *
C  *                                                                   *
C  *********************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  ESTIMATES THE MEAN AND VARIANCE OF THE INPUT DATA VECTOR AND
C           OUTPUTS THIS INFORMATION TO UNIT 6
C
C  ARGUMENTS ARE AS FOLLOWS;
C
C    ISTAT  UNIT NUMBER TO OUTPUT RUN-TIME STATS AND RESULTS. (INPUT)
C
C    ITERM  UNIT NUMBER CONNECTED TO THE USER'S TERMINAL. (INPUT)
C
C    JF     THE EXPECTED NUMBER OF SIMULATIONS PERFORMED. (INPUT)
C
C    NS     THE ACTUAL NUMBER OF SIMULATIONS PERFORMED (THIS MAY BE LESS THAN
C           JF SINCE SIMULATIONS WITH ERRORS ARE IGNORED). (INPUT)
C
C    NW     THE NUMBER OF WALLS IN THE PROBLEM (0, 1, OR 2). (INPUT)
C
C    TI     RUN TIME FOR THE RANDOM FIELD SET-UP PORTION OF THE ANALYSIS IN
C           SECONDS. (INPUT)
C
C    TS     RUN TIME FOR THE RANDOM FIELD SIMULATION PORTION OF THE ANALYSIS
C           IN SECONDS. (INPUT)
C
C    TF     RUN TIME FOR FINITE ELEMENT PORTION OF THE ANALYSIS IN SECONDS.
C           (INPUT)
C
C    KSEED  RANDOM NUMBER GENERATOR SEED USED FOR THIS RUN. (INPUT)
C
C    FLMN   ON INPUT FLMN CONTAINS THE FLOW RATE SUMMED OVER ALL THE
C           REALIZATIONS. ON OUTPUT FLMN WILL CONTAIN THE AVERAGE
C           FLOW RATE. (INPUT/OUTPUT)
C
C    FLSD   ON INPUT FLSD CONTAINS THE SQUARED FLOW RATES SUMMED OVER ALL
C           THE REALIZATIONS. ON OUTPUT, FLSD WILL CONTAIN THE STANDARD
C           DEVIATION OF THE FLOW RATE. (INPUT/OUTPUT)
C
C    EGMN   ON INPUT EGMN CONTAINS THE EXIT GRADIENTS SUMMED OVER ALL THE
C           REALIZATIONS. ON OUTPUT EGMN WILL CONTAIN THE AVERAGE
C           EXIT GRADIENT. (INPUT/OUTPUT)
C
C    EGSD   ON INPUT EGSD CONTAINS THE SQUARED EXIT GRADIENTS SUMMED OVER ALL
C           THE REALIZATIONS. ON OUTPUT, EGSD WILL CONTAIN THE STANDARD
C           DEVIATION OF THE EXIT GRADIENT. (INPUT/OUTPUT)
C
C    XKMN   ON INPUT, XKMN CONTAINS THE FIELD AVERAGE X-DIRECTION PERMEABILITY
C           SUMMED OVER ALL THE REALIZATIONS. ON OUTPUT, XKMN WILL CONTAIN
C           THE GLOBAL AVERAGE X-DIRECTION PERMEABILITY. (INPUT/OUTPUT)
C
C    XKSD   ON INPUT, XKSD CONTAINS THE FIELD AVERAGE X-DIRECTION SQUARED
C           PERMEABILITY SUMMED OVER ALL THE REALIZATIONS. ON OUTPUT, XKSD
C           WILL CONTAIN THE ESTIMATED GLOBAL STANDARD DEVIATION OF THE
C           X-DIRECTION PERMEABILITY. (INPUT/OUTPUT)
C
C    KMNX   (REAL) CONTAINS THE TARGET MEAN PERMEABILITY IN THE X DIRECTION.
C           (INPUT)
C
C    KSDX   (REAL) CONTAINS THE EXPECTED (TARGET) STANDARD DEVIATION, AFTER
C           LOCAL AVERAGING, OF THE PERMEABILITY IN THE X DIRECTION. (INPUT)
C
C    YKMN   ON INPUT, YKMN CONTAINS THE FIELD AVERAGE Y-DIRECTION PERMEABILITY
C           SUMMED OVER ALL THE REALIZATIONS. ON OUTPUT, YKMN WILL CONTAIN
C           THE GLOBAL AVERAGE Y-DIRECTION PERMEABILITY. (INPUT/OUTPUT)
C
C    YKSD   ON INPUT, YKSD CONTAINS THE FIELD AVERAGE Y-DIRECTION SQUARED
C           PERMEABILITY SUMMED OVER ALL THE REALIZATIONS. ON OUTPUT, YKSD
C           WILL CONTAIN THE ESTIMATED GLOBAL STANDARD DEVIATION OF THE
C           Y-DIRECTION PERMEABILITY. (INPUT/OUTPUT)
C
C    KMNY   (REAL) CONTAINS THE TARGET MEAN PERMEABILITY IN THE Y DIRECTION.
C           (INPUT)
C
C    KSDY   (REAL) CONTAINS THE EXPECTED (TARGET) STANDARD DEVIATION, AFTER
C           LOCAL AVERAGING, OF THE PERMEABILITY IN THE Y DIRECTION. (INPUT)
C
C    UPMN   ON INPUT, UPMN CONTAINS THE UPLIFT FORCE ACTING BETWEEN WALLS
C           (2-WALL CASE ONLY) SUMMED OVER ALL THE REALIZATIONS. ON OUTPUT
C           UPMN WILL CONTAIN THE AVERAGE UPLIFT FORCE. (INPUT/OUTPUT)
C
C    UPSD   ON INPUT, UPSD CONTAINS THE SUM OF THE SQUARES OF THE UPLIFT
C           FORCES CALCULATED FOR EACH REALIZATION. ON OUTPUT, UPSD WILL
C           CONTAIN THE STANDARD DEVIATION OF THE UPLIFT FORCE. (INPUT/OUTPUT)
C
C  VERBOS   LOGICAL FLAG WHICH IS TRUE IF RESULTS ARE TO BE ECHOED TO
C           STANDARD OUTPUT. (INPUT)
C
C   DEBUG   LOGICAL FLAG WHICH IS TRUE IF DEBUG DATA IS TO BE DUMPED TO
C           THE LOGICAL UNIT ISTAT. (INPUT)
C
C     NEL   INTEGER GIVING THE NUMBER OF ELEMENTS IN THE PROBLEM. (INPUT)
C---------------------------------------------------------------------------
      SUBROUTINE STATS( ISTAT, ITERM, JF, NS, NW, TI, TS, TF, KSEED,
     >                  FLMN, FLSD, EGMN, EGSD, XKMN, XKSD, KMNX, KSDX,
     >                  YKMN, YKSD, KMNY, KSDY, UPMN, UPSD, VERBOS,
     >                  DEBUG, NEL )
      LOGICAL VERBOS, DEBUG
      REAL KMNX, KSDX, KMNY, KSDY
      DATA ZERO/0.0/, ONE/1.0/, TWO/2.0/
      COMMON/DBGLAS/ ZMN, ZVR, ZCV, NN
C
   1  FORMAT(A)
   2  FORMAT()
   3  FORMAT(A,T20,E18.6,T50,E18.6)
   4  FORMAT(A,I10,A,I10,A)
   5  FORMAT(A,F16.2,A)
   6  FORMAT(A,T20,E11.4,' (',E11.4,')',T50,E11.4,' (',E11.4,')')
   7  FORMAT(A,F12.6)
 
      IF( NS .GT. 1 ) THEN
         FNS = FLOAT(NS)
         DNS = ONE/FNS
         DN1 = ONE/FLOAT(NS-1)
C						STATISTICS OF FLOW
         A = FLMN*DNS
         V = (FLSD - TWO*A*FLMN + FNS*A*A)*DN1
         FLMN = A
         FLSD = SQRT(V)
C						STATISTICS OF EXIT GRADIENT
         A = EGMN*DNS
         V = (EGSD - TWO*A*EGMN + FNS*A*A)*DN1
         EGMN = A
         EGSD = SQRT(V)
C						PERMEABILITY STATISTICS
         FNL = FLOAT(NS*NEL)
         DNL = ONE/FNL
         DNM = ONE/(FNL-ONE)
         A = XKMN*DNL
         V = (XKSD - TWO*A*XKMN + FNL*A*A)*DNM
         XKMN = A
         XKSD = SQRT(V)
         A = YKMN*DNL
         V = (YKSD - TWO*A*YKMN + FNL*A*A)*DNM
         YKMN = A
         YKSD = SQRT(V)
C						UPLIFT STATISTICS
         IF( NW .EQ. 2 ) THEN
            A = UPMN*DNS
            V = (UPSD - TWO*A*UPMN + FNS*A*A)*DN1
            UPMN = A
            UPSD = SQRT(V)
         ENDIF
C						RANDOM FIELD STATISTICS
         IF( DEBUG ) THEN
            F = FLOAT(NS*NN)
            A = ZMN/F
            V = (ZVR - TWO*A*ZMN + F*A*A)/(F-ONE)
            C = (ZCV - F*A*A)/((F-ONE)*V)
            WRITE(ISTAT,7)'Estimated Random Field Mean       = ',A
            WRITE(ISTAT,7)'Estimated Random Field Variance   = ',V
            WRITE(ISTAT,7)'Estimated Random Field Covariance = ',C
            WRITE(ISTAT,2)
         ENDIF
      ELSE
         FLSD = ZERO
         EGSD = ZERO
         XKSD = ZERO
         YKSD = ZERO
         UPSD = ZERO
      ENDIF
C						OUTPUT RESULTS
      K = ISTAT
  10  WRITE(K,1)'GLOBAL STATISTICS'
      WRITE(K,1)'================='
      WRITE(K,2)
      WRITE(K,4)'Generator seed used = ',KSEED
      WRITE(K,4)'Statistics based on ',NS,' of ',JF,' realizations.'
      WRITE(K,2)
      WRITE(K,1)'                           Average                  Sta
     >ndard Deviation'
      WRITE(K,1)'                     ===================           ====
     >================'
      WRITE(K,2)
      WRITE(K,3)'Normalized Flow',FLMN,FLSD
      WRITE(K,2)
      WRITE(K,3)'Exit Gradient',EGMN,EGSD
      IF( NW .EQ. 2 ) WRITE(K,3)'Uplift Force', UPMN, UPSD
      WRITE(K,2)
      WRITE(K,6)'Permeability (X)',XKMN,KMNX,XKSD,KSDX
      WRITE(K,6)'Permeability (Y)',YKMN,KMNY,YKSD,KSDY
      WRITE(K,2)
      WRITE(K,1)'(target statistics parenthesized)'
      WRITE(K,2)
      WRITE(K,1)'TIME LOG'
      WRITE(K,1)'========'
      WRITE(K,2)
      WRITE(K,5)'RF generator setup took ',TI,' seconds.'
      WRITE(K,5)'    RF simulations took ',TS,' seconds.'
      WRITE(K,5)'       FE analysis took ',TF,' seconds.'
 
      IF( VERBOS .AND. K .NE. ITERM ) THEN
         K = ITERM
         GO TO 10
      ENDIF
 
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                      SUBROUTINE TWOCUT                          *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY D. VAUGHAN GRIFFITHS, 1991
C
C  PURPOSE  FORMS THE COORDINATES AND STEERING VECTOR FOR 4-NODE
C           QUADRILATERALS COUNTING IN Y-DIRECTION (LAPLACES' EQUATION,
C           2-WALL PROBLEM)
C
C--------------------------------------------------------------------------
      SUBROUTINE TWOCUT(IP,IQ,NX1,NX2,NYE,IC1,IC2,AA,BB,COORD,G,NF,INF)
      REAL COORD(4,2)
      INTEGER NUM(4),G(*),NF(INF,*)
 
      NUM(1)=(IP-1)*(NYE+1)+IQ+1
      NUM(2)=NUM(1)-1
      NUM(3)=IP*(NYE+1)+IQ
      NUM(4)=NUM(3)+1
      IF(IP.LE.NX1)GOTO 10000
      IF(IP.GT.NX1+NX2)GOTO 2000
      IF(IP.GT.NX1+1)GOTO 3000
      IF(IQ.GE.IC1)GOTO 4000
      NUM(1)=NUM(1)+NYE+1
      NUM(2)=NUM(2)+NYE+1
      GOTO 5000
 4000 IF(IQ.GT.IC1)GOTO 5000
      NUM(2)=NUM(2)+NYE+1
      GOTO 5000
 3000 NUM(1)=NUM(1)+IC1
      NUM(2)=NUM(2)+IC1
 5000 NUM(3)=NUM(3)+IC1
      NUM(4)=NUM(4)+IC1
      GOTO 10000
 2000 IF(IP.GT.NX1+NX2+1)GOTO 6000
      IF(IQ.GE.IC2)GOTO 7000
      NUM(1)=NUM(1)+NYE+1+IC1
      NUM(2)=NUM(2)+NYE+1+IC1
      GOTO 9000
 7000 IF(IQ.GT.IC2)GOTO 6500
      NUM(1)=NUM(1)+IC1
      NUM(2)=NUM(2)+NYE+1+IC1
      GOTO 9000
 6500 NUM(1)=NUM(1)+IC1
      NUM(2)=NUM(2)+IC1
      GOTO 9000
 6000 NUM(1)=NUM(1)+IC1+IC2
      NUM(2)=NUM(2)+IC1+IC2
 9000 NUM(3)=NUM(3)+IC1+IC2
      NUM(4)=NUM(4)+IC1+IC2
10000 DO 1 I=1,4
    1 G(I)=NF(NUM(I),1)
      COORD(1,1)=(IP-1)*AA
      COORD(2,1)=(IP-1)*AA
      COORD(3,1)=IP*AA
      COORD(4,1)=IP*AA
      COORD(1,2)=-IQ*BB
      COORD(2,2)=-(IQ-1)*BB
      COORD(3,2)=-(IQ-1)*BB
      COORD(4,2)=-IQ*BB
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                      SUBROUTINE TYPE                            *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY D. VAUGHAN GRIFFITHS, 1991
C
C  PURPOSE  SETS UP PARAMETERS AND BOUNDARY CONDITIONS FOR 0-WALL, 1-WALL,
C           AND 2-WALL FEM SOIL FLOW PROBLEMS.
C
C--------------------------------------------------------------------------
      SUBROUTINE TYPE(NW,NX1,NX2,NX3,NYE,IC1,IC2,
     +                IFIX,NXE,NN,N,IW,IR,IG,NO,NF,INF,ITERM)
      INTEGER NO(*),NF(INF,*)
C
      IF(NW.EQ.0) THEN
         NXE = NX1
         NN  = (NXE+1)*(NYE+1)
         N   = NN - NYE - 1
         IW  = NYE + 2
         IG  = (NYE+1)*NXE
         NC  = IG
         DO 19 I = 1, NYE+1
            NF(NC+I,1) = 1
  19     CONTINUE
         IFIX = NYE + 1
         DO 20 I = 1, IFIX
            NO(I) = I
  20     CONTINUE
      ELSE IF(NW.EQ.1)THEN
         NXE = NX1 + NX2
         NN  = (NXE+1)*(NYE+1) + IC1
         N   = NN - NX2 - 1
         IW  = NYE + 2 + IC1
         IG  = (NYE+1)*(NX1+1)
         NC           = IG + 1
         NF(NC,1)     = 1
         NF(NC+IC1,1) = 1
         DO 21 I = 1, NX2-1
            NF(NC+IC1+I*(NYE+1),1) = 1
  21     CONTINUE
         IFIX = NX1 + 1
         DO 22 I = 1, IFIX
            NO(I) = 1+(I-1)*(NYE+1)
  22     CONTINUE
      ELSE IF( NW .EQ. 2 ) THEN
         NXE = NX1 + NX2 + NX3
         NN  = (NXE+1)*(NYE+1) + IC1 + IC2
         N   = NN - NX3 - 1
         IW  = NYE + 2 + MAX0(IC1,IC2)
         IG  = (NYE+1)*(NX1+NX2+1)+IC1
         NC           = IG + 1
         NF(NC,1)     = 1
         NF(NC+IC2,1) = 1
         DO 23 I = 1, NX3-1
            NF(NC+IC2+I*(NYE+1),1) = 1
  23     CONTINUE
         IFIX = NX1 + 1
         DO 24 I = 1, IFIX
            NO(I) = 1+(I-1)*(NYE+1)
  24     CONTINUE
      ELSE
         CALL ERROR( 6, NW, ITERM )
      END IF
 
      J = 0
      DO 25 I = 1, NN
         IF( NF(I,1) .EQ. 1 ) THEN
            NF(I,1) = 0
         ELSE
            J       = J + 1
            NF(I,1) = J
         END IF
  25  CONTINUE
      IR  = N*(IW+1)
 
      RETURN
      END
C  **********************************************************************
C  *                                                                    *
C  *                          FUNCTION UPLIFT                           *
C  *                                                                    *
C  **********************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY D. VAUGHAN GRIFFITHS, 1991
C
C  PURPOSE  RETURNS THE UPLIFT PRESSURE ACTING BETWEEN WALLS
C
C  ARGUMENTS ARE AS FOLLOWS;
C
C    LOADS    VECTOR OF NODAL PRESSURES. (INPUT)
C
C    NX1      NUMBER OF ELEMENTS IN THE X-DIRECTION TO THE LEFT OF THE
C             LEFT-MOST WALL. (INPUT)
C
C    NX2      NUMBER OF ELEMENTS IN THE X-DIRECTION BETWEEN THE TWO WALLS.
C             (INPUT)
C
C    NYE      TOTAL NUMBER OF ELEMENTS IN THE Y-DIRECTION. (INPUT)
C
C    IC1      NUMBER OF ELEMENTS IN THE Y-DIRECTION AGAINST THE LEFT-MOST
C             WALL. (INPUT)
C--------------------------------------------------------------------------
      REAL FUNCTION UPLIFT(LOADS,NX1,NX2,NYE,IC1)
      REAL LOADS(*)
      DATA HALF/0.5/
 
      NX = NX1 + 1
      NY = NYE + 1
      IF( IC1 .GE. 1 ) THEN
         S = LOADS(NX*NY+1)
      ELSE
         S = LOADS(NX1*NY+1)
      ENDIF
      NC = NX*NY + IC1 + 1
 
      UPLIFT = HALF*(S + LOADS(NC+(NX2-1)*NY))
      DO 10 I = 2, NX2
         UPLIFT = UPLIFT + LOADS(NC + (I-2)*NY)
  10  CONTINUE
 
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                     FUNCTION DAVEX2                             *
C  *                                                                 *
C  *******************************************************************
C  DOUBLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, TUNS, 1991
C
C  PURPOSE  RETURNS THE VARIANCE OF A LOCAL AVERAGE OF A 2-D PROCESS
C           HAVING AN EXPONENTIAL COVARIANCE FUNCTION.
C
C  RETURNS THE VARIANCE OF A LOCAL AREA AVERAGE OF THE PROCESS HAVING
C  THE EXPONENTIAL COVARIANCE FUNCTION
C
C         B(X,Y) = EXP{ -SQRT[(2*X)^2 + (2*Y)^2] }
C
C  THE VARIANCE FUNCTION CAN BE APPROXIMATED BY
C
C         V(X,Y) = 0.5*[V(X)*V(Y|X) + V(Y)*V(X|Y)]
C
C  WHERE INDIVIDUAL FUNCTIONS ARE BASED ON THE 1-D ANALYTICAL MODEL
C
C                        1
C         V(T) = ---------------------
C                [1 + (T/A)^1.5]^(2/3)
C
C  AND WHERE `A' IS THE DIRECTIONAL EFFECTIVE SCALE OF FLUCTUATION. SEE
C  PAGES 247-249 OF "RANDOM FIELDS" BY E. VANMARCKE AND MY THESIS FOR
C  MORE DETAILS. THE ARGUMENTS TO THIS ROUTINE ARE JUST THE DIMENSIONS
C  OF THE AVERAGING AREA.
C  THE PARAMETERS OF THE PROCESS ARE BROUGHT IN THROUGH THE COMMON DPARAM
C  IN WHICH PB IS THE RATIO OF THE CHARACTERISTIC
C  AREA TO THE PRODUCT OF THE DIRECTIONAL SCALES OF
C  FLUCTUATION (SEE VANMARCKE, PG 240) (PB = PI/2 FOR EXPONENTIAL PROCESSES).
C  `EMAX' IS THE MAXIMUM EXPONENT EXP() CAN HANDLE BEFORE OVER OR UNDERFLOW
C  OCCURS (IN DOUBLE PRECISION). (SINGLE PRECISION VALUE IS AROUND 85.)
C---------------------------------------------------------------------------
      REAL FUNCTION VARFNC( X, Y )
      COMMON/DPARAM/ PA, PB, PX, PY, PZ
      DATA HALF/0.50/, ONE/1.0/, ONEPT5/1.50/
      DATA TWOTHD/-0.66666666666666666666666666667/
      DATA EMAX/690.0/
C						AVOID UNDERFLOW REPORTS
      CALL XUFLOW(0)
      XX  = ABS(X)
      YY  = ABS(Y)
      PP  = PB*PB
      R1  = XX*XX/PP
      R2  = YY*YY/PP
C						AVOID UNDERFLOW
      IF( R1 .LT. EMAX ) THEN
         A1  = (PB + (ONE - PB)*EXP(-R1))
      ELSE
         A1  = PB
      ENDIF
      IF( R2 .LT. EMAX ) THEN
         A2  = (PB + (ONE - PB)*EXP(-R2))
      ELSE
         A2 = PB
      ENDIF
 
      VARFNC = HALF*(((ONE+YY**ONEPT5)*(ONE+(XX/A2)**ONEPT5))**TWOTHD
     >             + ((ONE+XX**ONEPT5)*(ONE+(YY/A1)**ONEPT5))**TWOTHD)
 
      RETURN
      END
C  *******************************************************************
C  *                                                                 *
C  *                    SUBROUTINE VNORM                             *
C  *                                                                 *
C  *******************************************************************
C  SINGLE PRECISION VERSION 1.0
C  WRITTEN BY GORDON A. FENTON, PRINCETON, 1989.
C
C  PURPOSE  RETURN A VECTOR OF N(0,1) RANDOM DISTRIBUTED REALIZATIONS
C
C  RETURNS A NORMALLY DISTRIBUTED, ZERO MEAN, UNIT VARIANCE RANDOM VECTOR
C  IN THE ARGUMENT `G'. `N' IS THE LENGTH OF THE DESIRED VECTOR.
C  ENSURE THAT THE RANDOM NUMBER GENERATOR IS INITIALIZED PRIOR TO CALLING
C  THIS ROUTINE.
C---------------------------------------------------------------------------
      SUBROUTINE VNORM( G, N )
      REAL G(1)
      DATA ONE/1.0/, TWO/2.0/, TWOPI/6.2831853071795864769/
C
      NN = N/2
      NN = 2*NN
      DO 10 I = 1, NN, 2
         A      = TWOPI*G05CAF(X)
         R      = G05CAF(X)
         R      = SQRT(-TWO*ALOG(ONE-R))
         G(I)   = R*COS(A)
         G(I+1) = R*SIN(A)
  10  CONTINUE
 
      IF( N .EQ. NN ) RETURN
C                                     FOR N ODD, SET THE LAST VALUE
      A = TWOPI*G05CAF(X)
      R = G05CAF(X)
      R = SQRT(-TWO*ALOG(ONE-R))
      G(N) = R*COS(A)
 
      RETURN
      END
C  *********************************************************************
C  *                                                                   *
C  *                    INTEGER FUNCTION LNBLNK                        *
C  *                                                                   *
C  *********************************************************************
C   VERSION 1.0
C   WRITTEN BY G. A. FENTON, PRINCETON, 1989.
C
C   PURPOSE  TO RETURN THE NON-BLANK LENGTH OF A STRING.
C
C   RETURNS THE INDEX OF THE LAST NON-BLANK CHARACTER IN THE ARGUMENT STRING.
C   RETURNS 0 IF THEY ARE ALL BLANK.
C   THIS ROUTINE IS FOR USE IF YOUR VERSION OF FORTRAN DOES NOT HAVE LNBLNK.
C   IT DOES, HOWEVER, DEPEND ON THE FUNCTION `LEN' WHICH RETURNS THE
C   DIMENSIONED LENGTH OF THE STRING `STR'.
C----------------------------------------------------------------------------
      INTEGER FUNCTION LNBLNK( STR )
      CHARACTER*(*) STR
 
      I = LEN( STR )
      DO 10 LNBLNK = I, 1, -1
         IF( STR(LNBLNK:LNBLNK) .NE. ' ' ) RETURN
  10  CONTINUE
      LNBLNK = 0
 
      RETURN
      END
