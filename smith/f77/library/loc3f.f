      SUBROUTINE LOC3F(LOCAL,GLOBAL,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE TRANSFORMS THE END REACTION AND MOMENTS
C      INTO THE ELEMENT'S LOCAL COORDINATE SYSTEM (3-D)
C
      REAL LOCAL(*),GLOBAL(*),COORD(ICOORD,*),R0(3,3),T(12,12)
      DO 1 I=1,12
      DO 1 J=1,12
    1 T(I,J)=0.
      X1=COORD(IP,1)
      Y1=COORD(IP,2)
      Z1=COORD(IP,3)
      X2=COORD(IP,4)
      Y2=COORD(IP,5)
      Z2=COORD(IP,6)
      PI=4.*ATAN(1.)
      GAMA=COORD(IP,7)*PI/180.
      CG=COS(GAMA)
      SG=SIN(GAMA)
      XL=X2-X1
      YL=Y2-Y1
      ZL=Z2-Z1
      ELL=SQRT(XL*XL+YL*YL+ZL*ZL)
      DEN=ELL*SQRT(XL*XL+ZL*ZL)
      IF(DEN.EQ.0.)GOTO 50
      R0(1,1)=XL/ELL
      R0(1,2)=YL/ELL
      R0(1,3)=ZL/ELL
      R0(2,1)=(-XL*YL*CG-ELL*ZL*SG)/DEN
      R0(2,2)=DEN*CG/(ELL*ELL)
      R0(2,3)=(-YL*ZL*CG+ELL*XL*SG)/DEN
      R0(3,1)=(XL*YL*SG-ELL*ZL*CG)/DEN
      R0(3,2)=-DEN*SG/(ELL*ELL)
      R0(3,3)=(YL*ZL*SG+ELL*XL*CG)/DEN
      GOTO 60
   50 R0(1,1)=0.
      R0(1,3)=0.
      R0(2,2)=0.
      R0(3,2)=0.
      R0(1,2)=1.
      R0(2,1)=-CG
      R0(3,3)=CG
      R0(2,3)=SG
      R0(3,1)=SG
   60 CONTINUE
      DO 2 I=1,3
      DO 2 J=1,3
      X=R0(I,J)
      DO 2 K=0,9,3
    2 T(I+K,J+K)=X
      DO 3 I=1,12
      SUM=0.
      DO 4 J=1,12
    4 SUM=SUM+T(I,J)*GLOBAL(J)
    3 LOCAL(I)=SUM
      RETURN
      END
