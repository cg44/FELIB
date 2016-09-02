      SUBROUTINE BMCOL3(KM,EA,EIY,EIZ,GJ,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX OF A
C      GENERAL 3-D BEAM-COLUMN ELEMENT
C
      REAL KM(12,12),COORD(ICOORD,*),T(12,12),TT(12,12),R0(3,3),C(12,12)
      PI=4.*ATAN(1.)
      GAMA=COORD(IP,7)*PI/180.
      X1=COORD(IP,1)
      Y1=COORD(IP,2)
      Z1=COORD(IP,3)
      X2=COORD(IP,4)
      Y2=COORD(IP,5)
      Z2=COORD(IP,6)
      XL=X2-X1
      YL=Y2-Y1
      ZL=Z2-Z1
      ELL=SQRT(XL*XL+YL*YL+ZL*ZL)
      CG=COS(GAMA)
      SG=SIN(GAMA)
      DEN=ELL*SQRT(XL*XL+ZL*ZL)
      DO 1 I=1,12
      DO 1 J=1,12
      KM(I,J)=0.
      T(I,J)=0.
    1 TT(I,J)=0.
      A1=EA/ELL
      A2=12.*EIZ/(ELL*ELL*ELL)
      A3=12.*EIY/(ELL*ELL*ELL)
      A4=6.*EIZ/(ELL*ELL)
      A5=6.*EIY/(ELL*ELL)
      A6=4.*EIZ/ELL
      A7=4.*EIY/ELL
      A8=GJ/ELL
      KM(1,1)=A1
      KM(7,7)=A1
      KM(1,7)=-A1
      KM(7,1)=-A1
      KM(2,2)=A2
      KM(8,8)=A2
      KM(2,8)=-A2
      KM(8,2)=-A2
      KM(3,3)=A3
      KM(9,9)=A3
      KM(3,9)=-A3
      KM(9,3)=-A3
      KM(4,4)=A8
      KM(10,10)=A8
      KM(4,10)=-A8
      KM(10,4)=-A8
      KM(5,5)=A7
      KM(11,11)=A7
      KM(5,11)=.5*A7
      KM(11,5)=.5*A7
      KM(6,6)=A6
      KM(12,12)=A6
      KM(6,12)=.5*A6
      KM(12,6)=.5*A6
      KM(2,6)=A4
      KM(6,2)=A4
      KM(2,12)=A4
      KM(12,2)=A4
      KM(6,8)=-A4
      KM(8,6)=-A4
      KM(8,12)=-A4
      KM(12,8)=-A4
      KM(5,9)=A5
      KM(9,5)=A5
      KM(9,11)=A5
      KM(11,9)=A5
      KM(3,5)=-A5
      KM(5,3)=-A5
      KM(3,11)=-A5
      KM(11,3)=-A5
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
      T(I+K,J+K)=X
    2 TT(J+K,I+K)=X
      DO 3 I=1,12
      DO 3 J=1,12
      SUM=0.
      DO 4 K=1,12
    4 SUM=SUM+KM(I,K)*T(K,J)
    3 C(I,J)=SUM
      DO 5 I=1,12
      DO 5 J=1,12
      SUM=0.
      DO 6 K=1,12
    6 SUM=SUM+TT(I,K)*C(K,J)
    5 KM(I,J)=SUM
      RETURN
      END
