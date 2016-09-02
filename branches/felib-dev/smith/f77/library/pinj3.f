      SUBROUTINE PINJ3(KM,EA,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX FOR A
C      GENERAL 3-D PIN-JOINTED ELEMENT
C
      REAL KM(6,6),COORD(ICOORD,*)
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
      XL=XL/ELL
      YL=YL/ELL
      ZL=ZL/ELL
      A=XL*XL
      B=YL*YL
      C=ZL*ZL
      D=XL*YL
      E=YL*ZL
      F=ZL*XL
      KM(1,1)=A
      KM(4,4)=A
      KM(2,2)=B
      KM(5,5)=B
      KM(3,3)=C
      KM(6,6)=C
      KM(1,2)=D
      KM(2,1)=D
      KM(4,5)=D
      KM(5,4)=D
      KM(2,3)=E
      KM(3,2)=E
      KM(5,6)=E
      KM(6,5)=E
      KM(1,3)=F
      KM(3,1)=F
      KM(4,6)=F
      KM(6,4)=F
      KM(1,4)=-A
      KM(4,1)=-A
      KM(2,5)=-B
      KM(5,2)=-B
      KM(3,6)=-C
      KM(6,3)=-C
      KM(1,5)=-D
      KM(5,1)=-D
      KM(2,4)=-D
      KM(4,2)=-D
      KM(2,6)=-E
      KM(6,2)=-E
      KM(3,5)=-E
      KM(5,3)=-E
      KM(1,6)=-F
      KM(6,1)=-F
      KM(3,4)=-F
      KM(4,3)=-F
      DO 1 I=1,6
      DO 1 J=1,6
    1 KM(I,J)=KM(I,J)*EA/ELL
      RETURN
      END
