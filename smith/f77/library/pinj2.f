      SUBROUTINE PINJ2(KM,EA,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX FOR AN
C      INCLINED 2-D PIN-JOINTED ELEMENT
C
      REAL KM(4,4),COORD(ICOORD,*)
      X1=COORD(IP,1)
      Y1=COORD(IP,2)
      X2=COORD(IP,3)
      Y2=COORD(IP,4)
      ELL=SQRT((Y2-Y1)**2+(X2-X1)**2)
      CS=(X2-X1)/ELL
      SN=(Y2-Y1)/ELL
      A=CS*CS
      B=SN*SN
      C=CS*SN
      KM(1,1)=A
      KM(3,3)=A
      KM(1,3)=-A
      KM(3,1)=-A
      KM(2,2)=B
      KM(4,4)=B
      KM(2,4)=-B
      KM(4,2)=-B
      KM(1,2)=C
      KM(2,1)=C
      KM(3,4)=C
      KM(4,3)=C
      KM(1,4)=-C
      KM(4,1)=-C
      KM(2,3)=-C
      KM(3,2)=-C
      DO 1 I=1,4
      DO 1 J=1,4
    1 KM(I,J)=KM(I,J)*EA/ELL
      RETURN
      END
