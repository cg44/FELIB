      SUBROUTINE FORMM(STRESS,M1,M2,M3)
C
C      THIS SUBROUTINE FORMS THE DERIVATIVES OF THE INVARIANTS
C      WITH RESPECT TO THE STRESSES
C
      REAL STRESS(*),M1(4,4),M2(4,4),M3(4,4)
      SX=STRESS(1)
      SY=STRESS(2)
      TXY=STRESS(3)
      SZ=STRESS(4)
      DX=(2.*SX-SY-SZ)/3.
      DY=(2.*SY-SZ-SX)/3.
      DZ=(2.*SZ-SX-SY)/3.
      SIGM=(SX+SY+SZ)/3.
      DO 1 I=1,4
      DO 1 J=1,4
      M1(I,J)=0.
      M2(I,J)=0.
    1 M3(I,J)=0.
      M1(1,1)=1.
      M1(1,2)=1.
      M1(2,1)=1.
      M1(1,4)=1.
      M1(4,1)=1.
      M1(2,2)=1.
      M1(2,4)=1.
      M1(4,2)=1.
      M1(4,4)=1.
      DO 2 I=1,4
      DO 2 J=1,4
    2 M1(I,J)=M1(I,J)/(9.*SIGM)
      M2(1,1)=.6666666666666666
      M2(2,2)=.6666666666666666
      M2(4,4)=.6666666666666666
      M2(3,3)=2.
      M2(2,4)=-.3333333333333333
      M2(4,2)=-.3333333333333333
      M2(1,2)=-.3333333333333333
      M2(2,1)=-.3333333333333333
      M2(1,4)=-.3333333333333333
      M2(4,1)=-.3333333333333333
      M3(1,1)=DX/3.
      M3(2,4)=DX/3.
      M3(4,2)=DX/3.
      M3(2,2)=DY/3.
      M3(1,4)=DY/3.
      M3(4,1)=DY/3.
      M3(4,4)=DZ/3.
      M3(1,2)=DZ/3.
      M3(2,1)=DZ/3.
      M3(3,3)=-DZ
      M3(3,4)=-2.*TXY/3.
      M3(4,3)=-2.*TXY/3.
      M3(1,3)=TXY/3.
      M3(3,1)=TXY/3.
      M3(2,3)=TXY/3.
      M3(3,2)=TXY/3.
      RETURN
      END
