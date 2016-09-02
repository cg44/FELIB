      SUBROUTINE FORMM3(STRESS,M1,M2,M3)
C
C      THIS SUBROUTINE FORMS THE DERIVATIVES OF THE INVARIANTS
C      WITH RESPECT TO THE STRESSES (3-D)
C
      REAL STRESS(*),M1(6,6),M2(6,6),M3(6,6)
      SX=STRESS(1)
      SY=STRESS(2)
      SZ=STRESS(3)
      TXY=STRESS(4)
      TYZ=STRESS(5)
      TZX=STRESS(6)
      SIGM=(SX+SY+SZ)/3.
      DX=SX-SIGM
      DY=SY-SIGM
      DZ=SZ-SIGM
      DO 1 I=1,6
      DO 1 J=I,6
      M1(I,J)=0.
    1 M2(I,J)=0.
      DO 2 I=1,3
      DO 2 J=1,3
    2 M1(I,J)=1./(3.*SIGM)
      DO 3 I=1,3
      M2(I,I)=2.
    3 M2(I+3,I+3)=6.
      M2(1,2)=-1.
      M2(1,3)=-1.
      M2(2,3)=-1.
      M3(1,1)=DX
      M3(1,2)=DZ
      M3(1,3)=DY
      M3(1,4)=TXY
      M3(1,5)=-2.*TYZ
      M3(1,6)=TZX
      M3(2,2)=DY
      M3(2,3)=DX
      M3(2,4)=TXY
      M3(2,5)=TYZ
      M3(2,6)=-2.*TZX
      M3(3,3)=DZ
      M3(3,4)=-2.*TXY
      M3(3,5)=TYZ
      M3(3,6)=TZX
      M3(4,4)=-3.*DZ
      M3(4,5)=3.*TZX
      M3(4,6)=3.*TYZ
      M3(5,5)=-3.*DX
      M3(5,6)=3.*TXY
      M3(6,6)=-3.*DY
      DO 4 I=1,6
      DO 4 J=I,6
      M1(I,J)=M1(I,J)/3.
      M1(J,I)=M1(I,J)
      M2(I,J)=M2(I,J)/3.
      M2(J,I)=M2(I,J)
      M3(I,J)=M3(I,J)/3.
    4 M3(J,I)=M3(I,J)
      RETURN
      END
