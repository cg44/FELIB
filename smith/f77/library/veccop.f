      SUBROUTINE VECCOP(A,B,N)
C
C      THIS SUBROUTINE COPIES VECTOR A INTO VECTOR B
C
      REAL A(*),B(*)
      DO 1 I=1,N
    1 B(I)=A(I)
      RETURN
      END
