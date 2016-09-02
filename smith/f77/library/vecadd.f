      SUBROUTINE VECADD(A,B,C,N)
C
C      THIS SUBROUTINE ADDS VECTORS  A+B=C
C
      REAL A(*),B(*),C(*)
      DO 1 I=1,N
    1 C(I)=A(I)+B(I)
      RETURN
      END
