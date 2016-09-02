      SUBROUTINE NULL3(A,IA1,IA2,L,M,N)
C
C      THIS SUBROUTINE NULLS A 3-D MATRIX
C
      REAL A(IA1,IA2,*)
      DO 1 I=1,L
      DO 1 J=1,M
      DO 1 K=1,N
    1 A(I,J,K)=0.
      RETURN
      END
