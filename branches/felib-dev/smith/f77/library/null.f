      SUBROUTINE NULL(A,IA,M,N)
C
C      THIS SUBROUTINE NULLS A 2-D ARRAY
C
      REAL A(IA,*)
      DO 1 I=1,M
      DO 1 J=1,N
    1 A(I,J)=0.0
      RETURN
      END
