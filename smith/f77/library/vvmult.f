      SUBROUTINE VVMULT(V1,V2,PROD,IPROD,M,N)
C
C      THIS SUBROUTINE FORMS A VECTOR PRODUCT
C
      REAL V1(*),V2(*),PROD(IPROD,*)
      DO 1 I=1,M
      DO 1 J=1,N
    1 PROD(I,J)=V1(I)*V2(J)
      RETURN
      END
