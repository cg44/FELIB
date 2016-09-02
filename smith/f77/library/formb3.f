      SUBROUTINE FORMB3(BEE,IBEE,DERIV,IDERIV,NOD)
C
C      THIS SUBROUTINE FORMS THE 3-D STRAIN-DISPLACEMENT MATRIX
C
      REAL BEE(IBEE,*),DERIV(IDERIV,*)
      DO 1 M=1,NOD
      N=3*M
      K=N-1
      L=K-1
      X=DERIV(1,M)
      BEE(1,L)=X
      BEE(4,K)=X
      BEE(6,N)=X
      Y=DERIV(2,M)
      BEE(2,K)=Y
      BEE(4,L)=Y
      BEE(5,N)=Y
      Z=DERIV(3,M)
      BEE(3,N)=Z
      BEE(5,K)=Z
      BEE(6,L)=Z
    1 CONTINUE
      RETURN
      END
