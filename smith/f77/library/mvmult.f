      SUBROUTINE MVMULT(M,IM,V,K,L,Y)
C
C      THIS SUBROUTINE MULTIPLIES A MATRIX BY A VECTOR
C
      REAL M(IM,*),V(*),Y(*)
      DO 1 I=1,K
      X=0.
      DO 2 J=1,L
    2 X=X+M(I,J)*V(J)
      Y(I)=X
    1 CONTINUE
      RETURN
      END
