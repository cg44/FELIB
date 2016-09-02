      SUBROUTINE FMTRI6(DER,IDER,FUN,SAMP,ISAMP,I)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 6-NODED TRIANGULAR ELEMENTS
C
      REAL DER(IDER,*),SAMP(ISAMP,*),FUN(*)
      C1=SAMP(I,1)
      C2=SAMP(I,2)
      C3=1.-C1-C2
      FUN(1)=(2.*C1-1.)*C1
      FUN(2)=4.*C1*C2
      FUN(3)=(2.*C2-1.)*C2
      FUN(4)=4.*C2*C3
      FUN(5)=(2.*C3-1.)*C3
      FUN(6)=4.*C3*C1
      DER(1,1)=4.*C1-1.
      DER(1,2)=4.*C2
      DER(1,3)=0.
      DER(1,4)=-4.*C2
      DER(1,5)=-(4.*C3-1.)
      DER(1,6)=4.*(C3-C1)
      DER(2,1)=0.
      DER(2,2)=4.*C1
      DER(2,3)=4.*C2-1.
      DER(2,4)=4.*(C3-C2)
      DER(2,5)=-(4.*C3-1.)
      DER(2,6)=-4.*C1
      RETURN
      END
