      SUBROUTINE FMDRAD(DEE,IDEE,E,V)
C
C      THIS SUBROUTINE FORMS THE ELASTIC AXISYMMETRIC
C      STRESS/STRAIN MATRIX
C
      REAL DEE(IDEE,*)
      V1=1.-V
      C=E/((1.+V)*(1.-2.*V))
      DEE(1,1)=V1*C
      DEE(2,2)=V1*C
      DEE(3,3)=.5*C*(1.-2.*V)
      DEE(4,4)=V1*C
      DEE(1,2)=V*C
      DEE(2,1)=V*C
      DEE(1,3)=0.
      DEE(3,1)=0.
      DEE(1,4)=V*C
      DEE(4,1)=V*C
      DEE(2,3)=0.
      DEE(3,2)=0.
      DEE(2,4)=V*C
      DEE(4,2)=V*C
      DEE(4,3)=0.
      DEE(3,4)=0.
      RETURN
      END
