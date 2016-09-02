      SUBROUTINE FORMD3(DEE,IDEE,E,V)
C
C      THIS SUBROUTINE FORMS THE 3-D STRAIN
C      STRESS/STRAIN MATRIX
C
      REAL DEE(IDEE,*)
      V1=V/(1.-V)
      VV=(1.-2.*V)/(1.-V)*.5
      DO 1 I=1,6
      DO 1 J=1,6
    1 DEE(I,J)=0.
      DEE(1,1)=1.
      DEE(2,2)=1.
      DEE(3,3)=1.
      DEE(1,2)=V1
      DEE(2,1)=V1
      DEE(1,3)=V1
      DEE(3,1)=V1
      DEE(2,3)=V1
      DEE(3,2)=V1
      DEE(4,4)=VV
      DEE(5,5)=VV
      DEE(6,6)=VV
      DO 2 I=1,6
      DO 2 J=1,6
    2 DEE(I,J)=DEE(I,J)*E/(2.*(1.+V)*VV)
      RETURN
      END
