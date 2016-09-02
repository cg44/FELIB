      SUBROUTINE FMLAG9(DER,IDER,FUN,SAMP,ISAMP,I,J)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND
C      THEIR DERIVATIVES FOR 9-NODED QUADRILATERAL ELEMENTS
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      ETA=SAMP(I,1)
      XI=SAMP(J,1)
      ETAM=ETA-1.
      XIM=XI-1.
      ETAP=ETA+1.
      XIP=XI+1.
      X2P1=2.*XI+1.
      X2M1=2.*XI-1.
      E2P1=2.*ETA+1.
      E2M1=2.*ETA-1.
      FUN(1)=.25*XI*XIM*ETA*ETAM
      FUN(2)=-.5*XI*XIM*ETAP*ETAM
      FUN(3)=.25*XI*XIM*ETA*ETAP
      FUN(4)=-.5*XIP*XIM*ETA*ETAP
      FUN(5)=.25*XI*XIP*ETA*ETAP
      FUN(6)=-.5*XI*XIP*ETAP*ETAM
      FUN(7)=.25*XI*XIP*ETA*ETAM
      FUN(8)=-.5*XIP*XIM*ETA*ETAM
      FUN(9)=XIP*XIM*ETAP*ETAM
      DER(1,1)=.25*X2M1*ETA*ETAM
      DER(1,2)=-.5*X2M1*ETAP*ETAM
      DER(1,3)=.25*X2M1*ETA*ETAP
      DER(1,4)=-XI*ETA*ETAP
      DER(1,5)=.25*X2P1*ETA*ETAP
      DER(1,6)=-.5*X2P1*ETAP*ETAM
      DER(1,7)=.25*X2P1*ETA*ETAM
      DER(1,8)=-XI*ETA*ETAM
      DER(1,9)=2.*XI*ETAP*ETAM
      DER(2,1)=.25*XI*XIM*E2M1
      DER(2,2)=-XI*XIM*ETA
      DER(2,3)=.25*XI*XIM*E2P1
      DER(2,4)=-.5*XIP*XIM*E2P1
      DER(2,5)=.25*XI*XIP*E2P1
      DER(2,6)=-XI*XIP*ETA
      DER(2,7)=.25*XI*XIP*E2M1
      DER(2,8)=-.5*XIP*XIM*E2M1
      DER(2,9)=2.*XIP*XIM*ETA
      RETURN
      END