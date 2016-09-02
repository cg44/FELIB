      SUBROUTINE FMLIN3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES FOR 8-NODED BRICK ELEMENTS
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      ETA=SAMP(I,1)
      XI=SAMP(J,1)
      ZETA=SAMP(K,1)
      ETAM=1.-ETA
      XIM=1.-XI
      ZETAM=1.-ZETA
      ETAP=ETA+1.
      XIP=XI+1.
      ZETAP=ZETA+1.
      FUN(1)=.125*XIM*ETAM*ZETAM
      FUN(2)=.125*XIM*ETAM*ZETAP
      FUN(3)=.125*XIP*ETAM*ZETAP
      FUN(4)=.125*XIP*ETAM*ZETAM
      FUN(5)=.125*XIM*ETAP*ZETAM
      FUN(6)=.125*XIM*ETAP*ZETAP
      FUN(7)=.125*XIP*ETAP*ZETAP
      FUN(8)=.125*XIP*ETAP*ZETAM
      DER(1,1)=-.125*ETAM*ZETAM
      DER(1,2)=-.125*ETAM*ZETAP
      DER(1,3)=.125*ETAM*ZETAP
      DER(1,4)=.125*ETAM*ZETAM
      DER(1,5)=-.125*ETAP*ZETAM
      DER(1,6)=-.125*ETAP*ZETAP
      DER(1,7)=.125*ETAP*ZETAP
      DER(1,8)=.125*ETAP*ZETAM
      DER(2,1)=-.125*XIM*ZETAM
      DER(2,2)=-.125*XIM*ZETAP
      DER(2,3)=-.125*XIP*ZETAP
      DER(2,4)=-.125*XIP*ZETAM
      DER(2,5)=.125*XIM*ZETAM
      DER(2,6)=.125*XIM*ZETAP
      DER(2,7)=.125*XIP*ZETAP
      DER(2,8)=.125*XIP*ZETAM
      DER(3,1)=-.125*XIM*ETAM
      DER(3,2)=.125*XIM*ETAM
      DER(3,3)=.125*XIP*ETAM
      DER(3,4)=-.125*XIP*ETAM
      DER(3,5)=-.125*XIM*ETAP
      DER(3,6)=.125*XIM*ETAP
      DER(3,7)=.125*XIP*ETAP
      DER(3,8)=-.125*XIP*ETAP
      RETURN
      END
