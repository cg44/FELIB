      SUBROUTINE STAB2D(KM,EA,EI,IP,COORD,ICOORD,PAX)
C
C      THIS SUBROUTINE FORMS THE STIFFNESS MATRIX OF AN
C      INCLINED 2-D BEAM-COLUMN ELEMENT TAKING ACCOUNT
C      OF THE EFFECTS OF AXIAL FORCES
C
      REAL KM(6,6),COORD(ICOORD,*)
      X1=COORD(IP,1)
      Y1=COORD(IP,2)
      X2=COORD(IP,3)
      Y2=COORD(IP,4)
      ELL=SQRT((X2-X1)**2+(Y2-Y1)**2)
      C=(X2-X1)/ELL
      S=(Y2-Y1)/ELL
      ALP=.5*ELL*SQRT(ABS(PAX)/EI)
      IF(PAX.GT.(5.E-5*EI/ELL**2))THEN
      SBAR=ALP*(1.-2.*ALP/TANH(2.*ALP))/(TANH(ALP)-ALP)
      CBAR=(2.*ALP-SINH(2.*ALP))/(SINH(2.*ALP)-2.*ALP*COSH(2.*ALP))
      ELSE IF(PAX.LT.(-5.E-5*EI/ELL**2))THEN
      SBAR=ALP*(1.-2.*ALP/TAN(2.*ALP))/(TAN(ALP)-ALP)
      CBAR=(2.*ALP-SIN(2.*ALP))/(SIN(2.*ALP)-2.*ALP*COS(2.*ALP))
      ELSE
      SBAR=4.
      CBAR=.5
      END IF
      BET1=2.*SBAR*(1.+CBAR)+PAX*ELL**2/EI
      BET2=SBAR*(1.+CBAR)
      E1=EA/ELL
      E2=BET1*EI/(ELL*ELL*ELL)
      E3=EI/ELL
      E4=BET2*EI/(ELL*ELL)
      KM(1,1)=C*C*E1+S*S*E2
      KM(4,4)=KM(1,1)
      KM(1,2)=S*C*(E1-E2)
      KM(2,1)=KM(1,2)
      KM(4,5)=KM(1,2)
      KM(5,4)=KM(4,5)
      KM(1,3)=-S*E4
      KM(3,1)=KM(1,3)
      KM(1,6)=KM(1,3)
      KM(6,1)=KM(1,6)
      KM(3,4)=S*E4
      KM(4,3)=KM(3,4)
      KM(4,6)=KM(3,4)
      KM(6,4)=KM(4,6)
      KM(1,4)=-KM(1,1)
      KM(4,1)=KM(1,4)
      KM(1,5)=S*C*(-E1+E2)
      KM(5,1)=KM(1,5)
      KM(2,4)=KM(1,5)
      KM(4,2)=KM(2,4)
      KM(2,2)=S*S*E1+C*C*E2
      KM(5,5)=KM(2,2)
      KM(2,5)=-KM(2,2)
      KM(5,2)=KM(2,5)
      KM(2,3)=C*E4
      KM(3,2)=KM(2,3)
      KM(2,6)=KM(2,3)
      KM(6,2)=KM(2,6)
      KM(3,3)=SBAR*E3
      KM(6,6)=KM(3,3)
      KM(3,5)=-C*E4
      KM(5,3)=KM(3,5)
      KM(5,6)=KM(3,5)
      KM(6,5)=KM(5,6)
      KM(3,6)=SBAR*CBAR*E3
      KM(6,3)=KM(3,6)
      RETURN
      END
