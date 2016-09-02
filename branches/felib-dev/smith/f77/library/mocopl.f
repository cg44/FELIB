      SUBROUTINE MOCOPL(PHI,PSI,E,V,STRESS,PL)
C
C      THIS SUBROUTINE FORMS THE PLASTIC STRESS/STRAIN MATRIX
C      FOR A MOHR-COULOMB MATERIAL  (PHI,PSI IN DEGREES)
C
      REAL STRESS(4),ROW(4),COL(4),PL(4,4)
      SX=STRESS(1)
      SY=STRESS(2)
      TXY=STRESS(3)
      SZ=STRESS(4)
      PI=4.*ATAN(1.)
      PHIR=PHI*PI/180.
      PSIR=PSI*PI/180.
      SNPH=SIN(PHIR)
      SNPS=SIN(PSIR)
      SQ3=SQRT(3.)
      CC=1.-2.*V
      DX=(2.*SX-SY-SZ)/3.
      DY=(2.*SY-SZ-SX)/3.
      DZ=(2.*SZ-SX-SY)/3.
      D2=SQRT(-DX*DY-DY*DZ-DZ*DX+TXY*TXY)
      D3=DX*DY*DZ-DZ*TXY*TXY
      TH=-3.*SQ3*D3/(2.*D2**3)
      IF(TH.GT.1.)TH=1.
      IF(TH.LT.-1.)TH=-1.
      TH=ASIN(TH)/3.
      SNTH=SIN(TH)
      IF(ABS(SNTH).GT..49)THEN
      SIG=-1.
      IF(SNTH.LT.0.)SIG=1.
      RPH=SNPH*(1.+V)/3.
      RPS=SNPS*(1.+V)/3.
      CPS=.25*SQ3/D2*(1.+SIG*SNPS/3.)
      CPH=.25*SQ3/D2*(1.+SIG*SNPH/3.)
      COL(1)=RPH+CPH*((1.-V)*DX+V*(DY+DZ))
      COL(2)=RPH+CPH*((1.-V)*DY+V*(DZ+DX))
      COL(3)=CPH*CC*TXY
      COL(4)=RPH+CPH*((1.-V)*DZ+V*(DX+DY))
      ROW(1)=RPS+CPS*((1.-V)*DX+V*(DY+DZ))
      ROW(2)=RPS+CPS*((1.-V)*DY+V*(DZ+DX))
      ROW(3)=CPS*CC*TXY
      ROW(4)=RPS+CPS*((1.-V)*DZ+V*(DX+DY))
      EE=E/((1.+V)*CC*(RPH*SNPS+2.*CPH*CPS*D2*D2*CC))
      ELSE
      ALP=ATAN(ABS((SX-SY)/(2.*TXY)))
      CA=COS(ALP)
      SA=SIN(ALP)
      DD=CC*SA
      S1=1.
      S2=1.
      IF((SX-SY).LT..0)S1=-1.
      IF(TXY.LT..0)S2=-1.
      COL(1)=SNPH+S1*DD
      COL(2)=SNPH-S1*DD
      COL(3)=S2*CC*CA
      COL(4)=2.*V*SNPH
      ROW(1)=SNPS+S1*DD
      ROW(2)=SNPS-S1*DD
      ROW(3)=S2*CC*CA
      ROW(4)=2.*V*SNPS
      EE=E/(2.*(1.+V)*CC*(SNPH*SNPS+CC))
      END IF
      DO 1 I=1,4
      DO 1 J=1,4
    1 PL(I,J)=EE*ROW(I)*COL(J)
      RETURN
      END
