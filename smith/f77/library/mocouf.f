      SUBROUTINE MOCOUF(PHI,C,SIGM,DSBAR,THETA,F)
C
C      THIS SUBROUTINE CALCULATES THE VALUE OF THE YIELD FUNCTION
C      FOR A MOHR-COULOMB MATERIAL (PHI IN DEGREES)
C
      PHIR=PHI*4.*ATAN(1.)/180.
      SNPH=SIN(PHIR)
      CSPH=COS(PHIR)
      CSTH=COS(THETA)
      SNTH=SIN(THETA)
      F=SNPH*SIGM+DSBAR*(CSTH/SQRT(3.)-SNTH*SNPH/3.)-C*CSPH
      RETURN
      END