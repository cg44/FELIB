      SUBROUTINE NUMIN3(SAMP,ISAMP,WT,NGP)
C
C      THIS SUBROUTINE FORMS THE SAMPLING POINTS AND
C      WEIGHTS FOR INTEGRATION OVER A TETRAHEDRON
C
      REAL SAMP(ISAMP,*),WT(*)
      IF(NGP.EQ.1)GOTO 10
      IF(NGP.EQ.4)GOTO 40
      IF(NGP.EQ.5)GOTO 50
   10 SAMP(1,1)=.25
      SAMP(1,2)=.25
      SAMP(1,3)=.25
      WT(1)=1.
      GOTO 99
   40 SAMP(1,1)=.58541020
      SAMP(1,2)=.13819660
      SAMP(1,3)=SAMP(1,2)
      SAMP(2,2)=SAMP(1,1)
      SAMP(2,3)=SAMP(1,2)
      SAMP(2,1)=SAMP(1,2)
      SAMP(3,3)=SAMP(1,1)
      SAMP(3,1)=SAMP(1,2)
      SAMP(3,2)=SAMP(1,2)
      SAMP(4,1)=SAMP(1,2)
      SAMP(4,2)=SAMP(1,2)
      SAMP(4,3)=SAMP(1,2)
      WT(1)=.25
      WT(2)=.25
      WT(3)=.25
      WT(4)=.25
      GOTO 99
   50 SAMP(1,1)=.25
      SAMP(1,2)=.25
      SAMP(1,3)=.25
      SAMP(2,1)=.5
      SAMP(2,2)=1./6.
      SAMP(2,3)=SAMP(2,2)
      SAMP(3,2)=.5
      SAMP(3,3)=1./6.
      SAMP(3,1)=SAMP(3,3)
      SAMP(4,3)=.5
      SAMP(4,1)=1./6.
      SAMP(4,2)=SAMP(4,1)
      SAMP(5,1)=1./6.
      SAMP(5,2)=SAMP(5,1)
      SAMP(5,3)=SAMP(5,1)
      WT(1)=-.8
      WT(2)=9./20.
      WT(3)=WT(2)
      WT(4)=WT(2)
      WT(5)=WT(2)
   99 CONTINUE
      RETURN
      END
