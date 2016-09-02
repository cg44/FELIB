      SUBROUTINE LOC2F(LOCAL,GLOBAL,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE TRANSFORMS THE END REACTIONS AND MOMENTS
C      INTO THE ELEMENT'S LOCAL COORDINATE SYSTEM (2-D)
C
      REAL COORD(ICOORD,*),LOCAL(*),GLOBAL(*)
      X1=COORD(IP,1)
      Y1=COORD(IP,2)
      X2=COORD(IP,3)
      Y2=COORD(IP,4)
      ELL=SQRT((X2-X1)**2+(Y2-Y1)**2)
      C=(X2-X1)/ELL
      S=(Y2-Y1)/ELL
      LOCAL(1)=C*GLOBAL(1)+S*GLOBAL(2)
      LOCAL(2)=C*GLOBAL(2)-S*GLOBAL(1)
      LOCAL(3)=GLOBAL(3)
      LOCAL(4)=C*GLOBAL(4)+S*GLOBAL(5)
      LOCAL(5)=C*GLOBAL(5)-S*GLOBAL(4)
      LOCAL(6)=GLOBAL(6)
      RETURN
      END
