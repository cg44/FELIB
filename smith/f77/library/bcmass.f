      SUBROUTINE BCMASS(MM,RHO,AREA,IP,COORD,ICOORD)
C
C      THIS SUBROUTINE FORMS THE CONSISTENT MASS MATRIX OF
C      AN INCLINED 2-D BEAM-COLUMN ELEMENT
C
      REAL MM(6,6),COORD(ICOORD,*)
      X1=COORD(IP,1)
      Y1=COORD(IP,2)
      X2=COORD(IP,3)
      Y2=COORD(IP,4)
      ELL=SQRT((Y2-Y1)**2+(X2-X1)**2)
      C=(X2-X1)/ELL
      S=(Y2-Y1)/ELL
      MM(1,1)=140.*C*C+156.*S*S
      MM(4,4)=MM(1,1)
      MM(2,2)=156.*C*C+140.*S*S
      MM(5,5)=MM(2,2)
      MM(3,3)=4.*ELL*ELL
      MM(6,6)=MM(3,3)
      MM(1,2)=-16.*C*S
      MM(4,5)=MM(1,2)
      MM(1,5)=-MM(1,2)
      MM(2,4)=-MM(1,2)
      MM(1,3)=-22.*ELL*S
      MM(4,6)=-MM(1,3)
      MM(2,3)=22.*ELL*C
      MM(5,6)=-MM(2,3)
      MM(1,4)=70.*C*C+54.*S*S
      MM(1,6)=13.*ELL*S
      MM(3,4)=-MM(1,6)
      MM(2,5)=54.*C*C+70.*S*S
      MM(2,6)=-13.*ELL*C
      MM(3,5)=-MM(2,6)
      MM(3,6)=-3.*ELL*ELL
      FAC=RHO*AREA*ELL/420.
      DO 1 I=1,6
      DO 1 J=I,6
      MM(I,J)=MM(I,J)*FAC
    1 MM(J,I)=MM(I,J)
      RETURN
      END
