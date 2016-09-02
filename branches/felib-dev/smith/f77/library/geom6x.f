      SUBROUTINE GEOM6X(IP,IQ,NXE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 6-NODE TRIANGLES COUNTING IN X-DIRECTION
C
      REAL COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(6)
      IF(MOD(IQ,2).EQ.0)GOTO 1
      NUM(1)=(IQ-1)*(2*NXE+1)+2*IP-1
      NUM(2)=(IQ-1)*(2*NXE+1)+2*NXE+2*IP
      NUM(3)=(IQ+1)*(2*NXE+1)+2*IP-1
      NUM(4)=NUM(2)+1
      NUM(5)=NUM(1)+2
      NUM(6)=NUM(1)+1
      COORD(1,1)=(IP-1)*AA
      COORD(1,2)=-(IQ-1)/2*BB
      COORD(3,1)=(IP-1)*AA
      COORD(3,2)=-(IQ+1)/2*BB
      COORD(5,1)=IP*AA
      COORD(5,2)=COORD(1,2)
      GOTO 2
    1 NUM(1)=IQ*(2*NXE+1)+2*IP+1
      NUM(2)=(IQ-2)*(2*NXE+1)+2*NXE+2*IP+2
      NUM(3)=(IQ-2)*(2*NXE+1)+2*IP+1
      NUM(4)=NUM(2)-1
      NUM(5)=NUM(1)-2
      NUM(6)=NUM(1)-1
      COORD(1,1)=IP*AA
      COORD(1,2)=-IQ/2*BB
      COORD(3,1)=IP*AA
      COORD(3,2)=-(IQ-2)/2*BB
      COORD(5,1)=(IP-1)*AA
      COORD(5,2)=COORD(1,2)
    2 DO 3 I=1,2
      COORD(2,I)=.5*(COORD(1,I)+COORD(3,I))
      COORD(4,I)=.5*(COORD(3,I)+COORD(5,I))
    3 COORD(6,I)=.5*(COORD(5,I)+COORD(1,I))
      INC=0
      DO 4 I=1,6
      DO 4 J=1,2
      INC=INC+1
    4 G(INC)=NF(NUM(I),J)
      RETURN
      END
