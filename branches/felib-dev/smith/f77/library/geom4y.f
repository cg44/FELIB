      SUBROUTINE GEOM4Y(IP,IQ,NYE,AA,BB,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 4-NODE QUADS COUNTING IN Y-DIRECTION
C
      REAL COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(4)
      NUM(1)=(IP-1)*(NYE+1)+IQ+1
      NUM(2)=NUM(1)-1
      NUM(3)=IP*(NYE+1)+IQ
      NUM(4)=NUM(3)+1
      INC=0
      DO 1 I=1,4
      DO 1 J=1,2
      INC=INC+1
    1 G(INC)=NF(NUM(I),J)
      COORD(1,1)=AA*(IP-1)
      COORD(2,1)=AA*(IP-1)
      COORD(3,1)=AA*IP
      COORD(4,1)=AA*IP
      COORD(1,2)=-BB*IQ
      COORD(2,2)=-BB*(IQ-1)
      COORD(3,2)=-BB*(IQ-1)
      COORD(4,2)=-BB*IQ
      RETURN
      END
