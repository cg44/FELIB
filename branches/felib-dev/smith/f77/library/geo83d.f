      SUBROUTINE GEO83D(IP,IQ,IS,NXE,NZE,AA,BB,CC,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE BRICK ELEMENTS COUNTING X-Z PLANES IN Y-DIRECTION
C
      REAL COORD(ICOORD,*)
      INTEGER G(*),NF(INF,*),NUM(8)
      NUM(1)=(IQ-1)*(NXE+1)*(NZE+1)+IS*(NXE+1)+IP
      NUM(2)=NUM(1)-NXE-1
      NUM(3)=NUM(2)+1
      NUM(4)=NUM(1)+1
      NUM(5)=NUM(1)+(NXE+1)*(NZE+1)
      NUM(6)=NUM(5)-NXE-1
      NUM(7)=NUM(6)+1
      NUM(8)=NUM(5)+1
      INC=0
      DO 1 I=1,8
      DO 1 J=1,3
      INC=INC+1
    1 G(INC)=NF(NUM(I),J)
      COORD(1,1)=(IP-1)*AA
      COORD(2,1)=(IP-1)*AA
      COORD(5,1)=(IP-1)*AA
      COORD(6,1)=(IP-1)*AA
      COORD(3,1)=IP*AA
      COORD(4,1)=IP*AA
      COORD(7,1)=IP*AA
      COORD(8,1)=IP*AA
      COORD(1,2)=(IQ-1)*BB
      COORD(2,2)=(IQ-1)*BB
      COORD(3,2)=(IQ-1)*BB
      COORD(4,2)=(IQ-1)*BB
      COORD(5,2)=IQ*BB
      COORD(6,2)=IQ*BB
      COORD(7,2)=IQ*BB
      COORD(8,2)=IQ*BB
      COORD(1,3)=-IS*CC
      COORD(4,3)=-IS*CC
      COORD(5,3)=-IS*CC
      COORD(8,3)=-IS*CC
      COORD(2,3)=-(IS-1)*CC
      COORD(3,3)=-(IS-1)*CC
      COORD(6,3)=-(IS-1)*CC
      COORD(7,3)=-(IS-1)*CC
      RETURN
      END
