      SUBROUTINE GEOV4Y(IP,IQ,NDE,RAD,DEP,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR EACH ELEMENT (NUMBERING IN THE Y-DIRECTION)
C
      REAL COORD(ICOORD,*),RAD(*),DEP(*)
      INTEGER G(*),NF(INF,*),NUM(4)
      NUM(1)=(IP-1)*(NDE+1)+IQ+1
      NUM(2)=NUM(1)-1
      NUM(3)=IP*(NDE+1)+IQ
      NUM(4)=NUM(3)+1
      INC=0
      DO 1 I=1,4
      DO 1 J=1,2
      INC=INC+1
    1 G(INC)=NF(NUM(I),J)
      COORD(1,1)=RAD(IP)
      COORD(2,1)=RAD(IP)
      COORD(3,1)=RAD(IP+1)
      COORD(4,1)=RAD(IP+1)
      COORD(1,2)=DEP(IQ+1)
      COORD(2,2)=DEP(IQ)
      COORD(3,2)=DEP(IQ)
      COORD(4,2)=DEP(IQ+1)
      RETURN
      END