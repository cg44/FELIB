      SUBROUTINE SLOGEO(IP,IQ,NYE,TOP,BOT,DEPTH,COORD,ICOORD,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE QUADRILATERALS IN A 'SLOPE' GEOMETRY
C      (NUMBERING IN THE Y-DIRECTION)
C
      REAL TOP(*),BOT(*),COORD(ICOORD,*),DEPTH(*)
      INTEGER G(*),NF(INF,*),NUM(8)
      NUM(1)=(IP-1)*(3*NYE+2)+2*IQ+1
      NUM(2)=NUM(1)-1
      NUM(3)=NUM(1)-2
      NUM(4)=(IP-1)*(3*NYE+2)+2*NYE+IQ+1
      NUM(5)=IP*(3*NYE+2)+2*IQ-1
      NUM(6)=NUM(5)+1
      NUM(7)=NUM(5)+2
      NUM(8)=NUM(4)+1
      INC=0
      DO 1 I=1,8
      DO 1 J=1,2
      INC=INC+1
    1 G(INC)=NF(NUM(I),J)
      FAC1=(BOT(IP)-TOP(IP))/(DEPTH(NYE+1)-DEPTH(1))
      FAC2=(BOT(IP+1)-TOP(IP+1))/(DEPTH(NYE+1)-DEPTH(1))
      COORD(1,1)=TOP(IP)+(DEPTH(IQ+1)-DEPTH(1))*FAC1
      COORD(3,1)=TOP(IP)+(DEPTH(IQ)-DEPTH(1))*FAC1
      COORD(5,1)=TOP(IP+1)+(DEPTH(IQ)-DEPTH(1))*FAC2
      COORD(7,1)=TOP(IP+1)+(DEPTH(IQ+1)-DEPTH(1))*FAC2
      COORD(2,1)=.5*(COORD(1,1)+COORD(3,1))
      COORD(6,1)=.5*(COORD(5,1)+COORD(7,1))
      COORD(4,1)=.5*(COORD(3,1)+COORD(5,1))
      COORD(8,1)=.5*(COORD(7,1)+COORD(1,1))
      COORD(1,2)=DEPTH(IQ+1)
      COORD(8,2)=DEPTH(IQ+1)
      COORD(7,2)=DEPTH(IQ+1)
      COORD(3,2)=DEPTH(IQ)
      COORD(4,2)=DEPTH(IQ)
      COORD(5,2)=DEPTH(IQ)
      COORD(2,2)=.5*(COORD(1,2)+COORD(3,2))
      COORD(6,2)=.5*(COORD(5,2)+COORD(7,2))
      RETURN
      END