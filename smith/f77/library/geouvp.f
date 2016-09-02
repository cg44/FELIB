      SUBROUTINE GEOUVP(IP,IQ,NXE,WIDTH,DEPTH,COORD,ICOORD,
     +                  COORDF,ICORDF,G,NF,INF)
C
C      THIS SUBROUTINE FORMS THE NODAL COORDINATES AND STEERING
C      VECTOR FOR A VARIABLE MESH OF 4-NODE/8-NODE
C      QUADRILATERAL ELEMENTS NUMBERING IN THE X-DIRECTION
C      (U,V,P  BIOT CONSOLIDATION)
C
      REAL COORD(ICOORD,*),COORDF(ICORDF,*),WIDTH(*),DEPTH(*)
      INTEGER NUM(8),G(*),NF(INF,*)
      NUM(1)=IQ*(3*NXE+2)+2*IP-1
      NUM(2)=IQ*(3*NXE+2)+IP-NXE-1
      NUM(3)=(IQ-1)*(3*NXE+2)+2*IP-1
      NUM(4)=NUM(3)+1
      NUM(5)=NUM(4)+1
      NUM(6)=NUM(2)+1
      NUM(7)=NUM(1)+2
      NUM(8)=NUM(1)+1
      INC=0
      DO 1 I=1,8
      DO 1 J=1,2
      INC=INC+1
    1 G(INC)=NF(NUM(I),J)
      DO 2 I=1,7,2
      INC=INC+1
    2 G(INC)=NF(NUM(I),3)
      COORD(1,1)=WIDTH(IP)
      COORD(2,1)=WIDTH(IP)
      COORD(3,1)=WIDTH(IP)
      COORDF(1,1)=WIDTH(IP)
      COORDF(2,1)=WIDTH(IP)
      COORD(5,1)=WIDTH(IP+1)
      COORD(6,1)=WIDTH(IP+1)
      COORD(7,1)=WIDTH(IP+1)
      COORDF(3,1)=WIDTH(IP+1)
      COORDF(4,1)=WIDTH(IP+1)
      COORD(4,1)=.5*(COORD(3,1)+COORD(5,1))
      COORD(8,1)=.5*(COORD(7,1)+COORD(1,1))
      COORD(1,2)=DEPTH(IQ+1)
      COORD(8,2)=DEPTH(IQ+1)
      COORD(7,2)=DEPTH(IQ+1)
      COORDF(1,2)=DEPTH(IQ+1)
      COORDF(4,2)=DEPTH(IQ+1)
      COORD(3,2)=DEPTH(IQ)
      COORD(4,2)=DEPTH(IQ)
      COORD(5,2)=DEPTH(IQ)
      COORDF(2,2)=DEPTH(IQ)
      COORDF(3,2)=DEPTH(IQ)
      COORD(2,2)=.5*(COORD(1,2)+COORD(3,2))
      COORD(6,2)=.5*(COORD(5,2)+COORD(7,2))
      RETURN
      END
