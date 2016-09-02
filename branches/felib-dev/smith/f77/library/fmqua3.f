      SUBROUTINE FMQUA3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C      THIS SUBROUTINE FORMS THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES FOR 20-NODED BRICK ELEMENTS
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      INTEGER XII(20),ETAI(20),ZETAI(20)
      XI=SAMP(I,1)
      ETA=SAMP(J,1)
      ZETA=SAMP(K,1)
      XII(1)=-1
      XII(2)=-1
      XII(3)=-1
      XII(9)=-1
      XII(10)=-1
      XII(13)=-1
      XII(14)=-1
      XII(15)=-1
      XII(4)=0
      XII(8)=0
      XII(16)=0
      XII(20)=0
      XII(5)=1
      XII(6)=1
      XII(7)=1
      XII(11)=1
      XII(12)=1
      XII(17)=1
      XII(18)=1
      XII(19)=1
      DO 1 L=1,8
    1 ETAI(L)=-1
      DO 2 L=9,12
    2 ETAI(L)=0
      DO 3 L=13,20
    3 ETAI(L)=1
      ZETAI(1)=-1
      ZETAI(7)=-1
      ZETAI(8)=-1
      ZETAI(9)=-1
      ZETAI(12)=-1
      ZETAI(13)=-1
      ZETAI(19)=-1
      ZETAI(20)=-1
      ZETAI(2)=0
      ZETAI(6)=0
      ZETAI(14)=0
      ZETAI(18)=0
      ZETAI(3)=1
      ZETAI(4)=1
      ZETAI(5)=1
      ZETAI(10)=1
      ZETAI(11)=1
      ZETAI(15)=1
      ZETAI(16)=1
      ZETAI(17)=1
      DO 4 L=1,20
      XIO=XI*XII(L)
      ETAO=ETA*ETAI(L)
      ZETAO=ZETA*ZETAI(L)
      IF(L.EQ.4.OR.L.EQ.8.OR.L.EQ.16.OR.L.EQ.20)THEN
      FUN(L)=.25*(1.-XI*XI)*(1.+ETAO)*(1.+ZETAO)
      DER(1,L)=-.5*XI*(1.+ETAO)*(1.+ZETAO)
      DER(2,L)=.25*ETAI(L)*(1.-XI*XI)*(1.+ZETAO)
      DER(3,L)=.25*ZETAI(L)*(1.-XI*XI)*(1.+ETAO)
      ELSE IF(L.GE.9.AND.L.LE.12)THEN
      FUN(L)=.25*(1.+XIO)*(1.-ETA*ETA)*(1.+ZETAO)
      DER(1,L)=.25*XII(L)*(1.-ETA*ETA)*(1.+ZETAO)
      DER(2,L)=-.5*ETA*(1.+XIO)*(1.+ZETAO)
      DER(3,L)=.25*ZETAI(L)*(1.+XIO)*(1.-ETA*ETA)
      ELSE IF(L.EQ.2.OR.L.EQ.6.OR.L.EQ.14.OR.L.EQ.18)THEN
      FUN(L)=.25*(1.+XIO)*(1.+ETAO)*(1.-ZETA*ZETA)
      DER(1,L)=.25*XII(L)*(1.+ETAO)*(1.-ZETA*ZETA)
      DER(2,L)=.25*ETAI(L)*(1.+XIO)*(1.-ZETA*ZETA)
      DER(3,L)=-.5*ZETA*(1.+XIO)*(1.+ETAO)
      ELSE
      FUN(L)=.125*(1.+XIO)*(1.+ETAO)*(1.+ZETAO)*(XIO+ETAO+ZETAO-2.)
      DER(1,L)=.125*XII(L)*(1.+ETAO)*(1.+ZETAO)*(2.*XIO+ETAO+ZETAO-1.)
      DER(2,L)=.125*ETAI(L)*(1.+XIO)*(1.+ZETAO)*(XIO+2.*ETAO+ZETAO-1.)
      DER(3,L)=.125*ZETAI(L)*(1.+XIO)*(1.+ETAO)*(XIO+ETAO+2.*ZETAO-1.)
      END IF
    4 CONTINUE
      RETURN
      END
