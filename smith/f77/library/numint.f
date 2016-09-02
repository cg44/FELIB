      SUBROUTINE NUMINT(S,IS,WT,NGP)
C
C      THIS SUBROUTINE FORMS THE SAMPLING POINTS AND
C      WEIGHTS FOR INTEGRATION OVER A TRIANGULAR AREA
C
      REAL S(IS,*),WT(*)
      GOTO(1,1,3,4,4,6,7,7,7,7,7,12,12,12,12,16),NGP
    1 S(1,1)=1./3.
      S(1,2)=1./3.
      WT(1)=1.
      GOTO 99
    3 S(1,1)=.5
      S(1,2)=.5
      S(2,1)=.5
      S(2,2)=0.
      S(3,1)=0.
      S(3,2)=.5
      WT(1)=1./3.
      WT(2)=WT(1)
      WT(3)=WT(1)
      GOTO 99
    4 S(1,1)=1./3.
      S(1,2)=1./3.
      S(2,1)=.6
      S(2,2)=.2
      S(3,1)=.2
      S(3,2)=.6
      S(4,1)=.2
      S(4,2)=.2
      WT(1)=-9./16.
      WT(2)=25./48.
      WT(3)=WT(2)
      WT(4)=WT(2)
      GOTO 99
    6 S(1,1)=.816847572980459
      S(1,2)=.091576213509771
      S(2,1)=S(1,2)
      S(2,2)=S(1,1)
      S(3,1)=S(1,2)
      S(3,2)=S(1,2)
      S(4,1)=.108103018168070
      S(4,2)=.445948490915965
      S(5,1)=S(4,2)
      S(5,2)=S(4,1)
      S(6,1)=S(4,2)
      S(6,2)=S(4,2)
      WT(1)=.109951743655322
      WT(2)=WT(1)
      WT(3)=WT(1)
      WT(4)=.223381589678011
      WT(5)=WT(4)
      WT(6)=WT(4)
      GOTO 99
    7 S(1,1)=1./3.
      S(1,2)=1./3.
      S(2,1)=.797426985353087
      S(2,2)=.101286507323456
      S(3,1)=S(2,2)
      S(3,2)=S(2,1)
      S(4,1)=S(2,2)
      S(4,2)=S(2,2)
      S(5,1)=.470142064105115
      S(5,2)=.059715871789770
      S(6,1)=S(5,2)
      S(6,2)=S(5,1)
      S(7,1)=S(5,1)
      S(7,2)=S(5,1)
      WT(1)=.225
      WT(2)=.125939180544827
      WT(3)=WT(2)
      WT(4)=WT(2)
      WT(5)=.132394152788506
      WT(6)=WT(5)
      WT(7)=WT(5)
      GOTO 99
   12 S(1,1)=.873821971016996
      S(1,2)=.063089014491502
      S(2,1)=S(1,2)
      S(2,2)=S(1,1)
      S(3,1)=S(1,2)
      S(3,2)=S(1,2)
      S(4,1)=.501426509658179
      S(4,2)=.249286745170910
      S(5,1)=S(4,2)
      S(5,2)=S(4,1)
      S(6,1)=S(4,2)
      S(6,2)=S(4,2)
      S(7,1)=.636502499121399
      S(7,2)=.310352451033785
      S(8,1)=S(7,1)
      S(8,2)=.053145049844816
      S(9,1)=S(7,2)
      S(9,2)=S(7,1)
      S(10,1)=S(7,2)
      S(10,2)=S(8,2)
      S(11,1)=S(8,2)
      S(11,2)=S(7,1)
      S(12,1)=S(8,2)
      S(12,2)=S(7,2)
      WT(1)=.050844906370207
      WT(2)=WT(1)
      WT(3)=WT(1)
      WT(4)=.116786275726379
      WT(5)=WT(4)
      WT(6)=WT(4)
      WT(7)=.082851075618374
      WT(8)=WT(7)
      WT(9)=WT(7)
      WT(10)=WT(7)
      WT(11)=WT(7)
      WT(12)=WT(7)
      GOTO 99
   16 S(1,1)=1./3.
      S(1,2)=1./3.
      S(2,1)=.658861384496478
      S(2,2)=.170569307751761
      S(3,1)=S(2,2)
      S(3,2)=S(2,1)
      S(4,1)=S(2,2)
      S(4,2)=S(2,2)
      S(5,1)=.898905543365938
      S(5,2)=.050547228317031
      S(6,1)=S(5,2)
      S(6,2)=S(5,1)
      S(7,1)=S(5,2)
      S(7,2)=S(5,2)
      S(8,1)=.081414823414554
      S(8,2)=.459292588292723
      S(9,1)=S(8,2)
      S(9,2)=S(8,1)
      S(10,1)=S(8,2)
      S(10,2)=S(8,2)
      S(11,1)=.008394777409958
      S(11,2)=.263112829634638
      S(12,1)=S(11,1)
      S(12,2)=.728492392955404
      S(13,1)=S(11,2)
      S(13,2)=S(11,1)
      S(14,1)=S(11,2)
      S(14,2)=S(12,2)
      S(15,1)=S(12,2)
      S(15,2)=S(11,1)
      S(16,1)=S(12,2)
      S(16,2)=S(11,2)
      WT(1)=.144315607677787
      WT(2)=.103217370534718
      WT(3)=WT(2)
      WT(4)=WT(2)
      WT(5)=.032458497623198
      WT(6)=WT(5)
      WT(7)=WT(5)
      WT(8)=.095091634267284
      WT(9)=WT(8)
      WT(10)=WT(8)
      WT(11)=.027230314174435
      WT(12)=WT(11)
      WT(13)=WT(11)
      WT(14)=WT(11)
      WT(15)=WT(11)
      WT(16)=WT(11)
   99 CONTINUE
      RETURN
      END