      PROGRAM P100                                                              
C                                                                               
C       PROGRAM 10.0 EIGENVALUES AND EIGENVECTORS OF A 2-D FRAME                
C       OF BEAM-COLUMN ELEMENTS,CONSISTENT MASS.                                
C                                                                               
C                                                                               
C       ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                  
C                                                                               
      PARAMETER(IKB1=100,IKB2=20,ICOORD=20,INF=30)                              
C                                                                               
      REAL KB(IKB1,IKB2),MB(IKB1,IKB2),KM(6,6),MM(6,6),COORD(ICOORD,4),         
     +BIGK(IKB1,IKB1),DIAG(IKB1),UDIAG(IKB1)                                    
      INTEGER G(6),NF(INF,3)                                                    
      DATA IKM,IDOF/2*6/,NODOF/3/                                               
C                                                                               
C       INPUT SECTION                                                           
C                                                                               
      READ(5,*)NXE,N,IW,NN,NR,NMODES                                            
      CALL READNF(NF,INF,NN,NODOF,NR)                                           
      IWP1=IW+1                                                                 
      CALL NULL(KB,IKB1,N,IWP1)                                                 
      CALL NULL(MB,IKB1,N,IWP1)                                                 
      CALL NULL(BIGK,IKB1,N,N)                                                  
C                                                                               
C       GLOBAL STIFFNESS MATRIX ASSEMBLY                                        
C                                                                               
      DO 10 IP=1,NXE                                                            
      READ(5,*)EA,EI,RHO,AREA,(COORD(IP,I),I=1,4)                               
      CALL BMCOL2(KM,EA,EI,IP,COORD,ICOORD)                                     
      CALL BCMASS(MM,RHO,AREA,IP,COORD,ICOORD)                                  
      CALL GSTRNG(IP,NODOF,NF,INF,G)                                            
      CALL FORMKB(KB,IKB1,KM,IKM,G,IW,IDOF)                                     
   10 CALL FORMKB(MB,IKB1,MM,IKM,G,IW,IDOF)                                     
C                                                                               
C      REDUCE TO STANDARD EIGENVALUE PROBLEM                                    
C                                                                               
      CALL CHOLIN(MB,IKB1,N,IW)                                                 
      CALL LBKBAN(MB,BIGK,KB,IKB1,IW,N)                                         
      DO 20 I=1,N                                                               
      DO 20 J=I+1,N                                                             
   20 BIGK(I,J)=BIGK(J,I)                                                       
      CALL LBBT(MB,BIGK,IKB1,IW,N)                                              
C                                                                               
C      EXTRACT EIGENVALUES AND EIGENVECTORS                                     
C                                                                               
      CALL TRIDIA(N,1.E-20,BIGK,DIAG,UDIAG,BIGK,IKB1)                           
      IFAIL=1                                                                   
      CALL EVECTS(N,1.E-20,DIAG,UDIAG,BIGK,IKB1,IFAIL)                          
      CALL PRINTV(DIAG,N)                                                       
      DO 30 J=1,NMODES                                                          
      DO 40 I=1,N                                                               
   40 UDIAG(I)=BIGK(I,J)                                                        
      CALL CHOBK2(MB,IKB1,UDIAG,N,IW)                                           
   30 CALL PRINTV(UDIAG,N)                                                      
      STOP                                                                      
      END                                                                       
