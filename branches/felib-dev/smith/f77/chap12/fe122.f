      PROGRAM P122                                                              
C                                                                               
C      PROGRAM 12.2 ANALYSIS OF PILE GROUPS USING 2-NODE LINE ELEMENTS          
C      FOR THE PILES, T-Z SPRINGS FOR THE SOIL AND MINDLIN'S EQUATION           
C      FOR THE INTERACTIONS                                                     
C                                                                               
C      ALTER NEXT LINE TO CHANGE PROBLEM SIZE                                   
C                                                                               
      PARAMETER(IKV=4950,ILOADS=300,INO=25,INX=31)                              
C                                                                               
      REAL KV(IKV),SM(IKV),LOADS(ILOADS),FORCE(ILOADS),KM(2,2),VAL(INO),        
     +COORD(INO,2),ZCOORD(INX),FMAX(INX),FSOIL(INO,INX),SUM(INO),KS(2,2)        
      INTEGER G(2),GP(INX),GQ(INX),NO(INO),NA(ILOADS)                           
      DATA IKM,IKS,IDOF/3*2/,RF/0.9/                                            
C                                                                               
C      INPUT AND INITIALISATION                                                 
C                                                                               
      READ(5,*)NXE,E,PILEN,R0,GSURF,GRATE,VSOIL,RHO,NPILE,INCS,IGROUP           
      READ(5,*)((COORD(I,J),J=1,2),I=1,NPILE)                                   
      READ(5,*)(FMAX(I),I=1,NXE+1)                                              
      READ(5,*)NL,(NO(I),VAL(I),I=1,NL)                                         
      NN=NXE+1                                                                  
      N=NN*NPILE                                                                
      ELL=PILEN/FLOAT(NXE)                                                      
      PI=ACOS(-1.)                                                              
      CSA=PI*R0*R0                                                              
      RM=2.5*RHO*(1.-VSOIL)*PILEN                                               
      DO 10 I=1,NN                                                              
   10 ZCOORD(I)=FLOAT(I-1)*ELL                                                  
      CALL NULL(FSOIL,INO,NPILE,NN)                                             
      FACT2=3.-4.*VSOIL                                                         
      FACT3=(1.-VSOIL)**2                                                       
C                                                                               
C      LOCATION OF DIAGONAL ELEMENTS OF STIFFNESS MATRIX                        
C                                                                               
      NSTO=0                                                                    
      DO 20 J=1,N                                                               
      NSTO=NSTO+J                                                               
   20 NA(J)=NSTO                                                                
      IR=NA(N)                                                                  
C                                                                               
C      LOAD/DISPLACEMENTS APPLIED IN INCREMENTS                                 
C                                                                               
      DO 30 II=1,INCS                                                           
      WRITE(6,2000)II                                                           
      CALL NULVEC(KV,IR)                                                        
      CALL NULVEC(SM,IR)                                                        
      IF(IGROUP.EQ.0)GOTO 40                                                    
C                                                                               
C      FORM SOIL FLEXIBILITY MATRIX USING MINDLIN'S EQUATION                    
C                                                                               
      DO 50 IP=1,NPILE                                                          
      CALL GEOM(IP,NN,NPILE,GP)                                                 
      IF(IP+1.GT.NPILE)GOTO 50                                                  
      DO 60 IQ=IP+1,NPILE                                                       
      CALL GEOM(IQ,NN,NPILE,GQ)                                                 
      RR=SQRT((COORD(IP,1)-COORD(IQ,1))**2 +                                    
     +        (COORD(IP,2)-COORD(IQ,2))**2)                                     
      DO 70 I=1,NN                                                              
      IF(FSOIL(IP,I)/FMAX(I).GT..9999)GOTO 70                                   
      ZZ=ZCOORD(I)                                                              
      GSOILI=GSURF+GRATE*ZZ                                                     
      IF(I.EQ.1)GSOILI=GSURF+GRATE*.25*ELL                                      
      IF(I.EQ.NXE)GSOILI=GSURF+GRATE*(ZZ+.25*ELL)                               
      DO 80 J=1,NN                                                              
      IF(FSOIL(IQ,J)/FMAX(J).GT..9999)GOTO 80                                   
      CC=ZCOORD(J)                                                              
      GSOILJ=GSURF+GRATE*CC                                                     
      IF(J.EQ.1)GSOILJ=GSURF+GRATE*.25*ELL                                      
      IF(J.EQ.NXE)GSOILJ=GSURF+GRATE*(CC+.25*ELL)                               
      GSOIL=.5*(GSOILI+GSOILJ)                                                  
      FACT1=1./(16.*PI*GSOIL*(1.-VSOIL))                                        
      ZMC2=(ZZ-CC)*(ZZ-CC)                                                      
      ZPC2=(ZZ+CC)*(ZZ+CC)                                                      
      R1=SQRT(RR*RR+ZMC2)                                                       
      R2=SQRT(RR*RR+ZPC2)                                                       
C                                                                               
C      SOIL DISPLACEMENT AT IP DUE TO UNIT LOAD AT IQ                           
C      OBTAINED USING MINDLIN'S EQUATION                                        
C                                                                               
      KOUNT=MAX0(NA(GQ(J))-GQ(J)+GP(I),NA(GP(I))-GP(I)+GQ(J))                   
      KV(KOUNT)=FACT1*(FACT2/R1+(8.*FACT3-FACT2)/R2+ZMC2/R1**3                  
     +        +(FACT2*ZPC2-2.*CC*ZZ)/R2**3+6.*ZZ*CC*ZPC2/R2**5)                 
   80 CONTINUE                                                                  
   70 CONTINUE                                                                  
   60 CONTINUE                                                                  
   50 CONTINUE                                                                  
C                                                                               
C      ADD IN FLEXIBILITY CONTRIBUTIONS FROM DISCRETE SOIL SPRINGS              
C                                                                               
   40 DO 90 IP=1,NPILE                                                          
      CALL GEOM(IP,NN,NPILE,GP)                                                 
      GSOIL=GSURF+GRATE*.25*ELL                                                 
      CALL FORMXI(FSOIL(IP,1),FMAX(1),RF,RM,R0,XI)                              
      FSHAFT=XI/(PI*GSOIL*ELL)                                                  
      IF(FSOIL(IP,1)/FMAX(1).GT.0.9999)FSHAFT=1.E12                             
      KV(NA(GP(1)))=KV(NA(GP(1)))+FSHAFT                                        
      GSOIL=GSURF+GRATE*(ZCOORD(NXE)+.25*ELL)                                   
      CALL FORMXI(FSOIL(IP,NXE),FMAX(NXE),RF,RM,R0,XI)                          
      FSHAFT=XI/(PI*GSOIL*3.*ELL)                                               
      IF(FSOIL(IP,NXE)/FMAX(NXE).GT..9999)FSHAFT=1.E12                          
      KV(NA(GP(NXE)))=KV(NA(GP(NXE)))+FSHAFT                                    
      DO 100 I=2,NN-2                                                           
      GSOIL=GSURF+GRATE*ZCOORD(I)                                               
      CALL FORMXI(FSOIL(IP,I),FMAX(I),RF,RM,R0,XI)                              
      FSHAFT=XI/(2.*PI*GSOIL*ELL)                                               
      IF(FSOIL(IP,I)/FMAX(I).GT..9999)FSHAFT=1.E12                              
  100 KV(NA(GP(I)))=KV(NA(GP(I)))+FSHAFT                                        
C                                                                               
C      ADD TIP FLEXIBILITY                                                      
C                                                                               
      GSOIL=GSURF+GRATE*ZCOORD(NN)                                              
      FTIP=(1.-VSOIL)/((4.*GSOIL*R0)*(1.-RF*FSOIL(IP,NN)/FMAX(NN))**2)          
      IF(FSOIL(IP,NN)/FMAX(NN).GT..9999)FTIP=1.E12                              
   90 KV(NA(GP(NN)))=KV(NA(GP(NN)))+FTIP                                        
C                                                                               
C      INVERT SOIL FLEXIBILITY MATRIX TO GIVE STIFFNESS MATRIX                  
C                                                                               
      CALL SPARIN(KV,N,NA)                                                      
      DO 110 I=1,N                                                              
      CALL NULVEC(LOADS,N)                                                      
      LOADS(I)=1.                                                               
      CALL SPABAC(KV,LOADS,N,NA)                                                
      DO 120 J=1,I                                                              
  120 SM(NA(I)-I+J)=LOADS(J)                                                    
  110 CONTINUE                                                                  
      CALL VECCOP(SM,KV,IR)                                                     
C                                                                               
C      FORM AND ASSEMBLE STIFFNESS MATRIX                                       
C                                                                               
      CALL AXIKM(KM,CSA,E,ELL)                                                  
      DO 130 IP=1,NPILE                                                         
      CALL GEOM(IP,NN,NPILE,GP)                                                 
      DO 140 I=1,NXE                                                            
      G(1)=GP(I)                                                                
      G(2)=GP(I+1)                                                              
  140 CALL FSPARV(SM,KM,IKM,G,NA,IDOF)                                          
  130 CONTINUE                                                                  
C                                                                               
C      PRESCRIBE DISPLACEMENTS : BIG SPRING TECHNIQUE                           
C                                                                               
      CALL NULVEC(LOADS,N)                                                      
      DO 150 I=1,NL                                                             
      IX=NA(NO(I))                                                              
      SM(IX)=SM(IX)+1.E20                                                       
  150 LOADS(NO(I))=VAL(I)*SM(IX)                                                
C                                                                               
C      SOLVE SM TO GIVE PILE DISPLACEMENTS                                      
C                                                                               
      CALL SPARIN(SM,N,NA)                                                      
      CALL SPABAC(SM,LOADS,N,NA)                                                
C                                                                               
C      RECOVER SOIL SPRING FORCES                                               
C                                                                               
      CALL LINMLS(KV,LOADS,FORCE,N,NA)                                          
      TOTSUM=0.                                                                 
      DO 160 IP=1,NPILE                                                         
      CALL GEOM(IP,NN,NPILE,GP)                                                 
      SUM(IP)=0.                                                                
      DO 170 I=1,NN                                                             
      FSOIL(IP,I)=FSOIL(IP,I)+FORCE(GP(I))                                      
      IF(FSOIL(IP,I).GT.FMAX(I))FSOIL(IP,I)=FMAX(I)                             
  170 SUM(IP)=SUM(IP)+FSOIL(IP,I)                                               
  160 TOTSUM=TOTSUM+SUM(IP)                                                     
      WRITE(6,1000)TOTSUM,VAL(1)*II                                             
      PAV=TOTSUM/FLOAT(NPILE)                                                   
      DO 180 IP=1,NPILE                                                         
  180 WRITE(6,2000)IP,SUM(IP)/PAV                                               
   30 CONTINUE                                                                  
 1000 FORMAT(2E12.4)                                                            
 2000 FORMAT(I10,E12.4)                                                         
      STOP                                                                      
      END                                                                       
