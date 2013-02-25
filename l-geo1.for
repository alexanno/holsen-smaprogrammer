      
      PROGRAM LGEO1
          
      
C*******************************************************************
C     WRITE(*,*) BEREGNING AV GEOGRAFISKE KOORDINATER FOR PUNKT 2
C     WRITE(*,*) OG ASIMUT FOR LINJEN I PUNKT 2.KOORDINATENE
C     WRITE(*,*) (B1,L1),GEODETISK LINJE S12 OG ASIMUT A12 GITT.
C     WRITE(*,*) GJELDER OGSè FOR MEGET LANGE AVSTANDER
C     WRITE(*,*)       J.H. AUGUST 90
C*******************************************************************
        
      IMPLICIT LOGICAL(A-Z)
      INTEGER P,F

      REAL*8 A,B,B1,L1,B2,L2,A1,A2,DS,S1,S2,S3,RB0,RB1,RB2,SI1,SI2
      REAL*8 DSI,T,N,E,PI,RO,W0,DS0,DL,DLA,N1,R,R1,R2,R3,B0
      REAL*8 C,D1,D2,D3,VINKELTR,K1,LA1,LA2                                     

      PI=DATAN(1)*4.D0
      RO=PI/180.D0
            
      WRITE(*,*)'LEGG INN B1,L1,S12,A12:'                      
      READ(*,*) B1,L1,DS,A1 
      B1=B1*RO
      L1=L1*RO
      A1=A1*RO

      WRITE(*,*) 'VALG AV ELLIPSOIDE.P SETTES LIK 1 FOR :'
      WRITE(*,*) 'NORSK BESSELS ,LIK 2 FOR INTERNASJONAL:'
      WRITE(*,*) 'OG LIK 3 FOR FOR WGS-84-ELLIPSOIDEN:'
      WRITE(*,*) 'LEGG INN P:'
      READ(*,*) P
      
                 CALL AKSER(P,A,B)

C****** BEREGNING AV SIGMA1,BETA0 OG LAMDA1

      RB1=VINKELTR(B,A,B1)
     
                
      T= -DCOS(A1)
      N= DTAN(RB1)
         CALL VINKEL(SI1,T,N)
                             
      T= 1.D0
      N=-DSIN(SI1)*DTAN(A1)
         CALL VINKEL(RB0,T,N)
      
      IF(A1.GT.PI) THEN
         RB0=PI-RB0
      ENDIF
                  
      T=DTAN(SI1)
      N=DCOS(RB0)
         CALL VINKEL(LA1,T,N)
      
      

C******BEREGNING AV S1, S2=S1+DS OG SIGMA2

      B0=VINKELTR(A,B,RB0)
      E=(A-B)*(A+B)/A**2.D0
      W0=SQRT(1-E*(DSIN(B0))**2.D0)
           
      K1=(1-W0)/(1+W0)
      C=B*(1.D0+K1**2.D0/4.D0)/(1.D0-K1)
      D1=(K1/2.D0-(3.D0*K1**3.D0)/16.D0)
      D2=(K1**2.D0)/16.D0
      D3=(K1**3.D0)/48.D0
           
             CALL AVSTAND(C,D1,D2,D3,S1,SI1)
                  
      S2=S1+DS
      SI2=0
      S3 =0
      F=0
    
   10 DS0=S2-S3
      SI2=SI2+DS0/C
      F=F+1
                            
      CALL AVSTAND(C,D1,D2,D3,S3,SI2)
          
      IF(F.LT.5) THEN    
          
          GOTO 10
      
      ENDIF
     
     
C******BEREGNING AV LAMDA2,ASIMUT 2 OG BREDDENE BETA2 OG B2

                                 
      T=DTAN(SI2)
      N=DCOS(RB0)
          CALL VINKEL(LA2,T,N)
      T= 1.D0
      N=DSIN(SI2)*DTAN(RB0)
          CALL VINKEL(A2,T,N)
      T=DCOS(A2)
      N=DTAN(SI2)
          CALL VINKEL(RB2,T,N)
                  
      B2=VINKELTR(A,B,RB2)
      IF(B2.LT.0) THEN
          LA2=LA2-PI
      ELSEIF(B1.LT.0) THEN
          LA2=LA2+PI
            
      ENDIF
      
          
      IF(LA2.LT.LA1) THEN
           LA2=LA2+2.D0*PI
      ENDIF

C******BEREGNING AV LENGDEFORSKJELLEN,DL,FRA 1 TIL 2

      DLA=LA2-LA1
      DSI=SI2-SI1
      N1=(A-B)/(A+B)      
      R=E*DCOS(RB0)/2.D0
      R1=(1.D0+N1-K1/2.D0-(K1**2.D0)/4.D0)
      R2=K1/4.D0
      R3=(K1**2.D0)/16.D0
      DL=DLA-R*(R1*DSI-R2*(DSIN(2.D0*SI2)-DSIN(2.D0*SI1)))
      DL=DL-R*R3*(DSIN(4*SI2)-DSIN(4*SI1))
      
      IF(DL.GT.PI*2.D0) THEN
            DL=DL-PI*2.D0
      ENDIF
            
      
      IF(A1.GT.PI) THEN
          L2=L1-DL
      ELSE
          L2=L1+DL
      ENDIF
          
      IF(L2.GE.2.D0*PI) THEN
          L2=L2-2.D0*PI
      ENDIF
   
      
      
      IF(A1.LT.PI)THEN
          A2=2.D0*PI-A2
      ELSE 
          A2=A2
      ENDIF                       
                 
      IF(B1.LT.0.AND.B2.LT.0) THEN
          IF(A1.LT.PI) THEN
          L2=L2-PI
          ENDIF
      ENDIF
      
      IF(B1.LT.0.AND.B2.LT.0) THEN
          IF(A1.GT.PI) THEN
          L2=L2+PI
          ENDIF
      ENDIF
      
      IF(L2.GT.PI) THEN
          L2=L2-PI-PI
      ELSEIF(L2.LT.-PI) THEN
          L2=L2+PI+PI
      ENDIF

      RB0=RB0/RO
      SI1=SI1/RO
      SI2=SI2/RO
      LA1=LA1/RO
      LA2=LA2/RO
      WRITE(*,*) RB0
      WRITE(*,*) SI1,SI2
      WRITE(*,*) LA1,LA2
      
      B2=B2/RO
      L2=L2/RO
      A2=A2/RO
      
      
      WRITE(*,*) 'GEOGRAFISKE KOORDINATER:'
      WRITE(*,110) 'B2:',B2
      WRITE(*,*) 'NEGATIV B2 BETYR SYDLIG BREDDE:'
      WRITE(*,*)
      
      WRITE(*,110) 'L2:',L2
      WRITE(*,*)
      WRITE(*,*) 'ASIMUT FRA PUNKT 2 TIL PUNKT 1:'
      WRITE(*,110) 'A2:',A2
      WRITE(*,*) 'PROGRAMMET REGNER A2 I FORHOLD TIL NORD,OGSè Pè'
      WRITE(*,*) 'Pè DEN SYDLIGE HALVELLIPSOIDE'
  
  100 FORMAT(1X,A,$)
  110 FORMAT(1X,A,F14.9)
      
      STOP

      END
      
      
      
      REAL*8 FUNCTION VINKELTR(A,B,V)
      REAL*8 A,B,V

             VINKELTR=DATAN(A*TAN(V)/B)
      
      RETURN
      
      END


                SUBROUTINE VINKEL(A1,T,N)
      REAL*8 A1,T,N,PI
      
           PI=DATAN(1)*4.D0
      IF(ABS(T).LT.5D-9.AND.N.GT.0)THEN
             A1=0
             RETURN
      ELSEIF(T.GT.0.AND.ABS(N).LT.5D-9)THEN
             A1=PI/2.D0
             RETURN
      ELSEIF(T.LT.0.AND.ABS(N).LT.5D-9)THEN
             A1=1.5D0*PI

             RETURN
      ELSEIF(ABS(T).LT.5D-9.AND.N.LT.0)THEN
             A1=PI
             RETURN
      ENDIF
             
      A1=DATAN(T/N)
      
      IF(T.GT.0.AND.N.GT.0)THEN
         A1=A1            
      RETURN
      ELSEIF(T.GT.0.AND.N.LT.0)THEN
      
         A1=A1+PI
      RETURN
      ELSEIF(T.LT.0.AND.N.LT.0)THEN
         A1=A1+PI
      RETURN
                     
      ELSE
         A1=A1+PI*2.D0
      RETURN
      
      ENDIF
               
      END


                    SUBROUTINE AKSER(P,A,B)


      INTEGER P
      REAL*8 A,B      
                           
      IF(P.EQ.1)THEN
        A=6377492.018D0
        B=6356173.509D0
                RETURN
      
                ELSEIF(P.EQ.2)THEN
        A =6378388.000D0
        B =6356911.946D0    
                RETURN
  
                ELSEIF(P.EQ.3)THEN
        A=6378137.000D0
        B=6356752.314D0
                RETURN
      
      ENDIF
    
      END


                SUBROUTINE AVSTAND(C,D1,D2,D3,S1,SI)

      REAL*8 C,D1,D2,D3,S1,SI

      S1=C*(SI+D1*DSIN(2.D0*SI)-D2*DSIN(4.D0*SI)+D3*DSIN(6*SI))

      RETURN
      END