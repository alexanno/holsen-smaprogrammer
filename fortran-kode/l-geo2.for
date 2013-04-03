                 
          
      PROGRAM LGEO2




C*******************************************************************
C     WRITE(*,*) GEOGRAFISKE KOORDINATER ER GITT FOR TO PUNKTER 
C     WRITE(*,*) 1 OG 2,(B1,L1) OG(B2,L2).BEREGN GEODETISK LINJE
C     WRITE(*,*) FRA 1 TIL 2 OG ASIMUT I PUNKT 1 OG I PUNKT 2.
C     WRITE(*,*) GJELDER OGS� FOR MEGET LANGE AVSTANDER
C     WRITE(*,*)       J.H. AUGUST 90
C*******************************************************************
        
      IMPLICIT LOGICAL(A-Z)
      INTEGER P,F

      REAL*8 A,B,B1,L1,B2,L2,A1,A2,DS,S1,S2,RB0,RB1,RB2,RB3,SI1
      REAL*8 SI2,DSI,E,PI,RO,W0,DL,DLA,N1,R,R1,R2,R3,B0
      REAL*8 C,D1,D2,D3,VINKELTR,K1,T,N,LA1,LA2,D                                     

      PI=DATAN(1.D0)*4.D0
      RO=PI/180.D0
            
      WRITE(*,*)'LEGG INN B1,L1,B2,L2 :'                    
      READ(*,*) B1,L1,B2,L2 
      B1=B1*RO
      L1=L1*RO
      B2=B2*RO
      L2=L2*RO 
                
      
      WRITE(*,*) 'VALG AV ELLIPSOIDE.P SETTES LIK 1 FOR :'
      WRITE(*,*) 'NORSK BESSELS ,LIK 2 FOR INTERNASJONAL:'
      WRITE(*,*) 'OG LIK 3 FOR FOR WGS-84-ELLIPSOIDEN:'
      WRITE(*,*) 'LEGG INN P:'
      READ(*,*) P
      
                 CALL AKSER(P,A,B)
         
      RB1=VINKELTR(B,A,B1)
      RB2=VINKELTR(B,A,B2)
      RB3=(RB1+RB2)/2.D0
      E=(A-B)*(A+B)/A**2.D0
      DL=L2-L1
     
      
      IF(ABS(DL).LT.2.E-8) THEN
             WRITE(*,*) 'BRUK HELLER ET MERIDIANBUEPROGRAM'
      
      D=PI-ABS(DL)
      GOTO 40      
      
      ENDIF
      
      IF(DL.LT.-PI) THEN
          DL=DL +PI*2.D0
      ELSEIF (DL.GT.PI) THEN
          DL=DL-PI*2.D0
      ENDIF                     
      
      DLA=DL/SQRT(1-E*(DCOS(RB3))**2.D0)
      
      F=0
      

C*****BEREGNING AV VERDIER FOR A1,RB0,SI1 OG SI2.

   11 T=DTAN(RB1)*DCOS(DLA)-DTAN(RB2)
      N=DTAN(RB1)*DSIN(DLA)
            CALL VINKEL(LA1,T,N)
      T=-DTAN(RB2)*COS(DLA)+DTAN(RB1) 
      N=DTAN(RB2)*SIN(DLA)    
            CALL VINKEL(LA2,T,N)
           
      
      T= DTAN(RB1)
      N=DCOS(LA1)
            CALL VINKEL(RB0,T,N)
              
      T=DCOS(RB0)*DTAN(LA1)
      N=1.D0
            CALL VINKEL(SI1,T,N)
      T= DCOS(RB0)*DTAN(LA2)
      N=1.D0
            CALL VINKEL(SI2,T,N)
      IF(SI2.LT.SI1) THEN
            SI2=SI2+PI*2.D0
      ENDIF
      
      IF(B2.LT.0) THEN
           SI2=SI2-PI
      ENDIF
      
      IF(B1.LT.0) THEN
           SI2=SI2-PI
      ENDIF
      
      F=F+1
                                 
      DSI=SI2-SI1      
      DLA=LA2-LA1
      IF(F.GT.4) THEN
           GOTO 12
      ENDIF
C******BEREGNING AV DLA,FRA 1 TIL 2

      B0=VINKELTR(A,B,RB0)
      W0=SQRT(1-E*(DSIN(B0))**2.D0)
      K1=(1-W0)/(1+W0)
      N1=(A-B)/(A+B)      
      R=E*DCOS(RB0)/2.D0
      R1=(1.D0+N1-K1/2.D0-(K1**2.D0)/4.D0)
      R2=K1/4.D0
      R3=(K1**2.D0)/16.D0
      DLA=DL+R*(R1*DSI-R2*(DSIN(2.D0*SI2)-DSIN(2.D0*SI1)))
      DLA=DLA+R*R3*(DSIN(4*SI2)-DSIN(4*SI1))  
            GOTO 11  


C******BEREGNING AV S1,S2 OG DS

         
   12 C=B*(1.D0+K1**2.D0/4.D0)/(1.D0-K1)
      D1=(K1/2.D0-(3.D0*K1**3.D0)/16.D0)
      D2=(K1**2.D0)/16.D0
      D3=(K1**3.D0)/48.D0
      
      CALL AVSTAND(C,D1,D2,D3,S1,SI1)
      CALL AVSTAND(C,D1,D2,D3,S2,SI2)
      DS=S2-S1
              
      
      T=1
      N=-DSIN(SI1)*DTAN(RB0)
             CALL VINKEL(A1,T,N)
      T=1
      N=-DSIN(SI2)*DTAN(RB0)
             CALL VINKEL(A2,T,N)
      
      IF(DS.LT.0.D0)THEN
             DS=-DS
      ENDIF
      
      RB0=RB0/RO
      LA1=LA1/RO
      LA2=LA2/RO
      SI1=SI1/RO
      SI2=SI2/RO
      A1=A1/RO
      A2=A2/RO
      WRITE(*,*) RB0
      WRITE(*,*) LA1,LA2
      WRITE(*,*) SI1,SI2
      WRITE(*,*) A1,A2  
      A1=A1*RO
      A2=A2*RO
 
      IF(DL.GT.0.D0) THEN
            A2=A2+PI
      ELSEIF(DL.LT.0.D0) THEN
            A1=A1+PI
      ENDIF
      
         
      A1=A1/RO
      A2=A2/RO

      WRITE(*,*)'ASIMUT FRA 1 TIL 2 OG FRA 2 TIL 1,I FORHOLD'
      WRITE(*,*)'TIL NORD OGS� P� DEN SYDLIGE HALVELLIPSOIDE'
      WRITE(*,110)'A1:',A1

      WRITE(*,110)'A2:',A2
      WRITE(*,*)'GEODETISK LINJE FRA 1 TIL 2:'
      WRITE(*,120)'S:',DS

  
  100 FORMAT(1X,A,$)
  110 FORMAT(1X,A,F14.9)
  120 FORMAT(1X,A,F12.3)    
  
   40 STOP

      END
               SUBROUTINE AKSER (P,A,B)
      
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
      
      
                SUBROUTINE VINKEL(A1,T,N)
      
      REAL*8 A1,T,N,PI
             PI=DATAN(1.D0)*4.D0
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

      REAL*8 FUNCTION VINKELTR(A,B,V)
      REAL*8 A,B,V

             VINKELTR=DATAN(A*TAN(V)/B)
      
      RETURN
      
      END