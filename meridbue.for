      PROGRAM MERIDBUE
C**************************************************************
C         BEREGNING AV MERIDIANBUEN MELLOM BREDDE B1 OG B2
C         ELLER B2 NèR B1 OG MERIDIANBUEN G ER GITT
C
C         MASKINEN SPùR ETTER B1 OG B2	
C         ELLER ETTER B1 OG G
C
C***************************************************************
      IMPLICIT LOGICAL (A-Z)
      
      INTEGER P,R
      
      REAL*8 RO,A,B,B1,B2,N,A1,A2,A3,A4,K,K1,K2,K3
      REAL*8 G,G1,PI
      PI=DATAN(1.D0)*4.D0
 
      RO=180.D0/PI
      WRITE(*,*)
      WRITE(*,*) '*********************************************'
      WRITE(*,*) '**                                         **'
      WRITE(*,*) '**                                         **'
      WRITE(*,*) '**          J. H. APRIL 1990               **'
      WRITE(*,*) '*********************************************'
      WRITE(*,*)
      WRITE(*,*) 'R SETTES LIK 0 NèR (B1,B2) ER GITT .I MOTSATT FALL'
      WRITE(*,*) 'SETTES R LIK 1 NèR B1 0G G ER GITT'
      READ(*,*) R
      
      WRITE(*,*)  'VALG AV ELLIPSOIDE.P SETTES LIK 1 FOR:'             
      WRITE(*,*)  'NORSK BESSELS,LIK 2 FOR INTERNASJONAL:'
      WRITE(*,*)  'OG LIK 3 FOR WGS-84-ELLIPSOIDEN:'         
      WRITE(*,*)  'LEGG INN P:'
      READ(*,*) P
            
           CALL AKSER(P,A,B)
      
      IF(R.EQ.1) THEN 
           GOTO 10
      ENDIF

      WRITE(*,*) 'LEGG INN BREDDENE B1 OG B2:'
      READ(*,*) B1,B2
      B1=B1/RO
      B2=B2/RO
      
C-----BEREGNING AV KONSTANTER      
      
      N=(A-B)/(A+B)
      K=A/(N+1.D0)

      A1=N
      A2=N**2
      A3=N**3
      A4=N**4
      K1=(1.D0+A2/4.D0+A4/64.D0)
      K2=( A1-A3*1.D0/8.D0)*3.D0
      K3=(A2-A4/4.D0)*15.D0/8.D0
      
      
C-----Meridianbue mellom B1 og B2
      
          CALL BUE(B2,B1,K,K1,K2,K3,G,G1,R,N)

      WRITE(*,110) 'MERIDIANBUE = ',G1
  
      STOP
      

C-----HER STARTER BEREGNINGEN AV B2 NèR B1 OG G ER GITT
      
   10 WRITE(*,*) 'LEGG INN B1 OG G:'
      READ(*,*) B1,G
      B1=B1/RO
      N=(A-B)/(A+B)
      K=A/(N+1.D0)
      K1=1.D0+N*N/4.D0+N**4.D0/64.D0
      K2=(N-N*N*N/8.D0)*3.D0
      K3=(N*N-(N**4.D0)*1/4.D0)*15.D0/8.D0
  
           CALL BUE(B2,B1,K,K1,K2,K3,G,G1,R,N)
      
      B2=B2*RO
     
      WRITE(*,120) 'BREDDE B2=',B2
      WRITE(*,110) 'MERIDIANBUE =',G

  
  100 FORMAT (1X,A,$)
  120 FORMAT (1X,A,F12.9)
  110 FORMAT (1X,A,F12.3)
      STOP
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
      
      SUBROUTINE BUE(BR,BO,K,K1,K2,K3,G,G1,R,N)
      
      INTEGER R
            
      REAL*8 BR,BO,K,K1,K2,K3,G1,DB,BM,G,N
      G1=0
      DB=0
      BM=0
      IF(R.EQ.0)THEN
        DB=BR-BO
        BM=BR-DB/2.D0
           
      ENDIF
        
          
   30 IF(R.EQ.1) THEN
      DB= DB+(G-G1)/(K*K1)
      BM=BO+DB/2.D0
      ENDIF
           
      G1 = K*(K1*DB-K2*DCOS(2.D0*BM)*DSIN(DB)+K3*DCOS(4.D0*BM)
     &     *DSIN(2.D0*DB))
      
      G1 = G1- K*(N**3*DCOS(6.D0*BM)*DSIN(3.D0*DB))*35.D0/24.D0
      
      G1 = G1+ K*(N**4.D0*DCOS(8*BM)*DSIN(4*DB)*315.D0)/256.D0
      
      IF(R.EQ.0) THEN
        G=G1
        RETURN
      ENDIF
      
      IF(R.EQ.1.AND.ABS(G-G1).GT.1D-4) THEN
        GOTO 30
      ELSE
        BR=BO+DB
               
         G=G1    

      RETURN
      ENDIF

      END