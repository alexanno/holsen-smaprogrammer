      PROGRAM BLXY
C***************************************************************************
C     BEREGING AV PLANE KOORDINATER (X,Y) AV GEOGRAFISKE(B,l)
C     OG OMVENDT.
C     (X,Y) ER GAUSS-KRUGERSKE ELLER UTM-KOORDINATER.
C***************************************************************************

      IMPLICIT LOGICAL (A-Z)

      REAL*8 A,B,BR,L,L1,BO,N,N1,T,A1,A2,A3,A4,A5
    
      REAL*8 X,Y,PI,RO,B1,B2,B3,B4,B5,A6,ET,ETF,NF,BF,G

      REAL*8 TF,K,K1,K2,K3,G1

      
      INTEGER P,R,S
   
      
      PI=DATAN(1.D0)*4.D0
      RO=180.D0/PI

      WRITE(*,*)
      WRITE(*,*)   '******************************************'
      WRITE(*,*)   '**BEREGNING AV (X,Y) NèR (B,L) ER KJENT **'             
      WRITE(*,*)   '** OG OMVENDT                           **'          
      WRITE(*,*)   '**                                      **'           
      WRITE(*,*)   '**          JON HOLSEN, APRIL 90        **'           
      WRITE(*,*)   '**                IKO, NTH              **'           
      WRITE(*,*)   '******************************************'
      
      WRITE(*,*) 'S SETTES LIK 0 FOR NGO-SYSTEM OG 1 FOR UTM'
      WRITE(*,100) 'LEGG INN S:'  
      READ(*,*) S
      
      WRITE(*,*)  'R SETTES LIK 0 NèR (B,L) ER GITT I MOTSATT FALL'
      WRITE(*,*)  'SETTES R LIK 1'                                
      WRITE(*,100) 'LEGG INN R:'
      READ(*,*) R   
      
      WRITE(*,*) 'VALG AV ELLIPSOIDE.P SETTES LIK 1 FOR'
      WRITE(*,*) 'NORSK BESSELS,LIK 2 FOR INTERNASJONAL'
      WRITE(*,*) 'OG LIK 3 FOR WGS-84-ELLIPSOIDEN'
      WRITE(*,100) 'LEGG INN P:'
      READ(*,*)   P
      
      IF(P.EQ.1)THEN
        A=6377492.018D0
        B=6356173.509D0
      ELSEIF(P.EQ.2)THEN
        A=6378388.000D0
        B=6356911.946D0
      ELSEIF(P.EQ.3)THEN
        A=6378137.000D0
        B=6356752.314D0
      ENDIF
          
           
      N=(A-B)/(A+B)
      K=A/(N+1.D0)
      K1=1.D0+N*N/4.D0+N**4.D0/64.D0
      K2=(N-N*N*N/8.D0)*3.D0
      K3=(N*N-(N**4.D0)*1/4.D0)*15.D0/8.D0
           
         
      WRITE(*,100) 'LEGG INN  BREDDEN FOR X-AKSENS NULLPUNKT:'
      READ(*,*) BO
      WRITE(*,100) 'LEGG INN GEOGRAFISK LENGDE FOR X-AKSEN:'
      READ(*,*) L1
      L1=L1/RO

      BO=BO/RO
                  
      IF(R.EQ.1) GOTO 10
      
      WRITE(*,100) 'LEGG INN GEOGRAFISK BREDDE:'
      READ(*,*) BR

      WRITE(*,100) 'LEGG INN GEOGRAFISK LENGDE:'
      READ (*,*) L
      L=L/RO

      BR=BR/RO
      
      L=L-L1
                                         
      
     
      WRITE(*,*) 'HER STARTER BEREGNINGEN MED B OG L SOM KJENT'

      ET=(A**2-B**2)/B**2
 
      ET=ET*DCOS(BR)**2
      N1=(A**2)/(DSQRT(1.D0+ET)*B)
      T =DTAN(BR)
      A1=N1*DCOS(BR)
      A2=-(N1*T*DCOS(BR)**2.D0)/2.D0
      A3=-(N1*DCOS(BR)**3.D0)*(1.D0-T**2.D0+ET)/6.D0
      A4=N1*T*DCOS(BR)**4.D0*(5.D0-T**2.D0+9*ET+4.D0*ET**2.D0)
     &/24.D0
      A6=N1*T*DCOS(BR)**6*(61.D0-58.D0*T**2+T**4+270*ET
     &-330*ET*T**2.D0)/720
                     
      A5=N1*DCOS(BR)**5.D0*(5.D0-18.D0*T**2.D0+T**4.D0
     &+14.D0*ET-58.D0*ET*T**2)/120  
    
      X=-A2*L**2.D0+A4*L**4.D0+A6*L**6.D0
      Y= A1*L-A3*L**3+A5*L**5
      
            
        CALL MERIDBUE (BR,BO,K,K1,K2,K3,G,G1,R,N)
      
     
      X=X+G
            
      IF(S.EQ.1) THEN
        X=X*0.9996
        Y=Y*0.9996+500000.D0
      ENDIF
                  
      WRITE(*,110) 'X:',X
      WRITE(*,*)
      WRITE(*,110) 'Y:',Y
      STOP  


   10 WRITE(*,*) 'HER STARTER BEREGNINGEN MED KJENT X OG Y'

      WRITE(*,100) 'LEGG INN X:'
      READ(*,*)   X
      WRITE(*,100) 'LEGG INN Y:'
      READ(*,*)    Y

      IF (S.EQ.1) THEN
      X=X/0.9996
      Y=(Y-500000.D0)/0.9996
      ENDIF

      G=X
      
           CALL MERIDBUE(BF,BO,K,K1,K2,K3,G,G1,R,N)
 
      
      ETF=(A-B)*(A+B)/B**2.D0
      ETF=ETF*COS(BF)**2.D0
      
      NF =A**2.D0/(DSQRT(1.D0 +ETF)*B)
      TF =DTAN(BF)
      B1=1.D0/(NF*DCOS(BF))
      B2=TF*(1.D0+ETF)/(2.D0*NF**2.D0) 
      B3=(1.D0+2.D0*TF**2.D0+ETF)/(6.D0*NF**3.D0*DCOS(BF))
      B4=TF*(5.D0+3.D0*TF**2.D0+6*ETF-6.D0*ETF*TF**2.D0)
      B4=B4/(24.D0*NF**4.D0)   
            
      B5=(5.D0+28.D0*TF**2.D0+24.D0*TF**4.D0)
      B5=B5/(120.D0*NF**5.D0*DCOS(BF))
                          
      BR=BF+(-B2*Y**2.D0+B4*Y**4.D0)
      L =(B1*Y-B3*Y**3.D0+B5*Y**5.D0)
      L=L+L1
      L=L*RO
      BR=BR*RO
      WRITE(*,100) 'GEOGRAFISKE KOORDINATER'
      WRITE(*,*)
      WRITE(*,120) 'B:',BR
      WRITE(*,*)
      WRITE(*,120) 'L:',L

  100 FORMAT(1X,A,$)
  120 FORMAT(1X,A,F13.9)
  110 FORMAT(1X,A,F13.4)

      STOP
      
      END



      SUBROUTINE MERIDBUE(BR,BO,K,K1,K2,K3,G,G1,R,N)
      
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
      