                   

                    PROGRAM KONVERG
C******************************************************************
C     WRITE(*,*) PROGRAMMET BEREGNER PLAN MERIDIANKONVERGENS I GAUSS 
C     WRITE(*,*) KONFORME PROJEKSJON OG UTM.(B,L) ELLER (X,Y) Mè VíRE 
C     WRITE(*,*) GITT I DET AKTUELLE PUNKT.TRE FORSKJELLIGE 
C     WRITE(*,*) ELLIPSOIDER KAN VELGES UNDER BRUK AV PROGRAMMET.                                        VED BRUK AV PROGRAMMET.
C                         
C                         FEBRUAR -91
C                            J.H.
C******************************************************************
       
      IMPLICIT LOGICAL (A-Z)

      INTEGER I,P,S
      REAL*8 A,B,B1,BO,L1,LO,PI,RO,DL,E,T,C,SI,CO,X,Y,N,NF,G,G1
      REAL*8 K,K1,K2,K3,BF

      PI=DATAN(1)*4.D0
      RO=PI/180.D0

      WRITE(*,*) 'I SETTES LIK 1 NèR (B,L) ER GITT OG 2 NèR (X,Y) '
      WRITE(*,*) 'ER GITT'
      READ(*,*) I
      WRITE(*,*)  'VALG AV ELLIPSOIDE.P SETTES LIK 1 FOR:'             
      WRITE(*,*)  'NORSK BESSELS,LIK 2 FOR INTERNASJONAL:'
      WRITE(*,*)  'OG LIK 3 FOR WGS-84-ELLIPSOIDEN:'         
      WRITE(*,*)  'LEGG INN P:'
      READ(*,*) P                  
                        CALL AKSER (P,A,B)

      WRITE(*,*)'LEGG INN LENGDEGRADEN FOR FOR X-AKSEN'
      READ(*,*) LO
      
      IF(I.EQ.2) THEN
          GOTO 12
      ENDIF
      WRITE(*,*)'LEGG INN BREDDE OG LENGDE'
      READ(*,*) B1,L1
      DL=L1-LO
      B1=B1*RO
      DL=DL*RO

      E=((DCOS(B1))**2.D0)*(A**2.D0-B**2.D0)/B**2.D0
      T=DTAN(B1)
      SI=DSIN(B1)
      CO=DCOS(B1) 
      
        
      C =DL*SI+(DL**3.D0)*(CO**2.D0)*(1+3*E+2*E**2)*SI/3
          
      C=C+(DL**5)*(CO**4)*(2-T**2)*SI/15
      
      C=C*1.11111111111/RO
      
      IF(I.EQ.1) THEN
           GOTO 50
      ENDIF


C**********HER BEREGNES C NèR X OG Y ER GITT
      
   12 WRITE(*,*)'LEGG INN BREDDE FOR X-AKSENS NULLPUNKT'
      READ(*,*) BO
      BO=BO*RO
      WRITE(*,*)'LEGG INN X,Y'
      READ(*,*) X,Y
      WRITE(*,*)'S SETTES LIK 1 FOR NGO-SYSTEM OG LIK 2 FOR UTM'
      READ(*,*) S

      IF(S.EQ.2) THEN
          X=X/.9996D0
          Y=(Y-500000.D0)/.9996D0
      ENDIF


C***********HER STARTER BEREGNINGEN AV BF OG NF
      
      N=(A-B)/(A+B)
      K=A/(N+1.D0)
      K1=1.D0+N*N/4.D0+N**4.D0/64.D0
      K2=(N-N*N*N/8.D0)*3.D0
      K3=(N*N-(N**4.D0)*1/4.D0)*15.D0/8.D0
      G=X

               CALL FOTBR(BF,BO,K,K1,K2,K3,G,G1,N)

      E=(A-B)*(A+B)/(B**2.D0)
      E=E*((DCOS(BF))**2.D0)
      NF=A**2.D0/(B*SQRT(1+E))
      T=DTAN(BF)
      
      C=Y*T/NF-(1+T**2-E-2*E**2)*T*Y**3.D0/(3*NF**3)
      C=C+(2+5*T**2+3*T**4)*T*Y**5/(15*NF**5)
      C=C/(RO*.9D0)
      
   50 WRITE(*,*)
      WRITE(*,100)'C ER MERIDIANKONVERGENS I GON'
      WRITE(*,*)
      WRITE(*,120)'C:',C
      
  100 FORMAT(1X,A,$)
  120 FORMAT(1X,A,F10.7)

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

                  SUBROUTINE FOTBR(BR,BO,K,K1,K2,K3,G,G1,N)
      
      REAL*8 BR,BO,K,K1,K2,K3,G,G1,N,DB,BM

      G1=0
      DB=0
      BM=0
      

     
   30 DB= DB+(G-G1)/(K*K1)
      BM=BO+DB/2.D0
      
      
      G1 = K*(K1*DB-K2*DCOS(2.D0*BM)*DSIN(DB)+K3*DCOS(4.D0*BM)
     &     *DSIN(2.D0*DB))
      
      G1 = G1- K*(N**3*DCOS(6.D0*BM)*DSIN(3.D0*DB))*35.D0/24.D0
      
      G1 = G1+ K*(N**4.D0*DCOS(8*BM)*DSIN(4*DB)*315.D0)/256.D0
         
          
      IF(ABS(G-G1).GT.1D-4) THEN
        GOTO 30
      ELSE
        BR=BO+DB
         G=G1    

      RETURN
      ENDIF

      END