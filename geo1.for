                       
                 PROGRAM GEO1
     
C*******************************************************************
C     WRITE(*,*) BEREGNING AV GEOGRAFISKE KOORDINATER FOR PUNKT 2
C     WRITE(*,*) OG ASIMUT FOR LINJEN I PUNKT 2.KOORDINATENE
C     WRITE(*,*) (B1,L1),GEODETISK LINJE S12 OG ASIMUT A12 GITT.
C     WRITE(*,*)       J.H. JUNI 90
C*******************************************************************
        
      IMPLICIT LOGICAL(A-Z)
      INTEGER P,R

      REAL*8 A,B,B1,L1,B2,L2,B3,BE1,BE2,BE3,M1,M2,M3,N1,N2,N3
      REAL*8 R1,R2,R3,X1,X2,X3,S,S1,S2,PI,RO,A1,A2,A3,GA
      REAL*8 K1,K2,K3,G1,G2,T1,T2,T3,C1,C2,C3,C4,W1,W2,D1,D2

      PI=DATAN(1)*4.D0
      RO=PI/180.D0
      R =0
      
      WRITE(*,*)'LEGG INN B1,L1,S12,A12:'                      
      READ(*,*) B1,L1,S,A1 
      B1=B1*RO
      L1=L1*RO
      A1=A1*RO

      WRITE(*,*) 'VALG AV ELLIPSOIDE.P SETTES LIK 1 FOR :'
      WRITE(*,*) 'NORSK BESSELS ,LIK 2 FOR INTERNASJONAL:'
      WRITE(*,*) 'OG LIK 3 FOR FOR WGS-84-ELLIPSOIDEN:'
      WRITE(*,*) 'LEGG INN P:'
      READ(*,*) P
      
                 CALL AKSER(P,A,B)

      
C********** BEREGNING AV FORELùPIGE VERDIER FOR B2,L2 OG A2
      
      BE1=DATAN(B*DTAN(B1)/A)
      M1= DCOS(BE1)/DCOS(B1)
      N1= A*M1
      M1=(M1**3.D0)*(B**2.D0)/A
      R1=M1*N1/(N1*(DCOS(A1))**2.D0+M1*(DSIN(A1))**2.D0)
      B2=B1+S*DCOS(A1)/M1
      L2=L1+S*DSIN(A1)/(N1*DCOS(B1))
              
   
   10 BE2=DATAN(B*DTAN(B2)/A)
      B3=(B1+B2)/2
      BE3=DATAN(B*DTAN(B3)/A)

C********** BEREGNING AV GEOSENTRISK DELTA X,DELTA Y,DELTA Z

      X1=A*(DCOS(BE2)*DCOS(L2)-DCOS(BE1)*DCOS(L1))
      X2=A*(DCOS(BE2)*DSIN(L2)-DCOS(BE1)*DSIN(L1))
      X3=B*(DSIN(BE2)-DSIN(BE1))

      S1=SQRT(X1**2.D0+X2**2.D0+X3**2.D0)
      

C******** KORDEN TIL GITT GEODETISK LINJE,S2

      CALL ASIMUT(-X1,-X2,-X3,B2,L2,A2)
      
           CALL KRRAD(A,B,B1,BE1,A1,M1,N1,R1)
           CALL KRRAD(A,B,B2,BE2,A2,M2,N2,R2)
           CALL KRRAD(A,B,B3,BE3,A3,M3,N3,R3)

           R3=(R1+R2+4*R3)/6
      GA=S/(2.D0*R3)
      S2=2.D0*R3*DSIN(GA)
      

C******** DIFFERENSIALKOEFFISIENTER
      
      K1=-A*DSIN(BE2)*DCOS(L2)
      K2=-A*DSIN(BE2)*DSIN(L2)
      K3= B*DCOS(BE2)
      G1=-A*DCOS(BE2)*DSIN(L2)
      G2= A*DCOS(BE2)*DCOS(L2)
      T1=DSIN(B1)*DCOS(L1)*DSIN(A1)-DSIN(L1)*DCOS(A1)
      T2=DSIN(B1)*DSIN(L1)*DSIN(A1)+DCOS(L1)*DCOS(A1)
      T3=-DCOS(B1)*DSIN(A1)

      C1=K1*T1+K2*T2+K3*T3
      C2=G1*T1+G2*T2
      C3=(K1*X1+K2*X2+K3*X3)/S2
      C4=(G1*X1+G2*X2)/S2


C******** KONSTANTLEDD

      W1=X1*T1+X2*T2+X3*T3
      W2=S1-S2

C******** LùSNING AV DIFFERENSIALLIKNINGENE

      D1=0
      D2=0
      
      D1=(-C4*W1+C2*W2)/(C1*C4-C2*C3)
      D2=(C3*W1-C1*W2)/(C1*C4-C2*C3)

      R=R+1

      BE2=BE2+D1
      L2 =L2 +D2
      B2=DATAN(A*DTAN(BE2)/B)
      IF(R.LT.6) GOTO 10

      B2=B2/RO
      L2=L2/RO
      A2=A2/RO

 
      WRITE(*,*) 'GEOGRAFISKE KOORDINATER:'
      WRITE(*,110) 'B2:',B2
      WRITE(*,110) 'L2:',L2
      WRITE(*,*)
      WRITE(*,*) 'ASIMUT FRA PUNKT 2 TIL PUNKT 1:'
      WRITE(*,110) 'A2:',A2
      
      
  100 FORMAT(1X,A,$)
  110 FORMAT(1X,A,F14.9)

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

                 SUBROUTINE ASIMUT(X1,X2,X3,B1,L1,A1)
       
           
      REAL*8 X1,X2,X3,B1,L1,A1,T,N,PI

      T=-X1*DSIN(L1)+X2*DCOS(L1)
      N=-X1*DSIN(B1)*DCOS(L1)-X2*DSIN(B1)*DSIN(L1)+X3*DCOS(B1)
              
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


                    SUBROUTINE KRRAD(A,B,B1,BE1,A1,M1,N1,R1)
      
      REAL*8 A,B,B1,BE1,A1,M1,N1,R1,PI

      PI=DATAN(1)*4.D0
         
      M1=DCOS(BE1)/DCOS(B1)
      N1=A*M1
      M1=M1**3.D0*B**2.D0/A
      
      R1=M1*N1/(N1*(DCOS(A1))**2.D0+M1*(DSIN(A1))**2.D0)
      
      IF(A1.EQ.0.OR.A1.EQ.2.D0*PI) THEN
            R1=M1
      ENDIF      
                           
         RETURN
      END            