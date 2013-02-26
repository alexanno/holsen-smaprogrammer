      PROGRAM KRRAD
C*****************************************************************
C     BEREGNING AV KRUMNINGSRADIER.M,N,MR,RA  
C     MASKINEN SPùR ETTER BR=BREDDE OG AS=ASIMUT  
C*****************************************************************
      IMPLICIT LOGICAL (A-Z)

      INTEGER P
      
      REAL*8  RO, A, B, AR,BR, E, M, N, MR, RA, PI

      PI = DATAN(1.D0)*4.D0

      RO=180.D0/PI

      WRITE(*,*)  
      WRITE(*,*)' ************************************************'
      WRITE(*,*)' **       BEREGNING AV KRUMNINGSRADIER         **'
      WRITE(*,*)' **                                            **'
      WRITE(*,*)' **                        J.H. April-90       **'
      WRITE(*,*)' ************************************************'
      WRITE(*,*)  

      WRITE(*,*) 'VALG AV ELLIPSOIDE.P SETTES LIK 1 FOR :'
      WRITE(*,*) 'NORSK BESSELS ,LIK 2 FOR INTERNASJONAL:'
      WRITE(*,*) 'OG LIK 3 FOR FOR WGS-84-ELLIPSOIDEN:'
      WRITE(*,*) 'LEGG INN P:'
      READ(*,*) P
          
      WRITE(*,100) ' Legg inn bredde for punkt :'
      READ(*,*) BR
      BR=BR/RO

      WRITE(*,100) ' Legg inn asimut for punkt :'
      READ(*,*) AR

C---- Beregning starter her
      
      CALL AKSER(P,A,B)
            
      AR=AR/RO
      E=(A**2.D0-B**2.D0)/A**2.D0
      M=(DSIN(BR))**2
      M =DSQRT(1.D0-E*M)
                             
      N =A/M
      M =(1.D0-E)*A/M**3.D0
                                  
             		MR = DSQRT(M*N)		 
      RA=N*M/( N*DCOS(AR)**2 + M*DSIN(AR)**2 )

      WRITE(*,*)  
      WRITE(*,110)  ' Meridiankrumningsradius M  = ', M
      WRITE(*,110)  ' Normalkrumningsradius   N  = ', N
      WRITE(*,110)  ' Midlere krumningsradius MR = ', MR
      WRITE(*,110)  ' Krumningsradius i retning asimut AR = ', RA

  100 FORMAT ( 1X,A,$ )
  110 FORMAT ( 1X,A,F12.3 )
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