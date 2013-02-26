
C*****************************************************************  

C     BEREGNING AV TOPOSENTRISKE KOORDINIATER AV GEOSENTRISKE OG

C     VICE VERSA

C*****************************************************************

      IMPLICIT LOGICAL (A-Z)
      REAL*8 B,L,X1,X2,X3,Y1,Y2,Y3,A,SH,S,ZD,PI,RO
      INTEGER T
      PI=DATAN(1.D0)*4.D0
      RO=180/PI
      WRITE(*,*) '**************************************************'
      WRITE(*,*) '**   BEREGNING AV TOPOSENTRISKE KOORDINATER AV  **' 
      WRITE(*,*) '**     GEOSENTRISKE KOORDINATER OG OMVENDT      **' 
      WRITE(*,*) '**                                              **'
      WRITE(*,*) '**                       J.H. APRIL 90          **'
      WRITE(*,*) '**************************************************'
      WRITE(*,100) 'LEGG INN BREDDE:'
      READ(*,*) B
      B=B/RO
      WRITE(*,100) 'LEGG INN LENGDE:'
      READ(*,*) L
      L=L/RO
      WRITE(*,100) 'LEGG INN  X,Y,Z:'
         
      READ(*,*) X1,X2,X3
             
               
      WRITE(*,100) 'T SETTES LIK 0 VED TRANSFORMASJON FRA GEOS',
     &   'TIL TOPOS. FRA TOPOS TIL GEOS SETTES T=1'
      WRITE(*,*)
      
      WRITE(*,100) 'LEGG INN T:'
      READ(*,*) T
      IF(T.EQ.1) THEN
      GOTO 10
      ENDIF

C- - - - BEREGNING AV TOPOSENTRISKE KOORDINATER
         
      Y1=-DSIN(B)*DCOS(L)*X1-DSIN(B)*DSIN(L)*X2+DCOS(B)*X3
      Y2=-DSIN(L)*X1+DCOS(L)*X2
      
      Y3=DCOS(B)*DCOS(L)*X1 +DCOS(B)*DSIN(L)*X2+DSIN(B)*X3
              
      SH=DSQRT(Y1**2+Y2**2)
      
      A=DASIN(Y2/SH)
     
      
      IF(Y1.LT.0.AND.Y2.GT.0) THEN
         A=PI-A
      ELSEIF(Y1.LT.0.AND.Y2.LT.0) THEN
         A=PI-A
      ELSEIF(Y1.GT.0.AND.Y2.LT.0) THEN
         A=2.D0*PI+A
      ENDIF
      
      A=A*RO
           
      S=DSQRT(Y1**2+Y2**2+Y3**2)
      ZD=DASIN(Y3/S)
      ZD=PI/2-ZD
      
      ZD=ZD*RO
      
      IF(T.EQ.0) THEN
      WRITE(*,100) 'ASIMUT,SENITDISTANSE,TOPOSENTRISKE KOORDIN.'
      ENDIF
      WRITE(*,*)                
      WRITE(*,*)    
      WRITE(*,100) 'ASIMUT (360-GRADERS DELING)'
      WRITE(*,*)
      WRITE(*,110) 'A:',A
      WRITE(*,*)
      
      WRITE(*,110) 'SENITDISTANSE (360-GRADER ):',ZD
      WRITE(*,*)
      WRITE(*,120) 'HORISONTALAVSTAND:',SH
   30 IF(T.EQ.1) THEN 
      WRITE(*,100) 'GEOSENTISKE KOORDINATER'
      ENDIF
      WRITE(*,*)         
            
      WRITE(*,120) 'X:',Y1
      WRITE(*,120) 'Y:',Y2
      WRITE(*,120) 'z:',Y3
      T=0
           
      IF(T.EQ.0) THEN
      GOTO 20
      ENDIF

C- - - -      HER REGNES GEOSENTRISKE KOORDINATER
      
      
   10 Y1=-DSIN(B)*DCOS(L)*X1-DSIN(L)*X2+DCOS(B)*DCOS(L)*X3
      Y2=-DSIN(B)*DSIN(L)*X1+DCOS(L)*X2+DCOS(B)*DSIN(L)*X3
      Y3=DCOS(B)*X1+DSIN(B)*X3
      
      IF(T.EQ.1)THEN
      GOTO 30
      ENDIF
  100 FORMAT(1X,A,$)
  110 FORMAT(1X,A,F12.7)
  120 FORMAT(1X,A,F12.3)
         
   20 END

      
            
      