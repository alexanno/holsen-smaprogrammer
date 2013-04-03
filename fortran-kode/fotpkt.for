                           PROGRAM FOTPKT

      INTEGER N3,I,T,S2,LUINN,NUM(1:50)
      REAL*8 KO(1:20),F1,F2,PI,RO
      REAL VK(1:10,1:10)
      CHARACTER*10 FILNAVN
      LUINN=10
      DO 1 I=1,10
      DO 2 T=1,10
      VK(I,T)=0
    2 CONTINUE
    1 CONTINUE
      PI=DATAN(1.D0)*4.D0
      RO=200.D0/PI
      DO 3 I=1,20
      KO(I)=0
    3 CONTINUE      
      WRITE(*,*)'LEGG INN F| OG F2'
      READ(*,*)F1,F2
      F1=F1/RO
      F2=F2/RO
      
      N3=6
      S2=6
      WRITE(*,*)'FILNAVN,GITTE KOORDINATER'
      READ(*,103) FILNAVN
  103 FORMAT(A) 
        
      
      OPEN(UNIT=LUINN,FILE=FILNAVN,FORM='FORMATTED',STATUS='OLD')

      DO 5 I=1,6

      T=2*I
      READ(LUINN,140) NUM(I),KO(T-1),KO(T)
  140 FORMAT(T5,I2,T7,F11.4,T18,F11.4)
      WRITE(*,310)NUM(I),KO(I-1),KO(I)
  310 FORMAT(T5,I2,T8,F11.4,T20,F11.4/)
        
    5 CONTINUE
      CLOSE (LUINN)
        
      LUINN=12      
      WRITE(*,*)'LEGG INN FILNAVN VEKTSKO'
      READ(*,103) FILNAVN

      OPEN(UNIT=LUINN,FILE=FILNAVN,FORM='FORMATTED',STATUS='OLD')
      T=1
      READ(LUINN,102)(VK(T,I),I=1,6)
      
  
      
  102 FORMAT(T4,F6.2,T11,F6.2,T18,F6.2,T25,F6.2,T32,F6.2,T37,F8.2)
            
      WRITE(*,104) (VK(T,I),I=1,6)
  104 FORMAT(T6,6(F6.2,TR2))

      CLOSE (LUINN)
      STOP
      END