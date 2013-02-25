PROGRAM Excent(input,output);

{$N+,E+}   { Links with the 8087 numeric coprocessor or if it is not   }
           { found emulate by the run-time-library. (Needed to perform }
           { all real-type calculations.)                              }

{*********************************************************}
{  Program f”r kontroll av excentricitetsber„kningar vid  }
{  trianguleringsm„tningar under felt”vningarna i Land-   }
{  m†lingskurset.                                         }
{                                                         }
{  F”rfattare :    Bobo Nordahl, IGF                      }
{                                                         }
{  Datum :         20/7 1990                              }
{                                                         }
{  F”r„ndringar :                                         }
{                                                         }
{*********************************************************}

USES   CRT;  { Kan ej anv„ndas p† alla pc, (en speciell processor beh”vs).}
             { Anv„nds f”r att f† b„ttre sk„rmhantering.                  }


CONST        maxant = 10;
                rho = 63.66197724;


TYPE          {*  ord = Array [1..15] of CHAR;  * Fel i turbo pascal *}
                ord = string[15];
            punkter = Array [1..maxant] of ord;
      Real_1_maxant = Array [1..maxant] of REAL;
       Ext_1_maxant = Array [1..maxant] of EXTENDED;

VAR      X,Y,avstand,retn,exc_vinkel   : Real_1_maxant;
         res,eps                       : Ext_1_maxant;
         totantal,ant_nypunkt          : INTEGER;
         alternativ                    : INTEGER;
         senter                        : ord;
         siktepunkt                    : punkter;
         X_sen,Y_sen,e                 : REAL;


{*********************************************************}
{*********************************************************}


PROCEDURE Start;

BEGIN

  WRITELN;
  WRITELN ('*****************************************************');
  WRITELN ('*****************************************************');
  WRITELN ('*                                                   *');
  WRITELN ('*      Program f”r excentricitetsber„kningar        *');
  WRITELN ('*                                                   *');
  WRITELN ('*****************************************************');
  WRITELN ('*****************************************************');
  WRITELN;
  WRITELN ('Programmet kan k”ras med tv† olika alternativ.');
  WRITELN;
  WRITE ('Alternativ 1: Excentrisk uppst„llning i Torshaug med sikt till');
  WRITELN (' Nidarvoll,');
  WRITELN ('              Havstein kirke samt en nypunkt vid NTH.');
  WRITELN ('Alternativ 2: Excentrisk uppst„llning i valfri punkt och sikt ');
  WRITELN ('              till max 10 valfria punkter.');
  WRITELN;
  WRITE ('V„lj alternativ ( 1/2 ) : ');
  READLN ( alternativ );

END; { Start }


{*********************************************************}


PROCEDURE Inles_alternativ1;

VAR       i                 : INTEGER;

BEGIN

  senter := 'TORSHAUG';
  totantal := 3;
  ant_nypunkt := 1;
  siktepunkt[1] := 'HAVSTEIN KIRKE';
  siktepunkt[2] := 'NIDARVOLL';
  siktepunkt[3] := 'NTH_(Nypunkt)';
  X_sen := -1711.82;
  Y_sen := -2935.94;
  X[1]  := -2453.57;
  Y[1]  := -1314.50;
  X[2]  := -3542.83;
  Y[2]  :=   141.18;

  WRITELN;
  WRITELN;
  WRITELN;
  WRITE ('Avst†nd ges i enheten meter och vinklar i gon med ');
  WRITE ('positiv riktning medsols.');
  WRITELN;
  WRITELN;
  WRITELN;
  WRITELN ('Ge ca-avst†nd fr†n ',senter,' till ',siktepunkt[3] );
  WRITE ('Avst†nd : ');
  READ ( avstand[3] );

  WRITELN;
  WRITELN ('Ge m„tta reducerade riktningar fr†n den excentriska uppst„llningen');
  FOR i:=1 TO totantal DO BEGIN
      WRITE ('till ',siktepunkt[i],' : ');
      READ ( retn[i] );
  END;

  WRITELN;
  WRITE ('Ge in excentricitetsvinklarna,');
  WRITELN ('(N.B! Riktningen 0 gon mot ',senter,')');
  FOR i:=1 TO totantal DO BEGIN
      WRITE ('till ',siktepunkt[i],' : ');
      READ ( exc_vinkel[i] );
  END;

  WRITELN;
  WRITELN ('Ge in avst†ndet fr†n instrumentet till ',senter);
  WRITE (' i meter : ');
  READLN ( e );

END; { Inles_alternativ1 }


{*********************************************************}


PROCEDURE Inles_alternativ2;

VAR      i,j                           : INTEGER;


BEGIN

  WRITELN;
  WRITELN;
  WRITELN;
  WRITE ('Avst†nd ges i enheten meter och vinklar i gon med ');
  WRITE ('positiv riktning medsols.');
  WRITELN;
  WRITELN;
  WRITELN;
  WRITELN;
  WRITE ('Ge punktnamn p† den centriska punkten : ');
  READ (senter);
  WRITELN;
  WRITELN ('Ange antal tillsiktningspunkter som m„ts fr†n den ');
  WRITE ('excentriska uppst„llningen.   Antal? : ');
  READ ( totantal );
  WRITELN;
  WRITE ('Hur m†nga av dessa „r nypunkter? : ');
  READLN ( ant_nypunkt );

  FOR i:=1 TO ant_nypunkt DO BEGIN
      j := totantal - ant_nypunkt +i ;
      WRITELN;
      WRITE ('Ge punktnamn p† nypunkt nr ',i,' :');
      READLN ( siktepunkt[j] );
  END;

  WRITELN;

  WRITELN ('Ge punktnamn p† ”vriga tillsiktningspunkter.');

  FOR i:=1 TO (totantal - ant_nypunkt) DO BEGIN
      j := i + ant_nypunkt;
      WRITELN;
      WRITE ('Tillsiktningspunkt nr ',j,' :');
      READLN ( siktepunkt[i] );
  END;

  WRITELN;
  WRITELN ('Ge koordinater f”r den centriska punkten ',senter);
  WRITE ('X= ');
  READ ( X_sen );
  WRITE ('Y= ');
  READ ( Y_sen );
  WRITELN;

  FOR i:=1 TO (totantal - ant_nypunkt) DO BEGIN
      WRITELN ('Ge koordinater f”r punkten ',siktepunkt[i] );
      WRITE ('X= ');
      READ ( X[i] );
      WRITE ('Y= ');
      READ ( Y[i] );
  END;

  WRITELN;

  FOR i:= (totantal - ant_nypunkt +1) TO totantal DO BEGIN
      WRITELN ('Ge ca-avst†nd fr†n ',senter,' till ',siktepunkt[i] );
      WRITE ('Avst†nd : ');
      READ ( avstand[i] );
  END;

  WRITELN;
  WRITELN ('Ge m„tta reducerade riktningar fr†n den excentriska uppst„llningen');
  FOR i:=1 TO totantal DO BEGIN
      WRITE ('till ',siktepunkt[i],' : ');
      READ ( retn[i] );
  END;

  WRITELN;
  WRITE ('Ge in excentricitetsvinklarna,');
  WRITELN ('(N.B! Riktningen 0 gon mot ',senter,')');
  FOR i:=1 TO totantal DO BEGIN
      WRITE ('till ',siktepunkt[i],' : ');
      READ ( exc_vinkel[i] );
  END;

  WRITELN;
  WRITELN ('Ge in avst†ndet fr†n instrumentet till ',senter);
  WRITE (' i meter : ');
  READLN ( e );

END; { Inles_alternativ2 }


{*********************************************************}


PROCEDURE Kontroll_inles;

  { Testar om excentervinklarna avviker "f”r" mycket fr†n }
  { de m„tta riktningsvinklarna. "F”r" mycket avg”rs vid  }
  { ett test med avviket i sidled mot centriska punkten   }
  { med konstanten maxavvik.                              }

LABEL     777;
CONST     maxavvik = 0.005;

VAR       i,j              : INTEGER;
          test1,test2      : Real_1_maxant;
          svar             : CHAR;

BEGIN

  777: { B”rja h„r n„r nya testv„rden skall ber„knas efter „ndring. }

  FOR i:=1 TO (totantal-1) DO BEGIN
      j := i+1;
      test1[i] := (retn[j]-retn[i]) - (exc_vinkel[j]-exc_vinkel[i]);
      test2[i] := e*Sin(Abs(test1[i]/rho));
  END;

  FOR i:=1 TO (totantal-1) DO BEGIN
      IF (test2[i] > maxavvik) OR (abs(test1[i]) > 1)
                              AND (abs(test1[i]) < 399) THEN BEGIN
         ClrScr; { Se kommentar vid ”ppnandet av "Crt Unit" ovan }
         WRITELN;
         WRITELN;
         WRITELN ('*************           N.B!        *************');
         WRITELN;
         WRITE ('Excentricitetsvinkeln mellan ',siktepunkt[i],' och ');
         WRITELN ( siktepunkt[i+1],', ');
         WRITELN ('avviker fr†n de m„tta riktningarna med ',test1[i]:6:4,'gon.');
         WRITELN;
         FOR j:=i TO i+1 DO BEGIN
             WRITELN;
             WRITE ('Den nuvarande excentricitetsvinkeln mot ');
             WRITELN ( siktepunkt[j],' = ',exc_vinkel[j]:6:4);
             WRITE ('Skall denna „ndras? (J/N) ');
             READLN ( svar );
             IF svar IN ['J','j','Y','y'] THEN BEGIN
                WRITE ('Ge ny excentricitetsvinkel f”r ',siktepunkt[j],' : ');
                READLN ( exc_vinkel[j] );
                GOTO 777;      { Ber„kna nya testv„rden efter en „ndring. }
             END;
         END;
      END;
  END;
END; { Kontroll_inles }


{*********************************************************}


PROCEDURE Beregn;

  { Ber„knar reducerade vinkler efter korrektion f”r  }
  { excentricitet p† en station.                      }

CONST     pi  =  3.141592654;

VAR       i,j               : INTEGER;
          delta_X2,delta_Y2 : REAL;
          sin_eps,cos_eps   : Ext_1_maxant;
          korr_retn         : Real_1_maxant;

BEGIN

  { Ber„knar avst†nden }
  FOR i:=1 TO (totantal-ant_nypunkt) DO BEGIN
      delta_X2 := Sqr(X_sen - X[i]);
      delta_Y2 := Sqr(Y_sen - Y[i]);
      avstand[i] := Sqrt(delta_X2 + delta_Y2);
  END;

  { Ber„knar korrektionsvinkeln "eps"  (excentricitetsvinkeln) }
  FOR i:=1 TO totantal DO BEGIN
      sin_eps[i] := ((e/avstand[i]) * Sin(exc_vinkel[i]/rho));
      cos_eps[i] := Sqrt(1 - Sqr(sin_eps[i]));
      IF ( cos_eps[i] > 0.00000005 ) THEN
         eps[i] := Arctan(sin_eps[i] / cos_eps[i])
      ELSE IF ( cos_eps[i] <= 0.00000005 ) AND ( sin_eps[i] < 0 ) THEN
         eps[i] := (3*pi)/2
      ELSE
         eps[i] := pi/2;
  END;

  { Ber„knar centrerade och reducerade vinklar }
  j := 1;
  FOR i:=1 TO totantal DO BEGIN
      korr_retn[i] := retn[i] + (eps[i]*rho);
      IF Abs(retn[i]) < 0.0005 THEN
         j:=i;
  END;
  FOR i:=1 TO totantal DO BEGIN
      res[i] := korr_retn[i] - korr_retn[j];
      IF ( res[i] < 0 ) THEN BEGIN
         res[i] := res[i] + 400;
      END;
  END;
END; { Beregn }


{*********************************************************}


PROCEDURE Resultat;

VAR       i,rad                 : INTEGER;


BEGIN

  ClrScr; { Se kommentar vid ”ppnandet av "Crt Unit" ovan }

  WRITELN;
  WRITELN;
  WRITELN ('****************************************************');
  WRITELN ('**                                                **');
  WRITELN ('**    Resultat av excentricitetsber„kningarna     **');
  WRITELN ('**                                                **');
  WRITELN ('****************************************************');
  WRITELN;
  WRITE ('Excentricitetskorrektion');
  rad := WhereY;
  FOR i:=1 TO totantal DO BEGIN
      GotoXY(2,(rad+i));
      WRITE ('Sikt till ',siktepunkt[i],' : ');
      GotoXY(30,(rad+i));
      WRITE ((rho*eps[i]):8:5);
  END;
  WRITELN;
  WRITELN;
  WRITE ('Reducerade och centrerade riktningar.');
  rad := WhereY;
  FOR i:=1 TO totantal DO BEGIN
      GotoXY(2,(rad+i));
      WRITE ('Sikt till ',siktepunkt[i],' : ');
      GotoXY(30,(rad+i));
      WRITE (res[i]:10:5);
  END;
END; { Resultat }


{*********************************************************}
{*********************************************************}


BEGIN     { HUVUDPROGRAM }

  Textbackground(Blue);   { Se kommentar vid ”ppnandet av "Crt Unit" ovan }
  ClrScr;

  Start;
  IF alternativ = 1 THEN
     Inles_alternativ1
  ELSE
     Inles_alternativ2;
  Kontroll_inles;
  Beregn;
  Resultat;

END.