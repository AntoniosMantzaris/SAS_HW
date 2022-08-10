LIBNAME SAS"/folders/myfolders";
DATA original;
	SET SAS.RCT;
	where ID<=75;
RUN;
PROC SORT DATA=original;
BY ID TRT TIME;
RUN;
PROC MIXED DATA=original METHOD=TYPE3 COVTEST CL;
CLASS ID TRT;
MODEL RESP=TRT/SOLUTION DDFM=SAT cl;		*DDFM=SAT for CI;
RANDOM ID(TRT);			*subject nested within oncologist;
RUN;
PROC MIXED DATA=original METHOD=TYPE3 CL;
CLASS  TRT;
MODEL RESP=TRT/SOLUTION cl;		
RUN;
PROC MIXED DATA=original METHOD=TYPE3 COVTEST CL;
CLASS ID TIME;
MODEL RESP=TIME/SOLUTION outp=pred;		
RANDOM ID;			
RUN;
PROC UNIVARIATE DATA=pred NORMAL;		*reject normality-->friedman;
VAR Resid;
PROBPLOT Resid /NORMAL(MU=est SIGMA=est);
HISTOGRAM Resid /NORMAL(MU=est SIGMA=est);
RUN;
proc freq data=pred;
table Resid;
run;
PROC FREQ DATA=original;
TABLES ID*TIME*RESP / CMH2 SCORES=RANK NOPRINT;
RUN;



*q7.2;
DATA original;
	SET SAS.IVF;
RUN;
PROC SORT DATA=original;
BY ID TRT PER;
RUN;
PROC MIXED DATA=original METHOD=TYPE3 COVTEST CL;
CLASS ID TRT PER ;
MODEL IMP=TRT PER TRT*PER/SOLUTION DDFM=SAT outp=pre CL;		*DDFM=SAT for CI;
RANDOM ID(TRT)/solution;			*subject nested within oncologist;
RUN;

PROC UNIVARIATE DATA=pre NORMAL;		*reject normality-->friedman;
VAR Resid;
PROBPLOT Resid /NORMAL(MU=est SIGMA=est);
HISTOGRAM Resid /NORMAL(MU=est SIGMA=est);
RUN;
proc freq data=pre;
table Resid;
run;

proc glm data=pre;
class TRT PER;
model resid =TRT*PER;			
means TRT*PER / hovtest=bartlett;			*or whatever the fixed the fixed effects are represented above;
run;

PROC MIXED DATA=original METHOD=TYPE3 COVTEST CL;
CLASS ID TRT PER ;
MODEL IMP=TRT PER TRT*PER/SOLUTION DDFM=SAT outp=pre CL;		*DDFM=SAT for CI;
RANDOM ID(TRT)/solution;			
lsmeans TRT/diff cl;
RUN;
PROC MIXED DATA=original METHOD=TYPE3 COVTEST CL;
CLASS ID TRT PER ;
MODEL IMP=TRT PER TRT*PER/SOLUTION DDFM=SAT outp=pre CL;		*DDFM=SAT for CI;
RANDOM ID(TRT)/solution;			
slice PER*TRT  / sliceby(per='10') diff cl; 
RUN;



*q7.3;
DATA original3;
	SET SAS.RCT;
RUN;
PROC SORT DATA=original3;
BY ID TRT TIME CENTER;
RUN;
PROC MIXED DATA=original3 METHOD=TYPE3 CL;
CLASS TRT TIME CENTER ID;
MODEL RESP = TRT TIME CENTER TRT*TIME TRT*CENTER TIME*CENTER TRT*TIME*CENTER/SOLUTION DDFM=SAT CL;
RANDOM ID(TRT*CENTER)/solution cl;
/*slice PER*TRT  / sliceby(per='10') diff cl; */     *instead of lsmeans when i want to include only measurments;
RUN;	

data new ;
set original3 ;
where CENTER =1;
if TIME < 4 then SENSOR = 1;
else SENSOR = 2;
keep ID RESP SENSOR ;
run;
PROC SORT DATA=new;
BY ID SENSOR ;
RUN;
PROC MIXED DATA=new METHOD=TYPE3 CL;
CLASS ID SENSOR;
MODEL RESP = ID /SOLUTION DDFM=SAT CL;
RANDOM SENSOR(ID)/solution cl;
RUN;

PROC MIXED DATA=new METHOD=TYPE3 CL;
CLASS ID SENSOR;
MODEL RESP = ID /SOLUTION DDFM=SAT ;
RANDOM SENSOR SENSOR*ID;
RUN;

PROC MIXED DATA=new METHOD=TYPE3 CL;
CLASS ID SENSOR;
MODEL RESP = SENSOR /SOLUTION DDFM=SAT ;
RANDOM ID SENSOR*ID;
RUN;




