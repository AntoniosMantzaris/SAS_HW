*Runs tests;
%MACRO RUNSCUC(data=,var=,alpha=);
PROC IML;
use &data;
read all var {&var};
close &data;

X=&var;
X=X[loc(cmiss(X)<1)];
n=nROW(X);
MED=median(X);

XC=X;
DO i=1 to n by 1;
	IF (XC[i] >= MED) then XC[i]=1;
	ELSE XC[i]=0;
END;

n1C=sum(XC);
n2C=n-n1C;

RC=1;
DO i=2 to n by 1;
	if(XC[i] ^= XC[i-1]) then RC=RC+1;
END;

MUC=1+(2*n1C*n2C)/(n1C+n2C);
VARC=2*n1C*n2C*(2*n1C*n2C-n1C-n2C)/((n1C+n2C-1)*(n1C+n2C)**2);

SC=(RC-MUC)/SQRT(VARC);
TC=QUANTILE('NORMAL',&alpha/2);
TCU=QUANTILE('NORMAL',1-&alpha/2);
PC=(1-CDF('NORMAL',abs(SC)))*2;

XUC=REPEAT(0,n-1,1);
TIES=0;
DO i=1 to (n-1) by 1;
	IF (X[i+1] > X[i]) then XUC[i]=1;
	IF (X[i+1] < X[i]) then XUC[i]=0;
	IF (X[i+1] = X[i]) then XUC[i]=XUC[i-1];
	IF (X[i+1] = X[i]) then TIES=TIES+1;
END;

RUC=1;
DO i=2 to (n-1) by 1;
	if(XUC[i] ^= XUC[i-1]) then RUC=RUC+1;
END;

MUUC=(2*(n-TIES)-1)/3;
VARUC=(16*(n-TIES)-29)/90;

SUC=(RUC-MUUC)/SQRT(VARUC);
TUC=QUANTILE('NORMAL',&alpha/2);
TUCU=QUANTILE('NORMAL',1-&alpha/2);
PUC=(1-CDF('NORMAL',abs(SUC)))*2;

PRINT("Median based (conditional) runs test");
PRINT(RC||MUC||sqrt(VARC)||PC||SC||TC||TCU||n);
PRINT("(unconditional) runst test for serial randomness");
PRINT(TIES);
PRINT(RUC||MUUC||sqrt(VARUC)||PUC||SUC||TUC||TUCU||(n-TIES));
quit;
%MEND;



data iq;
input supplement$ iq@@ ;
datalines ;
suppB 104.3 suppB 99.0 suppB 112.5 suppB 114.0
suppB 132.4 suppB 109.4 suppB 98.8 suppB 98.9
suppB 112.4 suppB 101.9 suppB 97.0 suppB 112.1
suppB 100.7 suppB 100.5 suppB 114.8 suppB 100.6
suppB 105.3 suppB 110.5 suppB 110.7 suppB 119.3
suppA 97.0 suppA 106.7 suppA 108.1 suppA 97.1
suppA 96.7 suppA 105.4 suppA 105.6 suppA 110.0
suppA 106.3 suppA 99.7 suppA 108.4 suppA 106.3
suppA 109.5 suppA 99.8 suppA 93.7 suppA 107.7
suppA 102.7 suppA 106.3 suppA 97.7 suppA 107.3
;
PROC MIXED DATA=iq METHOD=TYPE3 CL;
CLASS supplement;
MODEL iq = supplement /SOLUTION outp=pred CL;
lsmeans supplement/ cl;   
RUN;
proc ttest data=iq;
class supplement;
var iq;
run;
proc npar1way DATA=iq wilcoxon;
class supplement;
var iq;
run;
PROC UNIVARIATE DATA=PRED NORMAL;		*MODEL RESP = TRT/SOLUTION OUTP=PRED CL;
VAR RESID;
PROBPLOT RESID /NORMAL(MU=est SIGMA=est);
HISTOGRAM RESID /NORMAL(MU=est SIGMA=est);
RUN;

data iq2 ;
input supplement$ iq@@ ;
datalines ;
suppC 103.3 suppC 104.0 suppC 117.5 suppC 119.0
suppC 135.4 suppC 113.4 suppC 103.8 suppC 103.9
suppC 115.4 suppC 106.9 suppC 102.0 suppC 117.1
suppC 105.7 suppC 105.5 suppC 119.8 suppC 105.6
suppC 110.3 suppC 115.5 suppC 115.7 suppC 124.3
;
DATA iqf;
	SET iq iq2;
run;
PROC MIXED DATA=iqf METHOD=TYPE3 CL;
CLASS supplement;
MODEL iq = supplement /SOLUTION CL;
lsmeans supplement/diff cl;   
RUN;



*q4.4;
data COAG ;
input Patient C K@@;
datalines ;
1 120 132 8 145 133 15 117 123
2 114 116 9 120 123 16 125 108
3 129 135 10 129 116 17 136 131
4 128 115 11 126 127 18 151 119
5 155 134 12 136 140 19 130 129
6 105 56 13 135 140 20 136 124
7 114 114 14 125 114 21 113 112
;
run;
proc sort data=coag;
	by Patient;
run;
proc transpose data=COAG out=COAG2;
   by Patient;
run;
data COAG2;
	SET COAG2(rename=(col1=measurement));
run;

ods output SolutionR=EBLUPdata;	
proc mixed data=COAG2 method=TYPE3 cl;
class Patient;
model measurement = /solution outp=pred cl;
random Patient/solution;						*ICC=σ_G^2/(σ_G^2+σ_E^2);
run;

proc sort data=EBLUPdata;
by Patient;
run;
%RUNSCUC(data=EBLUPdata,var=Estimate,alpha=0.05);
proc sort data=pred;
by Patient;
run;
%RUNSCUC(data=pred,var=Resid,alpha=0.05);


*q4.6;
LIBNAME SAS"/folders/myfolders";
DATA original;
	SET SAS.RCT;
	where TIME=1;
RUN;
proc mixed data=original method=TYPE3 cl;
class ID;
model RESP = /solution outp=re cl;
random ID/solution;						*ICC=σ_G^2/(σ_G^2+σ_E^2);
run;
proc sort data=re;
by CENTER;
run;
%RUNSCUC(data=re,var=Resid,alpha=0.05);



*q4.9;
PROC IMPORT OUT= WORK.Gron
  DATAFILE= "/folders/myfolders/Groningen.csv"
    DBMS=CSV REPLACE;
    GETNAMES=YES;
RUN;
DATA q9;
SET Gron;
Ybar=0.027*(X**2)-0.275*X;
diff=Y-Ybar;
run;
%RUNSCUC(data=q9,var=diff,alpha=0.05);




