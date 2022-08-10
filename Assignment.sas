LIBNAME SAS"/folders/myfolders";
DATA original;
	SET SAS.bchc_2dmt00;
RUN;

proc sort data=original;
by Race_Ethnicity;
run;
proc mixed data=original method=TYPE3 cl;
	CLASS Race_Ethnicity;
	MODEL mortality = /SOLUTION cl OUTP=PRED;
	RANDOM Race_Ethnicity /SOLUTION;				*initial anova rejects null hypothesis;
run;
PROC UNIVARIATE DATA=PRED NORMAL;
VAR RESID;
PROBPLOT RESID /NORMAL(MU=est SIGMA=est);
HISTOGRAM RESID /NORMAL(MU=est SIGMA=est);
RUN;



*outlier test;

proc means data=PRED mean std n q1 q3;
var RESID;
output out=pred2 mean=mean median=median std=std n=n q1=p25 q3=p75 ;
run;
data PRED;
set PRED;
if _n_=1 then set pred2;
drop _TYPE_ _FREQ_;
run;

PROC BOXPLOT DATA=PRED;
  PLOT RESID*MEAN/BOXSTYLE=SCHEMATIC;
RUN;

DATA TUKEY;
SET PRED;
IQR=p75-p25;
LOWERT = p25 - 1.5*IQR;
UPPERT = p75 + 1.5*IQR;
RUN;

DATA TUKEY;
SET TUKEY;
T=(RESID>UPPERT OR RESID<LOWERT);
RUN;

PROC SORT data=TUKEY;
by descending T;
run;

PROC PRINT DATA=TUKEY;
RUN;

DATA PRED_new;
set TUKEY;
where T NE 1;
run;

proc freq data=PRED_new;    *check for ties-->ties are not significant-->Shapiro-Wilk.;
table Resid;
run;

PROC UNIVARIATE DATA=PRED_new NORMAL;
VAR RESID;
HISTOGRAM RESID /NORMAL;   *normality of residuals is not rejected;
RUN;

proc glm data=PRED_new;
class Race_Ethnicity;
model RESID = Race_Ethnicity;
means Race_Ethnicity / hovtest=BARTLETT;		*reject homogeneity;
run;

proc mixed data=PRED_new method=TYPE3 cl;
	CLASS Race_Ethnicity;
	MODEL mortality = /SOLUTION cl ;
	RANDOM Race_Ethnicity /SOLUTION;				*anova without the outliers rejects null hypothesis;
run;


*doornbos for outliers;
proc means data=PRED mean std n;
var RESID;
output out=door mean=mean median=median std=std n=n;
run;
data PRED;
set PRED;
if _n_=1 then set door;
drop _TYPE_ _FREQ_;
run;
data DOORNBOS;
	SET PRED;
	U = (RESID-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T', DOORNBOS_Y, N-2))*N,1);
RUN;

PROC SORT DATA =Doornbos;
	BY descending DOORNBOS_Y;
RUN;

proc print data=Doornbos;
run;

proc freq data=PRED;    *check for ties in mortality to exclude the outlier;
table mortality;
run;

*2nd outlier  ;
DATA PREDnew;
set PRED;
where  mortality NE 70.4;
run;

proc means data=PREDnew mean std n;
var RESID;
output out=door mean=mean median=median std=std n=n;
run;
data PREDnew;
set PREDnew;
if _n_=1 then set door;
drop _TYPE_ _FREQ_;
run;
data DOORNBOS;
	SET PREDnew;
	U = (RESID-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T', DOORNBOS_Y, N-2))*N,1);
RUN;

PROC SORT DATA =Doornbos;
	BY descending DOORNBOS_Y;
RUN;

proc print data=Doornbos;
run;
*3rd outlier ;
DATA PREDnew2;
set PREDnew;
where  mortality NE 80.7;
run;

proc means data=PREDnew2 mean std n;
var RESID;
output out=door mean=mean median=median std=std n=n;
run;
data PREDnew2;
set PREDnew2;
if _n_=1 then set door;
drop _TYPE_ _FREQ_;
run;
data DOORNBOS;
	SET PREDnew2;
	U = (RESID-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T', DOORNBOS_Y, N-2))*N,1);
RUN;

PROC SORT DATA =Doornbos;
	BY descending DOORNBOS_Y;
RUN;

proc print data=Doornbos;
run;   
*3rd;
DATA PREDnew3;
set PREDnew2;
where  mortality NE 71.3;
run;

proc means data=PREDnew3 mean std n;
var RESID;
output out=door mean=mean median=median std=std n=n;
run;
data PREDnew3;
set PREDnew3;
if _n_=1 then set door;
drop _TYPE_ _FREQ_;
run;
data DOORNBOS;
	SET PREDnew3;
	U = (RESID-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T', DOORNBOS_Y, N-2))*N,1);
RUN;

PROC SORT DATA =Doornbos;
	BY descending DOORNBOS_Y;
RUN;

proc print data=Doornbos;
run; 
*4th;
DATA PREDnew4;
set PREDnew3;
where  mortality NE 70.6;
run;

proc means data=PREDnew4 mean std n;
var RESID;
output out=door mean=mean median=median std=std n=n;
run;
data PREDnew4;
set PREDnew4;
if _n_=1 then set door;
drop _TYPE_ _FREQ_;
run;
data DOORNBOS;
	SET PREDnew4;
	U = (RESID-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T', DOORNBOS_Y, N-2))*N,1);
RUN;

PROC SORT DATA =Doornbos;
	BY descending DOORNBOS_Y;
RUN;

proc print data=Doornbos;
run;            *there is no outlier;



proc freq data=PREDnew4;    *check for ties-->ties are not significant-->Shapiro-Wilk.;
table Resid;
run;

PROC UNIVARIATE DATA=PREDnew4 NORMAL;
VAR RESID;
HISTOGRAM RESID /NORMAL;   *normality of residuals is rejected;
RUN;

*assumptions of anova rejected with 2 ways;

*transform data;
data original2;
set original;
y=log(mortality);
run;

proc mixed data=original2 method=TYPE3 cl;
	CLASS Race_Ethnicity;
	MODEL y = /SOLUTION cl OUTP=PRED_log;
	RANDOM Race_Ethnicity /SOLUTION;
run;
PROC UNIVARIATE DATA=PRED_log NORMAL;
VAR RESID;
PROBPLOT RESID /NORMAL(MU=est SIGMA=est);
HISTOGRAM RESID /NORMAL(MU=est SIGMA=est);
RUN;					*not normal;

proc glm data=PRED_log;
class Race_Ethnicity;
model RESID = Race_Ethnicity;
means Race_Ethnicity / hovtest=BARTLETT;		*homogeneity is not rejected;
run;



/*
DATA BoxCox;
	SET original;
	min2=(-1/2)*(mortality**(-2)-1);
	min12=(-2)*(mortality**(-1/2)-1);
	zero=log(mortality);
	plus12=2*(mortality**(1/2)-1);
	plus2=(0.5)*((mortality)**(2)-1);     bgainei to log kalytero transformation
RUN;
PROC UNIVARIATE data=BoxCox;
	histogram min2/normal;
	histogram min12/normal;
	histogram zero/normal;
	histogram plus12/normal;
	histogram plus2/normal;
RUN;
*/

proc means data=PRED_log mean std n;
var RESID;
output out=door mean=mean median=median std=std n=n;
run;
data PRED_log;
set PRED_log;
if _n_=1 then set door;
drop _TYPE_ _FREQ_;
run;
data DOORNBOS;
	SET PRED_log;
	U = (RESID-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T', DOORNBOS_Y, N-2))*N,1);
RUN;

PROC SORT DATA =Doornbos;
	BY descending DOORNBOS_Y;
RUN;

proc print data=Doornbos;
run;

*2nd;
DATA PREDlognew;
set PRED_log;
where  mortality NE 70.4;
run;

proc means data=PREDlognew mean std n;
var RESID;
output out=door mean=mean median=median std=std n=n;
run;
data PREDlognew;
set PREDlognew;
if _n_=1 then set door;
drop _TYPE_ _FREQ_;
run;
data DOORNBOS;
	SET PREDlognew;
	U = (RESID-MEAN)/STD;
	W = SQRT((N*(N-2)*U**2)/((N-1)**2-N*U**2));
	DOORNBOS_Y=ABS(W);
	CRITER= QUANTILE('T', 1-0.05/(2*N), N-2);
	P= MIN(2*MIN(CDF('T', DOORNBOS_Y, N-2),1-CDF('T', DOORNBOS_Y, N-2))*N,1);
RUN;

PROC SORT DATA =Doornbos;
	BY descending DOORNBOS_Y;
RUN;

proc print data=Doornbos;
run;

PROC UNIVARIATE DATA=PREDlognew NORMAL;
VAR RESID;
PROBPLOT RESID /NORMAL(MU=est SIGMA=est);
HISTOGRAM RESID /NORMAL(MU=est SIGMA=est);
RUN;														*normality rejected again;

proc npar1way DATA=original wilcoxon;
class Race_Ethnicity;						*kruskal-wallis rejects null hypothesis;
var mortality;
run;

*EBLUPS with CI; 
proc sort data=original;
by Race_Ethnicity;
run;
ods output SolutionR=EBLUPdata;
proc mixed data=original method=TYPE3 cl;
	CLASS Race_Ethnicity;
	MODEL mortality = /SOLUTION;
	RANDOM Race_Ethnicity /SOLUTION cl;				
run;


*****2;
proc sort data= original;
by  Place Race_Ethnicity;
run;
PROC MIXED DATA=original METHOD=TYPE3;
CLASS Place Race_Ethnicity;
MODEL mortality=/SOLUTION outp=re;
RANDOM Race_Ethnicity Place;

RUN;
PROC UNIVARIATE DATA=re NORMAL;
VAR RESID;
PROBPLOT RESID /NORMAL(MU=est SIGMA=est);
HISTOGRAM RESID /NORMAL(MU=est SIGMA=est);
RUN;
	
*2way non parametric-->fieldman;

/*
*******2';
proc sort data= original;
by  Place Race_Ethnicity;
run;													*ICC: Quantifies how much of the variability is explained by the random effect;
PROC MIXED DATA=original METHOD=TYPE3;					*ICC=0.55 after, ICC=0.49 before. For normal data-->assignment 2;
CLASS  Place Race_Ethnicity ;								*or 0.71;
MODEL mortality=Race_Ethnicity/SOLUTION;				*!!!NO;
RANDOM Place;											*without taking into account the variability of places;
LSMEANS Race_Ethnicity /DIFF=control adjust=tukey;		*only asian with white are the "same" while after taking;
RUN;													*into account places we add all-asian and all-hispanic;
*/



*outliers;
proc means data=re mean std n q1 q3;
var RESID;
output out=pred3 mean=mean median=median std=std n=n q1=p25 q3=p75 ;
run;
data re;
set re;
if _n_=1 then set pred3;
drop _TYPE_ _FREQ_;
run;

PROC BOXPLOT DATA=re;
  PLOT RESID*MEAN/BOXSTYLE=SCHEMATIC;
RUN;

DATA TUKEY;
SET re;
IQR=p75-p25;
LOWERT = p25 - 1.5*IQR;
UPPERT = p75 + 1.5*IQR;
RUN;

DATA TUKEY;
SET TUKEY;
T=(RESID>UPPERT OR RESID<LOWERT);
RUN;

PROC SORT data=TUKEY;
by descending T;
run;

PROC PRINT DATA=TUKEY;
RUN;

DATA re_new;
set TUKEY;
where T NE 1;
run;


PROC UNIVARIATE DATA=re_new NORMAL;
VAR RESID;					*skewness=0.06465698  kurtosis=0.02631656;
RUN;


DATA approx;

N=224;
G1=0.06465698;
G2=0.02631656;
b1=(N-2)*G1/(sqrt(N*(N-1)));
b2=G2*((N-2)*(N-3))/((N+1)*(N-1))+3*(N-1)/(N+1);
*JB=N*(b1**2/6+(G2-3)**2/24);
Cn=(3*(N**2+27*N-70)*(N+1)*(N+3))/((N-2)*(N+5)*(N+7)*(N+9));
Wn2=-1+SQRT(2*(Cn-1));
Alphan=SQRT(2/(Wn2-1));
Dn=1/sqrt(log(sqrt(Wn2)));
Bn=sqrt((N+1)*(N+3)/(6*(N-2)))*b1;
Ts=Dn*log(Bn/Alphan+sqrt(1+(Bn/Alphan)**2));
Mun=3*(N-1)/(N+1);
Sigman=sqrt((24*N*(N-2)*(N-3))/((N+3)*(N+5)*(N+1)**2));
Gamma1n=((6*(N**2-5*N+2))/((N+7)*(N+9)))*sqrt(6*(N+3)*(N+5)/(N*(N-2)*(N-3)));
An=6+(8/(Gamma1n))*(2/Gamma1n+sqrt(1+4/(Gamma1n**2)));
Un=(b2-Mun)/Sigman;
Tk=sqrt(9*An/2)*((9*An-2)/(9*An)-((1-2/An)/(1+Un*sqrt(2/(An-4))))**(1/3));
K2=Tk**2+Ts**2;   *combination of skewness and kurtosis tests-->chi square (but only for large n);
Ps=2*min(cdf('Normal',Ts,0,1),1-cdf('Normal',Ts,0,1));
Pk=2*min(cdf('Normal',Tk,0,1),1-cdf('Normal',Tk,0,1));
PK2=1-cdf('chisq',K2,2);
*PJB=1-cdf('chisq',JB,2);
run;

proc print data=approx;
var Ts Tk K2 Ps Pk PK2 ;    *normality of residuals not rejected;
run;

PROC MIXED DATA=re_new METHOD=TYPE3;
CLASS Place Race_Ethnicity;
MODEL mortality=/SOLUTION outp=re;
RANDOM Race_Ethnicity Place;

RUN;




/*
PROC MIXED DATA=original METHOD=TYPE3;    *fixed-fixed;
CLASS Place Race_Ethnicity ;
MODEL mortality=Place Race_Ethnicity/SOLUTION;

LSMEANS Place Race_Ethnicity/DIFF=control adjust=tukey;
RUN;
*/


******3;
proc sort data=original;
by Place;
run;
PROC MIXED DATA=original METHOD=TYPE3;    
CLASS Place  ;
MODEL mortality=Place /SOLUTION ;

LSMEANS Place /DIFF=control adjust=tukey;		*we want diversity-->keep Sesttle, Las Vegas, San Jose, according to estimates;
RUN;								*We have no reason to question the accurateness of the randomization(assignment2);
									*Cities should have data from the same years since there could be diifferences;
									*occured by the passing of years (however this assumption is rejected);
								



proc sort data=original;
by Year;
run;
PROC MIXED DATA=original METHOD=TYPE3;    
CLASS Year  ;								*do not reject null hypothesis, means are equal every year;
MODEL mortality=Year /SOLUTION;

LSMEANS Year /DIFF=control adjust=tukey;	
RUN;