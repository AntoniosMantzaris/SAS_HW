LIBNAME SAS"/folders/myfolders";
DATA original;
	SET SAS.IVF;
	where PER=4;
	keep ID AGEM;
RUN;
PROC UNIVARIATE DATA=original NORMAL;
VAR AGEM;
RUN;
proc freq data=original;
table AGEM;
run;

PROC UNIVARIATE DATA=original NORMAL;
	VAR AGEM;
RUN;
proc iml;
DATA approx;
*N=712;
*G1=0.43510643;
*G2=0.51594979;
N=253;
G1=-0.178018;
G2=	-0.1849355;
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
K2=Tk**2+Ts**2;
Ps=2*min(cdf('Normal',Ts,0,1),1-cdf('Normal',Ts,0,1));
Pk=2*min(cdf('Normal',Tk,0,1),1-cdf('Normal',Tk,0,1));
PK2=1-cdf('chisq',K2,2);
*PJB=1-cdf('chisq',JB,2);
run;

proc print data=approx;
var Ts Tk K2 Ps Pk PK2 ;
run;



*q5.2;
DATA Q13;
	set SAS.IVF;
	where PER=4;
	keep BW ID;
RUN;


proc means data=q13 mean std n q1 q3;
var BW;
output out=Q52 mean=mean median=median std=std n=n q1=p25 q3=p75 ;
run;
data Q13;
set Q13;
if _n_=1 then set Q52;
drop _TYPE_ _FREQ_;
run;


DATA Grubbs;
	SET Q13;
	U = (BW-MEAN)/STD;
	Grubbs=ABS(U);
	
	t = quantile("t", 0.05 / (2*N), N-2);
	C_twosided_approx = (n-1) * sqrt(t**2 / (n * (t**2 + n - 2)));
	u_inv = u*sqrt((n-2)) / sqrt(n-1-u**2);
	p_twosided_approx = min(2*n*MIN(1-cdf("t", u_inv, n-2),cdf("t", u_inv, n-2)), 1);
RUN;
PROC SORT DATA =Grubbs;
	BY descending Grubbs;
RUN;
proc print data=Grubbs;   *largest grubbs has p-value=0.33-->not outlier;
run;

PROC UNIVARIATE DATA=Q13 NORMAL;
VAR BW;
RUN;
proc freq data=Q13;
table BW;
run;

*Tukey's method;
PROC BOXPLOT DATA=Q13;
  PLOT BW*MEAN/BOXSTYLE=SCHEMATIC;
RUN;

DATA TUKEY;
SET Q13;
IQR=p75-p25;
LOWERT = p25 - 1.5*IQR;
UPPERT = p75 + 1.5*IQR;
RUN;

DATA TUKEY;
SET TUKEY;
T=(BW>UPPERT OR BW<LOWERT);
RUN;

PROC SORT data=TUKEY;
by descending T;
run;

PROC PRINT DATA=TUKEY;
RUN;

data Q13NEW;
SET TUKEY;
WHERE T NE 1;
run;
PROC UNIVARIATE DATA=Q13NEW NORMAL;
VAR BW;
RUN;
proc freq data=Q13NEW;
table BW;
run;

DATA Q13;
	set SAS.IVF;
	where PER=4;
	keep BW ID;
RUN;
data q13;
set q13;
plus2=(1/2)*(BW**2-1);
run;
proc means data=q13 mean std n q1 q3;
var plus2;
output out=Q52 mean=mean median=median std=std n=n q1=p25 q3=p75 ;
run;
data Q13;
set Q13;
if _n_=1 then set Q52;
drop _TYPE_ _FREQ_;
run;


DATA Grubbs;
	SET Q13;
	U = (plus2-MEAN)/STD;
	Grubbs=ABS(U);
	
	t = quantile("t", 0.05 / (2*N), N-2);
	C_twosided_approx = (n-1) * sqrt(t**2 / (n * (t**2 + n - 2)));
	u_inv = u*sqrt((n-2)) / sqrt(n-1-u**2);
	p_twosided_approx = min(2*n*MIN(1-cdf("t", u_inv, n-2),cdf("t", u_inv, n-2)), 1);
RUN;
PROC SORT DATA =Grubbs;
	BY descending Grubbs;
RUN;
proc print data=Grubbs;   *largest grubbs has p-value=0.33-->not outlier;
run;

PROC UNIVARIATE DATA=Q13 NORMAL;
VAR plus2;
RUN;
proc freq data=Q13;
table plus2;
run;



*q5.3;
DATA q16;
	SET SAS.IVF;
	where PER=4;
	keep ID GA;
RUN;
DATA Q16a;
	SET Q16;
	g=log(44-GA);
RUN;
*Doornbos test;
proc means data=q16 mean std n;
	var GA;
	output out=L10_sumstat mean=mean median=median std=std n=n;
run;

data q16;
set q16;
if _n_=1 then set L10_sumstat;
drop _TYPE_ _FREQ_;
run;

data DOORNBOS;
	SET q16;
	U = (GA-MEAN)/STD;
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

data L10_NEW;
SET q16;
WHERE ID NE 174;
run;
*Doornbos test;
proc means data=L10_NEW mean std n;
	var GA;
	output out=L10_sum mean=mean median=median std=std n=n;
run;

data L10_NEW;
set L10_NEW;
if _n_=1 then set L10_sum;
drop _TYPE_ _FREQ_;
run;

data DOORNBOS;
	SET L10_NEW;
	U = (GA-MEAN)/STD;
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

PROC UNIVARIATE DATA=L10_NEW NORMAL;
VAR GA;
RUN;
proc freq data=L10_NEW;
table GA;
run;


DATA q16;
	SET SAS.IVF;
	where PER=4;
	keep ID GA;
RUN;
*Hampel's rule;
proc means data=q16 mean std n;
	var GA;
	output out=L10_sumstat mean=mean median=median std=std n=n;
run;
data q16;
set q16;
if _n_=1 then set L10_sumstat;
drop _TYPE_ _FREQ_;
run;

DATA Hampel;
	SET q16;
	D = ABS(GA-MEDIAN);
RUN;

proc means data=Hampel;
	var D;
	output out=Hampel_med median=medianD;
run;

data Hampel;
set Hampel;
if _n_=1 then set Hampel_med;
drop _TYPE_ _FREQ_;
run;

DATA Hampel;
	SET Hampel;
	Z = ABS(GA-MEDIAN)/MEDIAND;
	H = (Z>3.5);
RUN;


Proc sort data=Hampel;
by descending H;
run;

PROC PRINT DATA=Hampel;
RUN;

data Q16NEW;
SET HAMPEL;
WHERE H NE 1;
run;

PROC UNIVARIATE DATA=Q16NEW NORMAL;
VAR GA;
RUN;
proc freq data=Q16NEW;
table GA;
run;


*Hampel's rule;
proc means data=q16a mean std n;
	var g;
	output out=L10_sumstat mean=mean median=median std=std n=n;
run;
data q16a;
set q16a;
if _n_=1 then set L10_sumstat;
drop _TYPE_ _FREQ_;
run;

DATA Hampel;
	SET q16a;
	D = ABS(g-MEDIAN);
RUN;

proc means data=Hampel;
	var D;
	output out=Hampel_med median=medianD;
run;

data Hampel;
set Hampel;
if _n_=1 then set Hampel_med;
drop _TYPE_ _FREQ_;
run;

DATA Hampel;
	SET Hampel;
	Z = ABS(g-MEDIAN)/MEDIAND;
	H = (Z>3.5);
RUN;


Proc sort data=Hampel;
by descending H;
run;

PROC PRINT DATA=Hampel;
RUN;

data Q16NEW;
SET HAMPEL;
WHERE H NE 1;
run;

PROC UNIVARIATE DATA=Q16NEW NORMAL;
VAR g;
RUN;
proc freq data=Q16NEW;
table g;
run;



*q5.4;
DATA original;
	SET SAS.IVF;
	where PER=4;
	keep ID IMP TRT;
RUN;
PROC MIXED DATA=original METHOD=TYPE3 CL;
CLASS TRT;
MODEL IMP = TRT /SOLUTION  outp=pred CL;
*lsmeans TRT/diff cl;    *lsmeans-->estimates summarized   diff-->estimates of differences between groups;
RUN;

PROC UNIVARIATE DATA=pred NORMAL;
VAR resid;
RUN;
proc freq data=pred;
table resid;
run;
proc glm data=pred;
class TRT;
model resid = TRT;
means TRT / hovtest= levene;
run;



*q5.5;
DATA original;
	SET SAS.RCT;
	where TIME=1;
RUN;
proc mixed data=original method=TYPE3 cl;
class ID;
model RESP = /solution outp=re cl;
random ID/solution;						*ICC=σ_G^2/(σ_G^2+σ_E^2);
run;
PROC UNIVARIATE DATA=re NORMAL;
VAR resid;
RUN;
proc freq data=re;
table resid;
run;



