LIBNAME SAS"/folders/myfolders";
DATA original;
	SET SAS.IVF;
	where PER=4;
RUN;
proc means data=original MEAN VAR STD N;
var AGEM;
run;
proc univariate data=original cibasic;
	var AGEM;
	histogram AGEM/ normal;
run;
PROC IML;
use original;
read all var{agem};
close WEEK1;
alpha=0.05;
Ybar=mean(agem);
s=var(agem);
n=nrow(agem);
qT=quantile('t',alpha/2,n-1);
UCL=Ybar-qT*sqrt(s/n);
LCL=Ybar+qT*sqrt(s/n);
A=Ybar||LCL||UCL;
create DATA from A[colname={'mean''LCL''UCL'}];
append from A;
close DATA;
quit;

PROC IML;
use original;
read all var{agem};
close WEEK1;
alpha=0.05;
Ybar=mean(agem);
s=var(agem);
n=nrow(agem);
qT=quantile('t',alpha/2,n-1);
UPL=Ybar-qT*sqrt((n+1)*s/n);
LPL=Ybar+qT*sqrt((n+1)*s/n);
A=Ybar||LPL||UPL;
create DATA from A[colname={'mean''LPL''UPL'}];
append from A;
close DATA;
quit;


DATA original2;
	set original;
	count= (AGEM>=40);
RUN;
proc means data=original2 sum;
var count;
run;


PROC IML;
use original;
read all var{agem};
close WEEK1;
alpha=0.05;
s=var(agem);
n=nrow(agem);
qCL=quantile('chisquare',alpha/2,n-1);
qCU=quantile('chisquare',1-alpha/2,n-1);
UCL=(n-1)*s/qCL;
LCL=(n-1)*s/qCU;
sd=sqrt(s);
UCLsd=sqrt((n-1)*s/qCL);
LCLsd=sqrt((n-1)*s/qCU);
A=(s||LCL||UCL)//(sd||LCLsd||UCLsd);
print(s||LCL||UCL);
create SD from A[colname={'statistic''LCL''UCL'}];
append from A;
close SD;
quit;


*q1.3;
proc means data=original MEAN STD N P25 P75 QRANGE;    *IQR;
var BW;
output out=DATAOUT MEAN=MEAN STD=STD P25=P25 P75=P75 Qrange=IQR;
run;

PROC IML;
use original;
read all var{BW};
close original;
alpha=0.05;
p1=0.25;
s1=p1*(1-p1);
p2=0.75;
s2=p2*(1-p2);
n=nrow(BW);
z=quantile('Normal',1-alpha/2);
pU1=p1+z*sqrt(s1/n);
pL1=p1-z*sqrt(s1/n);
nU1=min(floor(n*pU1)+1,n);
nL1=max(1,floor(n*pL1));
pU2=p2+z*sqrt(s2/n);
pL2=p2-z*sqrt(s2/n);
nU2=min(floor(n*pU2)+1,n);
nL2=max(1,floor(n*pL2));
call sort(BW);
call qntl(pct1, BW, p1);
call qntl(pct2, BW, p2);
LCL1=BW[nL1];
UCL1=BW[nU1];
LCL2=BW[nL2];
UCL2=BW[nU2];
A=(pct1||LCL1||UCL1||nL1||nU1||pct2||LCL2||UCL2||nL2||nU2);
create PCTL from A[colname={'pct1' 'LCL1' 'UCL1' 'LR1' 'UR1' 'pct2' 'LCL2' 'UCL2' 'LR2' 'UR2'}];
append from A;
close PCTL;
quit;

proc print data=DATAOUT;
	var MEAN IQR;
RUN;
data count2;
	set original;
	count= (BW<=4235 and BW>=2495);
run;
proc means data= count2 sum;
	var count;
RUN;


DATA WEEK1BOXCOX;
SET original;
	AGEMMINUS2= (-1/2)*(BW**-2 -1);
	AGEMMINUS1= (-1)*(BW**-1 -1);
	AGEMMINUS12= (-2)*(BW**-(0.5)-1);
	AGEM0= log(BW);
	AGEM13= (3)*(BW**(1/3) -1);
	AGEM12= (2)*(BW**(1/2) -1);
	AGEM2= (0.5)*(BW**(2) -1);
RUN;
proc univariate data=WEEK1BOXCOX;
	histogram AGEMMINUS2 /normal;
	histogram AGEMMINUS1 /normal;
	histogram AGEMMINUS12 /normal;
	histogram AGEM0 /normal;
	histogram AGEM13 /normal;
	histogram AGEM12 /normal;
	histogram AGEM2 /normal;
	histogram AGEM /normal;
run;

PROC IML;
use WEEK1BOXCOX;
read all var{AGEM2};
close WEEK1BOXCOX;
alpha=0.05;
Ybar=mean(AGEM2);
s=var(AGEM2);
n=nrow(AGEM2);
qT=quantile('t',alpha/2,n-1);
UPL=Ybar-qT*sqrt((n+1)*s/n);
LPL=Ybar+qT*sqrt((n+1)*s/n);
u=sqrt(2*UPL+1);
l=sqrt(2*LPL+1);
A=Ybar||l||u;
create DATA from A[colname={'mean' 'LPL' 'UPL'}];
append from A;
close DATA;
quit;


proc sort data=original;
	by descending BW;
RUN;


data boys;
	set original;
	where SEX=1;
run;
data girls;
	set original;
	where SEX=0;
run;
proc means data=boys N MEAN VAR SKEW KURT;   
var BW;
run;
proc means data=girls N MEAN VAR SKEW KURT;   
var BW;
run;


data boysBox;
	set boys;
	AGEM2=(0.5)*(BW**(2) -1);			*PREDICTION with BoxCox-->set the variable before;
RUN;
PROC IML;
use boysBox;
read all var{AGEM2};
close boys;
alpha=0.05;
Ybar=mean(AGEM2);
s=var(AGEM2);
n=nrow(AGEM2);
qT=quantile('t',alpha/2,n-1);
UPL=Ybar-qT*sqrt((n+1)*s/n);
LPL=Ybar+qT*sqrt((n+1)*s/n);
u=sqrt(2*UPL+1);
l=sqrt(2*LPL+1);
A=Ybar||l||u;
create DATA from A[colname={'mean' 'LPL' 'UPL'}];
append from A;
close DATA;
quit;



*q1.6;
DATA original;
	SET SAS.IVF;
	where PER=4;
	age=log(44-GA);
	keep age;
RUN;
PROC IML;
use original;
read all var{age};
close original;
alpha=0.05;
Ybar=mean(age);
s=var(age);
n=nrow(age);
qT=quantile('t',alpha/2,n-1);
UPL=Ybar-qT*sqrt((n+1)*s/n);
LPL=Ybar+qT*sqrt((n+1)*s/n);
u=44-exp(UPL);
l=44-exp(LPL);
A=Ybar||u||l;
create DATA from A[colname={'mean' 'LPL' 'UPL'}];
append from A;
close DATA;
quit;

PROC IML;
use original;
read all var{age};
close original;
alpha=0.05;
Ybar=mean(age);
s=var(age);
n=nrow(age);
qT=quantile('t',alpha/2,n-1);
UCL=Ybar-qT*sqrt(s/n);
LCL=Ybar+qT*sqrt(s/n);
A=Ybar||LCL||UCL;
create DATA from A[colname={'mean' 'LCL' 'UCL'}];
append from A;
close DATA;
quit;

PROC IML;
use original;
read all var{age};
close original;
alpha=0.05;
p=0.5;
s=p*(1-p);
n=nrow(age);
z=quantile('Normal',1-alpha/2);
pU=p+z*sqrt(s/n);
pL=p-z*sqrt(s/n);
nU=min(floor(n*pU)+1,n);
nL=max(1,floor(n*pL));
call sort(age);
call qntl(pct, age, p);
LCL=age[nL];
UCL=age[nU];
u=44-exp(UCL);							*we can return to initial values for CI of percentiles;
l=44-exp(LCL);
A=(pct||u||l||nL||nU);
create PCTL from A[colname={'pctl' 'LCL' 'UCL' 'LR' 'UR'}];
append from A;
close PCTL;
quit;

DATA original;
	SET SAS.IVF;
	where PER=4;
	age=log(44-GA);
	
	keep age GA FIS ;
RUN;
data original;
	set original;
	test=(GA<=38);
	count1=(FIS=1 AND test=1);
	count2=(FIS=1 AND test=0);
run;
proc means data=original n sum;
	var count1 count2;
run;

data testing;
	set original;
	where test=1;
RUN;
data stress1;
	set original;
	where FIS=1 AND test=1;
	
RUN;
data stress2;
	set original;
	where FIS=1 AND test=0;
	
RUN;

PROC IML;
use testing;
read all var{FIS};         *agemb=(agem<30);  *ποσοστό γυναικών κάτω από 30;
close testing; 
alpha=0.1;
Ybar=mean(FIS);
s=Ybar*(1-Ybar);
n=nrow(FIS);
z=quantile('Normal',1-alpha/2);
UCL=Ybar+z*sqrt(s/n);
LCL=Ybar-z*sqrt(s/n);
A=Ybar||LCL||UCL;
print(Ybar||LCL||UCL);
create PROP from A[colname={'p''LCL''UCL'}];
append from A;
close PROP;
quit;

data testing2;				*ξεκαθαρα ονοματα για να καταλαβαινω τι να χρησιμοποιω( οχι οπως εδω);
	set original;
	where test=0;
RUN;
PROC IML;
use testing2;
read all var{FIS};         
close testing; 
alpha=0.1;
Ybar=mean(FIS);
s=Ybar*(1-Ybar);
n=nrow(FIS);
z=quantile('Normal',1-alpha/2);
UCL=Ybar+z*sqrt(s/n);
LCL=Ybar-z*sqrt(s/n);
A=Ybar||LCL||UCL;
print(Ybar||LCL||UCL);
create PROP from A[colname={'p''LCL''UCL'}];
append from A;
close PROP;
quit;

proc freq data=testing2;
tables FIS /binomial(wald wilson exact level=2) alpha=0.1;		*exact method for small n;
run;


*q1.8;
DATA original;
	SET SAS.IVF;
	where PER=4;
RUN;

proc ttest data=original h0= 3200 sides=2 alpha=0.05;
	var BW;     *μ(f)=3200;
run;
proc ttest data=original h0= 3200 sides=U alpha=0.05;			*U for H1:μ>3200;
	var BW;     
run;

%macro samples(dataset=,ns=,n=);
                proc surveyselect data=&dataset NOPRINT
   			method=urs n=&n out=Final;
		run;
		
		Data Final;
		set final;
		sampleno=1;
		run;

               
               %do sn = 2 %to &ns;
                        proc surveyselect data=&dataset NOPRINT
   			method=urs n=&n out=SampleI;
			run;
			
			Data SampleI;
			Set SampleI;
			sampleno= &sn;
			run;

			Data Final;
			Set Final SampleI;
			run;
     		%end;
     		
     		proc datasets library=work NOPRINT;
     		delete SampleI;
     		run;
%MEND;

%samples(dataset=original,ns=1000,n=253);
proc means data = FINAL  mean NOPRINT ;
var BW ;
by sampleno ;
output out = MEANSBW mean = BW_MEAN ;
run ;

proc univariate data = MEANSBW ;
hist BW_MEAN / normal ;
run ;





