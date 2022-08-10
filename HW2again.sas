LIBNAME SAS"/folders/myfolders";
DATA original;
	SET SAS.IVF;
	where PER=4;
RUN;

proc ttest data=original;
class FIS;
var AGEM;
run;

proc glm data=original;
class FIS;
model AGEM = FIS;
means FIS / hovtest=BARTLETT;
run;
proc glm data=original;
class FIS;
model AGEM = FIS;
means FIS / hovtest=LEVENE;
run;

*Q2.2;
DATA original2;
	SET SAS.IVF;
	where PER=4 AND ID<=100;
RUN;

data w01;
set original2;
where TRT=0 OR TRT=1;
run;
data w02;
set original2;
where TRT=0 OR TRT=2;
run;

ods output WilcoxonScores=WRS (keep= Class N SumOfScores);
proc npar1way data=w01 correct=NO;
class TRT;
var BW;
exact wilcoxon    ;      
run;
proc npar1way data=w02 correct=NO;
class TRT;
var BW;
exact wilcoxon /mc   ;      
run;

PROC IML;
use WRS;
read all var{N SumOfScores};
close WRS;
G={1 , 2};
U=SumOfScores-N#(N+1)/2;
P=U/prod(N);
A=G||N||U||P;
create MWU from A[colname={'Group' 'N' 'U' 'P'}];			*if they ask for statistic give the 'U1';
append from A;												*p-values equal to wilcoxon's;
close MWU;
quit;


proc npar1way data=w01 correct=NO;
class TRT;
var BW;
exact KS /mc   ;      
run;

data transf;
set original2;
plus2=(0.5)*((BW)**(2)-1);
run;
data t01;
set transf;
where TRT=0 OR TRT=1;
RUN;
proc ttest data=t01;
class TRT;
var plus2;
run;


*q2.5;
proc format ;
   value TREATMENTFmt 0='Control'
                1='Treatment';
   value QUANTITYFmt 0='Low'
                1='High';
run;
data LTT;
   input TREATMENT QUANTITY COUNT;
   datalines;
0 0  77
0 1  23
1 0  81
1 1  19
;

proc freq data=LTT order=data;
   format TREATMENT TREATMENTFmt. QUANTITY QUANTITYFmt.;
   tables TREATMENT*QUANTITY / chisq;
   weight COUNT;
run;


*q2.7;
DATA original;
	SET SAS.IVF;
	where PER=4;
RUN;
data w01;
set original;
where TRT=0 OR TRT=1;
run;
data w02;
set original;
where TRT=0 OR TRT=2;
run;

Proc freq data=w01;
Table TRT*FIS / chisq;   
Exact chisq;
RUN;

Proc freq data=w01;
Table TRT*FIS / chisq;
Exact fisher;
RUN;


proc freq data=original;
tables TRT*FIS /chisq;
run;
