LIBNAME SAS"/folders/myfolders";
DATA original1;
	SET SAS.IVF;
RUN;
DATA original2;
	SET SAS.RCT;
RUN;


%Macro SIM_Gum(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

do i=1 to &nsim by 1;
U1=rand('Uniform');   *random values from uniform distribution in (a,b) with prob 1/|b-a|, here a=1,b=0;
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha);
return(Exp(-((-Log(x))**alpha + (-Log(U1))**alpha)**(1/alpha))*((-Log(x))**alpha + (-Log(U1))**alpha)**(-1 + 1/alpha)*((-Log(U1))**(alpha-1))/U1-U2);
finish;

intervals = {0.00001 1};        
U2C = froot("Func", intervals);  *finds zeros of the univariate function "func" 
-->converge to a root when given an interval in which the function changes signs;

X=X//U1;
Y=Y//U2C;
YI=YI//U2;
end;

Total=X||Y||YI;

create GumC from Total [colname={'X','Y','YI'}]; 
append from Total;       
close GumC;
quit;
%mend SIM_Gum;
%SIM_Gum(nsim=1000, alpha=5, seed=12345);


%Macro SIM_Frk(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

do i=1 to &nsim by 1;
U1=rand('Uniform');
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha);
return((Exp(alpha)*(-1 + Exp(alpha*x)))/(-Exp(alpha) + Exp(alpha*(1+x)) - Exp(alpha*(U1+x)) + Exp(alpha*(1 + U1)))-U2);
finish;

intervals = {0.00001 1};        
U2C = froot("Func", intervals);

X=X//U1;
Y=Y//U2C;
end;

Total=X||Y;

create FrkC from Total [colname={'X','Y'}]; 
append from Total;       
close FrkC;
quit;
%mend SIM_Frk;


%SIM_Frk(nsim=1000, alpha=-5, seed=12345);


*Spearman's rank correlation;
proc corr data=GumC spearman;
var X Y;
run;

*Kendall's tau;
proc corr data=GumC kendall;
var X Y;
run;



%Macro SpearmanRho(rho=);
proc iml;
pi = constant("pi");
tau=&rho;


start innerGum(y) global(alpha, x); 
   return(Exp(-((-Log(x))**alpha + (-Log(y))**alpha)**(1/alpha)));
finish; 

start outerGum(par) global(x,alpha); 
	x=par;
   yinterval = 0 || 1;
   /** evaluate inner integral for the parameter value, a=x **/ 
   call quad(w, "innerGum", yinterval);
   return (w);
finish; 

start finalGum(param) global(alpha, tau);
alpha=param;
xinterval= {0 1};
call quad(v, "outerGum", xinterval); /** outer integral **/ 
return(12*v-(3+tau));
finish;


/*
t = do(1, 100, 1);          
F = j(1, ncol(t));
do i = 1 to ncol(t);
   F[i] = finalGum(t[i]);      
end;
title "Integral - 3+ tau";
call Series(t, F) grid={x y} label={"x" "F(x)"}
                  other="refline 0 / axis=y";  
*/

intervalsGum = {1 100};        
SGum = froot("finalGum", intervalsGum);
print(SGum);

start innerClay(y) global(alpha, x); 
	return((x**(-alpha)+y**(-alpha)-1)**(-1/alpha));
finish; 

start outerClay(par) global(x, alpha); 
	x=par;

if(alpha>0) then yinterval= 0||1;
else yinterval= (1-x**(-alpha))**(-1/alpha)||1;
   /** evaluate inner integral for the parameter value, a=x **/ 
   call quad(w, "innerClay", yinterval);
   return (w);
finish; 

start finalClay(param) global(alpha, tau);
alpha=param;
xinterval= {0 1};
call quad(v, "outerClay", xinterval); /** outer integral **/ 
return(12*v-(3+tau));
finish;

/*
t = do(-1, 3, 0.11);          
F = j(1, ncol(t));
do i = 1 to ncol(t);
   F[i] = finalClay(t[i]);      
   print(finalClay(t[i]));
end;
title "Integral - 3+ tau";
call Series(t, F) grid={x y} label={"x" "F(x)"}
                  other="refline 0 / axis=y";  
*/
                 
intervalsClay = {-1 10};        
SClay = froot("finalClay", intervalsClay);
print(SClay);

SGau=2*sin(pi*tau/6);
print(SGau);

SFGM=3*tau;
print(SFGM);

start innerFrk(y) global(alpha, x); 
return(-(1/alpha)*Log(1+(Exp(-alpha*x)-1)*(Exp(-alpha*y)-1)/(Exp(-alpha)-1)));
finish; 

start outerFrk(par) global(x, alpha); 
	x=par;
   yinterval = 0 || 1;
   /** evaluate inner integral for the parameter value, a=x **/ 
   call quad(w, "innerFrk", yinterval);
   return (w);
finish; 

start finalFrk(param) global(alpha, tau);
alpha=param;
xinterval= {0 1};
call quad(v, "outerFrk", xinterval); /** outer integral **/ 
return(12*v-(3+tau));
finish;

/*

t = do(-10, 10, 1.1);          
F = j(1, ncol(t));
do i = 1 to ncol(t);
   F[i] = finalFrk(t[i]);      
  
end;
title "Integral - 3+ tau";
call Series(t, F) grid={x y} label={"x" "F(x)"}
                  other="refline 0 / axis=y";  
*/
                 
intervalsFrk = {-30 30};        
SFrk = froot("finalFrk", intervalsFrk);
print(SFrk);


CPAR=SGum||SClay||SFrk||SGau||SFGM;

create EstSpearman from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
append from CPAR;       
close EstSpearman;
quit;

%mend SpearmanRho;

%SpearmanRho(rho=0.95131);


***Estimate copula parameters using Kendall's Tau;
%Macro KendallTau(tau=);
proc iml;
pi = constant("pi");
tau=&tau;

SGum=1/(1-tau);
print(SGum);
SClay=2*tau/(1-tau);
print(SClay);
SGau=sin(pi*tau/2);
print(SGau);
SFGM=9*(tau/2);
print(SFGM);

start D(y);
return(y/(Exp(y)-1));
finish;

*IF alpha>0 / tau>0;
start FC(x) global(tau);
dinterval=0||x;
call quad(w, "D", dinterval);
return(1-(4/x)*(1-(1/x)*w)-tau);
finish;

intervals = {0.00001 20};        
SFrk = froot("FC", intervals);
print(SFrk);

/*

*IF alpha<0 / tau<0;
start FC(x) global(tau);
dinterval=0||-x;
call quad(w, "D", dinterval);
return(1-(4/x)*(1+(1/x)*w+0.5*x)-tau);
finish;

intervals = {-10 -0.00001};        
SFrk = froot("FC", intervals);
print(SFrk);

*/


CPAR=SGum||SClay||SFrk||SGau||SFGM;

create EstKendall from CPAR [colname={'Gumbel alpha','Clayton alpha','Frank alpha','Gaussian rho','FGM alpha'}]; 
append from CPAR;       
close EstKendall;
quit;

%mend KendallTau;

%KendallTau(tau=0.81311);


%SIM_Gum(nsim=1000, alpha=5, seed=12345);		*original;
proc sgplot data=GumC aspect=1;
scatter x=X y=Y;
run;
%SIM_Gum(nsim=1000, alpha=5.38, seed=12345);		*after using the kendall's tau computed from the original;
%SIM_Frk(nsim=1000, alpha=19, seed=12345);
proc sgplot data=GumC aspect=1;
scatter x=X y=Y;
run;
proc sgplot data=FrkC aspect=1;
scatter x=X y=Y;
run;


*q3.4;
DATA original2;
	SET SAS.RCT;
RUN;

PROC TRANSPOSE OUT=WIDE_RCT(DROP = _NAME_ _LABEL_) DATA=original2
PREFIX=RESP;
BY ID;
ID TIME;
VAR RESP;
RUN;

*Pearson's Correlation;
proc corr DATA=WIDE_RCT pearson;
var RESP1 RESP2;
run;
*Spearman's rank correlation;
proc corr data=WIDE_RCT spearman;
var RESP1 RESP2;
run;

*Kendall's tau;
proc corr data=WIDE_RCT kendall;
var RESP1 RESP2;
run;


proc corr data=WIDE_RCT pearson fisher(biasadj=no);
var RESP1 RESP2;
run;


%SpearmanRho(rho= 0.50477 );
%KendallTau(tau=0.35941);

*FGM Copula;
%Macro SIM_FGM(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;
do i=1 to &nsim by 1;
U1=rand('Uniform');
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha);
return(x*(1 + alpha*(1 - x)*(1 - U1)) - alpha*(1 - x)*x*U1-U2);
finish;

intervals = {0.00001 1};        
U2C = froot("Func", intervals);

X=X//U1;
Y=Y//U2C;
YI=YI//U2;
end;

Total=X||Y||YI;

create FGMC from Total [colname={'X','Y','YI'}]; 
append from Total;       
close FGMC;
quit;
%mend SIM_FGM;

%SIM_FGM(nsim=716, alpha=1.51431, seed=12345);		
%SIM_FGM(nsim=716, alpha=1.617345, seed=12345);		*ολα αυτα μπορουν να παραληφθουν αφου το α στο FGM;
proc sgplot data=WIDE_RCT aspect=1;					*πρεπει να ειναι μεταξυ -1,1;
scatter x=RESP1 y=RESP2;
run;
proc sgplot data=GumC aspect=1;
scatter x=X y=Y;
run;
proc sgplot data=FrkC aspect=1;
scatter x=X y=Y;
run;


*Clayton Copula;
%Macro SIM_Clay(alpha=, nsim=, seed=);
proc iml;
call streaminit(&seed);
alpha=&alpha;

do i=1 to &nsim by 1;
U1=rand('Uniform');
U2=rand('Uniform');

start Func(x) global(U1,U2,alpha);
return(U1**(-1 -alpha)*(x**(-alpha) + U1**(-alpha)-1)**(-1 - 1/alpha)-U2);
finish;

intervals = {0.001 1};        
U2C = froot("Func", intervals);

X=X//U1;
Y=Y//U2C;
end;

Total=X||Y;

create CC from Total [colname={'X','Y'}]; 
append from Total;       
close CC;
quit;
%mend SIM_Clay;


%KendallTau(tau=0.35941);
%SpearmanRho(rho= 0.50477 );




*ORIGINAL DATASET SIMULATION;
DATA original22;
SET original2;
RESP = RESP + (ranuni(1)-0.5);
KEEP ID TIME RESP;
run;

PROC TRANSPOSE OUT=WIDE32(DROP = _NAME_ _LABEL_) DATA=original22 PREFIX=RESP;
BY ID;
ID TIME;
VAR RESP;
RUN;
DATA WIDE32;
SET WIDE32;
if cmiss(of _all_) then delete;
run;
proc rank data=WIDE32 out=ranked_a;
var RESP1 RESP2;
ranks rank_RESP1 rank_RESP2;
run;
proc means data=ranked_a N;
var rank_RESP1 rank_RESP2;
run;
data marginals;
set ranked_a;
U_RESP1=rank_RESP1/695;
U_RESP2=rank_RESP2/695;
run;

proc sgplot data=marginals aspect=1;					
scatter x=U_RESP1 y=U_RESP2;
run;

%SIM_Clay(nsim=700, alpha=1.1221218, seed=12345);
%SIM_Frk(nsim=700, alpha=3.6264203, seed=12345);
proc sgplot data=CC;
scatter x=X y=Y;
run;
proc sgplot data=FrkC;
scatter x=X y=Y;
run;



*q3.6;
DATA original2;
	SET SAS.RCT;
	where CENTER=1 AND (TIME=5 OR TIME=6);
RUN;
PROC TRANSPOSE OUT=WIDE6(DROP = _NAME_ _LABEL_) DATA=original2 PREFIX=RESP;
BY ID;
ID TIME;
VAR RESP;
RUN;
DATA q6;
SET WIDE6;
dif=RESP6-RESP5;
logdif=log(RESP6)-log(RESP5);
rat=RESP6/RESP5;
RUN;

PROC UNIVARIATE DATA=q6 NORMAL;
	VAR dif logdif rat;
	HISTOGRAM dif logdif rat;
RUN; 



*q3.9;
DATA original1;
	SET SAS.IVF;
	where PER=10 OR PER=18;
	keep ID IMP ind PER;
RUN;
PROC TRANSPOSE OUT=WIDE9(DROP = _NAME_ _LABEL_) DATA=original1 PREFIX=IMP;
BY ID;
ID PER;
VAR IMP;
RUN;
data WIDE9 ;		*for mcnemar;
set WIDE9 ;
if cmiss (of IMP10 IMP18 ) then delete ;
run;

/*
DATA q9;
SET WIDE11;
dif=IMP18-IMP10;					t,sign,wlk-sign
logdif=log(IMP18)-log(IMP10);
rat=IMP18/IMP10;
RUN;

PROC UNIVARIATE DATA=q9 NORMAL;
	VAR dif logdif rat;
	HISTOGRAM dif logdif rat;
RUN; 
*/
DATA WIDE9;
set WIDE9;
IMP10=IMP10<85;
IMP18=IMP18<85;
RUN;
*McNemar;
PROC FREQ DATA=WIDE9;
	TABLES IMP10*IMP18/ AGREE;
	RUN;
*Exact McNemar;
PROC FREQ DATA=WIDE9;
	TABLES IMP10*IMP18; EXACT MCNEM;
RUN;



*q11;
proc format ;
   value CONTROLFmt 0='Low'
                1='High';
   value POSTFmt 0='Low'
                1='High';
run;
data LTT;
   input control post COUNT;
   datalines;
0 0  70
0 1  7
1 0  11
1 1  12
;
proc freq data=LTT order=data;
   format control CONTROLFmt. post POSTFmt.;
   tables control*post / agree;
   weight COUNT;
run;
