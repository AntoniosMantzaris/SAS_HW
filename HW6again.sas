data iq ;
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

PROC MIXED DATA=iqf METHOD=TYPE3;
	CLASS supplement;
	MODEL iq = supplement/SOLUTION;
	LSMEANS supplement / DIFF=control('suppC') ADJUST=DUNNETT CL;
RUN;


*q6.2;
LIBNAME SAS"/folders/myfolders";
DATA original;
	SET SAS.IVF;
	keep TRT IMP ID;
RUN;
/*
PROC SORT DATA=original;
BY ID TRT PER;
RUN;
PROC MEANS DATA=original NOPRINT;
VAR IMP;
BY ID TRT;
OUTPUT OUT=MEANS MEAN=VOLUME;
RUN;
*/	
PROC MIXED DATA=original METHOD=TYPE3;
CLASS TRT ID;
MODEL IMP = TRT/SOLUTION;			*NOT DDFM=SAT (quite different results)(no CI, no unbalanced, no higher);
RANDOM ID /SOLUTION; 
LSMEANS TRT /diff=control adjust=tukey ;
RUN;



*q6.3;
DATA original;
	SET SAS.RCT;
	where TIME=1;
	keep TRT RESP CENTER;
run;
PROC MIXED DATA=original METHOD=TYPE3;
CLASS TRT CENTER;
MODEL RESP = TRT CENTER/SOLUTION;		
LSMEANS TRT CENTER /diff adjust=tukey ;
RUN;
PROC MIXED DATA=original METHOD=TYPE3 CL;
CLASS TRT;
MODEL RESP = TRT /SOLUTION CL;
lsmeans TRT/diff cl;    *lsmeans-->estimates summarized   diff-->estimates of differences between groups;
RUN;

PROC MIXED DATA=original METHOD=TYPE3;
CLASS TRT CENTER;
MODEL RESP = TRT CENTER/SOLUTION;		
LSMEANS TRT CENTER /diff adjust=dunnett ;
RUN;



*q6.5;
DATA original;
	SET SAS.RCT;
	where TIME=6;
	keep ID;
run;


