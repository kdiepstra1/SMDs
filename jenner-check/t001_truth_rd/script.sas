/* Derived from: Diepstra.SMD exercise with estimating equations_11Mar26.sas
   (lines 11-40): the "truth" data-generation step, its X*Y cross-tab, and the
   identity-link binomial risk difference. &ntruth reduced so the step runs
   self-contained in the hosted engine; the row cap in autoexec.sas keeps the
   captured run small. Data-generation logic, model form, and PROC options are
   the author's, unchanged. */

%let ntruth=5000;                                    *sample size for generating truth;

/*Create Truth*/
data truth;
	call streaminit(9000);
	 do i=1 to &ntruth;
		X=rand("BERNOULLI",0.5);
		Z1=rand("BERNOULLI",0.5);
		Z2=rand("BERNOULLI",0.75);
		Z3=rand("BERNOULLI",0.90);
		Z4=rand("NORMAL", 0, 1); /* Standard normal distribution */ /*Currently does not affect Y*/
		Z5=rand("NORMAL", 0, 10); /* Normally distributed continuous variable with high variance */ /*Currently does not affect Y*/
		LnYprob = log(0.5) + Z1*log(1.4) + Z2*log(0.6) + Z3*log(0.3) + Z4*log(0.2) + Z5*log(0.02) ;
		Yprob = exp(LnYprob) / (1 + exp(LnYprob));
		Y = rand("BERNOULLI",Yprob);
		output;
	end;
run;

proc freq data=truth; /*percentages stayed the same to one decimal place using three different seeds*/
 table x*y ;
run;

*Risk difference;
proc genmod data=truth descending ;
 class x (param=ref ref="0") ;
 model y = x / dist=binomial link=identity ;
 ods output parameterestimates=truthrd ;
run;
