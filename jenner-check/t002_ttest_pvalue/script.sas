/* Derived from: Diepstra.SMD exercise with estimating equations_11Mar26.sas
   (lines 99-133): the %pvalue2 macro that computes, per simulated replicate,
   the pooled two-sample t-test p-value for each continuous covariate, then
   stacks the results. The macro body, its %scan/%countw variable loop, the
   BY-replicate PROC TTEST, the ODS OUTPUT of the TTESTS table, and the
   post-processing DATA steps are the author's, unchanged. Only adaptation: a
   small "simulate" dataset is generated inline (3 replicates) using the same
   call streaminit + rand draws as the upstream simulation step, so the macro
   runs self-contained; obs=100 cap in autoexec keeps the captured run small. */

%let convars=Z4 Z5 ;                                 *list of continuous covariates;

/*Simulate data (a few replicates, same generation as the upstream simulate step)*/
data simulate;
	call streaminit(654321);
	do replicate = 1 to 3;
	 do i=1 to 40;
		X=rand("BERNOULLI",0.5);
		Z4=rand("NORMAL", 0, 1); /* Standard normal distribution */
		Z5=rand("NORMAL", 0, 10); /* Normally distributed continuous variable with high variance */
		output;
	  end;
	end;
run;

/* P-values -- CONTINUOUS VARS */
%macro pvalue2(data, vars);

 %local var i;
 %do i=1 %to %sysfunc(countw(&vars));
  %let var=%scan(&vars,&i);

    proc sort data=&data ; by replicate ; run;
	ods exclude all;
    proc ttest data=&data ;
     by replicate ;
	 class x ;
     var &var ;
	 ods output TTESTS=ttest_results_&var ;
    run;
	ods select all;

	data ttest_results_&var ;
	 set ttest_results_&var (keep = replicate method probt) ;
	 where method = "Pooled" ;
     length covariate $ 3;
	 covariate = "&var" ;
	 rename probt = pvalue ;
	 drop method ;
	run;

 %end;

%mend;

%pvalue2(simulate, &convars);


data con_pvalue ;
 set ttest_: ;
run;

proc print data=con_pvalue; run;
