/*Karen Diepstra*/
/*Simulation: SMD estimates with 95% confidence intervals*/

/* For SMD estimating equations for categorical variables, see: "SMD (95% CIs) CALCULATION -- CATEGORICAL VARS" */
/* For SMD estimating equations for continuous variables, see: "SMD (95% CIs) CALCULATION -- CONTINUOUS VARS" */
/* See Ross et al. PMCID: PMC10904145 for PROC IML code (other than SMD EEs) */

options nofmterr;

/*Sample size and covariates*/
%let sampsize=/*100*/ 500 /*2000*/ /*10000*/ ;        *sample size for each replicate;          
%let rep=1000;                                        *number of simulated replicates;
%let catvars=Z1 Z2 Z3 ;                               *list of binary covariates;
%let convars=Z4 Z5 ;                                  *list of continuous covariates;
%let ntruth=20000000;                                 *sample size for generating truth;


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

/*Simulate data*/
data simulate;
	call streaminit(654321); 
	do replicate = 1 to &rep;
	 do i=1 to &sampsize;
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
	end;
run;


/* P-values -- CATEGORICAL VARS */
/* not reported */

%macro pvalue(data, vars);

 %local var i;                                            /*https://communities.sas.com/t5/SAS-Programming/How-can-I-loop-a-series-of-procedures-to-run-through-a-list-of/td-p/763397*/
 %do i=1 %to %sysfunc(countw(&vars));
  %let var=%scan(&vars,&i);

    proc sort data=&data ; by replicate ; run;
	ods exclude all;
    proc freq data=&data ;
     by replicate ;
     table x*(&var)/chisq warn=output ;
	 OUTPUT OUT = CHI_&VAR N PCHI;
    run;
	ods select all;

	data chi_&var ;
	 set chi_&var (keep = replicate p_pchi) ;
	 length variable $ 3;
	 variable = "&var" ;
	 rename p_pchi = pvalue ;
	run;

 %end;

%mend;

%pvalue(simulate, &catvars);


data cat_pvalue ;
 set chi_: ;
 rename variable = covariate ;
run;


/* P-values -- CONTINUOUS VARS */
/* not reported */

%macro pvalue2(data, vars);

 %local var i;                                            /*https://communities.sas.com/t5/SAS-Programming/How-can-I-loop-a-series-of-procedures-to-run-through-a-list-of/td-p/763397*/
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



/* SMD (95% CIs) CALCULATION -- CONTINUOUS VARS */

%macro replicate(rep=,vars=,output=) ;

%do r = 1 %to &rep;

  data replicate ;
   set simulate ;
   where replicate=&r ;
  run;

%local var i;                                            /*https://communities.sas.com/t5/SAS-Programming/How-can-I-loop-a-series-of-procedures-to-run-through-a-list-of/td-p/763397*/
 %do i=1 %to %sysfunc(countw(&vars));
  %let var=%scan(&vars,&i);

PROC IML;                            
	*Read data;
	use replicate;								
		read all var {x} into exposure;		
		read all var {&var} into z;
	close replicate;

	n = nrow(exposure);  

	
	/***********************************************
	Defining estimating equations */

	q = 5;								

	START efunc(theta) global(exposure, z, n);	


		ef_1 = exposure#(z - theta[1]);
		ef_2 = (1-exposure)#(z - theta[2]);
		ef_3 = exposure#((z - theta[1])##2 - theta[3]);
		ef_4 = (1-exposure)#((z - theta[2])##2 - theta[4]); 
		ef_5 = j(n, 1, (theta[1]-theta[2]) - (theta[5]*(sqrt(.5*(theta[3] + theta[4])))));


		ef_mat = ef_1||ef_2||ef_3||ef_4||ef_5; 
		RETURN(ef_mat);                         						
	FINISH efunc;        


	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {0.5,0.5,1,1,0.05};       * Initial parameter values; *value for each parameter;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"mean1","mean0","variance1","variance0","SMD"};  
	est = theta_hat`; 
    variance = vecdiag(sandwich); 
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;				

	PRINT variable est variance se lcl ucl;     	

	CREATE results VAR {variable est variance se lcl ucl};
		APPEND;                         		  	  
	CLOSE results;  
	QUIT;

	data results ;
	 set results ;
	 covariate = "&var" ;
	 replicate = &r ;
	RUN;


 	* CREATE OUTPUT DATASET FOR REPLICATE- *These will be editted depending on what output from you analysis looks like!;
	%if %sysfunc(exist(&output)) %then %do; 
		data &output;
		set &output results;
		run;
		/*POSSIBLE ALTERNATIVE DEPENDING ON RESULTS STRUCTURE;
		data &output;
		merge &output RESULTSFROMANALYSISMACRO;
		by interv ; 
		run;*/
	%end;
	%else %do;
		data &output; set results; run;
	%end;

%end;

%end;

%mend;

%replicate(rep=&rep,vars=&convars,output=smdscon) ;




/* SMD (95% CIs) CALCULATION -- CATEGORICAL VARS */

%macro replicate(rep=,vars=,output=) ;

%do r = 1 %to &rep;

  data replicate ;
   set simulate ;
   where replicate=&r ;
  run;


%local var i;                                            /*https://communities.sas.com/t5/SAS-Programming/How-can-I-loop-a-series-of-procedures-to-run-through-a-list-of/td-p/763397*/
 %do i=1 %to %sysfunc(countw(&vars));
  %let var=%scan(&vars,&i);

PROC IML;                            
	*Read data;
	use replicate;								
		read all var {x} into exposure;		
		read all var {&var} into z;
	close replicate;

	n = nrow(exposure);  

	
	/***********************************************
	Defining estimating equations */

	q = 3;								

	START efunc(theta) global(exposure, z, n);	


		ef_1 = exposure#(z - theta[1]);
		ef_2 = (1-exposure)#(z - theta[2]);
		ef_3 = j(n, 1, (theta[1]-theta[2]) - (theta[3]*sqrt((.5*theta[1]*(1-theta[1])) + (.5*theta[2]*(1-theta[2])))));


		ef_mat = ef_1||ef_2||ef_3; 
		RETURN(ef_mat);                         						
	FINISH efunc;        


	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {.5,.5,0};       * Initial parameter values; *value for each parameter;
	optn = q || 1;                      * Set options for nlplm, (3 - requests 3 roots,1 - printing summary output);
	tc = j(1, 12, .);                   * Create vector for Termination Criteria options, set all to default using .;
	tc[6] = 1e-9;                       * Replace 6th option in tc to change default tolerance;
	CALL nlplm(rc,                      /*Use the Levenberg-Marquardt root-finding method*/
			   theta_hat,                /*... name of output parameters that give the root*/
			   "eequat",                /*... function to find the roots of*/
			   initial,                 /*... starting values for root-finding*/
               optn, ,                  /*... optional arguments for root-finding procedure*/
               tc);                     /*... update convergence tolerance*/

	/***********************************************
	Baking the bread (approximate derivative) */
	par = q||.||.;                   	* Set options for nlpfdd, (3 - 3 parameters, . = default);
	CALL nlpfdd(func,                   /*Derivative approximation function*/
                deriv,                  /*... name of output matrix that gives the derivative*/
                na,                     /*... name of output matrix that gives the 2nd derivative - we do not need this*/
                "eequat",               /*... function to approximate the derivative of*/
                theta_hat,               /*... point where to find derivative*/
                par);                   /*... details for derivative approximation*/ 
	bread = - (deriv) / n;              * Negative derivative, averaged;

	/***********************************************
	Cooking the filling (matrix algebra) */
	residuals = efunc(theta_hat);		* Value of estimating functions at beta hat (n by q matrix);
	outerprod = residuals` * residuals; * Outerproduct of residuals (note transpose is flipped from slides);
	filling = outerprod / n; 				* Divide by n for filling;

	/***********************************************
	Assembling the sandwich (matrix algebra) */
	sandwich = ( inv(bread) * filling * inv(bread)` ) / n;

	/***********************************************
	Formatting results for output */
	variable = {"p1","p0","SMD"};  
	est = theta_hat`; 
    variance = vecdiag(sandwich); 
	se = sqrt(vecdiag(sandwich));       
	lcl = est - 1.96*se; 				
	ucl = est + 1.96*se;

	PRINT variable est variance se lcl ucl;     	

	CREATE results VAR {variable est variance se lcl ucl};
		APPEND;                         		  	  
	CLOSE results;  
	QUIT;

	data results ;
	 set results ;
	 covariate = "&var" ;
	 replicate = &r ;
	RUN;


 	* CREATE OUTPUT DATASET FOR REPLICATE- *These will be editted depending on what output from you analysis looks like!;
	%if %sysfunc(exist(&output)) %then %do; 
		data &output;
		set &output results;
		run;
		/*POSSIBLE ALTERNATIVE DEPENDING ON RESULTS STRUCTURE;
		data &output;
		merge &output RESULTSFROMANALYSISMACRO;
		by interv ; 
		run;*/
	%end;
	%else %do;
		data &output; set results; run;
	%end;

%end;

%end;

%mend;

%replicate(rep=&rep,vars=&catvars,output=smdscat) ;



/* CONCATENATE SMD DATASETS */

data smds ;
 set smdscat smdscon ;
 where variable = "SMD" ;

 *CLD;
 cld = ucl - lcl ;

 *CI coverage of null (0) ;
 if LCL<0 and UCL>0 then null_exc = 0 ;
  else null_exc = 1 ;

 *CI coverage of threshold (0.1) ;
 if LCL>0.1 or UCL<-0.1 then thres_exc = 1 ;
  else thres_exc = 0 ;

 *SMD >0.1 ;
 if est>=0.1 or est<=-0.1 then imbalance = 1 ;
  else imbalance = 0 ;
run;

proc sort data=smds ; by covariate replicate ; run;
proc sort data=cat_pvalue ; by covariate replicate ; run;
proc sort data=con_pvalue ; by covariate replicate ; run;

data results ;
 merge smds cat_pvalue con_pvalue ;
 by covariate replicate ;
run;

data results ;
 set results ;
 if pvalue <0.05 then imbalance_pvalue = 1 ;
 else imbalance_pvalue = 0 ;
run;

proc means data=results mean min max ;
 class covariate ;
 var est cld ;
run;

proc freq data=results ;
 by covariate ;
 table imbalance null_exc thres_exc ;
run;




/* CRUDE RD FOR EACH REPLICATE */

%macro cruderdrr(data=, replicate=);

*delete simrdrr if it exists;
%if %sysfunc(exist(cruderdrr)) %then %do; 
	proc datasets nolist; 
	delete cruderdrr; 
	quit; run;
%end;

 %do iter=1 %to &replicate;

     *Crude Risk difference ;
     proc genmod data=&data descending ;
	  where replicate=&iter ;
      class x (param=ref ref="0") ;
      model y = x  / dist=binomial link=identity ;
      ods output parameterestimates=rd_crude&iter ;
     run;
     
	 data rd_crude&iter;
	  set rd_crude&iter (keep = parameter estimate lowerwaldcl upperwaldcl);
	  where parameter="X";
	  rename estimate = simrd_crude;
	  rename lowerwaldcl = simlowerrd_crude;
	  rename upperwaldcl = simupperrd_crude;
	  replicate = &iter;
	  drop parameter;
	 run;


	* CREATE OUTPUT DATASET;
	%if %sysfunc(exist(cruderdrr)) %then %do; 
		data cruderdrr;
		set cruderdrr rd_crude&iter;
		run;
	%end;
	%else %do;
		data cruderdrr; set rd_crude&iter; run;
	%end;

	*Delete datasets used in analysis macro here so memory is clear before next rep;
	proc datasets nolist; 
	delete rd_crude&iter; 
	quit; run;

  %end;


%mend;

%cruderdrr(data=simulate, replicate=&rep);


/* RD ADJUSTED FOR ALL COVARIATES WITH SMD>0.1 FOR EACH REPLICATE */
/* adjustment done with IPTW */

%macro iptwrdrr(data=, replicate=);

*delete simrdrr if it exists;
%if %sysfunc(exist(iptwrdrr)) %then %do; 
	proc datasets nolist; 
	delete iptwrdrr; 
	quit; run;
%end;

 %do iter=1 %to &replicate;

     *Identify covariates for adjustment;
     data cov ;
	  set results ;
	  where replicate=&iter and (est>=.1 or est<=-.1);
	 run;



	 data list ;   
       length allvars $1000;   
       retain allvars ' ';   
       set cov end=eof;   
       allvars = trim(left(allvars))||' '||left(covariate);   
       if eof then call symput('varlist', allvars);   
     run; 

     %put &varlist;

    /* Check if dataset LIST has zero observations */
	%let dsid = %sysfunc(open(list));
	%let nobs = %sysfunc(attrn(&dsid, NOBS));
	%let rc   = %sysfunc(close(&dsid));

	%if &nobs = 0 %then %do;

	  *Risk difference;
     proc genmod data=&data descending ;
	  where replicate=&iter ;
      class x (param=ref ref="0") ;
      model y = x / dist=binomial link=identity ;
      ods output parameterestimates=rd_rep&iter ;
     run;
     
	 data rd_rep&iter;
	  set rd_rep&iter (keep = parameter estimate lowerwaldcl upperwaldcl);
	  where parameter="X";
	  rename estimate = simrd_iptw1;
	  rename lowerwaldcl = simlowerrd_iptw1;
	  rename upperwaldcl = simupperrd_iptw1;
	  replicate = &iter;
	  drop parameter;
	 run;

	%end;

	%else %do;
    *Probability of exposure;
     proc logistic data=&data descending ;
	  where replicate=&iter ;
      model x = &varlist ;
      output out=ipwdenoms p=pdenom ;
     run;

	 data weighteddata ;
	  set ipwdenoms ;
	  
	  if x = 1 then ipt_weight = 1/pdenom ;
	  if x = 0 then ipt_weight = 1/(1-pdenom) ;
     
	 run;

	 *Risk difference ;
	 proc genmod data=weighteddata descending ;
	 class i x (ref="0" param=ref) ;
	 model y = x / dist=binomial link=identity ;
	 weight ipt_weight ;
	 repeated subject=i/ corr=ind printmle;
	 ods output geeemppest=parms;
	 run;

	 data rd_rep&iter;
	  set parms (keep = parm estimate lowercl uppercl);
	  where parm="X";
	  rename estimate = simrd_iptw1;
	  rename lowercl = simlowerrd_iptw1;
	  rename uppercl = simupperrd_iptw1;
	  replicate = &iter;
	  drop parm;
	 run;

	%end; 



	* CREATE OUTPUT DATASET;
	%if %sysfunc(exist(iptwrdrr)) %then %do; 
		data iptwrdrr;
		set iptwrdrr rd_rep&iter;
		run;
	%end;
	%else %do;
		data iptwrdrr; set rd_rep&iter; run;
	%end;

	*Delete datasets used in analysis macro here so memory is clear before next rep;
	proc datasets nolist; 
	delete rd_rep&iter; 
	quit; run;

  %end;

%mend;

%iptwrdrr(data=simulate, replicate=&rep);



/* RD ADJUSTED FOR ALL COVARIATES WITH SMD 95%CI lower excluding zero FOR EACH REPLICATE */
/* adjustment done with IPTW */

%macro iptw2rdrr(data=, replicate=);

*delete simrdrr if it exists;
%if %sysfunc(exist(iptw2rdrr)) %then %do; 
	proc datasets nolist; 
	delete iptw2rdrr; 
	quit; run;
%end;

 %do iter=1 %to &replicate;

     *Identify covariates for adjustment;
     data cov ;
	  set results ;
	  where replicate=&iter and null_exc=1 ;
	 run;



	 data list ;   
       length allvars $1000;   
       retain allvars ' ';   
       set cov end=eof;   
       allvars = trim(left(allvars))||' '||left(covariate);   
       if eof then call symput('varlist', allvars);   
     run; 

     %put &varlist;

    /* Check if dataset LIST has zero observations */
	%let dsid = %sysfunc(open(list));
	%let nobs = %sysfunc(attrn(&dsid, NOBS));
	%let rc   = %sysfunc(close(&dsid));

	%if &nobs = 0 %then %do;

	  *Risk difference;
     proc genmod data=&data descending ;
	  where replicate=&iter ;
      class x (param=ref ref="0") ;
      model y = x / dist=binomial link=identity ;
      ods output parameterestimates=rd_rep&iter ;
     run;
     
	 data rd_rep&iter;
	  set rd_rep&iter (keep = parameter estimate lowerwaldcl upperwaldcl);
	  where parameter="X";
	  rename estimate = simrd_iptw2;
	  rename lowerwaldcl = simlowerrd_iptw2;
	  rename upperwaldcl = simupperrd_iptw2;
	  replicate = &iter;
	  drop parameter;
	 run;

	%end;

	%else %do;
    *Probability of exposure;
     proc logistic data=&data descending ;
	  where replicate=&iter ;
      model x = &varlist ;
      output out=ipwdenoms p=pdenom ;
     run;

	 data weighteddata ;
	  set ipwdenoms ;
	  
	  if x = 1 then ipt_weight = 1/pdenom ;
	  if x = 0 then ipt_weight = 1/(1-pdenom) ;
     
	 run;

	 *Risk difference ;
	 proc genmod data=weighteddata descending ;
	 class i x (ref="0" param=ref) ;
	 model y = x / dist=binomial link=identity ;
	 weight ipt_weight ;
	 repeated subject=i/ corr=ind printmle;
	 ods output geeemppest=parms;
	 run;

	 data rd_rep&iter;
	  set parms (keep = parm estimate lowercl uppercl);
	  where parm="X";
	  rename estimate = simrd_iptw2;
	  rename lowercl = simlowerrd_iptw2;
	  rename uppercl = simupperrd_iptw2;
	  replicate = &iter;
	  drop parm;
	 run;

	%end; 



	* CREATE OUTPUT DATASET;
	%if %sysfunc(exist(iptw2rdrr)) %then %do; 
		data iptw2rdrr;
		set iptw2rdrr rd_rep&iter;
		run;
	%end;
	%else %do;
		data iptw2rdrr; set rd_rep&iter; run;
	%end;

	*Delete datasets used in analysis macro here so memory is clear before next rep;
	proc datasets nolist; 
	delete rd_rep&iter; 
	quit; run;

  %end;

%mend;

%iptw2rdrr(data=simulate, replicate=&rep);


/* RD ADJUSTED FOR ALL COVARIATES WITH SMD 95%CI lower >0.1 or upper <-0.1 FOR EACH REPLICATE */
/* adjustment done with IPTW */
/* not reported */

%macro iptw3rdrr(data=, replicate=);

*delete simrdrr if it exists;
%if %sysfunc(exist(iptw3rdrr)) %then %do; 
	proc datasets nolist; 
	delete iptw3rdrr; 
	quit; run;
%end;

 %do iter=1 %to &replicate;

     *Identify covariates for adjustment;
     data cov ;
	  set results ;
	  where replicate=&iter and thres_exc=1 ;
	 run;



	 data list ;   
       length allvars $1000;   
       retain allvars ' ';   
       set cov end=eof;   
       allvars = trim(left(allvars))||' '||left(covariate);   
       if eof then call symput('varlist', allvars);   
     run; 

     %put &varlist;

    /* Check if dataset LIST has zero observations */
	%let dsid = %sysfunc(open(list));
	%let nobs = %sysfunc(attrn(&dsid, NOBS));
	%let rc   = %sysfunc(close(&dsid));

	%if &nobs = 0 %then %do;

	  *Risk difference;
     proc genmod data=&data descending ;
	  where replicate=&iter ;
      class x (param=ref ref="0") ;
      model y = x / dist=binomial link=identity ;
      ods output parameterestimates=rd_rep&iter ;
     run;
     
	 data rd_rep&iter;
	  set rd_rep&iter (keep = parameter estimate lowerwaldcl upperwaldcl);
	  where parameter="X";
	  rename estimate = simrd_iptw3;
	  rename lowerwaldcl = simlowerrd_iptw3;
	  rename upperwaldcl = simupperrd_iptw3;
	  replicate = &iter;
	  drop parameter;
	 run;

	%end;

	%else %do;
    *Probability of exposure;
     proc logistic data=&data descending ;
	  where replicate=&iter ;
      model x = &varlist ;
      output out=ipwdenoms p=pdenom ;
     run;

	 data weighteddata ;
	  set ipwdenoms ;
	  
	  if x = 1 then ipt_weight = 1/pdenom ;
	  if x = 0 then ipt_weight = 1/(1-pdenom) ;
     
	 run;

	 *Risk difference ;
	 proc genmod data=weighteddata descending ;
	 class i x (ref="0" param=ref) ;
	 model y = x / dist=binomial link=identity ;
	 weight ipt_weight ;
	 repeated subject=i/ corr=ind printmle;
	 ods output geeemppest=parms;
	 run;

	 data rd_rep&iter;
	  set parms (keep = parm estimate lowercl uppercl);
	  where parm="X";
	  rename estimate = simrd_iptw3;
	  rename lowercl = simlowerrd_iptw3;
	  rename uppercl = simupperrd_iptw3;
	  replicate = &iter;
	  drop parm;
	 run;

	%end; 



	* CREATE OUTPUT DATASET;
	%if %sysfunc(exist(iptw3rdrr)) %then %do; 
		data iptw3rdrr;
		set iptw3rdrr rd_rep&iter;
		run;
	%end;
	%else %do;
		data iptw3rdrr; set rd_rep&iter; run;
	%end;

	*Delete datasets used in analysis macro here so memory is clear before next rep;
	proc datasets nolist; 
	delete rd_rep&iter; 
	quit; run;

  %end;

%mend;

%iptw3rdrr(data=simulate, replicate=&rep);


/* Concatenate for bias estimation */

proc sort data=cruderdrr ; by replicate ; run;
proc sort data=iptwrdrr ; by replicate ; run;
proc sort data=iptw2rdrr ; by replicate ; run;
proc sort data=iptw3rdrr ; by replicate ; run;

data bias ;
merge cruderdrr iptwrdrr iptw2rdrr iptw3rdrr ;
by replicate ;
parameter="X" ;
run;

data truthrd2 ;
 set truthrd ;
 where parameter="X" ;
 rename estimate=truthrd ;
run;

data truthrd2 ;
 set truthrd2 (keep = parameter truthrd) ;
run;

proc sort data=bias ; by parameter ; run;
proc sort data=truthrd2 ; by parameter ; run;

data bias2 ;
 merge bias (in=a) truthrd2 ;
 by parameter ;
 if a ;
run;

/* calculate bias, absolute bias, and mse for each measure: crude, iptw1, and iptw2 */

data bias3 ;
 set bias2 ;

 crude_bias = simrd_crude - truthrd ;
 ab_crude_bias = abs(crude_bias) ;
 mse_crude = (truthrd - simrd_crude)**2 ;

 iptw1_bias = simrd_iptw1 - truthrd ;
 ab_iptw1_bias = abs(iptw1_bias) ;
 mse_iptw1 = (truthrd - simrd_iptw1)**2 ;

 iptw2_bias = simrd_iptw2 - truthrd ;
 ab_iptw2_bias = abs(iptw2_bias) ;
 mse_iptw2 = (truthrd - simrd_iptw2)**2 ;

 iptw3_bias = simrd_iptw3 - truthrd ;
 ab_iptw3_bias = abs(iptw3_bias) ;
 mse_iptw3 = (truthrd - simrd_iptw3)**2 ;

run;



/* mean bias */

proc univariate data=bias3 noprint;
    var crude_bias;
    output out=pctl_crude mean=crude_bias_est pctlpts=2.5 97.5 pctlpre=crude pctlname=_LL _UL;
run;

proc univariate data=bias3 noprint;
    var iptw1_bias;
    output out=pctl_iptw1 mean=iptw1_bias_est pctlpts=2.5 97.5 pctlpre=iptw1 pctlname=_LL _UL;
run;

proc univariate data=bias3 noprint;
    var iptw2_bias;
    output out=pctl_iptw2 mean=iptw2_bias_est pctlpts=2.5 97.5 pctlpre=iptw2 pctlname=_LL _UL;
run;

proc univariate data=bias3 noprint;
    var iptw3_bias;
    output out=pctl_iptw3 mean=iptw3_bias_est pctlpts=2.5 97.5 pctlpre=iptw3 pctlname=_LL _UL;
run;

data percentiles_output;
    merge pctl_crude pctl_iptw1 pctl_iptw2 pctl_iptw3;
run;

proc print data=percentiles_output noobs label;
    var crude_bias_est crude_LL crude_UL
        iptw1_bias_est iptw1_LL iptw1_UL
        iptw2_bias_est iptw2_LL iptw2_UL
        iptw3_bias_est iptw3_LL iptw3_UL;
    label crude_bias_est = "Crude Bias"   crude_LL = "Crude LL"   crude_UL = "Crude UL"
          iptw1_bias_est = "IPTW1 Bias"   iptw1_LL = "IPTW1 LL"   iptw1_UL = "IPTW1 UL"
          iptw2_bias_est = "IPTW2 Bias"   iptw2_LL = "IPTW2 LL"   iptw2_UL = "IPTW2 UL"
          iptw3_bias_est = "IPTW3 Bias"   iptw3_LL = "IPTW3 LL"   iptw3_UL = "IPTW3 UL";
run;


/* absolute bias */
proc univariate data=bias3 noprint;
    var ab_crude_bias;
    output out=pctl_ab_crude mean=ab_crude_bias_est pctlpts=2.5 97.5 pctlpre=ab_crude pctlname=_LL _UL;
run;
proc univariate data=bias3 noprint;
    var ab_iptw1_bias;
    output out=pctl_ab_iptw1 mean=ab_iptw1_bias_est pctlpts=2.5 97.5 pctlpre=ab_iptw1 pctlname=_LL _UL;
run;
proc univariate data=bias3 noprint;
    var ab_iptw2_bias;
    output out=pctl_ab_iptw2 mean=ab_iptw2_bias_est pctlpts=2.5 97.5 pctlpre=ab_iptw2 pctlname=_LL _UL;
run;
proc univariate data=bias3 noprint;
    var ab_iptw3_bias;
    output out=pctl_ab_iptw3 mean=ab_iptw3_bias_est pctlpts=2.5 97.5 pctlpre=ab_iptw3 pctlname=_LL _UL;
run;

data percentiles_output_ab;
    merge pctl_ab_crude pctl_ab_iptw1 pctl_ab_iptw2 pctl_ab_iptw3;
run;

proc print data=percentiles_output_ab noobs label;
    var ab_crude_bias_est ab_crude_LL ab_crude_UL
        ab_iptw1_bias_est ab_iptw1_LL ab_iptw1_UL
        ab_iptw2_bias_est ab_iptw2_LL ab_iptw2_UL
        ab_iptw3_bias_est ab_iptw3_LL ab_iptw3_UL;
    label ab_crude_bias_est = "Crude |Bias|"   ab_crude_LL = "Crude |Bias| LL"   ab_crude_UL = "Crude |Bias| UL"
          ab_iptw1_bias_est = "IPTW1 |Bias|"   ab_iptw1_LL = "IPTW1 |Bias| LL"   ab_iptw1_UL = "IPTW1 |Bias| UL"
          ab_iptw2_bias_est = "IPTW2 |Bias|"   ab_iptw2_LL = "IPTW2 |Bias| LL"   ab_iptw2_UL = "IPTW2 |Bias| UL"
          ab_iptw3_bias_est = "IPTW3 |Bias|"   ab_iptw3_LL = "IPTW3 |Bias| LL"   ab_iptw3_UL = "IPTW3 |Bias| UL";
run;




/* MSE and 95% CI */


/* Step 1: Calculate the mean (msehat_xx) of the squared error terms (mse_xx) */
proc means data=bias3 noprint;
    var mse_crude mse_iptw1 mse_iptw2 mse_iptw3;
    output out=mse_point (drop=_type_ _freq_)
        mean = MSEhat_crude MSEhat_iptw1 MSEhat_iptw2 MSEhat_iptw3
        n    = nsim_crude nsim_iptw1 nsim_iptw2 nsim_iptw3;
run;

/* Step 2: add MSEhat to every replicate and calculate squared deviation between squared error term (mse_xx) and MSE */
data bias4;
    if _n_=1 then set mse_point;
    set bias3;

    dev_crude = (mse_crude - MSEhat_crude)**2;
    dev_iptw1 = (mse_iptw1 - MSEhat_iptw1)**2;
    dev_iptw2 = (mse_iptw2 - MSEhat_iptw2)**2;
    dev_iptw3 = (mse_iptw3 - MSEhat_iptw3)**2;
run;

/* Step 3: Sum the squared deviations across replicates (output sumdev_xx) */
proc means data=bias4 noprint;
    var dev_crude dev_iptw1 dev_iptw2 dev_iptw3;
    output out=mcse_sums (drop=_type_ _freq_)
        sum = sumdev_crude sumdev_iptw1 sumdev_iptw2 sumdev_iptw3;
run;

/* Step 4: point estimate, MCSE, and 95% CI (see Morris et al. 2019) */
data mse_final;
    merge mse_point mcse_sums;

    MCSE_MSE_crude = sqrt(sumdev_crude / (nsim_crude*(nsim_crude-1)));
    LCL_MSE_crude  = MSEhat_crude - 1.96*MCSE_MSE_crude;
    UCL_MSE_crude  = MSEhat_crude + 1.96*MCSE_MSE_crude;
    RMSE_crude     = sqrt(MSEhat_crude);

    MCSE_MSE_iptw1 = sqrt(sumdev_iptw1 / (nsim_iptw1*(nsim_iptw1-1)));
    LCL_MSE_iptw1  = MSEhat_iptw1 - 1.96*MCSE_MSE_iptw1;
    UCL_MSE_iptw1  = MSEhat_iptw1 + 1.96*MCSE_MSE_iptw1;
    RMSE_iptw1     = sqrt(MSEhat_iptw1);

    MCSE_MSE_iptw2 = sqrt(sumdev_iptw2 / (nsim_iptw2*(nsim_iptw2-1)));
    LCL_MSE_iptw2  = MSEhat_iptw2 - 1.96*MCSE_MSE_iptw2;
    UCL_MSE_iptw2  = MSEhat_iptw2 + 1.96*MCSE_MSE_iptw2;
    RMSE_iptw2     = sqrt(MSEhat_iptw2);

    MCSE_MSE_iptw3 = sqrt(sumdev_iptw3 / (nsim_iptw3*(nsim_iptw3-1)));
    LCL_MSE_iptw3  = MSEhat_iptw3 - 1.96*MCSE_MSE_iptw3;
    UCL_MSE_iptw3  = MSEhat_iptw3 + 1.96*MCSE_MSE_iptw3;
    RMSE_iptw3     = sqrt(MSEhat_iptw3);
run;

proc print data=mse_final noobs;
    var MSEhat_crude LCL_MSE_crude UCL_MSE_crude
        MSEhat_iptw1 LCL_MSE_iptw1 UCL_MSE_iptw1
        MSEhat_iptw2 LCL_MSE_iptw2 UCL_MSE_iptw2
        MSEhat_iptw3 LCL_MSE_iptw3 UCL_MSE_iptw3;
run;

proc print data=mse_final noobs label;
    var RMSE_crude RMSE_iptw1 RMSE_iptw2 RMSE_iptw3;
    label RMSE_crude = "Crude Square Root of MSE"
          RMSE_iptw1 = "IPTW1 Square Root of MSE"
          RMSE_iptw2 = "IPTW2 Square Root of MSE"
          RMSE_iptw3 = "IPTW3 Square Root of MSE";
run;

