/*Karen Diepstra*/
/*Applied example: SMD estimates with 95% confidence intervals*/

/* Input dataset named ISO */
/* Input dataset not publicly available */

/* For unweighted SMD estimating equations for categorical variables, see: "SMD (95% CIs) CALCULATION -- CATEGORICAL VARS" */
/* For unweighted SMD estimating equations for continuous variables, see: "SMD (95% CIs) CALCULATION -- CONTINUOUS VARS" */
/* For weighted SMD estimating equations for categorical variables, see: "WEIGHTED SMD (95% CIs) CALCULATION -- CATEGORICAL VARS" */
/* For weighted SMD estimating equations for continuous variables, see: "WEIGHTED SMD (95% CIs) CALCULATION -- CONTINUOUS VARS" */
/* See Ross et al. PMCID: PMC10904145 for PROC IML code (other than SMD EEs) */


libname smd "xxxxx" ;

option nofmterr ;




*read ISO ;

data iso ;
 set smd.iso ;
run;




/*CREATE STANDARDIZED MEAN DIFFERENCES IN WEIGHTED AND UNWEIGHTED DATA FOR CATEGORICAL AND CONTINUOUS VARIABLES*/

/*DELETE PREVIOUS MEANS: AND WMEANS: DATASETS*/
proc datasets library=work ;
 delete means: wmeans: ;
quit;
run;


/* SPECIFY CATEGORICAL AND CONTINUOUS COVARIATES AND IPT WEIGHT */

%let catvars = hc artreg2 bf /*base_edu*/ primsomesec secondormore /*base_hivdx*/ hivonetofive hivfiveplus previso /*site*/ /*site1*/ site2 site3 site4 ;
%let convars = age /*visitnum agesq visitnumsq*/ ; 
%let weight = sipt_weight_trim ;


/* reformat artreg */
data iso ;
 set iso ;
 if artreg = 1 then artreg2 = 0 ;
 else if artreg = 2 then artreg2 = 1 ;
run;

/* restrict to complete cases */
data cc ;
 set iso ;
 where completecase = 1 ;
run;



/*************************************************************************************************************************
*************************CRUDE DATASET***********************************************************************************
**************************************************************************************************************************/

/* SMD (95% CIs) CALCULATION -- CATEGORICAL VARS */

%macro replicate(vars=,output=) ;

%local var i;                                            /*https://communities.sas.com/t5/SAS-Programming/How-can-I-loop-a-series-of-procedures-to-run-through-a-list-of/td-p/763397*/
 %do i=1 %to %sysfunc(countw(&vars));
  %let var=%scan(&vars,&i);

data rep_&var ;
 set cc ;
 where &var ne . ;
run;

PROC IML;                            
	*Read data;
	use rep_&var;								
		read all var {anyinhprev} into exposure;		
		read all var {&var} into z;
	close rep_&var;

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

  data results;
   set results;
   length covariate $ 20;  
   covariate = "&var";
run;


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

%mend;

%replicate(vars=&catvars,output=smdscat) ;




/* SMD (95% CIs) CALCULATION -- CONTINUOUS VARS */

%macro replicate(vars=,output=) ;

%local var i;                                            /*https://communities.sas.com/t5/SAS-Programming/How-can-I-loop-a-series-of-procedures-to-run-through-a-list-of/td-p/763397*/
 %do i=1 %to %sysfunc(countw(&vars));
  %let var=%scan(&vars,&i);

data rep_&var ;
 set cc ;
 where &var ne . ;
run;

PROC IML;                            
	*Read data;
	use rep_&var;								
		read all var {anyinhprev} into exposure;		
		read all var {&var} into z;
	close rep_&var;

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
	initial = {30,30,30,30,0.05};       * Initial parameter values; *value for each parameter;
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

%mend;

%replicate(vars=&convars,output=smdscon) ;


data crude ; 
 set smdscat
     smdscon ;
run;

data crude ;
 set crude ;
 where variable="SMD" ;
 
 rename est = crude ;
 rename lcl = crude_ll ;
 rename ucl = crude_ul ;

 drop variable variance se ;

run;








/*************************************************************************************************************************
*************************POST (WEIGHTED) DATASET***********************************************************************************
**************************************************************************************************************************/


/* WEIGHTED SMD (95% CIs) CALCULATION -- CATEGORICAL VARS */

%macro replicate(vars=,output=) ;

%local var i;                                            /*https://communities.sas.com/t5/SAS-Programming/How-can-I-loop-a-series-of-procedures-to-run-through-a-list-of/td-p/763397*/
 %do i=1 %to %sysfunc(countw(&vars));
  %let var=%scan(&vars,&i);

data rep_&var ;
 set cc ;
run;

PROC IML;                            
	*Read data;
	use rep_&var;								
		read all var {anyinhprev} into exposure;
        read all var {hc} into hc;	
		read all var {artreg2} into artreg;	
		read all var {bf} into bf;	
		read all var {secondormore} into edu2;	
		read all var {primsomesec} into edu1;
		read all var {hivfiveplus} into hiv2;
		read all var {hivonetofive} into hiv1;
		read all var {site2} into site2;
		read all var {site3} into site3;
		read all var {site4} into site4;
		read all var {previso} into previso;
		read all var {visitnum} into visitnum;
		read all var {age} into age;
		read all var {age1} into age1;
		read all var {age2} into age2;
		read all var {age3} into age3;
		read all var {age4} into age4;
		read all var {age5} into age5;
		read all var {&var} into z;
	close rep_&var;

	n = nrow(exposure);  

	
	/***********************************************
	Defining estimating equations */

	q = 23;								

	START efunc(theta) global(exposure, hc, artreg, bf, edu1, edu2, hiv1, hiv2, site2, site3, site4, previso, visitnum, age, age1, age2, age3, age4, age5, z, n);	


		pscore = 1 / (1 + exp(-(theta[1] + theta[2]*hc + theta[3]*artreg + theta[4]*bf + theta[5]*edu1 + theta[6]*edu2 + theta[7]*hiv1 + theta[8]*hiv2 + theta[9]*site2
                      + theta[10]*site3 + theta[11]*site4 + theta[12]*previso + theta[13]*visitnum + theta[14]*age + theta[15]*age1 + theta[16]*age2 + theta[17]*age3 
                      + theta[18]*age4 + theta[19]*age5)));	 /*Predicted propensity score */ 

		ef_1 = exposure - pscore;							 
		ef_2 = (exposure - pscore)#hc;						 
		ef_3 = (exposure - pscore)#artreg;						 
		ef_4 = (exposure - pscore)#bf;		
        ef_5 = (exposure - pscore)#edu1;		
        ef_6 = (exposure - pscore)#edu2;		
        ef_7 = (exposure - pscore)#hiv1;		
        ef_8 = (exposure - pscore)#hiv2;
		ef_9 = (exposure - pscore)#site2;	
		ef_10 = (exposure - pscore)#site3;	
		ef_11 = (exposure - pscore)#site4;
		ef_12 = (exposure - pscore)#previso;	
		ef_13 = (exposure - pscore)#visitnum;	
		ef_14 = (exposure - pscore)#age;	
		ef_15 = (exposure - pscore)#age1;
		ef_16 = (exposure - pscore)#age2;
		ef_17 = (exposure - pscore)#age3;
		ef_18 = (exposure - pscore)#age4;
		ef_19 = (exposure - pscore)#age5;
	

		ef_20 = exposure - theta[20];

	    wt = (exposure * theta[20] / pscore) + ((1 - exposure) * (1-theta[20]) / (1 - pscore));   /*Stabilized IPW*/

	    wt2 = wt;
		idx_low = loc(wt < 0.456275);
		idx_high = loc(wt > 5.7525184);

		if ncol(idx_low) > 0 then wt2[idx_low] = 0.456275;
		if ncol(idx_high) > 0 then wt2[idx_high] = 5.7525184;   /*Trimmed*/

		/* wt = exposure/pscore + (1-exposure)/(1-pscore); */  	 /*IPW */


        ef_21 = exposure#(z - theta[21])#wt2;
		ef_22 = (1-exposure)#(z - theta[22])#wt2;
		ef_23 = j(n, 1, (theta[21]-theta[22]) - (theta[23]*sqrt((.5*theta[21]*(1-theta[21])) + (.5*theta[22]*(1-theta[22]))))); 


		ef_mat = ef_1||ef_2||ef_3||ef_4||ef_5||ef_6||ef_7||ef_8||ef_9||ef_10||ef_11||ef_12||ef_13||ef_14||ef_15||ef_16||ef_17||ef_18||ef_19||ef_20||ef_21||ef_22||ef_23; 
		RETURN(ef_mat);                         						
	FINISH efunc;        


	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {2.3471054,-0.208242,-1.773456,-0.226322,-0.006354,-0.14653,0.0103589,0.0117231,-1.697518,-1.178631,-1.226237,2.3668119,-0.361178,-0.043765,0.2457653,0.2786596,-2.572694,7.6723784,-8.921592,.5,.5,.5,0};       * Initial parameter values; *value for each parameter;
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
	variable = {"int","hc","artreg","bf","edu1","edu2","hiv1","hiv2","site2","site3","site4","previso","visitnum","age","age1","age2","age3","age4","age5","meanexp","p1","p0","SMD"};  
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

  data results;
   set results;
   length covariate $ 20;  
   covariate = "&var";
run;


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

%mend;

%replicate(vars=&catvars,output=weightedsmdscat) ;



/* WEIGHTED SMD (95% CIs) CALCULATION -- CONTINUOUS VARS */


/* BOOTSTRAP TO GET 95% CONFIDENCE INTERVALS AROUND SMD */

*The bootstrap sandwich;

%macro theboot(seed=,data=,bsreps=,output=);

data seeds;
call streaminit(&seed);
	do j=0 to &bsreps;
		seed1 = ceil( (2**31 - 1)*rand("uniform") ); 
		seed2 = ceil( (2**31 - 1)*rand("uniform") ); 
		output;
	end;
run;


************************************;
* Loop through the BS 				;
************************************;
%do j = 0 %to &bsreps;

*Obtain seed;
proc sql noprint;
select seed1 format=16., seed2 format=16. into :seed1, :seed2 from seeds where j=&j;
quit;

	%if &j = 0 %then %do;
		*For rep 0 use all the data (point estimate);
		data rep0; set &data; run;
	%end;
	%else %do;
		*Obtain BS samples;
		proc surveyselect data = &data out = rep0 (drop = numberhits rename = (replicate = j) ) seed = &seed1 
			method = urs samprate = 1 outhits rep = 1 noprint; 
		run;
	%end;



************************************;
*CALL A MACRO WITH THE ANALYSIS PROGRAM HERE OR PASTE ANALYSIS PROGRAM. MUST USE THE REP DATASET FOR ANALYSIS;
*This analysis macro does every single step of the analysis: MI, g-comp, fitting p-score model, whatever your analysis entails, etc...;
************************************;


%macro replicate(vars=,output=) ;

%local var i;                                            /*https://communities.sas.com/t5/SAS-Programming/How-can-I-loop-a-series-of-procedures-to-run-through-a-list-of/td-p/763397*/
 %do i=1 %to %sysfunc(countw(&vars));
  %let var=%scan(&vars,&i);

data rep_&var ;
 set rep0 ; /*bootstrapped sample of cc*/
run;

PROC IML;                            
	*Read data;
	use rep_&var;								
		read all var {anyinhprev} into exposure;
        read all var {hc} into hc;	
		read all var {artreg2} into artreg;	
		read all var {bf} into bf;	
		read all var {secondormore} into edu2;	
		read all var {primsomesec} into edu1;
		read all var {hivfiveplus} into hiv2;
		read all var {hivonetofive} into hiv1;
		read all var {site2} into site2;
		read all var {site3} into site3;
		read all var {site4} into site4;
		read all var {previso} into previso;
		read all var {visitnum} into visitnum;
		read all var {age} into age;
		read all var {age1} into age1;
		read all var {age2} into age2;
		read all var {age3} into age3;
		read all var {age4} into age4;
		read all var {age5} into age5;
		read all var {&var} into z;
	close rep_&var;

	n = nrow(exposure);  

	
	/***********************************************
	Defining estimating equations */

		q = 22;								

	START efunc(theta) global(exposure, hc, artreg, bf, edu1, edu2, hiv1, hiv2, site2, site3, site4, previso, visitnum, age, age1, age2, age3, age4, age5, z, n);	


		pscore = 1 / (1 + exp(-(theta[1] + theta[2]*hc + theta[3]*artreg + theta[4]*bf + theta[5]*edu1 + theta[6]*edu2 + theta[7]*hiv1 + theta[8]*hiv2 + theta[9]*site2
                      + theta[10]*site3 + theta[11]*site4 + theta[12]*previso + theta[13]*visitnum + theta[14]*age + theta[15]*age1 + theta[16]*age2 + theta[17]*age3 
                      + theta[18]*age4 + theta[19]*age5)));	 /*Predicted propensity score */ 

		ef_1 = exposure - pscore;							 
		ef_2 = (exposure - pscore)#hc;						 
		ef_3 = (exposure - pscore)#artreg;						 
		ef_4 = (exposure - pscore)#bf;		
        ef_5 = (exposure - pscore)#edu1;		
        ef_6 = (exposure - pscore)#edu2;		
        ef_7 = (exposure - pscore)#hiv1;		
        ef_8 = (exposure - pscore)#hiv2;
		ef_9 = (exposure - pscore)#site2;	
		ef_10 = (exposure - pscore)#site3;	
		ef_11 = (exposure - pscore)#site4;
		ef_12 = (exposure - pscore)#previso;	
		ef_13 = (exposure - pscore)#visitnum;	
		ef_14 = (exposure - pscore)#age;	
		ef_15 = (exposure - pscore)#age1;
		ef_16 = (exposure - pscore)#age2;
		ef_17 = (exposure - pscore)#age3;
		ef_18 = (exposure - pscore)#age4;
		ef_19 = (exposure - pscore)#age5;
	

		ef_20 = exposure - theta[20];

	    wt = (exposure * theta[20] / pscore) + ((1 - exposure) * (1-theta[20]) / (1 - pscore));   /*Stabilized IPW*/

	    wt2 = wt;
		idx_low = loc(wt < 0.456275);
		idx_high = loc(wt > 5.7525184);

		if ncol(idx_low) > 0 then wt2[idx_low] = 0.456275;
		if ncol(idx_high) > 0 then wt2[idx_high] = 5.7525184;   /*Trimmed*/

		/* wt = exposure/pscore + (1-exposure)/(1-pscore); */  	 /*IPW */


        ef_21 = exposure#(z - theta[21])#wt2;
		ef_22 = (1-exposure)#(z - theta[22])#wt2;
		/*ef_23 = (exposure#((z - theta[21])##2 - theta[23]))#wt2;     */
		/*ef_24 = ((1-exposure)#((z - theta[22])##2 - theta[24]))#wt2; */
		/*ef_25 = j(n, 1, (theta[21]-theta[22]) - (theta[25]*(sqrt(.5*(theta[23] + theta[24])))));*/


		ef_mat = ef_1||ef_2||ef_3||ef_4||ef_5||ef_6||ef_7||ef_8||ef_9||ef_10||ef_11||ef_12||ef_13||ef_14||ef_15||ef_16||ef_17||ef_18||ef_19||ef_20||ef_21||ef_22; 
		RETURN(ef_mat);                         						
	FINISH efunc;        


	START eequat(theta);					 								
		ef = efunc(theta);
		RETURN(ef[+,]);                  								
	FINISH eequat;                       								

	/***********************************************
	Root-finding */
	initial = {2.3471054,-0.208242,-1.773456,-0.226322,-0.006354,-0.14653,0.0103589,0.0117231,-1.697518,-1.178631,-1.226237,2.3668119,-0.361178,-0.043765,0.2457653,0.2786596,-2.572694,7.6723784,-8.921592,.5,30,30};       * Initial parameter values; *value for each parameter;
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
	variable = {"int","hc","artreg","bf","edu1","edu2","hiv1","hiv2","site2","site3","site4","previso","visitnum","age","age1","age2","age3","age4","age5","meanexp","mean1","mean2"};  
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

data weightedsmdscon2 ;
 set weightedsmdscon ;
 where variable in ("mean1" "mean2") ;
 drop se lcl ucl ;
run;

/*long to wide*/
proc sort data=weightedsmdscon2;
    by covariate;
run;

proc transpose data=weightedsmdscon2 out=temp1(drop=_name_);
    by covariate;
    id variable;
    var est;
run;

proc transpose data=weightedsmdscon2 out=temp2(drop=_name_);
    by covariate;
    id variable;
    var variance;
run;

data weightedsmdscon3;
    merge temp1(rename=(mean1=est_mean1 mean2=est_mean2))
          temp2(rename=(mean1=variance_mean1 mean2=variance_mean2));
    by covariate;
    
    covariate_mean1 = covariate;
    covariate_mean2 = covariate;
run;

proc freq data=rep_age ;
 table anyinhprev / out=freq_counts (keep=anyinhprev count);
run;

proc transpose data=freq_counts out=temp3(drop=_name_ _label_);
    id anyinhprev;
    var count;
run;

data temp3 ;
 set temp3 ;
 covariate="&var" ;
run;

proc sort data=temp3 ; by covariate ; run;
proc sort data=weightedsmdscon3 ; by covariate ; run;

data weightedsmdscon3 ; 
 merge weightedsmdscon3 temp3 ;
 by covariate ;
run;

/*Calculate the SMD using the weighted means and their variances*/
data weightedsmdscon4 ;
 set weightedsmdscon3 ;

 var1 = variance_mean1 * _1 ; /*variance for the SD in the exposed*/
 var2 = variance_mean2 * _0 ; /*variance for the SD in the unexposed*/

 *SMD calculation;
 wSMD = (est_mean1 - est_mean2) / (sqrt((var1 + var2)/2)) ;

 bs=&j;   *bootstrap sample number;

 drop covariate est_mean1 est_mean2 variance_mean1 variance_mean2 covariate_mean1 covariate_mean2 ;

 run;

%end;

%mend;

%replicate(vars=&convars,output=weightedsmdscon) ;



	* CREATE OUTPUT DATASET - *These will be editted depending on what output from you analysis looks like!;
	%if %sysfunc(exist(&output)) %then %do; 
		data &output;
		set &output weightedsmdscon4;
		run;
		/*POSSIBLE ALTERNATIVE DEPENDING ON RESULTS STRUCTURE;
		data &output;
		merge &output RESULTSFROMANALYSISMACRO;
		by interv ; 
		run;*/
	%end;
	%else %do;
		data &output; set weightedsmdscon4; run;
	%end;

	*Delete datasets used in analysis macro here so memory is clear before next rep (do not delete seeds);
	proc datasets nolist; 
	delete weightedsmdscon weightedsmdscon2 weightedsmdscon3 weightedsmdscon4 ; 
	quit; run;
%end;

   proc datasets nolist; 
	delete seeds ; 
   quit; run;
%mend;


%theboot(seed=4448909,data=cc,bsreps=1000,output=weightedsmdscon5); 

data smd ;
 set weightedsmdscon5 ;
 where bs=0 ;

 covariate = "age" ;

 drop _0 _1 var1 var2 bs ;
run;

proc univariate data=weightedsmdscon5 noprint;
   var wSMD;
   output out=percentiles 
          pctlpts = 2.5 97.5 
          pctlpre = P
          std = wSMD_SD;
run;

data percentiles ;
 set percentiles ;
 covariate = "age" ;
run;

proc sort data=smd ; by covariate ; run;
proc sort data=percentiles ; by covariate ; run;

/*Create dataset with continuous SMD for age and 95% CI based on percentiles and 1.96*SD*/
/*variables: covariate, post, post_ll, post_ul, post_ll_sd, post_ul_sd*/
data weightedsmdscon6 ;
 merge smd percentiles ;
 by covariate ;
 
 rename wSMD = post ;
 rename p2_5 = post_ll ;
 rename p97_5 = post_ul ;

 post_ll_sd = wSMD - 1.96*wSMD_SD ; 
 post_ul_sd = wSMD + 1.96*wSMD_SD ;
run;

proc print data=weightedsmdscon6 ;
 var post_ll_sd post_ul_sd ;
run;


data weightedsmdscat2 ;
 set weightedsmdscat ;
 where variable="SMD" ;
 
 rename est = post ;
 rename lcl = post_ll ;
 rename ucl = post_ul ;

 drop variable variance se ;

run;


data post ;
 set weightedsmdscat2
     weightedsmdscon6 ;
run;







/*CONCATENATE CRUDE AND POST DATASETS*/
proc sort data=crude ; by covariate ; run;
proc sort data=post ; by covariate ; run;

data table ;
 merge crude post ;
 by covariate ;
run;

/*Rename covariates for Figure*/
data table ;
 length covariate $ 25;
 set table ;
 if covariate="age" then covariate="Age";
 if covariate="agesq" then covariate="Age squared";
 if covariate="artreg2" then covariate="Dolutegravir ART";
 if covariate="bf" then covariate="Breastfeeding";
 if covariate="hc" then covariate="DMPA";
 if covariate="hivfiveplus" then covariate="HIV dx 5+ years ago";
 if covariate="hivonetofive" then covariate="HIV dx 1-5 years ago";
 if covariate="previso" then covariate="Previous isoniazid";
 if covariate="primsomesec" then covariate="Education: some secondary";
 if covariate="secondormore" then covariate="Education: secondary+";
 if covariate="site1" then covariate="Site 1";
 if covariate="site2" then covariate="Site 2";
 if covariate="site3" then covariate="Site 3";
 if covariate="site4" then covariate="Site 4";
 if covariate="visitnum" then covariate="Visit number";
 if covariate="visitnumsq" then covariate="Visit number squared";

 /*scale all numbers by 100 for love plot*/
 crude = crude*100 ;
 post = post*100 ;
 crude_ll = crude_ll*100 ;
 crude_ul = crude_ul*100 ;
 post_ll = post_ll*100 ;
 post_ul = post_ul*100 ;

 crude_cld = crude_ul - crude_ll ;
 post_cld = post_ul - post_ll ;

run;

proc print data=table ; run;

data abstable;
 set table;
 smd=abs(post);
run;

proc means data=abstable n mean min max;
 var smd;
run;





/*SAVE DATASET IF MAKING LOVE PLOT*/

libname love 'xxxxx';


data love.table ;
 set table ;
run;




/*READ in TABLE */

data table ;
 set love.table ;
run;




/*************************************************************************************************************************
*************************CALL LOVE PLOT MACRO***********************************************************************************
**************************************************************************************************************************/

*run macro;
%include "xxxxx" ;


*call macro;
%love(working_directory = xxxx, 
			  sortchoice = original, 
			  method = w,
              threshold = 1); *my standardized mean differences are multiplied by 100, so threshold should be 0.1 x 100 = 10;
