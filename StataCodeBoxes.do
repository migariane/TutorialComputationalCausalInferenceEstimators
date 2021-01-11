/* 
Tutorial: causal inference methods made easy for applied resarchers/epidemiologists/statisticians 
=================================================================================================

ICON-LSHTM, LONDON, 12th January 2021

Miguel Angel Luque Fernandez, PhD
Assistant Professor of Epidemiology and Biostatistics
Camille Maringe, PhD
Assistant Professor

Inequalities in Cancer Outcomes Network, LSHTM, London, UK

Copyright (c) 2021 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Bug reports: miguel-angel.luque@lshtm.ac.uk	

The rhc dataset can be dowloaded at http://biostat.mc.vanderbilt.edu/wiki/Main/DataSets
*/
 
*** Preliminaries
			clear
			set more off
			cd "your path\Data" 					// this path should point to where the RHC data are
			use "rhc.dta", clear
			describe
			count
			* 83 variables and 5,735 observations
			
/* Box 1: Setting the data */
			* Define the outcome (Y), exposure (A), confounder (C), and confounders (W)
			global Y death_d30 
			global A rhc 
			global C gender  
			global W gender age edu race carcinoma 
	
/* Box 2: Naive estimate of the ATE */
			* Naive approach to estimate the causal effect
			regr $Y $A $C 	
			* The naive estimate of the causal effect is 0.07352

/* 3. G-formula */
/* 3.1 Non-parametric G-formula */
       
/* Box 3: Non-parametric G-Formula for the ATE */
            proportion $C
            matrix m=e(b)
            gen genderf = m[1,1]
            sum genderf
            gen genderm = m[1,2]
            sum genderm
			* you may need to install the command sumup, type:
			* ssc install sumup
            sumup $Y, by($A $C)
			* from sumup command extract the conditinal means by the given A and C levels i.e. zero and one
			* see matrix list y00: position subscript [3,1] is th one of interest
            matrix y00 = r(Stat1)
            matrix y01 = r(Stat2)
            matrix y10 = r(Stat3)
            matrix y11 = r(Stat4)
            gen EY1 = ((y11[3,1]-y01[3,1]))*genderm 
			gen EY0 = ((y10[3,1]-y00[3,1]))*genderf
            qui: mean EY1 EY0
			matrix ATE =  r(table) 
            display "The ATE is: "  ATE[1,1] + ATE[1,2]
            drop EY*
            * The ATE from non-parametric estimator is:  0.073692
			// Also one can try
			gen ATE = ((y11[3,1]-y01[3,1]))*genderm + ((y10[3,1]-y00[3,1]))*genderf
			qui sum ATE
			drop ATE
			
            * Check that Stata "teffects" command obtains the same estimate
            teffects ra ($Y $C) ($A)
            * The ATE from "teffects" implementation is:  0.073692
    	
/* Box 4: Bootstrap 95% Confidence Intervals (CI) for the ATE/ATT estimated using the Non-parametric G-Formula */
        
		* 1) For the ATE
			capture program drop ATE
            program define ATE, rclass
				capture drop y1
				capture drop y0
				capture drop ATE
				sumup $Y, by($A $C)
				matrix y00 = r(Stat1)
				matrix y01 = r(Stat2)
				matrix y10 = r(Stat3)
				matrix y11 = r(Stat4)
				gen ATE = ((y11[3,1]-y01[3,1]))*genderm + ((y10[3,1]-y00[3,1]))*genderf
				qui sum ATE
				return scalar ate = `r(mean)'
            end
            
            qui bootstrap r(ate), reps(1000): ATE
            estat boot, all	
    	
/* Box 5: Non-parametric G-Formula using a fully saturated regression model in Stata (A) */
			* method 1: conditional probabilities
            regress $Y ibn.$A ibn.$A#c.($C) , noconstant vce(robust) coeflegend
            predictnl ATE = (_b[1.rhc] + _b[1.rhc#c.gender]*gender) - (_b[0bn.rhc] + _b[0bn.rhc#c.gender]*gender)
			qui: sum ATE
            display "The ATE is:  " "`r(mean)'"
            drop ATE
			
/* Box 6: Non-parametric G-Formula using a fully saturated regression model in Stata (B) */
			* method 2: marginal probabilities
            regress $Y ibn.$A ibn.$A#c.($C) , noconstant vce(robust) coeflegend
			
			* Marginal probability in each treatment group
			margins $A , vce(unconditional) 
			
			* Difference in marginal probability between treatment groups
            margins r.$A , contrast(nowald) 
			
/* 3.2 PARAMETRIC G-FORMULA */			

* One confounder

/* Box 7: Parametric G-formula */		
			* Calculations by hand
			* Expected probability amongst treated
			regress $Y $C if $A==1      
			predict double y1hat
			
			* Expected probability amongst untreated
			regress $Y $C if $A==0      
			predict double y0hat
			mean y1hat y0hat            
			
			* Difference between expected probabilities (ATE) and biased confidence interval
			lincom _b[y1hat] - _b[y0hat]  
	
/* Box 8: Parametric regression adjustment using Stata's teffects (one confounder) */
			teffects ra ($Y $C) ($A) 

/* Box 9: Bootstrap for the parametric regression adjustment */
			capture program drop ATE
			program define ATE, rclass
			   capture drop y1
			   capture drop y0
			   reg $Y $C if $A==1 
			   predict double y1, xb
			   quiet sum y1
			   reg $Y $C if $A==0 
			   predict double y0, xb 
			   quiet sum y0
			   mean y1 y0 
			   lincom _b[y1]-_b[y0]
			   return scalar ace =`r(estimate)'
			end
			qui bootstrap r(ace), reps(1000): ATE
			estat boot, all

* More than one confounder

/* Box 10: Parametric multivariate regression adjustment implementation of the G-Formula */
			regress $Y $W if $A==1
			predict double y1hat 
			regress $Y $W if $A==0
			predict double y0hat
			mean y1hat y0hat
			lincom _b[y1hat] - _b[y0hat]
		
/* Box 11: Parametric multivariate regression adjustment using Stata’s teffects command */	
			teffects ra ($Y $W) ($A)

/* Box 12: Parametric multivariate regression adjustment using Stata’s margins command */
			regress $Y ibn.$A ibn.$A#c.($W) , noconstant vce(robust)
			margins $A, vce(unconditional)
			margins r.$A, contrast(nowald)
		
/* Box 13: Bootstrap for the multivariate parametric regression adjustment */
			capture program drop ATE
			program define ATE, rclass
			   capture drop y1
			   capture drop y0
			   reg $Y $W if $A==1 
			   predict double y1, xb
			   quiet sum y1
			   reg $Y $W if $A==0 
			   predict double y0, xb 
			   quiet sum y0
			   mean y1 y0 
			   lincom _b[y1]-_b[y0]
			   return scalar ace =`r(estimate)'
			end
			qui bootstrap r(ace), reps(1000): ATE dots
			estat boot, all

/* 4 Inverse probability of treatment weighting  */		
/* 4.1 Inverse probability of treatment weighting based on the propensity score plus regression adjustment */

/* Box 14: Computation of the IPTW estimator for the ATE */
			* propensity score model for the exposure
			logit $A $W, vce(robust) nolog
			
			* propensity score predictions
			predict double ps
			
			*  Sampling  weights  for  the  treated  group
			generate double ipw1 = ($A==1)/ps
			
			*  Weighted  outcome  probability  among  treated
			regress $Y [pw=ipw1]
			scalar Y1 = _b[_cons]
			
			*  Sampling  weights  for  the  non-treated  group
			generate double ipw0 = ($A==0)/(1-ps)
			regress $Y [pw=ipw0]
			scalar Y0 = _b[_cons]
			display "ATE =" Y1 - Y0
	
/* Box 15: Bootstrap computation for the IPTW estimator */
			* Bootstrap the confidence intervals
			capture program drop ATE
			program define ATE, rclass
			   capture drop y1
			   capture drop y0
			   regress $Y [pw=ipw1]
			   matrix y1 = e(b)
			   gen double y1 = y1[1,1]
			   regress $Y [pw=ipw0]
			   matrix y0 = e(b)
			   gen double y0 = y0[1,1]
			   mean y1 y0
			   lincom _b[y1]-_b[y0]
			   return scalar ace = `r(estimate)'
			end
			qui bootstrap r(ace), reps(1000): ATE
			estat boot, all
    
/* Box 16: Computation of the IPTW estimator for the ATE using Stata’s teffects command */ 
			teffects ipw ($Y) ($A $W, logit), nolog vsquish

/* Box 17: Assessing IPTW balance */
			* Stata teffects and tebalance commands
			qui teffects ipw ($Y) ($A $W)   
			tebalance summarize
			
			* By hand - with the example of gender
			egen  genderst = std(gender) 	//  Standardization
			logistic $A $W 				 	//  Propensity  score
			capture drop ps
			predict  double  ps
			gen  ipw = .
			replace  ipw=($A==1)/ps if $A==1
			replace  ipw=($A==0)/(1-ps) if $A==0
			regress  genderst  $A 		 	// Raw  difference
			regress  genderst  $A [pw=ipw]  // Standardized  difference
		
/* Box 18: Assessing IPTW overlap by hand */
			sort $A
			by $A: summarize ps
			kdensity ps if $A==1, generate(x1pointsa d1A) nograph n(10000) 
			kdensity ps if $A==0, generate(x0pointsa d0A) nograph n(10000) 
			label variable d1A "density for RHC=1"
			label variable d0A "density for RHC=0"
			twoway (line d0A x0pointsa , yaxis(1))(line d1A x1pointsa, yaxis(2))
		
/* Box 19: Assessing overlap using Stata's teffects overlap */
			qui: teffects ipw ($Y) ($A $W, logit), nolog vsquish
			teffects overlap

			
/* 4.2 Marginal structural model with stabilized weights */		
/* Box 20: Computation of the IPTW estimator for the ATE using a MSM */
			* Baseline treatment probabilities
			logit $A, vce(robust) nolog
			predict double nps, pr
			
			* propensity score model
			logit $A $W, vce(robust) nolog
			predict double dps, pr
			
			* Unstabilized weight
			cap drop ipw
			gen ipw = .
			replace ipw=($A==1)/dps if $A==1
			replace ipw=($A==0)/(1-dps) if $A==0
			sum ipw 
		
			* Stabilized weight 
			gen sws = .
			replace sws = nps/dps if $A==1
			replace sws = (1-nps)/(1-dps) if $A==0
			sum sws

			* MSM 
			reg $Y $A [pw=ipw], vce(robust) 	// MSM unstabilized weight
			reg $Y $A [pw=sws], vce(robust) 	// MSM stabilized weight


/* 4.3 IPTW with regression adjustment */
		
/* Box 21: Computation of the IPTW-RA estimator for the ATE and bootstrap for statistical inference */
			capture program drop ATE
			program define ATE, rclass
			   capture drop y1
			   capture drop y0
			   reg $Y $W if $A==1 [pw=sws]
			   predict double y1, xb
			   quiet sum y1
			   return scalar y1=`r(mean)'
			   reg $Y $W if $A==0 [pw=sws]
			   predict double y0, xb 
			   quiet sum y0
			   return scalar y0=`r(mean)'
			   mean y1 y0 
			   lincom _b[y1]-_b[y0]
			   return scalar ace =`r(estimate)'
			end
			qui bootstrap r(ace), reps(10): ATE
			estat boot, all
		
/* Box 22: Computation of the IPTW-RA estimator for the ATE using Stata’s teffects */
			teffects ipwra ($Y $W) ($A $W), nolog vsquish
    
/* 5. Augmented inverse probability weighting */	
/* Box 23: Computation of the AIPTW estimator for the ATE and bootstrap for statistical inference */
			* Step (i) prediction model for the outcome
			qui glm $Y $A $W, fam(bin) 
			predict double QAW, mu
			qui glm $Y $W if $A==1, fam(bin) 
			predict double Q1W, mu
			qui glm $Y $W if $A==0, fam(bin)
			predict double Q0W, mu
			
			* Step (ii): prediction model for the treatment
			cap drop dps nps sws y1 y0
			qui logit $A $W
			predict double dps, pr
			qui logit $A
			predict double nps, pr
			gen sws = .
			replace sws = nps/dps if $A==1
			replace sws = (1-nps)/(1-dps) if $A==0
		
			* Step (iii): Estimation equation
			gen double y1 = (sws*($Y-QAW) + (Q1W))
			quiet sum y1
			scalar y1=`r(mean)'
			gen double y0 = (sws*($Y-QAW) + (Q0W))
			quiet sum y0
			scalar y0=`r(mean)'
			mean y1 y0
			lincom _b[y1] - _b[y0]
		
			* Step (iv): Bootstrap confidence intervals
			capture program drop ATE
			program define ATE, rclass
			capture drop y1
			capture drop y0
			capture drop Q*
			qui glm $Y $A $W, fam(bin) 
			predict double QAW, mu
			qui glm $Y $W if $A==1, fam(bin) 
			predict double Q1W, mu
			qui glm $Y $W if $A==0, fam(bin)
			predict double Q0W, mu
			gen double y1 = (sws*($Y-QAW) + (Q1W))
			quiet sum y1
			return scalar y1=`r(mean)'
			gen double y0 = (sws*($Y-QAW) + (Q0W))
			quiet sum y0
			return scalar y0=`r(mean)'
			mean y1 y0
			lincom _b[y1] - _b[y0]
			return scalar ace =`r(estimate)'
			end
			qui bootstrap r(ace), reps(1000): ATE
			estat boot, all

/* Box 24: Computation of the AIPTW estimator for the ATE using Stata’s teffects */
			teffects aipw ($Y $W) ($A $W, logit), nolog vsquish

		
/* 6. DATA-ADAPTIVE ESTIMATION: ENSEMBLE LEARNING TARGETED MAXIMUMLIKELIHOOD ESTIMATION*/
/*Box 25: Computational implementation of TMLE by hand */

			* Step 1: prediction model for the outcome Q0 (g-computation)
			glm $Y $A $W, fam(binomial)
			predict double QAW_0, mu
			gen aa=$A
			replace $A = 0
			predict double Q0W_0, mu
			replace $A= 1
			predict double Q1W_0, mu
			replace $A = aa
			drop aa

			// Q to logit scale
			gen logQAW = log(QAW / (1 - QAW))
			gen logQ1W = log(Q1W / (1 - Q1W))
			gen logQ0W = log(Q0W / (1 - Q0W))

			* Step 2: prediction model for the treatment g0 (IPTW)
			glm $A $W, fam(binomial)
			predict gw, mu
			gen double H1W = $A  / gw
			gen double H0W = (1 - $A ) / (1 - gw)

			* Step 3: Computing the clever covariate H(A,W) and estimating the parameter (epsilon) (MLE)
			glm $Y H1W H0W, fam(binomial) offset(logQAW) noconstant
			mat a = e(b)
			gen eps1 = a[1,1]
			gen eps2 = a[1,2]

			* Step 4: update from Q0 to Q1
			gen double Q1W_1 = exp(eps1 / gw + logQ1W) / (1 + exp(eps1 / gw + logQ1W))
			gen double Q0W_1 = exp(eps2 / (1 - gw) + logQ0W) / (1 + exp(eps2 / (1 - gw) + logQ0W))

			* Step 5: Targeted estimate of the ATE 
			gen ATE = (Q1W_1 - Q0W_1)
			summ ATE
			global ATE = r(mean)
			drop ATE 
			
			* Step 6: Statistical inference (efficient influence curve)
			qui sum(Q1W_1)
			gen EY1tmle = r(mean)
			qui sum(Q0W_1)
			gen EY0tmle = r(mean)

			gen d1 = (($A  * ($Y - Q1W_1)/gw)) + Q1W_1 - EY1tmle
			gen d0 = ((1 - $A ) * ($Y - Q0W_1)/(1 - gw))  + Q0W_1 - EY0tmle

			gen IC = d1 - d0
			qui sum IC
			gen varIC = r(Var) / r(N)
			drop d1 d0 IC
			
			global LCI =  $ATE - 1.96*sqrt(varIC)
			global UCI =  $ATE + 1.96*sqrt(varIC)
			display "ATE:"  %05.4f  $ATE _col(15) "95%CI: " %05.4f  $LCI "," %05.4f  $UCI		

/* Box 26: TMLE with data-adaptive estimation using the Stata’s user writen eltmle */
			* if not already installed, type:
			* github install migariane/eltmle (last update January 2021)
			* alternatively: ssc install eltmle
			preserve
			eltmle $Y $A $W, tmle   
			restore

		
/* 7. Simulation */
/* Box 27: Data generation for the Monte Carlo experiment */

			* Data generation
			clear
			set obs 1000
			set seed 777
			gen w1    = round(runiform(1, 5)) //Quintiles of Socioeconomic Deprivation
			gen w2    = rbinomial(1, 0.45) //Binary: probability age >65 = 0.45
			gen w3    = round(runiform(0, 1) + 0.75*(w2) + 0.8*(w1)) //Stage 
			recode w3 (5/6=1)       //Stage (TNM): categorical 4 levels
			gen w4    = round(runiform(0, 1) + 0.75*(w2) + 0.2*(w1)) //Comorbidites: categorical four levels
			gen A     = (rbinomial(1,invlogit(-1 -  0.15*(w4) + 1.5*(w2) + 0.75*(w3) + 0.25*(w1) + 0.8*(w2)*(w4)))) //Binary treatment
			gen Y1    = (invlogit(-3 + 1 + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1))) // Potential outcome 1
			gen Y0    = (invlogit(-3 + 0 + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1))) // Potential outcome 2
			gen psi   = Y1-Y0  // Simulated ATE
			gen Y     = A*(Y1) + (1 - A)*Y0 //Binary outcome
	
	// Estimate the true simulated ATE
			mean psi
				
		//	ATE estimation
			* Regression adjustment
			teffects ra (Y w1 w2 w3 w4) (A)
			estimates store ra
				
			* IPTW
			teffects ipw (Y) (A w1 w2 w3 w4)
			estimates store ipw
				
			* IPTW-RA
			teffects ipwra (Y w1 w2 w3 w4) (A w1 w2 w3 w4)
			estimates store ipwra	
				
			* AIPTW
			teffects aipw (Y w1 w2 w3 w4) (A w1 w2 w3 w4)
			estimates store aipw
				
			* Results
			qui reg psi
			estimates store psi
			estout psi ra ipw ipwra aipw
		
		// Ensemble learning maximum likelihood estimation
			preserve
			eltmle Y A w1 w2 w3 w4, tmle
			restore
				
		// Relative bias of each ATE
			* Regression adjustment
			display abs(0.1787 - 0.203419)/0.1787 
			
			* IPTW
			display abs(0.1787 - 0.2776)/0.1787 
			
			* IPTW-RA
			display abs(0.1787 - .2052088)/0.1787
			
			* AIPTW
			display abs(0.1787 - 0.2030)/0.1787
			
			* ELTMLE
			display abs(0.1787 - 0.1784)/0.1787
				
