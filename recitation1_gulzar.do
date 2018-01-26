** Saad Gulzar
** Quant 2, recitation 1
** Feb 1, 2017

/*
*/

clear all
set seed 20160201
set maxvar 32767
set matsize 10000
set more off

** CLT example from Cyrus
* x ~ Gamma(1,2)
	local npop = 10000
	local k = 1
	local theta = 2

* drawing underlying X
	set obs `npop'
	gen x_sc = rgamma(`k', `theta')
	local sigma_x = sqrt(`k' * (`theta'^2))
	
	gen mu_x_sc = `k' * `theta'
	gen x = x_sc - mu_x_sc
	
	  * see what the distributionl looks like
	    kdensity x, scheme(s1color) 
	
* drawing 5 observations at random from x
	gen rnd = runiform()
	sort rnd
	gen x5 = x if _n<=5
	drop rnd
	
* calculating sum and then standardizing it
	egen x5_total = total(x5)
	gen x5_stansum = x5_total / (`sigma_x' * sqrt(5))
	
* now do this 1,000 times and plot the standardized sum	 
	drop x5*
	mat sums = .
		
	forvalues j = 1/1000 {
		gen rnd = runiform()
		sort rnd
		gen x5 = x if _n<=5
		drop rnd

		egen x5_total = total(x5)
		local x5_stansum = x5_total / (`sigma_x' * sqrt(5))
		
		*storing results
			if `j' == 1 mat sums = `x5_stansum'
			else mat sums = sums \ `x5_stansum'
		
		drop x5*
		}
	
* plot the distribution of sums
svmat sums
kdensity sums1

* next, try this out with larger samples: 5, 100	
drop sums1

	foreach i in 5 100  {
	mat sums`i' = .
	
		forvalues j = 1/1000 {
			gen rnd = runiform()
			sort rnd
			gen x`i' = x if _n<=`i'
			drop rnd

			egen x`i'_total = total(x`i')
			local x`i'_stansum = x`i'_total / (`sigma_x' * sqrt(`i'))
		
			*storing results
			if `j' == 1 mat sums`i' = `x`i'_stansum'
			else mat sums`i' = sums`i' \ `x`i'_stansum'
						
			drop x`i'*
			}
		}
svmat sums5
svmat sums100

* plotting results
twoway (kdensity sums51, lwidth(thick)) ///
	(kdensity sums1001, lwidth(thick)) ///
	///
	, scheme(s1color) legend(label(1 "50 obs") label(2 "100 obs") label(3 "1000 obs")) xtitle(Z_n)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
