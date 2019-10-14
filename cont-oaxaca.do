

/*******************************************************************************	
					
- Continuous Oaxaca Decomposition
	
*******************************************************************************/


** HH level dataset
********************************************* 
clear
use "$d3/HH_Merged_Ind.dta"

encode region, gen(nregion)
*drop if transp_time == 0

* -----------------------------------------------------------------------------
** Ulrick (2012) Continuous Oaxaca Decomposition
* -----------------------------------------------------------------------------

foreach z in 99 90 { // high percentile ranks
	foreach j in 50 75 { // low percentile ranks

		* Store variables in locals
		local Y `var'									// Outome variable
		local G `var'									// Group variable
		local X `vector'								// Control variables
		tokenize `X'
	
		* Percentiles
		local H p`z'											// High evaluation point
		local L p`j'											// Low evaluation point


		************************************************************************
		quietly { 	// Continuous Oaxaca Decomposition
	
			** Stage 1 
			* ------------------------------------------------------
			// Estimate E(X|i)
			** E(X|i_h)
			sum `G', d
			scalar _h = r(`H')

			forvalues i = 1/2 {
				reg ``i'' `G' 
				scalar Xhat`i'_h = _b[_cons] + _b[`G']*(_h)
				display Xhat`i'_h
			} 
			** E(X|i_l)
			sum `G', d
			scalar _l = r(`L')

			forvalues i = 1/2 {
				reg ``i'' `G' 
				scalar Xhat`i'_l = _b[_cons] + _b[`G']*(_l)
				display Xhat`i'_l
			} 
			scalar dx1 = Xhat1_h - Xhat1_l  
			scalar dx2 = Xhat2_h - Xhat2_l
	
		
			** Stage 2
			* ------------------------------------------------------
			// OLS regression
			reg `Y' `G' `X'
			return list
			matrix V = e(V)
			
			** Stage 3
			* ------------------------------------------------------
			// Generate estimates
			// Use delta method to generate SE's
			* Gap
			lincom _b[`1']*dx1 + _b[`2']*dx2 + _b[`3']*dx1*(_l) + _b[`4']*dx2*(_l) ///
					+ _b[`G']*(_h - _l) + _b[`3']*(_h - _l)*Xhat1_h + _b[`4']*(_h - _l)*Xhat2_h	
			scalar gap_`z'_`j' = r(estimate)
			scalar gap_se_`z'_`j' = r(se)
			scalar gap_p_`z'_`j' = r(p)
			* Explained
			lincom _b[`1']*dx1 + _b[`2']*dx2 + _b[`3']*dx1*(_l) + _b[`4']*dx2*(_l)	
			scalar ex_`z'_`j' = r(estimate)
			scalar ex_se_`z'_`j' = r(se)
			scalar ex_p_`z'_`j' = r(p)
			* Unexplained
			lincom _b[`G']*(_h - _l) + _b[`3']*(_h - _l)*Xhat1_h + _b[`4']*(_h - _l)*Xhat2_h	
			scalar un_`z'_`j' = r(estimate)
			scalar un_se_`z'_`j' = r(se)
			scalar un_p_`z'_`j' = r(p)
			
	
			** Store Results
			* ------------------------------------------------------
			*matrix Results = (gap, ex, un \ gap_se, ex_se, un_se \ gap_p, ex_p, un_p )
			

		************************************************************************
	
		* Results table
		frmttable using "$d2/results.doc", statmat(Results) sdec(3) ///
		title("Y ~ G : `z' - `j'") ///
		ctitle("","Gap","Explained","Unexplained") ///
		rtitle("Est."\"se"\"P-value") addtable 

		}
	}
}



