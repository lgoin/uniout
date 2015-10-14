* Author: Lucia Goin
* Purpose: Export univariate regressions with models in rows instead of columns
* Date: 100215
****************************************************************************
program uniout
	version 11
	
	syntax anything(equalok) using/, [replace] [se RObust CLuster(varlist) avg(string) r p n lab dec(real 2) PNOte SENote]
	
	ParseMod `anything'
	
	local rhs `r(rhs)'
	local lhs `r(lhs)'

	* matrix dimensions
	loc nrows: word count `rhs'
	loc ncols: word count `lhs'

	* error messages for cluster/robust/se options (must specify 1 only)
	if "`robust'" != "" & "`cluster'" != "" {
		di as err "cannot specify both robust and cluster options"
		exit 198
		}
	
	* error message for dec()
	if `dec' != 1 & `dec' != 2 & `dec' != 3 & `dec' != 4 & `dec' != 5{
		di as err "dec() incorrectly specified. must be an integer 1-5"
		exit 198
		}

qui {
		* run regressions
	local i 0
		foreach var1 of local lhs {
			foreach var2 of local rhs {
				local ++i
				
				* run regressions depending on which SE option was specified
					if "`robust'" != "" {
						qui reg `var1' `var2', vce(robust)
					}
					else if "`cluster'" != "" {
						qui reg `var1' `var2', vce(cluster `cluster')
						}
					else {
						qui reg `var1' `var2'
					}
					
				global beta = _b[`var2']
				
				* the two globals below (sterr/pval) are used in the stars macro after this set
				global sterr = _se[`var2']
				loc t $beta/$sterr
				global pval = 2*ttail(e(df_r), abs(`t'))

				* these globals are calculated based on which options are selected
				if "`se'" != "" | "`cluster'" != "" | "`robust'" != ""  {
					global se = _se[`var2']
					}
				if "`r'" != "" {
					global rsquared = e(r2)
					}
				if "`n'" != "" {
					global count = e(N)
					}
				if "`p'" != "" {
					global p = 2*ttail(e(df_r), abs(`t'))
					}
				if "`avg'" != "" {
						if "`avg'" == "lhs" {
							local n_avg 1
							qui summ `var1' if e(sample) == 1
							global lhs_avg `r(mean)'
							}
						else if "`avg'" == "rhs" {
							loc n_avg 1
							loc rhs_avg rhs
							qui summ `var2' if e(sample) == 1
							global rhs_avg `r(mean)'
							}
						else {
							di as err "option avg() incorrectly specified. must be either avg(lhs) or avg(rhs)"
							exit 198
						}
					}
				else if "`avg'" == "" {
					local n_avg 0
				}
				
				* determine a count for the number of options specified; used in setting # matrix rows
				local n_se = 1 - mi("`se'")
				local n_r = 1 - mi("`r'")
				local n_n = 1 - mi("`n'")
				local n_p = 1 - mi("`p'")
				
				loc opts "se rsquared N p avg"
				loc selected "`n_se' `n_r' `n_n' `n_p' `n_avg'"
				
				loc nopts =`n_se'+`n_r'+`n_n'+`n_p'+`n_avg'+1
				
				local optsselected "$beta"
				local optsselected_nobeta ""
				
				* this is 1-5 because there are 5 optional additions printed in each matrix:
					* SE(cluster/robust/none), p-val, R2, N, and mean(lhs/rhs)
					* note that beta is always included and already listed in optsselected above
				forvalues j = 1/5 {
					local optscount1: word `j' of `opts'
					local optscount2: word `j' of `selected'
						
						if `optscount2' == 1 {
						
							if "`optscount1'" == "se" | "`optscount1'" == "cluster" |"`optscount1'" == "robust" {
								loc optsselected "`optsselected' ,$se"
								loc optsselected_nobeta "`optsselected_nobeta' se"
								}
							else if "`optscount1'" == "rsquared" {
								loc optsselected "`optsselected' ,$rsquared"
								loc optsselected_nobeta "`optsselected_nobeta' r2"
								}
							else if "`optscount1'" == "N" {
								loc optsselected "`optsselected' ,$count"
								loc optsselected_nobeta "`optsselected_nobeta' N"
								}
							else if "`optscount1'" == "p" {
								loc optsselected "`optsselected' ,$p"
								loc optsselected_nobeta "`optsselected_nobeta' p-value"
								}
							else if "`optscount1'" == "avg" {
									if "`avg'" == "lhs" {
										loc optsselected "`optsselected' ,$lhs_avg"
										loc optsselected_nobeta "`optsselected_nobeta' mean(lhs)"
										}
									else if "`avg'" == "rhs" {
										loc optsselected "`optsselected' ,$rhs_avg"
										loc optsselected_nobeta "`optsselected_nobeta' mean(rhs)"
										}
							}
						}
				}
				
				* put regression results into a matrix, transpose
				matrix A`i' = (`optsselected')
				matrix A`i' = A`i''
				
				* To change the significance levels that the stars represent:
					if $pval <0.1 & $pval >= 0.05 {
						loc stars "`stars' one"
					}
					else if $pval <0.05 & $pval >= 0.01 {
						loc stars "`stars' two"
					}
					else if $pval < 0.01 {
						loc stars "`stars' three"
					}
					else {
						loc stars "`stars' none"
					}
			}
			
	}

		* define matrix dimensions. rows = number of RHS vars * number of options selected
		local numrows =`nrows' * `nopts'
		matrix H = J(`numrows',`ncols',.)
		local max = `nrows' * `ncols'
		numlist "1/`numrows'", integer ascending
		local range `r(numlist)'

		* this loop pulls the top-left hand corner of each submatrix that the loop below uses
		local iter 1
		forvalues s = 1/`nrows' {
			loc tmprows: word `iter' of `range'
			loc rows `rows' `tmprows'
			loc iter = `iter' + `nopts'
		}


		* this loop stacks the matrices together
		local iter 0
		forvalues j = 1/`ncols' {
			foreach row of local rows {
				local ++iter
				matrix H[`row', `j'] = A`iter'
				}
			}
			
		loc rownames ""
		* create local of rownames for matrix H
		if "`lab'" == ""{
			foreach name of local rhs {
				loc a "`name' `optsselected_nobeta'"
				loc rownames `rownames' `a'
			}
		}
		* this loop creates rownames if option lab specified
		else if "`lab'" != "" {
			foreach name of local rhs {
				loc lab1: var label `name'
				loc lab1 = trim("`lab1'")
				local lab1 = subinstr("`lab1'",".","",.)
				loc lab1 "`lab1'"
				loc rhs_labs `" `rhs_labs' "`lab1'" "'
				loc a `" "`lab1'" `optsselected_nobeta' "'
				loc rownames  "`rownames' `a'" 
			}
		}

		* set row and column names
		matrix rownames H =`rownames'
		matrix colnames H = `lhs'
		
		* create local w var labels to save 
		if "`lab'" != "" {
			foreach name of local lhs {
				loc lab1: var label `name'
				loc lab1 = trim("`lab1'")
				local lab1 = subinstr("`lab1'",".","",.)
				loc lab1 "`lab1'"
				loc col_labs `" `col_labs' "`lab1'" "' 
			}
		}	
		
		* rename variable names using first row
		noi estout matrix(H)
		esttab matrix(H) using regressions.csv, replace 
		
		preserve
		clear all
		insheet using regressions.csv
		drop in 1
		replace v1 = "models" in 1

		foreach var of varlist _all {
		  local newname  = `var'[1]
		  ren `var' `newname'
		}

		drop in 1
		destring, replace

		* string formatting
		ds, has (type numeric)
		loc vars `r(varlist)'

		* create local sigdigs for decimal place rounding
		if `dec' == 1 {
			loc sigdig 0.1
		}
		else if `dec' == 2 {
			loc sigdig 0.01
		}
		else if `dec' == 3 {
			loc sigdig 0.001
		}
		else if `dec' == 4 {
			loc sigdig 0.0001
		}
		else if `dec' == 5 {
			loc sigdig 0.00001
		}
		else if `dec' > 5 | `dec' < 1 {
            di as err "dec() must be an integer 1-5"
			exit 198
		}

		
	local i 0
		if "`lab'" == "" {
			foreach name of local vars {
				replace `name' = round(`name', `sigdig')
				tostring `name', replace force
				replace `name' = "(" + `name' + ")" if models == "se"
					foreach var of local rhs {
						local ++i
						local count: word `i' of `stars'
							if "`count'" == "one" {
								replace `name' = `name' + "*" if models == "`var'"
							}
							else if "`count'" == "two" {
								replace `name' = `name' + "**" if models == "`var'"
							}
							else if "`count'" == "three" {
								replace `name' = `name' + "***" if models == "`var'"
							}
					}
			}
		}
		else if "`lab'" != "" {
			foreach name of local vars {
				replace `name' = round(`name', `sigdig')
				tostring `name', replace force
				replace `name' = "(" + `name' + ")" if models == "se"
					foreach var of local rhs_labs {
						local ++i
						local count: word `i' of `stars'
							if "`count'" == "one" {
								replace `name' = `name' + "*" if models == "`var'"
							}
							else if "`count'" == "two" {
								replace `name' = `name' + "**" if models == "`var'"
							}
							else if "`count'" == "three" {
								replace `name' = `name' + "***" if models == "`var'"
							}
					}
			}
			expand 2 in 1
			egen list = seq()
			replace list = 0 if _N == list
			sort list
		
				foreach var of varlist _all {
					cap replace `var' ="" in 1
					cap replace `var' =. in 1
				}
		
				local n: word count `lhs'
				forvalues j = 1/`n' {
							local c: word `j' of `lhs'
							local d: word `j' of `col_labs'
							replace `c' = "`d'" in 1
						}
			
			
		}
		
		* include note at end of table
			if "`pnote'" !="" {
				loc newobs = _N+1
				set obs `newobs'
				replace models = "Note: * p < 0.1, ** p < 0.05, *** p < 0.01" in `newobs'
			}
			
			if "`senote'" != "" {
				if "`cluster'" != "" {
					loc newobs = _N+1
					set obs `newobs'
					replace models = "Note: Standard errors clustered by `cluster'" in `newobs'
				}
				else if "`robust'" != "" {
					loc newobs = _N+1
					set obs `newobs'
					replace models = "Note: Robust standard errors reported" in `newobs'
				}
			
			}
	
		
		rm regressions.csv
	}
	
	cap drop list
	xmlsave * using `using', doctype(excel) replace 
	mac drop pval lhs_avg rhs_avg count se rsquared sterr
	restore
end

* this program parses the syntax for the rhs and lhs variables
program ParseMod, rclass

        gettoken yvar 0 : 0 , parse(" /:=")  // no need for comma in list here
        if inlist(`"`yvar'"',":","=") {
                di as err "at least one dependent variable required"
                exit 198
        }
        if `"`yvar'"' == "/" {
                di as err "invalid use of /"
                exit 198
        }

        while !inlist(`"`yvar'"',":","=") {
                if `"`yvar'"' == "" || `"`yvar'"' == "," {
                        di as err ///
        "equal sign required between dependent variables and terms of the model"
                        exit 198
                }
                if `"`yvar'"' == "/" {
                        di as err ///
                      "invalid use of /; not allowed among dependent variables"
                        exit 198
                }
                local yvars `yvars' `yvar'
                
                gettoken yvar 0 : 0 , parse(" /:=,")
        }

        if `"`0'"' == "" { // no right-hand-side vars
                di as err "at least one independent variable required"
                exit 198
        }

        _anovaparse `0'

        return local rhs "`0'"
        return local lhs "`yvars'"
end
