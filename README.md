# uniout
Uniout is a stata program that exports univariate regressions in a matrix to an .xml file

#### Syntax
`uniout` *depvars* = *indepvars* `using` *filename* [, *options*]

If no options are specified, uniout will display the beta from each regression and stars for the significance of each coefficient.

##### options

    n                   prints the Ns

    p                   prints the p-values

    r                   prints the R-squared

    ncol                Adds column numbers under the variable names

    dec()               specifies the number of decimal places for N, p-value, R-squared, beta, standard errors and avg()

    lab                 prints the variable labels instead of the names in the rows, and under the names in the columns

    se                  prints the standard error, if cluster or robust are not specified, default is ols

    CLuster(clustvar)   prints clustered standard errors 

    RObust              prints robust standard errors

    avg(lhs|rhs)        prints the average of the sample of the lhs or the rhs variable

    PNOte               Adds a note with the significance for each p-value at the bottom of the table

    SENote              Adds a note at the bottom of the table with the type of standard error calculation used (if robust or cluster is specified, too)



#### Example 1

    sysuse auto, clear

    uniout displacement length weight = mpg rep78 using myfile, n lab dec(4) se cluster(turn)
