*********************************************************************************
SHEPHERD'S PI CORRELATION   (Version 2012-10-31)
*********************************************************************************

DISCLAIMER: We take no responsibility for any damage or data loss this software
may cause (not that it should). You may reuse and modify these scripts for your 
own purposes but please cite or acknowledge us if you do so. While this software
is not officially supported by us, if you have any questions or suggestions, 
please contact Sam Schwarzkopf (s.schwarzkopf@ucl.ac.uk).

These functions require the Statistics Toolbox from MATLAB.

*********************************************************************************
31 Oct 2012: Modified the corrci.m function so that it can also estimate the 
	     confidence intervals for Spearman's rho.
19 Oct 2012: Corrected a typo in the help section of corrci.m. This does not 
	     affect the calculation or the output of the function.
11 Oct 2012: Corrected a minor bug with ScatterOutliers.m preventing the title
	     of the graph to be displayed (Thanks for Carl Gaspar for pointing 
	     out this error). Added an explanation why we use bsmahal.m instead 
	     of the MATLAB bootstrp function. This does not affect the results.	      
09 Aug 2012: Added the corrci.m function. You can use this for calculating 
	     the nominal 95% confidence interval and for estimating the actual
	     confidence based on your data using a bootstrapping approach.  
08 Aug 2012: Added SetupRand.m script to initialise the random number generator 
  	     for bsmahal.m. This is necessary in order to prevent the boot-
	     strapping from always producing the same results every time you 
	     start MATLAB. You should run this initialisation before running 
	     bootstrapping, simulations or anything relying on randomization.

*********************************************************************************
Shepherd's pi is a robust test of statistical association between two 
variables. It can be used in lieu of Pearson's r or other tests (such as 
Spearman's rho or Kendall's tau) as these tests can be susceptible to the 
presence of influential outliers. It detects outliers by first bootstrapping
the Mahalanobis distance of each data point from the bivariate mean and then
excluding all observations whose distance is >=6. Shepherd's pi is simply
Spearman's rho after outlier removal. The p-value is doubled because the 
removal of outliers can inflate false positive rates. See also:

Schwarzkopf DS, de Haas B & Rees G (2012). Front. Hum. Neurosci.
(http://www.frontiersin.org/Human_Neuroscience/10.3389/fnhum.2012.00200/full)

*********************************************************************************
This archive contains four MATLAB functions. All these functions contain 
usage information for the MATLAB help feature:

SetupRand.m:	    Initializes the state of the random number generator.
		    	YOU SHOULD RUN THIS BEFORE ANYTHING ELSE!!!

Shepherd.m:         The actual Shepherd's pi correlation test.

ScatterOutliers.m:  A function to make scatter plots denoting the points 
                    that were removed as outliers and including a contour 
                    plot of Mahalanobis distances.
                    
bsmahal.m:          The function for bootstrapping Mahalanobis distances.
		    (You could also use the MATLAB bootstrp function here 
		     but this is buggy for large matrices so it is no good
		     for calculating the contour plot in ScatterOutliers.m).
        
round_decs.m:       Simple function for rounding to a decimal. 

corrci.m:	    Function for calculating the nominal 95% confidence interval
		    and for estimating a bootstrapped interval for a correlation.
    
*********************************************************************************
There is also an example data set containing the following variables:

x:          Data drawn from normal distribution
y0:         Data uncorrelated with x
y1:         Data weakly correlated with x
y2:         Data correlated with x
y3:         Data strongly correlated with x

xo1,yo1:    Two uncorrelated variables, which appear correlated under 
            Pearson's r due to the presence of a single outlier
            
xo3,yo3:    Two uncorrelated variables, which appear correlated under 
            Pearson's r due to the presence of three bivariate outliers

*********************************************************************************
Bootstrapping the Mahalanobis distance is more robust than simply using 
the raw Mahalanbois distance for outlier detection. This can be illustrated
by comparing ScatterOutliers(x,y0,1000) with ScatterOutliers(x,y0,0). 
The latter will use the raw Mahalanobis distance (because the number of 
bootstraps is 0). As you can see, the bootstrapped distances are larger, 
because they are less skewed by the outliers in the sample.

*********************************************************************************
