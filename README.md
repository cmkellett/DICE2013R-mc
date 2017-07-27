# DICE2013R-mc

DICE2013R-mc replicates the functionality of the DICE2013R code available for download from the website of Prof. W. Nordhaus.  DICE-2013R-mc is a MATLAB implementation of this model that makes use of the CasADi framework for algorithmic differentiation and numeric optimization.  Full details on the implementation and how to use it are in DICE2013R-mc.pdf.

Version 1.1 fixes a compatability issue with CasADi v3.2.1.

Version 2.0 adds the parameter set released by Nordhaus as DICE2016R.  In the main file, set the variable "Param_set" to either 2013 or 2016 for the desired parameter set.
