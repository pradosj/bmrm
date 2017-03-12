BMRM
===============
bmrm is an R package implementing a bundle method for minimization of convex and 
    non-convex risk 
    under L1 or L2 regularization. Implements the algorithm proposed by Teo et 
    al. (JMLR 2010) as well as the extension proposed by Do and Artieres (JMLR 
    2012). The package comes with lot of loss functions for machine learning 
    which make it powerfull for big data analysis. The applications includes:
    structured prediction, linear SVM, multiclass SVM, f-beta optimization, 
    ROC optimization, ordinal regression, quantile regression,
    epsilon insensitive regression, least mean square, logistic regression,
    least absolute deviation regression (see package examples), etc... all with
    L1 and L2 regularization.

Installation
---------------
This R package depends on lpSolve and LowRankQP packages, and requires them to 
be first installed on your system. The documentation is automatically generated
by roxygen2 package.

bmrm can be installed using devtools:

    devtools::install_github("pradosj/bmrm/bmrm")

