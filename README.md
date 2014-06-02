BMRM
===============
bmrm is an R package implementing "Bundle Methods for Regularized Risk Minimization" proposed by Teo et al., JMLR 2010. This universal data mining framework is particularly useful for structured prediction, classification and regression on big data. The package implements L1 and L2 regularization for: linear SVM, multiclass SVM, f-beta optimization, ROC optimization, ordinal regression, quantile regression, epsilon insensitive regression, least mean square, logistic regression, least absolute deviation (see package examples). Other extensions are possible very easily by implementing custom loss functions.


Installation
---------------
This R package depends on ClpAPI and kernlab packages, and require them to be first installed on your system.
Then bmrm can be build or installed with the commands:
	R CMD build bmrm
	R CMD INSTALL bmrm/

