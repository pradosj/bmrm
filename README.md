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


Usage
---------------

Here is an example to train a linear-multiclass-SVM on `iris`:

    # Prepare training set 
    # Allow the linear models to have a intercept by adding a constant input feature
    x <- cbind(intercept=10,data.matrix(iris[1:4]))
    y <- iris$Species


Train multiclass-SVM:
    w <- nrbm(ontologyLoss(x,y),LAMBDA=0.01)
    table(y,predict(w,x)) # Performance on training set

Train binary SVM to reconize viriginca species:

    w <- nrbm(hingeLoss(x,y=="virginica"),LAMBDA=0.01)
    table(y,predict(w,x)) # Performance on training set

Train binary max-margin linear classifier to reconize virginica and so that AUC is maximal, then display the ROC curve:
    w <- nrbm(rocLoss(x,y=="virginica"),LAMBDA=0.01)
    roc.stat(predict(w,x),y=="virginica") |>
    ggplot(aes(x=sensitivity,y=specificity)) + 
      geom_line() + coord_equal()



