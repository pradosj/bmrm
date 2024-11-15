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

    # Prepare the training set
    x <- cbind(intercept=10,data.matrix(iris[1:4]))
    y <- iris$Species

    # Train multiclass-SVM
    w <- nrbm(ontologyLoss(x,y),LAMBDA=0.01)

    # Performance on training set
    table(y,predict(w,x))


### Evaluate classifier performance with leave-one-out strategy

    svm_loo_pred <- function(x,y,...) {
        parallel::mclapply(mc.cores = 5,seq_along(y),function(i) {
                w <- nrbm(ontologyLoss(x[-i,,drop=FALSE],y[-i]),...)
                pred <- predict(w,x) # Predict all samples
                pred[-i] <- NA # Keep only the left out sample (Put NA everywhere else)
                pred
            }) |>
            simplify2array() |>
            diag()
    }
    loo_pred <- svm_loo_pred(x,y,LAMBDA=0.01)
    table(y,loo_pred) # Contingency matrix

### Sparse models with L1 regularisation

    # Train sparse multiclass-SVM with L1-regularization instead of the default L2-regularization
    w <- nrbmL1(ontologyLoss(x,y),LAMBDA=0.1)
    table(y,predict(w,x))


### Train other type of models

    # Train binary SVM to reconize viriginca species
    w <- nrbm(hingeLoss(x,y=="virginica"),LAMBDA=0.01)
    table(y=="virginica",predict(w,x)) # Performance on training set

    # Train max-margin classifier that maximize AUC, and display ROC curve
    w <- nrbm(rocLoss(x,y=="virginica"),LAMBDA=0.01)
    roc.stat(predict(w,x),y=="virginica") |>
      ggplot(aes(x=sensitivity,y=specificity)) + 
        geom_line() + coord_equal()

    # Train ordinal regression (assuming iris labels are ordered)
    w <- nrbm(ordinalRegressionLoss(x,y),LAMBDA=0.01)
    boxplot(predict(w,x)~y) 



