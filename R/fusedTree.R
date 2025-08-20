#' Construct design data used for fitting fusedTree models
#'
#' Prepares the full data design used to fit a fusedTree model, including
#' dummy-encoded clinical leaf node indicators, optional continuous clinical
#' variables, and a block-diagonal omics matrix structured per tree node.
#'
#' This function allows users to inspect the exact data structure used in
#' fusedTree model fitting. The \code{PenOpt()} and \code{fusedTreeFit()}
#' functions call this function internally so no need to call this function
#' to set-up the right data format. It is just meant for users to check what
#' is going on.
#'
#' @param Tree A fitted tree object, created using \pkg{rpart} or \pkg{partykit}.
#'   Must be an object of class \code{"rpart"} (from the \pkg{rpart} package) or
#'   \code{"constparty"} (from the \pkg{partykit} package).
#' @param X A numeric omics data matrix with dimensions
#'   (sample size × number of omics variables). Must be a \code{matrix}.
#' @param Z A \code{data.frame} of clinical covariates used in tree fitting.
#'   Must be the same data used to construct \code{Tree}.
#' @param LinVars Logical. Whether to include continuous clinical variables
#'   linearly in the model (in addition to tree clustering). Recommended,
#'   as trees may not capture linear effects well. Defaults to \code{TRUE}.
#'
#' @returns A list with the following components:
#' \describe{
#'   \item{Clinical}{A matrix encoding the clinical structure:
#'     \itemize{
#'       \item Dummy variables representing membership to leaf nodes of the tree,
#'       \item Continuous clinical covariates (if \code{LinVars = TRUE}).
#'     }
#'     Each row corresponds to a sample.
#'   }
#'   \item{Omics}{A matrix of omics data per leaf node. This matrix has
#'     dimensions: sample size × (number of leaf nodes × number of omics variables).
#'     For each observation, only the block of omics variables corresponding to its
#'     tree node is populated (other blocks are set to zero).
#'     If the matrix is high-dimensional, a sparse matrix is returned of class
#'     "dgCMatrix" for memory efficiency.}
#' }
#'
#' @export
#'
#' @examples
#' p = 5 # number of omics variables (low for illustration)
#' p_Clin = 5 # number of clinical variables
#' N = 100 # sample size
#' # simulate from Friedman-like function
#' g <- function(z) {
#'   15 * sin(pi * z[,1] * z[,2]) + 10 * (z[,3] - 0.5)^2 + 2 * exp(z[,4]) + 2 * z[,5]
#' }
#' Z <- as.data.frame(matrix(runif(N * p_Clin), nrow = N))
#' X <- matrix(rnorm(N * p), nrow = N)            # omics data
#' betas <- c(1,-1,3,4,2)                         # omics effects
#' Y <- g(Z) + X %*% betas + rnorm(N)             # continuous outcome
#' Y <- as.vector(Y)
#' dat = cbind.data.frame(Y, Z) #set-up data correctly for rpart
#' library(rpart)
#' rp <- rpart::rpart(Y ~ ., data = dat,
#'                    control = rpart::rpart.control(xval = 5, minbucket = 10),
#'                    model = TRUE)
#' cp = rp$cptable[,1][which.min(rp$cptable[,4])] # best model according to pruning
#' Treefit <- rpart::prune(rp, cp = cp)
#' plot(Treefit)
#' Dat_fusedTree <- Dat_Tree(Tree = Treefit, X = X, Z = Z, LinVars = FALSE)
#' Omics <- Dat_fusedTree$Omics
#' Clinical <- Dat_fusedTree$Clinical


Dat_Tree <- function(Tree, X, Z, LinVars = TRUE) {

  ## control statements ##
  if (!(inherits(Tree, "rpart") || inherits(Tree, "constparty"))) {
    stop("Tree must be an rpart (from rpart package) or constparty (from partykit package) object")}

  if (is.null(X) | is.null(Z)) {
    stop("Either X, or Z is unspecified")}

  if (ncol(X) == 0 || nrow(X) == 0) {
    stop("Omics covariates not specified")}

  if (ncol(Z) == 0 || nrow(Z) == 0) {
    stop("Clinical covariates  not specified")}

  if (!inherits(X, "matrix")) {
    stop("X should be specified as matrix object (so not a data.frame)")}

  if (!inherits(Z, "data.frame")) {
    stop("Z should be specified as a data.frame or matrix")}

  if (!is.logical(LinVars)) {
    stop("'LinVars' is not logical, specify as either TRUE or FALSE")}


  if (inherits(Tree, "rpart")) {
    nodes <- row.names(Tree$frame)[treeClust::rpart.predict.leaves(Tree, Z)]
    nodes <- as.numeric(nodes)
    NodeInds <- base::sort(base::unique(as.numeric(row.names(Tree$frame)[Tree$where])))
  } else if (inherits(Tree, "constparty")) {
    nodes <- partykit::predict.party(Tree, newdata = Z, type = "node")
    NodeInds <- base::sort(base::unique(nodes))
  }

  Nodenames <- base::paste0("N", NodeInds)

  p <- ncol(X)
  if (is.null(colnames(X))) {
    colnames(X) <- base::paste0("x",seq(1,ncol(X)))}

  if (is.null(colnames(Z))) {
    colnames(Z) <- base::paste0("z",seq(1,ncol(Z)))}

  namesZ <- colnames(Z)
  namesX <- colnames(X)

  NumNodes <- length(NodeInds)

  if (NumNodes < 2){
    message("Tree has single node, return design matrices")
    return(list(Clinical = stats::model.matrix(~., Z), Omics = X))
  } else {

    ClinIntercepts <- stats::model.matrix(~ 0 + factor(nodes, levels = NodeInds))

    X_tot <- Matrix::t(Matrix::KhatriRao(base::t(X), base::t(ClinIntercepts)))
    colnames(X_tot) <- c(sapply(namesX, function (x) base::paste0(x,base::paste0("_" ,Nodenames))))
    colnames(ClinIntercepts) <- Nodenames
    if (ncol(X) < nrow(X)) {
      X_tot <- as.matrix(X_tot)
    }

    if (isTRUE(LinVars)){
      idVars <- which(sapply(Z, is.numeric) & apply(Z, 2, function (x) length(unique(x)) > 2))
      nameZ <- colnames(Z)[idVars]
      Clinical <- as.matrix(cbind(ClinIntercepts, Z[,idVars]))
      colnames(Clinical)[- (1 : NumNodes)] <- nameZ
    } else {
      Clinical <- ClinIntercepts
    }
    return(list(Clinical = Clinical, Omics = X_tot))
  }
}

#' Create balanced cross-validation folds for hyperparameter tuning
#'
#' Constructs repeated K-fold cross-validation folds, balanced with respect to
#' the fitted tree structure and outcome (if applicable). The folds contain only
#' the test sample indices. This function is useful for tuning penalty parameters
#' in the fusedTree model.
#'
#' For binary and survival outcomes, the function ensures that the proportion
#' of cases vs. controls (or events vs. censored observations) remains
#' relatively constant across folds. In addition, samples are balanced across
#' the leaf nodes of the fitted tree to ensure consistency in node composition
#' between folds.
#'
#' @param Y The response variable. Should be one of:
#'   \itemize{
#'     \item Numeric (for linear regression),
#'     \item Binary (encoded as 0 and 1, for logistic regression),
#'     \item A survival object created using \code{Surv()} (for Cox regression).
#'   }
#'   Only right-censored survival data is currently supported.
#' @param Tree A fitted tree object, created using \pkg{rpart} or \pkg{partykit}.
#'   Must be an object of class \code{"rpart"} (from the \pkg{rpart} package) or
#'   \code{"constparty"} (from the \pkg{partykit} package).
#' @param Z A \code{data.frame} of clinical variables used to fit the tree.
#'   This is used to determine node membership for balancing folds.
#' @param model Character. Specifies the type of outcome model. Must be one of:
#'   \code{"linear"}, \code{"logistic"}, or \code{"cox"}.
#' @param kfold Integer. Number of folds K for cross-validation. Defaults to 5.
#' @param nrepeat Integer. Number of times the K-fold cross-validation is repeated.
#'   Defaults to 3.
#'
#' @returns A list of length \code{kfold × nrepeat}, where each element contains
#' the test indices for a specific fold. These indices can be used to
#' systematically split the data during cross-validation.
#' @export
#'
#' @examples
#' p = 5 # number of omics variables (low for illustration)
#' p_Clin = 5 # number of clinical variables
#' N = 100 # sample size
#' # simulate from Friedman-like function
#' g <- function(z) {
#'   15 * sin(pi * z[,1] * z[,2]) + 10 * (z[,3] - 0.5)^2 + 2 * exp(z[,4]) + 2 * z[,5]
#' }
#' Z <- as.data.frame(matrix(runif(N * p_Clin), nrow = N))
#' X <- matrix(rnorm(N * p), nrow = N)            # omics data
#' betas <- c(1,-1,3,4,2)                         # omics effects
#' Y <- g(Z) + X %*% betas + rnorm(N)             # continuous outcome
#' Y <- as.vector(Y)
#' dat = cbind.data.frame(Y, Z) #set-up data correctly for rpart
#' rp <- rpart::rpart(Y ~ ., data = dat,
#'                    control = rpart::rpart.control(xval = 5, minbucket = 10),
#'                    model = TRUE)
#' cp = rp$cptable[,1][which.min(rp$cptable[,4])] # best model according to pruning
#' Treefit <- rpart::prune(rp, cp = cp)
#' plot(Treefit)
#' folds <- CVfoldsTree(Y = Y, Tree = Treefit, Z = Z, model = "linear")


CVfoldsTree <- function(Y, Tree, Z, model = NULL, kfold = 5, nrepeat = 3) {

  ## ---------------------------------------------------------------------
  ## creates folds for cross-validation. Balance with respect to nodes
  ## and to response (for binary and survival)
  ## ---------------------------------------------------------------------


  #response: response vector, length N (required to balance CV)
  #model: "logistic", "linear", etc

  #kfold: scalar, the number of folds

  #nrepeat: number of repeats of the CV
  #Output: list object with kfold elements containing the sample indices of the left-out samples per fold

  if (!(inherits(Tree, "rpart") || inherits(Tree, "constparty"))) {
    stop("Tree must be an rpart (from rpart package) or constparty (from partykit package) object")
  }

  if (!(model %in% c("linear", "logistic", "cox"))) {
    stop("Model should be specified as linear, logistic, or cox")
  }

  if (ncol(Z) == 0 || nrow(Z) == 0) {
    stop("Clinical covariates not specified")
  }

  if (is.null(Y) | is.null(Z)) {
    stop("Either Y or Z is unspecified")}

  if (length(Y) == 0) {
    stop("Y vector is empty")
  }

  if (!(is.numeric(Y) || inherits(Y, "Surv"))) {
    stop("Y is not a numeric or Surv object. If Y is binary please specify it as numeric vector coded with 0 and 1")
  }

  if (model == "linear") {
    if (inherits(Y, "Surv")) {
      stop("Response is a Surv object, while model = linear is specified")
    }
    if (length(unique(Y)) < 3) {
      stop("Linear model, but Y has maximally 2 distinct values")
    }
  }

  if (model == "logistic") {
    if (inherits(Y, "Surv")) {
      stop("Response is a Surv object, while model = logistic is specified")
    }
    if (!all(Y == 1 | Y == 0)) {
      stop("Logistic model, specify binary response as numeric coded with 0 and 1")
    }
  }

  if (model == "cox") {
    if (!inherits(Y, "Surv")) {
      stop("Response is not Surv object, while model = cox is specified")
    }
  }

  if (!.is_single_positive_integer(kfold)) {
    stop("kfold should be a single positive integer")
  }

  if (kfold < 2) {
    stop("kfold should be at least 2")
  }

  if (!.is_single_positive_integer(nrepeat)) {
    stop("nrepeat should be a single positive integer")
  }

  ## Get node assignments
  if (inherits(Tree, "rpart")) {
    Nodes <- row.names(Tree$frame)[treeClust::rpart.predict.leaves(Tree, Z)]
  } else if (inherits(Tree, "constparty")) {
    Nodes <- partykit::predict.party(Tree, newdata = Z, type = "node")
  }
  Nodes <- factor(Nodes)


  if (model == "logistic"){
    StratDat <- cbind.data.frame(Resp = factor(Y), Node = Nodes)
    Strat <- splitTools::multi_strata(StratDat, strategy = "interaction")
    folds = splitTools::create_folds(Strat, k = kfold, invert = TRUE, m_rep = nrepeat)
  }

  if (model == "linear"){
    StratDat <- cbind.data.frame(Node = Nodes) # for continuous Y only stratify on tree cluster
    Strat <- splitTools::multi_strata(StratDat, strategy = "interaction")
    folds = splitTools::create_folds(Strat, k = kfold, invert = TRUE, m_rep = nrepeat)
  }

  if (model == "cox"){
    StratDat <- cbind.data.frame(Resp = factor(Y[,2]),Node = Nodes)
    Strat <- splitTools::multi_strata(StratDat, strategy = "interaction")
    folds <- splitTools::create_folds(Strat, k = kfold, invert = TRUE, m_rep = nrepeat)
  }
  return(folds)
}

#' Tuning of the penalty parameters of fusedTree using cross-validation
#'
#' Tuning is conducted by optimizing the cross-validated likelihood.
#' Users can either include the fusion penalty (by specifying `alphaInit > 0`),
#' or omit the fusion penalty (by specifying `alphaInit = 0`). If
#' `alphaInit = 0`, only the standard ridge penalty `lambda` is tuned.
#' Note that \code{Dat_Tree()} is called internally so please provide the
#' original data as input arguments.
#'
#' @param Tree A fitted tree object, created using \pkg{rpart} or \pkg{partykit}.
#'   Must be an object of class \code{"rpart"} (from the \pkg{rpart} package) or
#'   \code{"constparty"} (from the \pkg{partykit} package).
#' @param X The original omics data matrix. Has dimensions (sample size × number
#'   of omics variables). Should be a matrix.
#' @param Z The original clinical data matrix, which was used to fit the tree.
#'   Should be a \code{data.frame}.
#' @param Y The response; should be either numeric, binary (encoded by 0 and 1),
#'   or a survival object created by \code{Surv()} from the \code{survival} package.
#'   Only right-censored survival data is allowed.
#' @param model Character. Specifies the outcome model. One of \code{"linear"},
#'   \code{"logistic"}, or \code{"cox"}.
#' @param lambdaInit Numeric. Initial value for the standard ridge (L2) penalty
#'   \code{lambda}. Must be greater than zero. Defaults to 10.
#' @param alphaInit Numeric. Initial value for the fusion penalty \code{alpha}.
#'   If set to 0, fusion is omitted and only \code{lambda} is tuned. Must be zero
#'   or greater. Defaults to 10.
#' @param folds List. Each element contains the indices of the test samples for
#'   a fold. It is advisable to balance the samples with respect to the outcome
#'   (for binary and survival models) and the tree structure. If not provided,
#'   folds are generated internally.
#' @param loss Character. The loss function to optimize in cross-validation.
#'   For binary and survival outcomes, only \code{"loglik"} (cross-validated
#'   likelihood) is supported. For continuous outcomes, an alternative is
#'   \code{"sos"} (sum of squares loss). Defaults to \code{"loglik"}.
#' @param symFusion Logical. Whether fusion should be symmetric across nodes.
#' Setting this parameter to \code{FALSE} induces asymmetric fusion. That is,
#' nodes having more similar predictions (using only clinical variables)
#' of the response will be shrunk more to each other than nodes having more
#' distinct predictions. Setting to \code{FALSE} is not tested very well so
#' use on your own risk. Defaults to \code{TRUE}.
#' @param multistart Logical. Whether to initialize with different starting values when
#'   optimizing the cross-validated likelihood. Can help with stability when both
#'   \code{lambda} and \code{alpha} are tuned, at the cost of longer run time.
#'   Defaults to \code{FALSE}.
#' @param maxIter Integer. Maximum number of iterations for the IRLS (iterative
#'   reweighted least squares) algorithm. Used only for logistic and Cox models.
#'   Defaults to 30.
#' @param LinVars Logical. Whether to include continuous clinical variables
#'   linearly in the model (in addition to the tree structure). Can be helpful
#'   since trees may not capture linear effects well. Defaults to \code{TRUE}.
#'
#' @returns A numeric vector with the tuned values of the penalties:
#'   \itemize{
#'     \item \code{lambda}: standard ridge (L2) penalty.
#'     \item \code{alpha}: fusion penalty (only if \code{alphaInit > 0}).
#'   }
#'   If \code{alphaInit = 0}, only the tuned \code{lambda} is returned.
#'
#' @references
#' \CRANpkg{porridge}
#'
#' @details
#' The cross-validated likelihood is optimized using the \code{Nelder-Mead}
#' method from \code{stats::optim()}. When tuning both \code{lambda} and
#' \code{alpha}, the objective function can be noisy. Setting
#' \code{multistart = TRUE} performs optimization from several starting values
#' to improve robustness. This is only applicable when \code{alphaInit > 0}.
#'
#' @export
#'
#' @examples
#' p = 5 # number of omics variables (low for illustration)
#' p_Clin = 5 # number of clinical variables
#' N = 100 # sample size
#' # simulate from Friedman-like function
#' g <- function(z) {
#'   15 * sin(pi * z[,1] * z[,2]) + 10 * (z[,3] - 0.5)^2 + 2 * exp(z[,4]) + 2 * z[,5]
#' }
#' set.seed(11)
#' Z <- as.data.frame(matrix(runif(N * p_Clin), nrow = N))
#' X <- matrix(rnorm(N * p), nrow = N)            # omics data
#' betas <- c(1,-1,3,4,2)                         # omics effects
#' Y <- g(Z) + X %*% betas + rnorm(N)             # continuous outcome
#' Y <- as.vector(Y)
#' dat = cbind.data.frame(Y, Z) #set-up data correctly for rpart
#' rp <- rpart::rpart(Y ~ ., data = dat,
#'                    control = rpart::rpart.control(xval = 5, minbucket = 10),
#'                    model = TRUE)
#' cp = rp$cptable[,1][which.min(rp$cptable[,4])] # best model according to pruning
#' Treefit <- rpart::prune(rp, cp = cp)
#' plot(Treefit)
#' folds <- CVfoldsTree(Y = Y, Tree = Treefit, Z = Z, model = "linear")
#' optPenalties <- PenOpt(Tree = Treefit, X = X, Y = Y, Z = Z,
#'                        model = "linear", lambdaInit = 10, alphaInit = 10,
#'                        loss = "loglik",
#'                        LinVars = FALSE,
#'                        folds = folds, multistart = FALSE)
#' optPenalties


PenOpt <- function(Tree, X, Y, Z, model = NULL,
                   lambdaInit = 10, alphaInit = 10,
                   folds = CVfoldsTree(Y = Y, Tree = Tree, Z = Z, model = model),
                   loss = "loglik",
                   symFusion = TRUE,
                   multistart = FALSE,
                   maxIter = 30,
                   LinVars = FALSE) {

  ## control statements ##
  if (!(inherits(Tree, "rpart") || inherits(Tree, "constparty"))) {
    stop("Tree must be an rpart (from rpart package) or constparty (from partykit package) object")
  }

  if (!(model %in% c("linear", "logistic", "cox"))) {
    stop("Model should be specified as linear, logistic, or cox")
  }

  if (!(loss %in% c("loglik", "sos"))) {
    stop("Loss should be specified as loglik or sos")
  }

  if (ncol(X) == 0 || nrow(X) == 0) {
    stop("Omics covariates not specified")
  }

  if (ncol(Z) == 0 || nrow(Z) == 0) {
    stop("Clinical covariates not specified")
  }

  if (length(Y) == 0) {
    stop("Y vector is empty")
  }

  if (is.null(Y) | is.null(X) | is.null(Z)) {
    stop("Either Y, X, or Z is unspecified")}

  if (!(is.numeric(Y) || inherits(Y, "Surv"))) {
    stop("Y is not a numeric or Surv object. If Y is binary please specify it as numeric vector coded with 0 and 1")
  }

  if (!inherits(X, "matrix")) {
    stop("X should be specified as matrix object (so not a data.frame)")
  }

  if (inherits(Z, "matrix")) {
    Z <- data.frame(Z)
  }
  if (!inherits(Z, "data.frame")) {
    stop("Z should be specified as a data.frame or matrix")
  }

  if (!is.logical(LinVars)) {
    stop("LinVars is not logical, specify as either TRUE or FALSE")
  }

  if (!is.logical(symFusion)) {
    stop("symFusion is not logical, specify as either TRUE or FALSE")
  }

  if (!is.logical(multistart)) {
    stop("multistart is not logical, specify as either TRUE or FALSE")
  }

  if (!inherits(folds, "list")) {
    stop("Test fold ids should be collected in a list")
  }

  if (length(folds) < 3) {
    stop("Number of folds should be at least 3")
  }

  folds_classes <- unlist(lapply(folds, class))
  if (!all(folds_classes %in% c("integer", "numeric"))) {
    stop("Test fold ids should be specified as integer or numeric")
  }

  if (!all(unlist(folds) <= length(Y))) {
    stop("Fold id out of range")
  }

  if (model == "linear") {
    if (inherits(Y, "Surv")) {
      stop("Response is a Surv object, while model = linear is specified")
    }
    if (length(unique(Y)) < 3) {
      stop("Linear model, but Y has maximally 2 distinct values")
    }
  }

  if (model == "logistic") {
    if (inherits(Y, "Surv")) {
      stop("Response is a Surv object, while model = logistic is specified")
    }
    if (!all(Y == 1 | Y == 0)) {
      stop("Logistic model, specify binary response as numeric coded with 0 and 1")
    }
  }

  if (model == "cox") {
    if (!inherits(Y, "Surv")) {
      stop("Response is not Surv object")
    }
  }

  if (!is.numeric(lambdaInit) || length(lambdaInit) != 1) {
    stop("`lambdaInit` must be a single numeric value.")
  }

  if (!is.numeric(alphaInit) || length(alphaInit) != 1) {
    stop("`alphaInit` must be a single numeric value.")
  }

  if (lambdaInit <= 0) {
    stop("initial penalty lambdaInit should be larger than zero")
  }

  if (alphaInit < 0) {
    stop("initial penalty alphaInit should be zero or positive")
  }

  if (!.is_single_positive_integer(maxIter)) {
    stop("maxiter should be a single positive integer")
  }


  ## obtaining tree stats
  if (inherits(Tree, "rpart")) {
    nodes <- treeClust::rpart.predict.leaves(Tree, Z)
  } else if (inherits(Tree, "constparty")) {
    nodes <- partykit::predict.party(Tree, newdata = Z, type = "node")
  }

  NumNodes <- length(unique(nodes))

  if (NumNodes < 2){
    stop("Only single node present in the tree.
          Please consider a different model as the tree does not find
          any signal in the clinical variables")
  }


  if (NumNodes > 1){
    Dat <- Dat_Tree(Tree = Tree, X = X, Z = Z,
                   LinVars = LinVars)

    X1 <- Dat$Omics; U1 <- Dat$Clinical;

    remove(Dat)

    NumNod <- ncol(X1) / ncol(X)

    if (alphaInit == 0){
      Delta <- NULL
    }
    if (alphaInit > 0){
      Delta = .PenMatr(NumNodes = NumNod, p = ncol(X),
                       symmetric = symFusion, Tree = Tree)
      if (ncol(Delta) != ncol(X1)){
        stop("number of columns of penalty matrix does not equal number of
           columns of design matrix X")}
    }


    optPenalties <- .optPenaltyGLM.kCVauto2(Y = Y, X = X1, U = U1,
                                            lambdaInit = lambdaInit,
                                            lambdaGinit = alphaInit, Dg = Delta,
                                            model = model, folds = folds,
                                            symFusion = symFusion,
                                            Tree = Tree,
                                            loss = loss,
                                            multistart = multistart,
                                            maxIter = maxIter)
    if (length(optPenalties) == 1) {
      names(optPenalties)<-c("lambda")
    }
    if (length(optPenalties) == 2) {
      names(optPenalties)<-c("lambda","alpha")
    }

    return(optPenalties)
  }
}

#' Fit a fusedTree model with or without fusion penalty
#'
#' Fits a fusedTree model by solving a penalized regression problem using
#' either a linear, logistic, or Cox model. The model includes both a standard
#' ridge (L2) penalty and a fusion penalty to encourage similarity
#' between leaf node-specific omics effects. The fusion penalty can also be
#' omitted by specifying `alpha = 0`.
#'
#' @param Tree A fitted tree object, created using \pkg{rpart} or \pkg{partykit}.
#'   Must be an object of class \code{"rpart"} (from the \pkg{rpart} package) or
#'   \code{"constparty"} (from the \pkg{partykit} package).
#' @param X A matrix of omics data with dimensions (sample size × number of omics variables).
#' @param Z A data frame of clinical covariates used to fit the tree. Must be a
#'   \code{data.frame}, not a matrix.
#' @param Y The response variable. Can be:
#'   \itemize{
#'     \item numeric (for linear regression),
#'     \item binary (0/1, for logistic regression),
#'     \item a survival object created by \code{Surv()} (right-censored data only).
#'   }
#' @param LinVars Logical. Whether to include continuous clinical variables
#'   linearly in the model (in addition to the tree structure). Defaults to \code{TRUE}.
#' @param model Character. Specifies the type of outcome model to fit.
#'   One of: \code{"linear"}, \code{"logistic"}, or \code{"cox"}.
#' @param lambda Numeric. Value for the standard ridge (L2) penalty.
#' @param alpha Numeric. Value for the fusion penalty.
#' The fusion penalty is not incorporated when `alpha = 0` is specified.
#' @param symFusion Logical. Whether fusion should be symmetric across nodes.
#' Setting this parameter to \code{FALSE} induces asymmetric fusion. That is,
#' nodes having more similar predictions (using only clinical variables)
#' of the response will be shrunk more to each other than nodes having more
#' distinct predictions. Setting to \code{FALSE} is not tested very well so
#' use on your own risk. Defaults to \code{TRUE}.
#' @param maxIter Integer. Maximum number of iterations for the IRLS (iterative
#'   reweighted least squares) algorithm. Used only when \code{model = "logistic"} or \code{"cox"}.
#'   Defaults to 50.
#' @param minSuccDiff Numeric. The minimum difference in log-likelihood between
#'   successive iterations of IRLS to declare convergence. Only used when
#'   \code{model = "logistic"} or \code{"cox"}.
#' @param dat Logical. Whether to return the data used in model fitting
#'   (i.e., omics, clinical, and response). Defaults to \code{FALSE}.
#' @param verbose Logical. Whether to print progress updates from the IRLS algorithm.
#'   Only applies to \code{model = "logistic"} or \code{model ="cox"}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Tree}{The fitted tree object from 'rpart'.}
#'   \item{Effects}{A named numeric vector of estimated effect sizes, including:
#'     intercepts (tree leaf nodes), omics effects (per node), and linear
#'     clinical effects (if \code{LinVars = TRUE}).}
#'   \item{Breslow}{(Optional) The Breslow estimates of the baseline hazard \code{ht}
#'      and the cumulative baseline hazard \code{Ht} for each time point. Only returned for
#'      \code{model = "cox"}.}
#'   \item{Parameters}{A list of model parameters used in fitting (e.g.,
#'     \code{lambda}, \code{alpha}, \code{model}, etc.).}
#'   \item{Clinical}{(Optional) The clinical design matrix used in fitting, if \code{dat = TRUE}.}
#'   \item{Omics}{(Optional) The omics design matrix used in fitting, if \code{dat = TRUE}.}
#'   \item{Response}{(Optional) The response vector used in fitting, if \code{dat = TRUE}.}
#' }
#' The returned list object is of class S3 for which predict() is available
#'
#' @details
#' \strong{Linear model:} Estimated using a closed-form analytic solution.
#'
#' \strong{Logistic and Cox models:} Estimated using IRLS (iterative reweighted least
#' squares), equivalent to the Newton-Raphson algorithm.
#'
#' \strong{Cox model:} The full likelihood approach is used, following
#' van Houwelingen et al. (2005). See also van de Wiel et al. (2021) for
#' additional details on penalized regression for survival outcomes.
#'
#' @references
#' \CRANpkg{porridge}
#'
#' van Houwelingen, H. C., et al.. (2005). Cross-validated Cox regression on microarray gene expression data.
#' *Stad Med*
#'
#' van de Wiel, M. A., et al. (2021). Fast Cross-validation for Multi-penalty High-dimensional Ridge Regression.
#' *J Comput Graph Stat*
#'
#' @export
#'
#' @examples
#' p = 5 # number of omics variables (low for illustration)
#' p_Clin = 5 # number of clinical variables
#' N = 100 # sample size
#' # simulate from Friedman-like function
#' g <- function(z) {
#'   15 * sin(pi * z[,1] * z[,2]) + 10 * (z[,3] - 0.5)^2 + 2 * exp(z[,4]) + 2 * z[,5]
#' }
#' set.seed(11)
#' Z <- as.data.frame(matrix(runif(N * p_Clin), nrow = N))
#' X <- matrix(rnorm(N * p), nrow = N)            # omics data
#' betas <- c(1,-1,3,4,2)                         # omics effects
#' Y <- g(Z) + X %*% betas + rnorm(N)             # continuous outcome
#' Y <- as.vector(Y)
#' dat = cbind.data.frame(Y, Z) #set-up data correctly for rpart
#' rp <- rpart::rpart(Y ~ ., data = dat,
#'                    control = rpart::rpart.control(xval = 5, minbucket = 10),
#'                    model = TRUE)
#' cp = rp$cptable[,1][which.min(rp$cptable[,4])] # best model according to pruning
#' Treefit <- rpart::prune(rp, cp = cp)
#' plot(Treefit)
#' folds <- CVfoldsTree(Y = Y, Tree = Treefit, Z = Z, model = "linear")
#' optPenalties <- PenOpt(Tree = Treefit, X = X, Y = Y, Z = Z,
#'                        model = "linear", lambdaInit = 10, alphaInit = 10,
#'                        loss = "loglik",
#'                        LinVars = FALSE,
#'                        folds = folds, multistart = FALSE)
#' optPenalties
#'
#' # with fusion
#' fit <- fusedTree(Tree = Treefit, X = X, Y = Y, Z = Z,
#'                     LinVars = FALSE, model = "linear",
#'                     lambda = optPenalties[1],
#'                     alpha = optPenalties[2])
#' # without fusion
#' fit1 <- fusedTree(Tree = Treefit, X = X, Y = Y, Z = Z,
#'                      LinVars = FALSE, model = "linear",
#'                      lambda = optPenalties[1],
#'                      alpha = 0)
#' #compare effect estimates
#' fit$Effects
#' fit1$Effects


fusedTree <- function(Tree, X, Y, Z, LinVars = TRUE, model,
                      lambda, alpha, symFusion = TRUE,
                      maxIter = 50, minSuccDiff = 10^(-10),
                      dat = FALSE, verbose = TRUE) {

  ## control statements ##
  if (!(inherits(Tree, "rpart") || inherits(Tree, "constparty"))) {
    stop("Tree must be an rpart (from rpart package) or constparty (from partykit package) object")
  }

  if (!(model %in% c("linear", "logistic", "cox"))) {
    stop("Model should be specified as linear, logistic, or cox")
  }

  if (ncol(X) == 0 || nrow(X) == 0) {
    stop("Omics covariates not specified")}

  if (ncol(Z) == 0 || nrow(Z) == 0) {
    stop("Clinical covariates  not specified")}

  if (length(Y) == 0) {
    stop("Y vector has length zero")}

  if (is.null(Y) | is.null(X) | is.null(Z)) {
    stop("Either Y, X, or Z is unspecified")}

  if (!(is.numeric(Y) || inherits(Y, "Surv"))) {
    stop("Y is not a numeric or Surv object. If Y is binary please specify it as numeric vector coded with 0 and 1")}

  if (!inherits(X, "matrix")) {
    stop("X should be specified as matrix object (so not a data.frame)")}

  if (inherits(Z, "matrix")) {
    Z <- data.frame(Z)}

  if (!inherits(Z, "data.frame")) {
    stop("Z should be specified as a data.frame or matrix")}

  if (!is.logical(LinVars)) {
    stop("'LinVars' is not logical, specify as either TRUE or FALSE")}

  if (!is.logical(symFusion)) {
    stop("'symFusion' is not logical, specify as either TRUE or FALSE")}

  if (!is.logical(verbose)) {
    stop("'verbose' is not logical, specify as either TRUE or FALSE")}

  if (!is.logical(dat)) {
    stop("'dat' is not logical, specify as either TRUE or FALSE")}

  if (model == "linear"){
    if (inherits(Y, "Surv")) {
      stop("Response is a Surv object, while model = logistic is specified")}
    if (length(unique(Y)) < 3) {
      stop("Linear model, but Y has maximally 2 distinct values")}
  }

  if (model == "logistic"){
    if (inherits(Y, "Surv")) {
      stop("Response is a Surv object, while model = logistic is specified")}
    if(!all(Y==1 | Y==0)) {
      stop("Logistic model, specify binary response as numeric coded with 0 and 1")}
  }

  if (model == "cox"){
    if (!inherits(Y, "Surv")) {stop("Response is not Surv object")}
  }

  if (!is.numeric(lambda) || length(lambda) != 1) {
    stop("`lambda` must be a single numeric value.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1) {
    stop("`alpha` must be a single numeric value.")
  }

  if(lambda <= 0){
    stop("Ridge penalty 'lambda' should be larger than zero")}

  if(alpha < 0){
    stop("Fusion penalty 'alpha' should be zero or positive")}

  if (!is.numeric(minSuccDiff) || length(minSuccDiff) != 1) {
    stop("`minSuccDiff` must be a single numeric value.")
  }

  if (minSuccDiff <= 0 || minSuccDiff >= 1e-4) {
    stop("minSuccDiff should be larger than 0 and smaller than 10^-4")
  }

  if (!.is_single_positive_integer(maxIter)) {
    stop("maxiter should be a single positive integer")
  }


  ## obtaining tree stats
  ## obtaining tree stats
  if (inherits(Tree, "rpart")) {
    nodes <- treeClust::rpart.predict.leaves(Tree, Z)
  } else if (inherits(Tree, "constparty")) {
    nodes <- partykit::predict.party(Tree, newdata = Z, type = "node")
  }

  NumNodes <- length(unique(nodes))

  if (NumNodes < 2){
    stop("Only single node present in the tree.
          Please consider a different model as the tree does not find
          any signal in the clinical variables")
  }


  if (NumNodes > 1){

    Dat <- Dat_Tree(Tree = Tree, X = X, Z = Z,  LinVars = LinVars)

    X1 <- Dat$Omics; U1 <- Dat$Clinical
    remove(Dat)
    NumNod <- ncol(X1)/ncol(X)

    if (alpha == 0){
      Delta <- NULL
    }
    if (alpha > 0){
      Delta = .PenMatr(NumNodes = NumNod, p = ncol(X),
                       symmetric = symFusion, Tree = Tree)
      if (ncol(Delta) != ncol(X1)){
        stop("number of columns of penalty matrix does not equal number of
           columns of design matrix X")}
    }

    Fit = .ridgeGLM2(Y = Y, U = U1, X = X1,
                     lambda = lambda, lambdaG = alpha, Dg = Delta,
                     model = model, symFusion = symFusion, Tree = Tree,
                     maxIter = maxIter, minSuccDiff = minSuccDiff,
                     verbose = verbose)

    if (model == "linear" || model == "logistic") {
      names(Fit) <- c(colnames(U1),colnames(X1))
      Pars = cbind.data.frame("Model" = model, "LinVar" = LinVars, "Alpha" = alpha,"Lambda" = lambda)
      if (isTRUE(dat)){
        res <- list(Tree = Tree, Effects = Fit, Pars = Pars, Omics = X1, Clinical = U1, Response = Y)
      } else {
        res <- list(Tree = Tree, Effects = Fit, Pars = Pars)}
    }
    if (model == "cox") {
      Effects <- Fit$estimates
      Breslow <- Fit$Breslow
      names(Effects) <- c(colnames(U1),colnames(X1))
      Pars = cbind.data.frame("Model" = model, "LinVar" = LinVars, "Alpha" = alpha,"Lambda" = lambda)
      if (isTRUE(dat)){
        res <- list(Tree = Tree, Effects = Effects, Breslow = Breslow, Pars = Pars, Omics = X1, Clinical = U1, Response = Y)
      } else {
        res <- list(Tree = Tree, Effects = Effects, Breslow = Breslow, Pars = Pars)}
    }
  }
  class(res) <- "fusedTree"
  return(res)
}

#####################################
#### Methods for fusedTree class ####
#####################################
#' Predict Method for Fused Tree Models
#'
#' Generates predictions from a fitted `fusedTree` object using new clinical and omics data.
#'
#' @param object An object of class \code{"fusedTree"} returned by the \code{\link{fusedTree}} function.
#' @param newX A matrix of new omics covariates for prediction. Must have the same number of columns (variables) as used in the model fitting.
#' @param newZ A data frame of new clinical covariates for prediction.
#' @param newY Optional input that is only used for survival response.
#' Defaults to \code{NULL}. If provided, it should be a
#' \code{\link[survival]{Surv}} object. In that case, \code{newY} is used
#' to interpolate the baseline hazard to the event times of the test data.
#' @param ... Currently not used. Included for S3 method consistency.
#'
#' @return
#' A model-specific prediction object:
#' \itemize{
#'   \item For \code{"linear"} models: a \code{data.frame} with a single column \code{Ypred} (predicted values).
#'   \item For \code{"logistic"} models: a \code{data.frame} with columns
#'   \code{Probs} (predicted probabilities), and \code{LinPred} (linear predictor).
#'   \item For \code{"cox"} models: a \code{list} with two elements:
#'     \describe{
#'       \item{data}{A \code{data.frame} with a single column \code{LinPred} (linear predictor).}
#'       \item{Survival}{A matrix of predicted survival probabilities. Rows = test subjects, columns = unique times from \code{newY}.
#'       Only returned when \code{newY} is provided.}
#'     }
#' }
#'
#' @export
#'
#' @seealso \code{\link{fusedTree}} for model fitting.
#'
#' @examples
#' p = 5 # number of omics variables (low for illustration)
#' p_Clin = 5 # number of clinical variables
#' N = 100 # sample size
#' # simulate from Friedman-like function
#' g <- function(z) {
#'   15 * sin(pi * z[,1] * z[,2]) + 10 * (z[,3] - 0.5)^2 + 2 * exp(z[,4]) + 2 * z[,5]
#' }
#' set.seed(11)
#' Z <- as.data.frame(matrix(runif(N * p_Clin), nrow = N))
#' X <- matrix(rnorm(N * p), nrow = N)            # omics data
#' betas <- c(1,-1,3,4,2)                         # omics effects
#' Y <- g(Z) + X %*% betas + rnorm(N)             # continuous outcome
#' Y <- as.vector(Y)
#' dat = cbind.data.frame(Y, Z) #set-up data correctly for rpart
#' rp <- rpart::rpart(Y ~ ., data = dat,
#'                    control = rpart::rpart.control(xval = 5, minbucket = 10),
#'                    model = TRUE)
#' cp = rp$cptable[,1][which.min(rp$cptable[,4])] # best model according to pruning
#' Treefit <- rpart::prune(rp, cp = cp)
#' fit <- fusedTree(Tree = Treefit, X = X, Y = Y, Z = Z,
#'                     LinVars = FALSE, model = "linear",
#'                     lambda = 10,
#'                     alpha = 1000)
#' Preds <- predict(fit, newX = X, newZ = Z, newY = Y)


predict.fusedTree <- function(object, newX, newZ, newY = NULL, ...) {

  # ---- Check inputs ----
  if (!inherits(object, "fusedTree")) stop("Object must be of class 'fusedTree'.")
  if (is.null(newX) | is.null(newZ)) {stop("Either newY, newX, or newZ is unspecified")}
  if (!is.matrix(newX)) stop("newX must be a matrix.")
  if (!is.data.frame(newZ)) stop("newZ must be a data frame.")
  #if (!(is.numeric(newY) || inherits(newY, "Surv"))) stop("newY must be numeric or a Surv object.")
  if (ncol(newX) == 0 || nrow(newX) == 0) stop("newX has zero columns or rows.")
  if (ncol(newZ) == 0 || nrow(newZ) == 0) stop("newZ has zero columns or rows.")
  #if (length(newY) == 0) stop("newY is empty.")


  # ---- Extract components from model object ----
  Tree     <- object$Tree
  Ests     <- object$Effects
  Pars     <- object$Pars
  Linvars  <- Pars$LinVar
  model    <- Pars$Model
  Breslow  <- if (model == "cox") object$Breslow else NULL

  # ---- Generate new data in tree structure ----
  Dat_test <- Dat_Tree(Tree = Tree, X = newX, Z = newZ, LinVars = Linvars)
  Xtest    <- Dat_test$Omics
  Ztest    <- Dat_test$Clinical

  # ---- Align predictors with model coefficients ----
  keep_x   <- colnames(Xtest) %in% names(Ests)
  Xtest    <- Xtest[, keep_x, drop = FALSE]

  keep_ids <- names(Ests) %in% c(colnames(Ztest), colnames(Xtest))
  Ests     <- Ests[keep_ids]

  # ---- Compute linear predictor ----
  LP <- as.numeric(cbind(Ztest, Xtest) %*% Ests)

  # ---- Model-specific predictions ----
  if (model == "linear") {
    Ypred <- LP

    return(data.frame(Ypred = Ypred))

  } else if (model == "logistic") {
    if (!all(newY == 1 | newY == 0)) {
      stop("Logistic model, specify newY as numeric coded with 0 and 1")
    }
    Ypred <- 1 / (1 + exp(-LP))
    return(data.frame(Probs = Ypred, LinPred = LP))

  } else if (model == "cox") {

    if (!is.null(newY)){
      if (!inherits(newY, "Surv")) {
        stop("For Cox models, newY must be a Surv object.")
      }
      t_test <- sort(unique(newY[, 1]))
      #t_test <- seq(0, max(newY[,1]), length.out = 100)

      # Interpolate cumulative baseline hazard
      ord <- order(Breslow$time)
      time_sorted <- Breslow$time[ord]
      Ht_sorted   <- Breslow$Ht[ord]

      H0_fun <- stats::stepfun(time_sorted, c(0, Ht_sorted), right = TRUE)
      H0_test <- H0_fun(t_test)

      # Compute survival probabilities
      risk <- exp(LP)
      S_test <- outer(risk, H0_test, function(r, H) exp(-r * H))

      colnames(S_test) <- paste0("t=", sprintf("%.3f", t_test))

      rownames(S_test) <- paste0("Subject_", seq_len(nrow(S_test)))
      return(list(
        LinPred = data.frame(LinPred = LP),
        Survival = S_test
      ))
    } else {
      return(data.frame(LinPred = LP))
    }
  } else {
    stop("Unsupported model type: ", model)
  }
}

###############################################################################
################################### Auxiliary Functions #######################
###############################################################################

.PenMatr <- function(NumNodes, p, symmetric = TRUE, Tree = NULL){

  ## ---------------------------------------------------------------------
  ## Creates the right penalty matrix for given dimension of omics vars (p)
  ## and numer of leaf nodes (NumNodes)
  ## ---------------------------------------------------------------------

  if (isTRUE(symmetric)) {

    block <- base::diag(1,nrow = NumNodes, ncol = NumNodes) -
      1 / NumNodes * matrix(1, nrow = NumNodes,ncol = NumNodes)
  }

  if (isFALSE(symmetric)) {
    blocksqrt <- base::diag(1,nrow = NumNodes, ncol = NumNodes) -
      .compute_weight_matrix(Tree = Tree)
    block <- base::t(blocksqrt) %*% blocksqrt
  }

  Omega <- Matrix::kronecker(Matrix::Diagonal(p, 1), block)
  #Omega <- as.matrix(Omega)
  return(Omega)
}

.EigPenmatr <- function(NumNodes, p, symmetric = TRUE, Tree = NULL){

  ## ---------------------------------------------------------------------
  ## Computes eigen decomposition efficiently of the penalty matrix
  ## required for efficient cross-validation
  ## ---------------------------------------------------------------------

  Eig_base = base::eigen(.PenMatr(NumNodes = NumNodes, p = 1,
                                  symmetric = symmetric,
                                  Tree = Tree))$vectors

  Eigvecs = Matrix::kronecker(Matrix::Diagonal(p,1),Eig_base)
  diags = base::rep(c(base::rep(1, NumNodes - 1), 0), p)
  return(list(values = diags,vectors = as.matrix(Eigvecs)))
}

.compute_weight_matrix <- function(Tree) {
  # Extract predicted values for each terminal node
  if (!(inherits(Tree, "rpart") || inherits(Tree, "constparty"))) {
    stop("Tree must be an rpart (from rpart package) or constparty (from partykit package) object")}



  if (inherits(Tree, "rpart")) {
    method <- Tree$method

    # regression and survival
    if (method %in% c("anova", "poisson", "exp")) {
      terminal_preds <- Tree$frame$yval[Tree$frame$var == "<leaf>"]

    } else if (method == "class") {
      res <- Tree$frame$yval2
      terminal_preds <- res[Tree$frame$var == "<leaf>", ncol(res) - 1]
    } else {
      stop("Unsupported rpart method: ", method)
    }

  } else if (inherits(Tree, "constparty")) {
    terminal_nodes <- partykit::nodeids(Tree, terminal = TRUE)
    terminal_preds <- partykit::predict_party(Tree, terminal_nodes, type = "response")

    # distinguish for classification case
    if (is.factor(terminal_preds)) {
      terminal_preds <- partykit::predict_party(Tree, terminal_nodes, type = "prob")[,2]
    }
  }

  M <- base::length(terminal_preds)
  # Compute ranks
  ranks <- base::rank(terminal_preds)

  rank_diff <- matrix(0, nrow = M, ncol = M)



  # Fill upper triangle and copy to lower (symmetry)
  for (m in 1:(M - 1)) {
    for (m2 in (m + 1):M) {
      diff <- (ranks[m] - ranks[m2])^2
      rank_diff[m, m2] <- diff
      rank_diff[m2, m] <- diff
    }
  }

  # Set diagonal to 1
  # Convert to weight matrix

  base::diag(rank_diff) <- 0
  w_matrix <- 1 / (1 + rank_diff)

  # Normalize
  w_matrix <- w_matrix/base::rowSums(w_matrix)

  return(w_matrix)
}

.is_single_positive_integer <- function(x) {
  is.numeric(x) && length(x) == 1 && x == as.integer(x) && x > 0
}

###############################################################################
######## Auxiliary Functions for ridge estimation. Strong inspiration #########
######## from the porridge package developed by van Wieringen         #########
###############################################################################

.ridgeGLM2 <- function(Y, X, U,
                       lambda, lambdaG = 0,
                       Dg = NULL,
                       model, symFusion = TRUE, Tree = NULL,
                       minSuccDiff = 10^(-10), maxIter = 100, verbose = FALSE){

  ## ---------------------------------------------------------------------
  ## Ridge estimation of the regression parameter of the
  ## generalized linear model (GLM)
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y           : A numeric being the response vector.
  ## X           : A matrix, the design matrix. The number of rows should
  ##               match the number of elements of Y.
  ## lambda      : A numeric, the ridge penalty parameter.

  ## model       : A character indicating which generalized linear model
  ##               model instance is to be fitted.
  ## minSuccDiff : A numeric, the minimum distance between two successive
  ##               estimates to be achieved.
  ## maxIter     : A numeric specifying the maximum number of iterations.
  ## ---------------------------------------------------------------------
  ## Value:
  ## The ridge estimate of the regression parameter, a numeric.
  ## ---------------------------------------------------------------------
  ## Authors      : Wessel N. van Wieringen, modified by JM Goedhart
  ## ---------------------------------------------------------------------
  if (is.null(Y) | is.null(X) | is.null(U)) {stop("Either Y, X, or U is unspecified")}
  if (ncol(X) == 0 | nrow(X) == 0) {stop("Omics matrix has zero columns or rows")}
  if (ncol(U) == 0 | nrow(U) == 0) {stop("Clinical matrix has zero columns or rows")}
  if (length(Y) == 0) {stop("Length Y vector is zero")}



  if (model == "linear"){
    return(.ridgeLM(Y = Y, X = X, U = U, lambda = lambda,
                    lambdaG = lambdaG, Dg = Dg,
                    symFusion = symFusion,Tree = Tree,
                    verbose = verbose))
  }

  if (model == "logistic"){
    return(.ridgeBLM(Y = Y, X = X, U = U, lambda = lambda,
                     lambdaG = lambdaG, Dg = Dg,
                     symFusion = symFusion,Tree = Tree,
                     minSuccDiff = minSuccDiff, maxIter = maxIter,
                     verbose = verbose))
  }

  if (model == "cox"){
    return(.ridgeSurv(Y = Y, X = X, U = U, lambda = lambda,
                      lambdaG = lambdaG, Dg = Dg,
                      symFusion = symFusion,Tree = Tree,
                      minSuccDiff = minSuccDiff, maxIter = maxIter,
                      verbose = verbose))
  }
}


.optPenaltyGLM.kCVauto2 <- function(Y,
                                    X,
                                    U,
                                    lambdaInit,
                                    lambdaGinit,
                                    Dg,
                                    model = "linear",
                                    folds,
                                    symFusion = TRUE,
                                    Tree = NULL,
                                    loss = "loglik",
                                    multistart = FALSE,
                                    minSuccDiff = 10^(-5),
                                    maxIter = 30){

  ## ---------------------------------------------------------------------
  ## Function finds the optimal penalty parameter of the ridge
  ## regression estimator of the generalized linear model parameter. The
  ## optimum is defined as the minimizer of the cross-validated loss
  ## associated with the estimator.
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y           : response vector
  ## X           : design matrix
  ## lambdaInit  : initial guess of the optimal penalty parameter
  ## folds       : list with folds
  ## stratified  : boolean, should the data by stratified such that
  ##               range of the response is comparable among folds

  ## model       : A character indicating which generalized linear model
  ##               model instance is to be fitted.
  ## loss        : loss to be used in CV, either "sos" (sum-of-squares) or
  ##              "loglik" (loglikelihood)
  ## minSuccDiff : A numeric, the minimum distance between two successive
  ##               estimates to be achieved.
  ## maxIter     : A numeric specifying the maximum number of iterations.
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the optimal cross-validated penalty parameter of the
  ## ridge estimator of the generalized linear regression model parameter
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen, modified and extended by JM Goedhart
  ## ---------------------------------------------------------------------

  ####----------------------------------------------------------------###
  ####----------------------------------------------------------------###
  if (ncol(X) == 0 | nrow(X) == 0) {stop("Omics matrix has zero columns or rows")}
  if (ncol(U) == 0 | nrow(U) == 0) {stop("Clinical matrix has zero columns or rows")}
  if (length(Y) == 0){stop("Length Y vector is zero")}
  if (is.null(Y) | is.null(X) | is.null(U)) {stop("Either Y, X, or Z is unspecified")}

  if (lambdaGinit == 0){
    message("Tuning fusedTree without fusion penalty, as alphaInit = 0 is specified")


    if (model == "linear"){
      if (ncol(X) >= nrow(X)){
        lambdasOpt <- stats::optim(log(lambdaInit),
                            .kcvLMloss_PgeqN_UP_noGP,
                            Y = Y,
                            X = X,
                            U = U,
                            folds = folds,
                            loss = loss,
                            method = "Brent",
                            lower = log(0.1), upper = log(10^10),
                            control = list(reltol = 1e-5, trace = 1, parscale =log(2))
        )$par
        return(exp(lambdasOpt))
      }
      if (ncol(X) < nrow(X)){
        lambdasOpt <- stats::optim(log(lambdaInit),
                            .kcvLMloss_PleqN_UP_GPandNoGP,
                            Y = Y,
                            X = as.matrix(X),
                            U = U,
                            folds = folds,
                            loss = loss,
                            method = "Brent",
                            lower = log(0.1), upper = log(10^10),
                            control = list(reltol = 1e-5, trace = 1, parscale =log(2))
        )$par
        return(exp(lambdasOpt))
      }
    }

    if (model == "logistic"){

      if(loss != "loglik"){stop("For logistic model, only loglik loss is implemented")}

      lambdasOpt <- stats::optim(log(lambdaInit),
                          .kcvBLMloss_UP_noGP,
                          Y = Y,
                          X = X,
                          U = U,
                          folds = folds,
                          maxIter = maxIter,
                          minSuccDiff = minSuccDiff,
                          method = "Brent",
                          lower = log(0.1), upper = log(10^10),
                          control = list(reltol = 1e-5, trace = 1, parscale =log(2))
      )$par
      return(exp(lambdasOpt))
    }
    if (model == "cox"){
      if(loss != "loglik"){stop("For survival model, only loglik loss is implemented")}

      lambdasOpt <- stats::optim(log(lambdaInit),
                          .kcvSurvloss_UP_noGP,
                          Y = Y,
                          X = X,
                          U = U,
                          folds = folds,
                          maxIter = maxIter,
                          minSuccDiff = minSuccDiff,
                          method = "Brent",
                          lower = log(0.1), upper = log(10^10),
                          control = list(reltol=1e-5, trace = 1, parscale =log(2))
      )$par
      return(exp(lambdasOpt))
    }
  }

  if (lambdaGinit > 0) {
    message("Tuning fusedTree with fusion penalty")

    if (model == "linear") {
      if (ncol(X) >= nrow(X)){

        NumNodes = length(which(Dg[1,]!=0))
        Dg <- .EigPenmatr(NumNodes = NumNodes, p = ncol(X)/NumNodes,
                          symmetric = symFusion,
                          Tree = Tree)  #fast eigendecomposition

        if (isTRUE(multistart)) {
          starts <- list(
            log(c(lambdaInit, lambdaGinit)),
            log(0.1) * log(c(lambdaInit, lambdaGinit)),
            log(100) * log(c(lambdaInit, lambdaGinit)),
            log(c(lambdaInit, 10000*lambdaGinit)),
            log(c(1000*lambdaInit, lambdaGinit))
          )
          results <- lapply(starts, function(start) {
            stats::optim(start,
                  .kcvLMloss_PgeqN_UP_GP,
                  Y = Y,
                  X = X %*% Dg$vectors,
                  U = U,
                  Ds = Dg$values,
                  folds = folds,
                  loss = loss,
                  method = "Nelder-Mead",
                  control = list(reltol = 1e-5, maxit = 1000,
                                 parscale = c(log(2), log(2)))
            )
          })

          # Choose the best result
          best_res <- results[[which.min(sapply(results, `[[`, "value"))]]
          lambdasOpt <- exp(best_res$par)
          return(lambdasOpt)
        }
        if (isFALSE(multistart)) {
          lambdasOpt <- stats::optim(log(c(lambdaInit, lambdaGinit)),
                              .kcvLMloss_PgeqN_UP_GP,
                              Y = Y,
                              X = X %*% Dg$vectors,
                              U = U,
                              Ds = Dg$values,
                              folds = folds,
                              loss = loss,
                              method = "Nelder-Mead",
                              control = list(reltol = 1e-5,
                                             parscale = c(log(2), log(2)))
          )$par
          return(exp(lambdasOpt))
        }

      }
      if (ncol(X) < nrow(X)) {

        if (isTRUE(multistart)) {
          starts <- list(
            log(c(lambdaInit, lambdaGinit)),
            log(0.1) * log(c(lambdaInit, lambdaGinit)),
            log(1000) * log(c(lambdaInit, lambdaGinit)),
            log(c(lambdaInit, 10000*lambdaGinit)),
            log(c(1000*lambdaInit, lambdaGinit))
          )
          results <- lapply(starts, function(start) {
            stats::optim(start,
                  .kcvLMloss_PleqN_UP_GPandNoGP,
                  Y = Y,
                  X = as.matrix(X),
                  U = U,
                  Dg = as.matrix(Dg),
                  folds = folds,
                  loss = loss,
                  method = "Nelder-Mead",
                  control = list(reltol=1e-5,
                                 parscale = c(log(2), log(2)))
            )
          })
          best_res <- results[[which.min(sapply(results, `[[`, "value"))]]
          lambdasOpt <- exp(best_res$par)
          return(lambdasOpt)
        }

        if (isFALSE(multistart)) {
          lambdasOpt <- stats::optim(log(c(lambdaInit, lambdaGinit)),
                              .kcvLMloss_PleqN_UP_GPandNoGP,
                              Y = Y,
                              X = as.matrix(X),
                              U = U,
                              Dg = as.matrix(Dg),
                              folds = folds,
                              loss = loss,
                              method = "Nelder-Mead",
                              control = list(reltol=1e-5, parscale = c(log(2),log(2)))
          )$par
          return(exp(lambdasOpt))
        }
      }
    }

    if (model == "logistic"){
      if(loss != "loglik"){stop("For logistic model, only loglik loss is implemented")}

      if (isTRUE(multistart)) {
        starts <- list(
          log(c(lambdaInit, lambdaGinit)),
          log(0.1) * log(c(lambdaInit, lambdaGinit)),
          log(1000) * log(c(lambdaInit, lambdaGinit)),
          log(c(lambdaInit, 10000*lambdaGinit)),
          log(c(1000*lambdaInit, lambdaGinit))
        )
        results <- lapply(starts, function(start) {
          stats::optim(start,
                .kcvBLMloss_UP_GP,
                Y = Y,
                X = X,
                U = U,
                Dg = Dg,
                folds = folds,
                symFusion = symFusion,
                Tree = Tree,
                maxIter = maxIter,
                minSuccDiff = minSuccDiff,
                method = "Nelder-Mead",
                control = list(reltol = 1e-5,
                               parscale = c(log(2), log(2)))
          )
        })
        best_res <- results[[which.min(sapply(results, `[[`, "value"))]]
        lambdasOpt <- exp(best_res$par)
        return(lambdasOpt)
      }

      if (isFALSE(multistart)) {
        lambdasOpt <- stats::optim(log(c(lambdaInit, lambdaGinit)),
                            .kcvBLMloss_UP_GP,
                            Y = Y,
                            X = X,
                            U = U,
                            Dg = Dg,
                            folds = folds,
                            symFusion = symFusion,
                            Tree = Tree,
                            maxIter = maxIter,
                            minSuccDiff=minSuccDiff,
                            method = "Nelder-Mead",
                            control = list(reltol = 1e-5,
                                           parscale = c(log(2), log(2)))
        )$par
        return(exp(lambdasOpt))
      }
    }

    if (model == "cox"){
      if(loss != "loglik"){stop("For survival model, only loglik loss is implemented")}

      if (isTRUE(multistart)) {
        starts <- list(
          log(c(lambdaInit, lambdaGinit)),
          log(0.1) * log(c(lambdaInit, lambdaGinit)),
          log(1000) * log(c(lambdaInit, lambdaGinit)),
          log(c(lambdaInit, 10000*lambdaGinit)),
          log(c(1000*lambdaInit, lambdaGinit))
        )

        results <- lapply(starts, function(start) {
          stats::optim(start,
                .kcvSurvloss_UP_GP,
                Y = Y,
                X = X,
                U = U,
                Dg = Dg,
                folds = folds,
                symFusion = symFusion,
                Tree = Tree,
                maxIter = maxIter,
                minSuccDiff = minSuccDiff,
                method = "Nelder-Mead",
                control = list(reltol = 1e-5,
                               parscale = c(log(2), log(2)))
          )
        })
        best_res <- results[[which.min(sapply(results, `[[`, "value"))]]
        lambdasOpt <- exp(best_res$par)
        return(lambdasOpt)
      }

      if (isFALSE(multistart)) {
        lambdasOpt <- stats::optim(log(c(lambdaInit, lambdaGinit)),
                            .kcvSurvloss_UP_GP,
                            Y = Y,
                            X = X,
                            U = U,
                            Dg = Dg,
                            folds = folds,
                            symFusion = symFusion,
                            Tree = Tree,
                            maxIter = maxIter,
                            minSuccDiff = minSuccDiff,
                            method = "Nelder-Mead",
                            control = list(reltol = 1e-5,
                                           parscale = c(log(2), log(2)))
        )$par
        return(exp(lambdasOpt))
      }
    }
  }
}

.kcvLMloss_PgeqN_UP_noGP <- function(lambda, Y, X, U, folds, loss){

  ## ---------------------------------------------------------------------
  ## Internal function yields the cross-validated loglikelihood of the
  ## ridge regression estimator with a two-dimensional covariate layout.
  ## ---------------------------------------------------------------------
  ## Arguments
  ## lambdas : penalty parameter vector
  ## Y       : response vector
  ## X       : design matrix multiplied by the eigenvector matrix of the
  ##           nonnegative definite matrix that specifies the structure
  ##           of spatial fused ridge penalty.
  ## U       : design matrix of the unpenalized covariates

  ## Ds      : nonnegative eigenvalues of the nonnegative definite matrix
  ##           that specifies the structure of spatial fused ridge penalty.
  ## folds   : list with folds
  ## loss    : character, either 'loglik' of 'sos', specifying the loss
  ##           criterion to be used in the cross-validation
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the linear regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  # evaluate loss per fold
  cvLoss <- 0
  lambda <- exp(lambda) # positivity constraint

  # calculation prior to cv-loop
  XXT     <- Matrix::tcrossprod(X) / lambda


  for (k in 1:length(folds)){
    # evaluate the estimator of the unpenalized low-dimensional
    # regression parameter on the training samples
    if (all(dim(U) > 0)){
      gHat <- Matrix::solve(Matrix::crossprod(U[-folds[[k]],,drop=FALSE],
                            Matrix::solve(XXT[-folds[[k]], -folds[[k]]] +
                            Matrix::diag(nrow(XXT)-length(folds[[k]])),
                            U[-folds[[k]],,drop=FALSE])),
                            Matrix::crossprod(U[-folds[[k]],,drop=FALSE],
                            Matrix::solve(XXT[-folds[[k]], -folds[[k]]] +
                            Matrix::diag(nrow(XXT)-length(folds[[k]])),
                            Y[-folds[[k]]])))
      gHat <- as.numeric(gHat)
    }

    # evaluate the linear predictor on the left-out samples
    if (all(dim(U) > 0)){
      Yhat <- Matrix::tcrossprod(Matrix::solve(Matrix::diag(rep(1, nrow(XXT)-length(folds[[k]]))) +
                              XXT[-folds[[k]], -folds[[k]]],
                              Y[-folds[[k]]] -
                              as.numeric(Matrix::crossprod(Matrix::t(U[-folds[[k]],,drop=FALSE]), gHat))),
                        XXT[folds[[k]], -folds[[k]], drop=FALSE])
    } else {
      Yhat <- Matrix::tcrossprod(Matrix::solve(Matrix::diag(rep(1, nrow(XXT)-length(folds[[k]]))) +
                                XXT[-folds[[k]], -folds[[k]]],
                              Y[-folds[[k]]]),
                        XXT[folds[[k]], -folds[[k]], drop=FALSE])
    }

    if (loss == "loglik"){
      # evaluate the loglikelihood on the left-out samples
      cvLoss <- cvLoss + length(Y[folds[[k]]]) *
        (log(2 * pi * sum((Y[folds[[k]]] - Yhat)^2) /
               length(folds[[k]])) + 1) / 2
    }
    if (loss == "sos"){
      # evaluate the sum-of-sqaures on the left-out samples
      cvLoss <- cvLoss + sum((Y[folds[[k]]] - Yhat)^2)
    }
  }

  # average over the folds
  return(cvLoss / length(folds))
}

.kcvLMloss_PleqN_UP_GPandNoGP <- function(lambdas, Y, X, U,
                                          Dg = matrix(0, nrow = nrow(X), ncol = ncol(X)),
                                          folds, loss){

  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the linear regression
  ## model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds

  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  # evaluate loss per fold


  cvLoss <- 0
  lambda  <- exp(lambdas[1]) # positivity constraint
  lambdaG <- exp(lambdas[2]) # positivity constraint

  if (is.na(lambdaG)) {
    lambdaG <- 0
  }


  for (k in 1:length(folds)){
    # evaluate regression estimate
    bHat <- suppressMessages(.ridgeGLM2(Y[-folds[[k]]], X[-folds[[k]],,drop=FALSE],
                       U[-folds[[k]],,drop=FALSE], Dg=Dg,
                       lambda =lambda, lambdaG = lambdaG,
                       model="linear", verbose = FALSE))

    # evaluate linear predictor
    Yhat <- X[folds[[k]],,drop=FALSE] %*% bHat[-c(1:ncol(U))] +
      U[folds[[k]],,drop=FALSE] %*% bHat[c(1:ncol(U))]

    if (loss == "loglik"){
      # evaluate the loglikelihood on the left-out samples
      cvLoss <- cvLoss + length(Y[folds[[k]]]) *
        (log(2 * pi * sum((Y[folds[[k]]] - as.numeric(Yhat))^2) /
               length(folds[[k]])) + 1) / 2

    }
    if (loss == "sos"){
      # evaluate the sum-of-sqaures on the left-out samples
      cvLoss <- cvLoss + sum((Y[folds[[k]]] - as.numeric(Yhat))^2)
    }
  }

  # average over the folds
  return(cvLoss / length(folds))
}

.kcvBLMloss_UP_noGP <- function(lambda, Y, X, U, folds, minSuccDiff, maxIter){

  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the logistic regression
  ## model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds

  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen, modified by JM Goedhart
  ## ---------------------------------------------------------------------

  # initiate
  cvLoss  <- 0
  lambda <- exp(lambda) # positivity constraint

  # evaluate subexpressions of the estimator
  if (ncol(X) >= nrow(X)){
    XXT     <- Matrix::tcrossprod(X)/lambda
  }


  for (k in 1:length(folds)){


    # evaluate the initial linear predictor
    lpOld <- rep(log(mean(Y[-folds[[k]]])/(1-mean(Y[-folds[[k]]]))), length(Y[-folds[[k]]]))
    lpPrev <- rep(log(mean(Y)/(1-mean(Y))), length(Y))

    loglikPrev <- .loglikBLMlp(Y=Y[-folds[[k]]], lp = lpOld) #initial loglikelihood (penalty is zero)

    for (iter in 1:maxIter){
      # calculate the weights
      Ypred <- 1 / (1 + exp(-lpOld))
      W0    <- Ypred * (1 - Ypred)
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <-
          .Machine$double.eps
      }
      if (ncol(X) >= nrow(X)){
        # adjusted response
        Z <- lpOld + (Y[-folds[[k]]] - Ypred)/W0

        # evaluate unpenalized low-dim regression estimator
        Matrix::diag(XXT)[-folds[[k]]] <- Matrix::diag(XXT)[-folds[[k]]] + 1/W0
        slh        <- Matrix::solve(XXT[-folds[[k]], -folds[[k]]],
                            Z)
        gHat       <- Matrix::solve(Matrix::crossprod(U[-folds[[k]],],
                                    Matrix::solve(XXT[-folds[[k]], -folds[[k]]],
                                            U[-folds[[k]],])),
                                    Matrix::crossprod(U[-folds[[k]],], slh))
        Ug         <- Matrix::tcrossprod(U, Matrix::t(gHat))

        # evaluate linear predictor without evaluating the estimator
        slh        <- Matrix::solve(XXT[-folds[[k]], -folds[[k]], drop=FALSE],
                            Z - Ug[-folds[[k]],])
        Matrix::diag(XXT)[-folds[[k]]] <- Matrix::diag(XXT)[-folds[[k]]] - 1/W0
        lpAll      <- Matrix::crossprod(XXT[-folds[[k]], ,drop=FALSE], slh)
        penalty    <- sum(lpAll[-folds[[k]]] * slh) / 2
        lpAll      <- as.numeric(lpAll) + as.numeric(Ug)
      }
      if (ncol(X) < nrow(X)){
        # adjusted response
        X <- as.matrix(X)
        Z <- W0*lpOld + (Y[-folds[[k]]] - Ypred)

        # evaluate subexpressions of the estimator
        XTXpD <-crossprod(sweep(X[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*"))
        tUTX  <-crossprod(X[-folds[[k]], , drop=FALSE],
                          sweep(U[-folds[[k]], , drop=FALSE], 1, W0, FUN="*"))
        XTZ   <-crossprod(X[-folds[[k]], , drop=FALSE], Z)
        UTZ   <-crossprod(U[-folds[[k]], , drop=FALSE], Z)
        UTU   <-crossprod(sweep(U[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*"))
        diag(XTXpD) <- diag(XTXpD) + lambda

        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU -crossprod(solve(XTXpD, tUTX), tUTX),
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ),
                                            crossprod(X[-folds[[k]], , drop=FALSE],
                                                      sweep(U[-folds[[k]], , drop=FALSE],
                                                            1, W0, FUN="*")))))))
        Ug <- as.numeric(tcrossprod(gHat, U))

        # evaluate linear predictor
        XWZpT <- as.numeric(
          crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld +
                      (Y[-folds[[k]]] - Ypred) - W0*Ug[-folds[[k]]]))
        bHat    <- qr.solve(XTXpD, XWZpT,
                            tol=.Machine$double.eps)
        lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug
        penalty <- 0.5 * lambda * sum((bHat)^2)
      }

      # split linear predictor by fold
      lpOld      <- lpAll[-folds[[k]]]
      lpNew      <- lpAll[ folds[[k]]]


      # compute new likelihood
      loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld) - penalty

      # step-halving for stability
      if (loglik < loglikPrev){
        lpOld <- 0.5*lpOld + 0.5*lpPrev[-folds[[k]]]
        lpNew <- 0.5*lpNew + 0.5*lpPrev[folds[[k]]]
      }

      # assess convergence
      if (is.infinite(loglik) | is.na(loglik) | is.nan(loglik)){
        lpNew <- rep(0,nrow(Y[folds[[k]]]))
        break
      } else if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lpAll
      }
    }
    cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
  }
  # average over the folds
  return(-cvLoss / length(folds))
}

.kcvSurvloss_UP_noGP <- function(lambda, Y, X, U, folds, minSuccDiff, maxIter){

  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the logistic regression
  ## model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds

  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  # initiate
  RespTot <- c()
  LPTot <- c()
  lambda <- exp(lambda) # to ensure positive lambda

  # evaluate subexpressions of the estimator
  if (ncol(X) >= nrow(X)){
    XXT     <- Matrix::tcrossprod(X)/lambda
  }


  for (k in 1:length(folds)){

    # initialization
    penalty <- 0
    lpPrev <- rep(0, length(Y[,2]))
    lpOld <- rep(0, length(Y[-folds[[k]],2]))
    loglikPrev <- .loglikSurv(survival::Surv(time=Y[-folds[[k]],1],event = Y[-folds[[k]],2]),lpOld) - penalty
    Yev <- Y[-folds[[k]],2]

    for (iter in 1:maxIter){
      # calculate the weights
      H0 <- .breslow(survival::Surv(time=Y[-folds[[k]],1],event = Y[-folds[[k]],2]),lpOld)[,2]
      Ypred <- as.numeric(H0 * exp(lpOld))
      W0 <- Ypred
      if (min(W0) <= 2*.Machine$double.eps){
        W0[which(W0 < 2*.Machine$double.eps)] <-
          2*.Machine$double.eps
      }

      if (ncol(X) >= nrow(X)){
        # adjusted response
        Z <- lpOld + (Yev - Ypred)/W0

        # evaluate unpenalized low-dim regression estimator
        Matrix::diag(XXT)[-folds[[k]]] <- Matrix::diag(XXT)[-folds[[k]]] + 1/W0
        slh        <- Matrix::solve(XXT[-folds[[k]], -folds[[k]]],
                            Z)
        gHat       <- Matrix::solve(Matrix::crossprod(U[-folds[[k]],],
                                    Matrix::solve(XXT[-folds[[k]], -folds[[k]]],
                                            U[-folds[[k]],])),
                              Matrix::crossprod(U[-folds[[k]],], slh))
        Ug         <- Matrix::tcrossprod(U, Matrix::t(gHat))

        # evaluate linear predictor without evaluating the estimator
        slh        <- Matrix::solve(XXT[-folds[[k]], -folds[[k]], drop=FALSE],
                            Z - Ug[-folds[[k]],])
        Matrix::diag(XXT)[-folds[[k]]] <- Matrix::diag(XXT)[-folds[[k]]] - 1/W0
        lpAll      <- Matrix::crossprod(XXT[-folds[[k]], ,drop=FALSE], slh)
        penalty    <- sum(lpAll[-folds[[k]]] * slh) / 2
        lpAll      <- as.numeric(lpAll) + as.numeric(Ug)
      }
      if (ncol(X) < nrow(X)){
        # adjusted response
        X <- as.matrix(X)
        Z <- W0*lpOld + (Yev - Ypred)

        # evaluate subexpressions of the estimator
        XTXpD <-crossprod(sweep(X[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*"))
        tUTX  <-crossprod(X[-folds[[k]], , drop=FALSE],
                          sweep(U[-folds[[k]], , drop=FALSE], 1, W0, FUN="*"))
        XTZ   <-crossprod(X[-folds[[k]], , drop=FALSE], Z)
        UTZ   <-crossprod(U[-folds[[k]], , drop=FALSE], Z)
        UTU   <-crossprod(sweep(U[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*"))
        diag(XTXpD) <- diag(XTXpD) + lambda

        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU -crossprod(solve(XTXpD, tUTX), tUTX),
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ),
                                            crossprod(X[-folds[[k]], , drop=FALSE],
                                                      sweep(U[-folds[[k]], , drop=FALSE],
                                                            1, W0, FUN="*")))))))
        Ug <- as.numeric(tcrossprod(gHat, U))

        # evaluate linear predictor
        XWZpT <- as.numeric(
          crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld +
                      (Yev - Ypred) - W0*Ug[-folds[[k]]]))
        bHat <- solve(XTXpD,XWZpT)
        #bHat    <- qr.solve(XTXpD, XWZpT,
        #                    tol=.Machine$double.eps)
        lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug
        penalty <- 0.5 * lambda * sum((bHat)^2)
      }
      # split linear predictor by fold
      lpOld <- lpAll[-folds[[k]]]
      lpNew <- lpAll[folds[[k]]]

      # compute new likelihood
      loglik <- .loglikSurv(survival::Surv(time=Y[-folds[[k]],1],event = Y[-folds[[k]],2]), lpOld) - penalty

      # step-halving for stability
      if (loglik < loglikPrev){
        lpOld <- 0.5*lpOld + 0.5*lpPrev[-folds[[k]]]
        lpNew <- 0.5*lpNew + 0.5*lpPrev[folds[[k]]]
      }

      # assess convergence
      if (is.infinite(loglik) | is.na(loglik) | is.nan(loglik)){
        lpNew=rep(0,nrow(Y[folds[[k]],]))
        break
      } else if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        break
      } else {
        loglikPrev <- loglik
        lpPrev<-lpAll
      }

    }

    RespTot <- base::append(RespTot,Y[folds[[k]],])
    LPTot <- append(LPTot,as.numeric(lpNew))
  }

  RespTot1 <- survival::Surv(time=RespTot[,1], event = RespTot[,2])
  CVLoss = .loglikSurv(Y=RespTot1,lp=LPTot)

  # average over the folds
  return(-CVLoss/length(folds))
}

.kcvLMloss_PgeqN_UP_GP <- function(lambdas, Y, X, U, Ds, folds, loss){

  ## ---------------------------------------------------------------------
  ## Internal function yields the cross-validated loglikelihood of the
  ## ridge regression estimator with a two-dimensional covariate layout.
  ## ---------------------------------------------------------------------
  ## Arguments
  ## lambdas : penalty parameter vector
  ## Y       : response vector
  ## X       : design matrix multiplied by the eigenvector matrix of the
  ##           nonnegative definite matrix that specifies the structure
  ##           of spatial fused ridge penalty.
  ## U       : design matrix of the unpenalized covariates
  ##
  ## Ds      : nonnegative eigenvalues of the nonnegative definite matrix
  ##           that specifies the structure of spatial fused ridge penalty.
  ## folds   : list with folds
  ## loss    : character, either 'loglik' of 'sos', specifying the loss
  ##           criterion to be used in the cross-validation
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the linear regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  # evaluate loss per fold
  cvLoss <- 0
  lambda  <- exp(lambdas[1]) # positivity constraint
  lambdaG <- exp(lambdas[2]) # positivity constraint

  # calculation prior to cv-loop
  #X <- sweep(X, 2, sqrt(lambda + lambdaG * Ds), FUN="/")
  D <- Matrix::Diagonal(x = 1 / sqrt(lambda + lambdaG * Ds))
  X <- X %*% D
  XXT     <- Matrix::tcrossprod(X)


  for (k in 1:length(folds)){
    # evaluate the estimator of the unpenalized low-dimensional
    # regression parameter on the training samples
    if (all(dim(U) > 0)){
      gHat <- Matrix::solve(Matrix::crossprod(U[-folds[[k]],,drop=FALSE],
                            Matrix::solve(XXT[-folds[[k]], -folds[[k]]] +
                            Matrix::diag(nrow(XXT)-length(folds[[k]])),
                            U[-folds[[k]],,drop=FALSE])),
                            Matrix::crossprod(U[-folds[[k]],,drop=FALSE],
                            Matrix::solve(XXT[-folds[[k]], -folds[[k]]] +
                            Matrix::diag(nrow(XXT)-length(folds[[k]])),
                                    Y[-folds[[k]]])))
      gHat <- as.numeric(gHat)
    }
    if (all(dim(U) > 0)){
      Yhat   <- as.numeric(Matrix::crossprod(Matrix::t(XXT[folds[[k]], -folds[[k]]]),
                           Matrix::solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) +
                                             XXT[-folds[[k]], -folds[[k]]],
                          Y[-folds[[k]]] - as.numeric(Matrix::crossprod(Matrix::t(U[-folds[[k]],,drop=FALSE]), gHat)))))
      + as.numeric(Matrix::crossprod(Matrix::t(U[-folds[[k]],,drop=FALSE]),gHat))

    } else {
      Yhat   <-
        Matrix::tcrossprod(Matrix::solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) +
                           XXT[-folds[[k]], -folds[[k]]],
                         Y[-folds[[k]]]),
                   XXT[folds[[k]], -folds[[k]], drop=FALSE])
    }


    if (loss == "loglik"){
      # evaluate the loglikelihood on the left-out samples
      cvLoss <- cvLoss + length(Y[folds[[k]]]) *
        (log(2 * pi * sum((Y[folds[[k]]] - Yhat)^2) /
               length(folds[[k]])) + 1) / 2
    }
    if (loss == "sos"){
      # evaluate the sum-of-sqaures on the left-out samples
      cvLoss <- cvLoss + sum((Y[folds[[k]]] - Yhat)^2)
    }
  }

  # average over the folds
  return(cvLoss / length(folds))
}

.kcvBLMloss_UP_GP <- function(lambdas, Y, X, U, Dg, folds,
                              minSuccDiff, maxIter, symFusion, Tree){

  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the logistic regression
  ## model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds
  ##
  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen, modified by JM Goedhart
  ## ---------------------------------------------------------------------

  # initiate
  cvLoss  <- 0
  lambda  <- exp(lambdas[1]) # positivity constraint
  lambdaG <- exp(lambdas[2]) # positivity constraint

  # evaluate subexpressions of the estimator
  if (ncol(X) >= nrow(X)){
    NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
    Dg <-Matrix::tcrossprod(Matrix::kronecker(Matrix::diag(1,pX),Matrix::solve(lambdaG*.PenMatr(NumNodes,1, symmetric = symFusion, Tree = Tree) +
                            Matrix::diag(lambda,NumNodes)))
                    ,X)

    XDXT        <- Matrix::crossprod(Matrix::t(X), Dg)
    diagXDXTorg <- Matrix::diag(XDXT)
  }

  if (ncol(X) < nrow(X)){
    Dg <- as.matrix(Dg)
    Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
  }

  for (k in 1:length(folds)){
    # initialization
    penalty    <-  0
    lpOld <- rep(log(mean(Y[-folds[[k]]])/(1-mean(Y[-folds[[k]]]))), length(Y[-folds[[k]]]))
    lpPrev <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    loglikPrev <- .loglikBLMlp(Y = Y[-folds[[k]]], lp = lpOld) #initial likelihood, penalty equals zero

    for (iter in 1:maxIter){
      # calculate the weights
      Ypred <- 1 / (1 + exp(-lpOld))
      W0    <- Ypred * (1 - Ypred)
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <-
          .Machine$double.eps
      }
      if (ncol(X) >= nrow(X)){
        # adjusted response
        Z <- lpOld + (Y[-folds[[k]]] - Ypred)/W0

        # evaluate unpenalized low-dim regression estimator
        Matrix::diag(XDXT)[-folds[[k]]] <- Matrix::diag(XDXT)[-folds[[k]]] + 1/W0
        slh        <- Matrix::solve(XDXT[-folds[[k]], -folds[[k]]], Z)

        gHat       <- Matrix::solve(Matrix::crossprod(U[-folds[[k]],],
                                    Matrix::solve(XDXT[-folds[[k]], -folds[[k]]],
                                            U[-folds[[k]],])),
                                    Matrix::crossprod(U[-folds[[k]],], slh))
        Ug         <- Matrix::tcrossprod(U, Matrix::t(gHat))

        # evaluate linear predictor without evaluating the estimator
        slh        <- Matrix::solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE],
                            Z - Ug[-folds[[k]],])
        Matrix::diag(XDXT)[-folds[[k]]] <- diagXDXTorg[-folds[[k]]]
        lpAll      <- Matrix::crossprod(XDXT[-folds[[k]], ,drop=FALSE], slh)
        penalty    <- 0.5 * sum(lpAll[-folds[[k]]] * slh)
        lpAll      <- as.numeric(lpAll) + as.numeric(Ug)
      }

      if (ncol(X) < nrow(X)){
        # adjusted response

        Z <- W0*lpOld + (Y[-folds[[k]]] - Ypred)
        X <- as.matrix(X)


        # evaluate subexpressions of the estimator
        XTXpD <-crossprod(sweep(X[-folds[[k]], , drop=FALSE],
                                1, sqrt(W0), FUN="*"))
        tUTX  <-crossprod(X[-folds[[k]], , drop=FALSE],
                          sweep(U[-folds[[k]], , drop=FALSE],
                                1, W0, FUN="*"))
        XTZ   <-crossprod(X[-folds[[k]], , drop=FALSE], Z)
        UTZ   <-crossprod(U[-folds[[k]], , drop=FALSE], Z)
        UTU   <-crossprod(sweep(U[-folds[[k]], , drop=FALSE],
                                1, sqrt(W0), FUN="*"))
        XTXpD <- XTXpD + Dg

        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU -crossprod(solve(XTXpD, tUTX), tUTX),
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ),
                                            crossprod(X[-folds[[k]], , drop=FALSE],
                                                      sweep(U[-folds[[k]], , drop=FALSE],
                                                            1, W0, FUN="*")))))))
        Ug   <- as.numeric(tcrossprod(gHat, U))

        # evaluate linear predictor
        XWZpT <- as.numeric(
          crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld +
                      (Y[-folds[[k]]] - Ypred) - W0*Ug[-folds[[k]]]))
        bHat    <- qr.solve(XTXpD, XWZpT,
                            tol=.Machine$double.eps)
        lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug
        penalty <- 0.5 * sum(crossprod(Dg, bHat) *
                               (bHat))
      }

      # split linear predictor by fold
      lpOld      <- lpAll[-folds[[k]]]
      lpNew      <- lpAll[ folds[[k]]]

      # compute new likelihood
      loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld) - penalty

      # step-halving for stability
      if (loglik < loglikPrev){
        lpOld <- 0.5*lpOld + 0.5*lpPrev[-folds[[k]]]
        lpNew <- 0.5*lpNew + 0.5*lpPrev[folds[[k]]]
      }

      # assess convergence
      if (is.infinite(loglik) | is.na(loglik) | is.nan(loglik)){
        lpNew=rep(0,nrow(Y[folds[[k]]]))
        break
      } else if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lpAll
      }
    }
    cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
  }
  # average over the folds
  return(-cvLoss / length(folds))
}

.kcvSurvloss_UP_GP <- function(lambdas, Y, X, U, Dg, folds,
                               minSuccDiff, maxIter, symFusion, Tree){

  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the logistic regression
  ## model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds
  ##
  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  # initiate
  RespTot <- c()
  LPTot <- c()

  lambda  <- exp(lambdas[1]) #positivity constraint
  lambdaG <- exp(lambdas[2]) #positivity constraint

  # evaluate subexpressions of the estimator
  if (ncol(X) >= nrow(X)){
    NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
    Dg <- Matrix::tcrossprod(Matrix::kronecker(Matrix::diag(1,pX),Matrix::solve(lambdaG*.PenMatr(NumNodes,1, symmetric = symFusion, Tree = Tree)+
                                                    Matrix::diag(lambda,NumNodes)))
                    ,X)

    XDXT        <- Matrix::crossprod(Matrix::t(X), Dg)
    diagXDXTorg <- Matrix::diag(XDXT)
  }
  if (ncol(X) < nrow(X)){
    Dg <- as.matrix(Dg)
    Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * as.matrix(Dg)
  }


  for (k in 1:length(folds)){

    # initialization
    penalty    <-  0
    lpPrev <- rep(0, length(Y[,2]))
    lpOld <- rep(0, length(Y[-folds[[k]],2]))

    loglikPrev <- .loglikSurv(survival::Surv(time=Y[-folds[[k]],1],event = Y[-folds[[k]],2]), lp=lpOld)-penalty
    Yev = Y[-folds[[k]],2]


    for (iter in 1:maxIter){

      # calculate the weights
      H0 <- .breslow(survival::Surv(time = Y[-folds[[k]],1], event = Y[-folds[[k]],2]),lpOld)[,2]
      Ypred <- as.numeric(H0 * exp(lpOld)) + 10^(-15)
      W0 <- Ypred
      if (min(W0) <= 2*.Machine$double.eps){
        W0[which(W0 < 2*.Machine$double.eps)] <-
          2*.Machine$double.eps
      }

      if (ncol(X) >= nrow(X)){
        # adjusted response
        Z <- lpOld + (Yev - Ypred)/W0

        # evaluate unpenalized low-dim regression estimator
        Matrix::diag(XDXT)[-folds[[k]]] <- Matrix::diag(XDXT)[-folds[[k]]] + 1/W0
        slh        <- Matrix::solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE], Z)

        gHat       <- Matrix::solve(Matrix::crossprod(U[-folds[[k]], , drop=FALSE],
                                    Matrix::solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE],
                                            U[-folds[[k]], , drop=FALSE])),
                                    Matrix::crossprod(U[-folds[[k]], , drop=FALSE], slh))
        Ug         <- Matrix::tcrossprod(U, Matrix::t(gHat))

        # evaluate linear predictor without evaluating the estimator
        slh        <- Matrix::solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE],
                            Z - Ug[-folds[[k]],])
        Matrix::diag(XDXT)[-folds[[k]]] <- diagXDXTorg[-folds[[k]]]
        lpAll      <- Matrix::crossprod(XDXT[-folds[[k]], ,drop=FALSE], slh)
        penalty    <- 0.5 * sum(lpAll[-folds[[k]]] * slh)
        lpAll      <- as.numeric(lpAll) + as.numeric(Ug)
      }

      if (ncol(X) < nrow(X)){
        # adjusted response
        X <- as.matrix(X)
        Z <- W0*lpOld + (Yev - Ypred)


        # evaluate subexpressions of the estimator
        XTXpD <-crossprod(sweep(X[-folds[[k]], , drop=FALSE],
                                1, sqrt(W0), FUN="*"))
        tUTX  <-crossprod(X[-folds[[k]], , drop=FALSE],
                          sweep(U[-folds[[k]], , drop=FALSE],
                                1, W0, FUN="*"))
        XTZ   <-crossprod(X[-folds[[k]], , drop=FALSE], Z)
        UTZ   <-crossprod(U[-folds[[k]], , drop=FALSE], Z)
        UTU   <-crossprod(sweep(U[-folds[[k]], , drop=FALSE],
                                1, sqrt(W0), FUN="*"))
        XTXpD <- XTXpD + Dg

        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU -crossprod(solve(XTXpD, tUTX), tUTX),
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ),
                                            crossprod(X[-folds[[k]], , drop=FALSE],
                                                      sweep(U[-folds[[k]], , drop=FALSE],
                                                            1, W0, FUN="*")))))))
        Ug   <- as.numeric(tcrossprod(gHat, U))

        # evaluate linear predictor
        XWZpT <- as.numeric(
          crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld +
                      (Yev - Ypred) - W0*Ug[-folds[[k]]]))
        bHat    <- qr.solve(XTXpD, XWZpT,
                            tol=.Machine$double.eps)
        lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug
        penalty <- 0.5 * sum(crossprod(Dg, bHat) *
                               (bHat))
      }

      # split linear predictor by fold
      lpOld      <- lpAll[-folds[[k]]]
      lpNew      <- lpAll[ folds[[k]]]

      #compute likelihood
      loglik <- .loglikSurv(survival::Surv(time=Y[-folds[[k]],1],event = Y[-folds[[k]],2]), lpOld) - penalty

      # step-halving for stability
      if (loglik < loglikPrev){
        lpOld <- 0.5*lpOld + 0.5*lpPrev[-folds[[k]]]
        lpNew <- 0.5*lpNew + 0.5*lpPrev[folds[[k]]]
      }

      # assess convergence
      if (is.infinite(loglik) | is.na(loglik) | is.nan(loglik)){
        lpNew=rep(0,nrow(Y[folds[[k]],]))
        break
      } else if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lpAll
      }
    }

    RespTot <- base::append(RespTot,Y[folds[[k]],])
    LPTot <- c(LPTot,as.numeric(lpNew))

  }
  CVLoss = .loglikSurv(Y=survival::Surv(time=RespTot[,1],event = RespTot[,2]),lp=LPTot)
  # average over the folds
  return(-CVLoss)

}

.ridgeLM <- function(Y, X, U, lambda, lambdaG = 0,
                     symFusion, Tree,
                     Dg = NULL, verbose = TRUE){
  ## ---------------------------------------------------------------------
  ## Function that evaluates ridge regression estimator with a regular and
  ## generalized penalty. The fused penalty is specified by the matrix Dg.
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix multiplied by the eigenvector matrix of the
  ##           nonnegative definite matrix that specifies the structure
  ##           of spatial fused ridge penalty.
  ## lambda  : numeric, the regular ridge penalty parameter
  ## lambdaG : numeric, the fused ridge penalty parameter
  ## Dg      : nonnegative definite matrix that specifies the structure
  ##           of the generalized ridge penalty.

  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the generalized ridge regression estimate.
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen, modified by JM Goedhart
  ## ---------------------------------------------------------------------


  if (lambdaG == 0){

    # no fusion
    message("Fit fusedTree without fusion penalty, as alpha = 0 is specified")


    if (nrow(X) >= ncol(X)){
      X <- as.matrix(X)
      # efficient evaluation for n >= p

      # evaluate subexpressions of the estimator
      tUTX        <-crossprod(X, U)
      XTXpI       <-crossprod(X)
      XTY         <-crossprod(X, Y)
      UTY         <-crossprod(U, Y)
      UTU         <-crossprod(U)
      diag(XTXpI) <- diag(XTXpI) + lambda

      # evaluate unpenalized low-dim regression estimator
      gHat <- solve(UTU -crossprod(solve(XTXpI, tUTX), tUTX),
                    (UTY - as.numeric(crossprod(solve(XTXpI, XTY), tUTX))))

      # evaluate penalized high-dim regression estimator
      bHat <- solve(XTXpI, XTY -crossprod(t(tUTX), gHat))
    }

    if (nrow(X) < ncol(X)){
      # efficient evaluation for n < p

      # evaluate subexpressions of the estimator
      XXT       <- Matrix::tcrossprod(X) / lambda
      Matrix::diag(XXT) <- Matrix::diag(XXT) + 1
      Y         <- Y

      # evaluate unpenalized low-dim regression estimator
      gHat <- Matrix::solve(Matrix::crossprod(U, Matrix::solve(XXT, U)),
                            Matrix::crossprod(U, Matrix::solve(XXT, Y)))

      # evaluate penalized high-dim regression estimator
      Y    <- Y - Matrix::crossprod(Matrix::t(U), gHat)
      bHat <- Matrix::crossprod(X, Matrix::solve(XXT, Y)) / lambda
    }
  }

  if (lambdaG > 0){

    # with fusion
    message("Fit fusedTree with fusion penalty")

    if (nrow(X) >= ncol(X)){
      # efficient evaluation for n >= p

      # evaluate subexpressions of the estimator
      Dg <- as.matrix(Dg)
      X <- as.matrix(X)
      Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
      XTX  <-crossprod(X)
      tUTX <-crossprod(X, U)
      XTY  <-crossprod(X, Y)
      UTY  <-crossprod(U, Y)
      UTU  <-crossprod(U)


      # evaluate unpenalized low-dim regression estimator
      gHat <- solve(UTU - crossprod(solve(Dg + XTX, tUTX), tUTX),
                    (UTY -crossprod(solve(Dg + XTX, XTY), tUTX)[1,]))

      # evaluate penalized high-dim regression estimator
      bHat <- solve(XTX + Dg, XTY -crossprod(t(tUTX), gHat))
    }

    if (nrow(X) < ncol(X)){
      # efficient evaluation for n < p

      # evaluate subexpressions of the estimator
      NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
      Dg <-Matrix::tcrossprod(Matrix::kronecker(Matrix::diag(1,pX),Matrix::solve(lambdaG*.PenMatr(NumNodes,1, symmetric = symFusion, Tree = Tree) +
                                                Matrix::diag(lambda,NumNodes)))
                      ,X)

      XDXTpI <- Matrix::crossprod(Matrix::t(X), Dg) + Matrix::diag(1,nrow(X))

      # evaluate unpenalized low-dim regression estimator
      gHat <- Matrix::solve(Matrix::crossprod(U, Matrix::solve(XDXTpI, U)),
                            Matrix::crossprod(U, Matrix::solve(XDXTpI, Y)))

      # evaluate penalized high-dim regression estimator
      Y    <- Y - Matrix::crossprod(t(U), gHat)
      bHat <- Matrix::tcrossprod(Dg, Matrix::t(Matrix::solve(XDXTpI, Y)))[,1]
    }
  }
  return(c(as.numeric(gHat), as.numeric(bHat)))
}


.ridgeBLM <- function(Y, X, U, lambda, lambdaG, Dg,
                      symFusion, Tree,
                      minSuccDiff, maxIter, verbose = FALSE){

  if (lambdaG == 0){

    # no fusion
    message("Fit fusedTree without fusion penalty, as alpha = 0 is specified")

    # initiate

    if (ncol(X) >= nrow(X)){
      XXT <- Matrix::tcrossprod(X) / lambda
    }
    if (ncol(X) >= nrow(X)){
      tUTX <- Matrix::crossprod(X, U)
    }
    lpPrev <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    lp      <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    loglikPrev <- .loglikBLMlp(Y, lp)

    for (iter in 1:maxIter){
      # calculate the weights
      Ypred <- 1 / (1 + exp(-lp))
      W0    <- as.numeric(Ypred * (1 - Ypred))
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <-
          .Machine$double.eps
      }

      # distinguish between the high- and low-dimensional cases
      if (ncol(X) >= nrow(X)){
        # obtain part one of the estimator
        Z <- lp + (Y - Ypred)/W0

        # now obtain the IRWLS update efficiently
        Matrix::diag(XXT) <- Matrix::diag(XXT) + 1/W0
        slh       <- Matrix::solve(XXT, Z)
        gHat      <- Matrix::solve(Matrix::crossprod(U, Matrix::solve(XXT, U)),
                                   Matrix::crossprod(U, slh))
        slh       <- Matrix::solve(XXT, Z - U %*% gHat)
        Matrix::diag(XXT) <- Matrix::diag(XXT) - 1/W0
        lp        <- Matrix::crossprod(XXT, slh)
        penalty   <- sum(lp * slh) / 2
        lp        <- as.numeric(lp) + as.numeric(U %*% gHat)
        loglik    <- .loglikBLMlp(Y, lp) - penalty
      }

      if (ncol(X) < nrow(X)){
        # adjusted response
        X <- as.matrix(X)
        Z <- W0*lp + (Y - Ypred)

        # evaluate subexpressions of the estimator
        XTXpI       <-crossprod(sweep(X, 1, sqrt(W0), FUN="*"))
        tUTX        <-crossprod(X, sweep(U, 1, W0, FUN="*"))
        XTZ         <-crossprod(X, Z)
        UTZ         <-crossprod(U, Z)
        UTU         <-crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
        diag(XTXpI) <- diag(XTXpI) + lambda

        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU -crossprod(solve(XTXpI, tUTX), tUTX),
                (UTZ - as.numeric(crossprod(solve(XTXpI, XTZ),
                                            crossprod(X, sweep(U, 1, W0, FUN="*")))))))

        # obtain part one of the estimator
        XWZpT <- as.numeric(
          crossprod(X, W0*lp + (Y - Ypred) -
                      W0*as.numeric(tcrossprod(gHat, U))))

        # now obtain the IRWLS update efficiently
        bHat <- solve(XTXpI, XWZpT)
        lp   <- as.numeric(tcrossprod(X, t(bHat)) +
                             as.numeric(tcrossprod(gHat, U)))

        # evaluate the loglikelihood
        loglik <- .loglikBLMlp(Y, lp) -
          0.5 * lambda * sum((bHat)^2)
      }
      if (verbose) {
        cat(sprintf("Iteration %2d   log likelihood equals: %10.3f\n", iter, loglik))
      }
      #step-halving for stability
      if (loglik < loglikPrev){
        lp <- 0.5*lp + 0.5*lpPrev
      }

      # assess convergence
      if(is.nan(loglik) | is.infinite(loglik) | is.na(loglik)){
        stop("convergence error, please increase penalty")}

      if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        if (verbose) {
          cat("IRLS converged at iteration ",iter, "\n")
        }
        break
      } else {
        loglikPrev <- loglik
        lpPrev <-lp
      }
      if(iter == maxIter) {
        stop("Not converged yet, please increase maxIter")}
    }
    if (ncol(X) >= nrow(X)){
      bHat <- as.numeric(Matrix::crossprod(X, slh)) / lambda
    }
  }

  if (lambdaG > 0){

    # with fusion
    message("Fit fusedTree with fusion penalty")

    if (ncol(X) >= nrow(X)){
      # evaluate subexpressions of the estimator
      NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
      Dg <-Matrix::tcrossprod(Matrix::kronecker(Matrix::diag(1,pX),Matrix::solve(lambdaG*.PenMatr(NumNodes,1, symmetric = symFusion, Tree = Tree) +
                                                Matrix::diag(lambda,NumNodes)))
                      ,X)
      XDXT        <- Matrix::crossprod(Matrix::t(X), Dg)
      diagXDXTorg <- Matrix::diag(XDXT)
    } else {
      # evaluate subexpressions of the estimator
      Dg <- as.matrix(Dg)
      Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
    }

    # initiate
    lpPrev <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    lp      <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    loglikPrev <- .loglikBLMlp(Y, lp)


    for (iter in 1:maxIter){
      # calculate the weights
      Ypred <- 1 / (1 + exp(-lp))
      W0    <- as.numeric(Ypred * (1 - Ypred))
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <-
          .Machine$double.eps
      }

      # distinguish between the high- and low-dimensional cases
      if (ncol(X) >= nrow(X)){
        # obtain part one of the estimator
        Z <- lp + (Y - Ypred)/W0

        # now obtain the IRWLS update efficiently
        Matrix::diag(XDXT) <- Matrix::diag(XDXT) + 1/W0
        slh        <- Matrix::solve(XDXT, Z)
        gHat       <- Matrix::solve(Matrix::crossprod(U, Matrix::solve(XDXT, U)),
                                    Matrix::crossprod(U, slh))
        slh        <- Matrix::solve(XDXT, Z - U %*% gHat)
        Matrix::diag(XDXT) <- diagXDXTorg
        lp         <- Matrix::crossprod(XDXT, slh)
        penalty    <- sum(lp * slh) / 2
        lp         <- as.numeric(lp) + as.numeric(U %*% gHat)
        loglik    <- .loglikBLMlp(Y, lp) - penalty
      }
      if (ncol(X) < nrow(X)){
        # adjusted response
        X <- as.matrix(X)
        Z <- W0*lp + (Y - Ypred)

        # evaluate subexpressions of the estimator
        XTXpD <-crossprod(sweep(X, 1, sqrt(W0), FUN="*"))
        tUTX  <-crossprod(X, sweep(U, 1, W0, FUN="*"))
        XTZ   <-crossprod(X, Z)
        UTZ   <-crossprod(U, Z)
        UTU   <-crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
        XTXpD <- XTXpD + Dg

        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU -crossprod(solve(XTXpD, tUTX), tUTX),
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ),
                                            crossprod(X, sweep(U, 1, W0, FUN="*")))))))

        # obtain part one of the estimator
        XWZpT <-
          crossprod(X, W0*lp + (Y - Ypred) -
                      W0*tcrossprod(gHat, U)[1,])[,1]

        # now obtain the IRWLS update efficiently
        bHat    <- qr.solve(XTXpD, XWZpT,
                            tol=.Machine$double.eps)
        penalty <- sum(crossprod(Dg, bHat) *
                         (bHat))/2
        lp      <- as.numeric(tcrossprod(X, t(bHat)) +
                                as.numeric(tcrossprod(gHat, U)))

        # evaluate the loglikelihood
        loglik    <- .loglikBLMlp(Y, lp) - penalty
      }

      if (verbose) {
        cat(sprintf("Iteration %2d   log likelihood equals: %10.3f\n", iter, loglik))
      }

      # step-halving for stability
      if (loglik < loglikPrev){
        lp <- 0.5*lp + 0.5*lpPrev #step-halving for first iteration to guide search
      }

      # assess convergence
      if(is.nan(loglik) | is.infinite(loglik) | is.na(loglik)) {
        stop("convergence error, please increase penalty")}

      if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        if (verbose) {
          cat("IRLS converged at iteration ",iter, "\n")
        }
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lp
      }
      if(iter == maxIter) {
        stop("Not converged yet, please increase maxIter")}
    }
    if (ncol(X) >= nrow(X)){
      bHat <- as.numeric(Matrix::tcrossprod(as.numeric(slh), Dg))
    }
  }
  return(c(as.numeric(gHat), as.numeric(bHat)))
}

.ridgeSurv <- function(Y, X, U, lambda, lambdaG, Dg, symFusion, Tree,
                       minSuccDiff, maxIter, verbose = FALSE){


  if (lambdaG == 0){

    # no fusion
    message("Fit fusedTree without fusion penalty, as alpha = 0 is specified")

    # initiate

    if (ncol(X) >= nrow(X)){
      XXT <- Matrix::tcrossprod(X) /lambda
    }
    if (ncol(X) >= nrow(X)){
      tUTX <- Matrix::crossprod(X, U)
    }

    #initiate
    lpPrev <- rep(0, length(Y))
    lp <- rep(0, length(Y))
    loglikPrev <- .loglikSurv(Y, lp) #penalty is zero for initial betas
    Yev = Y[,2] # event part of response

    for (iter in 1:maxIter){
      # calculate the weights

      H0 <- .breslow(Y,lp)[,2]
      Ypred <- as.numeric(H0 * exp(lp))
      W0 <- Ypred

      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <-
          2*.Machine$double.eps}

      # distinguish between the high- and low-dimensional cases
      if (ncol(X) >= nrow(X)){
        # obtain part one of the estimator
        Z <- lp + (Yev - Ypred)/W0

        # now obtain the IRWLS update efficiently
        Matrix::diag(XXT) <- Matrix::diag(XXT) + 1/W0
        slh       <- Matrix::solve(XXT, Z)
        gHat      <- Matrix::solve(Matrix::crossprod(U, Matrix::solve(XXT, U)),
                                   Matrix::crossprod(U, slh))
        slh       <- Matrix::solve(XXT, Z - U %*% gHat)
        Matrix::diag(XXT) <- Matrix::diag(XXT) - 1/W0
        lp        <- Matrix::crossprod(XXT, slh)
        penalty   <- sum(lp * slh) / 2
        lp        <- as.numeric(lp) + as.numeric(U %*% gHat)
      }

      if (ncol(X) < nrow(X)){
        # adjusted response
        X <- as.matrix(X)
        Z <- W0*lp + (Yev - Ypred)

        # evaluate subexpressions of the estimator
        XTXpI       <-crossprod(sweep(X, 1, sqrt(W0), FUN="*"))
        tUTX        <-crossprod(X, sweep(U, 1, W0, FUN="*"))
        XTZ         <-crossprod(X, Z)
        UTZ         <-crossprod(U, Z)
        UTU         <-crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
        diag(XTXpI) <- diag(XTXpI) + lambda

        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU -crossprod(solve(XTXpI, tUTX), tUTX),
                (UTZ - as.numeric(crossprod(solve(XTXpI, XTZ),
                                            crossprod(X, sweep(U, 1, W0, FUN="*")))))))

        # obtain part one of the estimator
        XWZpT <- as.numeric(
          crossprod(X, W0*lp + (Yev - Ypred) -
                      W0*as.numeric(tcrossprod(gHat, U))))

        # now obtain the IRWLS update efficiently
        bHat <- qr.solve(XTXpI, XWZpT, tol = .Machine$double.eps)
        lp   <- as.numeric(tcrossprod(X, t(bHat)) +
                             as.numeric(tcrossprod(gHat, U)))
        penalty <- 0.5 * lambda * sum((bHat)^2)
      }

      loglik    <- .loglikSurv(Y, lp) - penalty

      if (verbose) {
        cat(sprintf("Iteration %2d   log likelihood equals: %10.3f\n", iter, loglik))
      }


      # step-halving for stability
      if (loglik < loglikPrev){
        lp <- 0.5*lp + 0.5*lpPrev
      }
      # assess convergence
      if(is.nan(loglik) | is.infinite(loglik) | is.na(loglik)) {
        stop("convergence error, please increase penalty")}

      if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        if (verbose) {
        cat("IRLS converged at iteration ",iter, "\n")
        }
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lp
      }
      if(iter == maxIter) {
        stop("Not converged yet, please increase maxIter")}
    }

    if (ncol(X) >= nrow(X)){
      bHat <- as.numeric(Matrix::crossprod(X, slh)) / lambda
    }
  }

  if (lambdaG > 0){
    # with fusion
    message("Fit fusedTree with fusion penalty")

    if (ncol(X) >= nrow(X)){
      # evaluate subexpressions of the estimator
      NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
      Dg <- Matrix::tcrossprod(Matrix::kronecker(Matrix::diag(1,pX),Matrix::solve(lambdaG*.PenMatr(NumNodes,1, symmetric = symFusion, Tree = Tree) +
                               Matrix::diag(lambda,NumNodes)))
                      ,X)
      XDXT        <- Matrix::crossprod(Matrix::t(X), Dg)
      diagXDXTorg <- Matrix::diag(XDXT)
    } else {
      # evaluate subexpressions of the estimator
      Dg <- as.matrix(Dg)
      Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
    }

    #initiate
    lpPrev <- rep(0, length(Y))
    lp <- rep(0, length(Y))
    loglikPrev <- .loglikSurv(Y, lp) #penalty is zero for initial betas
    Yev <- Y[,2]

    for (iter in 1:maxIter){
      # calculate the weights

      H0 <- .breslow(Y,lp)[,2]
      Ypred <- as.numeric(H0 * exp(lp))
      W0 <- Ypred

      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <-
          2*.Machine$double.eps}

      # distinguish between the high- and low-dimensional cases
      if (ncol(X) >= nrow(X)){
        # obtain part one of the estimator
        Z <- lp + (Yev - Ypred)/W0

        # now obtain the IRWLS update efficiently
        Matrix::diag(XDXT) <- Matrix::diag(XDXT) + 1/W0
        slh        <- Matrix::solve(XDXT, Z)
        gHat       <- Matrix::solve(Matrix::crossprod(U, Matrix::solve(XDXT, U)),
                                    Matrix::crossprod(U, slh))
        slh        <- Matrix::solve(XDXT, Z - U %*% gHat)
        Matrix::diag(XDXT) <- diagXDXTorg
        lp         <- Matrix::crossprod(XDXT, slh)
        penalty    <- sum(lp * slh) / 2
        lp         <- as.numeric(lp) + as.numeric(U %*% gHat)
      }

      if (ncol(X) < nrow(X)){
        # adjusted response
        Z <- W0*lp + (Yev - Ypred)
        X <- as.matrix(X)

        # evaluate subexpressions of the estimator
        XTXpD <-crossprod(sweep(X, 1, sqrt(W0), FUN="*"))
        tUTX  <-crossprod(X, sweep(U, 1, W0, FUN="*"))
        XTZ   <-crossprod(X, Z)
        UTZ   <-crossprod(U, Z)
        UTU   <-crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
        XTXpD <- XTXpD + Dg

        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU -crossprod(solve(XTXpD, tUTX), tUTX),
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ),
                                            crossprod(X, sweep(U, 1, W0, FUN="*")))))))

        # obtain part one of the estimator
        XWZpT <-
          crossprod(X, W0*lp + (Yev - Ypred) -
                      W0*tcrossprod(gHat, U)[1,])[,1]

        # now obtain the IRWLS update efficiently
        bHat    <- qr.solve(XTXpD, XWZpT,
                            tol=.Machine$double.eps)
        penalty <- sum(crossprod(Dg, bHat) *
                         (bHat))/2
        lp      <- as.numeric(tcrossprod(X, t(bHat)) +
                                as.numeric(tcrossprod(gHat, U)))
      }
      # compute new likelihood
      loglik    <- .loglikSurv(Y, lp) - penalty

      if (verbose) {
        cat(sprintf("Iteration %2d   log likelihood equals: %10.3f\n", iter, loglik))
      }

      # step-halving for stability
      if (loglik < loglikPrev){
        lp <- 0.5*lp + 0.5*lpPrev
      }

      # assess convergence
      if(is.nan(loglik) | is.infinite(loglik) | is.na(loglik)) {
        stop("convergence error, please increase penalty")}

      if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        if (verbose) {
        cat("IRLS converged at iteration ",iter, "\n")
        }
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lp
      }
      if(iter == maxIter){
        stop("Not converged yet, please increase maxIter")}
    }

    if (ncol(X) >= nrow(X)){
      bHat <- as.numeric(Matrix::tcrossprod(as.numeric(slh), Dg))
    }
  }
  effects <- c(as.numeric(gHat), as.numeric(bHat))
  baseline <- .breslow(Y = Y, lp = lp)
  return(list(estimates = effects, Breslow = baseline))
}



###### Auxiliary functions ######
#################################

.loglikBLM <- function(Y, X, betas){

  ## ---------------------------------------------------------------------
  ## Function calculates the loglikelihood of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## betas   : regression parameter
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the loglikelihood of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  # evaluate the linear predictor
  lp <- as.numeric(tcrossprod(betas, X))

  # evaluate the loglikelihood
  loglik1 <- exp(Y * lp) / (1 + exp(lp))
  loglik2 <- exp(((Y-1) * lp)) / (1 + exp(-lp))
  loglik1[!is.finite(log(loglik1))] <- NA
  loglik2[!is.finite(log(loglik2))] <- NA
  loglik  <- sum(log(apply(cbind(loglik1, loglik2), 1, mean, na.rm=TRUE)))
  return(loglik)
}

.loglikBLMlp <- function(Y, lp){

  ## ---------------------------------------------------------------------
  ## Function calculates the loglikelihood of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## lp      : linear predictor
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the loglikelihood of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  # evaluate the loglikelihood
  loglik1 <- exp(Y * lp) / (1 + exp(lp))
  loglik2 <- exp(((Y-1) * lp)) / (1 + exp(-lp))
  loglik1[!is.finite(log(loglik1))] <- NA
  loglik2[!is.finite(log(loglik2))] <- NA
  loglik  <- sum(log(apply(cbind(loglik1, loglik2), 1, mean, na.rm=TRUE)))
  return(loglik)
}

.loglikLM <- function(Y, X, betas){

  ## ---------------------------------------------------------------------
  ## Function calculates the loglikelihood of the linear regression model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## betas   : regression parameter
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the loglikelihood of the linear regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  return(-length(Y) * (log(2 * pi * sum((Y - as.numeric(tcrossprod(as.numeric(betas), X)))^2) / length(Y)) + 1) / 2)
}

.breslow <- function(Y,lp){
  #checked to coincide with basehaz from penalized (which coincides with survfit)
  #Y survival response; lp: linear predictor (length n)
  #Returns hazard and cumulative hazard at all timepoints
  # Author: Mark A van de Wiel

  if (!inherits(Y, "Surv")) {
    stop("Survival Model, please specify Y as survival object")}
  ord <- base::order(Y)
  invord <- base::order(ord) #to put back in original order
  Ys <- Y[ord]
  di <- Ys[,2] #event
  ti <- Ys[,1] #time
  lp <- lp[ord]
  htsort <- di/rev(base::cumsum(exp(rev(lp))))
  if (min(htsort) <= .Machine$double.eps){
    htsort[which(htsort < .Machine$double.eps)] <-
      2*.Machine$double.eps}
  #ht <- htsort
  #Ht <- cumsum(htsort)
  ht <- htsort[invord]
  Ht <- base::cumsum(htsort)[invord]
  return(data.frame(ht=ht,Ht=Ht,time=ti[invord]))
}

.loglikSurv <- function(Y,lp){


  if (!inherits(Y, "Surv")) {
    stop("Survival Model, please specify Y as survival object")}

  hazards <- .breslow(Y,lp)
  di <- Y[,2]
  ht <- hazards[,1]
  Ht <- hazards[,2]
  thescore <- sum(-Ht*exp(lp))
  di1 <- which(di==1)
  if(length(di1)>0) thescore <- thescore + sum(di[di1]*(log(ht[di1])+lp[di1]))
  return(thescore)
}

.CIndexSurv <- function(Y,lp){
  # Author: Mark A van de Wiel
  if (!inherits(Y, "Surv")) {
    stop("Survival Model, please specify Y as survival object")}
  lpmin <- -lp
  thescore <- survival::concordance(response ~ lpmin)$concordance
  if(is.na(thescore)) thescore <- 0.5
  return(thescore)
}

LogLik <- function(Y,X,betas,model){
  if (model=="linear"){
    return(.loglikLM(Y,X,betas))
  }
  if (model=="logistic"){
    return(.loglikBLM(Y,X,betas))
  }
}

.sosLM <- function(Y, X, betas){

  ## ---------------------------------------------------------------------
  ## Function calculates the sum-of-squares of the linear regression model
  ## ---------------------------------------------------------------------
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## betas   : regression parameter
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the sum-of-squares of the linear regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------

  return(sum((Y - as.numeric(tcrossprod(as.numeric(betas), X)))^2))
}

