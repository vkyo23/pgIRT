#' @title Polya-Gamma IRT with EM algorithm
#' @description \link{pgIRT} estimates IRT model with Polya-Gamma data augmentation and EM algorithm. The function is available for 
#' the data containing K >= 2 response category (binary ~ K-multinomial) and having different categories across items j.
#' 
#' @param data a matrix, data.frame or tbl_df of roll-call (i.e. individual-item) matrix (I x J). The elements should start from 1 and missing values should be NA.
#' @param model a string, one of "default" or "dynamic".
#' @param prior a list, containing prior distribution:
#' \itemize{
#'   \item \code{a0} a double or J x K matrix, prior mean of alpha. Default is 0.
#'   \item \code{A0} a double or J x K matrix, prior variance of alpha. Default is 25.
#'   \item \code{b0} a double or J x K matrix, prior mean of beta Default is 0.
#'   \item \code{B0} a double or J x K matrix, prior variance of beta. Default is 25.
#'   \item \code{theta0} a double or I length vector, prior mean of theta_i0 for dynamic model. Default is 0. 
#'   \item \code{Delta0} a double or I length vector, prior variance of theta_i0 for dynamic model. Default is 1.
#'   \item \code{Delta} a double, prior evolution variance of theta_it for dynamic model. Default is 0.01.
#' }
#' @param init a list, containing initial values. If not supplied, the function automatically set them from the data.
#' \itemize{
#'   \item \code{alpha} J x K-1 matrix of alpha.
#'   \item \code{beta} J x K-1 matrix of beta.
#'   \item \code{theta} I x T (sessions) matrix of theta. For default model, T = 1.
#' }
#' @param constraint an integer scalar or vector (for dynamic model, the same length as the number of sessions), 
#' index of the voter i (the location of i) whose ideal point is always set positive.
#' @param dyn_options a list, containing the options for dynamic model. If you choose "dynamic" for `model`, you must supply this argument. Using \link{make_dyn_options}() is strongly recommended:
#' \itemize{
#'   \item \code{session_individual} an integer matrix, the first column contains attending sessions (starts from 1) of individuals and 
#'   the second column contains the identifiers of individuals.
#'   \item \code{bill_session} an integer vector, indicating the sessions each bill belongs to (starts from 1).
#'   \item \code{matched_bill} (Optional). An integer vector which indicates matched bills (the location of matched bills, starts from 1) 
#'   and is the same length as the number of bills. If a bill does not have matches, please fill NA. Using \link{make_bill_match} is strongly recommended.
#' }
#' @param tol a double (< 1), convergence threshold. Default is 1e-6.
#' @param maxit am integer, maximum number of iterations for EM. Default is 500.
#' @param verbose a integer, the function prints the status every `verbose`. Default is NULL.
#' @param std a bool, whether ideal points are standardized or not. Default is FALSE.
#' 
#' @importFrom stats na.omit
#' @importFrom Rcpp sourceCpp
#' 
#' @useDynLib pgIRT, .registration = TRUE
#' @export
#' 
#' @examples 
#' \dontrun{
#' # Binary model (K = 2)
#' I <- 100
#' J <- 1000
#' theta_true <- seq(-2, 2, length = I)
#' alpha_true <- rnorm(J)
#' beta_true <- rnorm(J)
#' Y_star <- cbind(1, theta_true) %*% rbind(alpha_true, beta_true)
#' Y <- matrix(rbinom(I * J, 1, plogis(Y_star)), I, J)
#' Y[Y == 0] <- 2 # Elements should start from 1 for `pgIRT`
#' pgfit <- pgIRT(Y,
#'                model = "default",
#'                constraint = which.max(theta_true))
#' cor(pgfit$parameter$theta, theta_true)
#' 
#'
#' # Multinomial model (K >= 2)
#' data(m_data_dyn)
#' mat <- make_rollcall(m_data_dyn,
#'                      unit_id = "unit",
#'                      bill_id = "bill",
#'                      vote_col = "vote")
#' mat <- clean_rollcall(mat)
#' pgfit_mlt <- pgIRT(mat,
#'                    model = "default",
#'                    constraint = 1)
#' summary(pgfit_mlt, parameter = "theta")
#' 
#'
#' # Dynamic multinomial model 
#' ops <- make_dyn_ops(m_data_dyn,
#'                     unit_id = "unit",
#'                     bill_id = "bill",
#'                     time_id = "time",
#'                     vote_col = "vote",
#'                     clean = TRUE)
#' pgfit_mlt_dyn <- pgIRT(mat,
#'                        model = "dynamic",
#'                        dyn_options = ops,
#'                        constraint = 1)
#' summary(pgfit_mlt_dyn, parameter = "theta")
#' }


pgIRT <- function(data, 
                  model = c('default', 'dynamic'),
                  init = NULL,
                  prior = NULL,
                  constraint = NULL,
                  std = FALSE,
                  dyn_options = NULL,
                  maxit = 500,
                  tol = 1e-6,
                  verbose = NULL) {
  cat("=========================================================================\n")
  cat("Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm\n")
  cat("=========================================================================\n")
  # Input check
  ## Data
  if (!class(data)[1] %in% c("matrix", "tbl_df", "data.frame")) {
    stop("`data` must be 'matrix', 'data.frame' or 'tbl_df' object.")
  } else {
    data <- as.matrix(data)
  }
  if (is.null(rownames(data))) rownames(data) <- 1:nrow(data)
  if (is.null(colnames(data))) colnames(data) <- 1:ncol(data)
  ## Get the size of data
  I <- nrow(data)
  J <- ncol(data)
  ## Get categories 
  num_cat <- apply(data, 2, function(x) length(na.omit(unique(x))))
  max_cat <- apply(data, 2, max, na.rm = TRUE)
  K <- max(num_cat)
  ## Model
  if (length(model) > 1) stop("Please specify `model` argument. Only one of 'default' or 'dynamic' is allowed.")
  if (!model %in% c('default', 'dynamic')) stop("Please specify `model` argument. Only one of 'default' or 'dynamic' is allowed.")
  control <- list(maxit = maxit, 
                  tol = tol, 
                  verbose = verbose, 
                  std = std)
  if (model == "default") {
    input <- list(data = data, 
                  max_cat = max_cat,
                  prior = NULL, 
                  init = NULL,
                  constraint = NULL)
  } else {
    input <- list(data = data, 
                  max_cat = max_cat,
                  prior = NULL, 
                  init = NULL,
                  constraint = NULL,
                  dyn_options = dyn_options)
  }
  ## Check whether the data is binary matrix
  check_bin <- length(unique(na.omit(as.vector(data)))) == 2
  if (check_bin) {
    cat('= Format ------> Binary\n')
  } else {
    ## Stop if k = 0 detected
    if (any(as.vector(na.omit(data)) == 0)) {
      stop('`pgIRT` does not support k = 0 for multinomial model. The response categories must starts from 1 and missing values must be NA.')
    }
    cat('= Format ------> Multinomial')
    cat(' (# of categories: Min =', min(num_cat), '/ Max =', K, ')\n')
  }
  cat("= Model ------->", ifelse(model == 'default', 'Default pgIRT', 'Dynamic pgIRT'), "\n")
  ## Dynamic options
  if (model == 'dynamic') {
    if (is.null(dyn_options)) stop("Please supply `dyn_options`! See more detail: ?make_dyn_options.")
    session_individual <- dyn_options$session_individual
    if (min(session_individual[, 1]) != 1) stop("The first column of `dyn_options$session_individual` should start from 1.")
    session_individual[, 1] <- session_individual[, 1] - 1
    bill_session <- dyn_options$bill_session
    if (min(bill_session) != 1) stop("`dyn_options$bill_session` should start from 1.")
    bill_session <- bill_session - 1
    if (!is.null(dyn_options$matched_bill)) {
      matched_bill <- dyn_options$matched_bill
      matched_bill <- matched_bill - 1
    } else {
      matched_bill <- rep(NA, J)
    }
  }
  ## Constraint
  if (is.null(constraint)) {
    constraint <- 1
    if (model == 'dynamic') constraint <- rep(1, I)
  } else {
    if (model == 'default') {
      constraint <- constraint
      if (length(constraint) > 1) stop('`constraint` must be a scalar, not a vector.')
    } else {
      if (length(constraint) == 1) {
        constraint <- rep(constraint, length(unique(bill_session)))
      } else if (length(constraint) != length(unique(bill_session))) {
        stop("`constraint` must be the same length as the number of sessions.")
      }
    }
  }
  ## Initial values
  if (is.null(init)) {
    init <- make_init(data,
                      model = model,
                      T = ifelse(model == 'dynamic', length(unique(bill_session)), 1),
                      constraint = constraint)
  } else {
    inls <- c('alpha', 'beta', 'theta')
    for (l in 1:length(inls)) {
      if (!exists(inls[l], init)) {
        stop(paste0('A list `init` does not contain `', inls[l], '`.'))
      } else if ((inls[l] %in% c('alpha', 'beta') & nrow((init[inls[l]])[[1]]) != J) |
                 (inls[l] %in% c('alpha', 'beta') & ncol(init[inls[l]][[1]]) != K-1)) {
        stop(paste0('`', inls[l], '` must be J x K-1 matrix'))
      } else if (inls[l] == 'theta' & model == 'dynamic') {
        if (nrow(init[inls[l]][[1]]) != I | ncol(init[inls[l]][[1]]) != length(unique(bill_session))) {
          stop('`theta` must be I x T matrix for dynamic model.') 
        }
      } else if ((inls[l] == 'theta' & model == 'default' & nrow(init[inls[l]][[1]]) != I) |
                 (inls[l] == 'theta' & model == 'default' & ncol(init[inls[l]][[1]]) != 1)) {
        stop('`theta` must be I x 1 matrix for default model.')
      }
    }
  }
  ## Priors
  if (is.null(prior)) {
    if (model == 'default') {
      prior <- list(a0 = matrix(0, J, K-1),
                    A0 = matrix(25, J, K-1),
                    b0 = matrix(0, J, K-1),
                    B0 = matrix(25, J, K-1))
    } else {
      prior <- list(a0 = matrix(0, J, K-1),
                    A0 = matrix(25, J, K-1),
                    b0 = matrix(0, J, K-1),
                    B0 = matrix(25, J, K-1),
                    theta0 = rep(0, I),
                    Delta0 = rep(1, I),
                    Delta = matrix(0.01, I, length(unique(bill_session))))
    }
  } else {
    if (model == 'default') {
      prls <- c('a0', 'A0', 'b0', 'B0')
      for (l in 1:length(prls)) {
        if (!exists(prls[l], prior)) {
          stop(paste0('A list `prior` does not contain `', prls[l], '`.'))
        } else if (length(unlist(prior[prls[l]])) == 1) {
          prior[prls[l]][[1]] <- matrix(unlist(prior[prls[l]]), nrow(init$alpha), ncol(init$alpha))
        } else if (length(unlist(prior[prls[l]])) > 1 &
                   length(unlist(prior[prls[l]])) != nrow(init$alpha) * ncol(init$alpha)) {
          stop(paste0('A prior `', prls[l], '` must be the same length as the number of items.'))
        }
      }
    } else {
      prls <- c('a0', 'A0', 'b0', 'B0', 'theta0', 'Delta0', 'Delta')
      for (l in 1:length(prls)) {
        if (!exists(prls[l], prior)) {
          stop(paste0('Please supply a prior`', prls[l], '`.'))
        } 
        if (prls[l] %in% c('a0', 'A0', 'b0', 'B0')) {
          if (length(unlist(prior[prls[l]])) == 1) {
            prior[prls[l]][[1]] <- matrix(unlist(prior[prls[l]]), nrow(init$alpha), ncol(init$alpha))
          } else if (length(unlist(prior[prls[l]])) != nrow(init$alpha) * ncol(init$alpha)) {
            stop(paste0('A prior `', prls[l], '` must be a scalar or J x K matrix.'))
          }
        } else if (prls[l] %in% c('theta0', 'Delta0')) {
          if (length(unlist(prior[prls[l]])) == 1) {
            prior[prls[l]][[1]] <- rep(unlist(prior[prls[l]]), nrow(data))
          } else if (length(unlist(prior[prls[l]])) != nrow(data)) {
            stop(paste0('A prior `', prls[l], '` must be a scalar or I length vector.'))
          }
        } else if (prls[l] == "Delta") {
          if (length(unlist(prior[prls[l]])) != 1) {
            if (length(unlist(prior[prls[l]])) != nrow(data) * length(unique(bill_session))) {
              stop(paste0('A prior `Delta` must be a scalar or I x T matrix.'))
            }
          } else {
            prior[prls[l]][[1]] <- matrix(unlist(prior[prls[l]]), nrow(data), length(unique(bill_session)))
          }
        }
      }
    }
  }
  ## Tol
  if (tol >= 1) stop('`tol` must be lower than 1.')
  input$prior <- prior
  input$init <- init
  input$constraint <- constraint
  # Implementation
  stime <- proc.time()[3]
  cat("=\n---------- Implementing EM ----------\n")
  if (is.null(verbose)) verbose <- maxit + 1
  if (model == 'default') {
    em <- EMstep(Y = data,
                 alpha = init$alpha,
                 beta = init$beta,
                 theta = init$theta,
                 max_cat = max_cat,
                 constraint = constraint - 1,
                 a0 = prior$a0,
                 A0 = prior$A0,
                 b0 = prior$b0,
                 B0 = prior$B0,
                 maxit = maxit,
                 tol = tol,
                 verbose = verbose,
                 std = std)
    el <- round(proc.time()[3] - stime, 1)
    if (em$conv) {
      cat('Model converged at iteration', em$iter, ':', el, 'sec.\n')
    } else {
      warning('Model failed to converge :', el, 'sec.\n')
    }
    colnames(em$cor) <- c('alpha', 'beta', 'theta')
    rownames(em$theta) <- rownames(data)
    colnames(em$theta) <- paste0('theta_', 1:ncol(em$theta))
    rownames(em$alpha) <- rownames(em$beta) <- colnames(data)
    colnames(em$alpha) <- paste0('alpha_', 1:ncol(em$alpha))
    colnames(em$beta) <- paste0('beta_', 1:ncol(em$beta))
  } else {
    em <- dyn_EMstep(Y = data,
                     alpha = init$alpha,
                     beta = init$beta,
                     theta = init$theta,
                     max_cat = max_cat,
                     constraint = constraint - 1,
                     a0 = prior$a0,
                     A0 = prior$A0,
                     b0 = prior$b0,
                     B0 = prior$B0,
                     theta0 = prior$theta0,
                     Delta0 = prior$Delta0,
                     Delta = prior$Delta,
                     session_individual = session_individual,
                     bill_session = bill_session,
                     matched_bill = matched_bill,
                     maxit = maxit,
                     tol = tol,
                     verbose = verbose,
                     std = std)
    el <- round(proc.time()[3] - stime, 1)
    if (em$conv) {
      cat('Model converged at iteration', em$iter, ':', el, 'sec.\n')
    } else {
      warning('Model failed to converge :', el, 'sec.\n')
    }
    colnames(em$cor) <- c('alpha', 'beta', 'theta')
    rownames(em$theta) <- rownames(data)
    colnames(em$theta) <- paste0('theta_', 1:ncol(em$theta))
    rownames(em$alpha) <- rownames(em$beta) <- colnames(data)
    colnames(em$alpha) <- paste0('alpha_', 1:ncol(em$alpha))
    colnames(em$beta) <- paste0('beta_', 1:ncol(em$beta))
  }
  L <- list(model = model,
            parameter = list(alpha = em$alpha,
                             beta = em$beta,
                             theta = em$theta),
            iteration = em$iter,
            converge = em$conv,
            convstat = em$cor,
            input = input,
            control = control)
  class(L) <- 'pgIRT_fit'
  return(L)
}

