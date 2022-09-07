#' @title Polya-Gamma IRT with EM algorithm
#' @description \link{pgIRT} Estimates IRT model (binary, multinomial and dynamic) by Polya-Gamma data augmentation and EM algorithm.
#' 
#' @param data a matrix, roll-call matrix (I x J). For binary model, the matrix is allowed to contain only 1, 0 and NA. For multinomial, 1, 2, 3 and NA is allowed.
#' @param model a string, one of "bin" (binary), "bin_dyn" (dynamic binary), "multi" (multinomial) or "multi_dyn" (dynamic multinomial)
#' @param prior a list, containing prior distribution:
#' \itemize{
#'   \item \code{a0} a double (\code{bin}, \code{bin_dyn})  or 2-length vector (\code{multi}, \code{multi_dyn}), prior mean of alpha. Default is 0 (\code{bin}, \code{bin_dyn}), c(0, 0) (\code{multi}, \code{multi_dyn}).
#'   \item \code{A0} a double (\code{bin}, \code{bin_dyn}) or 2-length vector (\code{multi}, \code{multi_dyn}), prior variance of alpha. Default is 25 (\code{bin}, \code{bin_dyn}), c(25, 25) (\code{multi}, \code{multi_dyn}).
#'   \item \code{b0} a double (\code{bin}, \code{bin_dyn}) or 2-length vector (\code{multi}, \code{multi_dyn}), prior mean of beta Default is 0 (\code{bin}, \code{bin_dyn}), c(0, 0) (\code{multi}, \code{multi_dyn}).
#'   \item \code{B0} a double (\code{bin}, \code{bin_dyn}) or 2-length vector (\code{multi}, \code{multi_dyn}), prior variance of beta. Default is 25 (\code{bin}, \code{bin_dyn}), c(25, 25) (\code{multi}, \code{multi_dyn}).
#'   \item \code{theta0} a numeric vector (I length), prior mean of theta_i0 for dynamic model (\code{bin_dyn}, \code{multi_dyn}). Default is \code{rep(0, I)}. 
#'   \item \code{Delta0} a numeric vector (I length), prior variance of theta_i0 for dynamic model (\code{bin_dyn}, \code{multi_dyn}). Default is \code{rep(1, I)}.
#'   \item \code{Delta} a double, prior evolution variance of theta_it for dynamic model (\code{bin_dyn}, \code{multi_dyn}). Default is 0.01.
#' }
#' @param init a list, containing initial values (strongly recommended to use \link{make_init}):
#' \itemize{
#'   \item \code{alpha} J length vector of alpha for \code{bin} and \code{bin_dyn} model. J x 2 matrix for \code{multi} and \code{multi_dyn}.
#'   \item \code{beta} J length vector of beta for \code{bin} and \code{bin_dyn} model. J x 2 matrix for \code{multi} and \code{multi_dyn}.
#'   \item \code{theta} I length vector of theta for \code{bin} and \code{multi} model. I x T (sessions) matrix for \code{bin_dyn} and \code{multi_dyn}.
#' }
#' @param constraint an integer or integer vector (for dynamic model, the same length as the number of sessions), indicating the voter whose ideal point is always set positive.
#' @param dyn_options a list, containing the options for dynamic model. If you choose \code{bin_dyn} or \code{multi_dyn} for `model`, you must supply this argument. Using \link{make_dyn_options}() is strongly recommended:
#' \itemize{
#'   \item \code{session_individual} a matrix, the first column contains attending sessions (starts from 0) of individuals and 
#'   the second column contains the identifiers of individuals.
#'   \item \code{bill_session} an integer vector, indicating the sessions each bill belongs to (starts from 0).
#'   \item \code{bill_match} (Optional). An integer vector which indicates matched bills (the location of matched bills, starts from 0) 
#'   and is the same length as the number of bills. If a bill does not have matches, please fill NA. Using \link{make_bill_match} is strongly recommended.
#' }
#' @param tol a double (< 1), convergence threshold. Default is 1e-6.
#' @param maxit am integer, maximum number of iterations for EM. Default is 500.
#' @param verbose a integer, the function prints the status every `verbose`. Default is NULL.
#' @param std a bool, whether ideal points are standardized or not. Default is FALSE.
#' 
#' @importFrom crayon yellow
#' @importFrom stats cor na.omit
#' @useDynLib pgIRT, .registration = TRUE
#' @export

pgIRT <- function(data, 
                  model = c("bin", "bin_dyn", "multi", "multi_dyn"), 
                  prior = NULL, 
                  init = NULL, 
                  constraint = NULL, 
                  dyn_options = NULL, 
                  tol = 1e-6, 
                  maxit = 500, 
                  verbose = NULL, 
                  std = FALSE) {
  
  if (!class(data)[1] %in% c("matrix", "tbl_df", "data.frame")) {
    stop("`datamatrix` must be 'matrix', 'data.frame' or 'tbl_df' object.")
  } else {
    data <- as.matrix(data)
  }
  
  I <- nrow(data)
  J <- ncol(data)
  
  if (is.null(verbose)) verbose <- maxit + 1
  
  if (is.null(prior)) {
    prior <- list(a0 = 0,
                  A0 = 25,
                  b0 = 0,
                  B0 = 25,
                  Delta0 = rep(1, I),
                  Delta = .01,
                  theta0 = rep(0, I))
    if (model %in% c("multi", "multi_dyn")) {
      prior$a0 <- rep(0, 2)
      prior$A0 <- rep(25, 2)
      prior$b0 <- rep(0, 2)
      prior$B0 <- rep(25, 2)
    }
  }
  if (is.null(init)) {
    stop("Please supply `init`! See more detail: ?make_init.")
  }
  
  if (is.null(constraint)) {
    is_const <- FALSE
  } else {
    is_const <- TRUE
  }
  
  if (length(model) > 1) stop("Please specify `model` argument. Only one of 'bin', 'bin_dyn', 'multi' or 'multi_dyn' is allowed.")
  
  if (model %in% c("multi", "multi_dyn")) {
    md <- "Multinomial"
    ALLPOINT <- sum(is.na(data)) + sum(data == 1, na.rm = TRUE) + sum(data == 2, na.rm = TRUE) + sum(data == 3, na.rm = TRUE)
    if (ALLPOINT != I * J) stop("For multinomial model, the data is allowed to contain only NA, 1, 2 and 3.")
    Y1 <- Y2 <- data
    Y1[Y1 %in% c(2, 3)] <- 0
    Y2[Y2 == 3] <- 0; Y2[Y2 == 1] <- NA; Y2[Y2 == 2] <- 1
    category <- apply(data, 2, function(x) unique(x[!is.na(x)]))
    num_cat <- lapply(category, length) %>% unlist()
    max_cat <- lapply(category, max) %>% unlist()
  } else if (model %in% c("bin", "bin_dyn")) {
    md <- "Binomial"
    ALLPOINT <- sum(is.na(data)) + sum(data == 1, na.rm = TRUE) + sum(data == 0, na.rm = TRUE)
    if (ALLPOINT != I * J) stop("For binary model, the data is allowed to contain only NA, 1 and 0.")
    Y <- data
  }
  if (model %in% c("bin_dyn", "multi_dyn")) {
    md <- paste("Dynamic", md)
    if (is.null(dyn_options)) stop("Please supply `dyn_options`! See more detail: ?make_dyn_options.")
    session_individual <- dyn_options$session_individual
    bill_session <- dyn_options$bill_session
    if (!is.null(dyn_options$bill_match)) {
      bill_match <- dyn_options$bill_match
    } else {
      bill_match <- rep(NA, J)
    }
    if (!is.null(constraint)) {
      if (length(constraint) == 1) {
        constraint <- rep(constraint, length(unique(bill_session)))
      } else if (length(constraint) != length(unique(bill_session))) {
        stop("`constraint` must be the same length as the number of sessions.")
      }
    }
  }
  
  cat("=========================================================================\n")
  cat("Polya-Gamma data augmentation Item Response Theory Model via EM Algorithm\n")
  cat("Model =", crayon::yellow(md, "\n"))
  cat("=========================================================================\n")
  
  stime <- proc.time()[3]
  
  if (model == "bin") {
    alpha <- init$alpha
    beta <- init$beta
    theta <- init$theta
    if (is.null(constraint)) constraint <- which.max(theta)
    a0 <- prior$a0
    A0 <- prior$A0
    b0 <- prior$b0
    B0 <- prior$B0
    
    for (g in 1:maxit) {
      theta_old <- theta
      alpha_old <- alpha
      beta_old <- beta
      
      # Estep
      omega <- get_Eomega_bin(theta_old, alpha_old, beta_old)
      
      # Mstep
      theta <- update_theta_bin(Y, omega, alpha_old, beta_old, constraint, is_const)
      if (std) {
        theta <- scale(theta)[, 1]
      }
      alpha <- update_alpha_bin(Y, omega, beta_old, theta, a0, A0)
      beta <- update_beta_bin(Y, omega, alpha, theta, b0, B0)
      
      # Convergence check
      theta_cor <- stats::cor(stats::na.omit(as.numeric(theta_old)), stats::na.omit(as.numeric(theta)))
      alpha_cor <- stats::cor(alpha_old, alpha)
      beta_cor <- stats::cor(beta_old, beta)
      
      all_cor <- abs(c(theta_cor, alpha_cor, beta_cor))
      names(all_cor) <- c("theta", "alpha", "beta")
      
      if (g == 1) cor_store <- list()
      cor_store[[g]] <- all_cor
      
      if (g > 1 & max(1 - all_cor) < tol) {
        if (!is.null(rownames(Y))) names(theta) <- rownames(Y)
        if (!is.null(colnames(Y))) names(alpha) <- names(beta) <- colnames(Y)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "bin",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = TRUE)
        class(L) <- "pgIRT_fit"
        cat("Model converged at iteration", g, ":", crayon::yellow(round(proc.time()[3] - stime, 1), "sec\n"))
        return(L)
      }
      
      if (g == maxit) {
        if (!is.null(rownames(Y))) names(theta) <- rownames(Y)
        if (!is.null(colnames(Y))) names(alpha) <- names(beta) <- colnames(Y)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "bin",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = FALSE)
        class(L) <- "pgIRT_fit"
        warning("Model falied to converge!\n Return the result as it is, but it may be unreliable : total time ",
                round(proc.time()[3] - stime, 1), " sec\n")
        return(L)
      }
      
      if (g %% verbose == 0) {
        cat("Iteration", g, ": eval =", names(which.max(1 - all_cor)),
            max(1 - all_cor), crayon::yellow("elapsed", round(proc.time()[3] - stime, 1), "sec\n"))
      }
    }
  }
  
  if (model == "bin_dyn") {
    alpha <- init$alpha
    beta <- init$beta
    theta <- init$theta
    if (is.null(constraint)) constraint <- rep(which.max(theta[, 1]), I)
    a0 <- prior$a0
    A0 <- prior$A0
    b0 <- prior$b0
    B0 <- prior$B0
    Delta0 <- prior$Delta0
    Delta <- prior$Delta
    theta0 <- prior$theta0
    
    for (g in 1:maxit) {
      theta_old <- theta
      alpha_old <- alpha
      beta_old <- beta
      
      # Estep
      omega <- get_Eomega_bin_dyn(theta_old, alpha_old, beta_old, bill_session)
      
      # Mstep
      theta <- update_theta_bin_dyn(Y, omega, alpha_old, beta_old, theta0, Delta0, Delta,
                                    constraint, is_const, session_individual, bill_session)
      if (std) {
        theta <- apply(theta, 2, scale)
      }
      alpha <- update_alpha_bin_dyn(Y, omega, beta_old, theta, a0, A0, bill_session, bill_match)
      beta <- update_beta_bin_dyn(Y, omega, alpha, theta, b0, B0, bill_session)
      
      if (g == 1) theta_old[is.na(theta)] <- NA
      
      # Convergence check
      theta_cor <- stats::cor(stats::na.omit(as.numeric(theta_old)), stats::na.omit(as.numeric(theta)))
      alpha_cor <- stats::cor(alpha_old, alpha)
      beta_cor <- stats::cor(beta_old, beta)
      
      all_cor <- abs(c(theta_cor, alpha_cor, beta_cor))
      names(all_cor) <- c("theta", "alpha", "beta")
      
      if (g == 1) cor_store <- list()
      cor_store[[g]] <- all_cor
      
      if (g > 1 & max(1 - all_cor) < tol) {
        if (!is.null(rownames(Y))) rownames(theta) <- rownames(Y)
        if (!is.null(colnames(Y))) names(alpha) <- names(beta) <- colnames(Y)
        colnames(theta) <- 0:(ncol(theta) -1)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint, dyn_options = dyn_options)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "bin_dyn",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = TRUE)
        class(L) <- "pgIRT_fit"
        cat("Model converged at iteration", g, ":", crayon::yellow(round(proc.time()[3] - stime, 1), "sec\n"))
        return(L)
      }
      
      if (g == maxit) {
        if (!is.null(rownames(Y))) rownames(theta) <- rownames(Y)
        if (!is.null(colnames(Y))) names(alpha) <- names(beta) <- colnames(Y)
        colnames(theta) <- 0:(ncol(theta) -1)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint, dyn_options = dyn_options)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "bin_dyn",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = FALSE)
        class(L) <- "pgIRT_fit"
        warning("Model falied to converge!\n Return the result as it is, but it may be unreliable : total time ",
                round(proc.time()[3] - stime, 1), " sec\n")
        return(L)
      }
      
      if (g %% verbose == 0) {
        cat("Iteration", g, ": eval =", names(which.max(1 - all_cor)),
            max(1 - all_cor), crayon::yellow("elapsed", round(proc.time()[3] - stime, 1), "sec\n"))
      }
    }
  }
  
  if (model == "multi") {
    alpha <- init$alpha
    beta <- init$beta
    theta <- init$theta
    if (is.null(constraint)) constraint <- which.max(theta)
    a0 <- prior$a0
    A0 <- prior$A0
    b0 <- prior$b0
    B0 <- prior$B0
    
    for (g in 1:maxit) {
      theta_old <- theta
      alpha_old <- alpha
      beta_old <- beta
      
      # Estep
      omega <- get_Eomega_mlt(theta_old, alpha_old, beta_old)
      
      # Mstep
      theta <- update_theta_mlt(Y1, Y2, omega, alpha_old, beta_old, constraint, is_const, max_cat, num_cat)
      if (std) {
        theta <- scale(theta)[, 1]
      }
      alpha <- update_alpha_mlt(Y1, Y2, omega, beta_old, theta, a0, A0)
      beta <- update_beta_mlt(Y1, Y2, omega, alpha, theta, b0, B0)
      
      # Convergence check
      theta_cor <- stats::cor(theta_old, theta)
      alpha1_cor <- stats::cor(alpha_old[, 1], alpha[, 1])
      alpha2_cor <- stats::cor(stats::na.omit(alpha_old[, 2]), stats::na.omit(alpha[, 2]))
      beta1_cor <- stats::cor(beta_old[, 1], beta[, 1])
      beta2_cor <- stats::cor(stats::na.omit(beta_old[, 2]), stats::na.omit(beta[, 2]))
      
      all_cor <- abs(c(theta_cor, alpha1_cor, alpha2_cor, beta1_cor, beta2_cor))
      names(all_cor) <- c("theta", "alpha1", "alpha2", "beta1", "beta2")
      
      if (g == 1) cor_store <- list()
      cor_store[[g]] <- all_cor
      
      if (g > 1 & max(1 - all_cor) < tol) {
        if (!is.null(rownames(Y1))) names(theta) <- rownames(Y1)
        if (!is.null(colnames(Y1))) rownames(alpha) <- rownames(beta) <- colnames(Y1)
        colnames(alpha) <- c('alpha1', 'alpha2')
        colnames(beta) <- c('beta1', 'beta2')
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "multi",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = TRUE)
        class(L) <- "pgIRT_fit"
        cat("Model converged at iteration", g, ":", crayon::yellow(round(proc.time()[3] - stime, 1), "sec\n"))
        return(L)
      }
      
      if (g == maxit) {
        if (!is.null(rownames(Y1))) names(theta) <- rownames(Y1)
        if (!is.null(colnames(Y1))) rownames(alpha) <- rownames(beta) <- colnames(Y1)
        colnames(alpha) <- c('alpha1', 'alpha2')
        colnames(beta) <- c('beta1', 'beta2')
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "multi",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = FALSE)
        class(L) <- "pgIRT_fit"
        warning("Model falied to converge!\n Return the result as it is, but it may be unreliable : total time ",
                round(proc.time()[3] - stime, 1), " sec\n")
        return(L)
      }
      
      if (g %% verbose == 0) {
        cat("Iteration", g, ": eval =", names(which.max(1 - all_cor)),
            max(1 - all_cor), crayon::yellow("elapsed", round(proc.time()[3] - stime, 1), "sec\n"))
      }
    }
  }
  
  if (model == "multi_dyn") {
    alpha <- init$alpha
    beta <- init$beta
    theta <- init$theta
    if (is.null(constraint)) constraint <- rep(which.max(theta[, 1]), length(unique(bill_session)))
    a0 <- prior$a0
    A0 <- prior$A0
    b0 <- prior$b0
    B0 <- prior$B0
    theta0 <- prior$theta0
    Delta0 <- prior$Delta0
    Delta <- prior$Delta
    
    for (g in 1:maxit) {
      theta_old <- theta
      alpha_old <- alpha
      beta_old <- beta
      
      # Estep
      omega <- get_Eomega_mlt_dyn(theta_old, alpha_old, beta_old, bill_session)
      
      # Mstep
      theta <- update_theta_mlt_dyn(Y1, Y2, omega, alpha_old, beta_old, theta0, Delta0, Delta, constraint,
                                    is_const,
                                    session_individual, bill_session, max_cat, num_cat)
      if (std) {
        theta <- apply(theta, 2, scale)
      }
      alpha <- update_alpha_mlt_dyn(Y1, Y2, omega, beta_old, theta, a0, A0, bill_session, bill_match)
      beta <- update_beta_mlt_dyn(Y1, Y2, omega, alpha, theta, b0, B0, bill_session)
      
      # Convergence check
      if (g == 1) theta_old[is.na(theta)] <- NA
      theta_cor <- stats::cor(stats::na.omit(as.numeric(theta_old)), stats::na.omit(as.numeric(theta)))
      alpha1_cor <- stats::cor(alpha_old[, 1], alpha[, 1])
      alpha2_cor <- stats::cor(stats::na.omit(alpha_old[, 2]), stats::na.omit(alpha[, 2]))
      beta1_cor <- stats::cor(beta_old[, 1], beta[, 1])
      beta2_cor <- stats::cor(stats::na.omit(beta_old[, 2]), stats::na.omit(beta[, 2]))
      
      all_cor <- abs(c(theta_cor, alpha1_cor, alpha2_cor, beta1_cor, beta2_cor))
      names(all_cor) <- c("theta", "alpha1", "alpha2", "beta1", "beta2")
      
      if (g == 1) cor_store <- list()
      cor_store[[g]] <- all_cor
      
      if (g > 1 & max(1 - all_cor) < tol) {
        if (!is.null(rownames(Y1))) rownames(theta) <- rownames(Y1)
        if (!is.null(colnames(Y1))) rownames(alpha) <- rownames(beta) <- colnames(Y1)
        colnames(theta) <- 0:(ncol(theta) -1)
        colnames(alpha) <- c('alpha1', 'alpha2')
        colnames(beta) <- c('beta1', 'beta2')
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint, dyn_options = dyn_options)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "multi_dyn",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = TRUE)
        class(L) <- "pgIRT_fit"
        cat("Model converged at iteration", g, ":", crayon::yellow(round(proc.time()[3] - stime, 1), "sec\n"))
        return(L)
      }
      
      if (g == maxit) {
        if (!is.null(rownames(Y1))) rownames(theta) <- rownames(Y1)
        if (!is.null(colnames(Y1))) rownames(alpha) <- rownames(beta) <- colnames(Y1)
        colnames(theta) <- 0:(ncol(theta) -1)
        colnames(alpha) <- c('alpha1', 'alpha2')
        colnames(beta) <- c('beta1', 'beta2')
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint, dyn_options = dyn_options)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "multi_dyn",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = FALSE)
        class(L) <- "pgIRT_fit"
        warning("Model falied to converge!\n Return the result as it is, but it may be unreliable : total time ",
                round(proc.time()[3] - stime, 1), " sec\n")
        return(L)
      }
      
      if (g %% verbose == 0) {
        cat("Iteration", g, ": eval =", names(which.max(1 - all_cor)),
            max(1 - all_cor), crayon::yellow("elapsed", round(proc.time()[3] - stime, 1), "sec\n"))
      }
    }
  }
}

