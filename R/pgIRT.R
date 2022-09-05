#' @title Polya-Gamma IRT with EM algorithm
#' @description \link{pgIRT} Estimates IRT model (binary, multinomial and dynamic) by Polya-Gamma data augmentation and EM algorithm. 
#'
#' @param data a matrix, roll-call matrix. For binary model, the matrix is allowed to contain only 1, 0 and NA. For multinomial, 1, 2, 3 and NA is allowed.
#' @param prior a list, containing prior distribution.
#' @param init a list, containing initial values.
#' @param constraint a integer or vector, indicating the voter whose ideal point is always set positive
#' @param model string, one of "bin" (binary), "bin_dyn" (dynamic binary), "multi" (multinomial) or "multi_dyn" (dynamic multinomial)
#' @param dyn_options, a list, containing the options for dynamic model. If you choose "bin_dyn" or "multi_dyn" for `model`, you must supply this argument.
#' @param tol a double, convergence threshold. Default is 1e-6.
#' @param maxit a double, maximum number of iterations. Default is 500.
#' @param verbose a double, the function prints the status every `verbose`. 
#' @useDynLib pgIRT, .registration = TRUE
#' @export

pgIRT <- function(data, prior = NULL, init = NULL, constraint = NULL, model = c("bin", "bin_dyn", "multi", "multi_dyn"), 
                  dyn_options = NULL, tol = 1e-6, maxit = 500, verbose = NULL, std = FALSE) {
  
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
    md <- "Multinomial (Stick-Breaking)"
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
    if (ALLPOINT != I * J) stop("For binomial model, the data is allowed to contain only NA, 1 and 0.")
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
      if (length(constraint) != length(unique(bill_session))) {
        stop("`constraint` must be the same length as sessions.")
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
      theta_cor <- cor(na.omit(as.numeric(theta_old)), na.omit(as.numeric(theta)))
      alpha_cor <- cor(alpha_old, alpha)
      beta_cor <- cor(beta_old, beta)
      
      all_cor <- abs(c(theta_cor, alpha_cor, beta_cor))
      names(all_cor) <- c("theta", "alpha", "beta")
      
      if (g == 1) cor_store <- list()
      cor_store[[g]] <- all_cor
      
      if (g > 1 & max(1 - all_cor) < tol) {
        if (!is.null(rownames(Y))) names(theta) <- rownames(Y)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "bin",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = TRUE)
        class(L) <- "pgIRT"
        cat("Model converged at iteration", g, ":", crayon::yellow(round(proc.time()[3] - stime, 2), "sec\n"))
        return(L)
      }
      
      if (g == maxit) {
        if (!is.null(rownames(Y))) names(theta) <- rownames(Y)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "bin",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = FALSE)
        class(L) <- "pgIRT"
        warning("Model falied to converge!\n Return the result as it is, but it may be unreliable : total time ", 
                round(proc.time()[3] - stime, 2), " sec\n")
        return(L)
      }
      
      if (g %% verbose == 0) {
        cat("Iteration", g, ": eval =", names(which.max(1 - all_cor)), 
            max(1 - all_cor), crayon::yellow("elapsed", round(proc.time()[3] - stime, 2), "sec\n"))
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
      theta_cor <- cor(na.omit(as.numeric(theta_old)), na.omit(as.numeric(theta)))
      alpha_cor <- cor(alpha_old, alpha)
      beta_cor <- cor(beta_old, beta)
      
      all_cor <- abs(c(theta_cor, alpha_cor, beta_cor))
      names(all_cor) <- c("theta", "alpha", "beta")
      
      if (g == 1) cor_store <- list()
      cor_store[[g]] <- all_cor
      
      if (g > 1 & max(1 - all_cor) < tol) {
        if (!is.null(rownames(Y))) rownames(theta) <- rownames(Y)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint, dyn_options = dyn_options)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "bin_dyn", 
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = TRUE)
        class(L) <- "pgIRT"
        cat("Model converged at iteration", g, ":", crayon::yellow(round(proc.time()[3] - stime, 2), "sec\n"))
        return(L)
      }
      
      if (g == maxit) {
        if (!is.null(rownames(Y))) rownames(theta) <- rownames(Y)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint, dyn_options = dyn_options)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "bin_dyn",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = FALSE)
        class(L) <- "pgIRT"
        warning("Model falied to converge!\n Return the result as it is, but it may be unreliable : total time ", 
                round(proc.time()[3] - stime, 2), " sec\n")
        return(L)
      }
      
      if (g %% verbose == 0) {
        cat("Iteration", g, ": eval =", names(which.max(1 - all_cor)), 
            max(1 - all_cor), crayon::yellow("elapsed", round(proc.time()[3] - stime, 2), "sec\n"))
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
      theta_cor <- cor(theta_old, theta)
      alpha1_cor <- cor(alpha_old[, 1], alpha[, 1])
      alpha2_cor <- cor(na.omit(alpha_old[, 2]), na.omit(alpha[, 2]))
      beta1_cor <- cor(beta_old[, 1], beta[, 1])
      beta2_cor <- cor(na.omit(beta_old[, 2]), na.omit(beta[, 2]))
      
      all_cor <- abs(c(theta_cor, alpha1_cor, alpha2_cor, beta1_cor, beta2_cor))
      names(all_cor) <- c("theta", "alpha1", "alpha2", "beta1", "beta2")
      
      if (g == 1) cor_store <- list()
      cor_store[[g]] <- all_cor
      
      if (g > 1 & max(1 - all_cor) < tol) {
        if (!is.null(rownames(Y1))) names(theta) <- rownames(Y1)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "multi",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = TRUE)
        class(L) <- "pgIRT"
        cat("Model converged at iteration", g, ":", crayon::yellow(round(proc.time()[3] - stime, 2), "sec\n"))
        return(L)
      }
      
      if (g == maxit) {
        if (!is.null(rownames(Y1))) names(theta) <- rownames(Y1)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "multi",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = FALSE)
        class(L) <- "pgIRT"
        warning("Model falied to converge!\n Return the result as it is, but it may be unreliable : total time ", 
                round(proc.time()[3] - stime, 2), " sec\n")
        return(L)
      }
      
      if (g %% verbose == 0) {
        cat("Iteration", g, ": eval =", names(which.max(1 - all_cor)), 
            max(1 - all_cor), crayon::yellow("elapsed", round(proc.time()[3] - stime, 2), "sec\n"))
      }
    }
  }
  
  if (model == "multi_dyn") {
    alpha <- init$alpha
    beta <- init$beta
    theta <- init$theta
    if (is.null(constraint)) constraint <- rep(which.max(theta[, 1]), I)
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
      theta_cor <- cor(na.omit(as.numeric(theta_old)), na.omit(as.numeric(theta)))
      alpha1_cor <- cor(alpha_old[, 1], alpha[, 1])
      alpha2_cor <- cor(na.omit(alpha_old[, 2]), na.omit(alpha[, 2]))
      beta1_cor <- cor(beta_old[, 1], beta[, 1])
      beta2_cor <- cor(na.omit(beta_old[, 2]), na.omit(beta[, 2]))
      
      all_cor <- abs(c(theta_cor, alpha1_cor, alpha2_cor, beta1_cor, beta2_cor))
      names(all_cor) <- c("theta", "alpha1", "alpha2", "beta1", "beta2")
      
      if (g == 1) cor_store <- list()
      cor_store[[g]] <- all_cor
      
      if (g > 1 & max(1 - all_cor) < tol) {
        if (!is.null(rownames(Y1))) rownames(theta) <- rownames(Y1)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint, dyn_options = dyn_options)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "multi_dyn",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = TRUE)
        class(L) <- "pgIRT"
        cat("Model converged at iteration", g, ":", crayon::yellow(round(proc.time()[3] - stime, 2), "sec\n"))
        return(L)
      }
      
      if (g == maxit) {
        if (!is.null(rownames(Y1))) rownames(theta) <- rownames(Y1)
        parameter <- list(theta = theta, alpha = alpha, beta = beta)
        input <- list(data = data, constraint = constraint, dyn_options = dyn_options)
        control <- list(maxit = maxit, tol = tol, verbose = verbose, std = std)
        L <- list(parameter = parameter, prior = prior, init = init,
                  input = input, control = control, model = "multi_dyn",
                  correlation = dplyr::bind_rows(cor_store), iter = g, converge = FALSE)
        class(L) <- "pgIRT"
        warning("Model falied to converge!\n Return the result as it is, but it may be unreliable : total time ", 
                round(proc.time()[3] - stime, 2), " sec\n")
        return(L)
      }
      
      if (g %% verbose == 0) {
        cat("Iteration", g, ": eval =", names(which.max(1 - all_cor)), 
            max(1 - all_cor), crayon::yellow("elapsed", round(proc.time()[3] - stime, 2), "sec\n"))
      }
    }
  }
}

