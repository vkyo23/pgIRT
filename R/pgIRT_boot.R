#' @title Parametric bootstrap for pgIRT
#' @description \link{pgIRT_boot} implements paramtric bootstrap for pgIRT model to compute statistical uncertainty.
#' 
#' @param fit pgIRT object from `pgIRT()`.
#' @param boot an integer, number of bootstrap iterations. Default is 100.
#' @param verbose an integer, the function prints the status every `verbose`. Default is NULL.
#' 
#' @importFrom crayon red yellow
#' @useDynLib pgIRT, .registration = TRUE
#' @export

pgIRT_boot <- function(fit, 
                       boot = 100, 
                       verbose = NULL) {
  
  if (class(fit) != "pgIRT_fit") stop("Please supply the fit from `pgIRT`.")
  if (is.null(verbose)) verbose <- boot + 1

  # For suppressing the messages
  quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
  }

  md <- ifelse(fit$model == "bin", "Binomial",
               ifelse(fit$model == "bin_dyn", "Dynamic Binomial",
                      ifelse(fit$model == "multi", "Multinomial",
                             "Dynamic Multinomial")))

  cat("================================================================\n")
  cat("Parametric Bootstrap for pgIRT (", crayon::yellow(md), ")\n")
  cat("================================================================\n")

  stime <- proc.time()[3]

  if (fit$model == "bin") {
    theta_boot <- alpha_boot <- beta_boot <- list()
    for (b in 1:boot) {
      set.seed(b)
      predY <- Calc_predY_bin(fit$parameter$theta, fit$parameter$alpha, fit$parameter$beta)
      predY[is.na(fit$input$data)] <- NA
      fit_boot <- quiet(pgIRT(data = predY,
                              prior = fit$prior,
                              init = fit$init,
                              constraint = fit$input$constraint,
                              model = fit$model,
                              verbose = 10,
                              std = fit$control$std))
      if (!fit_boot$converge) {
        cat(crayon::red("Bootstrap", b, "failed to converge...\n"))
      }

      theta_boot[[b]] <- fit_boot$parameter$theta
      if (!is.null(names(fit$parameter$theta))) names(theta_boot[[b]]) <- names(fit$parameter$theta)
      if (!is.null(names(fit$parameter$alpha))) names(fit_boot$parameter$alpha) <- names(fit_boot$parameter$beta) <- names(fit$parameter$alpha)
      alpha_boot[[b]] <- fit_boot$parameter$alpha
      beta_boot[[b]] <- fit_boot$parameter$beta

      if (b %% verbose == 0) {
        cat("Boostrap", b, "DONE :", crayon::yellow(round(proc.time()[3] - stime, 1), "sec\n"))
      }
    }
    L <- list(theta = theta_boot, alpha = alpha_boot, beta = beta_boot, input = fit)
    class(L) <- c("pgIRT_boot")
    return(L)
  }
  if (fit$model == "bin_dyn") {
    theta_boot <- alpha_boot <- beta_boot <- list()
    for (b in 1:boot) {
      set.seed(b)
      predY <- Calc_predY_bin_dyn(fit$parameter$theta, fit$parameter$alpha, fit$parameter$beta,
                                  fit$input$dyn_options$bill_session)
      predY[is.na(fit$input$data)] <- NA
      fit_boot <- quiet(pgIRT(data = predY,
                              prior = fit$prior,
                              init = fit$init,
                              constraint = fit$input$constraint,
                              model = fit$model,
                              dyn_options = fit$input$dyn_options,
                              verbose = 10,
                              std = fit$control$std))
      if (!fit_boot$converge) {
        cat(crayon::red("Bootstrap", b, "failed to converge...\n"))
      }

      theta_boot[[b]] <- fit_boot$parameter$theta
      if (!is.null(rownames(fit$parameter$theta))) rownames(theta_boot[[b]]) <- rownames(fit$parameter$theta)
      if (!is.null(names(fit$parameter$alpha))) names(fit_boot$parameter$alpha) <- names(fit_boot$parameter$beta) <- names(fit$parameter$alpha)
      colnames(theta_boot[[b]]) <- 0:(ncol(fit$parameter$theta) - 1)
      alpha_boot[[b]] <- fit_boot$parameter$alpha
      beta_boot[[b]] <- fit_boot$parameter$beta

      if (b %% verbose == 0) {
        cat("Boostrap", b, "DONE :", crayon::yellow(round(proc.time()[3] - stime, 1), "sec\n"))
      }
    }
    L <- list(theta = theta_boot, alpha = alpha_boot, beta = beta_boot, input = fit)
    class(L) <- c("pgIRT_boot")
    return(L)
  }
  if (fit$model == "multi") {
    theta_boot <- alpha_boot <- beta_boot <- list()
    for (b in 1:boot) {
      set.seed(b)
      Y1 <- Y2 <- fit$input$data
      Y1[Y1 %in% c(2, 3)] <- 0
      Y2[Y2 == 3] <- 0; Y2[Y2 == 1] <- NA; Y2[Y2 == 2] <- 1
      predY <- Calc_predY_mlt(fit$parameter$theta, fit$parameter$alpha, fit$parameter$beta)
      y_boot <- organize_Y(predY$Y1, predY$Y2, predY$Y3)
      y_boot[is.na(fit$input$data)] <- NA
      fit_boot <- quiet(pgIRT(data = y_boot,
                              prior = fit$prior,
                              init = fit$init,
                              constraint = fit$input$constraint,
                              model = fit$model,
                              verbose = 10,
                              std = fit$control$std))
      if (!fit_boot$converge) {
        cat(crayon::red("Bootstrap", b, "failed to converge...\n"))
      }

      theta_boot[[b]] <- fit_boot$parameter$theta
      if (!is.null(rownames(fit$parameter$theta))) rownames(theta_boot[[b]]) <- rownames(fit$parameter$theta)
      if(!is.null(rownames(fit$parameter$alpha))) rownames(fit_boot$parameter$alpha) <- rownames(fit_boot$parameter$beta) <- rownames(fit$parameter$alpha)
      colnames(fit_boot$parameter$alpha) <- c('alpha1', 'alpha2')
      colnames(fit_boot$parameter$beta) <- c('beta1', 'beta2')
      alpha_boot[[b]] <- fit_boot$parameter$alpha
      beta_boot[[b]] <- fit_boot$parameter$beta

      if (b %% verbose == 0) {
        cat("Boostrap", b, "DONE :", crayon::yellow(round(proc.time()[3] - stime, 1), "sec\n"))
      }
    }
    L <- list(theta = theta_boot, alpha = alpha_boot, beta = beta_boot, input = fit)
    class(L) <- c("pgIRT_boot")
    return(L)
  }
  if (fit$model == "multi_dyn") {
    theta_boot <- alpha_boot <- beta_boot <- list()
    for (b in 1:boot) {
      set.seed(b)
      Y1 <- Y2 <- fit$input$data
      Y1[Y1 %in% c(2, 3)] <- 0
      Y2[Y2 == 3] <- 0; Y2[Y2 == 1] <- NA; Y2[Y2 == 2] <- 1
      predY <- Calc_predY_mlt_dyn(fit$parameter$theta, fit$parameter$alpha, fit$parameter$beta,
                                  fit$input$dyn_options$bill_session)
      y_boot <- organize_Y(predY$Y1, predY$Y2, predY$Y3)
      y_boot[is.na(fit$input$data)] <- NA
      fit_boot <- quiet(pgIRT(data = y_boot,
                              prior = fit$prior,
                              init = fit$init,
                              constraint = fit$input$constraint,
                              dyn_options = fit$input$dyn_options,
                              model = fit$model,
                              verbose = 10,
                              std = fit$control$std))
      if (!fit_boot$converge) {
        cat(crayon::red("Bootstrap", b, "failed to converge...\n"))
      }

      theta_boot[[b]] <- fit_boot$parameter$theta
      if (!is.null(rownames(fit$parameter$theta))) rownames(theta_boot[[b]]) <- rownames(fit$parameter$theta)
      colnames(theta_boot[[b]]) <- 0:(ncol(fit$parameter$theta) - 1)
      if(!is.null(rownames(fit$parameter$alpha))) rownames(fit_boot$parameter$alpha) <- rownames(fit_boot$parameter$beta) <- rownames(fit$parameter$alpha)
      colnames(fit_boot$parameter$alpha) <- c('alpha1', 'alpha2')
      colnames(fit_boot$parameter$beta) <- c('beta1', 'beta2')
      alpha_boot[[b]] <- fit_boot$parameter$alpha
      beta_boot[[b]] <- fit_boot$parameter$beta

      if (b %% verbose == 0) {
        cat("Boostrap", b, "DONE :", crayon::yellow(round(proc.time()[3] - stime, 1), "sec\n"))
      }
    }
    L <- list(theta = theta_boot, alpha = alpha_boot, beta = beta_boot, input = fit)
    class(L) <- c("pgIRT_boot")
    return(L)
  }
}
