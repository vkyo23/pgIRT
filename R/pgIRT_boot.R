#' @title Parametric bootstrap for pgIRT
#' @description \link{pgIRT_boot} implements paramtric bootstrap for pgIRT model to compute statistical uncertainty.
#' 
#' @param fit pgIRT object from `pgIRT()`.
#' @param boot an integer, number of bootstrap iterations. Default is 100.
#' @param verbose an integer, the function prints the status every `verbose`. Default is NULL.
#' 
#' @importFrom stats predict
#' @importFrom Rcpp sourceCpp
#' 
#' @useDynLib pgIRT, .registration = TRUE
#' @export
#' 
#' @examples 
#' \dontrun{
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
#' boot <- pgIRT_boot(pgfit_mlt, boot = 100, verbose = 10)
#' summary(boot)
#' }

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

  md <- ifelse(fit$model == 'default', 'Default model', 'Dynamic model')

  cat("================================================================\n")
  cat("Parametric Bootstrap for pgIRT (", md, ")\n")
  cat("================================================================\n")

  stime <- proc.time()[3]

  if (fit$model == "default") {
    theta_boot <- alpha_boot <- beta_boot <- list()
    for (b in 1:boot) {
      set.seed(b)
      predY <- predict(fit, type = 'response')
      predY[is.na(fit$input$data)] <- NA
      rownames(predY) <- rownames(fit$input$data)
      fit_boot <- quiet(pgIRT(data = predY, 
                              model = fit$model,
                              std = fit$control$std,
                              init = fit$input$init,
                              prior = fit$input$prior,
                              constraint = fit$input$constraint,
                              dyn_options = NULL,
                              maxit = fit$control$maxit,
                              tol = fit$control$tol,
                              verbose = NULL))
      if (!fit_boot$converge) {
        cat("NOTE: Bootstrap", b, "failed to converge...\n")
      }

      theta_boot[[b]] <- fit_boot$parameter$theta
      alpha_boot[[b]] <- fit_boot$parameter$alpha
      beta_boot[[b]] <- fit_boot$parameter$beta

      if (b %% verbose == 0) {
        cat("Boostrap", b, "DONE :", round(proc.time()[3] - stime, 1), "sec\n")
      }
    }
    L <- list(theta = theta_boot, alpha = alpha_boot, beta = beta_boot, input = fit)
    class(L) <- c("pgIRT_boot")
    return(L)
  } else {
    theta_boot <- alpha_boot <- beta_boot <- list()
    for (b in 1:boot) {
      set.seed(b)
      predY <- predict(fit, type = 'response')
      predY[is.na(fit$input$data)] <- NA
      rownames(predY) <- rownames(fit$input$data)
      fit_boot <- quiet(pgIRT(data = predY, 
                              model = fit$model,
                              std = fit$control$std,
                              init = fit$input$init,
                              prior = fit$input$prior,
                              constraint = fit$input$constraint,
                              dyn_options = fit$input$dyn_options,
                              maxit = fit$control$maxit,
                              tol = fit$control$tol,
                              verbose = NULL))
      if (!fit_boot$converge) {
        cat("NOTE: Bootstrap", b, "failed to converge...\n")
      }

      theta_boot[[b]] <- fit_boot$parameter$theta
      alpha_boot[[b]] <- fit_boot$parameter$alpha
      beta_boot[[b]] <- fit_boot$parameter$beta

      if (b %% verbose == 0) {
        cat("Boostrap", b, "DONE :", round(proc.time()[3] - stime, 1), "sec\n")
      }
    }
    L <- list(theta = theta_boot, alpha = alpha_boot, beta = beta_boot, input = fit)
    class(L) <- c("pgIRT_boot")
    return(L)
  }
}
