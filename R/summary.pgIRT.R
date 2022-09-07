#' Print function for EM
#' @param object An output of \code{summary.pgIRT_fit()} function.
#' @export
print.summary.pgIRT_fit <- function(object) {
  # print(obj)
  # print_tbl <- getFromNamespace("print.tbl", "tibble")
  # print_tbl(obj)
  div <- paste(rep('=', 15), collapse = '') 
  for (p in 1:length(object)) {
    cat(div, 'Parameter =', names(object)[p], div, '\n')
    print(object[[p]])
  }
  invisible(object)
}

#' Print function for EM boot
#' @param object An output of \code{summary.pgIRT_boot()} function.
#' @export
print.summary.pgIRT_boot <- function(object) {
  # print(obj)
  # print_tbl <- getFromNamespace("print.tbl", "tibble")
  # print_tbl(obj)
  div <- paste(rep('=', 20), collapse = '') 
  for (p in 1:length(object)) {
    cat(div, 'Parameter =', names(object)[p], div, '\n')
    print(object[[p]])
  }
  invisible(object)
}

#' Summary function for EM
#' @param fit an object class of \code{"pgIRT_fit"}
#' @param parameter a string or string vector, indicating what parameters you want to get. 
#' @importFrom dplyr %>% mutate select as_tibble rename tibble filter bind_rows
#' @importFrom tidyr pivot_longer
#' @export
summary.pgIRT_fit <- function(fit, 
                              parameter = c('alpha', 'beta', 'theta')
) {
  
  if (fit$model %in% c('bin', 'multi')) {
    item <- dplyr::tibble(variable = NA, bill_id = NA, estimate = NA)
    ind <- dplyr::tibble(variable = NA, unit_id = NA, estimate = NA)
  } else {
    item <- dplyr::tibble(variable = NA, bill_id = NA, session = NA, estimate = NA)
    ind <- dplyr::tibble(variable = NA, unit_id = NA, session = NA, estimate = NA)
  }
  is_item <- FALSE
  is_ind <- FALSE
  if (any(parameter == 'alpha')) {
    is_item <- TRUE
    al <- fit$parameter$alpha
    if (fit$model %in% c('multi', 'multi_dyn')) {
      bill_id <- rownames(al) %>% 
        as.numeric()
    } else {
      bill_id <- names(al) %>% 
        as.numeric()
    }
    al <- dplyr::as_tibble(al) %>% 
      suppressWarnings()
    if (fit$model %in% c('bin_dyn', 'multi_dyn')) {
      bill_session <- fit$input$dyn_options$bill_session
      al <- al %>% 
        dplyr::mutate(session = bill_session + 1)
    }
    
    if (fit$model == 'multi_dyn') {
      al1 <- al[, 1] %>% 
        dplyr::mutate(variable = 'alpha1', bill_id = bill_id, session = al$session, .before = alpha1) %>% 
        dplyr::rename(estimate = alpha1)
      al2 <- al[, 2] %>% 
        dplyr::mutate(variable = 'alpha2', bill_id = bill_id, session = al$session, .before = alpha2) %>% 
        dplyr::rename(estimate = alpha2)
      item <- dplyr::bind_rows(item, al1, al2)
    } else if (fit$model == 'bin_dyn') {
      al <- al %>% 
        dplyr::mutate(variable = 'alpha', bill_id = bill_id, session, .before = value) %>% 
        dplyr::rename(estimate = value)
      item <- dplyr::bind_rows(item, al)
    } else if (fit$model == 'multi') {
      al1 <- al[, 1] %>% 
        dplyr::mutate(variable = 'alpha1', bill_id = bill_id, .before = alpha1) %>% 
        dplyr::rename(estimate = alpha1)
      al2 <- al[, 2] %>% 
        dplyr::mutate(variable = 'alpha2', bill_id = bill_id, .before = alpha2) %>% 
        dplyr::rename(estimate = alpha2)
      item <- dplyr::bind_rows(item, al1, al2)
    } else if (fit$model == 'bin') {
      al <- al %>% 
        dplyr::mutate(variable = 'alpha', bill_id = bill_id, .before = value) %>% 
        dplyr::rename(estimate = value)
      item <- dplyr::bind_rows(item, al)
    }
  } 
  if (any(parameter == 'beta')) {
    is_item <- TRUE
    be <- fit$parameter$beta
    if (fit$model %in% c('multi', 'multi_dyn')) {
      bill_id <- rownames(be) %>% 
        as.numeric()
    } else {
      bill_id <- names(be) %>% 
        as.numeric()
    }
    be <- dplyr::as_tibble(be) %>% 
      suppressWarnings()
    if (fit$model %in% c('bin_dyn', 'multi_dyn')) {
      bill_session <- fit$input$dyn_options$bill_session
      be <- be %>% 
        dplyr::mutate(session = bill_session + 1)
    }
    
    if (fit$model == 'multi_dyn') {
      be1 <- be[, 1] %>% 
        dplyr::mutate(variable = 'beta1', bill_id = bill_id, session = be$session, .before = beta1) %>% 
        dplyr::rename(estimate = beta1)
      be2 <- be[, 2] %>% 
        dplyr::mutate(variable = 'beta2', bill_id = bill_id, session = be$session, .before = beta2) %>% 
        dplyr::rename(estimate = beta2)
      item <- dplyr::bind_rows(item, be1, be2)
    } else if (fit$model == 'bin_dyn') {
      be <- be %>% 
        dplyr::mutate(variable = 'beta', bill_id = bill_id, session, .before = value) %>% 
        dplyr::rename(estimate = value)
      item <- dplyr::bind_rows(item, be)
    } else if (fit$model == 'multi') {
      be1 <- be[, 1] %>% 
        dplyr::mutate(variable = 'beta1', bill_id = bill_id, .before = beta1) %>% 
        dplyr::rename(estimate = beta1)
      be2 <- be[, 2] %>% 
        dplyr::mutate(variable = 'beta2', bill_id = bill_id, .before = beta2) %>% 
        dplyr::rename(estimate = beta2)
      item <- dplyr::bind_rows(item, be1, be2)
    } else if (fit$model == 'bin') {
      be <- be %>% 
        dplyr::mutate(variable = 'beta', bill_id = bill_id, .before = value) %>% 
        dplyr::rename(estimate = value)
      item <- dplyr::bind_rows(item, be)
    }
  } 
  if (any(parameter == 'theta')) {
    is_ind <- TRUE
    th <- fit$parameter$theta
    if (fit$model %in% c('bin', 'multi')) {
      unit_id <- names(th)
      th <- dplyr::tibble(variable = 'theta', unit_id = as.numeric(unit_id), estimate = th)
    } else {
      unit_id <- rownames(th)
      th <- th %>% 
        dplyr::as_tibble() %>% 
        suppressWarnings() %>% 
        dplyr::mutate(unit_id = as.numeric(unit_id)) %>% 
        tidyr::pivot_longer(-unit_id) %>% 
        dplyr::rename(session = name, estimate = value) %>% 
        dplyr::mutate(variable = 'theta', session = as.numeric(session) + 1,
                      unit_id = as.numeric(unit_id)) %>% 
        dplyr::select(variable, unit_id, session, estimate)
    }
    ind <- dplyr::bind_rows(ind, th)
  }
  
  L <- list()
  if (is_item) {
    item <- item %>% 
      dplyr::filter(!is.na(variable))
    if (any(parameter == 'alpha')) L$alpha <- item[grep('alpha', item$variable), ]
    if (any(parameter == 'beta')) L$beta <- item[grep('beta', item$variable), ]
  }
  if (is_ind) {
    ind <- ind %>% 
      dplyr::filter(!is.na(variable))
    L$theta <- ind
  }
  class(L) <- c('summary.pgIRT_fit', class(fit))
  return(L)
}

#' Summary function for EM boot
#' @param boot_fit an object class of \code{"pgIRT_boot"}
#' @param parameter a string or string vector, indicating what parameters you want to get. 
#' @param ci a float (< 1). The function returns \code{ci} * 100 % confidence interval.
#' @importFrom dplyr %>% mutate select as_tibble rename tibble filter
#' @importFrom tidyr pivot_longer
#' @importFrom stats quantile
#' @export
summary.pgIRT_boot <- function(boot_fit, 
                               parameter = c('alpha', 'beta', 'theta'), 
                               ci = 0.95) {
  
  if (ci > 1) stop('CI should be smaller than 1.')
  lwr_q <- (1 - ci) / 2
  upr_q <- 1 - (1 - ci) / 2
  stat_name <- paste0(ci * 100, '%')
  if (boot_fit$input$model %in% c('bin', 'multi')) {
    item <- dplyr::tibble(variable = NA, bill_id = NA, ci = NA, estimate = NA, lwr = NA, upr = NA)
    ind <- dplyr::tibble(variable = NA, unit_id = NA, ci = NA, estimate = NA, lwr = NA, upr = NA)
  } else {
    item <- dplyr::tibble(variable = NA, bill_id = NA, session = NA, ci = NA, estimate = NA, lwr = NA, upr = NA)
    ind <- dplyr::tibble(variable = NA, unit_id = NA, session = NA, ci = NA, estimate = NA, lwr = NA, upr = NA)
  }
  is_item <- FALSE
  is_ind <- FALSE
  if (any(parameter == 'alpha')) {
    is_item <- TRUE
    if (boot_fit$input$model %in% c('bin', 'bin_dyn')) {
      al_store <- rep()
      for (b in 1:length(boot_fit$alpha)) {
        al_b <- boot_fit$alpha[[b]]
        al_store <- cbind(al_store, al_b)
      }
      al_bias <- boot_fit$input$parameter$alpha - apply(al_store, 1, mean, na.rm = TRUE)
      al_lwr <- apply(al_store, 1, stats::quantile, probs = lwr_q, na.rm = TRUE)
      al_upr <- apply(al_store, 1, stats::quantile, probs = upr_q, na.rm = TRUE)
      bill_id <- names(boot_fit$alpha[[1]]) %>% as.numeric()
      tmp <- dplyr::tibble(variable = 'alpha', bill_id, ci = stat_name, estimate = boot_fit$input$parameter$alpha, 
                           lwr = al_lwr, upr = al_upr)
      if (boot_fit$input$model == 'bin_dyn') {
        session <- boot_fit$input$input$dyn_options$bill_session + 1
        tmp <- tmp %>% 
          dplyr::mutate(session = session, .before = ci)
      }
      item <- dplyr::bind_rows(item, tmp)
    } else {
      al1_store <- rep()
      al2_store <- rep()
      for (b in 1:length(boot_fit$alpha)) {
        al_b <- boot_fit$alpha[[b]]
        al1_store <- cbind(al1_store, al_b[, 1])
        al2_store <- cbind(al2_store, al_b[, 2])
      }
      al1_bias <- boot_fit$input$parameter$alpha[, 1] - apply(al1_store, 1, mean, na.rm = TRUE)
      al2_bias <- boot_fit$input$parameter$alpha[, 2] - apply(al2_store, 1, mean, na.rm = TRUE)
      al1_lwr <- apply(al1_store, 1, stats::quantile, probs = lwr_q, na.rm = TRUE) + al1_bias
      al1_upr <- apply(al1_store, 1, stats::quantile, probs = upr_q, na.rm = TRUE) + al1_bias
      al2_lwr <- apply(al2_store, 1, stats::quantile, probs = lwr_q, na.rm = TRUE) + al2_bias
      al2_upr <- apply(al2_store, 1, stats::quantile, probs = upr_q, na.rm = TRUE) + al2_bias
      bill_id <- rownames(boot_fit$alpha[[1]]) %>% as.numeric()
      tmp <- dplyr::tibble(variable = 'alpha1', bill_id,  ci = stat_name, estimate = boot_fit$input$parameter$alpha[, 1], 
                           lwr = al1_lwr, upr = al1_upr)
      tmp2 <- dplyr::tibble(variable = 'alpha2', bill_id, ci = stat_name, estimate = boot_fit$input$parameter$alpha[, 2],
                            lwr = al2_lwr, upr = al2_upr)
      if (boot_fit$input$model == 'multi_dyn') {
        ses <- boot_fit$input$input$dyn_options$bill_session + 1
        tmp <- tmp %>% 
          dplyr::mutate(session = ses, .before = ci)
        tmp2 <- tmp2 %>% 
          dplyr::mutate(session = ses, .before = ci)
      }
      item <- dplyr::bind_rows(item, tmp, tmp2)
    }
  }
  
  if (any(parameter == 'beta')) {
    is_item <- TRUE
    if (boot_fit$input$model %in% c('bin', 'bin_dyn')) {
      be_store <- rep()
      for (b in 1:length(boot_fit$beta)) {
        be_b <- boot_fit$beta[[b]]
        be_store <- cbind(be_store, be_b)
      }
      be_bias <- boot_fit$input$parameter$beta - apply(be_store, 1, mean, na.rm = TRUE)
      be_lwr <- apply(be_store, 1, stats::quantile, probs = lwr_q, na.rm = TRUE)
      be_upr <- apply(be_store, 1, stats::quantile, probs = upr_q, na.rm = TRUE)
      bill_id <- names(boot_fit$beta[[1]]) %>% as.numeric()
      tmp <- dplyr::tibble(variable = 'beta', bill_id, ci = stat_name, estimate = boot_fit$input$parameter$beta, 
                           lwr = be_lwr, upr = be_upr)
      if (boot_fit$input$model == 'bin_dyn') {
        session <- boot_fit$input$input$dyn_options$bill_session + 1
        tmp <- tmp %>% 
          dplyr::mutate(session = session, .before = ci)
      }
      item <- dplyr::bind_rows(item, tmp)
    } else {
      be1_store <- rep()
      be2_store <- rep()
      for (b in 1:length(boot_fit$beta)) {
        be_b <- boot_fit$beta[[b]]
        be1_store <- cbind(be1_store, be_b[, 1])
        be2_store <- cbind(be2_store, be_b[, 2])
      }
      be1_bias <- boot_fit$input$parameter$beta[, 1] - apply(be1_store, 1, mean, na.rm = TRUE)
      be2_bias <- boot_fit$input$parameter$beta[, 2] - apply(be2_store, 1, mean, na.rm = TRUE)
      be1_lwr <- apply(be1_store, 1, stats::quantile, probs = lwr_q, na.rm = TRUE) + be1_bias
      be1_upr <- apply(be1_store, 1, stats::quantile, probs = upr_q, na.rm = TRUE) + be1_bias
      be2_lwr <- apply(be2_store, 1, stats::quantile, probs = lwr_q, na.rm = TRUE) + be2_bias
      be2_upr <- apply(be2_store, 1, stats::quantile, probs = upr_q, na.rm = TRUE) + be2_bias
      bill_id <- rownames(boot_fit$beta[[1]]) %>% as.numeric()
      tmp <- dplyr::tibble(variable = 'beta1', bill_id,  ci = stat_name, estimate = boot_fit$input$parameter$beta[, 1], 
                           lwr = be1_lwr, upr = be1_upr)
      tmp2 <- dplyr::tibble(variable = 'beta2', bill_id, ci = stat_name, estimate = boot_fit$input$parameter$beta[, 2],
                            lwr = be2_lwr, upr = be2_upr)
      if (boot_fit$input$model == 'multi_dyn') {
        ses <- boot_fit$input$input$dyn_options$bill_session + 1
        tmp <- tmp %>% 
          dplyr::mutate(session = ses, .before = ci)
        tmp2 <- tmp2 %>% 
          dplyr::mutate(session = ses, .before = ci)
      }
      item <- dplyr::bind_rows(item, tmp, tmp2)
    }
  }
  
  if (any(parameter == 'theta')) {
    is_ind <- TRUE
    if (boot_fit$input$model %in% c('bin', 'multi')) {
      th_store <- rep()
      uid <- names(boot_fit$input$parameter$theta) %>% as.numeric()
      for (b in 1:length(boot_fit$theta)) {
        th_store <- cbind(th_store, boot_fit$theta[[b]])
      }
      th_bias <- boot_fit$input$parameter$theta - apply(th_store, 1, mean, na.rm = TRUE)
      th_lwr <- apply(th_store, 1, stats::quantile, probs = lwr_q, na.rm = TRUE) + th_bias
      th_upr <- apply(th_store, 1, stats::quantile, probs = upr_q, na.rm = TRUE) + th_bias
      ind <- dplyr::bind_rows(ind,
                              dplyr::tibble(variable = 'theta', unit_id = uid, ci = stat_name,
                                            estimate = boot_fit$input$parameter$theta, lwr = th_lwr, upr = th_upr))
    } else {
      for (i in 1:nrow(boot_fit$input$parameter$theta)) {
        th_store <- rep()
        uid <- rownames(boot_fit$theta[[1]])[i] %>% as.numeric()
        for (b in 1:length(boot_fit$theta)) {
          th_store <- cbind(th_store, boot_fit$theta[[b]][i, ])
        }
        th_bias <- boot_fit$input$parameter$theta[i, ] - apply(th_store, 1, mean, na.rm = TRUE)
        th_lwr <- apply(th_store, 1, stats::quantile, probs = lwr_q, na.rm = TRUE) + th_bias
        th_upr <- apply(th_store, 1, stats::quantile, probs = upr_q, na.rm = TRUE) + th_bias
        ind <- dplyr::bind_rows(ind,
                                dplyr::tibble(variable = 'theta', unit_id = uid, session = 1:length(th_upr), ci = stat_name,
                                              estimate = boot_fit$input$parameter$theta[i, ], lwr = th_lwr, upr = th_upr))
      }
    }
  }
  
  L <- list()
  if (is_item) {
    item <- item %>% 
      dplyr::filter(!is.na(variable))
    if (any(parameter == 'alpha')) L$alpha <- item[grep('alpha', item$variable), ]
    if (any(parameter == 'beta')) L$beta <- item[grep('beta', item$variable), ]
  }
  if (is_ind) {
    ind <- ind %>% 
      dplyr::filter(!is.na(variable))
    L$theta <- ind
  }
  class(L) <- c('summary.pgIRT_boot', class(boot_fit), class(boot_fit$input))
  return(L)
}
