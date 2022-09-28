#' Print function for EM
#' @param x An output of \code{summary.pgIRT_fit()} function.
#' @param ... parameters to \code{print()}.
#' @rdname print.summary.pgIRT_fit
#' @export
print.summary.pgIRT_fit <- function(x, ...) {
  div <- paste(rep('=', 15), collapse = '') 
  for (p in 1:length(x)) {
    cat(div, 'Parameter =', names(x)[p], div, '\n')
    print(x[[p]])
  }
  invisible(x)
}

#' Print function for EM boot
#' @param x An output of \code{summary.pgIRT_boot()} function.
#' @param ... parameters to \code{print()}.
#' @rdname print.summary.pgIRT_fit
#' @export
print.summary.pgIRT_boot <- function(x, ...) {
  div <- paste(rep('=', 20), collapse = '') 
  for (p in 1:length(x)) {
    cat(div, 'Parameter =', names(x)[p], div, '\n')
    print(x[[p]])
  }
  invisible(x)
}

#' Summary function for EM
#' @param object an object class of \code{"pgIRT_fit"}
#' @param parameter a string or string vector, indicating what parameters you want to get. 
#' @param ... parameters to \code{summary()}.
#' @importFrom dplyr %>% mutate select as_tibble rename tibble filter bind_rows
#' @importFrom tidyr pivot_longer
#' @rdname summary.pgIRT_fit
#' @export
summary.pgIRT_fit <- function(object, 
                              parameter = c('alpha', 'beta', 'theta'),
                              ...) {
  
  if (object$model == 'default') {
    item <- dplyr::tibble(bill = NA, variable = NA, estimate = NA)
    ind <- dplyr::tibble(unit = NA, variable = NA, estimate = NA)
  } else {
    item <- dplyr::tibble(bill = NA, session = NA, variable = NA, estimate = NA)
    ind <- dplyr::tibble(unit = NA, session = NA, variable = NA, estimate = NA)
  }
  is_item <- FALSE
  is_ind <- FALSE
  if (any(parameter == 'alpha')) {
    is_item <- TRUE
    al <- object$parameter$alpha
    if (object$model == 'default') {
      tmp_al <- as.data.frame(al) %>% 
        dplyr::mutate(bill = rownames(al),
                      .before = 'alpha_1') %>% 
        tidyr::pivot_longer(-bill) %>% 
        dplyr::rename(variable = name,
                      estimate = value) %>% 
        dplyr::filter(!is.na(estimate))
    } else {
      tmp_al <- as.data.frame(al) %>% 
        dplyr::mutate(bill = rownames(al),
                      session = object$input$dyn_options$bill_session,
                      .before = 'alpha_1') %>% 
        tidyr::pivot_longer(-c(bill, session)) %>% 
        dplyr::rename(variable = name,
                      estimate = value) %>% 
        dplyr::filter(!is.na(estimate))
    }
    item <- bind_rows(item, tmp_al)
  } 
  if (any(parameter == 'beta')) {
    is_item <- TRUE
    be <- object$parameter$beta
    if (object$model == 'default') {
      tmp_be <- as.data.frame(be) %>% 
        dplyr::mutate(bill = rownames(be),
                      .before = 'beta_1') %>% 
        tidyr::pivot_longer(-bill) %>% 
        dplyr::rename(variable = name,
                      estimate = value) %>% 
        dplyr::filter(!is.na(estimate))
    } else {
      tmp_be <- as.data.frame(be) %>% 
        dplyr::mutate(bill = rownames(be),
                      session = object$input$dyn_options$bill_session,
                      .before = 'beta_1') %>% 
        tidyr::pivot_longer(-c(bill, session)) %>% 
        dplyr::rename(variable = name,
                      estimate = value) %>% 
        dplyr::filter(!is.na(estimate))
    }
    item <- bind_rows(item, tmp_be)
  } 
  if (any(parameter == 'theta')) {
    is_ind <- TRUE
    th <- object$parameter$theta
    if (object$model == 'default') {
      tmp_th <- as.data.frame(th) %>% 
        mutate(unit = rownames(th),
               .before = theta_1) %>% 
        pivot_longer(-unit) %>% 
        dplyr::rename(variable = name,
                      estimate = value) %>% 
        dplyr::mutate(variable = gsub('_\\d+', '', variable))
    } else {
      tmp_th <- as.data.frame(th) %>% 
        mutate(unit = rownames(th),
               .before = theta_1) %>% 
        pivot_longer(-unit) %>% 
        dplyr::rename(variable = name,
                      estimate = value) %>% 
        dplyr::mutate(session = rep(1:ncol(th), nrow(th)),
                      .after = unit) %>% 
        dplyr::mutate(variable = gsub('_\\d+', '', variable)) %>% 
        dplyr::filter(!is.na(estimate))
    }
    ind <- dplyr::bind_rows(ind, tmp_th)
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
  class(L) <- c('summary.pgIRT_object', class(object))
  return(L)
}

#' Summary function for EM boot
#' @param object an object class of \code{"pgIRT_boot"}
#' @param parameter a string or string vector, indicating what parameters you want to get. 
#' @param ci a float (< 1). The function returns \code{ci} * 100 % confidence interval.
#' @param ... parameters to \code{summary()}.
#' @importFrom dplyr %>% mutate select as_tibble rename tibble filter bind_rows
#' @importFrom tidyr pivot_longer
#' @importFrom stats quantile
#' @rdname summary.pgIRT_boot
#' @export
summary.pgIRT_boot <- function(object, 
                               parameter = c('alpha', 'beta', 'theta'), 
                               ci = 0.95,
                               ...) {
  # Auxiliary function
  get_q <- function(x, pr) apply(x, 1, quantile, probs = pr, na.rm = TRUE)
  if (ci > 1) stop('CI should be smaller than 1.')
  lwr_q <- (1 - ci) / 2
  upr_q <- 1 - (1 - ci) / 2
  if (object$input$model == 'default') {
    item <- dplyr::tibble(bill = NA, variable = NA, ci = NA, lwr = NA, estimate = NA, upr = NA)
    ind <- dplyr::tibble(unit = NA, variable = NA, ci = NA, lwr = NA, estimate = NA, upr = NA)
  } else {
    item <- dplyr::tibble(bill = NA, variable = NA, session = NA, ci = NA, lwr = NA, estimate = NA, upr = NA)
    ind <- dplyr::tibble(unit = NA, variable = NA, session = NA, ci = NA, lwr = NA, estimate = NA, upr = NA)
  }
  is_item <- FALSE
  is_ind <- FALSE
  if (any(parameter == 'alpha')) {
    is_item <- TRUE
    al <- object$alpha
    if (object$input$model == 'default') {
      K <- ncol(al[[1]])
      tmp_al <- do.call(cbind, al)
      for (k in 1:K) {
        alk <- tmp_al[, seq(k, ncol(tmp_al), K)]
        bias <- object$input$parameter$alpha[, k] - rowMeans(alk, na.rm = TRUE)
        alk <- (t(get_q(alk, pr = c(lwr_q, upr_q))) + bias) 
        colnames(alk) <- c('lwr', 'upr')
        alk <- alk %>% 
          dplyr::as_tibble() %>% 
          dplyr::mutate(estimate = object$input$parameter$alpha[, k],
                        .after = lwr) %>% 
          dplyr::mutate(bill = rownames(object$input$parameter$alpha),
                        variable = paste0('alpha_', k),
                        ci = paste0(ci * 100, '%'),
                        .before = lwr)
        if (k == 1) {
          alk_all <- alk 
        } else {
          alk_all <- dplyr::bind_rows(alk_all, alk)
        }
      }
    } else {
      K <- ncol(al[[1]])
      tmp_al <- do.call(cbind, al)
      for (k in 1:K) {
        alk <- tmp_al[, seq(k, ncol(tmp_al), K)]
        bias <- object$input$parameter$alpha[, k] - rowMeans(alk, na.rm = TRUE)
        alk <- (t(get_q(alk, pr = c(lwr_q, upr_q))) + bias) 
        colnames(alk) <- c('lwr', 'upr')
        alk <- alk %>% 
          dplyr::as_tibble() %>% 
          dplyr::mutate(estimate = object$input$parameter$alpha[, k],
                        .after = lwr) %>% 
          dplyr::mutate(bill = rownames(object$input$parameter$alpha),
                        variable = paste0('alpha_', k),
                        session = object$input$input$dyn_options$bill_session,
                        ci = paste0(ci * 100, '%'),
                        .before = lwr)
        if (k == 1) {
          alk_all <- alk 
        } else {
          alk_all <- dplyr::bind_rows(alk_all, alk)
        }
      }
    }
    item <- dplyr::bind_rows(item, alk_all) %>% 
      dplyr::filter(!is.na(estimate))
  }
  
  if (any(parameter == 'beta')) {
    is_item <- TRUE
    be <- object$beta
    if (object$input$model == 'default') {
      K <- ncol(al[[1]]) 
      tmp_be <- do.call(cbind, be)
      for (k in 1:K) {
        bek <- tmp_be[, seq(k, ncol(tmp_be), K)]
        bias <- object$input$parameter$beta[, k] - rowMeans(bek, na.rm = TRUE)
        bek <- (t(get_q(bek, pr = c(lwr_q, upr_q))) + bias) 
        colnames(bek) <- c('lwr', 'upr')
        bek <- bek %>% 
          dplyr::as_tibble() %>% 
          dplyr::mutate(estimate = object$input$parameter$beta[, k],
                        .after = lwr) %>% 
          dplyr::mutate(bill = rownames(object$input$parameter$beta),
                        variable = paste0('beta_', k),
                        ci = paste0(ci * 100, '%'),
                        .before = lwr)
        if (k == 1) {
          bek_all <- bek 
        } else {
          bek_all <- dplyr::bind_rows(bek_all, bek)
        }
      }
    } else {
      K <- ncol(be[[1]])
      tmp_be <- do.call(cbind, be)
      for (k in 1:K) {
        bek <- tmp_be[, seq(k, ncol(tmp_be), K)]
        bias <- object$input$parameter$beta[, k] - rowMeans(bek, na.rm = TRUE)
        bek <- (t(get_q(bek, pr = c(lwr_q, upr_q))) + bias) 
        colnames(bek) <- c('lwr', 'upr')
        bek <- bek %>% 
          dplyr::as_tibble() %>% 
          dplyr::mutate(estimate = object$input$parameter$beta[, k],
                        .after = lwr) %>% 
          dplyr::mutate(bill = rownames(object$input$parameter$beta),
                        variable = paste0('beta_', k),
                        session = object$input$input$dyn_options$bill_session,
                        ci = paste0(ci * 100, '%'),
                        .before = lwr)
        if (k == 1) {
          bek_all <- bek 
        } else {
          bek_all <- dplyr::bind_rows(bek_all, bek)
        }
      }
    }
    item <- dplyr::bind_rows(item, bek_all) %>% 
      dplyr::filter(!is.na(estimate))
  }
  
  if (any(parameter == 'theta')) {
    is_ind <- TRUE
    th <- object$th
    TT <- ncol(th[[1]])
    tmp_th <- do.call(cbind, th)
    for (t in 1:TT) {
      tht <- tmp_th[, seq(t, ncol(tmp_th), TT)]
      bias <- object$input$parameter$theta[, t] - rowMeans(tht, na.rm = TRUE)
      tht <- (t(get_q(tht, pr = c(lwr_q, upr_q))) + bias) 
      colnames(tht) <- c('lwr', 'upr')
      tht <- tht %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(estimate = object$input$parameter$theta[, t],
                      .after = lwr) %>% 
        dplyr::mutate(unit = rownames(object$input$parameter$theta),
                      variable = 'theta',
                      session = t,
                      ci = paste0(ci * 100, '%'),
                      .before = lwr)
      if (t == 1) {
        tht_all <- tht
      } else {
        tht_all <- dplyr::bind_rows(tht_all, tht)
      }
    }
    if (object$input$model == 'default') tht_all <- dplyr::select(tht_all, -session)
    ind <- dplyr::bind_rows(ind, tht_all) %>% 
      dplyr::filter(!is.na(estimate))
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
  class(L) <- c('summary.pgIRT_boot', class(object), class(object$input))
  return(L)
}
