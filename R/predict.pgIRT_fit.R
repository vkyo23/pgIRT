#' Predict function for `pgIRT_fit`
#' @param object an object class of \code{pgIRT_fit}
#' @param ... parameters to \code{predict()}.
#' @param type a string, one of "prob" or "response". "prob" returns predicted probabilities of choices 
#' and "response" returns predicted value of Y.
#' @rdname predict.pgIRT_fit
#' @importFrom Rcpp sourceCpp
#' @export

predict.pgIRT_fit <- function(object, 
                              type = 'prob',
                              ...) {
  if (length(type) > 1) stop('`type` must be one of "prob" or "response".')
  if (!type %in% c('prob', 'response')) stop('`type` must be one of "prob" or "response".')
  if (object$model == 'default') {
    bs <- 1
  } else {
    bs <- object$input$dyn_options$bill_session - 1
  }
  out <- prediction(alpha = object$parameter$alpha,
                    beta = object$parameter$beta,
                    theta = object$parameter$theta,
                    max_cat = object$input$max_cat,
                    bill_session = bs,
                    model = object$model,
                    type = type)
  if (type == 'prob' & dim(out[[1]])[3] == 1) {
    out <- out[[1]][, , 1]
  } else {
    out <- out[[1]]
  }
  return(out)
}
