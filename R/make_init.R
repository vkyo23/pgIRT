#' @title Generating initial values for pgIRT model
#' @description \code{make_init} is an auxiliary function for pgIRT model and generates optimal initial values based on roll-call matrix.
#'
#' @param data a matrix object, roll-call matrix
#' @param model string, one of "bin" (binary), "bin_dyn" (dynamic binary), "multi" (multinomial) or "multi_dyn" (dynamic multinomial)
#' @param T an integer, indicating the number of sessions. For dynamic model, you must specify this argument. 
#' @param constraint an integer or integer vector (for dynamic model, T-length), indicating the voter whose ideal point is always set positive.
#' 
#' @importFrom stats qlogis na.omit
#' @importFrom dplyr %>% as_tibble tibble mutate filter right_join select group_by summarise ungroup left_join bind_cols if_else n
#' @importFrom tidyr pivot_longer pivot_wider expand_grid
#' @export
#' @examples
#' \dontrun{
#' data(m_data_dyn)
#' mat <- make_rollcall(m_data_dyn,
#'                      unit = "unit",
#'                      bill = "bill",
#'                      vote_col = "vote")
#' mat <- clean_rollcall(mat)
#' init <- make_init(mat,
#'                   model = "default",
#'                   constraint = 1)
#'                   
#' }

make_init <- function(data, 
                      model = c('default', 'dynamic'),
                      T = NULL, 
                      constraint = NULL) {
  if (!class(data)[1] %in% c("matrix", "data.frame", "tbl_df")) {
    stop("`data` must be 'matrix', 'data.frame' or 'tbl_df' object.")
  } else {
    data <- as.matrix(data)
  }
  if (length(model) > 1) {
    stop("Please specify `model` argument. Only one of 'default' or 'dynamic' is allowed.")
  } else {
    if (!model %in% c('default', 'dynamic')) {
      stop("Please specify `model` argument. Only one of 'default' or 'dynamic' is allowed.")
    } 
  }
  if (is.null(T) & model == 'dynamic') stop('Please supply the number of sessions `T`. For default model, T = 1.')
  if (model == 'default') T <- 1
  if (is.null(constraint)) {
    constraint <- 1
    if (model == 'dynamic') constraint <- rep(1, I)
  } else {
    if (model == 'default') {
      constraint <- constraint
      if (length(constraint) > 1) stop('`constraint` must be a scalar, not a vector.')
    } else {
      if (length(constraint) == 1) {
        constraint <- rep(constraint, T)
      } else if (length(constraint) != T) {
        stop("`constraint` must be the same length as the number of sessions.")
      }
    }
  }
  if (is.null(rownames(data))) rownames(data) <- 1:nrow(data)
  if (is.null(colnames(data))) colnames(data) <- 1:ncol(data)
  idd <- dplyr::tibble(id_num = 1:ncol(data), name = colnames(data))
  temp_d <- data %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(id = rownames(data)) %>%
    tidyr::pivot_longer(-id) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::right_join(idd, by = "name")
  
  options(dplyr.summarise.inform = FALSE)
  if (length(unique(na.omit(as.vector(data)))) != 2) {
    sb_calc <- temp_d %>%
      dplyr::group_by(id_num, value) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(id_num) %>%
      dplyr::summarise(vote = value, count = count, all = sum(count),
                       ref_cat = max(value)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(id_num, vote) %>%
      dplyr::summarise(count = count, probs = count / all, all = all,
                       ref_cat = ref_cat) %>%
      dplyr::group_by(id_num) %>%
      dplyr::summarise(vote = vote, count = count,
                       probs = probs, cumsum = cumsum(count) / all,
                       sb = probs / (1 - (cumsum - probs)),
                       ref_cat = ref_cat) %>%
      dplyr::filter(vote != ref_cat) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sb_value = stats::qlogis(sb)) %>%
      dplyr::select(id_num, vote, sb_value, ref_cat)
  } else {
    sb_calc <- temp_d %>%
      dplyr::group_by(id_num, value) %>%
      dplyr::summarise(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(id_num) %>%
      dplyr::summarise(vote = value, count = count, all = sum(count),
                       ref_cat = max(value)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(id_num, vote) %>%
      dplyr::summarise(count = count, probs = count / all, all = all,
                       ref_cat = ref_cat) %>%
      dplyr::distinct(id_num, .keep_all = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sb_value = stats::qlogis(probs),
                    sb_value = if_else(is.na(sb_value) | is.infinite(sb_value), 
                                       0.1,
                                       sb_value),
                    vote = 1) %>%
      dplyr::select(id_num, vote, sb_value, ref_cat)
  }
  options(dplyr.summarise.inform = TRUE)
  item_rows <- max(sb_calc$vote)
  alpha_init <- tidyr::expand_grid(id_num = 1:ncol(data),
                                   vote = 1:item_rows) %>% 
    dplyr::left_join(sb_calc, by = c('id_num', 'vote')) %>% 
    dplyr::select(-ref_cat) %>% 
    tidyr::pivot_wider(names_from = vote, values_from = sb_value) %>% 
    dplyr::select(-id_num) %>% 
    as.matrix()
  colnames(alpha_init) <- paste0('alpha_', 1:ncol(alpha_init))
  theta_init <- scale(-rowMeans(data, na.rm = TRUE))[, 1]
  theta_init <- matrix(rep(theta_init, T), length(theta_init), T)
  theta_init <- ifelse(is.na(theta_init), 0.1, theta_init)
  if (model == 'default') {
    if (theta_init[constraint, 1] < 0) theta_init <- -theta_init
  } else {
    for (t in 1:T) {
      if (theta_init[constraint[t], t] < 0) theta_init[, t] <- -theta_init[, t]
    }
  }
  colnames(theta_init) <- paste0('theta_', 1:T)
  beta_init <- rep(rep(0.1, ncol(data)), ncol(alpha_init)) %>% 
    matrix(ncol(data), ncol(alpha_init))
  beta_init[is.na(alpha_init)] <- NA
  colnames(beta_init) <- gsub('alpha', 'beta', colnames(alpha_init))
  L <- list(alpha = alpha_init,
            beta = beta_init,
            theta = theta_init)
  return(L)
}
