#' @title Generating initial values for pgIRT model
#' @description \code{make_init} is an auxiliary function for pgIRT model and generates optimal initial values based on roll-call matrix.
#'
#' @param data a matrix object, roll-call matrix
#' @param model string, one of "bin" (binary), "bin_dyn" (dynamic binary), "multi" (multinomial) or "multi_dyn" (dynamic multinomial)
#' @param T an integer, indicating the number of sessions. If you choose "bin_dyn" or "multi_dyn" for `model`, you must supply this argument.
#' @param bill_match (Optional). (Optional). An integer vector which indicates matched bills (the location of matched bills, starts from 0) 
#' and is the same length as the number of bills. If a bill does not have matches, please fill NA. Using \link{make_bill_match} is strongly recommended.
#' @param constraint an integer or integer vector (for dynamic model, T-length), indicating the voter whose ideal point is always set positive.
#' 
#' @importFrom stats coef glm qlogis rnorm
#' @importFrom dplyr %>% as_tibble tibble mutate filter right_join select group_by summarise ungroup left_join bind_cols if_else n
#' @importFrom tidyr pivot_longer
#' @export

make_init <- function(data, 
                      model = c("bin", "bin_dyn", "multi", "multi_dyn"),
                      T = NULL, 
                      bill_match = NULL, 
                      constraint = NULL) {
  
  if (!class(data)[1] %in% c("matrix", "data.frame", "tbl_df")) {
    stop("`data` must be 'matrix', 'data.frame' or 'tbl_df' object.")
  } else {
    data <- as.matrix(data)
  }
  if (length(model) > 1) {
    stop("Please specify `model` argument. Only one of 'bin', 'bin_dyn', 'multi' or 'multi_dyn' is allowed.")
  } else {
    if (!model %in% c("bin", "bin_dyn", "multi", "multi_dyn")) {
      stop("Please specify `model` argument. Only one of 'bin', 'bin_dyn', 'multi' or 'multi_dyn' is allowed.")
    }
  }
  if (!is.null(bill_match)) {
    if (length(bill_match) != ncol(data)) stop('Size of columns of `data` does not match the size of `bill_match`.')
    if (all(colnames(data) != names(bill_match))) stop('Names of columns of `data` does not match the names of `bill_match`.')
  }
  
  I <- nrow(data)
  J <- ncol(data)
  
  if (model %in% c("bin", "bin_dyn")) {
    ALLPOINT <- sum(is.na(data)) + sum(data == 1, na.rm = TRUE) + sum(data == 0, na.rm = TRUE)
    if (ALLPOINT != I * J) stop("For binary model, the data is allowed to contain only NA, 1 and 0.")
    
    theta_init <- theta_ini <- scale(rowMeans(data, na.rm = TRUE))[, 1]
    if (model == "bin_dyn") {
      if (is.null(T)) stop("Please supply `T`!")
      theta_init <- matrix(rep(theta_init, T), ncol = T)
    }
    if (!is.null(constraint)) {
      if (model == 'bin') {
        if (theta_init[constraint] < 0) theta_init <- -theta_init
      } else if (model == 'bin_dyn') {
        if (length(constraint) == 1) {
          constraint <- rep(constraint, ncol(theta_init))
        } else if (length(constraint) != ncol(theta_init)) {
          stop('`constraint` must be the same length as the number of sessions.')
        }
        theta_init <- set_constraint(theta_init, constraint)
      }
    }
    
    alpha_init <- beta_init <- rep()
    for (j in 1:ncol(data)) {
      c <- stats::coef(stats::glm(data[, j] ~ theta_ini, family = "binomial")) %>%
        suppressWarnings()
      c <- ifelse(is.na(c) | abs(c) > 5, .1, c)
      alpha_init[j] <- c[1]
      beta_init[j] <- c[2]
    }
    tmp <- list(alpha = alpha_init, beta = beta_init, theta = theta_init)
  }
  
  
  if (model %in% c("multi", "multi_dyn")) {
    ALLPOINT <- sum(is.na(data)) + sum(data == 1, na.rm = TRUE) + sum(data == 2, na.rm = TRUE) + sum(data == 3, na.rm = TRUE)
    if (ALLPOINT != I * J) stop("For multinomial model, the data is allowed to contain only NA, 1, 2 and 3.")
    
    idd <- dplyr::tibble(id_num = 1:ncol(data), name = colnames(data))
    temp_d <- data %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(id = rownames(data)) %>%
      tidyr::pivot_longer(-id) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::right_join(idd, by = "name")
    
    options(dplyr.summarise.inform = FALSE)
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
      dplyr::filter(vote != max(vote)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(sb_value = stats::qlogis(sb)) %>%
      dplyr::select(id_num, vote, sb_value, ref_cat)
    options(dplyr.summarise.inform = TRUE)
    
    ## Assign same values to matched resolutions
    sb1 <- sb_calc %>%
      dplyr::filter(vote == 1) %>%
      dplyr::select(alpha1 = sb_value, id_num)
    sb2 <- sb_calc %>%
      dplyr::filter(vote == 2) %>%
      dplyr::select(alpha2 = sb_value, id_num)
    alpha_init <- idd %>%
      dplyr::left_join(sb1, by = "id_num") %>%
      dplyr::left_join(sb2, by = "id_num") %>%
      dplyr::select(name, alpha1, alpha2)
    nam <- alpha_init$name
    if (!is.null(bill_match)) {
      alpha_init <- alpha_init %>%
        dplyr::bind_cols(dplyr::tibble(match = as.character(bill_match))) %>%
        dplyr::left_join(alpha_init %>% 
                           dplyr::select(match = name,
                                         a1 = alpha1, a2 = alpha2),
                         by = "match") %>%
        dplyr::mutate(new_alpha1 = dplyr::if_else(!is.na(match), a1, alpha1),
                      new_alpha1 = dplyr::if_else(is.na(new_alpha1), alpha1, new_alpha1),
                      new_alpha2 = dplyr::if_else(!is.na(match), a2, alpha2),
                      new_alpha2 = dplyr::if_else(!is.na(match) & is.na(new_alpha2), alpha2, new_alpha2),
                      new_alpha2 = dplyr::if_else(is.na(alpha2), as.double(NA), new_alpha2)) %>%
        dplyr::select(alpha1 = new_alpha1, alpha2 = new_alpha2) %>%
        as.matrix()
      
    } else {
      alpha_init <- alpha_init %>%
        dplyr::select(-name) %>%
        as.matrix()
    }
    
    rownames(alpha_init) <- nam
    beta_init <- ifelse(!is.na(alpha_init), stats::rnorm(length(alpha_init), 0, 5), NA)
    colnames(beta_init) <- c("beta1", "beta2")
    rownames(beta_init) <- nam
    
    theta_init <- scale(rowMeans(data, na.rm = TRUE))[, 1]
    if (model == "multi_dyn") {
      if (is.null(T)) stop("Please supply `T`!")
      theta_init <- matrix(rep(theta_init, T), ncol = T)
    }
    if (!is.null(constraint)) {
      if (model == 'multi') {
        if (theta_init[constraint] < 0) theta_init <- -theta_init
      } else if (model == 'multi_dyn') {
        if (length(constraint) == 1) {
          constraint <- rep(constraint, ncol(theta_init))
        } else if (length(constraint) != ncol(theta_init)) {
          stop('`constraint` must be the same length as the number of sessions.')
        }
        theta_init <- set_constraint(theta_init, constraint)
      }
    }
    
    tmp <- list(alpha = alpha_init, beta = beta_init, theta = theta_init)
  }
  return(tmp)
}
