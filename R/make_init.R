#' @title Generating initial values for pgIRT model
#' @description \code{make_init} is an auxilary function for pgIRT model and generates optimal initial values based on roll-call matrix.
#' 
#' @param data a matrix object, roll-call matrix 
#' @param model string, one of "bin" (binary), "bin_dyn" (dynamic binary), "multi" (multinomial) or "multi_dyn" (dynamic multinomial)
#' @param T string, indicating the length of times. If you choose "bin_dyn" or "multi_dyn" for `model`, you must supply this argument.
#' @param bill_match Numeric vector which indicates matched votes.
#' @param constraint a double, indicating the voter whose ideal point is always set positive
#' @export

make_init <- function(data, model = c("bin", "bin_dyn", "multi", "multi_dyn"), 
                      T = NULL, bill_match = NULL, constraint = NULL) {
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
  
  I <- nrow(data)
  J <- ncol(data)
  
  if (model %in% c("bin", "bin_dyn")) {
    ALLPOINT <- sum(is.na(data)) + sum(data == 1, na.rm = TRUE) + sum(data == 0, na.rm = TRUE)
    if (ALLPOINT != I * J) stop("For binomial model, the data is allowed to contain only NA, 1 and 0.")
    
    theta_init <- scale(rowMeans(data, na.rm = TRUE))[, 1]
    if (!is.null(constraint)) {
      if (theta_init[constraint] < 0) theta_init <- -theta_init
    }
    alpha_init <- beta_init <- rep()
    for (j in 1:ncol(data)) {
      c <- coef(glm(data[, j] ~ theta_init, family = "binomial")) %>%
        suppressWarnings()
      c <- ifelse(is.na(c) | abs(c) > 5, .1, c)
      alpha_init[j] <- c[1]
      beta_init[j] <- c[2]
    }
    if (model == "bin_dyn") {
      if (is.null(T)) stop("Please supply `T`!")
      theta_init <- matrix(rep(theta_init, T), ncol = T)
    }
    tmp <- list(alpha = alpha_init, beta = beta_init, theta = theta_init)
  }
  
  
  if (model %in% c("multi", "multi_dyn")) {
    ALLPOINT <- sum(is.na(data)) + sum(data == 1, na.rm = TRUE) + sum(data == 2, na.rm = TRUE) + sum(data == 3, na.rm = TRUE) 
    if (ALLPOINT != I * J) stop("For multinomial model, the data is allowed to contain only NA, 1, 2 and 3.")
    
    idd <- tibble(id_num = 1:ncol(data), name = colnames(data))
    temp_d <- data %>% 
      as_tibble() %>% 
      mutate(id = rownames(data)) %>% 
      pivot_longer(-id) %>% 
      filter(!is.na(value)) %>% 
      right_join(idd, by = "name")
    
    sb_calc <- temp_d %>% 
      group_by(id_num, value) %>% 
      summarise(count = n()) %>% 
      ungroup() %>% 
      group_by(id_num) %>% 
      summarise(vote = value, count = count, all = sum(count),
                ref_cat = max(value)) %>% 
      ungroup() %>% 
      group_by(id_num, vote) %>% 
      summarise(count = count, probs = count / all, all = all,
                ref_cat = ref_cat) %>% 
      group_by(id_num) %>% 
      summarise(vote = vote, count = count, 
                probs = probs, cumsum = cumsum(count) / all,
                sb = probs / (1 - (cumsum - probs)),
                ref_cat = ref_cat) %>% 
      filter(vote != max(vote)) %>% 
      ungroup() %>% 
      mutate(sb_value = qlogis(sb)) %>% 
      select(id_num, vote, sb_value, ref_cat)
    
    ## Assign same values to matched resolutions
    sb1 <- sb_calc %>% 
      filter(vote == 1) %>%
      select(alpha1 = sb_value, id_num)
    sb2 <- sb_calc %>% 
      filter(vote == 2) %>% 
      select(alpha2 = sb_value, id_num)
    alpha_init <- idd %>% 
      left_join(sb1, by = "id_num") %>% 
      left_join(sb2, by = "id_num") %>% 
      select(name, alpha1, alpha2)
    nam <- alpha_init$name
    if (!is.null(bill_match)) {
      alpha_init <- alpha_init %>% 
        bind_cols(tibble(match = as.character(bill_match))) %>% 
        left_join(alpha_init %>% select(match = name, 
                                       a1 = alpha1, a2 = alpha2),
                  by = "match") %>% 
        mutate(new_alpha1 = if_else(!is.na(match), a1, alpha1),
               new_alpha1 = if_else(is.na(new_alpha1), alpha1, new_alpha1),
               new_alpha2 = if_else(!is.na(match), a2, alpha2),
               new_alpha2 = if_else(!is.na(match) & is.na(new_alpha2), alpha2, new_alpha2),
               new_alpha2 = if_else(is.na(alpha2), as.double(NA), new_alpha2)) %>% 
        select(alpha1 = new_alpha1, alpha2 = new_alpha2) %>% 
        as.matrix()
      
    } else {
      alpha_init <- alpha_init %>% 
        select(-name) %>% 
        as.matrix()
      
    }
    
    rownames(alpha_init) <- nam
    beta_init <- ifelse(!is.na(alpha_init), rnorm(length(alpha_init), 0, 5), NA)
    colnames(beta_init) <- c("beta1", "beta2")
    rownames(beta_init) <- nam
    
    theta_init <- scale(rowMeans(data, na.rm = TRUE))[, 1]
    if (!is.null(constraint)) {
      if (theta_init[constraint] < 0) theta_init <- -theta_init
    }
    if (model == "multi_dyn") {
      if (is.null(T)) stop("Please supply `T`!")
      theta_init <- matrix(rep(theta_init, T), ncol = T)
    }
    
    tmp <- list(alpha = alpha_init, beta = beta_init, theta = theta_init)
  }
  return(tmp)
}
