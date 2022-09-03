#' @title Generating options for dynamic pgIRT model
#' @description \code{make_dyn_options} is an auxilary function for dynamic pgIRT model and generates options for such models.
#' 
#' @param dataframe data.frame, tbl_df or matrix object; dataframe (long-format, i.e. voter-bill unit) of voting data. 
#' @param unit_id string, column name which indicates voters.
#' @param bill_id string, column name which indicates bills
#' @param time_id string, column name which indicates times.
#' @param vote_col string, column name which indicates the votes.
#' @param add_bill_match Optional. Numeric vector which indicates matched votes.
#' @param drop_unanimous bool, whether removing unanimous bills or not
#' @export

make_dyn_options <- function(dataframe, unit_id = NULL, bill_id = NULL, time_id = NULL, vote_col = NULL,
                             add_bill_match = NULL, drop_unanimous = FALSE) {
  if (!class(dataframe)[1] %in% c("data.frame", "tbl_df")) {
    stop("`dataframe` must be 'matrix', 'data.frame' or 'tbl_df' object.")
  } else {
    dataframe <- as_tibble(dataframe)
  }
  if (is.null(unit_id)) stop("Please specify `unit_id`.")
  if (is.null(bill_id)) stop("Please specify `bill_id`.")
  if (is.null(time_id)) stop("Please specify `time_id`.")
  if (is.null(vote_col)) stop("Please specify `vote_col`.")
  unit_id <- rlang::sym(unit_id)
  bill_id <- rlang::sym(bill_id)
  time_id <- rlang::sym(time_id)
  vote_col <- rlang::enquo(vote_col)
  
  `%!in%` <- function(x, y) !x %in% y
  
  temp <- dataframe %>% 
    select(!!unit_id, !!bill_id, !!vote_col) %>% 
    spread(key = !!bill_id, value = !!vote_col) %>% 
    as.matrix()
  
  ## Extract country name
  rname <- temp[, 1] 
  
  ## Convert into numeric matrix
  temp <- apply(temp[, -1], 2, as.numeric)
  
  ## Label rows
  rownames(temp) <- rname
  
  if (drop_unanimous) {
    ## Remove resolutions which are unanimous voting.
    ## First, get number of voting category by resolution and
    ## record unanimous resolutions
    category <- apply(temp, 2, function(x) unique(x[!is.na(x)]))
    num_cat <- lapply(category, length) %>% unlist()
    drop_res <- names(num_cat)[which(num_cat == 1)]
    
    ## Remove countries which have never voted
    no_attend <- which(apply(temp, 1, function(x) all(is.na(x))))
    
    if (length(drop_res) != 0) {
      cat("Remove some bills because they are unanimous votings:", drop_res, "\n")
      dataframe <- dataframe %>% 
        filter(!!bill_id %!in% drop_res)
    }
    if (length(no_attend) != 0) {
      cat("Remove some units who have no voting record:", no_attend, "\n")
      dataframe <- dataframe %>% 
        filter(!!unit_id %!in% no_attend)
    }
  }
  
  ## Time-individual map (matrix) -----
  tmp1 <- dataframe %>% 
    arrange(!!unit_id) %>% 
    select(!!unit_id, !!bill_id, !!time_id) %>% 
    mutate(id := paste0(!!unit_id, "-", !!time_id),
           time := !!time_id - min(!!time_id)) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    arrange(!!unit_id) %>% 
    select(!!time_id, !!unit_id) %>% 
    as.matrix()
  tmp1[, 1] <- tmp1[, 1] - min(tmp1[, 1])
  
  tmp2 <- dataframe %>% 
    mutate(time := !!time_id - min(!!time_id)) %>% 
    distinct(!!bill_id, .keep_all = TRUE) %>% 
    select(!!bill_id, !!time_id) %>% 
    as.matrix()
  tmp2 <- as.numeric(tmp2[, 2]) - as.numeric(min(tmp2[, 2]))
  
  L <- list(session_individual = tmp1,
            bill_session = tmp2)
  
  if (!is.null(add_bill_match)) {
    L$bill_match <- add_bill_match
  } else {
    L$bill_match <- NULL
  }
  return(L)
}


