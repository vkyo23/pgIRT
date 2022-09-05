#' @title Generating roll-call matrix for pgIRT model
#' @description \code{make_rolllcall} is an auxilary function for pgIRT model
#' and generates roll-call matrix from long-dataframe.
#' 
#' @param @param dataframe data.frame, tbl_df or matrix object; dataframe (long-format, i.e. voter-bill unit) of voting data. 
#' @param unit_id string, column name which indicates voters.
#' @param bill_id string, column name which indicates bills
#' @param vote_col string, column name which indicates the votes.' 
#' @param drop_unanimous bool, whether removing unanimous bills or not
#' @export

make_rollcall <- function(dataframe, unit_id = NULL, bill_id = NULL, vote_col = NULL, drop_unanimous = FALSE) {
  if (!class(dataframe)[1] %in% c("data.frame", "tbl_df")) {
    stop("`dataframe` must be 'matrix', 'data.frame' or 'tbl_df' object.")
  } else {
    dataframe <- as_tibble(dataframe)
  }
  if (is.null(unit_id)) stop("Please specify `unit_id`.")
  if (is.null(bill_id)) stop("Please specify `bill_id`.")
  if (is.null(vote_col)) stop("Please specify `vote_col`.")
  unit_id <- sym(unit_id)
  bill_id <- enquo(bill_id)
  vote_col <- enquo(vote_col)
  
  temp <- dataframe %>% 
    select(!!unit_id, !!bill_id, !!vote_col) %>% 
    pivot_wider(names_from = !!bill_id, values_from = !!vote_col) %>% 
    arrange(!!unit_id) %>% 
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
    max_cat <- lapply(category, max) %>% unlist()
    drop_res <- names(max_cat)[which(num_cat == 1)]
    
    ## Remove such resolutions
    if (length(drop_res) != 0) {
      cat("Remove some bills because they are unanimous votings:", drop_res, "\n")
      temp <- temp[, which(!colnames(temp) %in% drop_res)]
      num_cat <- num_cat[which(!names(num_cat) %in% drop_res)]
      max_cat <- max_cat[which(!names(max_cat) %in% drop_res)]
    }
    ## Remove countries which have never voted
    no_attend <- which(apply(temp, 1, function(x) all(is.na(x))))
    if (length(no_attend) != 0) {
      cat("Remove some units who have no voting record:", no_attend, "\n")
      temp <- temp[-no_attend, ]
      rname <- rname[-no_attend]
      rownames(temp) <- rname
    }
  }
  temp <- temp[, (apply(temp, 2, function(x) min(x, na.rm = TRUE) == 1) | apply(temp, 2, function(x) min(x, na.rm = TRUE) == 0))]
  return(temp)
}


