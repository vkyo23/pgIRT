#' @title Cleaning roll-call matrix for pgIRT model
#' @description \code{clean_rolllcall} is an auxiliary function for pgIRT model and cleans a roll-call matrix
#' by removing noisy items and individuals who have no response records.
#' 
#' @param data a matrix, data.frame or tbl_df of roll-call (i.e. individual-item) matrix (I x J).
#' @export
#' 
#' @examples
#' \dontrun{
#' data(m_data_dyn)
#' mat <- make_rollcall(m_data_dyn,
#'                      unit = "unit",
#'                      bill = "bill",
#'                      vote_col = "vote")
#' mat <- clean_rollcall(mat)
#' }

clean_rollcall <- function(data) {
  if (is.null(rownames(data))) rownames(data) <- 1:nrow(data)
  if (is.null(colnames(data))) colnames(data) <- 1:ncol(data)
  category <- apply(data, 2, function(x) unique(x[!is.na(x)]))
  num_cat <- unlist(lapply(category, length))
  drop_res <- which(num_cat == 1)
  
  no_attend <- which(apply(data, 1, function(x) all(is.na(x))))
  
  if (length(drop_res) != 0) {
    cat("Remove some bills because they are unanimous votings:", colnames(data)[drop_res], "\n")
    data <- data[, -drop_res]
  }
  dr <- which(!apply(data, 2, function(x) min(x, na.rm = TRUE) == 1))
  if (length(dr) != 0) {
    cat("Remove some bills which do not start from 1:", colnames(data)[dr], "\n")
    data <- data[, -dr]
  }
  if (length(no_attend) != 0) {
    cat("Remove some units who have no voting record:", rownames(data)[no_attend], "\n")
    data <- data[-no_attend, ]
  }
  
  return(data)
}
