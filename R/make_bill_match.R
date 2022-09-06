#' @title Generating matched bill indicator
#' @description \code{make_bill_match} is an auxiliary function for pgIRT model and generates matched bill indicator for dynamic pgIRT models.
#'
#' @param data a matrix object, roll-call matrix.
#' @param bill_pairs a data.frame, tbl_df or matrix object (any x 2), the first column contains ID (integer or numeric) of bills 
#' which must be matched with the colunm names of \code{data} and the second column contains matched bills' ID before the bill. 
#' 
#' @importFrom dplyr %>% mutate select left_join arrange if_else row_number
#' @export

make_bill_match <- function(data, 
                            bill_pairs) {
  
  if (!class(data)[1] %in% c('matrix', 'data.frame', 'tbl_df')) stop('`data` only allows "matrix", "data.frame" or "tbl_df".')
  if (class(bill_pairs)[1] == 'matrix') {
    colnames(bill_pairs) <- c('rcid1', 'rcid2')
  } else {
    names(bill_pairs) <- c('rcid1', 'rcid2')
  }
  rcids <- dplyr::tibble(rcid = colnames(data) %>% 
                           as.numeric())
  matched_res <- rcids %>% 
    dplyr::left_join(bill_pairs %>% 
                       dplyr::select(rcid = rcid1, match = rcid2), by = 'rcid') %>% 
    dplyr::mutate(number = dplyr::row_number() - 1, .before = rcid) %>% 
    dplyr::arrange(rcid)
  
  matched_res <- matched_res %>% 
    dplyr::left_join(matched_res %>% 
                       dplyr::select(match_num = number, match = rcid), by = "match") %>% 
    dplyr::mutate(nomatch = if_else(!is.na(match) == is.na(match_num), 1, 0),
                  match = dplyr::if_else(nomatch == 1, as.double(NA), match)) %>% 
    dplyr::select(-nomatch)
  ma <- matched_res$match_num
  names(ma) <- matched_res$rcid 
  return(ma)
}
