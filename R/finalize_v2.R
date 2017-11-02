#' Finalize PSI and varPSI files
#'
#' Finalizes PSI and varPSI files, to be handled by EVE. It changes PSI values to 0-1 scale and all NA's are replaced by -1. Also, any extra columns w/ exon names are removed.
#' @param 
#' data File created by filter_exons() to be formatted for EVE model.
#' sp Total number of columns w/ exon names.
#' PSI Is the file you are finalizing a PSI file? Defaults to TRUE. 
#' @export
#' @examples
#' finalize(data = PSI_brain, sp = 3, PSI = TRUE) -> PSI_brain
#' finalize(data = varPSI_brain, sp = 3, PSI = FALSE) -> varPSI_brain


finalize <- function (data, sp, PSI=TRUE) {
  
  if (PSI==TRUE) {
    # multiply all values by 0.01
    data[,(1+sp):ncol(data)]*0.01 -> data[,(1+sp):ncol(data)]
  
    # replace NA's w/ -1
    data[is.na(data)] <- -1
  
    # remove any exon names if there is more than 1
    if (sp !=0) {
    data[, c(1, (1+sp):ncol(data))] -> data
    }
  }
  else {
    # remove any exon names if there is more than 1
    if (sp !=0) {
      data[, c(1, (1+sp):ncol(data))] -> data
    }
    
    # replace NA's w/ -1
    data[is.na(data)] <- -1
    
  }
  return(data)
}