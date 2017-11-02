#' Filter Exons with Zero Variance
#'
#' Removes exons without evidence of alternative splicing
#' @param 
#' data The PSI file to be filtered
#' sp The total number of columns containing exon names
#' indivs A numerical vector containing the number of individuals for each tissue, in order from left to right
#' stringent Applies a stronger filter to return exons w/ at least 2 individuals per tissue/species (Required for EVE model). Defaults to FALSE.
#' @export
#' @examples
#' filter_exons(data = PSI_brain, sp = 3, indivs = c(6,6,3), stringent = TRUE) -> PSI_brain
#' 
#' # don't forget to subset the complementary varPSI file to match
#' varPSI_brain = varPSI_brain[varPSI_brain$V1 %in% PSI_brain$V1, ]


filter_exons <- function(data, sp, indivs, stringent = FALSE) {

  if (stringent==FALSE) {
    
    final_data <- data.frame()
    
    ## Drop exons w/ no variance (<2 unique values)
    for (i in 1:nrow(data)) { 
      
      # initialize empty vector that will contain data w/out NA's
      noNAvec <- c()
    
      # drop NA's
      for (j in (1+sp):ncol(data)) {
      
        if(!is.na(data[i,j])) {
        
          # assign non-NA values to noNAvec- noNAvec becomes list!
          noNAvec = c(noNAvec, data[i,j])
        
        }
      }
      
      # Count unique elements of noNAvec, drop row if < 2 (zero variance)
      if (length(unique(noNAvec)) >= 2) {
        #new_df = new_df[-c(i),]
        final_data = rbind(final_data, data[i,])
        
      }
    }
  }
  
  # 'stringent' requires at least 2 indivs per tissue/species
  else {
    
    final_data <- data.frame()
    
    ## Drop exons w/ no variance (<2 unique values)
    for (i in 1:nrow(data)) { 
      
      # initialize empty vector that will contain data w/out NA's
      noNAvec <- c()
      
      # drop NA's
      for (j in (1+sp):ncol(data)) {
        
        if(!is.na(data[i,j])) {
          
          # assign non-NA values to noNAvec- noNAvec becomes list!
          noNAvec = c(noNAvec, data[i,j])
          
        }
      }
      
      # Count unique elements of noNAvec, drop row if < 2 (zero variance)
      if (length(unique(noNAvec)) >= 2) {
        #new_df = new_df[-c(i),]
        final_data = rbind(final_data, data[i,])
        
      }
    }
    
    # initiate num to loop through each species
    num=sp+indivs[1]
    
    # at least 2 individuals per exon/species
    for (i in 1:length(indivs)) {
      
      if (i == 1) {
        final_data[rowSums(is.na(final_data[(1+sp):num])) <= indivs[i]-2, ] -> final_data
      }
      
      else {
        final_data[rowSums(is.na(final_data[(num+1):(num+indivs[i])])) <= indivs[i]-2,] -> final_data
        num=num+indivs[i]
      }
      
    }
  }
  return(final_data)
}