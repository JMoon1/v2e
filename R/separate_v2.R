#' Separate PSI and QS data
#'
#' Retrieves PSI data and calculates var(PSI) between all tissues.
#' @param
#' data The file created by 'align' function, containing all species' PSI and Quality Score information
#' sp The total number of columns containing exon names
#' tiss A character vector containing tissue names (not case-sensitive, but spelling counts)
#' @export
#' @examples
#' separate(data = PSI_QS_FULL, sp = 3, tiss = c("brain", "cerebellum", "heart", "kidney", "liver", "testis"))

separate <- function(data, sp, tiss) {
  
  # 1) separate data to get PSI's
  data[, c(1:sp, seq(sp+1, ncol(data), by = 2))] ->> PSI_FULL

  # 1) separate data to get QS
  data[, c(1:sp, seq(sp+2, ncol(data), by = 2))] ->> varPSI_FULL

  # calculate varPSI
  for(i in 1:nrow(varPSI_FULL)) {
    
    # 2) parse out QS to get raw reads
    sub(".*@", "", varPSI_FULL[i,(sp+1):ncol(varPSI_FULL)]) -> varPSI_FULL[i,(sp+1):ncol(varPSI_FULL)]
    
  }
  
  for(i in 1:nrow(varPSI_FULL)) {
    
    # 3) separate inclusion and exclusion reads
    sub(",.*", "", varPSI_FULL[i,(sp+1):ncol(varPSI_FULL)]) -> reads_i
    as.numeric(reads_i) -> reads_i
    
    sub(".*,", "", varPSI_FULL[i,(sp+1):ncol(varPSI_FULL)]) -> reads_e
    as.numeric(reads_e) -> reads_e
    
    # make total reads and PSI vector
    total= reads_i + reads_e
    PSI=(reads_i/total)
    
    # initialize PSI
    for (j in 1:length(total)) {
      
      # if total number of reads equals zero
      if (total[j] == 0 | is.na(total[j])) {
        
        varPSI_FULL[i,(sp+j)]  <- NA
      }
      
      # when PSI equals 100%
      if ((reads_i[j] != 0 && (!is.na(reads_i[j]))) && (reads_e[j] == 0 | is.na(reads_e[j]))) {
        
        varPSI_FULL[i,(sp+j)] <- 1/(2*(total[j]+1))
        
      }
      
      # when PSI does not equal zero or 100%
      if ((reads_i[j] != 0 && (!is.na(reads_i[j]))) && (reads_e[j] != 0 && (!is.na(reads_e[j])))) {
        
        
        varPSI_FULL[i,(sp+j)] <- PSI[j]*(1-PSI[j])/total[j]
        
      }
      
      #when inclusion reads are zero, but exclusion reads are not
      if ((reads_i[j] == 0 | is.na(reads_i[j])) && (reads_e[j] != 0 && (!is.na(reads_e[j])))) {
        
        varPSI_FULL[i,(sp+j)] <- 1/(2*(total[j]+1))
        
      }
      
    }
  }
  

  for (i in 1:length(tiss)) {
    
    # separate PSI_FULL by tissues
    paste("PSI", tiss[i], sep = "_") -> nam_PSI
    
    assign(nam_PSI, PSI_FULL[, c(1:sp, grep(pattern=tiss[i], names(PSI_FULL), ignore.case = TRUE))], envir = .GlobalEnv)

    # separate varPSI_FULL by tissues
    paste("varPSI", tiss[i], sep = "_") -> nam_varPSI
    
    assign(nam_varPSI, varPSI_FULL[, c(1:sp, grep(pattern=tiss[i], names(varPSI_FULL), ignore.case = TRUE))], envir = .GlobalEnv)
  
    }
}