#' Align Homologous Exons
#'
#' Parse and assemble AS data between species
#' @param 
#' key A file containing homologous exons between species. 
#' data A list containing INCLUSION files from vast-tools ouput. The order of species in list should match the exon 'key'
#' @export
#' @examples
#' assemble(map = AltEx.key, data = list(INCLUSION_LEVELS_FULL.Hsa, INCLUSION_LEVELS_FULL.Ptr, INCLUSION_LEVELS_FULL.Mma)) -> PSI_QS_FULL

assemble <- function(map, data) {
  
  for (i in 1:length(data)) {
  
    data[[i]] -> incl_table
  
    # only take columns w/ exon names and PSI/QS data
    incl_table[,c(2, 7:ncol(incl_table))] -> incl_table
  
    # subset data to align equivalent exons between species 
    incl_table=incl_table[incl_table[,1] %in% map[,i],]
  
    # remove exons w/ more than one equivalent exon in other species (duplicates)
    map[!(duplicated(map[,i])), ] -> map
  
    # sort based on species column
    map[ order(map[,i]),] -> map
    incl_table[ order(incl_table[,1]),] -> incl_table
  
    # combine species and map together
    cbind(map, incl_table[,2:ncol(incl_table)]) -> map
  
  }
  return(map)
}