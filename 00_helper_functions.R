#####----------------------------------------------------------------------#####
##### Helper functions
#####----------------------------------------------------------------------#####
`%ni%` = Negate(`%in%`)

#####----------------------------------------------------------------------#####
# function to create an interactive data table
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}


generate_variant_genotype_info <- function(mutation, peddf){
  gt <- subset(snp.table.1, snp.table.1$Name == mutation)$SNP
  allele1 <- str_replace_all(str_split(gt, "/")[[1]][[1]], "\\[", "")
  allele2 <- str_replace_all(str_split(gt, "/")[[1]][[2]], "\\]", "")
  
  convert_genotype_to_SNPCall <- function(allele1, allele2, genotype){
    if (genotype == "AA"){
      return(c(allele1, allele1))
    } else if(genotype == "BB"){
      return(c(allele2, allele2))
    } else if (genotype == "AB"){
      return(c(allele1, allele2))
    } else if (genotype == "NC"){
      return(c(0, 0))
    }
  }
  
  col1 <- unlist(lapply(peddf[[mutation]], function(x){
    return(convert_genotype_to_SNPCall(allele1, allele2, x)[[1]])
  }))
  
  col2 <- unlist(lapply(peddf[[mutation]], function(x){
    return(convert_genotype_to_SNPCall(allele1, allele2, x)[[2]])
  }))
  return(list(allele.col1 = col1, allele.col2 = col2))
}



