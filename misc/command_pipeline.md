# Introduction

Note on all commands used in conducting GWAS for the dataset from AG Bruns. 

# Main pipeline

## Preprocessing: From "full table" to "reduced table"

See `00_preprocessing_full_tables.R`. From the "full table" dataframe, keep only the column `GType` for GWAS. 

```R
gc()
rm(list = ls())

library(stringr)
library(tidyverse)
library(dplyr)
library(comprehenr)

path.to.maindir <- "/home/hieunguyen/CRC1382/storage/AGBruns_SNP_data"
all.full.tables <- Sys.glob(file.path(path.to.maindir, "*", "*Full_Data_Table.txt"))

for (file in all.full.tables){
  full.table <- read.csv(file, sep = "\t")
  cols.to.keep <- c("Name", to_vec(for (item in colnames(full.table)) if(grepl("GType", item) == TRUE) item))
  reduced.full.table <- full.table[, cols.to.keep] %>%#
    column_to_rownames("Name")
  write.csv(reduced.full.table, str_replace(file, ".txt", ".reduced.txt"))
}
```

## Split the "reduced table" to smaller chunks for faster parallel pre-processing
Since the current data table is in a `GenomeStudio` format, we need to convert it to a format that is compatible with `snpStats` and `PLINK`. 

Firstly, we split the full table into small chunks by the following function; see `01_split_batch_SBP_vs_nonSBP.R`. Depending on the condition we want to stratify the data, we change the input
variable: `stratum`

```r
gc()
rm(list = ls())

stratum <- "SBP"

library(vroom)
library(comprehenr)
library(readxl)
library(writexl)
library(stringr)
library(tidyverse)
library(dplyr)

#####----------------------------------------------------------------------#####
##### PATHS AND CONFIGURATIONS
#####----------------------------------------------------------------------#####
path.to.main.dir <- "/home/hieunguyen/CRC1382"

path.to.outdir <- file.path(path.to.main.dir, "outdir")

path.to.storage <- "/home/hieunguyen/CRC1382/storage"
main.data.dir <- file.path(path.to.storage, "AGBruns_SNP_data")

path.to.save.output <- file.path(path.to.outdir, "AGBruns_SNP_data_06032023", "OUTPUT")
path.to.01.output <- file.path(path.to.save.output, "01_output")

dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.parallel.workdir <- file.path(path.to.save.output, "parallel_processing")

#####----------------------------------------------------------------------#####
##### PREPROCESSING INPUT DATA TABLES
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "snp_table.map")) == FALSE){
  print("Generating snp_table.map")
  ##### reduce the SNP table
  all.snp.tables <- Sys.glob(file.path(main.data.dir, "*", "*SNP_Table.txt"))
  
  snp.table.1 <- read.csv(all.snp.tables[[1]], sep = "\t")
  reduced.snp.table1 <- subset(snp.table.1, select = c("Index", "Name", "Chr", "Position",
                                                       "ChiTest100", "Het.Excess", "AA.Freq", "AB.Freq", "BB.Freq", 
                                                       "Minor.Freq", "Call.Freq",  "SNP"))
  
  snp.table.1 <- subset(snp.table.1, grepl("D", snp.table.1$SNP) == FALSE)
  snp.table.1 <- subset(snp.table.1, grepl("I", snp.table.1$SNP) == FALSE)
  
  snp.table.2 <- read.csv(all.snp.tables[[2]], sep = "\t")
  reduced.snp.table2 <- subset(snp.table.2, select = c("Index", "Name", "Chr", "Position",
                                                       "ChiTest100", "Het.Excess", "AA.Freq", "AB.Freq", "BB.Freq", 
                                                       "Minor.Freq", "Call.Freq",  "SNP"))

  snp.table.2 <- subset(snp.table.2, grepl("D", snp.table.2$SNP) == FALSE)
  snp.table.2 <- subset(snp.table.2, grepl("I", snp.table.2$SNP) == FALSE)
  
  ##### Generate a .map file
  mapdf <- subset(snp.table.1, select = c(Chr, Name))
  mapdf$Pos.morgans <- 0
  mapdf$coord <- snp.table.1$Position
  
  write.csv(mapdf, file = file.path(path.to.01.output, "snp_table.map"))
  write.table(mapdf, file = file.path(path.to.01.output, "snp_table_for_GWAS.map"), sep = "\t", col.names = FALSE, row.names = FALSE)
  write.csv(snp.table.1, file = file.path(path.to.01.output, "snp_table_full_info_plate1.map"))
  write.csv(snp.table.2, file = file.path(path.to.01.output, "snp_table_full_info_plate2.map"))

} else {
  print("SNP_table.map existed!")
  mapdf <- read.csv(file.path(path.to.01.output, "snp_table.map"))
  if ("X" %in% colnames(mapdf)){
    mapdf <- subset(mapdf, select = -c(X))
  }
}

##### preprocessing sample tables
all.sample.tables <- Sys.glob(file.path(main.data.dir, "*", "*Samples_Table.txt"))
  
sample.table <- data.frame()
for (file in all.sample.tables){
  tmp <- read.csv(file, sep = "\t")
  sample.table <- rbind(sample.table, tmp)
}
  
sample.table <- subset(sample.table, select = c("Unique.sample.A.number", "Sample.ID", "Gender.Est"))
  
colnames(sample.table) <- c("unique_ID", "Sample.ID", "Gender")
  
meta.data <- read_excel(file.path(path.to.storage, "AGBruns_SNP_data", "Inkludierte_14022023.xlsx"))

if (stratum == "SBP"){
  meta.data <- subset(meta.data, select = c(`unique sample A-number`, Hieu_SBP))
} else if (stratum == "CDAD"){
  meta.data <- subset(meta.data, select = c(`unique sample A-number`, Hieu_CDAD))
} else if (stratum == "SAS"){
  meta.data <- subset(meta.data, select = c(`unique sample A-number`, Hieu_SAS))
} else if (stratum == "death") {
  meta.data <- subset(meta.data, select = c(`unique sample A-number`, Hieu_death))
} else if (stratum == "transplant"){
  meta.data <- subset(meta.data, select = c(`unique sample A-number`, Hieu_transplant))
}

colnames(meta.data) <- c("unique_ID", stratum)

meta.data <- subset(meta.data, meta.data$unique_ID %in% unique(sample.table$unique_ID))
merge.sample.table <- merge(sample.table, meta.data, by.x = "unique_ID", by.y = "unique_ID", all.x = TRUE)
  
merge.sample.table <- merge.sample.table[!duplicated(merge.sample.table$unique_ID),]
  
#  peddf <- subset(merge.sample.table, select = c(Sample.ID, Gender, SBP))
peddf <- merge.sample.table[,c("Sample.ID", "Gender", stratum)]

peddf <- replace(peddf, is.na(peddf), 0)
peddf$family_id <- to_vec(for (item in seq(1, nrow(merge.sample.table))) sprintf("fam_%s", item))
peddf$father_id <- 0
peddf$mother_id <- 0
# peddf <- subset(peddf, select = c(family_id, Sample.ID, father_id, mother_id, Gender, SBP))
peddf <- peddf[, c("family_id", "Sample.ID", "father_id", "mother_id", "Gender", stratum)]  

all.full.tables <- Sys.glob(file.path(main.data.dir, "*", "*Full_Data_Table.reduced.txt"))
  
full.table.df <- data.frame()
for (file in all.full.tables){
  tmp.full.table <- read.csv(file, sep = ",") %>% column_to_rownames("X")
  tmp.full.table <- t(tmp.full.table)
  row.names(tmp.full.table) <- to_vec(for (item in row.names(tmp.full.table)) 
    str_replace_all(item, ".GType", ""))
  full.table.df <- rbind(full.table.df, tmp.full.table)
}
  
full.table.df <- full.table.df[, mapdf$Name]

row.names(full.table.df) <- to_vec(for (item in row.names(full.table.df)) str_replace(item, "X", ""))
  
peddf <- merge(peddf, full.table.df %>% rownames_to_column("Sample.ID"), by.x = "Sample.ID", by.y = "Sample.ID")
  
peddf <- peddf[, c(c("family_id", "Sample.ID", "father_id", "mother_id", "Gender", stratum), colnames(peddf)[7:length(colnames(peddf))])]
  
#####----------------------------------------------------------------------#####
##### Split the current peddf file to smaller chunks for faster processing
#####----------------------------------------------------------------------#####
all.mutations <- colnames(peddf)[7:ncol(peddf)]

split.all.mutations <- split(all.mutations, ceiling(seq_along(all.mutations) / 10000))

path.to.input.chunks <- file.path(path.to.parallel.workdir, "input_chunks")
dir.create(path.to.input.chunks, showWarnings = FALSE, recursive = TRUE)

path.to.output.chunks <- file.path(path.to.parallel.workdir, "output_chunks")
dir.create(path.to.output.chunks, showWarnings = FALSE, recursive = TRUE)

for (group in names(split.all.mutations)){
  print(sprintf("Working on %s", group))
  selected.cols <- c(colnames(peddf)[1:6], split.all.mutations[[group]])
  tmp.peddf <- peddf[, selected.cols]
  
  write.csv(tmp.peddf, file.path(path.to.input.chunks, sprintf("chunk_%s.csv", group)))
}
```

## Convert the Genotype to PLINK/snpStats compatible format

### PLINK
Use the function `01_generate_ped_and_fam_files_for_PLINK.R`

```sh
for case in SBP CDAD death SAS transplant;do echo -e "working on the case " $case "\n" && parallel -j 20 Rscript 01_generate_ped_and_fam_files_for_PLINK.R --input_chunk {} --output /home/hieunguyen/CRC1382/outdir/AGBruns_SNP_data/OUTPUT/01_output/parallel_processing_${case}/output_chunks --stratum $case ::: /home/hieunguyen/CRC1382/outdir/AGBruns_SNP_data/OUTPUT/01_output/parallel_processing_${case}/input_chunks/*.csv;done

```

### snpStats
Use the function `01_generate_ped_and_fam_files.R`



## Final preprocessing: Cat and paste all chunks into one single file

Run the script `preprocess_output_chunks_to_PLINK_inputs.sh` 

```sh
path_to_plink_inputs="/home/hieunguyen/CRC1382/outdir/AGBruns_SNP_data/plink_inputs";

for group in SBP SAS CDAD death transplant;do \
    path_to_output_chunks=/home/hieunguyen/CRC1382/outdir/AGBruns_SNP_data/OUTPUT/01_output/parallel_processing_${group}/output_chunks;

    cd ${path_to_output_chunks};

    for file in chunk_{2..72}.csv;do cat $file | tail -n +2 | cut --complement -f1,2,3,4,5,6 > mutation_only_${file};done

    cat chunk_1.csv | tail -n +2 > mutation_only_chunk_1.csv 

    paste mutation_only_chunk_{1..72}.csv > all_mutations.ped

    mv all_mutations.ped ${group}.ped

    rsync -avh --progress ${group}.ped ${path_to_plink_inputs}

    rm -rf mutation_only_chunk*;done
```

- Could you please provide the minor allele frequencies of the candidate SNPs for the cases (SBP) and the controls (non-SBP) and give the univariate odds ratios associated.

- Could you check, whether there is a gene-dose-relation (SBP frequency in wild-type, heterozygous and homozygous mutants) of those candidate SNPs.

- Could you check the for the other comparisons to see whether these are less genetically determined.  

## Run snpStats 

## Run PLINK
