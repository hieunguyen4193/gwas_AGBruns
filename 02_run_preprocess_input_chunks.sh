# outdir="/media/hieunguyen/CRC1382H/outdir";
outdir="/home/hieunguyen/CRC1382/outdir";
PROJECT="AGBruns_SNP_data_full_corrected_20240312";
all_GWAS_cases=$(cat all_GWAS_cases.txt);

for case in ${all_GWAS_cases};do \
echo -e "working on the case " $case "\n" && parallel -j 20 Rscript 01_generate_ped_and_fam_files_for_PLINK.R --input_chunk {} \
--output ${outdir}/${PROJECT}/01_output/parallel_processing_${case}/output_chunks \
--stratum $case ::: ${outdir}/${PROJECT}/01_output/parallel_processing_${case}/input_chunks/*.csv;done