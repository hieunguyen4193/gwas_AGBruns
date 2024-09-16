##### INPUT ARGS

maf=$1;
stratum=$2;
chrom=$3;
region_start=$4;
region_end=$5;
model=$6;
output_dir=$7;

##### DIR
PAINTORdir="/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data_full_corrected_20240312/PAINTOR/PAINTOR_V3.0";
outdir=""/home/hieunguyen/CRC1382/outdir"";
PROJECT="AGBruns_SNP_data_full_corrected_20240312";

##### PROCESS INPUT FOR SUSIE
echo -e "#####--------------------------------------------------------##### \n"
echo -e "Generating input for SUSIE"
echo -e "#####--------------------------------------------------------##### \n"
Rscript 04_generate_ped_file_at_input_region_for_SUSIE.R --stratum $stratum --chrom $chrom --start $region_start --end $region_end --output_dir $output_dir;

##### CALCULATE LD MATRIX
echo -e "#####--------------------------------------------------------##### \n"
echo -e "CALCULATING LD MATRIX"
echo -e "#####--------------------------------------------------------##### \n"
bash 04_run_PLINK_LD_calculation.sh $stratum $chrom $region_start $region_end $output_dir;

##### PROCESS INPUT FOR PAINTOR
echo -e "#####--------------------------------------------------------##### \n"
echo -e "Generating input for PAINTOR"
echo -e "#####--------------------------------------------------------##### \n"
Rscript 04_PAINTOR_preprocess.R --stratum $stratum --chrom $chrom --start $region_start --end $region_end --maf $maf --output_dir $output_dir;

# sample_datadir=${outdir}/${PROJECT}/04_output/${stratum}/region_chr${chrom}_${region_start}_${region_end}/PAINTOR
sample_datadir=${output_dir}/PAINTOR;
${PAINTORdir}/PAINTOR -input ${sample_datadir}/input.file -Zhead ZSCORE.P1 -LDname LD -in $sample_datadir -out $sample_datadir -mcmc -annotations fake_annotations

##### RENDER FINAL HTML REPORTS
echo -e "#####--------------------------------------------------------##### \n"
echo -e "RENDERING REPORTS ..."
echo -e "#####--------------------------------------------------------##### \n"
echo -e $output_dir
Rscript 04_render_report.R --stratum $stratum --chrom $chrom --start $region_start --end $region_end --maf $maf --model $model --output_dir $output_dir;
