path_to_main_src="/home/hieunguyen/CRC1382/src_2023/AGBruns_SNP_data"
plink=${path_to_main_src}/plink-1.07-x86_64/plink

# check if the plink command line works! ${plink} --help

# Mofidy snp_table.map file from the current file by the following commands:
# cat snp_table.map | cut -d, -f2 > chrs.map
# cat snp_table.map | cut -d, -f3 > variant_names.map
# cat snp_table.map | cut -d, -f4 > morgans.map
# cat snp_table.map | cut -d, -f5 > positions.map
# paste chrs.map variant_names.map morgans.map positions.map | tail -n +2> snp_table.plink.map

# replace the SBP column by affected = 2, unaffected = 1 
# awk -F "\t" 'BEGIN{OFS="\t"}{if ($6==0) {$6=1} else {$6=2}}2' mytest.ped > mytest.modified.ped
# MAIN CMDS:
# awk -F "\t" 'BEGIN{OFS="\t"}{if ($6==0) {$6=1} else {$6=2}}2' '/home/hieunguyen/CRC1382/outdir/AGBruns_SNP_data/OUTPUT/parallel_processing/output_chunks/ready2input_snpStats.noheader.ped' > '/home/hieunguyen/CRC1382/outdir/AGBruns_SNP_data/OUTPUT/parallel_processing/plink_inputs/ready2input_PLINK.noheader.ped'

outdir="/home/hieunguyen/CRC1382/outdir";
PROJECT="AGBruns_SNP_data_full_corrected_20240312";

path_to_plink_inputs=${outdir}/${PROJECT};
groups=$(ls ${path_to_plink_inputs}/02_output);

for group in ${groups};do \
    echo -e "working on group " $group "\n";

    path_to_ped_file=${path_to_plink_inputs}/02_output/${group}/${group}.ped;
    path_to_map_file=${path_to_plink_inputs}/01_output/snp_table.plink.map;

    $plink --ped $path_to_ped_file --map $path_to_map_file --maf 0.01 --fast-epistasis --noweb --allow-no-sex --hwe 0.05

    mkdir -p ${outdir}/${PROJECT}/PLINK_results/maf_0.01_model_epistasis/${group}_with_epistasis_results_maf_0.01;

    mv plink.* ${outdir}/${PROJECT}/PLINK_results/maf_0.01_model_epistasis/${group}_with_epistasis_results_maf_0.01;
done
