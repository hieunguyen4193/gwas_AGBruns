cd /home/hieunguyen/CRC1382/outdir/AGBruns_SNP_data_full_corrected_20240312/01_output;
cat snp_table.map | cut -d, -f2 > chrs.map
cat snp_table.map | cut -d, -f3 > variant_names.map
cat snp_table.map | cut -d, -f4 > morgans.map
cat snp_table.map | cut -d, -f5 > positions.map
paste chrs.map variant_names.map morgans.map positions.map | tail -n +2> snp_table.plink.map
