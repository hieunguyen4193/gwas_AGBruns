# INPUT ARGS

stratum=$1;
chrom=$2;
region_start=$3;
region_end=$4;
path_to_04_output=$5;

echo -e "Working on region: chrom " $chrom " start = " $region_start " end = "$region_end " in condition "$stratum;

outdir="/home/hieunguyen/CRC1382/outdir";
PROJECT="AGBruns_SNP_data_full_corrected_20240312";

path_to_main_src=/home/hieunguyen/CRC1382/src_2023/${PROJECT};
plink=${path_to_main_src}/plink-1.07-x86_64/plink;

echo $path_to_04_output;

cd ${path_to_04_output};
cat region_chr${chrom}_${region_start}_${region_end}.raw.map | cut -d, -f2 > chrs.map
cat region_chr${chrom}_${region_start}_${region_end}.raw.map | cut -d, -f3 > variant_names.map
cat region_chr${chrom}_${region_start}_${region_end}.raw.map | cut -d, -f4 > morgans.map
cat region_chr${chrom}_${region_start}_${region_end}.raw.map | cut -d, -f5 > positions.map
paste chrs.map variant_names.map morgans.map positions.map | tail -n +2> region_chr${chrom}_${region_start}_${region_end}.map

path_to_ped_file=${path_to_04_output}/region_chr${chrom}_${region_start}_${region_end}

$plink --file $path_to_ped_file --r --noweb  --matrix
mv ${path_to_04_output}/plink.ld plink.r.matrix.ld

$plink --file $path_to_ped_file --r2 --noweb  --matrix
mv ${path_to_04_output}/plink.ld plink.r2.matrix.ld

$plink --file $path_to_ped_file --r --noweb --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0
mv ${path_to_04_output}/plink.ld plink.r.list.ld

$plink --file $path_to_ped_file --r2 --noweb --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0
mv ${path_to_04_output}/plink.ld plink.r2.list.ld