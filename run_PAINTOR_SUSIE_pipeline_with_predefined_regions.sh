# NOD2: Chromosome 16, NC_000016.10 (50693606..50733075)
# TLR2: Chromosome 4, NC_000004.12 (153684280..153710637)
# TRAF6: Chromosome 11, NC_000011.10 (36483769..36510272)
# CCL2 (MCP-1): Chromosome 17, NC_000017.11 (34255285..34257203)
# NR1H4 (FXR): Chromosome 12, NC_000012.12 (100473866..100564414)

outdir="/home/hieunguyen/CRC1382/outdir"
PROJECT="AGBruns_SNP_data_full_corrected_20240312"

maf=0.01;
stratum="SBP_ever";
model="basic_assoc"
path_to_04_output=${outdir}/${PROJECT}/04_output

#####
chrom=8;
region_start=8000000;
region_end=9900000;
bash 04_PIPELINE_PAINTOR_AND_SUSIE.sh $maf $stratum $chrom $region_start $region_end $model ${path_to_04_output}/${stratum}/MAF_${maf}/${model}/chr${chrom}_${region_start}_${region_end}

#####
chrom=16;
region_start=49693606;
region_end=51733075;
bash 04_PIPELINE_PAINTOR_AND_SUSIE.sh $maf $stratum $chrom $region_start $region_end $model ${path_to_04_output}/${stratum}/MAF_${maf}/${model}/chr${chrom}_${region_start}_${region_end}

#####
chrom=4;
region_start=152684280;
region_end=154710637;
bash 04_PIPELINE_PAINTOR_AND_SUSIE.sh $maf $stratum $chrom $region_start $region_end $model ${path_to_04_output}/${stratum}/MAF_${maf}/${model}/chr${chrom}_${region_start}_${region_end}

#####
chrom=11;
region_start=35483769;
region_end=37510272;
bash 04_PIPELINE_PAINTOR_AND_SUSIE.sh $maf $stratum $chrom $region_start $region_end $model ${path_to_04_output}/${stratum}/MAF_${maf}/${model}/chr${chrom}_${region_start}_${region_end}

#####
chrom=17;
region_start=33255285;
region_end=35257203;
bash 04_PIPELINE_PAINTOR_AND_SUSIE.sh $maf $stratum $chrom $region_start $region_end $model ${path_to_04_output}/${stratum}/MAF_${maf}/${model}/chr${chrom}_${region_start}_${region_end}


#####
chrom=12;
region_start=99473866;
region_end=101564414;
bash 04_PIPELINE_PAINTOR_AND_SUSIE.sh $maf $stratum $chrom $region_start $region_end $model ${path_to_04_output}/${stratum}/MAF_${maf}/${model}/chr${chrom}_${region_start}_${region_end}
