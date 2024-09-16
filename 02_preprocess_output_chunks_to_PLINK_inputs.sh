outdir="/home/hieunguyen/CRC1382/outdir";
PROJECT="AGBruns_SNP_data_full_corrected_20240312";
path_to_main_output=${outdir}/${PROJECT};

path_to_01_output=${path_to_main_output}/01_output;
path_to_02_output=${path_to_main_output}/02_output;
mkdir -p ${path_to_02_output};

all_GWAS_cases=$(cat all_GWAS_cases.txt);

for group in ${all_GWAS_cases};do \
    path_to_plink_inputs=${path_to_02_output}/${group};
    mkdir -p ${path_to_plink_inputs};

    echo -e "Working on the dataset " ${group} "\n";
    path_to_output_chunks=${path_to_01_output}/parallel_processing_${group}/output_chunks;

    cd ${path_to_output_chunks};

    for file in chunk_{2..68}.csv;do cat $file | tail -n +2 | cut --complement -f1,2,3,4,5,6 > mutation_only_${file};done

    cat chunk_1.csv | tail -n +2 > mutation_only_chunk_1.csv 

    paste mutation_only_chunk_{1..68}.csv > all_mutations.ped

    mv all_mutations.ped ${group}.ped

    rsync -avh --progress ${group}.ped ${path_to_plink_inputs}

    rm -rf mutation_only_chunk*;done