1. Run `00_preprocessing_full_tables.R`

2. Run `01_split_batch_all_cases.R`

3. Run `02_run_preprocess_input_chunks.sh`

4. Run `02_preprocess_output_chunks_to_PLINK_inputs.sh`

5. Run `02_generate_Rdata_for_SNPStats.R`

6. Run `00_modify_snp_table_to_PLINK_input.sh`

7. Run PLINK

8. Run `run_to_generate_GWAS_reports.R`

9. Run `run_PAINTOR_SUSIE_pipeline_with_predefined_regions.sh` and `run_PAINTOR_SUSIE_pipeline_with_candidate_regions.R` to generate all post-GWAS results.

10. Run `get_data_from_input_SNP.R` to get the output table. 

