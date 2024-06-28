This blog describes the outputs of cellranger count piepline!  

After run command line successfully:  
~~~
cellranger count --id run_count_1kpbmcs \
   --fastqs /mnt/home/user.name/yard/run_cellranger_count/pbmc_1k_v3_fastqs \
   --sample pbmc_1k_v3 \
   --transcriptome /mnt/home/user.name/yard/run_cellranger_count/refdata-gex-GRCh38-2024-A \
   --create-bam true \
   --localcores 16 \
   --localmem 100
~~~
The terminal will show you:
~~~
Outputs:
- Run summary HTML:                         /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/web_summary.html
- Run summary CSV:                          /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/metrics_summary.csv
- BAM:                                      /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/possorted_genome_bam.bam
- BAM BAI index:                            /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/possorted_genome_bam.bam.bai
- BAM CSI index:                            null
- Filtered feature-barcode matrices MEX:    /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/filtered_feature_bc_matrix
- Filtered feature-barcode matrices HDF5:   /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/filtered_feature_bc_matrix.h5
- Unfiltered feature-barcode matrices MEX:  /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/raw_feature_bc_matrix
- Unfiltered feature-barcode matrices HDF5: /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/raw_feature_bc_matrix.h5
- Secondary analysis output CSV:            /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/analysis
- Per-molecule read information:            /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/molecule_info.h5
- CRISPR-specific analysis:                 null
- Antibody aggregate barcodes:              null
- Loupe Browser file:                       /md01/anzy/PDJ/test/run_count_1kpbmcs/outs/cloupe.cloupe
- Feature Reference:                        null
- Target Panel File:                        null
- Probe Set File:                           null

Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

2024-06-28 23:23:44 Shutting down.
Saving pipestance info to "run_count_1kpbmcs/run_count_1kpbmcs.mri.tgz"
~~~
