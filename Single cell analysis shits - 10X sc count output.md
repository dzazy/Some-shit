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

2.output folder named by '--id' parameter(e.g. output folder named 'run_count_1kpbmcs')
The output folder structure:

      run_count_1kpbmcs/
         ├── _cmdline                    Cellranger count command
         ├── _filelist                   Detailed output files information,include size,permissions,owner,group,date,path...
         ├── _finalstate                 Detailed variable and values of the job
         ├── _invocation                 Detailed arguments of the job
         ├── _jobmode                    Run environment
         ├── _log                        Log file of the job
         ├── _mrosource                  Job information
         ├── outs                        Outputs folder
         ├── _perf                       Job information
         ├── _perf._truncated_           Job information
         ├── run_count_1kpbmcs.mri.tgz   Contains a number of logs (including error messages, status reports, and information about the system) that are useful for trouble-shooting failures
         ├── SC_RNA_COUNTER_CS           Folder of intermediate file
         ├── _sitecheck                  Detailed run environment parameters
         ├── _tags                       Blank
         ├── _timestamp                  Start and end timepoint
         ├── _uuid                       Uuid
         ├── _vdrkill                    Job information
         └── _versions                   Software and pipeline version

The result folder is 'run_count_1kpbmcs/outs'  
'outs' folder structure:  

      outs/
         ├── analysis
         │   ├── clustering
         │   │   ├── gene_expression_graphclust
         │   │   │   └── clusters.csv
         │   │   ├── gene_expression_kmeans_10_clusters
         │   │   │   └── clusters.csv
         |   |   |       ...
         │   │   └── gene_expression_kmeans_9_clusters
         │   │       └── clusters.csv
         │   ├── diffexp
         │   │   ├── gene_expression_graphclust
         │   │   │   └── differential_expression.csv
         │   │   ├── gene_expression_kmeans_10_clusters
         │   │   │   └── differential_expression.csv
         |   |   |       ...
         │   │   └── gene_expression_kmeans_9_clusters
         │   │       └── differential_expression.csv
         │   ├── pca
         │   │   └── gene_expression_10_components
         │   │       ├── components.csv
         │   │       ├── dispersion.csv
         │   │       ├── features_selected.csv
         │   │       ├── projection.csv
         │   │       └── variance.csv
         │   ├── tsne
         │   │   └── gene_expression_2_components
         │   │       └── projection.csv
         │   └── umap
         │       └── gene_expression_2_components
         │           └── projection.csv
         ├── cloupe.cloupe
         ├── filtered_feature_bc_matrix
         │   ├── barcodes.tsv.gz
         │   ├── features.tsv.gz
         │   └── matrix.mtx.gz
         ├── filtered_feature_bc_matrix.h5
         ├── metrics_summary.csv
         ├── molecule_info.h5
         ├── possorted_genome_bam.bam
         ├── possorted_genome_bam.bam.bai
         ├── raw_feature_bc_matrix
         │   ├── barcodes.tsv.gz
         │   ├── features.tsv.gz
         │   └── matrix.mtx.gz
         ├── raw_feature_bc_matrix.h5
         └── web_summary.html

      31 directories, 41 files



