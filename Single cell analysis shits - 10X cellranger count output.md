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

Each file or folder structure:
1.raw_feature_bc_matrix and filtered_feature_bc_matrix folder:single cell raw and filtered matrix

      Each folder contains 3 files:barcodes.tsv.gz,features.tsv.gz,matrix.mtx.gz,three files can represent one big MEX matrix like:
      
      feature_type     gene         feature_id      AAACCCAAGGAGAGTA-1  AAACGCTTCAGCCCAG-1  ...
      Gene Expression  MIR1302-2HG  ENSG00000243485 0                   0
      Gene Expression      FAM138A  ENSG00000237613 0                   0
      Gene Expression        OR4F5  ENSG00000186092 0                   0
      Gene Expression   AL627309.1  ENSG00000238009 0                   0
      Gene Expression   AL627309.3  ENSG00000239945 0                   0
      ...


      barcodes.tsv.gz:Barcode sequences correspond to column indices:
      The format of barcodes.tsv.gz file is similar to:
      AAACCCAAGAAACCCA-1
      AAACCCAAGAAACTCA-1
      AAACCCAAGAAATTCG-1
      AAACCCAAGAACTGAT-1
      AAACCCAAGAAGAAAC-1
      ...
      barcodes file has 1 column
      1.Droplet barcode sequence,each barcode sequence includes a suffix with a dash separator followed by a number(e.g.'AAACCCAAGAAACCCA-1' is barcode of one droplet)

      features.tsv.gz:Features correspond to row indices
      The format of features.tsv.gz file is similar to:
      ENSG00000290825 DDX11L2         Gene Expression
      ENSG00000243485 MIR1302-2HG     Gene Expression
      ENSG00000237613 FAM138A         Gene Expression
      ENSG00000290826 ENSG00000290826 Gene Expression
      ENSG00000186092 OR4F5           Gene Expression
      ...
      features file has 3 columns:
      1.Gene ensembl name(e.g.'ENSG00000290825',it corresponds to gene_id in the annotation field of the reference GTF)
      2.Gene symbol name(e.g.'DDX11L2',corresponds to gene_name in the annotation field of the reference GTF,if gene_name is omit,corresponds to gene_id)
      3.Type of feature(e.g.'Gene Expression',it can be one of 'Gene Expression','Antibody Capture','CRISPR Guide Capture','Multiplexing Capture','CUSTOM')

      matrix.mtx.gz:Expression matrix data
      The format of matrix.mtx.gz file is similar to:
      %%MatrixMarket matrix coordinate integer general
      %metadata_json: {"software_version": "Cell Ranger cellranger-8.0.1", "format_version": 2}
      38606 340411 5370658
      2057 2 1
      4622 2 1
      ...
      The first and second lines are information of data and software
      Third line is the number of rwos,columns and matrix data number except '0'
      Forth to last lines have 3 columns:
      1.row index
      2.column index
      3.data
      (e.g.'2057 2 1' means '1' is at '2057' row,'2' column of matrix)

We can use R,python,or cellranger mat2csv to convert MEX matrix into normal matrix
Firstly,we use R to read matrix
~~~{R}
library(Matrix)

matrix_dir = "./raw_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
~~~

output is :
~~~{R}
print(mat[1:6,1:6]) #show top 6 rows and 6 columns data
6 x 6 sparse Matrix of class "dgTMatrix"
                AAACCCAAGAAACCCA-1 AAACCCAAGAAACTCA-1 AAACCCAAGAAATTCG-1 AAACCCAAGAACTGAT-1 AAACCCAAGAAGAAAC-1 AAACCCAAGAAGAGTA-1
ENSG00000290825                  .                  .                  .                  .                  .                  .
ENSG00000243485                  .                  .                  .                  .                  .                  .
ENSG00000237613                  .                  .                  .                  .                  .                  .
ENSG00000290826                  .                  .                  .                  .                  .                  .
ENSG00000186092                  .                  .                  .                  .                  .                  .
ENSG00000238009                  .                  .                  .                  .                  .                  .
print(dim(mat)) #the dimsions of data
[1]  38606 340411

note:'.' means this position is '0'
~~~

      







