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
         ├── filtered_feature_bc_matrix  #single cell filtered matrix
         │   ├── barcodes.tsv.gz
         │   ├── features.tsv.gz
         │   └── matrix.mtx.gz
         ├── filtered_feature_bc_matrix.h5  #single cell filtered matrix(h5 format)
         ├── metrics_summary.csv
         ├── molecule_info.h5
         ├── possorted_genome_bam.bam
         ├── possorted_genome_bam.bam.bai
         ├── raw_feature_bc_matrix  #single cell raw matrix
         │   ├── barcodes.tsv.gz
         │   ├── features.tsv.gz
         │   └── matrix.mtx.gz
         ├── raw_feature_bc_matrix.h5  #single cell raw matrix(h5 format)
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

We can use R,python,or cellranger mat2csv to convert MEX matrix into normal matrix.(ref website:https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices#r-load-mat)  
Firstly,we use R to read matrix:
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

Second,we use python to read the matrix:
~~~{python}
import csv
import gzip
import os
import scipy.io
 
# define MEX directory
matrix_dir = "/your/path/to/raw_feature_bc_matrix"
# read in MEX format matrix as table
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))
 
# list of transcript ids, e.g. 'ENSG00000243485'
features_path = os.path.join(matrix_dir, "features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of gene names, e.g. 'MIR1302-2HG'
gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of feature_types, e.g. 'Gene Expression'
feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

#need to install pandas by using 'pip install pandas'
import pandas as pd

# transform table to pandas dataframe and label rows and columns
matrix = pd.DataFrame.sparse.from_spmatrix(mat)
matrix.columns = barcodes
matrix.insert(loc=0, column="feature_id", value=feature_ids)
matrix.insert(loc=0, column="gene", value=gene_names)
matrix.insert(loc=0, column="feature_type", value=feature_types)
# save the table as a CSV (note the CSV will be a very large file)
#matrix.to_csv("mex_matrix.csv", index=False)
~~~
output is:
~~~{python}
# display matrix
print(matrix)
          feature_type             gene       feature_id  AAACCCAAGAAACCCA-1  ...  TTTGTTGTCTTGGAAC-1  TTTGTTGTCTTTCTAG-1  TTTGTTGTCTTTGAGA-1  TTTGTTGTCTTTGCGC-1
0      Gene Expression          DDX11L2  ENSG00000290825                   0  ...                   0                   0                   0                   0
1      Gene Expression      MIR1302-2HG  ENSG00000243485                   0  ...                   0                   0                   0                   0
2      Gene Expression          FAM138A  ENSG00000237613                   0  ...                   0                   0                   0                   0
3      Gene Expression  ENSG00000290826  ENSG00000290826                   0  ...                   0                   0                   0                   0
4      Gene Expression            OR4F5  ENSG00000186092                   0  ...                   0                   0                   0                   0
...                ...              ...              ...                 ...  ...                 ...                 ...                 ...                 ...
38601  Gene Expression  ENSG00000277836  ENSG00000277836                   0  ...                   0                   0                   0                   0
38602  Gene Expression  ENSG00000278633  ENSG00000278633                   0  ...                   0                   0                   0                   0
38603  Gene Expression  ENSG00000276017  ENSG00000276017                   0  ...                   0                   0                   0                   0
38604  Gene Expression  ENSG00000278817  ENSG00000278817                   0  ...                   0                   0                   0                   0
38605  Gene Expression  ENSG00000277196  ENSG00000277196                   0  ...                   0                   0                   0                   0

[38606 rows x 340414 columns]
~~~
      
Third,we use cellranger mat2csv to read the matrix:
~~~{bash}
# convert from MEX
#cellranger mat2csv your/payh/to/raw_feature_bc_matrix raw_feature_bc_matrix.csv
# or, convert from HDF5
#cellranger mat2csv your/payh/to/raw_feature_bc_matrix.h5 raw_feature_bc_matrix.csv
~~~
Warning:Dense files can be very large and may cause Excel to crash, or even fail in mat2csv if your computer doesn't have enough memory. For example, a feature-barcode matrix from a human reference (~33k genes) with ~3k barcodes uses at least 200MB of disk space. Our 1.3 million mouse neuron dataset, if converted to this format, would use more than 60GB of disk space. Thus, while you can use mat2csv for small datasets, we strongly recommend using R or Python (as shown in the sections above) to examine these matrix files.  

By the way,we also can use bash code to convert MXE into long table format(omit all 'NA' data)(ref website:https://kb.10xgenomics.com/hc/en-us/articles/360023793031-How-can-I-convert-the-feature-barcode-matrix-from-Cell-Ranger-3-to-a-CSV-file)  
~~~{bash}
# Print line number along with contents of barcodes.tsv.gz and genes.tsv.gz 
zcat barcodes.tsv.gz | awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1}' | sort -t, -k 1b,1 > numbered_barcodes.csv
zcat features.tsv.gz | awk -F "\t" 'BEGIN { OFS = "," }; {print NR,$1,$2,$3}' | sort -t, -k 1b,1 > numbered_features.csv

# Skip the header lines and sort matrix.mtx.gz
zcat matrix.mtx.gz | tail -n +4 | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 1b,1 > feature_sorted_matrix.csv
zcat matrix.mtx.gz | tail -n +4 | awk -F " " 'BEGIN { OFS = "," }; {print $1,$2,$3}' | sort -t, -k 2b,2 > barcode_sorted_matrix.csv

# Use join to replace line number with barcodes and genes
join -t, -1 1 -2 1 numbered_features.csv feature_sorted_matrix.csv | cut -d, -f 2,3,4,5,6 | sort -t, -k 4b,4 | join -t, -1 1 -2 4 numbered_barcodes.csv - | cut -d, -f 2,3,4,5,6 > final_matrix.csv

# Remove temp files
rm -f barcode_sorted_matrix.csv feature_sorted_matrix.csv numbered_barcodes.csv numbered_features.csv
~~~

Here is a sample of what final_matrix.csv looks like:

      AAACCTGCACATTAGC-1,ENSG00000005075,POLR2J,Gene Expression,1
      AAACCTGCACATTAGC-1,ENSG00000006015,C19orf60,Gene Expression,1
      AAACCTGCACATTAGC-1,ENSG00000007944,MYLIP,Gene Expression,2
      AAACCTGCACATTAGC-1,ENSG00000008128,CDK11A,Gene Expression,1
      AAACCTGCACATTAGC-1,ENSG00000008952,SEC62,Gene Expression,2
      ...
      It has 5 columns:
      1.10x Genomics cellular barcode
      2.Feature ID
      3.Feature name
      4.Feature type
      5.UMI count
      e.g.'AAACCTGCACATTAGC-1,ENSG00000005075,POLR2J,Gene Expression,1' The POLR2J gene(ENSG00000005075) has one sequence in 'AAACCTGCACATTAGC-1' barcoded cell.

2.raw_feature_bc_matrix.h5 and filtered_feature_bc_matrix.h5 file:H5 format single cell raw and filtered matrix(binary format)  
Except cellranger mat2csv,We can use Seurat in R to load .h5 matrix file.
~~~{R}
BiocManager::install(Seurat)
library(Seurat)
Seurat_matrix = Read10X_h5("raw_feature_bc_matrix.h5")
~~~
output is:
~~~
str(Seurat_matrix)
Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  ..@ i       : int [1:5370658] 2056 4621 7667 10583 18339 20034 20099 20564 21162 23190 ...
  ..@ p       : int [1:340412] 0 0 13 14 15 16 20 21 23 25 ...
  ..@ Dim     : int [1:2] 38606 340411
  ..@ Dimnames:List of 2
  .. ..$ : chr [1:38606] "DDX11L2" "MIR1302-2HG" "FAM138A" "ENSG00000290826" ...
  .. ..$ : chr [1:340411] "AAACCCAAGAAACCCA-1" "AAACCCAAGAAACTCA-1" "AAACCCAAGAAATTCG-1" "AAACCCAAGAACTGAT-1" ...
  ..@ x       : num [1:5370658] 1 1 1 1 1 1 1 1 1 1 ...
  ..@ factors : list()
The data format is dgCMatrix,it has 6 solts:
1.i is row index of the data x
2.p is column index of the data x
3.Dim is dimsion of the data(e.g.the matrix has 38606 rows and 340411 columns)
4.Dimnames is the row names and column names of the data(e.g.row is symbol name of genes,column is barcodes)
5.x is the value
6.factors is a blank list

print(as.matrix(Seurat_matrix)[1:6,1:6])
                AAACCCAAGAAACCCA-1 AAACCCAAGAAACTCA-1 AAACCCAAGAAATTCG-1 AAACCCAAGAACTGAT-1 AAACCCAAGAAGAAAC-1 AAACCCAAGAAGAGTA-1
DDX11L2                          0                  0                  0                  0                  0                  0
MIR1302-2HG                      0                  0                  0                  0                  0                  0
FAM138A                          0                  0                  0                  0                  0                  0
ENSG00000290826                  0                  0                  0                  0                  0                  0
OR4F5                            0                  0                  0                  0                  0                  0
ENSG00000238009                  0                  0                  0                  0                  0                  0
...
Warning message:
In asMethod(object) :
  sparse->dense coercion: allocating vector of size 97.9 GiB
The results is same as previous.
note:The normal matrix uses a lot of memory,usually we don't use that format of matrix to process data,here just to show the data. 
~~~






