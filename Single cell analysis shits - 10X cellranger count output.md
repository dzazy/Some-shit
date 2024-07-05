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
         ├── filtered_feature_bc_matrix  #Contains only detected cell-associated barcodes in MEX format. Each element of the matrix is the number of UMIs associated with a feature (row) and a barcode (column)
         │   ├── barcodes.tsv.gz
         │   ├── features.tsv.gz
         │   └── matrix.mtx.gz
         ├── filtered_feature_bc_matrix.h5  #single cell filtered matrix(h5 format)
         ├── metrics_summary.csv  #The statstics information of the library,contains a number of key metrics about the barcoding and sequencing process.
         ├── molecule_info.h5  #Contains per-molecule information for all molecules that contain a valid barcode, valid UMI, and were assigned with high confidence to a gene or Feature Barcode. This file is a required input to run cellranger aggr.(h5 format)
         ├── possorted_genome_bam.bam  #Indexed BAM file containing position-sorted reads aligned to the genome and transcriptome, as well as unaligned reads, annotated with barcode information. 
         ├── possorted_genome_bam.bam.bai  #Index of possorted_genome_bam.bam file
         ├── raw_feature_bc_matrix  #Contains all detected barcodes in MEX format. Each element of the matrix is the number of UMIs associated with a feature (row) and a barcode (column)
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

3.metrics_summary.csv:The statstics information of the library.  

      The content of the data is:
      Estimated Number of Cells,Mean Reads per Cell,Median Genes per Cell,Number of Reads,Valid Barcodes,Sequencing Saturation,Q30 Bases in Barcode,Q30 Bases in RNA Read,Q30 Bases in UMI,Reads Mapped to Genome,Reads Mapped Confidently to Genome,Reads Mapped Confidently to Intergenic Regions,Reads Mapped Confidently to Intronic Regions,Reads Mapped Confidently to Exonic Regions,Reads Mapped Confidently to Transcriptome,Reads Mapped Antisense to Gene,Fraction Reads in Cells,Total Genes Detected,Median UMI Counts per Cell
      "1,230","54,148","3,286","66,601,887",97.4%,70.8%,94.1%,90.2%,92.7%,96.1%,93.7%,3.7%,31.1%,58.9%,81.4%,7.9%,95.6%,"25,864","9,975"
      We can reshape the data
      Estimated Number of Cells                          "1,230"
      Mean Reads per Cell                                "54,148"
      Median Genes per Cell                              "3,286"
      Number of Reads                                    "66,601,887"
      Valid Barcodes                                     97.4%
      Sequencing Saturation                              70.8%
      Q30 Bases in Barcode                               94.1%
      Q30 Bases in RNA Read                              90.2%
      Q30 Bases in UMI                                   92.7%
      Reads Mapped to Genome                             96.1%
      Reads Mapped Confidently to Genome                 93.7%
      Reads Mapped Confidently to Intergenic Regions     3.7%
      Reads Mapped Confidently to Intronic Regions       31.1%
      Reads Mapped Confidently to Exonic Regions         58.9%
      Reads Mapped Confidently to Transcriptome          81.4%
      Reads Mapped Antisense to Gene                     7.9%
      Fraction Reads in Cells                            95.6%
      Total Genes Detected                               "25,864"
      Median UMI Counts per Cell                         "9,975"
      The data has 19 columns:
      1.Estimated Number of Cells:The number of barcodes associated with cell-containing partitions, estimated from the barcode UMI count distribution.(e.g.The library contains 1230 cell-containing barcodes)
      2.Mean Reads per Cell:The total number of reads divided by the estimated number of cells.
      3.Median Genes per Cell:The median number of genes detected (with nonzero UMI counts) across all cell-associated barcodes.
      4.Number of Reads:Total number of sequenced reads.
      5.Valid Barcodes:Fraction of reads with cell-barcodes that match the whitelist.
      6.Sequencing Saturation:Fraction of reads originating from an already-observed UMI. This is a function of library complexity and sequencing depth. More specifically, this is a ratio where: the denominator is the number of confidently-mapped reads with a valid cell-barcode and valid UMI, and the numerator is the subset of those reads that had a non-unique combination of (cell-barcode, UMI, gene).
      7.Q30 Bases in Barcode:Fraction of bases with Q-score at least 30 in the cell barcode sequences. This is the i7 index (I1) read for the Single Cell 3' v1 chemistry and the R1 read for the Single Cell 3' v2 chemistry.
      8.Q30 Bases in RNA Read:Fraction of bases with Q-score at least 30 in the RNA read sequences. This is Illumina R1 for the Single Cell 3' v1 chemistry and Illumina R2 for the Single Cell 3' v2 chemistry.
      9.Q30 Bases in UMI:Fraction of bases with Q-score at least 30 in the UMI sequences. This is the R2 read for the Single Cell 3' v1 chemistry and the R1 read for the Single Cell 3' v2 chemistry.
      10.Reads Mapped to Genome:Fraction of reads that mapped to the genome.
      11.Reads Confidently Mapped to Genome:Fraction of reads that mapped uniquely to the genome. If a read mapped to exonic loci from a single gene and also to non-exonic loci, it is considered uniquely mapped to one of the exonic loci.
      12.Reads Mapped Confidently to Intergenic Regions:Fraction of reads that mapped to the intergenic regions of the genome with a high mapping quality score as reported by the aligner.
      13.Reads Mapped Confidently to Intronic Regions:Fraction of reads that mapped to the intronic regions of the genome with a high mapping quality score as reported by the aligner.
      14.Reads Mapped Confidently to Exonic Regions:Fraction of reads that mapped to the exonic regions of the genome with a high mapping quality score as reported by the aligner.
      15.Reads Mapped Confidently to Transcriptome:Fraction of reads that mapped to a unique gene in the transcriptome with a high mapping quality score as reported by the aligner. The read must be consistent with annotated splice junctions when include-introns=false. These reads are considered for UMI counting.
      16.Reads Confidently Mapped Antisense:Fraction of reads confidently mapped to the transcriptome, but on the opposite strand of their annotated gene. A read is counted as antisense if it has any alignments that are consistent with an exon of a transcript but antisense to it, and has no sense alignments.
      17.Fraction Reads in Cells:The fraction of cell-barcoded, confidently mapped reads with cell-associated barcodes.
      18.Total Genes Detected:The number of genes with at least one UMI count in any cell.
      19.Median UMI Counts per Cell:The median number of total UMI counts across all cell-associated barcodes.

4.possorted_genome_bam.bam:Indexed BAM file containing position-sorted reads aligned to the genome and transcriptome, as well as unaligned reads, annotated with barcode information. 
The bam file is compressed mapping results sam file.We can use samtools software to read the content.(ref website:https://samtools.github.io/hts-specs/SAMv1.pdf)
~~~{bash}
#conda install samtools  Need conda environment
samtools view possorted_genome_bam.bam | less -S

      A00228:279:HFWFVDMXX:2:1385:2085:18975  16      chr1    13402   0       91M     *       0       0       GAGCCTCCACCACCCCGAGATCACATTTCTCACTGCCTTTTGTCTGCCCAGTTTCACCAGAAGTAGGCCTCTTCCTGACAGGCAGCTGCAC     FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF     NH:i:6  HI:i:1  AS:i:89 nM:i:0  RG:Z:run_count_1kpbmcs:0:1:HFWFVDMXX:2  RE:A:E  AN:Z:ENST00000456328,-649,91M   xf:i:0  CR:Z:TCTGCATGTAACTTAA   CY:Z:FFFFFFFFFFFFFFFF   UR:Z:TACCACAACCAG       UY:Z:FFFFFFFFFFFF       UB:Z:TACCACAACCAG
      A00228:279:HFWFVDMXX:1:2388:30499:30827 1024    chr1    13421   255     91M     *       0       0       ATCACATTTCTCACTGCCTTTTGTCTGCCCAGTTTCACCAGAAGTAGGCCTCTTCCTGACAGGCAGCTGCACCACTGCCTGGCGCTGTGCC     FF,:,F,,,,FF:,FFF:FFFFFFFFFFFFFFF,FFFFFFF:FF:FFFFF,FFFFFFFFFFFF,:FFFF,,FFFFF,F,FFFF:,,FFFFF     NH:i:2  HI:i:2  AS:i:89 nM:i:0  RG:Z:run_count_1kpbmcs:0:1:HFWFVDMXX:1  TX:Z:ENST00000456328,+668,91M   GX:Z:ENSG00000290825    GN:Z:DDX11L2    fx:Z:ENSG00000290825    RE:A:E  mm:i:1  xf:i:17 CR:Z:CACTGAAGTCTTTCTA   CY:Z:FFF:FFFFFF:F,FFF   CB:Z:CACTGAAGTCTTTCTA-1 UR:Z:GCCTTACAGGTC       UY:Z::F,:,FF,FFFF       UB:Z:GCCTTACAGGTC
      A00228:279:HFWFVDMXX:1:2415:3640:8015   1024    chr1    13435   255     1S90M   *       0       0       ATGCCTTTTGTCTGCCCAGTTTCACCAGAAGTATGCATCTTCATGACAGGCAGCTGCACCACTGCCTGGCGCTGTGCCCTTCCTTTGCTCT     ,FFFFF:FFFFF:,F:FF,FFFFFF:F:FF:FF,,F,FFFFF,F:FFFF,FF:FFFFF:,F:F,,:FF,:FFFFF:F:FFF,:FFF,:F::     NH:i:2  HI:i:2  AS:i:82 nM:i:3  RG:Z:run_count_1kpbmcs:0:1:HFWFVDMXX:1  TX:Z:ENST00000456328,+682,1S90M GX:Z:ENSG00000290825    GN:Z:DDX11L2    fx:Z:ENSG00000290825    RE:A:E  mm:i:1  xf:i:17 CR:Z:CACTGAAGTCTTTCTA   CY:Z:FFFFFFFFF,::F,FF   CB:Z:CACTGAAGTCTTTCTA-1 UR:Z:GCCTTACAGGTC       UY:Z:,:FFFFFFF:F:       UB:Z:GCCTTACAGGTC
      A00228:279:HFWFVDMXX:1:2317:2790:31093  0       chr1    13464   255     89M2S   *       0       0       GTAGGCCTCTTCCTGACAGGCAGCTGCACCACTGCCTGGCGCTGTGCCCTTCCTTTGCTCTGCCCGCTGGAGACGGTGTTTGTAATGGGGA     FFFFFFF,,FF,F:FFFFFFFF,FF:FFFFFFFF,,FFFFFFFFFF,FFFFFFF:FF,F:FFFFFFFF:,F:FFFFF,FFF,F,::,::,,     NH:i:2  HI:i:1  AS:i:85 nM:i:1  RG:Z:run_count_1kpbmcs:0:1:HFWFVDMXX:1  TX:Z:ENST00000456328,+711,89M2S GX:Z:ENSG00000290825    GN:Z:DDX11L2    fx:Z:ENSG00000290825    RE:A:E  mm:i:1  xf:i:25 CR:Z:CACTGAAGTCTTTCTA   CY:Z:FFFFF:FFFF:FFF:F   CB:Z:CACTGAAGTCTTTCTA-1 UR:Z:GCCTTACAGGTC       UY:Z:FF,FF:FF:FFF       UB:Z:GCCTTACAGGTC
      A00228:279:HFWFVDMXX:2:2166:6289:9565   1024    chr1    13467   255     91M     *       0       0       GGCCTCTTCCTGACAGGCAGCTGCACCACTGCCTGGCGCTGTGCCCTTCCTTTGCTCTGCCCGCTGGAGACGGTGTTTGTGCTGGGCCTGG     FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FF:,,,:F:F,:,:F     NH:i:3  HI:i:1  AS:i:85 nM:i:2  RG:Z:run_count_1kpbmcs:0:1:HFWFVDMXX:2  TX:Z:ENST00000456328,+714,91M   GX:Z:ENSG00000290825    GN:Z:DDX11L2    fx:Z:ENSG00000290825    RE:A:E  mm:i:1  xf:i:17 CR:Z:CACTGAAGTCTTTCTA   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:CACTGAAGTCTTTCTA-1 UR:Z:GCCTTACAGGTC       UY:Z:FFFFFFFFFFFF       UB:Z:GCCTTACAGGTC
      ...
      The file has 12 columns:
      1.QNAME:Sequence name(e.g.'A00228:279:HFWFVDMXX:2:1385:2085:18975' is the read name)
      2.FLAG:FLAG number,it is a binary number convert to decimal number(e.g.'77' = 000001001101 = 1 + 4 + 8 +64,means 1.PE read,2.the read itself is unmapped,3.its mate is unmapped,4.this is read1)
      abstract the read is paired in sequencing, no matter whether it is mapped in a pair      1
      abstract the read is mapped in a proper pair                                             2
      abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR                  4
      abstract the mate is unmapped                                                            8
      abstract the read is mapped to the reverse strand                                       16
      abstract the mate is mapped to the reverse strand                                       32
      abstract this is read1                                                                  64
      abstract this is read2                                                                 128
      abstract not primary alignment                                                         256
      abstract QC failure                                                                    512
      abstract optical or PCR duplicate                                                     1024
      abstract supplementary alignment                                                      2048
      3.RNAME:Chromosome name(e.g. chr1 is the name of chromosome 1)
      4.POS:left mapped position(e.g.'13402' means this read start mapped to '13402' position of the genome)
      5.MAPQ:mapping quality(-10log10P[mapping position is wrong],'255' indicates that the mapping quality is not available)
      6.CIGAR:
~~~






