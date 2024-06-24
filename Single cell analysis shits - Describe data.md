This blogs describe the data we used for 10X genomics data analysis!  
本博客介绍10X genomics 数据分析中用到的文件  

The command line:(Reference website : https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct)  
~~~
cellranger count --id=run_count_1kpbmcs \
   --fastqs=/mnt/home/user.name/yard/run_cellranger_count/pbmc_1k_v3_fastqs \
   --sample=pbmc_1k_v3 \
   --transcriptome=/mnt/home/user.name/yard/run_cellranger_count/refdata-gex-GRCh38-2020-A
~~~
As you can see, the standard 'cellranger count' command have four parameters:  
标准的‘cellranger count’命令包含以下四个参数：  
   1.--id : The project name you named  用户命名的项目名称  
   2.--fastqs : The path of sequence fastq files  10X测序fastq文件路径    
   3.--sample : The sample name you named  用户命名样本名称  
   4.--transcriptome : The path of genome and annotation files  基因组与注释文件路径  
Genome and annotation folder structure (10X provided Human genome and annotation as an example):  
基因组与注释文件(以10X做好的人类基因组和注释文件为例)文件夹结构：  
Download reference data:  
文件下载命令：  
~~~
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz
~~~
