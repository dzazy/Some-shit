
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

Genome and annotation data (10X provided Human genome and annotation as an example):  
基因组与注释文件(以10X做好的人类基因组和注释文件为例)：  

1.Download reference data:  
1.参考基因组与注释文件文件下载：  
~~~
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -zxvf refdata-gex-GRCh38-2020-A.tar.gz
~~~

2.Reference data folder structure:
2.参考基因组与注释文件夹结构：
    refdata-gex-GRCh38-2020-A  
      ├── fasta  
          ├── genome.fa     Genome sequence fasta format file 参考基因组序列  
          └── genome.fa.fai Index of genome sequence fasta format file 参考基因组序列索引  
      ├── genes  
          └── genes.gtf     Genome annotation gtf format file 基因注释文件  
      ├── pickle  
          └── genes.pickle  Unknown  
      ├── reference.json    Inpot data and pipeline version 参考序列与注释文件信息及所用软件版本与参数  
      └── star  
          ├── Genome  
          ├── SA  
          ├── SAindex  
          ├── chrLength.txt  
          ├── chrName.txt  
          ├── chrNameLength.txt  
          ├── chrStart.txt  
          ├── exonGeTrInfo.tab  
          ├── exonInfo.tab  
          ├── geneInfo.tab  
          ├── genomeParameters.txt  
          ├── sjdbInfo.txt  
          ├── sjdbList.fromGTF.out.tab  
          ├── sjdbList.out.tab  
          └── transcriptInfo.tab  
    4 directories, 20 files  
