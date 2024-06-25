
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
```
    refdata-gex-GRCh38-2020-A  
      ├── fasta  
          ├── genome.fa     Genome sequence fasta format file 参考基因组序列  
          └── genome.fa.fai Index of genome sequence fasta format file 参考基因组序列索引  
      ├── genes  
          └── genes.gtf     Genome annotation gtf format file 基因注释文件  
      ├── pickle  
          └── genes.pickle  Unknown  
      ├── reference.json    Inpot data and pipeline version 参考序列与注释文件信息及所用软件版本与参数  
      └── star              files for STAR software  数据比对软件STAR所需文件 
          ├── Genome                    STAR builded Genome index file  STAR构建索引产出的Genome文件
          ├── SA                        STAR builded SA file  STAR构建索引产出的SA文件
          ├── SAindex                   STAR builded SA index file  STAR构建索引产出的SAindex文件
          ├── chrLength.txt             Every chromosomes length(chr1-chr22,X,Y,Mitochondrion and other small parts)  文件中含有每一条线粒体的长度（包括线粒体1-22,X,Y,线粒体及其他片段）  
          ├── chrName.txt               Every chromosomes name(chr1-chr22,X,Y,Mitochondrion and other small parts)  文件中含有每一条线粒体的名称（包括线粒体1-22,X,Y,线粒体及其他片段）
          ├── chrNameLength.txt         Every chromosomes name and length(chr1-chr22,X,Y,Mitochondrion and other small parts)  文件中含有每一条线粒体的名称和长度（包括线粒体1-22,X,Y,线粒体及其他片段）
          ├── chrStart.txt              Each chromosome start position(from 0 to final position)  文件中含有每一条染色体在基因组fasta文件中开始的字符的绝对位置
          ├── exonGeTrInfo.tab          Exon information file(First line is exon number,second to final lines are Exon strat position of gene,Exon stop position of gene,last 3 columns Unknown)
          ├── exonInfo.tab              Exon information file(First line is exon number,second to final lines Unknown)
          ├── geneInfo.tab              Gene information file(First line is gene number,second to final lines are gene ensembl name)  基因信息文件（第一行是总基因数目，第二行至最后是基因ensembl名称）
          ├── genomeParameters.txt      STAR genomeGenerate mode command and parameters  STAR genomeGenerate模式的命令和参数
          ├── sjdbInfo.txt  
          ├── sjdbList.fromGTF.out.tab  
          ├── sjdbList.out.tab  
          └── transcriptInfo.tab        Transcript information file(First line is transcript number,second to final lines are transcript ensembl name,trans start position,trans stop position,Unknown,Unknown,exon number,exon order)
    4 directories, 20 files  
```
Each file structure:  
每个文件的结构如下：

1.genome.fa:   

      The format of genome sequence fasta file is similar to:   
      >chr1 1
      NNN...NNN...TGGCGCAGGC...TTAGGGTTAG...NNN...NNN   
      ">chr1 1" is the name of chromosome 1;
      "NNN...NNN...TGGCGCAGGC...TTAGGGTTAG...NNN...NNN" is sequence of chromosome 1,it contains some 'N' at top and bottom of sequence and four kinds of nucleotides(A,T,C,G) in the middle.  
      This file contains ordered chromosome '1,10-19,2,20-22,3-9,MT,X,Y' and 'KI...,GL...' sequence.  

2.genome.fa.fai：  

      The format of genome sequence fasta fai file is similar to:  
      chr1    248956422       8               60      61  
      chr10   133797422       253105714       60      61  
      chr11   135086622       389133104       60      61  
      chr12   133275309       526471180       60      61  
      ...  
  
