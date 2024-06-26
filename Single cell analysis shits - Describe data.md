![image](https://github.com/dzazy/Some-shit/assets/45283502/8f87a94d-698a-47b1-8034-b41ab6e71d24)
This blogs describe the data we used for 10X genomics data analysis!  
本博客介绍10X genomics 数据分析中用到的文件  

The command line:(ref website : https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct)  
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
          ├── exonGeTrInfo.tab          Exon information file
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

1.genome.fa:Genome sequence fasta format file

      The format of genome sequence fasta file is similar to:   
      >chr1 1
      NNN...NNN...TGGCGCAGGC...TTAGGGTTAG...NNN...NNN   
      ">chr1 1" is the name of chromosome 1;
      "NNN...NNN...TGGCGCAGGC...TTAGGGTTAG...NNN...NNN" is sequence of chromosome 1,it contains some 'N' at top and bottom of sequence and four kinds of nucleotides(A,T,C,G) in the middle.  
      This file contains ordered chromosome '1,10-19,2,20-22,3-9,MT,X,Y' and 'KI...,GL...' sequence.  

2.genome.fa.fai:Index of genome sequence fasta format file

      The format of genome sequence fasta fai file is similar to:  
      chr1    248956422       8               60      61  
      chr10   133797422       253105714       60      61  
      chr11   135086622       389133104       60      61  
      chr12   133275309       526471180       60      61  
      ...  
      The fasta fai file have 5 columns:  
      1.Sequence name(e.g. chr1 is the name of first sequence)  
      2.Sequence length(e.g. 248956422 is the length of chr1)
      3.Absolute first nucleotide position of each sequence(e.g. 8=7+1[include blank letter number of '>chr1 1'+newline break],253105714=8[length of '>chr1 1']+248956422[chr1 letter counts]+4149274[line counts of chr1]+10[length of '>chr10 10'])
      4.Nucleotide counts of each line in genome.fa file without blank(e.g. 60 is nucleotide counts of first line of chr1 sequence[without blank])
      5.Nucleotide counts of each line in genome.fa file with blank(e.g. 61 is nucleotide counts of first line of chr1 sequence[with blank])

3.genes.gtf:Genome annotation gtf format file(ref website:https://www.gencodegenes.org/pages/data_format.html)

      The format of genome annotation gtf format file is similar to:
      chr1    HAVANA  gene    29554   31109   .       +       .       gene_id "ENSG00000243485"; gene_version "5"; gene_type"lncRNA"; gene_name "MIR1302-2HG"; level 2; hgnc_id "HGNC:52482"; tag "ncRNA_host"; havana_gene "OTTHUMG00000000959.2";  
      chr1    HAVANA  transcript      29554   31097   .       +       .       gene_id "ENSG00000243485"; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; gene_type "lncRNA"; gene_name "MIR1302-2HG"; transcript_type "lncRNA"; transcript_name "MIR1302-2HG-202"; level 2; transcript_support_level "5"; hgnc_id "HGNC:52482"; tag "not_best_in_genome_evidence"; tag "dotter_confirmed"; tag "basic"; havana_gene "OTTHUMG00000000959.2"; havana_transcript "OTTHUMT00000002840.1";  
      chr1    HAVANA  exon    29554   30039   .       +       .       gene_id "ENSG00000243485"; gene_version "5"; transcript_id "ENST00000473358"; transcript_version "1"; gene_type "lncRNA"; gene_name "MIR1302-2HG"; transcript_type "lncRNA";transcript_name "MIR1302-2HG-202"; exon_number 1; exon_id "ENSE00001947070"; exon_version "1"; level 2; transcript_support_level "5"; hgnc_id "HGNC:52482"; tag "not_best_in_genome_evidence"; tag "dotter_confirmed"; tag "basic"; havana_gene "OTTHUMG00000000959.2"; havana_transcript "OTTHUMT00000002840.1";  
      ...  
      The gtf file have 9 columns:
      1.Chromosome name(e.g. chr1 is the name of chromosome 1)
      2.Annotation source(e.g. HAVANA is an annotate team,HAVANA is HAVANA-Ensembl team)
      3.Annotation type(e.g. gene is a type of annotation,include 'gene','transcript','exon','CDS','Selenocysteine','UTR','start_codon','stop_codon',...)
          3.1 Gene-Gene
          3.2 Transcript-Single or alternative spliced transcripts of each gene
          3.3 exon-Sequence of DNA that is expressed (transcribed) into RNA and often translated into protein.
          3.4 CDS-coding sequences
          3.5 Selenocysteine-SELENON gene
          3.6 UTR-Untranslated regions at top or bottom transcripts
          3.7 start_codon-Initiation signal for translation that is found on a mRNA strand
          3.8 stop_codon-Single nucleotide triplet end protein synthesis
      4.Genomic start location(e.g. 29554 is start position of gene ENSG00000243485)
      5.Genomic end location(e.g. 31109 is stop position of gene ENSG00000243485)
      6.Score(not used)(e.g. usually replaced by '.')
      7.Genomic Strand(e.g. '+' is positive strand,'-' means reverse strand)
      8.Genomic phase(for CDS features)(e.g. usually replaced by '.',One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on.)
      9.Additional information attributes(Semicolon-separated list of tag-value pairs. e.g. gene_id "ENSG00000243485"; gene_version "5"; gene_type"lncRNA"; gene_name "MIR1302-2HG"; level 2; hgnc_id "HGNC:52482"; tag "ncRNA_host"; havana_gene "OTTHUMG00000000959.2";)
          9.1 gene_id-Ensembl gene ID
          9.2 gene_version-Gene version
          9.3 gene_type-Gene type of each feature,include IG_C_gene,IG_C_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_V_gene,IG_V_pseudogene,TR_C_gene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_V_gene,TR_V_pseudogene,lncRNA,protein_coding.(ref website:https://www.gencodegenes.org/pages/biotypes.html)
          9.4 gene_name-Gene symbol name
          9.5 level 2-1 (verified loci),2 (manually annotated loci),3 (automatically annotated loci)
          9.6 hgnc_id-HUGO gene nomenclature committee gene ID.ref website:https://www.genenames.org/
          9.7 tag-Tag(ref website:https://www.gencodegenes.org/pages/tags.html)
          9.8 havana_gene-Gene id in the HAVANA db
          9.9 transcript_id-Ensembl transcript ID
          9.10 transcript_version-Transcript version
          9.11 transcript_type-Transcript type of each feature,include IG_C_gene,IG_C_pseudogene,IG_D_gene,IG_J_gene,IG_J_pseudogene,IG_V_gene,IG_V_pseudogene,TEC,TR_C_gene,TR_D_gene,TR_J_gene,TR_J_pseudogene,TR_V_gene,TR_V_pseudogene,lncRNA,misc_RNA,non_stop_decay,nonsense_mediated_decay,processed_pseudogene,protein_coding,pseudogene,retained_intron,transcribed_unitary_pseudogene,transcribed_unprocessed_pseudogene.(ref website:https://www.gencodegenes.org/pages/biotypes.html)
          9.12 transcript_name-Transcript symbol name
          9.13 transcript_support_level-transcripts are scored according to how well mRNA and EST alignments match over its full length:1 (all splice junctions of the transcript are supported by at least one non-suspect mRNA),2 (the best supporting mRNA is flagged as suspect or the support is from multiple ESTs),3 (the only support is from a single EST),4 (the best supporting EST is flagged as suspect),5 (no single transcript supports the model structure),NA (the transcript was not analyzed)
          9.14 havana_transcript-Transcript id in the HAVANA db
          9.15 exon_number-exon position in the transcript from its 5' end
          9.16 exon_id-Ensembl exon ID
          9.17 exon_version-Exon version
          9.18 ccdsid-offical CCDS ID

4.genes.pickle:A binary file-Unknown  
5.reference.json:Inpot data and pipeline version json file

      The reference.json file is a json format file:
      {
          "fasta_hash": "b6f131840f9f337e7b858c3d1e89d7ce0321b243",
          "genomes": [
              "GRCh38"
          ],
          "gtf_hash": "78ce95ffc520688283c4fe050a27b25b2f45b605",
          "input_fasta_files": [
              "Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified"
          ],
          "input_gtf_files": [
              "gencode.v32.primary_assembly.annotation.gtf.filtered"
          ],
          "mem_gb": 16,
          "mkref_version": "4.0.0",
          "threads": 2,
          "version": "2020-A"
      }
      The genome version is 'GRCh38'
      Genome file is Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified,it is modified from 'Homo_sapiens.GRCh38.dna.primary_assembly.fa' file.
      Annotation file is gencode.v32.primary_assembly.annotation.gtf.filtered,it is filtered from 'gencode.v32.primary_assembly.annotation.gtf',Because 10X library is poly A+ RNA library,the annotation file only retains lncRNA,protein coding,IG and TR features,remove other features like miRNA,snRNA... .
      Software uses 16GB memory.
      mkref software version is '4.0.0'.
      Software uses 2 threads.
      The reference version is '2020-A'.

6.Genome:STAR builded Genome index file-binary file  
7.SA:STAR builded SA file-binary file  
8.SAindex:STAR builded SA index file-binary file   
9.chrLength.txt:Chromosomes length file  

      The format of chrLength.txt file is similar to:
      248956422
      133797422
      135086622
      ...
      The chrLength.txt file have 1 columns:
      1.Nucleotide counts of chromosome(e.g. 248956422 is nucleotide counts of chromosome 1)

10.chrName.txt:Chromosomes name file

      The format of chrName.txt file is similar to:
      chr1
      chr10
      chr11
      ...
      The chrName.txt file have 1 columns:
      1.Chromosome name(e.g. chr1 is name of chromosome 1)  

11.chrNameLength.txt:Chromosomes name and length file

      The format of chrNameLength.txt file is similar to:
      chr1    248956422
      chr10   133797422
      chr11   135086622
      ...
      The chrNameLength.txt file have 2 columns:
      1.Chromosome name(e.g. chr1 is name of chromosome 1)
      2.Nucleotide counts of chromosome(e.g. 248956422 is nucleotide counts of chromosome 1)   

12.chrStart.txt:Chromosome start position

      The format of chrStart.txt file is similar to:
      0
      249036800
      382992384
      ...
      The chrStart.txt file have 1 columns:
      1.Chromosome start position(from 0 to final position)[How to calculate it is unknown]  

13.exonGeTrInfo.tab:Exon information file

      The format of exonGeTrInfo.tab file is similar to:
      1305354
      29553   30038   1       0       0
      30266   30666   1       0       1
      30563   30666   1       0       0
      30975   31096   1       0       0
      30975   31108   1       0       1
      ...
      The first line is exon number
      second to final lines have 5 columns
      1.Exon strat position(gtf genomic start location-1)
      2.Exon stop position(gtf genomic stop location-1)
      3.Exon strand('1' is positive strand,'2' is reverse strand)
      4.Gene number contains the exon(e.g.'29553 30038 1 0 0' the forth column '0' represent gene0-ENSG00000243485;'34553 35173 2 1 2' the forth column '1' represent gene1-ENSG00000237613)
      5.Transcript number contains the exon(e.g.'29553 30038 1 0 0' the fifth column '0' represent transcript0-ENST00000473358;'30266 30666 1 0 1' the fifth column '1' represent transcript0-ENST00000469289)

14.
