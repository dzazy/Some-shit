This blogs describe the data we used for count mode of 10X genomics software cellranger !  
本博客介绍10X genomics cellranger count模式用到的文件 ！  

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
    2.--fastqs : The path of 10X single cell sequence fastq files  10X测序单细胞fastq文件路径    
    3.--sample : The sample name you named  用户命名样本名称  
    4.--transcriptome : The path of genome and annotation files  基因组与注释文件路径 
    
Sequence fastq data(10X provided 1,000 PBMC data set as an example):

1.Download 10X single cell sequence fastq files
~~~
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar
~~~
The size of this dataset is 5.17G and takes a few minutes to download.

2.10X single cell sequence fastq folder structure:
```
    pbmc_1k_v3_fastqs
      ├── pbmc_1k_v3_S1_L001_I1_001.fastq.gz
      ├── pbmc_1k_v3_S1_L001_R1_001.fastq.gz
      ├── pbmc_1k_v3_S1_L001_R2_001.fastq.gz
      ├── pbmc_1k_v3_S1_L002_I1_001.fastq.gz
      ├── pbmc_1k_v3_S1_L002_R1_001.fastq.gz
      └── pbmc_1k_v3_S1_L002_R2_001.fastq.gz
    0 directories, 6 files
```
The file name is [Sample Name]_S1_L00[Lane Number] _[Read Type]_001.fastq.gz,it includes 3 parts:(ref website:https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-specifying-fastqs)  
[Sample_Name] is sample name(e.g.'pbmc_1k_v3' is sample name)  
[Lane Number] is sequence lane number(e.g. '1' is lane 1)  
[Read Type] is read type of sequence(e.g.'I1' and 'I2' are Sample index read (optional);'R1' and 'R2' are Read 1 and Read 2)  

fastq file structure:

      The format of fastq file is similar to:
      @A00228:279:HFWFVDMXX:1:1101:8486:1000 2:N:0:NCATTACT
      NACAAAGTCCCCCCCATAATACAGGGGGAGCCACTTGGGCAGGAGGCAGGGAGGGGTCCATTCCCCCTGGTGGGGCTGGTGGGGAGCTGTA
      +
      #FFFFFFFFFFFFFFF:FFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF
      The first line is name of sequence(e.g.'@A00228:279:HFWFVDMXX:1:1101:8486:1000 2:N:0:NCATTACT' is sequence name,it begins with a '@' character and follows sequence identifer and optional description.)
      The second line is sequence
      The third line begin with a '+' character and sometimes followed by the same sequence name
      The forth line is the quality values for the sequence,it must contains the same number of letters with line 2

Sequence structure of different read type file:

      First line:
      I1:@A00228:279:HFWFVDMXX:1:1101:8486:1000 1:N:0:NCATTACT
      R1:@A00228:279:HFWFVDMXX:1:1101:8486:1000 1:N:0:NCATTACT
      R2:@A00228:279:HFWFVDMXX:1:1101:8486:1000 2:N:0:NCATTACT
      I1 and R1 sequence have same sequence names
      R2 sequence name is a little different,it is '2:N:0',I1 and R1 are '1:N:0'
      Second line:
      I1:NCATTACT
      R1:NGTGATTAGCTGTACTCGTATGTAAGGT
      R2:NACAAAGTCCCCCCCATAATACAGGGGGAGCCACTTGGGCAGGAGGCAGGGAGGGGTCCATTCCCCCTGGTGGGGCTGGTGGGGAGCTGTA
      I1 has 8 letters,R1 has 28 letters and R2 has 91 letters
      'I1' is sample index sequence
      'R1' is 12bp barcode and 16 bp UMI sequence(ref website:https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html)
      'R2' is transcript sequence
      Third lines are same
      Forth lines are quality score


      

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
          ├── chrStart.txt              Each chromosome start position  文件中含有每一条染色体在基因组fasta文件中开始的字符的绝对位置
          ├── exonGeTrInfo.tab          Exon information file
          ├── exonInfo.tab              Exon information file
          ├── geneInfo.tab              Gene information file  基因信息文件
          ├── genomeParameters.txt      STAR genomeGenerate mode command and parameters  STAR genomeGenerate模式的命令和参数
          ├── sjdbInfo.txt              Splice junction information file
          ├── sjdbList.fromGTF.out.tab  Splice junction information file
          ├── sjdbList.out.tab          Splice junction information file
          └── transcriptInfo.tab        Transcript information file
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
      second to final lines have 5 columns:
      1.Exon strat position(gtf genomic start location-1)
      2.Exon stop position(gtf genomic stop location-1)
      3.Exon strand('1' is positive strand,'2' is reverse strand)
      4.Gene number contains the exon(e.g.'29553 30038 1 0 0' the forth column '0' represent gene0-ENSG00000243485;'34553 35173 2 1 2' the forth column '1' represent gene1-ENSG00000237613)
      5.Transcript number contains the exon(e.g.'29553 30038 1 0 0' the fifth column '0' represent transcript0-ENST00000473358;'30266 30666 1 0 1' the fifth column '1' represent transcript0-ENST00000469289)

14.exonInfo.tab:Exon information file

      The format of exonInfo.tab file is similar to:
      1305354
      0       485     0
      1010    1113    486
      1422    1543    590
      0       400     0
      709     842     401
      ...
      The first line is exon number
      second to final lines have 3 columns:
      1.Exon start position - Transcript start position(e.g. '0' represents exon ENSE00001947070 start position(29554)-transcript ENST00000473358 start position(29554);'1010' represent exon ENSE00001922571 start position(30564)-transcript ENST00000473358 start position(29554))
      2.Exon end position - Transcript start position(e.g. '485' represents exon ENSE00001947070 end position(30039)-transcript ENST00000473358 start position(29554);'1113' represent exon ENSE00001922571 end position(30667)-transcript ENST00000473358 start position(29554))
      3.Exon-1 length(e.g. '0' represent '486-486';'486' represent '485-0+1';590 represent '486+1113-1010+1')

15.geneInfo.tab:Gene information file

      The format of geneInfo.tab file is similar to:
      36601
      ENSG00000243485
      ENSG00000237613
      ENSG00000186092
      ENSG00000238009
      ENSG00000239945
      ...
      The first line is gene number
      second to final lines are gene ensembl ID

16.transcriptInfo.tab:Transcript information file

      The format of transcriptInfo.tab file is similar to:
      199138
      ENST00000473358 29553   31096   31096   1       3       0
      ENST00000469289 30266   31108   31096   1       2       3
      ENST00000417324 34553   36080   31108   2       3       5
      ENST00000461467 35244   36072   36080   2       2       8
      ENST00000641515 65418   71584   36080   1       3       10
      ...
      The first line is transcript number
      second to final lines have 7 columns:
      1.Transcript ensembl ID
      2.Transcript start position(e.g. 29553 is [start position of transcript ENSG00000243485]-1)
      3.Transcript end position(e.g. 31096 is [end position of transcript ENSG00000243485]-1)
      4.The max gene end position(e.g. first number is max(0,31096)=31096,the second number is max(31096,31096)=31096,the third number is max(31096,31108)=31108,the forth number is max(31108,36080)=36080)
      5.Transcript strand('1' is positive strand,'2' is reverse strand)
      6.exon number(e.g. '3' is exon number of ENST00000473358)
      7.exon accumulate order(e.g. '0' is blank exon number;'3' is ENST00000473358 exon number;'5' is ENST00000473358 and ENST00000469289 sum number 5=3+2)

17.genomeParameters.txt:STAR genomeGenerate mode command and parameters file

      The content of genomeParameters.txt is:
      ### STAR --runMode genomeGenerate --runThreadN 1 --genomeDir /mnt/scratch2/spaceranger/references/GRCh38/star --genomeFastaFiles /mnt/scratch2/spaceranger/references/GRCh38/fasta/genome.fa --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 3 --limitGenomeGenerateRAM 17179869184 --sjdbGTFfile /mnt/scratch2/spaceranger/references/GRCh38/genes/genes.gtf
      versionGenome   20201
      genomeFastaFiles        /mnt/scratch2/spaceranger/references/GRCh38/fasta/genome.fa
      genomeSAindexNbases     14
      genomeChrBinNbits       18
      genomeSAsparseD 3
      sjdbOverhang    100
      sjdbFileChrStartEnd     -
      sjdbGTFfile     /mnt/scratch2/spaceranger/references/GRCh38/genes/genes.gtf
      sjdbGTFchrPrefix        -
      sjdbGTFfeatureExon      exon
      sjdbGTFtagExonParentTranscript  transcript_id
      sjdbGTFtagExonParentGene        gene_id
      sjdbInsertSave  Basic
      1.STAR genomeGenerate mode command
      2.STAR genomeGenerate mode parameters

18.sjdbInfo.txt:Splice junction information file

      The format of sjdbInfo.txt file is similar to:
      362117  100
      30039   30562   1       1       0       1
      30667   30974   1       2       0       1
      35174   35275   2       2       0       2
      35481   35719   2       0       5       2
      65433   65518   1       3       0       1
      ...
      The first line contains two number:
      1.sj number
      2.sjdbOverhang parameter number('100' is sjdbOverhang number)
      second to final lines have 6 columns:
      1.Intron start position(gtf genomic start location,e.g. 30039 is start of intron-1)
      2.Intron end position(gtf genomic stop location,e.g. 30562 is stop of intron-1)
      3.
      4.
      5.
      6.Strand('1' is positive strand,'2' is reverse strand)

19.sjdbList.fromGTF.out.tab:Splice junction information file

      chr1    30040   30563   +       1
      chr1    30668   30975   +       1
      chr1    35175   35276   -       2
      chr1    35482   35720   -       2
      chr1    65434   65519   +       3
      ...
      The sjdbList.fromGTF.out.tab file contains 5 columns:
      1.Chromosome name(e.g. chr1 is the name of chromosome 1)
      2.Intron start position(gtf genomic start location,e.g. 30040 is start of intron)
      3.Intron end position(gtf genomic stop location,e.g. 30563 is stop of intron)
      4.Strand('+' is positive strand,'-' is reverse strand)
      5.Transcript order contains the junction(e.g.The first and second lines belong to the first transcript,the third and forth lines belong to the second transcript)

20.sjdbList.out.tab:Splice junction information file

      chr1    30040   30563   +
      chr1    30668   30975   +
      chr1    35175   35276   -
      chr1    35482   35720   -
      chr1    65434   65519   +
      ...
      The sjdbList.out.tab file contains 4 columns:
      1.Chromosome name(e.g. chr1 is the name of chromosome 1)
      2.Intron start position(gtf genomic start location,e.g. 30040 is start of intron)
      3.Intron end position(gtf genomic stop location,e.g. 30563 is stop of intron)
      4.Strand('+' is positive strand,'-' is reverse strand)




      

      

