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
