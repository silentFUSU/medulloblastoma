#!/bin/bash
list=(S09-metastatic)

for i in ${list[*]}
do
    Rscript /storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library/numbat/bin/pileup_and_phase.R \
        --label ${i} \
        --samples ${i} \
        --bams /storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/data/rawdata/bam/${i}-rna.bam \
        --barcodes /storage/zhangyanxiaoLab/chaiguoshi/work/medulloblastoma/results/cellranger/${i}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
        --outdir /storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/data/rawdata/bam/${i} \
        --gmap /storage/zhangyanxiaoLab/suzhuojie/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz \
        --eagle /storage/zhangyanxiaoLab/suzhuojie/Eagle_v2.4.1/eagle \
        --snpvcf /storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
        --paneldir /storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/data/1000G_hg38 \
        --ncores 32
done