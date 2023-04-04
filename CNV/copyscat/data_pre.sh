#!/bin/bash
echo "begin"
cd /storage/zhangyanxiaoLab/suzhuojie/projects/medulloblastoma/
pwd
list=(S04-recurrent)
for i in ${list[*]}
do
    zcat /storage/zhangyanxiaoLab/chaiguoshi/work/medulloblastoma/results/cellranger_arc/${i}/outs/atac_fragments.tsv.gz | grep "chr" > data/rawdata/atac_fragments/${i}_fragments_chr2.tsv
    gzip data/rawdata/atac_fragments/${i}_fragments_chr2.tsv
    echo "gzip done"
    python3 code/old_data/copyscAT/process_fragment_file.py -i data/rawdata/atac_fragments/${i}_fragments_chr2.tsv.gz -o data/merge_data/copyscat/${i}_fragments_process2.tsv -b 1000000 -f 10 -g data/hg38_chrom_sizes.tsv
done
