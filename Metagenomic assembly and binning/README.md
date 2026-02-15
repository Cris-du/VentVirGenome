# Metagenomic assembly and binning

---
## Install dependencies  
### 为了reads质控  
`fastp v0.23.3`,相关配置方法可参照[fastp](https://github.com/OpenGene/fastp)  

### 为了拼接contigs与过滤≥ 1kb contigs  
`MEGAHIT v1.2.9`,相关配置方法可参照[megahit](https://github.com/voutcn/megahit)  
`seqkit v2.10.1`,相关配置方法可参照[seqkit](https://github.com/shenwei356/seqkit)  

### 为了binning  
`samtools 1.19.2`,相关配置方法可参照[samtools](https://github.com/samtools/samtools)  
`BBMap 38.18`,相关配置方法可参照[bbmap](https://github.com/BioInfoTools/BBMap?tab=readme-ov-file)  
`MetaBAT2 2.15`,相关配置方法可参照[metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)  
`MaxBin 2.2.7`,相关配置方法可参照[maxbin](https://sourceforge.net/projects/maxbin/)  
`metaWRAP v=1.3.2`,相关配置方法可参照[metawrap](https://github.com/bxlab/metaWRAP)  

### 为了微生物基因组分类  
`GTDB-Tk v2.4.0`,相关配置方法可参照[GTDBTk](https://github.com/Ecogenomics/GTDBTk)  

### 你需要可以运行以下命令  
`fastp`  
`megahit`  
`seqkit`  
`samtools`  
`bbmap.sh`  
`jgi_summarize_bam_contig_depths`  
`run_MaxBin.pl`  
`metabat2`  
`metaWRAP`  
`gtdbtk`  

## 执行操作  
### contigs拼接  
使用`fastp`进行reads质控（包括修剪, 去除接头）  
准备已经下载好的宏基因组测序reads文件  
单端测序：`run_id_single_reads.fq.gz`  
双端测序：`run_id_forward_reads.fq.gz`, `run_id_reverse_reads.fq.gz`  

合并相同样本的不同测序reads文件  
单端文件  
```
cat ./sample_id/*_single_reads.fq.gz > sample_id_merge_single_reads.fq.gz
```
双端文件  
```
cat ./sample_id/*_forward_reads.fq.gz > sample_id_merge_forward_reads.fq.gz, cat ./sample_id/*_reverse_reads.fq.gz > sample_id_merge_reverse_reads.fq.gz

```
单端  
```
fastp -i sample_id_merge_single_reads.fq.gz -o sample_id_merge_single_fastped_reads.fq.gz -h sample_id_merge_single_fastped_reads_report.html -Q --thread=20 --length_required=15 --n_base_limit=5 --compression=6
```
双端  
```
fastp -i sample_id_merge_forward_reads.fq.gz -I sample_id_merge_reverse_reads.fq.gz -o sample_id_merge_forward_fastped_reads.fq.gz -O sample_id_merge_reverse_fastped_reads.fq.gz --unpaired1 sample_id_merge_forward_unpaired_fastped_reads.fq.gz --unpaired2 sample_id_merge_reverse_unpaired_fastped_reads.fq.gz -h sample_id_merge_paired_fastped_reads_report.html -Q --thread=20 --length_required=15 --n_base_limit=5 --compression=6
```

使用`megahit`进行contig拼接  
单端测序文件拼接  
```
megahit --continue --min-count 2 --k-min 21 --k-max 255 --k-step 4 -t 20 -r sample_id_merge_single_fastped_reads.fq.gz -o ./sample_id_megahit
```
双端测序文件拼接  
```
megahit --continue --min-count 2 --k-min 21 --k-max 255 --k-step 4 -t 20 -1 sample_id_merge_forward_fastped_reads.fq.gz -2 sample_id_merge_reverse_fastped_reads.fq.gz -o ./sample_id_megahit
```
### binning分箱  
过滤出长度≥1kb的contigs  
```
seqkit seq -g -j 20 -m 1000 ./sample_id_megahit/final.contigs.fa > sample_id_filter_1kb_contigs.fa
```
将每个样品中长度≥1kb的contigs映射至相应样本测序文件,并且为了`maxbin`分箱,进行bbmap覆盖深度计算  
单端
```
bbmap.sh in=sample_id_merge_single_fastped_reads.fq.gz ref=sample_id_filter_1kb_contigs.fa nodisk k=15 minid=0.9 keepnames=t covstats=sample_id_depth_bbmap.txt minaveragequality=5 outm=sample_id_bbmap.sam threads=64
```
双端
```
bbmap.sh in=sample_id_merge_forward_fastped_reads.fq.gz in2=sample_id_merge_reverse_fastped_reads.fq.gz ref=sample_id_filter_1kb_contigs.fa nodisk k=15 minid=0.9 keepnames=t covstats=sample_id_depth_bbmap.txt minaveragequality=5 outm=sample_id_bbmap.sam threads=64
```

sam文件转换为bam文件,并进行排序与索引  
转换  
```
samtools view -bS -h -@ 4 sample_id_bbmap.sam -o sample_id_bbmap.bam
```
排序与索引  
```
samtools sort -@ 4 -o sample_id_sorted_bbmap.bam sample_id_bbmap.bam && samtools index sample_id_sorted_bbmap.bam
```
为了`metabat2`分箱,进行jgi_depth覆盖深度计算  
```
jgi_summarize_bam_contig_depths --outputDepth sample_id_jgi_depth.txt sample_id_sorted_bbmap.bam
```

`maxbin`与`metabat2`分箱  
maxbin  
```
run_MaxBin.pl -contig sample_id_filter_1kb_contigs.fa -abund sample_id_depth_bbmap.txt -min_contig_length 1000 -thread 64 -out ./sample_id_maxbin_bins
```
metabat2  
```
metabat2 -i sample_id_filter_1kb_contigs.fa -a sample_id_jgi_depth.txt -m 1500 -v --cvExt -o ./sample_id_metabat2_bins -t 64
```

使用`metawrap`对`maxbin`与`metabat2`分箱结果进行精炼  
```
metaWRAP bin_refinement -t 40 -c 50 -x 5 -o ./sample_id_metawrap_bins -A ./sample_id_maxbin_bins -B ./sample_id_metabat2_bins --keep-ambiguous
```

使用`metawrap`对仅在`maxbin`或`metabat2`成功分箱的sample_id的bin进行质控  
```
metaWRAP bin_refinement -t 40 -c 50 -x 5 -o ./sample_id_metabat2(maxbin)_bins_metawrap -A ./sample_id_metabat2(maxbin)_bins --keep-ambiguous
```

人工筛选完整度≥50%且污染度≤5%的bin,作为GOHMGD  

### 微生物MAG分类   
安装`gtdb-tk`数据库  
```
gtdbtk download-data --data-dir ./gtdbtk_db --batch 4
```

将GOHMGD的`bin_id.fna`存放于同一目录内,如`./GOHMGD/*bin_id.fna`以进行gtdbtk分类  
```
gtdbtk classify_wf --genome_dir ./GOHMGD --out_dir ./GOHMGD_gtdbtk_skip --data_dir ./gtdbtk_db --skip_ani_screen -x fa --cpus 10 --pplacer_cpus 10
```
