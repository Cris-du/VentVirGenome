# Metagenomic assembly and binning

---
## Install dependencies  
### For Reads Quality Control  
`fastp v0.23.3`: For configuration details, refer to [fastp](https://github.com/OpenGene/fastp)  

### for assembly contigs and filter ≥ 1kb contigs  
`MEGAHIT v1.2.9`: For configuration details, refer to [megahit](https://github.com/voutcn/megahit)  
`seqkit v2.10.1`: For configuration details, refer to [seqkit](https://github.com/shenwei356/seqkit)  

### For binning  
`samtools 1.19.2`: For configuration details, refer to [samtools](https://github.com/samtools/samtools)  
`BBMap 38.18`: For configuration details, refer to [bbmap](https://github.com/BioInfoTools/BBMap?tab=readme-ov-file)  
`MetaBAT2 2.15`: For configuration details, refer to [metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)  
`MaxBin 2.2.7`: For configuration details, refer to [maxbin](https://sourceforge.net/projects/maxbin/)  
`metaWRAP v=1.3.2`: For configuration details, refer to [metawrap](https://github.com/bxlab/metaWRAP)  

### For Microbial Genome Classification  
`GTDB-Tk v2.4.0`: For configuration details, refer to [GTDBTk](https://github.com/Ecogenomics/GTDBTk)  

### Ensure the following commands are in your `PATH`  
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

## Operational Procedures  
### Contig Assembly  
Reads Quality Control: Use `fastp` for QC (including trimming and adapter removal).  
Prepare Metagenomic Reads:  
Single-end：`run_id_single_reads.fq.gz`  
Paired-end：`run_id_forward_reads.fq.gz`, `run_id_reverse_reads.fq.gz`  

Merge Reads from the Same Sample  
Single-end:  
```
cat ./sample_id/*_single_reads.fq.gz > sample_id_merge_single_reads.fq.gz
```
Paired-end：  
```
cat ./sample_id/*_forward_reads.fq.gz > sample_id_merge_forward_reads.fq.gz, cat ./sample_id/*_reverse_reads.fq.gz > sample_id_merge_reverse_reads.fq.gz

```
Single-end:  
```
fastp -i sample_id_merge_single_reads.fq.gz -o sample_id_merge_single_fastped_reads.fq.gz -h sample_id_merge_single_fastped_reads_report.html -Q --thread=20 --length_required=15 --n_base_limit=5 --compression=6
```
Paired-end：  
```
fastp -i sample_id_merge_forward_reads.fq.gz -I sample_id_merge_reverse_reads.fq.gz -o sample_id_merge_forward_fastped_reads.fq.gz -O sample_id_merge_reverse_fastped_reads.fq.gz --unpaired1 sample_id_merge_forward_unpaired_fastped_reads.fq.gz --unpaired2 sample_id_merge_reverse_unpaired_fastped_reads.fq.gz -h sample_id_merge_paired_fastped_reads_report.html -Q --thread=20 --length_required=15 --n_base_limit=5 --compression=6
```

Run `megahit` for contig assembly  
Single-end assembly:  
```
megahit --continue --min-count 2 --k-min 21 --k-max 255 --k-step 4 -t 20 -r sample_id_merge_single_fastped_reads.fq.gz -o ./sample_id_megahit
```
Paired-end assembly:  
```
megahit --continue --min-count 2 --k-min 21 --k-max 255 --k-step 4 -t 20 -1 sample_id_merge_forward_fastped_reads.fq.gz -2 sample_id_merge_reverse_fastped_reads.fq.gz -o ./sample_id_megahit
```
### Binning Process  
Filter Contigs (Length ≥ 1kb)  
```
seqkit seq -g -j 20 -m 1000 ./sample_id_megahit/final.contigs.fa > sample_id_filter_1kb_contigs.fa
```
Map Reads to Contigs and Calculate Coverage (BBMap for MaxBin)  
Single-end:
```
bbmap.sh in=sample_id_merge_single_fastped_reads.fq.gz ref=sample_id_filter_1kb_contigs.fa nodisk k=15 minid=0.9 keepnames=t covstats=sample_id_depth_bbmap.txt minaveragequality=5 outm=sample_id_bbmap.sam threads=64
```
Paired-end:
```
bbmap.sh in=sample_id_merge_forward_fastped_reads.fq.gz in2=sample_id_merge_reverse_fastped_reads.fq.gz ref=sample_id_filter_1kb_contigs.fa nodisk k=15 minid=0.9 keepnames=t covstats=sample_id_depth_bbmap.txt minaveragequality=5 outm=sample_id_bbmap.sam threads=64
```

Convert SAM to BAM, Sort, and Index  
```
samtools view -bS -h -@ 4 sample_id_bbmap.sam -o sample_id_bbmap.bam
samtools sort -@ 4 -o sample_id_sorted_bbmap.bam sample_id_bbmap.bam && samtools index sample_id_sorted_bbmap.bam
```
Calculate Coverage Depth (for MetaBAT2)  
```
jgi_summarize_bam_contig_depths --outputDepth sample_id_jgi_depth.txt sample_id_sorted_bbmap.bam
```

Run MaxBin and MetaBAT2  
maxbin  
```
run_MaxBin.pl -contig sample_id_filter_1kb_contigs.fa -abund sample_id_depth_bbmap.txt -min_contig_length 1000 -thread 64 -out ./sample_id_maxbin_bins
```
metabat2  
```
metabat2 -i sample_id_filter_1kb_contigs.fa -a sample_id_jgi_depth.txt -m 1500 -v --cvExt -o ./sample_id_metabat2_bins -t 64
```

Bin Refinement and Quality Control with metaWRAP
Refine MaxBin and MetaBAT2 results:
```
metaWRAP bin_refinement -t 40 -c 50 -x 5 -o ./sample_id_metawrap_bins -A ./sample_id_maxbin_bins -B ./sample_id_metabat2_bins --keep-ambiguous
```

QC for samples where only one tool succeeded:  
```
metaWRAP bin_refinement -t 40 -c 50 -x 5 -o ./sample_id_metabat2(maxbin)_bins_metawrap -A ./sample_id_metabat2(maxbin)_bins --keep-ambiguous
```

Selection Criteria: Manually select bins with Completeness ≥ 50% and Contamination ≤ 5% to serve as VentProkGenome.  

### Microbial MAG Classification   
Install GTDB-Tk Database  
```
gtdbtk download-data --data-dir ./gtdbtk_db --batch 4
```

Run GTDB-Tk Classification
Store VentProkGenome `bin_id.fna` files in a single directory (e.g.,`./VentProkGenome/*bin_id.fna`):
```
gtdbtk classify_wf --genome_dir ./VentProkGenome --out_dir ./VentProkGenome_gtdbtk_skip --data_dir ./gtdbtk_db --skip_ani_screen -x fa --cpus 10 --pplacer_cpus 10
```
