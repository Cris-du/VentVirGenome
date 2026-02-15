# Viral protein prediction and functional annotation
---
## Install dependencies  
### 为了蛋白质预测  
`prodigal-gv v2.11.0`,相关配置方法可参照[prodigal-gv](https://github.com/apcamargo/prodigal-gv)  

### 为了蛋白质聚类成蛋白质簇  
`cd-hit v4.8.1`,相关配置方法可参照[cdhit](https://github.com/weizhongli/cdhit)  

### 为了进行病毒的蛋白功能注释  
`dram-bio v1.5.0`,相关配置方法可参照[DRAM](https://github.com/WrightonLabCSU/DRAM)  
`virsorter v2.2.3`,相关配置方法可参照[Virsorter2](https://github.com/jiarong/VirSorter2)  

### 为了召回VentVirGenome的原始蛋白质数据集  
`diamond v2.1.8`,相关配置方法可参照[diamond](https://github.com/bbuchfink/diamond)  

### 你需要可以运行以下命令  
`prodigal-gv`  
`cd-hit`  
`virsorter2.sif`  
`DRAM-v.py`  
`diamond`    

## 自定义脚本  
### 为了识别蛋白质簇的代表性蛋白  
`transformat_pcs_report.py`   

## 执行操作  
### 病毒ORF与蛋白质预测  
```
prodigal-gv -i VentVirGenome_all_viral_contigs.fna -o VentVirGenome_all_viral_contigs.gff -a VentVirGenome_all_viral_contigs_protein.faa -d VentVirGenome_all_viral_contigs_orf.fna -f gff -p meta
```  

### 蛋白质聚类为protein cluster
使用`cd-hit`对VentVirGenome所有viral contigs的蛋白质进行聚类  
```
cd-hit -i VentVirGenome_all_viral_contigs_protein.faa -o VentVirGenome_all_viral_contigs_protein_cluster.faa -c 0.6 -aS 0.8 -g 1 -n 4 -d 0 -T 64 -M 0
```

由于原始的聚类结果文件`VentVirGenome_all_viral_contigs_protein_cluster.faa.clstr`内容结构非传统表格文件(tsv,csv等),故提供`transformat_pcs_report.py`脚本以供文件格式转换,用户可酌情使用  
```
python ./transformat_pcs_report.py -i VentVirGenome_all_viral_contigs_protein_cluster.faa.clstr -o VentVirGenome_all_viral_contigs_protein_cluster.tsv
```

### 病毒蛋白功能多样性注释  
由于`DRAM-v`只能接受`virsorter2`的输出结果，故需要对原始VentVirGenome vOTU代表性contig数据集进行`virsorter2`病毒重新预测  
```
virsorter2.sif run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv  --rm-tmpdir --min-score 0 -i VentVirGenome_vOTU_precontigs.fna -w VentVirGenome_vOTU_precontigs_virsorter2_out -j 4 --min-length 0
```
DRAM-v进行病毒蛋白质注释  
```
DRAM-v.py annotate -i ./VentVirGenome_vOTU_precontigs_virsoter2_out/for-dramv/final-viral-combined-for-dramv.fa -v ./VentVirGenome_vOTU_precontigs_virsoter2_out/for-dramv/viral-affi-contigs-for-dramv.tab -o ./VentVirGenome_vOTU_precontigs_virsoter2_dramv_annot --skip_trnascan --threads 40 --min_contig_size 0
```

### 召回VentVirGenome vOTU代表性contig的原始蛋白质数据集功能多样性注释与AMG识别结果    
构建diamond blastp参考数据库  
```
diamond makedb --in ./VentVirGenome_vOTU_precontigs_virsoter2_dramv_annot/genes.faa --db ./VentVirGenome_vOTU_precontigs_dramv_protein_db --threads 4
```
进行diamond blastp召回比对  
```
diamond blastp --query VentVirGenome_vOTU_precontigs_protein.faa --db VentVirGenome_vOTU_precontigs_dramv_protein_db --out VentVirGenome_vOTU_precontigs_protein_map_dramv_out_blastp.txt --outfmt 6 --id 100 --query-cover 100 --subject-cover 100 --evalue 1e-3 --max-target-seqs 5000000 --threads 2
```
进行相同病毒蛋白质名称对应结果提取  
```
cut -f1,2 VentVirGenome_vOTU_precontigs_protein_map_dramv_out_blastp.txt > VentVirGenome_vOTU_precontigs_protein_map_dramv_protein_seqname_map.txt
```
根据`VentVirGenome_vOTU_precontigs_protein_map_dramv_protein_seqname_map.txt`对`./VentVirGenome_vOTU_precontigs_virsoter2_dramv_annot/annotations.tsv`与`VentVirGenome_vOTU_precontigs_raw_amg.txt`进行VentVirGenome结果召回，获得最终VentVirGenome的蛋白质功能注释结果以及AMG注释结果。  
