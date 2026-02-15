# Virus–Host infective relationship prediction  
---
## Install dependencies  
### 为了识别CRIPSR-SPACER序列  
`CRT-mod version 2.0rev1`,内置CRISPR-spacer识别工具为[CRT](https://www.room220.com/crt/),相关配置方法可参照[CRT-mod](https://github.com/caseyh/crt-mod)  

### 为了进行blastn  
`blastn 2.12.0+`,相关配置方法来源于[blastn+](https://ftp.ncbi.nlm.nih.gov/blastn/executables/blastn+/2.12.0/)  

### 你需要可以运行以下命令  
`blastn`  

## 自定义脚本  
### 标准化CRIPSR-SPACER序列结果  
`trans_format_report.py`  
`extract_spacer_seq.py`  
### 筛选包含≥3个spacers的CRISPR簇  
`filter_3_spacers.py`  
### 过滤blastn结果  
`filter_short_blastn_result.py`  
`filter_long_blastn_result.py`  
### 标准化病毒-宿主分配结果  
`standard_blastn_result.py`  

## 执行操作  
### 构建GOHMGD的CRIPSR-SPACER数据集  
识别GOHMGD中每个MAG基因组`bin_id.fna`的CRIPSR-SPACER序列  
```
python ./crt-mod.py -i bin_id.fna fasta -o bin_id_CRISPRs_raw_report.txt --threads=1
```
提取并标准化GOHMGD的CRIPSR-SPACER序列结果  
```
python ./trans_format_report.py -i bin_id_CRISPRs_raw_report.txt -o bin_id_trans_CRISPRs_report.txt
python ./standard_CRISPRs_raw_report.py -i bin_id_trans_CRISPRs_report.txt -o bin_id_standard_raw_CRISPRs.fna
```
筛选包含≥3个spacer序列的CRISPR簇  
```
python ./filter_3_spacers.py -i bin_id_standard_raw_CRISPRs.fna -o bin_id_standard_more_3_CRISPRs.fna
```

### GOHVGD-GOHMGD病毒-宿主识别  
构建GOHVGD的blastn比对数据库
```
makeblastndb -in GOHVGD_all_viral_contigs.fna -dbtype nucl -out ./GOHVGD_all_viral_contigs_db
```
#### 基于CRIPSR-SPACER blastn进行病毒-宿主识别  
GOHMGD的CRIPSR-SPACER序列与GOHVGD进行blastn比对  
```
blastn -query bin_id_standard_more_3_CRISPRs.fna -db ./GOHVGD_all_viral_contigs_db -task blastn-short -outfmt '6 std qlen slen' -max_target_seqs 50000000 -out bin_id_standard_more_3_CRISPRs_to_GOHVGD_blastn_out.txt -num_threads 1 -dust no
```
过滤CRIPSR-SPACER blastn结果(仅保留SNP≤1的匹配)
```
python ./filter_short_blastn_result.py -i bin_id_standard_more_3_CRISPRs_to_GOHVGD_blastn_out.txt -o filter_bin_id_crt_short_blastn.txt
```
标准化为基于CRIPSR-SPACER blastn的病毒-宿主分配结果  
```
python ./standard_blastn_result.py -i filter_bin_id_crt_short_blastn.txt -o bin_id_crt_virus-host_out.txt
```
#### 基于MAG-bin blastn进行病毒-宿主识别  
GOHMGD的MAG-bin与GOHVGD进行blastn比对  
```
blastn -query bin_id.fna -db ./GOHVGD_all_viral_contigs_db -outfmt '6 std qlen slen' -max_target_seqs 50000000 -out bin_id_long_blastn.txt -num_threads 1 -dust no
```
过滤MAG bin blastn结果(仅保留匹配长度≥2500 bp的匹配)
```
python ./filter_long_blastn_result.py -i bin_id_long_blastn.txt -o filter_bin_id_long_blastn.txt
```
标准化为基于MAG-bin blastn的病毒-宿主分配结果  
```
python ./standard_blastn_result.py -i filter_bin_id_long_blastn.txt -o bin_id_virus-host_out.txt
```
合并所有`bin_id_crt_virus-host_out.txt`与`bin_id_virus-host_out.txt`,去重得到GOHVGD与GOHMGD的最终病毒-宿主分配结果
