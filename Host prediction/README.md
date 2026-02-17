# Host prediction  
---
## Install dependencies  
### For identifying CRISPR-SPACER sequences  
`CRT-mod version 2.0rev1`, the built-in CRISPR-spacer identification tool is [CRT](https://www.room220.com/crt/). Refer to [CRT-mod](https://github.com/caseyh/crt-mod) for configuration.  

### For running blastn  
`blastn 2.12.0+`, refer to [blastn+](https://ftp.ncbi.nlm.nih.gov/blastn/executables/blastn+/2.12.0/) for configuration.  

### You need to be able to run the following command  
`blastn`  

## Custom scripts  
### Standardizing CRISPR-SPACER sequence results  
`trans_format_report.py`  
`extract_spacer_seq.py`  
### Filtering CRISPR clusters containing ≥3 spacers  
`filter_3_spacers.py`  
### Filtering blastn results  
`filter_short_blastn_result.py`  
`filter_long_blastn_result.py`  
### Standardizing Virus-Host assignment results  
`standard_blastn_result.py`  

## Execution  
### Building the CRISPR-SPACER dataset for VentProkGenome  
Identify CRISPR-SPACER sequences for each MAG genome `bin_id.fna` in VentProkGenome:  
```
python ./crt-mod.py -i bin_id.fna fasta -o bin_id_CRISPRs_raw_report.txt --threads=1
```
Extract and standardize CRISPR-SPACER sequence results for VentProkGenome:  
```
python ./trans_format_report.py -i bin_id_CRISPRs_raw_report.txt -o bin_id_trans_CRISPRs_report.txt
python ./standard_CRISPRs_raw_report.py -i bin_id_trans_CRISPRs_report.txt -o bin_id_standard_raw_CRISPRs.fna
```
Filter for CRISPR clusters containing ≥3 spacer sequences:  
```
python ./filter_3_spacers.py -i bin_id_standard_raw_CRISPRs.fna -o bin_id_standard_more_3_CRISPRs.fna
```

### VentVirGenome-VentProkGenome Virus-Host identification  
Build the blastn database for VentVirGenome:
```
makeblastndb -in VentVirGenome_all_viral_contigs.fna -dbtype nucl -out ./VentVirGenome_all_viral_contigs_db
```
#### Virus-Host identification based on CRISPR-SPACER blastn  
Perform blastn alignment between VentProkGenome CRISPR-SPACER sequences and VentVirGenome:  
```
blastn -query bin_id_standard_more_3_CRISPRs.fna -db ./VentVirGenome_all_viral_contigs_db -task blastn-short -outfmt '6 std qlen slen' -max_target_seqs 50000000 -out bin_id_standard_more_3_CRISPRs_to_VentVirGenome_blastn_out.txt -num_threads 1 -dust no
```
Filter CRISPR-SPACER blastn results (keep only matches with SNP ≤ 1):
```
python ./filter_short_blastn_result.py -i bin_id_standard_more_3_CRISPRs_to_VentVirGenome_blastn_out.txt -o filter_bin_id_crt_short_blastn.txt
```
Standardize into Virus-Host assignment results based on CRISPR-SPACER blastn:  
```
python ./standard_blastn_result.py -i filter_bin_id_crt_short_blastn.txt -o bin_id_crt_virus-host_out.txt
```
#### Virus-Host identification based on MAG-bin blastn  
Perform blastn alignment between VentProkGenome MAG-bins and VentVirGenome:  
```
blastn -query bin_id.fna -db ./VentVirGenome_all_viral_contigs_db -outfmt '6 std qlen slen' -max_target_seqs 50000000 -out bin_id_long_blastn.txt -num_threads 1 -dust no
```
Filter MAG-bin blastn results (keep only matches with alignment length ≥ 2500 bp):)
```
python ./filter_long_blastn_result.py -i bin_id_long_blastn.txt -o filter_bin_id_long_blastn.txt
```
Standardize into Virus-Host assignment results based on MAG-bin blastn:  
```
python ./standard_blastn_result.py -i filter_bin_id_long_blastn.txt -o bin_id_virus-host_out.txt
```
Merge all `bin_id_crt_virus-host_out.txt` and `bin_id_virus-host_out.txt` files and remove duplicates to obtain the final Virus-Host assignment results for VentVirGenome and VentProkGenome.
