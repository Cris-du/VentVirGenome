# Viral protein prediction and functional annotation
---
## Install dependencies  
### For protein prediction  
`prodigal-gv v2.11.0`, refer to [prodigal-gv](https://github.com/apcamargo/prodigal-gv)  

### For clustering proteins into protein clusters  
`cd-hit v4.8.1`, refer to [cdhit](https://github.com/weizhongli/cdhit)  

### For viral protein functional annotation  
`dram-bio v1.5.0`, refer to [DRAM](https://github.com/WrightonLabCSU/DRAM)  
`virsorter v2.2.3`, refer to [Virsorter2](https://github.com/jiarong/VirSorter2)  

### For recalling the original protein dataset of VentVirGenome  
`diamond v2.1.8`, refer to [diamond](https://github.com/bbuchfink/diamond)  

### You need to be able to run the following commands  
`prodigal-gv`  
`cd-hit`  
`virsorter2.sif`  
`DRAM-v.py`  
`diamond`    

## Custom scripts  
### For identifying representative proteins of protein clusters  
`transformat_pcs_report.py`   

## Execution  
### Viral ORF and protein prediction  
```
prodigal-gv -i VentVirGenome_all_viral_contigs.fna -o VentVirGenome_all_viral_contigs.gff -a VentVirGenome_all_viral_contigs_protein.faa -d VentVirGenome_all_viral_contigs_orf.fna -f gff -p meta
```  

### Protein clustering into protein clusters
Use `cd-hit` to cluster proteins from all VentVirGenome viral contigs:  
```
cd-hit -i VentVirGenome_all_viral_contigs_protein.faa -o VentVirGenome_all_viral_contigs_protein_cluster.faa -c 0.6 -aS 0.8 -g 1 -n 4 -d 0 -T 64 -M 0
```

Since the raw clustering result file `VentVirGenome_all_viral_contigs_protein_cluster.faa.clstr`is not in a traditional table format (tsv, csv, etc.), the `transformat_pcs_report.py`script is provided for format conversion. Users may use it as needed:  
```
python ./transformat_pcs_report.py -i VentVirGenome_all_viral_contigs_protein_cluster.faa.clstr -o VentVirGenome_all_viral_contigs_protein_cluster.tsv
```

### Viral protein functional diversity annotation  
Since `DRAM-v` only accepts output from `virsorter2`, viral re-prediction must be performed on the original VentVirGenome vOTU representative contig dataset:  
```
virsorter2.sif run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv  --rm-tmpdir --min-score 0 -i VentVirGenome_vOTU_precontigs.fna -w VentVirGenome_vOTU_precontigs_virsorter2_out -j 4 --min-length 0
```
DRAM-v viral protein annotation:  
```
DRAM-v.py annotate -i ./VentVirGenome_vOTU_precontigs_virsoter2_out/for-dramv/final-viral-combined-for-dramv.fa -v ./VentVirGenome_vOTU_precontigs_virsoter2_out/for-dramv/viral-affi-contigs-for-dramv.tab -o ./VentVirGenome_vOTU_precontigs_virsoter2_dramv_annot --skip_trnascan --threads 40 --min_contig_size 0
```

### Recalling functional diversity annotation and AMG identification results for VentVirGenome vOTU representative contig original protein dataset    
Build diamond blastp reference database:  
```
diamond makedb --in ./VentVirGenome_vOTU_precontigs_virsoter2_dramv_annot/genes.faa --db ./VentVirGenome_vOTU_precontigs_dramv_protein_db --threads 4
```
Perform diamond blastp recall alignment:  
```
diamond blastp --query VentVirGenome_vOTU_precontigs_protein.faa --db VentVirGenome_vOTU_precontigs_dramv_protein_db --out VentVirGenome_vOTU_precontigs_protein_map_dramv_out_blastp.txt --outfmt 6 --id 100 --query-cover 100 --subject-cover 100 --evalue 1e-3 --max-target-seqs 5000000 --threads 2
```
Extract corresponding results for identical viral protein names:  
```
cut -f1,2 VentVirGenome_vOTU_precontigs_protein_map_dramv_out_blastp.txt > VentVirGenome_vOTU_precontigs_protein_map_dramv_protein_seqname_map.txt
```
Based on `VentVirGenome_vOTU_precontigs_protein_map_dramv_protein_seqname_map.txt`, recall VentVirGenome results from `./VentVirGenome_vOTU_precontigs_virsoter2_dramv_annot/annotations.tsv` and `VentVirGenome_vOTU_precontigs_raw_amg.txt` to obtain the final protein functional annotation and AMG annotation results for VentVirGenome.  
