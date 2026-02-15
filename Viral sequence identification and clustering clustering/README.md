# Virus Prediction and vOTU Clustering
---
## Install dependencies  
### For Filtering Contigs (≥ 3kb)  
`seqkit v2.10.1`: For configuration details, refer to [seqkit](https://github.com/shenwei356/seqkit)  

### For Initial Virus Prediction  
`genomad v1.7.1`: For configuration details, refer to [genomad](https://github.com/apcamargo/genomad/tree/main)  

### For Viral Genome Quality Assessment  
`checkv v1.0.1`: For configuration details, refer to [checkv](https://bitbucket.org/berkeleylab/checkv/src/master/#markdown-header-checkv-database)  

### For vOTU Clustering  
`anicalc.py` and `aniclust.py`: These scripts are provided by `CheckV`. Refer to the [checkv](https://bitbucket.org/berkeleylab/checkv/src/master/#markdown-header-checkv-database)  

### Ensure the following commands are available in your environment:  
`seqkit`  
`genomad`  
`checkv`  
`makeblastdb`  
`blastn`  

## Custom Scripts  
### Provirus Boundary Identification and Validation  
`check_provirus.py`  

## Operational Procedures  
### Virus Identification  
Filter Contigs (Length ≥ 3kb)  
```
seqkit seq -g -j 20 -m 3000 ./sample_id_megahit/final.contigs.fa > sample_id_filter_3kb_contigs.fa
```

Initial Prediction with `genomad`  
Download the virus marker database `genomad_db`  
```
genomad download-database ./genomad_db
```  
Execute end-to-end prediction:  
```
genomad end-to-end --cleanup -t 8 --splits 8 sample_id_filter_3kb_contigs.fa ./sample_id_step1_genomad ./genomad_db
```

Quality Assessment with CheckV  
```
checkv end_to_end ./sample_id_step1_genomad/sample_id_step1_genomad_summary/sample_id_step1_genomad_virus.fna ./sample_id_step1_genomad_step1_checkv -t 4
```

Refine and Filter Viral Sequences  
Use `check_provirus.py` based on `quality_summary.tsv` to extract non-proviruses, unclipped proviruses, and proviruses with boundaries clipped by CheckV:  
```
python ./check_provirus.py -it ./sample_id_step1_genomad_step1_checkv/quality_summary.tsv -ifc ./sample_id_step1_genomad_step1_checkv/proviruses.fna -ifg ./sample_id_step1_genomad/sample_id_step1_genomad_summary/sample_id_step1_genomad_virus.fna -o sample_id_no_provirus_virus.fna -og sample_id_provirus_part1.fna -ogc sample_id_provirus_part2.fna
```
Secondary Validation for Clipped Proviruses  
Re-run geNomad on clipped provirus sequences to confirm viral signals:  
```
genomad end-to-end --cleanup -t 8 --splits 8 sample_id_provirus_part2.fna ./sample_id_provirus_part2_step2_genomad ./genomad_db
```

Perform a second round of CheckV on validated sequences:  
```
checkv end_to_end ./sample_id_provirus_part2_step2_genomad/sample_id_provirus_part2_step2_genomad_summary/sample_id_provirus_part2_step2_genomad_virus.fna sample_id_step2_genomad_provirus_step2_checkv -t 4
```
Merge `sample_id_no_provirus_virus.fna`,`sample_id_provirus_part1.fna`,`sample_id_provirus_part2_step2_genomad_virus.fna` into `sample_id_virus_final.fna`, The combined results from all samples constitute the `VentVirGenome`  

### vOTU Clustering  
Refer to [checkv](https://bitbucket.org/berkeleylab/checkv/src/master/#markdown-header-checkv-database): Nayfach, S., Camargo, A.P., Schulz, F. et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nat Biotechnol 39, 578–585 (2021).[https://doi.org/10.1038/s41587-020-00774-7](https://doi.org/10.1038/s41587-020-00774-7)  

Merge all viral contigs and create a database `VentVirGenome_all_viral_contigs.fna`:
create blast+ database  
```
makeblastdb -in VentVirGenome_all_viral_contigs.fna -dbtype nucl -out ./VentVirGenome_db
```
Run BLASTn per sample to optimize speed and memory usage. Set max_target_seqs to a sufficiently large value to ensure all hits are captured:  
```
blastn -query sample_id_virus_final.fna -db ./VentVirGenome_db -outfmt '6 std qlen slen' -max_target_seqs 100000000 -o sample_id_blastn.txt -num_threads 64
```
calculate pairwise ANI by combining local alignments between sequence pairs  
```
anicalc.py -i sample_id_blastn.txt -o sample_id_ani.txt
```
Merage all `sample_id_ani.tsv` into `VentVirGenome_ani.txt`  
```
(head -n 1 ./VentVirGenome_ani_dir/*.txt | head -n 1 && tail -n +2 -q ./VentVirGenome_ani_dir/*.txt) > merged.txt
```
Finally, perform CD-HIT-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):  
```
aniclust.py --fna VentVirGenome_all_viral_contigs.fna --ani VentVirGenome_ani.txt --out VentVirGenome_vOTU_clusters.txt --min_ani 95 --min_tcov 85 --min_qcov 0
```
