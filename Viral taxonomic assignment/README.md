# Viral taxonomic assignment  
---
## Install dependencies  
### For viral taxonomic classification  
`MMseqs2 v17.b804f`, refer to [MMseqs2](https://github.com/soedinglab/MMseqs2)  

### You need to be able to run the following command  
`mmseqs`  

Specific configurations are derived from [MetaVR](https://github.com/apcamargo/ictv-mmseqs2-protein-database/tree/master?tab=readme-ov-file#step-5-build-the-mmseqs2-database)  

## Execution  
After building the database according to [MetaVR](https://github.com/apcamargo/ictv-mmseqs2-protein-database/tree/master?tab=readme-ov-file#step-5-build-the-mmseqs2-database), you can perform viral taxonomic classification directly using the following command:  
```
mmseqs easy-taxonomy genomes.fna ictv_nr_db/ictv_nr_db result tmp --blacklist "" --tax-lineage 1
```

The classification result file is `result_lca.tsv`
