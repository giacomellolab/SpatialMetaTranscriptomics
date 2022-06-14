## Prepare taxonomic database
1. `mkdir <path>/accession2taxid/`
2. Download the [NCBI taxonomy database](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/) to: `<path>/accession2taxid`
3. Install [taxonkit](https://bioinf.shenwei.me/taxonkit/) and make sure it is availalbe on PATH
4. Download https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz and https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz to `<path>/accession2taxid/`
5. `cat <path>/accession2taxid/*.gz > <path>/accession2taxid/nucl_gb_wgs.accession2taxid.gz` # join files
6. `perl SpatialMetaTranscriptomics/tax_assignment/add_full_tax_rank_to_accession2taxid.pl <path>/accession2taxid/nucl_gb_wgs.accession2taxid <path>/accession2taxid/nucl_gb_wgs.accession2taxid.with_full_lineage <path>/accession2taxid/ 1> <path>/accession2taxid/nucl_gb_wgs.accession2taxid.with_full_lineage.std 2>&1`


## Assign taxnomy for each read, UMI filtering and create spatial tables for Microbes
* Install MMSEQS2 form: https://github.com/soedinglab/MMseqs2 together with its `NCBI NT` database
* install usearch from: https://www.drive5.com/usearch/


# taxonomy assignment
* For input reads file: `y.fastq` and spatial information file `x....` 
1. drep reads with usearch: `usearch11.0.667_i86linux32 -fastx_uniques <y.fastq> -fastaout <y.usearch_unique.fasta> -tabbedout <y.usearch_unique.tabbedout> -sizeout -relabel Uniq -strand both`
3. assign taxonomy per read: `mmseqs easy-search --search-type 3 --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen" <y.usearch_unique.fasta>  <nt.DB> <y_vs_NT.mmseq2.m8> /tmp/ --threads 32 --split-memory-limit 250G 
- gzip y_vs_NT.mmseq2.m8
- perl /ebio/abt6_projects7/small_projects/hashkenazy/Scripts/spatial_array_metatranscriptomics/extract_best_hit_mmseq_m8_with_multi_mem_efficient.pl <y_vs_NT_Jan2021.mmseq2.m8.gz> YES
- perl /ebio/abt6_projects7/small_projects/hashkenazy/Scripts/spatial_array_metatranscriptomics/Add_NCBI_tax_annotation_to_mmseq_m9.pl <y_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple.gz> <nucl_gb_wgs.accession2taxid.with_full_lineage.perlHash> <y_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax>

  
