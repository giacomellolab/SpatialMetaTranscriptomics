# Prepare taxonomic database
1. `mkdir <taxonomy_path>/accession2taxid/`
2. Download the [NCBI taxonomy database](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/) to: `<taxonomy_path>/accession2taxid/`
3. Install [taxonkit](https://bioinf.shenwei.me/taxonkit/) and make sure it is availalbe on PATH
4. Download https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz and https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz to `<taxonomy_path>/accession2taxid/`
5. `cat <taxonomy_path>/accession2taxid/*.gz > <taxonomy_path>/accession2taxid/nucl_gb_wgs.accession2taxid.gz` # join files
6. `perl SpatialMetaTranscriptomics/tax_assignment/add_full_tax_rank_to_accession2taxid.pl <taxonomy_path>/accession2taxid/nucl_gb_wgs.accession2taxid <taxonomy_path>/accession2taxid/nucl_gb_wgs.accession2taxid.with_full_lineage <taxonomy_path>/accession2taxid/ 1> <taxonomy_path>/accession2taxid/nucl_gb_wgs.accession2taxid.with_full_lineage.std 2>&1`


# Assign taxnomy for each read, UMI filtering, and creating spatial tables for Microbes
* Install MMSEQS2 form: https://github.com/soedinglab/MMseqs2 together with its `NCBI NT` database
* install usearch from: https://www.drive5.com/usearch/


## taxonomy assignment
* For input reads file: `y.fastq` and spatial information file `x....` 
0. `mkdir out_path`
1. drep reads with usearch:<br>`usearch11.0.667_i86linux32 -fastx_uniques y.fastq -fastaout <out_path>/y.usearch_unique.fasta -tabbedout <out_path>/y.usearch_unique.tabbedout -sizeout -relabel Uniq -strand both`<br><br>
3. assign taxonomy per read:<br> 
`mmseqs easy-search --search-type 3 --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen" <out_path>/y.usearch_unique.fasta <nt.DB> <out_path>/y_vs_NT.mmseq2.m8> /tmp/ --threads 32 --split-memory-limit 250G`<br><br>
`gzip <out_path>/y_vs_NT.mmseq2.m8`<br><br>
`perl SpatialMetaTranscriptomics/tax_assignment/extract_best_hit_mmseq_m8_with_multi_mem_efficient.pl <out_path>/y_vs_NT_Jan2021.mmseq2.m8.gz YES`<br><br>
`perl SpatialMetaTranscriptomics/tax_assignment/Add_NCBI_tax_annotation_to_mmseq_m9.pl <out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple.gz <taxonomy_path>/accession2taxid/nucl_gb_wgs.accession2taxid.with_full_lineage.perlHash <out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple_and_tax`
<br><br> `perl SpatialMetaTranscriptomics/tax_assignment/mmseq_m8_with_tax_LCA_prefer_classified.pl <out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple_and_tax <out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple_and_tax.query_taxid_for_LCA <out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple_and_tax.query_LCA <taxonomy_path>/accession2taxid/ <out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple_and_tax.LCA.uniq_lineage_count.txt 1><out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple_and_tax.LCA.std 2>&1`
<br><br>
`gzip <out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple_and_tax`
<br><br>
`gzip <out_path>/y_vs_NT.mmseq2.m8.best_hit_multiple_and_tax.query_taxid_for_LCA.hits.m8`


  
