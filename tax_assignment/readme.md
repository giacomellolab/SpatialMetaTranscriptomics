## Prepare taxonomic database
1. mkdir \<path\>/accession2taxid/
2. Download the NCBI taxonomy database to: \<path\>/accession2taxid
3. Install taxonkit [https://bioinf.shenwei.me/taxonkit/] and make sure it is availalbe on PATH
4. Download https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz and https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz to \<path\>/accession2taxid/
5. cat \<path\>\/accession2taxid/*.gz \> \<path\>\/accession2taxid\/nucl_gb_wgs.accession2taxid.gz # join files
6. perl SpatialMetaTranscriptomics\/tax_assignment\/add_full_tax_rank_to_accession2taxid.pl \<path\>/accession2taxid/nucl_gb_wgs.accession2taxid \<path\>/accession2taxid/nucl_gb_wgs.accession2taxid.with_full_lineage \<path\>/accession2taxid/ 1> \<path\>/accession2taxid/nucl_gb_wgs.accession2taxid.with_full_lineage.std 2>&1
7. 
  
