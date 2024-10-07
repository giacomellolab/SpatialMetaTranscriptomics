use strict;
use List::MoreUtils qw(uniq);

if (@ARGV<4) {die "USAGE: perl $0 <in_mmseq_with_taxonomy> <out_taxid_hits_per_Q> <out_LCA> <taxonomy_DB_dir> <?out_LCA_lineage_counts>"}

my $in=shift;
my $out_for_lca=shift;
my $LCA_out=shift;
my $taxonomy_DB_dir=shift;
my $LCA_lineage_counts=shift;


my $MMSEQS_M8_Filtered_unclassified_genus=$out_for_lca.".hits.m8";

open (my $IN,"<",$in) || die "Can't open IN '$in' $!";
my $prev_query="NA";
my %classified_hits=();
my %unclassified_hits=();

my %classified_hits_taxid=();
my %unclassified_hits_taxid=();

while (my $line=<$IN>)
{
	chomp $line;
	my ($query,$target,$fident,$alnlen,$mismatch,$gapopen,$qstart,$qend,$tstart,$tend,$evalue,$bits,$qlen,$q_cover,$tax_id,$lineage)=split(/\t/,$line);
	# Uniq66193;size=1;       KM159131.1      1.000   143     0       0       1       143     343     485     7.421E-71       258     143     1.000   77133   Bacteria;unclassified Bacteria phylum;unclassified Bacteria class;unclassified Bacteria order;unclassified Bacteria family;unclassified Bacteria genus;uncultured bacterium
	my @lineage=split(/;/,$lineage);
	if ($lineage[5]=~/unclassified/)
	{
		if (exists $unclassified_hits{$query})
		{
			push (@{$unclassified_hits{$query}},$line);
			push (@{$unclassified_hits_taxid{$query}},$tax_id);
		}
		else
		{
			$unclassified_hits{$query}=[($line)];
			$unclassified_hits_taxid{$query}=[($tax_id)];
		}
	}
	else
	{
		if (exists $classified_hits{$query})
		{
			push (@{$classified_hits{$query}},$line);
			push (@{$classified_hits_taxid{$query}},$tax_id);
		}
		else
		{
			$classified_hits{$query}=[($line)];
			$classified_hits_taxid{$query}=[($tax_id)];
		}
	}
	# Bacteria;unclassified Bacteria phylum;unclassified Bacteria class;unclassified Bacteria order;unclassified Bacteria family;unclassified Bacteria genus;uncultured bacterium
}
close ($IN);

# if genus is unclassified, ignore for LCA
open (my $OUT_MMSEQS_M8_FILTERED_UNCLASSIFIED_GENUS,">",$MMSEQS_M8_Filtered_unclassified_genus) || die "Can't open OUT_MMSEQS_M8_FILTERED_UNCLASSIFIED_GENUS '$MMSEQS_M8_Filtered_unclassified_genus' $!";


open (my $OUT_FOR_LCA,">",$out_for_lca) || die "Can't open OUT_FOR_LCA '$out_for_lca' $!";

foreach my $query (sort by_cluster_size uniq ((keys %classified_hits),(keys %unclassified_hits)))
{
	if (exists $classified_hits{$query})
	{
		print $OUT_MMSEQS_M8_FILTERED_UNCLASSIFIED_GENUS join ("\n",@{$classified_hits{$query}}),"\n";
		print $OUT_FOR_LCA "$query\t",join(",",uniq(@{$classified_hits_taxid{$query}})),"\n";
	}
	else
	{
		print $OUT_MMSEQS_M8_FILTERED_UNCLASSIFIED_GENUS $unclassified_hits{$query}->[0],"\n";
		print $OUT_FOR_LCA "$query\t",join(",",uniq(@{$unclassified_hits_taxid{$query}})),"\n";
	}
}
close ($OUT_FOR_LCA);
close ($OUT_MMSEQS_M8_FILTERED_UNCLASSIFIED_GENUS);

LCA_by_TaxonKit($out_for_lca,$LCA_out,$taxonomy_DB_dir);

if (defined $LCA_lineage_counts)
{
	my %lineage_counter=();
	my $unclassified_lineage_formated="unclassified superkingdom;unclassified phylum;unclassified class;unclassified order;unclassified family;unclassified genus;unclassified species";
	
	# Go over the LCA lineages and create the lineage count sum file
	open (my $LCA,"<",$LCA_out) || die "Can't open LCA '$LCA_out' $!";
	# Uniq10;size=284;        115915,72658,63678,81985,3719,3702      3700    cellular organisms;Eukaryota;Viridiplantae;Streptophyta;Streptophytina;Embryophyta;Tracheophyta;Euphyllophyta;Spermatophyta;Magnoliopsida;Mesangiospermae;eudicotyledons;Gunneridae;Pentapetalae;rosids;malvids;Brassicales;Brassicaceae        Eukaryota;Streptophyta;Magnoliopsida;Brassicales;Brassicaceae;unclassified Brassicaceae genus;unclassified Brassicaceae species
	while (my $line=<$LCA>)
	{
		chomp ($line);
		my ($read_id,$tax_ids,$lca_tax_id,$lineage,$formatted_lineage)=split(/\t/,$line);
		
		my $reads_count=1;
		if ($read_id=~/Uniq\d+;size=\d+/)
		{
			($reads_count)=$read_id=~/Uniq\d+;size=(\d+)/;
		}
		if ($formatted_lineage eq "NA")
		{
			$formatted_lineage=$unclassified_lineage_formated;
		}
		if (exists $lineage_counter{$formatted_lineage})
		{
			$lineage_counter{$formatted_lineage}=$lineage_counter{$formatted_lineage}+$reads_count
		}
		else
		{
			$lineage_counter{$formatted_lineage}=$reads_count;
		}
	}
	close ($LCA);
	open (my $LCA_COUNT,">",$LCA_lineage_counts)  || die "Can't open LCA_COUNT '$LCA_lineage_counts' $!";
	foreach my $lineage (sort {$lineage_counter{$b}<=>$lineage_counter{$a}} keys %lineage_counter)
	{
		print $LCA_COUNT "$lineage_counter{$lineage};$lineage\n";
	}
	close ($LCA_COUNT);
}

sub by_cluster_size
{
	my ($x)=$a=~/size\=(\d+);/;
	my ($y)=$b=~/size\=(\d+);/;
	return $y<=>$x;
}

sub LCA_by_TaxonKit
{
	my $TaxID_lists_file=shift; # /ebio/abt6_projects8/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/raw_data/New_probs_array_Derek_sample_Dec2020/mapped_data/multimapped_STAR_2.7.7e/splitted_by_probs/191029-069_A2_Ar_multimap_bacteria_R2.PROB_1265.unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.for_LCA
	my $out=shift;
	my $taaxonomy_DB=shift;
	my $TAXON_KIT_EXEC="/ebio/abt6_projects7/small_projects/hashkenazy/Programs/taxonkit_072/taxonkit";
	my $std_file=$out.".std";
	
	my $cmd="$TAXON_KIT_EXEC lca --data-dir $taaxonomy_DB -i 2 -s \",\" $TaxID_lists_file 2>$std_file | $TAXON_KIT_EXEC lineage --data-dir $taaxonomy_DB -i 3  2>>$std_file | $TAXON_KIT_EXEC reformat --data-dir $taaxonomy_DB -i 4 -F > $out 2>>$std_file";
	
	print STDERR "[RUN] TaxonKit_cmd: $cmd\n";
	system($cmd);
}