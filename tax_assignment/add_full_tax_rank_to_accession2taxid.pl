use strict;
use Storable;
 
if (@ARGV<3) {die "USAGE: perl $0 <accession2taxid_file_in> <Out> <taxonomic_DB_dir (which contain names.dmp and nodes.dmp)>";}
my $accession2taxid_file_in=shift;
my $Out=shift;
my $taxonomic_DB_dir=shift;
my $TAXON_KIT_EXEC="taxonkit";
my $TaxonKit_template="$TAXON_KIT_EXEC lineage --data-dir $taxonomic_DB_dir -i 3 | $TAXON_KIT_EXEC reformat --data-dir $taxonomic_DB_dir -F -i 4";

# out files
my $Out_hash=$Out.".perlHash";

# tmp files
my $Out_tmp_uniq_taxid=$Out.".uniq_taxid.tmp.txt";
my $Out_tmp_uniq_partial_lineages=$Out.".uniq_partial_lineages.tmp.txt";
my $taxon_kit_out_tmp=$Out.".uniq_lineages_formatted.tmp.txt";

# extract the seq_id and the assigned tax_id and write it to $Out.tmp
=cut
open (my $IN,"<",$accession2taxid_file_in) || die "Can't open IN '$accession2taxid_file_in' $!";
my %uniq_tax=(); # key $tax_id, value: $tax_name
while (my $line=<$IN>)
{
	chomp ($line);
	my ($accession,$accession_version,$taxid,$gi)=split (/\t/,$line);
	$uniq_tax{$taxid}=1;
}
close ($IN);

open (my $OUT_TMP,">",$Out_tmp_uniq_taxid) || die "Can't open OUT_TMP '$Out_tmp_uniq_taxid' $!";
foreach my $tax_id (keys %uniq_tax)
{
	print $OUT_TMP "$tax_id\n";	
}
close ($OUT_TMP);


# get the full lineage with taxon_kit
my $cmd="$TAXON_KIT_EXEC lineage $Out_tmp_uniq_taxid --data-dir $taxonomic_DB_dir -i 1 | awk -F \"\\t\" '\$2!=\"\" {printf (\"\%s\\t\%s\\n\",\$1,\$2)}' | sort | uniq> $Out_tmp_uniq_partial_lineages"; # the partial lineage
print STDERR "[INFO]\t",$cmd,"\n";
system ($cmd);
my $cmd1="$TAXON_KIT_EXEC reformat $Out_tmp_uniq_partial_lineages --data-dir $taxonomic_DB_dir -F -i 2 > $taxon_kit_out_tmp";
print STDERR "[INFO]\t",$cmd1,"\n";
system ($cmd1);

# my $Taxon_kit_cmd="cat $Out_tmp | $TaxonKit_template";
# print STDERR "[INFO]\t",$Taxon_kit_cmd,"\n";
#open (my $OUT_TMP_TAXON_KIT,">",$taxon_kit_out_tmp) || die "Can't open OUT_TMP_TAXON_KIT '$taxon_kit_out_tmp' $!";
#open (TAXON_KIT_RESULTS,$Taxon_kit_cmd."|") || die "Can't exec '$Taxon_kit_cmd' $!";
=cut
open (my $TAXON_KIT_RESULTS,"<",$taxon_kit_out_tmp) || die "Can't open TAXON_KIT_RESULTS '$taxon_kit_out_tmp' $!";
my %TaxId2Lineage=();
while (my $line=<$TAXON_KIT_RESULTS>)
{
	chomp ($line); # 
	my ($tax_id,$lineage,$lineage_formatted)=split(/\t/,$line);
	$lineage_formatted =~ s/\s+/ /g;
	$line=join("\t",($tax_id,$lineage,$lineage_formatted));
	$TaxId2Lineage{$tax_id}=$line;
#	print $OUT_TMP_TAXON_KIT $line,"\n";
	#print "TaxId2Lineage{$seq_id}=$line\n"; #QA
}
#close (TAXON_KIT_RESULTS);
#close ($OUT_TMP_TAXON_KIT);
print "[INFO]]\tFinish reading TaxonKit output '$taxon_kit_out_tmp'\n";

# Go over the input file again and add the species_name\ttax_id\tlineage\tformatted_lineage\n
# Also count the uniq lineages and print them to the file $out_prefix.lineage_only.uniq_count.txt
my %lineage_counter=();
my %uniq_taxid_not_found=(); # Will hold all the tax id not found and their count
my %tax_mapper=(); # an hash that would be stored for further use: key: accession_id; value: tqx_id\tlineage_formatted
open (my $IN,"<",$accession2taxid_file_in) || die "Can't open IN '$accession2taxid_file_in' $!";
my $header=<$IN>;
chomp ($header);
open (my $OUT,">",$Out) || die "Can't open OUT '$Out' $!";
print $OUT "$header\tlineage\tlineage_formated\n";

while (my $line=<$IN>)
{
	chomp ($line);
	my ($tax_id1,$lineage,$lineage_formatted)=("NA","NA","NA");
	
	my ($accession,$accession_version,$tax_id,$gi)=split (/\t/,$line);
	if (exists $TaxId2Lineage{$tax_id}){
		($tax_id1,$lineage,$lineage_formatted)=split(/\t/,$TaxId2Lineage{$tax_id});
	}
	else
	{
		if (exists $uniq_taxid_not_found{$tax_id})
		{
			$uniq_taxid_not_found{$tax_id}++;
		}
		else
		{
			$uniq_taxid_not_found{$tax_id}=1;
		}
	}
	print $OUT "$line\t$lineage\t$lineage_formatted\n";
	if ($lineage_formatted ne "NA")
	{
		$tax_mapper{$accession_version}=$tax_id."\t".$lineage_formatted;
	}
}

# report all the uniq taxid not found
foreach my $tax_id (sort {$uniq_taxid_not_found{$b}<=>$uniq_taxid_not_found{$a}} keys %uniq_taxid_not_found)
{
	print STDERR "[WARN] Could not find 'TaxId2Lineage{$tax_id} ($uniq_taxid_not_found{$tax_id} times)\n";
}
close ($IN);
close ($OUT);
store \%tax_mapper, $Out_hash;
print "[INFO] Save a perl hash with the accession->tax_id\\tformatted_lineage to file: '$Out_hash'\n";


