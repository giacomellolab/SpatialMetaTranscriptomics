use strict;
use List::MoreUtils qw(uniq);
use warnings;
use Storable;
use IO::Zlib;
use IO::Compress::Gzip qw(gzip $GzipError);

if (@ARGV<3) {die "USAGE: perl $0 <m8_file> <NCBI_accssion2lineage_hash> <out>\n";}
my $m8_file=shift;
my $NCBI_map=shift;
my $out_prefix=shift;

my $out_mmseq_with_tax=$out_prefix."_With_tax";
my $out_read_tax_profile=$out_prefix."_Tax_profile"; # read_id\t%id\t%coverage\tBacteria\tEuk\tTAIR10\tgenus_list
my $out_uniq_tax_profile; #TODO

my $NCBI_tax_Map_hashRef = retrieve($NCBI_map);

=cut
open (my $NCBI_TAX,"<",$NCBI_map) || die "Can't open NCBI_TAX '$NCBI_map' $!";
# primaryAccession        start   stop    path    organism_name   taxid
my $heaer=<$NCBI_TAX>;
while (my $line=<$NCBI_TAX>)
{
	chomp ($line);
	my ($accession,$accession_version,$taxid,$gi)=split(/\t/,$line);
	$NCBI_tax_Map{$accession_version}=$taxid;
}
close ($NCBI_TAX);
=cut
my $gzip="NO";
my $BLAST_M8;
my $OUT;
if ($m8_file=~/.gz$/)
{
	$BLAST_M8 = new IO::Zlib;
	$BLAST_M8->open($m8_file, "rb") || die "Can't open BLAST_M8 '$m8_file' $!";
	$gzip="YES";
	$out_mmseq_with_tax=$out_mmseq_with_tax.".gz";
	$OUT = IO::Compress::Gzip->new($out_mmseq_with_tax) or die "gzip failed to open OUT '$out_mmseq_with_tax': $GzipError\n";
}
else
{
	open ($BLAST_M8,"<",$m8_file) || die "Can't open BLAST_M8 '$m8_file' $!";
	open ($OUT,">",$out_mmseq_with_tax) || die "Can't open OUT '$out_mmseq_with_tax' $!";
}

while (my $line=<$BLAST_M8>)
{
	chomp ($line);
	my @line=split(/\t/,$line);
	my $hit_id=$line[1];
	my $read_id=$line[0];
	my $percent_id=$line[2];
	my $Qcover=abs($line[13]);
	my ($accession,$start,$end)=split(/\./,$hit_id);
	my $tax_id_and_lineage="NA\tNA";
	if (exists $NCBI_tax_Map_hashRef->{$hit_id})
	{
		$tax_id_and_lineage=$NCBI_tax_Map_hashRef->{$hit_id};
	}
	print $OUT "$line\t$tax_id_and_lineage\n";
}
if ($gzip eq "NO")
{
	close ($BLAST_M8);
	close ($OUT);
}
else
{
	$BLAST_M8->close;
	$OUT->close;	
}

	
	
# TILL HERE 
=cut
my $Ath_lineage="Eukaryota;Archaeplastida;Chloroplastida;Charophyta;Phragmoplastophyta;Streptophyta;Embryophyta;Tracheophyta;Spermatophyta;Magnoliophyta;Brassicales;Arabidopsis;Arabidopsis thaliana (thale cress)";

my %reads_profile=(); # key: read ID; key2: {Bacteria, Eukaryota ,Archaea, TAIR10, genus}; value: hits count for each or uniq list for genus
open (my $BLAST_M8,"<",$m8_file) || die "Can't open BLAST_M8 '$m8_file' $!";
open (my $OUT,">",$out_mmseq_with_tax) || die "Can't open OUT '$out_mmseq_with_tax' $!";
while (my $line=<$BLAST_M8>)
{
	chomp ($line);
	my @line=split(/\t/,$line);
	my $hit_id=$line[1];
	my $read_id=$line[0];
	my $percent_id=$line[2];
	my $Qcover=abs($line[13]);
	my ($accession,$start,$end)=split(/\./,$hit_id);
	if (!exists $reads_profile{$read_id}{Qcover})
	{
		$reads_profile{$read_id}{Qcover}=[($Qcover)];
	}
	else
	{
		push (@{$reads_profile{$read_id}{Qcover}},$Qcover);
		#if ($reads_profile{$read_id}{Qcover}!=$Qcover)
		#{
		#	print "[WARN] for read id '$read_id' there is an eaual hit with different Qcover ($reads_profile{$read_id}{Qcover}!=$Qcover) -- ignore the second one...\n";
		#}
	}
	if (!exists $reads_profile{$read_id}{Percent_ID})
	{
		$reads_profile{$read_id}{Percent_ID}=[($percent_id)];
	}
	else
	{
		push (@{$reads_profile{$read_id}{Percent_ID}},$percent_id);
		#if ($reads_profile{$read_id}{Percent_ID}!=$percent_id)
		#{
		#	print "[WARN] for read id '$read_id' there is an eaual hit with different percent identity ($reads_profile{$read_id}{Percent_ID}!=$percent_id) -- ignore the second one...\n";
		#}
	}
	
	
	if (exists $SILVA_Map{$accession})
	{
		my @lineage=split(/;/,$SILVA_Map{$accession});
		my $SuperKingdom=$lineage[0];
		my $Genus=$lineage[scalar(@lineage)-2];
		if (exists $reads_profile{$read_id}{$SuperKingdom})
		{
			$reads_profile{$read_id}{$SuperKingdom}++;
		}
		else
		{
			$reads_profile{$read_id}{$SuperKingdom}=1;
		}
		print $OUT "$line\t$SILVA_Map{$accession}\t$SuperKingdom\t$Genus\n";
		if (exists $reads_profile{$read_id}{"genus"}) {push (@{$reads_profile{$read_id}{"genus"}},$Genus)}
		else {$reads_profile{$read_id}{"genus"}=[($Genus)]}
	}
	elsif ($accession=~/^Chr[12345CM]$/) # TAIR10
	{
		print $OUT "$line\t$Ath_lineage\tEukaryota\tTAIR10\n";
		if (exists $reads_profile{$read_id}{"Eukaryota"})
		{
			$reads_profile{$read_id}{"Eukaryota"}++;
		}
		else
		{
			$reads_profile{$read_id}{"Eukaryota"}=1;
		}
		if (exists $reads_profile{$read_id}{"TAIR10"})
		{
			$reads_profile{$read_id}{"TAIR10"}++;
		}
		else
		{
			$reads_profile{$read_id}{"TAIR10"}=1;
		}
		if (exists $reads_profile{$read_id}{"genus"}) {push (@{$reads_profile{$read_id}{"genus"}},"TAIR10")}
		else {$reads_profile{$read_id}{"genus"}=[("TAIR10")]}
	}
	else
	{
		print $OUT "$line\tNA\tNA\n";
	}
}
close ($BLAST_M8);
close ($BLAST_M8);

my @tax_super_kingdoms=("Bacteria","Eukaryota" ,"Archaea", "TAIR10");
open (my $READS_PROF,">",$out_read_tax_profile) || die "Can't open READS_PROF '$out_read_tax_profile' $!";
# print the tax sum report
print $READS_PROF "read_id\tpercent_identity\tcoverage\t",join("\t",@tax_super_kingdoms),"\tgenus\n";
foreach my $read_id (keys %reads_profile)
{
	my $read_Qcover=join(";",sort {$b<=>$a} (uniq @{$reads_profile{$read_id}{Qcover}}));
	my $read_Percent_ID=join(";",sort {$b<=>$a} (uniq @{$reads_profile{$read_id}{Percent_ID}}));
	
	print $READS_PROF "$read_id\t$read_Percent_ID\t$read_Qcover";
	foreach my $kingdom (@tax_super_kingdoms)
	{
		if (exists $reads_profile{$read_id}{$kingdom})
		{
			print $READS_PROF "\t$reads_profile{$read_id}{$kingdom}";
		}
		else
		{
			print $READS_PROF "\t0";
		}
	}
	my @uniq_genus=sort(uniq(@{$reads_profile{$read_id}{"genus"}}));
	print $READS_PROF "\t",join(";",@uniq_genus),"\n";
}
close ($READS_PROF);
=cut 
