use strict;
use IO::Zlib;
use File::Basename;
use IO::Compress::Gzip qw(gzip $GzipError);

if (@ARGV<5) {die "USAGE: perl $0 <reads_spatial_table> <reads_tax_annotation> <tax_level> <out> <Type {[NA] | Bacteria | ITS>"}
my $reads_spatial_table=shift;
my $reads_tax_annotation=shift;
my $tax_level=shift;
my $out=shift;
my $Type=shift; # [NA], Bacteria, ITS

if (!defined $Type) {$Type="NA";}

# each 'Type' define a set of filters
my %Tax_Levels=(	superkingdom => 0,
	           	phylum => 1, 
    	           	class => 2, 
    	           	order => 3, 
	           	family => 4,
	           	genus => 5,
	           	species => 6);

my $Filter_of_Interest="ALL";
my @ToExclude=();
my @levels_to_account_for_exclude=();
if (uc ($Type) eq "ITS")
{
	$Filter_of_Interest="Eukaryota",

	@ToExclude=("unclassified","Chloroplast","Mitochondria","uncultured","Streptophyta","Chordata","Arthropoda");

	@levels_to_account_for_exclude=("genus","phylum");
}

if (uc ($Type) eq "BACTERIA")
{
	$Filter_of_Interest="Bacteria",
	@ToExclude=("unclassified","Chloroplast","Mitochondria","uncultured");

	@levels_to_account_for_exclude=("genus");
}

# read the reads_tax_annotation and filter un-wanted tax
my %reads_tax=();
my $READS_TAX;
if ($reads_tax_annotation=~/.gz$/)
{
	$READS_TAX = new IO::Zlib;
	$READS_TAX->open($reads_tax_annotation, "rb") || die "Can't open READS_TAX '$reads_tax_annotation' $!";
}
else
{
	open ($READS_TAX,"<",$reads_tax_annotation) || die "Can't open READS_TAX '$reads_tax_annotation' $!";
}
while (my $line=<$READS_TAX>)
{
	chomp ($line);
	# Uniq1;size=832; 1266557 1266557 cellular organisms;Eukaryota;Viridiplantae;Streptophyta;Streptophytina;Embryophyta;Tracheophyta;Euphyllophyta;Spermatophyta;Magnoliopsida;Mesangiospermae;eudicotyledons;Gunneridae;Pentapetalae;Santalales;Viscaceae;Viscum;Viscum minimum     Eukaryota;Streptophyta;Magnoliopsida;Santalales;Viscaceae;Viscum;Viscum minimum
	my ($read_id_and_count,$list_of_tax_id,$LCA_tax_id,$LCA_unformatted,$LCA_formatted)=split(/\t/,$line);
	my ($read_id,$orig_count)=split(/;/,$read_id_and_count);
	my @tax=split(";",$LCA_formatted);
	my $FILTER="NO";
	# start filtering
	if ($Filter_of_Interest ne "ALL")
	{
		foreach my $level (@levels_to_account_for_exclude)
		{
			my $level_field=$Tax_Levels{$level};
			my $tax_to_test=$tax[$level_field];
			foreach my $toExclude (@ToExclude)
			{
				if ($tax_to_test=~/$toExclude/i)
				{
					$FILTER="YES";
					# print "FILTER: $line\t$tax_to_test";<STDIN>;
					last;
				}
			}
		}
	}
	if ($LCA_formatted=~/$Filter_of_Interest/ and $Filter_of_Interest ne "ALL") # of interest
	{
		if ($FILTER eq "NO")
		{
			#print "TAKE: $line\treads_tax{$read_id}=$tax[$Tax_Levels{$tax_level}]";<STDIN>;
			$reads_tax{$read_id}=$tax[$Tax_Levels{$tax_level}];
		}
	}

}
close ($READS_TAX);

# read the spatial info of 
my %tax_spatial=();
my $READS_SPATIAL_TBL;
if ($reads_spatial_table=~/.gz$/)
{
	$READS_SPATIAL_TBL = new IO::Zlib;
	$READS_SPATIAL_TBL->open($reads_spatial_table, "rb") || die "Can't open READS_SPATIAL_TBL '$reads_spatial_table' $!";
}
else
{
	open ($READS_SPATIAL_TBL,"<",$reads_spatial_table) || die "Can't open READS_SPATIAL_TBL '$reads_spatial_table' $!";
}
my %SB_maps=();
my %Tax_count=();

while (my $line=<$READS_SPATIAL_TBL>)
{
	chomp ($line);
	my ($SB,$read_id,$count)=split('\t',$line);
	if (exists $reads_tax{$read_id})
	{
		my $tax=$reads_tax{$read_id};
		if (exists $SB_maps{$SB}{$tax})
		{
			$SB_maps{$SB}{$tax}=$SB_maps{$SB}{$tax}+$count;
		}
		else
		{
			$SB_maps{$SB}{$tax}=$count;
		}
		if (exists $Tax_count{$tax})
		{
			$Tax_count{$tax}=$Tax_count{$tax}+$count;
		}
		else
		{
			$Tax_count{$tax}=$count;
		}
	}
}
my @SB=sort keys %SB_maps;
=full_table_ver
my $line=<$READS_SPATIAL_TBL>;
chomp ($line);
my ($read_id,@SB)=split(';',$line);
while ($line=<$READS_SPATIAL_TBL>)
{
	chomp ($line);
	my ($read_id,@counts)=split(';',$line);
	for (my $i=1;$i<scalar @SB;$i++)
	{
		if (exists $reads_tax{$read_id})
		{
			my $SB=$SB[$i];
			my $tax=$reads_tax{$read_id};
			if (exists $SB_maps{$SB}{$tax})
			{
				$SB_maps{$SB}{$tax}=$SB_maps{$SB}{$tax}+$counts[$i];
			}
			else
			{
				$SB_maps{$SB}{$tax}=$counts[$i];
			}
			if (exists $Tax_count{$tax})
			{
				$Tax_count{$tax}=$Tax_count{$tax}+$counts[$i];
			}
			else
			{
				$Tax_count{$tax}=$counts[$i];
			}
		}
	}
}
=cut
open (my $OUT_SPATIAL_TBL,">",$out) || die "Can't open OUT_SPATIAL_TBL '$out' $!";

print $OUT_SPATIAL_TBL "tax;",join(";",@SB),"\n";
foreach my $tax (sort {$Tax_count{$b}<=>$Tax_count{$a}} keys %Tax_count)
{
	print $OUT_SPATIAL_TBL $tax;
	foreach my $SB (@SB)
	{
		if (exists $SB_maps{$SB}{$tax})
		{
			print $OUT_SPATIAL_TBL ";$SB_maps{$SB}{$tax}";
		}
		else
		{
			print $OUT_SPATIAL_TBL ";0";
		}
	}
	print $OUT_SPATIAL_TBL "\n";
}
close ($OUT_SPATIAL_TBL);
