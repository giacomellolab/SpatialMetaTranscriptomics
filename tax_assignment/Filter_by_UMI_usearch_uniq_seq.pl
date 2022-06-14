use strict;
use IO::Zlib;
use File::Basename;
use IO::Compress::Gzip qw(gzip $GzipError);

if (@ARGV<5) {die "USAGE: perl $0 <usearch_uniq_tblout_file> <reads_spatial_map> <uniq_seq_modified_count_out_file> <filtered_reads_out_file> <spatial_table_map_out_file>"}
my $usearch_uniq_tblout_file=shift;
my $reads_spatial_map=shift;
#my $uniq_seq_tax_assignment=shift;

my $uniq_seq_modified_count_out_file=shift;
my $filtered_reads_out_file=shift; #gz file
my $spatial_table_map_out_file=shift; # for each spatial pos print the usearhc_uniq_id --> col=SB,row=usearhc_uniq_id -> gz file

# $spatial_table_map_out_file=$spatial_table_map_out_file.".gz";
# SB1-UMI1-Pseudomonas vs. SB1-UMI1-Pseudomonas (both with the same usearch uniq) -> This is duplicate - filter out!
# SB1-UMI1-X vs. SB1-UMI2 -> This is OK
# SB1-UMI1-X vs. SB2-UMI1 -> This is OK

# TO DO: for the filtered, include also the prob information, filter only if the prob is also identical (Now done based on prob files)
# IMPORTANT NOTE: take into account only the reads for which usearch_uniq is availalbe! (assuming these are the only sequences analyzed)

if (!-e $reads_spatial_map) { die "Can't find reads_spatial_map: '$reads_spatial_map'\n"}
if (!-e $usearch_uniq_tblout_file) { die "Can't find usearch_uniq_tblout_file: '$usearch_uniq_tblout_file'\n"}
#if (!-e $uniq_seq_tax_assignment) { die "Can't find uniq_seq_tax_assignment: '$uniq_seq_tax_assignment'\n"}

my %cognates_counter=(); # key1: SB; key2: UMI; key3: usearch_uniq_identifier ; value: count of the cognate
# my %LCA_sign_to_lineage=(); # key1: Kraken_LCA_signature; value: lineage;


# C       VH00203:7:AAAHHTKM5:1:1101:48301:1000   Chloroplast (taxid 44561)       136     44561:56 0:23 46259:5 45332:2 3:2 0:14  Chloroplast     44561   Bacteria;Cyanobacteria;Cyanobacteriia;Chloroplast       Bacteria;Cyanobacteria;Cyanobacteriia;Chloroplast;unclassified Chloroplast family;unclassified Chloroplast genus;unclassified Chloroplast species
# my $KRAKEN_read_id_field=1;
# my $KRAKEN_lineage_field=8;
# my $KRAKEN_LCA_sign_field=4;

# my $genus_filed=5;

my %read_spatial_pos=();
my %read_UMI=();
my %read2uaerch_uniq=();
my %usearch_uniq_size=();
my %usearch_uniq_size_after_filter=();
my %filtered_counters=();
my %SB_maps=();

# read the uniq assignmet by usearch and map to each read
my $USEARCH_UNIQ_TBLOUT;
if ($usearch_uniq_tblout_file=~/.gz$/)
{
	$USEARCH_UNIQ_TBLOUT = new IO::Zlib;
	$USEARCH_UNIQ_TBLOUT->open($usearch_uniq_tblout_file, "rb") || die "Can't open USEARCH_UNIQ_TBLOUT '$usearch_uniq_tblout_file' $!";
	
}
else
{
	open ($USEARCH_UNIQ_TBLOUT,"<",$usearch_uniq_tblout_file) || die "Can't open USEARCH_UNIQ_TBLOUT '$usearch_uniq_tblout_file' $!";
}

while (my $line=<$USEARCH_UNIQ_TBLOUT>)
{
	chomp ($line);
	my ($read_id,$usearch_uniq_id,$cluster_number,$cluster_size)=split(/\t/,$line);
	
	$read2uaerch_uniq{$read_id}=$usearch_uniq_id;
	$usearch_uniq_size{$usearch_uniq_id}=$cluster_size+1;
}
close ($USEARCH_UNIQ_TBLOUT);

%usearch_uniq_size_after_filter=%usearch_uniq_size;

my $DUPLICATED_REMOVED = IO::Compress::Gzip->new($filtered_reads_out_file) or die "gzip failed to open DUPLICATED_REMOVED '$filtered_reads_out_file': $GzipError\n";

my $SPATIAL_MAP;
if ($reads_spatial_map=~/.gz$/)
{
	$SPATIAL_MAP = new IO::Zlib;
	$SPATIAL_MAP->open($reads_spatial_map, "rb") || die "Can't open SPATIAL_MAP '$reads_spatial_map' $!";
	
}
else
{
	open ($SPATIAL_MAP,"<",$reads_spatial_map) || die "Can't open SPATIAL_MAP '$reads_spatial_map' $!";
}
my $filterd_reads=0;
my $all_reads=0;
while (my $line=<$SPATIAL_MAP>)
{
	chomp ($line);
	my ($read_id,$UMI,$SB)=split(/\s+/,$line);
	$read_id=~s/^\@//;
	if (exists $read2uaerch_uniq{$read_id}) # included in the set got later tax assignment
	{
		my $usearch_uniq=$read2uaerch_uniq{$read_id};
		# print "CHECKING: read2uaerch_uniq{$read_id}\n"; # QA
		if (exists $cognates_counter{$SB}{$UMI}{$usearch_uniq})
		{
			# print "SEEN: cognates_counter{$SB}{$UMI}{$usearch_uniq}\n";<STDIN>; # QA
			# update the cognate_counter
			$cognates_counter{$SB}{$UMI}{$usearch_uniq}++;
			# filter and print to the removed file with the SB and UMI
			print $DUPLICATED_REMOVED "$read_id\t$SB\t$UMI\t$usearch_uniq\n";
			$usearch_uniq_size_after_filter{$usearch_uniq}--;
			$filterd_reads++;
			my $cognate=$UMI."\t".$SB."\t".$usearch_uniq;
			if (exists $filtered_counters{$cognate})
			{
				$filtered_counters{$cognate}++;
			}
			else
			{
				$filtered_counters{$cognate}=1;
			}
		}
		else # first seen
		{
			$cognates_counter{$SB}{$UMI}{$usearch_uniq}=1;
			if (!exists $SB_maps{$SB}{$usearch_uniq})
			{
				$SB_maps{$SB}{$usearch_uniq}=1;
			}
			else
			{
				$SB_maps{$SB}{$usearch_uniq}++;
			}
			# else # it could exists with different UMI
			# {
			# 	print "[WARN] SB_maps{$SB}{$usearch_uniq} was already seen...\n";
			# }
		}
		$all_reads++
	}
#	$read_spatial_pos{$read_id}=$pos;
#	$read_UMI{$read_id}=$UMI;
}
close ($SPATIAL_MAP);

print STDERR "[INFO] for $usearch_uniq_tblout_file and $reads_spatial_map total reads filtered by UMI+SB+usearch_uniq: $filterd_reads (",sprintf("%.3f",$filterd_reads/$all_reads),") out of $all_reads reads -> ",$all_reads-$filterd_reads," left.\n";
print "[INFO] for $usearch_uniq_tblout_file and $reads_spatial_map total reads filtered by UMI+SB+usearch_uniq: $filterd_reads (",sprintf("%.3f",$filterd_reads/$all_reads),") out of $all_reads reads -> ",$all_reads-$filterd_reads," left.\n";

# the final output is the updated uniq counts
open (my $UPDATED_UNIQ_COUNT,">",$uniq_seq_modified_count_out_file) || die "Can't open  UPDATED_UNIQ_COUNT '$uniq_seq_modified_count_out_file' $!";
print $UPDATED_UNIQ_COUNT "usearch_uniq_id\tunfiltered_size\tfiltered_size\n";
foreach my $usearch_uniq_id (sort {$usearch_uniq_size_after_filter{$b}<=>$usearch_uniq_size_after_filter{$a}} keys %usearch_uniq_size_after_filter)
{
	print $UPDATED_UNIQ_COUNT "$usearch_uniq_id\t$usearch_uniq_size{$usearch_uniq_id}\t$usearch_uniq_size_after_filter{$usearch_uniq_id}\n";
}
close ($UPDATED_UNIQ_COUNT);

open (my $OUT_SPATIAL_TBL,">",$spatial_table_map_out_file) || die "Can't open OUT_SPATIAL_TBL '$spatial_table_map_out_file' $!";
# my $OUT_SPATIAL_TBL = IO::Compress::Gzip->new($spatial_table_map_out_file) or die "gzip failed to open OUT_SPATIAL_TBL '$spatial_table_map_out_file': $GzipError\n";

# avoid the whole table just print for each SB the usearch_uniq and counts -> SB\tusearch_uniq\tcounts
my %QA_sums=();
foreach my $SB (sort keys %SB_maps)
{
	my $QA_Sum=0;
	foreach my $usearch_uniq (sort {$SB_maps{$SB}{$b}<=>$SB_maps{$SB}{$a}} keys %{$SB_maps{$SB}})
	{
		print $OUT_SPATIAL_TBL "$SB\t$usearch_uniq\t$SB_maps{$SB}{$usearch_uniq}\n";
		# For QA
		if (exists $QA_sums{$usearch_uniq})
		{
			$QA_sums{$usearch_uniq}=$QA_sums{$usearch_uniq}+$SB_maps{$SB}{$usearch_uniq};
		}
		else
		{
			$QA_sums{$usearch_uniq}=$SB_maps{$SB}{$usearch_uniq};
		}
	}	
}

# QA
foreach my $usearch_uniq (sort keys %usearch_uniq_size_after_filter)
{
	if ($QA_sums{$usearch_uniq} != $usearch_uniq_size_after_filter{$usearch_uniq})
	{
		print STDERR "[ERROR] usearch_uniq_size_after_filter{$usearch_uniq}: $QA_sums{$usearch_uniq} != $usearch_uniq_size_after_filter{$usearch_uniq}\n";
	}
}

=comment
# the full table
my @SB=sort(keys %SB_maps);
print $OUT_SPATIAL_TBL "usearch_id;",join(";",@SB),"\n";
foreach my $usearch_uniq (sort keys %usearch_uniq_size_after_filter)
{
	print $OUT_SPATIAL_TBL $usearch_uniq;
	my $QA_Sum=0;
	foreach my $SB (@SB)
	{
		if (exists $SB_maps{$SB}{$usearch_uniq})
		{
			# print $OUT_SPATIAL_TBL ";$usearch_uniq_size_after_filter{$usearch_uniq}";
			$QA_Sum=$QA_Sum+$SB_maps{$SB}{$usearch_uniq};
			print $OUT_SPATIAL_TBL ";$SB_maps{$SB}{$usearch_uniq}";
		}
		else
		{
			print $OUT_SPATIAL_TBL ";0";
		}
	}
	print $OUT_SPATIAL_TBL "\n";
	if ($QA_Sum != $usearch_uniq_size_after_filter{$usearch_uniq})
	{
		print STDERR "[ERROR] usearch_uniq_size_after_filter{$usearch_uniq}: $QA_Sum != $usearch_uniq_size_after_filter{$usearch_uniq}\n";
	}
}
=cut 
close ($OUT_SPATIAL_TBL);