use strict;
use IO::Zlib;
use IO::Compress::Gzip qw(gzip $GzipError);
# @ARGV=("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI6//191028-062_A1_S1_Ar_multimap_bacteria_R2.usearch_unique_vs_NT_Jan2021.mmseq2.m8","YES");

my $MMSeq2_file=shift;
my $multiple_best=shift;

if (!defined $multiple_best) {$multiple_best="NO";}
else {$multiple_best=uc($multiple_best);}

my $MMSeq2_best_hit=$MMSeq2_file.".best_hit";
if ($multiple_best eq "YES") {$MMSeq2_best_hit=$MMSeq2_best_hit."_multiple"}

my $bit_score_column=11;
my $q_start=6;
my $q_end=7;
my $q_length=12;

# The file is formatted as a tab-separated list with 12 columns: (1,2) identifiers for query and target sequences/profiles, 
# (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, 
# (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score (13) query length.
# --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen"

my $IN;
my $gzip_in="NO";
if ($MMSeq2_file=~/(.*?).gz$/)
{
	my $base=$1;
	$IN = new IO::Zlib;
	$IN->open($MMSeq2_file, "rb") || die "Can't open IN '$MMSeq2_file' $!";
	if ($multiple_best eq "YES") 
	{
		$MMSeq2_best_hit=$base.".best_hit_multiple.gz"
	} else 
	{
		$MMSeq2_best_hit=$base.".best_hit.gz";
	}
	$gzip_in="YES";
}
else
{
	open ($IN,"<",$MMSeq2_file) || die "Can't open IN '$MMSeq2_file' $!";
}
my $prev_q="NA";
my %hits=();

my $OUT;
if ($gzip_in eq "YES")
{
	$OUT = IO::Compress::Gzip->new($MMSeq2_best_hit) or die "gzip failed to open '$MMSeq2_best_hit': $GzipError\n";
}
else
{
	open ($OUT,">",$MMSeq2_best_hit) || die "Can't open OUT '$MMSeq2_best_hit' $!";
}

while (my $line=<$IN>)
{
	chomp ($line);
	my @line=split(/\t/,$line);
	my $q_id=$line[0];
	if ($prev_q eq "NA") {$prev_q=$q_id}
	if  ($prev_q eq $q_id) # same q -> select best hit
	{
		if (exists $hits{$q_id})
		{
			if ($hits{$q_id}[0][$bit_score_column]<$line[$bit_score_column])
			{
				my @hits=(\@line);
				$hits{$q_id}=[@hits];
			}
			elsif ($hits{$q_id}[0][$bit_score_column]==$line[$bit_score_column] and $multiple_best eq "YES") # the same, we take all the hits
			{
				my @hits=(@{$hits{$q_id}});
				push (@hits,\@line);
				$hits{$q_id}=[@hits];
			}
		}
		else
		{
			my @hits=(\@line);
			$hits{$q_id}=[@hits];
		}
		$prev_q=$q_id;
	}
	else # new q - > print the best of previous one and start aggregate
	{
		# print the best of previous one
		foreach my $hit_ArrRef (@{$hits{$prev_q}})
		{
			my $Q_cover=($hit_ArrRef->[$q_end]-$hit_ArrRef->[$q_start]+1)/$hit_ArrRef->[$q_length];
			$Q_cover=sprintf("%.3f",$Q_cover);
			print $OUT join("\t",(@{$hit_ArrRef},$Q_cover)),"\n";
		}
		# update the new one
		%hits=();
		my @hits=(\@line);
		$hits{$q_id}=[@hits];
		$prev_q=$q_id;
	}
	
	
}
if ($gzip_in eq "NO") {
	close ($IN);
}
else
{
	$IN->close;	
}
# print the last q hits
foreach my $hit_ArrRef (@{$hits{$prev_q}})
{
	my $Q_cover=($hit_ArrRef->[$q_end]-$hit_ArrRef->[$q_start]+1)/$hit_ArrRef->[$q_length];
	$Q_cover=sprintf("%.3f",$Q_cover);
	print $OUT join("\t",(@{$hit_ArrRef},$Q_cover)),"\n";
}
if ($gzip_in eq "NO") {
	close ($OUT);
}
else
{
	$OUT->close;	
}
# foreach my $q_id (sort {$hits{$b}[0][$bit_score_column]<=>$hits{$a}[0][$bit_score_column]} keys %hits)
# {
# 	foreach my $hit_ArrRef (@{$hits{$q_id}})
#	{
#		my $Q_cover=($hit_ArrRef->[$q_end]-$hit_ArrRef->[$q_start]+1)/$hit_ArrRef->[$q_length];
#		$Q_cover=sprintf("%.3f",$Q_cover);
#		print $OUT join("\t",(@{$hit_ArrRef},$Q_cover)),"\n";
#	}
#}
