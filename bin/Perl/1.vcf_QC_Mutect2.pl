#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opts = (n => 3, D => 10, h => 0.01 );
getopts( 'i:o:n:D:h:', \%opts );
print(
	qq/qulity control of vcf files
Options: 
	-i STRING	raw vcf from GATK [required]
	-o STRING	output file [required]
	-n INT		minimum number of reads support the var [$opts{n}]
	-D INT		minimum number of total depth [$opts{D}]
	-h DOUBLE	minimum heterozygous rate of supportive reads [$opts{h}]
\n/
);

die("Input vcf and output file names are needed!\n")
  if ( !$opts{i} || !$opts{o} );
my ( $fin, $fout, $ad2TH, $DpTH,$hetTH ) =
  ( $opts{i}, $opts{o}, $opts{n}, $opts{D}, $opts{h} );
print "Parameters are:\n\t -i $fin -o $fout -n $ad2TH -D $DpTH -h $hetTH!\n";

open( IN,  "$fin" )   || die("Could not open file $fin!\n");
open( OUT, ">$fout" ) || die("Could not create file $fout!\n");
my $title = join( "\t",
	qw(chrom. start refAllele varAllele end quality depth refDepth varDepth mapQual MAF ECNT MBQ MPOS POPAF)
);
print OUT '#' . $title, "\n";

while ( my $line = <IN> ) {
	next if ( $line =~ /^#/ );
	chomp $line;
	my ( $chr, $start, $ref, $alt, $qual, $dp, $ad1, $ad2, $mq, $het,
		$ECNT, $MBQ, $MPOS, $POPAF)
	  =
	  annoExtract($line);
	next unless ($het);
	next if ( length($chr) > 5 or $ad2 < $ad2TH or $dp < $DpTH or $het < $hetTH );
	my $end = $start + length($ref) - 1;
	print OUT join( "\t",
		( $chr, $start, $ref, $alt, $end, $qual, $dp, $ad1, $ad2, $mq, $het,    
		  $ECNT, $MBQ, $MPOS, $POPAF) ),
	  "\n";
}
close IN;
close OUT;

sub annoExtract {
	my ($line) = @_;
	my (
		$chr,    $start, $ID,     $ref,
		$alt,    $tmpQual,  $FILTER, $info,
		$FORMAT, @alleleInfoArray
	) = split( /\t/, $line );
	($alt) = split( ',', $alt );
	my ($mq, $qual) = $info =~ /MMQ=\d+,(\d+).*;TLOD=(.+)$/;
	$qual=~s/,.+//;
	my @formats = split( ':', $FORMAT );
	my ($adIdx) = grep { $formats[$_] eq "AD" } 0 .. $#formats;
	my ($dpIdx) = grep { $formats[$_] eq "DP" } 0 .. $#formats;
	my ( $ad1, $ad2 ) = ( 0, 0 );
	my $dp;
	for my $alleleInfo (@alleleInfoArray) {
		next unless ( $alleleInfo =~ /:/ );
		my @alleles = split( ':', $alleleInfo );
		$dp = $alleles[$dpIdx];
		my ( $ad1_i, $ad2_i ) = split( ',', $alleles[$adIdx] );
		if ( $ad1_i =~ /\d+/ && $ad2_i =~ /\d+/ ) {
			$ad1 = $ad1 + $ad1_i;
			$ad2 = $ad2 + $ad2_i;
		}
	}
	my $het = sprintf( "%0.4f", $ad2 / $dp ) if ($ad2);

	my ($ECNT, $MBQ, $MPOS, $POPAF)=('NA')x4;
	if($info =~ /ECNT=(\d+)/){
		$ECNT = $1;
	}
	if($info =~ /MBQ=\d+,(\d+)/){
		$MBQ = $1;
	}
	if($info =~ /MPOS=(\d+)/){
		$MPOS = $1;
	}
	if($info =~ /POPAF=(.+?)[,;]/){
		$POPAF = $1;
	}
	return ( $chr, $start, $ref, $alt, $qual, $dp, $ad1, $ad2, $mq, $het, 
	$ECNT, $MBQ, $MPOS, $POPAF);
}

