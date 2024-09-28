#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;

my %opts = ( b=> "/ref_genomes/dbSNP/human/153_hg38/dbSnp153.bb");
getopts( 'i:o:b:', \%opts );
print(
	qq/Filter common indels in dbsnp database. "bigBedToBed" tool is required.
	In indel file, the fields should be chr, poz, ...
Options: 
	-i STRING	input indel file name [required]
	-o STRING	output result file name [required]
	-b STRING	dbSNP bigBed file name [required]
\n/
);

die("input, output and dbSNP bigBed file names are needed!\n")
  if ( !$opts{i} || !$opts{o} || !$opts{b} );

my ( $fin, $fout, $fbb ) = ( $opts{i}, $opts{o}, $opts{b} );
print("Parameters are:\n\t -i $fin -o $fout -b $fbb!\n");

open( IN, "$fin" ) || die("Could not open file $fin!\n");
my $header = <IN>;
chomp $header;
open( OUT, ">$fout" ) || die("Could not create file $fout!\n");
print OUT $header,
  "\tRS_ID\tdbSNP_ref\tdbSNP_alt\taltCount\tfreqSourceCount\tclass\tucscNotes\n";
while (<IN>) {
	chomp;
	my ( $rs, $dbSNP_ref, $altCount, $dbSNP_alt, $freqSourceCount, $class, $ucscNotes ) = (0)x7;
	my ( $chr, $poz, $ref, $alt, $end ) = split;
	$chr = "chr" . $chr if ( $chr !~ /chr/ );
	my $poz_0 = $poz - 1;
	my $dbSNP_result =
	  `bigBedToBed -chrom=$chr -start=$poz_0 -end=$end $fbb stdout`;
	if ($dbSNP_result) {
		my @dbSNP_records = split( /\n/, $dbSNP_result );
		for (@dbSNP_records) {
			if (/\tdel|ins\t/) {
				my @parts = split( /\t/, $_ );
				( $rs, $dbSNP_ref, $altCount, $dbSNP_alt, $freqSourceCount, $class, $ucscNotes ) = @parts[ 3 .. 6, 8, 13, 14 ];
			}
		}
	}
	next if($ucscNotes =~/commonAll/);
	print OUT join("\t", ($_,  $rs, $dbSNP_ref, $dbSNP_alt, $altCount, $freqSourceCount, $class, $ucscNotes)), "\n";
}
close IN;
close OUT;
