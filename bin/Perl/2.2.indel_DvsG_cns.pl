#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opts = ( n => 6, x => 0.125, d => 0.5 );
getopts( 'f:a:c:o:n:x:d:', \%opts );
print(
	qq/
Description:
	Filter case indel with consensus information with control sample.
Options: 
	-f STRING	input reference fasta file [required]
	-a STRING	input case snv file after the first round of filter with control sample [required]
	-c STRING	input control bam file name [required]
	-o STRING	output result file name [required]
	-n INT		control depth >= n [$opts{n}]
	-x DOUBLE	control het. <= x [$opts{x}]
	-d DOUBLE	control het. <= d * case het. [$opts{d}]
	
\n/
);

die("Ref fasta, case, control and output file names are needed!\n")
  if ( !$opts{f} || !$opts{a} || !$opts{c} || !$opts{o} );

my ( $fastaRef, $caseIndel, $controlBam, $fout, $depthTH, $maxControlHet,
	$hetDiffR )
  = ( $opts{f}, $opts{a}, $opts{c}, $opts{o}, $opts{n}, $opts{x}, $opts{d} );

print
"Parameters are:\n\t -f $fastaRef -a $caseIndel -c $controlBam -o $fout -n $depthTH -x $maxControlHet -d $hetDiffR!\n";

die("File $caseIndel dose not exist!\n")
  unless ( -e $caseIndel );

unless ( -e "$caseIndel.cns" ) {
	my $sysReturn = system(
"samtools mpileup -ABQ0 -d10000 -l $caseIndel -f $fastaRef $controlBam  >  $caseIndel.cns"
	);
	die("Error in samtools mpileup!\n") if ( $sysReturn != 0 );
}

open( IN, "$caseIndel.cns" );
my %ctrlInfor = ();
while (<IN>) {
	chomp;
	my ( $chr, $start, $refBase, $depth, $seq, $qseq ) = split(/\t/);
	$ctrlInfor{$chr}{$start} = [ $depth, $seq, $qseq ];
}
close IN;

open( IN, "$caseIndel" ) || die("Could not open file $caseIndel!\n");
open( OUT, ">$fout" ) || die("Could not create file $fout!\n");
while (<IN>) {
	chomp;
	my (
		$chr,        $start,     $refBase,        $varBase,
		$end,        $gene,      $class,          $type,
	        $qual,       $depth,     $refDepth,       $varDepth,
		$mapQ,       $tumorHet,  $ECNT,	          $MQbase,
		$ReadPos,    $popfreq,   $ctrlDepth,      $ctrlVarDepth,
		$ctrlHet
	) = split(/\t/);
	if ( $ctrlDepth ne "NA" ) {
		print OUT "$_\n";
	}
	elsif ( exists $ctrlInfor{$chr}{$start} ) {
		my ( $ctrlDepthCns, $ctrlSeqCns, $ctrlQseqCns ) =
		  @{ $ctrlInfor{$chr}{$start} };
		next if ( $ctrlDepthCns < $depthTH );
		my $num        = $ctrlSeqCns =~ tr/-+*/-+*/;
		my $controlHet = $num / $ctrlDepthCns;
		next if ( $controlHet > $maxControlHet );
		next if ( $num >= 2 && $controlHet > $hetDiffR * $tumorHet );
		print OUT join(
			"\t",
			(
				$chr, $start, $refBase, $varBase,
				$end, $gene,  $class,   $type,
				$qual,$depth, $refDepth,$varDepth,
				$mapQ, $tumorHet, $ECNT,$MQbase,
				$ReadPos,$popfreq, $ctrlDepthCns,$num, 
				$controlHet
			)
		  ),
		  "\n";
	}
}
close OUT;

