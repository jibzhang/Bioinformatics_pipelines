#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

my %opts = ( n => 6, x => 0.9, d => 0.5 );
getopts( 'f:a:c:o:n:x:d:', \%opts );
print(
	qq/Filter case snv with consensus information with control sample.
Options: 
	-f STRING	input reference fasta file [required]
	-a STRING	input case snv file after the first round of filter with control sample [required]
	-c STRING	input control bam file name [required]
	-o STRING	output result file name [required]
	-n INT		min control depth (>= n) [$opts{n}]
	-x DOUBLE	max control het. (<= x) [$opts{x}]
	-d DOUBLE	max ratio'ed control het. (<= d*case_het) [$opts{d}]
	
\n/
);

die("Ref fasta, case, control and output file names are needed!\n")
  if ( !$opts{f} || !$opts{a} || !$opts{c} || !$opts{o} );

my ( $fastaRef, $caseSnv, $controlBam, $fout, $depthTH, $maxControlHet,
	$hetDiffR )
  = ( $opts{f}, $opts{a}, $opts{c}, $opts{o}, $opts{n}, $opts{x}, $opts{d} );

print
"Parameters are:\n\t -f $fastaRef -a $caseSnv -c $controlBam -o $fout -n $depthTH -x $maxControlHet -d $hetDiffR!\n";

die("File $caseSnv dose not exist!\n")
  unless ( -e $caseSnv );

unless ( -e "$caseSnv.cns" ) {
	my $sysReturn = system(
"samtools mpileup -ABQ0 -d10000 -l $caseSnv -f $fastaRef $controlBam  >  $caseSnv.cns"
	);
	die("Error in samtools mpileup!\n") if ( $sysReturn != 0 );
}

open( IN, "$caseSnv.cns" );
my %ctrlInfor = ();
while (<IN>) {
	chomp;
	my ( $chr, $start, $refBase, $depth, $seq, $qseq ) = split(/\t/);
	$ctrlInfor{$chr}{$start} = [ $depth, $seq, $qseq ];
}
close IN;

open( IN,  "$caseSnv" ) || die("Could not open file $caseSnv!\n");
open( OUT, ">$fout" )   || die("Could not create file $fout!\n");
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
		my ( $num, $qs ) = het( $varBase, $ctrlSeqCns, $ctrlQseqCns );
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

sub het {
	my ( $refB, $seq, $qseq ) = @_;
	$seq =~ s/\^.//g;
	$seq =~ s/[\$]//g;
	while ( $seq =~ /[+-](\d+)/ ) {
		$seq =~ s/[+-]\d+\w{$1}//;
	}
	$seq = uc($seq);
	my @seqs  = split( '', $seq );
	my @qseqs = split( '', $qseq );
	my @qs    = ();
	while ( $seqs[0] ) {
		my $b = shift(@seqs);
		my $q = shift(@qseqs);
		push( @qs, ord($q) - 33 ) if ( $b eq $refB );
	}
	@qs = sort { $a <=> $b } (@qs);
	my $n = $#qs + 1;
	$qs[0] = 0 unless ( $qs[0] );
	return ( $n, join( '|', @qs ) );
}
