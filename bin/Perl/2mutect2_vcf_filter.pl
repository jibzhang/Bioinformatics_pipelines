#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

## Author: Jibin Zhang <jibzhang@coh.org>

## work similarly as varscan filter

my $ipfi;
my ( $mcov, $mre2, $strf, $mvfr, $opfi, $help ) =
  ( 10, 4, 1, 0.2,, );
GetOptions(
    'input|i=s'       => \$ipfi,
    'min-coverage=i'  => \$mcov,
    'min-reads2=i'    => \$mre2,
    'strand-filter:1' => \$strf,
    'min-var-freq=f'  => \$mvfr,
    'output-file=s'   => \$opfi,
    'help!'           => \$help,
);
$ipfi or &usage;
&usage if $help;

if ($opfi) {    # redirect if output file is assigned
    open OT, ">$opfi";
    select(OT);
}

my $filterheader = 0;
open IN, $ipfi or die($!);
while (<IN>) {
    chomp;
    if (/^#/) {
        if ( $filterheader == 1 && $_ !~ /^##FILTER/ ) {
            print "##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">\n";
            print "##FILTER=<ID=rdp$mcov,Description=\"Quality read depth less than $mcov\">\n";
            print "##FILTER=<ID=vdp$mre2,Description=\"Variant supporting read depth less than $mre2\">\n";
            print "##FILTER=<ID=vaf$mvfr,Description=\"Variant allelic frequency less than $mvfr\">\n";
            $filterheader = 0;
        }
        $filterheader = 1 if /^##filtering_status/;
        print "$_\n";
        next;
    }
    my @a    = split;
    my $pass = 0;
    my ( $chro, $posi, $refa, $muta ) = @a[ 0, 1, 3, 4 ];
    foreach my $info ( @a[ 9 ] ) { $pass += &varfilter($info); }
    if ($pass) {
        print join "\t", @a;
        print "\n";
    }
}
close IN;

sub varfilter {
    my $r = 0;
    my @nums = split /[:,]/, $_[0];
    return $r unless $#nums == 12;
    return $r unless $nums[4] >= $mcov;    # min-coverage
    return $r unless $nums[2] >= $mre2;    # min-read2
    return $r
      unless (
        $strf == 0
        || (   $nums[11] >= 1
            && $nums[12] >= 1 )
      );                                   # strand-filter
    return $r unless $nums[3] >= $mvfr;    # min-var-freq
    $r = 1;
}

sub usage {
    die(
        qq/
USAGE: mutect2_filter.pl -i <*.vcf> [OPTIONS]

	variant file - A VCF file of SNPs or indels produced by VarScan mpileup2snp (v2.3.6 or maybe above)
	               There should be only one sample processed during the previous mpileup step

OPTIONS:
        --min-coverage  Minimum read depth at a position to make a call [10]
        --min-reads2    Minimum supporting reads at a position to call variants [4]
        --strand-filter Ignore variants with >90% support on one strand [1]
        --min-var-freq  Minimum variant allele frequency threshold [0.20]
        --output-file   File to contain variants passing filters [STDOUT]
        --help\/-h       Print this page

CONTACT: Jibin Zhang <jibzhang\@coh.org>

\n/
    );
}
