#!/usr/bin/perl
use strict;
use warnings;

my $hg38map = "/home/jibzhang/reference/hg38_region/ensembl_hg38.map";
my $file = shift or die("Usage: $0 <*_annotatr_annotation.txt>\n");

my %tss;
my %megagenes;
open IN, $hg38map or die($!);
while (<IN>) {
    chomp;
    my @a = split;
    my $t = 0;
    if ( $a[2] eq "+" ) {
        $t = $a[3];
    } else {
        $t = $a[4];
    }
    my $k = int( $t / 10**6 );
    $k = "$a[1]:$k";
    $megagenes{$k} .= "$a[0],";    # indexing genes by Mb
    $tss{ $a[0] } = "$a[1]:$t";
}
close IN;

open IN, $file or die($!);
while (<IN>) {
    chomp;
    my @a = split;
    if ( $a[0] eq "seqnames" ) {
        print join "\t", @a;
        print "\tNEARBY_GENE\n";
        next;
    }
    print join "\t", @a;
    my $ctr   = int( ( $a[2] + $a[3] ) / 2 );
    my $k1    = int( $ctr / 10**6 );
    my $phase = $ctr % 10**6;
    my $k2    = $k1 + 1;    # refer the current Mb and an adjcent Mb for closest Gene TSS;
    $k2 = $k1 - 1 if $phase < 5 * 10**5;
    $k1 = "$a[1]:$k1";
    $k2 = "$a[1]:$k2";
    my $ng = "NA";
    my @mg;
    foreach my $k ( $k1, $k2 ) {
        if ( exists $megagenes{$k} ) {
            $megagenes{$k} =~ s/,$//;
            @mg = ( @mg, ( split /,/, $megagenes{$k} ) );
        }
    }
    my $MinDist = 10**6;
    foreach my $g (@mg) {
        my ( $chrom, $pos ) = split /:/, $tss{$g};
        my $dist = abs( $pos - $ctr );
        if ( $dist < $MinDist ) {
            $MinDist = $dist;
            $ng      = $g;
        }
    }
    print "\t$ng\n";
}
close IN;

