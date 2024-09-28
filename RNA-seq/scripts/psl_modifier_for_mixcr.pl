#!/usr/bin/perl
use strict;
use warnings;

my $file1 = shift
  or die(
"Usage:$0 <*.clones.txt> <blat_output.psl (against hg19/hg38)> <chain(IGH|IGL|IGK)>\n"
  );
my $file2 = shift
  or die(
"Usage:$0 <*.clones.txt> <blat_output.psl (against hg19/hg38)> <chain(IGH|IGL|IGK)>\n"
  );
my $chain = shift
  or die(
"Usage:$0 <*.clones.txt> <blat_output.psl (against hg19/hg38)> <chain(IGH|IGL|IGK)>\n"
  );

my %black_list;
open IN, $file2 or die($!);
while (<IN>) {
    chomp;
    next unless $_;
    my @a = split;
    next unless $a[1] =~ /^\d+$/;
    my $outtarget = 1;
    if ( $chain eq "IGH" ) {
        $outtarget = 0
          if $a[13] eq "chr14"
              and $a[15] >= 105586437     # IGH locus start point on hg38
              and $a[16] <= 107288051;    # IGH locus end point on hg19
    }
    elsif ( $chain eq "IGL" ) {
        $outtarget = 0
          if $a[13] eq "chr22"
              and $a[15] >= 22026076     # IGL locus start position on hg19
              and $a[16] <= 23265085;    # IGL locus end position on hg38
    }
    elsif ( $chain eq "IGK" ) {
        $outtarget = 0
          if $a[13] eq "chr2"
              and $a[15] >= 88857361      # IGK locus start position on hg38
              and $a[16] <= 90274235;     # IGK locus end position on hg19
    }
    else {
        die(
"Unrecognizable Gene name, which should be one of IGH, IGL or IGK\n"
        );
    }
    my $matched_size       = $a[12] - $a[11];
    my $matched_proportion = $matched_size / $a[10];
    my $match_score        = $a[1] / $matched_proportion;
    $black_list{ $a[9] } = 1
      if ( $matched_proportion >= 0.9 and $match_score >= 0.9 )
      or $outtarget;                      # Filtering conditions
}
close IN;

my $n        = 0;
my $total_rc = 0;
my @cdr3s;
open IN, $file1 or die($!);
while (<IN>) {
    chomp;
    my @a = split /\t/;
    if ( $a[1] =~ /^\d+$/ ) {
        $n++;
        next if exists $black_list{$n};
        $total_rc += $a[1];
        push @cdr3s, $_;
    }
    else {
        print join "\t", @a;
        print "\n";
    }
}
close IN;

foreach my $cdr3 (@cdr3s) {
    my @a = split /\t/, $cdr3;
    $a[2] = $a[1] / $total_rc;
    print join "\t", @a;
    print "\n";
}

