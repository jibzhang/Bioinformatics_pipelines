#!/usr/bin/perl
use strict;
use warnings;

my $pfx = shift or die("Usage: $0 <prefix>\n");

my $file1 = "$pfx.blaclist_filtered.bed";
my $file2 = "$pfx.promoter5k.bed";
my $file3 = "$pfx.genebody2500plus.bed";
my $file4 = "$pfx.upstream50k.bed";
my $file5 = "$pfx.downstream50k.bed";

my %peakInfo;
my %anno;
open IN, $file1 or die($!);
while (<IN>) {
	chomp;
	my @a = split;
	$peakInfo{$a[3]} = join "\t", @a;
}
close IN;

open IN, $file2 or die($!);
while (<IN>) {
	chomp;
	my @a = split;
	next if $a[13] eq ".";
	$anno{$a[3]} = join "\t", $a[13], "Promoter";
}
close IN;

open IN, $file3 or die($!);
while (<IN>) {
	chomp;
	my @a = split;
	next if $a[13] eq ".";
	next if exists $anno{$a[3]};
	$anno{$a[3]} = join "\t", $a[13], "GeneBody";
}
close IN;

open IN, $file4 or die($!);
while (<IN>) {
	chomp;
	my @a = split;
	next if $a[13] eq ".";
	next if exists $anno{$a[3]};
	$anno{$a[3]} = join "\t", $a[13], "Upstream";
}
close IN;

open IN, $file5 or die($!);
while (<IN>) {
	chomp;
	my @a = split;
	next if $a[13] eq ".";
	next if exists $anno{$a[3]};
	$anno{$a[3]} = join "\t", $a[13], "Downstream";
}
close IN;

print join "\t", qw/chr start end peakName score . foldEnrichment log10P log10Q point_source gene funcRegion/;
print "\n";
foreach my $k (sort keys %peakInfo) {
	$anno{$k} = ".\tIntergenic" unless exists $anno{$k};
	print join "\t", $peakInfo{$k}, $anno{$k};
	print "\n";
}
