#!/bin/sh

name=$1
input=$2
output=$3

echo "freq	$name"	> $output
	cat $input |\
		grep -vE 'chrX|chrY|chrM' |\
		awk '$5+$6>9{print int($4)}' |\
		perl -e '
			my %fr;
			my $n;
			for(my $i=0;$i<=100;$i++){
				$fr{$i}=0
			}
			while(<>){
				chomp;
				$n+=1;
				$fr{$_}+=1;
			}
			for(my $i=0;$i<=100;$i++){
				print join "\t",$i,$fr{$i}/$n;
				print "\n";
			}' >> $output
