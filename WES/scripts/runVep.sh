#!/bin/sh
echo "vep annotation"
input=$1
output=$2
forkN=$3

PATH="/home/zgu_labs/anaconda3/envs/vep/bin:$PATH"

dir_cache=/ref_genomes/genome_anno/human/vep/

if [ -z "$3" ];then
   echo "input: inputFile outputFile forkNumber!"
   exit 0;
fi


perl -lane '$id=join("-",@F[0..3]);print join("\t", (@F[0,1],$id , @F[2..$#F]))' $input > $output.vcf

vep -cache -a GRCh38 --force --per_gene --pick -e --tab -i $output.vcf -o $output.vep.tmp.1 --fork $forkN --dir_cache $dir_cache --stats_file $output.vep.html

grep -Pv "^##" $output.vep.tmp.1 > $output.vep.tmp.2

paste $input $output.vep.tmp.2 > $output.vep.tmp.3
export output

perl -lane 'BEGIN{open IN, $ENV{output}.".vep.tmp.3"; $h=<IN>; @p=split(/\t/, $h); $i=0; map{ $i=$_ if ($p[$_] eq "MAX_AF")} 0..$#p;}print unless ($F[$i] >0.01)' $output.vep.tmp.3 > $output 

rm $output.vcf $output.vep.tmp.1 $output.vep.tmp.2 $output.vep.tmp.3
