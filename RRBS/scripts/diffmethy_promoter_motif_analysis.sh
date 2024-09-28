#!/bin/bash

##########################
FILE=$1
name=$(basename "$FILE" | sed 's/.txt//')
hyper_region="/coh_labs/jochan/CD4_TET2KO_RRBS/RnBeads/${name}_hyper.bed"
hypo_region="/coh_labs/jochan/CD4_TET2KO_RRBS/RnBeads/${name}_hypo.bed"
hyper="/coh_labs/jochan/CD4_TET2KO_RRBS/RnBeads/homer/${name}_hyper"
hypo="/coh_labs/jochan/CD4_TET2KO_RRBS/RnBeads/homer/${name}_hypo"

cat $FILE | awk '($15<0.05 && $13>0.25)' > $hyper_region
cat $FILE | awk '($15<0.05 && $13<-0.25)' > $hypo_region
findMotifsGenome.pl $hyper_region hg38 $hyper -size 50
findMotifsGenome.pl $hypo_region hg38 $hypo -size 50