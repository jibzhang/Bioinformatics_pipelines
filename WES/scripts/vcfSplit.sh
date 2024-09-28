#!/bin/sh
# This script segregates Single Nucleotide Variants (SNVs) and insertions/deletions (INDELs) from a provided VCF file

# Initialize variables for the options
VCF_FILE=""
SNV_FILE=""
INDEL_FILE=""

# Use getopts to handle named parameters (-v for VCF_FILE, -s for SNV_FILE, -i for INDEL_FILE)
while getopts "v:s:i:" opt; do
  case $opt in
    v) VCF_FILE="$OPTARG"
    ;;
    s) SNV_FILE="$OPTARG"
    ;;
    i) INDEL_FILE="$OPTARG"
    ;;
    \?) printf "Invalid option -$OPTARG" >&2
    ;;
  esac
done

# Check if VCF_FILE, SNV_FILE and INDEL_FILE are provided. If not, print error message and exit script.
if [ -z "$VCF_FILE" ] || [ -z "$SNV_FILE" ] || [ -z "$INDEL_FILE" ]; then
    printf "Please provide VCF, SNV, and INDEL file names with -v, -s, and -i options respectively!\n"
    exit 1
fi

# Use Perl to parse the VCF_FILE.
# Assign the first value of the split of the 5th field to alt_variant.
# If a line begins with '#' or if the length of the 4th and alt_variant fields is equal to 1 (indicating a Single Nucleotide Variant), print the line to SNV_FILE.
perl -lane '
  my ($alt_variant) = split(",", $F[4]);
  print $_ if ( (/^#/) || (length($F[3]) == 1 && length($alt_variant) == 1 && $alt_variant ne "*") )
' "$VCF_FILE" > "$SNV_FILE"

# Similarly, if a line begins with '#' or if the length of the 4th or alt_variant fields is greater than 1 (indicating an insertion or deletion), print the line to INDEL_FILE.
perl -lane '
  my ($alt_variant) = split(",", $F[4]);
  print $_ if ( (/^#/) || (length($F[3]) > 1 || length($alt_variant) > 1) && $alt_variant ne "*" )
' "$VCF_FILE" > "$INDEL_FILE"
