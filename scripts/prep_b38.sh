#!/bin/bash


#########################
### ALSPAC IMPUTATION ###
#########################

# DATE STARTED: 20/12/2021

## TASK 1: PRE-IMPUATION
## The raw bed/bim/fam files are located in scratch

## TOPMED VCF: https://bravo.sph.umich.edu/freeze5/hg38/download
## McCarthy Tools: https://www.well.ox.ac.uk/~wrayner/tools/
## TOPMED guidance: https://topmedimpute.readthedocs.io/en/latest/prepare-your-data.html
## Tom's script: https://github.com/tomdudding/0033_HNC_GWAS/commit/df46a7d36f5563f61bffad563c524b420d31b415

module load apps/plink/1.90 # --recode to ped/map does not work on plink 2.00
module load languages/perl/5.30.0
module load languages/anaconda3/2020-3.8.5


# It's good practice to avoid hard coded paths in your scripts
# So store the paths in a config file and then read them in
# Use jq to get paths from config.json
# Download jq from here to a folder called ~/bin/ in your home directory https://stedolan.github.io/jq/
# Then add ~/bin to your PATH in ~/.bash_profile 

# Get paths
duo_geno=`cat ../config.json | jq -r '.duo_geno'`
scratch=`cat ../config.json | jq -r '.scratch_dir'`

# Make scratch directory if necessary
mkdir -p ${scratch}

# Create map file from bim file
# This avoids needing to convert bed to ped
cut -f 1-4 ${duo_geno}.bim > ${scratch}/duos_b37.bim.map

# Liftover 
./liftOverPlink.py \
-m ${scratch}/duos_b37.bim.map \
-c ./hg19ToHg38.over_Xto23.chain \
-o ${scratch}/duos_b38 \
-e ./liftOver

# Have any SNPs been dropped?
wc -l ${scratch}/duos_b38.map
wc -l ${scratch}/duos_b37.bim.map

# Remove rsids that did not liftover
# We'll call this new file duos_b38 but it's currently in b37 - we're gonna update the coordinates next
cut -f 4 ${scratch}/duos_b38.bed.unlifted | sed "/^#/d" > ${scratch}/excluded_snps.dat
plink --bfile ${duo_geno} --exclude ${scratch}/excluded_snps.dat --make-bed --out ${scratch}/duos_b38


# Check that all the liftover rsids are the same as the original rsids
wc -l ${scratch}/duos_b38.bim
wc -l ${scratch}/duos_b38.map

cut -f 2 ${scratch}/duos_b38.bim > ${scratch}/rsid1
cut -f 2 ${scratch}/duos_b38.map > ${scratch}/rsid2

DIFF=$(diff ${scratch}/rsid1 ${scratch}/rsid2) 
if [ "$DIFF" != "" ] 
then
    echo "ERROR - RSIDS DON'T MATCH"
fi

# Can now reconstruct bim file using new chromosomes and coordinates
# Replace the old duos_b38.bim with the newly constructed one
paste ${scratch}/duos_b38.bim ${scratch}/duos_b38.map | head
paste ${scratch}/duos_b38.map <(cut -f 5-6 ${scratch}/duos_b38.bim) > temp && mv temp ${scratch}/duos_b38.bim

# Identify and remove "bad" SNPs
./rmBadLifts.py \
--map ${scratch}/duos_b38.map \
--out ${scratch}/duos_b38_good.map \
--log ${scratch}/duos_bad_lifted.dat
cut -f 2 ${scratch}/duos_bad_lifted.dat > ${scratch}/bad_snps.dat
plink --bfile ${scratch}/duos_b38 --exclude ${scratch}/bad_snps.dat --make-bed --out ${scratch}/duos_b38 --allow-extra-chr


# Get frequency using only the mothers
grep "M" ${scratch}/duos_b38.fam | cut -d " " -f 1-2 > ${scratch}/mums.txt
plink --bfile ${scratch}/duos_b38 --keep ${scratch}/mums.txt --freq --out ${scratch}/freq


## Note from Gib - I've stopped at this point as I don't have the ancilliary files below

######################

# Converting TOPMED VCF to HRC formatted reference legend
perl /mnt/storage/scratch/ds21941/ALSPAC/Scripts/CreateTOPMed.pl -i /mnt/storage/scratch/ds21941/ALSPAC/Scripts/Reference/ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz

# Running checks on bim files
perl /mnt/storage/scratch/ds21941/ALSPAC/Scripts/HRC-1000G-check-bim.pl -b /mnt/storage/scratch/ds21941/ALSPAC/Data/gwa_660_g0m/data.bim -f /mnt/storage/scratch/ds21941/ALSPAC/Data/gwa_660_g0m/freq.afreq -r /mnt/storage/scratch/ds21941/ALSPAC/Scripts/Reference/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz -h -o /mnt/storage/scratch/ds21941/ALSPAC/bimChecks/gwa_660_g0m/

# Creating files in plink
plink /mnt/storage/scratch/ds21941/ALSPAC/bimChecks/gwa_660_g0m/Run-plink.sh
