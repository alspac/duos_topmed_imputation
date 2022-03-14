# Script i have used to combine datasets (in this case alspac_old ONCOARRAY and HN5000_ALSPAC_2020)
# The script uses bcftools which you can load by loading the module samtools 
# Chr X is done seperately but again this is not required and could be combined with some small tweeks

# Files
# $combined/$study/chr$i.filtered.vcf.gz - this is the imputed files from each dataset to combine, i have pre-pruned to redcue the file size (but you do not need to do this)
# 'allPass_snpList.txt' - single column of IDs for the SNPs passing the thresholds set, fomat of the IDs will need to match that in the .vcf files (liekly chr:pos:REF:ALT)

# Filter on the AllSNPs list - using allPass_snpList

module load apps/samtools-1.9.1
for i in `seq 1 22`; do
    for study in  alspac_old ONCOARRAY HN5000_ALSPAC_2020 ; do
         bcftools view -i 'ID=@./data/combined_imputed1/allPass_snpList.txt' $combined/$study/chr$i.filtered.vcf.gz -Oz -o $clean/$study/imputed_round1/chr$i.AllSnp.filtered.vcf.gz &
    done
done

# Chr X
for i in X; do
    for study in  alspac_old ONCOARRAY HN5000_ALSPAC_2020 ; do
          bcftools view -i 'ID=@./data/combined_imputed1/allPass_chrX_snpList.txt' $combined/$study/chr$i.filtered.vcf.gz -Oz -o $combined/$study/chr$i.AllSnp.filtered.vcf.gz &
    done
done

# Index these new vcf files

for i in `seq 1 22`; do
    for study in  alspac_old ONCOARRAY HN5000_ALSPAC_2020 ; do
        bcftools index $combined/$study/chr$i.AllSnp.filtered.vcf.gz &
    done
done

# Chr X
for i in X; do
    for study in  alspac_old ONCOARRAY HN5000_ALSPAC_2020 ; do
        bcftools index $combined/$study/chr$i.AllSnp.filtered.vcf.gz &
    done
done


# Merge datasets within each chromosome

for i in `seq 1 22`; do
         bcftools merge --force-samples $combined/alspac_old/chr$i.AllSnp.filtered.vcf.gz $combined/ONCOARRAY/chr$i.AllSnp.filtered.vcf.gz $combined/HN5000_alspac_2020/chr$i.AllSnp.filtered.vcf.gz -Oz -o $combined/chr$i.vcf.gz &
done

# Chr X
for i in X; do
         bcftools merge --force-samples $combined/alspac_old/chr$i.AllSnp.filtered.vcf.gz $combined/ONCOARRAY/chr$i.AllSnp.filtered.vcf.gz $combined/HN5000_ALSPAC_2020/chr$i.AllSnp.filtered.vcf.gz -Oz -o $combined/chr$i.vcf.gz &
done

# Annotate the vcf file to force the label genotype

for i in `seq 1 22`; do
    bcftools annotate -x FORMAT $combined/chr$i.vcf.gz -Oz -o $combined/chr$i.GT.vcf.gz
done

# Chr X
for i in X; do
    bcftools annotate -x FORMAT $combined/chr$i.vcf.gz -Oz -o $combined/chr$i.GT.vcf.gz
done
