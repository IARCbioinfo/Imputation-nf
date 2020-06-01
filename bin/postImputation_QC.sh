#!/bin/bash
#SBATCH -A "MR_Signatures"
#SBATCH -c 1
#SBATCH --mem=40GB


chr=$1
refs_path="/data/gep/MR_Signatures/work/gabriela/imputation/January_2020/ref_data/"
pop=$2


if [ $pop = "ALL" ]
then

  # Generate the reference frequency file
  # Generate a tab-delimited header
  echo -e 'CHR\tSNP\tREF\tALT\tAF\tEUR_AF\tAFR_AF\tEAS_AF' > 1000GP_chr${chr}_imputation_all.frq
  # Query the required fields from the VCF file and append to the allele frequency file
  bcftools query -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\t%INFO/EUR_AF\t%INFO/AFR_AF\t%INFO/EAS_AF\n' ${refs_path}1000GP_chr${chr}.bcf >> 1000GP_chr${chr}_imputation_all.frq

  # Generate a header for the output file
  echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAF_GROUP' > INFO_group_chr${chr}.txt

  # Query the required fields and add frequency group (1, 2 or 3) as the last column  chr_${chr}_combined.vcf.gz
  bcftools +fill-tags chr_${chr}_combined.vcf.gz -- -t AF,AN,AC | bcftools  query -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\t%INFO/R2\t-\n' | \
  # Here $5 refers to AF values, $7 refers to AF group
  awk -v OFS="\t" \
      '{if ($5>=0.05 && $5<=0.95) $7=1; \
          else if(($5>=0.005 && $5<0.05) || \
          ($5<=0.995 && $5>0.95)) $7=2; else $7=3} \
          { print $1, $2, $3, $4, $5, $6, $7 }' \
      >> INFO_group_chr${chr}.txt

fi

if [ $pop != "ALL" ]
then
  #plink_path="/data/gep/MR_Signatures/work/gabriela/imputation/January_2020/plink_files/Genotyping_R2/combined/"
  paste -d "_" <(awk '{print $1}' out_pop_admixture/1000G_checking_${pop}/samples_${pop}.txt) <(awk '{print $2}' out_pop_admixture/1000G_checking_${pop}/samples_${pop}.txt)  > samples_${pop}_chr${chr}.txt

  # Generate a header for the output file
  echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAF_GROUP' > INFO_group_chr${chr}_${pop}.txt

  bcftools view -S samples_${pop}_chr${chr}.txt chr_${chr}_combined.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\t%INFO/R2\t-\n' | \
  # Here $5 refers to AF values, $7 refers to AF group
  awk -v OFS="\t" \
      '{if ($5>=0.05 && $5<=0.95) $7=1; \
          else if(($5>=0.005 && $5<0.05) || \
          ($5<=0.995 && $5>0.95)) $7=2; else $7=3} \
          { print $1, $2, $3, $4, $5, $6, $7 }' \
      >> INFO_group_chr${chr}_${pop}.txt

  # # Generate a header for the output file
  echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAF_GROUP' > tmp.txt
  cat INFO_group_chr${chr}_${pop}.txt >> tmp.txt
  mv tmp.txt INFO_group_chr${chr}_${pop}.txt
fi