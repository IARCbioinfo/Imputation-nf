#! /usr/bin/env nextflow

// Copyright (C) 2020 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  <PROGRAM_NAME> <VERSION>: <SHORT DESCRIPTION>         "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE : nextflow run IARCbioinfo/bin/Preparation.nf"
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run IARCbioinfo/bin/Preparation.nf"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--<OPTION>                      <TYPE>                      <DESCRIPTION>"
    log.info "--out                      directory                       Directory to use as environment"
    log.info ""
    log.info "Optional arguments:"
    log.info "--<OPTION>                      <TYPE>                      <DESCRIPTION>"
    log.info ""
    log.info "Flags:"
    log.info "--<FLAG>                                                    <DESCRIPTION>"
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}

params.out = '.'

process topmed_prep{
  publishDir params.out+'files/', mode: 'copy'

  output:
  file 'PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz' into legend_file
  file 'bravo-dbsnp-all.vcf.gz*' into vcf_file

  shell:
  '''
  curl 'https://bravo.sph.umich.edu/freeze5/hg38/download/all' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: _ga=GA1.2.1477458196.1624480500; _gid=GA1.2.800547376.1624480500; remember_token="gabriel.ag.aurelie@gmail.com|460bee66c4d92d96634bd5c507a5929490124b4d13c958db91aaabf91f91720372c8fbef1e59a85f865655f8df7008dc92c0376db6372d6b6d657aa086534fee"; _gat_gtag_UA_73910830_2=1' --compressed > bravo-dbsnp-all.vcf.gz
  curl -O https://www.well.ox.ac.uk/~wrayner/tools/CreateTOPMed.zip
  unzip CreateTOPMed.zip
  ./CreateTOPMed.pl -i bravo-dbsnp-all.vcf.gz -o PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz
  bcftools index bravo-dbsnp-all.vcf.gz
  '''
}

/*

process Download{
  publishDir params.out+'files/ref/vcf/', mode: 'copy'

  output:
  file '*recalibrated_variants*' into VCF_all
  file '*recalibrated_variants*' into VCF_X
  file '*recalibrated_variants*' into VCF_Y

  file '*recalibrated_variants*' into VCF_Legend_all
  file '*recalibrated_variants*' into VCF_Legend_X
  file '*recalibrated_variants*' into VCF_Legend_Y

  file '*GRCh38*' into fasta_files

  shell:
  """
  wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gzip -d > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
  samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
  #wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr{{1..22},X,Y}.recalibrated_variants.vcf.gz{,.tbi}
  curl -O "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr{{1..22},X,Y}.recalibrated_variants.vcf.gz{,.tbi}"
  """}

process BCF_ALL{
  publishDir params.out+'files/ref/bcf/', mode: 'copy'

  input:
  file data from VCF_all.collect()
  file fasta from fasta_files.collect()
  val chromosome from 1..22

  output:
  file '*bcf*' into Bcf_ALL

  shell:
  """
  chr=!{chromosome}
  (bcftools view --no-version -h 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz | grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\\([0-9XY]\\)/##contig=<ID=chr\\1/'; bcftools view --no-version -H -c 2 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz | sed 's/^/chr/') | bcftools norm --no-version -Ou -m -any | bcftools norm --no-version -Ob -o 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.bcf -d none -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && bcftools index -f 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.bcf
  """}


process Legend_ALL{
  input:
  file data from VCF_Legend_all.collect()
  val chromosome from 1..22

  output:
  file '*legend' into Legend_ALL

  shell:
  '''
  CHR=!{chromosome}
  bcftools query -f '%ID:%POS:%REF:%ALT\\ %CHROM\\ %POS\\ %REF\\ %ALT\\ %INFO/AFR_AF\\ %INFO/AMR_AF\\ %INFO/EAS_AF\\ %INFO/EUR_AF\\ %INFO/SAS_AF\\ %INFO/AF\\n' 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz > chr\${CHR}.legend
  '''}
process Legend_X{
  input:
  file data from VCF_Legend_X.collect()
  val chromosome from "X"

  output:
  file '*legend' into Legend_X

  shell:
  '''
  CHR=X
  bcftools query -f '%ID:%POS:%REF:%ALT\\ %CHROM\\ %POS\\ %REF\\ %ALT\\ %INFO/AFR_AF\\ %INFO/AMR_AF\\ %INFO/EAS_AF\\ %INFO/EUR_AF\\ %INFO/SAS_AF\\ %INFO/AF\\n' 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz > chr\${CHR}.legend
  sed -i 's/X/23/g' chr\${CHR}.legend
  '''}
process Legend_Y{
  input:
  file data from VCF_Legend_Y.collect()
  val chromosome from "Y"

  output:
  file '*legend' into Legend_Y

  shell:
  '''
  CHR=Y
  bcftools query -f '%ID:%POS:%REF:%ALT\\ %CHROM\\ %POS\\ %REF\\ %ALT\\ %INFO/AFR_AF\\ %INFO/AMR_AF\\ %INFO/EAS_AF\\ %INFO/EUR_AF\\ %INFO/SAS_AF\\ %INFO/AF\\n' 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz > chr\${CHR}.legend
  sed -i 's/Y/24/g' chr\${CHR}.legend
  '''}
process Legend_Concatenation{
  publishDir params.out+'files/', mode: 'copy'

  input:
  file data from Legend_ALL.collect()
  file data from Legend_X.collect()
  file data from Legend_Y.collect()

  output:
  file 'ALL.chr_GRCh38.genotypes.20170504.legend' into Concat_Legend

  shell:
  """
  echo id chr position a0 a1 AFR AMR EAS EUR SAS ALL > header.txt
  cat header.txt chr1.legend chr2.legend chr3.legend chr4.legend chr5.legend chr6.legend chr7.legend chr8.legend chr9.legend chr10.legend chr11.legend chr12.legend chr13.legend chr14.legend chr15.legend chr16.legend chr17.legend chr18.legend chr19.legend chr20.legend chr21.legend chr22.legend chrX.legend chrY.legend > 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05.legend

  """}
process M3VCF_ALL{
  validExitStatus 0,134,143,255
  publishDir params.out+'files/ref/m3vcf/', mode: 'copy'
  cpus=16

  input:
  file data from Bcf_ALL.collect()
  val chromosome from 1..22

  output:
  file '*m3vcf.gz' into M3VCF_All

  shell:
  """
  CHR=!{chromosome}
  bcftools index --csi --threads 16 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.bcf -o ALL.chr\${CHR}_GRCh38.bcf.csi
  bcftools view --threads 16 -g ^miss 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.bcf -Oz -o ALL.chr\${CHR}_GRCh38.bcf
  bcftools index --threads 16 ALL.chr\${CHR}_GRCh38.bcf
  Minimac3 --cpus 16 --refHaps ALL.chr\${CHR}_GRCh38.bcf --processReference --prefix ALL.chr\${CHR} --chr chr\${CHR}
  """}

process BCF_X{
  publishDir params.out+'files/ref/bcf/', mode: 'copy'

  input:
  file data from VCF_X.collect()
  val chromosome from "X"

  output:
  file '*X*' into Bcf_X

  shell:
  """
  chr=!{chromosome}
  (bcftools view --no-version -h 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz | grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\\([0-9XY]\\)/##contig=<ID=chr\\1/'; bcftools view --no-version -H -c 2 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz | sed 's/^/chr/') | bcftools norm --no-version -Ou -m -any | bcftools norm --no-version -Ob -o 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.bcf -d none -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && bcftools index -f 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.bcf
  """}
process BCF_Y{
  publishDir params.out+'files/ref/bcf/', mode: 'copy'

  input:
  file data from VCF_Y.collect()
  val chromosome from "Y"

  output:
  file '*Y*' into Bcf_Y

  shell:
  """
  chr=!{chromosome}
  (bcftools view --no-version -h 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz | grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\\([0-9XY]\\)/##contig=<ID=chr\\1/'; bcftools view --no-version -H -c 2 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz | sed 's/^/chr/') | bcftools norm --no-version -Ou -m -any | bcftools norm --no-version -Ob -o 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.bcf -d none -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && bcftools index -f 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.bcf
  """}
process M3VCF_X{
  validExitStatus 0,134,143,255
  publishDir params.out+'data/ref/m3vcf/', mode: 'copy'
  cpus=16

  input:
  val chromosome from "X"

  output:
  file '*m3vcf.gz' into M3VCF_X

  shell:
  """
  CHR=X
  bcftools index --tbi -f --threads 16 !{params.out}20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz -o 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz.tbi
  bcftools view --threads 16 -g ^miss !{params.out}20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz -Oz -o ALL.chr\${CHR}_GRCh38.vcf.gz
  bcftools index --threads 16 ALL.chr\${CHR}_GRCh38.vcf.gz
  /home/lipinskib/minimac3/Minimac3-master/bin/Minimac3 --cpus 16 --refHaps ALL.chr\${CHR}_GRCh38.vcf.gz --processReference --prefix ALL.chr\${CHR}
  """}

process M3VCF_Y{
  validExitStatus 0,134,143,255
  publishDir params.out+'data/ref/m3vcf/', mode: 'copy'
  cpus=16

  input:
  val chromosome from "Y"

  output:
  file '*m3vcf*' into M3VCF_Y

  shell:
  """
  CHR=Y
  bcftools index --tbi -f --threads 16 !{params.out}20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz -o 20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz.tbi
  bcftools view --threads 16 -g ^miss !{params.out}20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr\${chr}.recalibrated_variants.vcf.gz -Oz -o ALL.chr\${CHR}_GRCh38.vcf.gz
  bcftools index --threads 16 ALL.chr\${CHR}_GRCh38.vcf.gz
  /home/lipinskib/minimac3/Minimac3-master/bin/Minimac3 --cpus 16 --refHaps ALL.chr\${CHR}_GRCh38.vcf.gz --processReference --prefix ALL.chr\${CHR}
  """}
*/


process Download_files{
  publishDir params.out+'files/', mode: 'copy'

  output:
  file "{relationships_w_pops_121708.txt,checkVCF.py,GRCh38_full_analysis_set_plus_decoy_hla.fa,GRCh38_full_analysis_set_plus_decoy_hla.fa.fai}" into File_donwload

  shell:
  """
  wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/relationships_w_pops_121708.txt
  wget https://raw.githubusercontent.com/zhanxw/checkVCF/master/checkVCF.py
  wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
  wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

  """
}
process hapmap_peparation{
  publishDir params.out+'files/', mode: 'copy'

  input:
  file data from legend_file.collect()

  output:
  file "hapmap_r23a/" into Hapmap_file
  file "*.chain" into Chain_file

  shell:
  '''
  wget http://zzz.bwh.harvard.edu/plink/dist/hapmap_r23a.zip
  wget http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz
  wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
  wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2432498/bin/pone.0002551.s003.xls

  in2csv pone.0002551.s003.xls > pone.0002551.s003.csv

  unzip hapmap_r23a.zip
  gunzip *chain.gz

  awk '{print "chr" \$1, \$4 -1, \$4, \$2 }' hapmap_r23a.bim | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > dataset.tolift
  liftOver dataset.tolift hg18ToHg38.over.chain dataset1 dataset_NCBI36.unMapped

  awk '{print \$4}' dataset1 > dataset1.snps
  plink --bfile hapmap_r23a --extract dataset1.snps --make-bed --out dataset1

  awk '{print \$4, \$3}' dataset1  > dataset1.pos
  plink --bfile dataset1 --update-map dataset1.pos --make-bed --out dataset2

  #awk -F',' '{print \$1}' pone.0002551.s003.csv | grep -v SNPID > AIM_list.txt
  #plink --bfile dataset2 --extract AIM_list.txt --make-bed --out dataset3

  plink --freq --bfile dataset2 --make-bed --out dataset3
  perl !{baseDir}/bin/HRC-1000G-check-bim.pl -b dataset3.bim -f dataset3.frq -r PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz -h

  bash Run-plink.sh

  awk '{print \$2}' dataset3-updated.fam > ref_samples.txt

  mkdir hapmap_r23a/
  mv ref_samples.txt hapmap_r23a/ref_samples.txt
  #mv AIM_list.txt hapmap_r23a/AIM_list.txt
  mv dataset3-updated.bed hapmap_r23a/hapmap_r23a.bed
  mv dataset3-updated.bim hapmap_r23a/hapmap_r23a.bim
  mv dataset3-updated.fam hapmap_r23a/hapmap_r23a.fam
  '''
}
