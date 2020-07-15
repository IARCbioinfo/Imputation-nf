#! /usr/bin/env nextflow

// Copyright (C) 2017 IARC/WHO

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
log.info "  Imputation_genotypage 1.0 : Pipeline for the imputation of a target dataset against a reference "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE : nextflow run Imputation_genotypage.nf --ref hapmap_r23a --target HGDP /// nextflow run ../../Imputation-nf/Imputation_genotypage.nf --ref hapmap_r23a --target HGDP --input /data/gep/MR_Signatures/work/Boris/protocol_min/data/ --script ../../Imputation-nf/bin/ "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run iarcbioinfo/template-nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--target                      pattern                      pattern of the target dataset which do the link with the file .bed/.bim./fam for plink"
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

// ## -- Option :
params.ref = "hapmap_r23a"
params.target = null      //HGDP & TCGA
params.geno1 = 0.03 
params.geno2 = 0.03
params.maf = 0.01
params.pihat = 0.185
params.hwe = 1e-8 

// ## -- Path :
params.input = null // '/data/gep/MR_Signatures/work/Boris/protocol_min/data/'
params.folder = params.input+'files/'  // '/data/gep/MR_Signatures/work/Boris/protocol_min/data/files/'
params.refDir = params.input+params.ref+'/'
params.targetDir = params.input+params.target+'/'
params.script = null // '/data/gep/MR_Signatures/work/Boris/protocol_min/script/bin/'
params.out = params.input

/*
process UpdateMap{ 

  input:
  file data from Channel.fromPath(params.targetDir+"*").collect()
  file data from Channel.fromPath(params.folder+'HRC-1000G-check-bim-NoReadKey.pl').collect()
  file data from Channel.fromPath(params.folder+'ALL.chr_test.legend').collect()

  output:
  file ('*-updated.{bed,bim,fam}') into TargetMap

  shell:
  '''
  ############################################################################################
  ## -- 0 : Update map on the tergat dataset
  plink --freq --bfile !{params.target} --out !{params.target}_frep
  perl HRC-1000G-check-bim-NoReadKey.pl -b !{params.target}.bim -f !{params.target}_frep.frq -r ALL.chr_test.legend -g -x -n
  grep -v "real-ref-alleles" Run-plink.sh> Run-plink-update.sh 
  bash Run-plink-update.sh
  '''}*/
process Admixture{
  input:
  file data from Channel.fromPath(params.folder+'relationships_w_pops_121708.txt').collect()
  file data from Channel.fromPath(params.refDir+"*").collect()
  file data from Channel.fromPath(params.targetDir+"*").collect()
  //file data from TargetMap.collect()


  output:
  file ('target4.{bed,bim,fam}') into Merge
  file ('target4.{bed,bim,fam}') into Merge2
  file ('out_pop_admixture/') into Admixture
  file ('out_pop_admixture/') into Admixture2
  file ('admixture_results_withGroups.txt') into Target

  shell:
  '''
  ## -- 1 : Retrieve common SNPs between SNPs in the reference list and the target data.
  ## -- Ref -- ##
  awk '{print $2}' !{params.ref}.bim | sort > ref_SNPs.txt

  ## -- Target -- ##
  grep -Fwf AIM_list.txt !{params.target}.bim | awk '{print $2}' > target_SNPs.txt #-updated
  plink --bfile !{params.target} --extract target_SNPs.txt --make-bed --out target #-updated

  ## -- Get common SNPs -- ##
  grep -Fwf <(cat ref_SNPs.txt) <(awk '{print $2}' target.bim)  > target_common_SNPs.txt
  awk -F ":" '{print $1}' target_common_SNPs.txt > ref_common_SNPs.txt
  paste target_common_SNPs.txt ref_common_SNPs.txt > change_ID.txt


  ## -- 2 : Extract the genotypes associated to these common SNPs + merge of the dataset
  plink --bfile !{params.ref} --extract ref_common_SNPs.txt --make-bed --out ref1
  plink --bfile target --extract target_common_SNPs.txt --make-bed --out target1
  plink --bfile target1 --update-map change_ID.txt --update-name --make-bed --out target2
  plink --bfile target2 --bmerge ref1 --allow-no-sex --make-bed --out merge

  awk '{print $2}' merge.fam  >  all_samples.txt

  ## -- 3 : Associate each reference sample to the correct origin + run admixure
  sort -k2 relationships_w_pops_121708.txt > relationships_w_pops_121708_2.txt
  awk '{print $2}' !{params.ref}.fam | sort -k1,1 > ref_ind.txt
  join -11 -22 ref_ind.txt relationships_w_pops_121708_2.txt | awk '{print $1, $7}' | awk '{if( $2 == "CEU") print $0,1,1; else { if($2 == "YRI") print $0,2,1; else {print $0,3,1}}}' > ind_pop.txt
  Rscript !{baseDir}/bin/create_pop_file.r merge.fam ind_pop.txt merge.pop
  K=3
  admixture --cv merge.bed $K --supervised -j40 | tee log${K}.out
  Rscript !{baseDir}/bin/process_admixture.r !{params.target}

  ############################################################################################
  ## -- 4 : First filtering step
  plink --bfile !{params.target} --geno !{params.geno1} --make-bed --out target3
  plink --bfile target3 --maf !{params.maf} --make-bed --out target4
  '''}
process Filtering1{
  input:
  file data from Merge.collect()
  file data from Admixture.collect()
  val rspop from Channel.from("CEU","CHB_JPT","YRI")

  output:
  file ('*.{het,imiss.txt,genome,sexcheck}') into QC
  file ('target_*.bim') into Bim
  file ('target_*.{bed,bim,fam}') into TargetFilter

  shell:
  '''
  pop=!{rspop}
  subpop=out_pop_admixture/1000G_checking_$pop/samples_$pop.txt

  ## -- 5 : Keep only the samples in the ancerstry $pop
  plink --bfile target4 --keep ${subpop} --make-bed --out target5_${pop}

  ## -- 6 : Perform sex checking
  plink --bfile target5_${pop} --make-bed --merge-x no-fail --out target6_${pop}
  plink --bfile target6_${pop} --make-bed --split-x b37 no-fail --out target_${pop}
  plink --bfile target_${pop} --indep-pairwise 50 5 0.2 --out target_Independent_SNPs_${pop} # keep independent SNPs
  plink --bfile target_${pop} --extract target_Independent_SNPs_${pop}.prune.in --check-sex --out target_sexCheck_${pop}

  ## -- 7 : Test relatedness between samples
  plink --bfile target_${pop} --extract target_Independent_SNPs_${pop}.prune.in  --genome --out target_rel_${pop} --min !{params.pihat}

  ## -- 8 : Test heterozygocity and missing rates
  plink --bfile target_${pop} --het --chr 1-22  --out het_${pop}
  echo "FID IID obs_HOM N_SNPs prop_HET" > het_${pop}.txt
  awk 'NR>1 {print $1,$2,$3,$5,($5-$3)/$5}' het_${pop}.het >> het_${pop}.txt
  plink --bfile target_${pop} --missing  --out miss_${pop}
  awk 'NR==FNR {a[$1,$2]=$5;next}($1,$2) in a{print $1,$2,$6,a[$1,$2]}' het_${pop}.txt miss_${pop}.imiss > het_${pop}.imiss.txt
  '''}
process QC1{
  publishDir params.out+'../result/'+params.target+'/QC1/', mode: 'copy'

  input:
  file data from QC.collect()
  file target from Target.collect()

  output:
  file ('postGenotyping_samples_QCs.pdf') into FigureQC1
  file ('selected_samples_afterSamplesQCsChecking.txt') into TableQC1

  shell:
  '''
  ## -- 9 : Figure QC1
  Rscript !{baseDir}/bin/after_genotyping_qc.r
  '''}
process Filtering2{
  input:
  file data from TargetFilter.collect()
  file data from Channel.fromPath(params.folder+'HRC-1000G-check-bim-NoReadKey.pl').collect()
  file data from Channel.fromPath(params.folder+'ALL.chr_test.legend').collect()
  val rspop from Channel.from("CEU","CHB_JPT","YRI")

  output:
  file ('withFreqFiltering_*') into DirFiltering
  file ('target_hwe_*.bim') into HWresult
  file ('target_freq_*.frq') into FreqResult
  file ('ID-target_*-1000G.txt') into FreqResultId

  shell:
  '''
  pop=!{rspop}

  if [ $pop = "CEU" ]; then
    pop2="EUR"
  elif [ $pop = "CHB_JPT" ]; then
    pop2="EAS"
  else
    pop2="AFR"
  fi

  ## -- 10 : AF based filter
  plink --freq --bfile target_${pop} --output-chr chr26 --out target_freq_${pop}
  perl HRC-1000G-check-bim-NoReadKey.pl -b target_${pop}.bim -f target_freq_${pop}.frq -r ALL.chr_test.legend -g -p ${pop2} -x
  mkdir withFreqFiltering_${pop}
  cp *1000G* Run-plink.sh withFreqFiltering_${pop}

  ## -- 11 :  HWE filtering
  plink --bfile target_${pop} --geno !{params.geno2} --make-bed --out target_geno_${pop}
  plink --bfile target_geno_${pop} --hwe !{params.hwe} --make-bed --out target_hwe_${pop}
  '''}


process Make_SNP_Filtering{
  publishDir params.out+'../result/'+params.target+'/SNP_filtering/', mode: 'copy'
  input:
  file data from DirFiltering.collect()
  file data from HWresult.collect()
  file data from Bim.collect()

  output:
  file ('filtered_snps.txt') into SNPsFilter
  file ('filtered_snps.txt') into SNPsFilter2
  file ('*.pdf') into resultFigure

  shell:
  '''
  ## -- 12 : Create list of SNPs to filter
  Rscript !{baseDir}/bin/afterGenotyping_SNPs_filtering.r
  '''}
process Filtering3{
  input:
  file data from Merge2.collect()
  file data from Channel.fromPath(params.folder+'HRC-1000G-check-bim-NoReadKey.pl').collect()
  file data from Channel.fromPath(params.folder+'ALL.chr_test.legend').collect()
  file data from SNPsFilter.collect()

  output:
  file ('target5-updated-chr*') into TargetChr
  //file ('ID-target5-1000G.txt') into TargetID

  shell:
  '''
  ## -- 13 : Removing SNPs
  plink --bfile target4 --exclude filtered_snps.txt  --make-bed --out target5

  ## -- 14 : filter
  plink --freq --bfile target5  --out target6
  perl HRC-1000G-check-bim-NoReadKey.pl -b target5.bim -f target6.frq -r ALL.chr_test.legend -g -x -n
  bash Run-plink.sh
  '''}

process QC2{
  publishDir params.out+'../result/'+params.target+'/QC2/', mode: 'copy'

  input:
  file data from SNPsFilter2.collect()
  file data from FreqResult.collect()
  file data from FreqResultId.collect()
  //file data from TargetID.collect()
  val rspop from Channel.from("CEU","CHB_JPT","YRI")
  file data from Channel.fromPath(params.folder+'ALL.chr_test.legend').collect()
  file data from Channel.fromPath(params.targetDir+"*.bim").collect()

  output:
  file ('*.pdf') into FigureQC2

  shell:
  '''
  head -n1 ALL.chr_test.legend >> ref_freq_withHeader.txt
  grep -Fwf <(awk '{print $2}' !{params.target}.bim) <(cat ALL.chr_test.legend) > ref_freq.txt
  cat ref_freq.txt >> ref_freq_withHeader.txt

  pop=!{rspop}

  if [ $pop = "CEU" ]; then
    pop2="EUR"
  elif [ $pop = "CHB_JPT" ]; then
    pop2="EAS"
  else
    pop2="AFR"
  fi

  ## -- 15 : Figure QC2
  Rscript !{baseDir}/bin/preImputation_QC_plots.r ${pop} ${pop2}
  '''}
process Filtering4{
  input:
  file data from TargetChr.collect()
  file data from Channel.fromPath(params.folder+'checkVCF.py').collect()
  file data from Channel.fromPath(params.folder+'GRCh38_full_analysis_set_plus_decoy_hla.fa*').collect() //#human_g1k_v37.fasta
  val chromosome from 1..22 //23

  output:
  file ('*-REFfixed.vcf.gz') into FilterFinal
  file ('*-REFfixed.vcf.gz') into FilterFinal2

  shell:
  '''
  chr=!{chromosome}

  ## -- 16 : Remove ambiguous strand/unknown SNPs
  awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' target5-updated-chr${chr}.bim | grep -v ambig | grep -v -e --- | sort -u > NonAmbiguous${chr}.snplist.txt
  plink --bfile target5-updated-chr${chr} --extract NonAmbiguous${chr}.snplist.txt --output-chr chr26 --make-bed --out target6_chr${chr}

  ## -- 17 : Create VCF
  plink  --bfile target6_chr${chr}  --output-chr chr26 --recode vcf --out target6_chr${chr}_vcf
  bcftools sort target6_chr${chr}_vcf.vcf | bgzip -c  > chr${chr}.vcf.gz

  ## -- 18 : Check SNPs
  python2 checkVCF.py -r GRCh38_full_analysis_set_plus_decoy_hla.fa -o after_check_${chr} chr${chr}.vcf.gz #human_g1k_v37.fasta
  bcftools norm --check-ref ws -f GRCh38_full_analysis_set_plus_decoy_hla.fa chr${chr}.vcf.gz | bcftools view -m 2 -M 2  | bgzip -c > chr${chr}-REFfixed.vcf.gz #human_g1k_v37.fasta
  python2 checkVCF.py -r GRCh38_full_analysis_set_plus_decoy_hla.fa -o after_check2_${chr} chr${chr}-REFfixed.vcf.gz #human_g1k_v37.fasta
  '''}


process Make_Chunks{
  input:
  val chromosome from 1..22
  file data from FilterFinal.collect()
  file data from Channel.fromPath("/data/gep/MR_Signatures/work/Boris/protocol_min/data/files/ref/vcf/*").collect()

  output:
  env chunks into NbChunk
  env chunks into NbChunk2
  tuple env(chunks), val(chromosome) into InfoChrChunk
  file ('*.txt') into ChunkSplit

  shell:
  '''
  ## -- 19 : Create Chunks
  chr=!{chromosome}
  #ref_haps="/data/gep/MR_Signatures/work/Boris/protocol_min/data/files/ref/vcf/"   #"/data/gep/MR_Signatures/work/gabriela/imputation/January_2020/ref_data/"
  bcftools index -f chr${chr}-REFfixed.vcf.gz
  bcftools isec -n +2 chr${chr}-REFfixed.vcf.gz ALL.chr${chr}_GRCh38.genotypes.20170504.bcf| bgzip -c > isec_chr_${chr}.vcf.gz #1000GP_chr${chr}.bcf #${ref_haps}
  Rscript !{baseDir}/bin/create_chunks.r ${chr}
  chunks=$(wc -l chunk_split_chr${chr}.txt | awk '{print $1}')
  '''}
process Make_Multiprocessing{
  input:
  val merge from InfoChrChunk.toList()

  output:
  file 'multiprocess.txt' into NbChr
  file 'multiprocess.txt' into NbChr2

  shell:
  '''
  #!/usr/bin/env python3
  ## -- 20 : Preparation of multiprocessing for imputation

  file=open('multiprocess.txt','w')
  liste=!{merge}
  for tuples in liste:
     for i in range(0,tuples[0]):
        file.write(str(tuples[1]) + '\\n')

  file.close()
  '''}
process Imputation{
  input:
  val chunks from NbChunk.map{1.."$it".toInteger()}.flatten()
  val chromosomes from NbChr.splitText()

  file data from Channel.fromPath("/data/gep/MR_Signatures/work/Boris/protocol_min/data/files/ref/vcf/*").collect()
  file data from Channel.fromPath("/data/gep/MR_Signatures/work/Boris/protocol_min/data/files/ref/m3vcf/*").collect()


  file data from FilterFinal2.collect()
  file data from ChunkSplit.collect()

  output:
  file '*imputed.dose.vcf.gz*' into FileVCF

  shell:
  '''
  chr=!{chromosomes}
  chunk=!{chunks}
  echo "chr: ${chr} n_chunks: ${chunk}"
  cpu=1
  #ref_haps="/data/gep/MR_Signatures/work/Boris/protocol_min/data/files/ref/vcf/" #"/data/gep/MR_Signatures/work/gabriela/imputation/January_2020/ref_data/"
  start=$(awk '{print $1}' <(awk 'NR == c' c="${chunk}" chunk_split_chr${chr}.txt))
  end=$(awk '{print $2}' <(awk 'NR == c' c="${chunk}" chunk_split_chr${chr}.txt))

  ## -- 21 : Phasing
  bcftools index -f chr${chr}-REFfixed.vcf.gz
  eagle --vcfRef ALL.chr${chr}_GRCh38.genotypes.20170504.bcf --vcfTarget chr${chr}-REFfixed.vcf.gz --vcfOutFormat v --geneticMapFile /home/lipinskib/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz --outPrefix chr_${chr}_chunk${chunk}.phased --bpStart ${start} --bpEnd ${end} --bpFlanking 5000000 --chrom ${chr} --numThreads ${cpu}  > chr_${chr}_chunk${chunk}_phasing.logphase
  #${ref_haps}
  sed -i "s/chr${chr}/${chr}/g" chr_${chr}_chunk${chunk}.phased.vcf

  ## -- 22 : Imputation
  #ref_haps="/data/gep/MR_Signatures/work/Boris/protocol_min/data/files/ref/m3vcf/"    #="/data/references/Homo_sapiens/ref_haps_1000G_phase3/hg19/m3vcf/"
  minimac4 --refHaps ALL.chr${chr}.m3vcf.gz --haps chr_${chr}_chunk${chunk}.phased.vcf --prefix chr_${chr}_chunk${chunk}.imputed --allTypedSites --format GT,DS,GP --cpus ${cpu} --chr ${chr} --start $start --end $end --window 500000 > chr_${chr}_chunk${chunk}.logimpute
  #${ref_haps}
  bcftools index -f chr_${chr}_chunk${chunk}.imputed.dose.vcf.gz
  '''}
process Concatenation{
  publishDir params.out+'../result/'+params.target+'/Result_Imputation/', mode: 'copy'
  cpus=6

  input:
  val chromosomes from 1..22
  file data from FileVCF.collect()

  output:
  file '*_combined.vcf.gz' into Imputation

  shell:
  '''
  chr=!{chromosomes}
  mkdir chr_${chr}
  ls chr_${chr}_chunk*.imputed.dose.vcf.gz> chr_${chr}_imp_res.txt
  ## -- 23 : Concatenation
  bcftools concat --threads 6 -f chr_${chr}_imp_res.txt -Ou | bcftools sort --temp-dir chr_${chr} -Ou | bgzip -c > chr_${chr}_combined.vcf.gz 
  '''}
process QC3_sh{
  input:
  val population from('ALL','CEU','YRI','CHB_JPT')
  each chromosome from 1..22
  file data from Imputation.collect()
  file data from Admixture2.collect()
  file data from Channel.fromPath("/data/gep/MR_Signatures/work/Boris/protocol_min/data/files/ref/vcf/*.bcf*").collect()

  output:
  file '*.{txt,frq}' into PostImputation_QC_sh_result

  shell:
  '''
  pop=!{population}
  chr=!{chromosome}
  ## -- 24 : QC3
  bash !{baseDir}/bin/postImputation_QC.sh ${chr} ${pop}
  '''}
process QC3_R{
  publishDir params.out+'../result/'+params.target+'/QC3/', mode: 'copy'
  input:
  val population from('ALL','CEU','YRI','CHB_JPT')
  file data from PostImputation_QC_sh_result.collect()

  output:
  file '*v2.{txt,pdf}' into PostImputation_QC_R_result

  shell:
  '''
  pop=!{population}
  Rscript !{baseDir}/bin/postImputation_QC_plots.r ${pop} 0.3
  '''}



