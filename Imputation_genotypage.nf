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
    log.info "--ref                      pattern                      pattern of the reference dataset which do the link with the file .bed/.bim./fam for plink"
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

params.ref = null //hapmap_r23a
params.target = null //HGDP

params.input = null // '/data/gep/MR_Signatures/work/Boris/protocol_min/data/'
params.folder = params.input+'data_start/'  // '/data/gep/MR_Signatures/work/Boris/protocol_min/data/data_start/'
params.refDir = params.input+params.ref+'/'
params.targetDir = params.input+params.target+'/'

params.script = null // '/data/gep/MR_Signatures/work/Boris/protocol_min/script/bin/'
params.out = params.input // '/data/gep/MR_Signatures/work/Boris/protocol_min/data/'

process merge_admixture{
  input:
  file data_start from Channel.fromPath(params.folder+'relationships_w_pops_121708.txt').collect()
  file data_start from Channel.fromPath(params.refDir+"*").collect()
  file data_start from Channel.fromPath(params.targetDir+"*").collect()

  output:
  file ('target4.{bed,bim,fam}') into resultMerge
  file ('target4.{bed,bim,fam}') into resultMerge2
  file ('admixture_results_withGroups.txt') into resultTarget
  file ('out_pop_admixture/') into resultAdmixture

  shell:
  '''
  ############################################################################################
  ## -- 1 : Retrieve common SNPs between SNPs in the reference list and the target data.
  ## -- Ref -- ##
  awk '{print $2}' !{params.ref}.bim | sort > ref_SNPs.txt

  ## -- Target -- ##
  grep -Fwf AIM_list.txt !{params.target}.bim | awk '{print $2}' > target_SNPs.txt
  plink --bfile !{params.target} --extract target_SNPs.txt --make-bed --out target

  ## -- Get common SNPs -- ##
  grep -Fwf <(cat ref_SNPs.txt) <(awk '{print $2}' target.bim)  > target_common_SNPs.txt
  awk -F ":" '{print $1}' target_common_SNPs.txt > ref_common_SNPs.txt
  paste target_common_SNPs.txt  ref_common_SNPs.txt > change_ID.txt


  ## -- 2 : Extract the genotypes associated to these common SNPs + merge of the dataset
  plink --bfile !{params.ref} --extract ref_common_SNPs.txt --make-bed --out ref1
  plink --bfile target --extract target_common_SNPs.txt --make-bed --out target1
  plink --bfile target1 --update-map change_ID.txt --update-name --make-bed --out target2
  plink --bfile target2 --bmerge ref1 --make-bed --out merge

  awk '{print $2}' merge.fam  >  all_samples.txt


  ## -- 3 : Associate each reference sample to the correct origin + run admixure
  sort -k2 relationships_w_pops_121708.txt > relationships_w_pops_121708_2.txt
  awk '{print $2}' !{params.ref}.fam | sort -k1,1 > ref_ind.txt
  join -11 -22 ref_ind.txt relationships_w_pops_121708_2.txt | awk '{print $1, $7}' | awk '{if( $2 == "CEU") print $0,1,1; else { if($2 == "YRI") print $0,2,1; else {print $0,3,1}}}' > ind_pop.txt
  Rscript !{baseDir}/bin/create_pop_file.r merge.fam ind_pop.txt merge.pop
  K=3
  admixture --cv merge.bed $K --supervised -j40 | tee log${K}.out
  Rscript !{baseDir}/bin/process_admixture.r #### regarder les sorties pour chopper les id


  ############################################################################################
  ## -- 4 : First filtering step
  plink --bfile !{params.target}  --geno 0.03 --make-bed --out target3
  plink --bfile target3 --maf 0.01 --make-bed --out target4
  '''
}

process processing_filtering1{
  //publishDir params.out, mode: 'copy'

  input:
  file data from resultMerge.collect()
  file dir from resultAdmixture.collect()
  val rspop from Channel.from("CEU","CHB_JPT","YRI")

  output:
  file ('*.het') into resultHet
  file ('*.imiss.txt') into resultHetMiss
  file ('*.genome') into resultGenome
  file ('*.sexcheck') into resultSex
  file ('target_*.{bed,bim,fam}') into resultFiltering
  file ('target_*.bim') into resultBim

  shell:
  '''
  pop=!{rspop}
  #subpop_path=out_pop_admixture/1000G_checking_$pop/
  subpop=out_pop_admixture/1000G_checking_$pop/samples_$pop.txt

  ## -- 5 : Keep only the samples in the ancerstry $pop
  plink --bfile target4 --keep ${subpop} --make-bed --out target5_${pop}


  ## -- 6 : Perform sex checking
  #plink --bfile target5_${pop} --make-bed --merge-x --out target6_${pop}
  plink --bfile target5_${pop} --make-bed --split-x b37 no-fail --out target_${pop}
  plink --bfile target_${pop} --indep-pairwise 50 5 0.2 --out target_Independent_SNPs_${pop} # keep independent SNPs
  plink --bfile target_${pop} --extract target_Independent_SNPs_${pop}.prune.in --check-sex --out target_sexCheck_${pop}


  ## -- 7 : Test relatedness between samples
  plink --bfile target_${pop} --extract target_Independent_SNPs_${pop}.prune.in --genome --out target_rel_${pop} --min 0.185


  ## -- 8 : Test heterozygocity and missing rates
  plink --bfile target_${pop} --het --chr 1-22 --out het_${pop}
  echo "FID IID obs_HOM N_SNPs prop_HET" > het_${pop}.txt
  awk 'NR>1 {print $1,$2,$3,$5,($5-$3)/$5}' het_${pop}.het >> het_${pop}.txt
  plink --bfile target_${pop} --missing --out miss_${pop}
  awk 'NR==FNR {a[$1,$2]=$5;next}($1,$2) in a{print $1,$2,$6,a[$1,$2]}' het_${pop}.txt miss_${pop}.imiss > het_${pop}.imiss.txt
  '''
}


process QC_analyse{
  publishDir params.out, mode: 'copy'

  input:
  file het from resultHet.collect()
  file miss from resultHetMiss.collect()
  file genome from resultGenome.collect()
  file sexcheck from resultSex.collect()
  file target from resultTarget.collect()

  output:
  file ('postGenotyping_samples_QCs.pdf') into resultFigureQC
  file ('selected_samples_afterSamplesQCsChecking.txt') into resultTableQC

  shell:
  '''
  ## -- 9 : Figure time
  Rscript !{baseDir}/bin/after_genotyping_qc.r
  '''
}


process processing_filtering2{
  input:
  file data from resultFiltering.collect()
  file data2 from Channel.fromPath(params.folder+'HRC-1000G-check-bim.pl').collect()
  file data3 from Channel.fromPath(params.folder+'1000GP_Phase3_combined.legend').collect()
  val rspop from Channel.from("CEU","CHB_JPT","YRI")

  output:
  file ('withFreqFiltering_*') into resultDirFiltering
  file ('target_hwe_*.bim') into resultHW

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
  plink --freq --bfile target_${pop} --out target_freq_${pop}
  perl HRC-1000G-check-bim.pl -b target_${pop}.bim -f target_freq_${pop}.frq -r 1000GP_Phase3_combined.legend -g -p ${pop2} -x
  mkdir withFreqFiltering_${pop}
  mv *1000G* Run-plink.sh withFreqFiltering_${pop}

  ## -- 11 :  HWE filtering
  plink --bfile target_${pop} --geno 0.03 --make-bed --out target_geno_${pop}
  plink --bfile target_geno_${pop} --hwe 1e-8 --make-bed --out target_hwe_${pop}
  '''
}


process make_snp_filtering{
  publishDir params.out, mode: 'copy'
  input:
  file directories from resultDirFiltering.collect()
  file dataHW from resultHW.collect()
  file dataBim from resultBim.collect()

  output:
  file ('filtered_snps.txt') into fileSNPS
  file ('*.pdf') into resultFigure

  shell:
  '''
  ## -- 12 : filtering SNP
  Rscript !{baseDir}/bin/afterGenotyping_SNPs_filtering.r
  '''
}

process processing_filtering3{
  input:
  file data from resultMerge2.collect()
  file data2 from Channel.fromPath(params.folder+'HRC-1000G-check-bim-NoReadKey.pl').collect()
  file data3 from Channel.fromPath(params.folder+'1000GP_Phase3_combined.legend').collect()
  file data4 from fileSNPS.collect()

  output:
  file ('target5-updated-chr*') into resultChr

  shell:
  '''
  ## -- 13 : removing SNP
  plink --bfile target4 --exclude filtered_snps.txt --make-bed --out target5

  plink --freq --bfile target5 --out target6
  perl HRC-1000G-check-bim-NoReadKey.pl -b target5.bim -f target6.frq -r 1000GP_Phase3_combined.legend -g -x -n
  bash Run-plink.sh
  '''
}


process processing_filtering4{
  input:
  file dataChr from resultChr.collect()
  file data from Channel.fromPath(params.folder+'checkVCF.py').collect()
  file data2 from Channel.fromPath(params.folder+'human_g1k_v37.fasta').collect()
  val chromosome from Channel.from(1..23)

  output:
  file ('*-REFfixed.vcf.gz') into resultFilterFinal

  shell:
  '''
  chr=!{chromosome}

  #Remove ambiguous strand snps and SNPs with unknown SNPs ID
  awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' target5-updated-chr${chr}.bim | grep -v ambig | grep -v -e --- | sort -u > NonAmbiguous{chr}.snplist.txt
  plink --bfile target5-updated-chr${chr} --extract NonAmbiguous{chr}.snplist.txt --make-bed --out target6_chr${chr}

  #Create VCF / sort and zip
  plink  --bfile target6_chr${chr} --recode vcf --out target6_chr${chr}_vcf
  bcftools sort target6_chr${chr}_vcf.vcf  | bgzip -c  > chr${chr}.vcf.gz

  #Check SNPs using check_VCF.py
  python2 checkVCF.py -r human_g1k_v37.fasta -o after_check_${chr} chr${chr}.vcf.gz
  bcftools norm --check-ref ws -f human_g1k_v37.fasta chr${chr}.vcf.gz | bcftools view -m 2 -M 2  | bgzip -c > chr${chr}-REFfixed.vcf.gz
  python2 checkVCF.py -r human_g1k_v37.fasta -o after_check2_${chr} chr${chr}-REFfixed.vcf.gz
  '''
}

process phasing{
  publishDir params.out, mode: 'copy'

  input:
  val chromosome from Channel.from(1..22)
  file data from resultFilterFinal.collect()

  output:
  file ('*.logimpute') into imputationFinal
  
  shell:
  '''
  chr=!{chromosome}

  Rscript !{baseDir}/bin/create_chunks.r ${chr}

  n_chunks=$(wc -l chunk_split_chr${chr}.txt | awk '{print $1}')
  echo "chr: ${chr} n_chunks: ${n_chunks}"
  for chunk in `seq 1 ${n_chunks}`;
  do
    cpu=5
    #path_plink=${imputation_path}
    #path_chr=${path_plink}chr_${chr}/
    ref_haps="/data/gep/MR_Signatures/work/gabriela/imputation/January_2020/ref_data/"

    start=$(awk '{print $1}' <(awk 'NR == c' c="${chunk}" chunk_split_chr${chr}.txt))
    end=$(awk '{print $2}' <(awk 'NR == c' c="${chunk}" chunk_split_chr${chr}.txt))

    bcftools index -f chr${chr}-REFfixed.vcf.gz
    /home/lipinskib/Eagle_v2.4.1/./eagle --vcfRef ${ref_haps}1000GP_chr${chr}.bcf --vcfTarget chr${chr}-REFfixed.vcf.gz --vcfOutFormat v --geneticMapFile /home/lipinskib/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz --outPrefix chr_${chr}_chunk${chunk}.phased --bpStart ${start} --bpEnd ${end} --bpFlanking 5000000 --chrom ${chr} --numThreads ${cpu}  > chr_${chr}_chunk${chunk}_phasing.logphase

    ref_haps="/data/references/Homo_sapiens/ref_haps_1000G_phase3/hg19/m3vcf/"
    /home/lipinskib/Minimac4/release-build/./minimac4 --refHaps ${ref_haps}${chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz --haps chr_${chr}_chunk${chunk}.phased.vcf --prefix chr_${chr}_chunk${chunk}.imputed --allTypedSites --format GT,DS,GP --cpus ${cpu} --chr ${chr} --start $start --end $end --window 500000 > chr_${chr}_chunk${chunk}.logimpute
  done
  '''
}
