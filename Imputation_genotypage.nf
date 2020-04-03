#!/usr/bin/env nextflow

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
    log.info "  USAGE : nextflow run Imputation_genotypage.nf --ref hapmap_r23a --target HGDP  "
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
  val data_ref from Channel.fromPath(params.refDir+params.ref)
  val data_target from Channel.fromPath(params.targetDir+params.target)
  file data from Channel.fromPath(params.script+'*.r').collect()
  file data_start from Channel.fromPath(params.folder+'relationships_w_pops_121708.txt').collect()
  file data_start from Channel.fromPath(params.refDir+"*").collect()

  output:
  file ('target4.{bed,bim,fam}') into resultMerge
  file ('out_pop_admixture/') into resultAdmixture

  script:
  """
  ############################################################################################
  ## -- 1 : Retrieve common SNPs between SNPs in the reference list and the target data.
  ## -- Ref -- ##

  ref="$data_ref"
  target="$data_target"

  awk '{print \$2}' \${ref%%.*}.bim | sort > ref_SNPs.txt

  ## -- Target -- ##
  grep -Fwf AIM_list.txt \${target%%.*}.bim | awk '{print \$2}' > target_SNPs.txt
  plink --bfile \${target%%.*} --extract target_SNPs.txt --make-bed --out target

  ## -- Get common SNPs -- ##
  grep -Fwf <(cat ref_SNPs.txt) <(awk '{print \$2}' target.bim)  > target_common_SNPs.txt
  awk -F ":" '{print \$1}' target_common_SNPs.txt > ref_common_SNPs.txt
  paste target_common_SNPs.txt  ref_common_SNPs.txt > change_ID.txt


  ## -- 2 : Extract the genotypes associated to these common SNPs + merge of the dataset
  plink --bfile \${ref%%.*} --extract ref_common_SNPs.txt --make-bed --out ref1
  plink --bfile target --extract target_common_SNPs.txt --make-bed --out target1
  plink --bfile target1 --update-map change_ID.txt --update-name --make-bed --out target2
  plink --bfile target2 --bmerge ref1 --make-bed --out merge

  awk '{print \$2}' merge.fam  >  all_samples.txt


  ## -- 3 : Associate each reference sample to the correct origin + run admixure
  sort -k2 relationships_w_pops_121708.txt > relationships_w_pops_121708_2.txt
  awk '{print \$2}' \${ref%%.*}.fam | sort -k1,1 > ref_ind.txt
  join -11 -22 ref_ind.txt relationships_w_pops_121708_2.txt | awk '{print \$1,\$7}' | awk '{if(\$2 == "CEU") print \$0,1,1; else { if(\$2 == "YRI") print \$0,2,1; else {print \$0,3,1}}}' > ind_pop.txt
  Rscript create_pop_file.r merge.fam ind_pop.txt merge.pop
  K=3
  admixture --cv merge.bed \$K --supervised -j40 | tee log\${K}.out
  Rscript process_admixture.r


  ############################################################################################
  ## -- 4 : First filtration step
  plink --bfile \${target%%.*}  --geno 0.03 --make-bed --out target3
  plink --bfile target3 --maf 0.01 --make-bed --out target4
  """
}

process processing_filtration{
  publishDir params.out, mode: 'copy'

  input:
  file data from resultMerge.collect()
  file dir from resultAdmixture.collect()
  val rspop from Channel.from("CEU","CHB_JPT","YRI")

  output:
  file ('out_pop_admixture/') into resultFiltration

  script:
  """
  pop=$rspop
  subpop_path=out_pop_admixture/1000G_checking_\$pop/
  subpop=out_pop_admixture/1000G_checking_\$pop/samples_\$pop.txt

  ## -- 5 : Keep only the samples in the ancerstry \$pop
  plink --bfile target4 --keep \${subpop} --make-bed --out target5


  ## -- 6 : Perform sex checking
  #plink --bfile target5 --make-bed --merge-x --out target6
  plink --bfile target5 --make-bed --split-x b37 no-fail --out \${subpop_path}target_subpop
  plink --bfile \${subpop_path}target_subpop --indep-pairwise 50 5 0.2 --out target_Independent_SNPs # keep independent SNPs
  plink --bfile \${subpop_path}target_subpop --extract target_Independent_SNPs.prune.in --check-sex --out \${subpop_path}target_subpop_sexCheck


  ## -- 7 : Test relatedness between samples
  plink --bfile \${subpop_path}target_subpop --extract target_Independent_SNPs.prune.in --genome --out \${subpop_path}target_subpop_rel --min 0.185


  ## -- 8 : Test heterozygocity and missing rates
  plink --bfile \${subpop_path}target_subpop --het --chr 1-22 --out \${subpop_path}het
  echo "FID IID obs_HOM N_SNPs prop_HET" > \${subpop_path}het.txt
  awk 'NR>1 {print \$1,\$2,\$3,\$5,(\$5-\$3)/\$5}' \${subpop_path}het.het >> \${subpop_path}het2.txt
  plink --bfile \${subpop_path}target_subpop --missing --out \${subpop_path}miss
  awk 'NR==FNR {a[\$1,\$2]=\$5;next}(\$1,\$2) in a{print \$1,\$2,\$6,a[\$1,\$2]}' \${subpop_path}het2.txt \${subpop_path}miss.imiss > \${subpop_path}het.imiss.txt
  """
}
