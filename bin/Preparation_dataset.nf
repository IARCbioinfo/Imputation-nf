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
log.info "  Preparation_dataset.nf 1.0 : Dataset's preparation for the imputation's pipeline "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE : nextflow run Pipeline_preparation.nf --dataset target                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run Pipeline_preparation.nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--dataset                    pattern                      pattern of the target/reference dataset which link with the file .bed/.bim./fam for plink"
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

params.folder = '/data/gep/MR_Signatures/work/Boris/protocol_min/data/data_start/'
params.out = '/data/gep/MR_Signatures/work/Boris/protocol_min/data/'

params.dataset = null

data = Channel.fromPath(params.folder+params.dataset)

process liftOver_Coordinates{
  input:
  val pattern from data
  file chain_file from Channel.fromPath(params.folder+'hg18ToHg19.over.chain')

  output:
  file('dataset2*') into resultChannel

  script:
  """
  name="$pattern"
  awk '{print "chr" \$1, \$4 -1, \$4, \$2 }' \${name%%.*}.bim | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > dataset.tolift
  liftOver dataset.tolift hg18ToHg19.over.chain dataset1 dataset_NCBI36.unMapped

  awk '{print \$4}' dataset1 > dataset1.snps
  plink --bfile \${name%%.*} --extract dataset1.snps --make-bed --out dataset1

  awk '{print \$4, \$3}' dataset1  > dataset1.pos
  plink --bfile dataset1 --update-map dataset1.pos --make-bed --out dataset2

  """
}


process plink_processing{
  publishDir params.out, mode: 'copy'

  input:
  file data from resultChannel.collect()
  file csv_file from Channel.fromPath(params.folder+'pone.0002551.s003.csv')
  file data2 from Channel.fromPath(params.folder+'HRC-1000G-check-bim.pl').collect()
  file data3 from Channel.fromPath(params.folder+'1000GP_Phase3_combined.legend').collect()

  output:
  file(params.dataset+"/") into resultChannel2

  script:
  if( params.dataset == 'hapmap_r23a' )
    """
    awk -F',' '{print \$1}' $csv_file | grep -v SNPID > AIM_list.txt
    plink --bfile dataset2 --extract AIM_list.txt --make-bed --out dataset3

    plink --freq --bfile dataset3 --make-bed --out dataset4
    perl HRC-1000G-check-bim.pl -b dataset4.bim -f dataset4.frq -r 1000GP_Phase3_combined.legend -g -p ALL -x -n
    bash Run-plink.sh

    awk '{print \$2}' dataset4-updated.fam > ref_samples.txt

    mkdir ${params.dataset}/
    mv ref_samples.txt ${params.dataset}/ref_samples.txt
    mv AIM_list.txt ${params.dataset}/AIM_list.txt
    mv dataset4-updated.bed ${params.dataset}/${params.dataset}.bed
    mv dataset4-updated.bim ${params.dataset}/${params.dataset}.bim
    mv dataset4-updated.fam ${params.dataset}/${params.dataset}.fam
    """

  else
    """
    plink --freq --bfile dataset2 --make-bed --out dataset3
    perl HRC-1000G-check-bim.pl -b dataset3.bim -f dataset3.frq -r 1000GP_Phase3_combined.legend -g -p ALL -x -n
    bash Run-plink.sh

    mkdir ${params.dataset}/
    mv dataset3-updated.bed ${params.dataset}/${params.dataset}.bed
    mv dataset3-updated.bim ${params.dataset}/${params.dataset}.bim
    mv dataset3-updated.fam ${params.dataset}/${params.dataset}.fam
    """
}
