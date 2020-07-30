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
    log.info "  USAGE : nextflow run IARCbioinfo/Imputation-nf -r v1.0 -profile singularity --target TCGA --input ../data/ "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run IARCbioinfo/Imputation-nf [-r vX.X -profile singularity] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--target                      pattern                      pattern of the target dataset which do the link with the plink files (.bed/.bim./fam)"
    log.info ""
    log.info "Optional arguments:"
    log.info "--<OPTION>                      <TYPE>                      <DESCRIPTION>"
    log.info "--input                      FOLDER                      Folder where you can find your data to run the pipeline"
    log.info "--script                      FOLDER                      Folder where you can find the auxiliary scripts of the pipeline"
    log.info "--out                      FOLDER                      Folder where you want to put the results of the pipeline"
    log.info "--VCFref                      FOLDER                      Folder to use as VCF reference"
    log.info "--BCFref                      FOLDER                      Folder to use as BCF reference"
    log.info "--M3VCFref                      FOLDER                      Folder to use as M3VCF reference"
    log.info "--legend                      FILE                      File to use as .legend"
    log.info "--fasta                      FILE                      File to use as fasta reference"
    log.info "--chain                      FILE                      File to use as liftover conversion"
    log.info "--out                      FOLDER                      Folder to use for results"
    log.info "--geno1                      FLOAT                       Value for the first genotyping call rate plink option"
    log.info "--geno2                      FLOAT                     Value for the second genotyping call rate plink option"
    log.info "--maf                      FLOAT                      Minor Allele Frequencie thresold for the data filtering step"
    log.info "--pihat                      FLOAT                      PI_HAT thresold for the data filtering step"
    log.info "--hwe                      FLOAT                      Hardy-Weinberg Equilibrium thresold for the data filtering step"
    log.info "--conversion                      [hg38/hg18/hg19]                      Option to convert data to hg38 version of the genome. Choose 'hg18' to convert your data from hg18 to hg38 or 'hg19' to convert your data from hg19 to hg38. Standard value is 'hg38'."


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
params.target = null      //HGDP & TCGA
params.origin = "hapmap_r23a"
params.geno1 = 0.03 
params.geno2 = 0.03
params.maf = 0.01
params.pihat = 0.185
params.hwe = 1e-8 

// ## -- Path :
params.input = null //data/gep/MR_Signatures/work/Boris/protocol_min/data/
params.script = 'IARCbioinfo/Imputation-nf/bin/' 
params.targetDir = params.input+params.target+'/'

params.folder = params.input+'files/'
params.legend = params.folder+'ALL.chr_GRCh38.genotypes.20170504.legend'
legend_file = file( params.legend )
params.fasta = params.folder+'GRCh38_full_analysis_set_plus_decoy_hla.fa'
params.fasta_fai = params.folder+'GRCh38_full_analysis_set_plus_decoy_hla.fa.fai'

params.originDir = params.folder+params.origin+'/'

params.ref=params.folder+"ref/"
params.VCFref = params.ref+"vcf/*"
params.BCFref = params.ref+"bcf/"
params.M3VCFref = params.ref+"m3vcf/*"

params.conversion = "hg38"
params.chain = params.folder

params.out = params.input


process UpdateHG381{ 
  //input:
  //file data from Channel.fromPath(params.targetDir+'*').collect()
  //file data from Channel.fromPath(params.folder+'HRC-1000G-check-bim-NoReadKey.pl').collect()
  //file data from Channel.fromPath(params.folder+'params.legend').collect()
  //file data from Channel.fromPath(params.folder+'params.chain').collect()

  shell:
  '''
  head -n1 !{legend_file} >> ref_freq_withHeader.txt
  '''
  
}


process UpdateHG3813{ 
  input:
  //file data from Channel.fromPath(params.targetDir+'*').collect()
  file data from Channel.fromPath(params.folder+'HRC-1000G-check-bim-NoReadKey.pl').collect()
  //file data from Channel.fromPath(params.folder+'params.legend').collect()
  //file data from Channel.fromPath(params.folder+'params.chain').collect()

  shell:
  '''
  head -n1 !{legend_file} >> ref_freq_withHeader.txt
  '''
  
}