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
log.info "  Imputation 1.0 : Pipeline for the imputation of a target dataset against a reference "
log.info "--------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE : nextflow run IARCbioinfo/Imputation-nf -r v1.0 -profile singularity --target my_target --input data/ --output result/ "
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
    log.info "--output                      FOLDER                      Folder in which the pipeline outputs while be saved"
    log.info "--VCFref                      FOLDER                      Folder to use as VCF reference"
    log.info "--BCFref                      FOLDER                      Folder to use as BCF reference"
    log.info "--M3VCFref                      FOLDER                      Folder to use as M3VCF reference"
    log.info "--legend                      FILE                      File to use as .legend"
    log.info "--fasta                      FILE                      File to use as fasta reference"
    log.info "--chain                      FILE                      File to use as liftover conversion"
    log.info "--geno1                      FLOAT                       Value for the first genotyping call rate plink option"
    log.info "--geno2                      FLOAT                     Value for the second genotyping call rate plink option"
    log.info "--maf                      FLOAT                      Minor Allele Frequencie thresold for the data filtering step"
    log.info "--pihat                      FLOAT                      PI_HAT thresold for the data filtering step"
    log.info "--hwe                      FLOAT                      Hardy-Weinberg Equilibrium thresold for the data filtering step"
    log.info "--conversion                      [hg38/hg18/hg19]                      Option to convert data to hg38 version of the genome. Choose 'hg18' to convert your data from hg18 to hg38 or 'hg19' to convert your data from hg19 to hg38. Standard value is 'hg38'."
    log.info "--cloud                      [off/on]                      Option to run the imputation on cloud Michighan and/or TOPMed server. You have to give your Michighan and/or TOPMed token to run it (Add --token_Michighan and/or --token_TOPMed)."
    log.info "--token_Michighan                      FILE                      File where you can find your Michighan token as string"
    log.info "--token_TOPMed                      FILE                      File where you can find your TOPMed token as string"
    log.info "--QC_cloud                      FOLDER                      Folder where you can find the VCF file imputed from the Michighan or TOPMed server"

    log.info ""
    log.info "Flags:"
    log.info "--<FLAG>                                                    <DESCRIPTION>"
    log.info ""
    exit 0
} else {
    /* Software information */
    log.info "help:                               ${params.help}"
}

// -- Option :
params.target = null
params.origin = "hapmap_r23a"
params.geno1 = 0.03
params.geno2 = 0.03
params.maf = 0.01
params.pihat = 0.185
params.hwe = 1e-8

// -- Path :
params.input = null
params.output = null
params.targetDir = params.input+params.target+'/'

params.folder = params.input+'files/'
params.legend = params.folder+'PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz'
legend_file = file( params.legend )
params.fasta = params.folder+'GRCh38_full_analysis_set_plus_decoy_hla.fa'
params.fasta_fai = params.folder+'GRCh38_full_analysis_set_plus_decoy_hla.fa.fai'

params.originDir = params.folder+params.origin+'/'

params.ref=params.folder+"ref/"
params.VCFref = params.ref+"vcf/*"
params.BCFref = params.ref+"bcf/"
params.M3VCFref = params.ref+"m3vcf/"

params.conversion = "hg38"
params.chain = params.folder

params.cloud = "off"
params.token_Michighan = null
params.token_TOPMed = null
params.imputationbot_password = null


params.QC_cloud = null

// -- Pipeline :
//if(params.QC_cloud==null){
  process UpdateHG38{
    input:
    file data from Channel.fromPath(params.targetDir+'*').collect()

    output:
    file ('*-updated.{bed,bim,fam}') into TargetUpdate
    file ('*-updated.bim') into TargetQC2

    shell:
    '''
    ############################################################################################
    ## -- 0 : Update version of the tergat dataset : hg18/19 --> hg38
    if [ !{params.conversion} != "hg38" ] ; then
      awk '{print "chr" $1, $4 -1, $4, $2 }' !{params.target}.bim | sed 's/chr23/chrX/' | sed 's/chr24/chrY/' > dataset.tolift

      if [ !{params.conversion} == "hg18" ] ; then
        liftOver dataset.tolift !{params.chain}hg18ToHg38.over.chain dataset1 dataset_NCBI36.unMapped
      fi

      if [ !{params.conversion} == "hg19" ] ; then
        liftOver dataset.tolift !{params.chain}hg19ToHg38.over.chain dataset1 dataset_NCBI36.unMapped
      fi

      awk '{print $4}' dataset1 > dataset1.snps
      plink --bfile !{params.target} --extract dataset1.snps --make-bed --out dataset1

      mkdir old
      mv !{params.target}* old/

      awk '{print $4, $3}' dataset1  > dataset1.pos
      plink --bfile dataset1 --update-map dataset1.pos --make-bed --out !{params.target}
    fi

    plink --freq --bfile !{params.target} --allow-no-sex --make-bed --out dataset3
    perl !{baseDir}/bin/HRC-1000G-check-bim.pl -b !{params.target}.bim -f dataset3.frq -r !{params.legend} -h
    grep -v "real-ref-alleles" Run-plink.sh> Run-plink-update.sh
    bash Run-plink-update.sh
    '''}
  process Admixture{
    publishDir params.output+params.target+'/admixture/', mode: 'copy',
       saveAs: {filename ->
           if (filename == "target4.bim") null
           else if (filename == "target4.bed") null
           else if (filename == "target4.fam") null
           else filename
       }

    input:
    file data from Channel.fromPath(params.folder+'relationships_w_pops_121708.txt').collect()
    file data from Channel.fromPath(params.originDir+"*").collect()
    file data from TargetUpdate.collect()

    output:
    file ('target4.{bed,bim,fam}') into Merge
    file ('target4.{bed,bim,fam}') into Merge2
    file ('out_pop_admixture/') into Admixture
    file ('out_pop_admixture/') into Admixture2
    file ('admixture_results_withGroups.txt') into Target

    shell:
    '''
    ## -- 1 : Retrieve common SNPs between SNPs in the reference list and the target data.
    awk '{print $2}' !{params.origin}.bim | sort > ref_SNPs.txt
    grep -Fwf <(cat ref_SNPs.txt) <(awk '{print $2}' !{params.target}-updated.bim)  > target_common_SNPs.txt

    ## -- 2 : Extract the genotypes associated to these common SNPs + merge of the dataset
    plink --bfile !{params.origin} --extract target_common_SNPs.txt --make-bed --out ref1
    plink --bfile !{params.target}-updated --extract target_common_SNPs.txt --make-bed --out target1
    plink --bfile target1 --bmerge ref1 --allow-no-sex --make-bed --out merge1


    if [ -e merge1-merge.missnp ]
    then
      plink --bfile ref1 --exclude merge1-merge.missnp --make-bed --out ref2
      plink --bfile target1 --exclude merge1-merge.missnp --make-bed --out target2
      plink --bfile ref2 --bmerge target2 --allow-no-sex --make-bed --out merge
    else
      mv merge1.fam merge.fam
      mv merge1.bed merge.bed
      mv merge1.bim merge.bim
    fi

    # perform ld pruning here!!! look at ELN command
    #plink --noweb --bfile merge --r2 --out merge_pruned --ld-window-kb 10000 --ld-window 10000 --ld-window-r2 0.1
    plink --bfile merge --out merge_pruned --indep-pairwise 50 5 0.2
    plink --bfile merge --extract merge_pruned.prune.in --make-bed --out merge_pruned
    awk '{print $2}' merge_pruned.fam  >  all_samples.txt

    ## -- 3 : Associate each reference sample to the correct origin + run admixure
    sort -k2 relationships_w_pops_121708.txt > relationships_w_pops_121708_2.txt
    awk '{print $2}' !{params.origin}.fam | sort -k1,1 > ref_ind.txt
    join -11 -22 ref_ind.txt relationships_w_pops_121708_2.txt | awk '{print $1, $7}' | awk '{if( $2 == "CEU") print $0,1,1; else { if($2 == "YRI") print $0,2,1; else {print $0,3,1}}}' > ind_pop.txt
    Rscript !{baseDir}/bin/create_pop_file.r merge_pruned.fam ind_pop.txt merge_pruned.pop
    K=3
    admixture --cv merge_pruned.bed $K --supervised -j4 | tee log${K}.out
    Rscript !{baseDir}/bin/process_admixture.r !{params.target}-updated

    ############################################################################################
    ## -- 4 : First filtering step
    plink --bfile !{params.target}-updated --geno !{params.geno1} --make-bed --out target3
    plink --bfile target3 --maf !{params.maf} --make-bed --out target4
    '''
  }
  process Filtering1{
    input:
    file data from Merge.collect()
    file data from Admixture.collect()
    val rspop from Channel.from("ALL")

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
    '''
  }

  process QC1{
    publishDir params.output+params.target+'/QC1/', mode: 'copy'

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
    '''
  }
  process Filtering2{
    input:
    file data from TargetFilter.collect()
    val rspop from Channel.from("ALL")

    output:
    file ('withFreqFiltering_*') into DirFiltering
    file ('target_hwe_*.bim') into HWresult
    file ('target_freq_*.frq') into FreqResult
    file ('ID-target_*-HRC.txt') into FreqResultId

    shell:
    '''
    pop=!{rspop}

    ## -- 10 : AF based filter
    plink --freq --bfile target_${pop} --output-chr chr26 --out target_freq_${pop}
    perl !{baseDir}/bin/HRC-1000G-check-bim.pl -b target_${pop}.bim -f target_freq_${pop}.frq -r !{params.legend} -h
    mkdir withFreqFiltering_${pop}
    cp *-HRC* Run-plink.sh withFreqFiltering_${pop}

    ## -- 11 :  HWE filtering
    plink --bfile target_${pop} --geno !{params.geno2} --make-bed --out target_geno_${pop}
    plink --bfile target_geno_${pop} --hwe !{params.hwe} --make-bed --out target_hwe_${pop}
    '''
  }
  process Make_SNP_Filtering{
    publishDir params.output+params.target+'/SNP_filtering/', mode: 'copy'
    input:
    file data from DirFiltering.collect()
    file data from HWresult.collect()
    file data from Bim.collect()

    output:
    file ('filtered_snps.txt') into SNPsFilter
    file ('filtered_snps.txt') into SNPsFilter2

    shell:
    '''
    ## -- 12 : Create list of SNPs to filter
    Rscript !{baseDir}/bin/afterGenotyping_SNPs_filtering.r
    '''
  }
  process Filtering3{
    input:
    file data from Merge2.collect()
    file data from SNPsFilter.collect()

    output:
    file ('target5-updated-chr*') into TargetChr

    shell:
    '''
    ## -- 13 : Removing SNPs
    plink --bfile target4 --exclude filtered_snps.txt  --make-bed --out target5

    ## -- 14 : filter
    plink --freq --bfile target5  --out target6
    perl !{baseDir}/bin/HRC-1000G-check-bim.pl -b target5.bim -f target6.frq -r !{params.legend} -h
    bash Run-plink.sh
    '''
  }

  process QC2{
    publishDir params.output+params.target+'/QC2/', mode: 'copy'

    input:
    file data from SNPsFilter2.collect()
    file data from FreqResult.collect()
    file data from FreqResultId.collect()
    val rspop from Channel.from("ALL")
    file data from TargetQC2.collect()

    output:
    file ('*.pdf') into FigureQC2

    shell:
    '''
    zmore !{legend_file} | head -n1 > ref_freq_withHeader.txt
    grep -Fwf <(awk '{print $2}' !{params.target}-updated.bim) <(zcat !{legend_file}) > ref_freq.txt
    cat ref_freq.txt >> ref_freq_withHeader.txt

    pop=!{rspop}

    pop2="AF"


    ## -- 15 : Figure QC2
    Rscript !{baseDir}/bin/preImputation_QC_plots.r ${pop} ${pop2}
    '''
  }

  process Filtering4{
    input:
    file data from TargetChr.collect()
    file data from Channel.fromPath(params.folder+'checkVCF.py').collect()
    file data from Channel.fromPath(params.fasta+'*').collect()
    val chromosome from 21..22

    output:
    file ('*-REFfixed.vcf.gz') into FilterFinal
    file ('*-REFfixed.vcf.gz') into FilterFinal2
    file ('*-REFfixed.vcf.gz') into FilterFinal3

    shell:
    '''
    chr=!{chromosome}
    export TMPDIR=/tmp/

    ## -- 16 : Remove ambiguous strand/unknown SNPs
    awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' target5-updated-chr${chr}.bim | grep -v ambig | grep -v -e --- | sort -u > NonAmbiguous${chr}.snplist.txt
    plink --bfile target5-updated-chr${chr} --extract NonAmbiguous${chr}.snplist.txt --output-chr chr26 --make-bed --out target6_chr${chr}

    ## -- 17 : Create VCF
    plink --bfile target6_chr${chr} --output-chr chr26 --recode vcf --out target6_chr${chr}_vcf
    bcftools sort target6_chr${chr}_vcf.vcf | bgzip -c  > chr${chr}.vcf.gz

    ## -- 18 : Check SNPs
    python2 checkVCF.py -r !{params.fasta} -o after_check_${chr} chr${chr}.vcf.gz
    bcftools norm --check-ref ws -f !{params.fasta} chr${chr}.vcf.gz | bcftools view -m 2 -M 2  | bgzip -c > chr${chr}-REFfixed.vcf.gz
    python2 checkVCF.py -r !{params.fasta} -o after_check2_${chr} chr${chr}-REFfixed.vcf.gz
    '''
  }

//  }
if(params.cloud=="on"){
  if (params.token_Michighan){
    process Michighan_Imputation{
      publishDir params.output+params.target+'/Result_Imputation/', mode: 'move'
      input:
      file data from Channel.fromPath(params.imputationbot_password).collect()
      file data from FilterFinal2.collect()

      output:
      file 'job-*/local/*dose.vcf.gz' into Imputation_dose
      // file '*info.gz' into Imputation_info
      // file 'stat/' into Imputation_stat
      // file 'log/' into Imputation_log

      when:
      !params.token_Michighan!=null

      shell:
      '''
      pw=$(cat !{params.imputationbot_password})
      token=$(cat !{params.token_Michighan})
      imputationbot add-instance https://imputationserver.sph.umich.edu ${token}
      imputationbot impute --files chr*-REFfixed.vcf.gz --refpanel hrc-r1.1 --build hg38 --autoDownload --password ${pw} --population mixed
      '''
    }
  }
  if (params.token_TOPMed){
    process TOPMed_Imputation{
      publishDir params.output+params.target+'/Result_Imputation/', mode: 'move'
      input:
      file data from Channel.fromPath(params.imputationbot_password).collect()
      file data from FilterFinal2.collect()

      when:
      !params.token_TOPMed!=null

      output:
      file ('*.{dose.vcf,info}.gz') into Imputation_res

      shell:
      '''
      pw=$(cat !{params.imputationbot_password})
      token=$(cat !{params.token_TOPMed})
      imputationbot add-instance https://imputation.biodatacatalyst.nhlbi.nih.gov ${token}
      imputationbot impute --files chr*-REFfixed.vcf.gz --refpanel topmed-r2 --build hg38 --autoDownload --password ${pw} --population mixed
      mv job-*/local/*.gz .
      mv job-*/logfile/* .
      '''
    }
  }

  process QC3_sh{
    input:
    val population from('ALL')
    each chromosome from 21..22
    file data from Channel.fromPath(params.output+params.target+'/Result_Imputation/*dose.vcf.gz').collect()
    file data from Channel.fromPath(params.folder+'bravo-dbsnp-all.vcf.gz').collect()
    file data from Channel.fromPath(params.folder+'bravo-dbsnp-all.vcf.gz.csi').collect()
    file data from Admixture2.collect()

    output:
    file '*.{txt,frq}' into PostImputation_QC_sh_result

    shell:
    '''
    pop=!{population}
    chr=!{chromosome}

    ## -- 24 : QC3
    bash !{baseDir}/bin/postImputation_QC.sh ${chr} ${pop}
    '''
  }

  process QC3_R{
    if(params.QC_cloud==null){publishDir params.output+params.target+'/QC3/', mode: 'copy'}
    else{publishDir params.output+params.target+'/QC3_cloud/', mode: 'copy'}
    cpus 4

    input:
    val population from('ALL')
    file data from PostImputation_QC_sh_result.collect()

    output:
    file '*v2.{txt,pdf}' into PostImputation_QC_R_result
    file '*CHR.pdf' into PostImputation_QC_R_result2

    shell:
    '''
    pop=!{population}
    Rscript !{baseDir}/bin/postImputation_QC_plots.r ${pop} 0.3
    '''
  }
}
