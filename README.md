# Imputation genotypage : Pipeline V1.0
## A nextflow pipeline to realise a dataset's imputation genotypage

[![CircleCI](https://circleci.com/gh/IARCbioinfo/template-nf.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/Imputation-nf)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/imputation-nf/)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1404)
[![DOI](https://zenodo.org/badge/94193130.svg)](https://zenodo.org/badge/latestdoi/94193130)

![Workflow representation](template-nf.png)

## Description
The pipeline used to perform the imputation of several targets datasets processed with standard input.

Here is a summary of the method :
- Preprocessing of target dataset if necessary : by using the nextflow script Preparation_dataset.nf which take the dataset plink file (.bed/.bim/.fam) in input.
- First step : Origin estimation of sample from the target dataset by using admixture tools and the hapmap dataset as reference.
- Second step : Series of SNPs filters and quality checking from the target dataset before the imputation step.
- Third step : VCF production
- Last step : Phasing and imputation


See the Usage section to test the full pipeline with the hapmap dataset like admixture's reference and HGDP dataset as target.

## Dependencies
The pipeline works under Linux distributions.

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.

2. External software:
- LiftOver : conda install ucsc-liftover
- Plink (PLINK v1.90b6.12 64-bit (28 Oct 2019)) : conda install plink
- Admixture (ADMIXTURE Version 1.3.0) : conda install admixture
- Perl : conda install perl
- Term::ReadKey module : conda install perl-termreadkey
- BcfTools : conda install bcftools
- eagle 2.4.1 : [See instructions](https://data.broadinstitute.org/alkesgroup/Eagle/#x1-50002.2)
- minimac4 : conda install cmake ; pip install cget ; git clone https://github.com/statgen/Minimac4.git ; cd Minimac4 ; bash install.sh

3. File to download :
- [Hapmap Dataset](zzz.bwh.harvard.edu/plink/dist/hapmap_r23a.zip) : as reference's dataset for admixture
- [HGDP Dataset](http://www.hagsc.org/hgdp/data/hgdp.zip) : for the dataset's test, you have to use the toMap.py & toPed.py in the 'converstion' directory to convert files in the .map/.ped plink format. Next you have to convert this last output in the .bed/.bam/.fam plink format by using plink line command and run the imputation's pipeline.
- Perl tool : [HRC-1000G-check-bim-NoReadKey.pl](https://www.well.ox.ac.uk/~wrayner/tools/) & [1000GP_Phase3_combined.legend](https://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz)
- LiftOver tool : [hg19ToHg38.over.chain](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz) & [hg18ToHg38.over.chain](http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg38.over.chain.gz)
- Peparation dataset tool : [pone.0002551.s003.xls](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2432498/bin/pone.0002551.s003.xls) (Convert it in .csv format)
- Admixture tool : [relationships_w_pops_121708.txt](ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/relationships_w_pops_121708.txt)
- [CheckVCF](https://github.com/zhanxw/checkVCF/raw/master/checkVCF.py), [Fasta file in V37](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz) & [Fasta file in V38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/)
- [1000G Reference in Hg38](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/) with the [doc](https://data.broadinstitute.org/alkesgroup/Eagle/#x1-320005.3.2)

4. Other to know : 
- Create [legend](https://imputationserver.readthedocs.io/en/latest/create-reference-panels/#create-legend-files) & [m3vcf](https://imputationserver.readthedocs.io/en/latest/create-reference-panels/#create-m3vcf-files) files for the reference


You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.




## Input
  | Type      | Description     |
  |-----------|---------------|
  | Plink datasets | Corresponds to the admixture reference's dataset and the target dataset to be analysed, both already preprocessed by filters found in the Preparation_dataset.nf. Composed by the following files : bed, bim & fam |



  Specify the test files location






## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --target | HGDP | Pattern of the target dataset which do the link with the file .bed/.bim./fam for plink. The data have to be placed in the main directory and in a directory which have the same name that the parttern name of the dataset. (for example : 'user/main_data/HGDP/')|
|--input| user/main_data/ | The path of the main directory. WARNING : for now, you have to create a 'user/main_data/data_start' directory where we can find all the data to download before to run the script. |
|--script | my/directory/script/bin | The path of the bin script directory, to be able to run the annexe programme grom the pipeline |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --geno1   |          0.03 | First genotyping call rate plink threshold, apply in the full target dataset |
| --geno2   |          0.03 | Second genotyping call rate plink threshold, apply in the target dataset divide by population |
| --maf     |          0.01 | Minor allele frequencies plink threshold, apply in the full target dataset |
| --pihat   |         0.185 | Minimum pi_hat value use for the relatedness test, 0.185 is halfway between the expected IBD for third- and second-degree relatives |
| --hwe     |          1e-8 | Hardy-Weinberg Equilibrium plink p-value threshold |


  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |
| --flag2    |      .... |


## Usage
1. Download all the necessary file to run the full pipeline (See the dependencies section, part 3) and put all the file on /your/directory/path/data_start/.

2. Prepare the dataset to the imputation pipeline.

Here we take the HGDP dataset as target and the Hapmap dataset as reference to run admixture tools. To make sure that the data is available to use, you have to convert the .txt's HGDP file in the good plink format. To do this, use the script toMap.py & toPed.py that you can find in the directory bin/ (Make sure that more of 30gb of memory are free before to run the toPed.py script). Give in argument the path of your directory where the script can find files HGDP_Map.txt (for the toMap's script) and HGDP_FinalReport_Forward.txt (or the toPed's script)

  ```
  python3 toMap.py /directory/path/
  python3 toPed.py /directory/path/
  ```

After you have to convert the .ped/.map file in .bed/.bim/.bam file with the use of plink line and copy it in a directory.

  ```
  plink --map /your/directory/path/HGDP.map --ped /directory/path/HGDP.ped --make-bed --out /your/directory/path/HGDP
  ```

And then we can run the nextflow script to prepare the both dataset. Make sur that the files 'pone.0002551.s003.csv', 'HRC-1000G-check-bim.pl', '1000GP_Phase3_combined.legend' and 'hg18ToHg19.over.chain' can be found in /your/directory/path/

  ```
  nextflow run Preparation_dataset.nf --folder /directory/path/ --dataset HGDP
  nextflow run Preparation_dataset.nf --folder /directory/path/ --dataset hapmap_r23a
  ```

3. Run the imputation pipeline.


  ```
  nextflow run IARCbioinfo/Imputation_genotypage.nf -r v1.0 -profile singularity --target HGDP --input /directory/path/ --script /directory/path/script/bin/ 
  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | output1    | ...... |
  | output2    | ...... |


## Detailed description (optional section)
...

## Directed Acyclic Graph
[![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/Imputation-nf/blob/master/dag.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | contrib1*    |            xx | Developer to contact for support (link to specific gitter chatroom) |
  | contrib2    |            xx | Developer |
  | contrib3    |            xx | Tester |

## References (optional)

## FAQ (optional)
# test-pipeline
