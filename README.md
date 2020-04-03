# Imputation genotypage : Pipeline
## A nextflow pipeline to realise a dataset's imputation genotypage

[![CircleCI](https://circleci.com/gh/IARCbioinfo/template-nf.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/template-nf)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/template-nf/)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1404)
[![DOI](https://zenodo.org/badge/94193130.svg)](https://zenodo.org/badge/latestdoi/94193130)

![Workflow representation](template-nf.png)

## Description
The pipeline used to perform the imputation of several targets datasets processed with standard input with the use of a reference dataset (here, hapmap dataset).

Here is a summary of the method:
- Preprocessing of reference/target dataset if necessary : filtration using by liftOver, plink and perl command line
- First filtration's step of dataset with bash and plink command lines to keep only the SNPs of interest
- Merge of both datasets
- Run admixture tools to associate each reference sample to the correct origin
- Second filtration's step from the start target dataset to keep only the samples by ancerstry population and to perform a sex checking by population
-

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

3. File to download :
- [Hapmap Dataset](zzz.bwh.harvard.edu/plink/dist/hapmap_r23a.zip) : as reference's dataset
- [HGDP Dataset](http://www.hagsc.org/hgdp/data/hgdp.zip) : for the dataset's test, you have to use the toMap.py & toPed.py in the 'converstion' directory to convert files in the .map/.ped plink format. Next you have to convert this last output in the .bed/.bam/.fam plink format by using plink line command and run the imputation's pipeline.
- Perl tools : [1](https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.11.zip) & [2](https://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz)
- [LiftOver tools](http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain.gz)
- Pipeline tools : [1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2432498/bin/pone.0002551.s003.xls) (Convert it in .csv format) & [2](ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2009-01_phaseIII/plink_format/relationships_w_pops_121708.txt)
- [CheckVCF](https://github.com/zhanxw/checkVCF/raw/master/checkVCF.py)
- [HumanV37](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz)


You can avoid installing all the external software by only installing Docker. See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.


## Input
  | Type      | Description     |
  |-----------|---------------|
  | Reference/Target Plink dataset | Corresponds to the reference's dataset and the dataset to be analyzed, both already preprocessed by filters found in the Preparation_dataset.nf. Composed by the following files : bed, bim & fam |



  Specify the test files location

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
| --ref | hapmap_r23a | Pattern of the reference dataset which do the link with the file .bed/.bim./fam for plink |
| --target | HGDP | Pattern of the target dataset which do the link with the file .bed/.bim./fam for plink |

  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
| --param3   |            xx | ...... |
| --param4    |            xx | ...... |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
| --help    | Display help |
| --flag2    |      .... |


## Usage
  ```
  nextflow run Preparation_dataset.nf --dataset HGDP
  nextflow run Preparation_dataset.nf --dataset hapmap_r23a

  nextflow run Imputation_genotypage.nf --ref hapmap_r23a --target HGDP

  ```

## Output
  | Type      | Description     |
  |-----------|---------------|
  | output1    | ...... |
  | output2    | ...... |


## Detailed description (optional section)
...

## Directed Acyclic Graph
[![DAG](dag.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/template-nf/blob/master/dag.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | contrib1*    |            xx | Developer to contact for support (link to specific gitter chatroom) |
  | contrib2    |            xx | Developer |
  | contrib3    |            xx | Tester |

## References (optional)

## FAQ (optional)
# test-pipeline
