################## BASE IMAGE #####################
FROM continuumio/miniconda3:latest

################## METADATA #######################
LABEL base_image="continuumio/miniconda3"
LABEL version="4.8.3"
LABEL software="Imputation_genotypage"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for Imputation_genotypage"
LABEL about.home="http://github.com/IARCbioinfo/Imputation_genotypage"
LABEL about.documentation="http://github.com/IARCbioinfo/Imputation_genotypage/README.md"
LABEL about.license_file="http://github.com/IARCbioinfo/Imputation_genotypage/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -n Imputation_genotypage -f /environment.yml && conda clean -a

RUN apt-get install -y cmake python-pip python-dev
RUN pip install cget 
RUN git clone https://github.com/statgen/Minimac4.git ; cd Minimac4 ; bash install.sh

RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz
RUN tar xvzf Eagle_v2.4.1.tar.gz

ENV PATH /opt/conda/envs/Imputation_genotypage/bin:$PATH
ENV PATH="$PATH:/Eagle_v2.4.1"
ENV PATH="$PATH:/Minimac4/release-build" 