# Set the base image
FROM lucapinello/crispresso
# Dockerfile author / maintainer 
MAINTAINER Gerry Meixiong <gerry.meixiong@czbiohub.org> 

RUN conda update -n base conda && conda install -c bioconda trimmomatic && conda install -c conda-forge awscli 
 

