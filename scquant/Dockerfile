# Spectrum Count Quantification
# Docker Container: scquant
# Base container
FROM bioconductor/release_proteomics

# Maintainer
MAINTAINER Daniel Kristiyanto, daniel.kristiyanto@pnnl.gov

# Main Directory
WORKDIR /root

# Install Java
RUN apt-get update
RUN apt-get -q -y install python2.7 g++ build-essential libxml2-dev libcurl4-openssl-dev apt-utils libnetcdf-dev
RUN apt-get clean

# Add Entrypoint Script
ADD entry.py /root/
ADD _functions.py /root/
ADD labelfreequant.R /root/
ADD install.R /root/

RUN Rscript install.R

# Run on Entrypoint
CMD python /root/entry.py