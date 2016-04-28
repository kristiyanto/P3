# MSGF
# Base container
FROM bioconductor/release_proteomics

# Maintainer
MAINTAINER Daniel Kristiyanto, daniel.kristiyanto@pnnl.gov

# Main Directory
WORKDIR /root

# Install Java


RUN apt-get update

RUN apt-get --fix-missing -q -y install python3 build-essential python3-dev default-jdk unzip g++ libxml2-dev libcurl4-openssl-dev apt-utils libnetcdf-dev
RUN apt-get clean
ENV JAVA_HOME=/usr/lib/jvm/java-7-openjdk-amd64

# Add MSGF plus
ADD MSGFPlus.jar /root/

# Add Entrypoint Script
ADD entry.py /root/
ADD itraq.R /root/
ADD scquant.R /root/


# Run on Entrypoint
CMD python3 /root/entry.py