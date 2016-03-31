# Base container
FROM bioconductor/devel_proteomics

# Maintainer
MAINTAINER Daniel Kristiyanto, daniel.kristiyanto@pnnl.gov

# Main Directory
WORKDIR /root

# Install Java
RUN apt-get update
RUN apt-get -q -y install default-jdk unzip python2.7 python2.7-dev g++ build-essential libxml2-dev libcurl4-openssl-dev apt-utils libnetcdf-dev
RUN apt-get clean
ENV JAVA_HOME=/usr/lib/jvm/java-7-openjdk-amd64

# Add MSGF plus
#ADD https://omics.pnl.gov/sites/default/files/software/MSGFPlus.20140716.zip /root/MSGF/
#RUN cd /root/MSGF/ && \
#	unzip -q MSGFPlus.20140716.zip && \
#   rm MSGFPlus.20140716.zip 
ADD MSGFPlus.jar /root/

# Add Entrypoint Script
ADD entry.py /root/
ADD _functions.py /root/
ADD labelfreequant.R /root/
ADD install.R /root/

RUN Rscript install.R

# Run on Entrypoint
CMD python /root/entry.py