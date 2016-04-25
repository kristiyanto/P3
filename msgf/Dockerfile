# MSGF
# Base container
FROM python:3.5

# Maintainer
MAINTAINER Daniel Kristiyanto, daniel.kristiyanto@pnnl.gov

# Main Directory
WORKDIR /root

# Install Java
RUN apt-get update
RUN apt-get -q -y install default-jdk unzip 
RUN apt-get clean
ENV JAVA_HOME=/usr/lib/jvm/java-7-openjdk-amd64

# Add MSGF plus
ADD MSGFPlus.jar /root/

# Add Entrypoint Script
ADD entry.py /root/
ADD _functions.py /root/


# Run on Entrypoint
CMD python /root/entry.py