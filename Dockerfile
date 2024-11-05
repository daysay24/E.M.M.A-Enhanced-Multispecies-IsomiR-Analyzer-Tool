# Use an Ubuntu base image
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive

# Create a project folder in Docker container
WORKDIR /emma

# Copy source code from host to Docker container's project folder
COPY . /emma

# Install software
RUN apt-get update
RUN apt-get install -y python3
RUN apt-get install -y python3-pip
RUN apt-get install -y r-base 
RUN apt-get install bedtools

# Set python to point to python3 globally
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1

# Install Python packages
RUN pip3 install -r /emma/requirements.txt

# Install R packages
RUN R -e "install.packages(c('ggplot2', 'patchwork', 'dplyr', 'ggrepel', 'tidyr'), repos='http://cran.rstudio.com/')"
