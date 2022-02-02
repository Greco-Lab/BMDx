FROM rocker/r-base:latest

RUN apt-get update && apt-get install -y \
	openjdk-8-jdk \
	libcurl4-openssl-dev \
	libxml2-dev \
	libv8-dev 

RUN ["R", "CMD", "javareconf"]

USER docker
COPY install_dependencies.R /home/docker/install_dependencies.R
RUN ["Rscript", "-e",  "install.packages(\"trend\")"]
RUN ["Rscript", "-e",  "install.packages(\"doParallel\")"]
RUN ["Rscript", "/home/docker/install_dependencies.R"]
RUN ["mkdir", "/home/docker/bmdx"]
COPY ./ /home/docker/bmdx/

WORKDIR /home/docker/bmdx
ENTRYPOINT ["Rscript", "/home/docker/bmdx/shiny_start.R"]
