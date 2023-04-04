# Use base image with miniconda3 installed
FROM continuumio/miniconda3
MAINTAINER Yueyao Gao <tag@broadinstitute.org>

# Specify Workdir
WORKDIR /BaseImage

# Create the environment
COPY Griffin_env_explicit.yml .
RUN conda env create -f Griffin_env_explicit.yml

# Copy Python scripts from Griffin
RUN mkdir -p /BaseImage/Griffin/scripts
COPY scripts/*py /BaseImage/Griffin/scripts/

# Set the working directory to the Griffin
WORKDIR /BaseImage/Griffin