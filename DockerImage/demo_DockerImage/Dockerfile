# Use base image with miniconda3 installed
FROM continuumio/miniconda3
MAINTAINER Yueyao Gao <tag@broadinstitute.org>

# Specify Workdir
WORKDIR /BaseImage

# Create the environment
COPY Griffin_env_explicit.yml .
RUN conda env create -f Griffin_env_explicit.yml

# Clone the git repository
RUN git clone https://github.com/adoebley/Griffin.git

# Set the working directory to the cloned repository
WORKDIR /BaseImage/Griffin

# DEMO Specifc Commands:
# Create run-demo directory and move snakemakes workflow scripts to run_demo
RUN mkdir run_demo
RUN cp -r snakemakes/griffin_GC_and_mappability_correction run_demo
RUN cp -r snakemakes/griffin_nucleosome_profiling run_demo

# Edit the config yaml files for demo run
COPY demo_samples.yaml ./run_demo/griffin_GC_and_mappability_correction/config/samples.yaml
COPY demo_sites.yaml ./run_demo/griffin_nucleosome_profiling/config/sites.yaml
COPY demo_samples.GC.yaml ./run_demo/griffin_nucleosome_profiling/config/samples.GC.yaml

# Run the config with mappabilty task enabled
COPY mappabilty_enabled_config.yaml ./run_demo/griffin_GC_and_mappability_correction/config/config.yaml

# Copy the demo bash script
COPY run-demo.sh .