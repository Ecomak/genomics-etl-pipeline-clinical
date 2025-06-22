# Use Ubuntu 20.04 as base
FROM ubuntu:20.04

# Avoid interactive prompts during package install
ENV DEBIAN_FRONTEND=noninteractive

# Install system packages and ART simulator
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    build-essential \
    zlib1g-dev \
    art-nextgen-simulation-tools \
    bwa samtools bcftools \
    python3 \
    python3-pip \
    && apt-get clean

# Set working directory
WORKDIR /app

# Copy everything into the container
COPY . /app

# Install Python dependencies (if any)
COPY docker/requirements.txt .
RUN pip3 install -r requirements.txt

# Default shell
CMD ["/bin/bash"]
