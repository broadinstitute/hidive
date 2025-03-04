FROM python:3.9-slim AS builder
#FROM gurobi/optimizer:10.0.3 AS builder

MAINTAINER Kiran V Garimella

ARG branch

# Install dependencies
RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
    curl wget bzip2 make gcc cmake g++ keychain git build-essential zlib1g-dev libssl-dev lbzip2 \
    libbz2-dev libcurl4-gnutls-dev libncurses5-dev libstdc++6 libncursesw5-dev liblzma-dev clang \
    libclang-dev pkg-config libtcmalloc-minimal4 python3-dev python3-setuptools fontconfig && \
    rm -rf /var/lib/apt/lists/*

# Reduced gcloud installation (from framegrace: https://github.com/GoogleCloudPlatform/gsutil/issues/1732#issuecomment-2029591598)
ENV PATH /google-cloud-sdk/bin:$PATH
ENV CLOUD_SDK_VERSION=496.0.0
ENV ARCH=x86_64

RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${CLOUD_SDK_VERSION}-linux-${ARCH}.tar.gz && \
    tar xzf google-cloud-cli-${CLOUD_SDK_VERSION}-linux-${ARCH}.tar.gz && \
    rm google-cloud-cli-${CLOUD_SDK_VERSION}-linux-${ARCH}.tar.gz && \
    rm -rf /google-cloud-sdk/platform/bundledpythonunix && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    gcloud components remove -q bq && \
    rm -rf $(find google-cloud-sdk/ -regex ".*/__pycache__") && \
    rm -rf google-cloud-sdk/.install/.backup && \
    rm -rf google-cloud-sdk/bin/anthoscli && \
    gcloud --version

# Setup crcmodc for gsutil
RUN pip3 uninstall -y crcmod && \
    pip3 install --no-cache-dir -U crcmod

# Fix python malloc bug
RUN echo 'export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"' >> ~/.bashrc

# Install minimap2
RUN wget -O minimap2-2.28.tar.bz2 https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28.tar.bz2 && \
    tar -vxjf minimap2-2.28.tar.bz2 && \
    cd minimap2-2.28 && \
    make && \
    cp minimap2 /usr/local/bin/ && \
    cd .. && \
    rm -rf minimap2-2.28 minimap2-2.28.tar.bz2

# Install samtools
RUN wget -O samtools-1.21.tar.bz2 https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -vxjf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    make && \
    cp samtools /usr/local/bin/ && \
    cd .. && \
    rm -rf samtools-1.21 samtools-1.21.tar.bz2

# Install bcftools
RUN wget -O bcftools-1.21.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 && \
    tar -vxjf bcftools-1.21.tar.bz2 && \
    cd bcftools-1.21 && \
    make && \
    cp bcftools /usr/local/bin/ && \
    cd .. && \
    rm -rf bcftools-1.21 bcftools-1.21.tar.bz2

RUN wget -O htslib-1.21.tar.bz2 https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 && \
    tar -vxjf htslib-1.21.tar.bz2 && \
    cd htslib-1.21 && \
    make && \
    cp tabix bgzip /usr/local/bin/ && \
    cd .. && \
    rm -rf htslib-1.21 htslib-1.21.tar.bz2

# Final stage
FROM python:3.9-slim

# Install some libraries in final image
RUN apt-get update && \
    apt-get install -y --no-install-recommends fontconfig libcurl4-gnutls-dev \
    libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev liblzma-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy necessary files from builder stage
COPY --from=builder /usr/local/bin/minimap2 /usr/local/bin/minimap2
COPY --from=builder /usr/local/bin/samtools /usr/local/bin/samtools
COPY --from=builder /usr/local/bin/bcftools /usr/local/bin/bcftools
COPY --from=builder /usr/local/bin/tabix /usr/local/bin/tabix
COPY --from=builder /usr/local/bin/bgzip /usr/local/bin/bgzip
COPY --from=builder /usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4 /usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4

# Copy gcloud utilities from builder stage
COPY --from=builder /google-cloud-sdk /google-cloud-sdk

# Copy scripts from local directory
RUN mkdir -p /scripts
COPY *.py /scripts/
RUN pip3 install yattag

# Set environment variables
ENV PATH="/root/.cargo/bin:/google-cloud-sdk/bin:${PATH}"
ENV LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"

# We want this to be interactive
CMD ["/bin/bash"]

# Copy hidive binary from Github action artifact
COPY --from=binary /hidive /usr/local/bin/hidive
