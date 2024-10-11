FROM python:3.9-slim AS builder
#FROM gurobi/optimizer:10.0.3 AS builder

MAINTAINER Kiran V Garimella

ARG branch

# Install dependencies
RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
    curl keychain git-lfs time datamash wget bzip2 make gcc cmake g++ build-essential\
    zlib1g-dev libssl-dev lbzip2 libbz2-dev libncurses5-dev libstdc++6\
    libncursesw5-dev liblzma-dev clang libclang-dev pkg-config \
    libtcmalloc-minimal4 python3-dev python3-setuptools && \
    rm -rf /var/lib/apt/lists/*

# Install gcloud
#RUN curl https://sdk.cloud.google.com | bash

# Reduced gcloud installation (from framegrace: https://github.com/GoogleCloudPlatform/gsutil/issues/1732#issuecomment-2029591598)
ENV CLOUD_SDK_VERSION="496.0.0"
ENV ARCH="x86_64"

RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${CLOUD_SDK_VERSION}-linux-${ARCH}.tar.gz && \
    tar xzf google-cloud-cli-${CLOUD_SDK_VERSION}-linux-${ARCH}.tar.gz && \
    rm google-cloud-cli-${CLOUD_SDK_VERSION}-linux-${ARCH}.tar.gz && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    gcloud components remove bq app-engine-go appctl package-go-module \
        cbt bigtable cloud-datastore-emulator cloud-firestore-emulator \
        pubsub-emulator cloud-run-proxy cloud-sql-proxy kustomize \
        log-streaming managed-flink-client minikube nomos local-extract \
        skaffold terraform-tools anthos-auth config-connector enterprise-certificate-proxy \
        alpha beta app-engine-java app-engine-python app-engine-python-extras \
        gke-gcloud-auth-plugin istioctl kpt kubectl kubectl-oidc pkg && \
    gcloud components update && \
    rm -rf $(find google-cloud-sdk/ -regex ".*/__pycache__") && \
    rm -rf google-cloud-sdk/.install/.backup && \
    gcloud --version

RUN git config --system credential.'https://source.developers.google.com'.helper gcloud.sh

# Setup crcmodc for gsutil
RUN pip3 uninstall -y crcmod && \
    pip3 install --no-cache-dir -U crcmod

# Fix python malloc bug
RUN echo 'export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"' >> ~/.bashrc

# Create known hosts for github ssh info
RUN mkdir -p ~/.ssh && \
    touch ~/.ssh/known_hosts && \
    ssh-keyscan github.com >> ~/.ssh/known_hosts

# Install cargo
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install hidive
RUN git clone https://github.com/broadinstitute/hidive.git && \
    cd hidive && \
    git checkout $branch && \
    cargo build --release && \
    cp target/release/hidive /usr/local/bin/ && \
    cd .. && \
    rm -rf hidive

# Install minimap2
RUN git clone https://github.com/lh3/minimap2 && \
    cd minimap2 && \
    make && \
    cp minimap2 /usr/local/bin/ && \
    cd .. && \
    rm -rf minimap2

# Final stage
FROM python:3.9-slim

# Copy necessary files from builder stage
COPY --from=builder /usr/local/bin/hidive /usr/local/bin/hidive
COPY --from=builder /usr/local/bin/minimap2 /usr/local/bin/minimap2
COPY --from=builder /usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4 /usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4

# Copy gcloud utilities from builder stage
COPY --from=builder /root/google-cloud-sdk /root/google-cloud-sdk

# Set environment variables
ENV PATH="/root/.cargo/bin:/root/google-cloud-sdk/bin:${PATH}"
ENV LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"

# We want this to be interactive
CMD ["/bin/bash"]