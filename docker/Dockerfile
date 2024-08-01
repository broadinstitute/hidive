FROM python:3.9-slim AS builder

MAINTAINER Kiran V Garimella

ARG branch

# Install dependencies
RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
    curl keychain git-lfs time datamash wget bzip2 make gcc cmake \
    zlib1g-dev libssl-dev lbzip2 libbz2-dev libncurses5-dev \
    libncursesw5-dev liblzma-dev clang libclang-dev \
    libtcmalloc-minimal4 python3-dev python3-setuptools && \
    rm -rf /var/lib/apt/lists/*

# Install gsutil
RUN curl https://sdk.cloud.google.com | bash

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

# Final stage
FROM python:3.9-slim

# Copy necessary files from builder stage
COPY --from=builder /usr/local/bin/hidive /usr/local/bin/hidive
#COPY --from=builder /root/.cargo /root/.cargo
COPY --from=builder /usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4 /usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4

# Set environment variables
ENV PATH="/root/.cargo/bin:${PATH}"
ENV LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"

# We want this to be interactive
CMD ["/bin/bash"]