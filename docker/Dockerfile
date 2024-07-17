FROM python:3.9.19-slim-bullseye

MAINTAINER Kiran V Garimella

ARG branch

# install gsutil
RUN apt-get update --allow-releaseinfo-change
RUN apt-get update -m
RUN apt-get install -y curl git-lfs time datamash wget bzip2 zlib1g-dev libssl-dev lbzip2 make gcc cmake libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev clang libclang-dev
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# Fix python malloc bug
# Detailed here: https://github.com/tensorflow/tensorflow/issues/6968 
RUN apt-get install -y libtcmalloc-minimal4
RUN echo 'export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"' >> ~/.bashrc

# Create known hosts for github ssh info:
RUN mkdir -p ~/.ssh && touch ~/.ssh/known_hosts
RUN ssh-keyscan github.com >> ~/.ssh/known_hosts

# Install vim so we can edit files later:
RUN apt-get install -y vim

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 \
    && tar xjf samtools-1.20.tar.bz2 \
    && rm samtools-1.20.tar.bz2 \
    && cd samtools-1.20 \
    && ./configure \
    && make \
    && make install

# install abpoa
RUN wget https://github.com/yangao07/abPOA/releases/download/v1.5.1/abPOA-v1.5.1.tar.gz \
    && tar -zxvf abPOA-v1.5.1.tar.gz \
    && cd abPOA-v1.5.1 \
    && make \
    && cp bin/abpoa /usr/local/bin/

# install cargo
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# install hidive
RUN git clone https://github.com/broadinstitute/hidive.git
RUN cd hidive \
    && git checkout $branch \
    && cargo build --release --verbose \
    && cp target/release/hidive /usr/local/bin/

# We want this to be interactive:
CMD ["/bin/bash"]
