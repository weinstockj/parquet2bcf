# Base image
FROM ubuntu:22.04
ARG DUCKDB_VERSION=v1.1.3

ARG DUCKDB_ARCH=amd64

# Update package lists
RUN apt-get update && apt-get upgrade -y

# Install build dependencies for Rust + parallel
RUN apt-get install -y \
    curl \
    git \
    build-essential \
    pkg-config \
    parallel \
    zip unzip \
    libssl-dev 

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# Set PATH for Rust
ENV PATH="$PATH:$HOME/.cargo/bin"
ENV PATH="/root/.cargo/bin:${PATH}"
ENV PATH="/root/.cargo/bin:${PATH}"

# Install Bcftools
RUN apt-get install -y bcftools

RUN curl -L -o duckdb_cli.zip "https://github.com/duckdb/duckdb/releases/download/${DUCKDB_VERSION}/duckdb_cli-linux-${DUCKDB_ARCH}.zip" \
    && unzip duckdb_cli.zip \
    && rm duckdb_cli.zip

RUN mv duckdb /usr/local/bin/duckdb

RUN apt-get install libclang-dev -y


# Copy your Rust project directory
WORKDIR /app

# Copy your project here (replace with your actual copy command)
COPY . .

# Build your Rust project (replace with your actual build command)
RUN cargo build --release

ENV PATH="${PATH}:/app/target/release"
# Final image

ENV SHELL /bin/bash

CMD ["/bin/bash"]

