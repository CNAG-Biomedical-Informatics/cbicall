FROM ubuntu:22.04

# File Author / Maintainer
LABEL org.opencontainers.image.authors="Manuel Rueda <manuel.rueda@cnag.eu>"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
# Ensure that system sbin directories are in PATH
ENV PATH="/usr/local/sbin:/usr/sbin:/sbin:${PATH}"
ARG NEXTFLOW_VERSION=25.10.4
ARG CROMWELL_VERSION=92
ENV CROMWELL_DIR=/opt/cromwell
ENV CROMWELL_JAR=/opt/cromwell/cromwell-${CROMWELL_VERSION}.jar
ENV WOMTOOL_JAR=/opt/cromwell/womtool-${CROMWELL_VERSION}.jar

# Install Linux tools (use --no-install-recommends and headless JRE)
RUN apt-get update && \
    apt-get install -y \
      gcc make git vim sudo \
      libncurses5 libtinfo5 \
      python3-pip openjdk-17-jre-headless ca-certificates curl && \
    rm -rf /var/lib/apt/lists/*

# Install Nextflow. Ubuntu 22.04 does not provide a suitable default apt
# package, so the Docker image pins the official launcher version explicitly.
RUN curl -fsSL https://get.nextflow.io | NXF_VER="${NEXTFLOW_VERSION}" bash && \
    mv nextflow /usr/local/bin/nextflow && \
    chmod 0755 /usr/local/bin/nextflow && \
    nextflow -version

# Install Cromwell and WOMtool from the official Broad Institute release
# assets. CBIcall uses CROMWELL_JAR for execution and WOMTOOL_JAR for optional
# WDL validation in `validate-registry`.
RUN mkdir -p "${CROMWELL_DIR}" && \
    curl -fsSL \
      "https://github.com/broadinstitute/cromwell/releases/download/${CROMWELL_VERSION}/cromwell-${CROMWELL_VERSION}.jar" \
      -o "${CROMWELL_JAR}" && \
    curl -fsSL \
      "https://github.com/broadinstitute/cromwell/releases/download/${CROMWELL_VERSION}/womtool-${CROMWELL_VERSION}.jar" \
      -o "${WOMTOOL_JAR}" && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec java -jar "$CROMWELL_JAR" "$@"' > /usr/local/bin/cromwell && \
    printf '%s\n' '#!/usr/bin/env bash' 'exec java -jar "$WOMTOOL_JAR" "$@"' > /usr/local/bin/womtool && \
    chmod 0755 /usr/local/bin/cromwell /usr/local/bin/womtool && \
    cromwell --version && \
    test -s "${WOMTOOL_JAR}"

# Download CBICall
WORKDIR /usr/share/
RUN git clone https://github.com/mrueda/cbicall.git

# Remove Git metadata and CI configs to shrink image and avoid shipping them
RUN rm -rf cbicall/.git cbicall/.github cbicall/.gitattributes cbicall/.gitignore

# Remove documentation sources from the runtime image
RUN rm -rf cbicall/docs-site/

# Install Python packages, including Snakemake.
RUN pip3 install -r /usr/share/cbicall/requirements.txt

# Add user "dockeruser"
ARG UID=1000
ARG GID=1000

RUN groupadd -g "${GID}" dockeruser \
  && useradd --create-home --no-log-init -u "${UID}" -g "${GID}" dockeruser

# To change default user from root -> dockeruser
#USER dockeruser

# Set the entry dir
WORKDIR /usr/share/cbicall
