FROM ubuntu:22.04

ARG CBICALL_VERSION=development
ARG VCS_REF=unknown

LABEL org.opencontainers.image.authors="Manuel Rueda <manuel.rueda@cnag.eu>" \
      org.opencontainers.image.source="https://github.com/CNAG-Biomedical-Informatics/cbicall" \
      org.opencontainers.image.version="${CBICALL_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"

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
      bcftools \
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

# Install Python packages before copying the remaining source so dependency
# installation can be reused when only CBIcall code changes.
WORKDIR /usr/share/cbicall
COPY requirements.txt ./requirements.txt
RUN pip3 install --no-cache-dir -r requirements.txt

# The repository checkout is the Docker build context. In CI this is the exact
# commit selected by actions/checkout; local builds use the caller's checkout.
COPY . .

# Add user "dockeruser"
ARG UID=1000
ARG GID=1000

RUN groupadd -g "${GID}" dockeruser \
  && useradd --create-home --no-log-init -u "${UID}" -g "${GID}" dockeruser

# To change default user from root -> dockeruser
#USER dockeruser

# Set the entry dir
WORKDIR /usr/share/cbicall
