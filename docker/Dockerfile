FROM ubuntu:22.04

# File Author / Maintainer
LABEL org.opencontainers.image.authors="Manuel Rueda <manuel.rueda@cnag.eu>"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
# Ensure that system sbin directories are in PATH
ENV PATH="/usr/local/sbin:/usr/sbin:/sbin:${PATH}"

# Install Linux tools (use --no-install-recommends and headless JRE)
RUN apt-get update && \
    apt-get install -y \
      git vim sudo \
      libncurses5 libtinfo5 \
      python3-pip openjdk-17-jre-headless ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Download CBICall
WORKDIR /usr/share/
RUN git clone https://github.com/mrueda/cbicall.git

# Remove the .git folder to save space and avoid shipping VCS data
RUN rm -rf cbicall/.git

# Install Python packages
RUN pip3 install -r requirements.txt --break-system-packages

# Add user "dockeruser"
ARG UID=1000
ARG GID=1000

RUN groupadd -g "${GID}" dockeruser \
  && useradd --create-home --no-log-init -u "${UID}" -g "${GID}" dockeruser

# To change default user from root -> dockeruser
#USER dockeruser

# Get back to entry dir
WORKDIR /usr/share/cbicall
