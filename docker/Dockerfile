FROM ubuntu:22.04

# File Author / Maintainer
LABEL org.opencontainers.image.authors="Manuel Rueda <manuel.rueda@cnag.eu>"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
# Ensure that system sbin directories are in PATH
ENV PATH="/usr/local/sbin:/usr/sbin:/sbin:${PATH}"

# Install Linux tools
RUN apt-get update && \
    apt-get -y install gcc make git cpanminus perl-doc vim sudo libperl-dev \
    libncurses5 libtinfo \
    python3-pip openjdk-17-jre

# Download Convert-Pheno
WORKDIR /usr/share/
RUN git clone https://github.com/mrueda/cbicall.git

# Remove the .git folder to save space and avoid shipping VCS data
RUN rm -rf cbicall/.git

# Install Perl modules
WORKDIR /usr/share/cbicall
RUN cpanm --notest --installdeps .

# Add user "dockeruser"
ARG UID=1000
ARG GID=1000

RUN groupadd -g "${GID}" dockeruser \
  && useradd --create-home --no-log-init -u "${UID}" -g "${GID}" dockeruser

# To change default user from root -> dockeruser
#USER dockeruser

# Get back to entry dir
WORKDIR /usr/share/cbicall
