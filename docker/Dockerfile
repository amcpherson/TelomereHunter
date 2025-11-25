FROM python:2.7-slim

# --- Fix apt sources (Buster is EOL) ---
RUN sed -i 's|deb.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i 's|security.debian.org|archive.debian.org|g' /etc/apt/sources.list && \
    sed -i '/buster-updates/d' /etc/apt/sources.list && \
    sed -i '/deb-src/d' /etc/apt/sources.list

# --- Install system dependencies, R, samtools, and build deps for R packages ---
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gcc \
        g++ \
        gfortran \
        make \
        git \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libssl-dev \
        libffi-dev \
        libxml2-dev \
        libpng-dev \
        libfreetype6-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libcairo2-dev \
        libjpeg-dev \
        libtiff5-dev \
        libicu-dev \
        wget \
        ca-certificates \
        samtools \
        r-base \
        r-base-dev \
        && rm -rf /var/lib/apt/lists/*

# --- Upgrade pip (but keep Python 2 compatible) ---
RUN pip install --upgrade "pip<21.0"

# Install telomerehunter from PyPI
RUN pip install -e git+https://github.com/amcpherson/TelomereHunter.git#egg=telomerehunter

# --- Install required R packages from CRAN ---
RUN Rscript -e "install.packages('ggplot2', repos='https://packagemanager.posit.co/cran/2019-01-01'); library(ggplot2)"
RUN Rscript -e "install.packages('reshape2', repos='https://packagemanager.posit.co/cran/2019-01-01'); library(reshape2)"
RUN Rscript -e "install.packages('gridExtra', repos='https://packagemanager.posit.co/cran/2019-01-01'); library(gridExtra)"
RUN Rscript -e "install.packages('RColorBrewer', repos='https://packagemanager.posit.co/cran/2019-01-01'); library(RColorBrewer)"
RUN Rscript -e "install.packages('cowplot', repos='https://packagemanager.posit.co/cran/2019-01-01'); library(cowplot)"
RUN Rscript -e "install.packages('svglite', repos='https://packagemanager.posit.co/cran/2019-01-01'); library(svglite)"
RUN Rscript -e "install.packages('dplyr', repos='https://packagemanager.posit.co/cran/2019-01-01'); library(dplyr)"

# Optional: test that packages load
RUN Rscript -e "library(ggplot2); library(reshape2); library(gridExtra); library(RColorBrewer); library(cowplot); library(svglite); cat('R OK\n')"

ENTRYPOINT ["telomerehunter"]

