# ----------------------------------------------
# 🌟 K-CHOPORE 🌟
# A hearty and satisfying tool for Nanopore sequencing data analysis!
# This Dockerfile builds K-CHOPORE, a powerful pipeline for Nanopore data processing
# using ELIGOS2, Bonito, and other essential tools for bioinformatics analysis.
#
# Base image: Ubuntu 22.04
# Maintained by: pelayovic
# ----------------------------------------------

FROM ubuntu:22.04

# Evitar interacción durante las instalaciones
ENV DEBIAN_FRONTEND=noninteractive

# ----------------------------------------------
# Información del mantenedor
# ----------------------------------------------
LABEL maintainer="pelayovic"

# ----------------------------------------------
# Instalación de herramientas esenciales
# ----------------------------------------------
# - yq, curl, jq
# ----------------------------------------------
RUN echo "🔄 Installing essential tools: curl, jq, yq..." && \
    apt-get update && apt-get install -y curl jq && \
    curl -L https://github.com/mikefarah/yq/releases/download/v4.9.8/yq_linux_amd64 -o /usr/bin/yq && \
    chmod +x /usr/bin/yq && \
    echo "✅ Essential tools installed successfully!"

# ----------------------------------------------
# Instalar Java (default-jdk) y Picard
# ----------------------------------------------
RUN echo "🔄 Installing Java and Picard..." && \
    apt-get update && apt-get install -y default-jdk wget && \
    wget https://github.com/broadinstitute/picard/releases/download/2.25.7/picard.jar -P /usr/local/bin/ && \
    echo "✅ Java and Picard installed successfully!"

# ----------------------------------------------
# Instalar Python 3, pip y crear entornos virtuales
# ----------------------------------------------
RUN echo "🔄 Installing Python 3, pip, and setting up virtual environments..." && \
    apt-get update && apt-get install -y python3 python3-pip python3-venv && \
    rm -rf /var/lib/apt/lists/* && \
    echo "✅ Python 3, pip, and virtual environments installed!"

# ----------------------------------------------
# Instalación de Snakemake y Pulp
# ----------------------------------------------
RUN echo "🔄 Installing Snakemake and Pulp..." && \
    pip install --upgrade pip && \
    pip install pulp==2.7.0 && \
    pip install --upgrade snakemake && \
    sed -i 's/list_solvers/listSolvers/' /usr/local/lib/python3.10/dist-packages/snakemake/__init__.py && \
    echo "✅ Snakemake and Pulp installed!"

# ----------------------------------------------
# Instalación de Bonito para basecalling
# ----------------------------------------------
RUN echo "🔄 Installing Bonito for basecalling..." && \
    pip install ont-bonito && \
    echo "✅ Bonito installed!"

# ----------------------------------------------
# Instalación de herramientas de desarrollo y compilación
# ----------------------------------------------
RUN echo "🔄 Installing development tools and libraries..." && \
    apt-get update && apt-get install -y \
    build-essential \
    libbz2-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev \
    libssl-dev \
    libffi-dev \
    bedtools \
    git \
    r-base \
    default-jre && \
    rm -rf /var/lib/apt/lists/* && \
    echo "✅ Development tools and libraries installed!"

# ----------------------------------------------
# Instalar minimap2
# ----------------------------------------------
RUN echo "🔄 Installing Minimap2..." && \
    apt-get update && apt-get install -y minimap2 && \
    echo "✅ Minimap2 installed!"

# ----------------------------------------------
# Instalar Dorado para basecalling
# ----------------------------------------------
RUN echo "🔄 Installing Dorado for basecalling..." && \
    wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.0-linux-x64.tar.gz -P /opt/ \
    && tar -xzvf /opt/dorado-0.8.0-linux-x64.tar.gz -C /opt/ \
    && rm /opt/dorado-0.8.0-linux-x64.tar.gz \
    && ln -s /opt/dorado-0.8.0-linux-x64/bin/dorado /usr/local/bin/dorado && \
    echo "✅ Dorado installed!"

# ----------------------------------------------
# Instalar dependencias de R
# ----------------------------------------------
RUN echo "🔄 Installing R dependencies..." && \
    Rscript -e 'install.packages("samplesizeCMH", repos="https://cloud.r-project.org")' && \
    echo "✅ R dependencies installed!"

# ----------------------------------------------
# Descargar modelos de Bonito
# ----------------------------------------------
RUN echo "🔄 Downloading Bonito models..." && \
    bonito download --models && \
    echo "✅ Bonito models downloaded!"

# ----------------------------------------------
# Copiar archivo de requisitos de Python y instalar dependencias
# ----------------------------------------------
RUN pip install tabulate==0.9.0
# ----------------------------------------------
# Instalación de Guppy (versión CPU)
# ----------------------------------------------
RUN echo "🔄 Installing Guppy (CPU version)..." && \
    wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.1.5_linux64.tar.gz -O /opt/ont-guppy-cpu.tar.gz && \
    tar -xvzf /opt/ont-guppy-cpu.tar.gz -C /opt/ && \
    echo 'PATH=/opt/ont-guppy-cpu/bin:$PATH' >> ~/.bashrc && \
    rm /opt/ont-guppy-cpu.tar.gz && \
    echo "✅ Guppy (CPU version) installed!"

# ----------------------------------------------
# Instalar flair-brookslab
# ----------------------------------------------
RUN echo "🔄 Installing FLAIR for isoform analysis..." && \
    pip install flair-brookslab && \
    echo "✅ FLAIR installed!"

# ----------------------------------------------
# Instalar rpy2 para la integración de R con Python
# ----------------------------------------------
RUN echo "🔄 Installing rpy2 for R and Python integration..." && \
    pip install rpy2 && \
    echo "✅ rpy2 installed!"
# Las dependencias de Python necesarias para K-CHOPORE se especifican en este archivo,
# que será copiado al contenedor y las dependencias serán instaladas usando `pip`.

COPY requirements.txt /home/eligos2/
RUN pip3 install --no-cache-dir -r /home/eligos2/requirements.txt
RUN pip install numpy==1.23.5

# ----------------------------------------------
# Configuración del entorno de trabajo
# ----------------------------------------------
WORKDIR /workspace

# ----------------------------------------------
# Añadir directorios al PATH del sistema
# ----------------------------------------------
ENV PATH="/home/eligos2:/home/eligos2/Scripts:/opt/ont-guppy-cpu/bin:$PATH"

# ----------------------------------------------
# Mensaje de confirmación al terminar la instalación
# ----------------------------------------------
RUN echo "🌟 K-CHOPORE Docker image built successfully! Dive into Nanopore sequencing data analysis! 🌟"

# ----------------------------------------------
# (Opcional) Comando por defecto al iniciar el contenedor
# ----------------------------------------------
# ENTRYPOINT ["python3", "/home/eligos2/k_chopo.py"]
