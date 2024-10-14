# ----------------------------------------------
# 🌟 K-CHOPORE 🌟
# A hearty and satisfying tool for Nanopore sequencing data analysis!
# This Dockerfile builds K-CHOPORE, a powerful pipeline for Nanopore data processing
# using ELIGOS2, Bonito, and other essential tools for bioinformatics analysis.
#
# Base image: ELIGOS2 (piroonj/eligos2)
# Maintained by: pelayovic
# ----------------------------------------------

# Usar la imagen de ELIGOS2 desde DockerHub como base para empezar con un entorno especializado
# FROM piroonj/eligos2

# Cambiar a Ubuntu 22.04 como base para mayor compatibilidad con librerías y herramientas modernas
FROM ubuntu:22.04

# Evitar interacción durante las instalaciones, asegurando que no se soliciten confirmaciones al usuario
ENV DEBIAN_FRONTEND=noninteractive

# ----------------------------------------------
# Información del mantenedor del contenedor
# ----------------------------------------------
LABEL maintainer="pelayovic"

# ----------------------------------------------
# Instalación de herramientas esenciales
# ----------------------------------------------
# - yq: Utilizado para manipular archivos YAML, útil para el manejo de configuraciones en K-CHOPORE
# - curl: Necesario para descargar recursos desde la web
# - jq: Utilizado para manipular archivos JSON, necesario para algunos pasos del pipeline

RUN apt-get update && apt-get install -y curl jq && \
    curl -L https://github.com/mikefarah/yq/releases/download/v4.9.8/yq_linux_amd64 -o /usr/bin/yq && \
    chmod +x /usr/bin/yq

# ----------------------------------------------
# Instalar Java (default-jdk) y Picard
# ----------------------------------------------
# Java es necesario para ejecutar herramientas como Picard, que se utiliza para el manejo de archivos BAM
# Picard es una herramienta popular para manipulación y validación de archivos de secuenciación.

RUN apt-get update && apt-get install -y default-jdk wget && \
    wget https://github.com/broadinstitute/picard/releases/download/2.25.7/picard.jar -P /usr/local/bin/

# ----------------------------------------------
# Instalar Python 3, pip, y crear entornos virtuales
# ----------------------------------------------
# Python es el lenguaje base para el pipeline de K-CHOPORE.
# Aquí también se instala `pip`, el manejador de paquetes de Python, y `venv` para crear entornos virtuales.

RUN apt-get update && \
    apt-get install -y python3 python3-pip python3-venv && \
    rm -rf /var/lib/apt/lists/*

# ----------------------------------------------
# Instalación de Snakemake y Pulp
# ----------------------------------------------
# Snakemake es el motor de ejecución del pipeline de K-CHOPORE, mientras que Pulp es una biblioteca
# para resolver problemas de programación lineal que usa Snakemake para gestionar tareas.
# En esta sección instalamos las versiones compatibles y adecuadas de ambos.

RUN pip install --upgrade pip && \
    pulp==2.7.0
RUN pip install --upgrade snakemake
RUN sed -i 's/list_solvers/listSolvers/' /usr/local/lib/python3.10/dist-packages/snakemake/__init__.py

# ----------------------------------------------
# Instalación de Bonito para basecalling
# ----------------------------------------------
# Bonito es un basecaller especializado para Nanopore, que ofrece capacidades avanzadas como
# el procesamiento de bases modificadas y modelos basados en Transformers. Aquí se instala Bonito y
# sus dependencias principales, incluyendo `flash-attn` para modelos de transformación.

RUN pip install ont-bonito

# ----------------------------------------------
# Instalación de herramientas de desarrollo y compilación necesarias
# ----------------------------------------------
# Estas herramientas son fundamentales para compilar paquetes de Python y otros sistemas que se usarán
# en el pipeline de K-CHOPORE.

RUN apt-get update && apt-get install -y \
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
    rm -rf /var/lib/apt/lists/*

#
# ----------------------------------------------
# Instala minimap2
# ----------------------------------------------
# 

# Instalar minimap2
RUN apt-get update && apt-get install -y minimap2


# Download and install Dorado
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.0-linux-x64.tar.gz -P /opt/ \
    && tar -xzvf /opt/dorado-0.8.0-linux-x64.tar.gz -C /opt/ \
    && rm /opt/dorado-0.8.0-linux-x64.tar.gz \
    && ln -s /opt/dorado-0.8.0-linux-x64/bin/dorado /usr/local/bin/dorado

# Make sure the dorado executable is available in the $PATH
ENV PATH="/opt/dorado-0.8.0-linux-x64/bin:$PATH"

# ----------------------------------------------
# Instalación de dependencias de R
# ----------------------------------------------
# Algunos análisis estadísticos en K-CHOPORE se ejecutan en R. Aquí instalamos el paquete "samplesizeCMH"
# que es utilizado para realizar cálculos de tamaño de muestra en los análisis.

RUN Rscript -e 'install.packages("samplesizeCMH", repos="https://cloud.r-project.org")'
# ----------------------------------------------
# Instalación de flash-attn
# ----------------------------------------------
# Este paquete es utilizado por Bonito para acelerar los cálculos relacionados con modelos de trans># Clonamos el repositorio de `flash-attn` y lo instalamos manualmente para garantizar su integració>
# ----------------------------------------------
RUN pip install tabulate==0.9.0
# ----------------------------------------------
# Configuración de CUDA para evitar conflictos
# ----------------------------------------------
# Se establece CUDA_HOME para evitar conflictos de versiones de CUDA al usar modelos basados en GPU.ENV CUDA_HOME=/usr/local/cuda

# ----------------------------------------------
# Descargar modelos de Bonito
# ----------------------------------------------
# Bonito descargará automáticamente los modelos la primera vez que se utilicen,
# pero aquí se descargan todos los modelos disponibles para asegurar que estén listos para el pipel>
RUN bonito download --models

# Agregar el directorio de Bonito al PATH del sistema para poder ejecutarlo desde cualquier lugar.
ENV PATH="/usr/local/bin/bonito:$PATH"

# Verificar que Bonito está correctamente instalado y que los modelos están disponibles.
RUN bonito download --models --show
# ----------------------------------------------
# Copiar el archivo de requisitos de Python (requirements.txt)
# ----------------------------------------------
# Las dependencias de Python necesarias para K-CHOPORE se especifican en este archivo,
# que será copiado al contenedor y las dependencias serán instaladas usando `pip`.

# COPY requirements.txt /home/eligos2/
#RUN pip3 install --no-cache-dir -r /home/eligos2/requirements.txt
RUN pip install numpy==1.23.5

# ----------------------------------------------
# Configuración del entorno de trabajo
# ----------------------------------------------
# Se establece el directorio de trabajo `/workspace` donde se ejecutarán todos los comandos del pipeline.

WORKDIR /workspace

# ----------------------------------------------
# Añadir los directorios de eligos2 y Scripts al PATH del sistema
# ----------------------------------------------
# Esto permitirá ejecutar scripts de eligos2 desde cualquier lugar dentro del contenedor.

ENV PATH="/home/eligos2:/home/eligos2/Scripts:$PATH"

# ----------------------------------------------
# Confirmación de la correcta construcción de la imagen
# ----------------------------------------------
# Este comando imprime un mensaje para confirmar que la imagen de Docker se construyó con éxito.
# ----------------------------------------------
# Instalar Guppy
# ----------------------------------------------
# Elimina la versión anterior de Guppy si está instalada
RUN rm -rf /opt/ont-guppy

# Instala la versión de Guppy solo para CPU
RUN wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_6.1.5_linux64.tar.gz -O /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz && \
    tar -xvzf /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz -C /opt/ && \
    echo 'PATH=/opt/ont-guppy-cpu/bin:$PATH' >> ~/.bashrc && \
    rm /opt/ont-guppy-cpu_6.1.5_linux64.tar.gz

# Asegurarse de que Guppy está en el PATH
# Instalar CUDA
# RUN apt-get update && apt-get install -y nvidia-cuda-toolkit



ENV PATH="/opt/ont-guppy-cpu/bin:$PATH"
# Instala flair-brookslab directamente desde PyPI

RUN pip install flair-brookslab
# Instalar rpy2 para la integración de R con Python
RUN pip install rpy2
# Copiar el requirements.txt al contenedor
COPY requirements.txt .

# Instalar las dependencias

RUN pip install --no-cache-dir -r requirements.txt
# Copiar la carpeta scripts al contenedor
COPY scripts /workspace/scripts
COPY data /workspace/data
COPY config /workspace/config
COPY requirements.txt /workspace/requirements.txt
ENV PYTHONPATH="/usr/local/lib/python3.10/dist-packages"


RUN echo "🌟 K-CHOPORE Docker image built successfully! Dive into Nanopore sequencing data analysis! 🍴"

# ----------------------------------------------
# (Opcional) Comando por defecto al iniciar el contenedor
# ----------------------------------------------
# Si se desea un comando que se ejecute automáticamente al iniciar el contenedor, se puede descomentar
# el siguiente `ENTRYPOINT`. En este caso, podría ejecutar el pipeline de K-CHOPORE automáticamente.

# ENTRYPOINT ["python3", "/home/eligos2/k_chopo.py"]


