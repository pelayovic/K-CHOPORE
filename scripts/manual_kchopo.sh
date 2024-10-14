#!/bin/bash

# Definir colores
BLUE='\033[0;34m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # Sin color

# Función para imprimir mensajes en color
print_in_color() {
    local color="$1"
    shift
    echo -e "${color}$*${NC}"
}

# Imprimir mensaje de bienvenida
print_in_color $CYAN "====================================================================================="
print_in_color $CYAN "   🌟 K-CHOPORE (Keen Comprehensive High-throughput Omics Pipeline Organizer) 🌟"
print_in_color $CYAN "====================================================================================="
print_in_color $GREEN "K-CHOPORE es una herramienta automatizada para el análisis de datos de secuenciación de Nanopore, que integra varios pasos, desde el basecalling hasta la detección de modificaciones epitranscriptómicas."
print_in_color $CYAN "====================================================================================="
print_in_color $CYAN "🌟 Just like a classic Asturian cachopo, K-CHOPORE is a hearty and satisfying tool for Nanopore sequencing data analysis! 🧬🍴 Dive in and enjoy the data feast with K-CHOPORE! 📊💻"

# Sección de Visión General
print_in_color $YELLOW "\n📜 Overview"
print_in_color $GREEN "K-CHOPORE is an automated pipeline for comprehensive analysis of Nanopore sequencing data."
print_in_color $GREEN "It handles the entire process from basecalling to epitranscriptomic modification detection."
print_in_color $GREEN "The pipeline integrates various specialized tools to ensure robust and reproducible results."

# Componentes del Pipeline
print_in_color $YELLOW "\n🔧 Pipeline Components"

print_in_color $GREEN "1. Basecalling with Guppy"
print_in_color $CYAN "Description: Converts raw Nanopore sequencing signals into readable DNA sequences."
print_in_color $CYAN "Command:"
print_in_color $BLUE "guppy_basecaller -i input_directory -s output_directory -c dna_r9.4.1_450bps_hac.cfg"

print_in_color $GREEN "\n2. Sequence Quality Check with pycoQC"
print_in_color $CYAN "Description: Assesses the quality of generated sequences."
print_in_color $CYAN "Installation and Execution:"
print_in_color $BLUE "pip3 install pycoQC"
print_in_color $BLUE "pycoQC -f ./pycoQC/anac017-1_C_R1_sequencing_summary.txt -a ./pycoQC/anac017-1_C_R1_sorted.bam -o anac017-1_C_R1.html --report_title \"anac017-1_C_R1\""

print_in_color $GREEN "\n3. Transcriptome Construction with FLAIR"
print_in_color $CYAN "Description: Constructs and aligns transcriptomes using FLAIR. This step is optional and user-configurable."
print_in_color $CYAN "Installation and Execution:"
print_in_color $BLUE "conda env create -n FLAIR_original -f FLAIR.yml"
print_in_color $BLUE "cd Desktop/flair_master"
print_in_color $BLUE "python flair.py align -g TAIR10_chr_all.fas -r WT_C_R1.fq.gz"

print_in_color $CYAN "Transcriptome alignment with samtools:"
print_in_color $BLUE "samtools sort flair.aligned.unsorted.bam -o flair.aligned.sorted.bam"
print_in_color $BLUE "samtools index flair.aligned.sorted.bam"
print_in_color $BLUE "python ./bin/bam2Bed12.py -i flair.aligned.sorted.bam > flair.aligned.sorted.bed"

print_in_color $CYAN "\nReference Genome Normalization:"
print_in_color $BLUE "java -jar picard.jar NormalizeFasta I=TAIR_chr_all.fas O=TAIR_chr_all_norm.fas"

print_in_color $GREEN "\n4. Additional Analysis"
print_in_color $CYAN "Description: In-depth analysis of the results, including multivariable analysis."
print_in_color $CYAN "Proposed Commands:"
print_in_color $BLUE "Result Analysis with DEXSeq and DRIMSeq"
print_in_color $CYAN "Multivariable Analysis (PCA, sPLS):"
print_in_color $BLUE "Perform PCA, sPLS to determine isoforms with specific characteristics."

# Sección de Configuración
print_in_color $YELLOW "\n⚙ Configuration"
print_in_color $GREEN "Key configuration files include:"
print_in_color $BLUE "config.yaml: The main configuration file where you specify paths, parameters, and options."
print_in_color $BLUE "environment.yml: The Conda environment configuration file with required dependencies."
print_in_color $BLUE "Snakefile: Defines the workflow rules and dependencies."

# Sección de Estructura del Pipeline
print_in_color $YELLOW "\n🏗️ Pipeline Structure"
print_in_color $GREEN "K-CHOPORE uses Snakemake for managing the workflow and dependencies."
print_in_color $CYAN "The pipeline includes:"
print_in_color $BLUE "- Basecalling"
print_in_color $BLUE "- Normalization"
print_in_color $BLUE "- Mapping"
print_in_color $BLUE "- Modification Detection"
print_in_color $BLUE "- Isoform Analysis (optional)"
print_in_color $BLUE "- Quality Control"
print_in_color $BLUE "- Multivariable Analysis"

# Sección de Casos de Uso
print_in_color $YELLOW "\n📊 Examples and Use Cases"
print_in_color $GREEN "Case Study 1: Analysis of m6A modifications in RNA samples."
print_in_color $GREEN "Case Study 2: Isoform quantification and differential expression analysis."
print_in_color $CYAN "Detailed examples are available in the results directory."

# Sección de Troubleshooting
print_in_color $YELLOW "\n🛠️ Troubleshooting"
print_in_color $GREEN "If you encounter issues, check the following:"
print_in_color $BLUE "1. Ensure all dependencies are correctly installed."
print_in_color $BLUE "2. Verify paths and file names in configuration files."
print_in_color $BLUE "3. Consult the logs in the logs directory for detailed error messages."

# Sección de Citaciones
print_in_color $YELLOW "\n🧬 Citations and References"
print_in_color $GREEN "Please cite the tools and references if you use K-CHOPORE in your research."
print_in_color $CYAN "Tools Integrated in K-CHOPORE:"
print_in_color $BLUE "1. Guppy"
print_in_color $BLUE "2. pycoQC"
print_in_color $BLUE "3. FLAIR"
print_in_color $BLUE "4. Picard"
print_in_color $BLUE "5. IsoformSwitchAnalyzeR"

# Sección de Contribuciones
print_in_color $YELLOW "\n🤝 Contributing"
print_in_color $GREEN "We welcome contributions from the community!"
print_in_color $BLUE "1. Fork the repository"
print_in_color $BLUE "2. Create a feature branch"
print_in_color $BLUE "3. Make your changes and test them"
print_in_color $BLUE "4. Submit a pull request"

# Sección de Contacto
print_in_color $YELLOW "\n📞 Contact and Support"
print_in_color $GREEN "For support, reach out via:"
print_in_color $BLUE "1. Email: pelayo.gonzalez@ispasturias.es"
print_in_color $BLUE "2. GitHub Issues: Submit an issue"

# Mensaje de cierre
print_in_color $CYAN "====================================================================================="
print_in_color $CYAN "                               ¡Gracias por usar K-CHOPORE!"
print_in_color $CYAN "====================================================================================="

