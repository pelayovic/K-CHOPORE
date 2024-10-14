#!/bin/bash

# ---------------------------------------------------------
# K-CHOPORE (Keen Comprehensive High-throughput Omics Pipeline Organizer)
# Un sistema automatizado para el análisis de datos de secuenciación Nanopore.
#
# Integra múltiples pasos esenciales desde el basecalling hasta la detección de modificaciones epitranscriptómicas.
#
# Desarrollado por:
# Pelayo González de Lena Rodríguez, MSc
# FPI Severo Ochoa Fellow
# Cancer Epigenetics and Nanomedicine Lab | FINBA
# Systems Biology Lab | Universidad de Oviedo
#
# Contacto:
# LinkedIn: https://www.linkedin.com/in/biopelayo/
# GitLab: https://gitlab.com/bio.pelayo/
# ---------------------------------------------------------

# ---------------------------------------------------------
# Cargar el archivo de configuración desde config.yml
# ---------------------------------------------------------

#!/bin/bash

# ---------------------------------------------------------
# K-CHOPORE (Keen Comprehensive High-throughput Omics Pipeline Organizer)
# An automated tool for analyzing Nanopore sequencing data.
#
# This script handles basecalling, alignment, isoform analysis,
# and epitranscriptomic modifications detection.
#
# Developed by:
# Pelayo González de Lena Rodríguez, MSc
# FPI Severo Ochoa Fellow
# Cancer Epigenetics and Nanomedicine Lab | FINBA
# Systems Biology Lab | University of Oviedo
#
# Contact:
# LinkedIn: https://www.linkedin.com/in/biopelayo/
# GitLab: https://gitlab.com/bio.pelayo/
#
# ---------------------------------------------------------
# Config file path (adjust if located elsewhere)
CONFIG_FILE="/workspace/config/config.yml"

# ---------------------------------------------------------
# Function to pretty print headers
# ---------------------------------------------------------
pretty_print() {
    local message=$1
    local border_length=${#message}
    local border=$(printf '%*s' "$border_length" '' | tr ' ' '-')
    echo
    echo " $border"
    echo " $message"
    echo " $border"
    echo
}

# ---------------------------------------------------------
# Parse YAML using yq
# ---------------------------------------------------------
if ! command -v yq &> /dev/null; then
    echo "Error: 'yq' no está instalado. Instálalo antes de continuar."
    exit 1
fi

# Parse the configuration file using yq and verify all variables
if [ -f "$CONFIG_FILE" ]; then
    echo "DEBUG: Parsing configuration file: $CONFIG_FILE"

    # Leer sección de archivos de entrada
    config_input_files_fast5_dir=$(yq e '.input_files.fast5_dir' "$CONFIG_FILE")
    config_input_files_sequencing_summary=$(yq e '.input_files.sequencing_summary' "$CONFIG_FILE")
    config_input_files_reference_genome=$(yq e '.input_files.reference_genome' "$CONFIG_FILE")
    config_input_files_fastq_file=$(yq e '.input_files.fastq_file' "$CONFIG_FILE")
    config_input_files_bam_file=$(yq e '.input_files.bam_file' "$CONFIG_FILE")
    config_input_files_map_file=$(yq e '.input_files.map_file' "$CONFIG_FILE")
    config_input_files_flair_repo=$(yq e '.input_files.flair_repo' "$CONFIG_FILE")
    config_input_files_gtf_file=$(yq e '.input_files.gtf_file' "$CONFIG_FILE")
    config_input_files_bed_file=$(yq e '.input_files.bed_file' "$CONFIG_FILE")
    config_input_files_iso_usage_file=$(yq e '.input_files.iso_usage_file' "$CONFIG_FILE")
    config_input_files_counts_matrix=$(yq e '.input_files.counts_matrix' "$CONFIG_FILE")
    config_input_files_pca_output=$(yq e '.input_files.pca_output' "$CONFIG_FILE")
    config_input_files_spls_output=$(yq e '.input_files.spls_output' "$CONFIG_FILE")
    config_input_files_isoformswitch_analyzer_output=$(yq e '.input_files.isoformswitch_analyzer_output' "$CONFIG_FILE")

    # Leer sección de opciones
    config_options_run_flair=$(yq e '.options.run_flair' "$CONFIG_FILE")
    config_options_r_version_required=$(yq e '.options.r_version_required' "$CONFIG_FILE")

    # Leer sección de salida
    config_output_path=$(yq e '.output.path' "$CONFIG_FILE")
    config_output_logs=$(yq e '.output.logs' "$CONFIG_FILE")

    # Print configuration variables for debugging purposes
    echo "DEBUG: Configuration variables:"
    echo "  - FAST5 directory: $config_input_files_fast5_dir"
    echo "  - Sequencing summary: $config_input_files_sequencing_summary"
    echo "  - Reference genome: $config_input_files_reference_genome"
    echo "  - FASTQ output file: $config_input_files_fastq_file"
    echo "  - BAM file: $config_input_files_bam_file"
    echo "  - Mapping file: $config_input_files_map_file"
    echo "  - FLAIR repository: $config_input_files_flair_repo"
    echo "  - GTF file: $config_input_files_gtf_file"
    echo "  - BED file: $config_input_files_bed_file"
    echo "  - Isoform usage file: $config_input_files_iso_usage_file"
    echo "  - Counts matrix: $config_input_files_counts_matrix"
    echo "  - PCA output: $config_input_files_pca_output"
    echo "  - SPLS output: $config_input_files_spls_output"
    echo "  - IsoformSwitchAnalyzer output: $config_input_files_isoformswitch_analyzer_output"
else
    echo "Error: Archivo de configuración no encontrado en $CONFIG_FILE"
    exit 1
fi

# ------------------------------
# Instalar herramientas necesarias (MultiQC)
# ------------------------------
pretty_print "Instalando herramientas necesarias"

if ! command -v multiqc &> /dev/null; then
    echo "Instalando MultiQC..."
    pip install multiqc
fi

# ------------------------------
# Preguntar al usuario sobre qué basecaller utilizar
# ------------------------------
# ------------------------------


# ------------------------------
# Preguntar al usuario sobre qué basecaller utilizar
# ------------------------------
pretty_print "Selecciona el basecaller para convertir FAST5 a FASTQ"


# ------------------------------
# Preguntar al usuario sobre qué basecaller utilizar
# ------------------------------
pretty_print "Selecciona el basecaller para convertir FAST5 a FASTQ"

echo "Elige el software de basecalling:"
echo "1) Guppy"
echo "2) Albacore"
echo "3) Dorado"

# Bucle para obtener una selección válida
while true; do
    read -r -p "Ingresa el número de tu elección (1, 2 o 3): " basecaller_choice
    echo "DEBUG: Has ingresado '$basecaller_choice'"  # Depuración para ver lo que se ingresa
    if [[ "$basecaller_choice" =~ ^[123]$ ]]; then
        case "$basecaller_choice" in
            1)
                echo "Has seleccionado Guppy."
                basecaller="guppy_basecaller"
                ;;
            2)
                echo "Has seleccionado Albacore."
                basecaller="read_fast5_basecaller.py"
                ;;
            3)
                echo "Has seleccionado Dorado."
                basecaller="dorado_basecaller"
                ;;
        esac
        break
    else
        echo "Selección inválida. Por favor elige entre 1, 2 o 3."
    fi
done
# ------------------------------
# Convertir FAST5 a FASTQ + Sequencing Summary
# ------------------------------
pretty_print "Convirtiendo FAST5 a FASTQ + Sequencing Summary"
echo "Procesando archivos FAST5..."

# Lógica para ejecutar el basecaller correspondiente
case $basecaller_choice in
    1) # Guppy
        guppy_basecaller -i $config_input_files_fast5_dir -s $config_input_files_fastq_file --compress_fastq --flowcell FLO-MIN106 --kit SQK-LSK109
        ;;
    2) # Albacore
        read_fast5_basecaller.py --input $config_input_files_fast5_dir --save_path $config_input_files_fastq_file --recursive --output_format fastq
        ;;
    3) # Dorado
        dorado_basecaller $config_input_files_fast5_dir $config_input_files_reference_genome -o $config_input_files_fastq_file
        ;;
esac

# Verificar que se generó el archivo FASTQ
if [ ! -f "$config_input_files_fastq_file" ]; then
    echo "Error: El archivo FASTQ no se generó en $config_input_files_fastq_file"
    exit 1
fi

# Validar existencia de sequencing summary
if [ ! -f "$config_input_files_sequencing_summary" ]; then
    echo "Error: El archivo Sequencing Summary no se encontró en $config_input_files_sequencing_summary"
    exit 1
fi

# ------------------------------
# Mapeo con Minimap2
# ------------------------------
pretty_print "Mapeando FASTQ con Minimap2"
echo "DEBUG: Ejecutando Minimap2..."

minimap2 -ax map-ont "$config_input_files_reference_genome" "$config_input_files_fastq_file" > "$config_input_files_map_file"

# Verificar existencia del archivo SAM
if [ ! -f "$config_input_files_map_file" ]; then
    echo "Error: El archivo SAM no se generó correctamente"
    exit 1
fi

# ------------------------------
# Ordenar e Indexar con Pysam
# ------------------------------
pretty_print "Ordenando e Indexando BAM con Pysam"
python3 /workspace/scripts/pipelines/sort_and_index_bam.py "$config_input_files_map_file" "$config_input_files_bam_file"

# Verificar BAM generado
if [ ! -f "$config_input_files_bam_file" ]; then
    echo "Error: El archivo BAM no se generó correctamente"
    exit 1
fi

# ------------------------------
# Análisis de calidad con PycoQC
# ------------------------------
pretty_print "Corriendo PycoQC para análisis de calidad"

pycoQC -i "$config_input_files_sequencing_summary" -o /path/to/pycoQC_output

# ------------------------------
# Ejecución de MultiQC
# ------------------------------
pretty_print "Ejecutando MultiQC para resumen de calidad"

multiqc "$config_output_path" -o "$config_output_path/multiqc_report"

# ------------------------------
# FLAIR para construcción del transcriptoma
# ------------------------------
if [ "${config_options_run_flair}" == "true" ]; then
    pretty_print "Ejecutando FLAIR"
    
    cd "$config_input_files_flair_repo"
    echo "DEBUG: Cambiando al directorio del repo de FLAIR: $(pwd)"
    
    python flair.py align -g "$config_input_files_reference_genome" -r "$config_input_files_fastq_file"
    python flair.py collapse -g "$config_input_files_reference_genome" -r "$config_input_files_fastq_file" -f "$config_input_files_gtf_file"
    python flair.py quantify -g "$config_input_files_reference_genome" -r "$config_input_files_fastq_file"
fi

# ------------------------------
# ELIGOS2 para detección de modificaciones epitranscriptómicas
# ------------------------------
pretty_print "Ejecutando ELIGOS2 para detección de modificaciones epitranscriptómicas"

eligos2 rna_mod -i "$config_input_files_bam_file" -reg "$config_input_files_bed_file" -ref "$config_input_files_reference_genome" -o results --pval 0.05 --oddR 5 --esb 0.2

eligos2 multi_samples_test --test_mods ./results/*_baseExt0.txt --ctrl_mods *_baseExt0.txt --prefix output_prefix_

# ------------------------------
# Additional Isoform and Multivariable Analysis
# ------------------------------
pretty_print "Análisis Adicional de Isoformas y Multivariable"

echo "Ejecutando análisis adicional de isoformas..."
python /path/to/isoformswitch_analyzer/isoformswitch_analyzer.py -i "$config_input_files_iso_usage_file" -o "$config_input_files_isoformswitch_analyzer_output"
echo "DEBUG: Isoform switch analysis output directory: $config_input_files_isoformswitch_analyzer_output"

echo "Realizando análisis PCA..."
python /path/to/pca_analysis.py -i "$config_input_files_counts_matrix" -o "$config_input_files_pca_output"
echo "DEBUG: PCA analysis output directory: $config_input_files_pca_output"

echo "Realizando análisis SPLS..."
python /path/to/spls_analysis.py -i "$config_input_files_counts_matrix" -o "$config_input_files_spls_output"
echo "DEBUG: SPLS analysis output directory: $config_input_files_spls_output"

# ------------------------------
# Important Note on FASTA and GTF Name Matching
# ------------------------------
pretty_print "Nota Importante"

echo "Asegúrate de que los nombres en el archivo FASTA del genoma de referencia coincidan con los del archivo GTF para evitar problemas durante el análisis."

# ------------------------------
# Completion Message
# ------------------------------
pretty_print "Pipeline K-CHOPORE Completado"
echo "Todos los pasos se han ejecutado con éxito."

