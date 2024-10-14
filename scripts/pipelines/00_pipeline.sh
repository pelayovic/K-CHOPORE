#!/bin/bash

# ---------------------------------------------------------
# K-CHOPORE (Keen Comprehensive High-throughput Omics Pipeline Organizer)
# An automated tool for analyzing Nanopore sequencing data.
#
# Integrates various essential steps from basecalling to epitranscriptomic modification detection.
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

# Parse the configuration file using yq and display all values
if [ -f "$CONFIG_FILE" ]; then
    echo "DEBUG: Parsing configuration file: $CONFIG_FILE"

    # Reading the input files section
    config_input_files_sequencing_summary=$(yq e '.input_files.sequencing_summary' "$CONFIG_FILE")
    config_input_files_reference_genome=$(yq e '.input_files.reference_genome' "$CONFIG_FILE")
    config_input_files_bam_file=$(yq e '.input_files.bam_file' "$CONFIG_FILE")
    config_input_files_fastq_file=$(yq e '.input_files.fastq_file' "$CONFIG_FILE")
    config_input_files_map_file=$(yq e '.input_files.map_file' "$CONFIG_FILE")

    # Reading the options section
    config_options_run_flair=$(yq e '.options.run_flair' "$CONFIG_FILE")
    config_options_r_version_required=$(yq e '.options.r_version_required' "$CONFIG_FILE")

    # Reading the output section
    config_output_path=$(yq e '.output.path' "$CONFIG_FILE")
    config_output_logs=$(yq e '.output.logs' "$CONFIG_FILE")
else
    echo "Error: Configuration file not found at $CONFIG_FILE"
    exit 1
fi

# Print all parsed configuration variables
echo "DEBUG: Configuration Variables Loaded from config.yml"
echo "Sequencing summary file path: $config_input_files_sequencing_summary"
echo "Reference genome path: $config_input_files_reference_genome"
echo "BAM file path: $config_input_files_bam_file"
echo "FASTQ file path: $config_input_files_fastq_file"
echo "Mapping output file: $config_input_files_map_file"
echo "Run FLAIR option: $config_options_run_flair"
echo "R version required: $config_options_r_version_required"
echo "Output path: $config_output_path"
echo "Log directory: $config_output_logs"
echo "---------------------------------------------------------"

# ------------------------------
# Pipeline Execution Steps
# ------------------------------

# ------------------------------
# Welcome Message
# ------------------------------
pretty_print "K-CHOPORE Pipeline - Build 2024"
echo "Keen Comprehensive High-throughput Omics Pipeline Organizer"
echo "An automated tool for analyzing Nanopore sequencing data."
echo "Developed by: Pelayo González de Lena Rodríguez, MSc"
echo "---------------------------------------------------------"

# ------------------------------
# Check R Version
# ------------------------------
pretty_print "Checking R Version"
R_VERSION=$(R --version | head -n 1 | awk '{print $3}')
echo "DEBUG: Installed R version: $R_VERSION"

REQUIRED_R_VERSION=$config_options_r_version_required

if [[ $(echo -e "$R_VERSION\n$REQUIRED_R_VERSION" | sort -V | head -n 1) == "$REQUIRED_R_VERSION" ]]; then
    echo "R version is compatible. Proceeding with package installation..."
    if ! Rscript -e 'install.packages("samplesizeCMH", repos="https://cloud.r-project.org")'; then
        echo "Error installing R package 'samplesizeCMH'. Please check your R installation and internet connection."
        exit 1
    fi
else
    echo "Incompatible R version ($R_VERSION). R version $REQUIRED_R_VERSION or higher is required."
    exit 1
fi

# ------------------------------
# Install and Configure Necessary Tools
# ------------------------------
pretty_print "Installing and Configuring Tools"

echo "Configuring PATH..."
export PATH=$PWD:$PWD/Scripts:$PATH || { echo "Error configuring PATH. Please check the directory path."; exit 1; }
echo "DEBUG: PATH configured as: $PATH"

# ------------------------------
# Verification and Preparation for PycoQC
# ------------------------------
pretty_print "Verification for PycoQC"
echo "DEBUG: Checking if sequencing summary file exists at: $config_input_files_sequencing_summary"

if [ ! -f "$config_input_files_sequencing_summary" ]; then
    echo "Sequencing summary file ($config_input_files_sequencing_summary) not found. Ensure the file is in the specified location."
    exit 1
fi

# ------------------------------
# Normalize Reference Genome
# ------------------------------
pretty_print "Normalizing Reference Genome"
echo "Normalizing the reference genome..."
echo "DEBUG: Reference genome input: $config_input_files_reference_genome"
echo "DEBUG: Normalized genome output: $config_input_files_reference_genome_norm"
if ! command -v java &> /dev/null; then
    echo "Error: java is not installed. Please install Java to continue."
    exit 1
fi
java -jar /usr/local/bin/picard.jar NormalizeFasta \
    I=/workspace/data/reference/genome/TAIR10_chr_all.fas \
    O=/workspace/data/reference/genome/TAIR10_chr_all_norm.fas

# ------------------------------
# Mapping and Modifications using minimap2 and ELIGOS
# ------------------------------
pretty_print "Mapping and Epitranscriptomic Modifications"

# 1. Ejecutar minimap2 (mapping)
echo "DEBUG: Ejecutando minimap2 para el mapeo con el archivo FASTQ: $config_input_files_fastq_file y referencia: $config_input_files_reference_genome_norm"
python3 /home/pelamovic/K-CHOPORE/scripts/alignment/minimap2_mapping.py \
    -ax map-ont $config_input_files_reference_genome_norm $config_input_files_fastq_file \
    -o $config_input_files_map_file

# 2. Verificar que minimap2 generó el archivo SAM
if [ ! -f "$config_input_files_map_file" ]; then
    echo "Error: El archivo SAM no se encontró en $config_input_files_map_file"
    exit 1
fi

# 3. Ordenar e indexar el archivo BAM usando pysam (después del mapeo)
echo "DEBUG: Ordenando e indexando el archivo BAM"
python3 /workspace/scripts/pipelines/sort_and_index_bam.py

# 4. Ejecutar epitranscriptomic modifications con ELIGOS
echo "Ejecutando modificaciones epitranscriptómicas con ELIGOS..."
./eligos2 rna_mod -i ${config_input_files_map_file%.sam}_sorted.bam -reg $config_input_files_bed_file -ref $config_input_files_reference_genome -o results --pval 0.05 --oddR 5 --esb 0.2
echo "DEBUG: ELIGOS2 results directory: $(pwd)/results"

echo "Performing multi-sample analysis with ELIGOS..."
./eligos2 multi_samples_test --test_mods ./results/*_baseExt0.txt --ctrl_mods *_baseExt0.txt --prefix output_prefix_
echo "DEBUG: Multi-sample test output prefix: output_prefix_"

# ------------------------------
# Check if FLAIR Should Be Executed
# ------------------------------
if [ "${config_options_run_flair}" == "true" ]; then
    pretty_print "Running FLAIR for Transcriptome Construction"
    
    cd $config_input_files_flair_repo
    echo "DEBUG: Changed to FLAIR repository directory: $(pwd)"
    
    echo "Activating FLAIR virtual environment..."
    source k-chopo/bin/activate

    echo "Aligning samples with FLAIR..."
    python flair.py align -g $config_input_files_reference_genome_norm -r $config_input_files_fastq_file
    echo "DEBUG: FLAIR aligned samples"

    echo "Generating BED file..."
    samtools sort flair.aligned.unsorted.bam -o flair.aligned.sorted.bam
    samtools index flair.aligned.sorted.bam
    python ./bin/bam2Bed12.py -i flair.aligned.sorted.bam > flair.aligned.sorted.bed
    echo "DEBUG: Generated BED file: flair.aligned.sorted.bed"
    
    echo "Collapsing alignments with FLAIR..."
    python flair.py collapse -g $config_input_files_reference_genome_norm -r $config_input_files_fastq_file -q flair.aligned.sorted.bed -f $config_input_files_gtf_file
    echo "DEBUG: Collapsed alignments"

    echo "Quantifying isoforms..."
    python flair.py quantify -g $config_input_files_reference_genome_norm -r $config_input_files_fastq_file -q flair.aligned.sorted.bed -f $config_input_files_gtf_file --output $config_input_files_flair_output_dir
    echo "DEBUG: Isoform quantification output directory: $config_input_files_flair_output_dir"
    
    echo "Summarizing counts and performing Fisher test..."
    python flair.py merge --inputs $config_input_files_counts_matrix --output $config_input_files_summarized_counts_matrix
    python flair.py fisher --inputs $config_input_files_summarized_counts_matrix --output $config_input_files_diff_isos_file
    echo "DEBUG: Fisher test output file: $config_input_files_diff_isos_file"
fi

# ------------------------------
# Run pycoQC
# ------------------------------
pretty_print "Running pycoQC for Quality Analysis"

if [ -f "$config_input_files_sequencing_summary" ]; then
    pycoQC -i $config_input_files_sequencing_summary -o /path/to/pycoQC_output
    echo "DEBUG: pycoQC output directory: /path/to/pycoQC_output"
else
    echo "Sequencing summary file not found. Check the file and its path."
    exit 1
fi

# ------------------------------
# Additional Isoform and Multivariable Analysis
# ------------------------------
pretty_print "Additional Isoform and Multivariable Analysis"

echo "Executing additional isoform analysis..."
python /path/to/isoformswitch_analyzer/isoformswitch_analyzer.py -i $config_input_files_iso_usage_file -o $config_input_files_isoformswitch_analyzer_output
echo "DEBUG: Isoform switch analysis output directory: $config_input_files_isoformswitch_analyzer_output"

echo "Performing PCA..."
python /path/to/pca_analysis.py -i $config_input_files_counts_matrix -o $config_input_files_pca_output
echo "DEBUG: PCA analysis output directory: $config_input_files_pca_output"

echo "Performing SPLS analysis..."
python /path/to/spls_analysis.py -i $config_input_files_counts_matrix -o $config_input_files_spls_output
echo "DEBUG: SPLS analysis output directory: $config_input_files_spls_output"

# ------------------------------
# Important Note on FASTA and GTF Name Matching
# ------------------------------
pretty_print "Important Note"

echo "Ensure that the names in the reference genome FASTA file match those in the GTF file to avoid issues during analysis."

# ------------------------------
# Completion Message
# ------------------------------
pretty_print "K-CHOPORE Pipeline - Completed"
echo "All steps have been successfully executed."



