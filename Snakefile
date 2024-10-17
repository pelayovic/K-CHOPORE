# ---------------------------------------------------------
# Pipeline K-CHOPORE en Snakemake
# ---------------------------------------------------------
# Este pipeline de Snakemake automatiza el análisis de datos de secuenciación Nanopore,
# incluyendo basecalling, alineamiento, detección de modificaciones epitranscriptómicas,
# y análisis multivariable.
# Creado por Pelayo González de Lena Rodríguez, MSc
# FPI Severo Ochoa Fellow
# Laboratorio de Epigenética del Cáncer y Nanomedicina | FINBA
# Laboratorio de Biología de Sistemas | Universidad de Oviedo
# https://www.linkedin.com/in/biopelayo/
# https://gitlab.com/bio.pelayo/
# +34 660 74 11 39
# ---------------------------------------------------------

# ---------------------------------------------------------
# Pipeline K-CHOPORE en Snakemake
# ---------------------------------------------------------

import os

# Cargar el archivo de configuración
configfile: "config/config.yml"

# Definir variables desde el archivo de configuración
out = config["output"]["path"]

# Archivos FASTQ y resúmenes de secuenciación
FASTQ_FILES = config["input_files"]["fastq_files"]
SUMMARY_FILES = config["input_files"]["sequencing_summaries"]

# Directorios y archivos de referencia
FAST5_DIR = config["input_files"]["fast5_dir"]
reference_genome = config["input_files"]["reference_genome"]
reference_index = config["input_files"]["reference_genome_mmi"]
flair_output_dir = config["input_files"]["flair_output_dir"]
gtf_file = config["input_files"]["gtf_file"]
bed_file = config["input_files"]["bed_file"]

# Parámetros del pipeline
threads = config["params"]["threads"]

# Imprimir configuración para depuración
print(f"[INFO] Genoma de referencia: {reference_genome}")
print(f"[INFO] Archivos FASTQ: {FASTQ_FILES}")
print(f"[INFO] Resúmenes de secuenciación: {SUMMARY_FILES}")
print(f"[INFO] Directorio de salida: {out}")
print(f"[INFO] Directorio FAST5: {FAST5_DIR}")
print(f"[INFO] Índice del genoma de referencia: {reference_index}")
print(f"[INFO] Hilos: {threads}")

# ---------------------------------------------------------
# Reglas del Pipeline
# ---------------------------------------------------------

# Regla para configurar la estructura completa del proyecto
rule setup_complete_structure:
    output:
        "complete_structure_created.txt"
    run:
        # Directorios específicos del proyecto general
        base_dirs = [
            "config",
            "data",
            "docs",
            "envs",
            "logs",
            "notebooks",
            "publication",
            "results",
            "reviews",
            "scripts"
        ]

        # Directorios adicionales para procesamiento de datos dentro de "results"
        processing_dirs = [
            os.path.join("results", "basecalls"),
            os.path.join("results", "summaries"),
            os.path.join("results", "mapped"),
            os.path.join("results", "sorted_bam"),
            os.path.join("results", "quality_analysis"),
            os.path.join("results", "flair"),
            os.path.join("results", "eligos"),
            os.path.join("results", "multiqc")
        ]

        # Unir todos los directorios en una lista
        dirs = base_dirs + processing_dirs

        # Crear los directorios si no existen
        for d in dirs:
            if not os.path.exists(d):
                os.makedirs(d)
                print(f"[INFO] Directorio creado: {d}")
            else:
                print(f"[INFO] El directorio ya existe: {d}")

        # Marcar la creación de la estructura completa
        with open(output[0], 'w') as f:
            f.write("Estructura completa del proyecto creada.")
        print("[INFO] Estructura completa del proyecto creada exitosamente.")

# Regla para indexar el genoma de referencia usando Minimap2
rule index_genome:
    input:
        structure_created="complete_structure_created.txt",
        reference_genome=reference_genome
    output:
        reference_index=reference_index
    shell:
        """
        echo "[INFO] Indexando el genoma de referencia usando Minimap2..."
        minimap2 -d {output.reference_index} {input.reference_genome}
        echo "[INFO] Indexación del genoma completada."
        """

# Regla para mapear lecturas FASTQ usando Minimap2
rule map_with_minimap2:
    input:
        structure_created="complete_structure_created.txt",
        reference_index=reference_index,
        fastq_files=expand("data/raw/fastq/{sample}.fastq", sample=["WT_C_R1", "WT_C_R2"])
    output:
        sam_files=expand("results/mapped/{sample}.sam", sample=["WT_C_R1", "WT_C_R2"])
    params:
        threads=threads
    shell:
        """
        echo "[INFO] Mapeando archivos FASTQ con Minimap2..."
        for i in {input.fastq_files}; do
            sample=$(basename $i .fastq)
            minimap2 -ax map-ont {input.reference_index} $i -t {params.threads} > results/mapped/${{sample}}.sam
            echo "[INFO] Mapeo completado para la muestra $sample."
        done
        """

# Regla para ordenar e indexar archivos BAM
rule sort_and_index_bam:
    input:
        structure_created="complete_structure_created.txt",
        sam_file="results/mapped/{sample}.sam"
    output:
        bam_file="results/sorted_bam/{sample}_sorted.bam"
    shell:
        """
        echo "[INFO] Ordenando e indexando archivos BAM..."
        rm -f {output.bam_file}.tmp.*
        export PYTHONPATH=/usr/local/lib/python3.10/dist-packages/
        python3 /workspace/scripts/pipelines/sort_and_index_bam.py {input.sam_file} {output.bam_file}
        echo "[INFO] Ordenamiento e indexación completados."
        """

# Regla para análisis de calidad con pycoQC
rule quality_analysis_with_pycoQC:
    input:
        structure_created="complete_structure_created.txt",
        summary="data/raw/summaries/{sample}_sequencing_summary.txt",
        bam="results/sorted_bam/{sample}_sorted.bam"
    output:
        html="results/quality_analysis/pycoQC_output_{sample}.html"
    shell:
        """
        echo "[INFO] Ejecutando pycoQC para análisis de calidad..."
        pycoQC -f {input.summary} -b {input.bam} -o {output.html}
        echo "[INFO] Análisis de calidad completado para la muestra {wildcards.sample}."
        """

# Regla para ejecutar FLAIR para análisis de isoformas
# Regla para ejecutar FLAIR para análisis de isoformas

#rule run_flair:
#    input:
#        structure_created="complete_structure_created.txt",
#        fastq_files="data/raw/fastq/{sample}.fastq",
#        reference_genome=reference_genome,
#        gtf_file=gtf_file
#    output:
#        transcriptome_gtf="results/flair/flair.collapse.isoforms.gtf",
#        transcriptome_bed="results/flair/flair.collapse.isoforms.bed"
#    params:
#        flair_path="./scripts/pipelines/flair/src/flair"  # Ruta al directorio de FLAIR
#    shell:
#        """
#        echo "[INFO] Ejecutando FLAIR para análisis de isoformas..."
#        python {params.flair_path}/flair_align.py -g {input.reference_genome} -r {input.fastq_files} -o {output.transcriptome_gtf}
#        python {params.flair_path}/flair_collapse.py -g {input.reference_genome} -r {input.fastq_files} -f {input.gtf_file} -o {output.transcriptome_gtf}
#        gtfToGenePred {output.transcriptome_gtf} stdout | genePredToBed stdin {output.transcriptome_bed}
#        echo "[INFO] Análisis de isoformas con FLAIR completado."
#        """





# Regla para ejecutar FLAIR para análisis de isoformas
rule flair_quantify:
    input:
        isoforms_fa="results/flair/{sample}_flair.collapse.isoforms.fa",
        reads_manifest="data/raw/reads_manifest.tsv"
    output:
        tpm_output="results/flair/{sample}_quantify_tpm.tsv"
    params:
        threads=threads
    shell:
        """
        echo "[INFO] Cuantificando isoformas con FLAIR..."
        python3 /workspace/scripts/pipelines/flair/flair.py quantify -r {input.reads_manifest} -i {input.isoforms_fa} --salmon --tpm
        echo "[INFO] Cuantificación con FLAIR completada."

        """

# Regla para análisis diferencial de expresión con FLAIR
rule flair_diff_exp:
    input:
        counts_matrix="results/flair/counts_matrix.tsv"
    output:
        diff_exp_output="results/flair/diffExp_salmon/diffExp_results.tsv"
    params:
        out_dir="results/flair/diffExp_salmon",
        threads=threads
    shell:
        """
        echo "[INFO] Ejecutando análisis diferencial de expresión con FLAIR..."
        python3 /workspace/scripts/pipelines/flair/flair.py diffExp -q {input.counts_matrix} -o {params.out_dir}
        echo "[INFO] Análisis diferencial de expresión completado."
        """

# Regla para análisis de splicing diferencial con FLAIR
rule flair_diff_splice:
    input:
        isoforms_bed="results/flair/{sample}_flair.collapse.isoforms.bed",
        counts_matrix="results/flair/counts_matrix.tsv"
    output:
        diff_splice_output="results/flair/diffSplice_results.tsv"
    shell:
        """
        echo "[INFO] Ejecutando análisis de splicing diferencial con FLAIR..."
        python3 /workspace/scripts/pipelines/flair/flair.py diffSplice -i {input.isoforms_bed} -q {input.counts_matrix} --test
        echo "[INFO] Análisis de splicing diferencial completado."
        """




# regla eligos2 para epimarks
# Regla para ejecutar ELIGOS2
# Regla para ejecutar ELIGOS2 para análisis epitranscriptómico
rule run_eligos2:
    input:
        bam_file="results/sorted_bam/{sample}_sorted.bam",
        reference_genome="data/reference/genome/TAIR10_chr_all.fas.fasta",
        region_bed="/workspace/data/transcriptoma_FLAIR/flair.collapse.isoforms.bed"
    output:
        eligos_output="results/eligos/{sample}_eligos_output.txt"
    log:
        "logs/eligos2_{sample}.log"
    shell:
        """
        echo "[INFO] Ejecutando ELIGOS2 para análisis epitranscriptómico..."
        python3 /workspace/scripts/pipelines/eligos2/eligos2 rna_mod \
            -i {input.bam_file} \
            -reg {input.region_bed} \
            -ref {input.reference_genome} \
            -o results/eligos \
            --pval 0.05 \
            --oddR 5 \
            --esb 0.2 \
            > {log} 2>&1
        echo "[INFO] Análisis con ELIGOS2 completado para la muestra {wildcards.sample}."
        """



# Regla "all" para ejecutar todo el pipeline
# Regla "all" para ejecutar todo el pipeline
rule all:
    input:
        "complete_structure_created.txt",
        reference_index,
        expand("results/mapped/{sample}.sam", sample=["WT_C_R1", "WT_C_R2"]),
        expand("results/sorted_bam/{sample}_sorted.bam", sample=["WT_C_R1", "WT_C_R2"]),
        expand("results/quality_analysis/pycoQC_output_{sample}.html", sample=["WT_C_R1", "WT_C_R2"]),
        expand("results/flair/{sample}_flair.collapse.isoforms.gtf", sample=["WT_C_R1", "WT_C_R2"]),
        expand("results/flair/{sample}_flair.collapse.isoforms.bed", sample=["WT_C_R1", "WT_C_R2"]),
        expand("results/eligos/{sample}_eligos_output.txt", sample=["WT_C_R1", "WT_C_R2"])

