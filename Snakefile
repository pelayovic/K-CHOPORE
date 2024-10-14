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
counts_matrix = config["input_files"]["counts_matrix"]

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
        reference_genome = reference_genome
    output:
        reference_index = reference_index
    shell:
        """
        echo "[INFO] Indexando el genoma de referencia usando Minimap2..."
        minimap2 -d {output.reference_index} {input.reference_genome}
        echo "[INFO] Indexación del genoma completada."
        """


#Regla para indexar el genoma de referencia usando Minimap2
#rule basecalling:
#    input:
#        fast5_dir = FAST5_DIR
#    output:
#        fastq_dir = os.path.join(out, "basecalls")
#    params:
#        threads = threads,
#        basecaller = config["tools"]["basecaller"],
#        guppy_config = config["tools"]["guppy_config_file"],
#        bonito_model = config["tools"]["bonito_model"]
#    shell:
#        """
#        if [[ "{params.basecaller}" == "guppy" ]]; then
#            echo "[INFO] Ejecutando basecaller Guppy..."
#            guppy_basecaller --input_path {input.fast5_dir} --save_path {output.fastq_dir} \
#                             --config {params.guppy_config} --num_callers {params.threads}
#        elif [[ "{params.basecaller}" == "bonito" ]]; then
#            echo "[INFO] Ejecutando basecaller Bonito..."
#            bonito basecaller {params.bonito_model} {input.fast5_dir} > {output.fastq_dir}/bonito_output.fastq
#        fi
#        echo "[INFO] Basecalling completado."
#        """

# Regla para mapear lecturas FASTQ usando Minimap2


# Regla para mapear lecturas FASTQ usando Minimap2
rule map_with_minimap2:
    input:
        fastq_files=expand("data/raw/fastq/{sample}.fastq", sample=["WT_C_R1", "WT_C_R2"]),
        reference_index=reference_index
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
        sam_files = expand(os.path.join(out, "mapped", "{sample}.sam"), sample=["WT_C_R1", "WT_C_R2"])
    output:
        bam_files = expand(os.path.join(out, "sorted_bam", "{sample}_sorted.bam"), sample=["WT_C_R1", "WT_C_R2"])
    shell:
        """
        echo "[INFO] Ordenando e indexando archivos BAM..."
        export PYTHONPATH=/usr/local/lib/python3.10/dist-packages/
        python3 /workspace/scripts/pipelines/sort_and_index_bam.py {input.sam_files} {output.bam_files}
        echo "[INFO] Ordenamiento e indexación completados."
        """

# Regla para análisis de calidad con pycoQC
rule quality_analysis_with_pycoQC:
    input:
        summaries = SUMMARY_FILES,
        bam_files = expand(os.path.join(out, "sorted_bam", "{sample}_sorted.bam"), sample=["WT_C_R1", "WT_C_R2"])
    output:
        htmls = expand(os.path.join(out, "quality_analysis", "pycoQC_output_{sample}.html"), sample=["WT_C_R1", "WT_C_R2"])
    shell:
        """
        echo "[INFO] Ejecutando pycoQC para análisis de calidad..."
        for sample in WT_C_R1 WT_C_R2; do
            pycoQC -f {input.summaries} -b {input.bam_files} -o {output.htmls}
            echo "[INFO] Análisis de calidad completado para la muestra {sample}."
        done
        """

# Regla para ejecutar FLAIR para análisis de isoformas
rule run_flair:
    input:
        fastq_files = FASTQ_FILES,
        reference_genome = reference_genome,
        gtf_file = gtf_file
    output:
        transcriptome_gtf = os.path.join(out, "flair", "flair.collapse.isoforms.gtf"),
        transcriptome_bed = os.path.join(out, "flair", "flair.collapse.isoforms.bed")
    shell:
        """
        echo "[INFO] Ejecutando FLAIR para análisis de isoformas..."
        flair align -g {input.reference_genome} -r {input.fastq_files} -o {output.transcriptome_gtf}
        flair collapse -g {input.reference_genome} -r {input.fastq_files} -f {input.gtf_file} -o {output.transcriptome_gtf}
        gtfToGenePred {output.transcriptome_gtf} stdout | genePredToBed stdin {output.transcriptome_bed}
        echo "[INFO] Análisis de isoformas con FLAIR completado."
        """

# Regla para ejecutar ELIGOS2 para análisis epitranscriptómico
rule run_eligos2:
    input:
        bam_files = expand(os.path.join(out, "sorted_bam", "{sample}_sorted.bam"), sample=["WT_C_R1", "WT_C_R2"]),
        reference_genome = reference_genome,
        bed_file = bed_file
    output:
        eligos_output = expand(os.path.join(out, "eligos", "{sample}_eligos_output.txt"), sample=["WT_C_R1", "WT_C_R2"])
    shell:
        """
        echo "[INFO] Ejecutando ELIGOS2 para análisis epitranscriptómico..."
        for sample in WT_C_R1 WT_C_R2; do
            eligos2 rna_mod -i {input.bam_files} -reg {input.bed_file} -ref {input.reference_genome} \
            -o {output.eligos_output}
            echo "[INFO] Análisis con ELIGOS2 completado para la muestra {sample}."
        done
        """

# Regla para generar el informe de MultiQC
#rule run_multiqc:
#    input:
#        directories = [
#            os.path.join(out, dir) for dir in ["basecalls", "mapped", "sorted_bam", "quality_analysis", "flair", "eligos"]
#        ]
#    output:
#        multiqc_report = os.path.join(out, "multiqc", "multiqc_report.html")
#    shell:
#        """
#        echo "[INFO] Generando informe MultiQC..."
#        multiqc {input.directories} -o {output.multiqc_report}
#        echo "[INFO] Informe MultiQC generado en {output.multiqc_report}."
#        """

# Regla "all" para ejecutar todo el pipeline
rule all:
    input:
        "complete_structure_created.txt",
        expand(os.path.join(out, "mapped", "{sample}.sam"), sample=["WT_C_R1", "WT_C_R2"]),
        expand(os.path.join(out, "sorted_bam", "{sample}_sorted.bam"), sample=["WT_C_R1", "WT_C_R2"]),
        expand(os.path.join(out, "quality_analysis", "pycoQC_output_{sample}.html"), sample=["WT_C_R1", "WT_C_R2"]),
        os.path.join(out, "flair", "flair.collapse.isoforms.gtf"),
        os.path.join(out, "flair", "flair.collapse.isoforms.bed"),
        expand(os.path.join(out, "eligos", "{sample}_eligos_output.txt"), sample=["WT_C_R1", "WT_C_R2"]),
        #os.path.join(out, "multiqc", "multiqc_report.html")

