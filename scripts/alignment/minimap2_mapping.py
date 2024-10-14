import pysam
import subprocess

# Definir las rutas a los archivos y scripts
minimap2_script = "/home/pelamovic/K_CHOPO/k_chopo/bin/minimap2.py"
reference_genome = "/workspace/data/reference/genome/TAIR10_chr_all.fas"
fastq_file = "/workspace/data/raw/fastq/WT_C_R1.fq.gz"
output_bam = "/workspace/data/processed/alignments/mapped_reads.bam"

# Ruta temporal para el archivo SAM
sam_output = output_bam.replace(".bam", ".sam")

# 1. Ejecutar minimap2 para realizar el mapeo
cmd = [
    "python3", minimap2_script, "-ax", "map-ont", reference_genome, fastq_file, "-o", sam_output
]

# Ejecutar minimap2 y verificar que se ejecuta correctamente
subprocess.run(cmd, check=True)

# 2. Convertir el archivo SAM a BAM y ordenarlo usando pysam
pysam.sort("-o", output_bam, sam_output)

# 3. Indexar el archivo BAM resultante
pysam.index(output_bam)

# Imprimir el resultado final
print(f"Mapping completed. BAM file sorted and indexed at: {output_bam}")

