#!/usr/bin/env python3
import os
import subprocess
import sys

# Definir los archivos FASTQ, los nombres de muestras, y las rutas de salida
fastq_files = sys.argv[1].split(",")  # Separar los archivos FASTQ por comas
samples = sys.argv[2].split(",")      # Separar los nombres de las muestras por comas
reference_genome = sys.argv[3]
output_dir = sys.argv[4]

# Ruta completa a Minimap2
minimap2_path = "/usr/bin/minimap2"

# Verificar si el directorio de salida existe
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Iterar sobre los archivos FASTQ y sus respectivas muestras
for i, fastq in enumerate(fastq_files):
    sample = samples[i]
    sam_file = os.path.join(output_dir, f"{sample}.sam")
    
    print(f"[INFO] Running Minimap2 for sample {sample}...")
    
    # Ejecutar Minimap2 para cada archivo FASTQ por separado
    cmd = [minimap2_path, "-ax", "map-ont", reference_genome, fastq, "-o", sam_file]
    subprocess.run(cmd, check=True)
    
    print(f"[INFO] Mapping completed for sample {sample}, output saved to {sam_file}.")

