#!/usr/bin/env python3
# ---------------------------------------------------
# K-CHOPORE: RNA modifications and epitranscriptomics analysis
# Script to create multiple reduced FAST5 files from input with random sampling
# ---------------------------------------------------

import h5py
import random

# Ruta del archivo FAST5 original
input_fast5 = '/home/pelamovic/K-CHOPORE/FAR91957_pass_a56dafa5_2.fast5'

# Función para copiar lecturas seleccionadas a un nuevo archivo FAST5
def crear_fast5_reducido(input_file, output_file, selected_reads):
    try:
        with h5py.File(input_file, 'r') as archivo_entrada:
            with h5py.File(output_file, 'w') as archivo_salida:
                print(f"[INFO] Creando archivo FAST5 reducido: {output_file}")

                # Copiar las lecturas seleccionadas al nuevo archivo
                for read_key in selected_reads:
                    # Copiar la estructura del grupo
                    archivo_entrada.copy(read_key, archivo_salida)
                    print(f"[INFO] Copiando lectura: {read_key}")

                print(f"[INFO] Archivo FAST5 reducido creado con éxito con {len(selected_reads)} lecturas.")

    except Exception as e:
        print(f"[ERROR] Ocurrió un problema al crear el archivo FAST5 reducido: {e}")

# Obtener el tamaño total de lecturas del archivo original
with h5py.File(input_fast5, 'r') as archivo_entrada:
    total_reads = list(archivo_entrada.keys())
    total_count = len(total_reads)
    # Calcular el número de lecturas para el 33%
    n_lecturas = int(total_count * 0.33)

# Generar 10 archivos FAST5 reducidos
for i in range(1, 11):
    # Seleccionar lecturas aleatorias sin duplicados
    selected_reads = random.sample(total_reads, n_lecturas)
    output_fast5 = f'/home/pelamovic/K-CHOPORE/FAR91957_pass_a56dafa5_reducido_{i}.fast5'
    crear_fast5_reducido(input_fast5, output_fast5, selected_reads)

