import pysam
import sys
import os






import os
import pysam  # Asegura que pysam se cargue correctamente

# Añadir la ruta de pysam si es necesario
if '/usr/local/lib/python3.10/dist-packages/' not in os.sys.path:
    os.sys.path.append('/usr/local/lib/python3.10/dist-packages/')

# Aquí va el resto de tu código
# por ejemplo, la lógica para sortear e indexar archivos BAM



# Comprobar si se han proporcionado los archivos SAM y el genoma de referencia como argumentos
if len(sys.argv) < 3:
    raise ValueError("Por favor, proporciona los archivos SAM y el genoma de referencia como argumentos.")

sam_files = sys.argv[1:-1]  # Lista de archivos SAM
reference_genome = sys.argv[-1]  # Genoma de referencia

# Directorio de salida donde se escribirán los archivos BAM
output_dir = "/workspace/results/sorted_bam"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Procesar cada archivo SAM
for sam_file in sam_files:
    # Crear el nombre del archivo BAM añadiendo _sorted
    sample_name = os.path.basename(sam_file).replace(".sam", "")
    bam_file = os.path.join(output_dir, f"{sample_name}_sorted.bam")

    # Ordenar el archivo SAM y convertirlo a BAM
    print(f"Ordenando el archivo SAM: {sam_file}")
    pysam.sort("-o", bam_file, sam_file)
    print(f"Archivo BAM ordenado creado: {bam_file}")

    # Indexar el archivo BAM
    print(f"Indexando el archivo BAM: {bam_file}")
    pysam.index(bam_file)
    print(f"Archivo de índice creado para: {bam_file}")

