#!/bin/bash

# Detectar si el script está siendo ejecutado dentro de Docker
if grep -q docker /proc/1/cgroup; then
    FAST5_DIR="/workspace/data/raw/fast5"  # Ruta dentro del contenedor Docker
else
    FAST5_DIR="/home/pelamovic/K-CHOPORE/data/raw/fast5"  # Ruta local fuera de Docker
fi

SAMPLES_FILE="samples_list.txt"

# Verificar si el directorio de FAST5 existe
if [ ! -d "$FAST5_DIR" ]; then
    echo "Error: El directorio $FAST5_DIR no existe."
    exit 1
fi

# Listar archivos FAST5 y generar la lista de muestras
fast5_files=$(ls "$FAST5_DIR"/*.fast5 2> /dev/null)
if [ -z "$fast5_files" ]; then
    echo "Error: No se encontraron archivos FAST5 en $FAST5_DIR."
    exit 1
fi

# Crear o limpiar el archivo samples_list.txt
> "$SAMPLES_FILE"

# Extraer los nombres de las muestras y escribirlos en el archivo
for file in $fast5_files; do
    sample_name=$(basename "$file" .fast5)
    echo "$sample_name" >> "$SAMPLES_FILE"
done

echo "Lista de muestras generada en $SAMPLES_FILE."

