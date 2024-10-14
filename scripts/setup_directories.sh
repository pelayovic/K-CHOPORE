#!/bin/bash

# Función para imprimir mensajes con un borde decorativo
pretty_print() {
    local message=$1
    local border_length=${#message}
    local border=$(printf '%*s' "$border_length" '' | tr ' ' '=')
    echo
    echo " $border"
    echo " $message"
    echo " $border"
    echo
}

# Definir el directorio base (asumimos que es el directorio actual)
BASE_DIR="/home/pelamovic/K-CHOPORE"

# Crear la estructura de directorios si no existen
pretty_print "Creando la estructura de directorios para K-CHOPORE"

create_directory() {
    local dir=$1
    if [ ! -d "$dir" ]; then
        sudo mkdir -p "$dir"
        echo "Directorio creado: $dir"
    else
        echo "Directorio ya existe: $dir"
    fi
}

create_directory "$BASE_DIR/data/raw/fastq"
create_directory "$BASE_DIR/data/raw/bam"
create_directory "$BASE_DIR/data/raw/metadata"
create_directory "$BASE_DIR/data/processed/alignments"
create_directory "$BASE_DIR/data/processed/counts"
create_directory "$BASE_DIR/data/processed/quality_reports"
create_directory "$BASE_DIR/data/reference/genome"
create_directory "$BASE_DIR/data/reference/annotations"
create_directory "$BASE_DIR/scripts/preprocessing"
create_directory "$BASE_DIR/scripts/alignment"
create_directory "$BASE_DIR/scripts/postprocessing"
create_directory "$BASE_DIR/scripts/quality_control"
create_directory "$BASE_DIR/scripts/analysis"
create_directory "$BASE_DIR/scripts/utilities"
create_directory "$BASE_DIR/results/figures"
create_directory "$BASE_DIR/results/tables"
create_directory "$BASE_DIR/results/reports"
create_directory "$BASE_DIR/logs"
create_directory "$BASE_DIR/config"
create_directory "$BASE_DIR/notebooks"
create_directory "$BASE_DIR/docs"
create_directory "$BASE_DIR/envs"
create_directory "$BASE_DIR/publication/manuscript/drafts"
create_directory "$BASE_DIR/publication/manuscript/revisions"
create_directory "$BASE_DIR/publication/manuscript/submission"
create_directory "$BASE_DIR/publication/figures"
create_directory "$BASE_DIR/publication/tables"
create_directory "$BASE_DIR/publication/supplementary"
create_directory "$BASE_DIR/reviews/review1"
create_directory "$BASE_DIR/reviews/review2"
create_directory "$BASE_DIR/reviews/response"

pretty_print "Creación de directorios completada"

# Crear el archivo principal README.md si no existe
pretty_print "Creando el archivo README.md principal"

main_readme="$BASE_DIR/README.md"
if [ ! -f "$main_readme" ]; then
    sudo bash -c "cat > $main_readme" <<EOL
# K-CHOPORE
K-CHOPORE (Knowledgeable Comprehensive High-throughput Omics Pipeline Organizer) es una herramienta automatizada para el análisis de datos de secuenciación de Nanopore que integra varios pasos, desde el basecalling hasta la detección de modificaciones epitranscriptómicas.

## Componentes del Pipeline
- Basecalling con Guppy
- Calidad de secuencias con pycoQC
- Construcción de transcriptoma de referencia con FLAIR
- Mapeo y detección de modificaciones con Minimap2 y ELIGOS2

Este proyecto está organizado en varias carpetas para gestionar datos, scripts, resultados y otros recursos necesarios para el análisis.
EOL
    echo "Archivo README.md principal creado en: $main_readme"
else
    echo "Archivo README.md principal ya existe: $main_readme"
fi

# Crear README.md en los subdirectorios solo si no existen
pretty_print "Creando archivos README.md para los subdirectorios"

declare -A readme_texts=(
    ["data/raw/README.md"]="Este directorio contiene los datos en crudo, incluyendo archivos FASTQ, BAM y metadatos relacionados.\n\n## Estructura\n- fastq/: Archivos de secuencias en formato FASTQ.\n- bam/: Archivos de alineamientos en formato BAM.\n- metadata/: Archivos de metadatos asociados a los datos."
    ["data/processed/README.md"]="Este directorio contiene los datos procesados, incluyendo alineamientos, conteos y reportes de calidad.\n\n## Estructura\n- alignments/: Archivos de alineamientos procesados.\n- counts/: Archivos de conteos de expresión génica.\n- quality_reports/: Reportes de calidad generados durante el procesamiento."
    ["data/reference/README.md"]="Este directorio contiene los datos de referencia, incluyendo el genoma y las anotaciones.\n\n## Estructura\n- genome/: Archivos del genoma de referencia.\n- annotations/: Archivos de anotaciones genómicas."
    ["scripts/README.md"]="Este directorio contiene scripts para el procesamiento, alineación, post-procesamiento, control de calidad, análisis y utilidades.\n\n## Subdirectorios\n- preprocessing/: Scripts para la preparación y limpieza de datos.\n- alignment/: Scripts para el alineamiento de secuencias.\n- postprocessing/: Scripts para el post-procesamiento de datos.\n- quality_control/: Scripts para el control de calidad de los datos.\n- analysis/: Scripts para el análisis de datos.\n- utilities/: Scripts de utilidades generales y funciones auxiliares."
    ["results/README.md"]="Este directorio contiene los resultados del análisis, incluyendo figuras, tablas y reportes.\n\n## Estructura\n- figures/: Figuras generadas durante el análisis.\n- tables/: Tablas de resultados y datos resumidos.\n- reports/: Reportes finales y resúmenes de los análisis."
    ["logs/README.md"]="Este directorio contiene los archivos de registro generados durante el análisis.\n\n## Estructura\n- Archivos de log de ejecución del pipeline y de errores."
    ["config/README.md"]="Este directorio contiene archivos de configuración para el pipeline.\n\n## Archivos de Configuración\n- Archivos .conf o .json para configurar los parámetros del pipeline."
    ["notebooks/README.md"]="Este directorio contiene cuadernos (notebooks) para el análisis interactivo y visualización de datos.\n\n## Estructura\n- Cuadernos Jupyter o similares para análisis de datos interactivos."
    ["docs/README.md"]="Este directorio contiene la documentación adicional del proyecto.\n\n## Documentos\n- Documentación técnica, guías de usuario y manuales del proyecto."
    ["envs/README.md"]="Este directorio contiene archivos relacionados con entornos de ejecución y dependencias.\n\n## Archivos de Entorno\n- Archivos de configuración para entornos virtuales o contenedores."
    ["publication/README.md"]="Este directorio contiene materiales relacionados con la publicación, incluyendo borradores, revisiones y materiales suplementarios.\n\n## Subdirectorios\n- manuscript/: Documentos relacionados con el manuscrito.\n- figures/: Figuras para la publicación.\n- tables/: Tablas para la publicación.\n- supplementary/: Material suplementario para la publicación."
    ["reviews/README.md"]="Este directorio contiene archivos relacionados con las revisiones del manuscrito y las respuestas a los revisores.\n\n## Subdirectorios\n- review1/: Archivos relacionados con la primera revisión del manuscrito.\n- review2/: Archivos relacionados con la segunda revisión del manuscrito.\n- response/: Respuestas a los revisores."
)

for file in "${!readme_texts[@]}"; do
    target_file="$BASE_DIR/$file"
    if [ ! -f "$target_file" ]; then
        sudo bash -c "echo -e \"${readme_texts[$file]}\" > \"$target_file\""
        echo "Archivo README.md creado en: $target_file"
    else
        echo "Archivo README.md ya existe: $target_file"
    fi
done

pretty_print "Creación de archivos README.md completada con éxito"

