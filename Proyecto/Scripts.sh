#!/bin/bash

# Script completo para descargar secuencias, alinear y construir un árbol filogenético

# Crear carpetas necesarias
mkdir -p ../Data
mkdir -p ../Results

# Lista de especies y gen objetivo (mioglobina)
genes="myoglobin"
especies=(
  "Zalophus wollebaeki"
  "Zalophus californianus"
  "Arctocephalus gazella"
  "Otaria flavescens"
  "Callorhinus ursinus"
)

echo "Descargando secuencias de mioglobina desde NCBI..."

# Eliminar archivo previo si existe
rm -f ../Data/MB.fasta

# Descargar las secuencias desde NCBI
for sp in "${especies[@]}"; do
  echo "  $sp"
  esearch -db protein -query "$genes[Gene] AND $sp[Organism]" \
    | efetch -format fasta >> ../Data/MB.fasta
done

# Mostrar cuántas líneas se descargaron
echo "Secuencias descargadas:"
wc -l ../Data/MB.fasta

# Alinear las secuencias con MUSCLE
echo "Realizando alineamiento con MUSCLE..."
muscle -in ../Data/MB.fasta -out ../Results/MB_aligned.fasta

# Construir árbol filogenético con IQ-TREE
echo "Construyendo árbol con IQ-TREE..."
iqtree -s ../Results/MB_aligned.fasta -m MFP -bb 1000 -nt AUTO -pre ../Results/MB_tree

echo "Análisis completo. Archivos generados en la carpeta Results/"
