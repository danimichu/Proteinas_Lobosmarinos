

# Script completo para descargar secuencias, alinear y construir un árbol filogenético

# Crear carpetas necesarias
mkdir -p ../Data
mkdir -p ../Results

# Descargar las secuencias desde NCBI
./datasets download gene symbol MB --ortholog Otariidae --filename MB_otariidae.zip
./datasets download gene symbol UCP1 --ortholog Otariidae --filename UCP1_otariidae.zip
./datasets download gene symbol SOD1 --ortholog Otariidae --filename SOD1_otariidae.zip
./datasets download gene symbol CAT --ortholog Otariidae --filename CAT_otariidae.zip
./datasets download gene symbol HIF1A --ortholog Otariidae --filename HIF1A_otariidae.zip
./datasets download gene symbol COX4I1 --ortholog Otariidae --filename COX4I1_otariidae.zip
./datasets download gene symbol PPARGC1A --ortholog Otariidae --filename PPARGC1A_otariidae.zip
./datasets download gene symbol ATP5F1A --ortholog Otariidae --filename ATP5F1A_otariidae.zip
./datasets download gene symbol PPARGC1A --ortholog Otariidae --filename PPARGC1A_otariidae.zip
./datasets download gene symbol ATP5F1A --ortholog Otariidae --filename ATP5F1A_otariidae.zip

#Comandos para descomprimir y crear los archivos .fasta
# MB
unzip MB_otariidae.zip -d MB_otariidae_data
find MB_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > MB_otariidae.fasta

# UCP1
unzip UCP1_otariidae.zip -d UCP1_otariidae_data
find UCP1_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > UCP1_otariidae.fasta

# SOD1
unzip SOD1_otariidae.zip -d SOD1_otariidae_data
find SOD1_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > SOD1_otariidae.fasta

# CAT
unzip CAT_otariidae.zip -d CAT_otariidae_data
find CAT_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > CAT_otariidae.fasta

# HIF1A
unzip HIF1A_otariidae.zip -d HIF1A_otariidae_data
find HIF1A_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > HIF1A_otariidae.fasta

# COX4I1
unzip COX4I1_otariidae.zip -d COX4I1_otariidae_data
find COX4I1_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > COX4I1_otariidae.fasta

# PPARGC1A
unzip PPARGC1A_otariidae.zip -d PPARGC1A_otariidae_data
find PPARGC1A_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > PPARGC1A_otariidae.fasta

# ATP5F1A
unzip ATP5F1A_otariidae.zip -d ATP5F1A_otariidae_data
find ATP5F1A_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > ATP5F1A_otariidae.fasta

# GPX1
unzip GPX1_otariidae.zip -d GPX1_otariidae_data
find GPX1_otariidae_data/ncbi_dataset/data/ -name protein.faa -exec cat {} + > GPX1_otariidae.fasta

# Carga los módulos
module av blast
module load blast+/2.11.0

# Lista de tus archivos fasta
PROTEINS=(
  ATP5F1A_otariidae.fasta
  CAT_otariidae.fasta
  COX4I1_otariidae.fasta
  GPX1_otariidae.fasta
  HIF1A_otariidae.fasta
  MB_otariidae.fasta
  PPARGC1A_otariidae.fasta
  SOD1_otariidae.fasta
  UCP1_otariidae.fasta
)

# Itera sobre cada archivo y ejecuta makeblastdb + blastp
for PROT in "${PROTEINS[@]}"; do
  BASENAME="${PROT%.fasta}"
  echo "Procesando $BASENAME..."

  # Crear base de datos BLAST
  makeblastdb -in "$PROT" -dbtype prot -out "${BASENAME}_db" -parse_seqids -blastdb_version 4

  # Ejecutar BLASTp (contra sí mismo)
  blastp -query "$PROT" -db "${BASENAME}_db" -out "${BASENAME}_blast.out" -outfmt 6
done


# Alinear las secuencias con MUSCLE
echo "Realizando alineamiento con MUSCLE..."
muscle -in ../Data/MB.fasta -out ../Results/MB_aligned.fasta

# Construir árbol filogenético con IQ-TREE
echo "Construyendo árbol con IQ-TREE..."
iqtree -s ../Results/MB_aligned.fasta -m MFP -bb 1000 -nt AUTO -pre ../Results/MB_tree

echo "Análisis completo. Archivos generados en la carpeta Results/"
