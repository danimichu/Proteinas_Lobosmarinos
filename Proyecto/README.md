# Proyecto Bioinformático: Análisis de Genes Ortólogos en Otariidae  
**Autora:** Daniela Centeno

## Propósito del proyecto  
Este proyecto tiene como finalidad estudiar genes ortólogos en especies representativas de la familia **Otariidae**, con el fin de comprender la evolución molecular de proteínas relacionadas con adaptaciones marinas como el buceo, la termorregulación y el metabolismo energético. A través del análisis comparativo entre géneros de otáridos, se busca identificar patrones evolutivos funcionales relevantes para la vida marina.

## Justificación e Importancia  
Los lobos marinos y osos marinos constituyen un excelente modelo para investigar la evolución de adaptaciones fisiológicas a entornos marinos extremos en mamíferos. Comprender cómo han evolucionado ciertas proteínas funcionales puede revelar mecanismos moleculares relacionados con la tolerancia a la hipoxia durante el buceo, la eficiencia energética en condiciones frías y la protección celular frente a especies reactivas de oxígeno (ROS). Esta información es clave no solo para la biología evolutiva, sino también para futuras investigaciones en conservación y fisiología comparada de mamíferos marinos.

## Proteínas de Interés  

### Hipoxia:  
- Mioglobina (MB)  
- Citocromo c oxidasa subunidad 4 (COX4I1)  
- Factor inducible por hipoxia (HIF-1α)  

### Metabolismo energético:  
- PPARGC1A (coactivador del metabolismo mitocondrial)  
- ATP5F1A (subunidad alfa de la ATP sintasa mitocondrial)  

### Termorregulación:  
- UCP1 (proteína desacoplante 1)  

### Estrés oxidativo:  
- Superóxido dismutasa (SOD1)  
- Glutatión peroxidasa (GPX1)  
- Catalasa (CAT)

## Objetivos

### Objetivo general  
Estudiar la evolución molecular de genes ortólogos relacionados con adaptaciones fisiológicas marinas en especies de la familia Otariidae.

### Objetivos específicos  
1. Identificar y comparar secuencias ortólogas de proteínas funcionales en Otariidae.  
2. Analizar la conservación y divergencia en secuencias proteicas asociadas a hipoxia, metabolismo energético y estrés oxidativo.  
3. Construir árboles filogenéticos para explorar las relaciones evolutivas entre las especies.  
4. Interpretar los cambios moleculares en el contexto de adaptaciones funcionales a la vida marina.  
5. Generar una base bioinformática para estudios evolutivos y de conservación.

## Organismo  

**Grupo:** Mamíferos marinos de la familia *Otariidae*

**Géneros incluidos en el análisis:**  
- Zalophus (ej. *Zalophus californianus*)  
- Arctocephalus  
- Eumetopias  
- Callorhinus  

> Nota: Aunque no se incluye *Zalophus wollebaeki* por falta de datos disponibles, el análisis comparativo entre otros miembros de Otariidae permite inferencias funcionales y evolutivas que podrían ser relevantes para especies endémicas de interés ecológico.

## Metodología  
Se construirá un entorno bioinformático con **Conda**. Las secuencias proteicas serán recopiladas de bases de datos públicas como **NCBI** y **UniProt**. Se utilizará **BLASTp** para identificar ortólogos utilizando proteínas humanas o modelos como referencia cuando sea necesario.  
Las secuencias obtenidas serán alineadas mediante **MUSCLE** para evaluar su conservación. Posteriormente, se construirán árboles filogenéticos con **IQ-TREE**, seleccionando el mejor modelo de sustitución, y se visualizarán con **FigTree**. Se interpretarán los resultados en relación con la ecología y fisiología de los otáridos. El proyecto será documentado en **GitHub** para llevar control de versiones y asegurar la reproducibilidad.

## Herramientas y Software  
- Conda (gestión de entornos bioinformáticos)  
- NCBI / UniProt (fuentes de secuencias)  
- BLASTp (identificación de ortólogos)  
- MUSCLE (alineamiento múltiple de secuencias)  
- IQ-TREE (análisis filogenético por máxima verosimilitud)  
- FigTree (visualización de filogenias)  
- Atom / Nano (edición de archivos)  
- Git / GitHub (control de versiones y documentación del proyecto)



## Fotos

![Zalophus wollebaeki](https://datazone.darwinfoundation.org/images/checklist/cp_1194.jpg)  
*Zalophus wollebaeki*

![Zalophus californianus](https://www.racerocks.ca/wp-content/uploads/2015/06/rm2010calsl.jpg)  
*Zalophus californianus*

![Otaria flavescens](https://ecoregistros.org/site/images/dataimages/2019/09/01/348207/Leao-marinho--Otaria-flavescens--macho-.jpg)  
*Otaria flavescens*

![Eumetopias jubatus](https://sealion-world.com/wp-content/uploads/Steller-Sea-Lion_624.jpg)  
*Eumetopias jubatus*

![Arctocephalus townsendi](https://zooinstitutes.com/img/animals/24/24234.jpg)  
*Arctocephalus townsendi*


## Lista de comandos

* Descargar las secuencias desde NCBI

```
./datasets download gene symbol MB --ortholog Otariidae --filename MB_otariidae.zip
./datasets download gene symbol UCP1 --ortholog Otariidae --filename UCP1_otariidae.zip
./datasets download gene symbol SOD1 --ortholog Otariidae --filename SOD1_otariidae.zip
./datasets download gene symbol CAT --ortholog Otariidae --filename CAT_otariidae.zip
./datasets download gene symbol HIF1A --ortholog Otariidae --filename HIF1A_otariidae.zip
./datasets download gene symbol COX4I1 --ortholog Otariidae --filename COX4I1_otariidae.zip
./datasets download gene symbol PPARGC1A --ortholog Otariidae --filename PPARGC1A_otariidae.zip
./datasets download gene symbol ATP5F1A --ortholog Otariidae --filename ATP5F1A_otariidae.zip

```
* Carga los módulos

```
module av blast
module load blast+/2.11.0
```
* Lista de tus archivos fasta

```
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
```
* Itera sobre cada archivo y ejecuta makeblastdb + blastp
```
for PROT in "${PROTEINS[@]}"; do
  BASENAME="${PROT%.fasta}"
  echo "Procesando $BASENAME..."
```
* Crear base de datos BLAST
```
  makeblastdb -in "$PROT" -dbtype prot -out "${BASENAME}_db" -parse_seqids -blastdb_version 4

* Ejecutar BLASTp (contra sí mismo)
  blastp -query "$PROT" -db "${BASENAME}_db" -out "${BASENAME}_blast.out" -outfmt 6
done

```
* Alinear las secuencias con MUSCLE
```
cat *.fasta > combinado.fasta
./muscle3.8.31_i86linux64 -in combinado.fasta -out aligned_.fasta
```
* Limpiar y dar formato a todos los encabezados en terminal
```
awk '
BEGIN { OFS="\t" }
{
    if ($0 ~ /^>/) {
        match($0, /[A-Z0-9\-]+/, gen);
        match($0, /\[organism=[^]]+\]/);
        org = substr($0, RSTART+10, RLENGTH-11);
        gsub(" ", "_", org);
        clave = gen[0] "_" org;
        contador[clave]++;
        iso = "X" contador[clave];
        print ">" gen[0] "_" iso "_" org;
    } else {
        print $0
    }
}' aligned_.fasta > alineado.fasta

```
