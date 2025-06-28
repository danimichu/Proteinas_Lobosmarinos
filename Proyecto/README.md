# Proyecto Bioinformático: Análisis de Genes Ortólogos en Otariidae

## Autor: Daniela Centeno
## Propósito del programa de tu proyecto
Este proyecto tiene como finalidad estudiar genes ortólogos en especies representativas de la familia **Otariidae**, con énfasis en el *lobo marino de Galápagos* (*Zalophus wollebaeki*), con el fin de comprender la evolución molecular de proteínas relacionadas con adaptaciones marinas como el buceo, la termorregulación y el metabolismo energético.
## Justificación e Importancia

Los lobos marinos representan un excelente modelo de estudio para investigar la evolución de adaptaciones marinas en mamíferos. Comprender cómo han evolucionado ciertas proteínas funcionales puede revelar mecanismos moleculares detrás de la hipoxia durante el buceo, la conservación del calor en medios fríos, o la protección celular frente a especies reactivas de oxígeno. Esta información no solo es relevante para la biología evolutiva, sino también para estrategias de conservación de especies vulnerables, como *Zalophus wollebaeki*, endémica de Galápagos.
## Proteínas de Interés

Las proteínas que se analizarán están relacionadas con:

- **Hipoxia**: mioglobina (MB), citocromo c oxidasa (COX), HIF-1α (factor inducible por hipoxia).
- **Metabolismo energético**: PGC-1α (coactivador del metabolismo mitocondrial), ATP sintasa.
- **Termorregulación**: UCP1 (proteína desacoplante mitocondrial).
- **Estrés oxidativo**: superóxido dismutasa (SOD), glutatión peroxidasa (GPX), catalasa (CAT).
## Objetivos

1. **Identificar genes ortólogos** en diferentes géneros de la familia Otariidae.
2. **Comparar secuencias proteicas** relacionadas con adaptaciones fisiológicas clave.
3. **Analizar la evolución molecular** y la divergencia funcional de estas proteínas.
4. **Construir árboles filogenéticos** que reflejen relaciones evolutivas entre las especies.
5. **Visualizar y documentar** los resultados para generar una base bioinformática útil en conservación y estudios evolutivos.
## Organismo

- Grupo: Mamíferos marinos de la familia **Otariidae**.
- Géneros incluidos en el análisis:
  - Zalophus (ej. *Zalophus wollebaeki*, *Z. californianus*)
  - Arctocephalus
  - Eumetopias
  - Callorhinus
  - Otaria
- Foco especial: *Zalophus wollebaeki* (lobo marino de Galápagos), especie endémica de alto interés ecológico y evolutivo.
---
## Metodología

Se construirá un entorno bioinformático con Conda, y se recopilarán secuencias de proteínas desde bases de datos públicas como NCBI y UniProt. Se realizará la búsqueda de ortólogos mediante BLASTp, utilizando proteínas humanas como referencia inicial si fuera necesario. Las secuencias obtenidas serán alineadas con MUSCLE para evaluar su conservación. Posteriormente, se generarán árboles filogenéticos con IQ-TREE bajo modelos de máxima verosimilitud, y se visualizarán con FigTree para analizar relaciones evolutivas. Los resultados se interpretarán en función de la ecología de las especies, enfocándose en adaptaciones funcionales a la vida marina. Se usará GitHub para llevar el control de versiones y documentar el progreso del proyecto.

## Herramientas y Software

- **Conda** (gestión del entorno bioinformático)
- **NCBI / UniProt** (descarga de secuencias)
- **BLASTp** (identificación de ortólogos)
- **MUSCLE** (alineamiento múltiple)
- **IQ-TREE** (filogenia molecular)
- **FigTree** (visualización de árboles)
- **Atom / Nano** (edición de archivos)
- **Git / GitHub** (control de versiones)


## Fotos
![Lobo marino de Galápagos](https://datazone.darwinfoundation.org/images/checklist/dscn3722.jpg)  
*Lobo marino de Galápagos*

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




