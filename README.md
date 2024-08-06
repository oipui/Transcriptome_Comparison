### Aims of the Study

This study leverages a single-cell RNA-sequencing (scRNA-seq) dataset of *Ciona intestinalis* to explore the origins of the vertebrate brain and identify cell types in the last common ancestor of urochordates and vertebrates. The key aims include:

1. **Annotation of the *Ciona* Nervous System**
   - Re-analyze the data to annotate cell types across different developmental stages.
   - Expand cell type annotations by considering blastomere lineages using a sequential and manual integration approach.
   - Validate the expanded annotations using pseudotime trajectory analysis.

2. **Elucidating the Brain of the Last Common Ancestor of Urochordates and Vertebrates**
   - Compare the nervous system of *Ciona* larvae to adult vertebrate brains (mouse and sea lamprey).
   - Utilize the Self Assembling Manifold Mapping algorithm ([SAMap](https://github.com/atarashansky/SAMap)) for comparative cell type analysis across three chordate species.
   - Identify homologous cell types between *Ciona*, mouse, and lamprey to reveal shared gene expression patterns present before the emergence of vertebrates.

This investigation aims to uncover the deep origins of the vertebrate brain through a detailed comparison of cell types in urochordates and vertebrates, providing insights into the evolutionary development of neural structures.

### Tools

- Comparative analysis: [SAMap](https://github.com/atarashansky/SAMap)
- Single cell analysis: [Seurat](https://satijalab.org/seurat/)
- Data visualisation: [ggplot2](https://ggplot2.tidyverse.org/)

