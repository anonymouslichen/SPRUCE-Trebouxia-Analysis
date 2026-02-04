# Data Availability

## Overview

This directory would contain the data files used in the analysis. Due to repository size constraints and data archiving requirements, the full datasets are hosted externally.

## Available Data

### Dryad Repository
The following processed data files are available at [https://doi.org/10.5061/dryad.rbnzs7hcf](https://doi.org/10.5061/dryad.rbnzs7hcf):

1. **Trebouxia_Metadata.csv** - Sample metadata
   - `sample_id`: Unique identifier for each DNA extraction
   - `individual`: Lichen thallus ID
   - `thallus_type`: "core" or "tips" (sampling location within thallus)
   - `temp`: Temperature treatment (°C above ambient)
   - `CO2`: CO2 treatment level (450 or 900 ppm)

2. **Trebouxia_OTU_table.csv** - OTU abundance matrix
   - Rows: Samples
   - Columns: OTU identifiers (based on Muggia et al. 2020 reference sequences)
   - Values: Read counts per OTU per sample
   - Initial format before aggregation of duplicate OTUs

3. **Trebouxia_tax_assignments.csv** - Taxonomic classifications
   - `OTU_id`: Simplified OTU identifier
   - `taxonomy`: Full taxonomic string (Kingdom through Species)
   - Individual rank columns: Kingdom, Phylum, Class, Order, Family, Genus, Species

### GenBank
Representative sequences for novel *Trebouxia* OTUs identified in this study:
- **Accession numbers**: OP879470 - OP879475
- These represent the 6 unknown OTUs (UNK_01, UNK_03-UNK_07) that were phylogenetically placed and renamed

### Raw Sequencing Data
Raw FASTQ files from Illumina MiSeq sequencing are **available upon request** to the corresponding author:
- **Contact**: Abigail R. Meyer (meye2099@umn.edu)

Files include:
- Fall 2020 sequencing run
- Spring 2022 sequencing run
- 2×300bp paired-end reads
- ITS1-ITS2 amplicons

*Note: Raw data were not deposited to NCBI SRA but are archived locally and available for scientific use.*

## Additional Files Required

The following files are needed to run the complete pipeline but are not included in this repository:

### Reference Database
- **A_C_I_S_unaligned_ITS.qza**: *Trebouxia* reference database from Muggia et al. (2020)
  - Contains representative ITS sequences for all known *Trebouxia* OTUs
  - Used for open-reference clustering in QIIME2
  - Can be reconstructed from sequences in Muggia et al. (2020) publication

### Mapping Files
- **evernia_mappingfile_F20.txt**: QIIME2 manifest for Fall 2020 samples
- **spruce_evernia_mappingfile_S22.txt**: QIIME2 manifest for Spring 2022 samples
- Format: Tab-separated with columns for sample-id and absolute-filepath

### Phylogenetic Tree
- **RAxML_bestTree.legend.tree**: Phylogenetic tree with placed unknowns
  - Generated using RAxML evolutionary placement algorithm
  - Based on Muggia et al. (2020) multilocus phylogeny
  - Used for creating phylogram in Figure 8A

## Using This Data

To reproduce the analysis:

1. Download processed data from Dryad
2. Place files in this `data/` directory
3. Update file paths in the R script (`02_Trebouxia_analysis.R`)
4. For full QIIME2 pipeline, contact authors for raw FASTQ files

## Data Usage and Citation

If you use these data, please cite:

> Meyer, A.R., M. Valentin, L. Liulevicius, T.R. McDonald, M.P. Nelsen, J. Pengra, R.J. Smith, and D. Stanton. 2023. Climate warming causes photobiont degradation and carbon starvation in a boreal climate sentinel lichen. *American Journal of Botany* 110(2): e16114. https://doi.org/10.1002/ajb2.16114

## Questions?

For questions about data availability or access to raw sequencing files, please contact:

**Abigail R. Meyer**  
Email: meye2099@umn.edu  
ORCID: 0000-0001-6375-1897
