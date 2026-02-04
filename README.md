# SPRUCE Trebouxia Photobiont Analysis

Analysis pipeline for characterizing *Trebouxia* photobiont communities in the lichen *Evernia mesomorpha* under experimental climate warming conditions.

## Publication

**Meyer, A.R., M. Valentin, L. Liulevicius, T.R. McDonald, M.P. Nelsen, J. Pengra, R.J. Smith, and D. Stanton. 2023.** Climate warming causes photobiont degradation and carbon starvation in a boreal climate sentinel lichen. *American Journal of Botany* 110(2): e16114. [https://doi.org/10.1002/ajb2.16114](https://doi.org/10.1002/ajb2.16114)

## Project Overview

This repository contains the bioinformatics pipeline for analyzing algal photobiont community composition in *Evernia mesomorpha* lichens from the SPRUCE (Spruce and Peatland Responses Under Changing Environments) whole-ecosystem climate manipulation experiment. The study revealed that experimental warming (>2°C) caused photobiont loss and shifts in *Trebouxia* community composition, particularly an increase in clade 'I' photobionts with warming.

### Key Findings
- Photobiont communities shifted toward greater representation of *Trebouxia* clade 'I' with increasing temperature
- Core regions of thalli showed higher photobiont alpha diversity than tips
- Multiple *Trebouxia* genotypes coexist within individual lichen thalli

## Skills Demonstrated

This project showcases proficiency in:
- **High-throughput sequencing analysis**: Illumina MiSeq paired-end ITS metabarcoding
- **QIIME2**: Complete pipeline from raw reads to OTU tables
- **Bioinformatics**: Quality filtering, taxonomic assignment, open-reference clustering
- **R programming**: Data manipulation (dplyr, tidyr), statistical analysis (vegan, metacoder)
- **Phylogenetic analysis**: RAxML placement, tree manipulation (ggtree, ape)
- **Data visualization**: ggplot2, publication-quality figures
- **Unix/Bash**: Command-line bioinformatics workflows on HPC systems

## Repository Contents

```
SPRUCE_Trebouxia_Analysis/
├── README.md                           # This file
├── scripts/
│   ├── 01_QIIME2_pipeline.sh          # QIIME2 processing pipeline
│   └── 02_Trebouxia_analysis.R        # R analysis and figure generation
├── data/
│   └── README.md                       # Data availability information
└── docs/
    └── workflow_overview.md            # Detailed methods description
```

## Workflow Summary

### 1. Sequence Processing (QIIME2)
- Import Illumina MiSeq FASTQ files
- Remove ITS primers using cutadapt
- Quality control and denoising with DADA2
- Merge datasets from multiple sequencing runs
- Open-reference clustering (97% similarity) against *Trebouxia* reference database (Muggia et al. 2020)
- Filter low-abundance OTUs (<0.10% relative abundance threshold)

### 2. Community Analysis (R)
- Rename unknown OTUs based on phylogenetic placement
- Aggregate duplicate OTUs
- Calculate alpha diversity (Shannon index) with rarefaction
- Statistical testing (Wilcoxon signed-rank test for core vs. tip comparison)
- Generate publication figures (stacked bar plots, regression plots)
- Phylogenetic visualization with ggtree

## Software Requirements

### QIIME2
- QIIME2 version 2021.4
- Plugins: cutadapt, dada2, vsearch, diversity

### R Packages
```r
dplyr          # Data manipulation
tidyr          # Data tidying
ggplot2        # Data visualization
scales         # Scale functions for ggplot2
ggtree         # Phylogenetic tree visualization
treeio         # Tree data input/output
ape            # Phylogenetic analysis
metacoder      # Rarefaction and taxonomic analysis
vegan          # Ecological diversity analysis
```

## Input Data

### Required Files
1. **Raw sequencing data**: Illumina MiSeq paired-end FASTQ files
   - Fall 2020 samples
   - Spring 2022 samples
2. **Mapping files**: Sample metadata in QIIME2 manifest format
3. **Reference database**: *Trebouxia* ITS sequences (Muggia et al. 2020)
4. **Metadata**: Sample information (treatment, thallus type, etc.)

### Data Structure
**Metadata columns**:
- `sample_id`: Unique sample identifier
- `individual`: Lichen individual ID
- `thallus_type`: "core" or "tips" 
- `temp`: Temperature treatment (°C above ambient)
- `CO2`: CO2 treatment (450 or 900 ppm)

## Data Availability

The processed OTU table, metadata, and taxonomy assignments are available on Dryad:
[https://doi.org/10.5061/dryad.rbnzs7hcf](https://doi.org/10.5061/dryad.rbnzs7hcf)

Raw sequencing data (FASTQ files) are available upon request to the corresponding author.

Representative sequences for novel *Trebouxia* OTUs are deposited in GenBank (accession numbers OP879470-OP879475).

## Usage

### Running the QIIME2 Pipeline
```bash
# Adjust paths in the script to match your system
bash scripts/01_QIIME2_pipeline.sh
```

### Running the R Analysis
```bash
# From R console or RStudio
source("scripts/02_Trebouxia_analysis.R")
```

**Note**: You will need to update file paths in both scripts to match your local directory structure.

## Study Design

**Experimental Setup**: SPRUCE whole-ecosystem warming experiment
- Location: Marcell Experimental Forest, Minnesota, USA
- Treatments: 6 warming levels (0, +2.25, +4.5, +6.75, +9.0°C) × 2 CO2 levels (ambient, elevated)
- Sample types: Lichen thallus cores and tips
- Timepoint: August 2019 (4 years into experimental warming)

**Sampling Strategy**: 
- 5 *Evernia mesomorpha* thalli per enclosure
- DNA extracted separately from core and tip regions
- Targeted ITS1-ITS2 region for algal identification

## Key Results

- **OTU Discovery**: 12 *Trebouxia* OTUs identified (clades S, I, and A)
- **Temperature Effect**: Significant increase in clade 'I' proportion with warming (P = 0.0011, R² = 0.159)
- **Spatial Pattern**: Core regions had significantly higher alpha diversity than tips (P = 0.0151)
- **Diversity Metrics**: Shannon diversity ranged from 0.125 (median tips) to 0.334 (median cores)

## Citation

If you use this code or analysis approach, please cite:

Meyer, A.R., et al. 2023. Climate warming causes photobiont degradation and carbon starvation in a boreal climate sentinel lichen. *American Journal of Botany* 110(2): e16114.

## Author

**Abigail R. Meyer**
- Email: meye2099@umn.edu
- ORCID: [0000-0001-6375-1897](https://orcid.org/0000-0001-6375-1897)
- Affiliation: Ecology, Evolution and Behavior, University of Minnesota

## License

Code is provided as-is for academic and research purposes. Please cite the original publication if using these methods.

## Acknowledgments

This research was conducted at the USDA Forest Service SPRUCE experiment site. Sequencing was performed at the University of Minnesota Genomics Center.
