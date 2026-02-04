# Workflow Overview: Trebouxia Photobiont Metabarcoding Analysis

## Pipeline Summary

This document provides a detailed overview of the bioinformatics and statistical workflow used to analyze *Trebouxia* photobiont communities in *Evernia mesomorpha* lichens from the SPRUCE climate manipulation experiment.

## Workflow Diagram

```
Raw FASTQ Files (Illumina MiSeq 2×300bp)
          ↓
    [QIIME2 Import]
          ↓
    [Primer Removal] ← ITS1T/ITS2T primers
          ↓
    [Quality Control & DADA2 Denoising] ← Truncate at 240bp
          ↓
    [Merge Sequencing Runs] ← F20 + S22
          ↓
    [Open-Reference Clustering] ← Muggia et al. 2020 reference
          ↓                            97% similarity threshold
    [Filter Low-Abundance OTUs] ← 0.10% rule (1047 reads)
          ↓
    [Export OTU Table]
          ↓
═══════════════════════════════════════════════════════════
          ↓
    [R Analysis]
          ↓
    ┌─────┴─────────────────────┐
    ↓                            ↓
[Rename Unknown OTUs]    [Aggregate Duplicates]
  (RAxML placement)         (Sum by name)
    ↓                            ↓
    └─────┬─────────────────────┘
          ↓
    [Add Clade Assignments] ← I, S, or A based on name
          ↓
    ┌─────┴──────┐
    ↓            ↓
[Community      [Alpha Diversity]
 Analysis]       Analysis
    ↓                 ↓
[Figure 8A]      [Rarefaction]
    ↓                 ↓
[Figure 8B]      [Shannon Index]
                      ↓
                 [Wilcoxon Test]
                  Core vs. Tips
```

## Detailed Methods

### Phase 1: Sequencing and Quality Control (QIIME2)

#### 1.1 Amplicon Sequencing
- **Platform**: Illumina MiSeq
- **Chemistry**: 2×300bp paired-end
- **Target**: ITS1-ITS2 region (full length)
- **Primers**: 
  - Forward (ITS1T): `GGAAGGATCATTGAATCTATCGT`
  - Reverse (ITS2T): `GGTTCGCTCGCCGCTACTA`
- **Specificity**: Algal-specific (Kroken and Taylor 2000)

#### 1.2 Data Import
- QIIME2 artifact format (.qza)
- Manifest-based import (SingleEndFastqManifestPhred33V2)
- Separate imports for two sequencing runs

#### 1.3 Primer Trimming
- Tool: cutadapt (via QIIME2 plugin)
- Strategy: Remove forward primer from 5' end
- Error rate: 0 (no mismatches allowed)
- Purpose: Ensure amplicon uniformity

#### 1.4 Quality Filtering and Denoising
- Tool: DADA2 (Callahan et al. 2016)
- Truncation: 240bp (based on quality score inspection)
- Quality threshold: Median PHRED ≥20
- Features:
  - Error modeling and correction
  - Chimera detection and removal
  - Generation of exact Amplicon Sequence Variants (ASVs)

#### 1.5 Dataset Merging
- Combined ASV tables from both sequencing runs
- Merged representative sequences
- Maintained sample metadata throughout

#### 1.6 Taxonomic Assignment
**Method**: Open-reference clustering (vsearch)
- **Reference database**: Muggia et al. (2020) *Trebouxia* ITS sequences
- **Similarity threshold**: 97% (operational species-level)
- **Strategy**: 
  1. Match ASVs to reference sequences
  2. Cluster novel ASVs de novo
  3. Create expanded reference set
- **Advantage**: Captures both known and unknown diversity

#### 1.7 Quality Filtering
- **Threshold**: Minimum 1047 reads per OTU across all samples
- **Rationale**: 0.10% abundance rule (Reitmeier et al. 2021)
- **Purpose**: Remove spurious sequences from sequencing errors
- **Effect**: Conservative filtering retains biological diversity

### Phase 2: Phylogenetic Placement (External to Pipeline)

#### 2.1 Unknown OTU Identification
**Problem**: 7 OTUs (UNK_01 through UNK_07) didn't match reference at 97%

**Solution**: Phylogenetic placement in Muggia et al. (2020) tree
- Tool: RAxML v8.2.12 (Stamatakis 2014)
- Method: Evolutionary Placement Algorithm (EPA)
- Model: GTRGAMMA
- Result: Placement of unknowns in clades A, I, or S

**Outcomes**:
- UNK_01 → I02 (Clade I, impressa/gelatinosa)
- UNK_02 → Removed (failed to align)
- UNK_03 → A13_1 (Clade A, arboricola/gigantea)
- UNK_04 → S01 (Clade S, simplex/jamesii)
- UNK_05 → I18 (Clade I)
- UNK_06 → A13_2 (Clade A)
- UNK_07 → A04 (Clade A)

### Phase 3: Community Analysis (R)

#### 3.1 Data Preparation
**OTU Table Processing**:
1. Rename unknown OTUs based on phylogenetic placement
2. Aggregate duplicate OTUs by first 3 characters of name
3. Merge with sample metadata

**Data Transformation**:
- Wide to long format conversion (gather/pivot_longer)
- Calculate clade-level summaries
- Compute proportional abundances

#### 3.2 Community Composition Visualization (Figure 8A)

**Stacked Bar Plot**:
- X-axis: Thallus type (core vs. tips)
- Y-axis: Proportion of reads (0-100%)
- Facets: Temperature treatment (0, 2, 4, 6, 8, 10°C)
- Fill: OTU identity with phylogenetically-informed colors
  - Clade A: Brown/yellow gradient
  - Clade I: Blue gradient
  - Clade S: Green gradient

**Phylogram (Legend)**:
- Rooted phylogenetic tree
- Branch lengths ignored (cladogram)
- Tips labeled with OTU codes
- Purpose: Show evolutionary relationships among OTUs

#### 3.3 Temperature Effects Analysis (Figure 8B)

**Research Question**: Does clade composition shift with warming?

**Statistical Approach**:
- Response variable: Proportion of clade 'I' per sample
- Predictor: Temperature treatment (continuous)
- Method: Linear regression
- Formula: `clade_prop ~ temp`

**Results**:
- Slope: 0.0196 (±SE)
- Intercept: 0.0007
- R²: 0.159
- p-value: 0.0011 (significant)
- Interpretation: ~2% increase in clade I per °C warming

**Visualization**:
- Scatter plot with regression line
- X-axis: Temperature (categorical for display)
- Y-axis: Proportion clade 'I' (%)

#### 3.4 Alpha Diversity Analysis

**Research Question**: Do core and tip regions differ in photobiont diversity?

**Rarefaction**:
- Tool: metacoder package (Foster et al. 2017)
- Depth: 1790 reads (minimum sample total)
- Purpose: Normalize for sequencing depth differences
- Maintains taxonomic structure during rarefaction

**Diversity Metric**:
- Shannon diversity index (vegan package)
- Formula: H' = -Σ(pi × ln(pi))
  - pi = proportion of reads belonging to OTU i
- Accounts for both richness and evenness

**Statistical Test**:
- Method: Wilcoxon signed-rank test (non-parametric)
- Design: Paired (core and tip from same individual)
- Null hypothesis: No difference in median diversity
- Alternative: Two-sided test

**Results**:
- Median Shannon (core): 0.334
- Median Shannon (tips): 0.125
- p-value: 0.0151 (significant)
- Effect size: ~2.7× higher diversity in cores

**Interpretation**: 
- Tips experience more environmental stress (UV, desiccation)
- Cores maintain more stable conditions
- May relate to position-specific photobiont mortality

## Key Computational Decisions

### 1. Single-End vs. Paired-End Processing
**Decision**: Process as single-end (R1 only)  
**Rationale**: After quality trimming, R1 and R2 reads didn't overlap sufficiently for merging

### 2. Truncation Length (240bp)
**Decision**: Truncate at 240bp  
**Rationale**: Balance between:
- Retaining maximum read length
- Maintaining quality (median PHRED >20)
- Sufficient length for taxonomic resolution

### 3. Clustering Threshold (97%)
**Decision**: 97% similarity for OTU clustering  
**Rationale**:
- Standard for species-level diversity
- Balances oversplitting vs. lumping
- Consistent with reference database construction

### 4. Abundance Filtering (1047 reads)
**Decision**: Minimum 1047 total reads per OTU  
**Rationale**: 
- Represents 0.10% of total dataset
- Reduces false-positive rare taxa
- Conservative approach from Reitmeier et al. (2021)

### 5. Rarefaction Depth (1790 reads)
**Decision**: Rarefy to 1790 reads per sample  
**Rationale**:
- Minimum sample depth in dataset
- Maximizes sample retention
- Required for unbiased alpha diversity comparison

## Software Versions and Environment

### QIIME2 Environment
- **Core**: QIIME2 2021.4
- **Plugins**:
  - cutadapt
  - dada2
  - vsearch
  - diversity

### R Environment
Key packages (versions may vary):
- dplyr, tidyr, ggplot2 (tidyverse ecosystem)
- vegan (community ecology)
- metacoder (taxonomic analysis)
- ggtree, treeio, ape (phylogenetics)

## Quality Control Checkpoints

Throughout the pipeline, quality is assessed at multiple stages:

1. **Post-import**: Raw read quality visualizations
2. **Post-trimming**: Primer removal efficiency
3. **Post-DADA2**: Denoising statistics (reads retained, chimeras removed)
4. **Post-clustering**: Rarefaction curves (sampling depth assessment)
5. **Post-filtering**: OTU table summary statistics
6. **Pre-analysis**: Rarefaction curves for alpha diversity

## Reproducibility Notes

### File Paths
All file paths in scripts are relative to local directory structure. Users must update:
- Working directory (`setwd()` in R)
- Input file paths in QIIME2 commands
- Output directory locations

### Random Seeds
- DADA2: Uses internal error modeling (reproducible with same version)
- Rarefaction: Should set seed for reproducible subsampling
  - Can add: `set.seed(123)` before rarefaction in R

### Platform Considerations
- QIIME2 scripts designed for Unix/Linux/MacOS
- R scripts platform-independent
- Some file path separators may need adjustment for Windows

## References

Key methods papers cited in workflow:

- **Bolyen et al. 2019**: QIIME2 platform
- **Callahan et al. 2016**: DADA2 denoising algorithm  
- **Reitmeier et al. 2021**: Abundance filtering thresholds
- **Muggia et al. 2020**: *Trebouxia* reference phylogeny
- **Stamatakis 2014**: RAxML phylogenetic placement
- **Foster et al. 2017**: Metacoder package
- **Oksanen et al. 2020**: vegan package

Full references available in main README and manuscript.
