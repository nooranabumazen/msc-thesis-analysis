# Metatranscriptomic Analysis of Pediatric Acute Sinusitis

[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs13073--025--01447--3-blue)](https://doi.org/10.1186/s13073-025-01447-3)
[![PMID](https://img.shields.io/badge/PMID-40098147-green)](https://pubmed.ncbi.nlm.nih.gov/40098147/)
[![View Analysis](https://img.shields.io/badge/View-Analysis-orange)](https://nooranabumazen.github.io/msc-thesis-analysis/)

R code and visualizations from my MSc thesis, published in *Genome Medicine* (2025).

## 📄 Publication

**Metatranscriptomic profiling reveals pathogen and host response signatures of pediatric acute sinusitis and upper respiratory infection**

*Doxey A.C., Abu Mazen N. et al., Genome Medicine (2025)*

[Read the full paper →](https://doi.org/10.1186/s13073-025-01447-3)

## 🔬 Project Overview

This repository contains the R code and analysis for investigating pathogen detection and host immune responses in pediatric acute sinusitis using RNA-sequencing. The study analyzes nasopharyngeal samples from 221 children to distinguish bacterial from viral upper respiratory infections.

### Key Findings

- **RNA-seq pathogen detection:** 87%/81% sensitivity/specificity for 3 bacterial pathogens; 86%/92% for 12 viral pathogens
- **Novel pathogen discovery:** Identified 22 additional pathogens not in the clinical panel, including 4 human coronaviruses
- **Diagnostic gap filled:** Found pathogens in 58% (11/19) of cases where culture and qRT-PCR were negative
- **Host-response signatures:** Discovered 548 bacterial-specific and 273 viral-specific genes (p<0.001)
- **Predictive value:** Host gene expression predicts infection type with AUC 0.77 (bacterial) and 0.88 (viral)
- **Viral genomes reconstructed:** 196 viral genomes including novel coronavirus, RSV, and enterovirus strains

## 📊 Analysis Contents

### Section 1: Bacterial Pathogen Detection
- RNA-seq validation against culture tests for MCAT, SPN, HFLU
- ROC curve analysis with sensitivity/specificity calculations
- Heatmaps and boxplots showing concordance between methods

### Section 2: Viral Pathogen Detection & Discovery
- Detection of 12 respiratory viruses by RNA-seq vs. qRT-PCR
- **Section 2.5:** Exploratory analysis identifying 22 additional pathogens:
  - 4 human coronaviruses (NL63, HKU1, 229E, Betacoronavirus 1)
  - Atypical bacteria (Chlamydia pneumoniae, Mycoplasma pneumoniae)
  - Solved 58% of clinically "negative" cases

### Section 3: Host Gene Expression Analysis
- **Section 3.1-3.2:** Differential expression analysis (DESeq2) comparing bacterial vs. viral infections
- **Section 3.3:** PCA plots showing:
  - PC1 (70% variance) = overall immune activation/infection severity
  - PC2 (9% variance) = bacterial vs. viral infection type
- Volcano plots identifying 821 differentially expressed genes
- Expression patterns of key immune markers (CXCL10, IL1B, S100A9, IFI27, etc.)

### Section 4: Host Response Correlates with Pathogen Load
- Z-score-based response scores for bacterial and viral signatures
- Correlation analysis: bacterial response score vs. pathogen abundance (r=0.50)
- Heatmaps ordered by pathogen abundance showing gene expression gradients
- Permutation testing validating specificity of bacterial response to known pathogens

### Section 5: Diagnostic Potential
- ROC analysis for predicting infection type from host response alone
- Response score boxplots stratified by pathogen abundance categories
- Comparison across "Both low", "Viral high", "Bacterial high", "Both high" groups

## 🛠️ Technologies & Methods

**R packages:** DESeq2, tximport, pheatmap, pROC, ggplot2, EnhancedVolcano, reshape2, RColorBrewer, patchwork

**Methods:** 
- Kraken2 taxonomic classification (pathogen detection)
- Salmon transcript quantification (host gene expression)
- DESeq2 differential expression with batch correction
- Variance-stabilizing transformation for visualization
- ROC curve analysis for diagnostic performance
- Principal component analysis (PCA)


## 🚀 Viewing the Analysis

**Live analysis:** [https://nooranabumazen.github.io/msc-thesis-analysis/](https://nooranabumazen.github.io/msc-thesis-analysis/)

The rendered HTML includes:
- Complete R code for reproducibility
- Figures and visualizations
- Statistical outputs and performance metrics
- Workflow from data processing to publication figures

## 💡 Key Insights

### Why This Matters

**Clinical Problem:** Distinguishing bacterial sinusitis (needs antibiotics) from viral infections (antibiotics unnecessary) is challenging with current diagnostic methods.

**Our Solution:** 
1. **RNA-seq as comprehensive diagnostic:** Detects pathogens missed by culture/qRT-PCR (coronaviruses, atypical bacteria)
2. **Host response profiling:** Gene expression patterns differentiate bacterial from viral infections
3. **Dual approach:** Combining pathogen detection + host response provides robust classification

### Impact

- **Reduced diagnostic failures:** 58% of "negative" cases had identifiable pathogens by RNA-seq
- **Antibiotic stewardship:** Better bacterial/viral distinction could reduce unnecessary prescriptions
- **Emerging pathogen surveillance:** Identified novel viral strains and expanded pathogen repertoire

## 📖 Citation

If you use this work, please cite:

```
Doxey AC, Abu Mazen N, Homm M, et al. 
Metatranscriptomic profiling reveals pathogen and host response signatures of pediatric acute sinusitis and upper respiratory infection. 
Genome Med. 2025 Mar 17;17(1):22. 
doi: 10.1186/s13073-025-01447-3.
```

## 👤 Author

**Noora Abu Mazen**  
MSc Graduate, University of Waterloo  
[GitHub](https://github.com/nooranabumazen) | [LinkedIn](https://www.linkedin.com/in/nooran-abu-mazen-9a0036163/)
Email: nooranabumazen@gmail.com

## 📜 License & Ethics

This work was completed as part of my MSc thesis under the supervision of Dr. Andrew Doxey with ethics approval from the University of Waterloo Research Ethics Board (Application #45063) and the University of Pittsburgh IRB. Patient data is not included in this repository to protect privacy.

---

*Research conducted in the Doxey Lab at the University of Waterloo in collaboration with UPMC Children's Hospital of Pittsburgh*
