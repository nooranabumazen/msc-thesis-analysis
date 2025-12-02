# Python Analysis Scripts from MSc Thesis

Python scripts from my MSc thesis analyzing pediatric acute sinusitis and upper respiratory infections using RNA-sequencing data. These scripts complement the [main R analysis repository](https://github.com/nooranabumazen/msc-thesis-analysis).

**Publication:** Doxey A.C., Abu Mazen N, Homm M, et al. Metatranscriptomic profiling reveals pathogen and host response signatures of pediatric acute sinusitis and upper respiratory infection. *Genome Med.* 2025 Mar 17;17(1):22. [doi: 10.1186/s13073-025-01447-3](https://doi.org/10.1186/s13073-025-01447-3)

---

## Overview

This repository contains Python scripts for:
1. **Antimicrobial resistance (AMR) gene detection** - Parsing CARD RGI results to identify beta-lactam resistance genes and validate against clinical tests
2. **Machine learning classifier** - Predicting bacterial vs. non-bacterial infection from host gene expression (exploratory analysis, not included in final publication)

Both analyses used nasopharyngeal swab samples from 221 children with acute sinusitis or upper respiratory infections.

---

## Scripts

### 1. `parse_card_results_pandas.py`

**Purpose:** Parse CARD (Comprehensive Antibiotic Resistance Database) RGI output files to identify penam (beta-lactam) antibiotic resistance genes and compare with clinical beta-lactamase test results.

**Key Features:**
- Filters for high-confidence hits (Strict/Perfect matches only)
- Focuses on penam drug class (beta-lactam antibiotics)
- Calculates sensitivity and specificity against clinical tests
- Identifies genes appearing exclusively in clinically positive or negative samples
- Generates statistics on gene identity and reference sequence coverage

**Input:**
- CARD RGI output files (`*cardoutput.txt`)
- Clinical metadata with beta-lactamase test results
- Differentially expressed genes list

**Output:**
- Gene and gene family counts across samples
- Sensitivity/specificity metrics
- Lists of genes in clinically positive vs. negative samples
- Gene statistics (identity %, reference length)

**Dependencies:**
```python
pandas
os
```

**Usage:**
```bash
python parse_card_results_pandas.py
```

**Configuration:** Edit these variables at the top of the script:
```python
RESULT_DIR = "/path/to/card-results"
METADATA_FILE = "metadata.csv"
```

---

### 2. `bacterial_infection_classifier.py`

**Purpose:** Build a machine learning classifier to predict bacterial vs. non-bacterial infection from host gene expression in nasal swab samples.

**Method:** Bagged Support Vector Machine (BagSVM) with nested cross-validation
- **Outer CV:** Model evaluation (10 iterations)
- **Inner CV:** Hyperparameter tuning (GridSearchCV) and feature selection (RFECV)
- **Feature selection:** Recursive feature elimination keeping 2-100 genes
- **Evaluation:** AUC-ROC on training, test, and held-out sets

**Note:** This classifier was developed as exploratory analysis during my MSc but was not included in the final publication. The published work (Section 4.2) used simpler bacterial/viral response scores based on differentially expressed genes, which achieved AUC 0.77 (bacterial) and 0.88 (viral).

**Input:**
- Gene counts matrix (samples × genes)
- Metadata with infection type labels
- List of differentially expressed genes

**Output:**
- Trained classifier models (`.joblib` files)
- Selected feature lists for each CV iteration
- Performance summary table with AUC-ROC scores

**Dependencies:**
```python
numpy
pandas
scikit-learn
matplotlib
joblib
```

**Usage:**
```bash
python bacterial_infection_classifier.py
```

**Configuration:** Edit these variables at the top of the script:
```python
RESULTS_PATH = "/path/to/output"
OUTPUT_PREFIX = "bacterial_infection_classifier"
MODE = "create"  # or "load" to load existing models
NUM_CV = 5
TEST_PROPORTION = 0.25
NUM_OUTER_CV = 10
```

---

## Background

### Antimicrobial Resistance Analysis

*Haemophilus influenzae* beta-lactamase production is a key resistance mechanism against penam antibiotics. We used CARD RGI to detect resistance genes from RNA-seq data and validated against clinical beta-lactamase tests. This analysis achieved:
- **87% sensitivity** for detecting *H. influenzae* beta-lactamase
- **High specificity** in distinguishing resistant from susceptible strains

### Host Response Classification

Host immune responses differ between bacterial and viral infections. This classifier attempted to predict infection type using only host gene expression patterns, independent of pathogen detection. While the classifier showed promising results, the final publication used a simpler approach with bacterial/viral response scores that were easier to interpret clinically.

---

## Related Work

**Main Analysis Repository:** [msc-thesis-analysis](https://github.com/nooranabumazen/msc-thesis-analysis)
- R code for pathogen detection, differential expression, and visualization
- Live analysis: https://nooranabumazen.github.io/msc-thesis-analysis/

**Key Findings from Published Work:**
- RNA-seq detected pathogens with 87%/81% sensitivity/specificity (bacterial) and 86%/92% (viral)
- Identified 22 additional pathogens not in clinical panel, including 4 human coronaviruses
- Found pathogens in 58% of cases where culture/qRT-PCR were negative
- Discovered 548 bacterial-specific and 273 viral-specific host genes
- Host gene expression predicted infection type with AUC 0.77 (bacterial) and 0.88 (viral)

---

## Citation

If you use these scripts, please cite:

```
Doxey A.C., Abu Mazen N, Homm M, et al.
Metatranscriptomic profiling reveals pathogen and host response signatures of 
pediatric acute sinusitis and upper respiratory infection.
Genome Med. 2025 Mar 17;17(1):22.
doi: 10.1186/s13073-025-01447-3
```

---

## Contact

**Nooran Abu Mazen**  
MSc Bioinformatics, University of Waterloo  
[GitHub](https://github.com/nooranabumazen) | [LinkedIn](https://www.linkedin.com/in/nooran-abu-mazen-9a0036163/)  
Email: nooranabumazen@gmail.com

---

## Acknowledgments

This work was completed as part of my MSc thesis under the supervision of Dr. Andrew Doxey at the University of Waterloo, in collaboration with UPMC Children's Hospital of Pittsburgh.

**Ethics Approval:** University of Waterloo Research Ethics Board (Application #45063) and University of Pittsburgh IRB

**Funding:** Natural Sciences and Engineering Research Council of Canada (NSERC)

---

## License

This code is provided for research and educational purposes. Patient data is not included to protect privacy.
