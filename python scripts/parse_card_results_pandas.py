"""
Parse CARD RGI results files and compare with clinical metadata using pandas.

CARD (Comprehensive Antibiotic Resistance Database) is a curated resource providing
reference sequences and detection models for antimicrobial resistance (AMR) genes.
The Resistance Gene Identifier (RGI) tool predicts resistome(s) from protein or 
nucleotide data based on homology and SNP models.

This script:
1. Parses CARD RGI output files to identify penam (beta-lactam) antibiotic resistance genes
2. Filters for high-confidence hits (Strict/Perfect matches only)
3. Compares detected resistance genes against clinical beta-lactamase test results
4. Calculates sensitivity and specificity of the metagenomic approach
5. Identifies genes that appear exclusively in clinically positive or negative samples
"""
import os
import pandas as pd

# Configuration
RESULT_DIR = "/path/to/directory"
METADATA_FILE = "metadata.csv"

def parse_card_results(result_dir):
    """
    Parse CARD output files and return DataFrame.
    
    Extracts AMR gene hits from CARD RGI output files, filtering for:
    - Strict or Perfect hit types (high confidence matches)
    - Penam drug class (beta-lactam antibiotics)
    
    Returns DataFrame with columns: sample_id, hit_type, gene, identity, 
    ref_length, drug_class, resistance_mechanism, gene_family
    """
    records = []
    
    # Walk through directory tree to find all CARD output files
    for root, _, files in os.walk(result_dir):
        for file in files:
            if not file.endswith("cardoutput.txt"):
                continue
            
            # Extract sample ID from filename (format: XXXX_SXX_cardoutput.txt)
            sample_id = "_".join(file.split("_")[:2])
            
            with open(os.path.join(root, file), 'r') as f:
                next(f)  # Skip header line
                for line in f:
                    fields = line.strip().split("\t")
                    
                    # Filter for high-confidence hits (Strict/Perfect) with penam in drug class
                    # fields[5] = hit type, fields[14] = drug class
                    if fields[5] in ("Strict", "Perfect") and "penam" in fields[14].split("; "):
                        records.append({
                            'sample_id': sample_id,
                            'hit_type': fields[5],              # Strict or Perfect
                            'gene': fields[8],                  # Best hit ARO (gene name)
                            'identity': float(fields[9]),       # Percent identity to reference
                            'ref_length': float(fields[20]),    # Length of reference sequence
                            'drug_class': fields[14],           # Drug class
                            'resistance_mechanism': fields[15], # Resistance mechanism
                            'gene_family': fields[16],          # AMR gene family
                        })
    
    return pd.DataFrame(records)

def main():
    # Parse CARD results into DataFrame
    df = parse_card_results(RESULT_DIR)
    
    # ============================================================================
    # SECTION 1: Overall gene and gene family counts across all samples
    # ============================================================================
    
    # Count how many unique samples each gene appears in
    print("GENES--------------------")
    gene_counts = df.groupby('gene')['sample_id'].nunique().sort_values(ascending=False)
    for gene, count in gene_counts.items():
        print(f"{gene}: {count}")
    
    # Count how many unique samples each gene family appears in
    print("\nGENE FAMILIES ----------------")
    family_counts = df.groupby('gene_family')['sample_id'].nunique().sort_values(ascending=False)
    for family, count in family_counts.items():
        print(f"{family}: {count}")
    
    # Total number of samples with at least one resistance gene detected
    samples_with_resistance = df['sample_id'].nunique()
    print(f"\nTotal samples with 1 or more resistance genes: {samples_with_resistance}")
    
    # ============================================================================
    # SECTION 2: Compare mNGS predictions to clinical beta-lactamase test results
    # ============================================================================
    
    # Load clinical metadata (includes beta-lactamase test results)
    metadata = pd.read_csv(METADATA_FILE)
    
    # Create summary of which samples have resistance gene hits
    samples_with_hits = df.groupby('sample_id').size().reset_index(name='hit_count')
    
    # Merge with metadata (left join keeps all samples from metadata)
    merged = metadata.merge(samples_with_hits, left_on='Filename', right_on='sample_id', how='left')
    merged['hit_count'] = merged['hit_count'].fillna(0)  # Samples with no hits get 0
    merged['has_resistance'] = merged['hit_count'] > 0
    
    # Calculate sensitivity and specificity
    # Clinical status: 1 = positive for beta-lactamase, 2 = negative
    positive_samples = merged[merged['Hflu_beta_lactamase'] == 1]
    negative_samples = merged[merged['Hflu_beta_lactamase'] == 2]
    
    # Sensitivity: proportion of clinically positive samples correctly identified by mNGS
    true_pos = positive_samples['has_resistance'].sum()
    total_pos = len(positive_samples)
    
    # Specificity: proportion of clinically negative samples correctly identified by mNGS
    true_neg = (~negative_samples['has_resistance']).sum()
    total_neg = len(negative_samples)
    
    sensitivity = (true_pos / total_pos) * 100 if total_pos > 0 else 0
    specificity = (true_neg / total_neg) * 100 if total_neg > 0 else 0
    
    print(f"\nSensitivity = {sensitivity:.2f}% and Specificity = {specificity:.2f}%")
    
    # ============================================================================
    # SECTION 3: Analyze which genes appear in positive vs negative samples
    # ============================================================================
    
    # Add clinical status to main dataframe for filtering
    df = df.merge(metadata[['Filename', 'Hflu_beta_lactamase']], 
                  left_on='sample_id', right_on='Filename', how='left')
    df = df.rename(columns={'Hflu_beta_lactamase': 'clinical_status'})
    
    # Count genes in clinically positive samples (status = 1)
    print("\nGenes in clinically positive samples:\n-------------------------")
    pos_genes = df[df['clinical_status'] == 1].groupby('gene')['sample_id'].nunique().sort_values(ascending=False)
    for gene, count in pos_genes.items():
        print(f"{gene}: {count}")
    
    # Count genes in clinically negative samples (status = 2)
    print("\nGenes in clinically negative samples:\n-------------------------")
    neg_genes = df[df['clinical_status'] == 2].groupby('gene')['sample_id'].nunique().sort_values(ascending=False)
    for gene, count in neg_genes.items():
        print(f"{gene}: {count}")
    
    # Identify genes that appear exclusively in one group or the other
    positive_only = set(pos_genes.index) - set(neg_genes.index)
    negative_only = set(neg_genes.index) - set(pos_genes.index)
    
    print(f"\nNegative only genes: {list(negative_only)}")
    print(f"Positive only genes: {list(positive_only)}")
    
    # ============================================================================
    # SECTION 4: Calculate statistics for genes in clinically positive samples
    # ============================================================================
    
    # For each gene, calculate average, median, and minimum identity and reference length
    print("\nPositive gene statistics:")
    pos_stats = df[df['clinical_status'] == 1].groupby('gene').agg({
        'identity': ['mean', 'median', 'min'],
        'ref_length': ['mean', 'median', 'min']
    }).round(2)
    
    for gene in pos_stats.index:
        id_avg, id_med, id_min = pos_stats.loc[gene, 'identity']
        ref_avg, ref_med, ref_min = pos_stats.loc[gene, 'ref_length']
        print(f"{gene} - ID: avg={id_avg:.2f}, med={id_med:.2f}, min={id_min:.2f} | "
              f"Ref: avg={ref_avg:.2f}, med={ref_med:.2f}, min={ref_min:.2f}")

if __name__ == "__main__":
    main()
