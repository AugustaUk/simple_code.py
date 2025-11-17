def dna_to_protein(dna_sequence):
    """
    Translate a DNA sequence to a protein sequence using the standard genetic code.

    Args:
        dna_sequence (str): A string of DNA nucleotides (A, T, G, C)

    Returns:
        str: The translated protein sequence
    """
    # Genetic code dictionary (i got the genetic code from here, https://www.hgmd.cf.ac.uk/docs/cd_amino.html)
    genetic_code = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    # Convert to uppercase and remove any whitespace
    dna_sequence = dna_sequence.upper().replace(' ', '').replace('\n', '')

    # Ensure the sequence length is divisible by 3
    if len(dna_sequence) % 3 != 0:
        # Add 'N' nucleotides to make it divisible by 3
        padding_length = 3 - (len(dna_sequence) % 3)
        dna_sequence += 'N' * padding_length

    protein_sequence = ""

    # Translate each codon
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]

        # If codon is incomplete or not recognized, treat as stop codon
        if len(codon) < 3 or codon not in genetic_code:
            protein_sequence += '*'
        else:
            protein_sequence += genetic_code[codon]

    return protein_sequence


def hamming_distance(str1, str2):
    """
    Calculate the Hamming distance between two strings.
    The Hamming distance is the number of positions at which the corresponding symbols are different.

    Args:
        str1 (str): First string
        str2 (str): Second string

    Returns:
        int: The Hamming distance between the two strings
    """
    # Pad the shorter string with extra characters to make them equal length
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)
    str2 = str2.ljust(max_len)

    # Calculate the Hamming distance
    distance = 0
    for i in range(max_len):
        if str1[i] != str2[i]:
            distance += 1

    return distance


# Example:
if __name__ == "__main__":
    # Process team members  read DNA sequences from their respective FASTA files
    team_members = [
        {"name": "Teodora", "username": "Teodora", "fasta_file": "Teodora_sequence.fasta"},
        {"name": "Augusta Ilechukwu", "username": "Augusta Ilechukwu ", "fasta_file": "Augusta_sequence.fasta"},
        {"name": "Sultan Alkhazraji", "username": "Sultan Alkhazraji ", "fasta_file": "Sultan_sequence.fasta"},
        {"name": "Jovana Aleksic", "username": "Jovana Aleksic", "fasta_file": "Jovana_sequence.fasta"},
        {"name": "Auwal Ibrahim", "username": "Auwal Ibrahim", "fasta_file": "Auwal_sequence.fasta"}
    ]

    for member in team_members:
        print(f"\n--- {member['name']} ---")
        slack_username = member['username'].strip()
        # Generate bluesky handle if not in CSV (synthesize one)
        bluesky_handle = f"{slack_username.lower().replace(' ', '')}.bsky.social"
        hamming_dist = hamming_distance(slack_username, bluesky_handle)
        print(f"Slack username: {slack_username}")
        print(f"Bluesky handle: {bluesky_handle}")
        print(f"Hamming distance between '{slack_username}' and '{bluesky_handle}': {hamming_dist}")

        # Read DNA sequence from the FASTA file for this member
        try:
            with open(member['fasta_file'], 'r') as file:
                lines = file.readlines()
                # Skip the first line (because its a header) and connect the rest
                dna_sequence = ''.join(lines[1:]).replace('\n', '').replace(' ', '')

            protein = dna_to_protein(dna_sequence)
            print(f"DNA sequence length: {len(dna_sequence)}")
            print(f"Protein sequence: {protein}")

        except FileNotFoundError:
            print(f"DNA file for {member['name']} not found")



Gene Expression Analysis

# a. Heatmap
"""
    Plot a clustered heatmap of top differentially expressed genes.

    Parameters:
    - file_path: str, path to CSV file (genes x samples)
    - groups: dict, e.g., {"HBR": ["HBR_1","HBR_2"], "UHR": ["UHR_1","UHR_2"]}
    - top_n: int, number of top genes to select
    - cmap: str, color palette
    - normalize_rows: bool, whether to normalize rows (genes) to 0-1
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Set modern theme and color palette
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['font.size'] = 12

# Define color scheme
colors = {
    "Upregulated": "#2E8B57",   # Forest Green
    "Downregulated": "#DC143C",  # Crimson
    "Not significant": "#696969", # Dim Gray
    "Malignant": "#FF6B6B",     # Coral
    "Benign": "#4ECDC4"          # Turquoise
}

# Load normalized expression data
gene_expression_df = pd.read_csv('/content/normalized_counts_for_HBR_and_UHR.csv.txt', index_col=0)

# Define sample groups
sample_groups = {
    "HBR": ["HBR_1", "HBR_2", "HBR_3"],
    "UHR": ["UHR_1", "UHR_2", "UHR_3"]
}

# Compute mean expression per group
gene_expression_df["HBR_mean"] = gene_expression_df[sample_groups["HBR"]].mean(axis=1)
gene_expression_df["UHR_mean"] = gene_expression_df[sample_groups["UHR"]].mean(axis=1)

# Compute absolute fold change
gene_expression_df["fold_change"] = np.abs(gene_expression_df["HBR_mean"] - gene_expression_df["UHR_mean"])

# Select top 20 genes based on fold change
top_gene_subset = gene_expression_df.nlargest(20, "fold_change")
heatmap_matrix = top_gene_subset[sample_groups["HBR"] + sample_groups["UHR"]]

# Generate clustered heatmap
sns.set(font_scale=1)
clustered_heatmap = sns.clustermap(
    heatmap_matrix,
    cmap="RdYlBu_r",     # Red-Yellow-Blue reversed (better for fold changes)
    standard_scale=1,   # normalize rows to 0-1
    row_cluster=True,   # cluster genes
    col_cluster=True,   # cluster samples
    linewidths=0.5,
    figsize=(12, 10),
    cbar_pos=(0.02, 0.8, 0.03, 0.18)  # Adjust colorbar position
)

# Set axis labels
clustered_heatmap.ax_heatmap.set_xlabel("Samples", fontsize=14, fontweight='bold')
clustered_heatmap.ax_heatmap.set_ylabel("Genes", fontsize=14, fontweight='bold')
plt.show()

print('---------------------------------------------------------------------')
print('Volcano Plot of Differentially Expressed Genes (DEGs)')

# Load differential expression results
deg_results = pd.read_csv('/content/differential_expression_results.csv.txt')

# Categorize significance
def categorize_significance(row):
    if row["PAdj"] < 0.05 and row["log2FoldChange"] > 1:
        return "Upregulated"
    elif row["PAdj"] < 0.05 and row["log2FoldChange"] < -1:
        return "Downregulated"
    else:
        return "Not significant"

deg_results["category"] = deg_results.apply(categorize_significance, axis=1)

# Create volcano plot
plt.figure(figsize=(10, 6))
for category in ["Upregulated", "Downregulated", "Not significant"]:
    subset = deg_results[deg_results["category"] == category]
    plt.scatter(
        subset["log2FoldChange"],
        -np.log10(subset["PAdj"]),
        c=colors[category],
        label=category,
        alpha=0.7,
        edgecolors="white",
        linewidth=0.5,
        s=60
    )

# Add vertical lines
plt.axvline(x=1, color='black', linestyle='--', linewidth=1.2)
plt.axvline(x=-1, color='black', linestyle='--', linewidth=1.2)

# Labels and title
plt.xlabel("log2 Fold Change", fontsize=14, fontweight='bold')
plt.ylabel("-log10(adjusted p-value)", fontsize=14, fontweight='bold')
plt.title("Volcano Plot of DEGs (HBR vs UHR)", fontsize=16, fontweight='bold')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print('---------------------------------------------------------------------')
print('Breast Cancer Data Exploration - Scatter Plot of Radius vs Texture')

# Load breast cancer dataset
bc_data = pd.read_csv('/content/breast_cancer_wisconsin.csv.txt')

# Map diagnosis to colors
diagnosis_colors = {'M': colors["Malignant"], 'B': colors["Benign"]}

plt.figure(figsize=(10, 7))
plt.scatter(
    bc_data['radius_mean'],
    bc_data['texture_mean'],
    c=bc_data['diagnosis'].map(diagnosis_colors),
    alpha=0.7,
    edgecolor='white',
    linewidth=0.5
)

# Add legends
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Malignant', markerfacecolor=colors["Malignant"], markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Benign', markerfacecolor=colors["Benign"], markersize=8)
]
plt.legend(handles=legend_elements, title="Diagnosis", title_fontsize=12, fontsize=10)

# Labels and title
plt.xlabel('Radius Mean', fontsize=14, fontweight='bold')
plt.ylabel('Texture Mean', fontsize=14, fontweight='bold')
plt.title('Breast Cancer: Radius vs Texture', fontsize=16, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print('---------------------------------------------------------------------')
print('Breast Cancer Data Exploration - Correlation Heatmap of Key Features')

# Select key features
key_features = ['radius_mean', 'texture_mean', 'perimeter_mean',
                'area_mean', 'smoothness_mean', 'compactness_mean']

# Compute correlation matrix
feature_subset = bc_data[key_features]
correlation_matrix = feature_subset.corr()

plt.figure(figsize=(10, 8))
sns.heatmap(
    correlation_matrix,
    annot=True,
    cmap="RdBu_r",
    fmt=".2f",
    linewidths=0.5,
    square=True,
    cbar_kws={'shrink': 0.8}
)

plt.title("Correlation Heatmap of Key Breast Cancer Features", fontsize=16, fontweight='bold')
plt.tight_layout()
plt.show()

print('---------------------------------------------------------------------')
print('Breast Cancer Data Exploration - Scatter Plot of Compactness vs Smoothness')

# Map diagnosis to colors
point_colors = bc_data['diagnosis'].map(diagnosis_colors)

plt.figure(figsize=(10, 7))
plt.scatter(
    bc_data['smoothness_mean'],       # X-axis
    bc_data['compactness_mean'],      # Y-axis
    c=point_colors,
    alpha=0.7,
    edgecolor='white',
    linewidth=0.5
)

# Gridlines
plt.grid(True, linestyle='--', alpha=0.5)

# Axis labels and title
plt.xlabel('Smoothness Mean', fontsize=14, fontweight='bold')
plt.ylabel('Compactness Mean', fontsize=14, fontweight='bold')
plt.title('Breast Cancer: Compactness vs Smoothness', fontsize=16, fontweight='bold')

# Legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Malignant', markerfacecolor=colors["Malignant"], markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Benign', markerfacecolor=colors["Benign"], markersize=8)
]
plt.legend(handles=legend_elements, title="Diagnosis", title_fontsize=12, fontsize=10)
plt.tight_layout()
plt.show()

print('---------------------------------------------------------------------')
print('Breast Cancer Data Exploration - KDE of Area Mean by Diagnosis')

plt.figure(figsize=(10, 7))

# KDE plots
sns.kdeplot(
    data=bc_data[bc_data['diagnosis']=='M'],
    x='area_mean',
    fill=True,
    color=colors["Malignant"],
    alpha=0.6,
    label='Malignant',
    linewidth=2
)

sns.kdeplot(
    data=bc_data[bc_data['diagnosis']=='B'],
    x='area_mean',
    fill=True,
    color=colors["Benign"],
    alpha=0.6,
    label='Benign',
    linewidth=2
)

# Labels and title
plt.xlabel('Area Mean', fontsize=14, fontweight='bold')
plt.ylabel('Density', fontsize=14, fontweight='bold')
plt.title('Distribution of Area Mean by Diagnosis', fontsize=16, fontweight='bold')
plt.legend(title='Diagnosis', title_fontsize=12, fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
