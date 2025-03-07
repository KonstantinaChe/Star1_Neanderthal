import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import AlignIO

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_alignment.py alignment.fasta")
        sys.exit(1)

    alignment_file = sys.argv[1]

    # Read the alignment using Biopython (assuming FASTA format)
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except Exception as e:
        print(f"Error reading alignment: {e}")
        sys.exit(1)

    # Build a dictionary mapping each sequence ID to a list of uppercase bases
    seq_dict = {record.description: list(str(record.seq).upper()) for record in alignment}

    # Ensure 'Starosele' is present in the alignment
    starosele_keys = [key for key in seq_dict.keys() if "Starosele" in key]
    if not starosele_keys:
        print("Error: 'Starosele' sequence not found in the alignment.")
        sys.exit(1)
    starosele_key = starosele_keys[0]

    # Create a DataFrame where each column is a sequence and each row is an alignment position
    df = pd.DataFrame(seq_dict)

    # Prepare storage for:
    # 1. Filtered positions for CSV output
    # 2. Binary difference indicator relative to Starosele for valid positions
    filtered_rows = []
    diff_data = {seq_id: [] for seq_id in df.columns if seq_id != starosele_key}
    included_positions = []

    # Loop over each position in the alignment
    for pos in df.index:
        starosele_base = df.loc[pos, starosele_key]

        # Skip if Starosele's base is ambiguous
        if starosele_base == "N":
            continue

        # Get the bases for all other sequences
        others = df.loc[pos].drop(starosele_key)

        # Apply filtering rules:
        # Exclude positions where Starosele has 'T' and all others have 'C'
        if starosele_base == "T" and (others == "C").all():
            continue

        # Exclude positions where Starosele has 'A' and all others have 'G'
        if starosele_base == "A" and (others == "G").all():
            continue

        # If we reach here, this position will be included
        pos_index = pos + 1  # Convert to 1-indexed position
        included_positions.append(pos_index)

        # Save this row's alleles for CSV output
        row_data = {"Position": pos_index}
        row_data.update(df.loc[pos].to_dict())
        filtered_rows.append(row_data)

        # Record binary difference for each other sequence
        for seq_id, base in others.items():
            if base == "N":
                diff_data[seq_id].append(0)  # Count 'N' as no difference
            elif base != starosele_base:
                diff_data[seq_id].append(1)  # Count as difference
            else:
                diff_data[seq_id].append(0)  # No difference

    # Create a DataFrame for the filtered positions and save to CSV
    filtered_df = pd.DataFrame(filtered_rows)
    filtered_csv_file = "filtered_positions.csv"
    filtered_df.to_csv(filtered_csv_file, index=False)
    print(f"Filtered alignment positions saved to '{filtered_csv_file}'.")

    # Compute total differences per sample
    diff_counts = {seq_id: sum(diff_data[seq_id]) for seq_id in diff_data}
    diff_df = pd.DataFrame(list(diff_counts.items()), columns=["Sample", "Differences"])
    diff_df = diff_df[diff_df["Sample"] != "Pan troglodytes"]  # Exclude chimpanzee
    diff_df = diff_df.sort_values(by="Differences", ascending=False)

    # Generate bar plot
    plt.figure(figsize=(12, 6))
    cmap = sns.color_palette("coolwarm", as_cmap=True)
    colors = diff_df["Differences"].rank(method="min", pct=True)

    sns.barplot(x="Sample", y="Differences", data=diff_df, palette=cmap(colors))
    plt.xticks(rotation=90)
    plt.xlabel("Individuals")
    plt.ylabel("Number of nucleotide differences")
    plt.title("Filtered Mitochondrial Dissimilarity")

    plt.text(-0.5, max(diff_df["Differences"]) * 0.98, "Most dissimilar", fontsize=12, fontweight='bold', color='red')
    plt.text(len(diff_df) - 0.5, max(diff_df["Differences"]) * 0.98, "Most similar", fontsize=12, fontweight='bold', color='blue', ha='right')

    plt.tight_layout()
    barplot_file = "starosele_filtered_differences_barplot.png"
    plt.savefig(barplot_file, dpi=300)
    plt.close()
    print(f"Filtered bar plot saved to '{barplot_file}'.")

if __name__ == "__main__":
    main()
