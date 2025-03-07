import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio import AlignIO

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_alignment.py alignment.fasta")
        sys.exit(1)

    alignment_file = sys.argv[1]

    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except Exception as e:
        print(f"Error reading alignment: {e}")
        sys.exit(1)

    seq_dict = {}
    for record in alignment:
        seq_dict[record.description] = list(str(record.seq).upper())

    found_starosele = any("Starosele" in key for key in seq_dict.keys())
    if not found_starosele:
        print("Error: 'Starosele' sequence not found in the alignment.")
        sys.exit(1)

    df = pd.DataFrame(seq_dict)

    starosele_keys = [key for key in df.columns if "Starosele" in key]
    if not starosele_keys:
        print("Error: 'Starosele' sequence not found in the alignment columns.")
        sys.exit(1)
    starosele_key = starosele_keys[0]

    filtered_rows = []
    included_positions = []

    for pos in df.index:
        starosele_base = df.loc[pos, starosele_key]
        if starosele_base == "N":
            continue
        others = df.loc[pos].drop(starosele_key)
        if starosele_base == "T" and (others == "C").all():
            continue
        if starosele_base == "A" and (others == "G").all():
            continue
        pos_index = pos + 1
        included_positions.append(pos_index)
        row_data = {"Position": pos_index}
        row_data.update(df.loc[pos].to_dict())
        filtered_rows.append(row_data)

    filtered_df = pd.DataFrame(filtered_rows)
    filtered_df.to_csv("hopefully_last_filtered_alignment.csv", index=False)
    print(f"Filtered alignment saved to 'filtered_alignment.csv' ({len(filtered_df)} positions kept).")

    alignment_filtered = filtered_df.drop(columns=["Position"])
    seq_names = [name for name in alignment_filtered.columns if "pan" not in name.lower()]
    pairwise_diff = pd.DataFrame(index=seq_names, columns=seq_names)

    for seq1 in seq_names:
        for seq2 in seq_names:
            valid_positions = (alignment_filtered[seq1] != "N") & (alignment_filtered[seq2] != "N")
            differences = ((alignment_filtered[seq1] != alignment_filtered[seq2]) & valid_positions).sum()
            pairwise_diff.loc[seq1, seq2] = differences

    pairwise_diff = pairwise_diff.astype(int)
    mask = np.triu(np.ones_like(pairwise_diff, dtype=bool))

    plt.figure(figsize=(max(10, len(seq_names)), max(8, len(seq_names))))
    ax2 = sns.heatmap(pairwise_diff, annot=True, fmt="d", cmap="magma",
                      cbar_kws={"label": "Number of Differences"}, mask=mask)
    plt.xlabel("Sequence")
    plt.ylabel("Sequence")
    plt.title("Pairwise NUcleotide Differences")
    plt.tight_layout()
    pairwise_heatmap_file = "hopefully_the_last_heatmap.svg"
    plt.savefig(pairwise_heatmap_file, dpi=300)
    plt.close()
    print(f"Pairwise differences heatmap saved to '{pairwise_heatmap_file}'.")

if __name__ == "__main__":
    main()
