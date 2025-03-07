import sys
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main():
    if len(sys.argv) < 2:
        print("Usage: python filter_alignment.py alignment.fasta")
        sys.exit(1)

    alignment_file = sys.argv[1]

    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except Exception as e:
        print(f"Error reading alignment: {e}")
        sys.exit(1)

    seq_dict = {record.description: list(str(record.seq).upper()) for record in alignment}

    if not any("Starosele" in key for key in seq_dict.keys()):
        print("Error: 'Starosele' sequence not found in the alignment.")
        sys.exit(1)

    df = pd.DataFrame(seq_dict)
    starosele_key = next(key for key in df.columns if "Starosele" in key)

    filtered_positions = []
    for pos in df.index:
        starosele_base = df.loc[pos, starosele_key]
        if starosele_base == "N":
            continue
        others = df.loc[pos].drop(starosele_key)
        if starosele_base == "T" and (others == "C").all():
            continue
        if starosele_base == "A" and (others == "G").all():
            continue
        filtered_positions.append(pos)

    filtered_df = df.loc[filtered_positions]
    filtered_seq_dict = {col: "".join(filtered_df[col]) for col in filtered_df.columns}

    output_file = "filtered_for_tree.fasta"
    with open(output_file, "w") as f:
        SeqIO.write([SeqRecord(Seq(seq), id=key, description="") for key, seq in filtered_seq_dict.items()], f, "fasta")
    
    print(f"Filtered alignment saved to '{output_file}' ({len(filtered_positions)} positions kept).")

if __name__ == "__main__":
    main()

