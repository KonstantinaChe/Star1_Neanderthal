from Bio import SeqIO

# Input and output file paths
input_fasta = "MAFFT_input.fasta"
output_fasta = "filtered_no_Ns_thr1.fasta"
positions_removed_file = "removed_positions_thr1.txt"
starosele_name = "Starosele"

# Read sequences from the FASTA file
sequences = list(SeqIO.parse(input_fasta, "fasta"))

# Find the Starosele sequence
starosele_seq = None
for record in sequences:
    if record.id == starosele_name:
        starosele_seq = record.seq
        break

if not starosele_seq:
    raise ValueError(f"Sequence named '{starosele_name}' not found in the input file.")

# Identify positions to keep and positions to remove
positions_to_keep = []
positions_removed = []
for i, base in enumerate(starosele_seq):
    if base != "N":
        positions_to_keep.append(i)
    else:
        positions_removed.append(i + 1)  # Use 1-based indexing for human readability

# Filter all sequences to retain only positions where Starosele is not 'N'
filtered_sequences = []
for record in sequences:
    filtered_seq = "".join(record.seq[i] for i in positions_to_keep)
    record.seq = filtered_seq
    filtered_sequences.append(record)

# Write the filtered sequences to a new FASTA file
SeqIO.write(filtered_sequences, output_fasta, "fasta")

# Write the removed positions to a separate file
with open(positions_removed_file, "w") as f:
    f.write("Positions removed due to 'N' in Starosele:\n")
    f.write("\n".join(map(str, positions_removed)))

print(f"Filtered sequences saved to {output_fasta}.")
print(f"Positions removed saved to {positions_removed_file}.")
