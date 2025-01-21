from Bio import SeqIO
from collections import Counter

# Input and output file paths
input_fasta = "filtered_no_Ns_thr1.fasta"
output_fasta = "filtered_no_private_mutations_and_Ns_thr1.fasta"
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

# Identify positions where Starosele has non-private mutations
non_private_positions = []
for i, star_base in enumerate(starosele_seq):
    # Gather bases from all sequences at position i
    column_bases = [record.seq[i] for record in sequences]
    base_counts = Counter(column_bases)
    
    # Check if Starosele's base is private (appears only once)
    if base_counts[star_base] > 1:
        non_private_positions.append(i)

# Filter all sequences to retain only non-private positions
filtered_sequences = []
for record in sequences:
    filtered_seq = "".join(record.seq[i] for i in non_private_positions)
    record.seq = filtered_seq
    filtered_sequences.append(record)

# Write the filtered sequences to a new FASTA file
SeqIO.write(filtered_sequences, output_fasta, "fasta")

print(f"Filtered sequences saved to {output_fasta}.")
