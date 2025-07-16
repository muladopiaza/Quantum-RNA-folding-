import RNA

# Load sequence from FASTA file
def load_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        seq = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return seq

# Save dot-bracket and MFE to file
def save_dbn(seq, structure, mfe, file_path):
    with open(file_path, 'w') as f:
        f.write(f">{file_path}\n{seq}\n{structure} ({mfe:.2f} kcal/mol)\n")

# Save connectivity table (CT format)
def save_ct(seq, structure, file_path):
    pairs = RNA.ptable(structure)
    with open(file_path, 'w') as f:
        f.write(f"{len(seq)}\tENERGY = 0\n")
        for i in range(1, len(seq)+1):
            prev = i - 1 if i > 1 else 0
            next = i + 1 if i < len(seq) else 0
            pair = pairs[i]
            f.write(f"{i} {seq[i-1]} {prev} {next} {pair} {i}\n")

# Run ViennaRNA folding
def fold_rna(seq):
    fc = RNA.fold_compound(seq)
    structure, mfe = fc.mfe()
    return structure, mfe

# Extract base pairs (i, j)
def get_base_pairs(structure):
    return [(i, j) for i, j in enumerate(RNA.ptable(structure)) if i < j]
def print_rnapdbee_dotbracket(filepath):
    """Print sequence and dot-bracket notation from RNApdbee output file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
        if len(lines) >= 2:
            seq = lines[1].strip()
            struct = lines[2].strip()
            print(f"Length:{len(seq)}")
            print("\nExperimental Structure")
            print("Sequence : " + seq)
            print("Structure: " + struct)
        else:
            print("File format issue: expected 2 lines (sequence + structure).")

# === Main ===
fasta_path = "rcsb_pdb_1EHZ.fasta"  # Change if needed
seq = load_fasta(fasta_path)
structure, mfe = fold_rna(seq)
print_rnapdbee_dotbracket("1ehz-2D-dotbracket.txt")

base_pairs = get_base_pairs(structure)
print(f"ViennaRNA Prediction:")
print(f"Sequence: {seq}")
print(f"Structure: {structure}")
print(f"MFE: {mfe:.2f} kcal/mol")
print(f"Base pairs (i, j): {base_pairs}")

# Save outputs
save_dbn(seq, structure, mfe, "output.dbn")
save_ct(seq, structure, "output.ct")

def parse_ct(filepath):
    """Extract base pairs from CT file as (i, j) tuples (1-based)."""
    pairs = set()
    with open(filepath, 'r') as f:
        lines = f.readlines()[1:]  # Skip header
        for line in lines:
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            i = int(parts[0])
            j = int(parts[4])
            if j > 0 and i < j:  # Avoid duplicates (only one direction)
                pairs.add((i, j))
    return pairs

def compare_ct_files(predicted_ct, experimental_ct):
    pred_pairs = parse_ct(predicted_ct)
    exp_pairs = parse_ct(experimental_ct)

    true_positives = pred_pairs & exp_pairs
    false_positives = pred_pairs - exp_pairs
    false_negatives = exp_pairs - pred_pairs

    tp = len(true_positives)
    fp = len(false_positives)
    fn = len(false_negatives)

    precision = tp / (tp + fp) if (tp + fp) else 0
    recall = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) else 0

    print(f"\nComparison Metrics:")
    print(f"True Positives (Correct pairs): {tp}")
    print(f"False Positives (Incorrect extra pairs): {fp}")
    print(f"False Negatives (Missed real pairs): {fn}")
    print(f"Precision: {precision:.3f}")
    print(f"Recall: {recall:.3f}")
    print(f"F1 Score: {f1:.3f}")

    if false_positives:
        print("\n False Positives (Predicted but not in experimental):")
        for pair in sorted(false_positives):
            print(f"  Predicted extra: {pair}")

    if false_negatives:
        print("\n False Negatives (Missing from predicted structure):")
        for pair in sorted(false_negatives):
            print(f"  Missed true pair: {pair}")

# === Run this ===
predicted_ct = "output.ct"
experimental_ct = "1ehz-2D-ct.txt"

compare_ct_files(predicted_ct, experimental_ct)
