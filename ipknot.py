# print_experimental_from_file.py

def load_dot_bracket(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
        # Remove possible FASTA-like header
        lines = [line.strip() for line in lines if not line.startswith(">") and line.strip() != ""]
    if len(lines) < 2:
        raise ValueError("File does not contain both sequence and structure.")
    sequence = lines[0]
    structure = lines[1]
    return sequence, structure
def dotbracket_to_pairs(structure):
    """Convert dot-bracket notation to base pair list."""
    stack = []
    pairs = []
    pair_symbols = {'(': ')', '[': ']', '{': '}', '<': '>'}
    match = {v: k for k, v in pair_symbols.items()}

    for i, char in enumerate(structure):
        if char in pair_symbols:
            stack.append((char, i))
        elif char in match:
            if stack:
                opener, j = stack.pop()
                if pair_symbols[opener] == char:
                    pairs.append((j, i))
    return set(pairs)

def compare_structures(exp_structure, pred_structure):
    """Compare experimental and predicted dot-bracket structures and print metrics."""
    exp_pairs = dotbracket_to_pairs(exp_structure)
    pred_pairs = dotbracket_to_pairs(pred_structure)

    tp = exp_pairs & pred_pairs
    fp = pred_pairs - exp_pairs
    fn = exp_pairs - pred_pairs

    precision = len(tp) / (len(tp) + len(fp)) if (tp or fp) else 0.0
    recall = len(tp) / (len(tp) + len(fn)) if (tp or fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    print("\n IPknot vs Experimental Metrics")
    print(f"True Positives: {len(tp)}")
    print(f"False Positives: {len(fp)}")
    print(f"False Negatives: {len(fn)}")
    print(f"Precision: {precision:.3f}")
    print(f"Recall:    {recall:.3f}")
    print(f"F1 Score:  {f1:.3f}")

    if fn:
        print("\nFalse Negatives (missed base pairs):")
        for pair in sorted(fn):
            print(f"  Missed true pair: {pair}")
    if fp:
        print("\nFalse Positives (wrong extra pairs):")
        for pair in sorted(fp):
            print(f"  Incorrect predicted pair: {pair}")


# Replace with your actual file path
file_path_exp = "1ehz-2D-dotbracket.txt"
file_path = "ipknot.vienna"

# Load and print
sequence_exp, structure_exp = load_dot_bracket(file_path_exp)
print("Experimental Structure")
print("Sequence :", sequence_exp)
print("Structure:", structure_exp)
sequence , structure = load_dot_bracket(file_path)
print("Ipknot Result:")
print("Sequence :", sequence)
print("Structure:", structure)
compare_structures(structure_exp, structure)
# NUPACK-based IPknot result
ipknot_seq = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
ipknot_struct = "(({{{{{{...{{{{{((((.[[[[[))..}}}}}}}}}}}...{{{.{{{]]]]]....))}}}.....))}}}."
print("Experimental Structure")
print("Sequence :", sequence_exp)
print("Structure:", structure_exp)
print("\nNUPACK Model IPknot++ Result:")
print("Sequence :", ipknot_seq)
print("Structure:", ipknot_struct)

compare_structures(structure_exp, ipknot_struct)

