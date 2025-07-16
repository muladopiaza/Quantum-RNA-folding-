# qubo_build.py
import numpy as np
from itertools import combinations
import pickle

# Input RNA sequence
rna_sequence = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"

VALID_PAIRS = {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G')}
PAIR_SCORES = {
    ('G', 'C'): -2.0,
    ('C', 'G'): -2.0,
    ('A', 'U'): -1.5,
    ('U', 'A'): -1.5,
    ('G', 'U'): -1.0,
    ('U', 'G'): -1.0
}

def get_valid_pairs(seq, min_loop_len=3):
    n = len(seq)
    valid = []
    for i in range(n):
        for j in range(i + min_loop_len + 1, n):
            if (seq[i], seq[j]) in VALID_PAIRS:
                valid.append((i + 1, j + 1))  # 1-based indexing
    return valid

def score_pair(i, j, seq):
    base_pair = (seq[i - 1], seq[j - 1])
    return PAIR_SCORES.get(base_pair, -0.5)

def build_qubo(pairs, seq):
    n = len(pairs)
    Q = np.zeros((n, n))

    # Diagonal terms (pair energies)
    for idx, (i, j) in enumerate(pairs):
        Q[idx, idx] = score_pair(i, j, seq)

    # Off-diagonal terms: conflict + stacking
    for a, b in combinations(range(n), 2):
        i1, j1 = pairs[a]
        i2, j2 = pairs[b]

        if set((i1, j1)) & set((i2, j2)):
            Q[a, b] = Q[b, a] = 5.0  # hard conflict penalty
        elif (i1 + 1 == i2 and j1 - 1 == j2):
            Q[a, b] = Q[b, a] = -1.5  # stacking bonus

    return Q

def filter_top_percent(pairs, seq, percent=0.15):
    scored = [(pair, score_pair(pair[0], pair[1], seq)) for pair in pairs]
    scored.sort(key=lambda x: x[1])  # sort most negative (strong) first
    cutoff = int(len(scored) * percent)
    return [pair for pair, _ in scored[:cutoff]]

# ========================
# QUBO Construction Steps
# ========================

# Step 1: Find all biologically valid pairs
all_valid_pairs = get_valid_pairs(rna_sequence)
print("Initial valid base pairs (before scoring):", len(all_valid_pairs))

# Step 2: Keep only the strongest 15%
filtered_pairs = filter_top_percent(all_valid_pairs, rna_sequence, percent=0.15)
print("Filtered top 15% pairs:", len(filtered_pairs))

# Step 3: Build QUBO with stacking and hard conflicts
Q = build_qubo(filtered_pairs, rna_sequence)
print("QUBO shape:", Q.shape)

# Step 4: Print Diagonal and Example Off-Diagonal
diagonal_biases = [round(Q[i, i], 2) for i in range(min(10, len(Q)))]
print("Diagonal (biases):", diagonal_biases)

nonzero = [(i, j, round(Q[i, j], 2)) for i in range(len(Q)) for j in range(i + 1, len(Q)) if Q[i, j] != 0]
print("Number of QUBO terms (non-zero off-diagonal):", len(nonzero))
print("Example QUBO terms:", nonzero[:10])

# Step 5: Save to file
with open("qubo_data.pkl", "wb") as f:
    pickle.dump((Q, filtered_pairs), f)

print("Saved biologically-informed QUBO to qubo_data.pkl")
print("Total base pairs saved:", len(filtered_pairs))
