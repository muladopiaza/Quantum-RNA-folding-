# run_annealing.py
import pickle
from neal import SimulatedAnnealingSampler

# Load QUBO matrix and valid pairs
with open("qubo_data.pkl", "rb") as f:
    Q, valid_pairs = pickle.load(f)

# Convert Q to dict
qubo_dict = {(i, j): Q[i, j] for i in range(len(Q)) for j in range(i, len(Q)) if Q[i, j] != 0}

print("Diagonal (biases):", Q.diagonal()[:10])
print("Example QUBO terms:", [(i, j, Q[i,j]) for i in range(10) for j in range(i+1, 10) if Q[i,j] != 0])

# Run annealer
sampler = SimulatedAnnealingSampler()
response = sampler.sample_qubo(qubo_dict, num_reads=5000)

# Extract best result
best_sample = list(response.first.sample.values())
selected_pairs = [valid_pairs[i] for i, val in enumerate(best_sample) if val == 1]
print("Selected base pairs:", selected_pairs)

# -------------------------------
#  DOT-BRACKET CONVERTER
# -------------------------------

def base_pairs_to_dot_bracket(seq_len, base_pairs):
    structure = ['.'] * seq_len
    for i, j in base_pairs:
        i, j = i - 1, j - 1  # convert to 0-based
        if i < j:
            structure[i] = '('
            structure[j] = ')'
        else:
            structure[j] = '('
            structure[i] = ')'
    return ''.join(structure)

# Safely get sequence length
sequence_length = max(max(pair) for pair in selected_pairs)
dot_bracket = base_pairs_to_dot_bracket(sequence_length, selected_pairs)
print("\n QUBO Predicted Dot-Bracket Structure:")
print("Annealer base pairs:", selected_pairs)

print(dot_bracket)

# -------------------------------
#  Evaluation
# -------------------------------

# 0-based ViennaRNA base pairs
vienna_pairs = [
    (0, 76), (1, 72), (2, 71), (3, 70), (4, 69), (5, 68), (6, 67), (7, 66),
    (10, 25), (11, 24), (12, 23), (13, 22),
    (27, 43), (28, 42), (29, 41), (30, 40), (31, 39),
    (49, 65), (50, 64), (51, 63), (52, 62), (53, 61)
]

# 0-based Experimental base pairs (additional inferred)
experimental_pairs = vienna_pairs + [
    (8, 14), (15, 48), (16, 59), (18, 55), (19, 56),
    (26, 44), (32, 38), (33, 36), (54, 58)
]

# Convert selected pairs to 0-based
selected_pairs_0 = [(i - 1, j - 1) for i, j in selected_pairs]

# -------------------------------
#  Comparison Function
# -------------------------------

def evaluate_prediction(predicted, ground_truth):
    pred_set = set(tuple(sorted(p)) for p in predicted)
    true_set = set(tuple(sorted(p)) for p in ground_truth)

    tp = pred_set & true_set
    fp = pred_set - true_set
    fn = true_set - pred_set

    precision = len(tp) / (len(tp) + len(fp)) if tp or fp else 0.0
    recall = len(tp) / (len(tp) + len(fn)) if tp or fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall > 0 else 0.0

    return tp, fp, fn, precision, recall, f1

# -------------------------------
#  Evaluation Results
# -------------------------------

tp_v, fp_v, fn_v, prec_v, rec_v, f1_v = evaluate_prediction(selected_pairs_0, vienna_pairs)
print("\n Comparison with ViennaRNA Structure:")
print(f"   True Positives: {len(tp_v)}")
print(f"   False Positives: {len(fp_v)}")
print(f"   False Negatives: {len(fn_v)}")
print(f"   Precision: {prec_v:.3f}")
print(f"   Recall: {rec_v:.3f}")
print(f"   F1 Score: {f1_v:.3f}")

tp_e, fp_e, fn_e, prec_e, rec_e, f1_e = evaluate_prediction(selected_pairs_0, experimental_pairs)
print("\n Comparison with Experimental Structure:")
print(f"   True Positives: {len(tp_e)}")
print(f"   False Positives: {len(fp_e)}")
print(f"   False Negatives: {len(fn_e)}")
print(f"   Precision: {prec_e:.3f}")
print(f"   Recall: {rec_e:.3f}")
print(f"   F1 Score: {f1_e:.3f}")
