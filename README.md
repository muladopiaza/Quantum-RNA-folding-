# RNA Secondary Structure Prediction using QUBO and Simulated Annealing

This project explores RNA secondary structure prediction through the lens of **combinatorial optimization** by formulating the folding problem as a **Quadratic Unconstrained Binary Optimization (QUBO)** model and solving it using **Simulated Annealing**.

##  Objective

Predict the secondary structure of an RNA sequence (tRNA^Phe from *Saccharomyces cerevisiae*, 1EHZ) using:
- Biologically-informed QUBO formulation
- Simulated Annealing (via `neal`)
- Comparison with experimentally validated structure and classical methods (ViennaRNA, IPknot, IPknot++)

---

##  Methodology

### 1. QUBO Construction
- All **valid base pairs** (Watson-Crick and wobble) are extracted with a minimum loop length of 3.
- Each valid base pair is assigned a **score based on its stability**:
  - GC: -2.0
  - AU: -1.5
  - GU: -1.0
- Top 15% strongest-scoring pairs are retained.
- A QUBO matrix is built with:
  - **Diagonal terms** representing pair stability
  - **Off-diagonal terms** penalizing overlapping/conflicting pairs (+5.0)
  - **Bonus** for stacking interactions (-1.5)

### 2. Optimization
- Simulated Annealing is used to minimize the QUBO.
- The resulting configuration represents selected base pairs.
- These are converted into dot-bracket notation for structural comparison.

### 3. Evaluation
Predicted base pairs are evaluated against:
-  Experimentally derived structure (from PDB 1EHZ)
-  Predictions by ViennaRNA, IPknot, IPknot++

Metrics:
- **True Positives (TP)**
- **False Positives (FP)**
- **False Negatives (FN)**
- **Precision**, **Recall**, **F1 Score**

---

##  Results

| Method             | TP  | FP  | FN  | Precision | Recall | F1 Score |
|--------------------|-----|-----|-----|-----------|--------|----------|
| ViennaRNA (RNAfold)| 21  | 0   | 9   | 1.000     | 0.700  | 0.824    |
| IPknot             | 11  | 9   | 9   | 0.550     | 0.550  | 0.550    |
| IPknot++ (NUPACK)  | 0   | 8   | 20  | 0.000     | 0.000  | 0.000    |
| **QUBO + Annealing** | 0   | 9   | 30  | 0.000     | 0.000  | 0.000    |

> Despite biologically motivated constraints, the QUBO model failed to predict any true positives under simulated annealing, indicating limitations in the modeling assumptions.

---

##  Insights

- Even though constraints and biological priors were carefully encoded, the QUBO formulation may not sufficiently capture the **free energy landscape** or **pseudoknot complexity** of RNA folding.
- The classical methods (especially thermodynamic ones) still outperform combinatorial models like QUBO in realistic settings.


