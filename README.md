# Protein Conformational Dynamics: Perturbation, PCA & Clustering Pipeline

A modular Python framework designed to simulate structural variations resembling conformational heterogeneity observed in structural biology (e.g., cryo-EM ensembles), perform dimensionality reduction via Principal Component Analysis (PCA), partition conformational landscapes with K-Means clustering, and enable bidirectional coordinate reconstruction and unseen projection.

---

## Key Features

- **Conformational Ensemble Perturbation:** Introduces controlled atomic displacement (0–2 Å) on a starting PDB coordinate set.
- **PCA Dimensionality Reduction:** Reduces $3N$-dimensional coordinate trajectories to primary orthogonal modes of variance.
- **K-Means Clustering:** Discretizes the 2D conformational subspace into distinct structural clusters.
- **Inverse Coordinate Mapping:** Reconstructs full 3D Cartesian coordinates ($x, y, z$) from arbitrary PC space coordinates using inverse PCA transforms.
- **Novel Conformation Projection:** Projects unseen PDB structures directly onto the established reference eigen-space.

---

## Installation

Clone the repository and set up your Python environment:

```bash
git clone [https://github.com/](https://github.com/)<your-username>/protein-dynamics-pca-pipeline.git
cd protein-conformation-pca-pipeline
pip install -r requirements.txt
```


## Workflow & Execution

Run the analysis scripts sequentially in the following hierarchy:

| **Step** | **Script** | **Input** | **Output** | **Purpose** |
|---|---|---|---|---|
| **01** | `code_1_structure_building.py` | `data/1L2Y.pdb` | `new_structures/*.pdb` | Generates coordinate perturbation ensembles. |
| **02** | `code_2_PCA_analysis.py` | `new_structures/*.pdb` | `PCA_projection.csv`, `pca_model.pkl`, `log.txt` | Flattens coordinates and extracts principal components. |
| **03** | `code_3_clustering.py` | `PCA_projection.csv` | `PCA_clusters.csv`, `PCA_plot.png` | Partitions conformational space (`k=4`) and visualizes the conformational landscape. |
| **04** | `code_4_reverse_construction.py` | `pca_model.pkl`, template PDB | `predicted_structure.pdb` | Reconstructs 3D coordinates from custom PC1/PC2 values. |
| **05** | `code_5_projection.py` | `pca_model.pkl`, novel PDB | `new_structure_added_plot.png` | Projects a novel/unseen structure onto the established PCA conformational space. |
