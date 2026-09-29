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
git clone [https://github.com/](https://github.com/)<your-username>/protein-conformation-pca-pipeline.git
cd protein-conformation-pca-pipeline
pip install -r requirements.txt
