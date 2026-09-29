import os
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

PCA_MODEL_FILE = "pca_model.pkl"
CLUSTER_DATA_FILE = "PCA_clusters.csv"
OUTPUT_PLOT = "new_structure_added_plot.png"

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
}

def main():
    print("--- PCA Unseen Structure Projector ---\n")

    new_pdb_file = input("Please enter the name of the new .pdb structure (e.g., structure_X.pdb): ").strip()

    if not os.path.isfile(new_pdb_file):
        print(f"\nError: Could not find '{new_pdb_file}' in the current working directory.")
        return

    for required_file in [PCA_MODEL_FILE, CLUSTER_DATA_FILE]:
        if not os.path.isfile(required_file):
            print(f"\nError: Missing required file '{required_file}'. Ensure it is in the same folder.")
            return


#  LOAD MODEL AND CLUSTER DATA

    print(f"\nLoading PCA model from {PCA_MODEL_FILE}...")
    pca = joblib.load(PCA_MODEL_FILE)
    
    print(f"Loading cluster data from {CLUSTER_DATA_FILE}...")
    df = pd.read_csv(CLUSTER_DATA_FILE)


#  EXTRACT COORDINATES FROM THE NEW PDB-----------------------------------------------------

    print(f"\nExtracting coordinates from '{new_pdb_file}'...")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("unseen_structure", new_pdb_file)

    coords_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in AMINO_ACIDS:
                    for atom in residue:
                        coords_list.append(atom.get_coord())

    if len(coords_list) == 0:
        print("Error: No standard amino-acid atoms found in this PDB.")
        return

    unseen_flat = np.array(coords_list).flatten()

    expected_features = pca.mean_.shape[0]
    if len(unseen_flat) != expected_features:
        print(f"\nError: The PCA model expects {expected_features} coordinates, "
              f"but your new structure contains {len(unseen_flat)}.")
        return

#  PROJECT INTO PCA SPACE---------------------------------------------------

    print("Projecting structure into the 2D PCA space...")
    
    unseen_pc = pca.transform(unseen_flat.reshape(1, -1))[0]
    
    print(f" -> Projected Position: PC1 = {unseen_pc[0]:.4f}, PC2 = {unseen_pc[1]:.4f}")


#  GENERATE AND SAVE THE PLOT--------------------------------------

    print(f"\nGenerating cluster plot...")
    
    plt.figure(figsize=(10, 8))
    colors = ['red', 'blue', 'yellow', 'green']

    for cluster_id in range(4):
        cluster_data = df[df['Cluster'] == cluster_id]
        plt.scatter(
            cluster_data['PC1'], 
            cluster_data['PC2'], 
            c=colors[cluster_id], 
            label=f'Class {cluster_id}', 
            edgecolor='black', 
            s=60, 
            alpha=0.6
        )

    plt.scatter(
        unseen_pc[0], 
        unseen_pc[1], 
        c='black', 
        marker='^',   
        s=80,        
        label='New Structure'
    )

    plt.title('PCA Space: K-Means Clusters with New Structure', fontsize=14)
    plt.xlabel('Principal Component 1', fontsize=12)
    plt.ylabel('Principal Component 2', fontsize=12)
    
    plt.legend(title="Legend", loc="best")
    plt.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=300)
    print(f"Plot saved successfully as '{OUTPUT_PLOT}'.")
    
    plt.show()

if __name__ == "__main__":
    main()