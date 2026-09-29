import os
import joblib
import numpy as np
from Bio.PDB import PDBParser, PDBIO

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
}

PCA_MODEL_FILE = "pca_model.pkl"
OUTPUT_FILE = "predicted_structure.pdb"


def main():
    print("--- PCA 3D Structure Reconstructor ---\n")
    
    try:
        pc1_val = float(input("Enter the target coordinate for PC1: "))
        pc2_val = float(input("Enter the target coordinate for PC2: "))
    except ValueError:
        print("\nError: You must enter valid numerical values.")
        return

    template_file = input("Enter the name of the template PDB file (e.g., structure_001.pdb): ").strip()

    if not os.path.isfile(template_file):
        print(f"\nError: Could not find '{template_file}' in the current working directory.")
        return

    if not os.path.isfile(PCA_MODEL_FILE):
        print(f"\nError: Could not find '{PCA_MODEL_FILE}'. Make sure you saved it from the previous script.")
        return

#  RECONSTRUCT COORDINATES (INVERSE TRANSFORM)----------------------------------------

    print("\nLoading PCA model...")
    pca = joblib.load(PCA_MODEL_FILE)

    print(f"Calculating 3D coordinates for PC1={pc1_val}, PC2={pc2_val}...")
    
    reconstructed_flat = pca.inverse_transform([[pc1_val, pc2_val]])[0]
    
    reconstructed_coords = reconstructed_flat.reshape(-1, 3)

    print(f"Reading template framework from '{template_file}'...")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("template", template_file)

    atom_index = 0
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in AMINO_ACIDS:
                    for atom in residue:
                        if atom_index >= len(reconstructed_coords):
                            print("\nError: The template PDB has more atoms than the PCA model expects.")
                            return
                        
                        atom.set_coord(reconstructed_coords[atom_index])
                        atom_index += 1

    if atom_index != len(reconstructed_coords):
        print(f"\nWarning: Atom count mismatch. The PCA model generated {len(reconstructed_coords)} coordinates, "
              f"but only {atom_index} were applied to the template.")

    print(f"Saving rebuilt protein...")
    io = PDBIO()
    io.set_structure(structure)
    io.save(OUTPUT_FILE)
    
    print(f"\nSuccess! The predicted structure has been saved as '{OUTPUT_FILE}'.")


if __name__ == "__main__":
    main()