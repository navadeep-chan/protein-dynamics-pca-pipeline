import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from sklearn.decomposition import PCA
import joblib


INPUT_FOLDER = "new_structures"
LOG_FILE = "log.txt"
PROJECTION_FILE = "PCA_projection.csv"

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
}

pdb_files = sorted(
    [
        os.path.join(INPUT_FOLDER, file)
        for file in os.listdir(INPUT_FOLDER)
        if file.lower().endswith(".pdb")
    ]
)

if len(pdb_files) == 0:
    raise FileNotFoundError(
        f"No PDB files were found inside '{INPUT_FOLDER}'."
    )

print(f"Number of PDB structures found: {len(pdb_files)}")


all_coordinates = []
structure_names = []

reference_atom_count = None

parser = PDBParser(QUIET=True)

for pdb_file in pdb_files:

    structure_id = os.path.basename(pdb_file)
    structure = parser.get_structure(structure_id, pdb_file)

    coords_list = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in AMINO_ACIDS:
                    for atom in residue:
                        coords_list.append(atom.get_coord())

    if len(coords_list) == 0:
        raise ValueError(
            f"No protein atoms found in {pdb_file}"
        )
    coordinates = np.array(coords_list)

    atom_count = len(coordinates)

    if reference_atom_count is None:
        reference_atom_count = atom_count

    elif atom_count != reference_atom_count:
        raise ValueError(
            f"Atom count mismatch in {pdb_file}. "
            f"Expected {reference_atom_count}, "
            f"found {atom_count}."
        )

    flattened_coordinates = coordinates.flatten()

    all_coordinates.append(flattened_coordinates)

    structure_names.append(
        os.path.basename(pdb_file)
    )

#  BUILDING COORDINATE MATRIX---------------------------------------------


coordinate_matrix = np.array(all_coordinates)

number_of_structures = coordinate_matrix.shape[0] 
"""ROWS"""
number_of_coordinates = coordinate_matrix.shape[1] 
"""COLUMNS"""
print(f"Number of atoms per structure: {reference_atom_count}")
print(
    f"Coordinate matrix shape: "
    f"{coordinate_matrix.shape}"
)

print(
    f"Each structure contains "
    f"{number_of_coordinates} coordinate values."
)

#  PCA---------------------------------------------

print("\nPerforming PCA...")

pca = PCA(n_components=2, random_state=20)

pca_projection = pca.fit_transform(
    coordinate_matrix
)

eigenvalues = pca.explained_variance_

eigenvectors = pca.components_

explained_variance_ratio = (
    pca.explained_variance_ratio_
)

with open(LOG_FILE, "w") as log:

    log.write("PCA ANALYSIS LOG\n")
    log.write("================\n\n")

    log.write(
        f"Number of structures: "
        f"{number_of_structures}\n"
    )

    log.write(
        f"Number of atoms per structure: "
        f"{reference_atom_count}\n"
    )

    log.write(
        f"Number of coordinate variables: "
        f"{number_of_coordinates}\n"
    )

    log.write(
        f"Coordinate matrix shape: "
        f"{coordinate_matrix.shape}\n\n"
    )

    log.write("Mean coordinate vector:\n")
    log.write(
        np.array2string(
            pca.mean_,
            precision=6,
            separator=", "
        )
    )
    log.write("\n\n")

# Eigenvalues
    log.write("Eigenvalues:\n")
    log.write(
        f"PC1: {eigenvalues[0]:.6f}\n"
    )
    log.write(
        f"PC2: {eigenvalues[1]:.6f}\n\n"
    )

# Explained variance
    log.write("Explained variance ratio:\n")
    log.write(
        f"PC1: {explained_variance_ratio[0]:.6f} "
        f"({explained_variance_ratio[0] * 100:.2f}%)\n"
    )

    log.write(
        f"PC2: {explained_variance_ratio[1]:.6f} "
        f"({explained_variance_ratio[1] * 100:.2f}%)\n\n"
    )

    log.write(
        f"Total variance explained by PC1 + PC2: "
        f"{sum(explained_variance_ratio) * 100:.2f}%\n\n"
    )

# Eigenvectors
    log.write(
        "Eigenvectors / Principal Component directions:\n"
    )

    log.write(
        "\nPC1 eigenvector:\n"
    )

    log.write(
        np.array2string(
            eigenvectors[0],
            precision=6,
            separator=", "
        )
    )

    log.write("\n\nPC2 eigenvector:\n")

    log.write(
        np.array2string(
            eigenvectors[1],
            precision=6,
            separator=", "
        )
    )

    log.write("\n")


#  SAVING PCA PROJECTION-----------------------------------------

projection_df = pd.DataFrame(
    {
        "Structure": structure_names,
        "PC1": pca_projection[:, 0],
        "PC2": pca_projection[:, 1]
    }
)

projection_df.to_csv(
    PROJECTION_FILE,
    index=False
)

print("\nPCA analysis completed.")

print(f"PCA information saved to: {LOG_FILE}")

print(
    f"PCA projections saved to: "
    f"{PROJECTION_FILE}"
)

print("\nFirst five PCA projections:")
print(projection_df.head())

joblib.dump(pca, 'pca_model.pkl')
print("\nPCA model saved to 'pca_model.pkl' for later use.")