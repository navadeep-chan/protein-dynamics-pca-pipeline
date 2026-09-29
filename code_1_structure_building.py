import os
import sys
import random

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
}

def read_pdb(filename):
    atom_records = []

    with open(filename, "r") as file:
        for line in file:

            if not line.startswith("ATOM"):
                continue
            residue_name = line[17:20].strip()

            if residue_name not in AMINO_ACIDS:
                continue

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            atom_records.append({
                "line": line.rstrip("\n"),
                "x": x,
                "y": y,
                "z": z
            })

    return atom_records


def generate_structure(atom_records):

    new_lines = []

    for atom in atom_records:

        dx = random.uniform(0, 2)
        dy = random.uniform(0, 2)
        dz = random.uniform(0, 2)

        new_x = atom["x"] + dx
        new_y = atom["y"] + dy
        new_z = atom["z"] + dz

        line = atom["line"]

        new_line = (
            line[:30]
            + f"{new_x:8.3f}"
            + f"{new_y:8.3f}"
            + f"{new_z:8.3f}"
            + line[54:]
        )

        new_lines.append(new_line)

    return new_lines


def main():

    if len(sys.argv) != 2:
        print("Usage:")
        print("python generate_structures.py <input.pdb>")
        sys.exit(1)

    pdb_file = sys.argv[1]

    if not os.path.isfile(pdb_file):
        print(f"Error: File '{pdb_file}' not found.")
        sys.exit(1)

    atom_records = read_pdb(pdb_file)

    if len(atom_records) == 0:
        print("Error: No standard amino-acid ATOM records found.")
        sys.exit(1)

    print(f"Input PDB: {pdb_file}")
    print(f"Amino-acid atoms found: {len(atom_records)}")

    while True:
        try:
            number_of_structures = int(
                input("How many new structures do you want to generate?: ")
            )

            if 10 <= number_of_structures <= 1000:
                break

            print("Please enter a number between 10 and 1000.")

        except ValueError:
            print("Please enter a valid integer.")

    output_folder = "new_structures"
    os.makedirs(output_folder, exist_ok=True)

    print("\nGenerating structures...")

    for i in range(1, number_of_structures + 1):

        new_structure = generate_structure(atom_records)

        output_file = os.path.join(
            output_folder,
            f"structure_{i:03d}.pdb"
        )

        with open(output_file, "w") as file:

            for line in new_structure:
                file.write(line + "\n")

            file.write("END\n")

    print("\nDone!")
    print(
        f"{number_of_structures} structures were saved in "
        f"'{output_folder}/'"
    )


if __name__ == "__main__":
    main()