from collections import Counter
import mdtraj as md
import numpy as np
import plotly.graph_objects as go

# ============================================================================
# ============================================================================
# PARSE .itp FILE AND GET MOLECULES + ATOM RANGES
# ============================================================================
# ============================================================================


# Helper function to find the right section in the .itp file
def find_section(line: str, section_name: str) -> bool:
    """
    Returns True if the given line is a header for the section 'section_name'.
    For example, if line is "[ atoms ]" and section_name is "atoms", returns True.
    """
    return line.strip().startswith("[") and section_name.lower() in line.lower()


def parse_atoms(itp_file_path: str):
    """
    Parses an itp file looking for the [ atoms ] section.
    
    It returns two lists:
      - proteins: a list of dictionaries, each with one key "Protein X" mapping to
                  a list of atoms (each represented as a dict with keys: nr, resnr, residu, atom).
      - glycans: a list of dictionaries, each with one key "Glycan X" mapping to a
                 list of atoms (same structure as proteins).
    
    A new protein is started when resnr == "1" and atom == "CA" (case-insensitive),
    and a new glycan is started when resnr == "1" and atom == "C1" (case-insensitive).
    
    Lines that are comments (starting with ";"), header lines, or empty are skipped.
    """
    proteins = []
    glycans = []
    current_protein = None
    current_glycan = None
    current_section = None

    # List of three-letter amino acid codes
    AA = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL"]

    with open(itp_file_path, "r") as file:
        for line in file:
            stripped = line.strip()
            # Skip empty lines or comment lines starting with ";"
            if not stripped or stripped.startswith(";"):
                continue
            
            # Check for header lines; update current_section accordingly
            if stripped.startswith("["):
                if find_section(stripped, "atoms"):
                    current_section = "atoms"
                else:
                    # Leaving the [ atoms ] section if another header is encountered
                    current_section = None
                continue

            # Only process lines when inside the [ atoms ] section
            if current_section == "atoms":
                tokens = stripped.split()
                # Process the line only if it starts with a number (expected atom entries)
                if tokens and tokens[0].isdigit():
                    nr = tokens[0]
                    # According to the header, tokens are:
                    # nr, type, resnr, residu, atom, cgnr, charge, mass, ...
                    # We only need: nr, resnr, residu, atom
                    resnr = tokens[2]
                    residu = tokens[3]
                    atom = tokens[4]
                    
                    # Determine if this entry is part of a protein (if residu in AA)
                    if residu.upper() in AA:
                        # Start a new protein if resnr == "1" and atom is "N"
                        if resnr == "1" and atom.upper() == "N":
                            current_protein = []
                            proteins.append({f"Protein {len(proteins) + 1}": current_protein})
                        if current_protein is not None:
                            current_protein.append({
                                "nr": nr,
                                "resnr": resnr,
                                "residu": residu,
                                "atom": atom
                            })
                    else:
                        # Otherwise, treat as glycan
                        # Start a new glycan if resnr == "1" and atom is "C1"
                        if resnr == "1" and atom.upper() == "C1":
                            current_glycan = []
                            glycans.append({f"Glycan {len(glycans) + 1}": current_glycan})
                        if current_glycan is not None:
                            current_glycan.append({
                                "nr": nr,
                                "resnr": resnr,
                                "residu": residu,
                                "atom": atom
                            })
    return proteins, glycans


def atom_ranges_of_molecules(proteins: list, glycans: list) -> list:
    """
    Given lists of proteins and glycans (each a list of dictionaries, where each dictionary
    has one key ("Protein X" or "Glycan X") mapping to a list of atom dictionaries), this function
    returns two lists:
    
      - protein_ranges: A list of dictionaries, each of the form { "Protein X": [nr_start, nr_end] },
        where nr_start is the "nr" of the first atom and nr_end is the "nr" of the last atom in that protein.
        
      - glycan_ranges: A list of dictionaries, each of the form { "Glycan X": [nr_start, nr_end] },
        where nr_start is the "nr" of the first atom and nr_end is the "nr" of the last atom in that glycan.
    
    Note: It is assumed that the list of atoms for each molecule is sorted in the order they appear in the file.
    """
    protein_ranges = []
    for protein in proteins:
        for key, atoms in protein.items():
            if atoms:  # Ensure there is at least one atom
                start_nr = int(atoms[0]["nr"])
                end_nr = int(atoms[-1]["nr"])
                protein_ranges.append({key: [start_nr, end_nr]})
    
    glycan_ranges = []
    for glycan in glycans:
        for key, atoms in glycan.items():
            if atoms:
                start_nr = int(atoms[0]["nr"])
                end_nr = int(atoms[-1]["nr"])
                glycan_ranges.append({key: [start_nr, end_nr]})
    
    return protein_ranges, glycan_ranges

# ============================================================================
# ============================================================================
# PARSE .itp FILE AND COLLECT DIHEDRALS
# ============================================================================
# ============================================================================


def parse_dihedrals(itp_file_path: str, glycan_ranges: list) -> list:
    """
    Parses the itp file to extract dihedrals from the [dihedrals] section.
    For each dihedral line (ignoring comments, headers, and blank lines),
    it extracts the four tokens ai, aj, ak, al (atom numbers).

    Then, for each glycan range (provided as a list of dictionaries,
    e.g. [{"Glycan 1": [nr_start, nr_end]}, {"Glycan 2": [nr_start, nr_end]}, ...]),
    it checks:
      - If all four atom numbers fall within the glycan's range, the dihedral is
        added to the glycan–glycan dihedrals for that glycan.
      - If three out of the four atoms are within the range, the dihedral is
        added to the protein–glycan dihedrals for that glycan.

    The function returns two lists:
      - glycan_glycan_dihedrals: a list of dictionaries mapping each glycan name
        to a list of its dihedrals [ai, aj, ak, al] (all four atoms from the glycan).
      - protein_glycan_dihedrals: a list of dictionaries mapping each glycan name
        to a list of dihedrals that have three atoms within the glycan range.
    """
    # Initialize dictionaries to collect dihedrals for each glycan
    glycan_glycan = {}
    protein_glycan = {}

    # Initialize flags for the section
    inside_dihedrals = False

    with open(itp_file_path, "r") as file:
        for line in file:
            stripped = line.strip()
            # Skip empty lines and comments (lines starting with ';')
            if not stripped or stripped.startswith(";"):
                continue

            # Check for header lines
            if stripped.startswith("["):
                if find_section(stripped, "dihedrals"):
                    inside_dihedrals = True
                else:
                    inside_dihedrals = False
                continue

            # Process lines only when inside the [dihedrals] section
            if inside_dihedrals:
                tokens = stripped.split()
                if tokens and tokens[0].isdigit():
                    try:
                        ai, aj, ak, al = map(int, tokens[:4])
                    except Exception:
                        continue  # Skip lines that cannot be parsed correctly
                        
                    # For each glycan range, check how many atoms fall inside
                    for glycan_dict in glycan_ranges:
                        for glycan_name, (start, end) in glycan_dict.items():
                            count = 0
                            for atom_nr in (ai, aj, ak, al):
                                if start <= atom_nr <= end:
                                    count += 1
                            if count == 4:
                                # All four atoms are in this glycan range.
                                # Append the dihedral to glycan_glycan for this glycan.
                                glycan_glycan.setdefault(glycan_name, []).append([ai, aj, ak, al])
                            elif count >= 1 and count <= 3:
                                # At least one atom is in the range, but not all.
                                # Append the dihedral to protein_glycan for this glycan.
                                protein_glycan.setdefault(glycan_name, []).append([ai, aj, ak, al])
    
    # Convert the dictionaries into lists of single-key dictionaries for output.
    glycan_glycan_dihedrals = []
    for glycan_dict in glycan_ranges:
        for glycan_name in glycan_dict.keys():
            glycan_glycan_dihedrals.append({glycan_name: glycan_glycan.get(glycan_name, [])})

    protein_glycan_dihedrals = []
    for glycan_dict in glycan_ranges:
        for glycan_name in glycan_dict.keys():
            protein_glycan_dihedrals.append({glycan_name: protein_glycan.get(glycan_name, [])})

    return glycan_glycan_dihedrals, protein_glycan_dihedrals


# ============================================================================
# ============================================================================
# FILTER AND CLASSIFY DIHEDRALS
# ============================================================================
# ============================================================================


def convert_nr_to_atom(dihedral_nr_list: list, glycans: list, proteins: list) -> list:
    """
    Given a list of atom numbers (e.g., [ai, aj, ak, al]) and lists of
    glycan and protein molecules (as returned by parse_atoms), this function
    searches both molecules for each atom number and retrieves the corresponding
    "atom" and "resnr" fields.
    
    Returns:
      - dihedral_atom_list: a list of the atom names corresponding to the provided numbers.
      - dihedral_resnr_list: a list of the residue numbers (resnr) corresponding to the provided numbers.
      
    The function searches proteins first; if not found, it then searches glycans.
    """
    dihedral_atom_list = []
    dihedral_resnr_list = []
    
    def find_atom_by_nr(molecule_list, target_nr: int):
        for molecule in molecule_list:
            # Each molecule is a dictionary with a single key mapping to a list of atom dictionaries.
            for atoms in molecule.values():
                for atom in atoms:
                    if int(atom["nr"]) == target_nr:
                        return atom["atom"], atom["resnr"]
        return None, None

    for nr in dihedral_nr_list:
        target_nr = int(nr)
        atom_name, resnr = find_atom_by_nr(proteins, target_nr)
        if atom_name is None:
            atom_name, resnr = find_atom_by_nr(glycans, target_nr)
        dihedral_atom_list.append(atom_name)
        dihedral_resnr_list.append(int(resnr))
    
    return dihedral_atom_list, dihedral_resnr_list

# Example usage:
# Suppose we have a dihedral_nr_list = [11, 5, 13, 12]
# dihedral_atoms, dihedral_residues = convert_nr_to_atom([11, 5, 13, 12], glycans, proteins)
# print("Dihedral atoms:", dihedral_atoms)
# print("Dihedral residues:", dihedral_residues)


def is_glyco_phi(atoms: list) -> bool:
    """
    Glycosilic phi: The dihedral must include "O5" and "C1".
    After removing one occurrence each of "O5" and "C1", the remaining two atoms
    must include one that starts with "O" and one that starts with "C".
    """
    if len(atoms) != 4:
        return False
    if "O5" not in atoms or "C1" not in atoms:
        return False
    remaining = atoms.copy()
    remaining.remove("O5")
    remaining.remove("C1")
    if len(remaining) != 2:
        return False
    return any(a.startswith("O") for a in remaining) and any(a.startswith("C") for a in remaining)

def is_glyco_psi(atoms: list) -> bool:
    """
    Glycosilic psi: The dihedral must include "C1".
    After removing one "C1", among the remaining three atoms exactly one
    must start with "O" and the other two must start with "C".
    """
    if len(atoms) != 4:
        return False
    if "C1" not in atoms:
        return False
    remaining = atoms.copy()
    remaining.remove("C1")
    if len(remaining) != 3:
        return False
    count_O = sum(1 for a in remaining if a.startswith("O"))
    count_C = sum(1 for a in remaining if a.startswith("C"))
    return count_O == 1 and count_C == 2


def is_glyco_omega(atoms: list) -> bool:
    """
    Glycosidic omega: Expected atoms (in any order) must exactly be:
      {"O6", "C6", "C5", "O5"}
    """
    return set(atoms) == {"O6", "C6", "C5", "C1"}

def is_prot_gly_phi(atoms: list) -> bool:
    """
    Protein–glycan phi: Expected atoms (in any order) must exactly be:
      {"O5", "C1", "ND2", "CG"}
    """
    return set(atoms) == {"O5", "C1", "ND2", "CG"}

def is_prot_gly_psi(atoms: list) -> bool:
    """
    Protein–glycan psi: Expected atoms (in any order) must exactly be:
      {"C1", "ND2", "CG", "CB"}
    """
    return set(atoms) == {"C1", "ND2", "CG", "CB"}


def classify_dihedrals(dihedral_atom_list: list, dihedral_residues: list) -> str:
    """
    Classify the dihedral based on the atom names and the residue numbers.
    
    The classification uses two criteria:
    
      1. The atom-names pattern:
         - Glycosidic phi: O5 - C1 - O1 - C'x  
           (first three atoms are "O5", "C1", "O1"; fourth atom starts with "C")
         - Glycosidic psi: C1 - O1 - C'x - C'x-1  
           (first two atoms "C1" and "O1"; third and fourth atoms start with "C",
         - Glycosidic omega: O1 - C'6 - C'5 - O'5  
         - Protein–glycan phi: O5 - C1 - N'D2 - C'G  
         - Protein–glycan psi: C1 - N'D2 - C'G - C'B

         Note: The apostrophe (') indicates that the atom belongs to a different residue.
         
      2. The residue counts (provided in dihedral_residues):
         - For glycosidic phi, expect 3 atoms from the current glycan residue and 1 from another.
         - For glycosidic psi, expect a 2–2 split.
         - For glycosidic omega, expect 1 from the current and 3 from another.
         - For protein–glycan phi, expect a 2–2 split.
         - For protein–glycan psi, expect 1 current and 3 other.
         
    The function returns one of "phi", "psi", "omega", or "unknown".
    """
    # Normalize atom names: remove extra spaces and convert to uppercase.
    atoms = [atom.strip().upper() if atom is not None else "" for atom in dihedral_atom_list]
    
    # Count the frequency of each residue number.
    freq = Counter(dihedral_residues)
    if freq:
        current_count = max(freq.values())
        other_count = 4 - current_count
    else:
        current_count, other_count = 0, 0
    
    # Now check patterns along with residue frequency criteria.
        # For glycosidic omega, expect current_count == 3 and other_count == 1.
    if is_glyco_omega(atoms) and current_count == 3 and other_count == 1:
        return "omega"
    # For glycosidic phi, we need current_count == 3 (3 atoms from same glycan residue) and other_count == 1.
    if is_glyco_phi(atoms) and current_count == 3 and other_count == 1:
        return "phi"
    # For glycosidic psi, expect a 2-2 split.
    if is_glyco_psi(atoms) and current_count == 2 and other_count == 2:
        return "psi"
    # For protein-glycan phi, expect a 2-2 split.
    if is_prot_gly_phi(atoms) and current_count == 2 and other_count == 2:
        return "phi"
    # For protein-glycan psi, expect current_count == 3 and other_count == 1.
    if is_prot_gly_psi(atoms) and current_count == 3 and other_count == 1:
        return "psi"
    
    # If none of the criteria match, return "unknown".
    return "unknown"

# Example usage:
# For a glycosidic phi example:
# dihedral_atom_list = ["O5", "C1", "O1", "C8"]
# dihedral_residues  = [101, 101, 101, 102] # (three atoms from residue 101, one from 102)
# print(classify_dihedrals(dihedral_atom_list, dihedral_residues)) # should return "phi"

# For a protein–glycan psi example:
# dihedral_atom_list = ["C1", "ND2", "CG", "CB"]
# dihedral_residues  = [50, 51, 51, 51] # (1 from residue 50, 3 from residue 51 would be protein–glycan psi;
# print(classify_dihedrals(dihedral_atom_list, dihedral_residues)) # should return "psi"


def filter_dihedral_list(dihedral_list: list, proteins: list, glycans: list) -> list:
    """
    Filters a dihedral list (a list of dictionaries, each keyed by a glycan name)
    based on the following criteria:
    
      1. Convert each dihedral (a list of four atom numbers) into atom names and residue numbers
         using convert_nr_to_atom(dihedral_nr_list, glycans, proteins).
      2. Exclude dihedrals that include any hydrogen atoms (i.e. an atom name starting with "H").
      3. Ensure that the dihedral involves at least two different residues.
      4. Classify the dihedral using classify_dihedrals(dihedral_atom_list, dihedral_residues).
         Only keep those classified as "phi", "psi", or "omega".
    
    The function returns a list of dictionaries with the structure:
    
      [ { "Glycan 1": { "phi": [ [ai,aj,ak,al], ... ],
                        "psi": [ [ai,aj,ak,al], ... ],
                        "omega": [ [ai,aj,ak,al], ... ] } },
        { "Glycan 2": { ... } },
        ... ]
    """
    filtered = []
    for glycan_dict in dihedral_list:
        for glycan_name, dihedrals in glycan_dict.items():
            classified = {}
            for d_nr_list in dihedrals:
                # Convert dihedral numbers to atom names and residue numbers.
                dihedral_atoms, dihedral_resnrs = convert_nr_to_atom(d_nr_list, glycans, proteins)
                # Condition 1: Exclude if more than one atom is hydrogen.
                has_h = sum(1 for atom in dihedral_atoms if atom and atom.startswith("H")) > 1
                # Condition 2: Must involve at least two different residues.
                distinct_residues = len(set(dihedral_resnrs)) >= 2
                if (not has_h) and distinct_residues:
                    # Classify the dihedral.
                    classification = classify_dihedrals(dihedral_atoms, dihedral_resnrs)
                    if classification in ["phi", "psi", "omega"]:
                        if classification in classified:
                            classified[classification].append(d_nr_list)
                        else:
                            classified[classification] = [d_nr_list]
            filtered.append({glycan_name: classified})
    
    return filtered


def filter_dihedrals(
        protein_glycan_dihedrals: list,
        glycan_glycan_dihedrals: list,
        proteins: list,
        glycans: list) -> list:
    """
    Filters both protein–glycan and glycan–glycan dihedral lists.
    
    Each input list is assumed to have the structure:
      [ { "Glycan X": [ [ai,aj,ak,al], [ai,aj,ak,al], ... ] }, ... ]
      
    The function applies filter_dihedral_list() to each input and returns:
    
      - filtered_protein_glycan_dihedrals: filtered and classified protein–glycan dihedrals.
      - filtered_glycan_glycan_dihedrals: filtered and classified glycan–glycan dihedrals.
    """
    filtered_protein_glycan_dihedrals = filter_dihedral_list(protein_glycan_dihedrals, proteins, glycans)
    filtered_glycan_glycan_dihedrals = filter_dihedral_list(glycan_glycan_dihedrals, proteins, glycans)
    
    return filtered_protein_glycan_dihedrals, filtered_glycan_glycan_dihedrals


# ============================================================================
# ============================================================================
# PARSE .gro FILE AND GET DIHEDRAL COORDINATES
# ============================================================================
# ============================================================================


def parse_gro_coordinates(gro_file_path: str) -> dict:
    """
    Parse a .gro file and return a dictionary mapping atom number (int) to its coordinates (x, y, z) as floats.
    The .gro file typically has the following format for each atom line:
      residue_number  residue_name  atom_name  atom_number  x  y  z [velocities...]
    We ignore velocities.
    """
    coord_dict = {}
    with open(gro_file_path, "r") as f:
        lines = f.readlines()
    
    # The first line is a title and the last two lines are the number of atoms and the box vectors.
    # The atom lines are in between.
    # One common approach is to read the second line as the number of atoms:
    try:
        n_atoms = int(lines[1].strip())
    except Exception:
        raise ValueError("Could not parse the number of atoms from the gro file.")
    
    # Atom lines: from line index 2 up to 2+n_atoms
    atom_lines = lines[2:2+n_atoms]
    
    for line in atom_lines:
        # .gro file columns are fixed-width. We assume a standard formatting.
        # Format: residue number (5 chars), residue name (5 chars), atom name (5 chars), atom number (5 chars), x (8 chars), y (8 chars), z (8 chars)
        # We can use string slicing.
        try:
            atom_number = int(line[15:20].strip())
            x = float(line[20:28].strip())
            y = float(line[28:36].strip())
            z = float(line[36:44].strip())
            coord_dict[atom_number] = (x, y, z)
        except Exception as e:
            print("Error parsing line:", line)
            continue
    return coord_dict


def convert_dihedrals(filtered_dihedrals: list, coord_dict:dict) -> list:
        """
        Given a filtered dihedral list (for either protein–glycan or glycan–glycan),
        replace each dihedral's atom numbers with their coordinates.
        """
        converted = []
        for glycan_dict in filtered_dihedrals:
            for glycan_name, dihedral_types in glycan_dict.items():
                new_types = {}
                for dihedral_type, dihedrals in dihedral_types.items():
                    new_dihedrals = []
                    for d in dihedrals:
                        # For each dihedral (a list of four atom numbers)
                        coords = []
                        for atom_nr in d:
                            coord = coord_dict.get(int(atom_nr))
                            if coord is None:
                                print(f"Warning: Atom number {atom_nr} not found in gro file.")
                                coord = (None, None, None)
                            coords.append(coord)
                        new_dihedrals.append(coords)
                    new_types[dihedral_type] = new_dihedrals
                converted.append({glycan_name: new_types})
        return converted


def parse_gro_file(gro_file_path: str,
                   filtered_glycan_glycan_dihedrals: list,
                   filtered_protein_glycan_dihedrals: list) -> list:
    """
    Given a .gro file and filtered dihedral structures (for glycan–glycan and protein–glycan),
    this function retrieves the coordinates for each dihedral from the .gro file.
    
    The filtered structures have the format:
      filtered_dihedrals = [
         { "Glycan 1": { "phi": [ [ai, aj, ak, al], [ai, aj, ak, al], ... ],
                         "psi": [ [ai, aj, ak, al], ... ],
                         "omega": [ [ai, aj, ak, al], ... ]
                       }
         },
         { "Glycan 2": { ... } },
         ...
      ]
    
    The function returns two lists with the same structure but with each dihedral replaced by a list of coordinate quadruplets:
    
      glycosidic_coordinates_prot_glycan = [ { "Glycan 1": { "phi": [ [(x,y,z), ...], ... ],
                                                              "psi": [ ... ],
                                                              "omega": [ ... ]
                                                            }
                                             },
                                             ... ]
      glycosidic_coordinates_glycan_glycan = [ { "Glycan 1": { ... }, ... } ]
    """
    # First, parse the gro file into a dictionary mapping atom numbers to coordinates.
    coord_dict = parse_gro_coordinates(gro_file_path)

    glycosidic_coordinates_prot_glycan = convert_dihedrals(filtered_protein_glycan_dihedrals, coord_dict)
    glycosidic_coordinates_glycan_glycan = convert_dihedrals(filtered_glycan_glycan_dihedrals, coord_dict)
    
    return glycosidic_coordinates_prot_glycan, glycosidic_coordinates_glycan_glycan


# ============================================================================
# ============================================================================
# CALCULATE DIHEDRALS AND PLOT
# ============================================================================
# ============================================================================


def calculate_dihedral(traj_file_path: str, gro_file_path: str, atom_indices: tuple) -> np.ndarray:
    """
    Calculate the dihedral angle over time for the given four atom indices.
    
    Parameters:
      traj_file_path (str): Path to the trajectory file (which contains or is associated with topology).
      atom_indices (tuple): A tuple of four atom indices (0-based indices expected by MDTraj).
    
    Returns:
      np.ndarray: A 1D array of dihedral angle values in radians (one value per frame).
    """
    traj = md.load(traj_file_path, top=gro_file_path)
    dihedral_values = md.compute_dihedrals(traj, [atom_indices])
    # dihedral_values is a 2D array of shape (n_frames, 1); flatten to 1D.
    return dihedral_values.flatten()


def plot_dihedrals_for_glycan(glycan_name: str, dihedral_dict: dict, traj_file_path: str, gro_file_path: str, label: str) -> go.Figure:
    """
    For one glycan (with a given name and dihedral dictionary), compute the time series
    for each dihedral and return a Plotly Figure.

    Parameters:
      glycan_name (str): Name of the glycan (e.g., "Glycan 1").
      dihedral_dict (dict): Dict with keys 'phi', 'psi', 'omega' -> lists of [ai,aj,ak,al].
      traj_file_path (str): Path to the trajectory file.
      label (str): Label string ("Protein-Glycan" or "Glycan-Glycan") for the figure title.

    Returns:
      go.Figure: Plotly figure with dihedral time series (shifted to [0,360] degrees).
    """
    fig = go.Figure()

    for dihedral_type, dihedrals in dihedral_dict.items():
        for i, d in enumerate(dihedrals):
            # Convert atom numbers to a tuple of integers
            atom_indices = tuple(int(a) for a in d)

            # Calculate dihedral angles in radians
            angles_rad = calculate_dihedral(traj_file_path, gro_file_path, atom_indices)

            # Convert from radians to degrees
            angles_deg = np.degrees(angles_rad)

            # Shift from [-180,180) to [0,360)
            angles_deg = (angles_deg + 360) % 360

            frames = np.arange(len(angles_deg))
            trace_name = f"{dihedral_type} dihedral {i+1}"

            fig.add_trace(
                go.Scatter(
                    x=frames,
                    y=angles_deg,
                    mode='lines',
                    name=trace_name
                )
            )

    title = f"{label} Dihedrals for {glycan_name}"
    fig.update_layout(
        title=title,
        xaxis_title="Frame",
        yaxis_title="Dihedral Angle (degrees)",
        template="plotly_white"
    )

    figure_name = f"dihedrals_{label.lower().replace('-','_')}_{glycan_name.replace(' ', '_').lower()}"
    # Save the figure with specified width, height, and scale factor
    fig.write_image(
        f"results/IL3_GLY_R1/{figure_name}.png",
        width=1920,    # width in pixels
        height=1080,   # height in pixels
        scale=2        # scale factor (output image will be 2x the specified dimensions)
    )

    return fig


def plot_all_dihedrals(traj_file_path: str,
                       gro_file_path: str,
                       filtered_protein_glycan_dihedrals: list,
                       filtered_glycan_glycan_dihedrals: list):
    """
    For each glycan in the filtered protein–glycan and glycan–glycan dihedral lists,
    compute and show a Plotly figure with the dihedral time series.
    
    The function creates:
      - A page for each glycan in protein–glycan dihedrals.
      - A page for each glycan in glycan–glycan dihedrals.
    """
    # Protein-Glycan dihedrals pages.
    for glycan_dict in filtered_protein_glycan_dihedrals:
        for glycan_name, dihedral_dict in glycan_dict.items():
            fig = plot_dihedrals_for_glycan(glycan_name, dihedral_dict, traj_file_path, gro_file_path, "Protein-Glycan")
            # Display the figure (this will open it in your default browser or in a notebook cell).
            fig.show()
    
    # Glycan-Glycan dihedrals pages.
    for glycan_dict in filtered_glycan_glycan_dihedrals:
        for glycan_name, dihedral_dict in glycan_dict.items():
            fig = plot_dihedrals_for_glycan(glycan_name, dihedral_dict, traj_file_path, gro_file_path, "Glycan-Glycan")
            fig.show()


# ============================================================================
# ============================================================================
# REGROUP THE GLYCAN DIHEDRALS DATA
# ============================================================================
# ============================================================================


def regroup_glycan_data(glycosidic_coordinates_glycan_glycan: list, glycosidic_coordinates_prot_glycan: list) -> dict:
    system_glycans_and_dihedrals = {}
    
    # Process glycosidic_coordinates_glycan_glycan
    for glycan_entry in glycosidic_coordinates_glycan_glycan:
        for glycan_name, dihedrals in glycan_entry.items():
            # Initialize dictionary for this glycan if not present.
            if glycan_name not in system_glycans_and_dihedrals:
                system_glycans_and_dihedrals[glycan_name] = {"phi": [], "psi": [], "omega": []}
            # Loop through expected dihedral keys
            for key in ['phi', 'psi', 'omega']:
                if key in dihedrals:
                    # Extend the list with the list(s) of coordinates
                    system_glycans_and_dihedrals[glycan_name][key].extend(dihedrals[key])
                    
    # Process glycosidic_coordinates_prot_glycan
    for glycan_entry in glycosidic_coordinates_prot_glycan:
        for glycan_name, dihedrals in glycan_entry.items():
            if glycan_name not in system_glycans_and_dihedrals:
                system_glycans_and_dihedrals[glycan_name] = {"phi": [], "psi": [], "omega": []}
            for key in ['phi', 'psi', 'omega']:
                if key in dihedrals:
                    system_glycans_and_dihedrals[glycan_name][key].extend(dihedrals[key])
                    
    return system_glycans_and_dihedrals


# ============================================================================
# ============================================================================
# ============================================================================
