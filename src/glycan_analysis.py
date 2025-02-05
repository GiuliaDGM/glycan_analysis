import re
import math
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_dihedrals
import plotly.graph_objects as go

# ============================================================================
# Parsing the .itp file
# ============================================================================

def extract_section_lines(lines, section_name):
    """
    Extract the lines corresponding to a given section name (e.g., '[ atoms ]').
    
    Parameters:
        lines (list of str): All lines from the .itp file.
        section_name (str): The section name to look for (e.g., '[ atoms ]').
    
    Returns:
        list of str: The lines belonging to that section (excluding comments, blank lines, etc.).
    """
    section_lines = []
    in_section = False
    
    for line in lines:
        # Check if we are entering the section
        if line.lower().startswith(section_name.lower()):
            in_section = True
            continue  # Skip the line that has the section name itself
        
        # If we meet another bracket or empty line, we leave the section
        if in_section and (line.startswith("[") or line == ""):
            break
        
        if in_section and not line.startswith(";"):
            # Not a comment, so capture it
            section_lines.append(line)
    
    return section_lines


def parse_atoms_section(atom_lines, amino_acids):
    """
    Parse the [ atoms ] section lines to extract atom data and identify glycan atoms.
    
    Parameters:
        atom_lines (list of str): Lines from the [ atoms ] section.
        amino_acids (set of str): Set of standard amino acid residue names.
    
    Returns:
        (dict, list):
          - atom_data: dict keyed by atom_index -> {resnr, resname, atom_name}
          - glycan_atoms: list of dict {atom_index, resnr, resname, atom_name}
    """
    atom_data = {}
    glycan_atoms = []
    
    for line in atom_lines:
        tokens = re.split(r'\s+', line)
        atom_index = int(tokens[0])
        resnr = int(tokens[2])
        resname = tokens[3]
        atom_name = tokens[4]

        atom_data[atom_index] = {
            "resnr": resnr,
            "resname": resname,
            "atom_name": atom_name
        }

        # If residue name is not in standard amino acids, treat as glycan
        if resname not in amino_acids:
            glycan_atoms.append({
                "atom_index": atom_index,
                "resnr": resnr,
                "resname": resname,
                "atom_name": atom_name
            })
    
    return atom_data, glycan_atoms


def parse_dihedrals_section(dihedral_lines):
    """
    Parse the [ dihedrals ] section lines to extract dihedral definitions.
    
    Parameters:
        dihedral_lines (list of str): Lines from the [ dihedrals ] section.
    
    Returns:
        list of tuples: Each tuple is (a1, a2, a3, a4) for a dihedral definition.
    """
    dihedrals = []
    for line in dihedral_lines:
        tokens = re.split(r'\s+', line)
        # The first 4 columns are the atom indices in a dihedral
        dihedral_atoms = tuple(map(int, tokens[:4]))
        dihedrals.append(dihedral_atoms)
    return dihedrals


def filter_dihedrals(dihedrals, atom_data):
    """
    Filter out dihedrals that have more than two hydrogens or do not cross multiple residues.
    
    Parameters:
        dihedrals (list of tuple): List of 4-atom dihedrals.
        atom_data (dict): atom_data keyed by atom_index, containing resname and atom_name.
    
    Returns:
        list of tuple: Filtered dihedrals.
    """
    filtered = []
    for dihedral in dihedrals:
        atom_names = [atom_data[a]["atom_name"] for a in dihedral]
        resnames = [atom_data[a]["resname"] for a in dihedral]
        carbon_count = sum(1 for name in atom_names if name.startswith("C"))
        hydrogen_count = sum(1 for name in atom_names if name.startswith("H"))
        # Keep if <= 2 hydrogens and involves at least 2 different residues & less than 3 carbons
        if hydrogen_count <= 2 and len(set(resnames)) > 1 and carbon_count < 3 :
            filtered.append(dihedral)
    return filtered


def create_glycan_entry(current_glycan, current_residues, dihedrals):
    """
    Create a glycan dictionary entry with the atoms, glycan_dihedrals, and protein_glycan_dihedrals.
    """
    return {
        "atoms": current_glycan,
        "glycan_dihedrals": [
            d for d in dihedrals if all(a in current_residues for a in d)
        ],
        "protein_glycan_dihedrals": [
            d for d in dihedrals if any(a in current_residues for a in d) 
            and not all(a in current_residues for a in d)
        ]
    }


def attach_protein_atoms_to_glycans(glycans, dihedrals, atom_data):
    """
    For any dihedral that crosses protein and glycan atoms, ensure the glycan 
    also lists the involved protein atoms in its 'atoms' list.
    
    Parameters:
        glycans (list of dict): The glycan structures built so far.
        dihedrals (list of tuple): Filtered dihedrals.
        atom_data (dict): Information about each atom (resnr, resname, atom_name).
    """
    for dihedral in dihedrals:
        for glycan in glycans:
            glycan_residues = {a['atom_index'] for a in glycan['atoms']}
            
            glycan_side = [a for a in dihedral if a in glycan_residues]
            protein_side = [a for a in dihedral if a not in glycan_residues]
            
            # If the dihedral crosses glycan/protein boundaries
            if glycan_side and protein_side:
                # Ensure we add the protein atoms to the glycan's atom list if not present
                for prot_atom_idx in protein_side:
                    if not any(x['atom_index'] == prot_atom_idx for x in glycan['atoms']):
                        glycan['atoms'].append({
                            "atom_index": prot_atom_idx,
                            "resnr": atom_data[prot_atom_idx]["resnr"],
                            "resname": atom_data[prot_atom_idx]["resname"],
                            "atom_name": atom_data[prot_atom_idx]["atom_name"]
                        })


def group_glycans(glycan_atoms, dihedrals, atom_data):
    """
    Group glycan atoms into separate glycans.
    We start a new glycan whenever we find a residue #1 with atom name 'C1'.
    
    Parameters:
        glycan_atoms (list): List of glycan atoms (dicts).
        dihedrals (list of tuple): Filtered dihedrals.
        atom_data (dict): Dictionary with info about all atoms.
    
    Returns:
        list of dict: Each dict represents a glycan with 'atoms', 'glycan_dihedrals', and 'protein_glycan_dihedrals'.
    """
    glycans = []
    current_glycan = []
    current_residues = set()

    for atom in glycan_atoms:
        # If we detect the start of a new glycan
        if atom['resnr'] == 1 and atom['atom_name'] == 'C1':
            # Save the previous glycan if it exists
            if current_glycan:
                glycans.append(create_glycan_entry(current_glycan, current_residues, dihedrals))
            
            # Reset for the new glycan
            current_glycan = []
            current_residues = set()

        # Add this atom to the current glycan
        current_glycan.append(atom)
        current_residues.add(atom['atom_index'])

    # After the loop, add the last glycan if it exists
    if current_glycan:
        glycans.append(create_glycan_entry(current_glycan, current_residues, dihedrals))
    
    # Optionally, we can attach protein atoms that connect to the glycan:
    attach_protein_atoms_to_glycans(glycans, dihedrals, atom_data)
    
    return glycans


def classify_dihedral(dihedral, atom_data):
    """
    Classify a dihedral as:
      - 'phi' if it has O5, C1, plus at least one O != O5 and one C != C1
      - 'psi' if it has C1, plus at least one O, one C != C1, and one H
      - 'unknown' otherwise.

    We also keep the requirement that it must span at least two residues.
    """
    atom_names = [atom_data[a]["atom_name"] for a in dihedral]
    resnrs = [atom_data[a]["resnr"] for a in dihedral]
    
    # Must involve at least two different residues
    if len(set(resnrs)) <= 1:
        return "unknown"

    # Count certain categories
    has_O5 = ("O5" in atom_names)
    has_C1 = ("C1" in atom_names)
    # any O not named O5
    has_other_O = any(n.startswith("O") and n != "O5" for n in atom_names)
    # any C not named C1
    has_other_C = any(n.startswith("C") and n != "C1" for n in atom_names)
    # any hydrogen
    has_H = any(n.startswith("H") for n in atom_names)

    # ϕ (phi): O5–C1–Ox–Cx
    # "O5" and "C1", plus at least one more "O" (not O5) and one more "C" (not C1)
    if has_O5 and has_C1 and has_other_O and has_other_C:
        return "phi"

    # ψ (psi): C1–Ox–Cx–Hx
    # "C1" in the set, plus an "O" (could be O5 or O2 or O3..., your choice),
    # a "C" not "C1," and at least one "H"
    if has_C1 and has_other_C and has_H:
        # We said "at least one O" in the set. 
        # If you want to exclude "O5" from psi, use has_other_O. 
        # If you allow O5, then do:
        if any(n.startswith("O") for n in atom_names):
            return "psi"
    
        # Check presence
    has_O5 = ("O5" in atom_names)
    has_C1 = ("C1" in atom_names)
    has_O6 = ("O6" in atom_names)
    has_C6 = ("C6" in atom_names)
    has_C5 = ("C5" in atom_names)
        # 1) Omega = O6, C6, C5, O5 all in the dihedral
    # => "omega"
    if has_O6 and has_C6 and has_C5 and has_O5:
        torsion_type = "omega"

    # If neither condition matches, default:
    return "unknown"


def classify_dihedrals_in_glycans(glycans, atom_data):
    """
    Classify dihedrals in each glycan as phi, psi, omega, or unknown, 
    and store the results in 'classified_glycan_dihedrals' and 
    'classified_protein_glycan_dihedrals' keys.
    
    Parameters:
        glycans (list of dict): Glycans with 'glycan_dihedrals' and 'protein_glycan_dihedrals'.
        atom_data (dict): Info about all atoms.
    """
    for glycan in glycans:
        # For each of these dihedral categories
        for category in ["glycan_dihedrals", "protein_glycan_dihedrals"]:
            dihedral_list = glycan[category]
            classified_list = []
            
            for dihedral in dihedral_list:
                d_type = classify_dihedral(dihedral, atom_data)
                classified_list.append({
                    "dihedral": dihedral,
                    "type": d_type,
                    "atoms": [
                        {
                            "atom_index": atom_idx,
                            "resnr": atom_data[atom_idx]["resnr"],
                            "resname": atom_data[atom_idx]["resname"],
                            "atom_name": atom_data[atom_idx]["atom_name"]
                        }
                        for atom_idx in dihedral
                    ]
                })
            
            # Store in glycan dictionary with a new key
            glycan[f"classified_{category}"] = classified_list


def parse_itp(itp_file, amino_acids):
    """
    Top-level function that orchestrates reading the .itp file,
    parsing atoms/dihedrals, filtering, grouping into glycans, 
    classifying dihedrals, and returning the results.

    Parameters:
        itp_file (str): Path to the .itp file.
        amino_acids (set): Set of standard amino acid residue names.

    Returns:
        (list, int):
            - list of glycan dictionaries,
            - total number of glycans.
    """
    # 1. Read the file lines
    with open(itp_file, 'r') as f:
        lines = [line.strip() for line in f]

    # 2. Extract relevant sections
    atom_lines = extract_section_lines(lines, "[ atoms ]")
    dihedral_lines = extract_section_lines(lines, "[ dihedrals ]")

    # 3. Parse the sections
    atom_data, glycan_atoms = parse_atoms_section(atom_lines, amino_acids)
    dihedrals = parse_dihedrals_section(dihedral_lines)

    # 4. Filter dihedrals
    dihedrals = filter_dihedrals(dihedrals, atom_data)

    # 5. Group glycan atoms into glycans
    glycans = group_glycans(glycan_atoms, dihedrals, atom_data)

    # 6. Classify dihedrals in glycans
    classify_dihedrals_in_glycans(glycans, atom_data)

    # 7. Return the final data
    num_glycans = len(glycans)
    return glycans, num_glycans, atom_data


def write_glycan_ranges(output_file, glycans, num_glycans, amino_acids):
    """
    Write glycan atom ranges, dihedrals, and count to an output file.
    
    Parameters:
        output_file (str): Path to the output file.
        glycans (list): List of glycans with their atoms and dihedrals.
        num_glycans (int): Total number of glycans.
    """
    with open(output_file, 'w') as file:
        file.write(f"# Total Glycans: {num_glycans}\n")
        file.write("# Glycan Atom Ranges and Dihedrals\n")
        for i, glycan in enumerate(glycans, start=1):
            # Filter out bridging protein atoms
            glycan_only_atoms = [
                atom for atom in glycan['atoms']
                if atom['resname'] not in amino_acids
            ]
            
            if not glycan_only_atoms:
                # If for some reason there's nothing left, skip or print something special
                file.write(f"Glycan {i}: (No glycan residues found)\n")
                continue
            
            # Compute range from glycan-only atoms
            glycan_indices = [atom['atom_index'] for atom in glycan_only_atoms]
            start = min(glycan_indices)
            end = max(glycan_indices)
            file.write(f"Glycan {i}: {start} - {end}\n")
            
            # Raw dihedrals (these lines remain the same as before)
            file.write(f"  Glycan-only Dihedrals: {glycan['glycan_dihedrals']}\n")
            file.write(f"  Protein-Glycan Dihedrals: {glycan['protein_glycan_dihedrals']}\n")
            
            # Classified dihedrals
            file.write(f"  Classified Glycan-only Dihedrals:\n")
            for classified in glycan["classified_glycan_dihedrals"]:
                file.write(f"    Dihedral: {classified['dihedral']} (Type: {classified['type']})\n")
                for atom in classified["atoms"]:
                    file.write(f"      {atom}\n")
            
            file.write(f"  Classified Protein-Glycan Dihedrals:\n")
            for classified in glycan["classified_protein_glycan_dihedrals"]:
                file.write(f"    Dihedral: {classified['dihedral']} (Type: {classified['type']})\n")
                for atom in classified["atoms"]:
                    file.write(f"      {atom}\n")

# ============================================================================
# Diheral calculation over a trajectory
# ============================================================================

def compute_dihedral_time_series_vectorized(universe, dihedrals_of_interest, start=None, stop=None, step=None):
    """
    Use MDAnalysis to compute dihedral angles for a list of dihedrals
    across an entire trajectory in 'universe' (which has .gro + .xtc).
    
    Parameters:
        universe (MDAnalysis.Universe): Universe loaded with .gro and .xtc.
        dihedrals_of_interest (list of tuple): (a1,a2,a3,a4) in 1-based indices
        start, stop, step (int): optional slicing arguments for frames.
        
    Returns:
        times (list of float): the time (ps) for each frame
        angle_data (dict): { (a1,a2,a3,a4) : [angles_in_degrees_per_frame] }
    """
    # Convert dihedrals from 1-based to 0-based
    zero_based = np.array([[d[0]-1, d[1]-1, d[2]-1, d[3]-1] for d in dihedrals_of_interest])
    angle_data = {d: [] for d in dihedrals_of_interest}
    times = []

    for ts in universe.trajectory[start:stop:step]:
        times.append(ts.time)  # in ps
        # shape = (n_atoms, 3)
        positions = universe.atoms.positions
        # shape for dihedral_positions is (N_dihedrals, 4, 3)
        dihedral_positions = positions[zero_based]

        # returns angles in radians
        angles_radians = calc_dihedrals(
            dihedral_positions[:,0,:],
            dihedral_positions[:,1,:],
            dihedral_positions[:,2,:],
            dihedral_positions[:,3,:],
        )
        angles_degrees = np.degrees(angles_radians)
        
        # store in angle_data
        for i, dih in enumerate(dihedrals_of_interest):
            angle_data[dih].append(angles_degrees[i])
    
    return times, angle_data

# ============================================================================
# Plotting functionality
# ============================================================================


def dihedral_glycans(ai, aj, ak, al, glycan_atom_sets):
    """
    Return a list of glycan indices that this dihedral (ai,aj,ak,al) touches.
    ai, aj, ak, al are 1-based indices, and each glycan_atom_sets[i] is a set of
    1-based atom indices. If any of these 4 are in glycan_atom_sets[i], 
    that glycan is considered 'involved'.
    """
    involved = []
    for i, atom_set in enumerate(glycan_atom_sets):
        # If any of the 4 atoms is in that atom_set, we consider 
        # that glycan is 'involved'
        if (ai in atom_set) or (aj in atom_set) or (ak in atom_set) or (al in atom_set):
            involved.append(i)
    return involved

def plot_dihedrals_by_glycan_and_type_plotly(
    times, 
    angle_data, 
    protein_glycan_dihedrals, 
    glycan_glycan_dihedrals, 
    glycans, 
    glycan_atom_sets
):
    """
    Group dihedrals by glycan (protein-glycan vs. glycan-glycan) 
    and create interactive Plotly plots for each glycan.

    Parameters:
        times (list of float): Simulation time for each frame, e.g. [0.0, 2.0, 4.0, ...].
        angle_data (dict): { (a1,a2,a3,a4) : [angles_per_frame, ...], ... }
        protein_glycan_dihedrals (list of tuple): Dihedrals bridging protein and glycan.
        glycan_glycan_dihedrals (list of tuple): Dihedrals purely within the glycan.
        glycans (list): Your list of glycan dicts, e.g. from parse_itp().
        glycan_atom_sets (list of set): Each set is the atom indices for a single glycan.

    Returns:
        None. (Displays interactive Plotly figures.)
    """

    # Same grouping logic as your old code -- no changes here
    def dihedral_glycans(ai, aj, ak, al, glycan_sets):
        involved = []
        for i, atom_set in enumerate(glycan_sets):
            if (ai in atom_set) or (aj in atom_set) or (ak in atom_set) or (al in atom_set):
                involved.append(i)
        return involved

    # 1) Group Protein-Glycan dihedrals by glycan
    protein_glycan_by_glycan = {i: [] for i in range(len(glycans))}
    for (ai, aj, ak, al) in protein_glycan_dihedrals:
        involved = dihedral_glycans(ai, aj, ak, al, glycan_atom_sets)
        for g_ind in involved:
            protein_glycan_by_glycan[g_ind].append((ai, aj, ak, al))
    
    # 2) Group Glycan-Glycan dihedrals by glycan
    glycan_glycan_by_glycan = {i: [] for i in range(len(glycans))}
    for (ai, aj, ak, al) in glycan_glycan_dihedrals:
        involved = dihedral_glycans(ai, aj, ak, al, glycan_atom_sets)
        for g_ind in involved:
            glycan_glycan_by_glycan[g_ind].append((ai, aj, ak, al))
    
    # 3) For each glycan, create an INTERACTIVE Plotly figure for protein-glycan dihedrals
    for g_ind, dihedrals_list in protein_glycan_by_glycan.items():
        if dihedrals_list:
            # Create a new figure for this glycan
            fig = go.Figure()
            for dih in dihedrals_list:
                angles = angle_data[dih]  # the timeseries for this dihedral
                fig.add_trace(
                    go.Scatter(
                        x=times, 
                        y=angles,
                        mode='lines',
                        name=f"Dihedral {dih}"  # legend label
                    )
                )
            fig.update_layout(
                title=f"Protein-Glycan Dihedrals - Glycan {g_ind+1}",
                xaxis_title="Time (ps)",
                yaxis_title="Angle (degrees)",
            )
            fig.show()  # Interactive plot appears in a browser or notebook

    # 4) Same approach for Glycan-Glycan dihedrals
    for g_ind, dihedrals_list in glycan_glycan_by_glycan.items():
        if dihedrals_list:
            fig = go.Figure()
            for dih in dihedrals_list:
                angles = angle_data[dih]
                fig.add_trace(
                    go.Scatter(
                        x=times,
                        y=angles,
                        mode='lines',
                        name=f"Dihedral {dih}"
                    )
                )
            fig.update_layout(
                title=f"Glycan-Glycan Dihedrals - Glycan {g_ind+1}",
                xaxis_title="Time (ps)",
                yaxis_title="Angle (degrees)",
            )
            fig.show()

# ============================================================================
# ============================================================================
# ============================================================================