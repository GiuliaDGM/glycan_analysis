import numpy as np
import mdtraj as md
import math
import plotly.express as px

# ============================================================================
# ============================================================================
# COMPUTE DIHEDRAL ANGLES
# ============================================================================
# ============================================================================

def compute_dihedral_angles(traj_file_path: str, top_file_path: str, atom_indices: tuple) -> np.ndarray:
    """
    Calculate the dihedral angle over time for the given four atom indices
    using MDTraj's compute_dihedrals.

    Parameters:
      traj_file_path (str): Path to the trajectory file (e.g. .xtc).
      top_file_path (str): Path to the topology file (e.g. .gro or .pdb) containing connectivity.
      atom_indices (tuple): A tuple of four atom indices (0-based indices for MDTraj).

    Returns:
      np.ndarray: A 1D array of dihedral angle values in radians (one value per frame).
    """
    traj = md.load(traj_file_path, top=top_file_path)
    dihedrals = md.compute_dihedrals(traj, [atom_indices])
    return dihedrals.flatten()


# ============================================================================
# ============================================================================
# MAP DIHEDRAL ANGLES TO CATEGORICAL STATES
# ============================================================================
# ============================================================================


def assign_digit_state(angle:float, angle_type:str)-> str:
    """
    Assign a categorical state to a given dihedral angle based on IUPAC intervals.
    
    Parameters:
        angle (float): Dihedral angle in radians.
        angle_type (str): 'phi' 'psi'
        
    Returns:
        str: The state label (e.g., "C", "G+", "A+", "T", "A-", "G-").
    """
    # For phi and psi angles:
    if angle_type == 'phi' or angle_type == 'psi':
        if -0.52 <= angle < 0.52:
            return "C"
        elif 0.52 <= angle < 1.57:
            return "G+"
        elif 1.57 <= angle < 2.62:
            return "A+"
        elif (2.62 <= angle <= math.pi) or (-math.pi <= angle < -2.62):
            return "T"
        elif -2.62 <= angle < -1.57:
            return "A-"
        elif -1.57 <= angle < -0.52:
            return "G-"
    # For omega angles
    # elif angle_type == 'omega':
    #     if -2.62 <= angle < 0:
    #         return "gg"
    #     elif 0 <= angle < 2.62:
    #         return "gt"
    #     elif (2.62 <= angle <= math.pi) or (-math.pi <= angle < -2.62):
    #         return "tg"
    return "Undefined"


def compute_conformer_string(dihedral_angles:list, angle_types:list) -> str:
    """
    Build a conformer string for a glycan from its dihedral angles.
    
    Parameters:
        dihedral_angles (list of float): List of dihedral angles for one frame.
        angle_types (list of str): List indicating the type of each angle ('phi_psi' or 'omega').
        
    Returns:
        str: A string representing the conformer (e.g., "G-A+T-").
    """
    states = []
    for angle, angle_type in zip(dihedral_angles, angle_types):
        state = assign_digit_state(angle, angle_type)
        states.append(state)
    return " - ".join(states)


# ============================================================================
# ============================================================================
# PLOT PSEUDO RAMACHANDRAN PLOT
# ============================================================================
# ============================================================================


def plot_ramachandran(glycan_name:str, phi_angles:np.ndarray, psi_angles:np.ndarray, state_labels=None)-> None:
    """
    Plot a pseudo Ramachandran plot..
    
    Parameters:
        phi_angles (np.ndarray): Array of phi angles (in radians).
        psi_angles (np.ndarray): Array of psi angles (in radians).
        state_labels (list, optional): List of state labels to color-code the points.
    """
    # Convert angles from radians to degrees
    phi_deg = np.degrees(phi_angles)
    psi_deg = np.degrees(psi_angles)
    
    title = f"Pseudo Ramachandran Plot for {glycan_name}"
    # Create a scatter plot
    fig = px.scatter(x=phi_deg, y=psi_deg, color=state_labels,
                     labels={'x': 'Phi (φ)', 'y': 'Psi (ψ)'},
                     title=title)
    
    # Overlay IUPAC bin boundaries (example values in degrees)
    bin_boundaries = [math.degrees(0.52), -math.degrees(0.52),
                      math.degrees(1.57), -math.degrees(1.57),
                      math.degrees(2.62), -math.degrees(2.62),
                      math.degrees(math.pi), -math.degrees(math.pi)]
    
    for boundary in bin_boundaries:
        fig.add_vline(x=boundary, line_dash="dash", line_color="grey")
        fig.add_hline(y=boundary, line_dash="dash", line_color="grey")
    

        fig.update_layout(
        title=title,
        xaxis_title="Frame",
        yaxis_title="Dihedral Angle (degrees)",
        template="plotly_white"
    )

    figure_name = f"{title.replace(' ', '_').lower()}"
    # Save the figure with specified width, height, and scale factor
    fig.write_image(
        f"results/IL3_GLY_R1/{figure_name}.png",
        width=1920,    # width in pixels
        height=1080,   # height in pixels
        scale=2        # scale factor (output image will be 2x the specified dimensions)
    )
    fig.show()


def compute_and_plot_ramachandran(system_glycans_and_dihedrals:dict, top_file:str, traj_file:str)-> None:
    """
    High-level function to compute dihedral angles, assign states, build conformer strings,
    and plot pseudo Ramachandran plots for each glycan.
    
    Parameters:
        system_glycans_and_dihedrals (dict): Dictionary containing glycans and their associated atom index sets.
        top_file (str): Path to the topology file.
        traj_file (str): Path to the trajectory file.
    """
    for glycan, angles_data in system_glycans_and_dihedrals.items():
        print(f"Processing {glycan}...")
        phi_angles_list = []  # To collect phi values for plotting
        psi_angles_list = []  # To collect psi values for plotting
        conformer_strings = []  # To store conformer strings for each frame
        
        # angles_data has keys 'phi' and 'psi', each with a list of atom_indices per linkage.
        for linkage_idx, (phi_indices, psi_indices) in enumerate(zip(angles_data['phi'], angles_data['psi'])):
            # Compute dihedral angles for the current linkage
            phi_angles = compute_dihedral_angles(traj_file, top_file, phi_indices)
            psi_angles = compute_dihedral_angles(traj_file, top_file, psi_indices)
            
            # Append angles for the plot
            phi_angles_list.append(phi_angles)
            psi_angles_list.append(psi_angles)
            
            # For each frame, compute the conformer string for this linkage
            angle_types = ['phi', 'psi']
            for frame_phi, frame_psi in zip(phi_angles, psi_angles):
                conformer_str = compute_conformer_string([frame_phi, frame_psi], angle_types)
                conformer_strings.append(conformer_str)
        
        # Flatten the collected angle arrays (if necessary) for plotting
        phi_all = np.concatenate(phi_angles_list)
        psi_all = np.concatenate(psi_angles_list)
        
        # Plot the pseudo Ramachandran plot; pass conformer_strings if you wish to color-code points.
        plot_ramachandran(glycan, phi_all, psi_all, state_labels=conformer_strings)

# ============================================================================
# ============================================================================
# ============================================================================
