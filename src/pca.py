# Imports
import numpy as np
import mdtraj as md
from sklearn.decomposition import PCA
import plotly.graph_objects as go

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
    # Flatten the result from shape (n_frames, 1) to (n_frames,)
    return dihedrals.flatten()


# ============================================================================
# ============================================================================
# CREATE FEATURE MATRIX AND PERFORM PCA
# ============================================================================
# ============================================================================


def create_feature_matrix_from_system(system_data: dict, traj_file_path: str, top_file_path: str) -> np.ndarray:
    """
    Given a dictionary for one glycan (with keys "phi", "psi", "omega") where each
    key maps to a list of dihedrals (each dihedral is a list of four atom indices, which
    will be used with the trajectory to compute a time series of dihedral angles),
    this function constructs a feature matrix X.

    For each dihedral, we:
      - Compute the time series of dihedral angles (in radians) using compute_dihedral_angles.
      - Compute the sin and cos values.
    For each frame, we then concatenate these features (first all sin's, then all cos's).

    Returns:
      X (np.ndarray): Feature matrix of shape (n_frames, 2*m_total) where m_total is the total number of dihedrals.
    """
    feature_series = []  # list to hold 1D arrays (n_frames,) for each feature
    # Loop over each dihedral type
    for d_type in ["phi", "psi", "omega"]:
        if d_type not in system_data:
            continue
        # Each dihedral in this category is assumed to be a list of four atom indices.
        # (If you have multiple dihedrals, system_data[d_type] is a list of quadruplets.)
        for dihedral in system_data[d_type]:
            # Compute the time series for this dihedral.
            # dihedral is assumed to be a tuple/list of 4 atom indices.
            angles = compute_dihedral_angles(traj_file_path, top_file_path, tuple(dihedral))
            # Compute sin and cos features.
            sin_vals = np.sin(angles)
            cos_vals = np.cos(angles)
            feature_series.append(sin_vals)
            feature_series.append(cos_vals)
    if not feature_series:
        raise ValueError("No dihedral data found in system_data.")
    # Stack features column-wise. Each array in feature_series has shape (n_frames,)
    X = np.column_stack(feature_series)
    return X


def compute_pca(X: np.ndarray, n_components: int = 2) -> np.ndarray:
    """
    Perform PCA on the feature matrix X and return the projection T.
    
    Parameters:
      X (np.ndarray): Feature matrix.
      n_components (int): Number of components (default 2).
    
    Returns:
      T (np.ndarray): 2D projection of shape (n_frames, n_components).
    """
    pca = PCA(n_components=n_components)
    T = pca.fit_transform(X)
    return T


def plot_pca(T: np.ndarray, title="PCA of Glycan Conformations"):
    """
    Plot the 2D PCA projection T using Plotly.
    
    Parameters:
      T (np.ndarray): 2D array (n_frames, 2) from PCA.
      title (str): Title for the plot.
    """
    n_frames = T.shape[0]
    frames = np.arange(n_frames)
    fig = go.Figure(data=go.Scatter(
        x=T[:, 0],
        y=T[:, 1],
        mode='markers',
        marker=dict(
            size=5,
            color=frames,  # Optionally, color-code by frame index
            colorscale='Viridis',
            colorbar=dict(title="Frame")
        )
    ))
    fig.update_layout(
        title=title,
        xaxis_title="PC1",
        yaxis_title="PC2",
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
    
    fig.show()


# ============================================================================
# ============================================================================
# REGROUP EVERYTHING
# ============================================================================
# ============================================================================


def compute_and_plot_pca(system_glycans_and_dihedrals: dict, top_file: str, traj_file: str):
    """
    For each glycan in system_glycans_and_dihedrals (a dictionary where keys are glycan names and values
    are dictionaries with keys 'phi', 'psi', 'omega', each mapping to a list of dihedral quadruplets
    as atom index tuples), this function:
      - Computes the dihedral time series for each dihedral using MDTraj.
      - Creates a feature matrix from sin and cos values of all dihedrals.
      - Performs PCA on the feature matrix.
      - Plots the 2D PCA projection using Plotly.

    Parameters:
      system_glycans_and_dihedrals (dict): The input dihedral data.
      top_file (str): Topology file path (e.g., .gro or .pdb).
      traj_file (str): Trajectory file path (e.g., .xtc).
    """
    for glycan_name, dihedral_dict in system_glycans_and_dihedrals.items():
        print(f"Processing {glycan_name} ...")
        # Prepare feature matrix for this glycan.
        X = create_feature_matrix_from_system(dihedral_dict, traj_file, top_file)
        print(f"Feature matrix for {glycan_name} has shape: {X.shape}")
        # Perform PCA
        T = compute_pca(X, n_components=2)
        # Plot the PCA projection.
        plot_title = f"PCA of {glycan_name} Conformations"
        plot_pca(T, title=plot_title)


# ============================================================================
# ============================================================================
# ============================================================================
