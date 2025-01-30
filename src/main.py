# Description: Main script to run the glycan analysis.

# Imports
from glycan_analysis import *

if __name__ == "__main__":
    # Input files
    itp_file = "../data/IL3_GLY_R1/toppar/PROA.itp"
    gro_file = "../data/IL3_GLY_R1/IL3_GLY_R1_conv.gro"
    xtc_file = "../data/IL3_GLY_R1/IL3_GLY_R1_conv.xtc"

    # Set of standard amino acids
    AMINO_ACIDS = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
        "GLY", "HIS", "HSD", "HSP" "ILE", "LEU", "LYS",
        "MET", "PHE", "PRO", "SER", "THR",
        "TRP", "TYR", "VAL"
    }

    ##############
    # A) PARSING #
    ##############
    # Parse the .itp file
    glycans, num_glycans, atom_data = parse_itp(itp_file, AMINO_ACIDS)

    # Write glycan ranges to file
    output_file = "../results/IL3_GLY_R1/glycan_ranges.txt"
    write_glycan_ranges(output_file, glycans, num_glycans, AMINO_ACIDS)
    print(f"Glycan ranges and dihedrals written to {output_file}\n")

    ########################
    # B) BUILD ATOM SETS   #
    ########################
    # For grouping dihedrals by glycan. These are 1-based indices (matching .itp).
    glycan_atom_sets = []
    for glycan in glycans:
        atom_indices = set(a["atom_index"] for a in glycan["atoms"])
        glycan_atom_sets.append(atom_indices)

    # Collect all glycan-only and protein-glycan dihedrals
    protein_glycan_dihedrals = []
    glycan_glycan_dihedrals = []
    for g in glycans:
        protein_glycan_dihedrals.extend(g["protein_glycan_dihedrals"])
        glycan_glycan_dihedrals.extend(g["glycan_dihedrals"])
    # Remove any duplicates
    protein_glycan_dihedrals = list(set(protein_glycan_dihedrals))
    glycan_glycan_dihedrals = list(set(glycan_glycan_dihedrals))

    #########################
    # C) MDANALYSIS UNIVERSE#
    #########################
    print(f"Loading the trajectory from: {gro_file} and {xtc_file}")
    u = mda.Universe(gro_file, xtc_file)

    # We want angles for all these dihedrals
    dihedrals_of_interest = list(set(protein_glycan_dihedrals + glycan_glycan_dihedrals))

    ###################################################
    # D) COMPUTE DIHEDRALS ACROSS FRAMES WITH MDAnalysis
    ###################################################
    print("Computing dihedral angles over the trajectory...")
    times, angle_data = compute_dihedral_time_series_vectorized(
        u,
        dihedrals_of_interest,
        start=None,  # e.g., 0 to skip initial frames
        stop=None,   # e.g., 1000 to limit frames
        step=None    # e.g., 10 to skip frames
    )
    print("Done computing dihedral angles.\n")

    ###################################
    # E) USER PROMPT TO PLOT OR NOT   #
    ###################################
    want_plot = None

    while want_plot is None:
        try:
            # Try showing a Tkinter pop-up
            root = tk.Tk()
            root.withdraw()  # Hide the main Tk window
            answer = simpledialog.askstring(
                "Plot generation",
                "Do you want to generate dihedral angle plots ? (yes/no)",
            )
            root.destroy()  # Close the hidden root window

            if answer is None:
                # User closed pop-up or cancelled
                print("No input provided via pop-up. Will ask in terminal.")
                # Fallback to terminal approach below
                raise Exception("Pop-up closed.")
            else:
                response = answer.strip().lower()
                if response in ["yes", "y"]:
                    want_plot = True
                elif response in ["no", "n"]:
                    want_plot = False
                else:
                    print("Please enter 'yes' or 'no' in the pop-up.")
                    # loop again (won't set want_plot, so it stays None)
        except:
            # If pop-up fails or user closed it, ask in terminal
            while True:
                answer = input("Do you want to generate dihedral angle plots ? (yes/no): ")
                response = answer.strip().lower()
                if response in ["yes", "y"]:
                    want_plot = True
                    break
                elif response in ["no", "n"]:
                    want_plot = False
                    break
                else:
                    print("Please enter 'yes' or 'no'.")

    # ============================================
    # F) PLOT IF USER WANTS
    # ============================================
    if want_plot:
        plot_dihedrals_by_glycan_and_type_plotly(
            times,
            angle_data,
            protein_glycan_dihedrals,
            glycan_glycan_dihedrals,
            glycans,
            glycan_atom_sets
        )
        print("Plots have been generated.")
    else:
        print("Skipping plot generation.")

    print("Finished plotting. Script complete.")