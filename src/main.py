# Description: Main script to run the glycan analysis.

# Imports

from dihedral_retrieval import (
    parse_atoms,
    atom_ranges_of_molecules,
    parse_dihedrals,
    filter_dihedrals,
    parse_gro_file,
    plot_all_dihedrals,
    regroup_glycan_data,
)

from pca import (
    compute_and_plot_pca,
)

from ramachandran import (
    compute_and_plot_ramachandran,
)

def main():
    """
    Main function to run the glycan analysis.
    
    This script performs the following steps:
    1. Parse .itp file to get molecules and atom ranges.
    2. Collect dihedrals from the .itp file.
    3. Filter and classify dihedrals.
    4. Parse .gro file to get dihedral coordinates.
    5. Regroup glycan data.
    6. Calculate dihedrals and plot them.
    7. Compute and plot PCA.
    8. Plot Ramachandran plots.
    """

    # Input files
    itp_file = "data/IL3_GLY_R1/toppar/PROA.itp"
    gro_file = "data/IL3_GLY_R1/IL3_GLY_R1_conv.gro"
    xtc_file = "data/IL3_GLY_R1/IL3_GLY_R1_conv.xtc"


    # ========================================================================
    # 1) PARSE .itp FILE AND GET MOLECULES + ATOM RANGES
    # ========================================================================

    # Parse the .itp file in the atoms section
    # This will return the proteins and glycans as a list of dictionaries
    # e.g. proteins = [ { Protein 1 : [ {nr, resnr, residu, atom}, ... ] }, { Protein 2 : [ {nr, resnr, residu, atom}, ... ] }, ... ]
    # e.g. glycans = [ { Glycan 1 : [ {nr, resnr, residu, atom}, ... ] }, { Glycan 2 : [ {nr, resnr, residu, atom}, ... ] }, ... ]
    proteins, glycans = parse_atoms(itp_file)
    print("\nProteins in the system:", len(proteins))
    print("Glycans in the system:", len(glycans), "\n")

    # Get the atom ranges for proteins and glycans
    # e.g. protein_ranges = [ { Protein 1 : [nr_start, nr_end] }, { Protein 2 : [nr_start, nr_end] }, ... ]
    # e.g. glycan_ranges = [ { Glycan 1 : [nr_start, nr_end] }, { Glycan 2 : [nr_start, nr_end] }, ... ]
    protein_ranges, glycan_ranges = atom_ranges_of_molecules(proteins, glycans)
    print("\nProtein ranges:", protein_ranges)
    print("Glycan ranges:", glycan_ranges)


    # ========================================================================
    # 2) PARSE .itp FILE AND COLLECT DIHEDRALS
    # ========================================================================

    # Collect all glycan-only and protein-glycan dihedrals
    # e.g. glycan_glycan_dihedrals = [ { Glycan 1 : [ [ai, aj, ak, al], ... ] }, { Glycan 2 : [ [ai, aj, ak, al], ... ] }, ... ]
    # e.g. protein_glycan_dihedrals = [ { Protein 1 : [ [ai, aj, ak, al], ... ] }, { Protein 2 : [ [ai, aj, ak, al], ... ] }, ... ]
    gg_dihedrals, pg_dihedrals = parse_dihedrals(itp_file, glycan_ranges)

    # Count the total number of protein-glycan dihedrals
    total_pg_dihedrals = 0
    for glycan_dict in pg_dihedrals:
        for glycan_name, dihedrals in glycan_dict.items():
            total_pg_dihedrals += len(dihedrals)
    print("\nTotal protein-glycan dihedrals:", total_pg_dihedrals)

    # Count the total number of glycan-glycan dihedrals
    total_gg_dihedrals = 0
    for glycan_dict in gg_dihedrals:
        for glycan_name, dihedrals in glycan_dict.items():
            total_gg_dihedrals += len(dihedrals)
    print("Total glycan-glycan dihedrals:", total_gg_dihedrals)


    # ========================================================================
    # 3) FILTER AND CLASSIFY DIHEDRALS
    # ========================================================================

    # Here we filter out the dihedrals that are not considered glycosidic linkages
    # We have a set of criteria to fill out:
    # In general we are going to want to keep all the dihedrals that involve 2 residues
    # In general we want to exclude all the dihedrals that have more than 1 hydrogen H
    # phi angle:
    # - In the case of a protein-glycan phi: O5 - C1 - N'D2 - C'G
    # - In the case of a glycan-glycan phi: O5 - C1 - Ox - C'x
    # psi angle:
        # - In the case of a protein-glycan psi: C1 N'D2 - C'G - C'B
        # - In the case of a glycan-glycan psi: C1 - Ox - C'x - C'x
    # omega angle: O6 - C6 - C5 - O5

    # e.g.filtered_glycan_glycan_dihedrals = [ { Glycan 1 : [ { 'phi' : [ [ai, aj, ak, al], ... }, {'psi' : [ [ai, aj, ak, al], ... ] }, {'omega' : [ [ai, aj, ak, al], ...] } ] }, { Glycan 2 : ... }, ... ]
    # e.g. filtered_protein_glycan_dihedrals = [ { Protein 1 : [ { 'phi' : [ [ai, aj, ak, al], ... ] }, {'psi' : [ [ai, aj, ak, al], ...] }, {'omega' : [ [ai, aj, ak, al], ... ] } ] }, { Protein 2 : ... }, ... ]
    filtered_pg, filtered_gg = filter_dihedrals(pg_dihedrals, gg_dihedrals, proteins, glycans)
    print("\nFiltered Protein-Glycan Dihedrals:", filtered_pg)
    print("Filtered Glycan-Glycan Dihedrals:", filtered_gg)

    # ========================================================================
    # 4) PARSE .gro FILE AND GET DIHEDRAL COORDINATES
    # ========================================================================

    # Next we parse the .gro to go get the coordinates of the filtered dihedrals
    # e.g. glycosidic_coordinates_prot_glycan = [ { Protein 1 : [ { 'phi' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ... ] }, {'psi' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ...] }, {'omega' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ... ] } ] }, { Protein 2 : ... }, ... ]
    # e.g. glycosidic_coordinates_glycan_glycan = [ { Glycan 1 : [ { 'phi' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ... ] }, {'psi' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ...] }, {'omega' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ... ] } ] }, { Glycan 2 : ... }, ... ]
    prot_coords, glyco_coords = parse_gro_file(gro_file, filtered_gg, filtered_pg)
    print("\nProtein-Glycan Dihedral Coordinates:", prot_coords)


    # ========================================================================
    # 5) REGROUP GLYCAN DATA
    # ========================================================================

    # Here we regroup the glycan data to get the dihedrals for each glycan
    # e.g. system_glycans_and_dihedrals = { 'Glycan 1' : { 'phi' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ... ], 'psi' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ... ], 'omega' : [ [(x,y,z), (x,y,z), (x,y,z), (x,y,z)], ... ] }, 'Glycan 2' : [...], ... }
    system_glycans_and_dihedrals = regroup_glycan_data(filtered_pg, filtered_gg)
    print("\nGlycan Data:", system_glycans_and_dihedrals)

    # ========================================================================
    # 6) CALCULATE DIHEDRALS AND PLOT
    # ========================================================================
    
    # Prompt to ask the user if they want to plot the dihedral data
    user_choice = input("\nDo you want to plot the dihedral data? (y/n): ").strip().lower()
    if user_choice == "y":
        # Here we plot the dihedrals over the course of the trajectory
        print("Plotting dihedrals...")
        plot_all_dihedrals(xtc_file, gro_file, filtered_pg, filtered_gg)
    else:
        print("Plotting skipped.")


    # ========================================================================
    # COMPUTE AND PLOT PCAS
    # ========================================================================

    # Prompt to ask the user if they want to compute and plot the PCA plots
    user_choice = input("\nDo you want to plot the PCA? (y/n): ").strip().lower()
    if user_choice == "y":
        # Here we compute and plot the PCA for the glycan dihedral data
        print("Computing and plotting PCA...")
        compute_and_plot_pca(system_glycans_and_dihedrals, gro_file, xtc_file)
    else:
        print("PCA skipped.")
    

    # ========================================================================
    # PLOT RAMACHANDRAN PLOTS
    # ========================================================================


    # Prompt to ask the user if they want to compute and plot the Ramachandran plots
    user_choice = input("\nDo you want to plot the Ramachandran plots? (y/n): ").strip().lower()
    if user_choice == "y":
        # Here we compute and plot the PCA for the glycan dihedral data
        print("Computing and plotting Ramachandran plots...")
        compute_and_plot_ramachandran(system_glycans_and_dihedrals, gro_file, xtc_file)
    else:
        print("Ramachandran skipped.")


if __name__ == "__main__":
    main()
