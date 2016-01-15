# Protein-peptide-docking_1CZY

## About Protein-peptide-docking_1CZY
This repository contains all necessary information and scripts to do protein-peptide docking in HADDOCK, using the complex 1CZY as an example, which is the Supplementary Information for paper[reference]. For the usage of these files and scripts, please check the content of the paper[reference]

*Note that the provided scripts are in principle generic and can be used for other systems as well.*

This repository contains following folders and files:

# pdb_files

  - reference_complex.pdb, the crystal structure of the protein-peptide complex used as reference, which is chain A and D of 1CZY.
  - protein.pdb, the crystal structure of unbound protein used for docking, which is chain C of 1CZZ.
  - peptide_alpha.pdb, the alpha-helix conformation of peptide.
  - peptide_polypro.pdb, the polyproline-II conformation of peptide.
  - peptide_extended.pdb, the extended conformation of peptide.
  - protein_active_residues.dat, the list of active residues of protein.
  - peptide_passive_residues.dat, the list of passive residues of peptide.
  - protein_histidine_states.dat, the list of charged states of histidines in protein.
  - structures.list, the list of 33 conformations or structures of peptide used for docking.

  # MD_conformations
    - md_1.pdb ... md_30.pdb, the 30 structures of peptide from MD simulations.


# MD-Parameter_Files

  - vacuum.mdp, MD parameters for energy minimization in vacuum.
  - ions.mdp, MD parameters for energy minimization in the presence of solvent and ions.
  - nvt.mdp, MD parameters for energy minimization in constant volume conditions.
  - npt.mdp, MD parameters for energy minimization in constant pressure conditions.
  - unrestrained.mdp, MD parameters for energy minimization without any restraints.
  - production.mdp, MD parameters for final production MD simulation.


# docking
  - restraints.tbl, the AIRs file.
  - new.html, the file to start a new HADDOCK project.


# scripts-PyMOL
  -build_seq.py, PyMOL script to build protein structure from sequence.


# scripts-MD
  - automd.sh, bash script to prepare files and perform MD simulations using GROMACS.
  - dpca.sh, bash script to do dihedral PCA using GROMACS.


# scripts-HADDOCK
  - molprobity.py, python script to define charge states of histidines.


## How to Cite
C. Geng, S. Narasimhan, João P.G.L.M. Rodrigues and Alexandre M.J.J. Bonvin, Information-driven, ensemble flexible peptide docking using HADDOCK, Modeling of peptide-protein interactions, Methods in Molecular Biology series (Springer Press), 2016

