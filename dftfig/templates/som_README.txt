Site-of-metabolism weights for the mpo_metabolism figure, panels B-E
====================================================

Only needed when figures.som_source is 'csv' or 'both'. With the default
'fukui_radical' the pipeline computes the weights itself from the condensed
radical Fukui index and this directory can stay empty.

File name   <molecule-name>.csv   -- must match the key in molecules.definitions
Columns     atom_index,probability

  atom_index   1-based index into the HEAVY ATOMS of the molecule, in the order
               RDKit reads the SMILES in molecules.definitions. Hydrogens are
               excluded. Check the ordering before trusting a figure:

                   from rdkit import Chem
                   m = Chem.MolFromSmiles("COc1ccc2[nH]c3ccccc3c2c1Cl")
                   for a in m.GetAtoms():
                       print(a.GetIdx() + 1, a.GetSymbol())

  probability  0-1. Rows outside the heavy-atom range are ignored with a warning.

Example (CCB1.csv):

    atom_index,probability
    1,0.42
    2,0.08
    13,0.31

Exporting from XenoSite: download the per-atom table and renumber it to the
SMILES heavy-atom order above - XenoSite's own numbering follows the structure
file you uploaded, which need not match.
