NOTES ON QB EDITS TO DEFAULT CCG FILES

Generally, we do not like changing MOE SVL files. All QB-changes to CCG-provided files are 
submitted to CCG.

*ALL PATCHES ARE APPLIED TO MOST CURRENT PUBLIC VERSION OF MOE*

 * io_pdb.svl : The following changes were made:
                * We ran into problems with the PDB file reader/writer not properly
                maintaining symmetry for crystallography users. If the symmetry is in
                the PDB file then it should at least be read into MOE. Likewise, it should
                be preserved and written out to the file. Having to choose PBC or some
                other symmetry operation in order to preserve this setting was unintuitive
                to X-ray crystallographers.
                * When setting the valance, alternate atoms were not preserved. The
                getHeavyValence function is provided for this purpose.
                * We added an option to apply symmetry to H2O.
                
 * protonate3d.svl : There are times when a truncated residue will cause problems with
                protonate3D. When the residue is truncated down to a single atom, 
                protonate3D will attempt to "merge" the atom with a previous residue.
                For example, see LYS 184 Chain A in 2BOV. LYS184 has a single
                N atom. When you run Protonate3D, the N atom is flipped with the O atom of
                the SER183. This causes there to be 2 atoms named N within the SER183.

 * structprep_ui.svl : Over the years, numerous changes have been made/suggested for
                structprep_ui.svl and most of them have made their way back to CCG. The
                two key remaining changes include:
                * There is a bug in the TRS residue provided in the system MDB file used 
                in structure prep. The proton names are incorrect which phenix does not like.
                * The default is to cap a broken terminal residue which is not what
                crystallographers like to see happen. Unless a terminal residue is an
                actual terminal, it should simple be charged.

 * save.svl : When the cryst parameters are available, they should be saved.