# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_read_doscar():
    from matsdp.vasp import vasp_read

    poscar_file_path = './vasp/POSCAR'
    poscar_dict = vasp_read.read_poscar(poscar_file_path)  
    for atom_indx in range(0, len(poscar_dict['atom_species_arr']) + 1):
        vasp_read.read_doscar(
            doscar_file_path = './vasp/DOSCAR',
            atom_indx = atom_indx,
            save_dos_arr = True,
            )
