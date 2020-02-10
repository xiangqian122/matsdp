# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_create_multiple_dvm_jobs():
    from matsdp.dvm import dvm_build
    poscar_file_path_dict = {}
    origin_atom_name_list_dict = {}
    poscar_file_path_dict['dvm_example'] = './vasp/CONTCAR'
    origin_atom_name_list_dict['dvm_example'] = ['Re1', 'Ni5']
    retn_val = dvm_build.create_multiple_dvm_jobs(
        poscar_file_path_dict = poscar_file_path_dict,
        origin_atom_name_list_dict = origin_atom_name_list_dict,
        elmt_ind_file_dir = './dvm/ind_extended/',
        radius = 8,
        include_mirror_atoms = True
        )
    assert retn_val == 0

