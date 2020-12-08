# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_nn_map():
    from matsdp.vasp import vasp_analyze
    retn_val = vasp_analyze.nn_map(
        poscar_file_path = './vasp/POSCAR',
        a0 = 3.545,
        n_shell = 2,
        elmt_selector_list = ['Re'],
        atom_label = 'n_xyz',
        )
    assert isinstance(retn_val,dict) == True
