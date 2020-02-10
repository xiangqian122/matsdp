# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_selection_sphere():
    from matsdp.vasp import vasp_build
    retn_val = vasp_build.selection_sphere(
        poscar_file_path = './vasp/CONTCAR',
        origin_atom_name = 'Re1',
        radius = 7,
        include_mirror_atoms = False,
        output_file_name = 'example'
        )
    assert retn_val == 0

