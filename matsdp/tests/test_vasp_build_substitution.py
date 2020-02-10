# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_substitution():
    from matsdp.vasp import vasp_build

    retn_val = vasp_build.substitution(
        substitution_list_file = './vasp/example.subst',
        poscar_file_path = './vasp/POSCAR_NoDope',
        )
    assert retn_val == 0
