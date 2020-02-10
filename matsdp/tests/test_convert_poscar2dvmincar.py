# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_poscar2dvmincar():
    from matsdp import convert

    retn_val = convert.poscar2dvmincar(
        poscar_path = './vasp/CONTCAR'
        )
    assert retn_val == 0
