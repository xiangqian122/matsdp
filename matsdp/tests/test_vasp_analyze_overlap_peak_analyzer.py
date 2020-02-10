# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_overlap_peak_analyzer():
    from matsdp.vasp import vasp_analyze
    retn_val = vasp_analyze.overlap_peak_analyzer(
        doscar_file_path = './vasp/DOSCAR',
        sysname = 'DOS1',
        atom_indx_list = ['Ni1','Re1'],
        n_shell = 2,
        a0 = 3.52,
        dos_mode = {'Ni':['d'],'Re':['d']},
        fermi_shift_zero = True,
        )
    assert retn_val == 0
