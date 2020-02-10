# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_estruct():
    from matsdp.vasp import vasp_analyze
    from matsdp.vasp import vasp_plot
    retn_val1 = vasp_analyze.estruct(
        doscar_file_path = './vasp/DOSCAR',
        sysname = 'DOS1',
        )
    retn_val2 = vasp_plot.plot_poscar(
        poscar_file_path = './vasp/POSCAR_estruct_DOS1_Ef-7.0888.vasp',
        euler_angle_type = 'zyz',
        phi = -3,
        theta = 5,
        psi = 0,
        elmt_color = {'Ni':'red','Re':'blue'},
        draw_mirror_atom = True,
        box_on = True,
        axis_indicator =True,
        plot_cell_basis_vector_label = True,
        plot_atom_label = True,
        fig_format = 'png',
        fig_dpi = 100,
        draw_colormap = True,
        colormap_column_indx = 1,
        colormap_vmin = None,
        colormap_vmax = None,
        vmin_color = 'blue',
        vmax_color = 'red',
        colorbar_alignment = 'vertical'
        )
    retn_val3 = vasp_plot.plot_poscar(
        poscar_file_path = './vasp/POSCAR_estruct_DOS1_Ef-7.0888.vasp',
        euler_angle_type = 'zyz',
        phi = -3,
        theta = 5,
        psi = 0,
        elmt_color = {'Ni':'red','Re':'blue'},
        draw_mirror_atom = True,
        box_on = True,
        axis_indicator =True,
        plot_cell_basis_vector_label = True,
        plot_atom_label = True,
        fig_format = 'png',
        fig_dpi = 100,
        draw_colormap = True,
        colormap_column_indx = 1,
        colormap_vmin = -80,
        colormap_vmax = -40,
        vmin_color = 'blue',
        vmax_color = 'red',
        colorbar_alignment = 'vertical'
        )    
    assert retn_val1 == 0 and retn_val2 == 0 and retn_val3 ==0
