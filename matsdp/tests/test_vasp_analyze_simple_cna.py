# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_simple_cna():
    from matsdp.vasp import vasp_analyze
    from matsdp.vasp import vasp_plot

    retn_val1 = vasp_analyze.simple_cna(
        poscar_file_path = './vasp/POSCAR',
        a0 = 3.545,
        common_neighbor_elmt_list = ['Re', 'W', 'Ta','Ni']
        )
    retn_val2 = vasp_plot.plot_poscar(
        poscar_file_path = './vasp/POSCAR_simple_common_neighbor_pair_count_ReNi.vasp',
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
        colormap_column_indx = 2,
        colormap_vmin = None,
        colormap_vmax = None,
        vmin_color = 'blue',
        vmax_color = 'red',
        colorbar_alignment = 'vertical'
        )
    retn_val3 = vasp_plot.plot_poscar(
        poscar_file_path = './vasp/POSCAR_simple_common_neighbor_pair_count_ReNi.vasp',
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
        colormap_column_indx = 2,
        colormap_vmin = None,
        colormap_vmax = None,
        vmin_color = 'blue',
        vmax_color = 'red',
        colorbar_alignment = 'horizontal'
        )

    assert isinstance(retn_val1,dict) == True and retn_val2 == 0 and retn_val3 == 0
