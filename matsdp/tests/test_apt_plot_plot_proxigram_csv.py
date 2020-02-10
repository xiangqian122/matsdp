# -*- coding: utf-8 -*-
import os
import sys
package_path = '/user/specified/path/to/matsdp/'
sys.path.insert(0, os.path.abspath(package_path))

def test_plot_proxigram_csv():
    from matsdp.apt import apt_plot

    retn_val = apt_plot.plot_proxigram_csv(
        proxigram_csv_file_path = './apt/profile-interface0.csv',
        sysname = 'M2',
        visible_elmt_list = ['Ni','Al'],
        interplation_on = False,
        fig_width = 6,
        fig_height = 5,
        fig_dpi = 600,
        fig_format = 'png',
        )
    assert retn_val == 0

