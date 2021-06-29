# -*- coding: utf-8 -*-
def read_band(band_file_path):
    '''
    read the band data
    '''
    args_dict = locals()
    import os
    import numpy as np
    from .. import funcs
    import math

    band_dict = {}
    band_dict['eigs'] = None
    band_dict['num_kpoints'] = None
    band_dict['num_bands'] = None
    band_dict['kpath_len_list'] = None

    band_file_path = os.path.abspath(band_file_path)
    band_dict['num_bands'] = 0
    band_dict['num_kpoints'] = 0
    band_dict['kpath_len_list'] = []
    band_line_indx_list = []
    last_line_content = None
    with open(band_file_path, 'r') as f:
        lines = f.readlines()
        # Get number of kpoints and kpath length list
        data_found = False
        for i_line in lines:
            if i_line.strip().strip('\n') != '' and not i_line.startswith('#'):
                data_found = True
                band_dict['num_kpoints'] += 1
                band_dict['kpath_len_list'].append(float(funcs.split_line(i_line, separator = ' ')[0]))
            elif (i_line.strip().strip('\n') == '' or i_line.startswith('#')) and data_found == True:
                break
        # Get estimated number of bands
        estimated_num_bands = 0
        if (len(lines) % band_dict['num_kpoints']) < band_dict['num_kpoints']:
            estimated_num_bands = math.ceil(len(lines) / band_dict['num_kpoints'])
        else:
            estimated_num_bands = math.ceil(len(lines) / band_dict['num_kpoints']) + 1
        eigs_arr = np.array([None] * band_dict['num_kpoints'] * estimated_num_bands)
        eigs_arr.shape = (band_dict['num_kpoints'], estimated_num_bands) 

        first_line_data_found = lines[0].strip().strip('\n') == '' or i_line.startswith('#')
        if first_line_data_found:
            band_dict['num_bands'] += 1
            eigs_arr[0, 0] = float(funcs.split_line(i_line, separator = ' ')[1])
        # Get Eignevalue array
        i_kpt = 0
        i_band = 0
        for i_line_indx in range(1, len(lines)):
            data_found = (lines[i_line_indx].strip().strip('\n') != '' and not lines[i_line_indx].startswith('#'))
            last_line_data_found = (lines[i_line_indx - 1].strip().strip('\n') != '' and not lines[i_line_indx - 1].startswith('#'))
            if data_found and not last_line_data_found:
                band_dict['num_bands'] += 1
                i_kpt = 0
                i_band = band_dict['num_bands'] - 1
                eigs_arr[i_kpt, i_band] = float(funcs.split_line(lines[i_line_indx], separator = ' ')[1])
            elif data_found and last_line_data_found:
                i_kpt += 1
                i_band = band_dict['num_bands'] - 1
            
            #print(i_line_indx, i_kpt, i_band, lines[i_line_indx])
                eigs_arr[i_kpt, i_band] = float(funcs.split_line(lines[i_line_indx], separator = ' ')[1])
        band_dict['eigs'] = eigs_arr[:, 0: band_dict['num_bands']]
    return band_dict
