def read_ind(ind_file_path):
    '''
    read the *.ind file of the DVM program
    ind_file_path: the file path of the IND.DAT file
    '''
    import os
    import numpy as np
    from .. import funcs
    ind_file_path = os.path.abspath(ind_file_path)
    dvm_ind_dict = {}
    temp_list = funcs.split_line(line = funcs.grep(kwd = 'ATOM',file_dir = ind_file_path), separator = '\n')
    num_elmts = len(temp_list)
    dvm_ind_dict['elmt_species_arr'] = np.array([None] * num_elmts)
    for i_line in range(0, num_elmts):
        elmt_name = funcs.split_line(temp_list[i_line])[2]
        dvm_ind_dict['elmt_species_arr'][i_line] = elmt_name
    return dvm_ind_dict

def read_incar(incar_file_path):
    '''
    read the *.incar file of the DVM program
    incar_file_path: the file path of the *.incar file
    '''
    import os
    import numpy as np
    from .. import funcs
    incar_file_path = os.path.abspath(incar_file_path)
    # designate the working directory
    workdir, incar_file = funcs.file_path_name(incar_file_path)

    dvm_incar_dict = {}
    with open(incar_file_path, 'r') as f_incar:
        lines = f_incar.readlines()
        dvm_incar_dict['n_atoms'] = int(funcs.split_line(lines[0])[0])
        dvm_incar_dict['r_buff'] = float(funcs.split_line(lines[0])[1])
        dvm_incar_dict['pos_arr'] = np.array([None] * dvm_incar_dict['n_atoms'] * 3)
        dvm_incar_dict['pos_arr'].shape = dvm_incar_dict['n_atoms'], 3
        dvm_incar_dict['dvm_atom_elmt_indx_arr'] = np.array([None] * dvm_incar_dict['n_atoms'])
        dvm_incar_dict['dvm_elmt_indx_arr'] = np.array([])
        vacant_line = False
        keep_record = False
        atom_indx = 0
        dvm_incar_dict['atom_indx_arr'] = np.array([None] * dvm_incar_dict['n_atoms'])
        for i_line in range(0, len(lines)):
            if len(funcs.split_line(lines[i_line])) == 0:
                vacant_line = True
                keep_record = True
            else:
                vacant_line = False
            if vacant_line == False and keep_record == True:
                dvm_incar_dict['pos_start_line_indx'] = i_line
                dvm_incar_dict['pos_arr'][atom_indx, 0] = float(funcs.split_line(lines[i_line])[0])
                dvm_incar_dict['pos_arr'][atom_indx, 1] = float(funcs.split_line(lines[i_line])[1])
                dvm_incar_dict['pos_arr'][atom_indx, 2] = float(funcs.split_line(lines[i_line])[2])
                dvm_incar_dict['dvm_atom_elmt_indx_arr'][atom_indx] = int(funcs.split_line(lines[i_line])[3]) ## index starts from 1
                dvm_incar_dict['atom_indx_arr'][atom_indx] = atom_indx + 1   ## atom index starts from 1
                atom_indx += 1
        dvm_ind_dict = read_ind(workdir + '/IND.DAT')
        dvm_incar_dict['num_elmts'] = len(dvm_ind_dict['elmt_species_arr'])
        dvm_incar_dict['elmt_species_arr'] = dvm_ind_dict['elmt_species_arr']
        dvm_incar_dict['dvm_elmt_indx_arr'] = set(list(dvm_incar_dict['dvm_atom_elmt_indx_arr']))
        dvm_incar_dict['elmt_num_arr'] = [list(dvm_incar_dict['dvm_atom_elmt_indx_arr']).count(x) for x in set(list(dvm_incar_dict['dvm_atom_elmt_indx_arr']))]
        dvm_incar_dict['xmin'] = np.min(dvm_incar_dict['pos_arr'][:, 0])
        dvm_incar_dict['xmax'] = np.max(dvm_incar_dict['pos_arr'][:, 0])
        dvm_incar_dict['ymin'] = np.min(dvm_incar_dict['pos_arr'][:, 1])
        dvm_incar_dict['ymax'] = np.max(dvm_incar_dict['pos_arr'][:, 1])
        dvm_incar_dict['zmin'] = np.min(dvm_incar_dict['pos_arr'][:, 2])
        dvm_incar_dict['zmax'] = np.max(dvm_incar_dict['pos_arr'][:, 2])
    return dvm_incar_dict

def read_otput(otput_file_path):
    '''
    read the *.otput file of the DVM program
    otput_file_path: the file path of the *.otput file
    '''
    import os
    import numpy as np
    from .. import funcs
    otput_file_path = os.path.abspath(otput_file_path)
    # designate the working directory
    workdir, otput_file = funcs.file_path_name(otput_file_path)

    dvm_otput_dict = {}
    dvm_otput_dict['unitscale'] = float(funcs.split_line(funcs.grep('unitscale', otput_file_path))[2])
    dvm_otput_dict['nbasis'] = int(funcs.split_line(funcs.grep('nbasis', otput_file_path))[2])
    dvm_otput_dict['natom'] = int(funcs.split_line(funcs.grep('natom', otput_file_path))[2])
    dvm_otput_dict['spin'] = int(funcs.split_line(funcs.grep('spin', otput_file_path))[2])
    dvm_otput_dict['ie_arr'] = np.array([None] * dvm_otput_dict['natom'] * dvm_otput_dict['natom'])
    dvm_otput_dict['ie_arr'].shape = dvm_otput_dict['natom'], dvm_otput_dict['natom']
    dvm_otput_dict['ie_up_arr'] = dvm_otput_dict['ie_arr'].copy()
    dvm_otput_dict['ie_dw_arr'] = dvm_otput_dict['ie_arr'].copy()
    kwd_ie = 'The interatomic energy'
    kwd_line_indx = int(funcs.find_kwd_line_indx(kwd_ie, otput_file_path)[0])
    with open(otput_file_path, 'r') as f_otput:
        lines = f_otput.readlines()

        vacant_line = True
        last_state = True
        transition_times = 0  # count the number of transitions
        for i_line in range(kwd_line_indx + 1, len(lines)):
            # everytime the last_state changed, the transition_times adds one
            # this is useful for identifying finite number of string blocks
            if len(funcs.split_line(lines[i_line])) == 0:
                vacant_line = True
                if vacant_line == last_state:
                    pass
                else:
                    transition_times += 1
                last_state = True
            else:
                vacant_line = False
                if vacant_line == last_state:
                    pass
                else:
                    transition_times += 1
                last_state = False
            if transition_times == 1:
                # read the lines
                i_temp_indx = 0
                j_temp_indx = 1
                ie_au_temp_indx = 2
                line_list = funcs.split_line(lines[i_line].replace('[',' ').replace(']',' '))
                num_col_block = int(len(line_list)/3)
                for i_col_block in range(0,num_col_block):
                    i_dvm_atom_indx = int(line_list[i_temp_indx])
                    j_dvm_atom_indx = int(line_list[j_temp_indx])
                    ie_au = float(line_list[ie_au_temp_indx])
                    i_temp_indx += 3
                    j_temp_indx += 3
                    ie_au_temp_indx += 3
                    dvm_otput_dict['ie_arr'][i_dvm_atom_indx - 1, j_dvm_atom_indx - 1] = ie_au
                    dvm_otput_dict['ie_up_arr'][i_dvm_atom_indx - 1, j_dvm_atom_indx - 1] = ie_au
            if dvm_otput_dict['spin'] == 1 and transition_times == 3:
                # read the lines
                i_temp_indx = 0
                j_temp_indx = 1
                ie_au_temp_indx = 2
                line_list = funcs.split_line(lines[i_line].replace('[',' ').replace(']',' '))
                num_col_block = int(len(line_list)/3)
                for i_col_block in range(0,num_col_block):
                    i_dvm_atom_indx = int(line_list[i_temp_indx])
                    j_dvm_atom_indx = int(line_list[j_temp_indx])
                    ie_au = float(line_list[ie_au_temp_indx])
                    i_temp_indx += 3
                    j_temp_indx += 3
                    ie_au_temp_indx += 3
                    dvm_otput_dict['ie_dw_arr'][i_dvm_atom_indx - 1, j_dvm_atom_indx - 1] = ie_au
    return dvm_otput_dict

def read_input(input_file_path):
    '''
    read the *.incar file of the DVM program
    incar_file_path: the file path of the *.incar file
    '''
    import os
    import numpy as np
    from .. import funcs
    from . import dvm_default
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)

    input_file_path = os.path.abspath(input_file_path)
    # designate the working directory
    workdir, dvm_otput_file = funcs.file_path_name(input_file_path)

    dvm_input_dict = dvm_default.input_default()
    with open(input_file_path, 'r') as f:
        lines = f.readlines()
        for i_line in range(0, len(lines)):
            if len(funcs.split_line(lines[i_line])) == 0:
                continue
            else:
                temp_str = funcs.split_line(lines[i_line])[0]
                sepcial_input_key_list = [
                    'jobid',
                    'efld',
                    'grid3dseed',
                    'frozeninfo',
                    'grid3dinfo',
                    'plotgrid'
                    ]
                if temp_str[0] == '#':
                    continue
                else:
                    if temp_str in dvm_input_dict.keys() and temp_str not in sepcial_input_key_list:
                        dvm_input_dict[temp_str] = funcs.split_line(lines[i_line])[2]
                    if temp_str in ['jobid', 'efld', 'grid3dseed']:
                        dvm_input_dict[temp_str] = funcs.split_line(lines[i_line], '=')[1]
    return dvm_input_dict
