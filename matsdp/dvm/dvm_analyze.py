def ie_nn(dvm_otput_file_path, a0 = 3.545):
    '''
    Interatomic energy between the atoms and their nearest neighbor atoms.
    This module has been tested for the source_23oct05 version of the DVM program
    dvm_otput_file_path: the *.otput file of the DVM output
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import convert
    from ..vasp import vasp_read
    from ..vasp import vasp_analyze
    from . import dvm_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)

    dvm_otput_file_path = os.path.abspath(dvm_otput_file_path)
    # designate the working directory
    workdir, dvm_otput_file = funcs.file_path_name(dvm_otput_file_path)

    dvm_incar_file_path = workdir + '/' + dvm_otput_file[0:-6] + '.incar'
    convert.dvmincar2poscar(dvm_incar_file_path)
    dvmincar2poscar_poscar_file_path = workdir + '/' + dvm_otput_file[0:-6] + '_dvmincar2poscar.vasp'
    poscar_dict = vasp_read.read_poscar(dvmincar2poscar_poscar_file_path)
    atom_name_list = poscar_dict['atomname_list']
    atom_indx_arr = poscar_dict['atom_indx_arr'] 
    n_atoms = poscar_dict['n_atoms']
    added_atom_data_arr = poscar_dict['added_atom_data']
    dvm_otput_dict = dvm_read.read_otput(dvm_otput_file_path)

    if dvm_otput_dict['spin'] == 0:
        ie_au_arr = np.array([None] * n_atoms * n_atoms)
        ie_au_arr.shape = n_atoms, n_atoms
        for i_atom in range(0, n_atoms):
            for j_atom in range(0, n_atoms):
                ie_au_arr[i_atom, j_atom] = dvm_otput_dict['ie_arr'][int(added_atom_data_arr[i_atom, 3]) - 1, int(added_atom_data_arr[j_atom, 3]) - 1]
    if dvm_otput_dict['spin'] == 1:
        ie_au_up_arr = np.array([None] * n_atoms * n_atoms)
        ie_au_up_arr.shape = n_atoms, n_atoms
        ie_au_dw_arr = np.array([None] * n_atoms * n_atoms)
        ie_au_dw_arr.shape = n_atoms, n_atoms
        for i_atom in range(0, n_atoms):
            for j_atom in range(0, n_atoms):
                ie_au_up_arr[i_atom, j_atom] = dvm_otput_dict['ie_up_arr'][int(added_atom_data_arr[i_atom, 3]) - 1, int(added_atom_data_arr[j_atom, 3]) - 1]
                ie_au_dw_arr[i_atom, j_atom] = dvm_otput_dict['ie_dw_arr'][int(added_atom_data_arr[i_atom, 3]) - 1, int(added_atom_data_arr[j_atom, 3]) - 1]

#########################
# IE for 1NN results
#########################
    # rlated files
    nn_atomname_list_file_path = output_dir + '/nn_atomname_list_without_mirror_label_1NN.txt'
    ie_au_nn_file_path = workdir + '/' + dvm_otput_file[0:-6] + '_ie_au_nn.txt'
    ie_ev_nn_file_path = workdir + '/' + dvm_otput_file[0:-6] + '_ie_ev_nn.txt'
    # call nearest neighbor analysis
    vasp_analyze.nn_map(
    poscar_file_path = dvmincar2poscar_poscar_file_path,
    a0 = a0,
    n_shell = 1
    )
    # prepare the data
    first_nn_ie_au_arr = None
    first_nn_ie_au_up_arr = None
    first_nn_ie_au_dw_arr = None
    with open(nn_atomname_list_file_path, 'r') as f_1nn_atomname:
        lines = f_1nn_atomname.readlines()
        # initializetion
        max_num_first_nn = max([len(funcs.split_line(x, ',')[1:]) for x in lines])
        first_nn_atom_name_arr = np.array([None] * n_atoms * max_num_first_nn)
        first_nn_atom_name_arr.shape = n_atoms, max_num_first_nn
        first_nn_dvm_atom_indx_arr = np.array([None] * n_atoms * max_num_first_nn)
        first_nn_dvm_atom_indx_arr.shape = n_atoms, max_num_first_nn
        if dvm_otput_dict['spin'] == 0:
            first_nn_ie_au_arr = np.array([None] * n_atoms * max_num_first_nn)
            first_nn_ie_au_arr.shape = n_atoms, max_num_first_nn
        if dvm_otput_dict['spin'] == 1:
            first_nn_ie_au_up_arr = np.array([None] * n_atoms * max_num_first_nn)
            first_nn_ie_au_up_arr.shape = n_atoms, max_num_first_nn
            first_nn_ie_au_dw_arr = np.array([None] * n_atoms * max_num_first_nn)
            first_nn_ie_au_dw_arr.shape = n_atoms, max_num_first_nn
        # get values
        for i_atom in range(0, n_atoms):
            i_atom_name = atom_name_list[i_atom]
            i_dvm_atom_indx = added_atom_data_arr[i_atom, 3]
            ##print(i_dvm_atom_indx, i_atom_name)
            num_col = len(funcs.split_line(lines[i_atom], ',')[1:])
            first_nn_atom_name_list = funcs.split_line(lines[i_atom], ',')[1:]
            first_nn_atom_name_arr[i_atom, 0:num_col] = first_nn_atom_name_list[:]
            for j_atom in range(0, num_col):
                j_atom_name = first_nn_atom_name_list[j_atom]
                j_dvm_atom_indx = added_atom_data_arr[atom_name_list.index(j_atom_name), 3]
                ##print('  ', j_dvm_atom_indx, j_atom_name)
                first_nn_dvm_atom_indx_arr[i_atom, j_atom] = j_dvm_atom_indx
                if dvm_otput_dict['spin'] == 0:
                    first_nn_ie_au_arr[i_atom, j_atom] = ie_au_arr[i_atom, atom_indx_arr[atom_name_list.index(j_atom_name)] - 1]
                if dvm_otput_dict['spin'] == 1:
                    first_nn_ie_au_up_arr[i_atom, j_atom] = ie_au_up_arr[i_atom, atom_indx_arr[atom_name_list.index(j_atom_name)] - 1]
                    first_nn_ie_au_dw_arr[i_atom, j_atom] = ie_au_dw_arr[i_atom, atom_indx_arr[atom_name_list.index(j_atom_name)] - 1]
    # write the data into files
    with open(ie_au_nn_file_path, 'w') as f_ie_au_nn, open(ie_ev_nn_file_path, 'w') as f_ie_ev_nn:
        pass
    with open(ie_au_nn_file_path, 'a') as f_ie_au_nn, open(ie_ev_nn_file_path, 'a') as f_ie_ev_nn:
        for i_atom in range(0, n_atoms):
            i_atom_str = atom_name_list[i_atom] + ' (' + str(added_atom_data_arr[i_atom, 3]) + ')'
            temp_arr = first_nn_atom_name_arr[i_atom,:]
            temp_arr_remove_none = temp_arr[temp_arr != np.array(None)]
            first_nn_atom_name_str = ', '.join([str(x) + ' (' + str(added_atom_data_arr[atom_name_list.index(x), 3]) + ')' + ' ' * (13 - len(str(x) + ' (' + str(added_atom_data_arr[atom_name_list.index(x), 3]) + ')')) for x in temp_arr_remove_none]) + '\n'

            if dvm_otput_dict['spin'] == 0:
                temp_au_arr = first_nn_ie_au_arr[i_atom,:]
                temp_au_arr_remove_none = temp_au_arr[temp_au_arr != np.array(None)]
                temp_ev_arr_remove_none = temp_au_arr[temp_au_arr != np.array(None)] * convert.unitconvert('a.u.','eV')
                first_nn_ie_au_str = ', '.join([str(x) + ' ' * (13 - len(str(x))) for x in temp_au_arr_remove_none]) + '\n'
                first_nn_ie_ev_str = ', '.join(['{:.4f}'.format(x) + ' ' * (13 - len('{:.4f}'.format(x))) for x in temp_ev_arr_remove_none]) + '\n'
                f_ie_au_nn.write(
                    i_atom_str + ' ' * (13 - len(i_atom_str)) + ', ' + 
                    first_nn_atom_name_str +
                    ' IE (a.u.)' + ' ' * (13 - len(' IE (a.u.)')) + ', ' +
                    first_nn_ie_au_str
                    )
                f_ie_ev_nn.write(
                    i_atom_str + ' ' * (13 - len(i_atom_str)) + ', ' + 
                    first_nn_atom_name_str +
                    ' IE (eV)'  + ' ' * (13 - len(' IE (eV)')) + ', ' +
                    first_nn_ie_ev_str
                    )
            if dvm_otput_dict['spin'] == 1:
                temp_au_up_arr = first_nn_ie_au_up_arr[i_atom,:]
                temp_au_up_arr_remove_none = temp_au_up_arr[temp_au_up_arr != np.array(None)]
                temp_ev_up_arr_remove_none = temp_au_up_arr[temp_au_up_arr != np.array(None)] * convert.unitconvert('a.u.','eV')
                temp_au_dw_arr = first_nn_ie_au_dw_arr[i_atom,:]
                temp_au_dw_arr_remove_none = temp_au_dw_arr[temp_au_dw_arr != np.array(None)]
                temp_ev_dw_arr_remove_none = temp_au_dw_arr[temp_au_dw_arr != np.array(None)] * convert.unitconvert('a.u.','eV')
                first_nn_ie_au_up_str = ', '.join([str(x) + ' ' * (13 - len(str(x))) for x in temp_au_up_arr_remove_none]) + '\n'
                first_nn_ie_ev_up_str = ', '.join(['{:.4f}'.format(x) + ' ' * (13 - len('{:.4f}'.format(x))) for x in temp_ev_up_arr_remove_none]) + '\n'
                first_nn_ie_au_dw_str = ', '.join([str(x) + ' ' * (13 - len(str(x))) for x in temp_au_dw_arr_remove_none]) + '\n'
                first_nn_ie_ev_dw_str = ', '.join(['{:.4f}'.format(x) + ' ' * (13 - len('{:.4f}'.format(x))) for x in temp_ev_dw_arr_remove_none]) + '\n'
                f_ie_au_nn.write(
                    i_atom_str + ' ' * (13 - len(i_atom_str)) + ', ' + 
                    first_nn_atom_name_str +
                    ' IE (a.u.) up' + ' ' * (13 - len(' ' * len(' IE (a.u.) up'))) + ', ' +
                    first_nn_ie_au_up_str +
                    ' IE (a.u.) dw' + ' ' * (13 - len(' ' * len(' IE (a.u.) dw'))) + ', ' +
                    first_nn_ie_au_dw_str
                    )
                f_ie_ev_nn.write(
                    i_atom_str + ' ' * (13 - len(i_atom_str)) + ', ' + 
                    first_nn_atom_name_str +
                    ' IE (eV) up' + ' ' * (13 - len(' ' * len(' IE (eV) dw'))) + ', ' +
                    first_nn_ie_ev_up_str +
                    ' IE (eV) dw' + ' ' * (13 - len(' ' * len(' IE (eV) dw'))) + ', ' +
                    first_nn_ie_ev_dw_str
                    )
    funcs.write_log(
        logfile,
        'dvm_analyze.ie_nn(' + '\n' +
        '    dvm_otput_file_path = ' + 'r\'' + str(dvm_otput_file_path) + '\'' + ',\n' +
        '    a0 = ' + str(a0) + ')\n' +
        '###############################\n')
    return first_nn_ie_au_arr, first_nn_ie_au_up_arr, first_nn_ie_au_dw_arr

def job_status(job_parent_dir):
    '''
    Check the job status for multiple jobs
    job_parent_dir: This is the parent directory which contains the multiple DVM jobs.
    '''
    import os
    from .. import funcs
    from .. import default_params
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)

    def dvm_job_finished(otput_file_path):
        import os
        kwd = 'Job end'
        retn_val = False
        otput_file_path = os.path.abspath(otput_file_path)
        with open(otput_file_path, 'r') as f:
            lines = f.readlines()
            for i_line in lines:
                if kwd in i_line:
                    retn_val = True
        return retn_val

    job_parent_dir = os.path.abspath(job_parent_dir)
    job_dir_list = [job_parent_dir + '/' + x for x in os.listdir(job_parent_dir)]
    job_dir_name_list = os.listdir(job_parent_dir)
    max_length_job_dir_name = max([len(x) for x in job_dir_name_list])

    job_status_file_path = output_dir + '/dvm_job_status.txt'

    with open(job_status_file_path, 'w') as f:
        pass
    with open(job_status_file_path, 'a') as f:
        for i_job_indx in range(0, len(job_dir_list)):
            otput_exist = False
            # check job status
            job_status = 'unfinished'
            for file0 in os.listdir(job_dir_list[i_job_indx]):
                if '.otput' in file0:
                    otput_exist = True
                    otput_file_path = job_dir_list[i_job_indx] + '/' + file0
                    if dvm_job_finished(otput_file_path) == True:
                        job_status = 'finished'
                        # remove redundant files
                        for i_str in ['basfcn', 'bfnbak', 'grid3d', 'mixbak', 'moinfo', 'potent', 'scfbak']:
                            for file1 in os.listdir(job_dir_list[i_job_indx]):
                                if i_str in file1:
                                    redundant_file_path = job_dir_list[i_job_indx] + '/' + file1
                                    os.remove(redundant_file_path)
                    continue
            f.write(job_dir_name_list[i_job_indx] + ' ' * (max_length_job_dir_name - len(job_dir_name_list[i_job_indx])) + '   ' + job_status + ' ' * (10 - len(job_status)) + ' \n')
    funcs.write_log(logfile,
        'dvm_analyze.job_status(\n' +
        '    job_parent_dir = ' + 'r\'' + job_parent_dir + '\'' + '\n'
        '    )\n' +
        '#############\n'
        )

    return 0
