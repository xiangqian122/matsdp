def write_poscar_with_atom_property(output_poscar_file_path, poscar_dict, added_atom_data, added_atom_property_str = None, added_atom_property_columns_str = None):
    '''
    write the poscar file with atom property data
    
    added_atom_data can be a 1-dimensional list, it can also be a 1-dimensional or 2-dimensional np.array
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)

    output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    funcs.touch(output_poscar_file_path)
    n_atoms = np.sum(poscar_dict['elmt_num_arr'])

    # after conversion, the added_atom_data turns to an np.array with the shape: n_atoms * n_col 
    added_atom_data = np.array(added_atom_data)
    if len(added_atom_data.shape) == 1:
        added_atom_data.shape = added_atom_data.shape[0], 1

    formatted_len1 = 13 #For example len('-54.0812')=8, len('-54.081283453')=13   
    poscar_dict_header_temp = poscar_dict['header'].copy()
    added_atom_property_columns_str = str(added_atom_property_columns_str).strip('=').strip(':').strip('\n')
    if added_atom_property_str == None:
        added_atom_property_str = ''
    if added_atom_property_columns_str == None:
        added_atom_property_columns_str = ''
    if ',' in added_atom_property_columns_str:
        added_atom_property_columns_str = ' '.join(funcs.split_line(line = added_atom_property_columns_str, separator = ','))
    else:
        added_atom_property_columns_str = ' '.join(funcs.split_line(line = added_atom_property_columns_str, separator = ' '))

    with open(output_poscar_file_path,'w') as f1:
        poscar_dict['header'][0] = (poscar_dict['header'][0].strip('\n').rstrip() +
                                    ' added_atom_property=' + str(added_atom_property_str).strip('=').strip(':').strip('\n') + ': ' +
                                    added_atom_property_columns_str +
                                    '\n'
                                    )
        for i_line in range(0,len(poscar_dict['header'])):
            f1.write(poscar_dict['header'][i_line])
    with open(output_poscar_file_path,'a') as f1:
        for i_atom in range(n_atoms):
            temp_arr = added_atom_data[i_atom,:]
            added_atom_data_str = ' '.join([str(temp_arr[indx]) + (' '*(max([len(str(x)) for x in added_atom_data[:, indx]]) - len(str(temp_arr[indx])))) for indx in range(0, len(temp_arr))])
            if poscar_dict['slet_dyn_on'] == True:
                f1.write('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,0]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,0])))) + ' ' +
                         '{:.6f}'.format(poscar_dict['coord_arr'][i_atom,1]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,1])))) + ' ' +
                         '{:.6f}'.format(poscar_dict['coord_arr'][i_atom,2]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,2])))) + ' ' +
                         poscar_dict['fix_arr'][i_atom,0] + ' ' +
                         poscar_dict['fix_arr'][i_atom,1] + ' ' +
                         poscar_dict['fix_arr'][i_atom,2] + ' ' +
                         added_atom_data_str +
                         '\n')
            elif poscar_dict['slet_dyn_on'] == False:
                f1.write('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,0]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,0])))) + ' ' +
                         '{:.6f}'.format(poscar_dict['coord_arr'][i_atom,1]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,1])))) + ' ' +
                         '{:.6f}'.format(poscar_dict['coord_arr'][i_atom,2]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,2])))) + ' ' +
                         added_atom_data_str +
                         '\n')
    output_poscar_dict = {}
    output_poscar_dict['output_poscar_file_path'] = output_poscar_file_path
    output_poscar_dict['added_atom_data_shape'] = added_atom_data.shape
    # recover the poscar_dict['header']
    poscar_dict['header'] = poscar_dict_header_temp
    return output_poscar_dict

def write_poscar_with_force(outcar_file_path, ionic_step = 'last', output_poscar_file_name = None):
    '''
    write poscar file with ionic force
    '''
    import os
    import numpy as np
    from . import vasp_read
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    outcar_file_path = os.path.abspath(outcar_file_path)    
    incar_params_dict, outcar_params_dict = vasp_read.read_outcar(outcar_file_path)
    workdir, dos_file = funcs.file_path_name(outcar_file_path)
    poscar_file_path = workdir + '/POSCAR'
    poscar_dict = vasp_read.read_poscar(poscar_file_path)
    # write the POSCAR file with atom property data
    added_atom_property_str1 = 'force'
    added_atom_property_columns_str1 = 'fx fy fz'
    added_atom_property_str2 = 'absforce'
    added_atom_property_columns_str2 = 'abs_fx abs_fy abs_fz'
    if ionic_step == 'last':
        i_ionic_step = outcar_params_dict['num_ionic_step'] - 1
    elif ionic_step == 'first':
        i_ionic_step = 1 - 1
    else:
        i_ionic_step = ionic_step - 1
        if ionic_step > outcar_params_dict['num_ionic_step']:
            i_ionic_step = outcar_params_dict['num_ionic_step']
    if output_poscar_file_name == None:
        output_poscar_file_path = workdir + '/POSCAR_with_force_step_' + str(i_ionic_step + 1) + '.vasp'
    else:
        None
    write_poscar_with_atom_property(
        output_poscar_file_path = output_poscar_file_path,
        poscar_dict = poscar_dict,
        added_atom_data = outcar_params_dict['force'][i_ionic_step,:,3:],
        added_atom_property_str = added_atom_property_str1,
        added_atom_property_columns_str = added_atom_property_columns_str1)
    write_poscar_with_atom_property(
        output_poscar_file_path = output_poscar_file_path[0:-5] + '_absforce.vasp',
        poscar_dict = poscar_dict,
        added_atom_data = np.abs(outcar_params_dict['force'][i_ionic_step,:,3:]),
        added_atom_property_str = added_atom_property_str2,
        added_atom_property_columns_str = added_atom_property_columns_str2)
    ionic_step_list = [ionic_step]
    funcs.write_log(
        logfile,
        'vasp_write.write_poscar_with_force(' + '\n' +
        '    outcar_file_path=' + 'r\'' + str(outcar_file_path) + '\'' + ',\n' +
        '    ionic_step=' + str(ionic_step_list).strip('[').strip(']') + ',\n' +
        '    output_poscar_file_name=' + '\'' + str(output_poscar_file_name) + '\'' + ')\n' +
        '###############################\n')

def write_potcar(poscar_path, elmt_potcar_dir):
    '''
    Description:
    prepare the POTCAR for a calculation
    '''
    import os
    from . import vasp_read
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']

    # atom position file and element potential file
    poscar_path = os.path.abspath(poscar_path)
    elmt_potcar_dir = os.path.abspath(elmt_potcar_dir)
    # designate the working directory
    workdir, poscar_file = funcs.file_path_name(poscar_path)
    
    poscar_dict = vasp_read.read_poscar(poscar_path)  
    with open(poscar_path, 'r') as f_poscar:
        lines = f_poscar.readlines()
        elmt_species_arr = poscar_dict['elmt_species_arr']
    # concatenate potential files of elements
    potcar_file = 'POTCAR'
    potcar_file_path = workdir + '/' + potcar_file 
    with open(potcar_file_path, 'w') as f_pot:
        pass
    with open(potcar_file_path, 'a') as f_pot:
        for i_elmt_name in elmt_species_arr:
            if program_name == 'VASP' or program_name == 'vasp':
                elmt_potcar_file = 'POTCAR.Z'
                elmt_potcar_file_path = elmt_potcar_dir + '/' + i_elmt_name + '/' + elmt_potcar_file
            with open(elmt_potcar_file_path, 'r') as f_elmt_pot:
                lines = f_elmt_pot.readlines()
                for i_line in lines:
                    f_pot.write(i_line)
            f_pot.write('\n')
    return 0
