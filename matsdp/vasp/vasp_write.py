def write_poscar(output_poscar_file_path, poscar_dict, coord_system = 'Direct', added_atom_data = None, added_atom_property_str = None, added_atom_property_columns_str = None):
    '''
    write the poscar file with atom property data
    
    added_atom_data can be a 1-dimensional list, it can also be a 1-dimensional or 2-dimensional np.array
    '''
    args_dict = locals()
    import os
    import numpy as np
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    funcs.touch(output_poscar_file_path)
    n_atoms = np.sum(poscar_dict['elmt_num_arr'])

    formatted_len1 = 13 #For example len('-54.0812')=8, len('-54.081283453')=13   
    poscar_dict_header_temp = poscar_dict['header']
    if not added_atom_data is None:
        # after conversion, the added_atom_data turns to an np.array with the shape: n_atoms * n_col 
        added_atom_data = np.array(added_atom_data)
        if len(added_atom_data.shape) == 1:
            added_atom_data.shape = added_atom_data.shape[0], 1

        added_atom_property_columns_str = str(added_atom_property_columns_str).strip('=').strip(':').strip('\n')
        if added_atom_property_str in [None,'None','none']:
            added_atom_property_str = ''
        if added_atom_property_columns_str in [None,'None','none']:
            added_atom_property_columns_str = ''
        if ',' in added_atom_property_columns_str:
            added_atom_property_columns_str = ' '.join(funcs.split_line(line = added_atom_property_columns_str, separator = ','))
        else:
            added_atom_property_columns_str = ' '.join(funcs.split_line(line = added_atom_property_columns_str, separator = ' '))

    temp_str = ''
    coord_str = ''
    fix_str = ''
    added_atom_data_str = ''
    coord_str_arr = np.array([None] * poscar_dict['n_atoms'])
    fix_str_arr = np.array([None] * poscar_dict['n_atoms'])
    added_atom_data_str_arr = np.array([None] * poscar_dict['n_atoms'])

    if (coord_system.startswith('C') or coord_system.startswith('c')):
        column_span = slice(3,6)
    elif (coord_system.startswith('D') or coord_system.startswith('d')):
        column_span = slice(0,3)
    float_format1 = '{:17.9f}'

    temp_str = ''
    if added_atom_data is None:
        temp_str = temp_str + (poscar_dict['comment'] +
                                    '\n'
                                    )
    elif not added_atom_data is None:
        temp_str = temp_str + (poscar_dict['comment'] +
                                    ' added_atom_property=' + str(added_atom_property_str).strip('=').strip(':').strip('\n') + ': ' +
                                    added_atom_property_columns_str +
                                    '\n'
                                    )
    temp_str = temp_str + (
        str(poscar_dict['uni_scale_fac']) + '\n' + 
        str(' '.join(float_format1.format(x) for x in poscar_dict['box_len_arr'][0,:])) + '\n' + 
        str(' '.join(float_format1.format(x) for x in poscar_dict['box_len_arr'][1,:])) + '\n' + 
        str(' '.join(float_format1.format(x) for x in poscar_dict['box_len_arr'][2,:])) + '\n' + 
        str(' '.join(poscar_dict['elmt_species_arr'])) + '\n' + 
        str(' '.join(str(x) for x in poscar_dict['elmt_num_arr'])) + '\n' 
        )
    if poscar_dict['slet_dyn_on'] == True:
        temp_str = temp_str + ('Selective Dynamics\n' + 
            coord_system + '\n'
            )
    elif poscar_dict['slet_dyn_on'] == False:
        temp_str = temp_str + (
            coord_system + '\n'
            )
    for i_atom in range(n_atoms):
        if not added_atom_data is None:
            temp_arr = added_atom_data[i_atom,:]
            added_atom_data_str = ' '.join([str(temp_arr[indx]) + (' '*(max([len(str(x)) for x in added_atom_data[:, indx]]) - len(str(temp_arr[indx])))) for indx in range(0, len(temp_arr))]) + ' ' 
            added_atom_data_str_arr[i_atom] = added_atom_data_str
        coord_str = ' '.join(float_format1.format(x) for x in poscar_dict['pos_arr'][i_atom, column_span]) + ' '
        coord_str_arr[i_atom] = coord_str

        fix_str = (poscar_dict['fix_arr'][i_atom,0] + ' ' +
            poscar_dict['fix_arr'][i_atom,1] + ' ' +
            poscar_dict['fix_arr'][i_atom,2] + ' ')
        fix_str_arr[i_atom] = fix_str

    if added_atom_data is None:
        if poscar_dict['slet_dyn_on'] == True:
            for i_atom in range(n_atoms):
                temp_str = temp_str + (coord_str_arr[i_atom] +
                         fix_str_arr[i_atom] +
                         '\n')
        elif poscar_dict['slet_dyn_on'] == False:
            for i_atom in range(n_atoms):
                temp_str = temp_str + (coord_str_arr[i_atom] +
                         '\n')
    elif not added_atom_data is None:
        if poscar_dict['slet_dyn_on'] == True:
            for i_atom in range(n_atoms):
                temp_str = temp_str + (coord_str_arr[i_atom] +
                         fix_str_arr[i_atom] +
                         added_atom_data_str_arr[i_atom] +
                         '\n')
        elif poscar_dict['slet_dyn_on'] == False:
            for i_atom in range(n_atoms):
                temp_str = temp_str + (coord_str_arr[i_atom] +
                         added_atom_data_str_arr[i_atom] +
                         '\n')
    with open(output_poscar_file_path,'w') as f:
        f.write(temp_str)
    output_poscar_dict = {}
    output_poscar_dict['output_poscar_file_path'] = output_poscar_file_path
    if not added_atom_data is None:
        output_poscar_dict['added_atom_data_shape'] = added_atom_data.shape
    # recover the poscar_dict['header']
    poscar_dict['header'] = poscar_dict_header_temp
    return output_poscar_dict

def write_poscar_with_force(outcar_file_path, ionic_step = 'last', output_poscar_file_name = None):
    '''
    write poscar file with ionic force
    '''
    args_dict = locals()
    import os
    import numpy as np
    from . import vasp_read
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    outcar_file_path = os.path.abspath(outcar_file_path)    
    outcar_params_dict = vasp_read.read_outcar(outcar_file_path)
    workdir, dos_file = funcs.file_path_name(outcar_file_path)
    poscar_file_path = os.path.join(workdir, 'POSCAR')
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
    if output_poscar_file_name in [None,'None','none']:
        output_poscar_file_path = os.path.join(workdir, 'POSCAR_with_force_step_' + str(i_ionic_step + 1) + '.vasp')
    else:
        output_poscar_file_path = os.path.join(workdir, str(output_poscar_file_name) + '_step_' + str(i_ionic_step + 1) + '.vasp')
    write_poscar(
        output_poscar_file_path = output_poscar_file_path,
        poscar_dict = poscar_dict,
        added_atom_data = outcar_params_dict['force'][i_ionic_step,:,3:],
        added_atom_property_str = added_atom_property_str1,
        added_atom_property_columns_str = added_atom_property_columns_str1)
    write_poscar(
        output_poscar_file_path = output_poscar_file_path[0:-5] + '_absforce.vasp',
        poscar_dict = poscar_dict,
        added_atom_data = np.abs(outcar_params_dict['force'][i_ionic_step,:,3:]),
        added_atom_property_str = added_atom_property_str2,
        added_atom_property_columns_str = added_atom_property_columns_str2)
    ionic_step_list = [ionic_step]
    ##################################
    # Determine the args string
    ##################################
    log_str = ''
    func_name = 'vasp_write.write_poscar_with_force'
    args_str = func_name + '(' + '\n'
    for i_arg in args_dict.keys():
        arg_value = args_dict[i_arg]
        if isinstance(arg_value,str):
            arg_value_str = '\'' + arg_value + '\''
        else:
            arg_value_str = str(arg_value)

        if i_arg == 'outcar_file_path':
            arg_value_str = 'r\'' + str(outcar_file_path) + '\''
        if i_arg == 'ionic_step':
            arg_value_str = str(ionic_step_list).strip('[').strip(']')
        if i_arg == 'output_poscar_file_name':
            arg_value_str = '\'' + str(output_poscar_file_name) + '\''
        args_str += '    ' + i_arg + ' = ' + arg_value_str + ',\n'
    args_str += '    )\n'
    args_str += '################################################\n'
    log_str += args_str
    funcs.write_log(logfile, log_str)
    return 0

def write_potcar(poscar_path, elmt_potcar_dir):
    '''
    Description:
    prepare the POTCAR for a calculation
    /elmt_potcar_dir/
        /H/POTCAR.Z
        /He/POTCAR.Z
        ...
        /Ni/POTCAR.Z
        ...
        ...
    '''
    args_dict = locals()
    import os
    from . import vasp_read
    from .. import funcs
    from .. import default_params
    from .. import periodic_table

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    periodic_table_dict = periodic_table.periodic_tab()

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
    # untill now I can't find a way of reading POTCAR.Z file using Python, so we will call the system command 'zcat' in the Unix or the Linux
    potcar_file_path = os.path.join(workdir, 'POTCAR')
    funcs.touch(potcar_file_path)
    use_zcat = True
    use_gunzip = False
    use_python = False
    for i_elmt_name in elmt_species_arr:
        # use recommended setting of POTCAR of in the VASP manual (PAW)
        #i_elmt_potcar_dir = os.path.join(elmt_potcar_dir, i_elmt_name)
        recommended_potcar_str = i_elmt_name + periodic_table_dict['vasppot_paw'][i_elmt_name]
        i_elmt_potcar_dir = os.path.join(elmt_potcar_dir, recommended_potcar_str)
        #print(recommended_potcar_str)
        # if the recommended POTCAR directory does not exist, try other POTCAR settings
        if not os.path.exists(i_elmt_potcar_dir) and not os.path.isfile(i_elmt_potcar_dir):
            i_elmt_potcar_dir_temp = os.path.join(elmt_potcar_dir, i_elmt_name + '')
            if os.path.exists(i_elmt_potcar_dir_temp):
                print('# Recommended POTCAR ' + recommended_potcar_str + ' is not found. ' + i_elmt_name + ' is used for element ' + i_elmt_name)
                i_elmt_potcar_dir = i_elmt_potcar_dir_temp
            else:
                i_elmt_potcar_dir_temp = os.path.join(elmt_potcar_dir, i_elmt_name + '_sv')
                if os.path.exists(i_elmt_potcar_dir_temp):
                    print('# Recommended POTCAR ' + recommended_potcar_str + ' is not found. _sv is used for element ' + i_elmt_name)
                    i_elmt_potcar_dir = i_elmt_potcar_dir_temp
                else:
                    i_elmt_potcar_dir_temp = os.path.join(elmt_potcar_dir, i_elmt_name + '_pv')
                    if os.path.exists(i_elmt_potcar_dir_temp):
                        print('# Recommended POTCAR ' + recommended_potcar_str + ' is not found. _pv is used for element ' + i_elmt_name)
                        i_elmt_potcar_dir = i_elmt_potcar_dir_temp
                    else:
                        i_elmt_potcar_dir_temp = os.path.join(elmt_potcar_dir, i_elmt_name + '_s')
                        if os.path.exists(i_elmt_potcar_dir_temp):
                            print('# Recommended POTCAR ' + recommended_potcar_str + ' is not found. _s is used for element ' + i_elmt_name)
                            i_elmt_potcar_dir = i_elmt_potcar_dir_temp
                        else:
                            i_elmt_potcar_dir_temp = os.path.join(elmt_potcar_dir, i_elmt_name + '_h')
                            if os.path.exists(i_elmt_potcar_dir_temp):
                                print('# Recommended POTCAR ' + recommended_potcar_str + ' is not found. _h is used for element ' + i_elmt_name)
                                i_elmt_potcar_dir = i_elmt_potcar_dir_temp
                            else:
                                i_elmt_potcar_dir_temp = os.path.join(elmt_potcar_dir, i_elmt_name + '_d')
                                if os.path.exists(i_elmt_potcar_dir_temp):
                                    print('# Recommended POTCAR' + recommended_potcar_str + ' is not found. _d is used for element ' + i_elmt_name)
                                    i_elmt_potcar_dir = i_elmt_potcar_dir_temp
        if not os.path.exists(i_elmt_potcar_dir) and not os.path.isfile(i_elmt_potcar_dir):
            print('ERROR #2012271559 (from vasp_write): Please check POTCAR')
            exit() 
        elmt_potcar_file = 'POTCAR.Z'
        elmt_potcar_file_path = os.path.join(i_elmt_potcar_dir, elmt_potcar_file)
        #error handling
        if not os.path.exists(elmt_potcar_file_path):
            elmt_potcar_file = 'POTCAR'
            elmt_potcar_file_path_temp = os.path.join(i_elmt_potcar_dir, elmt_potcar_file)
            if os.path.exists(elmt_potcar_file_path_temp):
                elmt_potcar_file_path = elmt_potcar_file_path_temp
        # try to build POTCAR
        if use_zcat == True:
            try:
                os.system('zcat ' + elmt_potcar_file_path + ' >> ' + potcar_file_path + r'>/dev/null 2>&1')
            except:
                pass
        if i_elmt_name == elmt_species_arr[0] and os.path.exists(potcar_file_path) and os.path.getsize(potcar_file_path) == 0:
            use_zcat = False
            use_gunzip = True
            use_python = False
        if use_gunzip == True:
            try:
                os.system('gunzip -c ' + elmt_potcar_file_path + ' >> ' + potcar_file_path + r'>/dev/null 2>&1')
            except:
                print('vasp_write error: the call of the zcat and the gunzip -c command failed. Make sure you are in a Linux system and the zcat or the gunzip command can work successfully, or please check the files in the elemental POTCAR directory')
                funcs.write_log(
                    logfile,
                    '# vasp_write error: the call of the zcat and the gunzip -c command failedi. Make sure you are in a Linux system and the zcat or the gunzip command can work successfully, or please check the files in the elemental POTCAR directory')
                ##exit()
        if i_elmt_name == elmt_species_arr[0] and os.path.exists(potcar_file_path) and os.path.getsize(potcar_file_path) == 0:
            use_zcat = False
            use_gunzip = False
            use_python = True
        if use_python == True:
            try:
                with open(elmt_potcar_file_path, 'r') as f:
                    line = f.readlines()
                with open(potcar_file_path, 'a') as f1:
                    f1.writelines(line)
            except:
                print('WARNING (from vasp_write): Failed to build POTCAR')
    return 0

def write_incar(incar_file_path, incar_dict = None, mode = 'w'):
    '''
    Description:
    prepare the INCAR for a calculation
    mode: 'w', 'a', 's'. They denotes 'write', 'append', and 'substitute', respectively
    for the mode 's', if the input tag exists, the new tag will substitute the old tag. If the input tag is not in the INCAR, the new tag will be appended to the INCAR file.
    '''
    args_dict = locals()
    import os
    from collections import OrderedDict
    from . import vasp_read
    from . import vasp_default
    from . import vasp_tools
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])

    #INCAR file path
    incar_file_path = os.path.abspath(incar_file_path)
    # designate the working directory
    workdir, incar_file = funcs.file_path_name(incar_file_path)
    incar_dict = vasp_tools.ordered_incar_dict(incar_dict)
    #incar_dict = OrderedDict(incar_dict) 
    if mode in ['w', 'a']:
        with open(incar_file_path, mode) as f:
            if incar_dict in [None, 'None', 'none']:
                pass
            elif isinstance(incar_dict, dict):
                for key in incar_dict:
                    f.write(key + ' = ' + str(incar_dict[key]) + '\n')
    elif mode == 's':           
        with open(incar_file_path, 'r') as f:
            line = f.readlines()
            if incar_dict in [None, 'None', 'none']:
                pass
            elif isinstance(incar_dict, dict):
                found_key = False
                for key in incar_dict:
                    ##print(key, incar_dict[key])
                    for i_line_indx in range(0, len(line)):
                        if key in line[i_line_indx]:
                            found_key = True
                            line[i_line_indx] = (key + ' = ' + str(incar_dict[key]) + '\n')
                if found_key == False:
                    line.append('\n')
                    line.append(key + ' = ' + str(incar_dict[key]) + '\n')
            else:
                print('ERROR (from vasp_write): incar_dict must be dictionary')
                exit()
        with open(incar_file_path, 'w') as f:
            f.writelines(line)
    return 0 

def write_kpoints(kpoints_file_path, kpoints_dict = None, poscar_file_path = None, mode = 'w'):
    '''
    Description:
    prepare the KPOINTS for a calculation
    mode: 'w', 'a', 's'. They denotes 'write', 'append', and 'substitute', respectively
    '''
    args_dict = locals()
    import os
    from . import vasp_read
    from . import vasp_default
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])

    #KPOINTS file path
    kpoints_file_path = os.path.abspath(kpoints_file_path)
    if poscar_file_path not in [None, 'None', 'none']:
        poscar_file_path = os.path.abspath(poscar_file_path)
    # designate the working directory
    workdir, kpoints_file = funcs.file_path_name(kpoints_file_path)

    if mode in ['w', 'a']:
        with open(kpoints_file_path, mode) as f:
            if kpoints_dict in [None, 'None', 'none']:
                pass
            elif isinstance(kpoints_dict, dict):
                if kpoints_dict['scheme'] in ['a', 'A']:
                    f.write(kpoints_dict['comment'])
                    f.write(str(kpoints_dict['num_kpoints']))
                    f.write(str(kpoints_dict['scheme']))
                    f.write(str(kpoints_dict['length_param']))
                elif kpoints_dict['scheme'] in ['g', 'G', 'm', 'M']:
                    f.write(kpoints_dict['comment'])
                    f.write(str(kpoints_dict['num_kpoints']))
                    f.write(str(kpoints_dict['scheme']))
                    f.write(str(kpoints_dict['subdivisions_arr'][0]) + ' ' + str(kpoints_dict['subdivisions_arr'][1]) + ' ' + str(kpoints_dict['subdivisions_arr'][2]))
                    f.write(str(kpoints_dict['origin_shift_arr'][0]) + ' ' + str(kpoints_dict['origin_shift_arr'][1]) + ' ' + str(kpoints_dict['origin_shift_arr'][2]))
                elif kpoints_dict['scheme'] in ['c', 'C', 'k', 'K']:
                    pass
                elif kpoints_dict['scheme'] in ['l', 'L']:
                    # Line-mode
                    f.write(kpoints_dict['comment'])
                    f.write(str(kpoints_dict['num_intersections']))
                    f.write(str(kpoints_dict['scheme']))
                    f.write(str(kpoints_dict['coord_type']))
    elif mode == 's':           
        with open(kpoints_file_path, 'r') as f:
            line = f.readlines()
            if kpoints_dict in [None, 'None', 'none']:
                pass
            elif isinstance(kpoints_dict, dict):
                for key in kpoints_dict:
                    if key == 'comment':
                        line[0] = (str(kpoints_dict[key]) + '\n')
                    if key in ['num_kpoints', 'num_intersections']:
                        line[1] = (str(kpoints_dict[key]) + '\n')
                    if key == 'scheme':
                        line[2] = (str(kpoints_dict[key]) + '\n')
                    if key == 'subdivisions_arr':
                        line[3] = (str(kpoints_dict[key][0]) + ' ' + str(kpoints_dict[key][1]) + ' ' + str(kpoints_dict[key][2]) +  '\n')
                    if key == 'origin_shift_arr':
                        line[4] = (str(kpoints_dict[key][0]) + ' ' + str(kpoints_dict[key][1]) + ' ' + str(kpoints_dict[key][2]) +  '\n')
            else:
                print('ERROR (from vasp_write): kpoints_dict must be dictionary')
                exit()
        with open(kpoints_file_path, 'w') as f:
            f.writelines(line)
    return 0

def write_band_data(eigenval_or_procar_dict, e_fermi = 0, header = None):
    '''
    write band data into a file
    eigenval_or_procar_dict: the dictionary of EIGENVAL or PROCAR
    e_fermi: the value of E-fermi in OUTCAR
    '''
    args_dict = locals()
    import os
    import numpy as np
    file_path, filename = os.path.split(eigenval_or_procar_dict['file_path'])
    # Determine the output format of integers and floats
    int_len = len(str(eigenval_or_procar_dict['num_kpoints']))
    kpt_format = '{:' + str(int_len) + 'd}'
    eig_format = '{:15.6f}'
    #band_eigs_XXX.txt is the file containing information of the band structure
    kpt_number_arr = np.array(range(1, eigenval_or_procar_dict['num_kpoints'] + 1))
    if eigenval_or_procar_dict['ispin'] == 1:
        band_arr = np.concatenate((kpt_number_arr.reshape(-1,1), (eigenval_or_procar_dict['eigs'] - e_fermi)), axis = 1)
        band_data_file = os.path.join(file_path, 'band_eigs.txt')
        np.savetxt(band_data_file, band_arr, fmt = '%' + str(int_len) + 'i ' + '%15.6f ' * len(eigenval_or_procar_dict['eigs'][0,:]), delimiter = ',', header = header)
        band_data_file_1 = os.path.join(file_path, 'band.dat')
        with open(band_data_file_1, 'w') as f: 
            f.write('# ' + header + '\n')
        with open(band_data_file_1, 'a') as f: 
            for i_col in range(len(eigenval_or_procar_dict['eigs'][0,:])):
                temp_band_arr = np.concatenate((np.array(eigenval_or_procar_dict['kpath_len_list']).reshape(-1,1), (eigenval_or_procar_dict['eigs'] - e_fermi)[:, i_col].reshape(-1,1)), axis = 1)
                np.savetxt(f, temp_band_arr, fmt = '%11.6f %11.8e', delimiter = ',')
                f.write('\n')
    elif eigenval_or_procar_dict['ispin'] == 2:
        band_arr_up = np.concatenate((kpt_number_arr.reshape(-1,1), (eigenval_or_procar_dict['eigs_up'] - e_fermi)), axis = 1)
        band_arr_dn = np.concatenate((kpt_number_arr.reshape(-1,1), (eigenval_or_procar_dict['eigs_dn'] - e_fermi)), axis = 1)
        band_data_file_up = os.path.join(file_path, 'band_eigs_up.txt')
        band_data_file_dn = os.path.join(file_path, 'band_eigs_dn.txt')
        np.savetxt(band_data_file_up, band_arr_up, fmt = '%' + str(int_len) + 'i ' + '%15.6f ' * len(eigenval_or_procar_dict['eigs_up'][0,:]), delimiter = ',', header = header)
        np.savetxt(band_data_file_dn, band_arr_dn, fmt = '%' + str(int_len) + 'i ' + '%15.6f ' * len(eigenval_or_procar_dict['eigs_dn'][0,:]), delimiter = ',', header = header)
        band_data_file_1_up = os.path.join(file_path, 'band_up.dat')
        band_data_file_1_dn = os.path.join(file_path, 'band_dn.dat')
        with open(band_data_file_1_up, 'w') as f_up, open(band_data_file_1_dn, 'w') as f_dn: 
            f_up.write('# ' + header + '\n')
            f_dn.write('# ' + header + '\n')
        with open(band_data_file_1_up, 'a') as f_up,  open(band_data_file_1_dn, 'w') as f_dn: 
            for i_col in range(len(eigenval_or_procar_dict['eigs_up'][0,:])):
                temp_band_arr_up = np.concatenate((np.array(eigenval_or_procar_dict['kpath_len_list']).reshape(-1,1), (eigenval_or_procar_dict['eigs_up'] - e_fermi)[:, i_col].reshape(-1,1)), axis = 1)
                temp_band_arr_dn = np.concatenate((np.array(eigenval_or_procar_dict['kpath_len_list']).reshape(-1,1), (eigenval_or_procar_dict['eigs_dn'] - e_fermi)[:, i_col].reshape(-1,1)), axis = 1)
                np.savetxt(f_up, temp_band_arr_up, fmt = '%11.6f %11.8e', delimiter = ',')
                np.savetxt(f_dn, temp_band_arr_dn, fmt = '%11.6f %11.8e', delimiter = ',')
                f_up.write('\n')
                f_dn.write('\n')
    return 0
