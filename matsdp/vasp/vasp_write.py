def write_poscar_with_atom_property(output_poscar_file_path, poscar_dict, added_atom_data = None, added_atom_property_str = None, added_atom_property_columns_str = None):
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
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    funcs.touch(output_poscar_file_path)
    n_atoms = np.sum(poscar_dict['elmt_num_arr'])

    formatted_len1 = 13 #For example len('-54.0812')=8, len('-54.081283453')=13   
    poscar_dict_header_temp = poscar_dict['header']
    if added_atom_data not in [None, 'None', 'none']:
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

    with open(output_poscar_file_path,'w') as f1:
        if added_atom_data in [None, 'None', 'none']:
            poscar_dict['header'][0] = (poscar_dict['header'][0].strip('\n').rstrip() +
                                        '\n'
                                        )
        elif added_atom_data not in [None, 'None', 'none']:
            poscar_dict['header'][0] = (poscar_dict['header'][0].strip('\n').rstrip() +
                                        ' added_atom_property=' + str(added_atom_property_str).strip('=').strip(':').strip('\n') + ': ' +
                                        added_atom_property_columns_str +
                                        '\n'
                                        )
        for i_line in range(0,len(poscar_dict['header'])):
            f1.write(poscar_dict['header'][i_line])
    with open(output_poscar_file_path,'a') as f1:
        added_atom_data_str = ''
        coord_str = ''
        fix_str = ''
        for i_atom in range(n_atoms):
            if added_atom_data not in [None, 'None', 'none']:
                temp_arr = added_atom_data[i_atom,:]
                added_atom_data_str = added_atom_data_str + ' '.join([str(temp_arr[indx]) + (' '*(max([len(str(x)) for x in added_atom_data[:, indx]]) - len(str(temp_arr[indx])))) for indx in range(0, len(temp_arr))] + '\n')
            coord_str = coord_str + ('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,0]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,0])))) + ' ' + 
                '{:.6f}'.format(poscar_dict['coord_arr'][i_atom,1]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,1])))) + ' ' + 
                '{:.6f}'.format(poscar_dict['coord_arr'][i_atom,2]) + ' ' + (' '*(formatted_len1-len('{:.6f}'.format(poscar_dict['coord_arr'][i_atom,2])))) + ' ' + '\n')

            fix_str = fix_str + (poscar_dict['fix_arr'][i_atom,0] + ' ' +
                poscar_dict['fix_arr'][i_atom,1] + ' ' +
                poscar_dict['fix_arr'][i_atom,2] + ' ' + '\n')

        if added_atom_data in [None, 'None', 'none']:
            if poscar_dict['slet_dyn_on'] == True:
                f1.write(coord_str +
                         fix_str +
                         '\n')
            elif poscar_dict['slet_dyn_on'] == False:
                f1.write(coord_str +
                         '\n')
        elif added_atom_data not in [None, 'None', 'none']:
            if poscar_dict['slet_dyn_on'] == True:
                f1.write(coord_str +
                         fix_str +
                         added_atom_data_str +
                         '\n')
            elif poscar_dict['slet_dyn_on'] == False:
                f1.write(coord_str +
                         added_atom_data_str +
                         '\n')
    output_poscar_dict = {}
    output_poscar_dict['output_poscar_file_path'] = output_poscar_file_path
    if added_atom_data not in [None, 'None', 'none']:
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
    import os
    from . import vasp_read
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])

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
        i_elmt_potcar_dir = os.path.join(elmt_potcar_dir, i_elmt_name)
        if not os.path.exists(i_elmt_potcar_dir) and not os.path.isfile(i_elmt_potcar_dir):
            i_elmt_potcar_dir_temp = i_elmt_potcar_dir + '_sv'
            if os.path.exists(i_elmt_potcar_dir_temp):
                print('# _sv is used' + ' for element ' + i_elmt_name)
                i_elmt_potcar_dir = i_elmt_potcar_dir_temp
            else:
                i_elmt_potcar_dir_temp = i_elmt_potcar_dir + '_pv'
                if os.path.exists(i_elmt_potcar_dir_temp):
                    print('# _pv is used' + ' for element ' + i_elmt_name)
                    i_elmt_potcar_dir = i_elmt_potcar_dir_temp
                else:
                    i_elmt_potcar_dir_temp = i_elmt_potcar_dir + '_s'
                    if os.path.exists(i_elmt_potcar_dir_temp):
                        print('# _s is used' + ' for element ' + i_elmt_name)
                        i_elmt_potcar_dir = i_elmt_potcar_dir_temp
                    else:
                        i_elmt_potcar_dir_temp = i_elmt_potcar_dir + '_h'
                        if os.path.exists(i_elmt_potcar_dir_temp):
                            print('# _h is used' + ' for element ' + i_elmt_name)
                            i_elmt_potcar_dir = i_elmt_potcar_dir_temp
        if not os.path.exists(i_elmt_potcar_dir) and not os.path.isfile(i_elmt_potcar_dir):
            print('ERROR (from vasp_write): Please check POTCAR')
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
