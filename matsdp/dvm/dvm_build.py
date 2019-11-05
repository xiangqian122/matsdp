def sphere_poscar(poscar_dir, origin_atom_name, radius = 7.0, include_mirror_atoms = True, dvm_incar_file_name = 'dvm_example'):
    '''
    Descriptions:
    This module is used to generate model (*.incar file) for the DVM program
    select atoms in a sphere that is centered around an arbitrary atom from the POSCAR file
    radius: float type. The atoms within a distance 'radius' from the original atom are selected (units in Angstroms)
    include_mirror_atoms: Logical value. Whether to include the mirror atoms or not
    dvm_incar_file_name: user-defined dvm job name
    '''
    import os
    import sys
    import numpy as np
    from .. import funcs
    from .. import convert
    from .. vasp import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)

    poscar_dir = os.path.abspath(poscar_dir)

    # the *.incar file which contains the selected atoms
    dvm_incar_file = output_dir + '/' + dvm_incar_file_name + '.incar'
    dvm_incar_temp_file1 = dvm_incar_file_name + '_1.temp'
    dvm_incar_temp_file2 = dvm_incar_file_name + '_2.temp'
    # the *.vasp file which contains the selected atoms
    vasp_file = output_dir + '/' + dvm_incar_file_name + '.vasp'
    vasp_temp_file1 = dvm_incar_file_name + 'vasp_1.temp'
    vasp_temp_file2 = dvm_incar_file_name + 'vasp_2.temp'
    # the *.xyz file which contains the selected atoms
    xyz_file = output_dir + '/' + dvm_incar_file_name + '.xyz'
    xyz_temp_file1 = dvm_incar_file_name + '_xyz_1.temp'
    xyz_temp_file2 = dvm_incar_file_name + '_xyz_2.temp'

    # Extract information from the input POSCAR file
    poscar_dict = vasp_read.read_poscar(poscar_dir)
    n_atoms = poscar_dict['n_atoms']
    atom_indx = convert.atomname2indx(poscar_dir, origin_atom_name)
    origin_atom_pos_arr = np.array([0.000]*3, dtype = np.float)
    origin_atom_pos_arr[0] = poscar_dict['pos_arr'][atom_indx - 1, 3]
    origin_atom_pos_arr[1] = poscar_dict['pos_arr'][atom_indx - 1, 4]
    origin_atom_pos_arr[2] = poscar_dict['pos_arr'][atom_indx - 1, 5]
    # expand the model according to the periodic boundary condition
    nx_hi = 0
    nx_lo = 0
    ny_hi = 0
    ny_lo = 0
    nz_hi = 0
    nz_lo = 0
    if include_mirror_atoms == True:
        if radius >= np.abs(poscar_dict['xhi'] - origin_atom_pos_arr[0]):
            nx_hi = int(np.ceil((radius - np.abs(poscar_dict['xhi'] - origin_atom_pos_arr[0]))/(poscar_dict['xhi'] - poscar_dict['xlo'])))
        if radius >= np.abs(poscar_dict['xlo'] - origin_atom_pos_arr[0]):
            nx_lo = int(np.ceil((radius - np.abs(poscar_dict['xlo'] - origin_atom_pos_arr[0]))/(poscar_dict['xhi'] - poscar_dict['xlo'])))
        if radius >= np.abs(poscar_dict['yhi'] - origin_atom_pos_arr[1]):
            ny_hi = int(np.ceil((radius - np.abs(poscar_dict['yhi'] - origin_atom_pos_arr[1]))/(poscar_dict['yhi'] - poscar_dict['ylo'])))
        if radius >= np.abs(poscar_dict['ylo'] - origin_atom_pos_arr[1]):
            ny_lo = int(np.ceil((radius - np.abs(poscar_dict['ylo'] - origin_atom_pos_arr[1]))/(poscar_dict['yhi'] - poscar_dict['ylo'])))
        if radius >= np.abs(poscar_dict['zhi'] - origin_atom_pos_arr[2]):
            nz_hi = int(np.ceil((radius - np.abs(poscar_dict['zhi'] - origin_atom_pos_arr[2]))/(poscar_dict['zhi'] - poscar_dict['zlo'])))
        if radius >= np.abs(poscar_dict['zlo'] - origin_atom_pos_arr[2]):
            nz_lo = int(np.ceil((radius - np.abs(poscar_dict['zlo'] - origin_atom_pos_arr[2]))/(poscar_dict['zhi'] - poscar_dict['zlo'])))
        n_atoms = n_atoms * (len(range(-nx_lo, nx_hi + 1)) * len(range(-ny_lo, ny_hi + 1)) * len(range(-nz_lo, nz_hi + 1)))
    expanded_pos_arr = np.array([None] * n_atoms * 3, dtype = np.float)
    expanded_pos_arr.shape = n_atoms, 3
    expanded_atom_elmtindx_arr = np.array([None] * n_atoms)
    expanded_atom_species_arr = np.array([None] * n_atoms)
    expanded_atomname_arr = np.array([None] * n_atoms)

    temp_atom_indx = 0
    for nx in range(-nx_lo, nx_hi + 1):
        for ny in range(-ny_lo, ny_hi + 1):
            for nz in range(-nz_lo, nz_hi + 1):
                for i_atom in range(0, poscar_dict['n_atoms']):
                    expanded_pos_arr[temp_atom_indx, :] = poscar_dict['pos_arr'][i_atom, 3:6] + nx * poscar_dict['l_arr'][0,:] + ny * poscar_dict['l_arr'][1,:] + nz * poscar_dict['l_arr'][2,:]
                    expanded_atom_elmtindx_arr[temp_atom_indx] = poscar_dict['atom_elmtindx_arr'][i_atom] + 1
                    expanded_atom_species_arr[temp_atom_indx] = poscar_dict['atom_species_arr'][i_atom]
                    expanded_atomname_arr[temp_atom_indx] = poscar_dict['atomname_list'][i_atom]
                    temp_atom_indx += 1

    xmin = np.min(expanded_pos_arr[:,0])
    ymin = np.min(expanded_pos_arr[:,1])
    zmin = np.min(expanded_pos_arr[:,2])

    # write the selected atoms to designated files (*.incar, *.vasp, *.xyz)
    dvm_atom_num = 0
    dvm_elmt_species_list = []
    dvm_elmt_count_list = []
    # write the body of the *.incar file
    with open(dvm_incar_temp_file1, 'w') as f_temp1, open(vasp_temp_file1, 'w') as f_vasp_temp1, open(xyz_temp_file1, 'w') as f_xyz_temp1: 
        pass
    with open(dvm_incar_temp_file1, 'a') as f_temp1, open(vasp_temp_file1, 'a') as f_vasp_temp1, open(xyz_temp_file1, 'a') as f_xyz_temp1: 
        for i_atom in range(0, n_atoms):
            vector1 = origin_atom_pos_arr
            vector2 = np.array([expanded_pos_arr[i_atom, 0], expanded_pos_arr[i_atom, 1], expanded_pos_arr[i_atom, 2]])
            dist = np.linalg.norm(vector1 - vector2)
            if dist <= radius:
                dvm_atom_num += 1
                if expanded_atom_species_arr[i_atom] not in dvm_elmt_species_list and expanded_atom_species_arr[i_atom] in poscar_dict['elmt_species_arr']:
                    dvm_elmt_species_list.append(expanded_atom_species_arr[i_atom])
                    dvm_elmt_count_list.append(0)
                dvm_elmt_count_list[dvm_elmt_species_list.index(expanded_atom_species_arr[i_atom])] += 1
                # *.incar
                f_temp1.write(
                    str('{:.6f}'.format(expanded_pos_arr[i_atom, 0])) + ' ' +
                    str('{:.6f}'.format(expanded_pos_arr[i_atom, 1])) + ' ' +
                    str('{:.6f}'.format(expanded_pos_arr[i_atom, 2])) + ' ' +
                    str(expanded_atom_elmtindx_arr[i_atom])+ '\n'
                    )
                # *.vasp
                if include_mirror_atoms == True:
                    f_vasp_temp1.write(
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 0] - xmin)) + ' ' +
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 1] - ymin)) + ' ' +
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 2] - zmin)) + ' ' +
                        expanded_atom_species_arr[i_atom] + str(dvm_elmt_count_list[dvm_elmt_species_list.index(expanded_atom_species_arr[i_atom])]) + ' ' +
                        '    ' + str(expanded_atomname_arr[i_atom]) + ' ' +
                        str(dvm_atom_num) + '\n'
                        )
                else: 
                    f_vasp_temp1.write(
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 0])) + ' ' +
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 1])) + ' ' +
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 2])) + ' ' +
                        expanded_atom_species_arr[i_atom] + str(dvm_elmt_count_list[dvm_elmt_species_list.index(expanded_atom_species_arr[i_atom])]) + ' ' +
                        '    ' + str(expanded_atomname_arr[i_atom]) + ' ' +
                        str(dvm_atom_num) + '\n'
                        )
                # *.xyz
                f_xyz_temp1.write(
                    str(expanded_atom_species_arr[i_atom])+ ' ' +
                    str('{:.6f}'.format(expanded_pos_arr[i_atom, 0])) + ' ' +
                    str('{:.6f}'.format(expanded_pos_arr[i_atom, 1])) + ' ' +
                    str('{:.6f}'.format(expanded_pos_arr[i_atom, 2])) + ' ' +
                    str(expanded_atomname_arr[i_atom]) + ' ' +
                    str(dvm_atom_num) + '\n'
                    )

    # write the header of the files (*.incar, *.vasp, *.xyz)
    temp_str = ''
    for i_atom_indx in range(1, dvm_atom_num + 1):
        temp_str = temp_str + str(i_atom_indx) + ' '
        if i_atom_indx % 10 == 0:
            temp_str = temp_str + '\n'
    formatted_dvm_header = temp_str

    temp_str = ''
    for i_elmt_species in dvm_elmt_species_list:
        temp_str = temp_str + funcs.grep('    ' + i_elmt_species, vasp_temp_file1)
    vasp_selected_atom_position = temp_str
    with open(vasp_temp_file1, 'w') as f_vasp_temp1:
        f_vasp_temp1.write(vasp_selected_atom_position)
    
    with open(dvm_incar_temp_file2, 'w') as f_temp2, open(vasp_temp_file2, 'w') as f_vasp_temp2, open(xyz_temp_file2, 'w') as f_xyz_temp2: 
        f_temp2.write(
            str(dvm_atom_num) + ' 3.0\n' + formatted_dvm_header + '\n\n'
            )
        if include_mirror_atoms == True:
            l_arr_shifted = poscar_dict['l_arr'].copy()
            l_arr_shifted[0,:] = l_arr_shifted[0,:] * len(range(-nx_lo, nx_hi + 1)) 
            l_arr_shifted[1,:] = l_arr_shifted[1,:] * len(range(-ny_lo, ny_hi + 1)) 
            l_arr_shifted[2,:] = l_arr_shifted[2,:] * len(range(-nz_lo, nz_hi + 1)) 
            f_vasp_temp2.write(
                'generated by matsdp. The atom coordinates are shifted (x y z dvm_model_atomname original_atomname dvm_atom_indx)\n1\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in l_arr_shifted[0,:]) + '\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in l_arr_shifted[1,:]) + '\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in l_arr_shifted[2,:]) + '\n' +
                ' '.join(dvm_elmt_species_list) + '\n' +
                ' '.join(str(i) for i in dvm_elmt_count_list) + '\n' +
                'Cartesian\n' 
                )
        else:
            f_vasp_temp2.write(
                'generated by matsdp. (x y z dvm_model_atomname original_atomname dvm_atom_indx)\n1\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in poscar_dict['l_arr'][0,:]) + '\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in poscar_dict['l_arr'][1,:]) + '\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in poscar_dict['l_arr'][2,:]) + '\n' +
                ' '.join(dvm_elmt_species_list) + '\n' +
                ' '.join(str(i) for i in dvm_elmt_count_list) + '\n' +
                'Cartesian\n' 
                )
        f_xyz_temp2.write(
            str(dvm_atom_num) + '\n' +
            '#generated by matsdp. (elmt_species x y z original_atomname dvm_atom_indx)\n'
            )
    #merge files
    funcs.merge_files(dvm_incar_temp_file2, dvm_incar_temp_file1)
    os.remove(dvm_incar_temp_file1)
    funcs.mv(dvm_incar_temp_file2, dvm_incar_file)

    funcs.merge_files(vasp_temp_file2, vasp_temp_file1)
    os.remove(vasp_temp_file1)
    funcs.mv(vasp_temp_file2, vasp_file)

    funcs.merge_files(xyz_temp_file2, xyz_temp_file1)
    os.remove(xyz_temp_file1)
    funcs.mv(xyz_temp_file2, xyz_file)
    
    funcs.write_log(logfile,
        'dvm_build.sphere_vasp(' + '\n' +
        '    poscar_dir=' + 'r\'' + str(poscar_dir) + '\'' + ',\n' +
        '    origin_atom_name=\'' + str(origin_atom_name) + '\',\n' + 
        '    radius=' + str(radius) + ',\n' +
        '    include_mirror_atoms=' + str(include_mirror_atoms) + ',\n' +
        '    dvm_incar_file_name=\'' + str(dvm_incar_file_name) + '\')\n')
    return 0 
