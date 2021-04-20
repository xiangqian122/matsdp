# -*- coding: utf-8 -*-
def read_poscar(poscar_file_path, suppress_atom_number_ortho = True):
    '''
    Description:
        Read information from POSCAR file
    Args:
        @.poscar_file_path: String format. The directory of the POSCAR file. It can either be full path or relative path
    Return:
        poscar_dict
        l_arr: Numpy float array(np.float). Dimension is 3*3
        elmtindx_arr': Numpy integer array(np.int). Element index. Start from 0
        elmt_num_arr': Numpy integer array(np.int). Number of atoms for each element species
        elmt_species_arr': Numpy string array. element name in the periodic for each element species
        elmt_start_indx_arr': Numpy integer array(np.int). staring index for each element species in the atom_indxArr
        coord_system': String type. String value of either 'Direct' or 'Cartesian'
        uni_scale_fac': Float type. Universal scale factor(latice constant)
        slet_dyn_on': Logic type. Selective dynamics on or off? True or False
        num_header': Number of header lines. i.e. the number of lines before atomic positions
    '''
    import os
    import numpy as np
    import sys
    import time
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])

    poscar_file_path = os.path.abspath(poscar_file_path)
    poscar_dict = {}

    file_status = funcs.file_status(poscar_file_path)
    poscar_dict['file_status'] = file_status
    poscar_dict['file_path'] = poscar_file_path

    if file_status != 1:
        print('#WARNING #20120304 (from read_poscar): The file ' + poscar_file_path + ' does not exist or is empty, please check this file.')
    else:
        with open(poscar_file_path) as f:
            lines = f.readlines()
            poscar_dict['comment'] = lines[0].strip('\n').rstrip()
            poscar_dict['uni_scale_fac'] = float(funcs.split_line(lines[1])[0])
            box_len_arr = np.array([0.0] * 3 * 3,dtype = np.float)
            box_len_arr.shape = 3, 3
            l_arr = np.array([0.0] * 3 * 3,dtype = np.float)
            l_arr.shape = 3, 3
            box_len_arr[0,0] = float(funcs.split_line(lines[2])[0]);box_len_arr[0,1] = float(funcs.split_line(lines[2])[1]);box_len_arr[0,2] = float(funcs.split_line(lines[2])[2])
            box_len_arr[1,0] = float(funcs.split_line(lines[3])[0]);box_len_arr[1,1] = float(funcs.split_line(lines[3])[1]);box_len_arr[1,2] = float(funcs.split_line(lines[3])[2])
            box_len_arr[2,0] = float(funcs.split_line(lines[4])[0]);box_len_arr[2,1] = float(funcs.split_line(lines[4])[1]);box_len_arr[2,2] = float(funcs.split_line(lines[4])[2])
            poscar_dict['box_len_arr'] = box_len_arr
            l_arr = box_len_arr * poscar_dict['uni_scale_fac']
            # get the basis length and angle 
            vec_a = l_arr[0,:]
            vec_b = l_arr[1,:]
            vec_c = l_arr[2,:]
            poscar_dict['vec_a'] = vec_a
            poscar_dict['vec_b'] = vec_b
            poscar_dict['vec_c'] = vec_c
            basis_vector_dict = funcs.basis_vector_info(vec_a, vec_b, vec_c)
            poscar_dict['len_vec_a'] = basis_vector_dict['len_vec_a']
            poscar_dict['len_vec_b'] = basis_vector_dict['len_vec_b']
            poscar_dict['len_vec_c'] = basis_vector_dict['len_vec_c']
            poscar_dict['angle_alpha_radian'] = basis_vector_dict['angle_alpha_radian']
            poscar_dict['angle_beta_radian'] = basis_vector_dict['angle_beta_radian']
            poscar_dict['angle_gamma_radian'] = basis_vector_dict['angle_gamma_radian']
            poscar_dict['angle_alpha_degree'] = basis_vector_dict['angle_alpha_degree']
            poscar_dict['angle_beta_degree'] = basis_vector_dict['angle_beta_degree']
            poscar_dict['angle_gamma_degree'] = basis_vector_dict['angle_gamma_degree']
            poscar_dict['box_volume'] = basis_vector_dict['box_volume']
            poscar_dict['reciprocal_arr'] = basis_vector_dict['reciprocal_arr']
            poscar_dict['car2fra_matrix_arr'] = basis_vector_dict['car2fra_matrix_arr']
            poscar_dict['fra2car_matrix_arr'] = basis_vector_dict['fra2car_matrix_arr']
            poscar_dict['reciprocal_arr'] = basis_vector_dict['reciprocal_arr']

            elmt_species_info = False
            # extract string from the 6th line
            sixth_line_alpha_str_list = funcs.extract_alpha_str(lines[5])
            sixth_line_num_list = funcs.extract_num(lines[5])
            if len(sixth_line_alpha_str_list) == 0 and len(sixth_line_num_list) != 0:
                elmt_species_info = False
                temp_line_indx = 6
            elif len(sixth_line_alpha_str_list) != 0 and len(sixth_line_num_list) == 0:
                elmt_species_info = True
                temp_line_indx = 7
            ##print(sixth_line_alpha_str_list)
            ##print(sixth_line_num_list)
            first_char = lines[temp_line_indx][0:1]
            if first_char == "S" or first_char == "s":
                poscar_dict['slet_dyn_on'] = True 
                poscar_dict['num_header']  = temp_line_indx + 2
                first_char = lines[poscar_dict['num_header'] -1][0:1]
                if first_char == "D" or first_char == "d":
                    poscar_dict['coord_system'] = "Direct"
                elif first_char == "C" or first_char == "c":
                    poscar_dict['coord_system'] = "Cartesian"
            elif first_char == "D" or first_char == "d" or first_char == "C" or first_char == "c":
                poscar_dict['slet_dyn_on'] = False
                poscar_dict['num_header'] = temp_line_indx + 1
                if first_char == "D" or first_char == "d":
                    poscar_dict['coord_system'] = "Direct"
                elif first_char == "C" or first_char == "c":
                    poscar_dict['coord_system'] = "Cartesian"
            #-save the header of the POSCAR
            poscar_dict['header'] = lines[0:poscar_dict['num_header']]
            #-Get number of element species
            n_species = len(funcs.split_line(lines[5]))
            #-Number of atoms in the system
            poscar_dict['elmtindx_arr'] = np.array([0]*n_species,dtype = np.int)
            poscar_dict['elmt_species_arr'] = np.array(['--']*n_species,dtype = np.str)
            poscar_dict['elmt_num_arr'] = np.array([0]*n_species,dtype = np.int)
            for i in range(0,n_species):
                poscar_dict['elmtindx_arr'][i] = i
                if elmt_species_info == True:
                    poscar_dict['elmt_species_arr'][i] = funcs.split_line(lines[5])[i]
                    poscar_dict['elmt_num_arr'][i] = funcs.split_line(lines[6])[i]
                elif elmt_species_info == False:
                    #poscar_dict['elmt_species_arr'][i] = funcs.split_line(lines[5])[i]
                    poscar_dict['elmt_num_arr'][i] = funcs.split_line(lines[5])[i]
            n_atoms = sum(poscar_dict['elmt_num_arr'][:])
            #Element start index
            poscar_dict['elmt_start_indx_arr'] = np.array([0]*n_species,dtype=np.int)
            for i in range(0,n_species):
                if i == 0:
                    poscar_dict['elmt_start_indx_arr'][i] = 1
                else:
                    poscar_dict['elmt_start_indx_arr'][i] = sum(poscar_dict['elmt_num_arr'][range(0,i)])+1

            #-Array atom_species_arr and atom coordinates
            atom_species_arr = np.array(['--']*n_atoms,dtype=np.str)
            atom_elmtindx_arr = np.array([None]*n_atoms)
            atom_subindx_arr = np.array([0]*n_atoms,dtype=np.int)
            atomname_list = ['']*n_atoms
            indx = 0;
            for i in range(0,n_species):
                atom_species_arr[range(indx,indx+poscar_dict['elmt_num_arr'][i])] = [poscar_dict['elmt_species_arr'][i]]*poscar_dict['elmt_num_arr'][i]
                atom_elmtindx_arr[range(indx,indx+poscar_dict['elmt_num_arr'][i])] = [poscar_dict['elmtindx_arr'][i]]*poscar_dict['elmt_num_arr'][i]
                for j in range(0,poscar_dict['elmt_num_arr'][i]):
                    i_atom = poscar_dict['elmt_start_indx_arr'][i] + j
                    atom_subindx_arr[i_atom-1] = j + 1
                    atomname_list[i_atom-1] = str(poscar_dict['elmt_species_arr'][i]) + str(atom_subindx_arr[i_atom-1])
                indx += poscar_dict['elmt_num_arr'][i]
            
            # user added atom data in the POSCAR file
            num_atom_data_column = len(funcs.split_line(lines[poscar_dict['num_header']]))
            if poscar_dict['slet_dyn_on'] == True:
                num_added_atom_data_column = num_atom_data_column - 6
            elif poscar_dict['slet_dyn_on'] == False:
                num_added_atom_data_column = num_atom_data_column - 3
            added_atom_data_arr = np.array([None] * n_atoms * num_added_atom_data_column)
            added_atom_data_arr.shape = n_atoms, num_added_atom_data_column        
            #Atom coordinates
            #pos_arr is the coordinates of all the atoms, first three columns are in Direct coordinate, last three columns are in Cartesian coordinate
            #coord_arr is the x, y, z coordinate in POSCAR, without the influence of lattice vector or scale factor
            coord_arr = np.array([0.0]*n_atoms*3,dtype = np.float)
            coord_arr.shape = n_atoms,3
            fix_arr = np.array(['T']*n_atoms*3,dtype = np.str)
            fix_arr.shape = n_atoms,3
            pos_arr = np.array([0.0]*n_atoms*6,dtype = np.float)
            pos_arr.shape = n_atoms,6
            l_inv_arr = np.linalg.inv(l_arr)
            for i in range(0,n_atoms):
                temp = funcs.split_line(lines[poscar_dict['num_header'] +i])
                coord_arr[i,0] = float(temp[0])
                coord_arr[i,1] = float(temp[1])
                coord_arr[i,2] = float(temp[2])
                if poscar_dict['slet_dyn_on'] == True:
                    fix_arr[i,0] = temp[3]
                    fix_arr[i,1] = temp[4]
                    fix_arr[i,2] = temp[5]
                    added_atom_data_arr[i,0:num_added_atom_data_column] = temp[6:(6+num_added_atom_data_column)]
                elif poscar_dict['slet_dyn_on'] == False:
                    added_atom_data_arr[i,0:num_added_atom_data_column] = temp[3:(3+num_added_atom_data_column)]
                if poscar_dict['coord_system'] == "Direct":
                    # fractional coordinate
                    pos_arr[i,0] = coord_arr[i,0]
                    pos_arr[i,1] = coord_arr[i,1]
                    pos_arr[i,2] = coord_arr[i,2]
                    # cartesian coordinate
                    pos_arr[i,3] = np.dot(coord_arr[i,:], l_arr[:,[0]])
                    pos_arr[i,4] = np.dot(coord_arr[i,:], l_arr[:,[1]])
                    pos_arr[i,5] = np.dot(coord_arr[i,:], l_arr[:,[2]])
                elif poscar_dict['coord_system'] == "Cartesian":
                    #cartesian coordinate
                    pos_arr[i,3] = coord_arr[i,0] * poscar_dict['uni_scale_fac']
                    pos_arr[i,4] = coord_arr[i,1] * poscar_dict['uni_scale_fac']
                    pos_arr[i,5] = coord_arr[i,2] * poscar_dict['uni_scale_fac']
                    # fractional coordinate
                    pos_arr[i,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(pos_arr[i,3:6])
                    ##pos_arr[i,0] = np.dot(l_inv_arr[0,:], coord_arr[i,:])
                    ##pos_arr[i,1] = np.dot(l_inv_arr[1,:], coord_arr[i,:])
                    ##pos_arr[i,2] = np.dot(l_inv_arr[2,:], coord_arr[i,:])

            poscar_dict['n_atoms'] = int(np.sum(poscar_dict['elmt_num_arr']))
            poscar_dict['num_elmts'] = len(atom_species_arr)
            poscar_dict['atom_species_arr'] = atom_species_arr
            poscar_dict['atom_elmtindx_arr'] = atom_elmtindx_arr
            poscar_dict['atom_subindx_arr'] = atom_subindx_arr
            poscar_dict['atomname_list'] = atomname_list
            poscar_dict['atom_indx_arr'] = np.array([ x + 1 for x in range(0, poscar_dict['n_atoms'])])  # atom index start from 1
            poscar_dict['atom_key_arr'] = np.array([ x for x in range(0, poscar_dict['n_atoms'])])  # atom key start from 0
            poscar_dict['pos_arr'] = pos_arr
            poscar_dict['coord_arr'] = coord_arr
            poscar_dict['fix_arr'] = fix_arr
            poscar_dict['l_arr'] = l_arr
            poscar_dict['l_inv_arr'] = l_inv_arr
            poscar_dict['xlo'] = np.min([poscar_dict['l_arr'][0,0], poscar_dict['l_arr'][1,0], poscar_dict['l_arr'][2,0], 0])
            poscar_dict['xhi'] = np.max([poscar_dict['l_arr'][0,0], poscar_dict['l_arr'][1,0], poscar_dict['l_arr'][2,0], 0])
            poscar_dict['ylo'] = np.min([poscar_dict['l_arr'][0,1], poscar_dict['l_arr'][1,1], poscar_dict['l_arr'][2,1], 0])
            poscar_dict['yhi'] = np.max([poscar_dict['l_arr'][0,1], poscar_dict['l_arr'][1,1], poscar_dict['l_arr'][2,1], 0])
            poscar_dict['zlo'] = np.min([poscar_dict['l_arr'][0,2], poscar_dict['l_arr'][1,2], poscar_dict['l_arr'][2,2], 0])
            poscar_dict['zhi'] = np.max([poscar_dict['l_arr'][0,2], poscar_dict['l_arr'][1,2], poscar_dict['l_arr'][2,2], 0])       
            poscar_dict['added_atom_data'] = added_atom_data_arr
            poscar_dict['added_atom_property'] = None
            poscar_dict['added_atom_property_columns'] = None
            if suppress_atom_number_ortho == False:
                # atom numbering according to the order of the three orthogonal directions. This may be time consuming for large number of ions.
                for number_order_text in ['xyz','yzx','zxy','xzy','yxz','zyx']:
                    poscar_dict['atom_number_ortho_' + number_order_text] = funcs.atom_number_ortho(
                        atom_key_arr = poscar_dict['atom_key_arr'],
                        pos_arr_cartesian= poscar_dict['pos_arr'][:,3:6],
                        delta = 0.05,
                        number_order = number_order_text
                        )
            if 'added_atom_property=' in poscar_dict['header'][0]:
                poscar_dict['added_atom_property'] = funcs.split_line(line = (funcs.split_line(line = poscar_dict['header'][0], separator = '=')[-1]), separator = ':')[0]
                poscar_dict['added_atom_property_columns'] = funcs.split_line(line = poscar_dict['header'][0], separator = ':')[-1]
    #print(poscar_dict['atom_number_ortho_yxz'])
    return poscar_dict

def read_outcar(outcar_file_path):
    '''
    Description:
        Read information from OUTCAR file
    Args:
        outcar_file_path: String format. The directory of the OUTCAR file. It can either be full path or relative path
    Return:
        @.outcar_params_dict
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import convert
    from .. import default_params
    import copy
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    outcar_file_path = os.path.abspath(outcar_file_path)
    work_dir, outcar_file = funcs.file_path_name(outcar_file_path)
    # Initialize the outcar_params_dict
    outcar_params_dict = {}
    outcar_params_dict['file_path'] = None
    outcar_params_dict['work_dir'] = None
    outcar_params_dict['num_elmt_type'] = None
    outcar_params_dict['num_elmt_type_list'] = None
    outcar_params_dict['file_status'] = None
    outcar_params_dict['read_status'] = None
    outcar_params_dict['num_atoms'] = None
    outcar_params_dict['num_elmt_type'] = None
    outcar_params_dict['num_elmt_type_list'] = None
    outcar_params_dict['num_ionic_step'] = None
    outcar_params_dict['LORBIT'] = None
    outcar_params_dict['ISPIN'] = None
    outcar_params_dict['e_fermi'] = None
    outcar_params_dict['e_fermi_mod'] = None
    outcar_params_dict['XC(G=0)'] = None
    outcar_params_dict['alpha+bet'] = None
    outcar_params_dict['TOTEN'] = None
    outcar_params_dict['energy_without_entropy'] = None
    outcar_params_dict['energy(sigma->0)'] = None
    outcar_params_dict['force'] = None
    outcar_params_dict['NBANDS'] = None
    outcar_params_dict['elapsed_time'] = None
    outcar_params_dict['NELM'] = None
    outcar_params_dict['NSW'] = None
    outcar_params_dict['initial_value_dict'] = {}

    outcar_params_dict['initial_value_dict'] = copy.deepcopy(outcar_params_dict)

    outcar_params_dict['file_path'] = outcar_file_path
    outcar_params_dict['work_dir'] = work_dir
    file_status = funcs.file_status(outcar_file_path)
    outcar_params_dict['file_status'] = file_status

    if file_status != 1:
        print('WARNING #20120312 (from read_outcar): The file ' + outcar_file_path + ' does not exist or is empty, please check this file.')
        outcar_params_dict['read_status'] = 0
    else:
        with open(outcar_file_path,'r') as f:
            line = f.readlines()
            # get the elapsed time
            kwd = 'timing'
            num_lines = len(line)
            num_tail_lines = 300
            num_tail_lines = min(num_lines, num_tail_lines)
            outcar_params_dict['num_ionic_step'] = funcs.grep('TOTAL-FORCE',outcar_file_path).count('TOTAL-FORCE')
            for indx in range(num_lines - 1, num_lines - num_tail_lines -1, -1):
                i_line = line[indx]
                if kwd in i_line:
                    temp_line_indx = indx + 6
                    outcar_params_dict['elapsed_time'] = float(funcs.split_line(line = line[temp_line_indx])[-1])
           
            def outcar_extract_info(line, outcar_params_dict):
                from .. import funcs
                
                i_ionic_step = 0
                for i in range(len(line)):
                    if 'ions per type' in line[i]:
                        try:
                            outcar_params_dict['num_elmt_type'] = len(funcs.split_line(line = line[i], separator='=')[-1].split())
                        except:
                            pass
                        try:
                            outcar_params_dict['num_elmt_type_list'] = funcs.split_line(line = funcs.split_line(line = line[i], separator='=')[-1], separator = ' ')[:]
                        except:
                            pass
                        try:
                            outcar_params_dict['num_atoms'] = np.sum([int(item) for item in outcar_params_dict['num_elmt_type_list']])
                        except:
                            pass
                        try:
                            outcar_params_dict['force'] = np.array([None] * outcar_params_dict['num_ionic_step'] * outcar_params_dict['num_atoms'] * 6)
                        except:
                            pass
                        try:
                            outcar_params_dict['force'].shape = outcar_params_dict['num_ionic_step'] , outcar_params_dict['num_atoms'] , 6
                        except:
                            pass
                    if 'LORBIT' in line[i]:
                        try:
                            outcar_params_dict['LORBIT'] = int(funcs.split_line(line = line[i],separator = '=')[-1].split()[0]) 
                        except:
                            pass
                    if 'NELM' in line[i] and 'NELMIN' in line[i]:
                        try:
                            outcar_params_dict['NELM'] = int(funcs.split_line(line = line[i],separator = ';')[0].split('=')[-1]) 
                        except:
                            pass
                    if 'NSW' in line[i]:
                        try:
                            outcar_params_dict['NSW'] = int(funcs.split_line(line = line[i],separator = '=')[-1].split()[0]) 
                        except:
                            pass
                    if 'ISPIN' in line[i]:
                        try:
                            outcar_params_dict['ISPIN'] = int(funcs.split_line(line = line[i],separator = '=')[-1].split()[0]) 
                        except:
                            pass
                    if 'E-fermi' in line[i] and len(funcs.split_line(line = line[i],separator = ':')) > 1:
                        try:
                            outcar_params_dict['e_fermi'] = float(funcs.split_line(line = line[i],separator = ':')[1].split()[0]) 
                        except:
                            pass
                        try:
                            outcar_params_dict['XC(G=0)'] = float(funcs.split_line(line = line[i],separator = ':')[2].split()[0]) 
                        except:
                            pass
                        try:
                            outcar_params_dict['alpha+bet'] = float(funcs.split_line(line = line[i],separator = ':')[-1].strip()) 
                        except:
                            pass
                        try:
                            outcar_params_dict['e_fermi_mod'] = outcar_params_dict['e_fermi'] + outcar_params_dict['alpha+bet']
                        except:
                            pass
                    if 'TOTEN' in line[i] and 'time' not in line[i] and len(funcs.split_line(line[i], '=')) > 1:
                        try:
                            outcar_params_dict['TOTEN'] = float(funcs.split_line(line[i], '=')[1].split()[0])
                        except:
                            pass
                    if 'energy' in line[i] and 'without' in line[i] and 'entropy' in line[i] and 'time' not in line[i] and len(funcs.split_line(line[i], '=')) == 3 and 'energy(sigma->0)' in line[i]:
                        try:
                            outcar_params_dict['energy_without_entropy'] = float(funcs.split_line(line[i], '=')[1].split()[0])
                        except:
                            pass
                    if 'energy' in line[i] and 'sigma' in line[i] and 'entropy' in line[i] and 'time' not in line[i]:
                        try:
                        #print(outcar_file_path, line[i])
                            outcar_params_dict['energy(sigma->0)'] = float(funcs.split_line(line[i], '=')[-1].split()[0])
                        except:
                            pass
                    if 'NBANDS' in line[i] and 'NKPTS' in line[i] and 'NKDIM' in line[i]:
                        try:
                            outcar_params_dict['NBANDS'] = int(funcs.split_line(line = line[i],separator = '=')[-1].split()[0]) 
                        except:
                            pass
                    if 'TOTAL-FORCE' in line[i] and outcar_params_dict['num_elmt_type'] != None:
                        try:
                            for i_atom in range(0, outcar_params_dict['num_atoms']):
                                i_atom_line = i + i_atom + 2
                                if len(funcs.split_line(line[i_atom_line])) != 6:
                                    # in case of the following situation, this problem still exists until vasp5.3:
                                    #POSITION                                       TOTAL-FORCE (eV/Angst)
                                    #-----------------------------------------------------------------------------------
                                    #0.02082      0.00016     -0.02030    473254.270891  -6163.613158-493186.317939
                                    print('WARNING (from vasp_read): The OUTCAR file format is not correct for matsdp to read\n' + outcar_file_path + '\n' + 'line number: ' + str(i_atom_line + 1) + '\n' + line[i_atom_line] + '\n' + 'trying to solve this problem ...\n')
                                    funcs.write_log(logfile, '# WARNING (from vasp_read): Please check the format of the file: ' + outcar_file_path + '\n' + '#line number: ' + str(i_atom_line + 1) + '\n' + '#' + line[i_atom_line] + '\n' + '#trying to solve this problem ...\n')
                                    # this formatting problem is caused by too large force on an atom
                                    line[i_atom_line] = ' '.join(funcs.split_line(line[i_atom_line])[0:3]) + ' ' + line[i_atom_line].strip('\n').strip()[-42:-28] + ' ' + line[i_atom_line].strip('\n').strip()[-28:-14] + ' ' + line[i_atom_line].strip('\n').strip()[-14:] + '\n'
                                    all_are_digits = True
                                    for item in funcs.split_line(line[i_atom_line]):
                                        try:
                                            float(item)
                                        except ValueError:
                                            # not a float
                                            all_are_digits = False
                                    if all_are_digits == True:
                                        print('the line ' + str(i_atom_line + 1) + ' is converted to: \n' + line[i_atom_line])
                                        funcs.write_log(logfile, '# the line ' + str(i_atom_line + 1) + ' is converted to: \n' + '#' + line[i_atom_line])
                                    else:
                                        print('WARNING #20120201 (from vasp_read). The problem has not been solved. Please check the OUTCAR file ' + outcar_file_path + '\n')
                                        funcs.write_log(logfile, '# WARNING #20120201 (from vasp_read): The problem has not been solved. Please check the OUTCAR file ' + outcar_file_path + '\n')
                                else:
                                    outcar_params_dict['force'][i_ionic_step,i_atom,:] = [float(item) for item in funcs.split_line(line[i_atom_line])]
                        except:
                            pass
                        i_ionic_step += 1
                return outcar_params_dict
            try:
                outcar_params_dict = outcar_extract_info(line, outcar_params_dict)
                outcar_params_dict['read_status'] = 1
            except:
                print('WARNING #2103301011 (from vasp_read.read_outcar). Error in reading OUTCAR file. Please check the OUTCAR file ' + outcar_file_path + '\n')
                outcar_params_dict['read_status'] = 0
    def check_outcar_params_dict(outcar_params_dict):    
        read_status = 1
        for i_key in outcar_params_dict.keys():
            i_value = outcar_params_dict[i_key]
            j_value = outcar_params_dict['initial_value_dict'][i_key]
            ij_equal = funcs.variables_equal(i_value, j_value)
            if ij_equal == True:
                print('WARNING #2103301030 (from vasp_read.read_outcar). Error in reading outcar_params_dict[\'' + str(i_key) + '\']. Please check the OUTCAR file ' + outcar_file_path + '\n')
                read_status = 0
        return read_status
    outcar_params_dict['read_status'] = check_outcar_params_dict(outcar_params_dict)
    return outcar_params_dict

def read_doscar(doscar_file_path, atom_indx, save_dos_arr = False):
    '''
    Description:
        Read information from DOSCAR file.
        The energy is not shifed. i.e., the Fermi energy is not shifted.
        spin down DOS is mutiplied by -1
        atom_indx = 0 denotes extracting the total DOS
        atom_indx > 0 denotes extracting the partial DOS of each atom
        For The NCol columns in DOSCAR, each column denotes:
            num_col = 10. energy, s, py, pz, px, dxy, dyz, dz2, dxz, dx2 
            num_col = 19. energy, s(up), s(dw), py(up), py(dw), pz(up), pz(dw), px(up), px(dw), dxy(up), dxy(dw), dyz(up), dyz(dw), dz2(up), dz2(dw), dxz(up), dxz(dw), dx2(up), dx2(dw) 
            num_col = 17. energy, s, py, pz, px, dxy, dyz, dz2, dxz, dx2, f-3, f-2, f-1, f0, f1, f2, f3 
            num_col = 33. energy, s(up), s(dw), py(up), py(dw), pz(up), pz(dw), px(up), px(dw), dxy(up), dxy(dw), dyz(up), dyz(dw), dz2(up), dz2(dw), dxz(up), dxz(dw), dx2(up), dx2(dw) , f-3(up) , f-3(dw), f-2(up), f-2(dw), f-1(up), f-1(dw), f0(up), f0(dw), f1(up), f1(dw), f2(up), f2(dw), f3(up), f3(dw) 
    Args:
        @.doscar_file_path: String format. The directory of the DOSCAR file. It can either be full path or relative path
        @.atom_indx: Integer format. The real atom index in the POSCAR. If there are N atoms then the atom indices are frim 1 to N. Note that atom_indx = 0 means to extract TDOS inoformation
        @save_dos_arr: logical value. Determine whether to save the dos_arr to a file or not.
    Return:
        @.dos_arr: nedos * NCol array type. It contains the density of states
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    doscar_file_path = os.path.abspath(doscar_file_path)
    workdir, dos_file = funcs.file_path_name(doscar_file_path)
    
    doscar_dict = {}
    file_status = funcs.file_status(doscar_file_path)
    doscar_dict['file_status'] = file_status
    
    if file_status != 1:
        print('WARNING #20120313 (from read_doscar): The file ' + outcar_file_path + ' does not exist or is empty, please check this file.')
    else:
        #Also Required files are OUTCAR and POSCAR
        with open(doscar_file_path,'r') as f:
            lines = f.readlines()
            n_atoms = int(funcs.split_line(lines[0])[0])
            nedos = int(funcs.split_line(lines[5])[2])
            e_fermi = float(funcs.split_line(lines[5])[3])
##            funcs.write_log(logfile, 'n_atoms = ' + str(n_atoms) + ' nedos = ' + str(nedos) + ' e_fermi(Unmodified) = ' + str(e_fermi))
            if isinstance(atom_indx, int) and atom_indx == 0:
                #TDOS
                start_line = 6
                end_line = start_line + nedos - 1
                NC = len(funcs.split_line(lines[start_line]))
                dos_arr = np.array([0.000000]*nedos*NC,dtype = np.float)
                dos_arr.shape = nedos,NC
                for i in range(0,nedos):
                    for j in range(0,NC):
                        dos_arr[i][j] = float(funcs.split_line(lines[start_line + i])[j]) * funcs.logic_retn_val(j % 2 == 0 and j > 0, -1, 1)
                dos_file = os.path.join(workdir, 'TDOS.dat')
            elif isinstance(atom_indx, int)  and atom_indx != 0:
                #DOS
                start_line = 6 + (nedos + 1) * atom_indx
                end_line = start_line + nedos - 1
                NC = len(funcs.split_line(lines[start_line]))
                dos_arr = np.array([0.000000] * nedos * NC,dtype = np.float)
                dos_arr.shape = nedos,NC
                if NC == 10 or NC == 17:
                    for i in range(0,nedos):
                        for j in range(0,NC):
                            dos_arr[i][j] = float(funcs.split_line(lines[start_line + i])[j])
                elif NC == 19 or NC == 33:
                    for i in range(0,nedos):
                        for j in range(0,NC):
                            dos_arr[i][j] = float(funcs.split_line(lines[start_line + i])[j]) * funcs.logic_retn_val(j % 2 == 0 and j > 0, -1, 1)
                dos_file = os.path.join(workdir, 'DOS' + str(atom_indx) + '.dat')
            elif isinstance(atom_indx, list):
                # This is for the case of LDOS of multiple atoms.
                for i_atom_indx in atom_indx:
                    start_line = 6 + (nedos + 1) * i_atom_indx
                    end_line = start_line + nedos - 1
                    NC = len(funcs.split_line(lines[start_line]))
                    dos_arr = np.array([0.000000] * nedos * NC,dtype = np.float)
                    dos_arr.shape = nedos,NC
                    if NC == 10 or NC == 17:
                        for i in range(0,nedos):
                            for j in range(0,NC):
                                dos_arr[i][j] = dos_arr[i][j] + float(funcs.split_line(lines[start_line + i])[j])
                    elif NC == 19 or NC == 33:
                        for i in range(0,nedos):
                            for j in range(0,NC):
                                dos_arr[i][j] = dos_arr[i][j] + float(funcs.split_line(lines[start_line + i])[j]) * funcs.logic_retn_val(j % 2 == 0 and j > 0, -1, 1)
                    dos_file = os.path.join(workdir, 'DOS_' + '+'.join([str(x) for x in atom_indx]) + '.dat')
                
            if save_dos_arr == True:
                np.savetxt(dos_file, dos_arr)
            else:
                pass

        doscar_dict['atom_indx'] = atom_indx
        doscar_dict['num_col'] = NC
        doscar_dict['dos_arr'] = dos_arr
        doscar_dict['energy'] = dos_arr[:,[0]]
        if doscar_dict['num_col'] == 3:
            doscar_dict['TDOS'] = dos_arr[:,[1]]
            doscar_dict['IntDOS'] = dos_arr[:,[2]]
        elif doscar_dict['num_col'] == 5:
            doscar_dict['TDOS_up'] = dos_arr[:,[1]]
            doscar_dict['TDOS_dw'] = dos_arr[:,[2]]
            doscar_dict['TDOS'] = doscar_dict['TDOS_up'] - doscar_dict['TDOS_dw']
            doscar_dict['IntDOS_up'] = dos_arr[:,[3]]
            doscar_dict['IntDOS_dw'] = dos_arr[:,[4]]
            doscar_dict['IntDOS'] = doscar_dict['IntDOS_up'] - doscar_dict['IntDOS_dw']
        elif doscar_dict['num_col'] == 10 or doscar_dict['num_col'] == 17:
            doscar_dict['s'] = dos_arr[:,[1]]
            doscar_dict['py'] = dos_arr[:,[2]]
            doscar_dict['pz'] = dos_arr[:,[3]]
            doscar_dict['px'] = dos_arr[:,[4]]
            doscar_dict['dxy'] = dos_arr[:,[5]]
            doscar_dict['dyz'] = dos_arr[:,[6]]
            doscar_dict['dz2'] = dos_arr[:,[7]]
            doscar_dict['dxz'] = dos_arr[:,[8]]
            doscar_dict['dx2'] = dos_arr[:,[9]]
            doscar_dict['p'] = doscar_dict['py'] + doscar_dict['pz'] + doscar_dict['px']
            doscar_dict['d'] = doscar_dict['dxy'] + doscar_dict['dyz'] + doscar_dict['dz2']  + doscar_dict['dxz'] + doscar_dict['dx2']
            doscar_dict['LDOS'] = doscar_dict['s'] + doscar_dict['p'] + doscar_dict['d']
            if doscar_dict['num_col'] == 17:
                doscar_dict['fm3'] = dos_arr[:,[10]]
                doscar_dict['fm2'] = dos_arr[:,[11]]
                doscar_dict['fm1'] = dos_arr[:,[12]]
                doscar_dict['f0'] = dos_arr[:,[13]]
                doscar_dict['f1'] = dos_arr[:,[14]]
                doscar_dict['f2'] = dos_arr[:,[15]]
                doscar_dict['f3'] = dos_arr[:,[16]]
                doscar_dict['f'] = doscar_dict['fm3'] + doscar_dict['fm2'] + doscar_dict['fm1']  + doscar_dict['f0'] + doscar_dict['f1'] + doscar_dict['f2'] + doscar_dict['f3']
                doscar_dict['LDOS'] = doscar_dict['s'] + doscar_dict['p'] + doscar_dict['d'] + doscar_dict['f']
        elif doscar_dict['num_col'] == 19 or doscar_dict['num_col'] == 33:
            doscar_dict['s_up'] = dos_arr[:,[1]]
            doscar_dict['s_dw'] = dos_arr[:,[2]]
            doscar_dict['py_up'] = dos_arr[:,[3]]
            doscar_dict['py_dw'] = dos_arr[:,[4]]
            doscar_dict['pz_up'] = dos_arr[:,[5]]
            doscar_dict['pz_dw'] = dos_arr[:,[6]]
            doscar_dict['px_up'] = dos_arr[:,[7]]
            doscar_dict['px_dw'] = dos_arr[:,[8]]
            doscar_dict['dxy_up'] = dos_arr[:,[9]]
            doscar_dict['dxy_dw'] = dos_arr[:,[10]]
            doscar_dict['dyz_up'] = dos_arr[:,[11]]
            doscar_dict['dyz_dw'] = dos_arr[:,[12]]
            doscar_dict['dz2_up'] = dos_arr[:,[13]]
            doscar_dict['dz2_dw'] = dos_arr[:,[14]]
            doscar_dict['dxz_up'] = dos_arr[:,[15]]
            doscar_dict['dxz_dw'] = dos_arr[:,[16]]
            doscar_dict['dx2_up'] = dos_arr[:,[17]]
            doscar_dict['dx2_dw'] = dos_arr[:,[18]]
            doscar_dict['s'] = doscar_dict['s_up'] - doscar_dict['s_dw']
            doscar_dict['p_up'] = doscar_dict['py_up'] + doscar_dict['pz_up'] + doscar_dict['px_up']
            doscar_dict['p_dw'] = doscar_dict['py_dw'] + doscar_dict['pz_dw'] + doscar_dict['px_dw']
            doscar_dict['py'] = doscar_dict['py_up'] - doscar_dict['py_dw']
            doscar_dict['pz'] = doscar_dict['pz_up'] - doscar_dict['pz_dw']
            doscar_dict['px'] = doscar_dict['px_up'] - doscar_dict['px_dw']               
            doscar_dict['p'] = doscar_dict['p_up'] - doscar_dict['p_dw']
            doscar_dict['d_up'] = doscar_dict['dxy_up'] + doscar_dict['dyz_up'] + doscar_dict['dz2_up'] + doscar_dict['dxz_up'] + doscar_dict['dx2_up']
            doscar_dict['d_dw'] = doscar_dict['dxy_dw'] + doscar_dict['dyz_dw'] + doscar_dict['dz2_dw'] + doscar_dict['dxz_dw'] + doscar_dict['dx2_dw']
            doscar_dict['dxy'] = doscar_dict['dxy_up'] - doscar_dict['dxy_dw']
            doscar_dict['dyz'] = doscar_dict['dyz_up'] - doscar_dict['dyz_dw']
            doscar_dict['dz2'] = doscar_dict['dz2_up'] - doscar_dict['dz2_dw']
            doscar_dict['dxz'] = doscar_dict['dxz_up'] - doscar_dict['dx2_dw']
            doscar_dict['dx2'] = doscar_dict['dx2_up'] - doscar_dict['dx2_dw']
            doscar_dict['d'] = doscar_dict['d_up'] - doscar_dict['d_dw']
            doscar_dict['LDOS'] = doscar_dict['s'] + doscar_dict['p'] + doscar_dict['d']
            doscar_dict['LDOS_up'] = doscar_dict['s_up'] + doscar_dict['p_up'] + doscar_dict['d_up']
            doscar_dict['LDOS_dw'] = doscar_dict['s_dw'] + doscar_dict['p_dw'] + doscar_dict['d_dw']
            if doscar_dict['num_col'] == 33:
                doscar_dict['fm3_up'] = dos_arr[:,[19]]
                doscar_dict['fm3_dw'] = dos_arr[:,[20]]
                doscar_dict['fm2_up'] = dos_arr[:,[21]]
                doscar_dict['fm2_dw'] = dos_arr[:,[22]]
                doscar_dict['fm1_up'] = dos_arr[:,[23]]
                doscar_dict['fm1_dw'] = dos_arr[:,[24]]
                doscar_dict['f0_up'] = dos_arr[:,[25]]
                doscar_dict['f0_dw'] = dos_arr[:,[26]]
                doscar_dict['f1_up'] = dos_arr[:,[27]]
                doscar_dict['f1_dw'] = dos_arr[:,[28]]
                doscar_dict['f2_up'] = dos_arr[:,[29]]
                doscar_dict['f2_dw'] = dos_arr[:,[30]]
                doscar_dict['f3_up'] = dos_arr[:,[31]]
                doscar_dict['f3_dw'] = dos_arr[:,[32]]
                doscar_dict['f_up'] = doscar_dict['fm3_up'] + doscar_dict['fm2_up'] + doscar_dict['fm1_up'] + doscar_dict['f0_up'] + doscar_dict['f1_up'] + doscar_dict['f2_up'] + doscar_dict['f3_up']
                doscar_dict['f_dw'] = doscar_dict['fm3_dw'] + doscar_dict['fm2_dw'] + doscar_dict['fm1_dw'] + doscar_dict['f0_dw'] + doscar_dict['f1_dw'] + doscar_dict['f2_dw'] + doscar_dict['f3_dw']
                doscar_dict['f'] = doscar_dict['f_up'] - doscar_dict['f_dw']
                doscar_dict['LDOS'] = doscar_dict['s'] + doscar_dict['p'] + doscar_dict['d'] + doscar_dict['f']
                doscar_dict['LDOS_up'] = doscar_dict['s_up'] + doscar_dict['p_up'] + doscar_dict['d_up'] + doscar_dict['f_up']
                doscar_dict['LDOS_dw'] = doscar_dict['s_dw'] + doscar_dict['p_dw'] + doscar_dict['d_dw'] + doscar_dict['f_dw']
    return doscar_dict

def read_oszicar(oszicar_file_path, save_fig = False, dpi = 100):
    '''Read OSZICAR'''
    import os
    import numpy as np
    from .. import funcs
    import copy

    try:
        import matplotlib.pyplot as plt
    except:
        pass

    oszicar_file_path = os.path.abspath(oszicar_file_path)
    work_dir, oszicar_file = funcs.file_path_name(oszicar_file_path)
    #initializethe OSZICAR dictionary
    oszicar_dict = {}
    oszicar_dict['file_status'] = None
    oszicar_dict['file_path'] = None
    oszicar_dict['num_ionic_steps'] = None
    oszicar_dict['ionic_step_line_indx_list'] = []
    oszicar_dict['num_electronic_steps_list'] = []
    oszicar_dict['mag_list'] = []
    oszicar_dict['etot_arr'] = np.array([])
    oszicar_dict['initial_value_dict'] = {}

    oszicar_dict['initial_value_dict'] = copy.deepcopy(oszicar_dict)

    oszicar_dict['file_status'] = funcs.file_status(oszicar_file_path)
    oszicar_dict['file_path'] = oszicar_file_path

    if oszicar_dict['file_status'] != 1:
        print('WARNING #20120307 (from read_oszicar): File ' + oszicar_file_path + ' does not exist or is empty. Please check this file.')
    else:
        with open(oszicar_file_path,'r') as f:
            line = f.readlines()

        def oszicar_extract_info(line, oszicar_dict):
            from .. import funcs
            # get the number of ionic steps
            num_ionic_steps = 0
            for i_line in range(len(line)):
                if len(funcs.split_line(line[i_line])) > 1:
                    if funcs.split_line(line[i_line])[1] == 'F=':
                        num_ionic_steps += 1
                        oszicar_dict['ionic_step_line_indx_list'].append(i_line)
                        mag_str = funcs.split_line(line = line[i_line], separator = '=')[-1]
                        mag_list = mag_str.replace('\t', '    ').strip('\n').strip().split(' ')
                        oszicar_dict['mag_list'].append(funcs.split_line(line = line[i_line], separator = '=')[-1])
            oszicar_dict['num_ionic_steps'] = num_ionic_steps
            # get num_electronic_steps_list
            ionic_step = 0
            for i_line in range(len(line)):
                if len(funcs.split_line(line[i_line])) > 1:
                    if funcs.split_line(line[i_line])[1] == 'F=':
                        # get information about electronic steps
                        last_electronic_step_found = True
                        if ionic_step == 0:
                            temp_line = 0
                        else:
                            temp_line = oszicar_dict['ionic_step_line_indx_list'][ionic_step - 1]
                        #print('INI, FIN=', i_line, temp_line)
                        for j_line in range(i_line, temp_line, -1):
                            #print('ionic=', ionic_step, 'j_line=',j_line)
                            if line[j_line][0:5] in ['CG : ', 'DIA: ', 'NONE ', 'RMM: ', 'DAV: ']:
                                if last_electronic_step_found == True:
                                    oszicar_dict['num_electronic_steps_list'].append(int(line[j_line][5:8]))
                                    last_electronic_step_found = False
                        ionic_step = ionic_step + 1
            #Get etot_arr
            # solely read etot_arr will avoid interupting getting num_electronic)steps_list if etot_arr cannot be read.
            oszicar_dict['etot_arr'] = np.array([None] * num_ionic_steps)
            ionic_step = 0
            for i_line in range(len(line)):
                if len(funcs.split_line(line[i_line])) > 1:
                    if funcs.split_line(line[i_line])[1] == 'F=':
                        # get information about electronic steps
                        last_electronic_step_found = True
                        if ionic_step == 0:
                            temp_line = 0
                        else:
                            temp_line = oszicar_dict['ionic_step_line_indx_list'][ionic_step - 1]
                        for j_line in range(i_line, temp_line, -1):
                            if len(funcs.split_line(line[i_line])) > 2:
                                #print('------------ionic = ',ionic_step) 
                                 oszicar_dict['etot_arr'][ionic_step] = float(funcs.split_line(line[i_line])[2])
                        ionic_step = ionic_step + 1

            return oszicar_dict
        try:
            oszicar_dict = oszicar_extract_info(line, oszicar_dict)
        except:
            print('WARNING #2103311236 (from vasp_read.read_oszicar). Error in reading OSZICAR file. Please check the OSZICAR file ' + oszicar_file_path + '\n')

        if save_fig == True:
            try:
                fig, ax = plt.subplots()
                plt.plot(range(len(etot_arr)), etot_arr, marker = '.')
                ax.set(xlabel = 'Ionic steps')
                ax.set(ylabel = '$E_{tot}(eV)$')
                etot_oszicar_figfile = os.path.join(work_dir, 'etot_oszicar.pdf')
                plt.savefig(etot_oszicar_figfile,dpi = dpi)
                plt.close
            except:
                pass
    def check_oszicar_dict(oszicar_dict):    
        for i_key in oszicar_dict.keys():
            i_value = oszicar_dict[i_key]
            j_value = oszicar_dict['initial_value_dict'][i_key]
            #if isinstance(j_value,(list,pd.core.series.Series,np.ndarray)):
            ij_equal = funcs.variables_equal(i_value, j_value)
            if ij_equal == True:
                print('WARNING #2103302253 (from vasp_read.read_oszicar). Error in reading oszicar_dict[\'' + str(i_key) + '\']. Please check the OSZICAR file ' + oszicar_file_path + '\n')
        return 0
    check_oszicar_dict(oszicar_dict)
    return oszicar_dict

def read_kpoints(kpoints_file_path):
    '''Read KPOINTS'''
    import os
    import numpy as np
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    kpoints_file_path = os.path.abspath(kpoints_file_path)

    work_dir, kpt_file_name = funcs.file_path_name(kpoints_file_path)
    poscar_file_path = os.path.join(work_dir, 'POSCAR')
    poscar_dict = read_poscar(poscar_file_path)
    poscar_file_status = poscar_dict['file_status']

    kpoints_dict = {}
    kpoints_dict['file_path'] = kpoints_file_path
    
    kpoints_file_status = funcs.file_status(kpoints_file_path)
    kpoints_dict['file_status'] = kpoints_file_status

    if kpoints_file_status != 1 or poscar_file_status != 1:
        print('WARNING #20120308 (from read_kpoints): ' + kpoints_file_path + ' or ' + poscar_file_path + ' does not exist or is empty!')
    else: 
        subdivisions_arr = np.array([None] * 3)
        origin_shift_arr = np.array([None] * 3)
        with open(kpoints_file_path,'r') as f:
            line = f.readlines()
            num_lines = len(line)

            kpoints_dict['comment'] = line[0]
            kpoints_dict['scheme'] = line[2][0]

            if kpoints_dict['scheme'] in ['a', 'A']:
                kpoints_dict['num_kpoints'] = int(funcs.split_line(line[1])[0])
                kpoints_dict['length_param'] = float(funcs.split_line(line[3])[0])
            elif kpoints_dict['scheme'] in ['g', 'G', 'm', 'M']:
                num_kpoints_temp = funcs.split_line(line[1])[0]
                if isinstance(num_kpoints_temp, str):
                    kpoints_dict['num_kpoints'] = 0
                else:
                    kpoints_dict['num_kpoints'] = int(num_kpoints_temp)

                subdivisions_arr_0 = funcs.split_line(line[3])[0]
                subdivisions_arr_1 = funcs.split_line(line[3])[1]
                subdivisions_arr_2 = funcs.split_line(line[3])[2]
                ##if isinstance(subdivisions_arr_0, str):
                ##    subdivisions_arr[0] = 0
                ##else:
                ##    subdivisions_arr[0] = int(float(subdivisions_arr_0))
                ##if isinstance(subdivisions_arr_1, str):
                ##    subdivisions_arr[1] = 0
                ##else:
                ##    subdivisions_arr[1] = int(float(subdivisions_arr_1))
                ##if isinstance(subdivisions_arr_2, str):
                ##    subdivisions_arr[2] = 0
                ##else:
                ##    subdivisions_arr[2] = int(float(subdivisions_arr_2))
                try:
                    subdivisions_arr[0] = int(float(subdivisions_arr_0))
                    subdivisions_arr[1] = int(float(subdivisions_arr_1))
                    subdivisions_arr[2] = int(float(subdivisions_arr_2))
                except:
                    #if the type of value is string, we denote it as a meaningless negative number
                    subdivisions_arr[0] = -999999 
                    subdivisions_arr[1] = -999999 
                    subdivisions_arr[2] = -999999 
                kpoints_dict['subdivisions_arr'] = subdivisions_arr

                origin_shift_arr_0 = funcs.split_line(line[4])[0]
                origin_shift_arr_1 = funcs.split_line(line[4])[1]
                origin_shift_arr_2 = funcs.split_line(line[4])[2]
                if isinstance(origin_shift_arr_0, str):
                    origin_shift_arr[0] = 0
                else:
                    origin_shift_arr[0] = float(origin_shift_arr_0)  
                if isinstance(origin_shift_arr_1, str):
                    origin_shift_arr[1] = 0
                else:
                    origin_shift_arr[1] = float(origin_shift_arr_1)  
                if isinstance(origin_shift_arr_2, str):
                    origin_shift_arr[2] = 0
                else:
                    origin_shift_arr[2] = float(origin_shift_arr_2)  
                kpoints_dict['origin_shift_arr'] = origin_shift_arr
            elif kpoints_dict['scheme'] in ['c', 'C', 'k', 'K']:
                pass
            elif kpoints_dict['scheme'] in ['l', 'L']:
                # Line-mode
                kpoints_dict['coord_type'] = None
                kpoints_dict['kpath'] = None
                kpoints_dict['num_intersections'] = int(funcs.split_line(line[1])[0])
                kpoints_dict['coord_type'] = str(funcs.split_line(line[3])[0])[0]
                kpath_nodes_coord_x_list = []
                kpath_nodes_coord_y_list = []
                kpath_nodes_coord_z_list = []
                kpath_nodes_label_list = []
                last_high_symm_kpt_label = None
                kpt_label_list = []
                #get high symmetry k points (units in cartesian coordinate)
                for i_line in range(4, num_lines):
                    if len(funcs.split_line(line[i_line])) == 0:
                        continue
                    elif len(funcs.split_line(line[i_line])) != 0:
                        if funcs.split_line(line[i_line])[0][0] == '#':
                            continue
                    high_symm_kpt_x = float(funcs.split_line(line[i_line])[0])
                    high_symm_kpt_y = float(funcs.split_line(line[i_line])[1])
                    high_symm_kpt_z = float(funcs.split_line(line[i_line])[2])
                    temp_x = high_symm_kpt_x
                    temp_y = high_symm_kpt_y
                    temp_z = high_symm_kpt_z
                    high_symm_kpt_label = ''.join(funcs.split_line(line[i_line])[3:]).strip('!').strip('\\')
                    ##print('original',high_symm_kpt_x, high_symm_kpt_y, high_symm_kpt_z, high_symm_kpt_label)
                    #############################################################
                    # if the coordinate is reciprocal, change it to Cartesian (See VASP manual: VASP the Guide)
                    # \vec{k}=x_{1}\vec{b}_{1}+x_{2}\vec{b}_{2}+x_{3}\vec{b}_{3}
                    #############################################################
                    if kpoints_dict['coord_type'] in ['R', 'r']:
                        high_symm_kpt_x = temp_x * poscar_dict['reciprocal_arr'][0,0] + temp_y * poscar_dict['reciprocal_arr'][1,0] + temp_z * poscar_dict['reciprocal_arr'][2,0]
                        high_symm_kpt_y = temp_x * poscar_dict['reciprocal_arr'][0,1] + temp_y * poscar_dict['reciprocal_arr'][1,1] + temp_z * poscar_dict['reciprocal_arr'][2,1]
                        high_symm_kpt_z = temp_x * poscar_dict['reciprocal_arr'][0,2] + temp_y * poscar_dict['reciprocal_arr'][1,2] + temp_z * poscar_dict['reciprocal_arr'][2,2]
                        ##print('x,y,z',high_symm_kpt_x, high_symm_kpt_y, high_symm_kpt_z, high_symm_kpt_label)
                    # Change the labeling of the band to LaTeX format
                    if high_symm_kpt_label in ['GAMMA', 'gamma', 'Gamma']:
                        high_symm_kpt_label = 'Gamma'
                    for i_greek_capital_lett in defaults_dict['greek_capital_letter_list']:
                        if i_greek_capital_lett in high_symm_kpt_label:
                            high_symm_kpt_label = high_symm_kpt_label.replace(i_greek_capital_lett, '\\' + i_greek_capital_lett)
                    high_symm_kpt_label = '$' + high_symm_kpt_label + '$'
                    ##if high_symm_kpt_label in ['Gamma', 'gamma', 'g']:
                    ##    high_symm_kpt_label = '$\Gamma$'
                    kpt_label_list.append(high_symm_kpt_label)
                    if high_symm_kpt_label != last_high_symm_kpt_label:
                        kpath_nodes_coord_x_list.append(high_symm_kpt_x)
                        kpath_nodes_coord_y_list.append(high_symm_kpt_y)
                        kpath_nodes_coord_z_list.append(high_symm_kpt_z)
                        kpath_nodes_label_list.append(high_symm_kpt_label)
                    last_high_symm_kpt_label = high_symm_kpt_label
                #print(kpath_nodes_label_list)   #for debugging purpose
                #check the continuity of the kpath
                kpath_node_counter_list = [0] * len(kpath_nodes_label_list)
                kpath_continuous_list = [False] * len(kpath_nodes_label_list)
                last_kpt_label = None
                start_indx = 0
                for i_kpt in range(0, len(kpath_nodes_label_list)):
                    i_kpt_label = kpath_nodes_label_list[i_kpt]
                    for j_kpt in range(start_indx, len(kpt_label_list)):
                        j_kpt_label = kpt_label_list[j_kpt]
                        if i_kpt_label == j_kpt_label:
                            kpath_node_counter_list[i_kpt] += 1
                        elif i_kpt_label != j_kpt_label:
                            start_indx = j_kpt
                            break
                for i_kpt in range(0, len(kpath_nodes_label_list)):
                    if kpath_node_counter_list[i_kpt] != 1:
                        kpath_continuous_list[i_kpt] = True
                # based on the kpath_continuous_list, we determine the left_right_list
                kpath_left_right_list = ['L'] * len(kpath_nodes_label_list)
                kpath_left_right_list[0] = 'R'
                for i_kpt in range(0, len(kpath_nodes_label_list)):
                    if kpath_continuous_list[i_kpt] == True:
                        kpath_left_right_list[i_kpt] = 'LR'
                    if kpath_continuous_list[i_kpt] == False and i_kpt != 0:
                        #find the nearest point in the kpath_continuous_list which is True
                        for temp_i_kpt in reversed(range(0,i_kpt)): 
                            if kpath_continuous_list[temp_i_kpt] == True and (i_kpt - temp_i_kpt) % 2 != 0:
                                kpath_left_right_list[i_kpt] = 'R'
                                continue
                            elif kpath_continuous_list[temp_i_kpt] == True and (i_kpt - temp_i_kpt) % 2 == 0:
                                kpath_left_right_list[i_kpt] = 'L'
                                continue
                kpath_left_right_list[-1] = 'L'
                #build xaxis tick list for high symmetry k points
                kpath_nodes_xaxis_tick_list = []
                last_high_symm_kpt_x = kpath_nodes_coord_x_list[0]
                last_high_symm_kpt_y = kpath_nodes_coord_x_list[0]
                last_high_symm_kpt_z = kpath_nodes_coord_x_list[0]
                last_xaxis_tick = 0
                for i_kpt in range(0, len(kpath_nodes_label_list)):
                    kpts_dist = np.linalg.norm([kpath_nodes_coord_x_list[i_kpt] - last_high_symm_kpt_x,
                                                kpath_nodes_coord_y_list[i_kpt] - last_high_symm_kpt_y,
                                                kpath_nodes_coord_z_list[i_kpt] - last_high_symm_kpt_z])
                    new_xaxis_tick = last_xaxis_tick + kpts_dist 
                    ##if (i_kpt % 2) == 0 and kpath_continuous_list[i_kpt] == False and kpath_continuous_list[i_kpt - 1] == False:
                    ##if (i_kpt % 2) == 0 and kpath_continuous_list[i_kpt] == False and kpath_continuous_list[i_kpt - 1] == False and i_kpt != len(kpath_nodes_label_list) - 1:
                    ##if kpath_continuous_list[i_kpt] == False and kpath_continuous_list[i_kpt - 1] == False and i_kpt != len(kpath_nodes_label_list) - 1:
                    if kpath_continuous_list[i_kpt] == False and kpath_left_right_list[i_kpt] == 'R':
                        #print('i_kpt=',i_kpt, 'last=', last_xaxis_tick, 'new=',new_xaxis_tick)
                        new_xaxis_tick = last_xaxis_tick
                    kpath_nodes_xaxis_tick_list.append(new_xaxis_tick)
                    last_high_symm_kpt_x = kpath_nodes_coord_x_list[i_kpt]    
                    last_high_symm_kpt_y = kpath_nodes_coord_y_list[i_kpt]
                    last_high_symm_kpt_z = kpath_nodes_coord_z_list[i_kpt]
                    last_xaxis_tick = new_xaxis_tick
                kpoints_dict['kpath_nodes_xaxis_tick_list'] = kpath_nodes_xaxis_tick_list
                #build kpoints_dict
                num_kpath_nodes = len(kpath_nodes_label_list)
                kpoints_dict['num_kpath_nodes'] = num_kpath_nodes
                kpoints_dict['kpath_nodes_list'] = kpath_nodes_label_list
                kpoints_dict['kpath_left_right_list'] = kpath_left_right_list
                kpath_nodes_arr = np.array([None] * num_kpath_nodes * 8)
                kpath_nodes_arr.shape = num_kpath_nodes, 8
                #bs_xaxis_indx = 0
                bs_xaxis_label = None
                bs_xaxis_label_arr = [None] * num_kpath_nodes
                for i_kpt in range(0, num_kpath_nodes):
                    kpath_nodes_arr[i_kpt, 0] = i_kpt
                    kpath_nodes_arr[i_kpt, 1] = kpath_nodes_coord_x_list[i_kpt]  
                    kpath_nodes_arr[i_kpt, 2] = kpath_nodes_coord_y_list[i_kpt]
                    kpath_nodes_arr[i_kpt, 3] = kpath_nodes_coord_z_list[i_kpt]
                    kpath_nodes_arr[i_kpt, 4] = kpath_nodes_label_list[i_kpt]
                    kpath_nodes_arr[i_kpt, 5] = format(kpath_nodes_xaxis_tick_list[i_kpt], '.9f')
                    kpath_nodes_arr[i_kpt, 6] = kpath_continuous_list[i_kpt]
                    kpath_nodes_arr[i_kpt, 7] = kpath_left_right_list[i_kpt]
                kpoints_dict['kpath_nodes_arr'] = kpath_nodes_arr
                # determine the k points in the xaxis for the band structure plot
                intersections_end_kpt_list = [False] * len(kpath_continuous_list)
                intersections_end_kpt_list[0] = False
                num_intersections_interval = 0
                for i_kpt in range(1, len(kpath_continuous_list)):
                    if i_kpt != len(kpath_continuous_list) - 1:
                        #if not (kpath_continuous_list[i_kpt - 1] == False and kpath_continuous_list[i_kpt] == False):
                        if kpath_left_right_list[i_kpt] == 'LR' or kpath_left_right_list[i_kpt] == 'L':
                            num_intersections_interval = num_intersections_interval + 1
                            intersections_end_kpt_list[i_kpt] = True
                    elif i_kpt == len(kpath_continuous_list) - 1:
                        num_intersections_interval = num_intersections_interval + 1
                        intersections_end_kpt_list[i_kpt] = True
                kpoints_dict['num_intersections_interval'] = num_intersections_interval
                num_kpoints_xaxis = num_intersections_interval * kpoints_dict['num_intersections'] 
                kpoints_dict['num_kpoints_xaxis'] = num_kpoints_xaxis
                kpoints_dict['intersections_end_kpt_list'] = intersections_end_kpt_list
                #to avoid the following error: mask = np.isnan(self.x) TypeError: ufunc 'isnan' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''.    Use dtype = np.float64 instead of dtype = object.
                kpoints_xaxis_arr = np.array([None] * num_kpoints_xaxis, dtype=np.float64)
                start_kpt = kpath_nodes_xaxis_tick_list[0]
                i_intersection = 0
                for i_kpt in range(1, len(intersections_end_kpt_list)):
                    if intersections_end_kpt_list[i_kpt] == True:
                        i_intersection = i_intersection + 1
                        start_indx = (i_intersection - 1) * kpoints_dict['num_intersections']
                        end_indx = start_indx + kpoints_dict['num_intersections']
                        end_kpt =  kpath_nodes_xaxis_tick_list[i_kpt]
                        start_kpt = kpath_nodes_xaxis_tick_list[i_kpt - 1]
                        kpoints_xaxis_arr[start_indx:end_indx] = np.linspace(start_kpt, end_kpt, kpoints_dict['num_intersections'], True)
                kpoints_dict['kpoints_xaxis_arr'] = kpoints_xaxis_arr
                #get x axis tick for the band structure plot 
                bs_xaxis_label_list = []
                bs_xaxis_tick_list = []
                for i_kpt in range(1, len(intersections_end_kpt_list)):
                    if i_kpt != len(intersections_end_kpt_list) - 1:
                        if i_kpt == 1:
                            bs_xaxis_label_list.append(kpath_nodes_label_list[i_kpt - 1])
                            bs_xaxis_tick_list.append(kpath_nodes_xaxis_tick_list[i_kpt - 1])
                        if intersections_end_kpt_list[i_kpt] == True and intersections_end_kpt_list[i_kpt - 1] != False and i_kpt != 1:
                            bs_xaxis_label_list.append(kpath_nodes_label_list[i_kpt - 1])
                            bs_xaxis_tick_list.append(kpath_nodes_xaxis_tick_list[i_kpt - 1])
                        ##elif intersections_end_kpt_list[i_kpt] == True and intersections_end_kpt_list[i_kpt - 1] == False and i_kpt != 1:
                        ##    bs_xaxis_label_list.append('')
                        ##    bs_xaxis_tick_list.append(kpath_nodes_xaxis_tick_list[i_kpt - 1])
                        elif intersections_end_kpt_list[i_kpt] == False and kpath_left_right_list[i_kpt] == 'R':
                            bs_xaxis_label_list.append(kpath_nodes_label_list[i_kpt - 1] + '|' + kpath_nodes_label_list[i_kpt])
                            bs_xaxis_tick_list.append(kpath_nodes_xaxis_tick_list[i_kpt - 1])
                    elif i_kpt == len(intersections_end_kpt_list) - 1:
                            if intersections_end_kpt_list[i_kpt - 1] == True:
                                bs_xaxis_label_list.append(kpath_nodes_label_list[i_kpt - 1])
                                bs_xaxis_label_list.append(kpath_nodes_label_list[i_kpt])
                                bs_xaxis_tick_list.append(kpath_nodes_xaxis_tick_list[i_kpt - 1])
                                bs_xaxis_tick_list.append(kpath_nodes_xaxis_tick_list[i_kpt])
                            elif intersections_end_kpt_list[i_kpt - 1] == False:
                                bs_xaxis_label_list.append(kpath_nodes_label_list[i_kpt])
                                bs_xaxis_tick_list.append(kpath_nodes_xaxis_tick_list[i_kpt])
                kpoints_dict['bs_xaxis_label_list'] = bs_xaxis_label_list
                kpoints_dict['bs_xaxis_tick_list'] = bs_xaxis_tick_list
    return kpoints_dict

def read_eigenval(eigenval_file_path):
    '''
    Read EIGENVAL
    '''
    import os
    import numpy as np
    from .. import funcs
    from . import vasp_tools

    eigenval_file_path = os.path.abspath(eigenval_file_path)
    file_path, filename = os.path.split(eigenval_file_path)
    eigenval_dict = {}
    eigenval_dict['file_path'] = eigenval_file_path
    eigenval_dict['dict_type'] = 'eigenval'

    file_status = funcs.file_status(eigenval_file_path)
    eigenval_dict['file_status'] = file_status
    # read_status: 1: sucessfully read; 0: fail to read
    eigenval_dict['read_status'] = None
    eigenval_dict['num_lines'] = None
    eigenval_dict['num_ions'] = None
    eigenval_dict['loops_after_write_apcf_dos'] = None
    eigenval_dict['ispin'] = None
    eigenval_dict['cell_volume'] = None
    eigenval_dict['box_size_list'] = None
    eigenval_dict['num_valence_electrons'] = None
    eigenval_dict['num_kpoints'] = None
    eigenval_dict['num_bands'] = None
    eigenval_dict['kpt_coord_arr'] = None
    eigenval_dict['weights'] = None
    eigenval_dict['eigs'] = None
    eigenval_dict['eigs_up'] = None
    eigenval_dict['eigs_dw'] = None
    eigenval_dict['occupancy'] = None
    eigenval_dict['occupancy_up'] = None
    eigenval_dict['occupancy_dw'] = None

    if file_status != 1:
        print('WARNING #20120310 (from read_eigenval): The file ' + eigenval_dict['file_path'] + ' does not exist or is empty. Please check the EIGENVAL file.')
        eigenval_dict['read_status'] = 0
    else:
        file_type = vasp_tools.check_file_type(eigenval_dict['file_path'])
        if file_type != 'EIGENVAL':
            print('WARNING #20120301 (from read_eigenval): Incorrect EIGENVAL format. Please check the file' + eigenval_dict['file_path'])
            eigenval_dict['read_status'] = 0
        elif file_type == 'EIGENVAL':
            num_header = 7

            with open(eigenval_file_path,'r') as f:
                line = f.readlines()
                eigenval_dict['num_lines'] = len(line)
                eigenval_dict['num_ions'] = int(funcs.split_line(line[0])[0])
                eigenval_dict['loops_after_write_apcf_dos'] = int(funcs.split_line(line[0])[2]) # number of loops after writing the averaged pair correlation functions and DOS.
                #print(eigenval_file_path)
                ##print(line[0])
                eigenval_dict['ispin'] = int(funcs.split_line(line[0])[3])
                eigenval_dict['cell_volume'] = float(funcs.split_line(line[1])[0])
                eigenval_dict['box_size_list'] = [float(funcs.split_line(line[1])[1]), float(funcs.split_line(line[1])[2]), float(funcs.split_line(line[1])[3])]
                eigenval_dict['nnum_valence_electrons'] = int(funcs.split_line(line[5])[0])
                eigenval_dict['num_kpoints'] = int(funcs.split_line(line[5])[1])
                eigenval_dict['num_bands'] = int(funcs.split_line(line[5])[2])
                    
                #to avoid the following error: mask = np.isnan(self.x) TypeError: ufunc 'isnan' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''.    Use dtype = np.float64 instead of dtype = object.
                if eigenval_dict['ispin'] == 1:
                    eigenval_dict['eigs'] = np.array([None] * eigenval_dict['num_kpoints'] * eigenval_dict['num_bands'], dtype=np.float64)
                    eigenval_dict['eigs'].shape = eigenval_dict['num_kpoints'], eigenval_dict['num_bands']
                    eigenval_dict['occupancy'] = np.array([None] * eigenval_dict['num_kpoints'] * eigenval_dict['num_bands'], dtype=np.float64)
                    eigenval_dict['occupancy'].shape = eigenval_dict['num_kpoints'], eigenval_dict['num_bands']
                elif eigenval_dict['ispin'] == 2:
                    eigenval_dict['eigs_up'] = np.array([None] * eigenval_dict['num_kpoints'] * eigenval_dict['num_bands'], dtype=np.float64)
                    eigenval_dict['eigs_up'].shape = eigenval_dict['num_kpoints'], eigenval_dict['num_bands']
                    eigenval_dict['eigs_dw'] = np.array([None] * eigenval_dict['num_kpoints'] * eigenval_dict['num_bands'], dtype=np.float64)
                    eigenval_dict['eigs_dw'].shape = eigenval_dict['num_kpoints'], eigenval_dict['num_bands']
                    eigenval_dict['occupancy_up'] = np.array([None] * eigenval_dict['num_kpoints'] * eigenval_dict['num_bands'], dtype=np.float64)
                    eigenval_dict['occupancy_up'].shape = eigenval_dict['num_kpoints'], eigenval_dict['num_bands']
                    eigenval_dict['occupancy_dw'] = np.array([None] * eigenval_dict['num_kpoints'] * eigenval_dict['num_bands'], dtype=np.float64)
                    eigenval_dict['occupancy_dw'].shape = eigenval_dict['num_kpoints'], eigenval_dict['num_bands']

                eigenval_dict['kpt_coord_arr'] = np.array([None] * eigenval_dict['num_kpoints'] * 3)     
                eigenval_dict['kpt_coord_arr'].shape = eigenval_dict['num_kpoints'], 3
                eigenval_dict['weights'] = np.array([None] * eigenval_dict['num_kpoints'])

                i_line = num_header
                if eigenval_dict['ispin'] == 1:
                    for i_kpoint in range(0, eigenval_dict['num_kpoints']):
                        eigenval_dict['kpt_coord_arr'][i_kpoint, 0] = float(funcs.split_line(line[i_line])[0])
                        eigenval_dict['kpt_coord_arr'][i_kpoint, 1] = float(funcs.split_line(line[i_line])[1])
                        eigenval_dict['kpt_coord_arr'][i_kpoint, 2] = float(funcs.split_line(line[i_line])[2])
                        eigenval_dict['weights'][i_kpoint] = float(funcs.split_line(line[i_line])[3])
                        for i_band in range(0, eigenval_dict['num_bands']):
                            i_line = i_line + 1
                            try:
                                eigenval_dict['eigs'][i_kpoint,i_band] = float(funcs.split_line(line[i_line])[1])
                            except:
                                # This happens when the parameter cannot be converted to string, for example: sting: '***************'
                                eigenval_dict['eigs'][i_kpoint,i_band] = 999999
                            # Currently, the occupancy will not be read, because the format is unclear for now
                            ##occupancy[i_kpoint, i_band] = float(funcs.split_line(line[i_line])[2])
                            if i_band == (eigenval_dict['num_bands'] - 1):
                                i_line = i_line + 2
                elif eigenval_dict['ispin'] == 2:
                    for i_kpoint in range(0, eigenval_dict['num_kpoints']):
                        eigenval_dict['kpt_coord_arr'][i_kpoint, 0] = float(funcs.split_line(line[i_line])[0])
                        eigenval_dict['kpt_coord_arr'][i_kpoint, 1] = float(funcs.split_line(line[i_line])[1])
                        eigenval_dict['kpt_coord_arr'][i_kpoint, 2] = float(funcs.split_line(line[i_line])[2])
                        eigenval_dict['weights'][i_kpoint] = float(funcs.split_line(line[i_line])[3])
                        for i_band in range(0, eigenval_dict['num_bands']):
                            i_line = i_line + 1
                            try:
                                eigenval_dict['eigs_up'][i_kpoint,i_band] = float(funcs.split_line(line[i_line])[1])
                                eigenval_dict['eigs_dw'][i_kpoint,i_band] = float(funcs.split_line(line[i_line])[2])
                            except:
                                # This happens when the parameter cannot be converted to string, for example: sting: '***************'
                                eigenval_dict['eigs_up'][i_kpoint,i_band] = 999999
                                eigenval_dict['eigs_dw'][i_kpoint,i_band] = 999999
                            # Currently, the occupancy will not be read, because the format is unclear for now
                            #occupancy_up[i_kpoint, i_band] = float(funcs.split_line(line[i_line])[3])
                            #occupancy_dw[i_kpoint, i_band] = float(funcs.split_line(line[i_line])[4])
                            if i_band == (eigenval_dict['num_bands'] - 1):
                                i_line = i_line + 2
            eigenval_dict['read_status'] = 1
                
            band_data_txt = ''
            band_data_txt_up = ''
            band_data_txt_dw = ''
    return eigenval_dict

def read_procar(procar_file_path):
    '''Read PROCAR'''
    import os
    import numpy as np
    from . import vasp_read
    from . import vasp_tools
    from .. import funcs

    procar_file_path = os.path.abspath(procar_file_path)
    procar_dict = {}
    procar_dict['file_path'] = procar_file_path
    procar_dict['dict_type'] = 'procar'

    file_status = funcs.file_status(procar_file_path)
    procar_dict['file_status'] = file_status

    if file_status != 1:
        print('WARNING #20120311 (from read_procar): The file ' + procar_file_path + ' does not exist or is empty. PLease check the PROCAR file.')
    else:
        #num_header = 3
        num_header = 1

        #two kinds of PROCAR files: 'PROCAR_collinear' or 'PROCAR_noncollinear'
        file_type = vasp_tools.check_file_type(procar_file_path)
        with open(procar_file_path,'r') as f:
            line = f.readlines()
            #num_lines = len(line)
            num_lines = funcs.line_num(procar_file_path)
            num_kpoints = int(funcs.split_line(line = line[1], separator = ':')[1].split('#')[0]) 
            num_bands = int(funcs.split_line(line = line[1], separator = ':')[2].split('#')[0])
            num_ions = int(funcs.split_line(line = line[1], separator = ':')[3].split('#')[0])
            # the PROCAR file has different file format for the collilnear and noncollinear calculations. Care must be taken when dealing with the PROCAR file.
            if file_type == 'PROCAR_collinear':
                #estimated number of lines for PROCAR file in collinear calculations
                estimated_num_lines_ispin1 = num_header + 1 * ((num_kpoints * ((num_bands * (num_ions + 1 + 4)) + 3)) + 1)
                estimated_num_lines_ispin2 = num_header + 2 * ((num_kpoints * ((num_bands * (num_ions + 1 + 4)) + 3)) + 1)
            elif file_type == 'PROCAR_noncollinear':
                #estimated number of lines for PROCAR file in noncollinear calculations
                estimated_num_lines_ispin1 = num_header + 1 * ((num_kpoints * ((num_bands * ((num_ions + 1) * 4 + 4)) + 3)) + 1) 
                estimated_num_lines_ispin2 = num_header + 2 * ((num_kpoints * ((num_bands * ((num_ions + 1) * 4 + 4)) + 3)) + 1)

            if num_lines == estimated_num_lines_ispin1:
                ispin = 1 
            elif num_lines == estimated_num_lines_ispin2:
                ispin = 2
            else:
                print('ERROR: from vasp_read. Please check you PROCAR file')
            num_m = 4 #total partial charge, mx, my and mz (m = magnetization) contributions to that state 
            # orbitals      0     1      2      3     4      5      6      7      8      9
            # orbitals      s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot
            num_orbits = 9

            #to avoid the following error: mask = np.isnan(self.x) TypeError: ufunc 'isnan' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''.    Use dtype = np.float64 instead of dtype = object.
            if ispin == 1:
                eigs = np.array([None] * num_kpoints * num_bands, dtype=np.float64)
                eigs.shape = num_kpoints, num_bands
                occupancy = np.array([None] * num_kpoints * num_bands, dtype=np.float64)
                occupancy.shape = num_kpoints, num_bands
                if file_type == 'PROCAR_collinear':
                    projections = np.array([None] * num_kpoints * num_bands * (num_ions + 1) * (num_orbits + 1), dtype=np.float64)
                    projections.shape = num_kpoints, num_bands, (num_ions + 1), (num_orbits + 1)
                if file_type == 'PROCAR_noncollinear':
                    projections_noncollinear = np.array([None] * num_kpoints * num_bands * (num_ions + 1) * num_m * (num_orbits + 1), dtype=np.float64)
                    projections_noncollinear.shape = num_kpoints, num_bands, (num_ions + 1), num_m, (num_orbits + 1)
            elif ispin == 2:
                eigs_up = np.array([None] * num_kpoints * num_bands, dtype=np.float64)
                eigs_up.shape = num_kpoints, num_bands
                eigs_dw = np.array([None] * num_kpoints * num_bands, dtype=np.float64)
                eigs_dw.shape = num_kpoints, num_bands
                occupancy_up = np.array([None] * num_kpoints * num_bands, dtype=np.float64)
                occupancy_up.shape = num_kpoints, num_bands
                occupancy_dw = np.array([None] * num_kpoints * num_bands, dtype=np.float64)
                occupancy_dw.shape = num_kpoints, num_bands
                if file_type == 'PROCAR_collinear':
                    projections_up = np.array([None] * num_kpoints * num_bands * (num_ions + 1) * (num_orbits + 1), dtype=np.float64)
                    projections_up.shape = num_kpoints, num_bands, (num_ions + 1), (num_orbits + 1)
                    projections_dw = np.array([None] * num_kpoints * num_bands * (num_ions + 1) * (num_orbits + 1), dtype=np.float64)
                    projections_dw.shape = num_kpoints, num_bands, (num_ions + 1), (num_orbits + 1)
                elif file_type == 'PROCAR_noncollinear':
                    projections_up_noncollinear = np.array([None] * num_kpoints * num_bands * (num_ions + 1) * num_m * (num_orbits + 1), dtype=np.float64)
                    projections_up_noncollinear.shape = num_kpoints, num_bands, (num_ions + 1), num_m, (num_orbits + 1)
                    projections_dw_noncollinear = np.array([None] * num_kpoints * num_bands * (num_ions + 1) * num_m, (num_orbits + 1), dtype=np.float64)
                    projections_dw_noncollinear.shape = num_kpoints, num_bands, (num_ions + 1), num_m, (num_orbits + 1)
            
            kpt_coord_arr = np.array([None] * num_kpoints * 3)
            kpt_coord_arr.shape = num_kpoints, 3
            weights = np.array([None] * num_kpoints)

            i_line = num_header
            for i_ispin in range(0, ispin):
                i_line = i_line + 1
                for i_kpoint in range(0, num_kpoints):
                    i_line = i_line + 1
                    ###############################################################################
                    # the format in VASP source code:
                    # sphpro.F: 3201 FORMAT(/' k-point ',I4,' :',3X,3F11.8,'     weight = ',F10.8/)
                    ###############################################################################
                    temp_n = len(' k-point ') + 4 + len(' :') + 3 + 1
                    kpt_coord_arr[i_kpoint, 0] = float(line[i_line][temp_n:temp_n+11])
                    kpt_coord_arr[i_kpoint, 1] = float(line[i_line][temp_n+11:temp_n+11+11])
                    kpt_coord_arr[i_kpoint, 2] = float(line[i_line][temp_n+11+11:temp_n+11+11+11])
                    weights[i_kpoint] = float(funcs.split_line(line = line[i_line], separator = '=')[1])
                    i_line = i_line + 2
                    for i_band in range(0, num_bands):
                        if ispin == 1:
                            eigs[i_kpoint,i_band] = float(funcs.split_line(line = line[i_line], separator = '#')[1].split()[1])
                            occupancy[i_kpoint, i_band] = float(funcs.split_line(line = line[i_line], separator = '#')[2].split()[1])
                        elif ispin == 2 and i_ispin == 0:
                            eigs_up[i_kpoint,i_band] = float(funcs.split_line(line = line[i_line], separator = '#')[1].split()[1])
                            occupancy_up[i_kpoint, i_band] = float(funcs.split_line(line = line[i_line], separator = '#')[2].split()[1])
                        elif ispin == 2 and i_ispin == 1:
                            eigs_dw[i_kpoint,i_band] = float(funcs.split_line(line = line[i_line], separator = '#')[1].split()[1])
                            occupancy_dw[i_kpoint, i_band] = float(funcs.split_line(line = line[i_line], separator = '#')[2].split()[1])
                        i_line = i_line + 3
                        if file_type == 'PROCAR_collinear':
                            for i_ion in range(0, num_ions + 1):
                                for i_orbit in range(0, num_orbits + 1):
                                    if ispin == 1:
                                        projections[i_kpoint, i_band, i_ion, i_orbit] = float(funcs.split_line(line[i_line])[i_orbit + 1])
                                    elif ispin == 2 and i_ispin == 0:
                                        projections_up[i_kpoint, i_band, i_ion, i_orbit] = float(funcs.split_line(line[i_line])[i_orbit + 1])
                                    elif ispin == 2 and i_ispin == 1:
                                        projections_dw[i_kpoint, i_band, i_ion, i_orbit] = float(funcs.split_line(line[i_line])[i_orbit + 1])
                                i_line = i_line + 1
                            i_line = i_line + 1
                        elif file_type == 'PROCAR_noncollinear':
                            for i_m in range(0, num_m):
                                for i_ion in range(0, num_ions + 1):
                                    for i_orbit in range(0, num_orbits + 1):
                                        if ispin == 1:
                                            projections_noncollinear[i_kpoint, i_band, i_ion, i_m, i_orbit] = float(funcs.split_line(line[i_line])[i_orbit + 1])
                                        elif ispin == 2 and i_ispin == 0:
                                            projections_up_noncollinear[i_kpoint, i_band, i_ion, i_m, i_orbit] = float(funcs.split_line(line[i_line])[i_orbit + 1])
                                        elif ispin == 2 and i_ispin == 1:
                                            projections_dw_noncollinear[i_kpoint, i_band, i_ion, i_m, i_orbit] = float(funcs.split_line(line[i_line])[i_orbit + 1])
                                    i_line = i_line + 1
                            i_line = i_line + 1
        #build procar_dict
        procar_dict['file_type'] = file_type
        procar_dict['num_lines'] = num_ions
        procar_dict['num_kpoints'] = num_kpoints
        procar_dict['num_bands'] = num_bands
        procar_dict['num_ions'] = num_ions
        procar_dict['ispin'] = ispin
        procar_dict['kpt_coord_arr'] = kpt_coord_arr
        procar_dict['weights'] = weights 
        if ispin == 1:
            procar_dict['eigs'] = eigs
            procar_dict['occupancy'] = occupancy
            if file_type == 'PROCAR_collinear':
                procar_dict['projections'] = projections
            elif file_type == 'PROCAR_noncollinear':
                procar_dict['projections_noncollinear'] = projections_noncollinear
        elif ispin == 2:
            procar_dict['eigs_up'] = eigs_up
            procar_dict['eigs_dw'] = eigs_dw
            procar_dict['occupancy_up'] = occupancy_up
            procar_dict['occupancy_dw'] = occupancy_dw
            if file_type == 'PROCAR_collinear':
                procar_dict['projections_up'] = projections_up
                procar_dict['projections_dw'] = projections_dw
            elif file_type == 'PROCAR_noncollinear':
                procar_dict['projections_up_noncollinear'] = projections_up_noncollinear
                procar_dict['projections_dw_noncollinear'] = projections_dw_noncollinear
    return procar_dict

def read_incar(incar_file_path):
    '''
    Description:
        Read information from INCAR file
    Args:
        @.incar_file_path: String format. The file path of INCAR. It can either be full path or relative path
    Return:
        incar_dict: dictionary type.
    '''
    import os
    import numpy as np
    import sys
    import time
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])

    incar_file_path = os.path.abspath(incar_file_path)
    incar_dict = {}

    file_status = funcs.file_status(incar_file_path)
    incar_dict['file_status'] = file_status

    if file_status != 1:
        print('WARNING #20120313 (from read_incar): The file ' + incar_file_path + ' does not exist or is empty, please check this file.')
    else:
        with open(incar_file_path) as f:
            line = f.readlines()
            for i_line_indx in range(0, len(line)):
                if line[i_line_indx][0] == '#' or len(line[i_line_indx].strip()) == 0:
                    continue
                else:
                    incar_key = funcs.split_line(line = line[i_line_indx], separator = '=')[0]
                    incar_val = funcs.split_line(line = line[i_line_indx], separator = '=')[1]
                    try:
                        incar_val = float(incar_val)
                    except:
                        pass
                    try:
                        incar_val = int(incar_val)
                    except:
                        pass
                    incar_dict[incar_key] = incar_val
    return incar_dict

def read_potcar(potcar_file_path):
    '''
    Description:
        Read information from POTCAR file
    Args:
        @.potcar_file_path: String format. The file path of POTCAR. It can either be full path or relative path
    Return:
        potcar_dict: dictionary type.
    '''
    import os
    import numpy as np
    import sys
    import time
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])

    potcar_file_path = os.path.abspath(potcar_file_path)
    potcar_dict = {}

    file_status = funcs.file_status(potcar_file_path)
    potcar_dict['file_status'] = file_status

    potcar_dict['num_elmts'] = 0
    potcar_dict['elmt_list'] = []
    potcar_dict['zvalf_list'] = []
    if file_status != 1:
        print('WARNING #2103101917 (from read_potcar): The file ' + potcar_file_path + ' does not exist or is empty, please check this file.')
    else:
        with open(potcar_file_path) as f:
            line = f.readlines()
            eod_indx_list = funcs.find_kwd_line_indx(kwd = 'End of Dataset', file_path = potcar_file_path)
            num_elmts = len(eod_indx_list)
            elmt = funcs.split_line(line = line[0], separator = ' ')[1]
            zvalf = funcs.split_line(line = line[1], separator = ' ')[0]
            potcar_dict['elmt_list'].append(elmt)
            potcar_dict['zvalf_list'].append(zvalf)
            potcar_dict['num_elmts'] = num_elmts
            for i_elmt_indx in range(1, num_elmts):
                elmt = funcs.split_line(line = line[i_elmt_indx + 1], separator = ' ')[1]
                zvalf = funcs.split_line(line = line[i_elmt_indx + 2], separator = ' ')[0]
                potcar_dict['elmt_list'].append(elmt)
                potcar_dict['zvalf_list'].append(zvalf)
    return potcar_dict

##def read_vasprun(vasprun_file_path, screen_list = None):
##    '''
##    Read vasprun.xml file.
##    screen_list = ['incar', 'kpoints', 'parameters', 'structure', 'atominfo', 'calculation']
##    '''
##    import os
##    import xml.etree.ElementTree as ET
##
##    vasprun_dict = {}
##    vasprun_file_path = os.path.abspath(vasprun_file_path)
##    tree = ET.parse('vasprun.xml')
##    root = tree.getroot()
##
##    ##for child in root.iter():
##    ##    print(child.tag, child.attrib)
##    ##    #print(child.tag)
##
##    if screen_list in [None, 'None', 'none']:
##        screen_list = ['incar']
##        pass
##
##    for i_screen_item in screen_list:    
##        for child in root.find('.//' + i_screen_item):
##            print(child.tag, child.attrib)
##
####    for child in root.find('.//incar'):
####        print(child.tag, child.attrib)
####
####    for child in root.find('.//kpoints'):
####        print(child.tag, child.attrib)
####
####    for child in root.find('.//parameters'):
####        print(child.tag, child.attrib)
####
####    for child in root.find('.//atominfo'):
####        print(child.tag, child.attrib)
####
####    # name = initialpos
####    for child in root.find('.//structure'):
####        print(child.tag, child.attrib)
####
####    for child in root.find('.//calculation'):
####        print(child.tag, child.attrib)
####
####    # name = finalpos
####    for child in root.find('.//structure'):
####        print(child.tag, child.attrib)
##
##    print('--------------------------')
##    print('**************************')
##
##
##
##    ##for item in root:
##    ##    d = {}
##    ##    for elmt in item:
##    ##        print(elmt.tag)
##    ##        for k, v in elmt.items():
##    ##            print('    ', k,v) 
##
##
##    return vasprun_dict
