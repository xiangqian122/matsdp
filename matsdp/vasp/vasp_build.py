# -*- coding: utf-8 -*-
def substitution(substitution_list_file, poscar_file_path):
    '''
    Description:
        Generate new POSCAR based on the input POSCAR by alloying element susbstitution
        The substitution_list_file is required and should consists of system entries.
        A system corresponds to a specific model configuration.
        System entris specifies how atoms are substitued in different systems.
        A system entry is a block of sucessive lines without line breaks.
        Each system entris must be seperated by blank lines.
        A typical system entry has the following format:
            Nsubst
            ElementNameToBeSubstituted ElementSubindx NewElementName
            ElementNameToBeSubstituted ElementSubindx NewElementName
            ...
            (n_subst lines of ElementNameToBeSubstituted ElementSubindx NewElementName)
        where,
            ElementNameToBeSubstituted is he name of the element which is to be substituted.
            ElementSubindx is the subindex number of the element which is to be substituted in the POSCAR file;
            NewElementName is the name of the new element which take the place of the substituted atom. If NewElementName = Va, then a vacancy is added.
        As shown above, each system should start with a line which specifies a number: n_subst.
        n_subst is the number of atoms to be substituted in the system.
        Then each of the following n_subst lines specifies the element(s) to be substituted and the element(s) which take its/their place(s).
        A specific example SubstList.in file is as follow:
            1
            Ni244 W
            
            2
            Ni244 Re
            Al12 Re

            ...
        sysname is L$(line_number)_composition_D$(Duplicate)

    - Args
     * substitution_list_file: String format. It specifies the directory of the .subst file (substitution list file)
     * poscar_file_path: String format. The directory of the POSCAR file which is to be subsituted. It can either be full path or relative path.
    '''
    args_dict = locals()
       
    import os
    import sys
    import numpy as np
    from .. import funcs
    from . import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    poscar_file_path = os.path.abspath(poscar_file_path)
    substitution_list_file_abs_path = os.path.abspath(substitution_list_file)
    substitution_list_file_path,substitution_list_file = os.path.split(substitution_list_file_abs_path)
    if os.path.isfile(substitution_list_file_abs_path) == False:
        funcs.write_log(logfile,substitution_list_file + " doesn't exist, please check current folder. The file extension of substitution_list_file should be .subst")
    else:
        subst_file_name, subst_file_extension = os.path.splitext(substitution_list_file)
        if subst_file_extension != '.subst':
            funcs.write_log(logfile,substitution_list_file + " detected. But the file extension should be .subst ")
    if os.path.isfile(poscar_file_path) == False:
        funcs.write_log(logfile, '#' + poscar_file_path + " doesn't exist.")
    
    # Extract information from the input POSCAR file
    poscar_dict = vasp_read.read_poscar(poscar_file_path)
    n_atoms = np.sum(poscar_dict['elmt_num_arr']) 
    with open(substitution_list_file_abs_path) as s_lines:
        s_line = s_lines.readlines();
        s_n_lines = len(s_line)

    #Recognize atom index
    composition_arr = np.array([],dtype = np.str)
    composition_arr_count = np.array([],dtype = np.int)
    sysname_file = os.path.join(substitution_list_file_path, subst_file_name + '.sysname')
    if os.path.isfile(sysname_file):
        os.remove(sysname_file)
    i_line = 1
    model_number = 1
    while (i_line < s_n_lines):
        atom_species_name_subst_arr = poscar_dict['atom_species_arr'].copy()
        coord_subst_arr = poscar_dict['coord_arr'].copy()
        fix_subst_arr = poscar_dict['fix_arr'].copy()
        n_subst = int(funcs.split_line(s_line[i_line-1])[0])
        subst_elmt = np.array([],dtype=np.str)
        for i in range(i_line,i_line+n_subst):
            #Correspond the atom index with the line number
            substituted_elmt_species = ''.join(filter(str.isalpha, funcs.split_line(s_line[i])[0]))  
            substituted_atom_subindx = int(''.join(filter(str.isdigit, funcs.split_line(s_line[i])[0])))
            indx = np.argwhere(poscar_dict['elmt_species_arr']==substituted_elmt_species)
            atom_indx = poscar_dict['elmt_start_indx_arr'][indx] + substituted_atom_subindx - 1
            atom_species_name_subst_arr[atom_indx - 1] = funcs.split_line(s_line[i])[1]
            #Find unique substitution element species
            subst_elmt = np.append(subst_elmt,funcs.split_line(s_line[i])[1])
        unique_subst_elmt = np.unique(subst_elmt)
        elmt_species_mod = np.unique(np.append(poscar_dict['elmt_species_arr'], unique_subst_elmt))
        n_elmt_mod = np.array([0]*len(elmt_species_mod),dtype = np.int)
        for i in range(0,len(elmt_species_mod)):
            n_elmt_mod[i] = sum(atom_species_name_subst_arr==elmt_species_mod[i])
        # remove the element when the number of atoms for an element is zero
        val_indx = np.argwhere(n_elmt_mod == 0)
        n_elmt_mod = np.delete(n_elmt_mod, val_indx)
        elmt_species_mod = np.delete(elmt_species_mod, val_indx) 
        #Define the System Name
        composition = ""
        for i in range(len(elmt_species_mod)-1,-1,-1):
            composition = composition + str(elmt_species_mod[i]) + str(n_elmt_mod[i])
        if composition not in composition_arr:
            composition_arr = np.append(composition_arr,composition)
            composition_arr_count = np.append(composition_arr_count,1)
        else:
            indxtemp = np.argwhere(composition_arr == composition)
            composition_arr_count[indxtemp] = composition_arr_count[indxtemp] + 1
##        int_length = len(str(abs(s_n_lines)))
##        int_format_str = '{:0' + str(int_length) + 'd}'
        int_length = 7
        int_format_str = '{:0' + str(int_length) + 'd}' 
        sysname = 'L' + str(int_format_str.format(i_line)) + '_' + composition + "_D" + str(composition_arr_count[np.argwhere(composition_arr == composition)]).strip('[[').strip(']]')
        atom_species_name_sort_arr = [None] * len(atom_species_name_subst_arr)
        #Rearrange the Atoms according to AtomSpeciesName
        temp1 = 0
        for i in range(0,len(elmt_species_mod)):
            temp2 = 0
            for j in range(0,n_atoms):
                if elmt_species_mod[i] == atom_species_name_subst_arr[j]:
                    coord_subst_arr[temp1][0] = poscar_dict['coord_arr'][j][0]
                    coord_subst_arr[temp1][1] = poscar_dict['coord_arr'][j][1]
                    coord_subst_arr[temp1][2] = poscar_dict['coord_arr'][j][2]
                    fix_subst_arr[temp1][0] = poscar_dict['fix_arr'][j][0]
                    fix_subst_arr[temp1][1] = poscar_dict['fix_arr'][j][1]
                    fix_subst_arr[temp1][2] = poscar_dict['fix_arr'][j][2]
                    atom_species_name_sort_arr[temp1] = atom_species_name_subst_arr[j]
                    temp1 = temp1 + 1
                    temp2 = temp2 + 1
                    if temp2 >= n_elmt_mod[i]:
                        break
        # remove the vacancy
        val_indx = np.argwhere(elmt_species_mod == 'Va')
        elmt_species_mod_remove_va = np.delete(elmt_species_mod, val_indx)
        n_elmt_mod_remove_va = np.delete(n_elmt_mod, val_indx)
        #Export POSCAR file for each system
        models_path = os.path.join(output_dir, subst_file_name, sysname)
        funcs.mkdir(models_path)
        isExists = os.path.exists(models_path)
        if not isExists:
            funcs.write_log(logfile,'WARNING:' + models_path + ' does not exist!')
        poscar_out = models_path + '/POSCAR'
        with open(poscar_out, 'w') as fileobject, open(poscar_file_path,'r') as pline:
            pline = pline.readlines()
            fileobject.write(pline[0] +
                             pline[1] +
                             pline[2] +
                             pline[3] +
                             pline[4] +
                             str(" ".join(elmt_species_mod_remove_va)) + '\n' +
                             str(n_elmt_mod_remove_va).strip('[').strip(']') + '\n')
            if poscar_dict['slet_dyn_on'] == True:
                fileobject.write(pline[7] +
                                 pline[8])
            else:
                fileobject.write(pline[7])
        with open(poscar_out,'a') as fileobject:
            if poscar_dict['slet_dyn_on'] == True:
                for i in range(0, n_atoms):
                    if atom_species_name_sort_arr[i] == 'Va':
                        continue
                    fileobject.write('{:.9f}'.format(coord_subst_arr[i][0]) + ' ' + '{:.9f}'.format(coord_subst_arr[i][1]) + ' ' + '{:.9f}'.format(coord_subst_arr[i][2]) + ' ' + str(fix_subst_arr[i][0]) + ' ' + str(fix_subst_arr[i][1]) + ' ' + str(fix_subst_arr[i][2]) + '\n')
            else:
                for i in range(0, n_atoms):
                    if atom_species_name_sort_arr[i] == 'Va':
                        continue
                    fileobject.write('{:.9f}'.format(coord_subst_arr[i][0]) + ' ' + '{:.9f}'.format(coord_subst_arr[i][1]) + ' ' + '{:.9f}'.format(coord_subst_arr[i][2]) + '\n')
        with open(sysname_file,'a') as sysnamefileobject:
            sysnamefileobject.write(sysname + '\n')
        
        i_line = i_line + n_subst + 2
        model_number += 1
    
    ##################################
    # Determine the args string
    ##################################
    log_str = ''
    func_name = 'vasp_build.substitution'
    args_str = func_name + '(' + '\n'
    for i_arg in args_dict.keys():
        arg_value = args_dict[i_arg]
        if isinstance(arg_value,str):
            arg_value_str = '\'' + arg_value + '\''
        else:
            arg_value_str = str(arg_value)

        if i_arg == 'substitution_list_file':
            arg_value_str = 'r\'' + str(substitution_list_file_abs_path) + '\''
        args_str += '    ' + i_arg + ' = ' + arg_value_str + ',\n'
    args_str += '    )\n'
    args_str += '################################################\n'
    log_str += args_str
    funcs.write_log(logfile, log_str)
    return 0

def rep_elmt(substitution_list_file, poscar_file_path, old_elmt, elmt_group):
    '''
    - Descriptions
     * replace element in the .subst file by specific elements and generate corresponding models
    - Args
     * substitution_list_file: String format. It specifies the directory of the *.subst file (substitution list file).
     * poscar_file_path: String format. The directory of the POSCAR file. It can either be full path or relative path
     * old_elmt: String format. It specifies the element in the *.subst which you want to substute with.
     * elmt_group: list format. It specifies the elements which you want to subsitute for.
    '''
    args_dict = locals()
    import os
    import numpy as np
    import sys
    import time
    from .. import funcs, default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    time_start=time.time()
    poscar_file_path = os.path.abspath(poscar_file_path)
    substitution_list_file_abs_path = os.path.abspath(substitution_list_file)
    subst_file_name, subst_file_extension = os.path.splitext(substitution_list_file)
    sysname_file =subst_file_name + '.sysname'
    str1 = old_elmt
    if os.path.isfile('Tempsysname.dat'):
        os.remove('Tempsysname.dat')
    for i_elmt in elmt_group:
        funcs.replace_file_content(substitution_list_file,str1,i_elmt)
        str1 = i_elmt
        substitution(substitution_list_file, poscar_file_path)
        funcs.merge_files('Tempsysname.dat',sysname_file)
    if os.path.isfile(sysname_file):
        os.remove(sysname_file)
    os.rename('Tempsysname.dat',sysname_file)

    #restore the original SubstList.in
    funcs.replace_file_content(substitution_list_file,str1,old_elmt)
    time_end=time.time()
    funcs.write_log(logfile,'time=' + str(time_end-time_start) + '# unit in second \n')
    ##################################
    # Determine the args string
    ##################################
    log_str = ''
    func_name = 'vasp_build.rep_elmt'
    args_str = func_name + '(' + '\n'
    for i_arg in args_dict.keys():
        arg_value = args_dict[i_arg]
        if isinstance(arg_value,str):
            arg_value_str = '\'' + arg_value + '\''
        else:
            arg_value_str = str(arg_value)

        if i_arg == 'substitution_list_file':
            arg_value_str = 'r\'' + str(substitution_list_file_abs_path) + '\''
        if i_arg == 'poscar_file_path':
            arg_value_str = 'r\'' + str(poscar_file_path) + '\''
        if i_arg == 'old_elmt':
            arg_value_str = '\'' + str(old_elmt) + '\''
        args_str += '    ' + i_arg + ' = ' + arg_value_str + ',\n'
    args_str += '    )\n'
    args_str += '################################################\n'
    log_str += args_str
    funcs.write_log(logfile, log_str)
    return 0

def selection_sphere(poscar_file_path, origin_atom_name, radius = 7.0, include_mirror_atoms = True, output_file_name = 'example'):
    '''
    Descriptions:
    This module is used to select atoms in a sphere that is centered around an arbitrary atom from the POSCAR file.
    The selected atom coordinates can also be written to the *.incar file of the DVM program.
    poscar_file_path: String format. The directory of the POSCAR file. It can either be full path or relative path.
    origin_atom_name: String type. It defines the origin atom of the sphere
    radius: Float type. The atoms within a distance 'radius' from the original atom are selected (units in Angstroms);
    include_mirror_atoms: Logical value. Whether to include the mirror atoms or not;
    output_file_name: user-defined output file name.
    '''
    args_dict = locals()
    import os
    import sys
    import numpy as np
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)
    selection_dir_name = 'atom_selection_sphere'
    funcs.mkdir(os.path.join(output_dir, selection_dir_name, output_file_name))

    poscar_file_path = os.path.abspath(poscar_file_path)

    # the *.vasp file which contains the selected atoms
    vasp_file = os.path.join(output_dir, selection_dir_name, output_file_name, output_file_name + '.vasp')
    vasp_temp_file1 = output_file_name + 'vasp_1.temp'
    vasp_temp_file2 = output_file_name + 'vasp_2.temp'
    # the *.xyz file which contains the selected atoms
    xyz_file = os.path.join(output_dir, selection_dir_name, output_file_name, output_file_name + '.xyz')
    xyz_temp_file1 = output_file_name + '_xyz_1.temp'
    xyz_temp_file2 = output_file_name + '_xyz_2.temp'
    # the *.incar file which contains the selected atoms
    dvm_incar_file = os.path.join(output_dir, selection_dir_name, output_file_name, output_file_name + '.incar')
    dvm_incar_temp_file1 = output_file_name + '_1.temp'
    dvm_incar_temp_file2 = output_file_name + '_2.temp'

    # Extract information from the input POSCAR file
    poscar_dict = vasp_read.read_poscar(poscar_file_path)
    n_atoms = poscar_dict['n_atoms']
    atom_indx = convert.atomname2indx(poscar_file_path, origin_atom_name)
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

    # check whether negative value exist in the atom position
    negative_value_exist = False
    for item in np.sum(np.array(expanded_pos_arr) < 0, axis=0):
        if item > 0:
            negative_value_exist = True

    # write the selected atoms to designated files (*.vasp, *.xyz, *.incar)
    new_atom_num = 0
    new_elmt_species_list = []
    new_elmt_count_list = []
    # write the body of the *.vasp file
    with open(vasp_temp_file1, 'w') as f_vasp_temp1: 
        pass
    with open(vasp_temp_file1, 'a') as f_vasp_temp1: 
        for i_atom in range(0, n_atoms):
            vector1 = origin_atom_pos_arr
            vector2 = np.array([expanded_pos_arr[i_atom, 0], expanded_pos_arr[i_atom, 1], expanded_pos_arr[i_atom, 2]])
            dist = np.linalg.norm(vector1 - vector2)
            if dist <= radius:
                new_atom_num += 1
                if expanded_atom_species_arr[i_atom] not in new_elmt_species_list and expanded_atom_species_arr[i_atom] in poscar_dict['elmt_species_arr']:
                    new_elmt_species_list.append(expanded_atom_species_arr[i_atom])
                    new_elmt_count_list.append(0)
                new_elmt_count_list[new_elmt_species_list.index(expanded_atom_species_arr[i_atom])] += 1
                # *.vasp
                if include_mirror_atoms == True or (include_mirror_atoms == False and negative_value_exist == True):
                    # if the mirror atoms are included, then the atom coordinates will be shifted
                    x_value = expanded_pos_arr[i_atom, 0] - xmin
                    y_value = expanded_pos_arr[i_atom, 1] - ymin
                    z_value = expanded_pos_arr[i_atom, 2] - zmin
                else:
                    x_value = expanded_pos_arr[i_atom, 0]
                    y_value = expanded_pos_arr[i_atom, 1]
                    z_value = expanded_pos_arr[i_atom, 2]
                f_vasp_temp1.write(
                    str('{:.6f}'.format(x_value)) + ' ' +
                    str('{:.6f}'.format(y_value)) + ' ' +
                    str('{:.6f}'.format(z_value)) + ' ' +
                    str(expanded_atom_species_arr[i_atom])+ ' ' +
                    expanded_atom_species_arr[i_atom] + str(new_elmt_count_list[new_elmt_species_list.index(expanded_atom_species_arr[i_atom])]) + ' ' +
                    '    ' + str(expanded_atomname_arr[i_atom]) + ' ' + '\n'
                    )
    # rearrange the order of the elements in the list
    temp_elmt_species_list = []
    temp_elmt_count_list = []
    for i_elmt_species in poscar_dict['elmt_species_arr']:
        if i_elmt_species in new_elmt_species_list:
            temp_elmt_species_list.append(i_elmt_species)
            temp_elmt_count_list.append(new_elmt_count_list[new_elmt_species_list.index(i_elmt_species)])
    new_elmt_species_list = temp_elmt_species_list
    new_elmt_count_list = temp_elmt_count_list
    new_elmt_indx_list = list(range(1, len(new_elmt_species_list) + 1))

    # rearrange the order of the elements in the file
    temp_str = ''
    for i_elmt_species in new_elmt_species_list:
        temp_str = temp_str + funcs.grep('    ' + i_elmt_species, vasp_temp_file1)
    vasp_selected_atom_position = temp_str
    with open(vasp_temp_file1, 'w') as f_vasp_temp1:
        f_vasp_temp1.write(vasp_selected_atom_position)

    # add two columns (element index and atom index) to the file
    temp_str = ''
    temp_file1 = 'temp1.txt'
    with open(temp_file1, 'w') as f_temp1:
        pass
    with open(vasp_temp_file1, 'r') as f_vasp_temp1:
        lines = f_vasp_temp1.readlines()
        with open(temp_file1, 'a') as f_temp1:
            for i_line in range(len(lines)):
                temp_str = temp_str + ' '.join(funcs.split_line(lines[i_line])[0:3]) + ' ' + str(new_elmt_indx_list[new_elmt_species_list.index(funcs.split_line(lines[i_line])[3])]) + ' ' + ' '.join(funcs.split_line(lines[i_line])[3:6]) + ' ' + str(i_line + 1) + '\n'
    vasp_selected_atom_position = temp_str
    os.remove(temp_file1)
    
    with open(vasp_temp_file1, 'w') as f_vasp_temp1:
        f_vasp_temp1.write(vasp_selected_atom_position)
        
    with open(dvm_incar_temp_file1, 'w') as f_dvm_temp1, open(xyz_temp_file1, 'w') as f_xyz_temp1:
        pass
    with open(dvm_incar_temp_file1, 'a') as f_dvm_temp1, open(xyz_temp_file1, 'a') as f_xyz_temp1:
        with open(vasp_temp_file1, 'r') as f_vasp_temp1:
            lines = f_vasp_temp1.readlines()
            for i_line in range(len(lines)):
                f_dvm_temp1.write(' '.join(funcs.split_line(lines[i_line])[0:4]) + '\n')
                f_xyz_temp1.write(str(funcs.split_line(lines[i_line])[4]) + ' ' + ' '.join(funcs.split_line(lines[i_line])[0:4]) + ' ' + ' '.join(funcs.split_line(lines[i_line])[6:8]) + '\n')

    # write the header of the files (*.vasp, *.xyz, *.incar)    
    with open(vasp_temp_file2, 'w') as f_vasp_temp2, open(xyz_temp_file2, 'w') as f_xyz_temp2, open(dvm_incar_temp_file2, 'w') as f_dvm_temp2:
        # *.vasp
        if include_mirror_atoms == True or (include_mirror_atoms == False and negative_value_exist == True):
            l_arr_shifted = poscar_dict['l_arr'].copy()
            l_arr_shifted[0,:] = l_arr_shifted[0,:] * len(range(-nx_lo, nx_hi + 1)) 
            l_arr_shifted[1,:] = l_arr_shifted[1,:] * len(range(-ny_lo, ny_hi + 1)) 
            l_arr_shifted[2,:] = l_arr_shifted[2,:] * len(range(-nz_lo, nz_hi + 1)) 
            f_vasp_temp2.write(
                'Generated by matsdp. The atom coordinates are shifted. added_atom_property=selection: new_elmt_indx elmt_species new_model_atomname original_atomname new_atom_indx\n1\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in l_arr_shifted[0,:]) + '\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in l_arr_shifted[1,:]) + '\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in l_arr_shifted[2,:]) + '\n' +
                ' '.join(new_elmt_species_list) + '\n' +
                ' '.join(str(i) for i in new_elmt_count_list) + '\n' +
                'Cartesian\n' 
                )
        else:
            f_vasp_temp2.write(
                'Generated by matsdp. added_atom_property=selection: new_elmt_indx elmt_species new_model_atomname original_atomname new_atom_indx\n1\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in poscar_dict['l_arr'][0,:]) + '\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in poscar_dict['l_arr'][1,:]) + '\n' +
                ' '.join(str('{:.6f}'.format(i)) for i in poscar_dict['l_arr'][2,:]) + '\n' +
                ' '.join(new_elmt_species_list) + '\n' +
                ' '.join(str(i) for i in new_elmt_count_list) + '\n' +
                'Cartesian\n' 
                )
        # *.xyz
        f_xyz_temp2.write(
            str(new_atom_num) + '\n' +
            '#Generated by matsdp. added_atom_property=selection: new_elmt_indx original_atomname new_atom_indx\n'
            )
        # *.incar
        temp_str = ''
        for i_atom_indx in range(1, new_atom_num + 1):
            temp_str = temp_str + str(i_atom_indx) + ' '
            if i_atom_indx % 10 == 0:
                temp_str = temp_str + '\n'
        formatted_dvm_header = temp_str
        f_dvm_temp2.write(
            str(new_atom_num) + ' 3.0\n' + formatted_dvm_header + '\n\n'
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
    
    ##################################
    # Determine the args string
    ##################################
    log_str = ''
    func_name = 'vasp_build.selection_sphere'
    args_str = func_name + '(' + '\n'
    for i_arg in args_dict.keys():
        arg_value = args_dict[i_arg]
        if isinstance(arg_value,str):
            arg_value_str = '\'' + arg_value + '\''
        else:
            arg_value_str = str(arg_value)

        if i_arg == 'poscar_file_path':
            arg_value_str = 'r\'' + str(poscar_file_path) + '\''
        args_str += '    ' + i_arg + ' = ' + arg_value_str + ',\n'
    args_str += '    )\n'
    args_str += '################################################\n'
    log_str += args_str
    funcs.write_log(logfile, log_str)
    return 0 

def create_multiple_vasp_jobs(substitution_list_file, poscar_file_path, elmt_potcar_dir, incar_str, kpoints_str):
    '''
    create POSCAR, INCAR, KPOINTS, POTCAR files, this module is based on the substitution module
    The POSCAR is created by atom substitution
    '''
    args_dict = locals()
    import os
    from .. import funcs, default_params
    from . import vasp_write

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    substitution(
        substitution_list_file = os.path.abspath(substitution_list_file),
        poscar_file_path = os.path.abspath(poscar_file_path)
        )

    subst_work_dir, substitution_list_file_txt = funcs.file_path_name(substitution_list_file)
    sysname_file_path = os.path.join(subst_work_dir, substitution_list_file_txt[:-6] + '.sysname')
    with open(sysname_file_path, 'r') as f:
        lines = f.readlines()
        for indx in range(0, len(lines)):
            sysname = funcs.split_line(lines[indx])[0]
            if sysname != '':
                work_dir = os.path.join(output_dir, substitution_list_file[:-6], sysname)
                funcs.write_file(
                    input_str = incar_str,
                    dest_file_path = os.path.join(work_dir, 'INCAR')
                    )
                funcs.write_file(
                    input_str = kpoints_str,
                    dest_file_path = os.path.join(work_dir, 'KPOINTS')
                    )
                vasp_write.write_potcar(
                    poscar_path = os.path.join(work_dir, 'POSCAR'),
                    elmt_potcar_dir = os.path.abspath(elmt_potcar_dir)
                    )

    ##################################
    # Determine the args string
    ##################################
    log_str = ''
    func_name = 'vasp_build.create_multiple_vasp_jobs'
    args_str = func_name + '(' + '\n'
    for i_arg in args_dict.keys():
        arg_value = args_dict[i_arg]
        if isinstance(arg_value,str):
            arg_value_str = '\'' + arg_value + '\''
        else:
            arg_value_str = str(arg_value)

        if i_arg == 'substitution_list_file':
            arg_value_str = 'r\'' + str(substitution_list_file) + '\''
        if i_arg == 'poscar_file_path':
            arg_value_str = 'r\'' + str(poscar_file_path) + '\''
        if i_arg == 'elmt_potcar_dir':
            arg_value_str = 'r\'' + str(elmt_potcar_dir) + '\''
        args_str += '    ' + i_arg + ' = ' + arg_value_str + ',\n'
    args_str += '    )\n'
    args_str += '################################################\n'
    log_str += args_str
    funcs.write_log(logfile, log_str)
    return 0

def orientation(poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = None):
    '''
    Rotate the model, to reorientate the a, b, and c axes according to the x, y, and z axes
    The value of mode can be one of the following:
        'c along z, b in yz'   
        'b along z, a in yz'
        'a along z, c in yz'
 
        'a along x, b in xy'
        'c along x, a in xy'
        'b along x, c in xy'
        
        'b along y, a in xy'
        'a along y, c in xy'
        'c along y, b in xy'
        
        'c along z, a in xz'
        'b along z, c in xz'
        'a along z, b in xz'
        
        'a along x, c in xz'
        'c along x, b in xz'
        'b along x, a in xz'
        
        'b along z, c in yz'
        'a along z, b in yz'
        'c along y, a in yz'   
    '''
    args_dict = locals()
    import os
    import numpy as np
    import math
    from .. import funcs
    from . import vasp_read, vasp_write
    from copy import copy, deepcopy
    poscar_file_path = os.path.abspath(poscar_file_path)
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path is None:
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_orientation')
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)

    rot_mat_dict = rot_mat(mode = mode, l_arr = poscar_dict['l_arr'])
    poscar_dict['l_arr'] = rot_mat_dict['rot_mat'].dot(poscar_dict['l_arr'].T).T
    poscar_dict['box_len_arr'] = poscar_dict['l_arr'] / poscar_dict['uni_scale_fac']
    l_arr_temp = deepcopy(poscar_dict['l_arr'])
    box_len_arr_temp = deepcopy(poscar_dict['box_len_arr'])
    # After rotation, we must correctly rearrange the a, b, and c axes. This is very important.
    poscar_dict['l_arr'][rot_mat_dict['xyz_vec_indx_1'], :] = l_arr_temp[rot_mat_dict['abc_vec_indx_1'], :]
    poscar_dict['l_arr'][rot_mat_dict['xyz_vec_indx_2'], :] = l_arr_temp[rot_mat_dict['abc_vec_indx_2'], :]
    poscar_dict['l_arr'][rot_mat_dict['xyz_vec_indx_3'], :] = l_arr_temp[rot_mat_dict['abc_vec_indx_3'], :]
    poscar_dict['box_len_arr'][rot_mat_dict['xyz_vec_indx_1'], :] = box_len_arr_temp[rot_mat_dict['abc_vec_indx_1'], :]
    poscar_dict['box_len_arr'][rot_mat_dict['xyz_vec_indx_2'], :] = box_len_arr_temp[rot_mat_dict['abc_vec_indx_2'], :]
    poscar_dict['box_len_arr'][rot_mat_dict['xyz_vec_indx_3'], :] = box_len_arr_temp[rot_mat_dict['abc_vec_indx_3'], :]
    basis_vector_dict = funcs.basis_vector_info(poscar_dict['l_arr'])

    coord_system = 'Direct' 
    poscar_dict['pos_arr'][:,3:6] = rot_mat_dict['rot_mat'].dot(poscar_dict['pos_arr'][:,3:6].T).T
    if (coord_system.startswith('C') or coord_system.startswith('c')):
        pass
    elif (coord_system.startswith('D') or coord_system.startswith('d')):
        # Use the conventional cell orientation: 'a along x, b in xy'
        # basis_vector_dict['car2fra_matrix_arr'] uses the conventional cell orientation: 'a along x, b in xy
        rot_mat_dict = rot_mat(mode = 'a along x, b in xy', l_arr = poscar_dict['l_arr'])
        poscar_dict['l_arr'] = rot_mat_dict['rot_mat'].dot(poscar_dict['l_arr'].T).T
        poscar_dict['box_len_arr'] = poscar_dict['l_arr'] / poscar_dict['uni_scale_fac']
        l_arr_temp = deepcopy(poscar_dict['l_arr'])
        box_len_arr_temp = deepcopy(poscar_dict['box_len_arr'])
        # After rotation, we must correctly rearrange the a, b, and c axes. This is very important.
        poscar_dict['l_arr'][rot_mat_dict['xyz_vec_indx_1'], :] = l_arr_temp[rot_mat_dict['abc_vec_indx_1'], :]
        poscar_dict['l_arr'][rot_mat_dict['xyz_vec_indx_2'], :] = l_arr_temp[rot_mat_dict['abc_vec_indx_2'], :]
        poscar_dict['l_arr'][rot_mat_dict['xyz_vec_indx_3'], :] = l_arr_temp[rot_mat_dict['abc_vec_indx_3'], :]
        poscar_dict['box_len_arr'][rot_mat_dict['xyz_vec_indx_1'], :] = box_len_arr_temp[rot_mat_dict['abc_vec_indx_1'], :]
        poscar_dict['box_len_arr'][rot_mat_dict['xyz_vec_indx_2'], :] = box_len_arr_temp[rot_mat_dict['abc_vec_indx_2'], :]
        poscar_dict['box_len_arr'][rot_mat_dict['xyz_vec_indx_3'], :] = box_len_arr_temp[rot_mat_dict['abc_vec_indx_3'], :]
        basis_vector_dict = funcs.basis_vector_info(poscar_dict['l_arr'])
    poscar_dict['pos_arr'][:,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][:,3:6].T).T
    #print('before=',poscar_dict['pos_arr'][:,0:3])
    # adding 1 to np.modf(poscar_dict['pos_arr'][:,0:3])[0] will cause a problem of large number adding small number, this is only an approximate way
    poscar_dict['pos_arr'][:,0:3] = np.modf(np.modf(poscar_dict['pos_arr'][:,0:3])[0] + 1)[0]
    # Replace the value close to one with 0
    poscar_dict['pos_arr'][:,0:3][np.isclose(poscar_dict['pos_arr'][:,0:3], 1, atol = 1e-16)] = 0
    #print('after=',poscar_dict['pos_arr'][:,0:3])
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = coord_system)
    return 0

def rot_mat(mode, l_arr):
    '''
    get rotation matrix from the mode variable
    l_arr: [vec_a, vec_b, vec_c]

    Rotate the model, to reorientate the a, b, and c axes according to the x, y, and z axes
    The value of mode can be one of the following:
        'c along z, b in yz'   
        'b along z, a in yz'
        'a along z, c in yz'
 
        'a along x, b in xy'
        'c along x, a in xy'
        'b along x, c in xy'
        
        'b along y, a in xy'
        'a along y, c in xy'
        'c along y, b in xy'
        
        'c along z, a in xz'
        'b along z, c in xz'
        'a along z, b in xz'
        
        'a along x, c in xz'
        'c along x, b in xz'
        'b along x, a in xz'
        
        'b along z, c in yz'
        'a along z, b in yz'
        'c along y, a in yz'   

    Method:
        for m along x/y/z, m in xy/yz/xz
        first, find the vector p that is perpendicular to both m and n, then align p to the normal of xy/yz/xz.
        second, align the m vector to x/y/z
    '''
    args_dict = locals()
    import numpy as np
    from .. import funcs
    rot_mat_dict = {}
    rot_mat_dict['rot_mat'] = None
    rot_mat_dict['abc_vec_indx_1'] = None
    rot_mat_dict['abc_vec_indx_2'] = None
    rot_mat_dict['abc_vec_indx_3'] = None
    rot_mat_dict['xyz_vec_indx_1'] = None
    rot_mat_dict['xyz_vec_indx_2'] = None
    rot_mat_dict['xyz_vec_indx_3'] = None
    rot_mat_dict['xyz_plane_indx'] = None

    x_vector = np.array([1, 0, 0])    
    y_vector = np.array([0, 1, 0])    
    z_vector = np.array([0, 0, 1])    
    # abc_vec_1 along xyz_vec_1, abc_vec_2 in xyz_plane
    abc_vec_indx_list = [0, 1, 2]
    xyz_vec_indx_list = [0, 1, 2]

    abc_vec_str_1 = mode.split('along')[0].strip()
    if abc_vec_str_1 == 'a':
        abc_vec_indx_1 = 0
        abc_vec_1 = l_arr[abc_vec_indx_1, :]  
    elif abc_vec_str_1 == 'b':
        abc_vec_indx_1 = 1  
        abc_vec_1 = l_arr[abc_vec_indx_1, :]  
    elif abc_vec_str_1 == 'c':
        abc_vec_indx_1 = 2  
        abc_vec_1 = l_arr[abc_vec_indx_1, :]  

    abc_vec_str_2 = mode.split(',')[1].split('in')[0].strip()
    if abc_vec_str_2 == 'a':
        abc_vec_indx_2 = 0
        abc_vec_2 = l_arr[abc_vec_indx_2, :]  
    elif abc_vec_str_2 == 'b':
        abc_vec_indx_2 = 1  
        abc_vec_2 = l_arr[abc_vec_indx_2, :]  
    elif abc_vec_str_2 == 'c':
        abc_vec_indx_2 = 2  
        abc_vec_2 = l_arr[abc_vec_indx_2, :]  

    ab_vec_indx_list = [abc_vec_indx_1, abc_vec_indx_2]
    abc_vec_indx_3 = list(set(abc_vec_indx_list) - set(ab_vec_indx_list))[0]

    xyz_vec_str_1 = mode.split(',')[0].split()[-1].strip()
    if xyz_vec_str_1 == 'x':
        xyz_vec_indx_1 = 0
        xyz_vec_1 = x_vector
    elif xyz_vec_str_1 == 'y':
        xyz_vec_indx_1 = 1
        xyz_vec_1 = y_vector
    elif xyz_vec_str_1 == 'z':
        xyz_vec_indx_1 = 2
        xyz_vec_1 = z_vector

    xyz_plane = mode.split(',')[1].split()[-1].strip()
    if xyz_plane in ['xy', 'yx']:
        xyz_plane_indx = 2
        xyz_plane_vec = z_vector
    elif xyz_plane in ['yz', 'zy']:
        xyz_plane_indx = 0
        xyz_plane_vec = x_vector
    elif xyz_plane in ['xz' or 'zx']:
        xyz_plane_indx = 1
        xyz_plane_vec = y_vector
 
    xyz_vec_indx_temp_list = [xyz_vec_indx_1, xyz_plane_indx]
    xyz_vec_indx_2 = list(set(xyz_vec_indx_list) - set(xyz_vec_indx_temp_list))[0]
    xyz_vec_indx_3 = xyz_plane_indx
    
    if [abc_vec_indx_1, abc_vec_indx_2] in [[0, 1], [1, 2], [2, 0]]:
        abc_vec1_cross_vec2 = np.cross(abc_vec_1, abc_vec_2)
    elif [abc_vec_indx_1, abc_vec_indx_2] in [[1, 0], [2, 1], [0, 2]]:
        abc_vec1_cross_vec2 = np.cross(abc_vec_2, abc_vec_1)

    abc_vec_1 = (abc_vec_1 / np.linalg.norm(abc_vec_1)).reshape((3, 1))
    xyz_vec_1 = (xyz_vec_1 / np.linalg.norm(xyz_vec_1)).reshape((3, 1))
    abc_vec1_cross_vec2 = (abc_vec1_cross_vec2 / np.linalg.norm(abc_vec1_cross_vec2)).reshape((3, 1))
    xyz_plane_vec = (xyz_plane_vec / np.linalg.norm(xyz_plane_vec)).reshape((3, 1))

    ##matrix_arr_1 = funcs.rotation_matrix_from_vectors(abc_vec_1, xyz_vec_1)
    ##abc_vec1_cross_vec2 = matrix_arr_1.dot(abc_vec1_cross_vec2)
    ##matrix_arr_2 = funcs.rotation_matrix_from_vectors(abc_vec1_cross_vec2, xyz_plane_vec)
    ##rotation_matrix_arr = matrix_arr_2.dot(matrix_arr_1)

    matrix_arr_1 = funcs.rotation_matrix_from_vectors(abc_vec1_cross_vec2, xyz_plane_vec)
    abc_vec_1 = matrix_arr_1.dot(abc_vec_1)
    matrix_arr_2 = funcs.rotation_matrix_from_vectors(abc_vec_1, xyz_vec_1)
    rotation_matrix_arr = matrix_arr_2.dot(matrix_arr_1)

    rot_mat_dict['rot_mat'] = rotation_matrix_arr
    rot_mat_dict['abc_vec_indx_1'] = abc_vec_indx_1 
    rot_mat_dict['abc_vec_indx_2'] = abc_vec_indx_2
    rot_mat_dict['abc_vec_indx_3'] = abc_vec_indx_3
    rot_mat_dict['xyz_vec_indx_1'] = xyz_vec_indx_1
    rot_mat_dict['xyz_vec_indx_2'] = xyz_vec_indx_2
    rot_mat_dict['xyz_vec_indx_3'] = xyz_vec_indx_3
    rot_mat_dict['xyz_plane_indx'] = xyz_plane_indx
    return rot_mat_dict

def align(poscar_file_path, align_direction = 'z', mode = 'bottom', output_poscar_file_path = None, delta = 0.05, radius_scale_factor = 1.15, write_layer_info = False, suppress_warning = True):
    '''
    mode = 'bottom' or 'center' or 'upright_box'
    mode = 'bottom' or 'center': this is only applied to atoms. Align the atoms to the bottom or center of the box in a certain direction
    mode = 'upright_box': Make the tilted box an "upright" one. Align one of the basis vector (a, b, or c) to x, y, or z direction(whether its x, y, or z is desiginated by the align_direction tag), without changing the box volume. First rotate the model, to reorientate the a, b, and c axes according to the x, y, and z axes (align one of the a, b, or c to one of the x, y or z axis, then make two of the a, b, and c axes in the same plane of two of the x, y, and z axes), then make two of the a, b, and c axes penpendicular to each other.
    '''
    args_dict = locals()
    import os
    import numpy as np
    import math
    from .. import funcs
    from . import vasp_read, vasp_write, vasp_tools
    poscar_file_path = os.path.abspath(poscar_file_path)
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path in [None, 'None', 'none']:
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_align')
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)
    poscar_direction_dict = vasp_tools.poscar_direction_params(poscar_dict = poscar_dict, direction = align_direction)

    # if vacuum exist in the model, align the atoms to the bottom of the model to avoid atoms too close to the top boundary
    # Check the existence of the vacuum layer and readjust the position of the layer with the value of "-1" in the bond_status_list to the origin.
    if mode in ['bottom', 'center']:
        poscar_layer_dict = vasp_tools.poscar_layer_params(poscar_dict = poscar_dict, poscar_direction_dict = poscar_direction_dict, criteria = 'auto', delta = delta, layer_dist_tolerance = 'auto', radius_style = 'csd', radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning = suppress_warning)
        if -1 in poscar_layer_dict['bond_status_list']:
            i_layer_indx = int(np.argwhere(np.array(poscar_layer_dict['bond_status_list']) == -1)[0])
            atom_pos_list = []
            for i_atom_indx in poscar_layer_dict['atoms_indx_list'][i_layer_indx]:
                atom_pos_list.append(poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_direct_column']])
            close_to_one_list = np.isclose(atom_pos_list, 1, atol = delta / 5)
            if True in close_to_one_list: 
                pos_i = np.modf(poscar_dict['pos_arr'][:,poscar_direction_dict['pos_arr_direct_column']] + 1 - np.min(np.array(atom_pos_list)[close_to_one_list]))[0]
                poscar_dict['pos_arr'][:,poscar_direction_dict['pos_arr_direct_column']] = pos_i
                poscar_dict['pos_arr'][:,3:6] = np.dot(poscar_dict['pos_arr'][:,0:3], poscar_dict['l_arr'])
            else:
                pos_i = poscar_dict['pos_arr'][:,poscar_direction_dict['pos_arr_direct_column']] - np.min(atom_pos_list)
                base_i = np.floor(poscar_dict['pos_arr'][:,poscar_direction_dict['pos_arr_direct_column']] - np.min(atom_pos_list))
                poscar_dict['pos_arr'][:,poscar_direction_dict['pos_arr_direct_column']] = pos_i - base_i
                poscar_dict['pos_arr'][:,3:6] = np.dot(poscar_dict['pos_arr'][:,0:3], poscar_dict['l_arr'])
            vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')

    if mode == 'bottom':
        # align bottom
        lowest_atom_indx = int(np.argwhere(poscar_dict['pos_arr'][:, poscar_direction_dict['pos_arr_column']] == np.min(poscar_dict['pos_arr'][:, poscar_direction_dict['pos_arr_column']]))[0])
        displace_vector = poscar_dict['pos_arr'][lowest_atom_indx, 3:6].dot(poscar_direction_dict['ortho_vector']) / poscar_direction_dict['cos_angle'] * (-poscar_direction_dict['unit_vector'])
        poscar_dict['pos_arr'][:, 3:6] = poscar_dict['pos_arr'][:, 3:6] + displace_vector
    elif mode == 'center':
        # align bottom to the bottom first
        lowest_atom_indx = int(np.argwhere(poscar_dict['pos_arr'][:, poscar_direction_dict['pos_arr_column']] == np.min(poscar_dict['pos_arr'][:, poscar_direction_dict['pos_arr_column']]))[0])
        displace_vector = poscar_dict['pos_arr'][lowest_atom_indx, 3:6].dot(poscar_direction_dict['ortho_vector']) / poscar_direction_dict['cos_angle'] * (-poscar_direction_dict['unit_vector'])
        # Displace the atoms to the center of the axis
        crystal_len_ortho = np.max(poscar_dict['pos_arr'][:, poscar_direction_dict['pos_arr_column']]) - np.min(poscar_dict['pos_arr'][:, poscar_direction_dict['pos_arr_column']])
        displace_vector = displace_vector + 1 / 2 * (poscar_direction_dict['box_len_ortho'] - crystal_len_ortho) / poscar_direction_dict['cos_angle'] * poscar_direction_dict['unit_vector']
        poscar_dict['pos_arr'][:, 3:6] = poscar_dict['pos_arr'][:, 3:6] + displace_vector
    elif mode == 'upright_box':
        if align_direction == 'z':
            orientation_mode = 'a along x, b in xy' 
        if align_direction == 'y':
            orientation_mode = 'a along x, c in xz' 
        if align_direction == 'x':
            orientation_mode = 'b along y, c in yz'
        # reorientate the a, b, and c axes according to x, y, and z axes.
        orientation(poscar_file_path, mode = orientation_mode, output_poscar_file_path = output_poscar_file_path)
        poscar_dict = vasp_read.read_poscar(output_poscar_file_path)
        # get basic information
        vec_1 = poscar_direction_dict['vec_1']
        vec_2 = poscar_direction_dict['vec_2']
        vec_12 = vec_1 + vec_2
        poscar_direction_dict = vasp_tools.poscar_direction_params(poscar_dict = poscar_dict, direction = align_direction)
        projected_side_vector_len_1 = vec_1.dot(poscar_direction_dict['ortho_vector_1'])
        projected_side_vector_len_2 = vec_2.dot(poscar_direction_dict['ortho_vector_2'])
        projected_side_vector_1 = projected_side_vector_len_1 * poscar_direction_dict['ortho_vector_1']
        projected_side_vector_2 = projected_side_vector_len_2 * poscar_direction_dict['ortho_vector_2']
        projected_side_vector_len_12_1 = vec_12.dot(poscar_direction_dict['ortho_vector_1'])
        projected_side_vector_len_12_2 = vec_12.dot(poscar_direction_dict['ortho_vector_2'])
        area_parallelogram = np.linalg.norm(np.cross(vec_1, vec_2))
        # make the designated axis an upright one
        poscar_dict['box_len_arr'][poscar_direction_dict['l_arr_row'], poscar_direction_dict['l_arr_column_1']] = 0 
        poscar_dict['box_len_arr'][poscar_direction_dict['l_arr_row'], poscar_direction_dict['l_arr_column_2']] = 0
        poscar_dict['l_arr'][poscar_direction_dict['l_arr_row'], poscar_direction_dict['l_arr_column_1']] = 0 
        poscar_dict['l_arr'][poscar_direction_dict['l_arr_row'], poscar_direction_dict['l_arr_column_2']] = 0
        #Whether a point is inside the parallelogram composed of vec_1 and vec_2. point_1, point_2, point_3, point_4 are the four vertices of the parallelogram composed of vec_1 and vec_2. If the area of four triangles equals to the area of the parallelogram, then the point is inside the parallelogram..
        point_1 = np.array([0, 0])
        point_2 = np.array([vec_1[poscar_direction_dict['l_arr_column_1']], vec_1[poscar_direction_dict['l_arr_column_2']]])
        point_3 = np.array([vec_2[poscar_direction_dict['l_arr_column_1']], vec_2[poscar_direction_dict['l_arr_column_2']]])
        point_4 = np.array([vec_12[poscar_direction_dict['l_arr_column_1']], vec_12[poscar_direction_dict['l_arr_column_2']]])
        delta_area = 0.01
        test_projection = False
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            temp_pos_1 = poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_column_1']]
            temp_pos_2 = poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_column_2']]
            point_0 = np.array([temp_pos_1, temp_pos_2])
            area_four_triangle = (
                funcs.area_triangle_from_points(point_0, point_1, point_2) + 
                funcs.area_triangle_from_points(point_0, point_2, point_4) + 
                funcs.area_triangle_from_points(point_0, point_4, point_3) + 
                funcs.area_triangle_from_points(point_0, point_3, point_1))
            if test_projection == True:
                print(area_four_triangle)
                print(area_parallelogram)
                print('----------------------')
                import matplotlib.pyplot as plt
                if i_atom_indx != -1:
                    plt.plot(point_1[0], point_1[1], linestyle = '', marker = '.', markersize = 10, color = 'black')
                    plt.plot(point_2[0], point_2[1], linestyle = '', marker = '.', markersize = 10, color = 'black')
                    plt.plot(point_3[0], point_3[1], linestyle = '', marker = '.', markersize = 10, color = 'black')
                    plt.plot(point_4[0], point_4[1], linestyle = '', marker = '.', markersize = 10, color = 'black')
                    plt.plot(point_0[0], point_0[1], linestyle = '', marker = 's', markersize = 10, color = 'red')
                    plt.text(point_0[0], point_0[1], '{:.2f}'.format(area_four_triangle))
            # for points not in the parallelogram, translate them to a position inside the parallelogram.
            n1 = 1
            n2 = 1
            max_dim = 7
            n1_arr = np.array([None] * max_dim * 2)
            i = 0
            for i_indx in range(0, max_dim + 1):
                n1_arr[i] = i_indx
                if i_indx != max_dim:
                    n1_arr[i + 1] = -i_indx
                if i == 0:
                    i = i + 1
                else:
                    i = i + 2
            n2_arr = n1_arr.copy()
            if abs(area_four_triangle - area_parallelogram) > delta_area:
                shift_arr_1 = np.array([vec_1[poscar_direction_dict['l_arr_column_1']], vec_1[poscar_direction_dict['l_arr_column_2']]])
                shift_arr_2 = np.array([vec_2[poscar_direction_dict['l_arr_column_1']], vec_2[poscar_direction_dict['l_arr_column_2']]])
                for n1 in n1_arr:
                    for n2 in n2_arr:
                        if n1 == 0 and n2 == 0:
                            continue
                        point_1p = point_1 + n1 * shift_arr_1 + n2 * shift_arr_2  
                        point_2p = point_2 + n1 * shift_arr_1 + n2 * shift_arr_2
                        point_3p = point_3 + n1 * shift_arr_1 + n2 * shift_arr_2
                        point_4p = point_4 + n1 * shift_arr_1 + n2 * shift_arr_2
                        area_four_triangle_1 = (
                            funcs.area_triangle_from_points(point_0, point_1p, point_2p) + 
                            funcs.area_triangle_from_points(point_0, point_2p, point_4p) + 
                            funcs.area_triangle_from_points(point_0, point_4p, point_3p) + 
                            funcs.area_triangle_from_points(point_0, point_3p, point_1p))
                        if abs(area_four_triangle_1 - area_parallelogram) <= delta_area:
                            poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_column_1']] =  temp_pos_1 - n1 * shift_arr_1[0] - n2 * shift_arr_2[0]
                            poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_column_2']] =  temp_pos_2 - n1 * shift_arr_1[1] - n2 * shift_arr_2[1]
                        if test_projection == True:
                            if i_atom_indx != -1:
                                plt.plot(point_1p[0], point_1p[1], linestyle = '', marker = '.', markersize = 10, color = 'green')
                                plt.plot(point_2p[0], point_2p[1], linestyle = '', marker = '.', markersize = 10, color = 'green')
                                plt.plot(point_3p[0], point_3p[1], linestyle = '', marker = '.', markersize = 10, color = 'green')
                                plt.plot(point_4p[0], point_4p[1], linestyle = '', marker = '.', markersize = 10, color = 'green')
                                plt.plot(point_0[0], point_0[1], linestyle = '', marker = 's', markersize = 10, color = 'blue')
                            break
        if test_projection == True:
            plt.show()

    vec_a = poscar_dict['box_len_arr'][0, :] * poscar_dict['uni_scale_fac'] 
    vec_b = poscar_dict['box_len_arr'][1, :] * poscar_dict['uni_scale_fac']
    vec_c = poscar_dict['box_len_arr'][2, :] * poscar_dict['uni_scale_fac']
    cell_arr = np.vstack((vec_a, vec_b, vec_c))
    basis_vector_dict = funcs.basis_vector_info(cell_arr)
    poscar_dict['pos_arr'][:,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][:,3:6].T).T
    poscar_dict['pos_arr'][:,0:3] = np.modf(np.modf(poscar_dict['pos_arr'][:,0:3])[0] + 1)[0]
    poscar_dict['pos_arr'][:,0:3][np.isclose(poscar_dict['pos_arr'][:,0:3], 1, atol = 1e-16)] = 0
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    return 0

def transform(poscar_file_path, r_matrix_arr = None, t_matrix_arr = None, apply_to = 'atom', output_poscar_file_path = None):
    '''
    Transform by rotation and transformation operations.
    r_matrix_arr: Rotation matrix
    t_matrix_arr: translation matrix
    apply_to: 'box/coordinate_system/basis' or 'atom/object' or 'all/atom_and_box'
    If shear deformatoin is performed, please set apply_to == 'all'
    '''
    args_dict = locals()
    import os
    import numpy as np
    from .. import funcs
    from . import vasp_read, vasp_write
    poscar_file_path = os.path.abspath(poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)
    new_extension = '_transform.vasp'
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path is None:
        fpre, fext = os.path.splitext(poscar_file_path)
        output_poscar_file_path = fpre + new_extension
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    if r_matrix_arr is None:
        r_matrix_arr = np.identity(3)
    if t_matrix_arr is None:
        t_matrix_arr = np.zeros((3, 3))
    if isinstance(r_matrix_arr, list):
        r_matrix_arr = np.array(r_matrix_arr)
    if isinstance(t_matrix_arr, list):
        t_matrix_arr = np.array(t_matrix_arr)
    if apply_to in ['atom', 'object']:
        atom_coord_vec = poscar_dict['pos_arr'][:,3:6]
        atom_coord_vec_trans = r_matrix_arr.dot(atom_coord_vec.T).T + np.tile(np.sum(t_matrix_arr, axis = 0), (atom_coord_vec.shape[0],1))
        poscar_dict['pos_arr'][:,3:6] = atom_coord_vec_trans[:,0:3]
        poscar_dict['pos_arr'][:,0:3] = poscar_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][:,3:6].T).T
        poscar_dict['pos_arr'][:,0:3] = np.modf(np.modf(poscar_dict['pos_arr'][:,0:3])[0] + 1)[0]
        poscar_dict['pos_arr'][:,0:3][np.isclose(poscar_dict['pos_arr'][:,0:3], 1, atol = 1e-16)] = 0
        coord_system = 'Direct'
    elif apply_to in ['box', 'coordinate_system', 'basis']:
        poscar_dict['box_len_arr'] = r_matrix_arr.dot(poscar_dict['box_len_arr'].T).T + t_matrix_arr
        cell_arr = poscar_dict['box_len_arr'] * poscar_dict['uni_scale_fac']
        basis_vector_dict = funcs.basis_vector_info(cell_arr)
        poscar_dict['pos_arr'][:,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][:,3:6].T).T
        poscar_dict['pos_arr'][:,0:3] = np.modf(np.modf(poscar_dict['pos_arr'][:,0:3])[0] + 1)[0]
        poscar_dict['pos_arr'][:,0:3][np.isclose(poscar_dict['pos_arr'][:,0:3], 1, atol = 1e-16)] = 0
        # Direct coordinate will use the 'a along x, b in xy' convention, so it is recommended to write in the Cartesian coordinate
        coord_system = 'Cartesian'
    elif apply_to in ['all', 'atom_and_box']:
        # apply to box
        poscar_dict['box_len_arr'] = r_matrix_arr.dot(poscar_dict['box_len_arr'].T).T + t_matrix_arr
        # apply to atoms
        atom_coord_vec = poscar_dict['pos_arr'][:,3:6]
        atom_coord_vec_trans = r_matrix_arr.dot(atom_coord_vec.T).T + np.tile(np.sum(t_matrix_arr, axis = 0), (atom_coord_vec.shape[0],1))
        poscar_dict['pos_arr'][:,3:6] = atom_coord_vec_trans[:,0:3]
        # Direct coordinate will use the 'a along x, b in xy' convention, so it is recommended to write in the Cartesian coordinate
        coord_system = 'Cartesian'
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = coord_system)
    return 0

def create_supercell(poscar_file_path, num_xyz_list = [1,1,1], output_poscar_file_path = None):
    '''
    Create supercell by expanding the cell in the direction of a, b, and c basis vectors.
    num_xyz_list = [n_x, n_y, n_z], nx>=1, ny>=1, nz>=1
    '''    
    args_dict = locals()
    import os
    import numpy as np
    from .. import funcs
    from . import vasp_read, vasp_write
    poscar_file_path = os.path.abspath(poscar_file_path)
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path in [None, 'None', 'none']:
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_supercell')
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)

    nx = num_xyz_list[0]
    ny = num_xyz_list[1]
    nz = num_xyz_list[2]
    if nx <= 0 or ny <=0 or nz <=0:
        print('WARNING #2012291940 (from create_supercell): The setting of num_xyz_list is not correct. Please check num_xyz_list. num_xyz_list = [n_x, n_y, n_z], nx>=1, ny>=1, nz>=1. num_xyz_list is set to num_xyz_list = [1,1,1].')
        nx = 1
        ny = 1
        nz = 1

    orientation(poscar_file_path = poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(output_poscar_file_path)
    #poscar_dict = vasp_read.read_poscar(poscar_file_path)
    poscar_dict['n_atoms'] = poscar_dict['n_atoms'] * nx * ny * nz
    # manipulation of box
    n_species = len(poscar_dict['elmt_species_arr'])
    pos_arr_old = np.copy(poscar_dict['pos_arr'])
    fix_arr_old = np.copy(poscar_dict['fix_arr'])
    elmt_num_arr_old = np.copy(poscar_dict['elmt_num_arr'])
    elmt_start_indx_arr_old = np.array([0]*n_species,dtype=np.int)

    poscar_dict['box_len_arr'][0,:] = poscar_dict['box_len_arr'][0,:] * nx
    poscar_dict['box_len_arr'][1,:] = poscar_dict['box_len_arr'][1,:] * ny
    poscar_dict['box_len_arr'][2,:] = poscar_dict['box_len_arr'][2,:] * nz
    poscar_dict['elmt_num_arr'] = poscar_dict['elmt_num_arr'] * nx * ny * nz
    poscar_dict['l_arr'] = poscar_dict['box_len_arr'] * poscar_dict['uni_scale_fac']
    basis_vector_dict = funcs.basis_vector_info(poscar_dict['l_arr'])

    #Element start index
    poscar_dict['elmt_start_indx_arr'] = np.array([0]*n_species,dtype=np.int)
    for i in range(0, n_species):
        if i == 0:
            elmt_start_indx_arr_old[i] = 1
            poscar_dict['elmt_start_indx_arr'][i] = 1
        else:
            elmt_start_indx_arr_old[i] = sum(elmt_num_arr_old[range(0,i)]) + 1
            poscar_dict['elmt_start_indx_arr'][i] = sum(poscar_dict['elmt_num_arr'][range(0,i)]) + 1
    # manipulation of atoms
    start_indx = 0
    poscar_dict['pos_arr'] = np.array([0.0] * poscar_dict['n_atoms'] * 6,dtype = np.float)
    poscar_dict['pos_arr'].shape = poscar_dict['n_atoms'], 6
    poscar_dict['fix_arr'] = np.array(['T'] * poscar_dict['n_atoms'] * 3,dtype = np.str)
    poscar_dict['fix_arr'].shape = poscar_dict['n_atoms'], 3
    for i_elmt_species in range(0, len(poscar_dict['elmt_species_arr'])):
        start_indx = elmt_start_indx_arr_old[i_elmt_species] - 1
        end_indx = start_indx + elmt_num_arr_old[i_elmt_species]
        row_span = slice(start_indx, end_indx)

        start_indx_new = poscar_dict['elmt_start_indx_arr'][i_elmt_species] - 1  
        end_indx_new = sum(poscar_dict['elmt_num_arr'][range(0,i_elmt_species + 1)])
        #end_indx_new = start_indx + poscar_dict['elmt_num_arr'][i_elmt_species] * nx * ny * nz
        row_span_new = slice(start_indx_new, end_indx_new)
        temp_arr_coord = np.copy(pos_arr_old[row_span, 3:6])
        temp_arr_fix = np.copy(fix_arr_old[row_span, :])
        #temp_arr_coord_old = np.copy(pos_arr_old[row_span, 3:6])
        for n1 in range(0, nx):
            for n2 in range(0, ny):
                for n3 in range(0, nz):
                    if n1 == 0 and n2 == 0 and n3== 0:
                        continue
                    new_coord_arr = pos_arr_old[row_span, 3:6] + poscar_dict['vec_a'] * n1 + poscar_dict['vec_b'] * n2 + poscar_dict['vec_c'] * n3
                    temp_arr_coord = np.vstack((temp_arr_coord, new_coord_arr))

                    #temp_arr_fix = np.copy(poscar_dict['fix_arr'][row_span,:])
                    new_fix_arr = fix_arr_old[row_span, :]
                    temp_arr_fix = np.vstack((temp_arr_fix, new_fix_arr))
        poscar_dict['pos_arr'][row_span_new, 3:6] = temp_arr_coord
        poscar_dict['fix_arr'][row_span_new, :] = temp_arr_fix
    poscar_dict['pos_arr'][:,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][:,3:6].T).T
    poscar_dict['pos_arr'][:,0:3] = np.modf(np.modf(poscar_dict['pos_arr'][:,0:3])[0] + 1)[0]
    poscar_dict['pos_arr'][:,0:3][np.isclose(poscar_dict['pos_arr'][:,0:3], 1, atol = 1e-16)] = 0
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    return 0

def build_vacuum_layer(poscar_file_path, vacuum_layer_direction = 'z', vacuum_width = 20, threshold_vacuum_width = 5.98, output_poscar_file_path = None, delta = 0.05, radius_scale_factor = 1.15, write_layer_info = False, suppress_warning = True):
    '''
    Build vacuum layer
    NOTE: In principle, the normal of the layers should be parallel with vacuum_layer_direction.
    Please note that after building the vacuum layer, the crystal will be shifted to the bottom of the model.
    The vacuum_width is defined as the difference between the maximum and minimum porjected(onto x, y, or z axis) atom coordinates.
    vacuum_layer_direction: Designate the direction to add the vacuum layer
    threshold_vacuum_width: This is the threshold distance for the vacuum layer separation, if layer separation is larger than this value, then it is denoted as vacuum layer. Otherwise, it is not a vacuum layer. Default value of threshold_vacuum_width is (1.15 * 2 * (CSD_covalent_radius of Fr)) = 1.15 * 2 * 2.60 = 5.98 Angstrom
    '''
    args_dict = locals()
    import os
    import math
    import numpy as np
    from .. import funcs
    from . import vasp_read, vasp_write, vasp_tools

    poscar_file_path = os.path.abspath(poscar_file_path)
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path in [None, 'None', 'none']:
        #output_poscar_file_path = os.path.join(fpath, 'POSCAR_vacuum_' + str(vacuum_width))
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_vacuum')
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    #orientation(poscar_file_path = poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = output_poscar_file_path)
    #poscar_dict = vasp_read.read_poscar(output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)

    if vacuum_width <= threshold_vacuum_width:
        print('WARNING #2012301433 (from build_vacuum_layer): The vacuum layer width is too small, it is recommended to set it to a value larger than ' + str(threshold_vacuum_width) + ' Angstrom.')

    ##############################################################
    # initialization
    ##############################################################
    poscar_direction_dict = vasp_tools.poscar_direction_params(poscar_dict = poscar_dict, direction = vacuum_layer_direction)
    ##if vacuum_layer_direction in ['x', 'X']:
    ##    nx = math.ceil((vacuum_width + poscar_direction_dict['box_len_ortho']) / poscar_direction_dict['box_len_ortho'])
    ##    ny = 1
    ##    nz = 1
    ##if vacuum_layer_direction in ['y', 'Y']:
    ##    nx = 1
    ##    ny = math.ceil((vacuum_width + poscar_direction_dict['box_len_ortho']) / poscar_direction_dict['box_len_ortho'])
    ##    nz = 1
    ##if vacuum_layer_direction in ['z', 'Z']:
    ##    nx = 1
    ##    ny = 1
    ##    nz = math.ceil((vacuum_width + poscar_direction_dict['box_len_ortho']) / poscar_direction_dict['box_len_ortho'])

    ##############################################################
    # Displace the atoms to the bottom of the axis
    ##############################################################
    ##exit()
    align(poscar_file_path = poscar_file_path, align_direction = vacuum_layer_direction, mode = 'bottom', output_poscar_file_path = output_poscar_file_path, delta = delta, radius_scale_factor = radius_scale_factor, write_layer_info = False, suppress_warning = suppress_warning)
    #print('----build vacuum: align_bottom----');exit()
    ##############################################################
    # Build vacuum layer
    ##############################################################
    poscar_dict = vasp_read.read_poscar(output_poscar_file_path)

    crystal_len = np.max(poscar_dict['pos_arr'][:, poscar_direction_dict['pos_arr_column']]) - np.min(poscar_dict['pos_arr'][:,poscar_direction_dict['pos_arr_column']])
    layer_separation = poscar_direction_dict['box_len_ortho'] - crystal_len
    #Determine the added delta_l_arr
    ##poscar_dict = vasp_read.read_poscar(output_poscar_file_path)

    ##layer_separation = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['crystal_width']
    ##delta_len = poscar_layer_dict['vacuum_width'] - layer_separation

    delta_len = vacuum_width - layer_separation
    ##if vacuum_layer_exist == False:
    ##    delta_len = vacuum_width - layer_separation
    ##elif vacuum_layer_exist == True:    
    ##    if layer_separation <= vacuum_width:
    ##        delta_len = vacuum_width - layer_separation
    ##    elif layer_separation > vacuum_width:
    ##        delta_len = vacuum_width - layer_separation
    delta_l_arr = delta_len / poscar_direction_dict['cos_angle'] * poscar_direction_dict['unit_vector']
    ##print('crystal_len=',crystal_len)
    ##print('box_len_ortho=',poscar_direction_dict['box_len_ortho'])
    ##print('layer_separation=',layer_separation)
    ##print('delta_len=',delta_len)
    ##print('delta_l_arr=', delta_l_arr)
    if vacuum_layer_direction in ['x', 'X']:
        t_matrix_arr = np.array([delta_l_arr,[0,0,0],[0,0,0]])
    elif vacuum_layer_direction in ['y', 'Y']:
        t_matrix_arr = np.array([[0,0,0],delta_l_arr,[0,0,0]])
    elif vacuum_layer_direction in ['z', 'Z']:
        t_matrix_arr = np.array([[0,0,0],[0,0,0],delta_l_arr])
    transform(poscar_file_path = output_poscar_file_path,
        t_matrix_arr = t_matrix_arr,
        apply_to = 'box',
        output_poscar_file_path = output_poscar_file_path) 
    return 0

def exfoliation(poscar_file_path, num_periods = 1, layer_lim = None, vacuum_width = 20, exfoliation_direction = 'z', delta = 0.05, layer_dist_tolerance = 'auto', criteria = 'auto', radius_style = 'csd', radius_scale_factor = 1.15, align_center = True, align_upright = True, output_poscar_file_path = None, write_layer_info = False, overwrite = True, suppress_warning = True):
    '''
    Exfoliate 2d material
    NOTE: The inpout POSCAR file should not be started with 'exfoliation_'
    num_periods: the number of periods(collective layers grouped by periodicity) left after exfoliation.
    num_layers: the number of layers left after exfoloiation
    layer_lim: a 2-d list (layer_lim = [m, n]) specific the selected layer from m-th layer to n-th layer, m,n start form 1
    - layer_dist_tolerance: If layer_dist_tolerance = 'auto', then the program determines this value automatically. This is the spacial resolution for the layer distance, for example z=0.24 and z=0.27 are considered to be in the same layer if layer_dist_tolerance = 0.03. If the distance between two atoms in certain direction exceeds this value, then they belongto different layers, otherwise they are in the same layer. A reference value of layer_dist_tolerance is set to 115% of the He-He CSD covalent radius plus an arbitrary small distance (e.g. 0.02 Angstrom): (1.15 * 2 * (CSD_covalent_radius of He)) + arbitrary_delta = 1.15 * 2 * 0.28 + 0.02 = 0.664 Angstrom
    align_center: Align the atoms to the center of the model along the exfoliation direction
    align_upright: Align the exfoliation direction axis as an upright one (i.e. not a tilted one)
    criteria_list: ['auto'], ['bonding'] or ['periodicity'], the criteria in the list will be adopted. For ['auto'], this means a combination of ['bonding'] and ['periodicity'] (first check bonding, then check periodicity) will be adopted. If bonding criterion worked, then periodicity criterion will be skipped; If bonding criterion failed, then periodicity criterion will be adopted.
    If the following error occurs, please tune the value of delta smaller.
            poscar_dict['fix_arr'][i_atom,2] + ' ' + '\n')
            IndexError: index 60 is out of bounds for axis 0 with size 60
    '''
    args_dict = locals()
    import os
    import numpy as np
    import math
    from .. import funcs
    from . import vasp_read, vasp_write, vasp_tools
    from collections import Counter
    from copy import copy, deepcopy
    import itertools
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    poscar_file_path = os.path.abspath(poscar_file_path)
    # Check file type: check if the file is POSCAR format
    pass
    fpath, fname = os.path.split(poscar_file_path)
    if not output_poscar_file_path is None:
        fpath_o, fname_o = os.path.split(output_poscar_file_path)
        if fname == fname_o:
            return 0

    #if fname.startswith('exfoliation_'):
    new_extension = '_exfoliation.vasp'
    if fname.endswith(new_extension):
        return 0

    if output_poscar_file_path is None:
        #output_poscar_file_path = os.path.join(fpath, 'POSCAR_exfoliation.vasp')
        fpre, fext = os.path.splitext(poscar_file_path)
        output_poscar_file_path = fpre + new_extension
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    
    # Whether to overwirte previous output
    if overwrite == True:
        if os.path.isfile(output_poscar_file_path):
            os.remove(output_poscar_file_path)
    else:
        if os.path.isfile(output_poscar_file_path):
            return 0 

    # reorient the model
    #mode: ['a along x, b in xy', 'b along x, c in xy', 'c along x, a in xy']:
    orientation(poscar_file_path = poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(output_poscar_file_path)

    #poscar_dict = vasp_read.read_poscar(poscar_file_path)

    #######################################
    # Initialize
    #######################################
    #delta: This is the spacial resolution (in unit of Angstrom) for the atoms, for example z=0.24 and z=0.25 are considered to be in the same position if delta=0.01.
    delta = 0.05

    poscar_direction_dict = vasp_tools.poscar_direction_params(poscar_dict = poscar_dict, direction = exfoliation_direction)
    if exfoliation_direction in ['x', 'X']:
        nx = math.ceil((vacuum_width + poscar_direction_dict['box_len_ortho']) / poscar_direction_dict['box_len_ortho'])
        ny = 1
        nz = 1
        if num_periods >= nx:
            nx = num_periods + 1
        if nx > 100:
            nx = 1
    if exfoliation_direction in ['y', 'Y']:
        nx = 1
        ny = math.ceil((vacuum_width + poscar_direction_dict['box_len_ortho']) / poscar_direction_dict['box_len_ortho'])
        nz = 1
        if num_periods >= ny:
            ny = num_periods + 1
        if ny > 100:
            ny = 1
    if exfoliation_direction in ['z', 'Z']:
        nx = 1
        ny = 1
        nz = math.ceil((vacuum_width + poscar_direction_dict['box_len_ortho']) / poscar_direction_dict['box_len_ortho'])
        if num_periods >= nz:
            nz = num_periods + 1
        if nz > 100:
            nz = 1

    #######################################
    # Reduce to primitive cell
    #######################################
    pass
    ############################################################
    # Reorientate the model
    ############################################################
    orientation(poscar_file_path = poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = output_poscar_file_path)
    poscar_dict_reorientation = vasp_read.read_poscar(output_poscar_file_path)
    poscar_dict_old = deepcopy(poscar_dict_reorientation) 
    vec_a_old = poscar_dict['box_len_arr'][0, :] * poscar_dict['uni_scale_fac'] 
    vec_b_old = poscar_dict['box_len_arr'][1, :] * poscar_dict['uni_scale_fac']
    vec_c_old = poscar_dict['box_len_arr'][2, :] * poscar_dict['uni_scale_fac']
    cell_arr_old = np.vstack((vec_a_old, vec_b_old, vec_c_old))
    basis_vector_dict_old = funcs.basis_vector_info(cell_arr_old)
    #print('---orientation----');exit()
    ############################################################
    # Expand model in a certain direction by creating supercell
    ############################################################
    create_supercell(poscar_file_path = output_poscar_file_path,
        num_xyz_list = [nx, ny, nz],
        output_poscar_file_path = output_poscar_file_path,
        ) 
    #print('----supercell--------');exit()
    ##############################################################
    # Get layer parameters from POSCAR
    ##############################################################
    poscar_dict = vasp_read.read_poscar(output_poscar_file_path)
    poscar_layer_dict = vasp_tools.poscar_layer_params(poscar_dict, poscar_direction_dict, criteria = criteria, delta = delta, layer_dist_tolerance = layer_dist_tolerance, radius_style = radius_style, radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning = suppress_warning)

    ###########################################################
    # Screen the region of intrest (ROI) by certain criteria
    ###########################################################
    #For criteria = 'auto', this means a combination of bonding and periodicity (first check bonding, then check periodicity) will be adopted. If bonding criterion worked, then periodicity criterion will be skipped; If bonding criterion failed, then periodicity criterion will be adopted.
    #######################
    # bonding criterion
    #######################
    # If the criteria_list is 'bonding', then the layers will be distinguished by the bonding status. If all the atoms in a layer have no bonding with the atoms in the adjacent layer, then they are said to be the boundary layer. This option is suitable for vdW materials
    if poscar_layer_dict['criteria'] == 'bonding' and -1 in poscar_layer_dict['bond_status_list'] and 1 in poscar_layer_dict['bond_status_list']:
        print('# From vasp_build.exfoliation: USE BONDING CRITERION')
        num_elmt_list = []
        elmt_species_list = []
        atom_indx_list = []
        for i_elmt in poscar_dict['elmt_species_arr']:
            i_num_elmt = 0
            #stacked_period_number = 0
            if not layer_lim is None:
                if layer_lim[1] > poscar_layer_dict['num_layers']:
                    min_layer_indx = layer_lim[0] - 1
                    max_layer_indx = poscar_layer_dict['num_layers'] - 1
                else:
                    min_layer_indx = layer_lim[0] - 1
                    max_layer_indx = layer_lim[1] - 1
                num_periods = 999
            else:
                min_layer_indx = 0
                max_layer_indx = poscar_layer_dict['num_layers'] - 1
            for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
               if (poscar_layer_dict['cell_by_bond_list'][i_layer_indx] <= num_periods and poscar_layer_dict['cell_by_bond_list'][i_layer_indx] != 0) and ((i_layer_indx >= min_layer_indx) and (i_layer_indx <= max_layer_indx)):
                   for i_atom_indx in poscar_layer_dict['atoms_indx_list'][i_layer_indx]:
                       if (i_elmt == poscar_dict['atom_species_arr'][i_atom_indx]):
                           atom_indx_list.append(i_atom_indx)
                           i_num_elmt = i_num_elmt + 1
            if i_num_elmt != 0:
                num_elmt_list.append(i_num_elmt)
                elmt_species_list.append(i_elmt)
        poscar_dict['elmt_num_arr'] = np.array(num_elmt_list) 
        poscar_dict['elmt_species_arr'] = np.array(elmt_species_list) 
        pos_arr_old = np.copy(poscar_dict['pos_arr'])
        poscar_dict['n_atoms'] = len(atom_indx_list) 
        poscar_dict['pos_arr'] = np.array([0.0] * poscar_dict['n_atoms'] * 6,dtype = np.float)
        poscar_dict['pos_arr'].shape = poscar_dict['n_atoms'], 6
        ##new_atom_indx = 0
        ##for i_atom_indx in atom_indx_list:
        ##    poscar_dict['pos_arr'][new_atom_indx, 3:6] = pos_arr_old[i_atom_indx, 3:6]
        ##    new_atom_indx = new_atom_indx + 1
    elif poscar_layer_dict['criteria'] == 'bonding' and 0 in poscar_layer_dict['bond_status_list']:
        #There are one or more layers which do not form bonds with the adjacent layers. Besides the bonding criterion, periodicity criterion usually will not handle this situation perfectly.
        if suppress_warning == False:
            print('WARNING #2103221330 (from exfoliation): Cannot exfoliate layers according to the bonding state. There are one or more atomic layers which do not form bond with both sides of the adjacent layers. The bond_status_list is: ' + str(poscar_layer_dict['bond_status_list']) + '. Please check the model manually or enlarge the threshold value for atom bonding by setting a larger multiplying factor for the atomic radius. Current, the \'periodicity\' criterion will be used to exfoliate the layers : ' + poscar_file_path)
        poscar_layer_dict['criteria'] = 'periodicity'
    elif poscar_layer_dict['criteria'] == 'bonding' and 2 in poscar_layer_dict['bond_status_list'] and len(set(poscar_layer_dict['bond_status_list'])) == 1:
        if suppress_warning == False:
            print('WARNING #2012301623 (from exfoliation): Cannot exfoliate layers according to the bonding state. All the layers form bonds with their adjacent layers. The \'periodicity\' criterion will be used to exfoliate the layers: ' + poscar_file_path)
        poscar_layer_dict['criteria'] = 'periodicity'
    #######################
    # periodicity criterion
    #######################
    # If the criteria_list is 'periodicity', then the layers will be distinguished by the periodicity of layers, i.e. ABCABCABC. This option is also suitable for bulk materials.
    if poscar_layer_dict['criteria'] == 'periodicity' and None not in poscar_layer_dict['bond_status_list']:
        print('# From vasp_build.exfoliation: USE PERIODICITY CRITERION')
        num_elmt_list = []
        elmt_species_list = []
        atom_indx_list = []
        
        characteristic_periodicity_indx = 0
        for i_elmt in poscar_dict['elmt_species_arr']:
            #admit_new_period = True
            #stacked_period_number = 0
            i_num_elmt = 0
            # specify range of layers
            if not layer_lim is None:
                if layer_lim[1] > poscar_layer_dict['num_layers']:
                    min_layer_indx = layer_lim[0] - 1
                    max_layer_indx = poscar_layer_dict['num_layers'] - 1
                else:
                    min_layer_indx = layer_lim[0] - 1
                    max_layer_indx = layer_lim[1] - 1
                num_periods = 999
            else:
                min_layer_indx = 0
                max_layer_indx = poscar_layer_dict['num_layers'] - 1

            # Screen according to the periodicity and layer indices
            # For the ABCABCABC lattice, use the VACUUM + ABCA + VACUUM configuration instead of VACUUM + ABC + VACUUM configuration. This vacuum_abca_vacuum condition is suited for the 2-dimensional configuration. the a_layer_indx denotes the layer index of the last 'A' layer in the ABCA configuration.
            cell_by_periodicity_arr = np.array(poscar_layer_dict['cell_by_periodicity_list'])
            periodicity_status_arr = np.array(poscar_layer_dict['periodicity_status_list'])
            next_cell_layer_indx_list = np.argwhere(cell_by_periodicity_arr == num_periods + 1)
            periodicity_layer_characteristic_indx_list = np.argwhere(periodicity_status_arr == characteristic_periodicity_indx)
            temp_arr = next_cell_layer_indx_list[np.in1d(next_cell_layer_indx_list, periodicity_layer_characteristic_indx_list)]
            if len(temp_arr) != 0:
                a_layer_indx = temp_arr[0]
            else:
                a_layer_indx = None

            for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
               vacuum_abca_vacuum_condition = (i_layer_indx == a_layer_indx)
               if ((poscar_layer_dict['cell_by_periodicity_list'][i_layer_indx] <= num_periods and poscar_layer_dict['cell_by_periodicity_list'][i_layer_indx] != 0) and ((i_layer_indx >= min_layer_indx) and (i_layer_indx <= max_layer_indx))) or vacuum_abca_vacuum_condition:
                   for i_atom_indx in poscar_layer_dict['atoms_indx_list'][i_layer_indx]:
                       if (i_elmt == poscar_dict['atom_species_arr'][i_atom_indx]):
                           atom_indx_list.append(i_atom_indx)
                           i_num_elmt = i_num_elmt + 1

            if i_num_elmt != 0:
                num_elmt_list.append(i_num_elmt)
                elmt_species_list.append(i_elmt)
        ##print('atom_indx_list=',atom_indx_list)
        ##print([k for k,v in Counter(atom_indx_list).items() if v>1])
        ##print('len_atom_indx_list=',len(atom_indx_list))
        poscar_dict['elmt_num_arr'] = np.array(num_elmt_list) 
        poscar_dict['elmt_species_arr'] = np.array(elmt_species_list) 
        #pos_arr_old = np.copy(poscar_dict['pos_arr'])
        pos_arr_old = deepcopy(poscar_dict['pos_arr'])
        poscar_dict['n_atoms'] = len(atom_indx_list) 
        poscar_dict['pos_arr'] = np.array([0.0] * poscar_dict['n_atoms'] * 6,dtype = np.float)
        poscar_dict['pos_arr'].shape = (poscar_dict['n_atoms'], 6)
    ###################################################################################
    # find the coordinates of the selected atoms according to the selected atom index
    ###################################################################################
    new_atom_indx = 0
    basis_vector_dict = funcs.basis_vector_info(poscar_dict['l_arr'])
    for i_atom_indx in atom_indx_list:
        poscar_dict['pos_arr'][new_atom_indx, 3:6] = pos_arr_old[i_atom_indx, 3:6]
        # the direct coordinate of the atom relative to the original cell(no supercell built)
        poscar_dict['pos_arr'][new_atom_indx, 0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][new_atom_indx,3:6])
        poscar_dict['pos_arr'][new_atom_indx,0:3] = np.modf(np.modf(poscar_dict['pos_arr'][new_atom_indx,0:3])[0] + 1)[0]
        poscar_dict['pos_arr'][new_atom_indx,0:3][np.isclose(poscar_dict['pos_arr'][new_atom_indx,0:3], 1, atol = 1e-16)] = 0
        #poscar_dict['pos_arr'][new_atom_indx, 0:3] = basis_vector_dict_old['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][new_atom_indx,3:6])
        new_atom_indx = new_atom_indx + 1
    ##vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    #print('--periodicity-----');exit()

    #######################################################################################
    # Solve the problem of the topmost and the bottom atoms that are in the same layers
    #######################################################################################
    base_list = set(math.floor(poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['l_arr_row']]) for i_atom_indx in range(0, poscar_dict['n_atoms']))
    max_bound = None
    for i_z in range(min(base_list), nz):
        if (i_z + 1) not in base_list:
            max_bound = i_z
            break
    if not max_bound is None:
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            #poscar_dict['pos_arr'][i_atom_indx,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
            if math.floor(poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['l_arr_row']]) > max_bound:
                poscar_dict['pos_arr'][i_atom_indx,poscar_direction_dict['l_arr_row']] += (-math.floor(poscar_dict['pos_arr'][i_atom_indx,poscar_direction_dict['l_arr_row']]) + max_bound) 
                ##print(poscar_dict['pos_arr'][i_atom_indx,poscar_direction_dict['l_arr_row']])
    ##vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    ##print('--max bound-----');exit()

    ############################################################ 
    # Convert the atom coordinate to cartesian coordinate 
    ############################################################ 
    pass
    ############################################################ 
    # Use the direct coordinate to write the final configuration
    ############################################################ 
    poscar_dict['pos_arr'][:,0:3] = np.modf(np.modf(poscar_dict['pos_arr'][:,0:3])[0] + 1)[0]
    poscar_dict['pos_arr'][:,0:3][np.isclose(poscar_dict['pos_arr'][:,0:3], 1, atol = 1e-16)] = 0
    
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    #print('--periodicity-----');exit()
    #######################################
    # Adjust the vacuum layer width
    #######################################
    build_vacuum_layer(poscar_file_path = output_poscar_file_path, vacuum_layer_direction = exfoliation_direction, vacuum_width = vacuum_width, threshold_vacuum_width = 5.98, output_poscar_file_path = output_poscar_file_path, delta = delta, radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning = suppress_warning)
    #print('--build vacuum-----');exit()
    ##############################################################################
    # Align the atoms to the center of the model along the exfoliation direction
    ##############################################################################
    # Displace the atoms to the center of the axis
    if align_center == True:
        align(poscar_file_path = output_poscar_file_path, align_direction = exfoliation_direction, mode = 'center', output_poscar_file_path = output_poscar_file_path, delta = delta, radius_scale_factor = radius_scale_factor, write_layer_info = False, suppress_warning = suppress_warning)
    ##############################################################################
    # Align the exfoliation direction axis as an upright one (i.e. not a tilted one)
    ##############################################################################
    if align_upright == True:
        align(poscar_file_path = output_poscar_file_path, align_direction = exfoliation_direction, mode = 'upright_box', output_poscar_file_path = output_poscar_file_path)
    ##############################################################################
    # Check for the validity of the generated model, i.e. whether the layers 
    # have been distinguished from each other.
    # Specifically, We adaptively tune the layer_dist_tolerance to get a 
    # layer resolved model
    # There is another solution, if you don't want to check the validity of the 
    # model aftering exfoliating, you'd better wisely and automatically 
    # generate the the suitable value of layer_dist_tolerance, which might be a 
    # better solution.
    ##############################################################################
    # To be considered
    ##print('natoms=',poscar_dict['n_atoms'])
    #print('#########################################################################')
    log_str = ''
    ##################################
    # Determine the args string
    ##################################
    func_name = 'vasp_build.exfoliation'
    args_str = func_name + '(' + '\n'
    for i_arg in args_dict.keys():
        arg_value = args_dict[i_arg]
        if isinstance(arg_value,str):
            arg_value_str = '\'' + arg_value + '\''
        else:
            arg_value_str = str(arg_value)
        if i_arg == 'poscar_file_path':
            arg_value_str = 'r\'' + poscar_file_path + '\''
        args_str += '    ' + i_arg + ' = ' + arg_value_str + ',\n'
    args_str += '    )\n'
    args_str += '################################################\n'

    log_str += args_str
    funcs.write_log(logfile, log_str)
    return 0

def conventional_cell(poscar_file_path, delta = 0.05, radius_scale_factor = 1.15, align_center = False, align_upright = True, vacuum_width = None, output_poscar_file_path = None, write_layer_info = False, suppress_warning = True):
    '''
    generate conventional cell by rotating the model and other transformation.
    vacuum_width: This is used to readjust the vacuum layer width. If vacuum_width = None, the vacuum width is not adjusted.
    Reference:
    [Wahyu Setyawan and Stefano Curtarolo, Computational Materials Science, 2010, 49, 299-312]
    '''
    args_dict = locals()
    import os
    import numpy as np
    import math
    from .. import funcs
    from . import vasp_read, vasp_write, vasp_tools

    conv_cell_dict = {}
    conv_cell_dict['file_path'] = None
    conv_cell_dict['vacuum_exist'] = None
    conv_cell_dict['vacuum_direction'] = None
    conv_cell_dict['vacuum_width'] = None
    conv_cell_dict['angle_alpha_degree'] = None 
    conv_cell_dict['angle_beta_degree']  = None 
    conv_cell_dict['angle_gamma_degree'] = None 
    poscar_file_path = os.path.abspath(poscar_file_path)
    conv_cell_dict['file_path'] = poscar_file_path
    # Check file type: check if the file is POSCAR format
    pass
    new_extension = '_conventional.vasp'
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path is None:
        fpre, fext = os.path.splitext(poscar_file_path)
        output_poscar_file_path = fpre + new_extension
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
 
    vacuum_exist = False
    ############################################################
    # Check three situations:
    # a along x, b in xy
    # b along x, c in xy
    # c along x, a in xy
    ############################################################
    
    # Check whether there is vacuum layer in the model
    poscar_dict = vasp_read.read_poscar(poscar_file_path)
    conv_cell_dict['angle_alpha_degree'] = poscar_dict['angle_alpha_degree']
    conv_cell_dict['angle_beta_degree']  = poscar_dict['angle_beta_degree']
    conv_cell_dict['angle_gamma_degree'] = poscar_dict['angle_gamma_degree']
    conv_cell_dict['vacuum_exist'] = False
    for i_mode in ['a along x, b in xy', 'b along x, c in xy', 'c along x, a in xy']:
        orientation(poscar_file_path = poscar_file_path, mode = i_mode, output_poscar_file_path = output_poscar_file_path)
        poscar_dict = vasp_read.read_poscar(output_poscar_file_path)
        poscar_direction_dict = vasp_tools.poscar_direction_params(poscar_dict, direction = 'z')
        poscar_layer_dict = vasp_tools.poscar_layer_params(poscar_dict, poscar_direction_dict, criteria = 'bonding', delta = delta, layer_dist_tolerance = 'auto', radius_style = 'csd', radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning = suppress_warning)
        if poscar_layer_dict['vacuum_exist'] == True:
            ##print('Vacuum exist in orientation: ', i_mode)
            conv_cell_dict['vacuum_exist'] = True
            conv_cell_dict['vacuum_width'] =  poscar_layer_dict['vacuum_width'] 
            if i_mode == 'a along x, b in xy':
                conv_cell_dict['vacuum_direction'] = 'c'
            elif i_mode == 'b along x, c in xy':
                conv_cell_dict['vacuum_direction'] = 'a'
            elif i_mode == 'c along x, a in xy':
                conv_cell_dict['vacuum_direction'] = 'b'
            break

    #######################################
    # Adjust the vacuum layer width
    #######################################
    if conv_cell_dict['vacuum_exist'] == True:
        if vacuum_width not in [None, 'None', 'none']:
            build_vacuum_layer(poscar_file_path = output_poscar_file_path, vacuum_layer_direction = 'z', vacuum_width = vacuum_width, threshold_vacuum_width = 5.98, output_poscar_file_path = output_poscar_file_path, delta = delta, radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning = suppress_warning)
    elif vacuum_exist == False:
        orientation(poscar_file_path = poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = output_poscar_file_path)

    ##############################################################################
    # Align the atoms to the center of the model along the exfoliation direction
    ##############################################################################
    # Displace the atoms to the center of the axis
    if align_center == True:
        align(poscar_file_path = output_poscar_file_path, align_direction = 'z', mode = 'center', output_poscar_file_path = output_poscar_file_path, delta = delta, radius_scale_factor = radius_scale_factor, write_layer_info = False, suppress_warning = suppress_warning)
    ##############################################################################
    # Align the exfoliation direction axis as an upright one (i.e. not a tilted one)
    ##############################################################################
    if align_upright == True:
        align(poscar_file_path = output_poscar_file_path, align_direction = 'z', mode = 'upright_box', output_poscar_file_path = output_poscar_file_path, delta = delta, radius_scale_factor = radius_scale_factor, write_layer_info = False, suppress_warning = suppress_warning)
    return conv_cell_dict 

def car2frac(poscar_file_path, output_poscar_file_path = None, suppress_warning = True):
    '''
    convert Cartesian coordinate to Fractional coordinate
    '''
    args_dict = locals()
    import os
    from .. import funcs
    from . import vasp_read, vasp_write, vasp_tools
    poscar_file_path = os.path.abspath(poscar_file_path)
    new_extension = '_car2frac.vasp'
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path is None:
        fpre, fext = os.path.splitext(poscar_file_path)
        output_poscar_file_path = fpre + new_extension
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    orientation(poscar_file_path = poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = output_poscar_file_path)
    return 0
