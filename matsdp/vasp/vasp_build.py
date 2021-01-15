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
    
    funcs.write_log(
        logfile,
        'vasp_build.substitution(' + '\n' +
        '    substitution_list_file=' + 'r\'' + str(substitution_list_file_abs_path) + '\'' + ',\n' +
        '    poscar_file_path='  + 'r\'' + str(poscar_file_path) + '\'' + ')\n' +
        '###############################\n')
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
    import os
    import numpy as np
    import sys
    import time
    from .. import funcs
    from .. import default_params

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
    funcs.write_log(
        logfile,
        'vasp_build.rep_elmt(' + '\n' +
        '    substitution_list_file='  + 'r\'' + str(substitution_list_file) + '\'' + ',\n' +
        '    poscar_file_path='  + 'r\'' + str(poscar_file_path) + '\'' + ',\n' +
        '    old_elmt=' + '\'' + str(old_elmt) + '\'' + ',\n' +
        '    elmt_group=' + str(elmt_group) + ')\n' +
        '###############################\n')
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
    
    funcs.write_log(logfile,
        'vasp_build.selection_sphere(' + '\n' +
        '    poscar_file_path=' + 'r\'' + str(poscar_file_path) + '\'' + ',\n' +
        '    origin_atom_name=\'' + str(origin_atom_name) + '\',\n' + 
        '    radius=' + str(radius) + ',\n' +
        '    include_mirror_atoms=' + str(include_mirror_atoms) + ',\n' +
        '    output_file_name=\'' + str(output_file_name) + '\')\n')
    return 0 

def create_multiple_vasp_jobs(substitution_list_file, poscar_file_path, elmt_potcar_dir, incar_str, kpoints_str):
    '''
    create POSCAR, INCAR, KPOINTS, POTCAR files, this module is based on the substitution module
    The POSCAR is created by atom substitution
    '''
    import os
    from .. import funcs
    from . import vasp_write
    from .. import default_params

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
    funcs.write_log(logfile,
        'vasp_build.create_multiple_vasp_jobs(' + '\n' +
        '    substitution_list_file = ' + 'r\'' + str(substitution_list_file) + '\'' + ',\n' +
        '    poscar_file_path = ' + 'r\'' + str(poscar_file_path) + '\'' + ',\n' +
        '    elmt_potcar_dir = ' + 'r\'' + str(elmt_potcar_dir) + '\'' + ',\n' +
        '    incar_str = \'\'\'' + str(incar_str) + '\'\'\',\n' +
        '    kpoints_str = \'\'\'' + str(kpoints_str) + '\'\'\')\n' +
        '############################\n'
        )
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
    import os
    import numpy as np
    import math
    from .. import funcs
    from . import vasp_read
    from . import vasp_write
    poscar_file_path = os.path.abspath(poscar_file_path)
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path in [None, 'None', 'none']:
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_orientation')
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)

    def read_mode(mode, l_arr):
        '''
        get rotation from the mode variable
        l_arr: [vec_a, vec_b, vec_c]
        '''
        x_vector = np.array([1, 0, 0])    
        y_vector = np.array([0, 1, 0])    
        z_vector = np.array([0, 0, 1])    
        # abc_vec_1 along xyz_vec_1, abc_vec_2 in xyz_plane
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
        
        if [abc_vec_indx_1, abc_vec_indx_2] in [[0, 1], [1, 2], [2, 0]]:
            abc_vec1_cross_vec2 = np.cross(abc_vec_1, abc_vec_2)
        elif [abc_vec_indx_1, abc_vec_indx_2] in [[1, 0], [2, 1], [0, 2]]:
            abc_vec1_cross_vec2 = np.cross(abc_vec_2, abc_vec_1)

        abc_vec_1 = (abc_vec_1 / np.linalg.norm(abc_vec_1)).reshape((3, 1))
        xyz_vec_1 = (xyz_vec_1 / np.linalg.norm(xyz_vec_1)).reshape((3, 1))
        abc_vec1_cross_vec2 = (abc_vec1_cross_vec2 / np.linalg.norm(abc_vec1_cross_vec2)).reshape((3, 1))
        xyz_plane_vec = (xyz_plane_vec / np.linalg.norm(xyz_plane_vec)).reshape((3, 1))
        #print(np.shape(abc_vec_1[np.newaxis, :].T))
        #matrix_arr_1 = xyz_vec_1 * abc_vec_1[np.newaxis, :].T
        #matrix_arr_2 = xyz_plane_vec * abc_vec1_cross_vec2[np.newaxis, :].T
        #matrix_arr_1 = xyz_vec_1.dot(np.linalg.pinv(abc_vec_1))
        #matrix_arr_2 = xyz_plane_vec.dot(np.linalg.pinv(abc_vec1_cross_vec2))
        matrix_arr_1 = funcs.rotation_matrix_from_vectors(abc_vec_1, xyz_vec_1)
        abc_vec1_cross_vec2 = matrix_arr_1.dot(abc_vec1_cross_vec2)
        matrix_arr_2 = funcs.rotation_matrix_from_vectors(abc_vec1_cross_vec2, xyz_plane_vec)
        rotation_matrix_arr = matrix_arr_2.dot(matrix_arr_1)
        return rotation_matrix_arr

    rotation_matrix_arr = read_mode(mode = mode, l_arr = poscar_dict['l_arr'])
    # Only rotate the simulation box, do not rotate the atom coordinates
    vec_a_old = poscar_dict['l_arr'][0, :]
    vec_b_old = poscar_dict['l_arr'][1, :]
    vec_c_old = poscar_dict['l_arr'][2, :]
    poscar_dict['l_arr'][0, :] = rotation_matrix_arr.dot(poscar_dict['l_arr'][0, :])
    poscar_dict['l_arr'][1, :] = rotation_matrix_arr.dot(poscar_dict['l_arr'][1, :])
    poscar_dict['l_arr'][2, :] = rotation_matrix_arr.dot(poscar_dict['l_arr'][2, :])
    poscar_dict['box_len_arr'][0, :] = poscar_dict['l_arr'][0, :] / poscar_dict['uni_scale_fac']
    poscar_dict['box_len_arr'][1, :] = poscar_dict['l_arr'][1, :] / poscar_dict['uni_scale_fac']
    poscar_dict['box_len_arr'][2, :] = poscar_dict['l_arr'][2, :] / poscar_dict['uni_scale_fac']

    ### For a best solution, the atoms outside the simulation should be translated back into the box
    ##vec_1 = poscar_dict['l_arr'][0, :] 
    ##vec_2 = poscar_dict['l_arr'][1, :] 
    ##vec_3 = poscar_dict['l_arr'][2, :] 
    ##volume_parallelepiped = abs(np.cross(vec_1, vec_2).dot(vec_3))
    ### Find eight vertices of the parallelepiped. 
    ##point_1 = np.array([0, 0, 0])
    ##point_2 = vec_1
    ##point_3 = vec_2
    ##point_4 = vec_3
    ##point_5 = vec_1 + vec_2
    ##point_6 = vec_1 + vec_3
    ##point_7 = vec_2 + vec_3
    ##point_8 = vec_1 + vec_2 + vec_3
    ##delta_volume = 0.01
    ##for i_atom_indx in range(0, poscar_dict['n_atoms']):
    ##    point_0 = poscar_dict['pos_arr'][i_atom_indx, 3:6]
    ##    volume_six_pentahedra = (
    ##        funcs.volume_pentahedron_from_points(point_1, point_2, point_3, point_5, point_0) + 
    ##        funcs.volume_pentahedron_from_points(point_1, point_2, point_4, point_6, point_0) + 
    ##        funcs.volume_pentahedron_from_points(point_1, point_3, point_4, point_7, point_0) + 
    ##        funcs.volume_pentahedron_from_points(point_8, point_4, point_6, point_7, point_0) + 
    ##        funcs.volume_pentahedron_from_points(point_8, point_3, point_5, point_7, point_0) + 
    ##        funcs.volume_pentahedron_from_points(point_8, point_2, point_5, point_6, point_0))
    ##    # for points not in the parallelepiped, translate them to a position inside the parallelepiped.
    ##    n1 = 1
    ##    n2 = 1
    ##    n3 = 1
    ##    max_dim = 3
    ##    n1_arr = np.array([None] * max_dim * 2)
    ##    i = 0
    ##    for i_indx in range(0, max_dim + 1):
    ##        n1_arr[i] = i_indx
    ##        if i_indx != max_dim:
    ##            n1_arr[i + 1] = -i_indx
    ##        if i == 0:
    ##            i = i + 1
    ##        else:
    ##            i = i + 2
    ##    n2_arr = n1_arr.copy()
    ##    n3_arr = n1_arr.copy()
    ##    if abs(volume_six_pentahedra - volume_parallelepiped) > delta_volume:
    ##        #print('---',volume_six_pentahedra, volume_parallelepiped)
    ##        shift_arr_1 = vec_1 
    ##        shift_arr_2 = vec_2
    ##        shift_arr_3 = vec_3
    ##        for n1 in n1_arr:
    ##            for n2 in n2_arr:
    ##                for n3 in n3_arr:
    ##                    if n1 == 0 and n2 == 0 and n3 == 0:
    ##                        continue
    ##                    point_0p = point_0 - n1 * shift_arr_1 - n2 * shift_arr_2 - n3 * shift_arr_3  
    ##                    ##point_1p = point_1 + n1 * shift_arr_1 + n2 * shift_arr_2 + n3 * shift_arr_3  
    ##                    ##point_2p = point_2 + n1 * shift_arr_1 + n2 * shift_arr_2 + n3 * shift_arr_3
    ##                    ##point_3p = point_3 + n1 * shift_arr_1 + n2 * shift_arr_2 + n3 * shift_arr_3
    ##                    ##point_4p = point_4 + n1 * shift_arr_1 + n2 * shift_arr_2 + n3 * shift_arr_3
    ##                    ##point_5p = point_5 + n1 * shift_arr_1 + n2 * shift_arr_2 + n3 * shift_arr_3  
    ##                    ##point_6p = point_6 + n1 * shift_arr_1 + n2 * shift_arr_2 + n3 * shift_arr_3
    ##                    ##point_7p = point_7 + n1 * shift_arr_1 + n2 * shift_arr_2 + n3 * shift_arr_3
    ##                    ##point_8p = point_8 + n1 * shift_arr_1 + n2 * shift_arr_2 + n3 * shift_arr_3
    ##                    volume_six_pentahedra_1 = (
    ##                        funcs.volume_pentahedron_from_points(point_1, point_2, point_3, point_5, point_0p) + 
    ##                        funcs.volume_pentahedron_from_points(point_1, point_2, point_4, point_6, point_0p) + 
    ##                        funcs.volume_pentahedron_from_points(point_1, point_3, point_4, point_7, point_0p) + 
    ##                        funcs.volume_pentahedron_from_points(point_8, point_4, point_6, point_7, point_0p) + 
    ##                        funcs.volume_pentahedron_from_points(point_8, point_3, point_5, point_7, point_0p) + 
    ##                        funcs.volume_pentahedron_from_points(point_8, point_2, point_5, point_6, point_0p))
    ##                    if abs(volume_six_pentahedra_1 - volume_parallelepiped) <= delta_volume:
    ##                        #print('dddd',volume_six_pentahedra_1, volume_parallelepiped)
    ##                        ##poscar_dict['pos_arr'][i_atom_indx, 3:6] = poscar_dict['pos_arr'][i_atom_indx, 3:6] - n1 * shift_arr_1 - n2 * shift_arr_2 - n3 * shift_arr_3
    ##                        poscar_dict['pos_arr'][i_atom_indx, 3:6] = point_0p
    ##                        break

    vec_a = poscar_dict['box_len_arr'][0, :] * poscar_dict['uni_scale_fac'] 
    vec_b = poscar_dict['box_len_arr'][1, :] * poscar_dict['uni_scale_fac']
    vec_c = poscar_dict['box_len_arr'][2, :] * poscar_dict['uni_scale_fac']
    basis_vector_dict = funcs.basis_vector_info(vec_a, vec_b, vec_c)
   
    # The relations among orientation, Cartesian and Fractional coordinates are complex!
    # if in the future, the car2fra_matrix_arr function is generalized, the following codes may be changed too.
    write_mode = 3
    if write_mode == 1:
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            poscar_dict['pos_arr'][i_atom_indx,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
        vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    elif write_mode == 2:
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            poscar_dict['pos_arr'][i_atom_indx,3:6] = rotation_matrix_arr.dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
        vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Cartesian')
    elif write_mode == 3:
        # now is only valid for 'a along x, b in xy' case, the code will be generalized once the the car2fra_matrix_arr function is generalized
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            poscar_dict['pos_arr'][i_atom_indx,3:6] = rotation_matrix_arr.dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            poscar_dict['pos_arr'][i_atom_indx,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
        vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    return 0

def align(poscar_file_path, align_direction = 'z', mode = 'bottom', output_poscar_file_path = None):
    '''
    mode = 'bottom' or 'center' or 'upright_box'
    mode = 'bottom' or 'center': this is only applied to atoms. Align the atoms to the bottom or center of the box in a certain direction
    mode = 'upright_box': Make the tilted box an "upright" one. Align one of the basis vector (a, b, or c) to x, y, or z direction(whether its x, y, or z is desiginated by the align_direction tag), without changing the box volume. First rotate the model, to reorientate the a, b, and c axes according to the x, y, and z axes (align one of the a, b, or c to one of the x, y or z axis, then make two of the a, b, and c axes in the same plane of two of the x, y, and z axes), then make two of the a, b, and c axes penpendicular to each other.
    '''
    import os
    import numpy as np
    import math
    from .. import funcs
    from . import vasp_read
    from . import vasp_write
    poscar_file_path = os.path.abspath(poscar_file_path)
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path in [None, 'None', 'none']:
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_align')
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)

    if align_direction in ['x', 'X']:
        l_arr_row = 0
        side_vector = poscar_dict['vec_a']
        side_vector_len = poscar_dict['len_vec_a']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([1, 0, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 3
        pos_arr_column_1 = 4
        pos_arr_column_2 = 5
    if align_direction in ['y', 'Y']:
        l_arr_row = 1
        side_vector = poscar_dict['vec_b']
        side_vector_len = poscar_dict['len_vec_b']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 1, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 4
        pos_arr_column_1 = 3
        pos_arr_column_2 = 5
    if align_direction in ['z', 'Z']:
        l_arr_row = 2
        side_vector = poscar_dict['vec_c']
        side_vector_len = poscar_dict['len_vec_c']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 0, 1])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 5

    if mode == 'bottom':
        # align bottom
        lowest_atom_indx = int(np.argwhere(poscar_dict['pos_arr'][:, pos_arr_column] == np.min(poscar_dict['pos_arr'][:, pos_arr_column]))[0])
        displace_vector = poscar_dict['pos_arr'][lowest_atom_indx, 3:6].dot(ortho_vector) / cos_angle * (-unit_vector)
        poscar_dict['pos_arr'][:, 3:6] = poscar_dict['pos_arr'][:, 3:6] + displace_vector
    elif mode == 'center':
        # align bottom to the bottom first
        lowest_atom_indx = int(np.argwhere(poscar_dict['pos_arr'][:, pos_arr_column] == np.min(poscar_dict['pos_arr'][:, pos_arr_column]))[0])
        displace_vector = poscar_dict['pos_arr'][lowest_atom_indx, 3:6].dot(ortho_vector) / cos_angle * (-unit_vector)
        # Displace the atoms to the center of the axis
        crystal_len_ortho = np.max(poscar_dict['pos_arr'][:, pos_arr_column]) - np.min(poscar_dict['pos_arr'][:, pos_arr_column])
        displace_vector = displace_vector + 1 / 2 * (box_len_ortho - crystal_len_ortho) / cos_angle * unit_vector
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
        if align_direction in ['x', 'X']:
            l_arr_column_1 = 1
            l_arr_column_2 = 2
            pos_arr_column_1 = 4
            pos_arr_column_2 = 5
            vec_1 = poscar_dict['vec_b']
            vec_2 = poscar_dict['vec_c']
            ortho_vector_1 = np.array([0, 1, 0])
            ortho_vector_2 = np.array([0, 0, 1])
        elif align_direction in ['y', 'Y']:
            l_arr_column_1 = 0
            l_arr_column_2 = 2
            pos_arr_column_1 = 3
            pos_arr_column_2 = 5
            vec_1 = poscar_dict['vec_a']
            vec_2 = poscar_dict['vec_c']
            ortho_vector_1 = np.array([1, 0, 0])
            ortho_vector_2 = np.array([0, 0, 1])
        elif align_direction in ['z', 'Z']:
            l_arr_column_1 = 0
            l_arr_column_2 = 1
            pos_arr_column_1 = 3
            pos_arr_column_2 = 4
            vec_1 = poscar_dict['vec_a']
            vec_2 = poscar_dict['vec_b']
            ortho_vector_1 = np.array([1, 0, 0])
            ortho_vector_2 = np.array([0, 1, 0])
        projected_side_vector_len_1 = vec_1.dot(ortho_vector_1)
        projected_side_vector_len_2 = vec_2.dot(ortho_vector_2)
        projected_side_vector_1 = projected_side_vector_len_1 * ortho_vector_1
        projected_side_vector_2 = projected_side_vector_len_2 * ortho_vector_2
        vec_12 = vec_1 + vec_2
        projected_side_vector_len_12_1 = vec_12.dot(ortho_vector_1)
        projected_side_vector_len_12_2 = vec_12.dot(ortho_vector_2)
        area_parallelogram = np.linalg.norm(np.cross(vec_1, vec_2))
        # make the designated axis an upright one
        poscar_dict['box_len_arr'][l_arr_row, l_arr_column_1] = 0 
        poscar_dict['box_len_arr'][l_arr_row, l_arr_column_2] = 0
        poscar_dict['l_arr'][l_arr_row, l_arr_column_1] = 0 
        poscar_dict['l_arr'][l_arr_row, l_arr_column_2] = 0
        #Whether a point is inside the parallelogram composed of vec_1 and vec_2. point_1, point_2, point_3, point_4 are the four vertices of the parallelogram composed of vec_1 and vec_2. If the area of four triangles equals to the area of the parallelogram, then the point is inside the parallelogram..
        point_1 = np.array([0, 0])
        point_2 = np.array([vec_1[l_arr_column_1], vec_1[l_arr_column_2]])
        point_3 = np.array([vec_2[l_arr_column_1], vec_2[l_arr_column_2]])
        point_4 = np.array([vec_12[l_arr_column_1], vec_12[l_arr_column_2]])
        delta_area = 0.01
        test_projection = False
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            temp_pos_1 = poscar_dict['pos_arr'][i_atom_indx, pos_arr_column_1]
            temp_pos_2 = poscar_dict['pos_arr'][i_atom_indx, pos_arr_column_2]
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
                shift_arr_1 = np.array([vec_1[l_arr_column_1], vec_1[l_arr_column_2]])
                shift_arr_2 = np.array([vec_2[l_arr_column_1], vec_2[l_arr_column_2]])
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
                            poscar_dict['pos_arr'][i_atom_indx, pos_arr_column_1] =  temp_pos_1 - n1 * shift_arr_1[0] - n2 * shift_arr_2[0]
                            poscar_dict['pos_arr'][i_atom_indx, pos_arr_column_2] =  temp_pos_2 - n1 * shift_arr_1[1] - n2 * shift_arr_2[1]
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
    basis_vector_dict = funcs.basis_vector_info(vec_a, vec_b, vec_c)
    for i_atom_indx in range(0, poscar_dict['n_atoms']):
        poscar_dict['pos_arr'][i_atom_indx,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    return 0

def transform(poscar_file_path, r_matrix_arr = None, t_matrix_arr = None, apply_to = 'atom', output_poscar_file_path = None):
    '''
    Transform by rotation and transformation operations.
    r_matrix_arr: Rotation matrix
    t_matrix_arr: translation matrix
    apply_to: 'box' or 'atom'
    '''
    import os
    import numpy as np
    from .. import funcs
    from . import vasp_read
    from . import vasp_write
    
    poscar_file_path = os.path.abspath(poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)

    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path in [None, 'None', 'none']:
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_transform')
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
    if apply_to in ['atom']:
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            atom_coord_vec = poscar_dict['pos_arr'][i_atom_indx,3:6]
            atom_coord_vec_trans = r_matrix_arr.dot(atom_coord_vec) + np.sum(t_matrix_arr, axis = 0)
            poscar_dict['pos_arr'][i_atom_indx,3:6] = atom_coord_vec_trans[0:3]
            poscar_dict['pos_arr'][i_atom_indx,0:3] = poscar_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][i_atom_indx,3:6])

    elif apply_to in ['box']:
        poscar_dict['box_len_arr'][0,:] = r_matrix_arr.dot(poscar_dict['box_len_arr'][0,:]) + t_matrix_arr[0,:]    
        poscar_dict['box_len_arr'][1,:] = r_matrix_arr.dot(poscar_dict['box_len_arr'][1,:]) + t_matrix_arr[1,:]
        poscar_dict['box_len_arr'][2,:] = r_matrix_arr.dot(poscar_dict['box_len_arr'][2,:]) + t_matrix_arr[2,:]
        vec_a = poscar_dict['box_len_arr'][0, :] * poscar_dict['uni_scale_fac'] 
        vec_b = poscar_dict['box_len_arr'][1, :] * poscar_dict['uni_scale_fac']
        vec_c = poscar_dict['box_len_arr'][2, :] * poscar_dict['uni_scale_fac']
        basis_vector_dict = funcs.basis_vector_info(vec_a, vec_b, vec_c)
        for i_atom_indx in range(0, poscar_dict['n_atoms']):
            poscar_dict['pos_arr'][i_atom_indx,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    return 0

def create_supercell(poscar_file_path, num_xyz_list = [1,1,1], output_poscar_file_path = None):
    '''
    Create supercell by expanding the cell in the direction of a, b, and c basis vectors.
    num_xyz_list = [n_x, n_y, n_z], nx>=1, ny>=1, nz>=1
    '''    
    import os
    import numpy as np
    from .. import funcs
    from . import vasp_read
    from . import vasp_write
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
    vec_a = poscar_dict['box_len_arr'][0, :] * poscar_dict['uni_scale_fac'] 
    vec_b = poscar_dict['box_len_arr'][1, :] * poscar_dict['uni_scale_fac']
    vec_c = poscar_dict['box_len_arr'][2, :] * poscar_dict['uni_scale_fac']
    basis_vector_dict = funcs.basis_vector_info(vec_a, vec_b, vec_c)
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
    for i_atom_indx in range(0, poscar_dict['n_atoms']):
        poscar_dict['pos_arr'][i_atom_indx,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    return 0

def build_vacuum_layer(poscar_file_path, vacuum_layer_direction = 'z', vacuum_layer_width = 20, threshold_vacuum_layer_width = 5.98, output_poscar_file_path = None):
    '''
    Build vacuum layer
    The vacuum_layer_width is defined as the difference between the maximum and minimum porjected(onto x, y, or z axis) atom coordinates.
    vacuum_layer_direction: Designate the direction to add the vacuum layer
    threshold_vacuum_layer_width: This is the threshold distance for the vacuum layer separation, if layer separation is larger than this value, then it is denoted as vacuum layer. Otherwise, it is not a vacuum layer. Default value of threshold_vacuum_layer_width is (1.15 * 2 * (CSD_covalent_radius of Fr)) = 1.15 * 2 * 2.60 = 5.98 Angstrom
    '''
    import os
    import math
    import numpy as np
    from .. import funcs
    from . import vasp_read
    from . import vasp_write

    poscar_file_path = os.path.abspath(poscar_file_path)
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path in [None, 'None', 'none']:
        #output_poscar_file_path = os.path.join(fpath, 'POSCAR_vacuum_' + str(vacuum_layer_width))
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_vacuum')
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    #orientation(poscar_file_path = poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = output_poscar_file_path)
    #poscar_dict = vasp_read.read_poscar(output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)

    if vacuum_layer_width <= threshold_vacuum_layer_width:
        print('WARNING #2012301433 (from build_vacuum_layer): The vacuum layer width is too small, it is recommended to set it to a value larger than ' + str(threshold_vacuum_layer_width) + ' Angstrom.')

    ##############################################################
    # initialization
    ##############################################################
    if vacuum_layer_direction in ['x', 'X']:
        l_arr_row = 0
        side_vector = poscar_dict['vec_a']
        side_vector_len = poscar_dict['len_vec_a']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([1, 0, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 3
        nx = math.ceil((vacuum_layer_width + box_len_ortho) / box_len_ortho)
        ny = 1
        nz = 1
    if vacuum_layer_direction in ['y', 'Y']:
        l_arr_row = 1
        side_vector = poscar_dict['vec_b']
        side_vector_len = poscar_dict['len_vec_b']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 1, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 4
        nx = 1
        ny = math.ceil((vacuum_layer_width + box_len_ortho) / box_len_ortho)
        nz = 1
    if vacuum_layer_direction in ['z', 'Z']:
        l_arr_row = 2
        side_vector = poscar_dict['vec_c']
        side_vector_len = poscar_dict['len_vec_c']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 0, 1])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 5
        nx = 1
        ny = 1
        nz = math.ceil((vacuum_layer_width + box_len_ortho) / box_len_ortho)
    #get crytal length
    #crystal_len = abs(np.max((poscar_dict['pos_arr'][:,3:6] - side_vector)[:, l_arr_row]) - np.min(poscar_dict['pos_arr'][:,l_arr_row]))
    crystal_len = np.max(poscar_dict['pos_arr'][:, pos_arr_column]) - np.min(poscar_dict['pos_arr'][:,pos_arr_column])
    ##############################################################
    # Check whether there exist vacuum layer in the original model
    ##############################################################
    layer_separation = box_len_ortho - crystal_len
    if layer_separation <= threshold_vacuum_layer_width:
        vacuum_layer_exist = False
    elif layer_separation > threshold_vacuum_layer_width:
        vacuum_layer_exist = True

    ##############################################################
    # Displace the atoms to the bottom of the axis
    ##############################################################
    align(poscar_file_path = poscar_file_path, align_direction = vacuum_layer_direction, mode = 'bottom', output_poscar_file_path = output_poscar_file_path)
    ##############################################################
    # Build vacuum layer
    ##############################################################
    #Determine the added delta_l_arr
    ##poscar_dict = vasp_read.read_poscar(output_poscar_file_path)
    delta_len = vacuum_layer_width - layer_separation
    ##if vacuum_layer_exist == False:
    ##    delta_len = vacuum_layer_width - layer_separation
    ##elif vacuum_layer_exist == True:    
    ##    if layer_separation <= vacuum_layer_width:
    ##        delta_len = vacuum_layer_width - layer_separation
    ##    elif layer_separation > vacuum_layer_width:
    ##        delta_len = vacuum_layer_width - layer_separation
    delta_l_arr = delta_len / cos_angle * unit_vector
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

def exfoliation(poscar_file_path, num_layers = 1, vacuum_layer_width = 20, exfoliation_direction = 'z', layer_dist_tolerance = 0.664, criteria_list = ['bonding', 'periodicity'], align_center = False, align_upright = True, output_poscar_file_path = None):
    '''
    Exfoliate 2d material
    num_layers: the number of layers left after exfoliation.
    - layer_dist_tolerance: This is the spacial resolution for the layer distance, for example z=0.24 and z=0.25 are considered to be in the same position if layer_dist_tolerance=0.01. Default value of layer_dist_tolerance is (1.15 * 2 * (CSD_covalent_radius of He)) + arbitrary_delta = 1.15 * 2 * 0.28 + 0.02 = 0.664 Angstrom
    align_center: Align the atoms to the center of the model along the exfoliation direction
    align_upright: Align the exfoliation direction axis as an upright one (i.e. not a tilted one)
    criteria_list: ['bonding'], ['periodicity'], or ['bonding', 'periodicity'], the criteria in the list will be adopted. For ['bonding', 'periodicity'], this means a combination of bonding and periodicity (first check bonding, then check periodicity) will be adopted. If bonding criterion worked, then periodicity criterion will be skipped; If bonding criterion failed, then periodicity criterion will be adopted.
    If the following error occurs, please turn the value of delta smaller.
            poscar_dict['fix_arr'][i_atom,2] + ' ' + '\n')
            IndexError: index 60 is out of bounds for axis 0 with size 60
    '''
    import os
    import numpy as np
    import math
    from .. import funcs
    from . import vasp_read
    from . import vasp_write
    poscar_file_path = os.path.abspath(poscar_file_path)
    fpath, fname = os.path.split(poscar_file_path)
    if output_poscar_file_path in [None, 'None', 'none']:
        output_poscar_file_path = os.path.join(fpath, 'POSCAR_exfoliation')
    else:
        output_poscar_file_path = os.path.abspath(output_poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)

    recommended_layer_dist_tolerance = 0.664 
    if layer_dist_tolerance < recommended_layer_dist_tolerance:
        print('WARNING #2012301056 (from vasp_build.exfoliation): layer_dist_tolerance is smaller than ' + str(recommended_layer_dist_tolerance) + '(Default value of layer_dist_tolerance is (1.15 * 2 * (CSD_covalent_radius of He)) + arbitrary_delta = 1.15 * 2 * 0.28 + 0.02 = 0.664 Angstrom), It is recommended to be larger than this value.')
    #######################################
    # Initialize
    #######################################
    #delta: This is the spacial resolution, for example z=0.24 and z=0.25 are considered to be in the same position if delta=0.01.
    delta = 0.05

    if exfoliation_direction in ['x', 'X']:
        number_order_text = 'yzx'
        l_arr_row = 0
        side_vector = poscar_dict['vec_a']
        side_vector_len = poscar_dict['len_vec_a']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([1, 0, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 3
        nx = math.ceil((vacuum_layer_width + box_len_ortho) / box_len_ortho)
        ny = 1
        nz = 1
        if num_layers >= nx:
            nx = num_layers + 1
    if exfoliation_direction in ['y', 'Y']:
        number_order_text = 'zxy'
        l_arr_row = 1
        side_vector = poscar_dict['vec_b']
        side_vector_len = poscar_dict['len_vec_b']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 1, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 4
        nx = 1
        ny = math.ceil((vacuum_layer_width + box_len_ortho) / box_len_ortho)
        nz = 1
        if num_layers >= ny:
            ny = num_layers + 1
    if exfoliation_direction in ['z', 'Z']:
        number_order_text = 'xyz'
        l_arr_row = 2
        side_vector = poscar_dict['vec_c']
        side_vector_len = poscar_dict['len_vec_c']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 0, 1])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 5
        nx = 1
        ny = 1
        nz = math.ceil((vacuum_layer_width + box_len_ortho) / box_len_ortho)
        if num_layers >= nz:
            nz = num_layers + 1

    #######################################
    # Reduce to primitive cell
    #######################################
    pass
    ############################################################
    # Reorientate the model
    ############################################################
    orientation(poscar_file_path = poscar_file_path, mode = 'a along x, b in xy', output_poscar_file_path = output_poscar_file_path)
    ############################################################
    # Expand model in a certain direction by creating supercell
    ############################################################
    create_supercell(poscar_file_path = output_poscar_file_path,
        num_xyz_list = [nx, ny, nz],
        output_poscar_file_path = output_poscar_file_path,
        ) 
    ##create_supercell(poscar_file_path = poscar_file_path,
    ##    num_xyz_list = [nx, ny, nz],
    ##    output_poscar_file_path = output_poscar_file_path,
    ##    ) 
    #######################################
    # Creating layer property array
    #######################################
    # group layers according to different properties and adding tags
    layer_dict = {}
    layer_dict['num_layers'] = 0 
    layer_dict['layer_coord_list'] = [None]
    layer_dict['indx_list'] = [None]
    layer_dict['atoms_indx_list'] =  [list()]
    #bond_status = 0: no bond with the adjacent two layers
    #bond_status = -1: forms bond only with one side (left or bottom) of the two layers
    #bond_status = 1: forms bond only with one side (right or top) of the two layers
    #bond_status = 2: forms bond with both the adjacent two layers
    layer_dict['bond_status_list'] = [None]
    layer_dict['periodicity_status_list'] = [None]
      
    poscar_dict = vasp_read.read_poscar(output_poscar_file_path)
    poscar_dict['atom_number_ortho_' + number_order_text] = funcs.atom_number_ortho(
        atom_key_arr = poscar_dict['atom_key_arr'],
        pos_arr_cartesian= poscar_dict['pos_arr'][:,3:6],
        delta = delta,
        number_order = number_order_text
        )
    indx_0 = np.argwhere(poscar_dict['atom_number_ortho_' + number_order_text] == 1)
    # temp_pos can be used as a representitive position of each layer
    temp_pos_0 = poscar_dict['pos_arr'][indx_0, pos_arr_column]
    temp_pos = temp_pos_0
    i_layer_indx = 0
    bond_breaking = True
    shift_vec = poscar_dict['l_arr'][l_arr_row,:]
    for i in range(1, poscar_dict['n_atoms'] + 1):
        i_atom_indx = int(np.argwhere(poscar_dict['atom_number_ortho_' + number_order_text] == i))
        i_atom_pos = poscar_dict['pos_arr'][i_atom_indx, pos_arr_column]
        i_atom_pos_shift = (poscar_dict['pos_arr'][i_atom_indx, 3:6] - shift_vec)[pos_arr_column -3]
        if abs(i_atom_pos_shift - temp_pos_0) <= delta:
            # for the atoms near the topmost boundary, check whether it is in the same layer as the bottom layer
            i_layer_indx = 0
            layer_dict['atoms_indx_list'][i_layer_indx].append(i_atom_indx)
        elif  abs(i_atom_pos - temp_pos) <= layer_dist_tolerance:
            pass
        else:
            i_layer_indx = i_layer_indx + 1
            layer_dict['indx_list'].append(None)
            layer_dict['layer_coord_list'].append(None)
            layer_dict['atoms_indx_list'].append(list())
            layer_dict['bond_status_list'].append(None)
            layer_dict['periodicity_status_list'].append(None)
        layer_dict['indx_list'][i_layer_indx] = i_layer_indx
        layer_dict['layer_coord_list'][i_layer_indx] = i_atom_pos
        layer_dict['atoms_indx_list'][i_layer_indx].append(i_atom_indx)
        layer_dict['bond_status_list'][i_layer_indx] = 0
        layer_dict['periodicity_status_list'][i_layer_indx] = None
        temp_pos = i_atom_pos
    layer_dict['num_layers'] = len(layer_dict['indx_list'])
    # check bonding status of atoms among the layers
    for i_layer_indx in range(0, layer_dict['num_layers']):
        # j and k layers are the adjacent layers for layer i
        j_layer_indx = i_layer_indx - 1
        k_layer_indx = i_layer_indx + 1
        # bottom layer
        if i_layer_indx == 0:
            j_layer_indx = layer_dict['indx_list'][-1]
            k_layer_indx = i_layer_indx + 1
        # top layer
        if i_layer_indx == layer_dict['indx_list'][-1]:
            j_layer_indx = i_layer_indx - 1
            k_layer_indx = 0
        i_layer_atom_indx_list = layer_dict['atoms_indx_list'][i_layer_indx]
        j_layer_atom_indx_list = layer_dict['atoms_indx_list'][j_layer_indx]
        k_layer_atom_indx_list = layer_dict['atoms_indx_list'][k_layer_indx]
        ij_bonding_list = []
        ik_bonding_list = []
        for i_atom_indx in i_layer_atom_indx_list:
            i_elmt = poscar_dict['atom_species_arr'][i_atom_indx]
            i_atom_pos = poscar_dict['pos_arr'][i_atom_indx, 3:6]
            for j_atom_indx in j_layer_atom_indx_list:
                j_elmt = poscar_dict['atom_species_arr'][j_atom_indx]
                j_atom_pos = poscar_dict['pos_arr'][j_atom_indx, 3:6]
                if i_layer_indx == 0 and j_layer_indx == layer_dict['indx_list'][-1]:
                    # Due to periodic boundary condition, tanslate the boundary layer with the lattice vector
                    if (j_atom_pos - i_atom_pos).dot(shift_vec) > 0:
                        j_atom_pos = j_atom_pos - shift_vec
                    else:
                        j_atom_pos = j_atom_pos + shift_vec
                bonding = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos) 
                ij_bonding_list.append(bonding)
            for k_atom_indx in k_layer_atom_indx_list:
                k_elmt = poscar_dict['atom_species_arr'][k_atom_indx]
                k_atom_pos = poscar_dict['pos_arr'][k_atom_indx, 3:6]
                if i_layer_indx == layer_dict['indx_list'][-1] and k_layer_indx == 0:
                    # Due to periodic boundary condition, tanslate the boundary layer with the lattice vector
                    if (j_atom_pos - i_atom_pos).dot(shift_vec) > 0:
                        k_atom_pos = k_atom_pos - shift_vec
                    else:
                        k_atom_pos = k_atom_pos + shift_vec
                bonding = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos) 
                ik_bonding_list.append(bonding)

        if True in ij_bonding_list and True in ik_bonding_list:
            layer_dict['bond_status_list'][i_layer_indx] = 2
        elif True in ij_bonding_list and True not in ik_bonding_list:
            layer_dict['bond_status_list'][i_layer_indx] = 1
        elif True not in ij_bonding_list and True in ik_bonding_list:
            layer_dict['bond_status_list'][i_layer_indx] = -1
        elif True not in ij_bonding_list and True not in ik_bonding_list:
            layer_dict['bond_status_list'][i_layer_indx] = 0
    ##print('bons_status_list = ', layer_dict['bond_status_list'])
    # check periodicity of layers
    if exfoliation_direction in ['x', 'X']:
        l_arr_row = 0
        side_vector = poscar_dict['vec_a']
        side_vector_len = poscar_dict['len_vec_a']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([1, 0, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 3
    if exfoliation_direction in ['y', 'Y']:
        l_arr_row = 1
        side_vector = poscar_dict['vec_b']
        side_vector_len = poscar_dict['len_vec_b']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 1, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 4
    if exfoliation_direction in ['z', 'Z']:
        l_arr_row = 2
        side_vector = poscar_dict['vec_c']
        side_vector_len = poscar_dict['len_vec_c']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 0, 1])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
    atom_dist_tolerance = 0.05 
    periodicity_number = 0
    layer_dict['periodicity_status_list'][0] = periodicity_number 
    for i_layer_indx in range(0, layer_dict['num_layers'] - 1):
        if layer_dict['periodicity_status_list'][i_layer_indx] == None:
            periodicity_number = periodicity_number + 1
            layer_dict['periodicity_status_list'][i_layer_indx] = periodicity_number 
        original_periodicity_number = periodicity_number
        i_layer_atom_indx_list = layer_dict['atoms_indx_list'][i_layer_indx] 
        i_layer_coord = layer_dict['layer_coord_list'][i_layer_indx]
        for j_layer_indx in range(i_layer_indx + 1, layer_dict['num_layers']):
            if layer_dict['periodicity_status_list'][j_layer_indx] != None:
                continue
            j_layer_atom_indx_list = layer_dict['atoms_indx_list'][j_layer_indx] 
            j_layer_coord = layer_dict['layer_coord_list'][j_layer_indx]

            ortho_layer_dist = abs(j_layer_coord - i_layer_coord)
            displace_vec = ortho_layer_dist / cos_angle * unit_vector
            identical_atoms_list1 = []
            for i_atom_indx in i_layer_atom_indx_list:
                identical_atoms_list2 = []
                i_elmt = poscar_dict['atom_species_arr'][i_atom_indx]
                i_atom_pos = poscar_dict['pos_arr'][i_atom_indx, 3:6]
                for j_atom_indx in j_layer_atom_indx_list:
                    j_elmt = poscar_dict['atom_species_arr'][j_atom_indx]
                    j_atom_pos = poscar_dict['pos_arr'][j_atom_indx, 3:6]
                    if i_elmt == j_elmt and np.linalg.norm(i_atom_pos + displace_vec - j_atom_pos) <= atom_dist_tolerance:
                        identical_atoms_list2.append(True)
                    elif i_elmt != j_elmt and np.linalg.norm(i_atom_pos - j_atom_pos) > atom_dist_tolerance:
                        identical_atoms_list2.append(False)
                if True in identical_atoms_list2:
                    identical_atoms_list1.append(True)
                else:
                    identical_atoms_list1.append(False)
            if False in identical_atoms_list1:                 
                ##periodicity_number = periodicity_number + 1
                ##layer_dict['periodicity_status_list'][j_layer_indx] = periodicity_number 
                pass
            else:
                layer_dict['periodicity_status_list'][j_layer_indx] = original_periodicity_number 
    ##print('periodicity_status_list=',layer_dict['periodicity_status_list'])
    ###########################################################
    # Screen the region of intrest (ROI) by certain criteria_list
    ###########################################################
    #For ['bonding', 'periodicity'], this means a combination of bonding and periodicity (first check bonding, then check periodicity) will be adopted. If bonding criterion worked, then periodicity criterion will be skipped; If bonding criterion failed, then periodicity criterion will be adopted.
    #######################
    # bonding criterion
    #######################
    # If the criteria_list is 'bonding', then the layers will be distinguished by the bonding status. If all the atoms in a layer have no bonding with the atoms in the adjacent layer, then they are said to be the boundary layer. This option is suitable for vdW materials
    if 'bonding' in criteria_list and -1 in layer_dict['bond_status_list'] and 1 in layer_dict['bond_status_list']:
        num_elmt_list = []
        elmt_species_list = []
        atom_indx_list = []
        for i_elmt in poscar_dict['elmt_species_arr']:
            i_num_elmt = 0
            stacked_layer_number = 0
            for i_layer_indx in range(0, layer_dict['num_layers']):
               if layer_dict['bond_status_list'][i_layer_indx] == -1:
                   stacked_layer_number = stacked_layer_number + 1
               if stacked_layer_number <= num_layers and stacked_layer_number != 0:
                   for i_atom_indx in layer_dict['atoms_indx_list'][i_layer_indx]:
                       if i_elmt == poscar_dict['atom_species_arr'][i_atom_indx]:
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
        new_atom_indx = 0
        for i_atom_indx in atom_indx_list:
            poscar_dict['pos_arr'][new_atom_indx, 3:6] = pos_arr_old[i_atom_indx, 3:6]
            new_atom_indx = new_atom_indx + 1
    else:
        pass
        #print('WARNING #2012301623 (from exfoliation): Cannot exfoliate layers according to the bonding state. All the layers are bonded: ' + poscar_file_path)

        #######################
        # periodicity criterion
        #######################
        # If the criteria_list is 'periodicity', then the layers will be distinguished by the periodicity of layers, i.e. ABCABCABC. This option also suitable for bulk materials.
        if 'periodicity' in criteria_list and None not in layer_dict['bond_status_list']:
            num_elmt_list = []
            elmt_species_list = []
            atom_indx_list = []
            indices_0 = [i for i, x in enumerate(layer_dict['periodicity_status_list']) if x == 0]
            indices_1 = [i for i, x in enumerate(layer_dict['periodicity_status_list']) if x == 1]
            if len(indices_0) == 1:
                periodicity_indx = 0
            if len(indices_0) > 1 and len(indices_1) > 1:
                delta_indices_0 = indices_0[1] - indices_0[0]
                delta_indices_1 = indices_1[1] - indices_1[0]
                if delta_indices_0 > delta_indices_1:
                    periodicity_indx = 1
                else:
                    periodicity_indx = 0
                
            for i_elmt in poscar_dict['elmt_species_arr']:
                i_num_elmt = 0
                stacked_layer_number = 0
                for i_layer_indx in range(0, layer_dict['num_layers']):
                   if layer_dict['periodicity_status_list'][i_layer_indx] == periodicity_indx:
                       stacked_layer_number = stacked_layer_number + 1
                   if stacked_layer_number <= num_layers and stacked_layer_number != 0:
                       for i_atom_indx in layer_dict['atoms_indx_list'][i_layer_indx]:
                           if i_elmt == poscar_dict['atom_species_arr'][i_atom_indx]:
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
            new_atom_indx = 0
            for i_atom_indx in atom_indx_list:
                poscar_dict['pos_arr'][new_atom_indx, 3:6] = pos_arr_old[i_atom_indx, 3:6]
                new_atom_indx = new_atom_indx + 1
        else:
            if len(criteria_list) == 1 and 'bonding' in criteria_list:
                print('WARNING #2012301623 (from exfoliation): Cannot exfoliate layers according to the bonding state. All the layers are bonded: ' + poscar_file_path)
            elif len(criteria_list) == 1 and 'periodicity' in criteria_list:
                print('WARNING #2101041102 (from exfoliation): Cannot exfoliate layers according to periodicity of layers: ' + poscar_file_path)
            elif 'bonding' in criteria_list and 'periodicity' in criteria_list:
                print('WARNING #2101041311 (from exfoliation): Cannot exfoliate layers according to the bonding state and periodicity of layers: ' + poscar_file_path)
    vec_a = poscar_dict['box_len_arr'][0, :] * poscar_dict['uni_scale_fac'] 
    vec_b = poscar_dict['box_len_arr'][1, :] * poscar_dict['uni_scale_fac']
    vec_c = poscar_dict['box_len_arr'][2, :] * poscar_dict['uni_scale_fac']
    basis_vector_dict = funcs.basis_vector_info(vec_a, vec_b, vec_c)
    for i_atom_indx in range(0, poscar_dict['n_atoms']):
        poscar_dict['pos_arr'][i_atom_indx,0:3] = basis_vector_dict['car2fra_matrix_arr'].dot(poscar_dict['pos_arr'][i_atom_indx,3:6])
    ##print(output_poscar_file_path)
    vasp_write.write_poscar(output_poscar_file_path = output_poscar_file_path, poscar_dict = poscar_dict, coord_system = 'Direct')
    #######################################
    # Adjust the vacuum layer width
    #######################################
    build_vacuum_layer(poscar_file_path = output_poscar_file_path, vacuum_layer_direction = exfoliation_direction, vacuum_layer_width = vacuum_layer_width, threshold_vacuum_layer_width = 5.98, output_poscar_file_path = output_poscar_file_path)
    ##############################################################################
    # Align the atoms to the center of the model along the exfoliation direction
    ##############################################################################
    # Displace the atoms to the center of the axis
    if align_center == True:
        align(poscar_file_path = output_poscar_file_path, align_direction = exfoliation_direction, mode = 'center', output_poscar_file_path = output_poscar_file_path)
    ##############################################################################
    # Align the exfoliation direction axis as an upright one (i.e. not a tilted one)
    ##############################################################################
    if align_upright == True:
        align(poscar_file_path = output_poscar_file_path, align_direction = exfoliation_direction, mode = 'upright_box', output_poscar_file_path = output_poscar_file_path)
    return 0
