# -*- coding: utf-8 -*-
def substitution(substitution_list_file, poscar_dir):
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
            Ni 244 W
            
            2
            Ni 244 Re
            Al 12 Re

            ...
        sysname is L$(line_number)_composition_D$(Duplicate)

    - Args
     * substitution_list_file: String format. It specifies the directory of the .subst file (substitution list file)
     * poscar_dir: String format. The directory of the POSCAR file which is to be subsituted. It can either be full path or relative path.
    '''
       
    import os
    import sys
    import numpy as np
    from .. import funcs
    from . import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)

    poscar_dir = os.path.abspath(poscar_dir)
    substitution_list_file_abs_path = os.path.abspath(substitution_list_file)
    substitution_list_file_path,substitution_list_file = os.path.split(substitution_list_file_abs_path)
    if os.path.isfile(substitution_list_file_abs_path) == False:
        funcs.write_log(logfile,substitution_list_file + " doesn't exist, please check current folder. The file extension of substitution_list_file should be .subst")
    else:
        subst_file_name, subst_file_extension = os.path.splitext(substitution_list_file)
        if subst_file_extension != '.subst':
            funcs.write_log(logfile,substitution_list_file + " detected. But the file extension should be .subst ")
    if os.path.isfile(poscar_dir) == False:
        funcs.write_log(logfile,poscar_dir + " doesn't exist, please check current folder")
    
    # Extract information from the input POSCAR file
    poscar_dict = vasp_read.read_poscar(poscar_dir)
    n_atoms = np.sum(poscar_dict['elmt_num_arr']) 
    with open(substitution_list_file_abs_path) as s_lines:
        s_line = s_lines.readlines();
        s_n_lines = len(s_line)

    #Recognize atom index
    composition_arr = np.array([],dtype = np.str)
    composition_arr_count = np.array([],dtype = np.int)
    sysname_file = substitution_list_file_path + '/' + subst_file_name + '.sysname'
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
        val_indx = np.argwhere(elmt_species_mod == 'Va')
        elmt_species_mod_remove_va = np.delete(elmt_species_mod, val_indx)
        n_elmt_mod_RemoveVa = np.delete(n_elmt_mod, val_indx)
        #Export POSCAR file for each system
        models_path = output_dir + '/' + subst_file_name + '/' + sysname
        funcs.mkdir(models_path)
        isExists = os.path.exists(models_path)
        if not isExists:
            funcs.write_log(logfile,'WARNING:' + models_path + ' does not exist!')
        poscar_out = models_path + '/POSCAR'
        with open(poscar_out,'w') as fileobject, open(poscar_dir,'r') as pline:
            pline = pline.readlines()
            fileobject.write(pline[0] +
                             pline[1] +
                             pline[2] +
                             pline[3] +
                             pline[4] +
                             str(" ".join(elmt_species_mod_remove_va)) + '\n' +
                             str(n_elmt_mod_RemoveVa).strip('[').strip(']') + '\n')
            if poscar_dict['slet_dyn_on'] == True:
                fileobject.write(pline[7] +
                                 pline[8])
            else:
                fileobject.write(pline[7])
        with open(poscar_out,'a') as fileobject:
            if poscar_dict['slet_dyn_on'] == True:
                for i in range(0,n_atoms):
                    if atom_species_name_sort_arr[i] == 'Va':
                        continue
                    fileobject.write(str(coord_subst_arr[i][0])+' '+str(coord_subst_arr[i][1])+' '+str(coord_subst_arr[i][2])+' '+str(fix_subst_arr[i][0])+' '+str(fix_subst_arr[i][1])+' '+str(fix_subst_arr[i][2])+'\n')
            else:
                for i in range(0,n_atoms):
                    if atom_species_name_sort_arr[i] == 'Va':
                        continue
                    fileobject.write(str(coord_subst_arr[i][0])+' '+str(coord_subst_arr[i][1])+' '+str(coord_subst_arr[i][2])+'\n')
        with open(sysname_file,'a') as sysnamefileobject:
            sysnamefileobject.write(sysname + '\n')
        
        i_line = i_line + n_subst + 2
        model_number += 1
    
    funcs.write_log(
        logfile,
        'vasp_build.substitution(' + '\n' +
        '    substitution_list_file=' + 'r\'' + str(substitution_list_file_abs_path) + '\'' + ',\n' +
        '    poscar_dir='  + 'r\'' +str(poscar_dir) + '\'' + ')\n' +
        '###############################\n')
    return 0

def rep_elmt(substitution_list_file, poscar_dir, old_elmt, elmt_group):
    '''
    - Descriptions
     * replace element in the .subst file by specific elements and generate corresponding models
    - Args
     * substitution_list_file: String format. It specifies the directory of the *.subst file (substitution list file).
     * poscar_dir: String format. The directory of the POSCAR file. It can either be full path or relative path
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
    poscar_dir = os.path.abspath(poscar_dir)
    substitution_list_file_abs_path = os.path.abspath(substitution_list_file)
    subst_file_name, subst_file_extension = os.path.splitext(substitution_list_file)
    sysname_file =subst_file_name + '.sysname'
    str1 = old_elmt
    if os.path.isfile('Tempsysname.dat'):
        os.remove('Tempsysname.dat')
    for i_elmt in elmt_group:
        funcs.replace_file_content(substitution_list_file,str1,i_elmt)
        str1 = i_elmt
        substitution(substitution_list_file, poscar_dir)
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
        '    poscar_dir='  + 'r\'' + str(poscar_dir) + '\'' + ',\n' +
        '    old_elmt=' + '\'' + str(old_elmt) + '\'' + ',\n' +
        '    elmt_group=' + str(elmt_group) + ')\n' +
        '###############################\n')
    return 0

def selection_sphere(poscar_dir, origin_atom_name, radius = 7.0, include_mirror_atoms = True, output_file_name = 'example'):
    '''
    Descriptions:
    This module is used to select atoms in a sphere that is centered around an arbitrary atom from the POSCAR file.
    The selected atom coordinates can also be written to the *.incar file of the DVM program.
    poscar_dir: String format. The directory of the POSCAR file. It can either be full path or relative path.
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
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)
    selection_dir_name = 'atom_selection_sphere'
    funcs.mkdir(output_dir + '/' + selection_dir_name + '/' + output_file_name + '/')

    poscar_dir = os.path.abspath(poscar_dir)

    # the *.vasp file which contains the selected atoms
    vasp_file = output_dir + '/' + selection_dir_name + '/' + output_file_name + '/' + output_file_name + '.vasp'
    vasp_temp_file1 = output_file_name + 'vasp_1.temp'
    vasp_temp_file2 = output_file_name + 'vasp_2.temp'
    # the *.xyz file which contains the selected atoms
    xyz_file = output_dir + '/' + selection_dir_name + '/' + output_file_name + '/' + output_file_name + '.xyz'
    xyz_temp_file1 = output_file_name + '_xyz_1.temp'
    xyz_temp_file2 = output_file_name + '_xyz_2.temp'
    # the *.incar file which contains the selected atoms
    dvm_incar_file = output_dir + '/' + selection_dir_name + '/' + output_file_name + '/' + output_file_name + '.incar'
    dvm_incar_temp_file1 = output_file_name + '_1.temp'
    dvm_incar_temp_file2 = output_file_name + '_2.temp'

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
                if include_mirror_atoms == True:
                    # if the mirror atoms are included, then the atom coordinates will be shifted
                    f_vasp_temp1.write(
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 0] - xmin)) + ' ' +
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 1] - ymin)) + ' ' +
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 2] - zmin)) + ' ' +
                        str(expanded_atom_species_arr[i_atom])+ ' ' +
                        expanded_atom_species_arr[i_atom] + str(new_elmt_count_list[new_elmt_species_list.index(expanded_atom_species_arr[i_atom])]) + ' ' +
                        '    ' + str(expanded_atomname_arr[i_atom]) + ' ' + '\n'
                        )
                else: 
                    f_vasp_temp1.write(
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 0])) + ' ' +
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 1])) + ' ' +
                        str('{:.6f}'.format(expanded_pos_arr[i_atom, 2])) + ' ' +
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
        if include_mirror_atoms == True:
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
        '    poscar_dir=' + 'r\'' + str(poscar_dir) + '\'' + ',\n' +
        '    origin_atom_name=\'' + str(origin_atom_name) + '\',\n' + 
        '    radius=' + str(radius) + ',\n' +
        '    include_mirror_atoms=' + str(include_mirror_atoms) + ',\n' +
        '    output_file_name=\'' + str(output_file_name) + '\')\n')
    return 0 
