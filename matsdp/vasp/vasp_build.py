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
