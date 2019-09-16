# -*- coding: utf-8 -*-
def read_poscar(poscar_dir):
    '''
    Description:
        Read information from POSCAR file
    Args:
        @.poscar_dir: String format. The directory of the POSCAR file. It can either be full path or relative path
    Return:
        poscar_dict
        l_arr: Numpy float array(np.float). Dimension is 3*3
        elmtindxArr': Numpy integer array(np.int). Element index. Start from 0
        elmt_num_arr': Numpy integer array(np.int). Number of atoms for each element species
        ElmtSpeciesArr': Numpy string array. element name in the periodic for each element species
        ElmtStartindxArr': Numpy integer array(np.int). staring index for each element species in the atom_indxArr
        CoordSystem': String type. String value of either 'Direct' or 'Cartesian'
        UniScaleFac': Float type. Universal scale factor(latice constant)
        SletDynOn': Logic type. Selective dynamics on or off? True or False
        N_Header': Number of header lines. i.e. the number of lines before atomic positions
    '''
    import os
    import numpy as np
    import sys
    import time
    from .. import funcs

    poscar_dir = os.path.abspath(poscar_dir)
    poscar_dict = {}
    with open(poscar_dir) as f:
        lines = f.readlines()
        poscar_dict['UniScaleFac'] = float(funcs.split_line(lines[1])[0])
        box_len_arr = np.array([0.0]*3*3,dtype = np.float)
        box_len_arr.shape = 3,3
        l_arr = np.array([0.0]*3*3,dtype = np.float)
        l_arr.shape = 3,3
        box_len_arr[0,0] = float(funcs.split_line(lines[2])[0]);box_len_arr[0,1] = float(funcs.split_line(lines[2])[1]);box_len_arr[0,2] = float(funcs.split_line(lines[2])[2])
        box_len_arr[1,0] = float(funcs.split_line(lines[3])[0]);box_len_arr[1,1] = float(funcs.split_line(lines[3])[1]);box_len_arr[1,2] = float(funcs.split_line(lines[3])[2])
        box_len_arr[2,0] = float(funcs.split_line(lines[4])[0]);box_len_arr[2,1] = float(funcs.split_line(lines[4])[1]);box_len_arr[2,2] = float(funcs.split_line(lines[4])[2])
        l_arr = box_len_arr * poscar_dict['UniScaleFac']
        #-First POS_ParamsDic['N_Header']  lines are header of the POSCAR file
        first_char = lines[7][0:1]
        if first_char == "S" or first_char == "s":
            poscar_dict['SletDynOn'] = True 
            poscar_dict['N_Header']  = 9  #This is the real Header line number in the POSCAR, which means that the header has POS_ParamsDic['N_Header']  lines.
            first_char = lines[poscar_dict['N_Header'] -1][0:1]
            if first_char == "D" or first_char == "d":
                poscar_dict['CoordSystem'] = "Direct"
            elif first_char == "C" or first_char == "c":
                poscar_dict['CoordSystem'] = "Cartesian"
        elif first_char == "D" or first_char == "d" or first_char == "C" or first_char == "c":
            poscar_dict['SletDynOn'] = False
            poscar_dict['N_Header'] =8
            if first_char == "D" or first_char == "d":
                poscar_dict['CoordSystem'] = "Direct"
            elif first_char == "C" or first_char == "c":
                poscar_dict['CoordSystem'] = "Cartesian"
        #-Get number of element species
        n_species = len(funcs.split_line(lines[5]))
        #-Number of atoms in the system
        poscar_dict['elmtindxArr'] = np.array([0]*n_species,dtype = np.int)
        poscar_dict['ElmtSpeciesArr'] = np.array(['--']*n_species,dtype = np.str)
        poscar_dict['elmt_num_arr'] = np.array([0]*n_species,dtype = np.int)
        for i in range(0,n_species):
            poscar_dict['elmtindxArr'][i] = i
            poscar_dict['ElmtSpeciesArr'][i] = funcs.split_line(lines[5])[i]
            poscar_dict['elmt_num_arr'][i] = funcs.split_line(lines[6])[i]
        n_atoms = sum(poscar_dict['elmt_num_arr'][:])
        #Element start index
        poscar_dict['ElmtStartindxArr'] = np.array([0]*n_species,dtype=np.int)
        for i in range(0,n_species):
            if i == 0:
                poscar_dict['ElmtStartindxArr'][i] = 1
            else:
                poscar_dict['ElmtStartindxArr'][i] = sum(poscar_dict['elmt_num_arr'][range(0,i)])+1

        #-Array atom_species_arr and atom coordinates
        atom_species_arr = np.array(['--']*n_atoms,dtype=np.str)
        atom_subindx_arr = np.array([0]*n_atoms,dtype=np.int)
        atomname_list = ['']*n_atoms
        indx = 0;
        for i in range(0,n_species):
            atom_species_arr[range(indx,indx+poscar_dict['elmt_num_arr'][i])] = [poscar_dict['ElmtSpeciesArr'][i]]*poscar_dict['elmt_num_arr'][i]
            for j in range(0,poscar_dict['elmt_num_arr'][i]):
                i_atom = poscar_dict['ElmtStartindxArr'][i] + j
                atom_subindx_arr[i_atom-1] = j + 1
                atomname_list[i_atom-1] = str(poscar_dict['ElmtSpeciesArr'][i]) + str(atom_subindx_arr[i_atom-1])
            indx += poscar_dict['elmt_num_arr'][i]
        
        #Atom coordinates
        #pos_arr is the coordinates of all the atoms, first three columns are in Direct coordinate, last columns are in Cartesian coordinate
        #coord_arr is the x, y, z coordinate in POSCAR, without the influence of lattice vector or scale factor
        coord_arr = np.array([0.0]*n_atoms*3,dtype = np.float)
        coord_arr.shape = n_atoms,3
        fix_arr = np.array(['T']*n_atoms*3,dtype = np.str)
        fix_arr.shape = n_atoms,3
        pos_arr = np.array([0.0]*n_atoms*6,dtype = np.float)
        pos_arr.shape = n_atoms,6
        L_Inv_Arr = np.linalg.inv(l_arr)
        for i in range(0,n_atoms):
            temp = funcs.split_line(lines[poscar_dict['N_Header'] +i])
            coord_arr[i,0] = float(temp[0])
            coord_arr[i,1] = float(temp[1])
            coord_arr[i,2] = float(temp[2])
            if first_char == "S" or first_char == "s":
                fix_arr[i,0] = temp[3]
                fix_arr[i,1] = temp[4]
                fix_arr[i,2] = temp[5]
            if  poscar_dict['CoordSystem'] == "Direct":
                pos_arr[i,0] = coord_arr[i,0]
                pos_arr[i,1] = coord_arr[i,1]
                pos_arr[i,2] = coord_arr[i,2]
                pos_arr[i,3] = np.dot(coord_arr[i,:], l_arr[:,[0]])
                pos_arr[i,4] = np.dot(coord_arr[i,:], l_arr[:,[1]])
                pos_arr[i,5] = np.dot(coord_arr[i,:], l_arr[:,[2]])
            elif poscar_dict['CoordSystem'] == "Cartesian":
                pos_arr[i,0] = np.dot(L_Inv_Arr[0,:], coord_arr[i,:])
                pos_arr[i,1] = np.dot(L_Inv_Arr[1,:], coord_arr[i,:])
                pos_arr[i,2] = np.dot(L_Inv_Arr[2,:], coord_arr[i,:])
                pos_arr[i,3] = coord_arr[i,0] * poscar_dict['UniScaleFac']
                pos_arr[i,4] = coord_arr[i,1] * poscar_dict['UniScaleFac']
                pos_arr[i,5] = coord_arr[i,2] * poscar_dict['UniScaleFac']
        poscar_dict['atom_species_arr'] = atom_species_arr
        poscar_dict['atom_subindx_arr'] = atom_subindx_arr
        poscar_dict['atomname_list'] = atomname_list
        poscar_dict['pos_arr'] = pos_arr
        poscar_dict['coord_arr'] = coord_arr
        poscar_dict['fix_arr'] = fix_arr
        poscar_dict['l_arr'] = l_arr
    return poscar_dict

def read_outcar(outcar_dir):
    '''
    Description:
        Read information from OUTCAR file
    Args:
        outcar_dir: String format. The directory of the OUTCAR file. It can either be full path or relative path
    Return:
        @.incar_params_dict
        @.outcar_params_dict
    '''
    import os
    import numpy as np
    from .. import funcs

    outcar_dir = os.path.abspath(outcar_dir)
    incar_params_dict = {}
    outcar_params_dict = {}
    incar_params_dict['LORBIT'] = int(funcs.split_line(funcs.grep('LORBIT',outcar_dir))[2])

    outcar_params_dict['e_fermi'] = float(funcs.split_line(funcs.grep('E-fermi',outcar_dir))[2])
    outcar_params_dict['XC(G=0)'] = float(funcs.split_line(funcs.grep('E-fermi',outcar_dir))[4])
    outcar_params_dict['LORBIT'] = float(funcs.split_line(funcs.grep('LORBIT',outcar_dir))[2])
    outcar_params_dict['ISPIN'] = int(funcs.split_line(funcs.grep('ISPIN',outcar_dir))[2])
    if outcar_params_dict['ISPIN'] == 1:
        outcar_params_dict['alpha+bet'] = float(funcs.split_line(funcs.grep('E-fermi',outcar_dir))[-1].strip(':'))
    else:
        outcar_params_dict['alpha+bet'] = float(funcs.split_line(funcs.grep('E-fermi',outcar_dir))[6].split(':')[-1])


    num_elmt_type = len(funcs.split_line(funcs.grep('ions per type',outcar_dir))) - 4
    num_elmt_type_list = funcs.split_line(funcs.grep('ions per type',outcar_dir))[4:]
    n_atoms = np.sum([int(item) for item in num_elmt_type_list])
    num_ionic_step = funcs.grep('TOTaxis_indic_ref_length-FORCE',outcar_dir).count('TOTaxis_indic_ref_length-FORCE')
    force_arr = np.array([None] * num_ionic_step * n_atoms * 6)
    force_arr.shape = num_ionic_step , n_atoms , 6
    outcar_params_dict['num_ionic_step'] = num_ionic_step

    with open(outcar_dir,'r') as f:
        line = f.readlines()
        iionic_step = 0
        for i in range(len(line)):
            if 'TOTaxis_indic_ref_length-FORCE' in line[i]:
                for i_atom in range(0,n_atoms):
                    i_atomLine = i + i_atom + 2
                    force_arr[iionic_step,i_atom,:] = [float(item) for item in funcs.split_line(line[i_atomLine])]
                iionic_step += 1
    outcar_params_dict['force_arr'] = force_arr
    return incar_params_dict, outcar_params_dict

def read_doscar(doscar_dir, atom_indx, save_dos_arr = False):
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
        @.poscar_dir: String format. The directory of the DOSCAR file. It can either be full path or relative path
        @.atom_indx: Integer format. The real atom index in the POSCAR. If there are N atoms then the atom indices are frim 1 to N. Note that atom_indx = 0 means to extract TDOS inoformation
        @save_dos_arr: logical value. Determine whether to save the dos_arr to a file or not.
    Return:
        @.dos_arr: NEDOS * NCol array type. It contains the density of states
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    doscar_dir = os.path.abspath(doscar_dir)
    workdir, dos_file = funcs.file_path_name(doscar_dir)
    #Also Required files are OUTCAR and POSCAR
    with open(doscar_dir,'r') as f:
        lines = f.readlines()
        n_atoms = int(funcs.split_line(lines[0])[0])
        NEDOS = int(funcs.split_line(lines[5])[2])
        e_fermi = float(funcs.split_line(lines[5])[3])
        funcs.write_log(logfile, 'n_atoms = ' + str(n_atoms) + ' NEDOS = ' + str(NEDOS) + ' e_fermi(Unmodified) = ' + str(e_fermi))
        #TDOS
        if atom_indx == 0:
            StartLine = 6
            EndLine = StartLine + NEDOS - 1
            NC = len(funcs.split_line(lines[StartLine]))
            dos_arr = np.array([0.000000]*NEDOS*NC,dtype = np.float)
            dos_arr.shape = NEDOS,NC
            for i in range(0,NEDOS):
                for j in range(0,NC):
                    dos_arr[i][j] = float(funcs.split_line(lines[StartLine + i])[j]) * funcs.logic_retn_val(j % 2 == 0 and j > 0, -1, 1)
            dos_file = workdir + '/TDOS.dat'
        #DOS
        if atom_indx != 0:
            StartLine = 6 + (NEDOS + 1) * atom_indx
            EndLine = StartLine + NEDOS - 1
            NC = len(funcs.split_line(lines[StartLine]))
            dos_arr = np.array([0.000000]*NEDOS*NC,dtype = np.float)
            dos_arr.shape = NEDOS,NC
            if NC == 10 or NC == 17:
                for i in range(0,NEDOS):
                    for j in range(0,NC):
                        dos_arr[i][j] = float(funcs.split_line(lines[StartLine + i])[j])
            elif NC == 19 or NC == 33:
                for i in range(0,NEDOS):
                    for j in range(0,NC):
                        dos_arr[i][j] = float(funcs.split_line(lines[StartLine + i])[j]) * funcs.logic_retn_val(j % 2 == 0 and j > 0, -1, 1)
            dos_file = workdir + '/DOS' + str(atom_indx) + '.dat'
        if save_dos_arr == True:
            np.savetxt(dos_file, dos_arr)
        else:
            pass

    doscar_dict = {}
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

def read_oszicar(OSZICAR_Dir):
    '''Read OSZICAR'''
    import numpy as np
    import matplotlib.pyplot as plt
    from .. import funcs
    with open(OSZICAR_Dir,'r') as f:
        line = f.readlines()

    ionic_step_num = 0
    for i_line in range(len(line)):
            if funcs.split_line(line[i_line])[1] == 'F=':
                ionic_step_num += 1
    etot_arr = np.array([None]*ionic_step_num)
    ionic_step = 0
    for i_line in range(len(line)):
        if split_line(line[i_line])[1] == 'F=':
            etot_arr[ionic_step] = float(funcs.split_line(line[i_line])[2])
            ionic_step += 1
    fig = plt.figure('NEtot')
    ax = fig.add_subplot(111)
    plt.plot(range(len(etot_arr)), etot_arr, marker = '.')
    ax.set(xlabel = 'Ionic steps')
    ax.set(ylabel = '$E_{tot}(eV)$')
    etot_oszicar_figfile = 'EtotOSZICAR.pdf'
    plt.savefig(etot_oszicar_figfile,dpi = 100)
    plt.close
    return etot_arr

