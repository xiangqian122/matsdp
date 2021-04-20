# -*- coding: utf-8 -*-
def nn_map(poscar_file_path, a0, n_shell = 1, elmt_selector_list = None, atom_label = None):
    '''
    Description:
        Nearest neighbor map. Now only availavle for face-centered cubic structures
    - Args
     * poscar_file_path: String format. It specifies the directory of the POSCAR file
     * a0: Float format. The lattice constant of the model. Unit in Angstrom
     * n_shell: Integer format. It determines up to which crystallographic shell the nearest neighbour map calculates
     * elmt_selector_list. Determines which type of elements are printed out in the nearest neighbor list. if the value is None, then all the elements are selected.
     * atom_label: defines which kind of atom label is used, choices are: 'elmtname', 'atomname', 'n_xyz', 'n_yzx', 'n_zxy', 'n_xzy', 'n_yxz', 'n_zyx' ,None or 'None'.
    '''
    import os
    import sys
    import math
    import numpy as np
    from sklearn.neighbors import KDTree
    from .. import funcs
    from . import vasp_read
    from . import vasp_write
    from .. import default_params

    # Scikit-learn requires Python (>= 3.5)
    if sys.version_info < (3,5):
        sys.exit('Sorry, Python < 3.5 is not supported. This module uses the package Scikit-learn, which requires Python (>= 3.5)')

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    poscar_file_path = os.path.abspath(poscar_file_path)
    n_neighbor = 100  #number of nearest neighbors in kNNList
    tol = 0.1 * a0  #tol is the maximum allowed deviation of atom posotions from perfect lattice site.

    if atom_label not in ['n_xyz', 'n_yzx', 'n_zxy', 'n_xzy', 'n_yxz', 'n_zyx', 'elmt_n_xyz', 'elmt_n_yzx', 'elmt_n_zxy', 'elmt_n_xzy', 'elmt_n_yxz', 'elmt_n_zyx']:
        number_order_text = 'xyz'
    elif atom_label in ['n_xyz', 'n_yzx', 'n_zxy', 'n_xzy', 'n_yxz', 'n_zyx', 'elmt_n_xyz', 'elmt_n_yzx', 'elmt_n_zxy', 'elmt_n_xzy', 'elmt_n_yxz', 'elmt_n_zyx']:
        number_order_text = atom_label[-3::]
    #Define lattice structure
    structure = "FCC"
    r = np.array([0] * n_shell, dtype = np.float)
    if structure == "FCC":
        first_nn = math.sqrt(2.0) / 2 * a0;
        for i_shell in range(0,n_shell):
            r[i_shell] = first_nn * math.sqrt(i_shell + 1.0);
            
    workdir, poscar_file = funcs.file_path_name(poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path, suppress_atom_number_ortho = False)
    n_atoms = np.sum(poscar_dict['elmt_num_arr'])  
    '''Find atom neighbors over periodic boundaries:
    Eextend/duplicate the cell in three dimensions, Create duplicates in three dimensions.
    Define the duplicity in each direction(include positive and negative directions) as 'Extend'.
    In each dimension, there are (2*n_extend+1) number of cells. (2*n_extend) of the (2*n_extend+1) cells are duplicates.
    There are (2*n_extend+1)^3-1 duplicates created. The volume is (2*n_extend+1)^3 times larger'''
    n_extend = 1 # Define the cell duplicity in each direction.
    # atom key list (including the shifted atoms, with shift label)
    atom_key_extend_list = list(poscar_dict['atom_key_arr'].copy())
    # atom key list (including the shifted atoms, but without shift label)
    atom_key_extend_list_no_shift_label = list(poscar_dict['atom_key_arr'].copy())
    # atom name list (including the shifted atoms, with shift label)
    atomname_extend_list = poscar_dict['atomname_list'].copy()
    # atom name list (including the shifted atoms, but without shift label)
    atomname_extend_list_no_shift_label = poscar_dict['atomname_list'].copy()
    # atom name list (including the shifted atoms, but without shift label) with atom indices
    atomname_extend_list_no_shift_label_with_atom_indx = poscar_dict['atomname_list'].copy()
    # atom species list (including the shifted atoms, but without shift label) with atom number ortho (atom numbering according to the order of three orthogonal directions)
    atomspecies_extend_list_no_shift_label_with_atom_number_ortho = list(poscar_dict['atom_species_arr'].copy())
    atom_number_ortho_extend_list = list(poscar_dict['atom_number_ortho_' + number_order_text].copy())
    car_extend_arr = np.array([0.0] * n_atoms * pow(2 * n_extend+1,3) * 3,dtype = np.float)
    car_extend_arr.shape = n_atoms*pow(2 * n_extend + 1, 3), 3
    formatted_atomname_len = 6 #For example len('Ni324') = 5, len('Ni1324') = 6
    formatted_len = 13 #For example len('Ni324_-1-1-1')=12, len('Ni1324_-1-1-1')=13
    formatted_len1 = 6
    for i_atom in range(n_atoms):
        car_extend_arr[i_atom,0] = poscar_dict['pos_arr'][i_atom,3]
        car_extend_arr[i_atom,1] = poscar_dict['pos_arr'][i_atom,4]
        car_extend_arr[i_atom,2] = poscar_dict['pos_arr'][i_atom,5]
        # change the format of the entries in the list
        # atom_key_extend_list[i_atom] = str(atom_key_extend_list[i_atom]) + ' ' * (formatted_len - len(str(atom_key_extend_list[i_atom]))) #atom_key_extend_list shoud contain intergers, not strings
        atom_key_extend_list_no_shift_label[i_atom] = str(atom_key_extend_list_no_shift_label[i_atom]) + ' ' * (formatted_len - len(str(atom_key_extend_list_no_shift_label[i_atom])))
        atomname_extend_list[i_atom] = atomname_extend_list[i_atom] + ' ' * (formatted_len - len(atomname_extend_list[i_atom]))
        atomname_extend_list_no_shift_label[i_atom] = atomname_extend_list_no_shift_label[i_atom] + ' ' * (formatted_len - len(atomname_extend_list_no_shift_label[i_atom]))
        atomname_extend_list_no_shift_label_with_atom_indx[i_atom] = atomname_extend_list_no_shift_label_with_atom_indx[i_atom] + ' (' + str(poscar_dict['atom_indx_arr'][i_atom]) + ')' + ' '*(formatted_len - len(atomname_extend_list_no_shift_label_with_atom_indx[i_atom] + ' (' + str(poscar_dict['atom_indx_arr'][i_atom]) + ')'))
        atomspecies_extend_list_no_shift_label_with_atom_number_ortho[i_atom] = atomspecies_extend_list_no_shift_label_with_atom_number_ortho[i_atom] + ' (' + str(poscar_dict['atom_number_ortho_' + number_order_text][i_atom]) + ')' + ' '*(formatted_len - len(atomspecies_extend_list_no_shift_label_with_atom_number_ortho[i_atom] + ' (' + str(poscar_dict['atom_number_ortho_' + number_order_text][i_atom]) + ')'))
        atom_number_ortho_extend_list[i_atom] = str(atom_number_ortho_extend_list[i_atom]) + ' ' * (formatted_len - len(str(atom_number_ortho_extend_list[i_atom])))
    i_atom_extend = i_atom
    unit_vec = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
    for i in range(-n_extend, n_extend + 1):
        for j in range(-n_extend, n_extend + 1):
            for k in range(-n_extend, n_extend + 1):
                if i == 0 and j == 0 and k == 0:
                    continue
                for i_atom in range(n_atoms):
                    i_atom_extend += 1
                    car_extend_arr[i_atom_extend,0] = poscar_dict['pos_arr'][i_atom,3] + np.dot(poscar_dict['l_arr'][:,0], unit_vec[0,:]) * i
                    car_extend_arr[i_atom_extend,1] = poscar_dict['pos_arr'][i_atom,4] + np.dot(poscar_dict['l_arr'][:,1], unit_vec[1,:]) * j
                    car_extend_arr[i_atom_extend,2] = poscar_dict['pos_arr'][i_atom,5] + np.dot(poscar_dict['l_arr'][:,2], unit_vec[2,:]) * k
                    #extend_atom_key = str(poscar_dict['atom_key_arr'][i_atom]) + '_' + str(i) + str(j) + str(k)
                    extend_atom_key = poscar_dict['atom_key_arr'][i_atom]
                    extend_atomname = poscar_dict['atomname_list'][i_atom] + '_' + str(i) + str(j) + str(k)
                    #atom_key_extend_list.append(extend_atom_key + ' ' * (formatted_len - len(extend_atom_key)))
                    atom_key_extend_list.append(extend_atom_key)
                    atomname_extend_list.append(extend_atomname + ' ' * (formatted_len - len(extend_atomname)))
                    atom_key_extend_list_no_shift_label.append(str(poscar_dict['atom_key_arr'][i_atom]) + ' '*(formatted_len - len(str(poscar_dict['atom_key_arr'][i_atom]))))
                    atomname_extend_list_no_shift_label.append(poscar_dict['atomname_list'][i_atom] + ' '*(formatted_len - len(poscar_dict['atomname_list'][i_atom])))
                    atomname_extend_list_no_shift_label_with_atom_indx.append(poscar_dict['atomname_list'][i_atom] + ' (' + str(poscar_dict['atom_indx_arr'][i_atom]) + ')' + ' '*(formatted_len - len(poscar_dict['atomname_list'][i_atom] + ' (' + str(poscar_dict['atom_indx_arr'][i_atom]) + ')')))
                    atomspecies_extend_list_no_shift_label_with_atom_number_ortho.append(poscar_dict['atom_species_arr'][i_atom] + ' (' + str(poscar_dict['atom_number_ortho_' + number_order_text][i_atom]) + ')' + ' '*(formatted_len - len(poscar_dict['atom_species_arr'][i_atom] + ' (' + str(poscar_dict['atom_number_ortho_' + number_order_text][i_atom]) + ')')))
                    atom_number_ortho_extend_list.append(str(poscar_dict['atom_number_ortho_' + number_order_text][i_atom]) + ' ' * (formatted_len - len(str(poscar_dict['atom_number_ortho_' + number_order_text][i_atom]))))

    #Calculate nearest neighbor information
    tree = KDTree(car_extend_arr)
    #Prepare the output files
    nn_atomname_list_file = os.path.join(workdir, 'nn_atomname_list.txt')
    nn_atomname_list_file_without_shift_label = os.path.join(workdir, 'nn_atomname_list.txt')
    nn_dist_file = os.path.join(workdir, 'nn_dist.txt')
    nn_map_temp_file = os.path.join(workdir, 'nn_map_temp.txt')
    nn_map_temp_file_without_shift_label = os.path.join(workdir, 'nn_map_temp_without_shift_label.txt')
    nn_map_file = os.path.join(workdir, 'nn_map.txt')
    nn_map_file_without_shift_label = os.path.join(workdir, 'nn_map_without_shift_label.txt')
    nn_count_1nn_poscar_file = os.path.join(workdir, 'poscar_nn_count_1nn.txt')

    nn_atomname_list_i_shell_str_list = [''] * n_shell
    nn_atomname_list_i_shell_without_shift_label_str_list = [''] * n_shell
    nn_atomname_list_i_shell_without_shift_label_with_atom_indx_str_list = [''] * n_shell
    nn_atomspecies_list_i_shell_without_shift_label_with_atom_number_ortho_str_list = [''] * n_shell
    nn_atom_number_ortho_list_i_shell_str_list = [''] * n_shell
    cna_atomname_list_i_shell_str_list = [''] * n_shell
    cna_atomname_list_i_shell_without_shift_label_str_list = [''] * n_shell
    cna_atomname_list_i_shell_without_shift_label_with_atom_indx_str_list = [''] * n_shell
    cna_atomspecies_list_i_shell_without_shift_label_with_atom_number_ortho_str_list = [''] * n_shell
    cna_atom_number_ortho_list_i_shell_str_list = [''] * n_shell
    nn_dist_list_ishell_str_list = [''] * n_shell
    nn_count_i_shell_str_list = [''] * n_shell

    nn_atomname_list_i_shell_file = []
    nn_atomname_list_i_shell_file_without_shift_label = []
    nn_atomname_list_i_shell_file_without_shift_label_with_atom_indx = []
    nn_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho = []
    nn_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho_sorted = []
    nn_atom_number_ortho_list_i_shell_file = []
    nn_atom_number_ortho_sorted_list_i_shell_file = []
    nn_dist_list_ishell_file = []
    nn_count_i_shell_file = []
    cna_atomname_list_i_shell_file = []
    cna_atomname_list_i_shell_file_without_shift_label = []
    cna_atomname_list_i_shell_file_without_shift_label_with_atom_indx = []
    cna_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho = []
    cna_atom_number_ortho_list_i_shell_file = []
    max_neighbor_num = 50
    #nn_atomname_ishell_arr constains the atom names with shift atom label
    #nn_atomname_ishell_arr_no_shift_label constains the atom names without shift atom label    
    nn_atom_key_ishell_arr = np.array([None] * n_atoms * n_shell * max_neighbor_num)
    nn_atom_key_ishell_arr.shape = n_atoms, n_shell, max_neighbor_num
    nn_atomname_ishell_arr = np.array([None] * n_atoms * n_shell * max_neighbor_num)
    nn_atomname_ishell_arr.shape = n_atoms, n_shell, max_neighbor_num
    nn_atomname_ishell_arr_no_shift_label = np.array([None] * n_atoms * n_shell * max_neighbor_num)
    nn_atomname_ishell_arr_no_shift_label.shape = n_atoms, n_shell, max_neighbor_num
    nn_atomname_ishell_arr_no_shift_label_with_atom_indx = np.array([None] * n_atoms * n_shell * max_neighbor_num)
    nn_atomname_ishell_arr_no_shift_label_with_atom_indx.shape = n_atoms, n_shell, max_neighbor_num
    nn_atomspecies_ishell_arr_no_shift_label_with_atom_number_ortho = np.array([None] * n_atoms * n_shell * max_neighbor_num)
    nn_atomspecies_ishell_arr_no_shift_label_with_atom_number_ortho.shape = n_atoms, n_shell, max_neighbor_num
    nn_atom_number_ortho_ishell_arr = np.array([None] * n_atoms * n_shell * max_neighbor_num)
    nn_atom_number_ortho_ishell_arr.shape = n_atoms, n_shell, max_neighbor_num
    nn_count_str = ''
    for i_shell in range(n_shell):
        shell_indx = i_shell +1
        nn_atomname_list_i_shell_file.append(os.path.join(workdir, 'nn_atomname_list_' + str(shell_indx) + 'NN.txt'))
        nn_atomname_list_i_shell_file_without_shift_label.append(os.path.join(workdir, 'nn_atomname_list_without_shift_label_' + str(shell_indx) + 'NN.txt'))
        nn_atomname_list_i_shell_file_without_shift_label_with_atom_indx.append(os.path.join(workdir, 'nn_atomname_list_without_shift_label_with_atomindx_' + str(shell_indx) + 'NN.txt'))
        nn_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho.append(os.path.join(workdir, 'nn_atomspecies_list_without_shift_label_with_atom_number_ortho_' + str(shell_indx) + 'NN.txt'))
        nn_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho_sorted.append(os.path.join(workdir, 'nn_atomspecies_list_without_shift_label_with_atom_number_ortho_sorted_' + str(shell_indx) + 'NN.txt'))
        nn_atom_number_ortho_list_i_shell_file.append(os.path.join(workdir, 'nn_atom_number_ortho_list_' + str(shell_indx) + 'NN.txt'))
        nn_atom_number_ortho_sorted_list_i_shell_file.append(os.path.join(workdir, 'nn_atom_number_ortho_sorted_list_' + str(shell_indx) + 'NN.txt'))
        nn_dist_list_ishell_file.append(os.path.join(workdir, 'nn_dist_list_' + str(shell_indx) + 'NN.txt'))
        nn_count_i_shell_file.append(os.path.join(workdir, 'nn_count_' + str(shell_indx) + 'NN.txt'))
        cna_atomname_list_i_shell_file.append(os.path.join(workdir, 'cna_atomname_list_' + str(shell_indx) + 'NN.txt'))
        cna_atomname_list_i_shell_file_without_shift_label.append(os.path.join(workdir, 'cna_atomname_list_without_shift_label_' + str(shell_indx) + 'NN.txt'))
        cna_atomname_list_i_shell_file_without_shift_label_with_atom_indx.append(os.path.join(workdir, 'cna_atomname_list_without_shift_label_with_atomindx_' + str(shell_indx) + 'NN.txt'))
        cna_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho.append(os.path.join(workdir, 'cna_atomspecies_list_without_shift_label_with_atom_number_ortho_' + str(shell_indx) + 'NN.txt'))
        cna_atom_number_ortho_list_i_shell_file.append(os.path.join(workdir, 'cna_atom_number_ortho_list_' + str(shell_indx) + 'NN.txt'))
        nn_count_str = nn_count_str + 'atom_name' + ' ' * (formatted_len - len('atom_name')) + ', ' + ', '.join([x + ' ' * (formatted_len1 - len(x)) for x in poscar_dict['elmt_species_arr']]) + ', ' + 'tot_count' + ' '*(formatted_len1 - len('tot_count')) + '\n'
        nn_count_i_shell_str_list[i_shell] = nn_count_i_shell_str_list[i_shell] + nn_count_str

    nn_map_str = ''
    nn_map_without_shift_label_str = ''
    nn_dist_str = ''
    nn_atomname_list_str = ''
    nn_atomname_list_without_shift_label_str = ''

    split_text = ' ' * formatted_atomname_len + ' '
    nn_map_str = nn_map_str + 'atom ' + ' '*(formatted_atomname_len - len('atom')) + split_text.join([str(i+1) + 'NN    ' for i in range(n_shell)]) + '\n'
    nn_map_without_shift_label_str = nn_map_without_shift_label_str + 'atom ' + ' '*(formatted_atomname_len - len('atom')) + split_text.join([str(i+1) + 'NN    ' for i in range(n_shell)]) + '\n'
    nn_map_dict = {}
    nn_map_dict['nn_count_1NN'] = []
    for i_atom in range(n_atoms):
        #dists, indices = tree.query([car_extend_arr[i_atom_extend]], k=[x for x in range(1,n_neighbor+1)])
        dists, indices = tree.query([car_extend_arr[i_atom]], k = n_neighbor + 1)
        nn_atomname_list = []
        nn_dist_list = []
        for i_neighbor in range(n_neighbor + 1):
            nn_atomname_list.append(atomname_extend_list[indices[0, i_neighbor]])
            nn_dist_list.append(str(dists[0, i_neighbor]))

        # select elements in the nearest neighbor list
        if elmt_selector_list in [None, 'None', 'none']:
            pass
        else:
            #screening subset of a the list according to specific keyword in the elmt_selector_list
            nn_atomname_extract_name_list = [''.join(filter(str.isalpha, x)) for x in nn_atomname_list]
            selector_indx_list = [i for (i, temp_elmt_name) in enumerate(nn_atomname_extract_name_list) if temp_elmt_name in elmt_selector_list]
                
            nn_atomname_list = [nn_atomname_list[x] for x in selector_indx_list]
            nn_dist_list = [nn_dist_list[x] for x in selector_indx_list]
                
        nn_atomname_list_str = nn_atomname_list_str + str(', '.join(nn_atomname_list)) + '\n' 
        nn_atomname_list_without_shift_label_str = nn_atomname_list_without_shift_label_str + str(', '.join(nn_atomname_list)) + '\n'
        nn_dist_str = nn_dist_str + str(atomname_extend_list[i_atom]) + ', ' + str(', '.join(['{:.4f}'.format(float(item)) for item in nn_dist_list])) + '\n'

        for i_shell in range(n_shell):
            nn_atom_key_ishell_list = []
            nn_atomname_ishell_list = []
            nn_atomname_ishell_list_no_shift_label = []
            nn_atomname_ishell_list_no_shift_label_with_atom_indx = []
            nn_atomspecies_ishell_list_no_shift_label_with_atom_number_ortho = []
            nn_atom_number_ortho_ishell_list = []
            nn_dist_ishell_list = []
            for i_neighbor in range(n_neighbor+1):
                if abs(dists[0,i_neighbor]-r[i_shell]) < tol:
                    nn_atom_key_ishell_list.append(atom_key_extend_list[indices[0,i_neighbor]])
                    nn_atomname_ishell_list.append(atomname_extend_list[indices[0,i_neighbor]])
                    nn_atomname_ishell_list_no_shift_label.append(atomname_extend_list_no_shift_label[indices[0,i_neighbor]])
                    nn_atomname_ishell_list_no_shift_label_with_atom_indx.append(atomname_extend_list_no_shift_label_with_atom_indx[indices[0,i_neighbor]])
                    nn_atomspecies_ishell_list_no_shift_label_with_atom_number_ortho.append(atomspecies_extend_list_no_shift_label_with_atom_number_ortho[indices[0,i_neighbor]])
                    nn_atom_number_ortho_ishell_list.append(atom_number_ortho_extend_list[indices[0,i_neighbor]])
                    nn_dist_ishell_list.append(str('{:.4f}'.format(dists[0,i_neighbor])))

            # select elements in the nearest neighbor list
            if elmt_selector_list in [None, 'None', 'none']:
                pass
            else:
                #screening subset of the list according to specific keyword in the elmt_selector_list
                nn_atomname_ishell_extract_name_list = [''.join(filter(str.isalpha, x)) for x in nn_atomname_ishell_list]
                selector_indx_list = [i for (i, temp_elmt_name) in enumerate(nn_atomname_ishell_extract_name_list) if temp_elmt_name in elmt_selector_list]
                nn_atomname_ishell_list = [nn_atomname_ishell_list[x] for x in selector_indx_list]
                nn_atom_key_ishell_list = [nn_atom_key_ishell_list[x] for x in selector_indx_list]
                nn_atomname_ishell_list_no_shift_label = [nn_atomname_ishell_list_no_shift_label[x] for x in selector_indx_list]
                nn_atomname_ishell_list_no_shift_label_with_atom_indx = [nn_atomname_ishell_list_no_shift_label_with_atom_indx[x] for x in selector_indx_list]
                nn_atomspecies_ishell_list_no_shift_label_with_atom_number_ortho = [nn_atomspecies_ishell_list_no_shift_label_with_atom_number_ortho[x] for x in selector_indx_list]
                nn_atom_number_ortho_ishell_list = [nn_atom_number_ortho_ishell_list[x] for x in selector_indx_list]
                nn_dist_ishell_list = [nn_dist_ishell_list[x] for x in selector_indx_list]

            nn_atom_key_ishell_arr[i_atom,i_shell,0:len(nn_atomname_ishell_list)] =  nn_atom_key_ishell_list
            nn_atomname_ishell_arr[i_atom,i_shell,0:len(nn_atomname_ishell_list)] =  nn_atomname_ishell_list
            nn_atomname_ishell_arr_no_shift_label[i_atom,i_shell,0:len(nn_atomname_ishell_list)] =  nn_atomname_ishell_list_no_shift_label
            nn_atomname_ishell_arr_no_shift_label_with_atom_indx[i_atom,i_shell,0:len(nn_atomname_ishell_list)] =  nn_atomname_ishell_list_no_shift_label_with_atom_indx
            nn_atomspecies_ishell_arr_no_shift_label_with_atom_number_ortho[i_atom,i_shell,0:len(nn_atomname_ishell_list)] =  nn_atomspecies_ishell_list_no_shift_label_with_atom_number_ortho
            nn_atom_number_ortho_ishell_arr[i_atom,i_shell,0:len(nn_atomname_ishell_list)] =  nn_atom_number_ortho_ishell_list
            # Strings for each shell are put in a list
            nn_atomname_list_i_shell_str_list[i_shell] =  nn_atomname_list_i_shell_str_list[i_shell] + str(atomname_extend_list[i_atom]) + ', ' + str(', '.join(nn_atomname_ishell_list))+'\n'
            nn_atomname_list_i_shell_without_shift_label_str_list[i_shell] = nn_atomname_list_i_shell_without_shift_label_str_list[i_shell] + str(atomname_extend_list[i_atom]) + ', ' + str(', '.join(nn_atomname_ishell_list_no_shift_label)) +'\n'
            nn_atomname_list_i_shell_without_shift_label_with_atom_indx_str_list[i_shell] = nn_atomname_list_i_shell_without_shift_label_with_atom_indx_str_list[i_shell] + str(atomname_extend_list_no_shift_label_with_atom_indx[i_atom]) + ', ' + str(', '.join(nn_atomname_ishell_list_no_shift_label_with_atom_indx)) +'\n'
            nn_atomspecies_list_i_shell_without_shift_label_with_atom_number_ortho_str_list[i_shell] = nn_atomspecies_list_i_shell_without_shift_label_with_atom_number_ortho_str_list[i_shell] + str(atomspecies_extend_list_no_shift_label_with_atom_number_ortho[i_atom]) + ', ' + str(', '.join(nn_atomspecies_ishell_list_no_shift_label_with_atom_number_ortho)) +'\n'
            nn_atom_number_ortho_list_i_shell_str_list[i_shell] = nn_atom_number_ortho_list_i_shell_str_list[i_shell] + str(atom_number_ortho_extend_list[i_atom]) + ', ' + str(', '.join(nn_atom_number_ortho_ishell_list)) +'\n'
            nn_dist_list_ishell_str_list[i_shell] = nn_dist_list_ishell_str_list[i_shell] + str(atomname_extend_list[i_atom]) + ', ' + str(', '.join(nn_dist_ishell_list))+'\n'
            nn_atom_count_list = [sum([y == z for y in [''.join(filter(str.isalpha, x)) for x in nn_atomname_ishell_list_no_shift_label]]) for z in poscar_dict['elmt_species_arr']]
            nn_tot_count = sum(nn_atom_count_list)
            nn_map_dict['nn_count_1NN'].append(nn_tot_count)
            nn_count_i_shell_str_list[i_shell] = nn_count_i_shell_str_list[i_shell] + str(atomname_extend_list[i_atom]) + ', ' + ', '.join(str(x) + ' ' * (formatted_len1 - len(str(x))) for x in nn_atom_count_list) + ', ' + str(nn_tot_count) + '\n'

    for i_atom in range(n_atoms):
        for i_shell in range(n_shell):
            for i_neighbor in range(max_neighbor_num):
                if nn_atomname_ishell_arr[i_atom,i_shell,i_neighbor] in [None,'None','none']:
                    nn_atomname_ishell_arr[i_atom,i_shell,i_neighbor] = ' '*formatted_len                  
                if nn_atomname_ishell_arr_no_shift_label[i_atom,i_shell,i_neighbor] in [None,'None','none']:
                    nn_atomname_ishell_arr_no_shift_label[i_atom,i_shell,i_neighbor] = ' '*formatted_len                  
                if nn_atomname_ishell_arr_no_shift_label_with_atom_indx[i_atom,i_shell,i_neighbor] in [None,'None','none']:
                    nn_atomname_ishell_arr_no_shift_label_with_atom_indx[i_atom,i_shell,i_neighbor] = ' '*formatted_len                  
                if nn_atomspecies_ishell_arr_no_shift_label_with_atom_number_ortho[i_atom,i_shell,i_neighbor] in [None,'None','none']:
                    nn_atomspecies_ishell_arr_no_shift_label_with_atom_number_ortho[i_atom,i_shell,i_neighbor] = ' '*formatted_len                  
                if nn_atom_number_ortho_ishell_arr[i_atom,i_shell,i_neighbor] in [None,'None','none']:
                    nn_atom_number_ortho_ishell_arr[i_atom,i_shell,i_neighbor] = ' '*formatted_len                  

    for i_atom in range(n_atoms):
        split_text = ' ' * formatted_atomname_len + ' '
        nn_text = split_text.join([' '.join(i) + '\n' for i in [nn_atomname_ishell_arr[i_atom,0:n_shell,i_neighbor].tolist() for i_neighbor in range(max_neighbor_num)]])
        nn_text_without_shift_label = split_text.join([' '.join(i) + '\n' for i in [nn_atomname_ishell_arr_no_shift_label[i_atom,0:n_shell,i_neighbor].tolist() for i_neighbor in range(max_neighbor_num)]])
        nn_map_temp_str = str(poscar_dict['atomname_list'][i_atom])+ ' ' * (formatted_atomname_len - len(poscar_dict['atomname_list'][i_atom])) + ' ' + nn_text + '\n' 
        nn_map_without_shift_label_temp_str = str(poscar_dict['atomname_list'][i_atom])+ ' ' * (formatted_atomname_len - len(poscar_dict['atomname_list'][i_atom])) + ' ' + nn_text_without_shift_label + '\n'
        nn_map_str = nn_map_str + str(nn_map_temp_str).replace('\t', '    ').strip('\n').strip() + '\n' 
        nn_map_without_shift_label_str = nn_map_without_shift_label_str + str(nn_map_without_shift_label_temp_str).replace('\t', '    ').strip('\n').strip() + '\n'

    for i_atom in range(n_atoms):
        for i_shell in range(n_shell):
            for i_neighbor in range(max_neighbor_num):
                if nn_atomname_ishell_arr[i_atom,i_shell,i_neighbor] == ' ' * formatted_len:
                    nn_atomname_ishell_arr[i_atom,i_shell,i_neighbor] = None              
                if nn_atomname_ishell_arr_no_shift_label[i_atom,i_shell,i_neighbor] == ' ' * formatted_len:
                    nn_atomname_ishell_arr_no_shift_label[i_atom,i_shell,i_neighbor] = None
                if nn_atomname_ishell_arr_no_shift_label_with_atom_indx[i_atom,i_shell,i_neighbor] == ' ' * formatted_len:
                    nn_atomname_ishell_arr_no_shift_label_with_atom_indx[i_atom,i_shell,i_neighbor] = None                  
                if nn_atomspecies_ishell_arr_no_shift_label_with_atom_number_ortho[i_atom,i_shell,i_neighbor] == ' ' * formatted_len:
                    nn_atomspecies_ishell_arr_no_shift_label_with_atom_number_ortho[i_atom,i_shell,i_neighbor] = None                  
                if nn_atom_number_ortho_ishell_arr[i_atom,i_shell,i_neighbor] == ' ' * formatted_len:
                    nn_atom_number_ortho_ishell_arr[i_atom,i_shell,i_neighbor] = None                  

    nn_atom_ishell_dict = {}
    nn_atom_ishell_dict['key'] = nn_atom_key_ishell_arr
    nn_atom_ishell_dict['atomname'] = nn_atomname_ishell_arr
    nn_atom_ishell_dict['atomname_no_shift_label'] = nn_atomname_ishell_arr_no_shift_label
    nn_atom_ishell_dict['atomname_no_shift_label_with_atom_indx'] = nn_atomname_ishell_arr_no_shift_label_with_atom_indx
    nn_atom_ishell_dict['atomspecies_no_shift_label_with_atom_number_ortho'] = nn_atomspecies_ishell_arr_no_shift_label_with_atom_number_ortho
    nn_atom_ishell_dict['atom_number_ortho' + number_order_text] = nn_atom_number_ortho_ishell_arr

    # common neighbor analysis
    for i_shell in range(n_shell):
        for i_atom in range(n_atoms):
            for j_atom in range(n_atoms):
                if i_atom == j_atom:
                    continue
                i_atom_arr = list(filter(None, nn_atom_ishell_dict['atomname_no_shift_label'][i_atom, i_shell, :]))
                j_atom_arr = list(filter(None, nn_atom_ishell_dict['atomname_no_shift_label'][j_atom, i_shell, :]))
                ij_intersection_arr = np.array([x for x in i_atom_arr if x in j_atom_arr])
                if len(ij_intersection_arr) == 0:
                    continue
                # atomname without shift label style
                i_atomname = atomname_extend_list[i_atom]
                j_atomname = atomname_extend_list[j_atom]
                cna_str = str(i_atomname) + ' - ' + str(j_atomname) + ',' + ','.join(ij_intersection_arr)
                cna_atomname_list_i_shell_without_shift_label_str_list[i_shell] = cna_atomname_list_i_shell_without_shift_label_str_list[i_shell] + cna_str + '\n'
                # atomname without shift label with atom index style
                i_atom_arr = list(filter(None, nn_atom_ishell_dict['atomname_no_shift_label_with_atom_indx'][i_atom, i_shell, :]))
                j_atom_arr = list(filter(None, nn_atom_ishell_dict['atomname_no_shift_label_with_atom_indx'][j_atom, i_shell, :]))
                ij_intersection_arr = np.array([x for x in i_atom_arr if x in j_atom_arr])
                i_atomname = atomname_extend_list_no_shift_label_with_atom_indx[i_atom]
                j_atomname = atomname_extend_list_no_shift_label_with_atom_indx[j_atom]
                cna_str = str(i_atomname) + ' - ' + str(j_atomname) + ',' + ','.join(ij_intersection_arr)
                cna_atomname_list_i_shell_without_shift_label_with_atom_indx_str_list[i_shell] = cna_atomname_list_i_shell_without_shift_label_with_atom_indx_str_list[i_shell] + cna_str + '\n'
                # with atom species and atom number ortho style
                i_atom_arr = list(filter(None, nn_atom_ishell_dict['atomspecies_no_shift_label_with_atom_number_ortho'][i_atom, i_shell, :]))
                j_atom_arr = list(filter(None, nn_atom_ishell_dict['atomspecies_no_shift_label_with_atom_number_ortho'][j_atom, i_shell, :]))
                ij_intersection_arr = np.array([x for x in i_atom_arr if x in j_atom_arr])
                i_atom_species_with_atom_number_ortho = atomspecies_extend_list_no_shift_label_with_atom_number_ortho[i_atom]
                j_atom_species_with_atom_number_ortho = atomspecies_extend_list_no_shift_label_with_atom_number_ortho[j_atom]
                cna_str = str(i_atom_species_with_atom_number_ortho) + ' - ' + str(j_atom_species_with_atom_number_ortho) + ',' + ','.join(ij_intersection_arr)
                cna_atomspecies_list_i_shell_without_shift_label_with_atom_number_ortho_str_list[i_shell] = cna_atomspecies_list_i_shell_without_shift_label_with_atom_number_ortho_str_list[i_shell] + cna_str + '\n'
                # with atom number ortho style
                i_atom_arr = list(filter(None, nn_atom_ishell_dict['atom_number_ortho' + number_order_text][i_atom, i_shell, :]))
                j_atom_arr = list(filter(None, nn_atom_ishell_dict['atom_number_ortho' + number_order_text][j_atom, i_shell, :]))
                ij_intersection_arr = np.array([x for x in i_atom_arr if x in j_atom_arr])
                i_atom_number_ortho = poscar_dict['atom_number_ortho_' + number_order_text][i_atom]
                j_atom_number_ortho = poscar_dict['atom_number_ortho_' + number_order_text][j_atom]
                cna_str = str(i_atom_number_ortho) + ' - ' + str(j_atom_number_ortho) + ',' + ','.join(ij_intersection_arr)
                cna_atom_number_ortho_list_i_shell_str_list[i_shell] = cna_atom_number_ortho_list_i_shell_str_list[i_shell] + cna_str + '\n'
                # atomname with shift label style
                i_atom_arr = list(filter(None, nn_atom_ishell_dict['atomname'][i_atom, i_shell, :]))
                j_atom_arr = list(filter(None, nn_atom_ishell_dict['atomname'][j_atom, i_shell, :]))
                ij_intersection_arr = np.array([x for x in i_atom_arr if x in j_atom_arr])
                if len(ij_intersection_arr) == 0 :
                    continue
                i_atomname = atomname_extend_list[i_atom]
                j_atomname = atomname_extend_list[j_atom]
                cna_str = str(i_atomname) + ' - ' + str(j_atomname) + ',' + ','.join(ij_intersection_arr)
                cna_atomname_list_i_shell_str_list[i_shell] = cna_atomname_list_i_shell_str_list[i_shell] + cna_str + '\n'
        with open(cna_atomname_list_i_shell_file_without_shift_label[i_shell],'w') as f_cna_atomname_without_shift_label, open(cna_atomname_list_i_shell_file_without_shift_label_with_atom_indx[i_shell], 'w') as f_cna_atomname_without_shift_label_with_atom_indx, open(cna_atomname_list_i_shell_file[i_shell], 'w') as f_cna_atomname_with_shift_label, open(cna_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho[i_shell], 'w') as f_cna_atom_species_atom_number_ortho, open(cna_atom_number_ortho_list_i_shell_file[i_shell], 'w') as f_cna_atom_number_ortho:
            f_cna_atomname_without_shift_label.write(cna_atomname_list_i_shell_without_shift_label_str_list[i_shell])
            f_cna_atomname_without_shift_label_with_atom_indx.write(cna_atomname_list_i_shell_without_shift_label_with_atom_indx_str_list[i_shell])
            f_cna_atomname_with_shift_label.write(cna_atomname_list_i_shell_str_list[i_shell])
            f_cna_atom_species_atom_number_ortho.write(cna_atomspecies_list_i_shell_without_shift_label_with_atom_number_ortho_str_list[i_shell])
            f_cna_atom_number_ortho.write(cna_atom_number_ortho_list_i_shell_str_list[i_shell])

    #write results
    with open(nn_atomname_list_file,'w') as f1, open(nn_dist_file,'w') as f2, open(nn_atomname_list_file_without_shift_label,'w') as f11:
        f1.write(nn_atomname_list_str)
        f11.write(nn_atomname_list_without_shift_label_str)
        f2.write(nn_dist_str)
    for i_shell in range(n_shell):
        with open(nn_atomname_list_i_shell_file[i_shell],'w') as f3, open(nn_dist_list_ishell_file[i_shell],'w') as f4, open(nn_atomname_list_i_shell_file_without_shift_label[i_shell],'w') as f13, open(nn_count_i_shell_file[i_shell], 'w') as f_nn_count, open(nn_atomname_list_i_shell_file_without_shift_label_with_atom_indx[i_shell],'w') as f_with_atom_indx, open(nn_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho[i_shell],'w') as f_with_atom_species_atom_number_ortho, open(nn_atom_number_ortho_list_i_shell_file[i_shell],'w') as f_with_atom_number_ortho:
            f3.write(nn_atomname_list_i_shell_str_list[i_shell])
            f13.write(nn_atomname_list_i_shell_without_shift_label_str_list[i_shell])
            f_with_atom_indx.write(nn_atomname_list_i_shell_without_shift_label_with_atom_indx_str_list[i_shell])
            f_with_atom_species_atom_number_ortho.write(nn_atomspecies_list_i_shell_without_shift_label_with_atom_number_ortho_str_list[i_shell])
            f_with_atom_number_ortho.write(nn_atom_number_ortho_list_i_shell_str_list[i_shell])
            f4.write(nn_dist_list_ishell_str_list[i_shell])
            f_nn_count.write(nn_count_i_shell_str_list[i_shell])

    with open(nn_map_file,'w') as f6:
        f6.write(nn_map_str)
    with open(nn_map_file_without_shift_label,'w') as f16:
        f16.write(nn_map_without_shift_label_str)
    vasp_write.write_poscar(output_poscar_file_path = nn_count_1nn_poscar_file,
                                               poscar_dict = poscar_dict,
                                               added_atom_data = nn_map_dict['nn_count_1NN'],
                                               added_atom_property_str = 'nn_count_1nn',
                                               added_atom_property_columns_str = 'nn_count_1nn')
    

    for i_shell in range(n_shell):
        with open(nn_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho[i_shell], 'r') as f1, open(nn_atom_number_ortho_list_i_shell_file[i_shell], 'r') as f2:
            f1_lines = f1.readlines()
            f2_lines = f2.readlines()
            line_num = len(f1_lines)
            # the first column stores the atom_number_ortho, the second column contains the NN atoms information for each atom
            atom_number_ortho_nn_f1_text_arr = np.array([None] * line_num * 2)
            atom_number_ortho_nn_f1_text_arr.shape = line_num, 2
            atom_number_ortho_nn_f2_text_arr = np.array([None] * line_num * 2)
            atom_number_ortho_nn_f2_text_arr.shape = line_num, 2
            for i in range(0, line_num):
                temp_list_1 = funcs.split_line(f1_lines[i], separator = ',')
                temp_list_2 = funcs.split_line(f2_lines[i], separator = ',')
                atom_number_ortho_nn_f1_text_arr[i, 0] = int(str(funcs.extract_num(temp_list_1[0])[0]))
                atom_number_ortho_nn_f1_text_arr[i, 1] = f1_lines[i]
                atom_number_ortho_nn_f2_text_arr[i, 0] = int(str(funcs.extract_num(temp_list_2[0])[0]))
                atom_number_ortho_nn_f2_text_arr[i, 1] = f2_lines[i]
        f1_argsort_indx_arr = np.argsort(atom_number_ortho_nn_f1_text_arr[:,0])
        f2_argsort_indx_arr = np.argsort(atom_number_ortho_nn_f2_text_arr[:,0])
        temp_text_1 = ''
        temp_text_2 = ''
        for i_atom in range(n_atoms):
            temp_text_1 = temp_text_1 + atom_number_ortho_nn_f1_text_arr[f1_argsort_indx_arr][i_atom,1]
            temp_text_2 = temp_text_2 + atom_number_ortho_nn_f2_text_arr[f2_argsort_indx_arr][i_atom,1]
        with open(nn_atomspecies_list_i_shell_file_without_shift_label_with_atom_number_ortho_sorted[i_shell], 'w') as f1, open(nn_atom_number_ortho_sorted_list_i_shell_file[i_shell], 'w') as f2:
            f1.write(temp_text_1)
            f2.write(temp_text_2)

#    print(atom_number_ortho_nn_text_arr[argsort_indx_arr])

    funcs.write_log(
        logfile,
        'vasp_analyze.nn_map(' + '\n' +
        '    poscar_file_path = ' + 'r\'' + str(poscar_file_path) + '\'' + ',\n' +
        '    a0 = ' + str(a0) + ',\n' +
        '    n_shell = ' + str(n_shell) + ',\n' +
        '    elmt_selector_list = ' + str(elmt_selector_list) + ',\n' +
        '    atom_label = ' + str(atom_label) + ')\n' +
        '###############################\n')
    return nn_atom_ishell_dict

def simple_cna(poscar_file_path, a0, common_neighbor_elmt_list = [], atom_label = None):
    '''
    simple common neighbor analysis (CNA)
    '''
    import os
    import re
    import itertools
    import numpy as np
    from .. import funcs
    from .. import default_params
    from . import vasp_read
    from . import vasp_write

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    ##For test purpose, the next two lines are commented.
    ##common_neighbor_elmt_list = ['Re','W','Ta','Ru']
    ##first_nn_list = ['Ni1', 'Ni1', 'Ni2', 'Ni38', 'Re1', 'Re1', 'Re5', 'W2', 'W4', 'Ta5',None,None]
    poscar_file_path = os.path.abspath(poscar_file_path)
    workdir, poscar_file = funcs.file_path_name(poscar_file_path)
    nn_atom_ishell_dict = nn_map(poscar_file_path, a0, n_shell = 1, elmt_selector_list = None, atom_label = None)

    poscar_dict = vasp_read.read_poscar(poscar_file_path, suppress_atom_number_ortho = False)
    n_atoms = np.sum(poscar_dict['elmt_num_arr'])

    if atom_label not in ['n_xyz', 'n_yzx', 'n_zxy', 'n_xzy', 'n_yxz', 'n_zyx']:
        number_order_text = 'xyz'
    elif atom_label in ['n_xyz', 'n_yzx', 'n_zxy', 'n_xzy', 'n_yxz', 'n_zyx']:
        number_order_text = atom_label[-3::]
    # define the element pairs which need to be anlyzed
    elmt_pair_list = [item for item in common_neighbor_elmt_list if item in poscar_dict['elmt_species_arr']]
    num_first_nn_comb = len(list(itertools.combinations([0] * 12, 2)))
    num_elmt_comb = len(list(itertools.combinations_with_replacement(elmt_pair_list, 2)))
    common_neighbor_pair_list = [None] * num_elmt_comb
    common_neighbor_pair_count_arr = np.array([None] * n_atoms*num_elmt_comb)
    common_neighbor_pair_count_arr.shape = n_atoms, num_elmt_comb
    common_neighbor_pair_content_arr = np.array([None] * n_atoms*num_elmt_comb*num_first_nn_comb)
    common_neighbor_pair_content_arr.shape = n_atoms,num_elmt_comb, num_first_nn_comb
    common_neighbor_pair_content_arr_with_atom_indx = np.array([None] * n_atoms * num_elmt_comb * num_first_nn_comb)
    common_neighbor_pair_content_arr_with_atom_indx.shape = n_atoms, num_elmt_comb, num_first_nn_comb
    common_neighbor_pair_content_arr_with_atomspecies_atom_number_ortho = np.array([None] * n_atoms * num_elmt_comb * num_first_nn_comb)
    common_neighbor_pair_content_arr_with_atomspecies_atom_number_ortho.shape = n_atoms, num_elmt_comb, num_first_nn_comb

    indx = 0
    for i_pair in itertools.combinations_with_replacement(elmt_pair_list, 2):
        common_neighbor_pair_list[indx] = i_pair[0] + '-' + i_pair[1]
        common_neighbor_pair_count_arr[:, indx] = 0
        indx += 1

    formatted_atomname_len = 6 #For example len('Ni324') = 5, len('Ni1324') = 6
    formatted_atomname_len1 = 13 #For example len('Ni324(324)') = 10, len('Ni1324(1324)') = 12
    formatted_len1 = 5 #For example len('Ni-Re')=5
    formatted_len2 = 13 #For example len('Ni329-Re123')=11, len('Ni1329-Re1123')=13
    formatted_len3 = 12 #For example len('-54.0812')=8, len('-54.08128345')=12    
    formatted_len4 = 24 #For example len('Ni329(329)-Re123(421)')=21, len('Ni1329(329)-Re1123(1234)')=24

    for i_atom in range(0, n_atoms):
        i_atom_name = poscar_dict['atomname_list'][i_atom]

        first_nn_list = list(filter(None, nn_atom_ishell_dict['atomname_no_shift_label'][i_atom, 0, :]))
        # for the first_nn_key_list, we can't use filter function, because the filter function will remove the value 0
        first_nn_key_list = [x for x in nn_atom_ishell_dict['key'][i_atom,0,:] if x is not None]

        first_nn_list_unique = list(set(first_nn_list))
        num_first_nn = len(first_nn_list)
        ##print(first_nn_list)

        # element type list
        first_nn_elmt_list = []
        for i in range(0,len(first_nn_list)):
            first_nn_elmt_list.append(''.join(re.findall(r'[A-Za-z]', first_nn_list[i])))
        # unique element type list
        first_nn_elmt_list_unique = list(set(first_nn_elmt_list))

        if len(list(set(first_nn_elmt_list_unique).intersection(set(elmt_pair_list)))) == 0:
            continue
        first_nn_list_temp = first_nn_list.copy()
        for i in range(0,len(first_nn_list_temp)):
            if ''.join(re.findall(r'[A-Za-z]', first_nn_list_temp[i])) not in elmt_pair_list:
                first_nn_list.remove(first_nn_list_temp[i])

##        for i_pair in itertools.combinations(first_nn_list, 2):
        for item in itertools.combinations(first_nn_key_list, 2):
            i_pair = [poscar_dict['atomname_list'][item[0]], poscar_dict['atomname_list'][item[1]]]
            if ''.join(re.findall(r'[A-Za-z]', i_pair[0])) not in common_neighbor_elmt_list or ''.join(re.findall(r'[A-Za-z]', i_pair[1])) not in common_neighbor_elmt_list:
                continue
            pair_str1 = ''.join(re.findall(r'[A-Za-z]', i_pair[0])) + '-' + ''.join(re.findall(r'[A-Za-z]', i_pair[1]))
            pair_str2 = ''.join(re.findall(r'[A-Za-z]', i_pair[1])) + '-' + ''.join(re.findall(r'[A-Za-z]', i_pair[0]))  
            if pair_str1 in common_neighbor_pair_list:
                indx = common_neighbor_pair_list.index(pair_str1)
                pair_str = i_pair[0].strip() + '-' + i_pair[1].strip()
                pair_str_with_atom_indx = i_pair[0].strip() + '(' + str(poscar_dict['atom_indx_arr'][item[0]]) + ')' + '-' + i_pair[1].strip() + '(' + str(poscar_dict['atom_indx_arr'][item[1]]) + ')'
                pair_str_with_atomspecies_atom_number_ortho = funcs.extract_alpha_str(i_pair[0].strip())[0] + '(' + str(poscar_dict['atom_number_ortho_' + number_order_text][item[0]]) + ')' + '-' + funcs.extract_alpha_str(i_pair[1].strip())[0] + '(' + str(poscar_dict['atom_number_ortho_' + number_order_text][item[1]]) + ')'
                pair_full_str = pair_str + (' ' * (formatted_len2-len(pair_str)))
                pair_full_str_with_atom_indx = pair_str_with_atom_indx + (' ' * (formatted_len4-len(pair_str_with_atom_indx)))
                pair_full_str_with_atomspecies_atom_number_ortho = pair_str_with_atomspecies_atom_number_ortho + (' ' * (formatted_len4-len(pair_str_with_atomspecies_atom_number_ortho)))
            elif pair_str2 in common_neighbor_pair_list:
                indx = common_neighbor_pair_list.index(pair_str2)
                pair_str = i_pair[1].strip() + '-' + i_pair[0].strip()
                pair_str_with_atom_indx = i_pair[1].strip() + '(' + str(poscar_dict['atom_indx_arr'][item[1]]) + ')' + '-' + i_pair[0].strip() + '(' + str(poscar_dict['atom_indx_arr'][item[0]]) + ')'
                pair_str_with_atomspecies_atom_number_ortho = funcs.extract_alpha_str(i_pair[1].strip())[0] + '(' + str(poscar_dict['atom_number_ortho_' + number_order_text][item[1]]) + ')' + '-' + funcs.extract_alpha_str(i_pair[0].strip())[0] + '(' + str(poscar_dict['atom_number_ortho_' + number_order_text][item[0]]) + ')'
                pair_full_str = pair_str + (' ' * (formatted_len2-len(pair_str)))
                pair_full_str_with_atom_indx = pair_str_with_atom_indx + (' ' * (formatted_len4-len(pair_str_with_atom_indx)))
                pair_full_str_with_atomspecies_atom_number_ortho = pair_str_with_atomspecies_atom_number_ortho + (' ' * (formatted_len4-len(pair_str_with_atomspecies_atom_number_ortho)))
            common_neighbor_pair_content_arr[i_atom, indx, common_neighbor_pair_count_arr[i_atom, indx]] =  pair_full_str
            common_neighbor_pair_content_arr_with_atom_indx[i_atom, indx, common_neighbor_pair_count_arr[i_atom, indx]] =  pair_full_str_with_atom_indx
            common_neighbor_pair_content_arr_with_atomspecies_atom_number_ortho[i_atom, indx, common_neighbor_pair_count_arr[i_atom, indx]] =  pair_full_str_with_atomspecies_atom_number_ortho
            common_neighbor_pair_count_arr[i_atom, indx] += 1

        common_neighbor_dict={}
        common_neighbor_dict['pair_name'] = common_neighbor_pair_list
        common_neighbor_dict['pair_count'] = common_neighbor_pair_count_arr
        common_neighbor_dict['pair_content'] = common_neighbor_pair_content_arr
        common_neighbor_dict['pair_content_with_atom_indx'] = common_neighbor_pair_content_arr_with_atom_indx
        common_neighbor_dict['pair_content_with_atomspecies_atom_number_ortho'] = common_neighbor_pair_content_arr_with_atomspecies_atom_number_ortho

    elmt_pair_list_str = ''.join(elmt_pair_list)
    simple_common_neighbor_pair_count_file = os.path.join(workdir, 'simple_common_neighbor_pair_count_' + elmt_pair_list_str + '.txt')
    simple_common_neighbor_pair_count_file_with_atom_indx = os.path.join(workdir, 'simple_common_neighbor_pair_count_with_atom_indx_' + elmt_pair_list_str + '.txt')
    simple_common_neighbor_pair_count_file_with_atomspecies_atom_number_ortho = os.path.join(workdir, 'simple_common_neighbor_pair_count_with_atomspecies_atom_number_ortho_' + elmt_pair_list_str + '.txt')
    simple_common_neighbor_file = os.path.join(workdir, 'simple_common_neighbor_' + elmt_pair_list_str + '.txt')
    simple_common_neighbor_file_with_atom_indx = os.path.join(workdir, 'simple_common_neighbor_with_atom_indx_' + elmt_pair_list_str + '.txt')
    simple_common_neighbor_file_with_atomspecies_atom_number_ortho = os.path.join(workdir, 'simple_common_neighbor_with_atomspecies_atom_number_ortho_' + elmt_pair_list_str + '.txt')

    for i_atom in range(n_atoms):
        for i_elmt_comb in range(num_elmt_comb):
            for i_pair_indx in range(num_first_nn_comb):
                if common_neighbor_dict['pair_content'][i_atom, i_elmt_comb, i_pair_indx] == None:
                    common_neighbor_dict['pair_content'][i_atom, i_elmt_comb, i_pair_indx] = ' '*formatted_len2
                if common_neighbor_dict['pair_content_with_atom_indx'][i_atom, i_elmt_comb, i_pair_indx] == None:
                    common_neighbor_dict['pair_content_with_atom_indx'][i_atom, i_elmt_comb, i_pair_indx] = ' '*formatted_len4
                if common_neighbor_dict['pair_content_with_atomspecies_atom_number_ortho'][i_atom, i_elmt_comb, i_pair_indx] == None:
                    common_neighbor_dict['pair_content_with_atomspecies_atom_number_ortho'][i_atom, i_elmt_comb, i_pair_indx] = ' '*formatted_len4
    
    simple_common_neighbor_pair_count_str = ''
    simple_common_neighbor_pair_count_str_with_atom_indx = ''
    simple_common_neighbor_pair_count_str_with_atomspecies_atom_number_ortho = ''
    simple_common_neighbor_str = ''
    simple_common_neighbor_str_with_atom_indx = ''
    simple_common_neighbor_str_with_atomspecies_atom_number_ortho = ''
##    simple_common_neighbor_pair_count_str = simple_common_neighbor_pair_count_str + 'atom   ' + str(' '.join(common_neighbor_dict['pair_name']))+'\n'
    simple_common_neighbor_pair_count_str = simple_common_neighbor_pair_count_str + 'atom   ' + ' '.join([x + ' '*(formatted_len1 - len(x)) for x in common_neighbor_dict['pair_name']]) + '\n'
    simple_common_neighbor_pair_count_str_with_atom_indx = simple_common_neighbor_pair_count_str_with_atom_indx + 'atom          ' + ' '.join([x + ' '*(formatted_len1 - len(x)) for x in common_neighbor_dict['pair_name']]) + '\n'
    simple_common_neighbor_pair_count_str_with_atomspecies_atom_number_ortho = simple_common_neighbor_pair_count_str_with_atomspecies_atom_number_ortho + 'atom          ' + ' '.join([x + ' '*(formatted_len1 - len(x)) for x in common_neighbor_dict['pair_name']]) + '\n'
##    simple_common_neighbor_pair_count_str_with_atom_indx = simple_common_neighbor_pair_count_str_with_atom_indx + 'atom          ' + str(' '.join(common_neighbor_dict['pair_name']))+'\n'
    simple_common_neighbor_str = simple_common_neighbor_str + 'atom   ' + str(' '.join([item + (' '*(formatted_len2 - len(str(item)))) for item in common_neighbor_dict['pair_name']]))+'\n'
    simple_common_neighbor_str_with_atom_indx = simple_common_neighbor_str_with_atom_indx + 'atom          ' + str(' '.join([item + (' '*(formatted_len4 - len(str(item)))) for item in common_neighbor_dict['pair_name']]))+'\n'
    simple_common_neighbor_str_with_atomspecies_atom_number_ortho = simple_common_neighbor_str_with_atomspecies_atom_number_ortho + 'atom          ' + str(' '.join([item + (' '*(formatted_len4 - len(str(item)))) for item in common_neighbor_dict['pair_name']]))+'\n'
    for i_atom in range(n_atoms):
        if len(list(set(first_nn_elmt_list_unique).intersection(set(elmt_pair_list)))) == 0:
            continue
        i_atom_name = poscar_dict['atomname_list'][i_atom]
        i_atom_name_with_atom_indx = poscar_dict['atomname_list'][i_atom] + '(' + str(poscar_dict['atom_indx_arr'][i_atom]) + ')' 
        i_atom_name_with_atomspecies_atom_number_ortho = poscar_dict['atom_species_arr'][i_atom] + '(' + str(poscar_dict['atom_number_ortho_' + number_order_text][i_atom]) + ')' 
        simple_common_neighbor_pair_count_str = simple_common_neighbor_pair_count_str + i_atom_name + ' '*(formatted_atomname_len - len(i_atom_name)) + ' ' + ' '.join([str(item) + (' '*(formatted_len1 - len(str(item)))) for item in common_neighbor_dict['pair_count'][i_atom,:]])+'\n'
        simple_common_neighbor_pair_count_str_with_atom_indx = simple_common_neighbor_pair_count_str_with_atom_indx + i_atom_name_with_atom_indx + ' '*(formatted_atomname_len1 - len(i_atom_name_with_atom_indx)) + ' ' + ' '.join([str(item) + (' '*(formatted_len1 - len(str(item)))) for item in common_neighbor_dict['pair_count'][i_atom,:]])+'\n'
        simple_common_neighbor_pair_count_str_with_atomspecies_atom_number_ortho = simple_common_neighbor_pair_count_str_with_atomspecies_atom_number_ortho + i_atom_name_with_atomspecies_atom_number_ortho + ' '*(formatted_atomname_len1 - len(i_atom_name_with_atomspecies_atom_number_ortho)) + ' ' + ' '.join([str(item) + (' '*(formatted_len1 - len(str(item)))) for item in common_neighbor_dict['pair_count'][i_atom,:]])+'\n'

        simple_common_neighbor_str = simple_common_neighbor_str + i_atom_name + ' '*(formatted_atomname_len - len(i_atom_name)) + ' ' + ' '.join([str(item) + (' '*(formatted_len2 - len(str(item)))) for item in common_neighbor_dict['pair_count'][i_atom,:]])+'\n'
        simple_common_neighbor_str_with_atom_indx = simple_common_neighbor_str_with_atom_indx + i_atom_name_with_atom_indx + ' '*(formatted_atomname_len1 - len(i_atom_name_with_atom_indx)) + ' ' + ' '.join([str(item) + (' '*(formatted_len4 - len(str(item)))) for item in common_neighbor_dict['pair_count'][i_atom,:]])+'\n'
        simple_common_neighbor_str_with_atomspecies_atom_number_ortho = simple_common_neighbor_str_with_atomspecies_atom_number_ortho + i_atom_name_with_atomspecies_atom_number_ortho + ' '*(formatted_atomname_len1 - len(i_atom_name_with_atomspecies_atom_number_ortho)) + ' ' + ' '.join([str(item) + (' '*(formatted_len4 - len(str(item)))) for item in common_neighbor_dict['pair_count'][i_atom,:]])+'\n'
        split_text = ' '*formatted_atomname_len + ' '
        split_text_with_atom_indx = ' '*formatted_atomname_len1 + ' '
        split_text_with_atomspecies_atom_number_ortho = ' '*formatted_atomname_len1 + ' '
        nn_text = split_text.join([' '.join(i) + '\n' for i in [common_neighbor_dict['pair_content'][i_atom, 0:num_elmt_comb, i_pair_indx].tolist() for i_pair_indx in range(num_first_nn_comb)]])
        nn_text_with_atom_indx = split_text_with_atom_indx.join([' '.join(i) + '\n' for i in [common_neighbor_dict['pair_content_with_atom_indx'][i_atom, 0:num_elmt_comb, i_pair_indx].tolist() for i_pair_indx in range(num_first_nn_comb)]])
        nn_text_with_atomspecies_atom_number_ortho = split_text_with_atomspecies_atom_number_ortho.join([' '.join(i) + '\n' for i in [common_neighbor_dict['pair_content_with_atomspecies_atom_number_ortho'][i_atom, 0:num_elmt_comb, i_pair_indx].tolist() for i_pair_indx in range(num_first_nn_comb)]])
        simple_common_neighbor_str = simple_common_neighbor_str + (' '*formatted_atomname_len) + ' ' + nn_text + '\n'
        simple_common_neighbor_str_with_atom_indx = simple_common_neighbor_str_with_atom_indx + (' '*formatted_atomname_len1) + ' ' + nn_text_with_atom_indx + '\n'
        simple_common_neighbor_str_with_atomspecies_atom_number_ortho = simple_common_neighbor_str_with_atomspecies_atom_number_ortho + (' '*formatted_atomname_len1) + ' ' + nn_text_with_atomspecies_atom_number_ortho + '\n'

    # write file
    with open(simple_common_neighbor_pair_count_file,'w') as f1, open(simple_common_neighbor_file,'w') as f2, open(simple_common_neighbor_file_with_atom_indx,'w') as f3, open(simple_common_neighbor_pair_count_file_with_atom_indx,'w') as f4, open(simple_common_neighbor_file_with_atomspecies_atom_number_ortho,'w') as f5, open(simple_common_neighbor_pair_count_file_with_atomspecies_atom_number_ortho,'w') as f6:
        f1.write(simple_common_neighbor_pair_count_str)
        f2.write(simple_common_neighbor_str)
        f3.write(simple_common_neighbor_str_with_atom_indx)
        f4.write(simple_common_neighbor_pair_count_str_with_atom_indx)
        f5.write(simple_common_neighbor_str_with_atomspecies_atom_number_ortho)
        f6.write(simple_common_neighbor_pair_count_str_with_atomspecies_atom_number_ortho)
        
    # write the POSCAR file with atom property data
    simple_cna_pair_count_poscar_file_path = os.path.join(workdir, 'POSCAR_simple_common_neighbor_pair_count_' + elmt_pair_list_str + '.vasp')
    added_atom_property_str = 'simple_cna'
    added_atom_property_columns_str = str(' '.join([item for item in common_neighbor_dict['pair_name']]))
    vasp_write.write_poscar(output_poscar_file_path = simple_cna_pair_count_poscar_file_path,
                                               poscar_dict = poscar_dict,
                                               added_atom_data = common_neighbor_dict['pair_count'],
                                               added_atom_property_str = added_atom_property_str,
                                               added_atom_property_columns_str = added_atom_property_columns_str)
    
    # delete empty lines of the file
    funcs.rm_empty_lines_in_file(simple_common_neighbor_file)
    funcs.rm_empty_lines_in_file(simple_common_neighbor_file_with_atom_indx)
    funcs.rm_empty_lines_in_file(simple_common_neighbor_file_with_atomspecies_atom_number_ortho)
    funcs.write_log(
        logfile,
        'vasp_analyze.simple_cna(' + '\n' +
        '    poscar_file_path = ' + 'r\'' + str(poscar_file_path) + '\'' + ',\n' +
        '    a0 = ' + str(a0) + ',\n' +
        '    common_neighbor_elmt_list = ' + str(common_neighbor_elmt_list) + ',\n' +
        '    atom_label = ' + str(atom_label) + ')\n' +
        '###############################\n')
    return common_neighbor_dict

def estruct(doscar_file_path, sysname = ''):
    '''
    - Descriptions
     * Calculates the structural energies at each atomic site
     * Reference: Chongyu Wang, Feng An, Gubing Lin, Fusui Liu, Ying Chen, Physical Review B, 1988

    - Args
     * doscar_file_path: String format. It specifies the directory of the DOSCAR
     * sysname: String format. User defined system name
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import convert
    from .. import default_params
    from . import vasp_read
    from . import vasp_write

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    valence_electron_num_dict = {'Al':3,
                             'Ti':4,'V':5,'Cr':6,'Mn':7,'Fe':8,'Co':9,'Ni':10,'Cu':11,'Zn':12,
                             'Zr':4,'Nb':5,'Mo':6,'Tc':7,'Ru':8,'Rh':9,'Pd':10,'Ag':11,'Cd':12,
                             'Hf':4,'Ta':5,'W':6,'Re':7,'Os':8,'Ir':9,'Pt':10,'Au':11,'Hg':12,
                             }

    d_s_valence_elmtname = ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
                           'Y','Zr','Nb','Mo','Tc','Ru','Rh','Ag','Cd',
                           'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg'
                           ]
    s_p_valence_elmtname = ['B','C','N','O','F','Ne',
                           'Al','Si','P','S','Cl','Ar',
                           ]
    s_valence_elmtname = ['H','Li','Be','Na','Mg','K','Ca','Rb','Sr','Cs','Ba']

    workdir, dos_file = funcs.file_path_name(doscar_file_path)
    outcar_file_path = os.path.join(workdir, 'OUTCAR')
    doscar_file_path = os.path.join(workdir, 'DOSCAR')
    incar_file_path = os.path.join(workdir, 'INCAR')
    poscar_file_path = os.path.join(workdir, 'POSCAR')
    poscar_dict = vasp_read.read_poscar(poscar_file_path)
    outcar_params_dict = vasp_read.read_outcar(outcar_file_path)
    LORBIT = outcar_params_dict['LORBIT']
    e_fermi = outcar_params_dict['e_fermi']
    if not outcar_params_dict['e_fermi'] is None:
        e_fermi_mod = outcar_params_dict['e_fermi'] + outcar_params_dict['alpha+bet']

    atom_num = len(poscar_dict['atom_species_arr'])

    estrut_arr = np.array([0.00000000] * atom_num)
    n1 = np.array([0.00000000] * atom_num)
    n2 = np.array([0.00000000] * atom_num)
    ne1 = 0
    ne2 = 0
    valence_electron_num_sum = 0
    eband = 0
    
    estruct_file = os.path.join(workdir, 'estruct_' + str(sysname) + '_Ef' + str('{:.4f}'.format(e_fermi_mod)) + '.txt')
    estruct_with_pos_file = os.path.join(workdir, 'estruct_with_position_' + str(sysname) + '_Ef' + str('{:.4f}'.format(e_fermi_mod)) + '.txt')

    formatted_atomname_len = 6 #For example len('Ni324') = 5, len('Ni1324') = 6
    formatted_len1 = 12 #For example len('-54.0812')=8, len('-54.08128345')=12    

    estruct_str = ''
    estruct_with_pos_str = ''
    estruct_str = estruct_str + ('atom' + ' ' + (' '*(formatted_atomname_len-len('atom'))) +
            'estruct\n')
    estruct_with_pos_str = estruct_with_pos_str + ('x' + ' ' + (' '*(formatted_len1-len('x'))) +
             'y' + ' ' + (' '*(formatted_len1-len('y'))) +
             'z' + ' ' + (' '*(formatted_len1-len('z'))) +
             'atom' + ' ' + (' '*(formatted_atomname_len-len('atom'))) +
             'estruct\n')

    doscar_dict = vasp_read.read_doscar(doscar_file_path, 1, False)
    NEDOS = len(doscar_dict['energy'])
    energy = doscar_dict['energy'] + outcar_params_dict['alpha+bet']
    e_grid = (np.max(energy)-np.min(energy)) / (NEDOS - 1)
    funcs.write_log(logfile, '## Modified Ef=' + str(e_fermi_mod) + ' eV' + '\n')
    funcs.write_log(logfile, '## NEDOS=' + str(NEDOS) + '\n')
    funcs.write_log(logfile, '## Emin=' + '{:.2f}'.format(np.min(energy)) + ' eV' +' Emax=' + '{:.2f}'.format(np.max(energy)) + ' eV' + '\n')
    funcs.write_log(logfile, '## dE=' + '{:.4f}'.format(e_grid) + ' eV' + '\n')
    dos_sum = np.array([0.0000000000]*NEDOS*atom_num)
    dos_sum.shape = NEDOS, atom_num

    for i_atom in range(0, atom_num):
        doscar_dict = vasp_read.read_doscar(doscar_file_path, i_atom + 1, False)
        ValenceElectronNum = valence_electron_num_dict[poscar_dict['atom_species_arr'][i_atom]]
        valence_electron_num_sum += ValenceElectronNum
        if doscar_dict['num_col'] == 10:
            #For transition elements:
            if poscar_dict['atom_species_arr'][i_atom] in d_s_valence_elmtname:
                dos_sum[:,i_atom] += np.abs(doscar_dict['dos_arr'][:,1])
                for i_col in range(5,10):
                    dos_sum[:,i_atom] += np.abs(doscar_dict['dos_arr'][:,i_col])
            elif poscar_dict['atom_species_arr'][i_atom] in s_p_valence_elmtname:
                dos_sum[:,i_atom] += np.abs(doscar_dict['dos_arr'][:,1])
                for i_col in range(2,5):
                    dos_sum[:,i_atom] += np.abs(doscar_dict['dos_arr'][:,i_col])
        if doscar_dict['num_col'] == 19:
            #For transition elements:
            if poscar_dict['atom_species_arr'][i_atom] in d_s_valence_elmtname:
                dos_sum[:,i_atom] += np.abs(doscar_dict['dos_arr'][:,1])
                for i_col in range(9,19):
                    dos_sum[:,i_atom] += np.abs(doscar_dict['dos_arr'][:,i_col])
            elif poscar_dict['atom_species_arr'][i_atom] in s_p_valence_elmtname:
                dos_sum[:,i_atom] += np.abs(doscar_dict['dos_arr'][:,1])
                for i_col in range(3,9):
                    dos_sum[:,i_atom] += np.abs(doscar_dict['dos_arr'][:,i_col])
        #Integrate over energy spectrum
        for j in range(0,NEDOS):
            if energy[j] < e_fermi_mod:
                estrut_arr[i_atom] += dos_sum[j,i_atom] * energy[j] * e_grid
                n1[i_atom] += dos_sum[j,i_atom] * e_grid
        ne1 += n1[i_atom] * 2
        eband += estrut_arr[i_atom]
        estruct_str = estruct_str + (poscar_dict['atomname_list'][i_atom] + (' '*(formatted_atomname_len-len(poscar_dict['atomname_list'][i_atom]))) +
                ' ' + '{:.4f}'.format(estrut_arr[i_atom]) + '\n')
        estruct_with_pos_str = estruct_with_pos_str + ('{:.4f}'.format(poscar_dict['pos_arr'][i_atom,3]) + ' ' + (' '*(formatted_len1-len('{:.4f}'.format(poscar_dict['pos_arr'][i_atom,3])))) +
                 '{:.4f}'.format(poscar_dict['pos_arr'][i_atom,4]) + ' ' + (' '*(formatted_len1-len('{:.4f}'.format(poscar_dict['pos_arr'][i_atom,4])))) +
                 '{:.4f}'.format(poscar_dict['pos_arr'][i_atom,5]) + ' ' + (' '*(formatted_len1-len('{:.4f}'.format(poscar_dict['pos_arr'][i_atom,5])))) +
                 poscar_dict['atomname_list'][i_atom] + ' ' + (' '*(formatted_atomname_len-len(poscar_dict['atomname_list'][i_atom]))) +
                 '{:.4f}'.format(estrut_arr[i_atom]) + '\n')
##        f.write(poscar_dict['atomname_list'][i_atom] + ' ' + '{:.4f}'.format(estrut_arr[i_atom]) + ' ' + '{:.4f}'.format(2*n1[i_atom]) + '\n')

##        #Integrate over energy spectrum: trapezoidal integration
##        NEDOS_BelowFermi = 0
##        for j in range(0,NEDOS):
##            if energy[j] < e_fermi_mod:
##                NEDOS_BelowFermi += 1
##        n_Mult_E = np.array([None]*NEDOS_BelowFermi)
##        for j in range(0,NEDOS):
##            if energy[j] < e_fermi_mod:
##                n_Mult_E[j] = dos_sum[j,i]         
##        n2[i] = np.trapz(n_Mult_E, dx = e_grid)
##        ne2 += n2[i] * 2

    with open(estruct_file,'w') as f, open(estruct_with_pos_file,'w') as f1:
        f.write(estruct_str)
        f1.write(estruct_with_pos_str)
    # write the POSCAR file with atom property data
    estruct_poscar_file_path = os.path.join(workdir, 'POSCAR_estruct_' + str(sysname) + '_Ef' + str('{:.4f}'.format(e_fermi_mod)) + '.vasp')
    added_atom_property_str = 'estruct'
    added_atom_property_columns_str = 'estruct'
    vasp_write.write_poscar(output_poscar_file_path = estruct_poscar_file_path,
                                               poscar_dict = poscar_dict,
                                               added_atom_data = estrut_arr,
                                               added_atom_property_str = added_atom_property_str,
                                               added_atom_property_columns_str = added_atom_property_columns_str)

    funcs.write_log(logfile, '## Calculated number of valence electrons=' + '{:.2f}'.format(ne1) + '\n')
    funcs.write_log(logfile, '## Actual number of valence electrons=' + str(valence_electron_num_sum) + '\n')
    funcs.write_log(logfile, '## eband=' + '{:.2f}'.format(eband) + ' eV' + '\n')
    funcs.write_log(
        logfile,
        'vasp_analyze.estruct(' + '\n' +
        '    doscar_file_path=' + 'r\'' + str(doscar_file_path) + '\'' + ',\n' +
        '    sysname=' + '\'' + str(sysname) + '\'' + ')\n' +
        '###############################\n')
##        NonZeroEiFound = False
##        for j in range(0,NEDOS):
##            for i in range(0, atom_num):
##              if dos_sum[j,i] != 0:
##                  E0 = energy[j]
##                  print('Lowest occupied energy state=', E0, 'eV')
##                  NonZeroEiFound = True
##                  break
##            if NonZeroEiFound == True:
##                break
##            
##        Nl_Temp = 0
##        for j in range(0,NEDOS):
##            for i in range(0, atom_num):
##                Nl_Temp += dos_sum[j,i] * e_grid * 2
##            if Nl_Temp >= valence_electron_num_sum:
##                print('Nl_Temp=',Nl_Temp)
##                print('iNEDOS=',j)
##                print('Estimated Efermi1=',energy[j], 'eV')
##                print('Estimated Efermi2=',E0+j*e_grid, 'eV')
##                break

    return 0

##########################################
#Below are modules for analyzation of DOS
##########################################

'''
Description:
DOS peak analyses for selected atoms with their neighboring atoms.
Mainly function is to find the overlapped orbitals and their corresponding energy levels.
'''

def label_peaks(energy, sample_dos_arr):
    '''
    Description:
    Roughly label the peaks
    sample_dos_arr is the DOS array which is to be labeled, for example, Sampledos_arr=d, or Sampledos_arr=dxy, or Sampledos_arr=s, etc.
    '''
    import numpy as np
    from scipy.signal import find_peaks
    def find_noise_baseline(interval_states_num, dos_arr_1d):
        '''
        description:
        Find the "noise baseline of the DOS curve, this is for the purpose of peak finding in DOS curve"
        interval_states_num: number of states included in the check interval.
        dos_arr_1d: One dimensional DOS array
        Return:
        noise_baseline: The values of the DOS wihch is larger than the noise_baseline will be used for peak finding
        '''
        dos_arr_1d_nonzero = dos_arr_1d[np.nonzero(dos_arr_1d)]
        if len(dos_arr_1d_nonzero) == 0:
            print('ERROR: vasp.analyze error. Nonzero array dos_arr_1d is []. This may be caused by the dos_mode setting. If an orbit does not exist for an element, this error may raise. For example, the d orbit does not exist for the Al atom. So the user should not designate the d orbital to the Al atom. In this case, please the orbits (dos_mode) for each element')
            exit()
        interval_num = int(np.ceil(np.max(dos_arr_1d_nonzero) / interval_states_num))
        #print(interval_num)
        interval_states_num_redefine = np.max(dos_arr_1d_nonzero) / interval_num
        states_num_list = []
        for i_interval in range(0,interval_num):
            interval_lo_val = i_interval * interval_states_num_redefine
            interval_hi_val = interval_lo_val + interval_states_num_redefine
            i_states_num = 0
            for i_state in range(0,len(dos_arr_1d_nonzero)):
                if dos_arr_1d_nonzero[i_state] >= interval_lo_val and dos_arr_1d_nonzero[i_state] < interval_hi_val:
                    i_states_num += 1
            states_num_list.append(i_states_num)
        #print(states_num_list)
        #print(np.argmax(states_num_list))
        #print(interval_states_num_redefine)
        
        #Sort the states num, if the neighboring value of states are within a range, then lower the states number
        temp_statex_indx = -1
        sorted_states_num = np.sort(states_num_list)[::-1]
        sorted_states_indx = np.argsort(states_num_list)[::-1]
        #print(sorted_states_num)
        for i_states_bin in range(0,len(sorted_states_num)):
            if i_states_bin == len(sorted_states_num)-2:
                continue
            a1 = sorted_states_num[i_states_bin]
            a2 = sorted_states_num[i_states_bin + 1]
            a3 = sorted_states_num[i_states_bin + 2]
            #print(np.abs((a1 - a2)/a1), np.abs((a2 - a3)/a2))
            if np.abs((a1 - a2)/a1) > 0.20:
                temp_statex_indx = sorted_states_indx[i_states_bin]
                break
            else:
                temp_statex_indx = sorted_states_indx[i_states_bin + 1]
        #print(temp_statex_indx, states_num_list[temp_statex_indx])
        #print('-----')
        noise_baseline = (temp_statex_indx * interval_states_num_redefine + interval_states_num_redefine)
        #noise_baseline = (np.argmax(states_num_list) * interval_states_num_redefine + 0.5* interval_states_num_redefine)
        return noise_baseline
    sample_dos_arr_1d = sample_dos_arr.reshape((1,len(sample_dos_arr)))[0]
    #the np.abs() function is used to turn the values of spin down occupation(intentionally set to negative occupation for better visulization and comparison) to positive values
    characteristic_states_num_sample_dos_arr_1d = 0.06 * np.max(np.abs(sample_dos_arr_1d))
    noise_baseline_sample_dos_arr = find_noise_baseline(characteristic_states_num_sample_dos_arr_1d, np.abs(sample_dos_arr_1d))
    peaks, properties = find_peaks(np.abs(sample_dos_arr_1d), height = noise_baseline_sample_dos_arr, prominence = characteristic_states_num_sample_dos_arr_1d)
    return energy[peaks], np.abs(sample_dos_arr[peaks])

def get_orbit_dos_peak_list(dos_arr, orbits, energy):
    import numpy as np
    energy_peak_list = []
    peak_val_list = []
    orbitname_list = []
    dos_temp_arr = np.array([None] * len(dos_arr[:,[0]]))
    dos_temp_arr.shape = len(dos_arr[:,[0]]), 1
    if dos_arr.shape[1] == 10 or dos_arr.shape[1] == 17:
        #Spin restricted calculation
        s = dos_arr[:,[1]] 
        py = dos_arr[:,[2]]
        pz = dos_arr[:,[3]]
        px = dos_arr[:,[4]]
        dxy = dos_arr[:,[5]] 
        dyz = dos_arr[:,[6]]
        dz2 = dos_arr[:,[7]]
        dxz = dos_arr[:,[8]]
        dx2 = dos_arr[:,[9]]
        p = py + pz + px
        d = dxy + dyz + dz2 + dxz + dx2
        ldos = s + p + d
        if 'd' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, d),axis = 0)
        if 'dxy' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dxy),axis = 0)
        if 'dyz' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dyz),axis = 0)
        if 'dz2' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dz2),axis = 0)
        if 'dxz' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dxz),axis = 0)
        if 'dx2' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dx2),axis = 0)
        if 'p' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, p),axis = 0)
        if 'py' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, py),axis = 0)
        if 'pz' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, pz),axis = 0)
        if 'px' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, px),axis = 0)
        if 's' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, s),axis = 0)
        if dos_arr.shape[1] == 17:
            fm3 = dos_arr[:,[10]]
            fm2 = dos_arr[:,[11]]
            fm1 = dos_arr[:,[12]]
            f0 = dos_arr[:,[13]]
            f1 = dos_arr[:,[14]]
            f2 = dos_arr[:,[15]]
            f3 = dos_arr[:,[16]]
            f = fm3 + fm2 + fm1 + f0 + f1 + f2 + f3
            ldos = ldos + f
            if 'f' in orbits:
                dos_temp_arr = np.concatenate((dos_temp_arr, f),axis = 0)
        if 'LDOS' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, ldos),axis = 0)
        dos_temp_arr = dos_temp_arr[len(dos_arr[:,[0]]):,0]

        if 'LDOS' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, ldos)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('LDOS')
        if 'f' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, f)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('f')
        if 'd' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, d)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('d')
        if 'p' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, p)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('p')
        if 's' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, s)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('s')
            
    elif dos_arr.shape[1] == 19 or dos_arr.shape[1] == 33:
        #Spin polarized calculation
        s_up = dos_arr[:,[1]] 
        s_dw = dos_arr[:,[2]]
        s = s_up - s_dw
        py_up = dos_arr[:,[3]] 
        py_dw = dos_arr[:,[4]]
        pz_up = dos_arr[:,[5]]
        pz_dw = dos_arr[:,[6]]
        px_up = dos_arr[:,[7]]
        px_dw = dos_arr[:,[8]]
        dxy_up = dos_arr[:,[9]]
        dxy_dw = dos_arr[:,[10]] 
        dyz_up = dos_arr[:,[11]] 
        dyz_dw = dos_arr[:,[12]]
        dz2_up = dos_arr[:,[13]]
        dz2_dw = dos_arr[:,[14]]
        dxz_up = dos_arr[:,[15]]
        dxz_dw = dos_arr[:,[16]]
        dx2_up = dos_arr[:,[17]]
        dx2_dw = dos_arr[:,[18]]
        p_up = py_up + pz_up + px_up
        p_dw = py_dw + pz_dw + px_dw
        p = p_up - p_dw
        d_up = dxy_up + dyz_up + dz2_up + dxz_up + dx2_up
        d_dw = dxy_dw + dyz_dw + dz2_dw + dxz_dw + dx2_dw
        d = d_up - d_dw
        ldos_up = s_up + p_up + d_up
        ldos_dw = s_dw + p_dw + d_dw
        ldos = s + p + d
        if 'd' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, d_up, d_dw),axis = 0)
        if 'dxy' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dxy_up, dxy_dw),axis = 0)
        if 'dyz' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dyz_up, dyz_dw),axis = 0)
        if 'dz2' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dz2_up, dz2_dw),axis = 0)
        if 'dxz' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dxz_up, dxz_dw),axis = 0)
        if 'dx2' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, dx2_up, dx2_dw),axis = 0)
        if 'p' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, p_up, p_dw),axis = 0)
        if 'py' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, py_up, py_dw),axis = 0)
        if 'pz' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, pz_up, pz_dw),axis = 0)
        if 'px' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, px_up, px_dw),axis = 0)
        if 's' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, s_up, s_dw),axis = 0)
        if dos_arr.shape[1] == 17:
            fm3_up = dos_arr[:,[19]]
            fm3_dw = dos_arr[:,[20]]
            fm2_up = dos_arr[:,[21]]
            fm2_dw = dos_arr[:,[22]]
            fm1_up = dos_arr[:,[23]]
            fm1_dw = dos_arr[:,[24]]
            f0_up = dos_arr[:,[25]]
            f0_dw = dos_arr[:,[26]]
            f1_up = dos_arr[:,[27]]
            f1_dw = dos_arr[:,[28]]
            f2_up = dos_arr[:,[29]] 
            f2_dw = dos_arr[:,[30]]
            f3_up = dos_arr[:,[31]]
            f3_dw = dos_arr[:,[32]]
            f_up = fm3_up + fm2_up + fm1_up + f0_up + f1_up + f2_up + f3_up
            f_dw = fm3_dw + fm2_dw + fm1_dw + f0_dw + f1_dw + f2_dw + f3_dw
            f = f_up - f_dw
            ldos_up = ldos_up + f_up
            ldos_dw = ldos_dw + f_dw
            ldos = ldos + f
            if 'f' in orbits:
                dos_temp_arr = np.concatenate((dos_temp_arr, f_up, f_dw),axis = 0)
        if 'LDOS' in orbits:
            dos_temp_arr = np.concatenate((dos_temp_arr, ldos_up, ldos_dw),axis = 0)
        dos_temp_arr = dos_temp_arr[len(dos_arr[:,[0]]):,0]

        if 'LDOS' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, ldos_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('ldos_up')
            energy_peaks, peaks_val = label_peaks(energy, ldos_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('ldos_dw')
        if 'f' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, f_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('f_up')
            energy_peaks, peaks_val = label_peaks(energy, f_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('f_dw')
        if 'd' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, d_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('d_up')
            energy_peaks, peaks_val = label_peaks(energy, d_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('d_dw')
        if 'dxy' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, dxy_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dxy_up')
            energy_peaks, peaks_val = label_peaks(energy, dxy_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dxy_dw')
        if 'dyz' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, dyz_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dyz_up')
            energy_peaks, peaks_val = label_peaks(energy, dyz_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dyz_dw')
        if 'dz2' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, dz2_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dz2_up')
            energy_peaks, peaks_val = label_peaks(energy, dz2_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dz2_dw')
        if 'dxz' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, dxz_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dxz_up')
            energy_peaks, peaks_val = label_peaks(energy, dxz_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dxz_dw')
        if 'dx2' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, dx2_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dx2_up')
            energy_peaks, peaks_val = label_peaks(energy, dx2_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('dx2_dw')
        if 'p' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, p_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('p_up')
            energy_peaks, peaks_val = label_peaks(energy, p_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('p_dw')
        if 'py' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, py_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('py_up')
            energy_peaks, peaks_val = label_peaks(energy, py_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('py_dw')
        if 'pz' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, pz_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('pz_up')
            energy_peaks, peaks_val = label_peaks(energy, pz_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('pz_dw')
        if 'px' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, px_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('px_up')
            energy_peaks, peaks_val = label_peaks(energy, px_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('px_dw')
        if 's' in orbits:
            energy_peaks, peaks_val = label_peaks(energy, s_up)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('s_up')
            energy_peaks, peaks_val = label_peaks(energy, s_dw)
            energy_peak_list.append(energy_peaks.reshape((1, len(energy_peaks)))[0])
            peak_val_list.append(peaks_val.reshape((1, len(peaks_val)))[0])
            orbitname_list.append('s_dw')
    return energy_peak_list, peak_val_list, orbitname_list

def overlap_peak_analyzer(doscar_file_path, sysname, atom_indx_list,n_shell,a0, dos_mode = None, fermi_shift_zero = True):
    '''
    procedure:
    Specify the atoms to be checked through
    Find the NN atoms of the selected atoms
    get dos of the atoms
    analyze the peaks
    find the overlapped peaks using a combinatorial cross check of the orbitals
    Determine the effective overlap energy(EOE):$E_{overlap}=\sum_{i}E_{i}=\sum_i{}n_{i}^{eff}E_{i} \delta(E−E_{i})$
    - Descriptions
     * Finding the overlapped orbitals of two neighboring atoms in the DOS analyis.
     * DOS peak analyses for selected atoms with their neighboring atoms.
     * Find the overlapped orbitals and their corresponding energy levels.	
     
    - Args
     * doscar_file_path: String format. The directory which contains the DOSCAR file, abstract path can be accepted
     * sysname: String format. A string character which specifies the name of the system, this string will be used as part of the output file name
     * atom_indx_list: List of strings. Specifies the list of selected atoms.
     * n_shell: float format. Up to which crystallographic shell(up to which nearest neighbor) of the selected atom will be considered
     * a0: float format. The approximate lattice constant of the crystal
     * dos_mode: dictionary format. Determins which orbital will be considered, f, d, p, s, dxy, dyz, ... can be used
     * fermi_shift_zero: A logical value determining whether the energy levels of the DOS will be shifted to zero
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    #Default values
    if dos_mode == None:
        dos_mode = defaults_dict['dos_mode']
    
    round_decimal = 2    #Rounding of the energy level
    max_neighbor_num = 50   #Maximum number of neighbors, please do not modify this variable unless you know what you are doing
    formatted_len = 12 #For example len('Ni324_-1-1-1')=12
    
    doscar_file_path = os.path.abspath(doscar_file_path)
    if sysname == None:
        sysname = ''

    workdir, dos_file = funcs.file_path_name(doscar_file_path)
    outcar_file_path = os.path.join(workdir, 'OUTCAR')
    doscar_file_path = os.path.join(workdir, 'DOSCAR')
    incar_file_path = os.path.join(workdir, 'INCAR')
    poscar_file_path = os.path.join(workdir, 'POSCAR')
    poscar_dict = vasp_read.read_poscar(poscar_file_path)
    outcar_params_dict = vasp_read.read_outcar(outcar_file_path)
    LORBIT = outcar_params_dict['LORBIT']
    e_fermi = outcar_params_dict['e_fermi']
    e_fermi_mod = e_fermi + outcar_params_dict['alpha+bet']

    overlap_peak_file = os.path.join(workdir, 'overlap_peak_' + str(sysname) + '.txt')
    eff_overlap_energy_file = os.path.join(workdir, 'effective_overlap_energy_' + str(sysname) + '.txt')
    overlap_peak_str = ''
    eff_overlap_energy_str = ''
    if fermi_shift_zero == True:
        overlap_peak_str = overlap_peak_str + 'e_fermi = 0\n'
    else:
        overlap_peak_str = overlap_peak_str + 'e_fermi = '+str(e_fermi_mod)+'\n'
    for i in range(0,len(atom_indx_list)):
        if isinstance(atom_indx_list[i],str) and atom_indx_list[i] != 'TDOS':
            atom_indx = convert.atomname2indx(os.path.join(workdir, 'POSCAR'),atom_indx_list[i])
            atomname_temp = str(atom_indx_list[i])
        elif isinstance(atom_indx_list[i],str) and atom_indx_list[i] == 'TDOS':
            atom_indx = 0
        elif isinstance(atom_indx_list[i],int):
            atom_indx = atom_indx_list[i]
            atomname_temp = str(poscar_dict['atom_species_arr'][atom_indx - 1])
        overlap_peak_str = overlap_peak_str + ('################################\n')
        overlap_peak_str = overlap_peak_str + (atomname_temp + '\n')
        eff_overlap_energy_str = eff_overlap_energy_str + ('################################\n')
        eff_overlap_energy_str = eff_overlap_energy_str + (atomname_temp + '\n')            
        orbits = dos_mode[str(poscar_dict['atom_species_arr'][atom_indx - 1])]
        doscar_dict = vasp_read.read_doscar(doscar_file_path, atom_indx, False)
        energy = doscar_dict['energy'] + outcar_params_dict['alpha+bet'] - e_fermi_mod*funcs.logic_retn_val(fermi_shift_zero,1,0)

        energy_peak_list, peak_val_list, orbitname_list = get_orbit_dos_peak_list(doscar_dict['dos_arr'], orbits, energy)
        eff_overlap_energy = np.array([0.000000000]*len(orbitname_list))
        #######################################
        # For the neighbors of the selected atom
        #######################################
        nn_atom_ishell_dict = nn_map(poscar_file_path, a0, n_shell)
        nn_atomname_ishell_arr = nn_atom_ishell_dict['atomname']
        nn_atomname_ishell_arr_no_shift_label = nn_atom_ishell_dict['atomname_no_shift_label']
        for i_shell in range(n_shell):
            overlap_peak_str = overlap_peak_str + ('\n' + str(i_shell+1) + 'NN Orbit(' + atomname_temp + ') Orbit(X) OverlappedPeaks...' + '\n')
            eff_overlap_energy_str = eff_overlap_energy_str + ('\n' + str(i_shell+1) + 'NN Orbit(' + atomname_temp + ')' + '\n')
            for i_neighbor in range(max_neighbor_num):
                if nn_atomname_ishell_arr_no_shift_label[atom_indx-1,i_shell,i_neighbor] == None:
                    continue
                i_nn_atomname = nn_atomname_ishell_arr_no_shift_label[atom_indx-1,i_shell,i_neighbor].strip()
                if isinstance(i_nn_atomname,str) and i_nn_atomname != 'TDOS':
                    i_nn_atom_indx = convert.atomname2indx(os.path.join(workdir, 'POSCAR'),i_nn_atomname)
                    iNN_atomname_temp = i_nn_atomname
                elif isinstance(i_nn_atomname,str) and i_nn_atomname == 'TDOS':
                    i_nn_atom_indx = 0
                elif isinstance(i_nn_atomname,int):
                    i_nn_atom_indx = i_nn_atomname
                    iNN_atomname_temp = str(poscar_dict['atom_species_arr'][atom_indx - 1])

                iNN_orbits = dos_mode[str(poscar_dict['atom_species_arr'][i_nn_atom_indx - 1])]
                
                i_nn_doscar_dict = vasp_read.read_doscar(doscar_file_path, i_nn_atom_indx, False)
                i_nn_dos_arr = i_nn_doscar_dict['dos_arr']
    
                i_nn_energy_peak_list, i_nn_peak_val_list, iNN_orbitname_list = get_orbit_dos_peak_list(i_nn_dos_arr, iNN_orbits, energy)

                for i_orbit in range(0,len(orbitname_list)):
                    for j_orbit in range(0,len(orbitname_list)):
                        energy_peaks1 = np.round(energy_peak_list[i_orbit],round_decimal)
                        energy_peaks2 = np.round(i_nn_energy_peak_list[j_orbit],round_decimal)
                        #intersection of two arrays
                        energy_peaks3,comm1,comm2 = np.intersect1d(energy_peaks1,energy_peaks2,return_indices=True)
                        sorted_energy_peaks3 = np.sort(energy_peaks3)
                        sorted_energy_peaks3_arg = np.argsort(energy_peaks3)
                        overlap_peak_val_list1 = peak_val_list[i_orbit][comm1]
                        overlap_peak_val_list2 = i_nn_peak_val_list[j_orbit][comm2]
                        Sortedoverlap_peak_val_list1 = overlap_peak_val_list1[sorted_energy_peaks3_arg]
                        Sortedoverlap_peak_val_list2 = overlap_peak_val_list2[sorted_energy_peaks3_arg]
                        overlap_peak_str = overlap_peak_str + (iNN_atomname_temp + ' ' + str(orbitname_list[i_orbit]) + ' ' + str(orbitname_list[j_orbit]) +
                                ' ' + ' '.join([str('{:.2f}'.format(i)) for i in sorted_energy_peaks3]) + '\n')
                        #Calculate effective overlap energy
                        if fermi_shift_zero == True:
                            sorted_energy_peaks3 = sorted_energy_peaks3 + e_fermi_mod*funcs.logic_retn_val(fermi_shift_zero,1,0)
                        else:
                            pass
                        for iPeak in range(0,len(sorted_energy_peaks3)):
                            if sorted_energy_peaks3[iPeak] <= e_fermi_mod:
##                                eff_overlap_energy[i_orbit] += sorted_energy_peaks3[iPeak] * np.min([Sortedoverlap_peak_val_list1[iPeak],Sortedoverlap_peak_val_list2[iPeak]])
                                eff_overlap_energy[i_orbit] += sorted_energy_peaks3[iPeak] * Sortedoverlap_peak_val_list1[iPeak]
        for i_orbit in range(0,len(orbitname_list)):
            eff_overlap_energy_str = eff_overlap_energy_str + str(orbitname_list[i_orbit]) + ' = '
            eff_overlap_energy_str = eff_overlap_energy_str + str(eff_overlap_energy[i_orbit]) + '\n'
        eff_overlap_energy_str = eff_overlap_energy_str + 'sum = ' + str(np.sum(eff_overlap_energy)) + '\n'
    with open(overlap_peak_file,'w') as f, open(eff_overlap_energy_file,'w') as f1:
        f.write(overlap_peak_str)
        f1.write(eff_overlap_energy_str)
    funcs.write_log(
        logfile,
        'dos_mode=' + str(dos_mode) + '\n' +
        'vasp_analyze.overlap_peak_analyzer(' + '\n' +
        '    doscar_file_path=' + 'r\'' + str(doscar_file_path) + '\'' + ',\n' +
        '    sysname=' + '\'' + str(sysname) + '\'' + ',\n' +
        '    atom_indx_list=' + str(atom_indx_list) + ',\n' +
        '    n_shell=' + str(n_shell) + ',\n' +
        '    a0=' + str(a0) + ',\n' +
        '    dos_mode=dos_mode' + ',\n' +
        '    fermi_shift_zero=' + str(fermi_shift_zero) + ')\n' +
        '###############################\n')
    return 0

def get_band_gap(eigenval_or_procar_dict, outcar_params_dict, kpoints_dict, fermi_shift_zero = True, suppress_warning = True):
    '''
    Get the band_gap, CBM and CVM
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import default_params
    from scipy import optimize
    from scipy.interpolate import interp1d

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    log_str = ''
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    #print(kpoints_dict)

    band_gap_dict = {}
    band_gap_dict['file_path'] = None
    band_gap_dict['band_gap'] = None
    band_gap_dict['gap_type'] = None
    band_gap_dict['CBM'] = None
    band_gap_dict['VBM'] = None

    band_gap_dict['band_gap_up'] = None
    band_gap_dict['gap_type_up'] = None
    band_gap_dict['CBM_up'] = None
    band_gap_dict['VBM_up'] = None

    band_gap_dict['band_gap_dw'] = None
    band_gap_dict['gap_type_dw'] = None
    band_gap_dict['CBM_dw'] = None
    band_gap_dict['VBM_dw'] = None

    band_gap_dict['file_path'] = eigenval_or_procar_dict['file_path'] 
    #get E_fermi
    e_fermi = outcar_params_dict['e_fermi']
    e_fermi_mod = outcar_params_dict['e_fermi_mod']
    ##if e_fermi not in [None, 'None', 'none'] and outcar_params_dict['alpha+bet'] not in [None, 'None', 'none']:
    ##    e_fermi_mod = e_fermi + outcar_params_dict['alpha+bet']
    #get KPOINTS 
    kpoints_scheme = kpoints_dict['scheme']
    #Check OUTCAR parameters
    valid_outcar_par = True
    temp_outcar_par_list = ['file_path', 'e_fermi', 'alpha+bet']
    for i_par in temp_outcar_par_list:
        if funcs.variables_equal(outcar_params_dict[i_par], outcar_params_dict['initial_value_dict'][i_par]):
            valid_outcar_par = False
            print('WARNING #2103311558 (from vasp_analyze). outcar_params_dict[\'' + str(i_par) + '\'] has an invalid value, please check the OUTCAR file ' + outcar_params_dict['file_path'] + '\n')

    if kpoints_scheme not in ['L', 'l']:
        if suppress_warning == True:
            pass
        elif suppress_warning ==False:
            print('WARNING #20120302 (from get_band_gap): KPOINTS file is not in line mode, this may not be a band structure calculation. File directory is : ' + kpoints_dict['file_path'])
    elif kpoints_scheme in ['L', 'l'] and valid_outcar_par:
        kpoints_arr = kpoints_dict['kpoints_xaxis_arr']
        kpath_num_intersections_interval = kpoints_dict['num_intersections_interval']
        kpath_num_intersections = kpoints_dict['num_intersections']

        general_params_dict = {}
        general_params_dict['fermi_shift_zero'] = fermi_shift_zero

        if eigenval_or_procar_dict['dict_type'] == 'eigenval':
            #get information from EIGENVAL
            ispin = eigenval_or_procar_dict['ispin']
            num_kpoints = eigenval_or_procar_dict['num_kpoints']
            num_bands = eigenval_or_procar_dict['num_bands'] 
            num_ions = eigenval_or_procar_dict['num_ions']
            band_gap_dict['num_bands'] = num_bands
            if ispin == 1:
                eigs = eigenval_or_procar_dict['eigs'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
                band_gap_dict['eigs'] = eigs
            elif ispin == 2:
                eigs_up = eigenval_or_procar_dict['eigs_up'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
                eigs_dw = eigenval_or_procar_dict['eigs_dw'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0)
                band_gap_dict['eigs_up'] = eigs_up
                band_gap_dict['eigs_dw'] = eigs_dw
        elif eigenval_or_procar_dict['dict_type'] == 'procar':
            ispin = eigenval_or_procar_dict['ispin']
            num_kpoints = eigenval_or_procar_dict['num_kpoints']
            num_bands = eigenval_or_procar_dict['num_bands']
            num_ions = eigenval_or_procar_dict['num_ions']
            band_gap_dict['num_bands'] = num_bands
            if ispin == 1:
                eigs = eigenval_or_procar_dict['eigs'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
                band_gap_dict['eigs'] = eigs
            elif ispin == 2:
                eigs_up = eigenval_or_procar_dict['eigs_up'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
                eigs_dw = eigenval_or_procar_dict['eigs_dw'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0)
                band_gap_dict['eigs_up'] = eigs_up
                band_gap_dict['eigs_dw'] = eigs_dw

        # band gap, CBM, VBM parameter initialization
        cbm_optimum = 9999999      #ispin = 1
        vbm_optimum = -9999999     #ispin = 1
        cbm = 9999999
        vbm = -9999999
        cbm_up = 9999999
        cbm_dw = 9999999
        vbm_up = -9999999
        vbm_dw = -9999999
        cbm_up_optimum = 9999999
        cbm_dw_optimum = 9999999
        vbm_up_optimum = -9999999
        vbm_dw_optimum = -9999999

        cbm_k_optimum = None   #ispin =1
        vbm_k_optimum = None   #ispin =1
        cbm_k = None
        vbm_k = None 
        cbm_k_up = None
        cbm_k_dw = None
        vbm_k_up = None
        vbm_k_dw = None
        cbm_k_up_optimum = None
        cbm_k_dw_optimum = None
        vbm_k_up_optimum = None
        vbm_k_dw_optimum = None

        band_gap = None
        band_gap_up = None
        band_gap_dw = None
        # get the current Fermi energy
        current_ef = e_fermi_mod - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0)
        # band indices of the valence band and the conduction band that are closest to the Fermi level
        hi_band_indx = 0
        lo_band_indx = 0
        hi_band_indx_up = 0
        lo_band_indx_up = 0
        hi_band_indx_dw = 0
        lo_band_indx_dw = 0

        max_eig = -9999999
        min_eig = 9999999
        max_eig_up = -9999999
        min_eig_up = 9999999
        max_eig_dw = -9999999
        min_eig_dw = 9999999

        for i_band in range(0, num_bands):
            if ispin ==1:
                band_arr = eigs[:, i_band]
                max_e = np.max(band_arr) 
                min_e = np.min(band_arr)
                if max_e >= current_ef and min_e <= current_ef:
                    hi_band_indx = i_band
                    lo_band_indx = i_band
                if max_e > current_ef and min_e > current_ef:
                    if min_e < min_eig:
                        min_eig = min_e
                        hi_band_indx = i_band
                if max_e < current_ef and min_e < current_ef:
                    if max_e > max_eig:
                        max_eig = max_e
                        lo_band_indx = i_band

                if len(kpoints_arr) != len(band_arr):
                    print('WARNING #2103311540 (from vasp_analyze.get_band_gap): length of kpoints_arr and band_arr are not identical. Please check EIGENVAL or PROCAR file in the directory: ' + outcar_params_dict['work_dir'])
                    break
                for i_interval in range(0, kpath_num_intersections_interval):
                    ##print('---------------------------')
                    ##print('i_interval=',i_interval)
                    ##print('kpath_num_intersections=', kpath_num_intersections)
                    start_indx = i_interval * kpath_num_intersections
                    end_indx = start_indx + kpath_num_intersections
                    kpoints_slice_arr = kpoints_arr[start_indx:end_indx]
                    eigs_slice_arr = band_arr[start_indx:end_indx]
                    ##print('len_band_arr=',len(band_arr),'len_kpoints_arr=',len(kpoints_arr))
                    ##print('start_indx,end_indx=',start_indx,end_indx)
                    ##print('kpoints_arr=',kpoints_arr)
                    ##print('band_arr=',band_arr)
                    ##print('eigs_slice_arr=',eigs_slice_arr)
                    kpoints_dense_arr = np.linspace(min(kpoints_slice_arr),max(kpoints_slice_arr),num = num_kpoints * 5,endpoint = True)
                    if len(kpoints_slice_arr) != len(eigs_slice_arr):
                        #print(kpoints_slice_arr, eigs_slice_arr)
                        #print('WARNING #2012131929 (from vasp_analyze.get_band_gap): length of kpoints_slice_arr and eigs_slice_arr do not conform. Please check EIGENVAL or PROCAR file in the directory: ' + outcar_params_dict['work_dir'])
                        break
                    interp = interp1d(kpoints_slice_arr, eigs_slice_arr, kind='cubic')

                    max_val = np.max(interp(kpoints_dense_arr))
                    min_val = np.min(interp(kpoints_dense_arr))

                    max_k = kpoints_dense_arr[np.argmax(interp(kpoints_dense_arr))]  
                    min_k = kpoints_dense_arr[np.argmin(interp(kpoints_dense_arr))]

                    # check if the band crosses the Fermi level
                    if max_val >= current_ef and min_val <= current_ef:
                        band_gap = 0
                    if max_val > current_ef and min_val > current_ef:
                        #bands above the Fermi level
                        cbm = min_val
                        cbm_k = min_k
                        #band_indx_valence_near_ef
                    elif max_val < current_ef and min_val < current_ef:
                        #bands below the Fermi level
                        vbm = max_val
                        vbm_k = max_k

                    cbm_optimum = np.min([cbm, cbm_optimum])
                    vbm_optimum = np.max([vbm, vbm_optimum])

                    temp_arg = np.argmin([cbm, cbm_optimum])
                    if temp_arg == 0:
                        cbm_k_optimum = cbm_k
                    temp_arg = np.argmin([vbm, cbm_optimum])
                    if temp_arg == 0:
                        vbm_k_optimum = vbm_k
            elif ispin == 2:
                band_arr_up = eigs_up[:, i_band]
                band_arr_dw = eigs_dw[:, i_band]
                max_up = np.max(band_arr_up) 
                min_up = np.min(band_arr_up)
                max_dw = np.max(band_arr_dw) 
                min_dw = np.min(band_arr_dw)
                if max_up >= current_ef and min_up <= current_ef:
                    hi_band_indx_up = i_band
                    lo_band_indx_up = i_band
                if max_dw >= current_ef and min_dw <= current_ef:
                    hi_band_indx_dw = i_band
                    lo_band_indx_dw = i_band
                if max_up > current_ef and min_up > current_ef:
                    if min_up < min_eig_up:
                        min_eig_up = min_up
                        hi_band_indx_up = i_band
                if max_up < current_ef and min_up < current_ef:
                    if max_up > max_eig_up:
                        max_eig_up = max_up
                        lo_band_indx_up = i_band
                if max_dw > current_ef and min_dw > current_ef:
                    if min_dw < min_eig_dw:
                        min_eig_dw = min_dw
                        hi_band_indx_dw = i_band
                if max_dw < current_ef and min_dw < current_ef:
                    if max_dw > max_eig_dw:
                        max_eig_dw = max_dw
                        lo_band_indx_dw = i_band
                
                for i_interval in range(0, kpath_num_intersections_interval):
                    start_indx = i_interval * kpath_num_intersections
                    end_indx = start_indx + kpath_num_intersections
                    kpoints_slice_arr = kpoints_arr[start_indx:end_indx]
                    eigs_slice_arr_up = band_arr_up[start_indx:end_indx]
                    eigs_slice_arr_dw = band_arr_dw[start_indx:end_indx]
                    kpoints_dense_arr = np.linspace(min(kpoints_slice_arr),max(kpoints_slice_arr),num = num_kpoints * 5,endpoint = True)
                    if len(kpoints_slice_arr) != len(eigs_slice_arr_up):
                        #print('WARNING #2012131927 (from vasp_analyze.get_band_gap): length of kpoints_slice_arr and eigs_slice_arr_up do not conform. Please check EIGENVAL or PROCAR file.')
                        break
                    interp_up = interp1d(kpoints_slice_arr, eigs_slice_arr_up, kind='cubic')
                    interp_dw = interp1d(kpoints_slice_arr, eigs_slice_arr_dw, kind='cubic')

                    max_val_up = np.max(interp_up(kpoints_dense_arr))
                    min_val_up = np.min(interp_up(kpoints_dense_arr))
                    max_val_dw = np.max(interp_dw(kpoints_dense_arr))
                    min_val_dw = np.min(interp_dw(kpoints_dense_arr))

                    max_k_up = kpoints_dense_arr[np.argmax(interp_up(kpoints_dense_arr))]  
                    min_k_up = kpoints_dense_arr[np.argmin(interp_up(kpoints_dense_arr))]
                    max_k_dw = kpoints_dense_arr[np.argmax(interp_dw(kpoints_dense_arr))]
                    min_k_dw = kpoints_dense_arr[np.argmin(interp_dw(kpoints_dense_arr))]

                    # check if the band crosses the Fermi level
                    if max_val_up >= current_ef and min_val_up <= current_ef:
                        band_gap_up = 0                        
                        band_gap = 0
                    if max_val_dw >= current_ef and min_val_dw <= current_ef:
                        band_gap_dw = 0                        
                        band_gap = 0
                    # for spin up channel
                    if max_val_up > current_ef and min_val_up > current_ef:
                        #bands above the Fermi level
                        cbm_up = min_val_up
                        cbm_k_up = min_k_up
                    elif max_val_up < current_ef and min_val_up < current_ef:
                        #bands below the Fermi level
                        vbm_up = max_val_up
                        vbm_k_up = max_k_up
                    # for spin down channel
                    if max_val_dw > current_ef and min_val_dw > current_ef:
                        #bands above the Fermi level
                        cbm_dw = min_val_dw
                        cbm_k_dw = min_k_dw
                    elif max_val_dw < current_ef and min_val_dw < current_ef:
                        #bands below the Fermi level
                        vbm_dw = max_val_dw
                        vbm_k_dw = max_k_dw

                    cbm = np.min([cbm_up, cbm_dw, cbm])
                    vbm = np.max([vbm_up, vbm_dw, vbm])
                    cbm_up_optimum = np.min([cbm_up, cbm_up_optimum])
                    vbm_up_optimum = np.max([vbm_up, vbm_up_optimum])
                    cbm_dw_optimum = np.min([cbm_dw, cbm_dw_optimum])
                    vbm_dw_optimum = np.max([vbm_dw, vbm_dw_optimum])

                    # consider both spin up and spin down channel
                    temp_arg = np.argmin([cbm_up, cbm_dw, cbm])
                    if temp_arg == 0:
                        cbm_k = cbm_k_up
                    elif temp_arg == 1:
                        cbm_k = cbm_k_dw
                    temp_arg = np.argmax([vbm_up, vbm_dw, vbm])
                    if temp_arg == 0:
                        vbm_k = vbm_k_up
                    elif temp_arg == 1:
                        vbm_k = vbm_k_dw
                    #only consider spin up channel 
                    temp_arg_up = np.argmin([cbm_up, cbm_up_optimum])
                    if temp_arg_up == 0:
                        cbm_k_up_optimum = cbm_k_up
                    temp_arg_up = np.argmax([vbm_up, vbm_up_optimum])
                    if temp_arg_up == 0:
                        vbm_k_up_optimum = vbm_k_up
                    #only consider spin down channel 
                    temp_arg_dw = np.argmin([cbm_dw, cbm_dw_optimum])
                    if temp_arg_dw == 0:
                        cbm_k_dw_optimum = cbm_k_dw
                    temp_arg_dw = np.argmax([vbm_dw, vbm_dw_optimum])
                    if temp_arg_dw == 0:
                        vbm_k_dw_optimum = vbm_k_dw
        #band gap information
        if ispin == 1:
            if band_gap == 0:
                # band gap does not exists
                gap_type = None
                cbm_optimum = None
                vbm_optimum = None
            if band_gap != 0:
                # band gap exists
                band_gap = cbm_optimum - vbm_optimum
                if cbm_k_optimum == vbm_k_optimum:
                    gap_type = 'direct'
                elif cbm_k_optimum != vbm_k_optimum:
                    gap_type = 'indirect'
            band_gap_dict['band_gap'] = band_gap
            band_gap_dict['gap_type'] = gap_type
            band_gap_dict['CBM'] = cbm_optimum
            band_gap_dict['VBM'] = vbm_optimum
            band_gap_dict['kpoint_CBM'] = cbm_k_optimum
            band_gap_dict['kpoint_VBM'] = vbm_k_optimum
            band_gap_dict['CBM_band_indx'] = hi_band_indx
            band_gap_dict['VBM_band_indx'] = lo_band_indx
            #print('---',band_gap_dict['VBM_band_indx'], band_gap_dict['CBM_band_indx'])
        elif ispin == 2:
            if band_gap == 0:
                # band gap does not exists
                gap_type = None
                cbm = None
                vbm = None
            if band_gap_up == 0:
                # band gap for spin up channel does not exists
                gap_type_up = None
                cbm_up_optimum = None
                vbm_up_optimum = None
            if band_gap_dw == 0:
                # band gap for spin down channel does not exists
                gap_type_dw = None
                cbm_dw_optimum = None
                vbm_dw_optimum = None
            if band_gap != 0:
                # band gap exists
                band_gap = cbm - vbm
                if cbm_k == vbm_k:
                    gap_type = 'direct'
                elif cbm_k != vbm_k:
                    gap_type = 'indirect'
            if band_gap_up != 0:
                # band gap exists for spin up channel
                band_gap_up = cbm_up_optimum - vbm_up_optimum
                if cbm_k_up_optimum == vbm_k_up_optimum:
                    gap_type_up = 'direct'
                elif cbm_k_up_optimum != vbm_k_up_optimum:
                    gap_type_up = 'indirect'
            if band_gap_dw != 0:
                # band gap exists for spin down channel
                band_gap_dw = cbm_dw_optimum - vbm_dw_optimum
                if cbm_k_dw_optimum == vbm_k_dw_optimum:
                    gap_type_dw = 'direct'
                elif cbm_k_dw_optimum != vbm_k_dw_optimum:
                    gap_type_dw = 'indirect'
            band_gap_dict['band_gap'] = band_gap
            band_gap_dict['gap_type'] = gap_type
            band_gap_dict['CBM'] = cbm
            band_gap_dict['VBM'] = vbm
            band_gap_dict['kpoint_CBM'] = cbm_k
            band_gap_dict['kpoint_VBM'] = vbm_k
            band_gap_dict['band_gap_up'] = band_gap_up
            band_gap_dict['band_gap_dw'] = band_gap_dw
            band_gap_dict['gap_type'] = gap_type
            band_gap_dict['gap_type_up'] = gap_type_up
            band_gap_dict['gap_type_dw'] = gap_type_dw
            band_gap_dict['CBM_up'] = cbm_up_optimum
            band_gap_dict['VBM_up'] = vbm_up_optimum
            band_gap_dict['CBM_dw'] = cbm_dw_optimum
            band_gap_dict['VBM_dw'] = vbm_dw_optimum
            band_gap_dict['kpoint_CBM_up'] = cbm_k_up_optimum
            band_gap_dict['kpoint_VBM_up'] = vbm_k_up_optimum
            band_gap_dict['kpoint_CBM_dw'] = cbm_k_dw_optimum
            band_gap_dict['kpoint_VBM_dw'] = vbm_k_dw_optimum
            band_gap_dict['CBM_band_indx_up'] = hi_band_indx_up
            band_gap_dict['VBM_band_indx_up'] = lo_band_indx_up
            band_gap_dict['CBM_band_indx_dw'] = hi_band_indx_dw
            band_gap_dict['VBM_band_indx_dw'] = lo_band_indx_dw
            #print('---',band_gap_dict['VBM_band_indx_up'], band_gap_dict['CBM_band_indx_up'])
    return band_gap_dict
