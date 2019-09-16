# -*- coding: utf-8 -*-
def nn_map(poscar_dir, a0, n_shell = 1):
    '''
    Description:
        Nearest neighbor map. Now only availavle for face-centered based structures
    - Args
     * poscar_dir: String format. It specifies the directory of the POSCAR file
     * a0: Float format. The lattice constant of the model. Unit in Angstrom
     * n_shell: Integer format. It determines up to which crystallographic shell the nearest neighbour map calculates
    '''
    import os
    import math
    import numpy as np
    from sklearn.neighbors import KDTree
    from .. import funcs
    from . import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    poscar_dir = os.path.abspath(poscar_dir)
    n_neighbor = 100  #number of nearest neighbors in kNNList
    tol = 0.1 * a0  #tol is the maximum allowed deviation of atom posotions from perfect lattice site.

    #Define lattice structure
    structure = "FCC"
    r = np.array([0]*n_shell,dtype=np.float)
    if structure == "FCC":
        first_nn = math.sqrt(2.0) / 2 * a0;
        for i_shell in range(0,n_shell):
            r[i_shell] = first_nn * math.sqrt(i_shell + 1.0);
            
    workdir, poscar_file = funcs.file_path_name(poscar_dir)
    poscar_dict = vasp_read.read_poscar(poscar_dir)
    n_atoms = np.sum(poscar_dict['elmt_num_arr'])  
    '''Find atom neighbors over periodic boundaries:
    Eextend/duplicate the cell in three dimensions, Create duplicates in three dimensions.
    Define the duplicity in each direction(include positive and negative directions) as 'Extend'.
    In each dimension, there are (2*n_extend+1) number of cells. (2*n_extend) of the (2*n_extend+1) cells are duplicates.
    There are (2*n_extend+1)^3-1 duplicates created. The volume is (2*n_extend+1)^3 times larger'''
    n_extend = 1 # Define the cell duplicity in each direction.
    atomname_extend_list = poscar_dict['atomname_list'].copy()
    atomname_extend_list_no_extend_atomname = poscar_dict['atomname_list'].copy()
    car_extend_arr = np.array([0.0]*n_atoms*pow(2*n_extend+1,3)*3,dtype = np.float)
    car_extend_arr.shape = n_atoms*pow(2*n_extend+1,3),3
    formatted_atomname_len = 5 #For example len('Ni324') = 5
    formatted_len = 12 #For example len('Ni324_-1-1-1')=12
    for i_atom in range(n_atoms):
        car_extend_arr[i_atom,0] = poscar_dict['pos_arr'][i_atom,3]
        car_extend_arr[i_atom,1] = poscar_dict['pos_arr'][i_atom,4]
        car_extend_arr[i_atom,2] = poscar_dict['pos_arr'][i_atom,5]
        atomname_extend_list[i_atom] = atomname_extend_list[i_atom] + ' '*(formatted_len - len(atomname_extend_list[i_atom]))    #This is for the purpose of formatted
    i_atom_extend = i_atom
    unit_vec = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
    for i in range(-n_extend, n_extend+1):
        for j in range(-n_extend, n_extend+1):
            for k in range(-n_extend, n_extend+1):
                if i == 0 and j == 0 and k == 0:
                    continue
                for i_atom in range(n_atoms):
                    i_atom_extend += 1
                    car_extend_arr[i_atom_extend,0] = poscar_dict['pos_arr'][i_atom,3] + np.dot(poscar_dict['l_arr'][:,0], unit_vec[0,:]) * i
                    car_extend_arr[i_atom_extend,1] = poscar_dict['pos_arr'][i_atom,4] + np.dot(poscar_dict['l_arr'][:,1], unit_vec[1,:]) * j
                    car_extend_arr[i_atom_extend,2] = poscar_dict['pos_arr'][i_atom,5] + np.dot(poscar_dict['l_arr'][:,2], unit_vec[2,:]) * k
                    extend_atomname = poscar_dict['atomname_list'][i_atom] + '_' + str(i) + str(j) + str(k)
                    atomname_extend_list.append(extend_atomname + ' '*(formatted_len - len(extend_atomname)))
                    atomname_extend_list_no_extend_atomname.append(poscar_dict['atomname_list'][i_atom])
    #Calculate nearest neighbor inofrmation
    tree = KDTree(car_extend_arr)
    output_dir = os.getcwd() + '/outputs'
    funcs.mkdir(output_dir)    
    nn_atomname_list_file = output_dir + '/nn_atomname_list.dat'
    nn_atomname_list_file_without_mirror_label = output_dir + '/nn_atomname_list.dat'
    nn_dist_file = output_dir + '/nn_dist.dat'
    nn_map_temp_file = output_dir + '/nn_map_temp.dat'
    nn_map_temp_file_without_mirror_label = output_dir + '/nn_map_temp_without_mirror_label.dat'
    nn_map_file = output_dir + '/nn_map.dat'
    nn_map_file_without_mirror_label = output_dir + '/nn_map_without_mirror_label.dat'
    with open(nn_atomname_list_file,'w') as f1, open(nn_dist_file,'w') as f2:
        pass
    with open(nn_atomname_list_file_without_mirror_label,'w') as f11:
        pass
    nn_atomname_list_i_shell_file = []
    nn_atomname_list_i_shell_file_without_mirror_label = []
    nn_dist_list_ishell_file = []
    max_neighbor_num = 50
    #nn_atomname_ishell_arr constains the atoms names with mirror atom label
    #nn_atomname_ishell_arr_no_extend_atomname constains the atoms names without mirror atom label    
    nn_atomname_ishell_arr = np.array([None]*n_atoms*n_shell*max_neighbor_num)
    nn_atomname_ishell_arr.shape = n_atoms,n_shell,max_neighbor_num
    nn_atomname_ishell_arr_no_extend_atomname = np.array([None]*n_atoms*n_shell*max_neighbor_num)
    nn_atomname_ishell_arr_no_extend_atomname.shape = n_atoms,n_shell,max_neighbor_num
    for i_shell in range(n_shell):
        shell_indx = i_shell +1
        nn_atomname_list_i_shell_file.append(output_dir + '/nn_atomname_list_' + str(shell_indx) + 'NN.dat')
        nn_atomname_list_i_shell_file_without_mirror_label.append(output_dir + '/nn_atomname_list_without_mirror_label_' + str(shell_indx) + 'NN.dat')
        nn_dist_list_ishell_file.append(output_dir + '/nn_dist_list_' + str(shell_indx) + 'NN.dat')
        with open(nn_atomname_list_i_shell_file[i_shell],'w') as f3, open(nn_dist_list_ishell_file[i_shell],'w') as f4:
            pass
        with open(nn_atomname_list_i_shell_file_without_mirror_label[i_shell],'w') as f13:
            pass
    with open(nn_map_file,'w') as f:
        f.write('atom ' + ' '.join([str(i+1) + 'NN' for i in range(n_shell)]) + '\n')
    with open(nn_map_file_without_mirror_label,'w') as f10:
        f10.write('atom ' + ' '.join([str(i+1) + 'NN' for i in range(n_shell)]) + '\n')

    with open(nn_atomname_list_file,'a') as f1, open(nn_dist_file,'a') as f2, open(nn_atomname_list_file_without_mirror_label,'a') as f11:
        for i_atom in range(n_atoms):
            #dists, indices = tree.query([car_extend_arr[i_atom_extend]], k=[x for x in range(1,n_neighbor+1)])
            dists, indices = tree.query([car_extend_arr[i_atom]], k=n_neighbor+1)
            nn_atomname_list = []
            nn_dist_list = []
            for i_neighbor in range(n_neighbor+1):
                nn_atomname_list.append(atomname_extend_list[indices[0,i_neighbor]])
                nn_dist_list.append(str(dists[0,i_neighbor]))
            f1.write(str(', '.join(nn_atomname_list))+'\n')
            f11.write(str(', '.join(nn_atomname_list))+'\n')
            f2.write(str(atomname_extend_list[i_atom]) + ', ' + str(", ".join(nn_dist_list))+'\n')
            for i_shell in range(n_shell):
                nn_atomname_ishell_list = []
                nn_atomname_ishell_list_no_extend_atomname = []
                nn_dist_ishell_list = []
                for i_neighbor in range(n_neighbor+1):
                    if abs(dists[0,i_neighbor]-r[i_shell]) < tol:
                        nn_atomname_ishell_list.append(atomname_extend_list[indices[0,i_neighbor]])
                        nn_atomname_ishell_list_no_extend_atomname.append(atomname_extend_list_no_extend_atomname[indices[0,i_neighbor]])
                        nn_dist_ishell_list.append(str(dists[0,i_neighbor]))
                nn_atomname_ishell_arr[i_atom,i_shell,0:len(nn_atomname_ishell_list)] =  nn_atomname_ishell_list
                nn_atomname_ishell_arr_no_extend_atomname[i_atom,i_shell,0:len(nn_atomname_ishell_list)] =  nn_atomname_ishell_list_no_extend_atomname
                with open(nn_atomname_list_i_shell_file[i_shell],'a') as f3, open(nn_dist_list_ishell_file[i_shell],'a') as f4, open(nn_atomname_list_i_shell_file_without_mirror_label[i_shell],'a') as f13:
                    f3.write(str(atomname_extend_list[i_atom]) + ', ' + str(', '.join(nn_atomname_ishell_list))+'\n')
                    f13.write(str(atomname_extend_list[i_atom]) + ', ' + str(', '.join(nn_atomname_ishell_list_no_extend_atomname))+'\n')
                    f4.write(str(atomname_extend_list[i_atom]) + ', ' + str(', '.join(nn_dist_ishell_list))+'\n')
    for i_atom in range(n_atoms):
        for i_shell in range(n_shell):
            for i_neighbor in range(max_neighbor_num):
                if nn_atomname_ishell_arr[i_atom,i_shell,i_neighbor] == None:
                    nn_atomname_ishell_arr[i_atom,i_shell,i_neighbor] = ' '*formatted_len                  
                if nn_atomname_ishell_arr_no_extend_atomname[i_atom,i_shell,i_neighbor] == None:
                    nn_atomname_ishell_arr_no_extend_atomname[i_atom,i_shell,i_neighbor] = ' '*formatted_len                  
    with open(nn_map_temp_file,'w') as f5:
        f5.write('atom  ' + ' '.join([str(i+1) + 'NN         ' for i in range(n_shell)]) + '\n')
    with open(nn_map_temp_file_without_mirror_label,'w') as f15:
        f15.write('atom  ' + ' '.join([str(i+1) + 'NN         ' for i in range(n_shell)]) + '\n')
    with open(nn_map_temp_file,'a') as f5, open(nn_map_temp_file_without_mirror_label,'a') as f15:
        for i_atom in range(n_atoms):
            split_text = ' '*formatted_atomname_len + ' '
            nn_text = split_text.join([' '.join(i) + '\n' for i in [nn_atomname_ishell_arr[i_atom,0:n_shell,i_neighbor].tolist() for i_neighbor in range(max_neighbor_num)]])
            nn_text_without_mirror_label = split_text.join([' '.join(i) + '\n' for i in [nn_atomname_ishell_arr_no_extend_atomname[i_atom,0:n_shell,i_neighbor].tolist() for i_neighbor in range(max_neighbor_num)]])
            f5.write(str(poscar_dict['atomname_list'][i_atom])+ ' '*(formatted_atomname_len - len(poscar_dict['atomname_list'][i_atom])) + ' ' + nn_text + '\n')
            f15.write(str(poscar_dict['atomname_list'][i_atom])+ ' '*(formatted_atomname_len - len(poscar_dict['atomname_list'][i_atom])) + ' ' + nn_text_without_mirror_label + '\n')

    with open(nn_map_file,'w') as f6:
        pass
    with open(nn_map_file_without_mirror_label,'w') as f16:
        pass
    with open(nn_map_temp_file,'r') as f5, open(nn_map_file,'a') as f6, open(nn_map_temp_file_without_mirror_label,'r') as f15, open(nn_map_file_without_mirror_label,'a') as f16:
        line = f5.readlines()
        line1 = f15.readlines()
        for i in range(len(line)):
            if str(line[i]).replace('\t', '    ').strip('\n').strip() != '':
                f6.write(str(line[i]))
            if str(line1[i]).replace('\t', '    ').strip('\n').strip() != '':
                f16.write(str(line1[i]))
        
    for i_atom in range(n_atoms):
        for i_shell in range(n_shell):
            for i_neighbor in range(max_neighbor_num):
                if nn_atomname_ishell_arr[i_atom,i_shell,i_neighbor] == ' '*formatted_len:
                    nn_atomname_ishell_arr[i_atom,i_shell,i_neighbor] = None              
                if nn_atomname_ishell_arr_no_extend_atomname[i_atom,i_shell,i_neighbor] == ' '*formatted_len:
                    nn_atomname_ishell_arr_no_extend_atomname[i_atom,i_shell,i_neighbor] = None               

    funcs.write_log(logfile, 'nn_map:' + '\n' + 'n_shell=' + str(n_shell) + '\n' + 'a0=' + str(a0) + '\n' + 'POSCARDir=' + poscar_dir + '\n' + '#############################################')
    return nn_atomname_ishell_arr, nn_atomname_ishell_arr_no_extend_atomname


def estruct(doscar_dir, sysname = ''):
    '''
    - Descriptions
     * Calculates the structural energies at each atomic site
     * Reference: Chongyu Wang, Feng An, Gubing Lin, Fusui Liu, Ying Chen, Physical Review B, 1988

    - Args
     * doscar_dir: String format. It specifies the directory of the DOSCAR
     * sysname: String format. User defined system name
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
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

    workdir, dos_file = funcs.file_path_name(doscar_dir)
    outcar_dir = workdir + '/OUTCAR'
    doscar_dir = workdir + '/DOSCAR'
    incar_dir = workdir + '/INCAR'
    poscar_dir = workdir + '/POSCAR'
    poscar_dict = vasp_read.read_poscar(poscar_dir)
    incar_params_dict, outcar_params_dict = vasp_read.read_outcar(outcar_dir)
    LORBIT = incar_params_dict['LORBIT']
    e_fermi = outcar_params_dict['e_fermi']
    e_fermi_mod = e_fermi + outcar_params_dict['alpha+bet']

    atom_num = len(poscar_dict['atom_species_arr'])

    estru = np.array([0.00000000] * atom_num)
    n1 = np.array([0.00000000] * atom_num)
    n2 = np.array([0.00000000] * atom_num)
    ne1 = 0
    ne2 = 0
    valence_electron_num_sum = 0
    eband = 0

    output_dir = os.getcwd() + '/outputs'
    funcs.mkdir(output_dir)
    
    estruct_file = output_dir + '/estruct_' + str(sysname) + '_Ef' + str('{:.4f}'.format(e_fermi_mod)) + '.txt'
    estruct_with_pos_file = output_dir + '/estruct_with_position_' + str(sysname) + '_Ef' + str('{:.4f}'.format(e_fermi_mod)) + '.txt'

    with open(estruct_file,'w') as f, open(estruct_with_pos_file,'w') as f1:
        f.write('atom_name estruct\n')
        f1.write('x y z atom_name estruct\n')
    with open(estruct_file,'a') as f, open(estruct_with_pos_file,'a') as f1:
        doscar_dict = vasp_read.read_doscar(doscar_dir, 1, False)
        NEDOS = len(doscar_dict['energy'])
        energy = doscar_dict['energy'] + outcar_params_dict['alpha+bet']
        e_grid = (np.max(energy)-np.min(energy)) / (NEDOS - 1)
        funcs.write_log(logfile, '#################################################\n' + 'estruct\n' + '#################################################\n')
        funcs.write_log(logfile, 'Modified Ef=' + str(e_fermi_mod) + ' eV' + '\n')
        funcs.write_log(logfile, 'NEDOS=' + str(NEDOS) + '\n')
        funcs.write_log(logfile, 'Emin=' + '{:.2f}'.format(np.min(energy)) + ' eV' +' Emax=' + '{:.2f}'.format(np.max(energy)) + ' eV' + '\n')
        funcs.write_log(logfile, 'dE=' + '{:.4f}'.format(e_grid) + ' eV' + '\n')
        dos_sum = np.array([0.0000000000]*NEDOS*atom_num)
        dos_sum.shape = NEDOS, atom_num

        for i in range(0, atom_num):
            doscar_dict = vasp_read.read_doscar(doscar_dir, i + 1, False)
            ValenceElectronNum = valence_electron_num_dict[poscar_dict['atom_species_arr'][i]]
            valence_electron_num_sum += ValenceElectronNum
            if doscar_dict['num_col'] == 10:
                #For transition elements:
                if poscar_dict['atom_species_arr'][i] in d_s_valence_elmtname:
                    dos_sum[:,i] += np.abs(doscar_dict['dos_arr'][:,1])
                    for i_col in range(5,10):
                        dos_sum[:,i] += np.abs(doscar_dict['dos_arr'][:,i_col])
                elif poscar_dict['atom_species_arr'][i] in s_p_valence_elmtname:
                    dos_sum[:,i] += np.abs(doscar_dict['dos_arr'][:,1])
                    for i_col in range(2,5):
                        dos_sum[:,i] += np.abs(doscar_dict['dos_arr'][:,i_col])
            if doscar_dict['num_col'] == 19:
                #For transition elements:
                if poscar_dict['atom_species_arr'][i] in d_s_valence_elmtname:
                    dos_sum[:,i] += np.abs(doscar_dict['dos_arr'][:,1])
                    for i_col in range(9,19):
                        dos_sum[:,i] += np.abs(doscar_dict['dos_arr'][:,i_col])
                elif poscar_dict['atom_species_arr'][i] in s_p_valence_elmtname:
                    dos_sum[:,i] += np.abs(doscar_dict['dos_arr'][:,1])
                    for i_col in range(3,9):
                        dos_sum[:,i] += np.abs(doscar_dict['dos_arr'][:,i_col])
            #Integrate over energy spectrum
            for j in range(0,NEDOS):
                if energy[j] < e_fermi_mod:
                    estru[i] += dos_sum[j,i] * energy[j] * e_grid
                    n1[i] += dos_sum[j,i] * e_grid
            ne1 += n1[i] * 2
            eband += estru[i]
            f.write(poscar_dict['atomname_list'][i] + ' ' + '{:.4f}'.format(estru[i]) + '\n')
            f1.write('{:.4f}'.format(poscar_dict['pos_arr'][i,3]) + ' ' + '{:.4f}'.format(poscar_dict['pos_arr'][i,4]) + ' ' + '{:.4f}'.format(poscar_dict['pos_arr'][i,5]) + ' ' + poscar_dict['atomname_list'][i] + ' ' + '{:.4f}'.format(estru[i]) + '\n')
##            f.write(poscar_dict['atomname_list'][i] + ' ' + '{:.4f}'.format(estru[i]) + ' ' + '{:.4f}'.format(2*n1[i]) + '\n')

##            #Integrate over energy spectrum: trapezoidal integration
##            NEDOS_BelowFermi = 0
##            for j in range(0,NEDOS):
##                if energy[j] < e_fermi_mod:
##                    NEDOS_BelowFermi += 1
##            n_Mult_E = np.array([None]*NEDOS_BelowFermi)
##            for j in range(0,NEDOS):
##                if energy[j] < e_fermi_mod:
##                    n_Mult_E[j] = dos_sum[j,i]         
##            n2[i] = np.trapz(n_Mult_E, dx = e_grid)
##            ne2 += n2[i] * 2

        funcs.write_log(logfile, 'Calculated number of valence electrons=' + '{:.2f}'.format(ne1) + '\n')
        funcs.write_log(logfile, 'Actual number of valence electrons=' + str(valence_electron_num_sum) + '\n')
        funcs.write_log(logfile, 'eband=' + '{:.2f}'.format(eband) + ' eV' + '\n')

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

    return eband

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
    sample_dos_arr is the DOS array which is to be labeled, for example, Sampltdos_arr=d, or Sampltdos_arr=dxy, or Sampltdos_arr=s, etc.
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
        #Spin unpolarized calculation
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
        dos_temp_arr = dos_temp_arr[len(dos_arr[:,[0]]):,0]
        if dos_arr.shape[1] == 17:
            fm3 = dos_arr[:,[10]]
            fm2 = dos_arr[:,[11]]
            fm1 = dos_arr[:,[12]]
            f0 = dos_arr[:,[13]]
            f1 = dos_arr[:,[14]]
            f2 = dos_arr[:,[15]]
            f3 = dos_arr[:,[16]]
            f = fm3 + fm2 + fm1 + f0 + f1 + f2 + f3

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
        dos_temp_arr = dos_temp_arr[len(dos_arr[:,[0]]):,0]
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
            if 'f' in orbits:
                dos_temp_arr = np.concatenate((dos_temp_arr, f_up, f_dw),axis = 0)
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

def overlap_peak_analyzer(doscar_dir, sysname, atom_indx_list,n_shell,a0, dos_mode = None, fermi_shift_zero = True):
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
     * doscar_dir: String format. The directory which contains the DOSCAR file, abstract path can be accepted
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

    #Default values
    if dos_mode == None:
        dos_mode = defaults_dict['dos_mode']
    
    round_decimal = 2    #Rounding of the energy level
    max_neighbor_num = 50   #Maximum number of neighbors, please do not modify this variable unless you know what you are doing
    formatted_len = 12 #For example len('Ni324_-1-1-1')=12
    
    doscar_dir = os.path.abspath(doscar_dir)
    if sysname == None:
        sysname = ''

    workdir, dos_file = funcs.file_path_name(doscar_dir)
    outcar_dir = workdir + '/OUTCAR'
    doscar_dir = workdir + '/DOSCAR'
    incar_dir = workdir + '/INCAR'
    poscar_dir = workdir + '/POSCAR'
    poscar_dict = vasp_read.read_poscar(poscar_dir)
    incar_params_dict, outcar_params_dict = vasp_read.read_outcar(outcar_dir)
    LORBIT = incar_params_dict['LORBIT']
    e_fermi = outcar_params_dict['e_fermi']
    e_fermi_mod = e_fermi + outcar_params_dict['alpha+bet']

    output_dir = os.getcwd() + '/outputs'
    funcs.mkdir(output_dir)

    overlap_peak_file = output_dir + '/overlap_paek_' + str(sysname) + '.dat'
    eff_overlap_energy_file = output_dir + '/EffectiveOverlapEnergy_' + str(sysname) + '.dat'
    with open(overlap_peak_file,'w') as f, open(eff_overlap_energy_file,'w') as f1:
        if fermi_shift_zero == True:
            f.write('e_fermi = 0\n')
        else:
            f.write('e_fermi = '+str(e_fermi_mod)+'\n')
    with open(overlap_peak_file,'a') as f, open(eff_overlap_energy_file,'a') as f1:
        for i in range(0,len(atom_indx_list)):
            if isinstance(atom_indx_list[i],str) and atom_indx_list[i] != 'TDOS':
                atom_indx = convert.atomname2indx(workdir + '/POSCAR',atom_indx_list[i])
                atomname_temp = str(atom_indx_list[i])
            elif isinstance(atom_indx_list[i],str) and atom_indx_list[i] == 'TDOS':
                atom_indx = 0
            elif isinstance(atom_indx_list[i],int):
                atom_indx = atom_indx_list[i]
                atomname_temp = str(poscar_dict['atom_species_arr'][atom_indx - 1])
            f.write('################################\n')
            f.write(atomname_temp + '\n')
            f1.write('################################\n')
            f1.write(atomname_temp + '\n')            
            orbits = dos_mode[str(poscar_dict['atom_species_arr'][atom_indx - 1])]
            doscar_dict = vasp_read.read_doscar(doscar_dir, atom_indx, False)
            energy = doscar_dict['energy'] + outcar_params_dict['alpha+bet'] - e_fermi_mod*funcs.logic_retn_val(fermi_shift_zero,1,0)

            energy_peak_list, peak_val_list, orbitname_list = get_orbit_dos_peak_list(doscar_dict['dos_arr'], orbits, energy)
            eff_overlap_energy = np.array([0.000000000]*len(orbitname_list))
            #######################################
            # For the neighbors of the selected atom
            #######################################
            nn_atomname_ishell_arr, nn_atomname_ishell_arr_no_extend_atomname = nn_map(poscar_dir, a0, n_shell)
            for i_shell in range(n_shell):
                f.write('\n' + str(i_shell+1) + 'NN Orbit(' + atomname_temp + ') Orbit(X) OverlappedPeaks...' + '\n')
                f1.write('\n' + str(i_shell+1) + 'NN Orbit(' + atomname_temp + ')' + '\n')
                for i_neighbor in range(max_neighbor_num):
                    if nn_atomname_ishell_arr_no_extend_atomname[atom_indx-1,i_shell,i_neighbor] == None:
                        continue
                    i_nn_atomname = nn_atomname_ishell_arr_no_extend_atomname[atom_indx-1,i_shell,i_neighbor].strip()
                    if isinstance(i_nn_atomname,str) and i_nn_atomname != 'TDOS':
                        i_nn_atom_indx = convert.atomname2indx(workdir + '/POSCAR',i_nn_atomname)
                        iNN_atomname_temp = i_nn_atomname
                    elif isinstance(i_nn_atomname,str) and i_nn_atomname == 'TDOS':
                        i_nn_atom_indx = 0
                    elif isinstance(i_nn_atomname,int):
                        i_nn_atom_indx = i_nn_atomname
                        iNN_atomname_temp = str(poscar_dict['atom_species_arr'][atom_indx - 1])

                    iNN_orbits = dos_mode[str(poscar_dict['atom_species_arr'][i_nn_atom_indx - 1])]
                    
                    i_nn_doscar_dict = vasp_read.read_doscar(doscar_dir, i_nn_atom_indx, False)
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
                            f.write(iNN_atomname_temp + ' ' + str(orbitname_list[i_orbit]) + ' ' + str(orbitname_list[j_orbit]) +
                                    ' ' + ' '.join([str('{:.2f}'.format(i)) for i in sorted_energy_peaks3]) + '\n')
                            #Calculate effective overlap energy
                            if fermi_shift_zero == True:
                                sorted_energy_peaks3 = sorted_energy_peaks3 + e_fermi_mod*funcs.logic_retn_val(fermi_shift_zero,1,0)
                            else:
                                pass
                            for iPeak in range(0,len(sorted_energy_peaks3)):
                                if sorted_energy_peaks3[iPeak] <= e_fermi_mod:
##                                    eff_overlap_energy[i_orbit] += sorted_energy_peaks3[iPeak] * np.min([Sortedoverlap_peak_val_list1[iPeak],Sortedoverlap_peak_val_list2[iPeak]])
                                    eff_overlap_energy[i_orbit] += sorted_energy_peaks3[iPeak] * Sortedoverlap_peak_val_list1[iPeak]
            for i_orbit in range(0,len(orbitname_list)):
                f1.write(str(orbitname_list[i_orbit]) + ' = ')
                f1.write(str(eff_overlap_energy[i_orbit]) + '\n')
            f1.write('sum = ' + str(np.sum(eff_overlap_energy)) + '\n')
    return 0


            
