# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use("Agg")

def plot_dos(atom_doscar_file_path_list, atom_sysname_list, atom_indx_list, atom_palette_list, atom_subplot_arg_list,
             subplot_arg_list, subplot_xlo_list, subplot_xhi_list, subplot_ylo_list, subplot_yhi_list,
             subplot_xtick_list, subplot_ytick_list, subplot_xlabel_list, subplot_ylabel_list, subplot_share_xy_list = [False, False] , mainplot_axis_label_list = [False, False],
             dos_mode = None, fermi_shift_zero = True, peak_analyzer = False, fig_format = 'png', fig_size = [15,10], fig_dpi = 600):
    '''
    - Descriptions
     * Plot PDOS, LDOS, TDOS, now only available for LORBIT = 11.
     * There are three types of input arguments: atom related input arguments, subplot related input arguments, and others

    - Atom related Args
     * atom_doscar_file_path_list: list format. Conatains DOSCAR files for each atom. The directory of DOSCAR files can either be full path or relative path
     * atom_sysname_list: system name for each atom, it corresponds to the atoms in the atom_doscar_file_path_list. This is for the purpose of labeling the DOS curves in the legend.
               If sysnameList = None, then the label of system name will not shown in the legend
               For example, sysnameList = ['System1','System1','System2']
     * atom_indx_list: list format. Atom index list, it corresponding to the atoms in  atom_doscar_file_path_list. If it is integer type then it denotes the atom index, if it is string type then it denotes the atom name
               atom_indx_list = [1,2,45] denotes the 1st, 2nd, and the 45th atoms in the POSCAR
               atom_indx_list = ['Ni1','Al3','Re3'] denotes Ni1, Al3, and Re3 in the POSCAR
               atom_indx_list = ['TDOS'] and atom_indx_list = [0] donotes the total dos
     * atom_palette_list: list format. Color for DOS curves of each atom.
     * atom_subplot_arg_list: list format. Defines the DOS curves of the atom are in which subplot. For example, atom_subplot_arg_list = [221, 222] denotes that the DOS curves of the first and the second atoms are in the subplot(221) and subplot(222) subplots, respectively. 

    - Subplot related Args
     * subplot_arg_list: list format. The subplot argument list, for example subplot_arg_list=[221,222] corresponds to subplot(221) and subplot(222)
     * subplot_xlo_list: list format. Low boundary of the x axis for each subplots. If None value is given, the low boundary of x axis in the data set will be chosen.
     * subplot_xhi_list: list format. High boundary of the x axis for each subplots. If None value is given, the high boundary of x axis in the data set will be chosen.
     * subplot_ylo_list: list format. Low boundary of the y axis for each subplots. If None value is given, the low boundary of y axis in the data set will be chosen.
     * subplot_yhi_list: list format. High boundary of the y axis for each subplots. If None value is given, the high boundary of y axis in the data set will be chosen.
     * subplot_xtick_list: list of logical values. If the list element is True (False), then the tick of the x axis will be shown (removed). 
     * subplot_ytick_list: list of logical values. If the list element is True (False), then the tick of the x axis will be shown (removed).
     * subplot_xlabel_list: list of logical values. Defines whether the x-label of each subplots are shown, it won't work for subplot=(111) figure.
     * subplot_xlabel_list: list of logical values. Defines whether the y-label of each subplots are shown, it won't work for subplot=(111) figure.
     * subplot_share_xy_list: list of logical values of length two. Defines whether the x or y axis will be shared or not. [False, False] denotes both x and y axes will not be shared. 
    '''
    import os
    import sys
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from .. import default_params
    from .. import periodic_table

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)
    periodic_table_dict = periodic_table.periodic_tab()
    #Default values
    if dos_mode == None:
        dos_mode = periodic_table_dict['dos_mode']

    #Possible further modifications: check the consistency of lengths of lists: atom_doscar_file_path_list and atom_indx_list; subplot_arg_list and SubplotXYLabelLogicList
    ###############################################################################################################
    #User defined region, change the params when necessary. Usually we don't change these parameters unless needed.
    ###############################################################################################################
    if atom_palette_list == None or atom_palette_list == 'None':
        atom_palette_list = ['black','magenta','cyan','darkviolet','red','gray','darkorange','darkcyan','palegreen',
                   'goldenrod','lightpink','tan','mediumpurple','steelblue','tomato','mediumturquoise',
                   'mediumslateblue','brown','darkseagreen','wheat','seagreen','maroon','deeppink','chocolate',
                   'bisque','plum']
    dos_lty = ['-','--','-.',':','-','--','-.',':']    #Line style for DOS orbital
    line_width = 1.0            # Recommended value 0.5 -- 2.0 (from thin to fat)
    upper_spine_scale = 1.02     #Scale the coordinate of the upper axis spine
    lower_spine_scale = 1.02     #Scale the coordinate of the lower axis spine
    #scale the DOS's: for purpose of plotting arbitrary unit DOS, default value = 1
    orbit_scaling_dict = {'s':1, 'p':1, 'd':1,'f':1,'LDOS':1}
    font_size = 18
    ##################################
    #Codes maintained by the developer 
    ##################################

    def label_peaks(energy, sample_dos_arr, subplot_indx, subplot_dict):
        '''
        Description:
        Roughly label the peaks
        sample_dos_arr is the DOS array which is to be labeled, for example, Sampltdos_arr=d, or Sampltdos_arr=dxy, or Sampltdos_arr=s, etc.
        '''
        import numpy as np
        import matplotlib.pyplot as plt
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
        for ipeak in range(0,len(peaks)):
            if energy[peaks][ipeak] >= subplot_dict['xlo'][subplot_indx] and energy[peaks][ipeak] <= subplot_dict['xhi'][subplot_indx]:
                plt.text(energy[peaks][ipeak], sample_dos_arr_1d[peaks][ipeak], str(energy[peaks][ipeak]),color='red',size=3)
        plt.plot(energy[peaks], sample_dos_arr_1d[peaks], "x",color='red')
        return energy[peaks]
    for i in range(0,len(atom_doscar_file_path_list)):
        atom_doscar_file_path_list[i] = os.path.abspath(atom_doscar_file_path_list[i])
    if atom_sysname_list == None:
        atom_sysname_list = [''] * len(atom_doscar_file_path_list)
    else:
        for i in range(0,len(atom_sysname_list)):
            if atom_sysname_list[i] == None or atom_sysname_list[i] == 'None' or atom_sysname_list[i] == 'none' or atom_sysname_list[i] == '':
                atom_sysname_list[i] = ''
            else:
                atom_sysname_list[i] = str(atom_sysname_list[i]) + ','
    subplot_arg_unique_list = list(set(atom_subplot_arg_list))
    #Initialize Subplot dictionary
    subplot_dict = {}
    subplot_dict['Arg'] = subplot_arg_unique_list
    subplot_dict['indx'] = list(range(len(subplot_arg_unique_list)))
    subplot_dict['ShareXY'] = [False,False]
    subplot_dict['xlo'] = [999999.9999] * len(subplot_arg_unique_list)
    subplot_dict['xhi'] = [-999999.9999] * len(subplot_arg_unique_list)
    subplot_dict['ylo'] = [999999.9999] * len(subplot_arg_unique_list)
    subplot_dict['yhi'] = [-999999.9999] * len(subplot_arg_unique_list)
    subplot_dict['xtick'] = [False] * len(subplot_arg_unique_list)
    subplot_dict['ytick'] = [False] * len(subplot_arg_unique_list)
    subplot_dict['xlabel'] = [None] * len(subplot_arg_unique_list)
    subplot_dict['ylabel'] = [None] * len(subplot_arg_unique_list)
    #Extract data from the user input into the subplot_dict
    subplot_dict['ShareXY'] = subplot_share_xy_list
    for i_subplot in range(0,len(subplot_arg_list)):
        if subplot_arg_list[i_subplot] in subplot_dict['Arg']:
            subplot_dict_indx = subplot_dict['Arg'].index(subplot_arg_list[i_subplot])
            if subplot_xlo_list[i_subplot] != None:
                subplot_dict['xlo'][subplot_dict_indx] = subplot_xlo_list[i_subplot]
            if subplot_xhi_list[i_subplot] != None:
                subplot_dict['xhi'][subplot_dict_indx] = subplot_xhi_list[i_subplot]
            if subplot_ylo_list[i_subplot] != None:
                subplot_dict['ylo'][subplot_dict_indx] = subplot_ylo_list[i_subplot]
            if subplot_yhi_list[i_subplot] != None:
                subplot_dict['yhi'][subplot_dict_indx] = subplot_yhi_list[i_subplot]
            subplot_dict['xtick'][subplot_dict_indx] = subplot_xtick_list[i_subplot]
            subplot_dict['ytick'][subplot_dict_indx] = subplot_ytick_list[i_subplot]
            subplot_dict['xlabel'][subplot_dict_indx] = subplot_xlabel_list[i_subplot]
            subplot_dict['ylabel'][subplot_dict_indx] = subplot_ylabel_list[i_subplot]
    subplot_arg_unique_count_list = []
    counter_list = []
    for i in range(len(subplot_arg_unique_list)):
        subplot_arg_unique_count_list.append(atom_subplot_arg_list.count(subplot_arg_unique_list[i]))
        counter_list.append(0)
    sysnameList = []
    dos_atomname_list = []
    e_fermi_list = []
    plt.rcParams['figure.figsize'] = (fig_size[0],fig_size[1])
    if len(atom_doscar_file_path_list) == 1:
        doscar_file_path = atom_doscar_file_path_list[0]
        workdir, dos_file = funcs.file_path_name(doscar_file_path)
        if isinstance(atom_indx_list[0],str) and atom_indx_list[0] != 'TDOS':
            atom_indx = convert.atomname2indx(workdir + '/POSCAR',atom_indx_list[0])
        elif isinstance(atom_indx_list[0],str) and atom_indx_list[0] == 'TDOS':
            atom_indx = 0
        elif isinstance(atom_indx_list[0],int):
            atom_indx = atom_indx_list[0]
        fig_dos = plt.figure('Fig'+dos_file)
    else:
        fig_dos = plt.figure('Fig_DOS_Cmpr')
    ax0 = fig_dos.add_subplot(111)
    #Below is the axis setting for the main figure, and the axis settings only work for main figure with multiple subplots(number of subplots >=2):
    #set interspacing between subplots
    if subplot_dict['ShareXY'][0] == True:
        hspace_val = 0
    else:
        hspace_val = None
    if subplot_dict['ShareXY'][1] == True:
        wspace_val = 0
    else:
        wspace_val = None
    plt.subplots_adjust(bottom =0.070, left = 0.090, top = 0.99, right = 0.98, wspace = wspace_val, hspace = hspace_val)
    subplot_last_arg = atom_subplot_arg_list[len(atom_subplot_arg_list)-1]
    mult_plot_logic = isinstance(subplot_last_arg,int) and subplot_last_arg != 111 or isinstance(subplot_last_arg,list) and subplot_last_arg != [1,1,1]
    xlabel_str = 'Energy (eV)'
    if '0' in atom_indx_list or 'TDOS' in atom_indx_list:
        ylabel_str = 'TDOS (arb. units)'
    else:
        if ['LDOS'] not in dos_mode.values():
            ylabel_str = 'PDOS (arb. units)'
        elif ['LDOS'] in dos_mode.values():
            ylabel_str = 'LDOS (arb. units)'
    if mult_plot_logic:
        ax0.set(xticks = [])
        ax0.set(yticks = [])
        ax0.spines['top'].set_visible(False)
        ax0.spines['bottom'].set_visible(False)
        ax0.spines['left'].set_visible(False)
        ax0.spines['right'].set_visible(False)
        ax0.spines['bottom'].set_position(('outward',20))
        ax0.spines['left'].set_position(('outward',30))
        plt.xlabel(funcs.logic_retn_val(mainplot_axis_label_list[0],xlabel_str,''), fontsize = font_size)
        plt.ylabel(funcs.logic_retn_val(mainplot_axis_label_list[1],ylabel_str,''), fontsize = font_size, labelpad=font_size*1.3)
    def plot_lines(orbits, doscar_dict, energy, atom_palette_list, dos_lty, line_width, peak_analyzer, subplot_indx, subplot_dict, atom_sysname_list, dos_file_indx, i_elmt_name, orbit_scaling_dict):
        import matplotlib.pyplot as plt
        orbitname_list = ['f','d','dxy','dyz','dz2','dxz','dx2','p','py','pz','px','s','LDOS']
        orbit_label_list = ['-$f$',
                          '-$d$','-$d_{xy}$','-$d_{yz}$','-$d_{z2}$','-$d_{xz}$','-$d_{x^{2}-y^{2}}$',
                          '-$p$','-$p_{y}$','-$p_{z}$','-$p_{x}$',
                          '-$s$',
                          '']
        if 'LDOS' in orbits and len(orbits) > 1:
            # when LDOS and PDOS are plotted in one figure
            orbit_label_list[-1] = '-LDOS'
        dos_lty_indx = 0
        for orbitname_indx in range(0,len(orbitname_list)):
            if orbitname_list[orbitname_indx] in orbits:
                orbit_scaling_factor = orbit_scaling_dict[orbitname_list[orbitname_indx]]
                if doscar_dict['num_col'] == 5 or doscar_dict['num_col'] == 19 or doscar_dict['num_col'] == 33:
                    plt.plot(energy,doscar_dict[orbitname_list[orbitname_indx]+'_up']*orbit_scaling_factor,color=atom_palette_list[dos_file_indx],linestyle = dos_lty[dos_lty_indx],
                             marker='',linewidth = line_width, label=atom_sysname_list[dos_file_indx]+i_elmt_name + orbit_label_list[orbitname_indx])
                    plt.plot(energy,doscar_dict[orbitname_list[orbitname_indx]+'_dw']*orbit_scaling_factor,color=atom_palette_list[dos_file_indx],linestyle = dos_lty[dos_lty_indx],
                             marker='',linewidth = line_width)
                    dos_lty_indx += 1                
                    if peak_analyzer == True:
                        label_peaks(energy, doscar_dict[orbitname_list[orbitname_indx]+'_up']*orbit_scaling_factor, subplot_indx, subplot_dict['xlo'][subplot_indx], subplot_dict['xhi'][subplot_indx])
                        label_peaks(energy, doscar_dict[orbitname_list[orbitname_indx]+'_dw']*orbit_scaling_factor, subplot_indx, subplot_dict['xlo'][subplot_indx], subplot_dict['xhi'][subplot_indx])
                elif doscar_dict['num_col'] == 3 or doscar_dict['num_col'] == 10 or doscar_dict['num_col'] == 17:
                    plt.plot(energy,doscar_dict[orbitname_list[orbitname_indx]]*orbit_scaling_factor,color=atom_palette_list[dos_file_indx],linestyle = dos_lty[dos_lty_indx],
                             marker='',linewidth = line_width, label=atom_sysname_list[dos_file_indx]+i_elmt_name + orbit_label_list[orbitname_indx])
                    dos_lty_indx += 1                
                    if peak_analyzer == True:
                        label_peaks(energy, doscar_dict[orbitname_list[orbitname_indx]]*orbit_scaling_factor, subplot_indx, subplot_dict)                        
    #Below plot DOS in each subplot
    for dos_file_indx in range(0,len(atom_doscar_file_path_list)):
        doscar_file_path = atom_doscar_file_path_list[dos_file_indx]
        workdir, dos_file = funcs.file_path_name(doscar_file_path)
        if isinstance(atom_indx_list[dos_file_indx],str) and atom_indx_list[dos_file_indx] != 'TDOS':
            atom_indx = convert.atomname2indx(workdir + '/POSCAR',atom_indx_list[dos_file_indx])
        elif isinstance(atom_indx_list[dos_file_indx],str) and atom_indx_list[dos_file_indx] == 'TDOS':
            atom_indx = 0
        elif isinstance(atom_indx_list[dos_file_indx],int):
            atom_indx = atom_indx_list[dos_file_indx]
        doscar_dict = vasp_read.read_doscar(doscar_file_path, atom_indx, False)
        outcar_file_path = workdir + '/OUTCAR'
        doscar_file_path = workdir + '/DOSCAR'
        incar_file_path = workdir + '/INCAR'
        poscar_file_path = workdir + '/POSCAR'
        poscar_dict = vasp_read.read_poscar(poscar_file_path)
        incar_params_dict, outcar_params_dict = vasp_read.read_outcar(outcar_file_path)
        LORBIT = incar_params_dict['LORBIT']
        e_fermi = outcar_params_dict['e_fermi']
        e_fermi_mod = e_fermi + outcar_params_dict['alpha+bet']
        subplot_indx = subplot_arg_unique_list.index(atom_subplot_arg_list[dos_file_indx])
        if atom_indx == 0:
            #TDOS
            '''If the subplot's parameters coincide with the previous subplot's parameters, a MatplotlibDeprecationWarning will occur. So we take the following action:
               If new added subplot's parameters is the same as the previous subplot's parameters, we don't add the new subplot, we use the previous subplot instead
               If new added subplot's parameters is different from the previous subplot's parameters, we add the new subplot.'''
            if dos_file_indx == 0:
                previous_subplot_arg = 111
            else:
                previous_subplot_arg = atom_subplot_arg_list[dos_file_indx - 1]
            if atom_subplot_arg_list[dos_file_indx] == previous_subplot_arg:
                pass
            else:
                fig_dos.add_subplot(atom_subplot_arg_list[dos_file_indx])       
            
            counter_list[subplot_indx] += 1
            sysnameList.append(''.join(str(poscar_dict['elmt_species_arr'][i])+str(poscar_dict['elmt_num_arr'][i]) for i in range(len(poscar_dict['elmt_species_arr'])-1,-1,-1)))
            dos_atomname_list.append('TDOS')
            e_fermi_list.append(e_fermi_mod)
            energy = doscar_dict['energy'] + outcar_params_dict['alpha+bet'] - e_fermi_mod*funcs.logic_retn_val(fermi_shift_zero,1,0)
            if subplot_xlo_list[subplot_arg_list.index(subplot_dict['Arg'][subplot_indx])] == None or subplot_xhi_list[subplot_arg_list.index(subplot_dict['Arg'][subplot_indx])] == None:
                subplot_dict['xlo'][subplot_indx] = min(subplot_dict['xlo'][subplot_indx],energy.min())
                subplot_dict['xhi'][subplot_indx] = max(subplot_dict['xhi'][subplot_indx],energy.max())
            if subplot_ylo_list[subplot_arg_list.index(subplot_dict['Arg'][subplot_indx])] == None or subplot_yhi_list[subplot_arg_list.index(subplot_dict['Arg'][subplot_indx])] == None:
                subplot_dict['ylo'][subplot_indx] = subplot_dict['ylo'][subplot_indx] * lower_spine_scale
                subplot_dict['yhi'][subplot_indx] = subplot_dict['yhi'][subplot_indx] * upper_spine_scale
            if doscar_dict['num_col'] == 3:
                plt.plot(energy,doscar_dict['TDOS'],color = atom_palette_list[dos_file_indx],marker='',linewidth = line_width,
                         label = atom_sysname_list[dos_file_indx]+'TDOS')
                if peak_analyzer == True:  
                    label_peaks(energy, doscar_dict['TDOS'], subplot_indx, subplot_dict)
            elif doscar_dict['num_col'] == 5:
                plt.plot(energy,doscar_dict['TDOS_up'],color = atom_palette_list[dos_file_indx],marker='',linewidth = line_width,
                         label = atom_sysname_list[dos_file_indx]+'TDOS')
                plt.plot(energy,doscar_dict['TDOS_dw'],color = atom_palette_list[dos_file_indx],marker='',linewidth = line_width)
                if peak_analyzer == True:  
                    label_peaks(doscar_dict['energy'], doscar_dict['TDOS_up'], subplot_indx, subplot_dict)
                    label_peaks(doscar_dict['energy'], doscar_dict['TDOS_dw'], subplot_indx, subplot_dict)                  
        elif atom_indx != 0:
            #PDOS and LDOS
            if dos_file_indx == 0:
                previous_subplot_arg = 111
            else:
                previous_subplot_arg = atom_subplot_arg_list[dos_file_indx - 1]
            if atom_subplot_arg_list[dos_file_indx] == previous_subplot_arg:
                pass
            else:
                fig_dos.add_subplot(atom_subplot_arg_list[dos_file_indx])     
            counter_list[subplot_indx] += 1
            i_elmt_name = poscar_dict['atomname_list'][atom_indx - 1]
            sysnameList.append(''.join(str(poscar_dict['elmt_species_arr'][i])+str(poscar_dict['elmt_num_arr'][i]) for i in range(len(poscar_dict['elmt_species_arr'])-1,-1,-1)))
            dos_atomname_list.append(i_elmt_name)
            e_fermi_list.append(e_fermi_mod)
            energy = doscar_dict['energy']  + outcar_params_dict['alpha+bet'] - e_fermi_mod*funcs.logic_retn_val(fermi_shift_zero,1,0)
            orbits = dos_mode[str(poscar_dict['atom_species_arr'][atom_indx - 1])]
            if len(orbits) == 0:
                print('## Error: The value of the dos_mode is None. The dos_mode must be defined!')                
            orbitname_list = ['f','d','dxy','dyz','dz2','dxz','dx2','p','py','pz','px','s','LDOS']
            #Subplot x axis limit
            if subplot_xlo_list[subplot_arg_list.index(subplot_dict['Arg'][subplot_indx])] == None or subplot_xhi_list[subplot_arg_list.index(subplot_dict['Arg'][subplot_indx])] == None:
                subplot_dict['xlo'][subplot_indx] = min(subplot_dict['xlo'][subplot_indx],energy.min())
                subplot_dict['xhi'][subplot_indx] = max(subplot_dict['xhi'][subplot_indx],energy.max())
            #Subplot y axis limit
            if doscar_dict['num_col'] == 10 or doscar_dict['num_col'] == 17:
                #Spin unpolarized calculation
                for orbit_indx in range(0,len(orbitname_list)):
                    if orbitname_list[orbit_indx] in orbits:
                        subplot_dict['ylo'][subplot_indx] = min(subplot_dict['ylo'][subplot_indx],min(doscar_dict[orbitname_list[orbit_indx]])[0])
                        subplot_dict['yhi'][subplot_indx] = max(subplot_dict['yhi'][subplot_indx],max(doscar_dict[orbitname_list[orbit_indx]])[0])
            elif doscar_dict['num_col'] == 19 or doscar_dict['num_col'] == 33:
                #Spin polarized calculation
                for orbit_indx in range(0,len(orbitname_list)):
                    if orbitname_list[orbit_indx] in orbits:
                        subplot_dict['ylo'][subplot_indx] = min(subplot_dict['ylo'][subplot_indx],min(doscar_dict[orbitname_list[orbit_indx]+'_up'])[0],min(doscar_dict[orbitname_list[orbit_indx]+'_dw'])[0])
                        subplot_dict['yhi'][subplot_indx] = max(subplot_dict['yhi'][subplot_indx],max(doscar_dict[orbitname_list[orbit_indx]+'_up'])[0],max(doscar_dict[orbitname_list[orbit_indx]+'_dw'])[0])
            if subplot_ylo_list[subplot_arg_list.index(subplot_dict['Arg'][subplot_indx])] == None or subplot_yhi_list[subplot_arg_list.index(subplot_dict['Arg'][subplot_indx])] == None:
                subplot_dict['ylo'][subplot_indx] = subplot_dict['ylo'][subplot_indx] * lower_spine_scale
                subplot_dict['yhi'][subplot_indx] = subplot_dict['yhi'][subplot_indx] * upper_spine_scale
            plot_lines(orbits, doscar_dict, energy, atom_palette_list, dos_lty, line_width, peak_analyzer, subplot_indx, subplot_dict, atom_sysname_list, dos_file_indx, i_elmt_name, orbit_scaling_dict)
        #Plot Fermi energy indicator line:
        if len(atom_doscar_file_path_list) == 1:
            #Figure only contain one subplot, single DOS is plotted or multiple DOS files are plotted onto one subplot        
            line1 = [(e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['ylo'][subplot_indx]), (e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['yhi'][subplot_indx])]
            (line1_xs, line1_ys) = zip(*line1)
            plt.plot(line1_xs, line1_ys, linestyle = '--',linewidth = 0.5, color = 'black')
        elif len(atom_doscar_file_path_list) != 1 and fermi_shift_zero == True and counter_list[subplot_indx] == subplot_arg_unique_count_list[subplot_indx]:
            #Figure contain multiple subplots, Fermi energy is shifted to zero. Each subplot only plot the Fermi energy level once.
            line1 = [(e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['ylo'][subplot_indx]), (e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['yhi'][subplot_indx])]
            (line1_xs, line1_ys) = zip(*line1)
            plt.plot(line1_xs, line1_ys, linestyle = '--',linewidth = 0.5, color = 'black')
        elif len(atom_doscar_file_path_list) != 1 and fermi_shift_zero == False:
            #Figure contain multiple subplots, Fermi energy is not shifted to zero.
            e_fermi = e_fermi_list[dos_file_indx]
            line1 = [(e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['ylo'][subplot_indx]),(e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['yhi'][subplot_indx])]
            (line1_xs, line1_ys) = zip(*line1)
            plt.plot(line1_xs, line1_ys, linestyle = '--',linewidth = 0.5,color = atom_palette_list[dos_file_indx])
        #axes, ticks, labels and legends
        ax = plt.gca()
        ax.set(title='')
        plt.xlim(subplot_dict['xlo'][subplot_indx], subplot_dict['xhi'][subplot_indx])
        plt.ylim(subplot_dict['ylo'][subplot_indx], subplot_dict['yhi'][subplot_indx])
        ax.margins(x=0.0, y=0.0)
##        ax.xaxis.set_major_formatter(FormatStrFormatter( '%.1f'))
##        ax.yaxis.set_major_formatter(FormatStrFormatter( '%.1f'))
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        #if the hspace and wspace ==0, then remove the endpoints of the axis label and axis ticks.
        if subplot_dict['ShareXY'][0] == True:
            yticks = ax.yaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)
##            yticks[-1].label1.set_visible(False)
        if subplot_dict['ShareXY'][1] == True:
            xticks = ax.xaxis.get_major_ticks()
            xticks[0].label1.set_visible(False)
##            xticks[-1].label1.set_visible(False)
        if subplot_dict['xtick'][subplot_indx] == False:
            plt.xticks([])
        if subplot_dict['ytick'][subplot_indx] == False:
            plt.yticks([])
        xlabel_str = 'Energy (eV)'
        if atom_indx == 0:
            ylabel_str = 'TDOS (arb. units)'
        else:
            if 'LDOS' not in orbits:
                ylabel_str = 'PDOS (arb. units)'
            elif 'LDOS' in orbits and len(orbits) == 1:
                ylabel_str = 'LDOS (arb. units)'
            elif 'LDOS' in orbits and len(orbits) > 1:
                ylabel_str = 'DOS (arb. units)'
            elif len(orbits) == 0:
                print('## Error: The value of the dos_mode is None. The dos_mode must be defined!')
        subplot_xlabel_on_off_str = funcs.logic_retn_val(subplot_dict['xlabel'][subplot_indx], xlabel_str,'')
        subplot_ylabel_on_off_str = funcs.logic_retn_val(subplot_dict['ylabel'][subplot_indx], ylabel_str,'')
        ax.set_xlabel(funcs.logic_retn_val(mult_plot_logic, subplot_xlabel_on_off_str, xlabel_str),fontsize=font_size)
        ax.set_ylabel(funcs.logic_retn_val(mult_plot_logic, subplot_ylabel_on_off_str, ylabel_str),fontsize=font_size)
        ax.legend(loc='best',fontsize=font_size)
    #Export Figures
    #sysindx_list denotes the systems as integers. Identical doscar_file_path have the same integer value.
    #As the items of atom_doscar_file_path_list goes from the first one to the last one, The corresponding position of sysindx_list is assigned a unique interger value(from 1 to N).
    sysindx_list = []
    Uniquesysindx_list = []
    unique_atom_doscar_file_path_list = []
    Uniquee_fermi_list = []
    sysindx = 1
    for idoscar_file_path in atom_doscar_file_path_list:
        if idoscar_file_path in unique_atom_doscar_file_path_list:
            sysindx_list.append(Uniquesysindx_list[unique_atom_doscar_file_path_list.index(idoscar_file_path)])
        else:
            Uniquesysindx_list.append(sysindx)
            unique_atom_doscar_file_path_list.append(idoscar_file_path)
            sysindx_list.append(sysindx)
            Uniquee_fermi_list.append(e_fermi_list[sysindx - 1])
            sysindx += 1
    formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time.time()))
    fig_file = output_dir + '\Fig_' + str(formatted_time) + '_sys' + ''.join([str(x) for x in sysindx_list]) + '_' + ''.join(dos_atomname_list) + '.' + fig_format
    plt.savefig(fig_file,dpi = fig_dpi)
    plt.close()
    #Write Figure information into logfile
    sys_list = []
    for iSys in range(0,len(atom_indx_list)):
        sys_list.append('sys' + str(sysindx_list[iSys]))
    funcs.write_log(logfile, '## The system directorie(s) are given below:')
    for i_unique_sys in range(0,len(unique_atom_doscar_file_path_list)):
        funcs.write_log(logfile,
                        'sys' + str(Uniquesysindx_list[i_unique_sys]) + ' = r\'' + str(unique_atom_doscar_file_path_list[i_unique_sys]) + '\'' + '\n'
                        'e_fermi_mod = ' + str(e_fermi_list[i_unique_sys]))
    funcs.write_log(logfile,
                    'dos_mode = ' + str(dos_mode) + '\n' +
                   'vasp_plot.plot_dos(' + '\n' +
                   '    atom_doscar_file_path_list=' + ('[' + ','.join(sys_list) + ']') + ',\n' +
                   '    atom_sysname_list=' + str(atom_sysname_list) + ',\n' +
                   '    atom_indx_list=' + str(atom_indx_list) + ',\n' +
                   '    atom_palette_list=' + str(atom_palette_list) + ',\n' +
                   '    atom_subplot_arg_list=' + str(atom_subplot_arg_list) + ',\n' +
                   '    subplot_arg_list=' + str(subplot_arg_list) + ',\n' +
                   '    subplot_xlo_list=' + str(subplot_xlo_list) + ',\n' +
                   '    subplot_xhi_list=' + str(subplot_xhi_list) + ',\n' +
                   '    subplot_ylo_list=' + str(subplot_ylo_list) + ',\n' +
                   '    subplot_yhi_list=' + str(subplot_yhi_list) + ',\n' +
                   '    subplot_xtick_list=' + str(subplot_xtick_list) + ',\n' +             
                   '    subplot_ytick_list=' + str(subplot_ytick_list) + ',\n' +
                   '    subplot_xlabel_list=' + str(subplot_xlabel_list) + ',\n' +              
                   '    subplot_ylabel_list=' + str(subplot_ylabel_list) + ',\n' +
                   '    subplot_share_xy_list=' + str(subplot_share_xy_list) + ',\n' +
                   '    mainplot_axis_label_list=' + str(mainplot_axis_label_list) + ',\n' +
                   '    dos_mode=dos_mode' + ',\n' +
                   '    fermi_shift_zero=' + str(fermi_shift_zero) + ',\n' +
                   '    peak_analyzer=' + str(peak_analyzer) + ',\n' +
                   '    fig_format=' + '\'' + str(fig_format) + '\'' + ',\n' +
                   '    fig_size=' + str(fig_size) + ',\n' +
                   '    fig_dpi=' + str(fig_dpi) + ')\n' +
                   'fig_file= ' + 'r\'' + fig_file + '\'' + '\n' +
                   '################################################\n')
    return subplot_dict


def plot_poscar(poscar_file_path, euler_angle_type = 'zyx', phi = -3, theta = 5, psi = 0, elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
                plot_cell_basis_vector_label = False, plot_atom_label = True, fig_format = 'png', fig_dpi = 100,
                draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical'):
    '''
    - Descriptions
     * Plot POSCAR model. Euler angles are used to rotate the view of the model.
     * Viewer direction is in x direction. The original orientation: x direction is perpendicular to the paper, z direction is in the paper and point to upper direction
     * Reference for Eulerian angles: Herbert Goldstein and Charles P. Poole Jr. and John L. Safko,Classical Mechanics (3rd Edition),2001.

    - Args:
     * poscar_file_path: String format. Directory of the POSCAR file which you want to plot
     * euler_angle_type: string of length 3. It specify the type of rotations based on Eulerian angles. Choices are 'zyz', 'zxz', 'zyx', etc.. Usually the 'zyz' type is used.

                'zyz' : proper Euler angle, y-convention. Performe consecutive rotations at axes counter-clockwisely. z-y-z rotation.
                        First rotate the z axes of atoms by an angle phi, then rotate the intermidiate y axis of atoms by an angle theta, finally rotate the final z axis of atoms by an angle psi

                'zxz' : proper Euler angle, x-convention. Performe consecutive rotations at axes counter-clockwisely. z-x-z rotation.
                        First rotate the z axes of atoms by an angle phi, then rotate the intermidiate x axis of atoms by an angle theta, finally rotate the final z axis of atoms by an angle psi

                'zyx' : Tait-Bryan angles. z-y-x rotation. Performe consecutive rotations at axes counter-clockwisely. z-y-x rotation.
                        First rotate the z axes of atoms by an angle phi, then rotate the intermidiate y axis of atoms by an angle theta, finally rotate the final x axis of atoms by an angle psi

     * phi, theta, psi: float formats. The first, second, and third rotation Eulerian angles, units in degrees.
     * draw_mirror_atom: Logical value. Whether to plot the mirror atoms at the periodic boundary
     * box_on: Logical value. Whether to plot the box or not
     * axis_indicator: Logical value. Whethether to plot the axis indicator
     * plot_cell_basis_vector_label: Logical value. Whether to plot the cell basis vector labels( i.e., to label the three basis vectors of the cell as a, b, and c)
     * plot_atom_label: Logical value. If true, then plot the atom name of each atom.
     * fig_format: String format. Figformat is a string that defines output figure format. Supported fig_format: 'png', 'eps', 'pdf', 'tif', 'tiff', 'jpg', 'jpeg', 'svg', 'svgz', 'pgf', 'ps', 'raw', 'rgba'     * fig_dpi: float format. The DPI for non-vector graphics.
    '''
    import os
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn.neighbors import KDTree
    from .. import funcs
    from .. import periodic_table
    from . import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)
    ##########################
    #User defined params
    ##########################
    #Don't recommend to modify these parameters unless needed
    #If axis_indicator is TURE, please specify IndicatorLen, which controls the magnitude of the axis indicator
    #indicator_len_scale donotes the axis indicator reference length, the value of reference length equals to indicator_len_scale*max(projected BoxSideLength in three direction after rotation) and the reference length gives the maximum length for the axis indicator in three directions
    #Recommended value is indicator_len_scale=0.08;
    indicator_len_scale = 0.08
    pos_precision = 0.5  #Unit in Angstrom
    #Define reference width of the TikZ picture, units in cm, a reference value is around ref_width=8.3 for double column paper and ref_width=17.6 for one column paper.
    ref_width=8.3
    radius_scale_fac = 3.0  # Default is radius_scale_fac = 3.0. This is also the recommand value.
    model_basis_elmt = 'Ni' # Define basis element. The relative atomic radii are relative to the atomic radius of this basis element.
    label_size = 14
    elmt_color_tikz = ['gray','magenta','blue','red','green','yellow']
    periodic_table_dict = periodic_table.periodic_tab()
    if elmt_color == None:
        pass
    else:
        periodic_table_dict['elmt_color'].update(elmt_color)

    ###########################################################################
    #Developer maintained codes.
    #Please don't modify the following codes unless you know what you are doing
    ###########################################################################
    view_axis = 'x' #view_axis: Decide viewing from which direction: 'x' or 'y' or 'z'. Default is view from x axis
    poscar_file_path = os.path.abspath(poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)    
    def disp_and_rotation(poscar_file_path, phi, theta, psi):
        '''
        Description:
            Rotation center is at automatically set as its geometrical center.
            Displace the origin of the model to the center rot_center_arr and then rotate the model according to Eulerian angle
        Args:
            poscar_file_path: POSCAR of the model
            rot_center_arr: rotation center of the model, units in Angstroms
            phi, theta, psi: Eulerian angles, units in degrees
        Return:
            pos_shift_arr, t, t_pos_shift_arr, d_pos_shift_arr, t_vertex, t_vertex_min, t_vertex_max, Lo, Hi 
        '''
        import numpy as np
        from .. import funcs
        angle_unit = 'deg' # define the unit of the Eulerian angles. 'deg' or 'rad'
        poscar_dict = vasp_read.read_poscar(poscar_file_path)
        # Treat each point as a vector. Think of the problem use the knowledge of vector analysis 
        origin_arr = np.array([0,0,0]) # The default origin is (0,0,0). It is a vector
        diagonal_arr = poscar_dict['l_arr'][0,:] + poscar_dict['l_arr'][1,:] + poscar_dict['l_arr'][2,:] # diagonal_arr is the furthest vertex from the origin, i.e. the diagonal relative to the origin. it is a vector    
        rot_center_arr = (poscar_dict['l_arr'][0,:] + poscar_dict['l_arr'][1,:] + poscar_dict['l_arr'][2,:]) / 2 #rot_center_arr is the defined rotation center. Here it is taken as the geometrical center. When tha model is rotated, the rotation center is taken as the origin
        shifted_origin_arr = origin_arr - rot_center_arr  # shifted_origin_arr is the shifted origin of the model. It is a vector.
        shifted_diagonal_arr = diagonal_arr - rot_center_arr # shifted_diagonal_arr is the shifted furthest vertex point from the origin, it is a vector
        #Define 8 vertices of the simulation box. The vertices are taken as the shifted vertices relative to the rotation center(i.e., the new origin)
        # If the unit cell vector is a, b, and c. Then the unshifted eight vertices are:
        # 1: (0,0,0); 2: /vec{a}; 3: /vec{a} + /vec{c}; 4: \vec{c}; 5: \vec{b}; 6: \vec{a} + \vec{b}; 7: \vec{a}+\vec{b}+\vec{c}; 8: \vec{b} + \vec{c}
        # The shifted eight vertices are each of the above vertices minus the rotation center array.
        vertex = np.array([0.0]*8*3, dtype=np.float)
        vertex.shape = 8,3
        vertex[0,:] = shifted_origin_arr
        vertex[1,:] = poscar_dict['l_arr'][0,:] - rot_center_arr  #\vec{a}
        vertex[2,:] = poscar_dict['l_arr'][0,:] + poscar_dict['l_arr'][2,:] - rot_center_arr  
        vertex[3,:] = poscar_dict['l_arr'][2,:] - rot_center_arr  #\vec{c}
        vertex[4,:] = poscar_dict['l_arr'][1,:] - rot_center_arr  #\vec{b}
        vertex[5,:] = poscar_dict['l_arr'][0,:] + poscar_dict['l_arr'][1,:] - rot_center_arr
        vertex[6,:] = poscar_dict['l_arr'][0,:] + poscar_dict['l_arr'][1,:] + poscar_dict['l_arr'][2,:] - rot_center_arr
        vertex[7,:] = poscar_dict['l_arr'][1,:] + poscar_dict['l_arr'][2,:] - rot_center_arr
        #Define Rotation matrix. Passive rotation(rotate the axis, not the atoms)
        t = funcs.t_euler(euler_angle_type, angle_unit, phi, theta, psi)
        #Inverse of the transformation matrix t_inv
        #The tranformation matrix t acts on the axes, making the axes to rotate in a counter-clockwise manner. This means that atom coordinates rotate in a clockwise manner
        #The inverse of t_inv should act on coordinates of atoms and cause coordinates of atoms to rotate in a counter-clockwise manner. 
        t_inv = np.linalg.inv(t)
        #Coordinate rotation based on Eulerian angles
        pos_shift_arr = poscar_dict['pos_arr'][:,3:6].copy()
        t_pos_shift_arr = poscar_dict['pos_arr'][:,3:6].copy()
        d_pos_shift_arr = poscar_dict['pos_arr'][:,3:6].copy()
        for i in range(0,len(poscar_dict['atom_species_arr'])):
            pos_shift_arr[i,0:3] = poscar_dict['pos_arr'][i,3:6] - rot_center_arr[:]
            t_pos_shift_arr[i,0:3] = pos_shift_arr[i,0:3].copy()
            d_pos_shift_arr[i,0:3] = pos_shift_arr[i,0:3].copy()
            #Atom Coordinate rotation(counter-clockwise) based on Eulerian angles
            t_pos_shift_arr[i,:] = np.dot(t_inv,pos_shift_arr[i,:])
            #Find the atoms which are at the boundaries, if found then displace the boundary atom to the mirror boundary
            for j in range(0,3):                        
                if np.abs(pos_shift_arr[i,:] - shifted_origin_arr)[j] <= pos_precision:
                    d_pos_shift_arr[i,j] = pos_shift_arr[i,j] + (poscar_dict['l_arr'][0,:] + poscar_dict['l_arr'][1,:] + poscar_dict['l_arr'][2,:])[j]
                if np.abs(pos_shift_arr[i,:] - shifted_diagonal_arr)[j] <= pos_precision:
                    d_pos_shift_arr[i,j] = pos_shift_arr[i,j] - (poscar_dict['l_arr'][0,:] + poscar_dict['l_arr'][1,:] + poscar_dict['l_arr'][2,:])[j]
                        
        #maximum and mininum projected position of the rotated vertices in x[1,00],y[0,1,0],z[0,0,1] directions. This is calculated for the purpose of plotting.
        t_vertex_min = np.array([999.99, 999.99, 999.99])
        t_vertex_max = np.array([-999.99, -999.99, -999.99])
        t_vertex = np.array([0.0]*8*3, dtype=np.float)
        t_vertex.shape = 8,3
        #Rotation of verteces in counter-clockwise mannner
        for i in range(0, 8):
            t_vertex[i:] = np.dot(t_inv, vertex[i,:])
        t_vertex_min[0] = min(t_vertex[:,0])
        t_vertex_max[0] = max(t_vertex[:,0])
        t_vertex_min[1] = min(t_vertex[:,1])
        t_vertex_max[1] = max(t_vertex[:,1])
        t_vertex_min[2] = min(t_vertex[:,2])
        t_vertex_max[2] = max(t_vertex[:,2])
        #########################
        # Rotate axis indicator
        #########################
        #The axis indicator is initially located at the origin of the model.
        #Define axis indicator initial reference length.
        #axis_indic_ref_length is the length of the unit axis indicator vector in three dimensions. Initially, axis indicator in three directions have the same lengths.
        #axis indicator initial reference length = max(Projected BoxLength after rotation) * $s(the axis scale factor)
        axis_indic_ref_length = max(t_vertex_max[:] - t_vertex_min[:]) * indicator_len_scale
        axisIndicorigin_arr = np.array([0,0,0]) # origin of the axis indicator 
        #vectors of the axis indicator
        indic_vec_arr = poscar_dict['l_arr'].copy()
        indic_vec_arr[0,:] = poscar_dict['l_arr'][0,:] / np.linalg.norm(poscar_dict['l_arr'][0,:]) * axis_indic_ref_length
        indic_vec_arr[1,:] = poscar_dict['l_arr'][1,:] / np.linalg.norm(poscar_dict['l_arr'][1,:]) * axis_indic_ref_length
        indic_vec_arr[2,:] = poscar_dict['l_arr'][2,:] / np.linalg.norm(poscar_dict['l_arr'][2,:]) * axis_indic_ref_length
        #Shift the axis indicator origin and three axis indicator basis vectors
        axisIndicOriginShiftArr = axisIndicorigin_arr - rot_center_arr
        indic_vec_shift_arr = indic_vec_arr.copy()
        for i in range(0,3):
            indic_vec_shift_arr[i,:] = indic_vec_arr[i,:] - rot_center_arr
        #rotate the axis indicator origin and three axis indicator basis vectors according to Eulerian angles
        t_axis_indic_origin_shift_arr = axisIndicOriginShiftArr.copy()
        t_indic_vec_shift_arr = indic_vec_shift_arr.copy()
        t_axis_indic_origin_shift_arr = np.dot(t_inv, axisIndicOriginShiftArr)
        #Rotation of axis indicator point in counter-clockwise manner
        for i in range(0,3):
            t_indic_vec_shift_arr[i,:] = np.dot(t_inv, indic_vec_shift_arr[i,:])
        write_log_logic = False
        if write_log_logic == True:
            funcs.write_log(
                logfile,
                'poscar_file_path = ' + str(poscar_file_path) + '\n' +
                'plot_poscar(poscar_file_path, view_axis=' + str(view_axis) + ', phi=' + str(phi) + ', theta=' + str(theta) + ', psi=' + str(psi) +
                ', draw_mirror_atom=' + str(draw_mirror_atom) + ', box_on=' + str(box_on) + ', axis_indicator=' + str(axis_indicator) +
                ', plot_cell_basis_vector_label=' + str(plot_cell_basis_vector_label) + ', plot_atom_label=' + str(plot_atom_label) + ', fig_format=' + str(fig_format) +
                ', logfile=' + logfile + '\n' +
                'Original Origin vector = ' + str(origin_arr) + '\n' +
                'Original Diagonal vector = ' + str(diagonal_arr) + '\n' +
                'Shifted Origin vector = ' + ' '.join([str(x) for x in shifted_origin_arr]) +  '\n' +
                'Shifted Diagonal vector = ' + ' '.join([str(x) for x in shifted_diagonal_arr]) + '\n' +
                'Original RotCenter =  ' + ' '.join([str(x) for x in rot_center_arr]) + '\n' +
                'Shifted 8 Vertices =  ' + '\n' + str(vertex) + '\n' +
                'T = ' + '\n' + str(t) + '\n' +
                'Shifted and rotated 8 Vertices = ' + '\n' + str(t_vertex)
                )
        else:
            pass
        return pos_shift_arr, t, t_pos_shift_arr, d_pos_shift_arr, t_vertex, t_vertex_min, t_vertex_max, shifted_origin_arr, shifted_diagonal_arr, t_axis_indic_origin_shift_arr, t_indic_vec_shift_arr, axis_indic_ref_length

    pos_shift_arr, t, t_pos_shift_arr, d_pos_shift_arr, t_vertex, t_vertex_min, t_vertex_max, shifted_origin_arr, shifted_diagonal_arr, t_axis_indic_origin_shift_arr, t_indic_vec_shift_arr, axis_indic_ref_length = disp_and_rotation(poscar_file_path, phi, theta, psi)
    #Inverse of the transformation matrix t_inv
    #The tranformation matrix t acts on the axes, making the axes to rotate in a counter-clockwise manner. This means that atom coordinates rotate in a clockwise manner
    #The inverse of t_inv should act on coordinates of atoms and cause coordinates of atoms to rotate in a counter-clockwise manner. 
    t_inv = np.linalg.inv(t)    

    ###############################################
    #Find equivalent mirror atoms at the boundary
    ###############################################
    # Find mirror atoms at periodic boundary. d_ donotes 'displaced'
    # Find equivalent mirror atoms at the boundary. One atom can have 1, 3, or 7 equivalent atoms at the boundary due to periodic boundary conditions          
    t_pos_shift_add_mirror_atoms_arr = t_pos_shift_arr.copy()
    atom_species_add_mirror_atoms_arr = poscar_dict['atom_species_arr'].copy()
    atomname_arr = np.array(poscar_dict['atomname_list'])
    atomname_add_mirror_atom_arr = atomname_arr.copy()
    if len(poscar_dict['added_atom_data'][0,:]) != 0:
        if colormap_column_indx > len(poscar_dict['added_atom_data'][0,:]) or colormap_column_indx <= 0:
            funcs.write_log(logfile, '## Error: colormap_column_indx out of range.' +
                            ' The value of colormap_column_indx should be within the range: 0--' + str(len(poscar_dict['added_atom_data'][0,:])) +
                            '. We now set colormap_column_indx = 1. The data in the first column of the added atom data will be plotted.'
                            )
            colormap_column_indx = 1
        selected_added_atom_data_arr = poscar_dict['added_atom_data'][:, colormap_column_indx - 1].copy()
        if min(selected_added_atom_data_arr) == max(selected_added_atom_data_arr):
            zero_color_gradient = True
        else:
            zero_color_gradient = False
    for i in range(0,len(poscar_dict['atom_species_arr'])):
        if draw_mirror_atom == True:
            pos_dup1_arr = d_pos_shift_arr[i,:].copy()        
            pos_dup1_arr = d_pos_shift_arr[i,:].copy()
            pos_dup2_arr = d_pos_shift_arr[i,:].copy()
            pos_dup3_arr = d_pos_shift_arr[i,:].copy()
            pos_dup4_arr = d_pos_shift_arr[i,:].copy()
            pos_dup5_arr = d_pos_shift_arr[i,:].copy()
            pos_dup6_arr = d_pos_shift_arr[i,:].copy()
            pos_dup7_arr = d_pos_shift_arr[i,:].copy()
            #Find equivalent position of atoms at periodic boundary
            if d_pos_shift_arr[i,0] == pos_shift_arr[i,0] and d_pos_shift_arr[i,1] == pos_shift_arr[i,1] and d_pos_shift_arr[i,2] == pos_shift_arr[i,2]:
                n_dup = 0
                continue
            elif d_pos_shift_arr[i,0] != pos_shift_arr[i,0] and d_pos_shift_arr[i,1] == pos_shift_arr[i,1] and d_pos_shift_arr[i,2] == pos_shift_arr[i,2]:
                n_dup = 1
            elif d_pos_shift_arr[i,0] == pos_shift_arr[i,0] and d_pos_shift_arr[i,1] != pos_shift_arr[i,1] and d_pos_shift_arr[i,2] == pos_shift_arr[i,2]:
                n_dup = 1   
            elif d_pos_shift_arr[i,0] == pos_shift_arr[i,0] and d_pos_shift_arr[i,1] == pos_shift_arr[i,1] and d_pos_shift_arr[i,2] != pos_shift_arr[i,2]:
                n_dup = 1
            elif d_pos_shift_arr[i,0] != pos_shift_arr[i,0] and d_pos_shift_arr[i,1] != pos_shift_arr[i,1] and d_pos_shift_arr[i,2] == pos_shift_arr[i,2]:
                n_dup = 3
                pos_dup1_arr[0] = pos_shift_arr[i,0]
                pos_dup2_arr[1] = pos_shift_arr[i,1]
            elif d_pos_shift_arr[i,0] == pos_shift_arr[i,0] and d_pos_shift_arr[i,1] != pos_shift_arr[i,1] and d_pos_shift_arr[i,2] != pos_shift_arr[i,2]:
                n_dup = 3
                pos_dup1_arr[1] = pos_shift_arr[i,1]
                pos_dup2_arr[2] = pos_shift_arr[i,2]
            elif d_pos_shift_arr[i,0] != pos_shift_arr[i,0] and d_pos_shift_arr[i,1] == pos_shift_arr[i,1] and d_pos_shift_arr[i,2] != pos_shift_arr[i,2]:
                n_dup = 3
                pos_dup1_arr[0] = pos_shift_arr[i,0]
                pos_dup2_arr[2] = pos_shift_arr[i,2]
            elif d_pos_shift_arr[i,0] != pos_shift_arr[i,0] and d_pos_shift_arr[i,1] != pos_shift_arr[i,1] and d_pos_shift_arr[i,2] != pos_shift_arr[i,2]:  
                n_dup = 7
                pos_dup1_arr[0] = pos_shift_arr[i,0]
                pos_dup2_arr[1] = pos_shift_arr[i,1]                
                pos_dup3_arr[2] = pos_shift_arr[i,2]     
                pos_dup4_arr[0] = pos_shift_arr[i,0]
                pos_dup4_arr[1] = pos_shift_arr[i,1]   
                pos_dup5_arr[1] = pos_shift_arr[i,1]                
                pos_dup5_arr[2] = pos_shift_arr[i,2] 
                pos_dup6_arr[0] = pos_shift_arr[i,0]
                pos_dup6_arr[2] = pos_shift_arr[i,2]
            if n_dup == 1:
                t_pos_dup1_arr = np.array([np.dot(t_inv,pos_dup1_arr[:])])
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup1_arr),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                if len(poscar_dict['added_atom_data'][0,:]) != 0:
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
            elif n_dup == 3:
                t_pos_dup1_arr = np.array([np.dot(t_inv,pos_dup1_arr[:])])
                t_pos_dup2_arr = np.array([np.dot(t_inv,pos_dup2_arr[:])])
                t_pos_dup3_arr = np.array([np.dot(t_inv,pos_dup3_arr[:])])
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup1_arr),axis = 0)
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup2_arr),axis = 0)
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup3_arr),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                if len(poscar_dict['added_atom_data'][0,:]) != 0:
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
            elif n_dup == 7:
                t_pos_dup1_arr = np.array([np.dot(t_inv,pos_dup1_arr[:])])
                t_pos_dup2_arr = np.array([np.dot(t_inv,pos_dup2_arr[:])])
                t_pos_dup3_arr = np.array([np.dot(t_inv,pos_dup3_arr[:])])
                t_pos_dup4_arr = np.array([np.dot(t_inv,pos_dup4_arr[:])])
                t_pos_dup5_arr = np.array([np.dot(t_inv,pos_dup5_arr[:])])
                t_pos_dup6_arr = np.array([np.dot(t_inv,pos_dup6_arr[:])])
                t_pos_dup7_arr = np.array([np.dot(t_inv,pos_dup7_arr[:])])
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup1_arr),axis = 0)
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup2_arr),axis = 0)
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup3_arr),axis = 0)
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup4_arr),axis = 0)
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup5_arr),axis = 0)
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup6_arr),axis = 0)
                t_pos_shift_add_mirror_atoms_arr = np.concatenate((t_pos_shift_add_mirror_atoms_arr,t_pos_dup7_arr),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array([poscar_dict['atom_species_arr'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array([poscar_dict['atomname_list'][i]])),axis = 0)
                if len(poscar_dict['added_atom_data'][0,:]) != 0:
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)
                    selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([poscar_dict['added_atom_data'][i, colormap_column_indx - 1]])),axis = 0)

    #Build KDTree:
    # Define the point where the viewer(our eyes) is
    if view_axis == 'x':
        viewer_point_arr = np.array([[t_vertex_max[0], 0, 0]]) * 2.0
    if view_axis == 'y':
        viewer_point_arr = np.array([[0, t_vertex_max[0], 0]]) * 2.0
    if view_axis == 'z':
        viewer_point_arr = np.array([[0, 0, t_vertex_max[0]]]) * 2.0
    #Build rotated atom list with one viwer point added. Find the points which is closest to the viewer point.
    t_pos_shift_add_viewer_point_arr = t_pos_shift_arr.copy()
    t_pos_shift_add_viewer_point_arr = np.concatenate((t_pos_shift_add_viewer_point_arr, viewer_point_arr), axis =0)
    #Calculate nearest neighbor inofrmation
    tree = KDTree(t_pos_shift_add_viewer_point_arr)
    dists, indices = tree.query([t_pos_shift_add_viewer_point_arr[len(poscar_dict['atom_species_arr'])]], k=len(poscar_dict['atom_species_arr']))
    view_atom_indx_arr = indices[0][:]

    #Build rotated atom list with mirror atoms at the boundary and one viwer point added. Find the points which is closest to the viewer point.
    t_pos_shift_add_mirror_atoms_and_viewer_point_arr = t_pos_shift_add_mirror_atoms_arr.copy()
    t_pos_shift_add_mirror_atoms_and_viewer_point_arr = np.concatenate((t_pos_shift_add_mirror_atoms_and_viewer_point_arr, viewer_point_arr), axis =0)
    #Calculate nearest neighbor atoms from the viwer point(i.e., from our eyes). KDTree method is adopted to find the nearest neighbor list
    tree = KDTree(t_pos_shift_add_mirror_atoms_and_viewer_point_arr)
    dists, indices = tree.query([t_pos_shift_add_mirror_atoms_and_viewer_point_arr[len(t_pos_shift_add_mirror_atoms_and_viewer_point_arr)-1]], k=len(t_pos_shift_add_mirror_atoms_and_viewer_point_arr))
    view_atom_indx_arr = indices[0][:]
    atom_species_add_mirror_atoms_arr = np.concatenate((atom_species_add_mirror_atoms_arr,np.array(['Vi','Vi','Vi'])),axis = 0) #Add dimension to atom_species_add_mirror_atoms_arr due to the viwer point is added.
    atomname_add_mirror_atom_arr = np.concatenate((atomname_add_mirror_atom_arr,np.array(['Vi','Vi','Vi'])),axis = 0) #Add dimension to atomname_add_mirror_atom_arr due to the viwer point is added.
##    if len(poscar_dict['added_atom_data'][0,:]) != 0:
##        selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([0,0,0])),axis = 0)

    #################################
    #Plotting
    #################################
    #Calculate the maximum projected box length in three dimensions
    box_length_vec_list = t_vertex_max[:] - t_vertex_min[:]

    if view_axis == 'x':
        hori_indx = 1
        vert_indx = 2
        fig_length = box_length_vec_list[hori_indx]
        fig_height = box_length_vec_list[vert_indx]
    elif view_axis == 'y':
        hori_indx = 0
        vert_indx = 2
        fig_length = box_length_vec_list[hori_indx]
        fig_height = box_length_vec_list[vert_indx]
    elif view_axis == 'z':
        hori_indx = 0
        vert_indx = 1
        fig_length = box_length_vec_list[hori_indx]
        fig_height = box_length_vec_list[vert_indx]

    if axis_indicator == True:
        #Define the multiplier of axis indicator length, This factor 'mA' places the origin of the axis indicator a distance mA*(axis length) away from the model.
        #axis_indic_ref_length is the length of the unit axis indicator vector in three directions. The length of the unit axis indicator vector in three directions are the same
        mA = 1.5
    max_projected_box_length = max(box_length_vec_list)
    if axis_indicator == True:      
        max_angstrom_width = box_length_vec_list[hori_indx] + axis_indic_ref_length * mA + axis_indic_ref_length
        fig_aspect_ratio = box_length_vec_list[vert_indx] / max_angstrom_width 
    else:
        max_angstrom_width = box_length_vec_list[hori_indx]
        fig_aspect_ratio = box_length_vec_list[vert_indx] / max_angstrom_width

    ##########################
    #Initialize pyplot figure
    ##########################
    #the scale factor 'scale' is useful for scaling the size of pyplot figure. The tested optimum scaling factor is 0.45
    scale = 0.45
    fig_size = [max_angstrom_width *scale,max_angstrom_width* fig_aspect_ratio*scale]
    plt.rcParams['figure.figsize'] = (fig_size[0],fig_size[1])
    fig_poscar = plt.figure('fig_poscar')
    fig_poscar.add_subplot(111)
    plt.subplots_adjust(bottom = 0, left = 0, top = 1, right = 1, wspace = 0, hspace = 0)
    ax = plt.gca()
    ##########################
    #Initialize TikZ figure
    ##########################
    #Define TikZ file directory
    tikz_dir = output_dir + '/poscar_tikz_latex_files'
    funcs.mkdir(tikz_dir)
    formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time.time()))
    system_folder = str(str(os.path.abspath(poscar_file_path)).split('/')).split('\\\\')[-2]
    tikz_texfile = tikz_dir + '/TikZ_' + str(formatted_time) + '_' + system_folder + '_' + str(euler_angle_type) + '_phi' + str(phi) + '_theta' + str(theta) + '_psi' + str(psi) + '_' + funcs.logic_retn_val(box_on,'BoxOn','BoxOff') + '_' + funcs.logic_retn_val(axis_indicator,'AxisOn','AxisOff') + '.tex' 
    
    #Define scaled element radius according to Atomic Radius. This is for the purpose of plotting. 
    reference_atom_radius = periodic_table_dict['atomic_radius'][model_basis_elmt]
    scaled_atomic_radius_for_tikz_dict = {}
    scaled_atomic_radius_for_pyplot_dict = {}
    reduc_fac = 80
    enlargefac_pyplot = 13
    enlargefac_tikz = 0.15
    for i in periodic_table_dict['atomic_radius']:
        if periodic_table_dict['atomic_radius'][i] != None:
            scaled_atomic_radius_for_tikz_dict[i] = (periodic_table_dict['atomic_radius'][i] - reduc_fac) / (reference_atom_radius - reduc_fac) * enlargefac_tikz * radius_scale_fac
            scaled_atomic_radius_for_pyplot_dict[i] = (periodic_table_dict['atomic_radius'][i] - reduc_fac) / (reference_atom_radius - reduc_fac) * enlargefac_pyplot * radius_scale_fac 
        else:
            pass
    #Define the number of round off digits of numbers in LaTeX TikZ picture.
    round_digits = 4

    #the scale factor 'scale' is useful for scaling the picture size in TikZ figure. The tested optimum scaling factor is 0.5
    scale = 0.5
    #Write header of Tikz picture
    with open(tikz_texfile,'w') as f1:
        f1.write('\documentclass[a4paper,12]{article}' + '\n' +
                 r'\usepackage{graphicx}' + '\n' +
                 r'\usepackage{tikz}' + '\n' +
                 r'\begin{document}' + '\n'
                 )
    for i in range(0, len(poscar_dict['elmt_species_arr'])):
        j = elmt_color_tikz[i]
        with open(tikz_texfile,'a') as f1:
            f1.write('\pgfdeclareradialshading{ballshading' +
                     str(j) +'}{\\pgfpoint{-10bp}{10bp}}{color(0bp)=(' +
                     str(j) +
                     '!15!white); color(9bp)=(' + str(j) + '!75!white);color(18bp)=(' + str(j) + '!70!black); color(25bp)=(' + str(j) + '!50!black);color(50bp)=(black)}' + '\n'
                     )
    with open(tikz_texfile,'a') as f1:
        f1.write(r'\begin{figure}' + '\n' +
                 '\centering' + '\n' +
                 r'\begin{tikzpicture}[scale=' + str(round(scale,round_digits)) + ']' + '\n')

    #######################
    #Draw the bounding box
    #######################
    if box_on == True:
        #draw bounding box for pyplot figure
        line1 = [(t_vertex[0,hori_indx], t_vertex[0,vert_indx]), (t_vertex[1,hori_indx], t_vertex[1,vert_indx])]
        line2 = [(t_vertex[1,hori_indx], t_vertex[1,vert_indx]), (t_vertex[2,hori_indx], t_vertex[2,vert_indx])]
        line3 = [(t_vertex[2,hori_indx], t_vertex[2,vert_indx]), (t_vertex[3,hori_indx], t_vertex[3,vert_indx])]
        line4 = [(t_vertex[3,hori_indx], t_vertex[3,vert_indx]), (t_vertex[0,hori_indx], t_vertex[0,vert_indx])]
        line5 = [(t_vertex[0,hori_indx], t_vertex[0,vert_indx]), (t_vertex[4,hori_indx], t_vertex[4,vert_indx])]
        line6 = [(t_vertex[4,hori_indx], t_vertex[4,vert_indx]), (t_vertex[5,hori_indx], t_vertex[5,vert_indx])]
        line7 = [(t_vertex[5,hori_indx], t_vertex[5,vert_indx]), (t_vertex[6,hori_indx], t_vertex[6,vert_indx])]
        line8 = [(t_vertex[6,hori_indx], t_vertex[6,vert_indx]), (t_vertex[7,hori_indx], t_vertex[7,vert_indx])]
        line9 = [(t_vertex[7,hori_indx], t_vertex[7,vert_indx]), (t_vertex[4,hori_indx], t_vertex[4,vert_indx])]
        line10 = [(t_vertex[1,hori_indx], t_vertex[1,vert_indx]), (t_vertex[5,hori_indx], t_vertex[5,vert_indx])]
        line11 = [(t_vertex[2,hori_indx], t_vertex[2,vert_indx]), (t_vertex[6,hori_indx], t_vertex[6,vert_indx])]
        line12 = [(t_vertex[3,hori_indx], t_vertex[3,vert_indx]), (t_vertex[7,hori_indx], t_vertex[7,vert_indx])]        
        (line1_xs, line1_ys) = zip(*line1)
        (line2_xs, line2_ys) = zip(*line2)
        (line3_xs, line3_ys) = zip(*line3)
        (line4_xs, line4_ys) = zip(*line4)
        (line5_xs, line5_ys) = zip(*line5)
        (line6_xs, line6_ys) = zip(*line6)
        (line7_xs, line7_ys) = zip(*line7)
        (line8_xs, line8_ys) = zip(*line8)
        (line9_xs, line9_ys) = zip(*line9)
        (line10_xs, line10_ys) = zip(*line10)
        (line11_xs, line11_ys) = zip(*line11)
        (line12_xs, line12_ys) = zip(*line12)
        line_xs = np.array([line1_xs,line2_xs,line3_xs,line4_xs,line5_xs,line6_xs,line7_xs,line8_xs,line9_xs])
        line_ys = np.array([line1_ys,line2_ys,line3_ys,line4_ys,line5_ys,line6_ys,line7_ys,line8_ys,line9_ys])
        plt.plot(line_xs, line_ys, linestyle = '-',linewidth = 0.5, color = 'black')
        plt.plot(line10_xs, line10_ys, linestyle = '-',linewidth = 0.5, color = 'black')
        plt.plot(line11_xs, line11_ys, linestyle = '-',linewidth = 0.5, color = 'black')
        plt.plot(line12_xs, line12_ys, linestyle = '-',linewidth = 0.5, color = 'black')
        ax.set(aspect='equal')
        #draw bounding box for TikZ figure
        with open(tikz_texfile,'a') as f1:
            f1.write('\draw [help lines] ' +
                     '(' + str(round(t_vertex[0,hori_indx],round_digits)) + ',' + str(round(t_vertex[0,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[1,hori_indx],round_digits)) + ',' + str(round(t_vertex[1,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[2,hori_indx],round_digits)) + ',' + str(round(t_vertex[2,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[3,hori_indx],round_digits)) + ',' + str(round(t_vertex[3,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[0,hori_indx],round_digits)) + ',' + str(round(t_vertex[0,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[4,hori_indx],round_digits)) + ',' + str(round(t_vertex[4,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[5,hori_indx],round_digits)) + ',' + str(round(t_vertex[5,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[6,hori_indx],round_digits)) + ',' + str(round(t_vertex[6,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[7,hori_indx],round_digits)) + ',' + str(round(t_vertex[7,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[4,hori_indx],round_digits)) + ',' + str(round(t_vertex[4,vert_indx],round_digits)) + ');' + '\n' +
                     '\draw [help lines] ' +
                     '(' + str(round(t_vertex[1,hori_indx],round_digits)) + ',' + str(round(t_vertex[1,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[5,hori_indx],round_digits)) + ',' + str(round(t_vertex[5,vert_indx],round_digits)) + ');' + '\n' +
                     '\draw [help lines] ' +
                     '(' + str(round(t_vertex[2,hori_indx],round_digits)) + ',' + str(round(t_vertex[2,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[6,hori_indx],round_digits)) + ',' + str(round(t_vertex[6,vert_indx],round_digits)) + ');' + '\n' +
                     '\draw [help lines] ' +
                     '(' + str(round(t_vertex[3,hori_indx],round_digits)) + ',' + str(round(t_vertex[3,vert_indx],round_digits)) + ') -- ' +
                     '(' + str(round(t_vertex[7,hori_indx],round_digits)) + ',' + str(round(t_vertex[7,vert_indx],round_digits)) + ');' + '\n'
                     )

    ##################
    #Draw atoms
    ##################
    plt.axis('off')
##    # the plt.scatter only aims at plotting the color bar 
##    # Use plt.scatter won't ensure the atom colors plotted according to the order of the view_atom_indx_arr
##    plt.scatter(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[view_atom_indx_arr[1:][::-1],hori_indx],
##                t_pos_shift_add_mirror_atoms_and_viewer_point_arr[view_atom_indx_arr[1:][::-1],vert_indx],
##                c = selected_added_atom_data_arr[view_atom_indx_arr[1:][::-1]],
##                cmap = colormap_mode,
##                vmin = colormap_vmin,
##                vmax = colormap_vmax,
##                s = 0,
##                alpha = 1,
##                )
##    plt.colorbar()
    if vmin_color == None or vmin_color == 'None' or vmax_color == None or vmax_color == 'None':
        vmin_color = 'blue'
        vmax_color = 'red'
    with open(tikz_texfile,'a') as f1:
        for i in view_atom_indx_arr[1:][::-1]:
            elmtindx = int(np.argwhere(poscar_dict['elmt_species_arr'] == atom_species_add_mirror_atoms_arr[i]))
##              plt.scatter(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],
##                      t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],
##                      color=periodic_table_dict['elmt_color'][atom_species_add_mirror_atoms_arr[i]],
##                      s=scaled_atomic_radius_for_pyplot_dict[atom_species_add_mirror_atoms_arr[i]])
            #Draw atoms for pyplot figure (with shiny shpere effect)
            original_ball_size = scaled_atomic_radius_for_pyplot_dict[atom_species_add_mirror_atoms_arr[i]]
            if len(poscar_dict['added_atom_data'][0,:]) != 0:
                colormap_normalized_data, colormap_vmin, colormap_vmax = funcs.data_normalize(input_data_value = float(selected_added_atom_data_arr[i]),
                                                                                              data_list = selected_added_atom_data_arr.astype(float),
                                                                                              colormap_vmin = colormap_vmin,
                                                                                              colormap_vmax = colormap_vmax)
            #Real color background
            if draw_colormap == False or (draw_colormap == True and len(poscar_dict['added_atom_data'][0,:]) == 0):
                plt.plot(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],
                         t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],
                         marker = '.',
                         markersize = original_ball_size,
                         markeredgecolor = 'None',
                         color = periodic_table_dict['elmt_color'][atom_species_add_mirror_atoms_arr[i]],
                         alpha = 1,
                         linestyle = '',
                         )
            if draw_colormap == True and len(poscar_dict['added_atom_data'][0,:]) != 0:
                # higher spectrum of the colormap
                plt.plot(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],
                         t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],
                         marker = '.',
                         markersize = original_ball_size,
                         markeredgecolor = 'None',
                         color = vmax_color,
                         alpha = colormap_normalized_data,
                         linestyle = '',
                         )
                # lower spectrum of the colormap
                plt.plot(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],
                         t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],
                         marker = '.',
                         markersize = original_ball_size,
                         markeredgecolor = 'None',
                         color = vmin_color,
                         alpha = (1 - colormap_normalized_data),
                         linestyle = '',
                         )

            if fig_format == 'eps' or fig_format == 'ps' :
                # the transparent white foreground does not work with the eps and the ps figures. Thus the shiny sphere effect won't work for eps and ps format figures.
                pass
            else:
                #white foreground
                for istep in np.arange(0, 1, 0.1):
                    scaled_ball_size = (1 - istep) * original_ball_size
                    plt.plot(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],
                             t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],
                             marker = '.',
                             markersize = scaled_ball_size,
                             markeredgecolor = 'None',
                             color = 'white',
                             alpha = istep,
                             linestyle = ''
                             )
            #Draw atoms for TikZ figure
            f1.write('\pgfpathcircle{\pgfpoint{' +
                     str(round(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],round_digits)) +
                     'cm}{' + str(round(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],round_digits)) +
                     'cm}}{' + str(round(scaled_atomic_radius_for_tikz_dict[atom_species_add_mirror_atoms_arr[i]],round_digits)) +
                     'cm}' + '\n' +
                     '\pgfshadepath{ballshading' + str(elmt_color_tikz[elmtindx]) + '}{0}' + '\n' +
                     '\pgfusepath{}' + '\n')
            if plot_atom_label == True:
                #Draw atoms for pyplot figure
                plt.text(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],
                         t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],
                         atomname_add_mirror_atom_arr[i],
                         size = label_size)   
                #Draw atoms for TikZ figure
                f1.write(r'\node [above] at (' + str(round(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],round_digits)) + ','
                         + str(round(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],round_digits)) + r') {\scriptsize ' + atomname_add_mirror_atom_arr[i] + '};' + '\n')
    # add the colorbar of the colormap
    if draw_colormap == True and len(poscar_dict['added_atom_data'][0,:]) != 0:
        if colorbar_alignment == 'vertical':
            plt.subplots_adjust(bottom = 0, left = 0, top = 1, right = 0.97, wspace = 0, hspace = 0)
        elif colorbar_alignment == 'horizontal':
            plt.subplots_adjust(bottom = 0.09, left = 0, top = 1, right = 1, wspace = 0, hspace = 0)
        funcs.userdefined_colormap(selected_added_atom_data_arr , colormap_vmin = colormap_vmin, colormap_vmax = colormap_vmax, vmax_color = vmax_color, vmin_color = vmin_color,
                 alignment = colorbar_alignment, left = None, bottom = None, width = None, height = None)

    ####################
    #Add axis indicator
    ####################
    #First put axis indicator in the rotation center, then rotate the axis indicator according to Eulerian angles
    #Then displace the axis indicator to desiginated position according to a specific displacement, usually it is placed in the bottom-right of the figure.
    if axis_indicator == True: 
        #Adjust the axis indicator to specific position
        with open(tikz_texfile,'a') as f1: 
            t_axis_indic_origin_shift_adjust_arr = t_axis_indic_origin_shift_arr.copy()
            t_indic_vec_shift_adjust_arr = t_indic_vec_shift_arr.copy()
            relative_pos_str = np.array([])
            adjust_dist = np.array([0.0]*2, dtype = float)
            adjust_dist[0] = t_vertex_max[hori_indx] - t_axis_indic_origin_shift_arr[hori_indx]
            adjust_dist[1] = t_vertex_min[vert_indx] - t_axis_indic_origin_shift_arr[vert_indx]
            t_axis_indic_origin_shift_adjust_arr[hori_indx] = t_axis_indic_origin_shift_arr[hori_indx] + adjust_dist[0] + axis_indic_ref_length * mA
            t_axis_indic_origin_shift_adjust_arr[vert_indx] = t_axis_indic_origin_shift_arr[vert_indx] + adjust_dist[1]
            t_indic_vec_shift_adjust_arr[0,hori_indx] = t_indic_vec_shift_arr[0,hori_indx] + adjust_dist[0] + axis_indic_ref_length * mA
            t_indic_vec_shift_adjust_arr[0,vert_indx] = t_indic_vec_shift_arr[0,vert_indx] + adjust_dist[1]
            t_indic_vec_shift_adjust_arr[1,hori_indx] = t_indic_vec_shift_arr[1,hori_indx] + adjust_dist[0] + axis_indic_ref_length * mA
            t_indic_vec_shift_adjust_arr[1,vert_indx] = t_indic_vec_shift_arr[1,vert_indx] + adjust_dist[1]
            t_indic_vec_shift_adjust_arr[2,hori_indx] = t_indic_vec_shift_arr[2,hori_indx] + adjust_dist[0] + axis_indic_ref_length * mA
            t_indic_vec_shift_adjust_arr[2,vert_indx] = t_indic_vec_shift_arr[2,vert_indx] + adjust_dist[1]
            line1 = [(t_axis_indic_origin_shift_adjust_arr[hori_indx], t_axis_indic_origin_shift_adjust_arr[vert_indx]), (t_indic_vec_shift_adjust_arr[0,hori_indx], t_indic_vec_shift_adjust_arr[0,vert_indx])]
            line2 = [(t_axis_indic_origin_shift_adjust_arr[hori_indx], t_axis_indic_origin_shift_adjust_arr[vert_indx]), (t_indic_vec_shift_adjust_arr[1,hori_indx], t_indic_vec_shift_adjust_arr[1,vert_indx])]
            line3 = [(t_axis_indic_origin_shift_adjust_arr[hori_indx], t_axis_indic_origin_shift_adjust_arr[vert_indx]), (t_indic_vec_shift_adjust_arr[2,hori_indx], t_indic_vec_shift_adjust_arr[2,vert_indx])]
            (line1_xs, line1_ys) = zip(*line1)
            (line2_xs, line2_ys) = zip(*line2)
            (line3_xs, line3_ys) = zip(*line3)
            for i in range(0,3):
                if t_indic_vec_shift_adjust_arr[i,hori_indx] != 0 and t_indic_vec_shift_adjust_arr[i,vert_indx] != 0:
                    relative_pos_str = np.concatenate((relative_pos_str, ['right']),axis = 0)
                elif t_indic_vec_shift_adjust_arr[i,hori_indx] < 0 and t_indic_vec_shift_adjust_arr[i,vert_indx] == 0:
                    relative_pos_str = np.concatenate((relative_pos_str, ['left']),axis = 0)
                elif t_indic_vec_shift_adjust_arr[i,hori_indx] > 0 and t_indic_vec_shift_adjust_arr[i,vert_indx] == 0:                    
                    relative_pos_str = np.concatenate((relative_pos_str, ['right']),axis = 0)
                elif t_indic_vec_shift_adjust_arr[i,hori_indx] == 0 and t_indic_vec_shift_adjust_arr[i,vert_indx] < 0:                    
                    relative_pos_str = np.concatenate((relative_pos_str, ['below']),axis = 0)
                elif t_indic_vec_shift_adjust_arr[i,hori_indx] == 0 and t_indic_vec_shift_adjust_arr[i,vert_indx] > 0:                    
                    relative_pos_str = np.concatenate((relative_pos_str, ['above']),axis = 0)
                elif t_indic_vec_shift_adjust_arr[i,hori_indx] == 0 and t_indic_vec_shift_adjust_arr[i,vert_indx] == 0:                    
                    relative_pos_str = np.concatenate((relative_pos_str, ['below']),axis = 0)
            #Add axis Indicator for pyplot figure
            plt.plot(line1_xs, line1_ys, linestyle = '-',linewidth = 0.5, color = 'black')
            plt.plot(line2_xs, line2_ys, linestyle = '-',linewidth = 0.5, color = 'black')
            plt.plot(line3_xs, line3_ys, linestyle = '-',linewidth = 0.5, color = 'black')
            plt.text(line1_xs[1], line1_ys[1],'x',size = label_size)
            plt.text(line2_xs[1], line2_ys[1],'y',size = label_size)
            plt.text(line3_xs[1], line3_ys[1],'z',size = label_size)
            #Add axis Indicator for Latex TixZ figure
            f1.write('\draw [->] ' +
                     '(' + str(line1_xs[0]) + ',' + str(line1_ys[0]) + ')--' +
                     '(' + str(line1_xs[1]) + ',' + str(line1_ys[1]) + ');' + '\n' +
                     r'\node [' + relative_pos_str[0] + '] at (' + str(line1_xs[1]) + ',' + str(line1_ys[1]) + ') {x};'
                     )
            f1.write('\draw [->] ' +
                     '(' + str(line2_xs[0]) + ',' + str(line2_ys[0]) + ')--' +
                     '(' + str(line2_xs[1]) + ',' + str(line2_ys[1]) + ');' + '\n' +
                     r'\node [' + relative_pos_str[1] + '] at (' + str(line2_xs[1]) + ',' + str(line2_ys[1]) + ') {y};'
                     )
            f1.write('\draw [->] ' +
                     '(' + str(line3_xs[0]) + ',' + str(line3_ys[0]) + ')--' +
                     '(' + str(line3_xs[1]) + ',' + str(line3_ys[1]) + ');' + '\n' +
                     r'\node [' + relative_pos_str[2] + '] at (' + str(line3_xs[1]) + ',' + str(line3_ys[1]) + ') {z};'
                     )
            #Plot labels of the cell basis vector
            if plot_cell_basis_vector_label == True:
                #Add labels of the cell basis vector for pyplot figure
                plt.text(t_vertex[0,hori_indx], t_vertex[0,vert_indx],'o',size = label_size)
                plt.text(t_vertex[1,hori_indx], t_vertex[1,vert_indx],'a',size = label_size)
                plt.text(t_vertex[4,hori_indx], t_vertex[4,vert_indx],'b',size = label_size)
                plt.text(t_vertex[3,hori_indx], t_vertex[3,vert_indx],'c',size = label_size)
                #Add labels of the cell basis vector for TikZ figure
                f1.write(r'\node [] at (' + str(round(t_vertex[0,hori_indx],round_digits)) + ',' + str(round(t_vertex[0,vert_indx],round_digits)) + ') {o};')
                f1.write(r'\node [] at (' + str(round(t_vertex[1,hori_indx],round_digits)) + ',' + str(round(t_vertex[1,vert_indx],round_digits)) + ') {a};')
                f1.write(r'\node [] at (' + str(round(t_vertex[4,hori_indx],round_digits)) + ',' + str(round(t_vertex[4,vert_indx],round_digits)) + ') {b};')
                f1.write(r'\node [] at (' + str(round(t_vertex[3,hori_indx],round_digits)) + ',' + str(round(t_vertex[3,vert_indx],round_digits)) + ') {c};')
    ######################################
    #Write end lines of TikZ figure code
    ######################################
    with open(tikz_texfile,'a') as f1:
        f1.write('\end{tikzpicture}' + '\n' +
                 '\end{figure}' + '\n' + '\n' +
                 '\end{document}'
                 )
    ###############
    # Dump figures
    ###############
    (tikz_parent_path , tikz_full_filename) = os.path.split(tikz_texfile)
    (tikz_filename , tikz_file_extension) = os.path.splitext(tikz_full_filename)
    if poscar_dict['added_atom_property'] == None:
        fig_file = output_dir + '\Fig_' + tikz_filename.strip('TikZ_')  + '.' + fig_format
    elif poscar_dict['added_atom_property'] != None and draw_colormap == False:
        fig_file = output_dir + '\Fig_' + tikz_filename.strip('TikZ_')  + '.' + fig_format
    elif poscar_dict['added_atom_property'] != None and draw_colormap == True:
        added_atom_property_columns_str = funcs.split_line(line = poscar_dict['added_atom_property_columns'], separator = ' ')[colormap_column_indx - 1]
        fig_file = output_dir + '\Fig_' + tikz_filename.strip('TikZ_')  + '_' + poscar_dict['added_atom_property'] + '_' + added_atom_property_columns_str + '.' + fig_format
    plt.savefig(fig_file,dpi = fig_dpi)
    plt.close()

    funcs.write_log(
        logfile,
        'elmt_color = ' + str(elmt_color) + '\n' +
        'vasp_plot.plot_poscar(' + '\n' +
        '    poscar_file_path=' + 'r\'' + str(poscar_file_path) + '\'' + ',\n' +
        '    euler_angle_type=' + '\'' + str(euler_angle_type) + '\'' + ',\n' +
        '    phi=' + str(phi) + ',\n' +
        '    theta=' + str(theta) + ',\n' +
        '    psi=' + str(psi) + ',\n' +
        '    elmt_color=elmt_color' + ',\n' +
        '    draw_mirror_atom=' + str(draw_mirror_atom) + ',\n' +
        '    box_on=' + str(box_on) + ',\n' +
        '    axis_indicator=' + str(axis_indicator) + ',\n' +
        '    plot_cell_basis_vector_label=' + str(plot_cell_basis_vector_label) + ',\n' +
        '    plot_atom_label=' + str(plot_atom_label) + ',\n' +                 
        '    fig_format=' + '\'' + str(fig_format) + '\'' + ',\n' +
        '    fig_dpi=' + str(fig_dpi) + ',\n' +                
        '    draw_colormap=' + str(draw_colormap) + ',\n' +
        '    colormap_column_indx=' + str(colormap_column_indx) + ',\n' +
        '    colormap_vmin=' + str(colormap_vmin) + ',\n' +
        '    colormap_vmax=' + str(colormap_vmax) + ',\n' +
        '    vmin_color=' + '\'' + str(vmin_color) + '\'' + ',\n' +
        '    vmax_color=' + '\'' + str(vmax_color) + '\'' + ',\n' +
        '    colorbar_alignment=' + '\'' + str(colorbar_alignment) + '\'' + ')\n' +
        'fig_file= ' + 'r\'' + fig_file + '\'' + '\n')
    funcs.write_log(
        logfile,
        'fig_filePyplot = ' + 'r\'' + str(fig_file) + '\'' + '\n' +
        'fig_fileTikZ = ' + 'r\'' + str(tikz_texfile) + '\'' + '\n' +
        '##################################\n')
##    os.system('pdflatex ' + tikz_texfile)
    return 0

def plot_poscar_for_workdir(workdir, euler_angle_type, phi, theta, psi, elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
                            plot_cell_basis_vector_label = False, plot_atom_label = True, poscar_or_contcar = 'POSCAR', fig_format = 'png', fig_dpi = 100,
                            draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical'):

    '''
    - Descriptions
     * Visualization of POSCARs. 
     * The mother folder needs to be specified which contains the folders with POSCARs
     * Euler angles are used to rotate the view of the model
     
    - Args:
     * workdir: String format. The mother folder which contains the folders with POSCARs
     * euler_angle_type: string of length 3. It specify the type of rotations based on Eulerian angles. Choices are 'zyz', 'zxz', 'zyx', etc.. Usually the 'zyz' type is used.

                'zyz' : proper Euler angle, y-convention. Performe consecutive rotations at axes counter-clockwisely. z-y-z rotation.
                        First rotate the z axes of atoms by an angle phi, then rotate the intermidiate y axis of atoms by an angle theta, finally rotate the final z axis of atoms by an angle psi

                'zxz' : proper Euler angle, x-convention. Performe consecutive rotations at axes counter-clockwisely. z-x-z rotation.
                        First rotate the z axes of atoms by an angle phi, then rotate the intermidiate x axis of atoms by an angle theta, finally rotate the final z axis of atoms by an angle psi

                'zyx' : Tait-Bryan angles. z-y-x rotation. Performe consecutive rotations at axes counter-clockwisely. z-y-x rotation.
                        First rotate the z axes of atoms by an angle phi, then rotate the intermidiate y axis of atoms by an angle theta, finally rotate the final x axis of atoms by an angle psi

     * phi, theta, psi: float formats. The first, second, and third rotation Eulerian angles, units in degrees.
     * draw_mirror_atom: Logical value. Whether to plot the mirror atoms at the periodic boundary
     * box_on: Logical value. Whether to plot the box or not
     * axis_indicator: Logic value. Whethether to plot the axis indicator
     * plot_cell_basis_vector_label: Logical value. Whether to plot the cell basis vector labels( i.e., to label the three basis vectors of the cell as a, b, and c)
     * plot_atom_label: Logical value. If true, then plot the atom name of each atom.
     * poscar_or_contcar: String format. Determine whether to plot POSCAR or CONTCAR. Either 'POSCAR' or 'CONTCAR' can be used. 
     * fig_format: String format. Figformat is a string that defines output figure format. Supported fig_format: 'png','eps','pdf','tif','tiff','jpg','jpeg','svg','svgz','pgf','ps','raw','rgba'
     * fig_dpi: float format. The DPI for non-vector graphics.
    '''
    import os
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    workdir = os.path.abspath(workdir)
    if os.path.exists(workdir) == False:
        funcs.write_log(logfile, '## ERROR:' + workdir +
                       ' does not exist.')
    workdir_list = [x[0] for x in os.walk(workdir)]
    
    if poscar_or_contcar != 'POSCAR' and poscar_or_contcar != 'CONTCAR':
        funcs.write_log(logfile,'## Error: The value of the poscar_or_contcar is not correct. Whether POSCAR or CONTCAR will be plotted must be defined.')
    for i in range(0,len(workdir_list)):
        iposcar_file = workdir_list[i] + '/' + poscar_or_contcar
        if os.path.exists(iposcar_file) == False:
            funcs.write_log(logfile, '## WARNING:' + iposcar_file +
                           ' does not exist.')
        else:
            plot_poscar(iposcar_file, euler_angle_type, phi, theta, psi, elmt_color, draw_mirror_atom, box_on, axis_indicator,
                        plot_cell_basis_vector_label, plot_atom_label, fig_format, fig_dpi,
                        draw_colormap, colormap_column_indx, colormap_vmin, colormap_vmax, vmax_color, vmin_color, colorbar_alignment)
    funcs.write_log(
        logfile,
        'elmt_color = ' + str(elmt_color) + '\n' +
        'vasp_plot.plot_poscar_for_workdir(' + '\n' +
        '    workdir=' + 'r\'' + str(workdir) + '\'' + ',\n' +
        '    euler_angle_type=' + '\'' + str(euler_angle_type) + '\'' + ',\n' +
        '    phi=' + str(phi) + ',\n' +
        '    theta=' + str(theta) + ',\n' +
        '    psi=' + str(psi) + ',\n' +
        '    elmt_color=elmt_color' + ',\n' +
        '    draw_mirror_atom=' + str(draw_mirror_atom) + ',\n' +
        '    box_on=' + str(box_on) + ',\n' +
        '    axis_indicator=' + str(axis_indicator) + ',\n' +
        '    plot_cell_basis_vector_label=' + str(plot_cell_basis_vector_label) + ',\n' +
        '    plot_atom_label=' + str(plot_atom_label) + ',\n' +
        '    poscar_or_contcar=' + '\'' + str(poscar_or_contcar) + '\'' + ',\n' +
        '    fig_format=' + '\'' + str(fig_format) + '\'' + ',\n' +
        '    fig_dpi=' + str(fig_dpi) + ',\n' +               
        '    draw_colormap=' + str(draw_colormap) + ',\n' +
        '    colormap_column_indx=' + str(colormap_column_indx) + ',\n' +
        '    colormap_vmin=' + str(colormap_vmin) + ',\n' +
        '    colormap_vmax=' + str(colormap_vmax) + ',\n' +
        '    vmin_color=' + '\'' + str(vmin_color) + '\'' + ',\n' +
        '    vmax_color=' + '\'' + str(vmax_color) + '\'' + ',\n' +
        '    colorbar_alignment=' + '\'' + str(colorbar_alignment) + '\'' + ')\n' +
        '###############################\n')
    return workdir_list


def plot_poscar_for_sysname(sysname_file, euler_angle_type = 'zyx', phi = -3, theta = 4, psi = 0, elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
                            plot_cell_basis_vector_label = False, plot_atom_label = True, fig_format = 'png', fig_dpi = 100,
                            draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical'):
    '''
    Description:
        Visualization of POSCAR models for systems in .sysname file. Euler angles are used to rotate the view of the model
        When lots of models with substituted atoms are generated, it is difficult and time-consuming to plot each model configuration one by one manually.
        The current function is used to plot the generated models for systems in .sysname file
        Special Note: first run ST.substitution() or RE.rep_elmt() to generate models with substituted atoms, then run current module. If not, an error may occur.
        In the same folder as *.sysname file, there must exist a folder with the same name as *.sysname file. For example, in the folder which contains a.sysname,
        there must exist the folder "a" which contains all the job folders containing all the POSCARs
    Args:
        sysname_file: String format. It specifies the system name(file folder) that constains the POSCAR for each sysetm. Each system occupies a line
        view_axis: Decide viewing from which direction: 'x' or 'y' or 'z'
        phi, theta, kappa:Define rotation angles, units in degrees. Reference:J. Hirth and J. Lothe, Theory of dislocations, 1982
        draw_mirror_atom: Logic value. Whether to plot the mirror atoms at the periodic boundary
        box_on: Logic value. Whether to plot the box or not
        axis_indicator: Logic value. Whethether to plot the axis indicator
        plot_cell_basis_vector_label: Logical value. Whether to plot the cell basis vector labels( i.e., to label the three basis vectors of the cell as a, b, and c)
        plot_atom_label: Logical value. If true, then plot the atom name of each atom.
        fig_format: Figformat is a string that defines output figure format. Supported fig_format: 'png','eps','pdf','tif','tiff','jpg','jpeg','svg','svgz','pgf','ps','raw','rgba'
        Fig_dpi: dpi for non-vector graphics
        logfile: Runtime output are written to logfile
    '''
    import os
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']

    (filename, file_extension) = os.path.splitext(sysname_file)
    subst_sysname_file = sysname_file
    with open(subst_sysname_file) as f:
        lines = f.readlines()
        n_lines = len(lines)
    if os.path.exists(filename) == False:
        funcs.write_log(logfile, '## ERROR: sysname file folder ' + str(filename) +
                       ' does not exist. Please check whether you have generated models according to the input .subst file ' + str(filename) + '.subst' + '\n' +
                       'Hints: Please first run ST.substitution() or RE.rep_elmt() to generate models with substituted atoms, then run plot_poscar_for_sysname()' + '\n' +
                       '####################################')
    else:
        for i_line in range(0,n_lines):
            sysname = funcs.split_line(lines[i_line])[0]
            workdir = os.getcwd() + '/' + str(filename) + '/' + sysname
            poscar_file_path = workdir + '/POSCAR'
            plot_poscar(poscar_file_path, euler_angle_type, phi, theta, psi, elmt_color, draw_mirror_atom, box_on, axis_indicator,
                        plot_cell_basis_vector_label, plot_atom_label, fig_format, fig_dpi)
    funcs.write_log(
        logfile,
        'elmt_color = ' + str(elmt_color) + '\n' +
        'vasp_plot.plot_poscar_for_workdir(' + '\n' +
        '    sysname_file=' + 'r\'' + str(sysname_file) + '\'' + ',\n' +
        '    euler_angle_type=' + '\'' + str(euler_angle_type) + '\'' + ',\n' +
        '    phi=' + str(phi) + ',\n' +
        '    theta=' + str(theta) + ',\n' +
        '    psi=' + str(psi) + ',\n' +
        '    elmt_color=elmt_color' + ',\n' +
        '    draw_mirror_atom=' + str(draw_mirror_atom) + ',\n' +
        '    box_on=' + str(box_on) + ',\n' +
        '    axis_indicator=' + str(axis_indicator) + ',\n' +
        '    plot_cell_basis_vector_label=' + str(plot_cell_basis_vector_label) + ',\n' +
        '    plot_atom_label=' + str(plot_atom_label) + ',\n' +
        '    fig_format=' + '\'' + str(fig_format) + '\'' + ',\n' +
        '    fig_dpi=' + str(fig_dpi) + ',\n' +              
        '    draw_colormap=' + str(draw_colormap) + ',\n' +
        '    colormap_column_indx=' + str(colormap_column_indx) + ',\n' +
        '    colormap_vmin=' + str(colormap_vmin) + ',\n' +
        '    colormap_vmax=' + str(colormap_vmax) + ',\n' +
        '    vmin_color=' + '\'' + str(vmin_color) + '\'' + ',\n' +
        '    vmax_color=' + '\'' + str(vmax_color) + '\'' + ',\n' +
        '    colorbar_alignment=' + '\'' + str(colorbar_alignment) + '\'' + ')\n' +
        '###############################\n')
    return 0
