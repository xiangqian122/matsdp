# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use("Agg")

def plot_dos(atom_doscar_file_path_list, atom_sysname_list = None, atom_indx_list = None, atom_palette_list = None, atom_subplot_arg_list = None,
             subplot_arg_list = None, subplot_xlo_list = None, subplot_xhi_list = None, subplot_ylo_list = None, subplot_yhi_list = None,
             subplot_xtick_list = None, subplot_ytick_list = None, subplot_xlabel_list = None, subplot_ylabel_list = None, subplot_share_xy_list = [False, False] , mainplot_axis_label_list = [False, False], xtick_direction = 'out', ytick_direction = 'out',
             dos_mode_dict = None, fermi_shift_zero = True, peak_analyzer = False, peak_analyzer_factor = 0.02, smoothing = False, smoothing_factor = 0.05, line_width = 2.0, font_size = 18, fig_format = 'png', fig_size = [15,10], fig_dpi = 600):
    '''
    - Descriptions
     * Plot PDOS, LDOS, TDOS, now only availabel for LORBIT = 11.
     * There are three types of input arguments: atom related input arguments, subplot related input arguments, and others

    - Atom related args
     * atom_doscar_file_path_list: list format. Conatains DOSCAR files for each atom. The directory of DOSCAR files can either be full path or relative path
     * atom_sysname_list: system name for each atom, it corresponds to the atoms in the atom_doscar_file_path_list. This is for the purpose of labeling the DOS curves in the legend.
               If sysname_list = None, then the label of system name will not shown in the legend
               For example, sysname_list = ['System1','System1','System2']
     * atom_indx_list: list format. Atom index list, it corresponding to the atoms in  atom_doscar_file_path_list. If it is integer type then it denotes the atom index, if it is string type then it denotes the atom name
               atom_indx_list = [1,2,45] denotes the 1st, 2nd, and the 45th atoms in the POSCAR
               atom_indx_list = ['Ni1','Al3','Re3'] denotes Ni1, Al3, and Re3 in the POSCAR
               atom_indx_list = ['TDOS'] and atom_indx_list = [0] donotes the total dos
               atom_indx_list = ['LDOS'] donotes the local dos
     * atom_palette_list: list format. Color for DOS curves of each atom.
     * atom_subplot_arg_list: list format. Defines the DOS curves of the atom are in which subplot. For example, atom_subplot_arg_list = [221, 222] denotes that the DOS curves of the first and the second atoms are in the subplot(221) and subplot(222) subplots, respectively. 

    - Subplot related args
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
     * xtick_direction = 'in' or 'out', if None value is taken, then xtick_direction = 'out'
     * ytick_direction = 'in' or 'out', if None value is taken, then ytick_direction = 'out'
     * smoothing: smoothing of the DOS. Currently only the Lorentzian broadening scheme is availabel. The value can be: True, False, 'Lorentzian'
     * peak_analyzer_factor: recommended value temp_scale >= 0.01 and temp_scale <= 0.06, usually 0.02 is selected. The smaller this value, the more fine peaks can be found.
     * smoothing_factor: float type. This defines the smoothing factor. In the case of the Lorentzian broadening scheme, the smoothing factor is the broadening width
     * line_width = 2.0            # Recommended value 0.5 -- 3.0 (from thin to fat)
     * font_size = 18            # Recommended value 18. This value designate the font size for the axis label font size and the legend font size
    - fig_size: set this value to avoid overlapping of subfigures
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
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)
    periodic_table_dict = periodic_table.periodic_tab()

    #Possible further modifications: check the consistency of lengths of lists: atom_doscar_file_path_list and atom_indx_list; subplot_arg_list and SubplotXYLabelLogicList
    ###############################################################################################################
    #User defined region, change the params when necessary. Usually we don't change these parameters unless needed.
    ###############################################################################################################
    default_palette_list = ['black','magenta','cyan','darkviolet','red','gray','darkorange','darkcyan','palegreen',
                   'goldenrod','lightpink','tan','mediumpurple','steelblue','tomato','mediumturquoise',
                   'mediumslateblue','brown','darkseagreen','wheat','seagreen','maroon','deeppink','chocolate',
                   'bisque','plum']
    default_lty_list = ['-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':']    #Line style for DOS orbital
    upper_spine_scale = 1.02     #Scale the coordinate of the upper axis spine. Recommended value 1.02
    lower_spine_scale = 1.02     #Scale the coordinate of the lower axis spine. Recommended value 1.02
    #scale the DOS's: for purpose of plotting arbitrary unit DOS, default value = 1
    orbit_scaling_dict = {'s':1, 'p':1, 'd':1,'f':1,'LDOS':1}
    #Set the string format of the x and the y axis. for example xaxis_major_format = '%.2f'. These parameters are useful for the generalized format of the axes of the subplots
    xaxis_major_format = None
    yaxis_major_format = None
    ##################################
    #Codes maintained by the developer 
    ##################################

    def label_peaks(axes_i, energy, sample_dos_arr, subplot_indx, subplot_dict, peak_analyzer_factor):
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
                temp_delta = 0.2  # recommended value temp_data = 0.20
                if np.abs((a1 - a2)/a1) > temp_delta:
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
        # peak_analyzer_factor: recommended value temp_scale >= 0.01 and temp_scale <= 0.06, usually 0.02 is selected. The smaller this value, the more fine peaks can be found.
        characteristic_states_num_sample_dos_arr_1d = peak_analyzer_factor * np.max(np.abs(sample_dos_arr_1d))
        noise_baseline_sample_dos_arr = find_noise_baseline(characteristic_states_num_sample_dos_arr_1d, np.abs(sample_dos_arr_1d))
        peaks, properties = find_peaks(np.abs(sample_dos_arr_1d), height = noise_baseline_sample_dos_arr, prominence = characteristic_states_num_sample_dos_arr_1d)
        for ipeak in range(0,len(peaks)):
            if energy[peaks][ipeak] >= subplot_dict['xlo_list'][subplot_indx] and energy[peaks][ipeak] <= subplot_dict['xhi_list'][subplot_indx]:
                axes_i.text(energy[peaks][ipeak], sample_dos_arr_1d[peaks][ipeak], str(energy[peaks][ipeak]),color='red',size=4)
        axes_i.plot(energy[peaks], sample_dos_arr_1d[peaks], "x",color='red')
        return energy[peaks]

    # initialization
    for i in range(0,len(atom_doscar_file_path_list)):
        i_atom_doscar_file_path = os.path.abspath(str(atom_doscar_file_path_list[i]))
        if os.path.exists(i_atom_doscar_file_path) and os.path.getsize(i_atom_doscar_file_path) > 0:
            pass
        else:
            print('ERROR: vasp_plot error. The file ' + str(i_atom_doscar_file_path) + ' does not exist or is empty')
            funcs.write_log(logfile, '#' + i_atom_doscar_file_path + " doesn't exist or is empty")
            quit()
        atom_doscar_file_path_list[i] = i_atom_doscar_file_path
    if atom_sysname_list in [None,'None','none']:
        atom_sysname_list = [''] * len(atom_doscar_file_path_list)
    else:
        for i in range(0,len(atom_sysname_list)):
            if atom_sysname_list[i] == None or atom_sysname_list[i] == 'None' or atom_sysname_list[i] == 'none' or atom_sysname_list[i] == '':
                atom_sysname_list[i] = ''
            else:
                atom_sysname_list[i] = str(atom_sysname_list[i])
    if atom_indx_list in [None,'None','none']:
        atom_indx_list = [0] * len(atom_doscar_file_path_list)

    if atom_palette_list in [None,'None','none']:
        atom_palette_list = default_palette_list[0:len(atom_doscar_file_path_list)]

    if atom_subplot_arg_list in [None,'None','none']:
        atom_subplot_arg_list = ['1,1,1'] * len(atom_doscar_file_path_list)


    for i_atom_subplot in range(0,len(atom_subplot_arg_list)):
        #Single argument to subplot must be a 3-digit integer, we convert all the subplot arguent to the format of 'xx,xx,xx'
        if type(atom_subplot_arg_list[i_atom_subplot]) == str and ',' in atom_subplot_arg_list[i_atom_subplot]:
            #for example '4,4,10'
            pass
        if type(atom_subplot_arg_list[i_atom_subplot]) == str and ',' not in atom_subplot_arg_list[i_atom_subplot]:
            #for example '442'
            temp_str = atom_subplot_arg_list[i_atom_subplot]
            atom_subplot_arg_list[i_atom_subplot] = ','.join([temp_str[0],temp_str[1],temp_str[2]])
        else:
            #for example 442
            temp_str = str(atom_subplot_arg_list[i_atom_subplot])
            atom_subplot_arg_list[i_atom_subplot] = ','.join([temp_str[0],temp_str[1],temp_str[2]])

    subplot_arg_unique_list = list(set(atom_subplot_arg_list))
    if subplot_arg_list in [None,'None','none']:
        subplot_arg_list = ['1,1,1'] * len(subplot_arg_unique_list) 

    for i_subplot in range(0,len(subplot_arg_list)):
        #Single argument to subplot must be a 3-digit integer, we convert all the subplot arguent to the format of 'xx,xx,xx'
        if type(subplot_arg_list[i_subplot]) == str and ',' in subplot_arg_list[i_subplot]:
            #for example '4,4,10'
            pass
        if type(subplot_arg_list[i_subplot]) == str and ',' not in subplot_arg_list[i_subplot]:
            #for example '442'
            temp_str = subplot_arg_list[i_subplot]
            subplot_arg_list[i_subplot] = ','.join([temp_str[0],temp_str[1],temp_str[2]])
        else:
            #for example 442
            temp_str = str(subplot_arg_list[i_subplot])
            subplot_arg_list[i_subplot] = ','.join([temp_str[0],temp_str[1],temp_str[2]])

    if subplot_xlo_list in [None,'None','none']:
        subplot_xlo_list = [None] * len(subplot_arg_unique_list) 
    if subplot_xhi_list in [None,'None','none']:
        subplot_xhi_list = [None] * len(subplot_arg_unique_list) 
    if subplot_ylo_list in [None,'None','none']:
        subplot_ylo_list = [None] * len(subplot_arg_unique_list) 
    if subplot_yhi_list in [None,'None','none']:
        subplot_yhi_list = [None] * len(subplot_arg_unique_list) 
    if subplot_xtick_list in [None,'None','none']:
        subplot_xtick_list = [False] * len(subplot_arg_unique_list) 
    if subplot_ytick_list in [None,'None','none']:
        subplot_ytick_list = [False] * len(subplot_arg_unique_list) 
    if subplot_xlabel_list in [None,'None','none']:
        subplot_xlabel_list = [False] * len(subplot_arg_unique_list) 
    if subplot_ylabel_list in [None,'None','none']:
        subplot_ylabel_list = [False] * len(subplot_arg_unique_list) 
    if dos_mode_dict in [None,'None','none']:
        dos_mode_dict = periodic_table_dict['dos_mode']

    #Initialize the atom dictionary
    atom_dict = {}
    atom_dict['doscar_file_path_list'] = atom_doscar_file_path_list 
    atom_dict['sysname_list'] = atom_sysname_list
    atom_dict['indx_list'] = atom_indx_list
    atom_dict['palette_list'] = atom_palette_list
    atom_dict['subplot_arg_list'] = atom_subplot_arg_list

    #Initialize Subplot dictionary
    subplot_dict = {}
    subplot_dict['arg_list'] = subplot_arg_unique_list
    subplot_dict['indx_list'] = list(range(len(subplot_arg_unique_list)))
    subplot_dict['xlo_list'] = [999999.9999] * len(subplot_arg_unique_list)
    subplot_dict['xhi_list'] = [-999999.9999] * len(subplot_arg_unique_list)
    subplot_dict['ylo_list'] = [999999.9999] * len(subplot_arg_unique_list)
    subplot_dict['yhi_list'] = [-999999.9999] * len(subplot_arg_unique_list)
    subplot_dict['xtick_list'] = [False] * len(subplot_arg_unique_list)
    subplot_dict['ytick_list'] = [False] * len(subplot_arg_unique_list)
    subplot_dict['xlabel_list'] = [None] * len(subplot_arg_unique_list)
    subplot_dict['ylabel_list'] = [None] * len(subplot_arg_unique_list)
    subplot_dict['share_xy'] = [False,False]

    #Initialize the gereral paramters list
    general_params_dict = {}
    general_params_dict['mainplot_axis_label_list'] = mainplot_axis_label_list 
    general_params_dict['xtick_direction'] = xtick_direction
    general_params_dict['ytick_direction'] = ytick_direction
    general_params_dict['fermi_shift_zero'] = fermi_shift_zero 
    general_params_dict['peak_analyzer'] = peak_analyzer
    general_params_dict['peak_analyzer_factor'] = peak_analyzer_factor
    general_params_dict['smoothing'] = smoothing
    general_params_dict['smoothing_factor'] = smoothing_factor
    general_params_dict['line_width'] = line_width
    general_params_dict['font_size'] = font_size
    general_params_dict['fig_format'] = fig_format
    general_params_dict['fig_size'] = fig_size
    general_params_dict['fig_dpi'] = fig_dpi 

    #Extract data from the user input into the subplot_dict
    subplot_dict['share_xy'] = subplot_share_xy_list
    for i_subplot in range(0,len(subplot_arg_unique_list)):
        if subplot_arg_list[i_subplot] in subplot_dict['arg_list']:
            subplot_dict_indx = subplot_dict['arg_list'].index(subplot_arg_list[i_subplot])
            if subplot_xlo_list[i_subplot] not in [None,'None','none']:
                subplot_dict['xlo_list'][subplot_dict_indx] = subplot_xlo_list[i_subplot]
            if subplot_xhi_list[i_subplot] not in [None,'None','none']:
                subplot_dict['xhi_list'][subplot_dict_indx] = subplot_xhi_list[i_subplot]
            if subplot_ylo_list[i_subplot] not in [None,'None','none']:
                subplot_dict['ylo_list'][subplot_dict_indx] = subplot_ylo_list[i_subplot]
            if subplot_yhi_list[i_subplot] not in [None,'None','none']:
                subplot_dict['yhi_list'][subplot_dict_indx] = subplot_yhi_list[i_subplot]
            subplot_dict['xtick_list'][subplot_dict_indx] = subplot_xtick_list[i_subplot]
            subplot_dict['ytick_list'][subplot_dict_indx] = subplot_ytick_list[i_subplot]
            subplot_dict['xlabel_list'][subplot_dict_indx] = subplot_xlabel_list[i_subplot]
            subplot_dict['ylabel_list'][subplot_dict_indx] = subplot_ylabel_list[i_subplot]
    if general_params_dict['xtick_direction'] in [None,'None','none']:
        general_params_dict['xtick_direction'] = 'out'
    if general_params_dict['ytick_direction'] in [None,'None','none']:
        general_params_dict['ytick_direction'] = 'out'
    subplot_arg_unique_count_list = []
    counter_list = []
    for i in range(len(subplot_dict['arg_list'])):
        subplot_arg_unique_count_list.append(atom_dict['subplot_arg_list'].count(subplot_arg_unique_list[i]))
        counter_list.append(0)
    sysname_list = []
    dos_atomname_list = []
    e_fermi_list = []
    initial_peak_analyzer_factor = peak_analyzer_factor
    if general_params_dict['peak_analyzer'] not in [True, False]:
        print('ERROR: vasp_plot input value error. The value of peak_analyzer should be logic value True or False')
        funcs.write_log(logfile, 'ERROR: vasp_plot input value error. The value of peak_analyzer should be logic value True or False')
        quit()
    if general_params_dict['peak_analyzer_factor'] in [None,'None','none']:
        general_params_dict['peak_analyzer_factor'] = 0.02
    initial_line_width = general_params_dict['line_width']
    if general_params_dict['line_width'] in [None,'None','none']:
        general_params_dict['line_width'] = 1.0
    initial_font_size = general_params_dict['font_size']
    if general_params_dict['font_size'] in [None,'None','none']:
        general_params_dict['font_size'] = 18
    plt.rcParams['figure.figsize'] = (general_params_dict['fig_size'][0],general_params_dict['fig_size'][1])
    # make use of the 'golden ratio' when treating the sizes of the x-ticks or y-ticks
    golden_ratio = 0.618
    axes_line_width = 2.0
    plt.rcParams['xtick.direction'] = general_params_dict['xtick_direction']
    plt.rcParams['ytick.direction'] = general_params_dict['ytick_direction']
    plt.rcParams['xtick.major.width'] = axes_line_width
    plt.rcParams['ytick.major.width'] = axes_line_width
    plt.rcParams['xtick.minor.width'] = axes_line_width * golden_ratio
    plt.rcParams['ytick.minor.width'] = axes_line_width * golden_ratio
    plt.rcParams['xtick.major.size'] = general_params_dict['font_size'] * (1 - golden_ratio)
    plt.rcParams['ytick.major.size'] = general_params_dict['font_size'] * (1 - golden_ratio)
    plt.rcParams['xtick.minor.size'] = general_params_dict['font_size'] * (1 - golden_ratio) * golden_ratio 
    plt.rcParams['ytick.minor.size'] = general_params_dict['font_size'] * (1 - golden_ratio) * golden_ratio
    plt.rcParams['axes.linewidth'] = axes_line_width
    if len(atom_dict['doscar_file_path_list']) == 1:
        doscar_file_path = atom_dict['doscar_file_path_list'][0]
        workdir, dos_file = funcs.file_path_name(doscar_file_path)
        if isinstance(atom_dict['indx_list'][0],str) and atom_dict['indx_list'][0] != 'TDOS':
            atom_indx = convert.atomname2indx(os.path.join(workdir, 'POSCAR'),atom_dict['indx_list'][0])
        elif isinstance(atom_dict['indx_list'][0],str) and atom_dict['indx_list'][0] == 'TDOS':
            atom_indx = 0
        elif isinstance(atom_dict['indx_list'][0],int):
            atom_indx = atom_dict['indx_list'][0]
        fig_dos = plt.figure('fig' + dos_file)
    else:
        fig_dos = plt.figure('fig_dos_cmpr')
    axes_list = [None]
    axes_list[0] = plt.subplot(1,1,1)
    active_axes = axes_list[0]
    axes_arg_list = [None]
    axes_arg_list[0] = '1,1,1'
    #Below is the axis setting for the main figure, and the axis settings only work for main figure with multiple subplots(number of subplots >=2):
    #set interspacing between subplots
    if subplot_dict['share_xy'][0] == True:
        hspace_val = 0
    else:
        hspace_val = None
    if subplot_dict['share_xy'][1] == True:
        wspace_val = 0
    else:
        wspace_val = None
    plt.subplots_adjust(bottom =0.095, left = 0.090, top = 0.99, right = 0.98, wspace = wspace_val, hspace = hspace_val)
    subplot_last_arg = atom_subplot_arg_list[len(atom_dict['subplot_arg_list'])-1]
##    mult_plot_logic = (isinstance(subplot_last_arg,int) and subplot_last_arg != 111) or (isinstance(subplot_last_arg,list) and subplot_last_arg != [1,1,1])
    mult_plot_logic = (subplot_last_arg not in [111,'1,1,1']) or (isinstance(subplot_last_arg,list) and subplot_last_arg != [1,1,1])
    xlabel_str = 'Energy (eV)'
    dos_mode_dict_str_list = []
    for item in dos_mode_dict.values():
        dos_mode_dict_str_list = dos_mode_dict_str_list + item
    if '0' in atom_dict['indx_list'] or 'TDOS' in atom_dict['indx_list']:
        ylabel_str = 'TDOS (arb. units)'
    else:
        if 'LDOS' not in dos_mode_dict_str_list:
            ylabel_str = 'PDOS (arb. units)'
        elif 'LDOS' in dos_mode_dict_str_list and len(set(dos_mode_dict_str_list)) == 1:
            ylabel_str = 'LDOS (arb. units)'
        elif 'LDOS' in dos_mode_dict_str_list and len(set(dos_mode_dict_str_list)) > 1:
            ylabel_str = 'DOS (arb. units)'


    def plot_lines(axes_i, orbits, doscar_dict, energy, atom_dict, dos_lty, subplot_indx, subplot_dict, dos_file_indx, i_elmt_name, orbit_scaling_dict, general_params_dict):
        from .. import funcs
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
            if len(atom_dict['sysname_list'][dos_file_indx].strip()) == 0:
                label_text = i_elmt_name + orbit_label_list[orbitname_indx]
            elif ',' == atom_dict['sysname_list'][dos_file_indx].strip()[-1]:
                label_text = atom_dict['sysname_list'][dos_file_indx].strip() + ' ' + i_elmt_name + orbit_label_list[orbitname_indx]
            else:
                label_text = atom_dict['sysname_list'][dos_file_indx].strip() + ', ' + i_elmt_name + orbit_label_list[orbitname_indx]
            if orbitname_list[orbitname_indx] in orbits:
                orbit_scaling_factor = orbit_scaling_dict[orbitname_list[orbitname_indx]]
                if doscar_dict['num_col'] == 5 or doscar_dict['num_col'] == 19 or doscar_dict['num_col'] == 33:
                    i_dos_arr_original_up = doscar_dict[orbitname_list[orbitname_indx]+'_up']
                    i_dos_arr_original_dw = doscar_dict[orbitname_list[orbitname_indx]+'_dw']
                    if general_params_dict['smoothing'] == True or general_params_dict['smoothing'] == 'Lorentzian':
                        i_energy_arr, i_dos_arr_up = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original_up, delta = general_params_dict['smoothing_factor']) 
                        i_energy_arr, i_dos_arr_dw = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original_dw, delta = general_params_dict['smoothing_factor'])
                    else:
                        i_energy_arr = energy
                        i_dos_arr_up = i_dos_arr_original_up 
                        i_dos_arr_dw = i_dos_arr_original_dw
                        
                    axes_i.plot(i_energy_arr, i_dos_arr_up * orbit_scaling_factor, color = atom_dict['palette_list'][dos_file_indx], linestyle = dos_lty[dos_lty_indx],
                             marker = '', linewidth = general_params_dict['line_width'], label = label_text)
                    axes_i.plot(i_energy_arr, i_dos_arr_dw * orbit_scaling_factor, color = atom_dict['palette_list'][dos_file_indx], linestyle = dos_lty[dos_lty_indx],
                             marker = '', linewidth = general_params_dict['line_width'])
                    dos_lty_indx += 1                
                    if general_params_dict['peak_analyzer'] == True:
                        label_peaks(axes_i, i_energy_arr, i_dos_arr_up * orbit_scaling_factor, subplot_indx, subplot_dict, general_params_dict['peak_analyzer_factor'])
                        label_peaks(axes_i, i_energy_arr, i_dos_arr_dw * orbit_scaling_factor, subplot_indx, subplot_dict, general_params_dict['peak_analyzer_factor'])
                elif doscar_dict['num_col'] == 3 or doscar_dict['num_col'] == 10 or doscar_dict['num_col'] == 17:
                    i_dos_arr_original = doscar_dict[orbitname_list[orbitname_indx]]
                    if general_params_dict['smoothing'] == True or general_params_dict['smoothing'] == 'Lorentzian':
                        i_energy_arr, i_dos_arr = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original, delta = general_params_dict['smoothing_factor']) 
                    else:
                        i_energy_arr = energy
                        i_dos_arr = i_dos_arr_original

                    axes_i.plot(i_energy_arr, i_dos_arr * orbit_scaling_factor, color = atom_dict['palette_list'][dos_file_indx], linestyle = dos_lty[dos_lty_indx],
                             marker = '', linewidth = general_params_dict['line_width'], label = label_text)
                    dos_lty_indx += 1                
                    if general_params_dict['peak_analyzer'] == True:
                        label_peaks(axes_i, i_energy_arr, i_dos_arr * orbit_scaling_factor, subplot_indx, subplot_dict, general_params_dict['peak_analyzer_factor'])
        return 0 

    #Below plot DOS in each subplot
    for dos_file_indx in range(0,len(atom_dict['doscar_file_path_list'])):
        doscar_file_path = atom_dict['doscar_file_path_list'][dos_file_indx]
        #print(doscar_file_path)  # this is for the purpose of debugging
        workdir, dos_file = funcs.file_path_name(doscar_file_path)
        if isinstance(atom_dict['indx_list'][dos_file_indx], str) and atom_dict['indx_list'][dos_file_indx] != 'TDOS' and '+' not in atom_dict['indx_list'][dos_file_indx]:
            atom_indx = convert.atomname2indx(os.path.join(workdir, 'POSCAR'), atom_dict['indx_list'][dos_file_indx])
        elif isinstance(atom_dict['indx_list'][dos_file_indx], str) and atom_dict['indx_list'][dos_file_indx] != 'TDOS' and '+' in atom_dict['indx_list'][dos_file_indx]:
            # This is for the case of LDOS of multiple atoms. For example, 'Ni123+Ni124+Re2', or '23+33+1'. In this case, the atom_indx must be string type.
            atom_indx = [x.strip() for x in atom_dict['indx_list'][dos_file_indx].split('+')]
            for i_temp in range(0, len(atom_indx)):
                if isinstance(atom_indx[i_temp], str):
                    atom_indx[i_temp] = convert.atomname2indx(os.path.join(workdir, 'POSCAR'), atom_indx[i_temp]) 
                elif isinstance(atom_indx[i_temp], int):
                    atom_indx[i_temp] = int(atom_indx[i_temp])
        elif isinstance(atom_dict['indx_list'][dos_file_indx], str) and atom_dict['indx_list'][dos_file_indx] == 'TDOS':
            atom_indx = 0
        elif isinstance(atom_dict['indx_list'][dos_file_indx], int):
            atom_indx = atom_dict['indx_list'][dos_file_indx]
        doscar_dict = vasp_read.read_doscar(doscar_file_path, atom_indx, False)
        outcar_file_path = os.path.join(workdir, 'OUTCAR')
        doscar_file_path = os.path.join(workdir, 'DOSCAR')
        incar_file_path = os.path.join(workdir, 'INCAR')
        poscar_file_path = os.path.join(workdir, 'POSCAR')
        poscar_dict = vasp_read.read_poscar(poscar_file_path)
        outcar_params_dict = vasp_read.read_outcar(outcar_file_path)
        LORBIT = outcar_params_dict['LORBIT']
        e_fermi = outcar_params_dict['e_fermi']
        e_fermi_mod = e_fermi + outcar_params_dict['alpha+bet']
        subplot_indx = subplot_dict['arg_list'].index(atom_dict['subplot_arg_list'][dos_file_indx])
        if atom_indx == 0:
            #TDOS
            '''If the subplot's parameters coincide with the previous subplot's parameters, a MatplotlibDeprecationWarning will occur. So we take the following action:
               If new added subplot's parameters is the same as the previous subplot's parameters, we don't add the new subplot, we use the previous subplot instead
               If new added subplot's parameters is different from the previous subplot's parameters, we add the new subplot.'''
            if dos_file_indx == 0:
                previous_subplot_arg = '1,1,1'
            else:
                previous_subplot_arg = atom_dict['subplot_arg_list'][dos_file_indx - 1]
            if atom_dict['subplot_arg_list'][dos_file_indx] == previous_subplot_arg:
                pass
            else:
                temp_arg = atom_dict['subplot_arg_list'][dos_file_indx]
                if temp_arg not in axes_arg_list:
                    axes_list.append(plt.subplot(int(temp_arg.split(',')[0]), int(temp_arg.split(',')[1]), int(temp_arg.split(',')[2])))
                    axes_arg_list.append(temp_arg)
                active_axes = axes_list[axes_arg_list.index(temp_arg)]

            counter_list[subplot_indx] += 1
            sysname_list.append(''.join(str(poscar_dict['elmt_species_arr'][i]) + str(poscar_dict['elmt_num_arr'][i]) for i in range(len(poscar_dict['elmt_species_arr'])-1,-1,-1)))
            dos_atomname_list.append('TDOS')
            e_fermi_list.append(e_fermi_mod)
            energy = doscar_dict['energy'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0)
            #Subplot x axis limit
            if subplot_xlo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                subplot_dict['xlo_list'][subplot_indx] = min(subplot_dict['xlo_list'][subplot_indx],energy.min())
            if subplot_xhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                subplot_dict['xhi_list'][subplot_indx] = max(subplot_dict['xhi_list'][subplot_indx],energy.max())
            if subplot_xlo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_xhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_xlo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] >= subplot_xhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])]:
                subplot_dict['xlo_list'][subplot_indx] = min(subplot_dict['xlo_list'][subplot_indx],energy.min())
                subplot_dict['xhi_list'][subplot_indx] = max(subplot_dict['xhi_list'][subplot_indx],energy.max())
            # find indices of the the lowest and highest values in the energy spectrum
            delta_xlo_minus_energy_list = abs(subplot_dict['xlo_list'][subplot_indx] - energy[:,0])
            delta_xhi_minus_energy_list = abs(subplot_dict['xhi_list'][subplot_indx] - energy[:,0])
            xlo_energy_indx = np.argwhere(delta_xlo_minus_energy_list == min(delta_xlo_minus_energy_list)).squeeze()
            xhi_energy_indx = np.argwhere(delta_xhi_minus_energy_list == min(delta_xhi_minus_energy_list)).squeeze()
            # lower index = (the index of the value closest to xlo) - 1  (except the lowest index in the energy spectrum)
            # higher index = (the index of the value closest to xhi) + 1 (except the highest index in the energy spectrum)
            if xlo_energy_indx != 0:
                xlo_energy_indx = xlo_energy_indx - 1
            if xhi_energy_indx != (len(energy[:,0]) - 1):
                xhi_energy_indx = xhi_energy_indx + 1

            if doscar_dict['num_col'] == 3:
                i_dos_arr_original = doscar_dict['TDOS']
                if general_params_dict['smoothing'] == True or general_params_dict['smoothing'] == 'Lorentzian':
                    i_energy_arr, i_dos_arr = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original, delta = general_params_dict['smoothing_factor'])
                else:
                    i_energy_arr = energy
                    i_dos_arr = i_dos_arr_original
                active_axes.plot(i_energy_arr, i_dos_arr, color = atom_dict['palette_list'][dos_file_indx],marker='',linewidth = general_params_dict['line_width'],
                         label = atom_dict['sysname_list'][dos_file_indx]+'TDOS')
                if peak_analyzer == True:  
                    label_peaks(axes_i, energy, i_dos_arr, subplot_indx, subplot_dict, general_params_dict['peak_analyzer_factor'])
                #Subplot y axis limit
                if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                    subplot_dict['ylo_list'][subplot_indx] = min(subplot_dict['ylo_list'][subplot_indx], min(i_dos_arr[xlo_energy_indx:xhi_energy_indx]))
                if subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                    subplot_dict['yhi_list'][subplot_indx] = max(subplot_dict['yhi_list'][subplot_indx], max(i_dos_arr[xlo_energy_indx:xhi_energy_indx]))
                if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] >= subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])]:
                    subplot_dict['ylo_list'][subplot_indx] = min(subplot_dict['ylo_list'][subplot_indx], min(i_dos_arr[xlo_energy_indx:xhi_energy_indx]))
                    subplot_dict['yhi_list'][subplot_indx] = max(subplot_dict['yhi_list'][subplot_indx], max(i_dos_arr[xlo_energy_indx:xhi_energy_indx]))
       
            elif doscar_dict['num_col'] == 5:
                i_dos_arr_original_up = doscar_dict['TDOS_up']
                i_dos_arr_original_dw = doscar_dict['TDOS_dw']
                if general_params_dict['smoothing'] == True or general_params_dict['smoothing'] == 'Lorentzian':
                    i_energy_arr, i_dos_arr_up = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original_up, delta = general_params_dict['smoothing_factor'])
                    i_energy_arr, i_dos_arr_dw = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original_dw, delta = general_params_dict['smoothing_factor'])
                else:
                    i_energy_arr = energy
                    i_dos_arr_up = i_dos_arr_original_up
                    i_dos_arr_dw = i_dos_arr_original_dw
                active_axes.plot(i_energy_arr,i_dos_arr_up, color = atom_dict['palette_list'][dos_file_indx], marker = '', linewidth = general_params_dict['line_width'],
                         label = atom_dict['sysname_list'][dos_file_indx] + 'TDOS')
                active_axes.plot(i_energy_arr,i_dos_arr_dw, color = atom_dict['palette_list'][dos_file_indx], marker = '', linewidth = general_params_dict['line_width'])
                if general_params_dict['peak_analyzer'] == True:  
                    label_peaks(axes_i, energy, i_dos_arr_up, subplot_indx, subplot_dict, general_params_dict['peak_analyzer_factor'])
                    label_peaks(axes_i, energy, i_dos_arr_dw, subplot_indx, subplot_dict, general_params_dict['peak_analyzer_factor'])                  
                #Subplot y axis limit
                if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                    subplot_dict['ylo_list'][subplot_indx] = min(subplot_dict['ylo_list'][subplot_indx], min(i_dos_arr_up[xlo_energy_indx:xhi_energy_indx]), min(i_dos_arr_dw[xlo_energy_indx:xhi_energy_indx]))
                if subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                    subplot_dict['yhi_list'][subplot_indx] = max(subplot_dict['yhi_list'][subplot_indx], max(i_dos_arr_up[xlo_energy_indx:xhi_energy_indx]), min(i_dos_arr_dw[xlo_energy_indx:xhi_energy_indx]))
                if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] >= subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])]:
                    subplot_dict['ylo_list'][subplot_indx] = min(subplot_dict['ylo_list'][subplot_indx], min(i_dos_arr_up[xlo_energy_indx:xhi_energy_indx]), min(i_dos_arr_dw[xlo_energy_indx:xhi_energy_indx]))
                    subplot_dict['yhi_list'][subplot_indx] = max(subplot_dict['yhi_list'][subplot_indx], max(i_dos_arr_up[xlo_energy_indx:xhi_energy_indx]), min(i_dos_arr_dw[xlo_energy_indx:xhi_energy_indx]))

            #Subplot y axis limit
            subplot_dict['ylo_list'][subplot_indx] = subplot_dict['ylo_list'][subplot_indx] * lower_spine_scale
            subplot_dict['yhi_list'][subplot_indx] = subplot_dict['yhi_list'][subplot_indx] * upper_spine_scale
        elif atom_indx != 0:
            #PDOS and LDOS
            if dos_file_indx == 0:
                previous_subplot_arg = '1,1,1'
            else:
                previous_subplot_arg = atom_dict['subplot_arg_list'][dos_file_indx - 1]
            if atom_dict['subplot_arg_list'][dos_file_indx] == previous_subplot_arg:
                pass
            else:
                temp_arg = atom_dict['subplot_arg_list'][dos_file_indx]
                if temp_arg not in axes_arg_list:
                    axes_list.append(plt.subplot(int(temp_arg.split(',')[0]), int(temp_arg.split(',')[1]), int(temp_arg.split(',')[2])))
                    axes_arg_list.append(temp_arg)
                active_axes = axes_list[axes_arg_list.index(temp_arg)]

            counter_list[subplot_indx] += 1
            if isinstance(atom_indx, int):
                i_elmt_name = poscar_dict['atomname_list'][atom_indx - 1]
            elif isinstance(atom_indx, list):
                # LDOS for multiple atoms
                #naming convension 1
                #i_elmt_name = ';'.join([poscar_dict['atomname_list'][x - 1] for x in atom_indx])
                #naming convension 2
                #i_elmt_name = ';'.join(list(set([poscar_dict['atom_species_arr'][x - 1] for x in atom_indx])))
                #naming convension 3
                #unique_elmt_name_temp_list = list(set([poscar_dict['atom_species_arr'][x - 1] for x in atom_indx]))
                #unique_elmt_count_temp_list = [0.000000 for i in unique_elmt_name_temp_list]
                #for i_unique_elmt in range(0, len(unique_elmt_name_temp_list)):
                #    unique_elmt_count_temp_list[i_unique_elmt] = unique_elmt_name_temp_list.count(unique_elmt_name_temp_list[i_unique_elmt])
                #i_elmt_name = ''.join([str(unique_elmt_count_temp_list[i_unique_elmt]) + unique_elmt_name_temp_list[i_unique_elmt] for i_unique_elmt in range(0, len(unique_elmt_name_temp_list))])
                #naming convension 4
                i_elmt_name = 'G' + str(dos_file_indx + 1) 
            sysname_list.append(''.join(str(poscar_dict['elmt_species_arr'][i])+str(poscar_dict['elmt_num_arr'][i]) for i in range(len(poscar_dict['elmt_species_arr'])-1,-1,-1)))
            if isinstance(atom_indx, int):
                dos_atomname_list.append(i_elmt_name)
            elif isinstance(atom_indx, list):
                #LDOS for multiple atoms
                dos_atomname_list.append(i_elmt_name)

            e_fermi_list.append(e_fermi_mod)
            energy = doscar_dict['energy']  + outcar_params_dict['alpha+bet'] - e_fermi_mod*funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0)
            if isinstance(atom_indx, int):
                orbits = dos_mode_dict[str(poscar_dict['atom_species_arr'][atom_indx - 1])]
            elif isinstance(atom_indx, list):
                #LDOS for multiple atoms
                orbits = 'LDOS'
            if len(orbits) == 0:
                print('## Error: The value of the dos_mode_dict is None. The dos_mode_dict must be defined!')                
            orbitname_list = ['f','d','dxy','dyz','dz2','dxz','dx2','p','py','pz','px','s','LDOS']
            #Subplot x axis limit
            if subplot_xlo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                subplot_dict['xlo_list'][subplot_indx] = min(subplot_dict['xlo_list'][subplot_indx],energy.min())
            if subplot_xhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                subplot_dict['xhi_list'][subplot_indx] = max(subplot_dict['xhi_list'][subplot_indx],energy.max())
            if subplot_xlo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_xhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_xlo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] >= subplot_xhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])]:
                subplot_dict['xlo_list'][subplot_indx] = min(subplot_dict['xlo_list'][subplot_indx],energy.min())
                subplot_dict['xhi_list'][subplot_indx] = max(subplot_dict['xhi_list'][subplot_indx],energy.max())
            # find indices of the the lowest and highest values in the energy spectrum
            delta_xlo_minus_energy_list = abs(subplot_dict['xlo_list'][subplot_indx] - energy[:,0])
            delta_xhi_minus_energy_list = abs(subplot_dict['xhi_list'][subplot_indx] - energy[:,0])
            xlo_energy_indx = np.argwhere(delta_xlo_minus_energy_list == min(delta_xlo_minus_energy_list)).squeeze()
            xhi_energy_indx = np.argwhere(delta_xhi_minus_energy_list == min(delta_xhi_minus_energy_list)).squeeze()
            # lower index = (the index of the value closest to xlo) - 1  (except the lowest index in the energy spectrum)
            # higher index = (the index of the value closest to xhi) + 1 (except the highest index in the energy spectrum)
            if xlo_energy_indx != 0:
                xlo_energy_indx = xlo_energy_indx - 1
            if xhi_energy_indx != (len(energy[:,0]) - 1):
                xhi_energy_indx = xhi_energy_indx + 1
            #Subplot y axis limit
            #if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none'] or subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
            if doscar_dict['num_col'] == 10 or doscar_dict['num_col'] == 17:
                #Spin restricted calculation
                for orbit_indx in range(0,len(orbitname_list)):
                    if orbitname_list[orbit_indx] in orbits:
                        i_dos_arr_original = doscar_dict[orbitname_list[orbit_indx]]
                        if general_params_dict['smoothing'] == True or general_params_dict['smoothing'] == 'Lorentzian':
                            i_energy_arr, i_dos_arr = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original, delta = general_params_dict['smoothing_factor'])
                        else:
                            i_energy_arr = energy
                            i_dos_arr = i_dos_arr_original
                        if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                            subplot_dict['ylo_list'][subplot_indx] = min(subplot_dict['ylo_list'][subplot_indx], min(i_dos_arr[xlo_energy_indx:xhi_energy_indx]))
                        if subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                            subplot_dict['yhi_list'][subplot_indx] = max(subplot_dict['yhi_list'][subplot_indx], max(i_dos_arr[xlo_energy_indx:xhi_energy_indx]))
                        if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] >= subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])]:
                            subplot_dict['ylo_list'][subplot_indx] = min(subplot_dict['ylo_list'][subplot_indx], min(i_dos_arr[xlo_energy_indx:xhi_energy_indx]))
                            subplot_dict['yhi_list'][subplot_indx] = max(subplot_dict['yhi_list'][subplot_indx], max(i_dos_arr[xlo_energy_indx:xhi_energy_indx]))
            elif doscar_dict['num_col'] == 19 or doscar_dict['num_col'] == 33:
                #Spin polarized calculation
                for orbit_indx in range(0,len(orbitname_list)):
                    if orbitname_list[orbit_indx] in orbits:
                        i_dos_arr_original_up = doscar_dict[orbitname_list[orbit_indx]+'_up'] 
                        i_dos_arr_original_dw = doscar_dict[orbitname_list[orbit_indx]+'_dw']
                        if smoothing == True or smoothing == 'Lorentzian':
                            i_energy_arr, i_dos_arr_up = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original_up, delta = general_params_dict['smoothing_factor'])
                            i_energy_arr, i_dos_arr_dw = funcs.lorentzian_broadening(x_arr = energy, y_arr = i_dos_arr_original_dw, delta = general_params_dict['smoothing_factor'])
                        else:
                            i_energy_arr = energy
                            i_dos_arr_up = i_dos_arr_original_up
                            i_dos_arr_dw = i_dos_arr_original_dw
                        if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                            subplot_dict['ylo_list'][subplot_indx] = min(subplot_dict['ylo_list'][subplot_indx], min(i_dos_arr_up[xlo_energy_indx:xhi_energy_indx]), min(i_dos_arr_dw[xlo_energy_indx:xhi_energy_indx]))
                        if subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] in [None,'None','none']:
                            subplot_dict['yhi_list'][subplot_indx] = max(subplot_dict['yhi_list'][subplot_indx], max(i_dos_arr_up[xlo_energy_indx:xhi_energy_indx]), max(i_dos_arr_dw[xlo_energy_indx:xhi_energy_indx]))
                        if subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] not in [None,'None','none'] and subplot_ylo_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])] >= subplot_yhi_list[subplot_dict['arg_list'].index(subplot_dict['arg_list'][subplot_indx])]:
                            subplot_dict['ylo_list'][subplot_indx] = min(subplot_dict['ylo_list'][subplot_indx], min(i_dos_arr_up[xlo_energy_indx:xhi_energy_indx]), min(i_dos_arr_dw[xlo_energy_indx:xhi_energy_indx]))
                            subplot_dict['yhi_list'][subplot_indx] = max(subplot_dict['yhi_list'][subplot_indx], max(i_dos_arr_up[xlo_energy_indx:xhi_energy_indx]), max(i_dos_arr_dw[xlo_energy_indx:xhi_energy_indx]))
                            
            subplot_dict['ylo_list'][subplot_indx] = subplot_dict['ylo_list'][subplot_indx] * lower_spine_scale
            subplot_dict['yhi_list'][subplot_indx] = subplot_dict['yhi_list'][subplot_indx] * upper_spine_scale
            dos_lty = default_lty_list
            plot_lines(active_axes, orbits, doscar_dict, energy, atom_dict, dos_lty, subplot_indx, subplot_dict, dos_file_indx, i_elmt_name, orbit_scaling_dict, general_params_dict)
        #Plot Fermi energy indicator line:
        fermi_level_linewidth = 1.0
        if len(atom_dict['doscar_file_path_list']) == 1:
            #Figure only contain one subplot, single DOS is plotted or multiple DOS files are plotted onto one subplot        
            line1 = [(e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['ylo_list'][subplot_indx]), (e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['yhi_list'][subplot_indx])]
            (line1_xs, line1_ys) = zip(*line1)
            active_axes.plot(line1_xs, line1_ys, linestyle = '--',linewidth = fermi_level_linewidth, color = 'black')
        elif len(atom_dict['doscar_file_path_list']) != 1 and general_params_dict['fermi_shift_zero'] == True and counter_list[subplot_indx] == subplot_arg_unique_count_list[subplot_indx]:
            #Figure contain multiple subplots, Fermi energy is shifted to zero. Each subplot only plot the Fermi energy level once.
            line1 = [(e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],0,1), subplot_dict['ylo_list'][subplot_indx]), (e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],0,1), subplot_dict['yhi_list'][subplot_indx])]
            (line1_xs, line1_ys) = zip(*line1)
            active_axes.plot(line1_xs, line1_ys, linestyle = '--',linewidth = fermi_level_linewidth, color = 'black')
        elif len(atom_dict['doscar_file_path_list']) != 1 and general_params_dict['fermi_shift_zero'] == False:
            #Figure contain multiple subplots, Fermi energy is not shifted to zero.
            e_fermi = e_fermi_list[dos_file_indx]
            line1 = [(e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1), subplot_dict['ylo_list'][subplot_indx]),(e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],0,1), subplot_dict['yhi_list'][subplot_indx])]
            (line1_xs, line1_ys) = zip(*line1)
            active_axes.plot(line1_xs, line1_ys, linestyle = '--',linewidth = fermi_level_linewidth,color = atom_dict['palette_list'][dos_file_indx])
        #axes, ticks, labels and legends
        active_axes.set(title='')
        active_axes.set_xlim(subplot_dict['xlo_list'][subplot_indx], subplot_dict['xhi_list'][subplot_indx])
        active_axes.set_ylim(subplot_dict['ylo_list'][subplot_indx], subplot_dict['yhi_list'][subplot_indx])
        active_axes.margins(x=0.0, y=0.0)
        active_axes.xaxis.set_major_locator(plt.MaxNLocator(5))
        active_axes.yaxis.set_major_locator(plt.MaxNLocator(5))
        active_axes.xaxis.set_minor_locator(AutoMinorLocator(5))
        active_axes.yaxis.set_minor_locator(AutoMinorLocator(5))

        # set the axis major format
        if xaxis_major_format in [None,'None','none']:
            pass
        else:
            # e.g. x_major_format = '%.2f'
            active_axes.xaxis.set_major_formatter(FormatStrFormatter(xaxis_major_format))
        if yaxis_major_format in [None,'None','none']:
            pass
        else:
            active_axes.yaxis.set_major_formatter(FormatStrFormatter(yaxis_major_format))
        plt.xticks(fontsize=general_params_dict['font_size'])
        plt.yticks(fontsize=general_params_dict['font_size'])
        #if the hspace and wspace ==0, then remove the endpoints of the axis label and axis ticks.
        if subplot_dict['share_xy'][0] == True:
            yticks = active_axes.yaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)
##            yticks[-1].label1.set_visible(False)
        if subplot_dict['share_xy'][1] == True:
            xticks = active_axes.xaxis.get_major_ticks()
            xticks[0].label1.set_visible(False)
##            xticks[-1].label1.set_visible(False)
        if subplot_dict['xtick_list'][subplot_indx] == False:
            plt.xticks([])
        if subplot_dict['ytick_list'][subplot_indx] == False:
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
                print('## Error: The value of the dos_mode_dict is None. The dos_mode_dict must be defined!')
        subplot_xlabel_on_off_str = funcs.logic_retn_val(subplot_dict['xlabel_list'][subplot_indx], xlabel_str,'')
        subplot_ylabel_on_off_str = funcs.logic_retn_val(subplot_dict['ylabel_list'][subplot_indx], ylabel_str,'')
        active_axes.set_xlabel(funcs.logic_retn_val(mult_plot_logic, subplot_xlabel_on_off_str, xlabel_str),fontsize=general_params_dict['font_size'])
        active_axes.set_ylabel(funcs.logic_retn_val(mult_plot_logic, subplot_ylabel_on_off_str, ylabel_str),fontsize=general_params_dict['font_size'])
        active_axes.legend(loc='best', frameon = False, fontsize=general_params_dict['font_size'])

    if mult_plot_logic:
        # if the subplot(1,1,1) is added before the other subplots, the subplot(1,1,1) will not be shown.
        # by adding the transparent subplot(1,1,1) after the other subplots can the mainplot labels be shown, or the subplot(111) will be overlapped.
        axes_list[0] = fig_dos.add_subplot(1,1,1)
        axes_list[0].patch.set_alpha(0)
        # control the subplot(1,1,1)
        axes_list[0].set(xticks = [])
        axes_list[0].set(yticks = [])
        axes_list[0].spines['top'].set_visible(False)
        axes_list[0].spines['bottom'].set_visible(False)
        axes_list[0].spines['left'].set_visible(False)
        axes_list[0].spines['right'].set_visible(False)
        axes_list[0].spines['bottom'].set_position(('outward',20))
        if 'TDOS' in atom_indx_list or 0 in atom_indx_list:
            #for the TDOS the value of the y axix is relatively large
            axes_list[0].spines['left'].set_position(('outward',53))
        else:
            axes_list[0].spines['left'].set_position(('outward',30))
        axes_list[0].set_xlabel(funcs.logic_retn_val(mainplot_axis_label_list[0],xlabel_str,''), fontsize = general_params_dict['font_size'])
        axes_list[0].set_ylabel(funcs.logic_retn_val(mainplot_axis_label_list[1],ylabel_str,''), fontsize = general_params_dict['font_size'], labelpad=general_params_dict['font_size']*1.3)

    #Export Figures
    #sysindx_list denotes the systems as integers. Identical doscar_file_path have the same integer value.
    #As the items of atom_doscar_file_path_list goes from the first one to the last one, The corresponding position of sysindx_list is assigned a unique interger value(from 1 to N).
    sysindx_list = []
    unique_sysindx_list = []
    unique_atom_doscar_file_path_list = []
    unique_fermi_list = []
    sysindx = 1
    for idoscar_file_path in atom_dict['doscar_file_path_list']:
        if idoscar_file_path in unique_atom_doscar_file_path_list:
            sysindx_list.append(unique_sysindx_list[unique_atom_doscar_file_path_list.index(idoscar_file_path)])
        else:
            unique_sysindx_list.append(sysindx)
            unique_atom_doscar_file_path_list.append(idoscar_file_path)
            sysindx_list.append(sysindx)
            unique_fermi_list.append(e_fermi_list[sysindx - 1])
            sysindx += 1
    formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time.time()))
    if general_params_dict['smoothing'] == True or general_params_dict['smoothing'] == 'Lorentzian':
        fig_file = os.path.join(output_dir, 'fig_' + str(formatted_time) + '_sys' + ''.join([str(x) for x in sysindx_list]) + '_' + ''.join(dos_atomname_list) + '_lorentz' + str(general_params_dict['smoothing_factor']).replace('.','p') + '.' + general_params_dict['fig_format'])
    else:
        fig_file = os.path.join(output_dir, 'fig_' + str(formatted_time) + '_sys' + ''.join([str(x) for x in sysindx_list]) + '_' + ''.join(dos_atomname_list) + '.' + general_params_dict['fig_format'])
    plt.savefig(fig_file,dpi = general_params_dict['fig_dpi'])
    plt.close()
    #Write Figure information into logfile
    fig_log_file = os.path.splitext(fig_file)[0] + '.log'
    sys_list = []
    for i_sys in range(0,len(atom_dict['indx_list'])):
        sys_list.append('sys' + str(sysindx_list[i_sys]))
    funcs.write_log(logfile, '## The system directorie(s) are given below:')
    log_str = ''
    for i_unique_sys in range(0,len(unique_atom_doscar_file_path_list)):
        log_str = log_str + ('sys' + str(unique_sysindx_list[i_unique_sys]) + ' = r\'' + str(unique_atom_doscar_file_path_list[i_unique_sys]) + '\'' + '\n' +
                  'e_fermi_mod = ' + str(e_fermi_list[i_unique_sys]) + '\n')
    log_str = log_str + (
        'dos_mode_dict = ' + str(dos_mode_dict) + '\n' +
        'vasp_plot.plot_dos(' + '\n' +
        '    atom_doscar_file_path_list = ' + ('[' + ','.join(sys_list) + ']') + ',\n' +
        '    atom_sysname_list = ' + str(atom_sysname_list) + ',\n' +
        '    atom_indx_list = ' + str(atom_indx_list) + ',\n' +
        '    atom_palette_list = ' + str(atom_palette_list) + ',\n' +
        '    atom_subplot_arg_list = ' + str(atom_subplot_arg_list) + ',\n' +
        '    subplot_arg_list = ' + str(subplot_arg_list) + ',\n' +
        '    subplot_xlo_list = ' + str(subplot_xlo_list) + ',\n' +
        '    subplot_xhi_list = ' + str(subplot_xhi_list) + ',\n' +
        '    subplot_ylo_list = ' + str(subplot_ylo_list) + ',\n' +
        '    subplot_yhi_list = ' + str(subplot_yhi_list) + ',\n' +
        '    subplot_xtick_list = ' + str(subplot_xtick_list) + ',\n' +             
        '    subplot_ytick_list = ' + str(subplot_ytick_list) + ',\n' +
        '    subplot_xlabel_list = ' + str(subplot_xlabel_list) + ',\n' +              
        '    subplot_ylabel_list = ' + str(subplot_ylabel_list) + ',\n' +
        '    subplot_share_xy_list = ' + str(subplot_share_xy_list) + ',\n' +
        '    mainplot_axis_label_list = ' + str(mainplot_axis_label_list) + ',\n' +
        '    xtick_direction = ' + '\'' + str(xtick_direction) + '\'' + ',\n' +
        '    ytick_direction = ' + '\'' + str(ytick_direction) + '\'' + ',\n' +
        '    dos_mode_dict = dos_mode_dict' + ',\n' +
        '    fermi_shift_zero = ' + str(fermi_shift_zero) + ',\n' +
        '    peak_analyzer = ' + str(peak_analyzer) + ',\n' +
        '    peak_analyzer_factor = ' + str(initial_peak_analyzer_factor) + ',\n' +
        '    smoothing = ' + str(smoothing) + ',\n' +
        '    smoothing_factor = ' + str(smoothing_factor) + ',\n' +
        '    line_width = ' + str(initial_line_width) + ',\n' +
        '    font_size = ' + str(initial_font_size) + ',\n' +
        '    fig_format = ' + '\'' + str(fig_format) + '\'' + ',\n' +
        '    fig_size = ' + str(fig_size) + ',\n' +
        '    fig_dpi = ' + str(fig_dpi) + ')\n' +
        'fig_file = ' + 'r\'' + fig_file + '\'' + '\n' +
        '################################################\n')
    funcs.write_log(fig_log_file, log_str)
    funcs.write_log(logfile, log_str)
    #return subplot_dict
    return fig_file, subplot_dict

def plot_poscar(poscar_file_path, euler_angle_type = 'zyx', phi = -3, theta = 5, psi = 0, elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
                plot_cell_basis_vector_label = False, plot_atom_label = None, label_size = 16, fig_format = 'png', fig_dpi = 100,
                draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical'):
    '''
    - Descriptions
     * Plot POSCAR model. Euler angles are used to rotate the view of the model.
     * Viewer direction is in x direction. The original orientation: x direction is perpendicular to the paper, z direction is in the paper and point to upper direction
     * Reference for Eulerian angles: Herbert Goldstein and Charles P. Poole Jr. and John L. Safko,Classical Mechanics (3rd Edition),2001.

    - args:
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
     * axis_indicator: Logical value. Whethether to plot the axis indicator or not
     * plot_cell_basis_vector_label: Logical value. Whether to plot the cell basis vector labels( i.e., to label the three basis vectors of the cell as a, b, and c)
     * plot_atom_label: String value. values: 'atom_name', 'atom_index', 'atom_species', 'fix_info', 'position_direct', 'position_cartesian', 'added_atom_data', 'n_xyz', 'n_yzx', 'n_zxy', 'n_xzy', 'n_yxz', 'n_zyx', 'elmt_n_xyz', 'elmt_n_yzx', 'elmt_n_zxy', 'elmt_n_xzy', 'elmt_n_yxz', 'elmt_n_zyx', or None/'None'. It plot atom label to each atom
     * label_size: label size for the label of each atom, the aixs indicator label, default is label_size=16.
     * fig_format: String format. Figformat is a string that defines output figure format. Supported fig_format: 'png', 'eps', 'pdf', 'tif', 'tiff', 'jpg', 'jpeg', 'svg', 'svgz', 'pgf', 'ps', 'raw', 'rgba'     * fig_dpi: float format. The DPI for non-vector graphics.
    '''
    import os
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn.neighbors import KDTree
    import re
    from .. import funcs
    from .. import periodic_table
    from . import vasp_read
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
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
    label_size = float(label_size)
    label_size1 = float(label_size) / 2  # This is for labeling the atom position (the x, y, z coordinates), recommended value label_size1=6.5
    elmt_color_tikz = ['gray','magenta','blue','red','green','yellow','cyan','gray','magenta','blue','red','green','yellow','cyan','gray','magenta','blue','red','green','yellow','cyan','gray','magenta','blue','red','green','yellow','cyan']
    periodic_table_dict = periodic_table.periodic_tab()
    periodic_table_dict['elmt_color']['--'] = 'black'
    periodic_table_dict['atomic_radius']['--'] = 100
    if elmt_color in [None,'None','none']:
        pass
    else:
        periodic_table_dict['elmt_color'].update(elmt_color)

    ###########################################################################
    #Developer maintained codes.
    #Please don't modify the following codes unless you know what you are doing
    ###########################################################################
    fig_file = None
    view_axis = 'x' #view_axis: Decide viewing from which direction: 'x' or 'y' or 'z'. Default is view from x axis
    poscar_file_path = os.path.abspath(poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path)    

    if poscar_dict['file_status'] != 1:
        print('WARNING #2012032033 (from plot_poscar): ' + poscar_file_path + ' does not exist or is empty!')
    else:    
        def disp_and_rotation(poscar_file_path, phi, theta, psi):
            '''
            Description:
                Rotation center is at automatically set as its geometrical center.
                Displace the origin of the model to the center rot_center_arr and then rotate the model according to Eulerian angle
            args:
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
##        if len(poscar_dict['added_atom_data'][0,:]) != 0:
##            selected_added_atom_data_arr = np.concatenate((selected_added_atom_data_arr,np.array([0,0,0])),axis = 0)

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
            #Define the multiplier of axis indicator length, This factor 'mA' places the origin of the axis indicator a distance mA*(axis length) away from the model. The tested optimum value of mA is 1.5
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
        #the scale factor 'scale' is useful for scaling the size of pyplot figure. The tested optimum scaling factor is 0.55
        #if the plotted model are 'truncated' by the boundary of the picture, please enlarge this factor 'scale'
        #scale = 0.45
        scale = 0.55
        fig_size = [max_angstrom_width * scale, max_angstrom_width * fig_aspect_ratio * scale]
        plt.rcParams['figure.figsize'] = (fig_size[0],fig_size[1])
        fig_poscar = plt.figure('fig_poscar')
        fig_poscar.add_subplot(111)
        plt.subplots_adjust(bottom = 0, left = 0, top = 1, right = 1, wspace = 0, hspace = 0)
        ax = plt.gca()
        ##########################
        #Initialize TikZ figure
        ##########################
        #Define TikZ file directory
        tikz_str = ''
        tikz_dir = os.path.join(output_dir, 'poscar_tikz_latex_files')
        funcs.mkdir(tikz_dir)
        formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time.time()))
        fpath_temp, fname_temp = os.path.split(os.path.abspath(poscar_file_path))
        ##print(str(str(os.path.abspath(poscar_file_path)).split('/')).split('\\\\'))
        ##print(str(str(os.path.abspath(poscar_file_path)).split('/')).split('\\\\')[-2])
        system_folder = os.path.split(fpath_temp)[-1]
        ##print(system_folder)
        ##system_folder = str(str(os.path.abspath(poscar_file_path)).split('/')).split('\\\\')[-2]
        tikz_texfile = os.path.join(tikz_dir, 'tikz_' + str(formatted_time) + '_' + system_folder.replace('.','p') + '_' + str(euler_angle_type) + '_' + str(phi).replace('.','p') + '_' + str(theta).replace('.','p') + '_' + str(psi).replace('.','p') + '.tex') 
        
        #Define scaled element radius according to Atomic Radius. This is for the purpose of plotting. 
        reference_atom_radius = periodic_table_dict['atomic_radius'][model_basis_elmt]
        scaled_atomic_radius_for_tikz_dict = {}
        scaled_atomic_radius_for_pyplot_dict = {}
        reduc_fac = 80
        enlargefac_pyplot = 13
        enlargefac_tikz = 0.15
        for i in periodic_table_dict['atomic_radius']:
            if periodic_table_dict['atomic_radius'][i] not in [None,'None','none']:
                scaled_atomic_radius_for_tikz_dict[i] = (periodic_table_dict['atomic_radius'][i] - reduc_fac) / (reference_atom_radius - reduc_fac) * enlargefac_tikz * radius_scale_fac
                scaled_atomic_radius_for_pyplot_dict[i] = (periodic_table_dict['atomic_radius'][i] - reduc_fac) / (reference_atom_radius - reduc_fac) * enlargefac_pyplot * radius_scale_fac 
            else:
                pass
        #Define the number of round off digits of numbers in LaTeX TikZ picture.
        round_digits = 4

        #the scale factor 'scale' is useful for scaling the picture size in TikZ figure. The tested optimum scaling factor is 0.5
        scale = 0.5
        #Write header of Tikz picture
        tikz_str = tikz_str + (r'\documentclass[a4paper,12]{article}' + '\n' +
            r'\usepackage{graphicx}' + '\n' +
            r'\usepackage{tikz}' + '\n' +
            r'\begin{document}' + '\n'
            )
        for i in range(0, len(poscar_dict['elmt_species_arr'])):
            j = elmt_color_tikz[i]
            tikz_str = tikz_str + (r'\pgfdeclareradialshading{ballshading' +
                str(j) + r'}{\pgfpoint{-10bp}{10bp}}{color(0bp)=(' +
                str(j) +
                r'!15!white); color(9bp)=(' + str(j) + '!75!white);color(18bp)=(' + str(j) + r'!70!black); color(25bp)=(' + str(j) + r'!50!black);color(50bp)=(black)}' + '\n')
        tikz_str = tikz_str + (r'\begin{figure}' + '\n' +
            r'\centering' + '\n' +
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
            tikz_str = tikz_str + (r'\draw [help lines] ' +
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
                r'\draw [help lines] ' +
                '(' + str(round(t_vertex[1,hori_indx],round_digits)) + ',' + str(round(t_vertex[1,vert_indx],round_digits)) + ') -- ' +
                '(' + str(round(t_vertex[5,hori_indx],round_digits)) + ',' + str(round(t_vertex[5,vert_indx],round_digits)) + ');' + '\n' +
                r'\draw [help lines] ' +
                '(' + str(round(t_vertex[2,hori_indx],round_digits)) + ',' + str(round(t_vertex[2,vert_indx],round_digits)) + ') -- ' +
                '(' + str(round(t_vertex[6,hori_indx],round_digits)) + ',' + str(round(t_vertex[6,vert_indx],round_digits)) + ');' + '\n' +
                r'\draw [help lines] ' +
                '(' + str(round(t_vertex[3,hori_indx],round_digits)) + ',' + str(round(t_vertex[3,vert_indx],round_digits)) + ') -- ' +
                '(' + str(round(t_vertex[7,hori_indx],round_digits)) + ',' + str(round(t_vertex[7,vert_indx],round_digits)) + ');' + '\n'
                )

        ##################
        #Draw atoms
        ##################
        plt.axis('off')
##        # the plt.scatter only aims at plotting the color bar 
##        # Use plt.scatter won't ensure the atom colors plotted according to the order of the view_atom_indx_arr
##        plt.scatter(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[view_atom_indx_arr[1:][::-1],hori_indx],
##                    t_pos_shift_add_mirror_atoms_and_viewer_point_arr[view_atom_indx_arr[1:][::-1],vert_indx],
##                    c = selected_added_atom_data_arr[view_atom_indx_arr[1:][::-1]],
##                    cmap = colormap_mode,
##                    vmin = colormap_vmin,
##                    vmax = colormap_vmax,
##                    s = 0,
##                    alpha = 1,
##                    )
##        plt.colorbar()
        if vmin_color in [None,'None','none'] or vmax_color in [None,'None','none']:
            vmin_color = 'blue'
            vmax_color = 'red'
        for i in view_atom_indx_arr[1:][::-1]:
            elmtindx = int(np.argwhere(poscar_dict['elmt_species_arr'] == atom_species_add_mirror_atoms_arr[i]))
##              plt.scatter(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],
##                      t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],
##                      color=periodic_table_dict['elmt_color'][atom_species_add_mirror_atoms_arr[i]],
##                      s=scaled_atomic_radius_for_pyplot_dict[atom_species_add_mirror_atoms_arr[i]])
            #Draw atoms for pyplot figure (with shiny shpere effect)
            if atom_species_add_mirror_atoms_arr[i] not in ['--', None, 'None', 'none']:
                original_ball_size = scaled_atomic_radius_for_pyplot_dict[atom_species_add_mirror_atoms_arr[i]]
            else: 
                original_ball_size = 100
            if len(poscar_dict['added_atom_data'][0,:]) != 0 and not isinstance(poscar_dict['added_atom_data'][0,:], str):
                if not isinstance(selected_added_atom_data_arr[i], str):
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
            tikz_str = tikz_str + (r'\pgfpathcircle{\pgfpoint{' +
                str(round(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],round_digits)) +
                'cm}{' + str(round(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],round_digits)) +
                'cm}}{' + str(round(scaled_atomic_radius_for_tikz_dict[atom_species_add_mirror_atoms_arr[i]],round_digits)) +
                'cm}' + '\n' +
                r'\pgfshadepath{ballshading' + str(elmt_color_tikz[elmtindx]) + '}{0}' + '\n' +
                r'\pgfusepath{}' + '\n')
            if plot_atom_label not in [None,'None','none',False,'False']:
                #atom label for pyplot figure
                i_label_size = label_size
                atom_indx_add_mirror_atom_arr = np.array([None] * len(atomname_add_mirror_atom_arr))
                atom_species_add_mirror_atom_arr = np.array([None] * len(atomname_add_mirror_atom_arr))
                atom_fix_info_add_mirror_atom_arr = np.array([None] * len(atomname_add_mirror_atom_arr))
                atom_position_direct_add_mirror_atom_arr = np.array([None] * len(atomname_add_mirror_atom_arr))
                atom_position_cartesian_add_mirror_atom_arr = np.array([None] * len(atomname_add_mirror_atom_arr))
                atom_added_atom_data_add_mirror_atom_arr = np.array([None] * len(atomname_add_mirror_atom_arr))
                atom_number_ortho_add_mirror_atom_arr = np.array([None] * len(atomname_add_mirror_atom_arr))
                atom_species_atom_number_ortho_add_mirror_atom_arr = np.array([None] * len(atomname_add_mirror_atom_arr))
                for indx in range(0, len(atomname_add_mirror_atom_arr)):
                    if atomname_add_mirror_atom_arr[indx] == 'Vi':
                        pass
                    else:
                        atom_indx_add_mirror_atom_arr[indx] = poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx]) + 1
                        atom_species_add_mirror_atom_arr[indx] = poscar_dict['atom_species_arr'][poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx])]
                        atom_fix_info_add_mirror_atom_arr[indx] = ''.join(poscar_dict['fix_arr'][poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx]),:])
                        atom_position_direct_add_mirror_atom_arr[indx] = '(' + ','.join(['{:.2f}'.format(x) for x in poscar_dict['pos_arr'][poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx])][0:3]]) + ')'
                        atom_position_cartesian_add_mirror_atom_arr[indx] = '(' + ','.join(['{:.2f}'.format(x) for x in poscar_dict['pos_arr'][poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx])][3:6]]) + ')'
                        atom_added_atom_data_add_mirror_atom_arr[indx] = poscar_dict['added_atom_data'][poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx]),:]
                        if plot_atom_label in ['n_xyz','n_yzx','n_zxy','n_xzy','n_yxz','n_zyx','elmt_n_xyz','elmt_n_yzx','elmt_n_zxy','elmt_n_xzy','elmt_n_yxz','elmt_n_zyx']:
                            atom_number_ortho_add_mirror_atom_arr[indx] = poscar_dict['atom_number_ortho_' + plot_atom_label[-3::]][poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx])]
                            atom_species_atom_number_ortho_add_mirror_atom_arr[indx] = poscar_dict['atom_species_arr'][poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx])] + '(' + str(poscar_dict['atom_number_ortho_' + plot_atom_label[-3::]][poscar_dict['atomname_list'].index(atomname_add_mirror_atom_arr[indx])]) + ')'

                if plot_atom_label == 'atom_name':
                    i_atom_label = atomname_add_mirror_atom_arr[i]
                elif plot_atom_label == 'atom_index':
                    i_atom_label = str(atom_indx_add_mirror_atom_arr[i])
                elif plot_atom_label == 'atom_species':
                    i_atom_label = atom_species_add_mirror_atom_arr[i]
                elif plot_atom_label == 'fix_info':
                    i_atom_label = atom_fix_info_add_mirror_atom_arr[i]
                elif plot_atom_label == 'position_direct':
                    i_atom_label = atom_position_direct_add_mirror_atom_arr[i]
                    i_label_size = label_size1
                elif plot_atom_label == 'position_cartesian':
                    i_atom_label = atom_position_cartesian_add_mirror_atom_arr[i]
                    i_label_size = label_size1
                elif plot_atom_label in ['n_xyz','n_yzx','n_zxy','n_xzy','n_yxz','n_zyx']:
                    i_atom_label = atom_number_ortho_add_mirror_atom_arr[i]
                elif plot_atom_label in ['elmt_n_xyz','elmt_n_yzx','elmt_n_zxy','elmt_n_xzy','elmt_n_yxz','elmt_n_zyx']:
                    i_atom_label = atom_species_atom_number_ortho_add_mirror_atom_arr[i]
                elif 'added_atom_data' in str(plot_atom_label):
                    if plot_atom_label == 'added_atom_data':
                        if len(atom_added_atom_data_add_mirror_atom_arr[i]) == 0:
                            i_atom_label = 'None'
                            i_label_size = label_size
                        else:
                            i_atom_label = atom_added_atom_data_add_mirror_atom_arr[i][0]
                            i_label_size = label_size1
                    else:
                        temp_column = int(re.search(r'\(([A-Za-z0-9_]+)\)', plot_atom_label).group(1))               
                        i_atom_label = atom_added_atom_data_add_mirror_atom_arr[i][temp_column - 1]
                        i_label_size = label_size1
                elif plot_atom_label == True or plot_atom_label == 'True':
                    i_atom_label = atomname_add_mirror_atom_arr[i]
                else:
                    print('ERROR: vasp_plot error. The value of the plot_atom_label is incorrect, it should be one of the following values: atom_name, atom_index, atom_species, fix_info, position_direct, position_cartesian, added_atom_data, added_atom_data(i_column), True, False, and None. Please check your input')
                    funcs.write_log(logfile, '# ERROR: vasp_plot error. The value of the plot_atom_label is incorrect, it should be one of the following values: atom_name, atom_index, atom_species, fix_info, position_direct, position_cartesian, added_atom_data, added_atom_data(i_column), True, False, and None. Please check your input')
                    exit()
                plt.text(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],
                         t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],
                         i_atom_label,
                         size = i_label_size)   
                #Draw atoms for TikZ figure
                tikz_str = tikz_str + (r'\node [above] at (' + str(round(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,hori_indx],round_digits)) + ','
                    + str(round(t_pos_shift_add_mirror_atoms_and_viewer_point_arr[i,vert_indx],round_digits)) + r') {\scriptsize ' + str(i_atom_label) + '};' + '\n')
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
            tikz_str = tikz_str + (r'\draw [->] ' +
                '(' + str(line1_xs[0]) + ',' + str(line1_ys[0]) + ')--' +
                '(' + str(line1_xs[1]) + ',' + str(line1_ys[1]) + ');' + '\n' +
                r'\node [' + relative_pos_str[0] + '] at (' + str(line1_xs[1]) + ',' + str(line1_ys[1]) + ') {x};'
                )
            tikz_str = tikz_str + (r'\draw [->] ' +
                '(' + str(line2_xs[0]) + ',' + str(line2_ys[0]) + ')--' +
                '(' + str(line2_xs[1]) + ',' + str(line2_ys[1]) + ');' + '\n' +
                r'\node [' + relative_pos_str[1] + '] at (' + str(line2_xs[1]) + ',' + str(line2_ys[1]) + ') {y};'
                )
            tikz_str = tikz_str + (r'\draw [->] ' +
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
                tikz_str = tikz_str + (r'\node [] at (' + str(round(t_vertex[0,hori_indx],round_digits)) + ',' + str(round(t_vertex[0,vert_indx],round_digits)) + ') {o};')
                tikz_str = tikz_str + (r'\node [] at (' + str(round(t_vertex[1,hori_indx],round_digits)) + ',' + str(round(t_vertex[1,vert_indx],round_digits)) + ') {a};')
                tikz_str = tikz_str + (r'\node [] at (' + str(round(t_vertex[4,hori_indx],round_digits)) + ',' + str(round(t_vertex[4,vert_indx],round_digits)) + ') {b};')
                tikz_str = tikz_str + (r'\node [] at (' + str(round(t_vertex[3,hori_indx],round_digits)) + ',' + str(round(t_vertex[3,vert_indx],round_digits)) + ') {c};')
        ######################################
        #Write end lines of TikZ figure code
        ######################################
        tikz_str = tikz_str + (r'\end{tikzpicture}' + '\n' +
            r'\end{figure}' + '\n' + '\n' +
            r'\end{document}'
            )
        ###############
        # Dump figures
        ###############
        with open(tikz_texfile,'w') as f1:
            f1.write(tikz_str)

        (tikz_parent_path , tikz_full_filename) = os.path.split(tikz_texfile)
        (tikz_filename , tikz_file_extension) = os.path.splitext(tikz_full_filename)
        if poscar_dict['added_atom_property'] == None:
            fig_file = os.path.join(output_dir, 'fig_' + tikz_filename.strip('tikz_')  + '.' + fig_format)
        elif poscar_dict['added_atom_property'] != None and draw_colormap == False:
            fig_file = os.path.join(output_dir, 'fig_' + tikz_filename.strip('tikz_')  + '.' + fig_format)
        elif poscar_dict['added_atom_property'] != None and draw_colormap == True:
            added_atom_property_columns_str = funcs.split_line(line = poscar_dict['added_atom_property_columns'], separator = ' ')[colormap_column_indx - 1]
            fig_file = os.path.join(output_dir, 'fig_' + tikz_filename.strip('tikz_')  + '_' + poscar_dict['added_atom_property'] + '_' + added_atom_property_columns_str + '.' + fig_format)
        plt.savefig(fig_file,dpi = fig_dpi)
        plt.close()

        fig_log_file = os.path.splitext(fig_file)[0] + '.log'
        log_str = ''
        log_str = log_str + (
            'elmt_color = ' + str(elmt_color) + '\n' +
            'vasp_plot.plot_poscar(' + '\n' +
            '    poscar_file_path = ' + 'r\'' + str(poscar_file_path) + '\'' + ',\n' +
            '    euler_angle_type = ' + '\'' + str(euler_angle_type) + '\'' + ',\n' +
            '    phi = ' + str(phi) + ',\n' +
            '    theta = ' + str(theta) + ',\n' +
            '    psi = ' + str(psi) + ',\n' +
            '    elmt_color = elmt_color' + ',\n' +
            '    draw_mirror_atom = ' + str(draw_mirror_atom) + ',\n' +
            '    box_on = ' + str(box_on) + ',\n' +
            '    axis_indicator = ' + str(axis_indicator) + ',\n' +
            '    plot_cell_basis_vector_label = ' + str(plot_cell_basis_vector_label) + ',\n' +
            '    plot_atom_label = \'' + str(plot_atom_label) + '\',\n' +                 
            '    label_size = ' + str(label_size) + ',\n' +
            '    fig_format = ' + '\'' + str(fig_format) + '\'' + ',\n' +
            '    fig_dpi = ' + str(fig_dpi) + ',\n' +                
            '    draw_colormap = ' + str(draw_colormap) + ',\n' +
            '    colormap_column_indx = ' + str(colormap_column_indx) + ',\n' +
            '    colormap_vmin = ' + str(colormap_vmin) + ',\n' +
            '    colormap_vmax = ' + str(colormap_vmax) + ',\n' +
            '    vmin_color = ' + '\'' + str(vmin_color) + '\'' + ',\n' +
            '    vmax_color = ' + '\'' + str(vmax_color) + '\'' + ',\n' +
            '    colorbar_alignment = ' + '\'' + str(colorbar_alignment) + '\'' + ')\n' +
            'fig_file = ' + 'r\'' + fig_file + '\'' + '\n')
        log_str = log_str + ( 
            'fig_fileTikZ = ' + 'r\'' + str(tikz_texfile) + '\'' + '\n' +
            '##################################\n')
        #funcs.write_log(fig_log_file, log_str)
        funcs.write_log(logfile, log_str)
##            os.system('pdflatex ' + tikz_texfile)
    return fig_file

def plot_poscar_for_workdir(workdir, euler_angle_type, phi, theta, psi, elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
                            plot_cell_basis_vector_label = False, plot_atom_label = None, poscar_or_contcar = 'POSCAR', fig_format = 'png', fig_dpi = 100,
                            draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical'):

    '''
    - Descriptions
     * Visualization of POSCARs. 
     * The mother folder needs to be specified which contains the folders with POSCARs
     * Euler angles are used to rotate the view of the model
     
    - args:
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
     * plot_atom_label: String value. Atom label
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
        iposcar_file = os.path.join(workdir_list[i], poscar_or_contcar)
        if os.path.exists(iposcar_file) == False:
            funcs.write_log(logfile, '## WARNING:' + iposcar_file +
                           ' does not exist.')
        else:
            plot_poscar(poscar_file_path = iposcar_file, euler_angle_type = euler_angle_type, phi = phi, theta = theta, psi = psi, 
                        elmt_color = elmt_color, draw_mirror_atom = draw_mirror_atom, box_on = box_on, axis_indicator = axis_indicator,
                        plot_cell_basis_vector_label = plot_cell_basis_vector_label, plot_atom_label = plot_atom_label, 
                        fig_format = fig_format, fig_dpi = fig_dpi, draw_colormap = draw_colormap, colormap_column_indx = colormap_column_indx, 
                        colormap_vmin = colormap_vmin, colormap_vmax = colormap_vmax, vmin_color = vmin_color, vmax_color = vmax_color, 
                        colorbar_alignment = colorbar_alignment)
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
    return 0


def plot_poscar_for_sysname(sysname_file, euler_angle_type = 'zyx', phi = -3, theta = 4, psi = 0, elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
                            plot_cell_basis_vector_label = False, plot_atom_label = None, fig_format = 'png', fig_dpi = 100,
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
        plot_atom_label: String value. Atom label
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
            workdir = os.path.join(os.getcwd(), str(filename), sysname)
            poscar_file_path = os.path.join(workdir, 'POSCAR')
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

def plot_bs(infile_path_list, xlim = None, ylim = None, fermi_shift_zero = True, band_list = None,
            interp_on = True, show_band_data_point = False,
            band_gap_label = False,
            band_palette_dict = None, band_lty_dict = None, #for single system
            system_color_list = None, system_lty_list = None, #multiple systems
            spd_and_site_projections_file_path_list = None, projections_point_size_factor = 1,
            legend_on = True, plot_fermi_level = False,
            xtick_direction = 'out', ytick_direction = 'out',
            line_width = 2.0, font_size = 23, fig_format = 'png', fig_size = [15,10], fig_dpi = 600,
            write_band_data = True,
           ):
    '''
    Functions: Plot band structure
    - infile_path_list: A list of input files. The input file can either be EIGENVAL or PROCAR.
    - xlim: the limit of the x axis for the band structure plot.
    - ylim: the limit of the y axis for the band structure plot.
    - fermi_shift_zero: whether the energies are shifted to zero or not.
    - band_list: Defines which bands are to be plotted. Band number starts from 1 (numbered from 1).
    - interp_on: whether to interpolate the band data pionts or not.
    - show_band_data_point: this is only valid when interp_on = True. if show_band_data_point = True, the raw data of the data points will be shown.
    - band_gap_label: Logic value. If band_gap_label = True, then the band gap, CBM, CVM will be labeled.
    - band_palette_dict: this is used to define the color of specific band. Dictionary key is the band index (numbered from 1).
    - band_lty_dict: this is used to define the linestyle of specific band. Dictionary key is the band index (numbered from 1).
    - system_color_list: if multiple infile are to be plotted. the label only label the system index. The line color of different systems, this is used when the band structure for different systems are to be compared.
    - system_lty_list: if multiple infile are to be plotted. the band linestyles are designated for each infile. The line style of different systems, this is used when the band structure for different systems are to be compared.
    - spd_and_site_projections_file_path_list: The parameter denotes the file which can be used to designate the spd- and site projected wave function character or each orbit. This parameter is only valid when the input file has the PROCAR file format. This file contains five columns: spin, mag, ion, orbit, color.
        - spin: spin status of the projected contribution.
        - mag: Magnetization_density. The total and local magnetizations, i.e. rho (orbital-projected contributions to the wavefunctions for each ion) and magnetization density (orbital-projected contributions to the magnetization) in the x(mx), y(my) and z(mz) directions. Please choose from the quadruplet of texts: ['rho', 'mx', 'my', 'mz'] or ['tot', 'mx', 'my', 'mz']. The 'rho' represents 'tot'. The default is ['rho'].
        - ion: the ion name. for example, 'Ni1', 'Al3', 'Te2'.
        - orbit: It contains one of the following orbits: 's', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2'(or 'dx2-y2'), 'tot'. 
        - color: the color for each plotted projected contribution.
        One example of the spd_and_site_projections_file is as follows:
            spin, mag, ions_list  , orbit, color   , legend 
            None, tot, ['Li1',2]  , s    , black   ,  None
            None, tot, ['tot']    , px   , red     ,  Auto
            None, tot, ['Li']     , py   , pink    ,  'Li'
            None, tot, ['Li3']    , pz   , blue    ,  None
            None, mx , [1,3]      , py   , green   ,  Li1+Li3 $p_{y}$ $m_x$
            ...
        in which, the first line 'spin, mag, ion, orbit, color' is the header line.
        The parameter choices for the spd_and_site_projections_file is as follows (In " x / y ", x is the parameter, y is the internal label of the parameters in the program.):
               spin     ,         mag      ,             ions_list     ,      orbit      ,    color    ,     legend
            None /  0   ,      None / 0    ,    ['tot']       / [4 ]   ,    s   /  0     ,    black    ,     None
            up   /  1   ,      tot  / 0    ,    ['Te2','Li5'] / [2, 7] ,    py  /  1     ,    red      ,     Auto
            dw   /  2   ,      mx   / 1    ,    ['Se5']       / [4 ]   ,    pz  /  2     ,    blue     ,     'Li'
                        ,      my   / 2    ,    ['Al3']       / [3 ]   ,    px  /  3     ,    pink     ,     'Li3 py $m_{x}$'
                        ,      mz   / 3    ,    ['Bi3']       / [12]   ,    dxy /  4     ,    green    ,
                        ,                  ,                           ,    dyz /  5     ,             ,   
                        ,                  ,                           ,    dz2 /  6     ,             ,   
                        ,                  ,                           ,    dxz /  7     ,             ,  
                        ,                  ,                           ,    dx2 /  8     ,             ,     
                        ,                  ,                           ,    tot /  9     ,             , 
    - projections_point_size_factor: the scaling factor for the fat band point size. Default value is 1.
    - legend_on: if True, the legend will be shown.
    - plot_fermi_level: if True, the Fermi level will be shown.
    - write_band_data: write band data file or not
    '''
    import os
    import sys
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
    from scipy import optimize
    from scipy.interpolate import interp1d
    import pandas as pd
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from . import vasp_write
    from . import vasp_tools
    from . import vasp_analyze
    from .. import default_params
    from .. import periodic_table

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    log_str = ''
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    #initialization
    default_palette_list = ['black','red','cyan','darkviolet','magenta','gray','darkorange','darkcyan','palegreen',
                   'goldenrod','lightpink','tan','mediumpurple','steelblue','tomato','mediumturquoise',
                   'mediumslateblue','brown','darkseagreen','wheat','seagreen','maroon','deeppink','chocolate',
                   'bisque','plum']
    default_lty_list = ['-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':']

    initial_band_palette_dict = band_palette_dict 
    initial_band_lty_dict = band_lty_dict
    initial_system_color_list = system_color_list
    initial_system_lty_list = system_lty_list

    if system_color_list in [None,'None','none']:
        system_color_list = default_palette_list
    if system_lty_list in [None,'None','none']:
        system_lty_list = default_lty_list

    num_sys = len(infile_path_list)
    for i_sys in range(0, num_sys):        
        infile_path_list[i_sys] = os.path.abspath(str(infile_path_list[i_sys]))
        if os.path.exists(infile_path_list[i_sys]) and os.path.getsize(infile_path_list[i_sys]) > 0:
            pass
        else:
            print('ERROR: vasp_plot error. The file ' + str(infile_path_list[i_sys]) + ' does not exist or is empty')
            funcs.write_log(logfile, '#' + infile_path_list[i_sys] + " doesn't exist or is empty")

    #Initialize the gereral paramters list
    general_params_dict = {}
    general_params_dict['xlim'] = xlim
    general_params_dict['ylim'] = ylim
    general_params_dict['fermi_shift_zero'] = fermi_shift_zero
    general_params_dict['xtick_direction'] = xtick_direction
    general_params_dict['ytick_direction'] = ytick_direction
    #general_params_dict['smoothing'] = smoothing
    #general_params_dict['smoothing_factor'] = smoothing_factor
    general_params_dict['line_width'] = line_width
    general_params_dict['font_size'] = font_size
    general_params_dict['fig_format'] = fig_format
    general_params_dict['fig_size'] = fig_size
    general_params_dict['fig_dpi'] = fig_dpi 

    initial_line_width = general_params_dict['line_width']
    if general_params_dict['line_width'] in [None,'None','none']:
        general_params_dict['line_width'] = 1.0
    initial_font_size = general_params_dict['font_size']
    if general_params_dict['font_size'] in [None,'None','none']:
        general_params_dict['font_size'] = 18
    plt.rcParams['figure.figsize'] = (general_params_dict['fig_size'][0],general_params_dict['fig_size'][1])
    # make use of the 'golden ratio' when treating the sizes of the x-ticks or y-ticks
    golden_ratio = 0.618
    axes_line_width = 2.0
    plt.rcParams['xtick.direction'] = general_params_dict['xtick_direction']
    plt.rcParams['ytick.direction'] = general_params_dict['ytick_direction']
    plt.rcParams['xtick.major.width'] = axes_line_width
    plt.rcParams['ytick.major.width'] = axes_line_width
    plt.rcParams['xtick.minor.width'] = axes_line_width * golden_ratio
    plt.rcParams['ytick.minor.width'] = axes_line_width * golden_ratio
    plt.rcParams['xtick.major.size'] = general_params_dict['font_size'] * (1 - golden_ratio)
    plt.rcParams['ytick.major.size'] = general_params_dict['font_size'] * (1 - golden_ratio)
    plt.rcParams['xtick.minor.size'] = general_params_dict['font_size'] * (1 - golden_ratio) * golden_ratio 
    plt.rcParams['ytick.minor.size'] = general_params_dict['font_size'] * (1 - golden_ratio) * golden_ratio
    plt.rcParams['axes.linewidth'] = axes_line_width

    fig_bs = plt.figure('fig_bs')
    active_axes = plt.subplot(1,1,1)
    active_axes = plt.gca()
    plt.subplots_adjust(bottom =0.13, left = 0.12, top = 0.97, right = 0.98)
    #plt.subplots_adjust(bottom =0.13, left = 0.16, top = 0.98, right = 0.98)

    if general_params_dict['xlim'] in [None, 'None', 'none']:
        xlo = 9999999
        xhi = -9999999
    else:
        if general_params_dict['xlim'][0] in [None, 'None', 'none']:
            xlo = 9999999
        else:
            xlo = general_params_dict['xlim'][0]
        if general_params_dict['xlim'][1] in [None, 'None', 'none']:
            xhi = -9999999
        else:
            xhi = general_params_dict['xlim'][1]
    if general_params_dict['ylim'] in [None, 'None', 'none']:
        ylo = 9999999
        yhi = -9999999
    else:
        if general_params_dict['ylim'][0] in [None, 'None', 'none']:
            ylo = 9999999
        else:
            ylo = general_params_dict['ylim'][0]
        if general_params_dict['ylim'][1] in [None, 'None', 'none']:
            yhi = -9999999
        else:
            yhi = general_params_dict['ylim'][1]
    #plot band structure for each system.
    handle_list = []
    label_list = []

    for i_sys in range(0, num_sys):
        workdir, infile = funcs.file_path_name(infile_path_list[i_sys])
        #automatically check the file type. Choices are: EIGENVAL or PROCAR
        file_type = vasp_tools.check_file_type(infile_path_list[i_sys])

        kpoints_file_path = os.path.join(workdir, 'KPOINTS')
        kpoints_dict = vasp_read.read_kpoints(kpoints_file_path)
        outcar_file_path = os.path.join(workdir, 'OUTCAR')
        outcar_params_dict = vasp_read.read_outcar(outcar_file_path)
        poscar_file_path = os.path.join(workdir, 'POSCAR')
        poscar_dict = vasp_read.read_poscar(poscar_file_path)

        #get KPOINTS information
        if kpoints_dict['file_status'] != 1:
            print('WARNING #20120307 (from read_kpoints): ' + kpoints_file_path + ' does not exist or is empty!')
        else:
            if 'kpoints_xaxis_arr' not in kpoints_dict.keys():
                print('Error (from vasp_plot): Please check your KPOINTS file')
                exit()
            kpoints_arr = kpoints_dict['kpoints_xaxis_arr']
            bs_xaxis_label_list = kpoints_dict['bs_xaxis_label_list']
            bs_xaxis_tick_list = kpoints_dict['bs_xaxis_tick_list']
            kpath_num_intersections = kpoints_dict['num_intersections']
            kpath_num_intersections_interval = kpoints_dict['num_intersections_interval']

        if kpoints_dict['file_status'] != 1 or poscar_dict['file_status'] != 1 or outcar_params_dict['file_status'] != 1:
            print('WARNING #20120308 (from read_kpoints): ' + kpoints_file_path + ' or ' + poscar_file_path + ' or ' + outcar_file_path + ' does not exist or is empty!')
            continue

        if file_type != 'EIGENVAL' and 'PROCAR' not in str(file_type):
            print('WARNING #20120303 (from vasp_plot): Please check the format of the file ' + infile_path_list[i_sys])
            continue

        #get E_fermi
        e_fermi = outcar_params_dict['e_fermi']
        e_fermi_mod = e_fermi + outcar_params_dict['alpha+bet']
        
        if file_type == 'EIGENVAL': 
            #get information from EIGENVAL
            eigenval_dict = vasp_read.read_eigenval(infile_path_list[i_sys])
            eigenval_or_procar_dict = eigenval_dict
            ispin = eigenval_dict['ispin']
            #ispin = outcar_params_dict['ISPIN']
            num_kpoints = eigenval_dict['num_kpoints']
            num_bands = eigenval_dict['num_bands'] 
            num_ions = eigenval_dict['num_ions']
            if ispin == 1:
                #eigs = eigenval_dict['eigs']     
                eigs = eigenval_dict['eigs'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
            elif ispin == 2:
                #eigs_up = eigenval_dict['eigs_up']     
                #eigs_dw = eigenval_dict['eigs_dw']
                eigs_up = eigenval_dict['eigs_up'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
                eigs_dw = eigenval_dict['eigs_dw'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0)
        elif file_type == 'PROCAR_collinear' or file_type == 'PROCAR_noncollinear':
            procar_dict = vasp_read.read_procar(infile_path_list[i_sys])
            eigenval_or_procar_dict = procar_dict
            ispin = procar_dict['ispin']
            num_kpoints = procar_dict['num_kpoints']
            num_bands = procar_dict['num_bands']
            num_ions = procar_dict['num_ions']
            if ispin == 1:
                #eigs = procar_dict['eigs'] 
                eigs = procar_dict['eigs'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
                if file_type == 'PROCAR_collinear':
                    projections = procar_dict['projections']
                elif file_type == 'PROCAR_noncollinear':
                    projections_noncollinear = procar_dict['projections_noncollinear']
            elif ispin == 2:
                #eigs_up = procar_dict['eigs_up']     
                #eigs_dw = procar_dict['eigs_dw']
                eigs_up = procar_dict['eigs_up'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
                eigs_dw = procar_dict['eigs_dw'] + outcar_params_dict['alpha+bet'] - e_fermi_mod * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0)
                if file_type == 'PROCAR_collinear':
                    projections_up = procar_dict['projections_up'] 
                    projections_dw = procar_dict['projections_dw'] 
                elif file_type == 'PROCAR_noncollinear':
                    projections_up_noncollinear = procar_dict['projections_up_noncollinear']
                    projections_dw_noncollinear = procar_dict['projections_dw_noncollinear']
        #write band data to a specific file
        if write_band_data == True:
            vasp_write.write_band_data(eigenval_or_procar_dict, e_fermi = e_fermi)

        # get information form spd_and_site_projections_file_path
        if spd_and_site_projections_file_path_list not in [None, 'None', 'none'] and isinstance(spd_and_site_projections_file_path_list, list):
            if spd_and_site_projections_file_path_list[i_sys] not in [None, 'None', 'none']:
                spd_and_site_projections_file_path = os.path.abspath(spd_and_site_projections_file_path_list[i_sys])
                with open(spd_and_site_projections_file_path, 'r') as f:
                    line = f.readlines()
                    line_0 = funcs.split_line(line = line[0], separator = ',')
                    num_lines_projections_file = len(line)
                    num_columns_projections_file = len(line_0)
                    projections_raw_arr = np.array([None] * (num_lines_projections_file - 1) * num_columns_projections_file)
                    projections_raw_arr.shape = (num_lines_projections_file - 1), num_columns_projections_file
                    projections_arr = np.array([None] * (num_lines_projections_file - 1) * num_columns_projections_file)
                    projections_arr.shape = (num_lines_projections_file - 1), num_columns_projections_file
                    projections_file_header_list = [line_0[0], line_0[1], line_0[2], line_0[3], line_0[4]]
                    for i_line in range(1, num_lines_projections_file):
                        temp_indx = i_line - 1
                        spin = funcs.split_line(line = line[i_line], separator = ',')[0]
                        if len(funcs.extract_alpha_str(spin)) == 0:
                            spin = int(spin)
                        projections_raw_arr[temp_indx, 0] = spin  # spin

                        mag = funcs.split_line(line = line[i_line], separator = ',')[1]
                        if len(funcs.extract_alpha_str(mag)) == 0:
                            mag = int(mag)
                        projections_raw_arr[temp_indx, 1] = mag  # mag

                        ions_list = []
                        ions_temp_list = (funcs.split_line(line = line[i_line], separator = '[')[1]).split(']')[0].split(',')
                        for item in ions_temp_list:
                            try:
                                item = int(item)
                            except:
                                item = item.strip('\'').strip('\"')
                            ions_list.append(item)
                        projections_raw_arr[temp_indx, 2] = ions_list  # ions_list

                        orbit = funcs.split_line(line = line[i_line], separator = ',')[-3] 
                        if len(funcs.extract_alpha_str(orbit)) == 0:
                            orbit = int(orbit)
                        projections_raw_arr[temp_indx, 3] = orbit  # orbit
                        
                        color = funcs.split_line(line = line[i_line], separator = ',')[-2]
                        projections_raw_arr[temp_indx, 4] = color  # color

                        legend = funcs.split_line(line = line[i_line], separator = ',')[-1]
                        projections_raw_arr[temp_indx, 5] = legend  # legend

                    spin_name_list = []
                    mag_name_list = []
                    orbit_name_list = []
                    label_list = []
                    for i_indx in range(0, num_lines_projections_file - 1):
                        spin = projections_raw_arr[i_indx, 0]
                        if spin in [None, 'None', 'none']: # spin unpolarized
                            projections_arr[i_indx, 0] = 0
                            spin_name_list.append('')
                        elif spin in ['up']: # up
                            projections_arr[i_indx, 0] = 1
                            mag_name_list.append('up')
                        elif spin in ['dw', 'down']: # down
                            projections_arr[i_indx, 0] = 2
                            mag_name_list.append('dw')

                        mag = projections_raw_arr[i_indx, 1]
                        if mag in [None, 'None', 'none']:   # tot
                            projections_arr[i_indx, 1] = 0
                            mag_name_list.append('')
                        elif mag in ['tot', 'rho']: # tot
                            projections_arr[i_indx, 1] = 0
                            mag_name_list.append('')
                        elif mag in ['mx', 'Mx']:   # mx
                            projections_arr[i_indx, 1] = 1
                            mag_name_list.append('$m_{x}$')
                        elif mag in ['my', 'My']:   # my
                            projections_arr[i_indx, 1] = 2
                            mag_name_list.append('$m_{y}$')
                        elif mag in ['mz', 'Mz']:   # mz
                            projections_arr[i_indx, 1] = 3
                            mag_name_list.append('$m_{z}$')
                        
                        ions_list = projections_raw_arr[i_indx, 2]
                        atom_indx_list = []
                        for item in ions_list:
                            string_len = len(''.join(filter(str.isalpha, str(item))))
                            number_len = len(''.join(filter(str.isdigit, str(item))))
                            if string_len != 0 and number_len != 0:
                                item_identifier = 'atomname'
                            elif string_len != 0 and number_len == 0:
                                item_identifier = 'elmt_species'
                            elif string_len == 0 and number_len != 0:
                                item_identifier = 'atom_indx'
                            elif string_len == 0 and number_len == 0:
                                print('ERROR (from vasp.plot): Please check your projections file.')
                                exit()
                            
                            if isinstance(item, str) and item_identifier == 'atomname':
                                atom_indx = convert.atomname2indx(os.path.join(workdir, 'POSCAR'), item)
                                atom_indx_list.append(atom_indx)
                            elif isinstance(item, int):
                                atom_indx = item
                                atom_indx_list.append(atom_indx)
                            elif isinstance(item, str) and item_identifier == 'elmt_species' and item == 'tot':
                                atom_indx = num_ions + 1
                                atom_indx_list.append(atom_indx)
                            elif isinstance(item, str) and item_identifier == 'elmt_species' and item != 'tot':
                                atom_species_arr = poscar_dict['atom_species_arr']
                                for i_atom in range(0, len(atom_species_arr)):
                                    atom_species = atom_species_arr[i_atom]
                                    if item == atom_species:
                                        atom_indx_list.append(i_atom + 1)

                        projections_arr[i_indx, 2] = atom_indx_list
                        for i_item in range(0, len(atom_indx_list)): 
                            projections_arr[i_indx, 2][i_item] = projections_arr[i_indx, 2][i_item] - 1

                        orbit = projections_raw_arr[i_indx, -3]
                        #orbit names: s     py     pz     px    dxy    dyz    dz2    dxz    dx2(or dx2-y2)    tot
                        if orbit in ['s']:
                            projections_arr[i_indx, 3] = 0
                            orbit_name_list.append('$s$')
                        elif orbit in ['py']:
                            projections_arr[i_indx, 3] = 1
                            orbit_name_list.append('$p_{y}$')
                        elif orbit in ['pz']:
                            projections_arr[i_indx, 3] = 2
                            orbit_name_list.append('$p_{z}$')
                        elif orbit in ['px']:
                            projections_arr[i_indx, 3] = 3
                            orbit_name_list.append('$p_{x}$')
                        elif orbit in ['dxy']:
                            projections_arr[i_indx, 3] = 4
                            orbit_name_list.append('$d_{xy}$')
                        elif orbit in ['dyz']:
                            projections_arr[i_indx, 3] = 5
                            orbit_name_list.append('$d_{yz}$')
                        elif orbit in ['dz2']:
                            projections_arr[i_indx, 3] = 6
                            orbit_name_list.append('$d_{z^{2}}$')
                        elif orbit in ['dxz']:
                            projections_arr[i_indx, 3] = 7
                            orbit_name_list.append('$d_{xz}$')
                        elif orbit in ['dx2', 'dx2-y2']:
                            projections_arr[i_indx, 3] = 8
                            orbit_name_list.append('$d_{x^{2}-y^{2}}$')
                        elif orbit in ['tot']:   # tot
                            projections_arr[i_indx, 3] = 9
                            orbit_name_list.append('')
                        else:                    # tot
                            projections_arr[i_indx, 3] = 9
                            orbit_name_list.append('')
                        
                        color = projections_raw_arr[i_indx, -2]
                        projections_arr[i_indx, 4] = color

                        legend = projections_raw_arr[i_indx, -1]
                        projections_arr[i_indx, 5] = legend

                    # convert projections_arr to string
                    def arr2str(input_arr):
                        import numpy as np
                        '''
                           convert numpy array to string
                        '''
                        temp_arr = input_arr
                        temp_str = 'np.array(['
                        for i_x in range(0, temp_arr.shape[0]):
                            temp_str = temp_str + '['
                            for i_y in range(0, temp_arr.shape[1]):
                                if isinstance(temp_arr[i_x, i_y], int):
                                    temp_str = temp_str + str(temp_arr[i_x, i_y]) + ',' 
                                elif isinstance(temp_arr[i_x, i_y], str):
                                    temp_str = temp_str + '\'' + temp_arr[i_x, i_y] + '\'' + ','
                                elif isinstance(temp_arr[i_x, i_y], list):
                                    temp_str = temp_str + str(temp_arr[i_x, i_y]) + ','
                            temp_str = temp_str + '],'
                        temp_str = temp_str + '])'
                        return temp_str
                    #log_str = log_str + ('#The columns of the projections_arr are: spin, mag, ions_list, orbit, color, legend.\n') 
                    log_str = log_str + ('import numpy as np' + '\n')
                    log_str = log_str + ('import pandas as pd' + '\n')
                    log_str = log_str + ('projections_raw_arr' + str(i_sys + 1) + ' = ' + arr2str(projections_raw_arr) + '\n')
                    log_str = log_str + ('projections_arr' + str(i_sys + 1) + ' = ' + arr2str(projections_arr) + '\n')
                    # write projections array to file
                    proj_arr = 'projections_raw_arr' + str(i_sys + 1)
                    csv_file = os.path.join(output_dir, 'projections_raw_arr' + str(i_sys + 1) + '.csv')
                    #pd.DataFrame(proj_arr).to_csv(csv_file,header=['spin', 'mag', 'ions_list'      , 'orbit', 'color' , 'legend'],index=False,doublequote=False)
                    log_str = log_str + ('proj_arr' + str(i_sys + 1) + ' = ' + proj_arr + '\n') 
                    log_str = log_str + ('csv_file' + str(i_sys + 1) + ' = \'' + proj_arr + '.csv' + '\'\n') 
                    log_str = log_str + ("pd.DataFrame(proj_arr" + str(i_sys + 1) + ").to_csv(csv_file" + str(i_sys + 1) + ",header=['spin', 'mag', 'ions_list'      , 'orbit', 'color' , 'legend'],index=False,doublequote=False)" + '\n') 


        band_palette_customize = False
        band_lty_customize = False
        if band_palette_dict not in [None,'None','none']:
            band_palette_customize = True
            temp_band_palette_dict = band_palette_dict.copy()
        if band_lty_dict not in [None,'None','none']:
            band_lty_customize = True
            temp_band_lty_dict = band_lty_dict.copy()
        band_palette_dict = {}
        band_lty_dict = {}
        for i_band in range(0, num_bands):
            i_band_indx = i_band + 1
            if band_list in [None,'None','none']:
                pass 
            else:
                if i_band_indx not in band_list:
                    continue
            # line color for each band
            #band_palette_dict[str(i_band_indx)] = 'black'
            # line style for each band
            if ispin == 1:
                band_palette_dict[str(i_band_indx)] = [default_palette_list[0]]
                band_lty_dict[str(i_band_indx)] = ['-']
            elif ispin == 2:
                #band_palette_dict[str(i_band_indx)] = ['black', 'red'] # one for spin up, one for spin down
                band_palette_dict[str(i_band_indx)] = [default_palette_list[0], default_palette_list[1]] # one for spin up, one for spin down
                #band_lty_dict[str(i_band_indx)] = ['-', '--'] # one for spin up, one for spin down
                band_lty_dict[str(i_band_indx)] = ['-', '-'] # one for spin up, one for spin down

            if band_palette_customize == True:
                if ispin == 1:
                    band_palette_dict[str(i_band_indx)] = [temp_band_palette_dict[str(i_band_indx)]]
                elif ispin == 2:
                    band_palette_dict[str(i_band_indx)] = [temp_band_palette_dict[str(i_band_indx)][0], temp_band_palette_dict[str(i_band_indx)][1]]
            if band_lty_customize == True:
                if ispin == 1:
                    band_lty_dict[str(i_band_indx)] = [temp_band_lty_dict[str(i_band_indx)]]
                elif ispin == 2:
                    band_lty_dict[str(i_band_indx)] = [temp_band_lty_dict[str(i_band_indx)][0], temp_band_lty_dict[str(i_band_indx)][1]]

        num_added_legend_item = 0
        for i_band in range(0, num_bands):
            i_band_indx = i_band + 1
            if band_list in [None,'None','none']:
                pass
            else:
                if i_band_indx not in band_list:
                    continue
            if ispin == 1:
                band_arr = eigs[:, i_band]
                xlo = min(xlo, np.min(kpoints_arr))
                xhi = max(xhi, np.max(kpoints_arr))
                ylo = min(ylo, np.min(band_arr))
                yhi = max(yhi, np.max(band_arr))
                band_label = 'system ' + str(i_sys + 1)
                if num_sys == 1: #for single system
                    band_color = band_palette_dict[str(i_band_indx)][0]
                    band_lty = band_lty_dict[str(i_band_indx)][0]
                    band_label = ''
                elif num_sys > 1: #for multiple system
                    band_color = system_color_list[i_sys] 
                    #band_lty = system_lty_list[i_sys] 
                    band_lty = system_lty_list[0] 
                if show_band_data_point == False:
                    band_marker = ''
                elif show_band_data_point == True:
                    band_marker = '.'
                if interp_on == True:
                    #interpolation of bands
                    for i_interval in range(0, kpath_num_intersections_interval):
                        start_indx = i_interval * kpath_num_intersections
                        end_indx = start_indx + kpath_num_intersections
                        kpoints_slice_arr = kpoints_arr[start_indx:end_indx]
                        eigs_slice_arr = band_arr[start_indx:end_indx]
                        interp = interp1d(kpoints_slice_arr, eigs_slice_arr, kind='cubic')
                        kpoints_dense_arr = np.linspace(min(kpoints_slice_arr),max(kpoints_slice_arr),num = num_kpoints * 5,endpoint = True)
                        plot_line, = active_axes.plot(kpoints_dense_arr, interp(kpoints_dense_arr), color = band_color, linestyle = band_lty, linewidth = general_params_dict['line_width'])
                        if spd_and_site_projections_file_path_list in [None, 'None', 'none'] or spd_and_site_projections_file_path_list[i_sys] in [None, 'None', 'none']:
                            plot_dot, = plt.plot(kpoints_slice_arr, eigs_slice_arr, color = band_color, linestyle = '', marker = band_marker, label = band_label)
                elif interp_on == False:
                    if spd_and_site_projections_file_path_list in [None, 'None', 'none'] or spd_and_site_projections_file_path_list[i_sys] in [None, 'None', 'none']:
                        plot_dot, = plt.plot(kpoints_arr, band_arr, color = band_color, linestyle = '', marker = '.', label = band_label)
            elif ispin == 2:
                #band plotting
                band_arr_up = eigs_up[:, i_band]
                band_arr_dw = eigs_dw[:, i_band]
                xlo = min(xlo, np.min(kpoints_arr))
                xhi = max(xhi, np.max(kpoints_arr))
                ylo = min(ylo, np.min(band_arr_up), np.min(band_arr_dw))
                yhi = max(yhi, np.max(band_arr_up), np.max(band_arr_dw))
                if num_sys == 1: #for single system
                    band_color_up = band_palette_dict[str(i_band_indx)][0]
                    band_color_dw = band_palette_dict[str(i_band_indx)][1]
                    band_lty_up = band_lty_dict[str(i_band_indx)][0]
                    band_lty_dw = band_lty_dict[str(i_band_indx)][1]
                elif num_sys > 1: #for multiple system
                    band_color_up = system_color_list[i_sys] 
                    band_color_dw = system_color_list[i_sys]
                    band_lty_up = band_lty_dict[str(i_band_indx)][0]
                    band_lty_dw = band_lty_dict[str(i_band_indx)][1]
                band_label_up = 'spin up'
                band_label_dw = 'spin dw'
                if show_band_data_point == False:
                    band_marker = ''
                elif show_band_data_point == True:
                    band_marker = '.'
                if interp_on == True:
                    #interpolation of bands
                    for i_interval in range(0, kpath_num_intersections_interval):
                        start_indx = i_interval * kpath_num_intersections
                        end_indx = start_indx + kpath_num_intersections
                        kpoints_slice_arr = kpoints_arr[start_indx:end_indx]
                        eigs_slice_arr_up = band_arr_up[start_indx:end_indx]
                        eigs_slice_arr_dw = band_arr_dw[start_indx:end_indx]
                        kpoints_dense_arr = np.linspace(min(kpoints_slice_arr),max(kpoints_slice_arr),num = num_kpoints * 5,endpoint = True)
                        interp_up = interp1d(kpoints_slice_arr, eigs_slice_arr_up, kind='cubic')
                        interp_dw = interp1d(kpoints_slice_arr, eigs_slice_arr_dw, kind='cubic')
                        plot_line_up, = active_axes.plot(kpoints_dense_arr, interp_up(kpoints_dense_arr), color = band_color_up, linestyle = band_lty_up, linewidth = general_params_dict['line_width'])
                        plot_line_dw, = active_axes.plot(kpoints_dense_arr, interp_dw(kpoints_dense_arr), color = band_color_dw, linestyle = band_lty_dw, linewidth = general_params_dict['line_width'])
                        #plot_line_up, = plt.plot(kpoints_slice_arr, eigs_slice_arr_up, color = band_color_up, linestyle = '-', marker = band_marker, label = band_label_up)
                        #plot_line_dw, = plt.plot(kpoints_slice_arr, eigs_slice_arr_dw, color = band_color_dw, linestyle = '-', marker = band_marker, label = band_label_dw)
                        plot_dot_up, = plt.plot(kpoints_slice_arr, eigs_slice_arr_up, color = band_color_up, linestyle = '', marker = band_marker, label = band_label_up)
                        plot_dot_dw, = plt.plot(kpoints_slice_arr, eigs_slice_arr_dw, color = band_color_dw, linestyle = '', marker = band_marker, label = band_label_dw)
                        
                elif interp_on == False:
                    if spd_and_site_projections_file_path_list in [None, 'None', 'none'] or spd_and_site_projections_file_path_list[i_sys] in [None, 'None', 'none']:
                        plot_dot_up, = plt.plot(kpoints_arr, band_arr_up, color = band_color_up, linestyle = '', marker = '.')
                        plot_dot_dw, = plt.plot(kpoints_arr, band_arr_dw, color = band_color_dw, linestyle = '', marker = '.', mfc='none')
            # plot the spd- and site projected wave function character of each band (fat band)
            if spd_and_site_projections_file_path_list not in [None, 'None', 'none'] and spd_and_site_projections_file_path_list[i_sys] not in [None, 'None', 'none']:
                # first, the file in infile_path_list should be the PROCAR file, not EIGENVAL. Codes needs to be added to avoid loading the EIGNEVAL file here!!!!! Or an ERROR would occur.
                markersize_scale = 200
                markersize_scale = markersize_scale * projections_point_size_factor
                poscar_dict = vasp_read.read_poscar(os.path.join(workdir, 'POSCAR'))
                for i_item in range(0, projections_arr.shape[0]):
                    if ispin == 1:
                        num_added_legend_item = num_added_legend_item + 1
                    elif ispin == 2:
                        num_added_legend_item = num_added_legend_item + 2
                    spin = projections_arr[i_item, 0]
                    mag = projections_arr[i_item, 1]
                    ions_list = projections_arr[i_item, 2]
                    ions_raw_list = projections_raw_arr[i_item, 2]
                    orbit = projections_arr[i_item, 3]
                    color = projections_arr[i_item, 4]
                    legend = projections_arr[i_item, 5]
                    if spin == 0:  #spin unpolarized
                        if ispin == 2:
                            print('ERROR (from vasp_plot.plot_bs): This calculation is spin polarized. Please check the spin tag (the first column) of the file ' + str(spd_and_site_projections_file_path_list[i_sys]))
                            exit()
                        band_color = color
                    elif spin == 1: #up
                        if ispin == 1:
                            print('ERROR (from vasp_plot.plot_bs): This calculation is spin unpolarized. Please check the spin tag (the first column) of the file ' + str(spd_and_site_projections_file_path_list[i_sys]))
                            exit()
                        band_color_up = color
                    elif spin == 2: #dw
                        if ispin == 1:
                            print('ERROR (from vasp_plot.plot_bs): This calculation is spin unpolarized. Please check the spin tag (the first column) of the file ' + str(spd_and_site_projections_file_path_list[i_sys]))
                            exit()
                        band_color_dw = color
                    # define the ionic informations of the projected band structure
                    if len(ions_list) == 1:
                        string_len = len(''.join(filter(str.isalpha, str(ions_raw_list[0]))))
                        number_len = len(''.join(filter(str.isdigit, str(ions_raw_list[0]))))
                        if string_len != 0 and number_len != 0:
                            item_identifier = 'atomname'
                        elif string_len != 0 and number_len == 0:
                            item_identifier = 'elmt_species'
                        elif string_len == 0 and number_len != 0:
                            item_identifier = 'atom_indx'
                        if ions_raw_list[0] == 'tot':
                            ions_txt = 'tot'
                        elif ions_raw_list[0] != 'tot' and item_identifier == 'atomname':
                            ions_txt = poscar_dict['atomname_list'][ions_list[0]]
                        elif ions_raw_list[0] != 'tot' and item_identifier == 'elmt_species':
                            ions_txt = ions_list[0]
                    else:
                        if len(ions_raw_list) == 1:
                            string_len = len(''.join(filter(str.isalpha, str(ions_raw_list[0]))))
                            number_len = len(''.join(filter(str.isdigit, str(ions_raw_list[0]))))
                            if string_len != 0 and number_len != 0:
                                item_identifier = 'atomname'
                            elif string_len != 0 and number_len == 0:
                                item_identifier = 'elmt_species'
                            elif string_len == 0 and number_len != 0:
                                item_identifier = 'atom_indx'
                            if ions_raw_list[0] == 'tot':
                                ions_txt = 'tot'
                            elif ions_raw_list[0] != 'tot' and item_identifier == 'elmt_species':
                                ions_txt = ions_raw_list[0]
                            else:
                                ions_txt = str(i_item + 1)
                        else:
                            ions_txt = str(i_item + 1)
                    if file_type == 'PROCAR_collinear':
                        if ispin == 1:
                            if legend in [None, 'None', 'none']:
                                band_label = ''
                            elif legend in ['Auto', 'auto']:
                                band_label = ions_txt + ' ' + str(orbit_name_list[i_item])
                            else:
                                band_label = legend 
                            projections = procar_dict['projections']
                            for i_kpoint in range(0, num_kpoints):
                                markersize = 0
                                for i_ion in ions_list: #sum of the projection contributions for designated atoms
                                    markersize = markersize + projections[i_kpoint, i_band, i_ion, orbit] 
                                plot_dot, = plt.plot(kpoints_arr[i_kpoint], band_arr[i_kpoint], color = band_color, linestyle = '', marker = '.', markersize = markersize * markersize_scale, markeredgecolor = 'None', label = band_label)
                        elif ispin == 2:
                            if legend in [None, 'None', 'none']:
                                band_label = ''
                            elif legend in ['Auto', 'auto']:
                                if spin == 1:
                                    band_label_up = ions_txt + ' ' + str(orbit_name_list[i_item]) + ' up'
                                elif spin == 2:
                                    band_label_dw = ions_txt + ' ' + str(orbit_name_list[i_item]) + ' dw'
                            else:
                                if spin == 1:
                                    band_label_up = legend
                                elif spin == 2:
                                    band_label_dw = legend
                            projections_up = procar_dict['projections_up']
                            projections_dw = procar_dict['projections_dw']
                            for i_kpoint in range(0, num_kpoints):
                                markersize_up = 0
                                markersize_dw = 0
                                for i_ion in ions_list:  #sum of the projection contributions for designated atoms
                                    if spin == 1:
                                        markersize_up = markersize_up + projections_up[i_kpoint, i_band, i_ion, orbit]
                                    elif spin == 2:
                                        markersize_dw = markersize_dw + projections_dw[i_kpoint, i_band, i_ion, orbit]
                                if spin == 1:
                                    plot_dot_up, = plt.plot(kpoints_arr[i_kpoint], band_arr_up[i_kpoint], color = band_color_up, linestyle = '', marker = '.', markersize = markersize_up * markersize_scale, markeredgecolor = 'None', label = band_label_up)
                                elif spin == 2:
                                    plot_dot_dw, = plt.plot(kpoints_arr[i_kpoint], band_arr_dw[i_kpoint], color = band_color_dw, linestyle = '', marker = '.', markersize = markersize_dw * markersize_scale, markeredgecolor = 'None', label = band_label_dw)
                    elif file_type == 'PROCAR_noncollinear':
                        if ispin == 1:
                            if legend in [None, 'None', 'none']:
                                band_label = ''
                            elif legend in ['Auto', 'auto']:
                                band_label = ions_txt + ' ' + str(orbit_name_list[i_item]) + ' ' + str(mag_name_list[mag])
                            else:
                                band_label = legend 
                            projections_noncollinear = procar_dict['projections_noncollinear']
                            for i_kpoint in range(0, num_kpoints):
                                markersize = 0
                                for i_ion in ions_list:  #sum of the projection contributions for designated atoms
                                    markersize = markersize + projections_noncollinear[i_kpoint, i_band, i_ion, mag, orbit]
                                plot_dot, = plt.plot(kpoints_arr[i_kpoint], band_arr[i_kpoint], color = band_color, linestyle = '', marker = '.', markersize = markersize * markersize_scale, markeredgecolor = 'None', label = band_label)
                        elif ispin == 2:
                            if legend in [None, 'None', 'none']:
                                band_label = ''
                            elif legend in ['Auto', 'auto']:
                                if spin == 1:
                                    band_label_up = ions_txt + ' ' + str(orbit_name_list[i_item]) + ' ' + str(mag_name_list[mag]) + ' up'
                                elif spin == 2:
                                    band_label_dw = ions_txt + ' ' + str(orbit_name_list[i_item]) + ' ' + str(mag_name_list[mag]) + ' dw'
                            else:
                                if spin == 1:
                                    band_label_up = legend
                                elif spin == 2:
                                    band_label_dw = legend
                            projections_up_noncollinear = procar_dict['projections_up_noncollinear']
                            projections_dw_noncollinear = procar_dict['projections_dw_noncollinear']
                            for i_kpoint in range(0, num_kpoints):
                                markersize_up = 0
                                markersize_dw = 0
                                for i_ion in ions_list:  #sum of the projection contributions for designated atoms
                                    if spin == 1:
                                        markersize_up = markersize_up + projections_noncollinear[i_kpoint, i_band, i_ion, mag, orbit]
                                    elif spin == 2:
                                        markersize_dw = markersize_dw + projections_noncollinear[i_kpoint, i_band, i_ion, mag, orbit]
                                if spin == 1:
                                    plot_dot_up, = plt.plot(kpoints_arr[i_kpoint], band_arr_up[i_kpoint], color = band_color_up, linestyle = '', marker = '.', markersize = markersize_up * markersize_scale, markeredgecolor = 'None', label = band_label_up)
                                elif spin == 2:
                                    plot_dot_dw, = plt.plot(kpoints_arr[i_kpoint], band_arr_dw[i_kpoint], color = band_color_dw, linestyle = '', marker = '.', markersize = markersize_dw * markersize_scale, markeredgecolor = 'None', label = band_label_dw)
                    if ispin == 1:
                        if num_added_legend_item <= projections_arr.shape[0]:
                            handle_list.append(plot_dot)
                            label_list.append(band_label)
                    elif ispin == 2:
                        if num_added_legend_item <= (projections_arr.shape[0] * 2):
                            if spin == 1:
                                handle_list.append(plot_dot_up)
                                label_list.append(band_label_up)
                            elif spin == 2:
                                handle_list.append(plot_dot_dw)
                                label_list.append(band_label_dw)
        #############################
        # get band gap information
        #############################
        band_gap_dict = vasp_analyze.get_band_gap(
            eigenval_or_procar_dict = eigenval_or_procar_dict,
            outcar_params_dict = outcar_params_dict, 
            kpoints_dict = kpoints_dict,
            )

        #plot and export band gap information
        log_str = log_str + ('# Egap= ' + '{:.6f}'.format(band_gap_dict['band_gap']) + ' (eV)' + '\n') 
        log_str = log_str + ('# Gap Type= ' + str(band_gap_dict['gap_type']) + '\n') 
        if band_gap_dict['CBM'] not in [None, 'None', 'none']:   
            log_str = log_str + ('# CBM = ' + '{:.6f}'.format(band_gap_dict['CBM']) + ' (eV)' + '\n')
        if band_gap_dict['VBM'] not in [None, 'None', 'none']:   
            log_str = log_str + ('# VBM = ' + '{:.6f}'.format(band_gap_dict['VBM']) + ' (eV)' + '\n')
        if ispin == 2:
            # spin up channel
            log_str = log_str + ('# Egap(up)= ' + '{:.6f}'.format(band_gap_dict['band_gap_up']) + ' (eV)' + '\n') 
            log_str = log_str + ('# Gap Type(up)= ' + str(band_gap_dict['gap_type_up']) + '\n') 
            if band_gap_dict['CBM_up'] not in [None, 'None', 'none']:   
                log_str = log_str + ('# CBM(up) = ' + '{:.6f}'.format(band_gap_dict['CBM_up']) + ' (eV)' + '\n')
            if band_gap_dict['VBM_up'] not in [None, 'None', 'none']:   
                log_str = log_str + ('# VBM(up) = ' + '{:.6f}'.format(band_gap_dict['VBM_up']) + ' (eV)' + '\n')
            # spin down channel
            log_str = log_str + ('# Egap(dw)= ' + '{:.6f}'.format(band_gap_dict['band_gap_dw']) + ' (eV)' + '\n') 
            log_str = log_str + ('# Gap Type(dw)= ' + str(band_gap_dict['gap_type_dw']) + '\n') 
            if band_gap_dict['CBM_dw'] not in [None, 'None', 'none']:   
                log_str = log_str + ('# CBM(dw) = ' + '{:.6f}'.format(band_gap_dict['CBM_dw']) + ' (eV)' + '\n')
            if band_gap_dict['VBM_dw'] not in [None, 'None', 'none']:   
                log_str = log_str + ('# VBM(dw) = ' + '{:.6f}'.format(band_gap_dict['VBM_dw']) + ' (eV)' + '\n')

        cbm_vbm_text_label = False
        cbm_vbm_markersize = 20
        if band_gap_label == True:
            if band_gap_dict['band_gap'] != 0:
                plt.plot(band_gap_dict['kpoint_CBM'], band_gap_dict['CBM'], marker = 'o', color = 'blue', markersize = cbm_vbm_markersize)
                plt.plot(band_gap_dict['kpoint_VBM'], band_gap_dict['VBM'], marker = 'o', color = 'blue', markersize = cbm_vbm_markersize)
                if cbm_vbm_text_label == True:
                    plt.text(band_gap_dict['kpoint_CBM'], band_gap_dict['CBM'], '{:.6f}'.format(band_gap_dict['CBM']) + ' eV', color = 'blue', fontsize = font_size * golden_ratio)
                    plt.text(band_gap_dict['kpoint_VBM'], band_gap_dict['VBM'], '{:.6f}'.format(band_gap_dict['VBM']) + ' eV', color = 'blue', fontsize = font_size * golden_ratio)
            if ispin == 2:
                if band_gap_dict['band_gap_up'] != 0:
                    plt.plot(band_gap_dict['kpoint_CBM_up'], band_gap_dict['CBM_up'], marker = '^', color = 'black', markersize = cbm_vbm_markersize)
                    plt.plot(band_gap_dict['kpoint_VBM_up'], band_gap_dict['VBM_up'], marker = '^', color = 'black', markersize = cbm_vbm_markersize)
                    if cbm_vbm_text_label == True:
                        plt.text(band_gap_dict['kpoint_CBM_up'], band_gap_dict['CBM_up'], '{:.6f}'.format(band_gap_dict['CBM_up']) + ' eV', color = 'black', fontsize = font_size * golden_ratio)
                        plt.text(band_gap_dict['kpoint_VBM_up'], band_gap_dict['VBM_up'], '{:.6f}'.format(band_gap_dict['VBM_up']) + ' eV', color = 'black', fontsize = font_size * golden_ratio)
                if band_gap_dict['band_gap_dw'] != 0:
                    plt.plot(band_gap_dict['kpoint_CBM_dw'], band_gap_dict['CBM_dw'], marker = 'v', color = 'red', markersize = cbm_vbm_markersize)
                    plt.plot(band_gap_dict['kpoint_VBM_dw'], band_gap_dict['VBM_dw'], marker = 'v', color = 'red', markersize = cbm_vbm_markersize)
                    if cbm_vbm_text_label == True:
                        plt.text(band_gap_dict['kpoint_CBM_dw'], band_gap_dict['CBM_dw'], '{:.6f}'.format(band_gap_dict['CBM_dw']) + ' eV', color = 'red', fontsize = font_size * golden_ratio)
                        plt.text(band_gap_dict['kpoint_VBM_dw'], band_gap_dict['VBM_dw'], '{:.6f}'.format(band_gap_dict['VBM_dw']) + ' eV', color = 'red', fontsize = font_size * golden_ratio)

        ###################
        #legend setting
        ###################
        if ispin == 1:
            if spd_and_site_projections_file_path_list in [None, 'None', 'none'] or spd_and_site_projections_file_path_list[i_sys] in [None, 'None', 'none']:
                if interp_on == True:
                    handle_list.append((plot_line, plot_dot))
                    label_list.append(band_label)
                elif interp_on == False:
                    handle_list.append((plot_dot))
                    label_list.append(band_label)
        elif ispin == 2:
            if spd_and_site_projections_file_path_list in [None, 'None', 'none'] or spd_and_site_projections_file_path_list[i_sys] in [None, 'None', 'none']:
                if interp_on == True:
                    handle_list.append((plot_line_up, plot_dot_up))
                    handle_list.append((plot_line_dw, plot_dot_dw))
                    label_list.append(band_label_up)
                    label_list.append(band_label_dw)
                elif interp_on == False:
                    handle_list.append((plot_dot_up))
                    handle_list.append((plot_dot_dw))
                    label_list.append(band_label_up)
                    label_list.append(band_label_dw)
        # reinitialize the parameters to avoid them affecting the next loop
        band_palette_dict = None
        band_lty_dict = None
        #plot Fermi level
        if plot_fermi_level == True:
            line1 = [(xlo, (e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1))), (xhi,(e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1)))]
            (line1_xs, line1_ys) = zip(*line1)
            horizontal_linewidth_scale = golden_ratio
            plt.plot(line1_xs, line1_ys, linestyle = '--',linewidth = general_params_dict['line_width'] * horizontal_linewidth_scale, color = 'black')
        else:
            pass

    #plot vertical lines for high symmetry k points
    for i_tick in bs_xaxis_tick_list:     
        line1 = [(i_tick, ylo), (i_tick, yhi)]
        (line1_xs, line1_ys) = zip(*line1)
        vertical_linewidth_scale = golden_ratio
        plt.plot(line1_xs, line1_ys, linestyle = '-',linewidth = general_params_dict['line_width'] * vertical_linewidth_scale, color = 'black')
    ###plot Fermi level
    ##if plot_fermi_level == True:
    ##    line1 = [(xlo, (e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1))), (xhi,(e_fermi_mod * funcs.logic_retn_val(fermi_shift_zero,0,1)))]
    ##    (line1_xs, line1_ys) = zip(*line1)
    ##    horizontal_linewidth_scale = golden_ratio
    ##    plt.plot(line1_xs, line1_ys, linestyle = '--',linewidth = general_params_dict['line_width'] * horizontal_linewidth_scale, color = 'black')
    ##else:
    ##    pass
    
    #axes, ticks, labels and legends
    plt.xticks(bs_xaxis_tick_list, bs_xaxis_label_list, fontsize = general_params_dict['font_size'])
    plt.yticks(fontsize = general_params_dict['font_size'])
    #plt.xlabel('Wave vector', fontsize = general_params_dict['font_size'])
    plt.ylabel('E-E$_{f}$ (eV)', fontsize = general_params_dict['font_size'])
    if legend_on == True:
        #add legend
        ncol = 1
        legend_fontsize_scale = golden_ratio
        if len(handle_list) <= 10:
            ncol = 1
        elif len(handle_list) > 10:
            ncol = 2
        elif len(handle_list) > 20:
            ncol = 3
        elif len(handle_list) > 30:
            ncol = 4
        elif len(handle_list) > 40:
            ncol = 5
        else:
            ncol = 6
        if num_sys != 1:
            lgnd = plt.legend(handles = handle_list, labels = label_list, loc = 'best', frameon = True, fontsize = general_params_dict['font_size'] * legend_fontsize_scale, ncol = ncol)
            # rescale the legend marker sizes to ensure all the marker sizes are the same.
            if spd_and_site_projections_file_path_list not in [None, 'None', 'none']:
                uniform_legend_marker_size = font_size
                for lgnd_handle in lgnd.legendHandles:
                    #lgnd_handle._legmarker.set_markersize(120)
                    lgnd_handle._legmarker.set_markersize(uniform_legend_marker_size)
    active_axes.set(title='')
    #plt.set_xlim(xlim[0], xlim[1])
    if xlim not in [None,'None','none']:
        plt.xlim(xlim[0], xlim[1])
    if ylim not in [None,'None','none']:
        plt.ylim(ylim[0], ylim[1])
    plt.margins(x=0.0, y=0.0)
    #plt.xaxis.set_major_locator(plt.MaxNLocator(5))
    #plt.yaxis.set_major_locator(plt.MaxNLocator(5))
    #plt.xaxis.set_minor_locator(AutoMinorLocator(5))
    #plt.yaxis.set_minor_locator(AutoMinorLocator(5))
    formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time.time()))

    fig_file = os.path.join(output_dir, 'fig_' + str(formatted_time) + '_band.' + general_params_dict['fig_format'])
    plt.savefig(fig_file, dpi = general_params_dict['fig_dpi'])
    plt.close()
    #Write Figure information into logfile
    fig_log_file = os.path.splitext(fig_file)[0] + '.log'
    infile_list = []
    projections_file_list = []
    for i_sys in range(0,len(infile_path_list)):
        infile_list.append('infile' + str(i_sys + 1))
        projections_file_list.append('projections_file' + str(i_sys + 1))
    #funcs.write_log(logfile, '## The infile directories and parameters are given below:')
    funcs.write_log(logfile, '## plot_bs() ')
    for i_sys in range(0,len(infile_path_list)):
        log_str = log_str + ('infile' + str(i_sys + 1) + ' = r\'' + str(infile_path_list[i_sys]) + '\'' + '\n')
    if spd_and_site_projections_file_path_list not in [None, 'None', 'none']:
        for i_sys in range(0,len(infile_path_list)):
            if spd_and_site_projections_file_path_list[i_sys] not in [None, 'None', 'none']:
                log_str = log_str + ('projections_file' + str(i_sys + 1) + ' = r\'' + str(spd_and_site_projections_file_path_list[i_sys]) + '\'' + '\n')
                log_str = log_str + ('##projections_file' + str(i_sys + 1) + ' = ' + 'csv_file' + str(i_sys + 1) + ' #Please uncomment this line when you are running without a projections file' + '\n')
            elif spd_and_site_projections_file_path_list[i_sys] in [None, 'None', 'none']:
                log_str = log_str + ('projections_file' + str(i_sys + 1) + ' = ' + str(spd_and_site_projections_file_path_list[i_sys]) + '\n')
        log_str = log_str + ('spd_and_site_projections_file_path_list = ' + ('[' + ','.join(projections_file_list) + ']') + '\n')
    else:
            log_str = log_str + ('spd_and_site_projections_file_path_list = None' + '\n')
    log_str = log_str + ('band_list = ' + str(band_list) + '\n' +
                         'band_palette_dict = ' + str(initial_band_palette_dict) + '\n' +
                         'band_lty_dict = ' + str(initial_band_lty_dict) + '\n' +
                         'system_color_list = ' + str(initial_system_color_list) + '\n' +
                         'system_lty_list = ' + str(initial_system_lty_list) + '\n'
                         )
    log_str = log_str + (
        'vasp_plot.plot_bs(' + '\n' +
        '    infile_path_list = ' + ('[' + ','.join(infile_list) + ']') + ',\n' +
        '    xlim = ' + str(xlim) + ',\n' + 
        '    ylim = ' + str(ylim) + ',\n' + 
        '    fermi_shift_zero = ' + str(fermi_shift_zero) + ',\n' +
        '    band_list = ' + str(band_list) + ',\n' + 
        '    interp_on = ' + str(interp_on) + ',\n' + 
        '    show_band_data_point = ' + str(show_band_data_point) + ',\n' + 
        '    band_gap_label = ' + str(band_gap_label) + ',\n' + 
        '    band_palette_dict = band_palette_dict' + ',\n' + 
        '    band_lty_dict = band_lty_dict' + ',\n' + 
        '    system_color_list = system_color_list' + ',\n' + 
        '    system_lty_list = system_lty_list' + ',\n' + 
        '    spd_and_site_projections_file_path_list = spd_and_site_projections_file_path_list' + ',\n' +
        '    projections_point_size_factor = ' + str(projections_point_size_factor) + ',\n' +
        '    legend_on = ' + str(legend_on) + ',\n' +
        '    plot_fermi_level = ' + str(plot_fermi_level) + ',\n' +
        '    xtick_direction = ' + '\'' + str(xtick_direction) + '\'' + ',\n' +
        '    ytick_direction = ' + '\'' + str(ytick_direction) + '\'' + ',\n' +
        '    line_width = ' + str(initial_line_width) + ',\n' +
        '    font_size = ' + str(initial_font_size) + ',\n' +
        '    fig_format = ' + '\'' + str(fig_format) + '\'' + ',\n' +
        '    fig_size = ' + str(fig_size) + ',\n' +
        '    fig_dpi = ' + str(fig_dpi) + ',\n' +
        '    write_band_data = ' + str(write_band_data) + ')\n' +
        '################################################\n')
    funcs.write_log(fig_log_file, log_str)
    funcs.write_log(logfile, log_str)
    return fig_file
