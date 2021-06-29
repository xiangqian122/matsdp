# -*- coding: utf-8 -*-
def plot_bs(band_dict, 
            xlim = None, ylim = None, 
            fermi_shift_zero = True, e_fermi = 0, 
            band_list = None,
            interp_on = True, plot_band_data_point = False, plot_line = True,
            band_gap_label = False,
            band_palette_dict = None, band_lty_dict = None, #for single system
            system_color_list = None, system_lty_list = None, #multiple systems
            spd_and_site_projections_file_path_list = None, projections_point_size_factor = 1,
            alpha = None,
            legend_on = True, plot_fermi_level = False,
            xtick_direction = 'out', ytick_direction = 'out',
            line_width = 2.0, font_size = 23, fig_format = 'png', fig_size = [15,10], fig_dpi = 600,
            write_band_data = True,
            band_dat_file = None,
            output_fig_file_path = None,
            ):
    '''
    plot the band structure
    '''
    args_dict = locals()
    import os
    import matplotlib.pyplot as plt
    from .. import funcs
    from .. import default_params
    import time
    import numpy as np

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    log_str = ''
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)
    
    default_palette_list = ['black','red','cyan','darkviolet','magenta','gray','darkorange','darkcyan','palegreen',
                   'goldenrod','lightpink','tan','mediumpurple','steelblue','tomato','mediumturquoise',
                   'mediumslateblue','brown','darkseagreen','wheat','seagreen','maroon','deeppink','chocolate',
                   'bisque','plum']
    default_lty_list = ['-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':','-','--','-.',':']

    plot_bs_dict = {}
    plot_bs_dict['fig_file'] = None
    
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
        xlo = 999999
        xhi = -999999
    else:
        if general_params_dict['xlim'][0] in [None, 'None', 'none']:
            xlo = 999999
        else:
            xlo = general_params_dict['xlim'][0]
        if general_params_dict['xlim'][1] in [None, 'None', 'none']:
            xhi = -999999
        else:
            xhi = general_params_dict['xlim'][1]
    if general_params_dict['ylim'] in [None, 'None', 'none']:
        ylo = 999999
        yhi = -999999
    else:
        if general_params_dict['ylim'][0] in [None, 'None', 'none']:
            ylo = 999999
        else:
            ylo = general_params_dict['ylim'][0]
        if general_params_dict['ylim'][1] in [None, 'None', 'none']:
            yhi = -999999
        else:
            yhi = general_params_dict['ylim'][1]


    ispin = 1
    if ispin == 1:
        #eigs = eigenval_dict['eigs']     
        eigs = band_dict['eigs'] - e_fermi * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
    elif ispin == 2:
        eigs_up = band_dict['eigs_up'] - e_fermi * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0) 
        eigs_dn = band_dict['eigs_dn'] - e_fermi * funcs.logic_retn_val(general_params_dict['fermi_shift_zero'],1,0)

    band_color = 'black'
    kpoints_arr = band_dict['kpath_len_list']
    for i_band in range(0, band_dict['num_bands']):
        i_band_indx = i_band + 1
        if band_list is None:
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
            band_label = 'band ' + str(i_band_indx)
            if plot_band_data_point == False:
                band_marker = ''
            elif plot_band_data_point == True:
                band_marker = '.'
            if plot_line == False:
                band_lty = ''
            elif plot_line == True:
                band_lty = default_lty_list[0]
            plot_dot, = plt.plot(kpoints_arr, band_arr, color = band_color, linestyle = band_lty, marker = band_marker, label = band_label, alpha = alpha)
    #plt.set_xlim(xlim[0], xlim[1])
    if xlim not in [None,'None','none']:
        plt.xlim(xlim[0], xlim[1])
    if ylim not in [None,'None','none']:
        plt.ylim(ylim[0], ylim[1])
    plt.margins(x=0.0, y=0.0)

    formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time.time()))
    plot_bs_dict['fig_file'] = os.path.join(output_dir, 'fig_' + str(formatted_time) + '_wannier_band.' + general_params_dict['fig_format'])
    if output_fig_file_path not in [None, 'None', 'none']:
        plot_bs_dict['fig_file'] = output_fig_file_path
    plt.savefig(plot_bs_dict['fig_file'], dpi = general_params_dict['fig_dpi'])
    plt.close()
    return 0

def plot_arc(arc_file_path):
    '''
    Plot the Fermi Arc
    '''
    args_dict = locals()
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib.colors import LogNorm
    import scipy.interpolate

    arc_file_path = os.path.abspath(arc_file_path)
    #data = np.loadtxt(arc_file_path)
    x, y, z = np.loadtxt(arc_file_path, unpack=True)
    ### Set up a regular grid of interpolation points
    ##xi, yi = np.linspace(np.amin(x), np.amax(x), 100), np.linspace(np.amin(y), np.amax(y), 100)
    ##xi, yi = np.meshgrid(xi, yi)
    ### Interpolate
    ##rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    ##zi = rbf(xi, yi)
    ##plt.imshow(zi, vmin = np.amin(z), vmax = np.amax(z), origin='lower',
    ##           extent = (np.amin(x), np.amax(x), np.amin(y), np.amax(y)))

    N = int(len(z)**.5)
    z = z.reshape(N, N)
    plt.imshow(z+10, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),
            cmap=cm.hot, norm=LogNorm())
    plt.contourf(x,y,z, 20, cmap='RdGy')
    plt.colorbar()
    plt.show()
    
    ##x0 = 0.47854750 
    ##y0 = 0.45928792
    ##pts_list = []
    ##circum_point_arr = points_in_circum(x0, y0, r = 0.2)
    ##arc_arr = np.concatenate((x.reshape(-1,1),y.reshape(-1,1)),axis = 1)
    ##for i in range(circum_point_arr.shape[0]):
    ##    point, indx = find_nearest_point(arc_arr, circum_point_arr[i,0], circum_point_arr[i,1])
    ##    pts_list.append(point.tolist() + [z[indx]])
    ##pts_arr = np.array(pts_list)
    ##length = 0
    ##length_list = []
    ##for i in range(1, pts_arr.shape[0]):
    ##    length += np.linalg.norm(pts_arr[i][0:2] - pts_arr[i-1][0:2])
    ##    length_list.append(length)
    ##plt.plot(length_list, pts_arr[1:,2])
    ##plt.show() 
    return 0
