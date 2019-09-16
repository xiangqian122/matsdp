# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use("Agg")

def plot_proxigram_csv(proxigram_csv_dir, sysname, visible_elmt_list, interplation_on = False, fig_width=6, fig_height=5, fig_dpi = 600, fig_format = 'png'):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.optimize as opt
    import numpy as np
    from .. import funcs
    from . import apt_read
    
    font_size = 15
    label_size = 15
    line_width = 0.5
    fig_size = (fig_width,fig_height)

    data_set, elmtname_list = apt_read.read_proxigram_csv(os.path.abspath(proxigram_csv_dir))
    plt.rcParams['figure.figsize'] = fig_size
    distance = np.array(data_set[:,0], dtype = np.float32)
    fig = plt.figure('concentration' + str(sysname))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(bottom = 0.10, left = 0.12, top = 0.94, right = 0.99, wspace=0.0, hspace=0.0)
    plt.xticks([])
    plt.yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    a, c = np.random.exponential(size=2)
    b, d = np.random.randn(2)

    ax1 = np.array([None]*len(visible_elmt_list))
    for i_plot in range(0, len(visible_elmt_list)):
        elmt_indx = elmtname_list.index(visible_elmt_list[i_plot])
        ax1[i_plot] = fig.add_subplot(str(len(visible_elmt_list))+'1'+str(i_plot+1))
        concentration = np.array(data_set[:,2*elmt_indx+1], dtype = np.float32)
        
        handle_list = []
        label_list = []
        p1, = plt.plot(distance, concentration,color='black',marker='o', mfc='none', linestyle = '')

        if interplation_on == True:
            (a_, b_, c_, d_), _ = opt.curve_fit(funcs.sigmoid, distance, concentration)
            xdense = np.linspace(min(distance),max(distance),num = 100,endpoint = True)
            y_fit = funcs.sigmoid(xdense, a_, b_, c_, d_)
            p1_fit, = plt.plot(xdense,y_fit,'--',color='black',linewidth=line_width)
            handle_list.append((p1,p1_fit))
            label_list.append(elmtname_list[elmt_indx])
        else:
            handle_list.append((p1))
            label_list.append(elmtname_list[elmt_indx])
        ax1[i_plot].legend(handles = handle_list, labels = label_list, loc = 'best', fontsize = font_size)

        xmin = min(distance)
        xmax = max(distance)
        MarginFraction = 0.09
        plt.vlines(x = 0,
                   ymin = min(concentration)-(max(concentration)-min(concentration))*MarginFraction,
                   ymax = max(concentration)+(max(concentration)-min(concentration))*MarginFraction,
                   linestyle='--',
                   linewidth=0.5)
        plt.ylim(min(concentration)-(max(concentration)-min(concentration))*MarginFraction, max(concentration)+(max(concentration)-min(concentration))*MarginFraction)

        if i_plot < len(visible_elmt_list)-1:
            plt.xticks([])
    plt.sca(ax)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
##    plt.text(xlim[1]*0.03, ylim[1]*1.02, '$\gamma$ phase', fontsize = font_size)
##    plt.text(xlim[1]*0.8, ylim[1]*1.02, '$\gamma^{\prime}$ phase', fontsize = font_size)  

    plt.xlabel('distance from the interface ($nm$)', size = label_size, labelpad=font_size*1.4)
    plt.ylabel('concentration ($at.\%$)', size = label_size, labelpad=font_size*2.3)

    fig_dir = os.getcwd() + './outputs'
    funcs.mkdir(fig_dir)
    fig_file = fig_dir + '/' + 'apt_concentration_profile_' + str(sysname) + '.' + fig_format
    plt.savefig(fig_file, dpi = fig_dpi)
    plt.close()
    return 0

