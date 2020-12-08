# -*- coding: utf-8 -*-
'''
This is the Task Management Tool (TMT)
'''

def write_task_summary(task_dir, task_type = 'VASP', band_gap_label = True, plot_fatband = False, band_struct_ylim = None, num_added_bands = 3):
    '''
    write task summary
    band_struct_ylim: A list of length two. The limit of the y axis (i.e. the energy range around the Fermi level). This is used to check the band structure near the Fermi level. If the value is None, then the program will determine the ylim automatically.
    num_added_bands: This is used to automatically determine the added bands around the Fermi level. 
    '''
    import time
    import os
    import numpy as np
    from .. import funcs
    from .. import convert
    from .. import default_params
    from ..vasp import vasp_plot
    from ..vasp import vasp_read
    from ..vasp import vasp_tools
    from ..vasp import vasp_analyze

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    projects_summary_dir = os.path.join(output_dir, defaults_dict['projects_summary_dir_name'])
    task_summary_dir = os.path.join(output_dir, defaults_dict['task_summary_dir_name'])
    funcs.mkdir(output_dir)
    #funcs.mkdir(projects_summary_dir)
    funcs.mkdir(task_summary_dir)
    
    default_palette_list = ['black','red','cyan','darkviolet','magenta','gray','darkorange','darkcyan','palegreen',
                   'goldenrod','lightpink','tan','mediumpurple','steelblue','tomato','mediumturquoise',
                   'mediumslateblue','brown','darkseagreen','wheat','seagreen','maroon','deeppink','chocolate',
                   'bisque','plum']

    task_dir = os.path.abspath(task_dir)
    task_dir_name = os.path.split(task_dir)[-1]
    #project_dir_name = os.path.split(os.path.abspath(os.path.join(task_dir, '../')))[-1]
    #project_summary_dir = os.path.join(projects_summary_dir, 'summary_' + project_dir_name)
    task_summary_dir = os.path.join(task_summary_dir, 'summary_' + task_dir_name)
    #task_summary_dir = os.path.join(project_summary_dir, task_dir_name)
    #funcs.mkdir(project_summary_dir)
    funcs.mkdir(task_summary_dir)

    opt_dir = os.path.join(task_dir, 'opt/')
    scf_dir = os.path.join(task_dir, 'scf/') 
    dos_dir = os.path.join(task_dir, 'dos/')
    bs_dir = os.path.join(task_dir, 'bs/')
    bs_soc_dir = os.path.join(task_dir, 'bs_soc/')

    task_summary_ctext_file_path = os.path.join(task_summary_dir, task_dir_name + '.tex')
    task_summary_ctex_str = ''
    
    current_time = time.time()
    formatted_time = time.strftime('## %Y%m%d %H:%M:%S',time.localtime(current_time))
    #################
    #Write header
    #################
    task_summary_ctex_str = task_summary_ctex_str + (r'\documentclass[a4paper,12]{ctexart}' + '\n' +
        r'\usepackage{amsmath,bm}' + '\n' +
        r'\usepackage[hidelinks]{hyperref}' + '\n' +
        r'\usepackage{longtable}' + '\n' +
        r'\usepackage{graphicx}' + '\n' +
        r'\usepackage{tikz}' + '\n' +
        r'\usepackage{lscape}' + '\n' +
        r'\usepackage{listings}' + '\n')
    task_summary_ctex_str = task_summary_ctex_str + (
r'\usepackage[                                                  ' + '\n' +   
r'    a4paper,                                                  ' + '\n' +
r'    left=3.17cm,right=3.17cm,top=4cm,bottom=3cm]{geometry}    ' + '\n' +       
r'\usepackage{fancyhdr}                                         ' + '\n' +    
r'                                                              ' + '\n' +  
r'\renewcommand{\headrulewidth}{0.6pt}%                         ' + '\n' +  
r'\pagestyle{fancy}                                             ' + '\n' +  
r'\lhead{}                                                      ' + '\n' +        
r'\chead{}                                                      ' + '\n' +            
r'\rhead{}                                                      ' + '\n' +          
r'\lfoot{}                                                      ' + '\n' +          
r'\cfoot{}                                                      ' + '\n' +          
r'\rfoot{}                                                      ' + '\n' +          
r'\fancyhead[CE,CO]{\zihao{-5}\nouppercase{\leftmark}}          ' + '\n' +         
r'\fancyfoot[CE,CO]{\zihao{-5}\thepage}                         ' + '\n')      

    task_summary_ctex_str = task_summary_ctex_str + r'\graphicspath{{' + r'./' + defaults_dict['output_dir_name'] + r'/}{' +  r'../../../' + r'}{./}}' + '\n'

    task_summary_ctex_str = task_summary_ctex_str + r'\renewcommand{\figurename}{Figure}' + '\n'
    task_summary_ctex_str = task_summary_ctex_str + r'\renewcommand{\tablename}{Table}' + '\n'
    task_summary_ctex_str = task_summary_ctex_str + r'\renewcommand{\contentsname}{Table of Contents}' + '\n'

    task_summary_ctex_str = task_summary_ctex_str +(
r'\definecolor{light-gray}{gray}{0.99}                               ' + '\n' +   
r'\definecolor{lbcolor}{rgb}{0.99,0.99,0.99}                         ' + '\n' +       
r'\lstset{%                                                          ' + '\n' + 
r'frame=shadowbox,                                                   ' + '\n' + 
r'language=python,                                                   ' + '\n' + 
r'basicstyle=\footnotesize,       %\scriptsize                       ' + '\n' + 
r'numbers=left,                                                      ' + '\n' + 
r'rulecolor=,                                                        ' + '\n' + 
r'upquote=true,                                                      ' + '\n' + 
r'aboveskip={1.5\baselineskip},                                      ' + '\n' + 
r'numberstyle=\footnotesize,                                         ' + '\n' + 
r'stepnumber=1,                                                      ' + '\n' + 
r'numbersep=5pt,                                                     ' + '\n' + 
r'backgroundcolor=\color{lbcolor},        %\color{white},            ' + '\n' + 
r'showspaces=false,                                                  ' + '\n' + 
r'showstringspaces=false,                                            ' + '\n' + 
r'showtabs=false,                                                    ' + '\n' + 
r'identifierstyle=\ttfamily,                                         ' + '\n' + 
r'keywordstyle=\color[rgb]{0,0,1},                                   ' + '\n' + 
r'%commentstyle=\color[rgb]{0.133,0.545,0.133}                       ' + '\n' + 
r'%stringstyle=\color[rgb]{0.627,0.126,0.941}                        ' + '\n' + 
r'%frame=single,                                                     ' + '\n' + 
r'tabsize=2,                                                         ' + '\n' + 
r'captionpos=b,                                                      ' + '\n' + 
r'prebreak=\raisebox{0ex}[0ex][0ex]{\ensuremath{\hookleftarrow}},    ' + '\n' + 
r'breaklines=true,                                                   ' + '\n' + 
r'extendedchars=true,                                                ' + '\n' + 
r'columns=fixed,                                                     ' + '\n' + 
r'breakatwhitespace=false,                                           ' + '\n' + 
r'xleftmargin=2em,                                                   ' + '\n' +    
r'xrightmargin=2em,                                                  ' + '\n' +  
r'aboveskip=1em                                                      ' + '\n' + 
r'}                                                                  ' + '\n') 

    task_summary_ctex_str = task_summary_ctex_str + r'\newcommand{\ucite}[1]{\textsuperscript{\cite{#1}}}' + '\n' 
    task_summary_ctex_str = task_summary_ctex_str + (
r'\newcommand{\wideunderline}[2][2em]{%                           ' + '\n' +   
r'  \underline{\makebox[\ifdim\width>#1\width\else#1\fi]{#2}}%    ' + '\n' +        
r'}                                                               ' + '\n')          

    task_summary_ctex_str = task_summary_ctex_str + '\n'
    task_summary_ctex_str = task_summary_ctex_str + (r'\begin{document}' + '\n')
    task_summary_ctex_str = task_summary_ctex_str + (
r'\zihao{4}                                                            ' + '\n' +
r'\begin{titlepage}                                                    ' + '\n' +
r'\begin{center}                                                       ' + '\n' +
r'							               ' + '\n' +
r'%\begin{figure}                                                      ' + '\n' +
r'%\centering                                                          ' + '\n' +
r'%\includegraphics[width=0.2\linewidth]{logo.png}                     ' + '\n' +
r'%\end{figure}                                                        ' + '\n' +
r'%\LARGE \textbf{MATSDP}                                              ' + '\n' +
r'								       ' + '\n' +
r'\vspace*{2cm}                                                        ' + '\n' +
r'								       ' + '\n' +
r'\Huge\textbf{Task Summary}                                           ' + '\n' +
r'								       ' + '\n' +
r'								       ' + '\n' +
r'\vspace{3.5cm}                                                       ' + '\n' +
r'						        	       ' + '\n' +
r'\begin{table}[h]                                                     ' + '\n' +
r'\Large                                                               ' + '\n' +
r'\centering                                                           ' + '\n' +
r'\begin{tabular}{rc}                                                  ' + '\n' +
r'Generated with: &\wideunderline[10em]{MATSDP}\\[1cm]                 ' + '\n' +
#r'Project: &\wideunderline[10em]{' + project_dir_name.replace("_", "\_") + r'}\\[1cm]' + '\n' +
r'Task: &\wideunderline[10em]{' + task_dir_name.replace("_", "\_") + r'}\\[1cm]' + '\n' +
r'Time: &\wideunderline[10em]{' + formatted_time.strip('#') + r'}\\[1cm]          ' + '\n' +
r'\end{tabular}                                                        ' + '\n' +
r'\end{table}                                                          ' + '\n' +
r'	             						       ' + '\n' +
r'%\Large xxxx年x月                                                    ' + '\n' +
r'							               ' + '\n' +
r'\end{center}                                                         ' + '\n' +
r'\end{titlepage}                                                      ' + '\n' +
r'				                                       ' + '\n' +
r'\newcommand\frontmatter{%                                            ' + '\n' +
r'    \cleardoublepage                                                 ' + '\n' +
r'  \thispagestyle{empty}                                              ' + '\n' +
r'  \pagenumbering{roman}}                                             ' + '\n' +
r'							               ' + '\n' +
r'\frontmatter                                                         ' + '\n' +
r'						                       ' + '\n' +
r'\tableofcontents                                                     ' + '\n' +
r'\newpage                                                             ' + '\n')

    task_summary_ctex_str = task_summary_ctex_str + (
r'\pagenumbering{arabic}  ' + '\n' +
r'\setcounter{page}{1}    ' + '\n')

    # materials structure and input symmetry

    #path_recorder_list = []
    path1, name1 = os.path.split(task_dir)
    path_recorder_list = [path1]
    path_level_list = [None] * 99
    temp_indx = 0
##    for root, dirs, files in os.walk(task_dir):
##        dirs.sort()
##        temp_indx = temp_indx + 1
##        level = root.replace(task_dir, '').count(os.sep)
##        path_level_list[level] = os.path.basename(root)
##        if level >= 1:
##            path_recorder_list = path_level_list[:level]
##            ##print('****',path_recorder_list)
##        path_recorder_list.append(os.path.basename(root))
##        i_dir_path = os.path.join(path1, '/'.join(path_recorder_list))

    for root, dirs, files in os.walk(task_dir):
        dirs.sort()
        temp_indx = temp_indx + 1
        level = root.replace(task_dir, '').count(os.sep)
        path_level_list[level] = os.path.basename(root)
        if level >= 1:
            path_recorder_list = [path1] + path_level_list[:level]
        path_recorder_list.append(os.path.basename(root))
        path_recorder_list_temp = path_recorder_list.copy()
        path_recorder_list_temp.pop(0)
        i_dir_path = os.path.join('', *path_recorder_list)
        
        i_task_str = '/'.join(path_recorder_list_temp)

        ##print(level, '/'.join(path_recorder_list)) 
        ##print('i_dir_path=', i_dir_path) 
        i_dir_name = os.path.basename(root).strip('\\').strip('#')
        i_task_tex_str = i_task_str.replace("_", "\_")
        i_task_tex_str_wrapped = ''
        wrap_length = 40
        for i_char_indx in range(0, len(i_task_tex_str)):
            i_char = i_task_tex_str[i_char_indx]
            if i_char_indx % wrap_length == 0 and i_char_indx != 0:
                i_task_tex_str_wrapped = i_task_tex_str_wrapped + i_char + '\\\\'
            else:
                i_task_tex_str_wrapped = i_task_tex_str_wrapped + i_char
        if check_job_type(i_dir_path) != 'VASP':
            task_summary_ctex_str = task_summary_ctex_str + (r'\section{' + i_task_tex_str_wrapped + '}' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + ('This is not a valid VASP job (one or more of the following files are missing: INCAR, POSCAR, KPOINTS, POTCAR).' + '\n')
            continue
        task_summary_ctex_str = task_summary_ctex_str + (r'\section{' + i_task_tex_str + '}' + '\n')
        #################################
        #Basic properties
        #################################
        task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{Job Summary}' + '\n')
        outcar_file_path = os.path.join(i_dir_path, 'OUTCAR')
        kpoints_file_path = os.path.join(i_dir_path, 'KPOINTS')
        eigenval_file_path = os.path.join(i_dir_path, 'EIGENVAL')
        poscar_file_path = os.path.join(i_dir_path, 'POSCAR')
        contcar_file_path = os.path.join(i_dir_path, 'CONTCAR')
        job_finished = vasp_tools.vasp_job_finished(outcar_file_path, suppress_warning = True) 
        if job_finished == False:
            job_status = 'Unfinished'
            #task_summary_ctex_str = task_summary_ctex_str + ('Job Status = ' + str(job_finished) + '\n')
            table_label = 'tab:' + i_dir_name + str(temp_indx)
            task_summary_ctex_str = task_summary_ctex_str + (r'The job summary is shown in Table \ref{' + table_label + '}.' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\begin{table}[!h]' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\centering' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\caption{\label{' + table_label + '}Job summary}' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\begin{tabular}{ll}' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'Parameters & Value' + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Job Status & ' + str(job_status) + '\\\\\n')

            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\end{tabular}' + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\end{table}' + '\\\\\n')
        elif job_finished == True:
            job_status = 'Finished'
            outcar_params_dict = vasp_read.read_outcar(outcar_file_path)
            energy_without_entropy = '{:.4f}'.format(outcar_params_dict['energy_without_entropy'])
            toten = '{:.4f}'.format(outcar_params_dict['TOTEN'])
            energy_sigma_0 = '{:.4f}'.format(outcar_params_dict['energy(sigma->0)'])
            elapsed_time_hour = '{:.4f}'.format(convert.time_converter(hour = 0, min = 0, sec = outcar_params_dict['elapsed_time'], unit = 'hour'))
            e_fermi = outcar_params_dict['e_fermi']
            e_fermi_mod = outcar_params_dict['e_fermi_mod']
            kpoints_dict = vasp_read.read_kpoints(kpoints_file_path)
            eigenval_dict = vasp_read.read_eigenval(eigenval_file_path)
            band_gap_dict = vasp_analyze.get_band_gap(
                eigenval_or_procar_dict = eigenval_dict, 
                outcar_params_dict = outcar_params_dict, 
                kpoints_dict = kpoints_dict, 
                fermi_shift_zero = True,
                )
            band_gap = band_gap_dict['band_gap']
            gap_type = band_gap_dict['gap_type']
            cbm = band_gap_dict['CBM']
            vbm = band_gap_dict['VBM']

            band_gap_up = band_gap_dict['band_gap_up']
            gap_type_up = band_gap_dict['gap_type_up']
            cbm_up = band_gap_dict['CBM_up']
            vbm_up = band_gap_dict['VBM_up']

            band_gap_dw = band_gap_dict['band_gap_dw']
            gap_type_dw = band_gap_dict['gap_type_dw']
            cbm_dw = band_gap_dict['CBM_dw']
            vbm_dw = band_gap_dict['VBM_dw']

            if band_gap in [None, 'None', 'none']:
                band_gap = 'None'
            else:
                band_gap = '{:.6f}'.format(band_gap)
            if cbm in [None, 'None', 'none']:
                cbm = 'None'
            else:
                cbm = '{:.6f}'.format(cbm)
            if vbm in [None, 'None', 'none']:
                vbm = 'None'
            else:
                vbm = '{:.6f}'.format(vbm)

            if band_gap_up in [None, 'None', 'none']:
                band_gap_up = 'None'
            else:
                band_gap_up = '{:.6f}'.format(band_gap_up)
            if cbm_up in [None, 'None', 'none']:
                cbm_up = 'None'
            else:
                cbm_up = '{:.6f}'.format(cbm_up)
            if vbm_up in [None, 'None', 'none']:
                vbm_up = 'None'
            else:
                vbm_up = '{:.6f}'.format(vbm_up)

            if band_gap_dw in [None, 'None', 'none']:
                band_gap_dw = 'None'
            else:
                band_gap_dw = '{:.6f}'.format(band_gap_dw)
            if cbm_dw in [None, 'None', 'none']:
                cbm_dw = 'None'
            else:
                cbm_dw = '{:.6f}'.format(cbm_dw)
            if vbm_dw in [None, 'None', 'none']:
                vbm_dw = 'None'
            else:
                vbm_dw = '{:.6f}'.format(vbm_dw)

            poscar_dict = vasp_read.read_poscar(poscar_file_path)
            contcar_dict = vasp_read.read_poscar(contcar_file_path)

            poscar_len_vec_a = poscar_dict['len_vec_a']
            poscar_len_vec_b = poscar_dict['len_vec_b']
            poscar_len_vec_c = poscar_dict['len_vec_c']
            poscar_angle_alpha_degree = poscar_dict['angle_alpha_degree']
            poscar_angle_beta_degree = poscar_dict['angle_beta_degree']
            poscar_angle_gamma_degree = poscar_dict['angle_gamma_degree']
            poscar_box_volume = poscar_dict['box_volume']
            poscar_len_vec_a = '{:.4f}'.format(poscar_len_vec_a)
            poscar_len_vec_b = '{:.4f}'.format(poscar_len_vec_b)
            poscar_len_vec_c = '{:.4f}'.format(poscar_len_vec_c)
            poscar_angle_alpha_degree = '{:.4f}'.format(poscar_angle_alpha_degree)
            poscar_angle_beta_degree = '{:.4f}'.format(poscar_angle_beta_degree)
            poscar_angle_gamma_degree = '{:.4f}'.format(poscar_angle_gamma_degree)
            poscar_box_volume = '{:.4f}'.format(poscar_box_volume)

            contcar_len_vec_a = contcar_dict['len_vec_a']
            contcar_len_vec_b = contcar_dict['len_vec_b']
            contcar_len_vec_c = contcar_dict['len_vec_c']
            contcar_angle_alpha_degree = contcar_dict['angle_alpha_degree']
            contcar_angle_beta_degree = contcar_dict['angle_beta_degree']
            contcar_angle_gamma_degree = contcar_dict['angle_gamma_degree']
            contcar_box_volume = contcar_dict['box_volume']
            contcar_len_vec_a = '{:.4f}'.format(contcar_len_vec_a)
            contcar_len_vec_b = '{:.4f}'.format(contcar_len_vec_b)
            contcar_len_vec_c = '{:.4f}'.format(contcar_len_vec_c)
            contcar_angle_alpha_degree = '{:.4f}'.format(contcar_angle_alpha_degree)
            contcar_angle_beta_degree = '{:.4f}'.format(contcar_angle_beta_degree)
            contcar_angle_gamma_degree = '{:.4f}'.format(contcar_angle_gamma_degree)
            contcar_box_volume = '{:.4f}'.format(contcar_box_volume)

            table_label = 'tab:' + i_dir_name + str(temp_indx)
            task_summary_ctex_str = task_summary_ctex_str + (r'The job summary is shown in Table \ref{' + table_label + '}.' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\begin{table}[!h]' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\centering' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\caption{\label{' + table_label + '}Job summary}' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\begin{tabular}{ll}' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'Parameters & Value' + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Job Status & ' + str(job_status) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('elapsed time(h) & ' + str(elapsed_time_hour) + '\\\\\n')

            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Energy without entropy &' + str(energy_without_entropy) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('TOTEN & ' + str(toten) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('energy(sigma->0) & ' + str(energy_sigma_0) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('E-fermi & ' + str(e_fermi) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('E-fermi (corrected) & ' + str(e_fermi_mod) + '\\\\\n')

            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$|\\vec_{a}|$ ($\mathring{A}$) (POSCAR)& ' + str(poscar_len_vec_a) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$|\\vec_{b}|$ ($\mathring{A}$) (POSCAR)& ' + str(poscar_len_vec_b) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$|\\vec_{c}|$ ($\mathring{A}$) (POSCAR)& ' + str(poscar_len_vec_c) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$\\alpha$ (degree) (POSCAR)& ' + str(poscar_angle_alpha_degree) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$\\beta$ (degree) (POSCAR)& ' + str(poscar_angle_beta_degree) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$\\gamma$ (degree) (POSCAR)& ' + str(poscar_angle_gamma_degree) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Volume ($\mathring{A}^{3}$) (POSCAR)& ' + str(poscar_box_volume) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$|\\vec_{a}|$ ($\mathring{A}$) (CONTCAR)& ' + str(contcar_len_vec_a) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$|\\vec_{b}|$ ($\mathring{A}$) (CONTCAR)& ' + str(contcar_len_vec_b) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$|\\vec_{c}|$ ($\mathring{A}$) (CONTCAR)& ' + str(contcar_len_vec_c) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$\\alpha$ (degree) (CONTCAR)& ' + str(contcar_angle_alpha_degree) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$\\beta$ (degree) (CONTCAR)& ' + str(contcar_angle_beta_degree) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('$\\gamma$ (degree) (CONTCAR)& ' + str(contcar_angle_gamma_degree) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Volume ($\mathring{A}^{3}$) (CONTCAR)& ' + str(contcar_box_volume) + '\\\\\n')

            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + ('E_{gap} (eV) & ' + str(band_gap) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Gap Type & ' + str(gap_type) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('CBM (eV) & ' + str(cbm) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('VBM (eV) & ' + str(vbm) + '\\\\\n')

            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + ('E_{gap}(up) (eV) & ' + str(band_gap_up) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Gap Type(up) & ' + str(gap_type_up) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('CBM(up) (eV) & ' + str(cbm_up) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('VBM(up) (eV) & ' + str(vbm_up) + '\\\\\n')

            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Egap(down) (eV) & ' + str(band_gap_dw) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('Gap Type(down) & ' + str(gap_type_dw) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('CBM(down) (eV) & ' + str(cbm_dw) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + ('VBM(down) (eV) & ' + str(vbm_dw) + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\hline' + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\end{tabular}' + '\\\\\n')
            task_summary_ctex_str = task_summary_ctex_str + (r'\end{table}' + '\\\\\n')
        ###############################
        # Model
        ###############################
        task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{Model}' + '\n')
        poscar_file_path = os.path.join(i_dir_path, 'POSCAR')
        if os.path.exists(poscar_file_path) and os.path.isfile(poscar_file_path):
            ###front view
            ##fig_poscar1 = vasp_plot.plot_poscar(poscar_file_path, euler_angle_type = 'zyx', phi = -3, theta = 5, psi = 0, 
            ##    elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
            ##    plot_cell_basis_vector_label = False, plot_atom_label = 'atom_name', label_size = 16, fig_format = 'png', fig_dpi = 150,
            ##    draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical')
            ###side view
            ##fig_poscar2 = vasp_plot.plot_poscar(poscar_file_path, euler_angle_type = 'zxy', phi = -87, theta = 5, psi = 0, 
            ##    elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
            ##    plot_cell_basis_vector_label = False, plot_atom_label = 'atom_name', label_size = 16, fig_format = 'png', fig_dpi = 150,
            ##    draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical')
            ###top view
            ##fig_poscar3 = vasp_plot.plot_poscar(poscar_file_path, euler_angle_type = 'yxz', phi = 93, theta = 5, psi = 0, 
            ##    elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
            ##    plot_cell_basis_vector_label = False, plot_atom_label = 'atom_name', label_size = 16, fig_format = 'png', fig_dpi = 150,
            ##    draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical')

            ##fig_poscar_list = [fig_poscar1, fig_poscar2, fig_poscar3]
            viewing_direction_list = ['front view', 'side view', 'top view']
            #for i_fig_poscar_indx in range(0, 3):
            euler_angle_type_phi_theta_psi_list_dict= {}
            euler_angle_type_phi_theta_psi_list_dict['front view'] = ['zyx',-3, 5, 0]
            euler_angle_type_phi_theta_psi_list_dict['side view'] = ['zxy',-87, 5, 0]
            euler_angle_type_phi_theta_psi_list_dict['psi view'] = ['yxz',93, 5, 0]
           
            ##for i_fig_poscar_indx in [0, 1, 2]:
            for i_fig_poscar_indx in [0]:
                viewing_direction = viewing_direction_list[i_fig_poscar_indx]
                i_fig_poscar = vasp_plot.plot_poscar(
                    poscar_file_path,
                    euler_angle_type = euler_angle_type_phi_theta_psi_list_dict[viewing_direction][0], 
                    phi = euler_angle_type_phi_theta_psi_list_dict[viewing_direction][1], 
                    theta = euler_angle_type_phi_theta_psi_list_dict[viewing_direction][2], 
                    psi = euler_angle_type_phi_theta_psi_list_dict[viewing_direction][3], 
                    elmt_color = None, draw_mirror_atom = True, box_on = True, axis_indicator = True,
                    plot_cell_basis_vector_label = False, plot_atom_label = 'atom_name', label_size = 16, fig_format = 'png', fig_dpi = 150,
                    draw_colormap = False, colormap_column_indx = 1, colormap_vmin = None, colormap_vmax = None, vmin_color = 'blue', vmax_color = 'red', colorbar_alignment = 'vertical')
                #i_fig_poscar = fig_poscar_list[i_fig_poscar_indx]
                if i_fig_poscar in [None, 'None', 'none']:
                    continue
                parent_path, filename1 = funcs.file_path_name(i_fig_poscar)
                tex_fig_label1 = 'fig:model' + str(i_fig_poscar_indx) +'_' + str(i_dir_name) 
                task_summary_ctex_str = task_summary_ctex_str + (r'The model (' + viewing_direction + ') is shown in Figure ' +  r'\ref{' + tex_fig_label1 +  '}.' + '\n\n')
                task_summary_ctex_str = task_summary_ctex_str + (r'\begin{figure}[h!]' + '\n')
                task_summary_ctex_str = task_summary_ctex_str + (r'\centering' + '\n')
                task_summary_ctex_str = task_summary_ctex_str + (r'\includegraphics[width=0.8\linewidth, height=10cm, keepaspectratio]{' + filename1 + '}' + '\n')
                task_summary_ctex_str = task_summary_ctex_str + (r'\caption{\label{' + tex_fig_label1 + '}Model (front view)' + '}' + '\n')
                task_summary_ctex_str = task_summary_ctex_str + (r'\end{figure}' + '\n')

        #################################
        #inputs of the calculations
        #################################
        task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{Inputs Summary}' + '\n')
        # materials structure
        ##for i_subtask in task_full_item_list: 
        ##    if i_subtask in task_subdir_list:
        ##        i_dir_path = os.path.join(task_dir, i_subtask)
        ##i_task_tex_str = i_subtask.replace("_", "\_")
        #task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{' + i_task_tex_str + '}' + '\n')
        ##if i_dir_path == 'opt':
        task_summary_ctex_str = task_summary_ctex_str + (r'\subsubsection{POSCAR}' + '\n')
        poscar_str = funcs.get_file_str(os.path.join(i_dir_path, 'POSCAR'))
        task_summary_ctex_str = task_summary_ctex_str + (r'\begin{lstlisting}' + '\n')
        task_summary_ctex_str = task_summary_ctex_str + (poscar_str) + '\n'
        task_summary_ctex_str = task_summary_ctex_str + (r'\end{lstlisting}' + '\n')

        task_summary_ctex_str = task_summary_ctex_str + (r'\subsubsection{INCAR}' + '\n')
        incar_file_path = os.path.join(i_dir_path, 'INCAR')
        if os.path.exists(incar_file_path) and os.path.isfile(incar_file_path): 
            incar_str = funcs.get_file_str(incar_file_path)
            task_summary_ctex_str = task_summary_ctex_str + (r'\begin{lstlisting}' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (incar_str) + '\n'
            task_summary_ctex_str = task_summary_ctex_str + (r'\end{lstlisting}' + '\n')
 
        kpoints_file_path = os.path.join(i_dir_path, 'KPOINTS')
        if os.path.exists(kpoints_file_path) and os.path.isfile(kpoints_file_path): 
            task_summary_ctex_str = task_summary_ctex_str + (r'\subsubsection{KPOINTS}' + '\n')
            kpoints_str = funcs.get_file_str(kpoints_file_path)
            task_summary_ctex_str = task_summary_ctex_str + (r'\begin{lstlisting}' + '\n')
            task_summary_ctex_str = task_summary_ctex_str + (kpoints_str) + '\n'
            task_summary_ctex_str = task_summary_ctex_str + (r'\end{lstlisting}' + '\n')
        #################################
        #Electronic structure
        #################################
        task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{Electronic Structure}' + '\n')
        # electronic structure - band structure
        #task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{Band Structure}' + '\n')
        ##for i_bs_name in ['bs', 'bs_soc']:
        kpoints_dict = vasp_read.read_kpoints(kpoints_file_path)

        ##if kpoints_dict['file_status'] != 1:
        ##    print('WARNING #2012032054 (from task_manager): ' + kpoints_file_path + ' does not exist or is empty!')
        ##else: 
        ##    print(kpoints_dict['file_status'])   
        if 'scheme' in kpoints_dict.keys() and job_finished == True:
            if kpoints_dict['scheme'] in ['l', 'L']:
                #i_task_tex_str = funcs.logic_retn_val(i_bs_name == 'bs_soc',' (SOC)','(w/o SOC)')
                ##if i_bs_name not in task_subdir_list:
                ##    task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{' + 'Band Structure' + i_task_tex_str + '}' + '\n')
                ##    task_summary_ctex_str = task_summary_ctex_str + ('None\n')
                ##    task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{Fat Band' + i_task_tex_str + '}' + '\n')
                ##    task_summary_ctex_str = task_summary_ctex_str + ('None\n')
                ##if i_bs_name in task_subdir_list:
                ########################
                #band structure
                ########################
                task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{' + 'Band Structure}' + '\n')
                eigenval_file_path = os.path.join(i_dir_path, 'EIGENVAL')
                if  os.path.exists(eigenval_file_path) and os.path.isfile(eigenval_file_path):
                    if eigenval_dict['ispin'] == 1:
                        band_indx_max = band_gap_dict['CBM_band_indx'] + num_added_bands 
                        band_indx_min = band_gap_dict['VBM_band_indx'] - num_added_bands 
                        if band_indx_max > band_gap_dict['num_bands']:
                            band_indx_max = band_gap_dict['num_bands'] - 1
                        if band_indx_min < 0:
                            band_indx_min = 0
                        e_max = np.max(band_gap_dict['eigs'][:, band_indx_max])
                        e_min = np.min(band_gap_dict['eigs'][:, band_indx_min])
                    elif eigenval_dict['ispin'] == 2: 
                        band_indx_max_up = band_gap_dict['CBM_band_indx_up'] + num_added_bands 
                        band_indx_min_up = band_gap_dict['VBM_band_indx_up'] - num_added_bands 
                        if band_indx_max_up > band_gap_dict['num_bands']:
                            band_indx_max_up = band_gap_dict['num_bands'] - 1
                        if band_indx_min_up < 0:
                            band_indx_min_up = 0
                        band_indx_max_dw = band_gap_dict['CBM_band_indx_dw'] + num_added_bands 
                        band_indx_min_dw = band_gap_dict['VBM_band_indx_dw'] - num_added_bands 
                        if band_indx_max_dw > band_gap_dict['num_bands']:
                            band_indx_max_dw = band_gap_dict['num_bands'] - 1
                        if band_indx_min_dw < 0:
                            band_indx_min_dw = 0
                        e_max_up = np.max(band_gap_dict['eigs_up'][:, band_indx_max_up])
                        e_min_up = np.min(band_gap_dict['eigs_up'][:, band_indx_min_up])
                        e_max_dw = np.max(band_gap_dict['eigs_dw'][:, band_indx_max_dw])
                        e_min_dw = np.min(band_gap_dict['eigs_dw'][:, band_indx_min_dw])
                        e_max = max(e_max_up, e_max_dw)
                        e_min = min(e_min_up, e_min_dw)
                    if band_struct_ylim in [None, 'None', 'none']:
                        band_struct_ylim = [e_min, e_max]
                    fig_bs = vasp_plot.plot_bs(infile_path_list = [eigenval_file_path], xlim = None, ylim = band_struct_ylim, fermi_shift_zero = True, band_list = None,
                        interp_on = True, show_band_data_point = False,
                        band_gap_label = band_gap_label,
                        band_palette_dict = None, band_lty_dict = None, #for single system
                        system_color_list = None, system_lty_list = None, #multiple systems
                        spd_and_site_projections_file_path_list = None, projections_point_size_factor = 1,
                        legend_on = True, plot_fermi_level = True,
                        xtick_direction = 'out', ytick_direction = 'out',
                        line_width = 2.0, font_size = 18, fig_format = 'png', fig_size = [15,10], fig_dpi = 150)
                    ##fig_bs1 = vasp_plot.plot_bs(infile_path_list = [eigenval_file_path], xlim = None, ylim = band_struct_ylim, fermi_shift_zero = True, band_list = None,
                    ##    interp_on = True, show_band_data_point = False,
                    ##    band_gap_label = band_gap_label,
                    ##    band_palette_dict = None, band_lty_dict = None, #for single system
                    ##    system_color_list = None, system_lty_list = None, #multiple systems
                    ##    spd_and_site_projections_file_path_list = None, projections_point_size_factor = 1,
                    ##    legend_on = True, plot_fermi_level = True,
                    ##    xtick_direction = 'out', ytick_direction = 'out',
                    ##    line_width = 2.0, font_size = 18, fig_format = 'png', fig_size = [15,10], fig_dpi = 150)
                    parent_path, filename = funcs.file_path_name(fig_bs)
                    # band structure over the whole energy range
                    tex_fig_label = 'fig:bs_' + i_task_str.strip('\\').strip('#')
                    task_summary_ctex_str = task_summary_ctex_str + (r'The band structure is shown in Figure ' +  r'\ref{' + tex_fig_label +  '}.' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\begin{figure}[h!]' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\centering' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\includegraphics[width=0.8\linewidth]{' + filename + '}' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\caption{\label{' + tex_fig_label + '}Band structure}' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\end{figure}' + '\n')
                    ##parent_path, filename1 = funcs.file_path_name(fig_bs1)
                    ### band structure in a specific energy range (usually near the Fermi level)
                    ##tex_fig_label1 = 'fig:bs1_' + i_task_str.strip('\\').strip('#')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'The band structure (in a specific energy range) of is shown in Figure ' +  r'\ref{' + tex_fig_label1 +  '}.' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\begin{figure}[h!]' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\centering' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\includegraphics[width=0.8\linewidth]{' + filename1 + '}' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\caption{\label{' + tex_fig_label1 + '}Band structure (in a specific energy range)}' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\end{figure}' + '\n')
                ################################
                #band structure (fat band)
                ################################
                poscar_file_path = os.path.join(i_dir_path, 'POSCAR')
                procar_file_path = os.path.join(i_dir_path, 'PROCAR')
                if plot_fatband == True and os.path.exists(poscar_file_path) and os.path.isfile(poscar_file_path) and os.path.exists(procar_file_path) and os.path.isfile(procar_file_path):
                    task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{Fat Band}' + '\n')
                    #prepare projections file
                    projections_txt = ''
                    projections_txt = projections_txt +(
r'spin, mag, ions_list  , orbit    , color      , legend' + '\n'
)                   
                    poscar_dict = vasp_read.read_poscar(poscar_file_path)
                    procar_dict = vasp_read.read_procar(procar_file_path)
                    elmt_species_arr = poscar_dict['elmt_species_arr']
                    ispin = procar_dict['ispin']
                    item_counter = 0
                    for i_elmt_species in elmt_species_arr:
                        if ispin == 1:
                            projections_txt = projections_txt +(
r"None, tot, " + r"['" + i_elmt_species + r"']     , tot      , " + default_palette_list[item_counter] + r"      , Auto" + "\n"
)
                        elif ispin == 2:
                            for i_spin in range(0, ispin):
                                spin_label = funcs.logic_retn_val(i_spin == 0,'up','dw') 
                                projections_txt = projections_txt +(
spin_label + r"  , tot, " + r"['" + i_elmt_species + r"']     , tot      , " + default_palette_list[item_counter] + r"      , Auto" + "\n"
)
                                item_counter = item_counter + 1
                    #write projections text into projections file
                    projections_file_path = os.path.join(task_summary_dir, 'projections_elements.txt')
                    funcs.write_file(projections_txt, projections_file_path, mode = 'w')
                    ####################################################################################################
                    # the fatband plotting is time consuming, turn on the switch if you want to plot the fatband
                    ####################################################################################################
                    #plot band structure (SOC)
                    if procar_dict['ispin'] == 1:
                        band_indx_max = band_gap_dict['CBM_band_indx'] + num_added_bands 
                        band_indx_min = band_gap_dict['VBM_band_indx'] - num_added_bands 
                        if band_indx_max > band_gap_dict['num_bands']:
                            band_indx_max = band_gap_dict['num_bands'] - 1
                        if band_indx_min < 0:
                            band_indx_min = 0
                        e_max = np.max(band_gap_dict['eigs'][:, band_indx_max])
                        e_min = np.min(band_gap_dict['eigs'][:, band_indx_min])
                    elif procar_dict['ispin'] == 2: 
                        band_indx_max_up = band_gap_dict['CBM_band_indx_up'] + num_added_bands 
                        band_indx_min_up = band_gap_dict['VBM_band_indx_up'] - num_added_bands 
                        if band_indx_max_up > band_gap_dict['num_bands']:
                            band_indx_max_up = band_gap_dict['num_bands'] - 1
                        if band_indx_min_up < 0:
                            band_indx_min_up = 0
                        band_indx_max_dw = band_gap_dict['CBM_band_indx_dw'] + num_added_bands 
                        band_indx_min_dw = band_gap_dict['VBM_band_indx_dw'] - num_added_bands 
                        if band_indx_max_dw > band_gap_dict['num_bands']:
                            band_indx_max_dw = band_gap_dict['num_bands'] - 1
                        if band_indx_min_dw < 0:
                            band_indx_min_dw = 0
                        e_max_up = np.max(band_gap_dict['eigs_up'][:, band_indx_max_up])
                        e_min_up = np.min(band_gap_dict['eigs_up'][:, band_indx_min_up])
                        e_max_dw = np.max(band_gap_dict['eigs_dw'][:, band_indx_max_dw])
                        e_min_dw = np.min(band_gap_dict['eigs_dw'][:, band_indx_min_dw])
                        e_max = max(e_max_up, e_max_dw)
                        e_min = min(e_min_up, e_min_dw)
                    if band_struct_ylim in [None, 'None', 'none']:
                        band_struct_ylim = [e_min, e_max]
                    fig_bs_fatband = vasp_plot.plot_bs(infile_path_list = [procar_file_path], xlim = None, ylim = band_struct_ylim, fermi_shift_zero = True, band_list = None,
                        interp_on = False, show_band_data_point = False,
                        band_gap_label = band_gap_label,
                        band_palette_dict = None, band_lty_dict = None, #for single system
                        system_color_list = None, system_lty_list = None, #multiple systems
                        spd_and_site_projections_file_path_list = [projections_file_path], projections_point_size_factor = 1,
                        legend_on = True, plot_fermi_level = True,
                        xtick_direction = 'out', ytick_direction = 'out',
                        line_width = 2.0, font_size = 18, fig_format = 'png', fig_size = [15,10], fig_dpi = 150)
                    ##fig_bs_fatband1 = vasp_plot.plot_bs(infile_path_list = [procar_file_path], xlim = None, ylim = band_struct_ylim, fermi_shift_zero = True, band_list = None,
                    ##    interp_on = False, show_band_data_point = False,
                    ##    band_gap_label = band_gap_label,
                    ##    band_palette_dict = None, band_lty_dict = None, #for single system
                    ##    system_color_list = None, system_lty_list = None, #multiple systems
                    ##    spd_and_site_projections_file_path_list = [projections_file_path], projections_point_size_factor = 1,
                    ##    legend_on = True, plot_fermi_level = True,
                    ##    xtick_direction = 'out', ytick_direction = 'out',
                    ##    line_width = 2.0, font_size = 18, fig_format = 'png', fig_size = [15,10], fig_dpi = 150)
                    parent_path, filename = funcs.file_path_name(fig_bs_fatband)
                    # band structure over the whole energy range
                    tex_fig_label = 'fig:fatband_' + i_task_str.strip('\\').strip('#')
                    task_summary_ctex_str = task_summary_ctex_str + (r'The fat band is shown in Figure ' +  r'\ref{' + tex_fig_label +  '}.' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\begin{figure}[h!]' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\centering' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\includegraphics[width=0.8\linewidth]{' + filename + '}' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\caption{\label{' + tex_fig_label + '}Fat band}' + '\n')
                    task_summary_ctex_str = task_summary_ctex_str + (r'\end{figure}' + '\n')
                    ##parent_path, filename1 = funcs.file_path_name(fig_bs_fatband1)
                    ### band structure in a specific energy range (usually near the Fermi level)
                    ##tex_fig_label1 = 'fig:fatband1_' + i_task_str.strip('\\').strip('#')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'The fat band (in a specific energy range) is shown in Figure ' +  r'\ref{' + tex_fig_label1 +  '}.' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\begin{figure}[h!]' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\centering' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\includegraphics[width=0.8\linewidth]{' + filename1 + '}' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\caption{\label{' + tex_fig_label1 + '}Fat band (in a specific energy range)}' + '\n')
                    ##task_summary_ctex_str = task_summary_ctex_str + (r'\end{figure}' + '\n')
                task_summary_ctex_str = task_summary_ctex_str + (r'\cleardoublepage' + '\n')
            else:
                ##task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{' + 'Band Structure}' + '\n')
                ##task_summary_ctex_str = task_summary_ctex_str + ('None\n')
                ##task_summary_ctex_str = task_summary_ctex_str + (r'\subsection{Fat Band}' + '\n')
                ##task_summary_ctex_str = task_summary_ctex_str + ('None\n')
                task_summary_ctex_str = task_summary_ctex_str + (r'\cleardoublepage' + '\n')
                pass

    # job_status
    #dir_tree_str = funcs.dir_tree(task_dir)
    job_status_str = vasp_tools.job_status(job_parent_dir = task_dir)
    task_summary_ctex_str = task_summary_ctex_str + (r'\section{APPENDIX: Job Status Summary}' + '\n')
    task_summary_ctex_str = task_summary_ctex_str + (r'\begin{landscape}' + '\n')
    task_summary_ctex_str = task_summary_ctex_str + (r'\begin{lstlisting}' + '\n')
    task_summary_ctex_str = task_summary_ctex_str + (job_status_str + '\n')
    task_summary_ctex_str = task_summary_ctex_str + (r'\end{lstlisting}' + '\n')
    task_summary_ctex_str = task_summary_ctex_str + (r'\end{landscape}' + '\n')

        # comparision of recommended and human input
        # comparision of recommended job and human input job
    ########################
    #Write .tex file
    ########################
    #task_summary_ctex_file_name = 'task_summary' + task_dir_name + '.tex'
    task_summary_ctex_file_name = task_dir_name + '.tex'
    task_summary_ctex_file_path = os.path.join(task_summary_dir, task_summary_ctex_file_name)
    task_summary_ctex_str = task_summary_ctex_str + (r'\end{document}' + '\n')
    funcs.write_file(task_summary_ctex_str, task_summary_ctex_file_path, mode = 'w')

    ########################
    # compile tex file
    ########################
    tex_debug_mode = False
    debug_str = funcs.logic_retn_val(tex_debug_mode, '', ' --interaction=batchmode')
    try:
        if tex_debug_mode == False:
            os.system('xelatex' + debug_str + ' -aux-directory=' + task_summary_dir + ' -output-directory=' + task_summary_dir + ' ' + task_summary_ctex_file_path)
            os.system('bibtex ' + debug_str + ' -aux-directory=' + task_summary_dir + ' -output-directory=' + task_summary_dir + ' ' + task_summary_ctex_file_path[:-4] + '.aux')
            os.system('xelatex' + debug_str + ' -aux-directory=' + task_summary_dir + ' -output-directory=' + task_summary_dir + ' ' + task_summary_ctex_file_path)
            os.system('xelatex' + debug_str + ' -aux-directory=' + task_summary_dir + ' -output-directory=' + task_summary_dir + ' ' + task_summary_ctex_file_path)
    except:
        print('WARNING: .tex file ' + task_summary_ctex_file_path + ' compilation failed.')

    return task_summary_ctex_file_path

def gen_inputs(poscar_file_path_list, project_name = 'Untitled', task_type = 'VASP', elmt_potcar_dir = None,
               apply_to_existed_task = True,
              ):
    '''
    generate inputs for a specified structure. The generated inputs are used for further structural optimization.
    the POSCAR file should be provided beforehand.
    For a given structure, to keep the uniqueness of the task, this module should only run once.
    '''
    import time
    import math
    import os
    import random
    import numpy as np
    from .. import funcs
    from .. import default_params
    from ..vasp import vasp_write
    from ..vasp import vasp_read
    from ..vasp import vasp_analyze
    #if params_test == True:
    #get params from params_test_output
    #modify the params according to the params test results

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    projects_dir = os.path.join(output_dir, defaults_dict['projects_dir_name'])
    funcs.mkdir(output_dir)
    funcs.mkdir(projects_dir)

    elmt_potcar_dir = os.path.abspath(elmt_potcar_dir)    

    # build project
    project_dir = os.path.join(projects_dir, project_name)
    funcs.mkdir(project_dir)
    project_dir_name = os.path.split(project_dir)[-1]   

    job_id_file_name = 'job_id.txt'
    # loop over each task 
    for poscar_file_path in poscar_file_path_list:
        poscar_file_path = os.path.abspath(poscar_file_path)
        poscar_dict = vasp_read.read_poscar(poscar_file_path)
        chem_formula_str = ''
        for i_elmt_indx in range(0, len(poscar_dict['elmt_species_arr'])):
            chem_formula_str = chem_formula_str + poscar_dict['elmt_species_arr'][i_elmt_indx] + str(poscar_dict['elmt_num_arr'][i_elmt_indx])

        ##current_time = time.time()
        ##formatted_time = time.strftime('## %Y%m%d%H%M',time.localtime(current_time))
        ##int1 = random.randint(0, 99).zfill(2)
        ##int2 = random.randint(0, 99).zfill(2)
        ##task_id = formatte_time[2:] + str(int1) + str(int2)
        ##task_dir_name = task_id + '_' + chem_formula_str
        ##task_dir = os.path.join(project_dir, task_dir_name)
        task_dir_name = chem_formula_str
        task_dir = os.path.join(project_dir, task_dir_name)
        existed_task_dir_list = [os.path.join(project_dir, x) for x in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, x))]
        if task_dir in existed_task_dir_list:
            if apply_to_existed_task == True:
                pass
            elif apply_to_existed_task == False:
                print('WARNING (from task_manager.gen_inputs): Duplicate chemical composition warning: the task for ' + poscar_file_path + ' is omitted.')
                continue
            ##continue
            #task_dir_name = task_dir_name + '_d1'
        funcs.mkdir(task_dir)

        opt_dir = os.path.join(task_dir, 'opt/')
        scf_dir = os.path.join(task_dir, 'scf/') 
        dos_dir = os.path.join(task_dir, 'dos/')
        bs_dir = os.path.join(task_dir, 'bs/')
        bs_soc_dir = os.path.join(task_dir, 'bs_soc/')
        nscf_dir = os.path.join(task_dir, 'nscf/')
        test_encut_dir = os.path.join(task_dir, 'test_encut/')
        test_kpoints_dir = os.path.join(task_dir, 'test_kpoints/')
        funcs.mkdir(opt_dir)
        funcs.mkdir(scf_dir)
        funcs.mkdir(dos_dir)
        funcs.mkdir(bs_dir)
        funcs.mkdir(bs_soc_dir)
        funcs.mkdir(nscf_dir)
        funcs.mkdir(test_encut_dir)
        funcs.mkdir(test_kpoints_dir)

        task_dir_list = [opt_dir, scf_dir, dos_dir, bs_dir, bs_soc_dir, nscf_dir, test_encut_dir, test_kpoints_dir]
        for i_task_dir in task_dir_list:
            if os.path.exists(i_task_dir) and not os.path.isfile(i_task_dir):
                pass
            else:
                print('WARNING (from task_manager.gen_inputs): ' + i_task_dir + ' does not exist! Cannot be created.')
        ############################
        # opt
        ############################
        outcar_path = os.path.join(opt_dir, 'OUTCAR')
        # when the job is running, or it has been finished, it will be skipped
        job_id_file_path = os.path.join(opt_dir, job_id_file_name)
        if os.path.exists(job_id_file_path):
            with open(job_id_file_path, 'r') as f:
                line = f.readlines()
                try:
                    job_id = int(line[0])
                except:
                    job_id = None
        else:
            job_id = None
        if job_id_exists(job_id) == True:
            # if the job is running, then pass
            pass
        elif check_job_type(opt_dir) == 'VASP' and os.path.exists(outcar_path) and os.path.isfile(outcar_path) and vasp_tools.vasp_job_finished(outcar_path, suppress_warning = True) == True:
            pass
        else:
            # prepare POSCAR
            opt_poscar_file_path = os.path.join(opt_dir, 'POSCAR')
            if not os.path.exists(opt_poscar_file_path) and not os.path.isfile(opt_poscar_file_path):
                funcs.cp(poscar_file_path, opt_poscar_file_path)
            # prepare INCAR
            incar_dict = {}
            incar_dict['SYSTEM']    = 'opt'
            incar_dict['ISTART']    = 0
            incar_dict['ENCUT']     = 500
            incar_dict['EDIFF']     = 1E-06
            incar_dict['EDIFFG']    = -1E-02
            incar_dict['PREC']      = 'Accurate'
            incar_dict['IBRION']    = 2
            incar_dict['NSW']       = 200
            incar_dict['ISIF']      = 3
            incar_dict['POTIM']     = 0.5
            incar_dict['ALGO']      = 'Fast'
            incar_dict['ISMEAR']    = -5
            incar_dict['SIGMA']     = 0.05
            incar_dict['LREAL']     = '.FALSE.'
            incar_dict['ISPIN']     = 2 
            incar_dict['NELM']      = 200
            incar_dict['ICHARG']    = 2
            incar_dict['LWAVE']     = '.FALSE.'
            incar_dict['LCHARG']    = '.FALSE.'
            incar_dict['NPAR']      = 4
            incar_file_path = os.path.join(opt_dir, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 'w')
            # prepare KPOINTS
            vec1_len = np.linalg.norm(poscar_dict['l_arr'][0,:])
            vec2_len = np.linalg.norm(poscar_dict['l_arr'][1,:])
            vec3_len = np.linalg.norm(poscar_dict['l_arr'][2,:])
            num_k1 = int(math.ceil(60 / vec1_len))
            num_k2 = int(math.ceil(60 / vec2_len))
            num_k3 = int(math.ceil(60 / vec3_len))
            # odd number of kpoints
            if (num_k1 % 2) == 0:
                num_k1 = num_k1 + 1
            if (num_k2 % 2) == 0:
                num_k2 = num_k2 + 1
            if (num_k3 % 2) == 0:
                num_k3 = num_k3 + 1
            kpoints_str = ''
            kpoints_str = kpoints_str + 'Gamma centered Monkhorst-Pack grids' + '\n'
            kpoints_str = kpoints_str + '0' + '\n'
            kpoints_str = kpoints_str + 'Gamma' + '\n'
            kpoints_str = kpoints_str + str(num_k1) + ' ' + str(num_k2) + ' ' + str(num_k3) + '\n'
            kpoints_str = kpoints_str + '0 0 0' + '\n'
            kpoints_file_path = os.path.join(opt_dir, 'KPOINTS')
            funcs.write_file(kpoints_str, kpoints_file_path, mode = 'w')
            # prepare POTCAR
            try:
                vasp_write.write_potcar(opt_poscar_file_path, elmt_potcar_dir)
            except:
                print('Cannot write POTCAR file')

            # set default ENCUT according to 1.3 * ENMAX (Value of ENMAX can be found in POTCAR)
            opt_potcar_file_path = os.path.join(task_dir, 'opt', 'POTCAR')
            enmax_val = -9999
            if os.path.exists(opt_potcar_file_path) and os.path.isfile(opt_potcar_file_path):
                if funcs.file_status(opt_potcar_file_path, suppress_warning = True) == 1: 
                    grep_str = funcs.grep('ENMAX', opt_potcar_file_path) 
                    pot_list = grep_str.split('\n')
                    for i_elmt in range(0, len(poscar_dict['elmt_species_arr'])):
                        enmax_temp = float(pot_list[i_elmt].split(';')[0].split('=')[-1])
                        enmax_val = max(enmax_val, enmax_temp)
                    recommended_encut_val = math.ceil(enmax_val * 1.3 / 10) * 10
                    encut_val = recommended_encut_val
                elif funcs.file_status(opt_potcar_file_path, suppress_warning = True) == 2:
                    print('WARNING (from task_manager.gen_inputs): ' + opt_potcar_file_path + ' is empty. Please check your POTCAR file. This problem may caused by the absense of the elemental POTCAR file (e.g. POTCAR.Z), i.e. the designation of the elemental POTCAR directory is not correct or it actually does not exist. Due to this problem, we set ENCUT = 350 (an arbitrary value) for the INCAR of the job directory ' + opt_dir)
                    recommended_encut_val = 350   # an arbitrary ENCUT value
                    encut_val = recommended_encut_val
                elif funcs.file_status(opt_potcar_file_path, suppress_warning = True) == 0:
                    pass

            vasp_write.write_incar(incar_file_path, incar_dict = {'ENCUT': recommended_encut_val}, mode = 's')
            # set default LREAL according to pseudo.F in the VASP source code         
            if incar_dict['LREAL'] in ['.FALSE.', '.F.'] and poscar_dict['n_atoms'] > 16:
                recommended_lreal_val = 'Auto'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            elif incar_dict['LREAL'] in ['.TRUE.', '.T.', 'On', 'O', 'Auto', 'A'] and poscar_dict['n_atoms'] < 8:
                recommended_lreal_val = '.FALSE.'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')

            # set ENCUT and KPOINTS according to the ENCUT and KPOINT test results (if the tests has been done)
            if os.path.exists(test_encut_dir) and not os.path.isfile(test_encut_dir):
                opt_outcar_file_path = os.path.join(task_dir, 'opt', 'OUTCAR')
                i_job_finished = vasp_tools.vasp_job_finished(opt_outcar_file_path, suppress_warning = True)
                i_job_type = check_job_type(os.path.join(task_dir, 'opt'))
                if i_job_type == 'VASP' and i_job_finished != True:
                    recommended_encut_test = get_recommended_test_param(test_encut_dir)
                    funcs.touch(os.path.join(test_encut_dir, 'test.log'))
                    if recommended_encut_test <= encut_val:
                        recommended_encut_test = encut_val
                    with open(os.path.join(test_encut_dir, 'test.log'), 'w') as f:
                        f.write(str(recommended_encut_test))
                    vasp_write.write_incar(incar_file_path, incar_dict = {'ENCUT': recommended_encut_test}, mode = 's')

            if os.path.exists(test_kpoints_dir) and not os.path.isfile(test_kpoints_dir):
                opt_outcar_file_path = os.path.join(task_dir, 'opt', 'OUTCAR')
                i_job_finished = vasp_tools.vasp_job_finished(opt_outcar_file_path, suppress_warning = True)
                i_job_type = check_job_type(os.path.join(task_dir, 'opt'))
                if i_job_type == 'VASP' and i_job_finished != True:
                    recommended_kpoints_test = int(get_recommended_test_param(test_kpoints_dir))
                    funcs.touch(os.path.join(test_kpoints_dir, 'test.log'))
                    kpoints_dict_opt = vasp_read.read_kpoints(os.path.join(opt_dir, 'KPOINTS'))
                    kpoints_val = kpoints_dict_opt['subdivisions_arr'][0]
                    if recommended_kpoints_test <= kpoints_val:
                        recommended_kpoints_test = kpoints_val
                    with open(os.path.join(test_kpoints_dir, 'test.log'), 'w') as f:
                        f.write(str(recommended_kpoints_test))
                    kpoints_dict = vasp_read.read_kpoints(os.path.join(opt_dir, 'KPOINTS'))
                    if kpoints_dict['scheme'] in ['g', 'G', 'm', 'M']:
                        subdivisions_arr = (kpoints_dict['subdivisions_arr'])
                        subdivisions_arr[0] = recommended_kpoints_test   
                        subdivisions_arr[1] = recommended_kpoints_test
                        subdivisions_arr[2] = recommended_kpoints_test
                        kpoints_dict_temp = {}
                        kpoints_dict_temp['subdivisions_arr'] = subdivisions_arr
                        vasp_write.write_kpoints(os.path.join(opt_dir, 'KPOINTS'), kpoints_dict = kpoints_dict_temp, mode = 's')

        # when opt job is running, or it has been finished, the test will be skipped
        job_id_file_path = os.path.join(opt_dir, job_id_file_name)
        if os.path.exists(job_id_file_path):
            with open(job_id_file_path, 'r') as f:
                line = f.readlines()
                try:
                    job_id = int(line[0])
                except:
                    job_id = None
        else:
            job_id = None
        if job_id_exists(job_id) == True:
            # if the job is running, then pass
            pass
        elif check_job_type(opt_dir) == 'VASP' and os.path.exists(outcar_path) and os.path.isfile(outcar_path) and vasp_tools.vasp_job_finished(outcar_path, suppress_warning = True) == True:
            pass
        else:
            ############################
            # ENCUT test
            ############################
            # copy A B
            src_dir = opt_dir
            dst_dir = test_encut_dir
            for i_file_name in ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR']:
                src_file_name = i_file_name
                dst_file_name = i_file_name
                dst_file_path = os.path.join(dst_dir, dst_file_name)
                if os.path.exists(src_dir) and not os.path.isfile(src_dir):
                    src_file_path = os.path.join(opt_dir, src_file_name)
                    if os.path.exists(src_file_path) and os.path.isfile(src_file_path):
                        funcs.cp(src_file_path, dst_file_path)
            # prepare INCAR
            incar_dict_temp = {}
            incar_dict_temp['ENCUT'] = '$i'
            vasp_write.write_incar(os.path.join(test_encut_dir, 'INCAR'), incar_dict = incar_dict_temp, mode = 's')
            #generate ENCUT test list
            incar_dict_opt = vasp_read.read_incar(os.path.join(opt_dir, 'INCAR'))
            encut_val = incar_dict_opt['ENCUT']
            encut_test_val_list = [encut_val]
            encut_val_temp = encut_val
            for indx in range(0, 4):
                encut_val_temp = encut_val_temp - 50
                if encut_val_temp >= 170:
                    encut_test_val_list.append(encut_val_temp)
            encut_val_temp = encut_val
            for indx in range(0, 5):
                encut_val_temp = encut_val_temp + 50
                if encut_val_temp >= 170:
                    encut_test_val_list.append(encut_val_temp)
            encut_test_val_list = sorted(encut_test_val_list)

            params_test(test_encut_dir, variant_name = '$i', variant_value_list = encut_test_val_list)
            # build log file
            #funcs.touch(os.path.join(test_encut_dir, 'test.log'))

            ############################
            # KPOINTS test
            ############################
            # copy A B
            src_dir = opt_dir
            dst_dir = test_kpoints_dir
            for i_file_name in ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR']:
                src_file_name = i_file_name
                dst_file_name = i_file_name
                dst_file_path = os.path.join(dst_dir, dst_file_name)
                if os.path.exists(src_dir) and not os.path.isfile(src_dir):
                    src_file_path = os.path.join(src_dir, src_file_name)
                    if os.path.exists(src_file_path) and os.path.isfile(src_file_path):
                        funcs.cp(src_file_path, dst_file_path)
            # prepare INCAR: get ENCUT from the test_encut job
            # get ENCUT from test ENCUT task
            recommended_encut_test = get_recommended_test_param(test_encut_dir)
            funcs.touch(os.path.join(test_encut_dir, 'test.log'))
            if recommended_encut_test <= encut_val:
                recommended_encut_test = encut_val
            with open(os.path.join(test_encut_dir, 'test.log'), 'w') as f:
                f.write(str(recommended_encut_test))

            test_log_path = os.path.join(test_encut_dir, 'test.log')
            if os.path.exists(test_log_path) and os.path.isfile(test_log_path) and os.path.getsize(test_log_path) > 0:
                with open(test_log_path, 'r') as f:
                    line = f.readlines()
                    recommended_encut = float(funcs.split_line(line = line[0], separator = '=')[-1])

                incar_dict_test_encut = vasp_read.read_incar(os.path.join(test_encut_dir, 'INCAR'))
                incar_dict_temp = {}
                incar_dict_temp['ENCUT'] = recommended_encut
                vasp_write.write_incar(os.path.join(test_kpoints_dir, 'INCAR'), incar_dict = incar_dict_temp, mode = 's')

            # modify the KPOINTS
            poscar_dict = vasp_read.read_poscar(os.path.join(opt_dir, 'POSCAR'))
            vec1_len = np.linalg.norm(poscar_dict['l_arr'][0,:])
            vec2_len = np.linalg.norm(poscar_dict['l_arr'][1,:])
            vec3_len = np.linalg.norm(poscar_dict['l_arr'][2,:])
            num_k1 = int(math.ceil(60 / vec1_len))  
            num_k2 = int(math.ceil(60 / vec2_len))
            num_k3 = int(math.ceil(60 / vec3_len))
            # odd number of kpoints
            if (num_k1 % 2) == 0:
                num_k1 = num_k1 + 1
            if (num_k2 % 2) == 0:
                num_k2 = num_k2 + 1
            if (num_k3 % 2) == 0:
                num_k3 = num_k3 + 1
            kpoints_dict = vasp_read.read_kpoints(os.path.join(opt_dir, 'KPOINTS'))
            if kpoints_dict['scheme'] in ['g', 'G', 'm', 'M']:
                subdivisions_arr = (kpoints_dict['subdivisions_arr'])
                subdivisions_arr[0] = '$i'
                subdivisions_arr[1] = '$i'
                subdivisions_arr[2] = '$i'
                kpoints_dict_temp = {}
                kpoints_dict_temp['subdivisions_arr'] = subdivisions_arr
                vasp_write.write_kpoints(os.path.join(test_kpoints_dir, 'KPOINTS'), kpoints_dict = kpoints_dict_temp, mode = 's')
            #generate test list
            kpoints_dict_opt = vasp_read.read_kpoints(os.path.join(opt_dir, 'KPOINTS'))
            #kpoints_val = kpoints_dict_opt['subdivisions_arr'][0]
            kpoints_val = num_k1
            kpoints_test_val_list = [kpoints_val]
            kpoints_val_temp = kpoints_val
            for indx in range(0, 4):
                kpoints_val_temp = kpoints_val_temp - 2
                if kpoints_val_temp >= 0:
                    kpoints_test_val_list.append(kpoints_val_temp)
            kpoints_val_temp = kpoints_val
            for indx in range(0, 5):
                kpoints_val_temp = kpoints_val_temp + 2
                if kpoints_val_temp >= 0:
                    kpoints_test_val_list.append(kpoints_val_temp)
            kpoints_test_val_list = sorted(kpoints_test_val_list)

            params_test(test_kpoints_dir, variant_name = '$i', variant_value_list = kpoints_test_val_list)
            # build log file
            #funcs.touch(os.path.join(test_kpoints_dir, 'test.log'))

        ############################
        # scf
        ############################
        outcar_path = os.path.join(scf_dir, 'OUTCAR')
        # when the job is running, or it has been finished, it will be skipped
        job_id_file_path = os.path.join(scf_dir, job_id_file_name)
        if os.path.exists(job_id_file_path):
            with open(job_id_file_path, 'r') as f:
                line = f.readlines()
                try:
                    job_id = int(line[0])
                except:
                    job_id = None
        else:
            job_id = None
        if job_id_exists(job_id) == True:
            # if the job is running, then pass
            pass
        elif check_job_type(scf_dir) == 'VASP' and os.path.exists(outcar_path) and os.path.isfile(outcar_path) and vasp_tools.vasp_job_finished(outcar_path, suppress_warning = True) == True:
            pass
        else:
            # copy A B
            src_dir = opt_dir
            dst_dir = scf_dir
            for i_file_name in ['CHGCAR', 'WAVECAR', 'CONTCAR', 'KPOINTS', 'POTCAR']:
                src_file_name = i_file_name
                dst_file_name = i_file_name
                if i_file_name == 'CONTCAR':
                    src_file_name = 'CONTCAR'
                    dst_file_name = 'POSCAR'
                dst_file_path = os.path.join(dst_dir, dst_file_name)
                if os.path.exists(src_dir) and not os.path.isfile(src_dir):
                    src_file_path = os.path.join(src_dir, src_file_name)
                    if os.path.exists(src_file_path) and os.path.isfile(src_file_path):
                        funcs.cp(src_file_path, dst_file_path)
            # prepare INCAR
            incar_dict = {}
            incar_dict['SYSTEM']    = 'scf'
            incar_dict['ISTART']    = 0
            incar_dict['ENCUT']     = 500
            incar_dict['EDIFF']     = 1E-06
            incar_dict['EDIFFG']    = -1E-02
            incar_dict['PREC']      = 'Accurate'
            incar_dict['IBRION']    = -1
            incar_dict['NSW']       = 0
            incar_dict['ISIF']      = 3
            incar_dict['POTIM']     = 0.5
            incar_dict['ALGO']      = 'Fast'
            incar_dict['ISMEAR']    = -5
            incar_dict['SIGMA']     = 0.05
            incar_dict['LREAL']     = '.FALSE.'
            incar_dict['ISPIN']     = 2 
            incar_dict['NELM']      = 200
            incar_dict['ICHARG']    = 2
            incar_dict['LWAVE']     = '.FALSE.'
            incar_dict['LCHARG']    = '.TRUE.'
            incar_dict['NPAR']      = 4
            incar_file_path = os.path.join(scf_dir, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 'w')
       
            # modify ENCUT 
            incar_dict_opt = vasp_read.read_incar(os.path.join(opt_dir, 'INCAR'))
            incar_dict_temp = {}
            incar_dict_temp['ENCUT'] = incar_dict_opt['ENCUT']
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict_temp, mode = 's')
            # set default LREAL according to pseudo.F in the VASP source code         
            if incar_dict['LREAL'] in ['.FALSE.', '.F.'] and poscar_dict['n_atoms'] > 16:
                recommended_lreal_val = 'Auto'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            elif incar_dict['LREAL'] in ['.TRUE.', '.T.', 'On', 'O', 'Auto', 'A'] and poscar_dict['n_atoms'] < 8:
                recommended_lreal_val = '.FALSE.'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            
            # modify the KPOINTS
            kpoints_dict = vasp_read.read_kpoints(os.path.join(opt_dir, 'KPOINTS'))
            if kpoints_dict['scheme'] in ['g', 'G', 'm', 'M']:
                subdivisions_arr = (kpoints_dict['subdivisions_arr'] * 2) - 1 
                kpoints_dict_temp = {}
                kpoints_dict_temp['subdivisions_arr'] = subdivisions_arr
                vasp_write.write_kpoints(os.path.join(scf_dir, 'KPOINTS'), kpoints_dict = kpoints_dict_temp, mode = 's')
        ############################
        # bs
        ############################
        outcar_path = os.path.join(bs_dir, 'OUTCAR')
        # when the job is running, or it has been finished, it will be skipped
        job_id_file_path = os.path.join(bs_dir, job_id_file_name)
        if os.path.exists(job_id_file_path):
            with open(job_id_file_path, 'r') as f:
                line = f.readlines()
                try:
                    job_id = int(line[0])
                except:
                    job_id = None
        else:
            job_id = None
        if job_id_exists(job_id) == True:
            # if the job is running, then pass
            pass
        elif check_job_type(bs_dir) == 'VASP' and os.path.exists(outcar_path) and os.path.isfile(outcar_path) and vasp_tools.vasp_job_finished(outcar_path, suppress_warning = True) == True:
            pass
        else:
            # copy A B
            src_dir = scf_dir
            dst_dir = bs_dir
            for i_file_name in ['CHGCAR', 'WAVECAR', 'CONTCAR', 'KPOINTS', 'POTCAR']:
                src_file_name = i_file_name
                dst_file_name = i_file_name
                if i_file_name == 'CONTCAR':
                    src_file_name = 'CONTCAR'
                    dst_file_name = 'POSCAR'
                dst_file_path = os.path.join(dst_dir, dst_file_name)
                if os.path.exists(src_dir) and not os.path.isfile(src_dir):
                    src_file_path = os.path.join(src_dir, src_file_name)
                    if os.path.exists(src_file_path) and os.path.isfile(src_file_path):
                        funcs.cp(src_file_path, dst_file_path)
            # prepare INCAR
            incar_dict = {}
            incar_dict['SYSTEM']    = 'bs'
            incar_dict['ISTART']    = 0
            incar_dict['ENCUT']     = 500
            incar_dict['EDIFF']     = 1E-06
            incar_dict['EDIFFG']    = -1E-02
            incar_dict['PREC']      = 'Accurate'
            incar_dict['IBRION']    = -1
            incar_dict['NSW']       = 0
            incar_dict['ISIF']      = 3
            incar_dict['POTIM']     = 0.5
            incar_dict['ALGO']      = 'Fast'
            incar_dict['ISMEAR']    = 0
            incar_dict['SIGMA']     = 0.05
            incar_dict['LREAL']     = '.FALSE.'
            incar_dict['ISPIN']     = 2 
            incar_dict['NELM']      = 200
            incar_dict['ICHARG']    = 11
            incar_dict['LWAVE']     = '.FALSE.'
            incar_dict['LCHARG']    = '.TRUE.'
            incar_dict['NPAR']      = 4
            incar_dict['LORBIT']    = 11
            incar_dict['NBANDS']    = 16
            incar_file_path = os.path.join(bs_dir, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 'w')

            # modify ENCUT 
            incar_dict_opt = vasp_read.read_incar(os.path.join(opt_dir, 'INCAR'))
            incar_dict_temp = {}
            incar_dict_temp['ENCUT'] = incar_dict_opt['ENCUT']
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict_temp, mode = 's')
            # set default LREAL according to pseudo.F in the VASP source code         
            if incar_dict['LREAL'] in ['.FALSE.', '.F.'] and poscar_dict['n_atoms'] > 16:
                recommended_lreal_val = 'Auto'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            elif incar_dict['LREAL'] in ['.TRUE.', '.T.', 'On', 'O', 'Auto', 'A'] and poscar_dict['n_atoms'] < 8:
                recommended_lreal_val = '.FALSE.'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            # modify NBANDS
            outcar_path_scf = os.path.join(scf_dir, 'OUTCAR')
            nbands_scf = None
            if os.path.exists(outcar_path_scf) and os.path.isfile(outcar_path_scf): 
                outcar_params_dict_scf = vasp_read.read_outcar(outcar_path_scf)
                nbands_scf = outcar_params_dict_scf['NBANDS']
            if nbands_scf in [None, 'None', 'none']:
                pass
            else:
                incar_dict_temp['NBANDS'] = nbands_scf
                vasp_write.write_incar(incar_file_path, incar_dict = incar_dict_temp, mode = 's')
            # generate KPOINTS (currently, we use pymatgen to generate k-path)
            temp_str = ''
            temp_str = temp_str + "from pymatgen.io.vasp.inputs import Kpoints" + '\n'
            temp_str = temp_str + "from pymatgen.core import Structure" + '\n'
            temp_str = temp_str + "from pymatgen.symmetry.bandstructure import HighSymmKpath" + '\n'
            temp_str = temp_str + "struct = Structure.from_file('POSCAR')" + '\n'
            temp_str = temp_str + "kpath = HighSymmKpath(struct)" + '\n'
            temp_str = temp_str + "kpts = Kpoints.automatic_linemode(divisions = 20, ibz = kpath)" + '\n'
            temp_str = temp_str + "kpts.write_file('KPOINTS')" + '\n'
            with open('gen_kpoints_line_mode_pymatgen.py', 'w') as f:
                f.write(temp_str)

            try:
                poscar_file_path_bs = os.path.join(bs_dir, 'POSCAR')
                kpoints_file_path_bs = os.path.join(bs_dir, 'KPOINTS')
                from pymatgen.io.vasp.inputs import Kpoints
                from pymatgen.core import Structure
                from pymatgen.symmetry.bandstructure import HighSymmKpath
                #from pymatgen import Structure
                from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
                ###################
                # prepare KPOINTS line-mode
                ###################
                struct = Structure.from_file(poscar_file_path_bs)
                kpath = HighSymmKpath(struct)
                kpts = Kpoints.automatic_linemode(divisions = 20, ibz = kpath)
                kpts.write_file(kpoints_file_path_bs)
                ###################
                # prepare primitive cell for band structure calculation
                ###################
                analyzer = SpacegroupAnalyzer(struct)
                prim_cell = analyzer.find_primitive()
                #print(prim_cell)
                prim_cell.to(filename = poscar_file_path_bs)
            except:
                pass

        ############################
        # bs_soc
        ############################
        outcar_path = os.path.join(bs_soc_dir, 'OUTCAR')
        # when the job is running, or it has been finished, it will be skipped
        job_id_file_path = os.path.join(bs_soc_dir, job_id_file_name)
        if os.path.exists(job_id_file_path):
            with open(job_id_file_path, 'r') as f:
                line = f.readlines()
                try:
                    job_id = int(line[0])
                except:
                    job_id = None
        else:
            job_id = None
        if job_id_exists(job_id) == True:
            # if the job is running, then pass
            pass
        elif check_job_type(bs_soc_dir) == 'VASP' and os.path.exists(outcar_path) and os.path.isfile(outcar_path) and vasp_tools.vasp_job_finished(outcar_path, suppress_warning = True) == True:
            pass
        else:
            # copy A B
            src_dir = scf_dir
            dst_dir = bs_soc_dir
            for i_file_name in ['CHGCAR', 'WAVECAR', 'CONTCAR', 'KPOINTS', 'POTCAR']:
                src_file_name = i_file_name
                dst_file_name = i_file_name
                if i_file_name == 'CONTCAR':
                    src_file_name = 'CONTCAR'
                    dst_file_name = 'POSCAR'
                dst_file_path = os.path.join(dst_dir, dst_file_name)
                if os.path.exists(src_dir) and not os.path.isfile(src_dir):
                    src_file_path = os.path.join(src_dir, src_file_name)
                    if os.path.exists(src_file_path) and os.path.isfile(src_file_path):
                        funcs.cp(src_file_path, dst_file_path)
            # prepare INCAR
            incar_dict = {}
            incar_dict['SYSTEM']     = 'bs_soc'
            incar_dict['ISTART']     = 0
            incar_dict['ENCUT']      = 500
            incar_dict['EDIFF']      = 1E-06
            incar_dict['EDIFFG']     = -1E-02
            incar_dict['PREC']       = 'Accurate'
            incar_dict['IBRION']     = -1
            incar_dict['NSW']        = 0
            incar_dict['ISIF']       = 3
            incar_dict['POTIM']      = 0.5
            incar_dict['ALGO']       = 'Fast'
            incar_dict['ISMEAR']     = 0
            incar_dict['SIGMA']      = 0.05
            incar_dict['LREAL']      = '.FALSE.'
            incar_dict['ISPIN']      = 2 
            incar_dict['NELM']       = 200
            incar_dict['ICHARG']     = 11
            incar_dict['LWAVE']      = '.FALSE.'
            incar_dict['LCHARG']     = '.TRUE.'
            incar_dict['NPAR']       = 4
            incar_dict['LORBIT']     = 11
            incar_dict['LSORBIT']    = '.TRUE.'
            incar_dict['LMAXMIX']    = 2
            incar_dict['SAXIS']      = '0 0 1'
            incar_dict['NBANDS']     = 32
            incar_dict['GGA_COMPAT'] = '.FALSE.'
            incar_file_path = os.path.join(bs_soc_dir, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 'w')

            # modify ENCUT 
            incar_dict_opt = vasp_read.read_incar(os.path.join(opt_dir, 'INCAR'))
            incar_dict_temp = {}
            incar_dict_temp['ENCUT'] = incar_dict_opt['ENCUT']
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict_temp, mode = 's')
            # set default LREAL according to pseudo.F in the VASP source code         
            if incar_dict['LREAL'] in ['.FALSE.', '.F.'] and poscar_dict['n_atoms'] > 16:
                recommended_lreal_val = 'Auto'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            elif incar_dict['LREAL'] in ['.TRUE.', '.T.', 'On', 'O', 'Auto', 'A'] and poscar_dict['n_atoms'] < 8:
                recommended_lreal_val = '.FALSE.'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            # modify NBANDS
            outcar_path_bs = os.path.join(bs_dir, 'OUTCAR')
            nbands_bs = None
            if os.path.exists(outcar_path_bs) and os.path.isfile(outcar_path_bs): 
                outcar_params_dict_bs = vasp_read.read_outcar(outcar_path_bs)
                nbands_bs = outcar_params_dict_bs['NBANDS']
            if nbands_bs in [None, 'None', 'none']:
                pass
            else:
                incar_dict_temp['NBANDS'] = int(nbands_bs * 2)
                vasp_write.write_incar(incar_file_path, incar_dict = incar_dict_temp, mode = 's')

        ############################
        # dos
        ############################
        outcar_path = os.path.join(dos_dir, 'OUTCAR')
        # when the job is running, or it has been finished, it will be skipped
        job_id_file_path = os.path.join(dos_dir, job_id_file_name)
        if os.path.exists(job_id_file_path):
            with open(job_id_file_path, 'r') as f:
                line = f.readlines()
                try:
                    job_id = int(line[0])
                except:
                    job_id = None
        else:
            job_id = None
        if job_id_exists(job_id) == True:
            # if the job is running, then pass
            pass
        elif check_job_type(dos_dir) == 'VASP' and os.path.exists(outcar_path) and os.path.isfile(outcar_path) and vasp_tools.vasp_job_finished(outcar_path, suppress_warning = True) == True:
            pass
        else:
            # copy A B
            src_dir = scf_dir
            dst_dir = dos_dir
            for i_file_name in ['CHGCAR', 'WAVECAR', 'CONTCAR', 'KPOINTS', 'POTCAR']:
                src_file_name = i_file_name
                dst_file_name = i_file_name
                if i_file_name == 'CONTCAR':
                    src_file_name = 'CONTCAR'
                    dst_file_name = 'POSCAR'
                dst_file_path = os.path.join(dst_dir, dst_file_name)
                if os.path.exists(src_dir) and not os.path.isfile(src_dir):
                    src_file_path = os.path.join(src_dir, src_file_name)
                    if os.path.exists(src_file_path) and os.path.isfile(src_file_path):
                        funcs.cp(src_file_path, dst_file_path)
            # prepare INCAR
            incar_dict = {}
            incar_dict['SYSTEM']    = 'dos'
            incar_dict['ISTART']    = 1
            incar_dict['ENCUT']     = 500
            incar_dict['EDIFF']     = 1E-06
            incar_dict['EDIFFG']    = -1E-02
            incar_dict['PREC']      = 'Accurate'
            incar_dict['IBRION']    = -1
            incar_dict['NSW']       = 0
            incar_dict['ISIF']      = 3
            incar_dict['POTIM']     = 0.5
            incar_dict['ALGO']      = 'Fast'
            incar_dict['ISMEAR']    = -5
            incar_dict['SIGMA']     = 0.05
            incar_dict['LREAL']     = '.FALSE.'
            incar_dict['ISPIN']     = 2 
            incar_dict['NELM']      = 200
            incar_dict['ICHARG']    = 11
            incar_dict['LWAVE']     = '.FALSE.'
            incar_dict['LCHARG']    = '.FALSE.'
            incar_dict['NPAR']      = 4
            incar_dict['LORBIT']    = 11
            incar_dict['NEDOS']     = 501
            incar_file_path = os.path.join(dos_dir, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 'w')

            # modify ENCUT 
            incar_dict_opt = vasp_read.read_incar(os.path.join(opt_dir, 'INCAR'))
            incar_dict_temp = {}
            incar_dict_temp['ENCUT'] = incar_dict_opt['ENCUT']
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict_temp, mode = 's')
            # set default LREAL according to pseudo.F in the VASP source code         
            if incar_dict['LREAL'] in ['.FALSE.', '.F.'] and poscar_dict['n_atoms'] > 16:
                recommended_lreal_val = 'Auto'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            elif incar_dict['LREAL'] in ['.TRUE.', '.T.', 'On', 'O', 'Auto', 'A'] and poscar_dict['n_atoms'] < 8:
                recommended_lreal_val = '.FALSE.'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')

            # modify the KPOINTS
            kpoints_dict = vasp_read.read_kpoints(os.path.join(opt_dir, 'KPOINTS'))
            if kpoints_dict['scheme'] in ['g', 'G', 'm', 'M']:
                subdivisions_arr = (kpoints_dict['subdivisions_arr'] * 2) + 3 
                kpoints_dict_temp = {}
                kpoints_dict_temp['subdivisions_arr'] = subdivisions_arr
                vasp_write.write_kpoints(os.path.join(dos_dir, 'KPOINTS'), kpoints_dict = kpoints_dict_temp, mode = 's')

        ############################
        # nscf
        ############################
        outcar_path = os.path.join(nscf_dir, 'OUTCAR')
        # when the job is running, or it has been finished, it will be skipped
        job_id_file_path = os.path.join(nscf_dir, job_id_file_name)
        if os.path.exists(job_id_file_path):
            with open(job_id_file_path, 'r') as f:
                line = f.readlines()
                try:
                    job_id = int(line[0])
                except:
                    job_id = None
        else:
            job_id = None
        if job_id_exists(job_id) == True:
            # if the job is running, then pass
            pass
        elif check_job_type(nscf_dir) == 'VASP' and os.path.exists(outcar_path) and os.path.isfile(outcar_path) and vasp_tools.vasp_job_finished(outcar_path, suppress_warning = True) == True:
            pass
        else:
            # copy A B
            src_dir = scf_dir
            dst_dir = nscf_dir
            for i_file_name in ['CHGCAR', 'WAVECAR', 'CONTCAR', 'KPOINTS', 'POTCAR']:
                src_file_name = i_file_name
                dst_file_name = i_file_name
                if i_file_name == 'CONTCAR':
                    src_file_name = 'CONTCAR'
                    dst_file_name = 'POSCAR'
                dst_file_path = os.path.join(dst_dir, dst_file_name)
                if os.path.exists(src_dir) and not os.path.isfile(src_dir):
                    src_file_path = os.path.join(src_dir, src_file_name)
                    if os.path.exists(src_file_path) and os.path.isfile(src_file_path):
                        funcs.cp(src_file_path, dst_file_path)
            # prepare INCAR
            incar_dict = {}
            incar_dict['SYSTEM']    = 'nscf'
            incar_dict['ISTART']    = 1
            incar_dict['ENCUT']     = 500
            incar_dict['EDIFF']     = 1E-06
            incar_dict['EDIFFG']    = -1E-02
            incar_dict['PREC']      = 'Accurate'
            incar_dict['IBRION']    = -1
            incar_dict['NSW']       = 0
            incar_dict['ISIF']      = 3
            incar_dict['POTIM']     = 0.5
            incar_dict['ALGO']      = 'Fast'
            incar_dict['ISMEAR']    = -5
            incar_dict['SIGMA']     = 0.05
            incar_dict['LREAL']     = '.FALSE.'
            incar_dict['ISPIN']     = 2 
            incar_dict['NELM']      = 0
            incar_dict['ICHARG']    = 12
            incar_dict['LWAVE']     = '.FALSE.'
            incar_dict['LCHARG']    = '.TRUE.'
            incar_dict['NPAR']      = 4
            incar_file_path = os.path.join(dos_dir, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 'w')

            # modify ENCUT 
            incar_dict_opt = vasp_read.read_incar(os.path.join(opt_dir, 'INCAR'))
            incar_dict_temp = {}
            incar_dict_temp['ENCUT'] = incar_dict_opt['ENCUT']
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict_temp, mode = 's')
            # set default LREAL according to pseudo.F in the VASP source code         
            if incar_dict['LREAL'] in ['.FALSE.', '.F.'] and poscar_dict['n_atoms'] > 16:
                recommended_lreal_val = 'Auto'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
            elif incar_dict['LREAL'] in ['.TRUE.', '.T.', 'On', 'O', 'Auto', 'A'] and poscar_dict['n_atoms'] < 8:
                recommended_lreal_val = '.FALSE.'
                vasp_write.write_incar(incar_file_path, incar_dict = {'LREAL': recommended_lreal_val}, mode = 's')
    return 0

def params_test(job_dir, variant_name = '$i', variant_value_list = None):
    '''
    For the files in the job_dir, this function substitute the keyword "$i" with the values of the "variant_value_list" and generate test jobs.
    job_dir: the directory of the calculation job.
    variant_value_list: a list of the variant values. 
    '''
    import os
    from .. import funcs
    from .. import default_params
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    ##projects_dir = os.path.join(output_dir, defaults_dict['projects_dir_name'])
    ##funcs.mkdir(projects_dir)
    funcs.mkdir(output_dir)

    job_dir = os.path.abspath(job_dir)
    files_list, file_name_list = funcs.get_files(job_dir, extension = None)

##    temp_list = []
##    for i_file in files_list:
##        if i_file.split('.')[-1] in ['png', 'eps', 'pdf', 'tif', 'tiff', 'jpg', 'jpeg', 'svg', 'svgz', 'pgf', 'ps', 'raw', 'rgba']:
##            # ignore the files that is not ASCII format to avoid file reading error
##            pass
##        else:
##            temp_list.append(i_file)
##    files_list = temp_list

    #test_dir_name = defaults_dict['test_dir_name']
    #funcs.mkdir(os.path.join(output_dir,test_dir_name))
    for i in range(0, len(variant_value_list)):
        #new_dir = os.path.join(output_dir, test_dir_name, str(variant_value_list[i]))
        #funcs.mkdir(new_dir)
        new_dir = os.path.join(job_dir, str(variant_value_list[i]))
        funcs.mkdir(new_dir)
        for j in range(0, len(files_list)):
            dest_file = os.path.join(new_dir, file_name_list[j])
            try:
                funcs.cp(files_list[j], dest_file)
                funcs.replace_file_content(dest_file, variant_name, str(variant_value_list[i]))
            except:
                pass
    funcs.write_log(logfile,
        'task_manager.params_test(\n' +
        '    job_dir = ' + 'r\'' + job_dir + '\'' + ',\n' +
        '    variant_name = ' + '\'' + variant_name + '\'' + ',\n' + 
        '    variant_value_list = ' + str(variant_value_list) + ',\n' 
        '    )\n' +
        '#############\n'
        )
    return 0

def get_recommended_test_param(test_dir, fit_order = 6, plot_fitted_curve = True, suppress_warning = True):
    '''
    test_dir: the directory of the parameter test jobs.
    '''
    import os
    from .. import funcs
    from .. import default_params
    import scipy.optimize as opt
    import numpy as np
    from ..vasp import vasp_read
    from ..vasp import vasp_analyze
    use_plt = True
    try:
        import matplotlib.pyplot as plt
    except:
        use_plt = False

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    ##projects_dir = os.path.join(output_dir, defaults_dict['projects_dir_name'])
    ##funcs.mkdir(projects_dir)
    funcs.mkdir(output_dir)

    general_params_dict = {}
    if use_plt == True:
        general_params_dict['marker_size'] = 10.0
        general_params_dict['line_width'] = 2.0
        general_params_dict['fig_format'] = 'png'
        general_params_dict['fig_dpi'] = 300
        general_params_dict['fig_size'] = [15, 10]
        general_params_dict['font_size'] = 18
        plt.rcParams['figure.figsize'] = (general_params_dict['fig_size'][0],general_params_dict['fig_size'][1])

    test_dir = os.path.abspath(test_dir)
    dir_list = [os.path.join(test_dir, x) for x in sorted([x for x in os.listdir(test_dir) if os.path.isdir(os.path.join(test_dir, x))])]
    dir_name_list = sorted([x for x in os.listdir(test_dir) if os.path.isdir(os.path.join(test_dir, x))])
    energy_list = []
    params_list = []
    for i in range(0, len(dir_list)):
        outcar_path = os.path.join(dir_list[i], 'OUTCAR')
        if funcs.file_status(outcar_path, suppress_warning = True) == 1 and vasp_tools.vasp_job_finished(outcar_path, suppress_warning = True):
            outcar_params_dict = vasp_read.read_outcar(outcar_path)
            energy_list.append(outcar_params_dict['energy(sigma->0)'])
            params_list.append(dir_name_list[i])
        else:
            if suppress_warning == False:
                print('WARNING (from get_recommended_test_param): ' + outcar_path + ' dose not exist or is empty or the VASP job has not been finished!')

    energy_arr = np.array(energy_list, dtype = float)
    params_arr = np.array(params_list, dtype = float)

    recommended_param_val = None 
    if len(params_arr) == 0:
        recommended_param_val = 0 
    else:
        #############################################################
        # determine the optimum parameter
        # fit the data and get the point with zero second derivative
        #############################################################
        #fit_order = 6
        coeffs = np.polyfit(params_arr, energy_arr, fit_order)
        ffit = np.poly1d(coeffs)
        params_dense_arr = np.linspace(min(params_arr), max(params_arr), num = 100, endpoint = True)

        ffit_dense = ffit(params_dense_arr)
        if use_plt == True:
            fig_file = os.path.join(test_dir, 'E_param_relation.' + general_params_dict['fig_format'])
            plt.title('E-param relation', fontsize = general_params_dict['font_size'])
            if plot_fitted_curve == True:
                plt.plot(params_dense_arr, ffit_dense, linestyle = '-', color = 'black', linewidth = general_params_dict['line_width'], markersize = general_params_dict['marker_size'])
            plt.plot(params_arr, energy_arr, linestyle = '', marker = '.', color = 'black', linewidth = general_params_dict['line_width'], markersize = general_params_dict['marker_size'])
            plt.xlabel('param', fontsize = general_params_dict['font_size'])
            plt.ylabel('E (eV)', fontsize = general_params_dict['font_size'])
            plt.xticks(fontsize=general_params_dict['font_size'])
            plt.yticks(fontsize=general_params_dict['font_size'])
            ##plt.show()
            plt.savefig(fig_file, dpi = general_params_dict['fig_dpi'])
            plt.close()
        else:
            pass
        #derivative
        fderiv = ffit.deriv()    #first derivative
        fderiv_n = ffit.deriv(1)  #nth derivative, either n=1 or n=2 would give good results, I choose n=1 here

        fderiv_arr = fderiv(params_arr)
        fderiv_dense = fderiv(params_dense_arr)
        if use_plt == True:
            fig_file = os.path.join(test_dir, 'E_derivative_param_relation.' + general_params_dict['fig_format'])
            plt.title('dE/d(param)-param relation', fontsize = general_params_dict['font_size'])
            if plot_fitted_curve == True:
                plt.plot(params_dense_arr, fderiv_dense, linestyle = '-', color = 'black', linewidth = general_params_dict['line_width'], markersize = general_params_dict['marker_size'])
            plt.plot(params_arr, fderiv_arr, linestyle = '', marker = '.', color = 'black', linewidth = general_params_dict['line_width'], markersize = general_params_dict['marker_size'])
            plt.xlabel('param', fontsize = general_params_dict['font_size'])
            plt.ylabel('E\' (eV)', fontsize = general_params_dict['font_size'])
            plt.xticks(fontsize=general_params_dict['font_size'])
            plt.yticks(fontsize=general_params_dict['font_size'])
            ##plt.show()
            plt.savefig(fig_file, dpi = general_params_dict['fig_dpi'])
            plt.close()
        else:
            pass

        #https://stackoverflow.com/questions/29634217/get-minimum-points-of-numpy-poly1d-curve
        crit = fderiv_n.r
        r_crit = crit[crit.imag==0].real
        ##print('r_crit=', r_crit)
        ##fderiv_n_val = fderiv_n(r_crit)
        ##print('fderiv_n_val=', test)

        #determine the optimum paramter
        sorted_r_crit = sorted(r_crit)
        num_r_crit = len(sorted_r_crit)
        if num_r_crit == 0:
            print('WARNING (from get_recommended_test_param): No optimum parameter is found.')
            recommended_param_val = params_arr[-1]
        elif num_r_crit == 1:
            recommended_param_val = sorted_r_crit[0]
        elif num_r_crit > 1 and num_r_crit < 3:
            # Among the zero second derivative pionts, choose the second smallest piont
            recommended_param_val = sorted_r_crit[1]
        elif num_r_crit >= 3:
            # Among the zero second derivative pionts, choose the second smallest piont
            recommended_param_val = sorted_r_crit[2]
        # find the cloest parameter point to the critical point
        if len(params_arr) > 1:
            optimum_indx = int(np.argsort(np.abs(params_arr - recommended_param_val))[1])
            recommended_param_val = params_arr[optimum_indx]
        elif len(params_arr) == 1:
            optimum_indx = int(np.argsort(np.abs(params_arr - recommended_param_val))[0])
            recommended_param_val = params_arr[optimum_indx]
        ##print('recommended parameter = ', recommended_param_val)

    # write test.log
    funcs.touch(os.path.join(test_dir, 'test.log'))
    with open(os.path.join(test_dir, 'test.log'), 'w') as f:
        f.write(str(recommended_param_val))

    return recommended_param_val

def latt_const_test(job_dir, model_dimension = 2, increment = 0.005, num_points = 5, latt_const_list = None):
    '''
    For the files in the job_dir, this function modifies the POSCAR generates test jobs.
    job_dir: the directory of the calculation job.
    model_dimension: dimension of the model. 1 for chain, 2 for slab, 3 for bulk.
    increment: the increment of the lattice constant, units in percentage/100.
    num_points: the number of scaled systems (including the original one).
    latt_const_list: a list of the lattice constant. If this list is provided, the increment tag is omitted.
    '''
    import os
    import numpy as np
    from .. import funcs
    from .. import default_params
    from ..vasp import vasp_read
    from ..vasp import vasp_write
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    #output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    ##projects_dir = os.path.join(output_dir, defaults_dict['projects_dir_name'])
    ##funcs.mkdir(projects_dir)
    #funcs.mkdir(output_dir)

    job_dir = os.path.abspath(job_dir)
    files_list, file_name_list = funcs.get_files(job_dir, extension = None)

    work_last_dir = os.path.join(job_dir, '..')
    latt_const_test_dir = os.path.join(work_last_dir, 'latt_const_test')
    funcs.mkdir(latt_const_test_dir)

    poscar_file_path = os.path.join(job_dir, 'POSCAR')
    poscar_dict = vasp_read.read_poscar(poscar_file_path)

    if model_dimension == 2:
        indx_list = [0]
        scale_fac_list = []
        if (num_points % 2) != 0:
            up_bound = (num_points - 1) / 2
        elif (num_points % 2) == 0:
            up_bound = num_points / 2
        for i in range(0, up_bound):
            indx_list.append(i+1)
            indx_list.append(-(i+1))
        for i_point in range(0, num_points):
            scale_fac_list.append(1 + increment * indx_list[i_point])
        for i_scale_fac in scale_fac_list:
            box_len_arr = poscar_dict['box_len_arr']
            len_vec1 = np.linalg.norm([box_len_arr[0,0], box_len_arr[0,1], box_len_arr[0,2]]) 
            len_vec2 = np.linalg.norm([box_len_arr[1,0], box_len_arr[1,1], box_len_arr[1,2]]) 
            len_vec3 = np.linalg.norm([box_len_arr[2,0], box_len_arr[2,1], box_len_arr[2,2]])
            #find the axis with the longest length
            len_vec_list = [len_vec1, len_vec2, len_vec3]
            max_len_axis_indx = len_vec_list.index(max(len_vec_list))
            if max_len_axis_indx == 0:
                box_len_arr[1,:] = box_len_arr[1,:] * i_scale_fac
                box_len_arr[2,:] = box_len_arr[2,:] * i_scale_fac
            elif max_len_axis_indx == 1:
                box_len_arr[0,:] = box_len_arr[0,:] * i_scale_fac
                box_len_arr[2,:] = box_len_arr[2,:] * i_scale_fac
            elif max_len_axis_indx == 2:
                box_len_arr[0,:] = box_len_arr[0,:] * i_scale_fac
                box_len_arr[1,:] = box_len_arr[1,:] * i_scale_fac
            vec1_str = ' '.join([str(x) for x in box_len_arr[0,:]]) + '\n'
            vec2_str = ' '.join([str(x) for x in box_len_arr[1,:]]) + '\n'
            vec3_str = ' '.join([str(x) for x in box_len_arr[2,:]]) + '\n'
            poscar_dict['header'][2] = vec1_str
            poscar_dict['header'][3] = vec2_str
            poscar_dict['header'][4] = vec3_str

            if max_len_axis_indx == 0:
                poscar_dict['coord_arr'][:,0] = poscar_dict['coord_arr'][:,0] * i_scale_fac 
            elif max_len_axis_indx == 1:
                poscar_dict['coord_arr'][:,1] = poscar_dict['coord_arr'][:,1] * i_scale_fac 
            elif max_len_axis_indx == 2:
                poscar_dict['coord_arr'][:,2] = poscar_dict['coord_arr'][:,2] * i_scale_fac 

            incar_dict= {}
            incar_dict['ISIF'] = 4
            new_dir = os.path.join(latt_const_test_dir, str(i_scale_fac))
            new_poscar_file_path = os.path.join(new_dir, 'POSCAR')
            new_incar_file_path = os.path.join(new_dir, 'INCAR')
            funcs.mkdir(new_dir)
            for j in range(0, len(files_list)):
                dest_file = os.path.join(new_dir, file_name_list[j])
                try:
                    funcs.cp(files_list[j], dest_file)
                    #funcs.replace_file_content(dest_file, variant_name, str(variant_value_list[i]))
                    vasp_write.write_poscar_with_atom_property(output_poscar_file_path = new_poscar_file_path, poscar_dict = poscar_dict)
                    vasp_write.write_incar(incar_file_path = new_incar_file_path, incar_dict = incar_dict, mode = 's')
                except:
                    pass
    elif model_dimension ==3:
        for i in range(0, len(latt_const_list)):
            #new_dir = os.path.join(output_dir, test_dir_name, str(variant_value_list[i]))
            #funcs.mkdir(new_dir)
            new_dir = os.path.join(job_dir, str(latt_const_list[i]))
            funcs.mkdir(new_dir)
            for j in range(0, len(files_list)):
                dest_file = os.path.join(new_dir, file_name_list[j])
                try:
                    funcs.cp(files_list[j], dest_file)
                    #funcs.replace_file_content(dest_file, variant_name, str(variant_value_list[i]))
                except:
                    pass
    funcs.write_log(logfile,
        'task_manager.latt_const_test(\n' +
        '    job_dir = ' + 'r\'' + job_dir + '\'' + ',\n' +
        '    model_dimension = ' + str(model_dimension) + ',\n' +
        '    increment = ' + str(increment) + ',\n' +
        '    num_points = ' + str(num_points) + ',\n' +
        '    latt_const_list = ' + str(latt_const_list) + ',\n' 
        '    )\n' +
        '#############\n'
        )
    return 0

def gen_submit_script(task_dir, submit_script_path_dict = None, queue_system = 'PBS'):
    '''
    generate submit script for the specified task
    submit_script_path_dict: the keyword 'universal' must exist in submit_script_path_dict, if no submit files are provided for other calculations, the submit file in submit_script_path_dict['universal'] will be used. 
    '''
    import os
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    projects_dir = os.path.join(output_dir, defaults_dict['projects_dir_name'])

    task_flow_list = ['opt', 'test_encut', 'test_kpoints', 'scf', 'bs', 'bs_soc', 'dos', 'nscf']

    if submit_script_path_dict in [None, 'None', 'none']:
        submit_script_path_dict = {}
        submit_script_path_dict['universal'] = os.path.join(project_dir, 'submit.sh')

    if 'universal' not in submit_script_path_dict.keys():
        print('ERROR (from gen_submit_script): The keyword universal must exist in the submit_script_path_dict')
        exit()
    else:
        universal_submit_script_path = os.path.abspath(submit_script_path_dict['universal'])

    for i_task in task_flow_list:
        if i_task not in submit_script_path_dict.keys():
            submit_script_path_dict[i_task] = None
        elif i_task in submit_script_path_dict.keys() and submit_script_path_dict[i_task] in [None, 'None', 'none']:
            submit_script_path_dict[i_task] = os.path.abspath(submit_script_path_dict['universal'])
    
    line_dict = {}    
    # read example submit script
    with open(submit_script_path_dict['universal'], 'r') as f:
        line_dict['universal'] = f.readlines()
    for i_task in task_flow_list:
        if submit_script_path_dict[i_task] in [None, 'None', 'none']:
            line_dict[i_task] = None
        else:
            with open(submit_script_path_dict[i_task], 'r') as f:
                line_dict[i_task] = f.readlines()
        
    # write submit script to each job
    submit_script_name = 'submit.sh'
    run_submit_script_name = 'run_submit.sh'
    i_task_dir = os.path.abspath(task_dir)
    for i_job_indx in range(0, len(task_flow_list)):
        i_job_dir = os.path.join(i_task_dir, task_flow_list[i_job_indx])
        if 'test' in task_flow_list[i_job_indx]:
            subjob_name_list = [x for x in os.listdir(i_job_dir) if os.path.isdir(os.path.join(i_job_dir, x))]
            for i_subjob_indx in range(0, len(subjob_name_list)):
                i_subjob_dir = os.path.join(i_job_dir, subjob_name_list[i_subjob_indx])
                submit_script_file_path = os.path.join(i_subjob_dir, submit_script_name)
                if os.path.exists(submit_script_file_path) == False:
                    fpath, fname = funcs.file_path_name(submit_script_file_path)
                    if os.path.exists(fpath) and os.path.isfile(fpath):
                        os.remove(fpath)
                    funcs.mkdir(fpath)
                    with open(submit_script_file_path, 'w') as f:
                        if line_dict[task_flow_list[i_job_indx]] in [None, 'None', 'none']:
                            f.writelines(line_dict['universal'])
                        else:
                            f.writelines(line_dict[task_flow_list[i_job_indx]])
        else:
            submit_script_file_path = os.path.join(i_job_dir, submit_script_name)
            if os.path.exists(submit_script_file_path) == False:
                fpath, fname = funcs.file_path_name(submit_script_file_path)
                if os.path.exists(fpath) and os.path.isfile(fpath):
                    os.remove(fpat)
                funcs.mkdir(fpath)
                with open(submit_script_file_path, 'w') as f:
                    if line_dict[task_flow_list[i_job_indx]] in [None, 'None', 'none']:
                        f.writelines(line_dict['universal'])
                    else:
                        f.writelines(line_dict[task_flow_list[i_job_indx]])
    return 0 

def job_id_exists(job_id, queue_system = 'PBS'):
    import os
    import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    if isinstance(job_id, int) == True:
        job_id = int(job_id)
        job_id_log_file_name = 'job_id.log'
        job_id_log_file_path = os.path.join(output_dir, job_id_log_file_name)
        funcs.touch(job_id_log_file_path)
        if queue_system == 'PBS':
            os.system('qstat -r ' + ' | grep ' + str(job_id) + ' > ' + job_id_log_file_path)
            ##print('qstat -r ' + ' | grep ' + str(job_id) + ' > ' + job_id_log_file_path)
            with open(job_id_log_file_path, 'r') as f:
                line = f.readlines() 
                num_lines = len(line)
        if num_lines == 0:
            id_exists = False
        elif num_lines == 1:
            id_exists =  True
    else: 
        id_exists = False
    ##print('job ' + str(job_id) + ' exists? ' + str(id_exists))
    return id_exists

def check_job_type(job_dir):
    '''
    decide the task type of the directory.
    '''
    import os
    job_type = 'VASP'

    job_dir = os.path.abspath(job_dir)
    file_list = os.listdir(job_dir)
    if job_type in ['VASP', 'Vasp', 'vasp']:
        job_input_file_name_list = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR']
        for i_file_str in job_input_file_name_list:
            if i_file_str in file_list:
                pass
            else:
                job_type = None
    else:
        job_type = None
    return job_type

