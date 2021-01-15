# -*- coding: utf-8 -*-
def ordered_incar_dict(incar_dict):
    '''
    rearrange the order of the keys of the incar_dict according to vasp_default.incar_default()
    '''  
    from collections import OrderedDict
    from . import vasp_default
    incar_default_dict = vasp_default.incar_default()
    temp_dict = {}
    temp_dict = OrderedDict(temp_dict)
    if incar_dict in [None, 'None', 'none']:
        pass
    elif isinstance(incar_dict, dict):
        for i_key in incar_default_dict.keys():
            for j_key in incar_dict.keys():
                if j_key == i_key:
                    temp_dict[i_key] = incar_dict[i_key]
        incar_dict = temp_dict
    return incar_dict

def check_warning(job_dir):
    '''
    check WARNING of the vasp job
    '''

    #########################################################################
    # WARNING: aliasing errors must be expected set NGX to  34 to avoid them
    #########################################################################


    # -----------------------------------------------------------------------------
    #|                                                                             |
    #|           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           |
    #|           W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!           |
    #|           W    W  A    A  R    R  N N  N  II  N N  N  G       !!!           |
    #|           W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !            |
    #|           WW  WW  A    A  R   R   N   NN  II  N   NN  G    G                |
    #|           W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!           |
    #|                                                                             |
    #|      Your generating k-point grid is not commensurate to the symmetry       |
    #|      of the lattice.  This can cause   slow convergence with respect        |
    #|      to k-points for HF type calculations                                   |
    #|      suggested SOLUTIONS:                                                   |
    #|       ) if not already the case, use automatic k-point generation           |
    #|       ) shift your grid to Gamma (G) (e.g. required for hex or fcc lattice) |
    #|                                                                             |
    # -----------------------------------------------------------------------------

    ############################################
    #WARNING: stress and forces are not correct
    ############################################

    ##################################################################
    #WARNING in EDDRMM: call to ZHEGV failed, returncode =   8  4  6
    ##################################################################

    ##################################################################
    #WARNING: small aliasing (wrap around) errors must be expected
    ##################################################################

    ##################################################################
    #WARNING: random wavefunctions but no delay for mixing, default for NELMDL
    ##################################################################

    pass
    return 0

def check_error(job_dir):
    '''
    check ERROR of the vasp job
    it returns a error list
    NOTE: errors are those prompt that cannot be ignored, and the job has to stop due to an error, which is different from warnings
    '''
    import os
    from .. import funcs
    from . import vasp_analyze

    job_dir = os.path.abspath(job_dir)
    outcar_file_path = os.path.join(job_dir, 'OUTCAR')
    fname_list = os.listdir(job_dir)
    unrelated_file_name_list = ['INCAR', 'POSCAR', 'CONTCAR', 'POTCAR', 'KPOINTS',
                                'DOSCAR', 'EIGNEVAL', 'IBZKPT', 'PCDAT', 'WAVECAR',
                                'PROCAR', 'CHGCAR', 'CHG', 'XDATCAR',
                               ]
    ##for i_item in unrelated_file_name_list:
    ##    if i_item in fname_list:
    ##        fname_list.remove(i_item)
    
    # keywords of warning or error, note that the keywords  must be in the same line
    error_kwds_list_dict = {}

    ########################################################
    # ZBRENT: fatal error in bracketing
    #     please rerun with smaller EDIFF, or copy CONTCAR
    #     to POSCAR and continue
    ########################################################
    error_kwds_list_dict['ZBRENT'] = ['error', 'ZBRENT', 'bracketing', 'rerun', 'copy CONTCAR']

    ###############################################################################
    # ERROR: charge density could not be read from file CHGCAR for ICHARG>10
    ###############################################################################
    error_kwds_list_dict['charge density could not be read ICHARG'] = ['charge density could not be read from file CHGCAR for ICHARG>10']

    #===================================================================================
    #=   BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES
    #=   PID 12196 RUNNING AT e16
    #=   EXIT CODE: 9
    #=   CLEANING UP REMAINING PROCESSES
    #=   YOU CAN IGNORE THE BELOW CLEANUP MESSAGES
    #===================================================================================
    error_kwds_list_dict['BAD TERMINATION CODE: 9'] = ['BAD', 'TERMINATION', 'APPLICATION', 'EXIT CODE: 9']

    ## -----------------------------------------------------------------------------
    ##|                                                                             |
    ##|     EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###     |
    ##|     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
    ##|     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
    ##|     EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       #      |
    ##|     E        R   R    R   R    O     O  R   R                               |
    ##|     E        R    R   R    R   O     O  R    R      ###     ###     ###     |
    ##|     EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###     |
    ##|                                                                             |
    ##|      The linear tetrahedron method can not  be used with the KPOINTS file   |
    ##|      (generation of strings of k-points)                                    |
    ##|                                                                             |
    ## -----------------------------------------------------------------------------
    error_kwds_list_dict['linear tetrahedron strings k-points'] = ['linear tetrahedron method', 'can not', 'strings of k-points']

    ## -----------------------------------------------------------------------------
    ##|                                                                             |
    ##|     EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###     |
    ##|     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
    ##|     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
    ##|     EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       #      |
    ##|     E        R   R    R   R    O     O  R   R                               |
    ##|     E        R    R   R    R   O     O  R    R      ###     ###     ###     |
    ##|     EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###     |
    ##|                                                                             |
    ##|      The linear tetrahedron method can not  be used with the KPOINTS file   |
    ##|      (generation of strings of k-points)                                    |
    ##|                                                                             |
    ## -----------------------------------------------------------------------------
    error_kwds_list_dict['The linear tetrahedron method can not  be used with the KPOINTS file'] = ['The linear tetrahedron method can not  be used with the KPOINTS file']

    ###############################################################################
    # ERROR: charge density could not be read from file CHGCAR for ICHARG>10
    ###############################################################################
    error_kwds_list_dict['charge density could not be read from file CHGCAR for ICHARG>10'] = ['charge density could not be read from file CHGCAR for ICHARG>10']

    ###############################################################################
    #internal error in GENERATE_KPOINTS_TRANS: number of G-vector changed in star
    ###############################################################################
    error_kwds_list_dict['error GENERATE_KPOINTS_TRANS changed in star'] = ['error', 'GENERATE_KPOINTS_TRANS', 'G-vector', 'star']

    ###############################################################################
    #internal error in GENERATE_KPOINTS_TRANS: G vector not found   232    0    0    1    0    0   -1    1 mkpoints_change.F
    ###############################################################################
    error_kwds_list_dict['error GENERATE_KPOINTS_TRANS vector not found'] = ['error', 'GENERATE_KPOINTS_TRANS', 'vector not found']

    ###############################################################################
    # VERY BAD NEWS! internal error in subroutine IBZKPT:
    # Volume of generating cell for k-mesh is zero!       3
    ###############################################################################
    error_kwds_list_dict['Volume of generating cell for k-mesh is zero'] = ['Volume of generating cell for k-mesh is zero']

    # Below are the errors from WANNIER90

    #############################################################
    #wannier90 error: examine the output/error file for details
    #############################################################
    error_kwds_list_dict['wannier90 error: examine the output'] = ['wannier90 error', 'examine the output']

    #################################################################################################
    #param_get_projections: Fewer projections defined than the number of Wannier functions requested
    #################################################################################################
    error_kwds_list_dict['param_get_projections: Fewer projections'] = ['param_get_projections: Fewer projections defined than the number of Wannier functions requested']

    #################################################################################################
    #Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid
    #################################################################################################
    error_kwds_list_dict['Monkhorst-Pack mp_grid'] = ['Monkhorst-Pack', 'setting mp_grid']

    #################################################################################################
    # VERY BAD NEWS! internal error in subroutine INVGRP:
    # inverse of rotation matrix was not found (increase SYMPREC)       3
    #################################################################################################
    error_kwds_list_dict['error in subroutine INVGRP'] = ['error in subroutine INVGRP',]

    #################################################################################################
    # VERY BAD NEWS! internal error in subroutine SGRCON:
    # Found some non-integer element in rotation matrix       3
    #################################################################################################
    error_kwds_list_dict['error in subroutine SGRCON'] = ['error in subroutine SGRCON',]
 
    error_type = None
    if vasp_job_finished(outcar_file_path, suppress_warning = True) == False:
        file_path_list = [os.path.join(job_dir, x) for x in fname_list]
        eo_file_list = []
        latest_eo_file = None
        for i_file in file_path_list:
            if os.path.split(i_file)[-1] == 'OUTCAR':
                continue
            # use keywords to identify the error or output file
            distrk_str = funcs.grep('distrk', i_file, suppress_warning = True)
            distr_str = funcs.grep('distr', i_file, suppress_warning = True)
            if len(distrk_str) != 0 and len(distr_str) != 0:
                eo_file_list.append(i_file)
        # find the latest (most recently generated) error or output file
        if len(eo_file_list) != 0:
            latest_eo_file = max(eo_file_list, key = os.path.getctime)
        # find warning or error from the error or output file
        if latest_eo_file not in [None, 'None', 'none']:
            with open(latest_eo_file, 'r') as f:
                line = f.readlines()
                for i_error_type in error_kwds_list_dict.keys():
                    num_kwd = len(error_kwds_list_dict[i_error_type])
                    ##counter = 0
                    kwd_check_list = [False] * len(error_kwds_list_dict[i_error_type])
                    for i_line_indx in range(0, len(line)):
                        i_line = line[i_line_indx]
                        for i_kwd_indx in range(0, len(error_kwds_list_dict[i_error_type])):
                            i_kwd = error_kwds_list_dict[i_error_type][i_kwd_indx]
                            if i_kwd in i_line:
                                kwd_check_list[i_kwd_indx] = True
                                ##counter = counter + 1
                        ##if counter == num_kwd:
                        if kwd_check_list.count(True) == len(error_kwds_list_dict[i_error_type]):
                            error_type = i_error_type
                            break
    return error_type

def error_handler(task_dir):
    '''
    provide solution for the specific error or warning
    '''
    import os
    from .. import funcs
    from . import vasp_write 

    dir_list = funcs.get_dirs(task_dir)
    for i_dir_path in dir_list: 

        error_type = check_error(i_dir_path)

        if error_type in [None, 'None', 'none', '']:
            pass
        if error_type == 'ZBRENT':
            # SOLUTION: cp CONTCAR POSCAR
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' + 
                  'SOLUTION: cp CONTCAR POSCAR and continue.\n')
            src_file_path = os.path.join(i_dir_path, 'CONTCAR')
            dst_file_path = os.path.join(i_dir_path, 'POSCAR')
            funcs.cp(src_file_path, dst_file_path)
        if error_type == 'BAD TERMINATION CODE: 9':
            # SOLUTION: check the machine memory and computational cost of your job
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' + 
                  'SOLUTION: check the machine memory and computational cost of your job.\n')
            pass
        if error_type == 'error in subroutine SGRCON':
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' + 
                  'SOLUTION: increase SYMPREC (e.g. SYMPREC = 1E-04) and continue.\n')
            incar_dict = {}
            incar_dict['SYMPREC'] = 1E-04
            incar_file_path = os.path.join(i_dir_path, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 's')
        if error_type == 'error in subroutine INVGRP':
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' + 
                  'SOLUTION: increase SYMPREC (e.g. SYMPREC = 1E-04) and continue.\n')
            incar_dict = {}
            incar_dict['SYMPREC'] = 1E-04
            incar_file_path = os.path.join(i_dir_path, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 's')
        if error_type == 'The linear tetrahedron method can not  be used with the KPOINTS file':
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' + 
                  'SOLUTION: set ISMEAR to ISMEAR = 0 and continue.\n')
            incar_dict = {}
            incar_dict['ISMEAR'] = 0
            incar_file_path = os.path.join(i_dir_path, 'INCAR')
            vasp_write.write_incar(incar_file_path, incar_dict = incar_dict, mode = 's')
    return 0

def vasp_job_finished(outcar_file_path, suppress_warning = True):
    '''
    check whether the VASP job has been finished or not
    '''
    import os
    from .. import funcs
    kwd = 'timing'
    retn_val = False
    outcar_file_path = os.path.abspath(outcar_file_path)
    if funcs.file_status(outcar_file_path, suppress_warning = True) == 0:
        if suppress_warning == False:
            print('# WARNING #2012041414 (from vasp_tools.vasp_job_finished): The file ' + str(outcar_file_path) + ' does not exist')
        else:
            pass
    # find the keyword from the bottom of the OUTCAR
    try:
        with open(outcar_file_path, 'r') as f:
            lines = f.readlines()
            num_lines = len(lines)
            num_tail_lines = 300
            num_tail_lines = min(num_lines, num_tail_lines)
            for indx in range(num_lines - 1, num_lines - num_tail_lines -1, -1):
                i_line = lines[indx]
                if kwd in i_line:
                    retn_val = True
                    break
    except:
        pass
    return retn_val

def job_status(job_parent_dir = './'):
    '''
    Check the job status for multiple jobs
    job_parent_dir: This is the parent directory which contains multiple VASP jobs
    '''
    import os
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from . import vasp_analyze
    from .. import default_params
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    job_parent_dir = os.path.abspath(job_parent_dir)

    if funcs.dir_status(job_parent_dir) == 0:
        quit()

    folder_level_list = []
    folder_path_list = []
    folder_path_str_list = []
    for root, dirs, files in os.walk(job_parent_dir):
        dirs.sort()
        level = root.replace(job_parent_dir, '').count(os.sep)
        indent = ' ' * 4 * (level)
        folder_level_list.append(level)
        folder_path_list.append(os.path.abspath(root))
        folder_path_str_list.append('{}{}/'.format(indent, os.path.basename(root)))
     
    max_length_folder_path_str = max([len(x) for x in folder_path_str_list])

    job_status_file_path = os.path.join(output_dir, 'job_status_vasp_' + os.path.split(job_parent_dir)[-1] + '.txt')
    
    temp_str = ''
    temp_str = temp_str + (
        'job_dir' + ' ' * (max_length_folder_path_str - len('job_dir')) + 
        '   status' + ' ' * (11 - len('status')) + 
        'energy_without' + ' ' * (16 - len('energy_without')) + 
        'TOTEN' + ' ' * (13 - len('TOTEN')) + 
        'energy(sigma->0)' + ' ' * (17 - len('energy(sigma->0)')) + 
        'elapsed_time' + ' ' * (14 - len('elapsed_time')) + 
##        'E-fermi' + ' ' * (8 - len('E-fermi')) + '\n' + 
        'Egap' + ' ' * (8 - len('Egap')) + '\n' + 

        ' ' * len('job_dir') + ' ' * (max_length_folder_path_str - len('job_dir')) + 
        ' ' * len('   status') + ' ' * (11 - len('status')) + 
        '_entropy' + ' ' * (16 - len('_entropy')) + 
        ' ' * len('TOTEN') + ' ' * (13 - len('TOTEN')) + 
        ' ' * len('energy(sigma->0)') + ' ' * (17 - len('energy(sigma->0)')) + 
        '(hour)' + ' ' * (14 - len('(hour)')) + 
##        ' ' * len('E-fermi') + ' ' * (8 - len('E-fermi')) + '\n' + 
        ' ' * len('Egap') + ' ' * (8 - len('Egap')) + '\n' + 

        '-' * len('job_dir') + ' ' * (max_length_folder_path_str - len('job_dir')) + 
        '-' * len('   status') + ' ' * (11 - len('status')) + 
        '-' * len('energy_without') + ' ' * (16 - len('energy_without')) + 
        '-' * len('TOTEN') + ' ' * (13 - len('TOTEN')) + 
        '-' * len('energy(sigma->0)') + ' ' * (17 - len('energy(sigma->0)')) + 
        '-' * len('elasped_time') + ' ' * (14 - len('elasped_time')) + 
##        '-' * len('E-fermi') + ' ' * (8 - len('E-fermi')) + '\n' 
        '-' * len('Egap') + ' ' * (8 - len('Egap')) + '\n' 

        )
    for i_folder_indx in range(0, len(folder_path_list)):
        output_exist = False
        # check job status
        job_status = 'unfinished'
        job_status_str = '--'
        energy_without_entropy = None
        for file0 in os.listdir(folder_path_list[i_folder_indx]):
            if 'OUTCAR' == file0:
                output_exist = True
                output_file_path = os.path.join(folder_path_list[i_folder_indx], 'OUTCAR')
                kpoints_file_path = os.path.join(folder_path_list[i_folder_indx], 'KPOINTS')
                eigenval_file_path = os.path.join(folder_path_list[i_folder_indx], 'EIGENVAL')
                if vasp_job_finished(output_file_path) == True:
                    job_status = 'finished'
                    job_status_str = 'finished'
                    outcar_params_dict = vasp_read.read_outcar(output_file_path)
                    energy_without_entropy = outcar_params_dict['energy_without_entropy']
                    toten = outcar_params_dict['TOTEN']
                    energy_sigma_0 = outcar_params_dict['energy(sigma->0)']
                    elapsed_time_hour = convert.time_converter(hour = 0, minute = 0, second = outcar_params_dict['elapsed_time'], unit = 'hour')
                    e_fermi = outcar_params_dict['e_fermi']
                    kpoints_dict = vasp_read.read_kpoints(kpoints_file_path)
                    eigenval_dict = vasp_read.read_eigenval(eigenval_file_path)
                    band_gap_dict = vasp_analyze.get_band_gap(
                        eigenval_or_procar_dict = eigenval_dict, 
                        outcar_params_dict = outcar_params_dict, 
                        kpoints_dict = kpoints_dict, 
                        fermi_shift_zero = True,
                        )
                    band_gap = band_gap_dict['band_gap']
                    if band_gap in [None, 'None', 'none']:
                        band_gap = '--'
                    else:
                        band_gap = '{:.4f}'.format(band_gap)
                continue
        if job_status == 'finished':
            temp_str = temp_str + (folder_path_str_list[i_folder_indx] + ' ' * (max_length_folder_path_str - len(folder_path_str_list[i_folder_indx])) + '   ' + 
                job_status_str + ' ' * (11 - len(job_status_str)) + 
                '{:.4f}'.format(energy_without_entropy) + ' ' * (16 - len('{:.4f}'.format(energy_without_entropy))) + 
                '{:.4f}'.format(toten) + ' ' * (13 - len('{:.4f}'.format(toten))) + 
                '{:.4f}'.format(energy_sigma_0) + ' ' * (17 - len('{:.4f}'.format(energy_sigma_0))) + 
                '{:.4f}'.format(elapsed_time_hour) + ' ' * (14 - len('{:.4f}'.format(elapsed_time_hour))) + 
##                '{:.4f}'.format(e_fermi) + ' ' * (8 - len('{:.4f}'.format(e_fermi))) + 
                band_gap + ' ' * (8 - len(band_gap)) + 
                ' \n')
        else:
            temp_str = temp_str + (folder_path_str_list[i_folder_indx] + ' ' * (max_length_folder_path_str - len(folder_path_str_list[i_folder_indx])) + '   ' + 
                job_status_str + ' ' * (11 - len(job_status_str)) + 
                '--' + ' ' * (16 - len('--')) + 
                '--' + ' ' * (13 - len('--')) + 
                '--' + ' ' * (17 - len('--')) + 
                '--' + ' ' * (14 - len('--')) + 
##                '--' + ' ' * (8  - len('--')) + 
                '--' + ' ' * (8  - len('--')) + 
                ' \n')
    with open(job_status_file_path, 'w') as f:
        f.write(temp_str)
    funcs.write_log(logfile,
        'vasp_tools.job_status(\n' +
        '    job_parent_dir = ' + 'r\'' + job_parent_dir + '\'' + '\n'
        '    )\n' +
        '#############\n'
        )

    return temp_str

def check_file_type(file_path):
    '''
    check the file type
    file_type: string format. It defines the file type, e.g. 'EIGENVAL', 'PROCAR'.
    the separator of each blocks in a line is white space.
    '''
    import os
    from .. import funcs

    file_type = None
    #check the number of string blocks of the header lines
    file_path = os.path.abspath(file_path)
    num_str_blocks_list_eigenval = [4, 5, 1, 1, None, 3, 0, 4]
    num_str_blocks_list_procar = [None, 12, 0, 9, 0, 8, 0]
    def check_num_str_blocks_of_header(file_path, num_str_blocks_list):
        with open(file_path,'r') as f:
            line = f.readlines()
            num_lines = len(line)
            #example of num_str_blocks_list: num_str_blocks_list = [4, 5, 1, 1, 1, 3, 0, 4]
            num_header_lines = len(num_str_blocks_list)
            temp_list = []
            all_passed = True
            if num_lines >= num_header_lines:
                for i_line in range(0, num_header_lines):
                    ##print('file path = ', file_path)
                    ##print('i_line = ', i_line)
                    num_str_blocks = len(funcs.split_line(line = line[i_line], separator = ' '))
                    temp_list.append(num_str_blocks)
                    if num_str_blocks_list[i_line] == None:
                        continue
                    else:
                        if num_str_blocks_list[i_line] != num_str_blocks:
                            all_passed = False
            elif num_lines < num_header_lines:
                all_passed = False
        return all_passed

    str_blocks_logic_val_eigenval = check_num_str_blocks_of_header(file_path, num_str_blocks_list_eigenval)
    str_blocks_logic_val_procar = check_num_str_blocks_of_header(file_path, num_str_blocks_list_procar)
    #check whether the file has the EIGENVAL format or not
    if str_blocks_logic_val_eigenval == True:
        file_type = 'EIGENVAL'                
    #check whether the file has the PROCAR format or not
    if str_blocks_logic_val_procar == True:
        file_type = 'PROCAR'                
        #check the PROCAR file is collinear or noncollinear calculation
        num_header = 3
        with open(file_path,'r') as f:
            line = f.readlines()
            num_lines = len(line)
            num_ions = int(funcs.split_line(line = line[1], separator = ':')[3].split('#')[0])
            temp_len = len(funcs.split_line(line = line[8 + num_ions + 1], separator = ' '))
            if temp_len == 0:
                file_type = 'PROCAR_collinear'
            elif temp_len != 0:
                file_type = 'PROCAR_noncollinear'
    return file_type

def get_latt_param(job_parent_dir = './'):
    '''
    Check the lattice properties (lattice size and angle) for multiple jobs
    job_parent_dir: This is the parent directory which contains multiple VASP jobs
    '''
    import os
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from . import vasp_analyze
    from .. import default_params
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    job_parent_dir = os.path.abspath(job_parent_dir)

    if funcs.dir_status(job_parent_dir) == 0:
        quit()

    folder_level_list = []
    folder_path_list = []
    folder_path_str_list = []
    for root, dirs, files in os.walk(job_parent_dir):
        dirs.sort()
        level = root.replace(job_parent_dir, '').count(os.sep)
        indent = ' ' * 4 * (level)
        folder_level_list.append(level)
        folder_path_list.append(os.path.abspath(root))
        folder_path_str_list.append('{}{}/'.format(indent, os.path.basename(root)))
     
    max_length_folder_path_str = max([len(x) for x in folder_path_str_list])

    job_status_file_path = os.path.join(output_dir, 'latt_param_vasp_' + os.path.split(job_parent_dir)[-1] + '.txt')
    
    temp_str = ''
    temp_str = temp_str + (
        'job_dir' + ' ' * (max_length_folder_path_str - len('job_dir')) + 
        'a(POS)' + ' ' * (13 - len('a(POS)')) + 
        'b(POS)' + ' ' * (13 - len('b(POS)')) + 
        'c(POS)' + ' ' * (13 - len('c(POS)')) + 
        'a(CON)' + ' ' * (13 - len('a(CON)')) + 
        'b(CON)' + ' ' * (13 - len('b(CON)')) + 
        'c(CON)' + ' ' * (13 - len('c(CON)')) + 

        'alpha(POS)' + ' ' * (11 - len('alpha(POS)')) + 
        'beta(POS)' + ' ' * (11 - len('beta(POS)')) + 
        'gamma(POS)' + ' ' * (11 - len('gamma(POS)')) + 
        'alpha(CON)' + ' ' * (11 - len('alpha(CON)')) + 
        'beta(CON)' + ' ' * (11 - len('beta(CON)')) + 
        'gamma(CON)' + ' ' * (11 - len('gamma(CON)')) + 

        'volume(POS)' + ' ' * (13 - len('volume(POS)')) + 
        'volume(CON)' + ' ' * (13 - len('volume(CON)')) + '\n' +  

        '-' * len('job_dir') + ' ' * (max_length_folder_path_str - len('job_dir')) + 
        '-' * len('a(POS)') + ' ' * (13 - len('a(POS)')) + 
        '-' * len('b(POS)') + ' ' * (13 - len('b(POS)')) + 
        '-' * len('c(POS)') + ' ' * (13 - len('c(POS)')) + 
        '-' * len('a(CON)') + ' ' * (13 - len('a(CON)')) + 
        '-' * len('b(CON)') + ' ' * (13 - len('b(CON)')) + 
        '-' * len('c(CON)') + ' ' * (13 - len('c(CON)')) + 

        '-' * len('alpha(POS)') + ' ' * (11 - len('alpha(POS)')) + 
        '-' * len('beta(POS)') + ' ' * (11 - len('beta(POS)')) + 
        '-' * len('gamma(POS)') + ' ' * (11 - len('gamma(POS)')) + 
        '-' * len('alpha(CON)') + ' ' * (11 - len('alpha(CON)')) + 
        '-' * len('beta(CON)') + ' ' * (11 - len('beta(CON)')) + 
        '-' * len('gamma(CON)') + ' ' * (11 - len('gamma(CON)')) + 

        '-' * len('volume(POS)') + ' ' * (13 - len('volume(POS)')) + 
        '-' * len('volume(CON)') + ' ' * (13 - len('volume(CON)')) + '\n' 

        )
    for i_folder_indx in range(0, len(folder_path_list)):
        output_exist = False
        # check job status
        job_status = 'unfinished'
        job_status_str = '--'
        energy_without_entropy = None
        for file0 in os.listdir(folder_path_list[i_folder_indx]):
            if 'OUTCAR' == file0:
                output_exist = True
                output_file_path = os.path.join(folder_path_list[i_folder_indx], 'OUTCAR')
                poscar_file_path = os.path.join(folder_path_list[i_folder_indx], 'POSCAR')
                contcar_file_path = os.path.join(folder_path_list[i_folder_indx], 'CONTCAR')
                if vasp_job_finished(output_file_path) == True:
                    job_status = 'finished'
                    job_status_str = 'finished'
                    poscar_dict = vasp_read.read_poscar(poscar_file_path)
                    contcar_dict = vasp_read.read_poscar(contcar_file_path)
                    len_a_poscar = poscar_dict['len_vec_a']
                    len_b_poscar = poscar_dict['len_vec_b']
                    len_c_poscar = poscar_dict['len_vec_c']
                    len_a_contcar = contcar_dict['len_vec_a']
                    len_b_contcar = contcar_dict['len_vec_b']
                    len_c_contcar = contcar_dict['len_vec_c']
                    angle_alpha_degree_poscar = poscar_dict['angle_alpha_degree']
                    angle_beta_degree_poscar = poscar_dict['angle_beta_degree']
                    angle_gamma_degree_poscar = poscar_dict['angle_gamma_degree']

                    angle_alpha_degree_contcar = contcar_dict['angle_alpha_degree']
                    angle_beta_degree_contcar = contcar_dict['angle_beta_degree']
                    angle_gamma_degree_contcar = contcar_dict['angle_gamma_degree']

                    box_volume_poscar = poscar_dict['box_volume']
                    box_volume_contcar = contcar_dict['box_volume']
                continue
        if job_status == 'finished':
            temp_str = temp_str + (folder_path_str_list[i_folder_indx] + ' ' * (max_length_folder_path_str - len(folder_path_str_list[i_folder_indx])) + '   ' + 
                '{:.4f}'.format(len_a_poscar) + ' ' * (13 - len('{:.4f}'.format(len_a_poscar))) + 
                '{:.4f}'.format(len_b_poscar) + ' ' * (13 - len('{:.4f}'.format(len_b_poscar))) + 
                '{:.4f}'.format(len_c_poscar) + ' ' * (13 - len('{:.4f}'.format(len_c_poscar))) + 
                '{:.4f}'.format(len_a_contcar) + ' ' * (13 - len('{:.4f}'.format(len_a_contcar))) + 
                '{:.4f}'.format(len_b_contcar) + ' ' * (13 - len('{:.4f}'.format(len_b_contcar))) + 
                '{:.4f}'.format(len_c_contcar) + ' ' * (13 - len('{:.4f}'.format(len_c_contcar))) + 

                '{:.2f}'.format(angle_alpha_degree_poscar) + ' ' * (11 - len('{:.2f}'.format(angle_alpha_degree_poscar))) + 
                '{:.2f}'.format(angle_beta_degree_poscar) + ' ' * (11 - len('{:.2f}'.format(angle_beta_degree_poscar))) + 
                '{:.2f}'.format(angle_gamma_degree_poscar) + ' ' * (11 - len('{:.2f}'.format(angle_gamma_degree_poscar))) + 
                '{:.2f}'.format(angle_alpha_degree_contcar) + ' ' * (11 - len('{:.2f}'.format(angle_alpha_degree_contcar))) + 
                '{:.2f}'.format(angle_beta_degree_contcar) + ' ' * (11 - len('{:.2f}'.format(angle_beta_degree_contcar))) + 
                '{:.2f}'.format(angle_gamma_degree_contcar) + ' ' * (11 - len('{:.2f}'.format(angle_gamma_degree_contcar))) + 

                '{:.4f}'.format(box_volume_poscar) + ' ' * (13 - len('{:.4f}'.format(box_volume_poscar))) + 
                '{:.4f}'.format(box_volume_contcar) + ' ' * (13 - len('{:.4f}'.format(box_volume_contcar))) + 

                ' \n')
        else:
            temp_str = temp_str + (folder_path_str_list[i_folder_indx] + ' ' * (max_length_folder_path_str - len(folder_path_str_list[i_folder_indx])) + '   ' + 
                '--' + ' ' * (13 - len('--')) + 
                '--' + ' ' * (13 - len('--')) + 
                '--' + ' ' * (13 - len('--')) + 
                '--' + ' ' * (13 - len('--')) + 
                '--' + ' ' * (13 - len('--')) + 
                '--' + ' ' * (13 - len('--')) + 

                '--' + ' ' * (11 - len('--')) + 
                '--' + ' ' * (11 - len('--')) + 
                '--' + ' ' * (11 - len('--')) + 
                '--' + ' ' * (11 - len('--')) + 
                '--' + ' ' * (11 - len('--')) + 
                '--' + ' ' * (11 - len('--')) + 

                '--' + ' ' * (13 - len('--')) + 
                '--' + ' ' * (13 - len('--')) + 
                ' \n')
    with open(job_status_file_path, 'w') as f:
        f.write(temp_str)
    funcs.write_log(logfile,
        'vasp_tools.get_latt_param(\n' +
        '    job_parent_dir = ' + 'r\'' + job_parent_dir + '\'' + '\n'
        '    )\n' +
        '#############\n'
        )

    return temp_str

def convert_coord_system(poscar_file_path, mode = 'direct2cartesian'):
    '''
    Convert from Direct coordinate to Cartesian coordinate
    mode: 'direct2cartesian' or 'cartesian2direct'
    '''
    import os
    from . import vasp_read
    from . import vasp_write
    from .. import default_params

    defaults_dict = default_params.default_params()
    poscar_file_path = os.path.abspath(poscar_file_path)
    poscar_dict = vasp_read.read_poscar(poscar_file_path) 
    fpath, fname = os.path.split(poscar_file_path)

    if mode == 'direct2cartesian':
        coord_system = 'Cartesian'
        poscar_suffix = '_car'
    elif mode == 'cartesian2direct':
        coord_system = 'Direct'
        poscar_suffix = '_dir'
    poscar_file_path_cartesian = os.path.join(fpath, fname + poscar_suffix)

    vasp_write.write_poscar(output_poscar_file_path = poscar_file_path_cartesian, poscar_dict = poscar_dict, coord_system = coord_system)
    return 0
