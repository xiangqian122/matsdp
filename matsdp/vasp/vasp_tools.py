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

    ##################################################################
    #WARNING in EDDRMM: call to ZHEGV failed, returncode =   6  3 27
    ##################################################################

    #################################################################################################
    # RNING in EDDRMM: call to ZHEGV failed, returncode =   6  3 31
    # num prob  num prob  num prob  num prob  num prob  num prob  num prob  num prob  num prob  num prob
    #################################################################################################

    #################################################################################################
    #  -----------------------------------------------------------------------------
    # |                                                                             |
    # |           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           |
    # |           W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!           |
    # |           W    W  A    A  R    R  N N  N  II  N N  N  G       !!!           |
    # |           W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !            |
    # |           WW  WW  A    A  R   R   N   NN  II  N   NN  G    G                |
    # |           W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!           |
    # |                                                                             |
    # |      For optimal performance we recommend to set                            |
    # |        NCORE= 4 - approx SQRT( number of cores)                             |
    # |      NCORE specifies how many cores store one orbital (NPAR=cpu/NCORE).     |
    # |      This setting can  greatly improve the performance of VASP for DFT.     |
    # |      The default, NPAR=number of cores might be grossly inefficient         |
    # |      on modern multi-core architectures or massively parallel machines.     |
    # |      Do your own testing !!!!                                               |
    # |      Unfortunately you need to use the default for GW and RPA calculations. |
    # |      (for HF NCORE is supported but not extensively tested yet)             |
    # |                                                                             |
    #  -----------------------------------------------------------------------------
    #################################################################################################


    #################################################################################################
    #  -----------------------------------------------------------------------------
    # |                                                                             |
    # |           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           |
    # |           W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!           |
    # |           W    W  A    A  R    R  N N  N  II  N N  N  G       !!!           |
    # |           W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !            |
    # |           WW  WW  A    A  R   R   N   NN  II  N   NN  G    G                |
    # |           W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!           |
    # |                                                                             |
    # |      Your generating k-point grid is not commensurate to the symmetry       |
    # |      of the lattice.  This can cause   slow convergence with respect        |
    # |      to k-points for HF type calculations                                   |
    # |      suggested SOLUTIONS:                                                   |
    # |       ) if not already the case, use automatic k-point generation           |
    # |       ) shift your grid to Gamma (G) (e.g. required for hex or fcc lattice) |
    # |                                                                             |
    #  -----------------------------------------------------------------------------
    #################################################################################################

    #################################################################################################
    #  -----------------------------------------------------------------------------
    # |                                                                             |
    # |           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           |
    # |           W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!           |
    # |           W    W  A    A  R    R  N N  N  II  N N  N  G       !!!           |
    # |           W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !            |
    # |           WW  WW  A    A  R   R   N   NN  II  N   NN  G    G                |
    # |           W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!           |
    # |                                                                             |
    # |      Your highest band is occupied at some k-points! Unless you are         |
    # |      performing a calculation for an insulator or semiconductor, without    |
    # |      unoccupied bands, you have included TOO FEW BANDS!! Please increase    |
    # |      the parameter NBANDS in file INCAR to ensure that the highest band     |
    # |      is unoccupied at all k-points. It is always recommended to             |
    # |      include a few unoccupied bands to accelerate the convergence of        |
    # |      molecular dynamics runs (even for insulators or semiconductors).       |
    # |      Because the presence of unoccupied bands improves wavefunction         |
    # |      prediction, and helps to suppress 'band-crossings.'                    |
    # |                                                                             |
    #  -----------------------------------------------------------------------------
    #################################################################################################



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
    error_kwds_list_dict['The linear tetrahedron method can not  be used with the KPOINTS file'] = ['The linear tetrahedron method can not  be used with the KPOINTS file']

    ###############################################################################
    #  WARNING: dimensions on CHGCAR file are different
    #  ERROR: charge density could not be read from file CHGCAR for ICHARG>10
    #  ERROR: charge density could not be read from file CHGCAR for ICHARG>10
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

    ###############################################################################
    # ERROR: charge density could not be read from file CHGCAR for ICHARG>10
    ###############################################################################
    error_kwds_list_dict['ERROR: charge density could not be read from file CHGCAR for ICHARG>10'] = ['ERROR: charge density could not be read from file CHGCAR for ICHARG>10']

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

    #################################################################################################
    # Error reading item 'MAGMOM' from file INCAR.
    # Error code was IERR=0 ... . Found N=  114 data.
    #################################################################################################
    error_kwds_list_dict['Error reading item MAGMOM from file INCAR'] = ['Error code was IERR=0', 'Found N=', 'data']

    #################################################################################################
    #  -----------------------------------------------------------------------------
    # |                                                                             |
    # |     EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###     |
    # |     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
    # |     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
    # |     EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       #      |
    # |     E        R   R    R   R    O     O  R   R                               |
    # |     E        R    R   R    R   O     O  R    R      ###     ###     ###     |
    # |     EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###     |
    # |                                                                             |
    # |      The number of bands is not sufficient to hold all electrons.           |
    # |      I am sorry, but this calculation does  not make sense.                 |
    # |      If you want to calculate only selected states you need to set          |
    # |      either EFERMI oder EREF. EREF specifies around which energy states     |
    # |      are supposed to be calculated.                                         |
    # |      I found NBANDS    =      264       NELECT  =     180                   |
    # |                                                                             |
    #  -----------------------------------------------------------------------------
    # or
    #  -----------------------------------------------------------------------------
    # |                                                                             |
    # |     EEEEEEE  RRRRRR   RRRRRR   OOOOOOO  RRRRRR      ###     ###     ###     |
    # |     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
    # |     E        R     R  R     R  O     O  R     R     ###     ###     ###     |
    # |     EEEEE    RRRRRR   RRRRRR   O     O  RRRRRR       #       #       #      |
    # |     E        R   R    R   R    O     O  R   R                               |
    # |     E        R    R   R    R   O     O  R    R      ###     ###     ###     |
    # |     EEEEEEE  R     R  R     R  OOOOOOO  R     R     ###     ###     ###     |
    # |                                                                             |
    # |      The number of bands is not sufficient to hold all electrons.           |
    # |      I am sorry, but this calculation does  not make sense.                 |
    # |      If you want to calculate only selected states you need to set          |
    # |      either EFERMI oder EREF. EREF specifies around which energy states     |
    # |      are supposed to be calculated.                                         |
    # |      I found NBANDS    =      160       NELECT  =      32                   |
    # |                                                                             |
    # |      ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ..., BYE!!! <----       |
    # |                                                                             |
    #  -----------------------------------------------------------------------------
    #################################################################################################
    error_kwds_list_dict['Number of electrons'] = ['The number of bands is not sufficient to hold all electrons',]

    #################################################################################################
    # ERROR: non collinear calculations require that VASP is compiled
    # without the flag -DNGXhalf and -DNGZhalf
    #################################################################################################
    error_kwds_list_dict['non collinear calculations require that VASP is compiled'] = ['non collinear calculations require that VASP is compiled',]

    ## The following errors results from the machine, not VASP program

    #################################################################################################
    # forrtl: severe (174): SIGSEGV, segmentation fault occurred
    #################################################################################################

    #################################################################################################
    # error in IBZKPT_HF: two k-points are equivalent           1          60
    # this will cause problems in the HF routine
    #################################################################################################

    #################################################################################################
    #################################################################################################

    # read errors 
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
                  'SOLUTION: Manually check the machine memory and computational cost of your job.\n')
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
        if error_type == 'ERROR: charge density could not be read from file CHGCAR for ICHARG>10':
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' + 
                  'This may be caused by an empty or invalid CHGCAR file\n' +
                  'SOLUTION: Manually copy the CHGCAR file to the directory and continue.\n')
        if error_type == 'Error reading item MAGMOM from file INCAR':
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' + 
                  'This may be caused by an incorrect setting of the MAGMOM tag in the INCAR file (N<NMAGMOM).\n' +
                  'SOLUTION: Manually check the MAGMOM tag of INCAR file and continue.\n')
        if error_type == 'Number of electrons':
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' +
                  'This may be caused by an incorrect setting of the NBANDS tag in the INCAR file.\n' +
                  'e.g. NELECT>NBANDS*2 for collinear calculations or NELECT*2>NBANDS*2 for noncollinear calculations.\n' + 
                  'SOLUTION: Manually check the NBANDS tag of INCAR file and continue.\n' +
                  'Please set NBANDS>NELECT/2 for collinear calculations or NBANDS*2>NELECT*2 for noncollinear calculations.\n'
                  )
                  
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

def job_status(job_parent_dir = './', write_abspath = False, write_band_gap = False):
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
    import multiprocessing
    from multiprocessing import Pool
    import time
    import math

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
        if write_abspath == False:
            folder_path_str_list.append('{}{}/'.format(indent, os.path.basename(root)))
        elif write_abspath == True:
            folder_path_str_list.append(root.replace(os.path.split(job_parent_dir)[0], ''))
     
    max_length_folder_path_str = max([len(x) for x in folder_path_str_list])

    job_status_file_path = os.path.join(output_dir, 'job_status_vasp_' + os.path.split(job_parent_dir)[-1] + '.txt')
    
    temp_str = ''
    # prepare header
    temp_str = temp_str + (
        'job_dir' + ' ' * (max_length_folder_path_str - len('job_dir')) + 
        '   status' + ' ' * (13 - len('status')) + 
        'energy_without' + ' ' * (16 - len('energy_without')) + 
        'TOTEN' + ' ' * (13 - len('TOTEN')) + 
        'energy(sigma->0)' + ' ' * (17 - len('energy(sigma->0)')) + 
        'elapsed_time' + ' ' * (14 - len('elapsed_time'))
        )
##        'E-fermi' + ' ' * (8 - len('E-fermi')) + '\n' + 
    if write_band_gap == True:
        temp_str = temp_str +  'Egap' + ' ' * (8 - len('Egap'))
    temp_str = temp_str + '\n'
    temp_str = temp_str + (
        ' ' * len('job_dir') + ' ' * (max_length_folder_path_str - len('job_dir')) + 
        ' ' * len('   status') + ' ' * (13 - len('status')) + 
        '_entropy' + ' ' * (16 - len('_entropy')) + 
        ' ' * len('TOTEN') + ' ' * (13 - len('TOTEN')) + 
        ' ' * len('energy(sigma->0)') + ' ' * (17 - len('energy(sigma->0)')) + 
        '(hour)' + ' ' * (14 - len('(hour)')) + 
##        ' ' * len('E-fermi') + ' ' * (8 - len('E-fermi')) + '\n' + 
        ' ' * len('Egap') + ' ' * (8 - len('Egap')) + '\n' + 

        '-' * len('job_dir') + ' ' * (max_length_folder_path_str - len('job_dir')) + 
        '-' * len('   status') + ' ' * (13 - len('status')) + 
        '-' * len('energy_without') + ' ' * (16 - len('energy_without')) + 
        '-' * len('TOTEN') + ' ' * (13 - len('TOTEN')) + 
        '-' * len('energy(sigma->0)') + ' ' * (17 - len('energy(sigma->0)')) + 
        '-' * len('elasped_time') + ' ' * (14 - len('elasped_time')) 
##        '-' * len('E-fermi') + ' ' * (8 - len('E-fermi')) + '\n' 
        )
    if write_band_gap == True:
        temp_str = temp_str + '-' * len('Egap') + ' ' * (8 - len('Egap')) 
    temp_str = temp_str + '\n'

    # check each job
    args = []
    for i_indx in range(0, len(folder_path_list)):
        #primitive_structure = data[i_indx]['cif']
        #material_id_list.append(data[i_indx]['material_id'])
        args.append((folder_path_list, folder_path_str_list, i_indx, job_status_file_path, write_band_gap))
    cores = math.ceil(multiprocessing.cpu_count() * 3 / 4)

    parallel_run = True
    # Check whether there is a "Permission denied" error when running in parallel 
    # PermissionError: [Errno 13] Permission denied
    try:
        with Pool(cores) as p:
            pass
    except:
        print('# WARNING (2104082043): Running in serial mode.')
        parallel_run = False

    results = ''
    if parallel_run == True:
        with Pool(cores) as p:
            results = p.starmap(single_job_status, args)
    elif parallel_run == False:
        for i_indx in range(0, len(folder_path_list)):
            # use unpacked arguments in the tuple
            results = results + single_job_status(*args[i_indx])

    for i_str in results:
        temp_str = temp_str + i_str


    with open(job_status_file_path, 'w') as f:
        f.write(temp_str)
    funcs.write_log(logfile,
        'vasp_tools.job_status(\n' +
        '    job_parent_dir = ' + 'r\'' + job_parent_dir + '\'' + '\n'
        '    )\n' +
        '#############\n'
        )

    return temp_str

def single_job_status(folder_path_list, folder_path_str_list, i_folder_indx, job_status_file_path, write_band_gap):
    import os
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from . import vasp_analyze
    import numpy as np

    data_arr = np.array([None] * len(folder_path_list) * 6)
    data_arr.shape = len(folder_path_list), 6

    max_length_folder_path_str = max([len(x) for x in folder_path_str_list])
#for i_folder_indx in range(0, len(folder_path_list)):
    output_exist = False
    # check job status
    job_status = 'unfinished'
    job_status_str = 'undone'
    if funcs.check_job_type(folder_path_list[i_folder_indx])['VASP'] == False:
        job_status = 'not_vasp_job'
        job_status_str = 'not_VASP_job'
        #job_status_str = '--'
    energy_without_entropy = None
    temp_str = ''
    for file0 in os.listdir(folder_path_list[i_folder_indx]):
        if 'OUTCAR' == file0:
            output_exist = True
            output_file_path = os.path.join(folder_path_list[i_folder_indx], 'OUTCAR')
            oszicar_file_path = os.path.join(folder_path_list[i_folder_indx], 'OSZICAR')
            kpoints_file_path = os.path.join(folder_path_list[i_folder_indx], 'KPOINTS')
            eigenval_file_path = os.path.join(folder_path_list[i_folder_indx], 'EIGENVAL')
            if vasp_job_finished(output_file_path) == True:
                job_status = 'finished'
                job_status_str = 'finished'
                outcar_params_dict = vasp_read.read_outcar(output_file_path)
                elapsed_time_hour = convert.time_converter(hour = 0, minute = 0, second = outcar_params_dict['elapsed_time'], unit = 'hour')
                kpoints_dict = vasp_read.read_kpoints(kpoints_file_path)
                oszicar_dict = vasp_read.read_oszicar(oszicar_file_path)

                # check if the parameters in the outcar_params_dict are valid
                valid_outcar_par = True
                temp_outcar_par_list = ['file_path', 'NELM', 'energy_without_entropy', 'TOTEN', 'energy(sigma->0)', 'elapsed_time']
                for i_par in temp_outcar_par_list:
                    if funcs.variables_equal(outcar_params_dict[i_par], outcar_params_dict['initial_value_dict'][i_par]):
                        valid_outcar_par = False
                        print('WARNING #2103301239 (from vasp_tools). outcar_params_dict[\'' + str(i_par) + '\'] has an invalid value, please check the OUTCAR file ' + outcar_params_dict['file_path'] + '\n')

                valid_oszicar_par = True
                temp_oszicar_par_list = ['file_path', 'num_electronic_steps_list', 'num_ionic_steps']
                for i_par in temp_oszicar_par_list:
                    if funcs.variables_equal(oszicar_dict[i_par], oszicar_dict['initial_value_dict'][i_par]):
                        valid_oszicar_par = False
                        print('WARNING #2103311429 (from vasp_tools). oszicar_dict[\'' + str(i_par) + '\'] has an invalid value, please check the OSZICAR file ' + oszicar_dict['file_path'] + '\n')

                ##if outcar_params_dict['NELM'] in oszicar_dict['num_electronic_steps_list'] or oszicar_dict['num_ionic_steps'] == outcar_params_dict['NSW']:
                ##    job_status_str = 'finished*'
                if outcar_params_dict['NELM'] in oszicar_dict['num_electronic_steps_list']:
                    job_status_str = 'NELM_reached'
                if outcar_params_dict['NSW'] == oszicar_dict['num_ionic_steps']:
                    job_status_str = 'NSW_reached'

                valid_band_gap = True
                if write_band_gap == True:
                    eigenval_dict = vasp_read.read_eigenval(eigenval_file_path)
                    if eigenval_dict['read_status'] == 1:
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

                    # check if the parameters in the band_gaps_dict are valid
                    #if band_gap_dict['band_gap'] is None:
                    #    valid_band_gap = False
                    #    print('WARNING #2103301306 (from vasp_tools). band_gap_dict[\'band_gap\'] has an invalid value, please check the EIGENVAL/PROCAR file ' + band_gap_dict['file_path'] + '\n')
            continue
  
    if job_status == 'finished' and valid_outcar_par and valid_oszicar_par:
        temp_str = temp_str + (
            #(folder_path_str_list[i_folder_indx] + ' ' * (max_length_folder_path_str - len(folder_path_str_list[i_folder_indx])) + '   ' + 
            #job_status_str + ' ' * (11 - len(job_status_str)) + 
            #'{:.4f}'.format(outcar_params_dict['energy_without_entropy']) + ' ' * (16 - len('{:.4f}'.format(outcar_params_dict['energy_without_entropy']))) + 
            #'{:.4f}'.format(outcar_params_dict['TOTEN']) + ' ' * (13 - len('{:.4f}'.format(outcar_params_dict['TOTEN']))) + 
            #'{:.4f}'.format(outcar_params_dict['energy(sigma->0)']) + ' ' * (17 - len('{:.4f}'.format(outcar_params_dict['energy(sigma->0)']))) + 
            #'{:.4f}'.format(elapsed_time_hour) + ' ' * (14 - len('{:.4f}'.format(elapsed_time_hour))) 
            ## '{:.4f}'.format(outcar_params_dict['e_fermi']) + ' ' * (8 - len('{:.4f}'.format(outcar_params_dict['e_fermi']))) + 
            funcs.str_format(folder_path_str_list[i_folder_indx], max_len = max_length_folder_path_str) + '   '
            + funcs.str_format(job_status_str, max_len = 13) 
            + funcs.float_format(outcar_params_dict['energy_without_entropy'], max_len = 16, format_str = '{:.4f}')
            + funcs.float_format(outcar_params_dict['TOTEN'], max_len = 13, format_str = '{:.4f}')
            + funcs.float_format(outcar_params_dict['energy(sigma->0)'], max_len = 17, format_str = '{:.4f}')
            + funcs.float_format(elapsed_time_hour, max_len = 14, format_str = '{:.4f}')
            )
        if write_band_gap == True:
            temp_str = temp_str + funcs.str_format(band_gap, max_len = 8)
        temp_str = temp_str + ' \n'
    else:
        temp_str = temp_str + (
            funcs.str_format(folder_path_str_list[i_folder_indx], max_len = max_length_folder_path_str) + '   '
            + funcs.str_format(job_status_str, max_len = 13)
            + funcs.str_format('--', max_len = 16)
            + funcs.str_format('--', max_len = 13)
            + funcs.str_format('--', max_len = 17)
            + funcs.str_format('--', max_len = 14)
            )
        if write_band_gap == True:
            temp_str = temp_str + funcs.str_format('--', max_len = 8)
        temp_str = temp_str + ' \n'
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

def rm_outputs(dst_dir = None, clean_subdir = True):
    '''
    remove the output files of VASP jobs
    clean_subdir = True: also clean the subdirectories
    '''
    import os
    from .. import default_params
    from .. import funcs

    defaults_dict = default_params.default_params()
    vasp_output_list = ['CHGCAR', 'CHG', 'CONTCAR', 'DOSCAR', 'EIGENVAL', 'IBZKPT', 'OSZICAR', 'OUTCAR', 'PCDAT', 'PROCAR', 'WAVECAR', 'XDATCAR', 'vasprun.xml'] 
    aux_file_list = ['e.*', 'o.*', 'error.*', 'output.*', 'vaspjob.*', 'REPORT']
    
    dst_dir = os.path.abspath(dst_dir)
    if clean_subdir == True:
        dir_list =  funcs.get_dirs(dst_dir)
    elif clean_subdir == False:
        dir_list = [dst_dir]
    for i_dir in dir_list:
        for item in vasp_output_list:
            try:
                os.remove(os.path.join(i_dir, item)) 
            except OSError:
                pass
    return 0

def poscar_layer_dist_tolerance(poscar_dict, radius_style = 'csd', radius_scale_factor = 1.15):
    ###############################
    # Set layer_dist_tolerance
    ###############################
    from .. import periodic_table

    periodic_table_dict = periodic_table.periodic_tab()

    if radius_style in ['csd', 'CSD']:
        # Find the smallest two CSD covalent radii of the elements in the system. The layer_dist_tolerance is 115% of the sum of two smallest CSD covalent radii.
        csd_covalent_radius_list = [periodic_table_dict['csd_covalent_radius'][i_elmt] for i_elmt in poscar_dict['elmt_species_arr']]
        csd_covalent_radius_list.sort()
        if len(csd_covalent_radius_list) == 1:
            smallest_atom_radius = csd_covalent_radius_list[0]
            #second_smallest_atom_radius = csd_covalent_radius_list[0]
            largest_atom_radius = csd_covalent_radius_list[-1]
            #second_largest_atom_radius = csd_covalent_radius_list[-1]
        else:
            smallest_atom_radius = csd_covalent_radius_list[0]
            #second_smallest_atom_radius = csd_covalent_radius_list[1]
            largest_atom_radius = csd_covalent_radius_list[-1]
            #second_largest_atom_radius = csd_covalent_radius_list[-2]
        min_atom_bonding_len = smallest_atom_radius * 2 / 100 * radius_scale_factor
        max_atom_bonding_len = largest_atom_radius * 2 / 100 * radius_scale_factor
        layer_dist_tolerance = min_atom_bonding_len
        ##print('automatically generated layer_dist_tolerance=', layer_dist_tolerance)
    return layer_dist_tolerance

def poscar_direction_params(poscar_dict, direction = 'z'):
    '''
    get parameters related with a specified direction in a POSCAR file.
    '''
    import numpy as np
    from .. import funcs
    poscar_direction_dict = {}
    poscar_direction_dict['direction'] = {}
    poscar_direction_dict['number_order_text'] = None
    poscar_direction_dict['l_arr_row'] = None
    poscar_direction_dict['l_arr_column_1'] = None
    poscar_direction_dict['l_arr_column_2'] = None
    poscar_direction_dict['side_vector'] = None
    poscar_direction_dict['side_vector_len'] = None
    poscar_direction_dict['vec_1'] = None
    poscar_direction_dict['vec_2'] = None
    poscar_direction_dict['unit_vector'] = None
    poscar_direction_dict['ortho_vector'] = None
    poscar_direction_dict['ortho_vector_1'] = None 
    poscar_direction_dict['ortho_vector_2'] = None
    poscar_direction_dict['box_len_ortho'] = None
    poscar_direction_dict['cos_angle'] = None
    poscar_direction_dict['pos_arr_column'] = None
    poscar_direction_dict['pos_arr_column_1'] = None
    poscar_direction_dict['pos_arr_column_2'] = None
    poscar_direction_dict['pos_arr_direct_column'] = None
    poscar_direction_dict['pos_arr_direct_column_1'] = None
    poscar_direction_dict['pos_arr_direct_column_2'] = None

    poscar_direction_dict['direction'] = direction
    if direction in ['x', 'X']:
        #number_order_text = 'yzx'
        number_order_text = 'x'
        l_arr_row = 0
        l_arr_column_1 = 1
        l_arr_column_2 = 2
        side_vector = poscar_dict['vec_a']
        side_vector_len = poscar_dict['len_vec_a']
        vec_1 = poscar_dict['vec_b']
        vec_2 = poscar_dict['vec_c']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([1, 0, 0])
        ortho_vector_1 = np.array([0, 1, 0])
        ortho_vector_2 = np.array([0, 0, 1])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 3
        pos_arr_column_1 = 4
        pos_arr_column_2 = 5
        pos_arr_direct_column = 0
        pos_arr_direct_column_1 = 1
        pos_arr_direct_column_2 = 2
    if direction in ['y', 'Y']:
        #number_order_text = 'zxy'
        number_order_text = 'y'
        l_arr_row = 1
        l_arr_column_1 = 0
        l_arr_column_2 = 2
        side_vector = poscar_dict['vec_b']
        side_vector_len = poscar_dict['len_vec_b']
        vec_1 = poscar_dict['vec_a']
        vec_2 = poscar_dict['vec_c']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 1, 0])
        ortho_vector_1 = np.array([1, 0, 0])
        ortho_vector_2 = np.array([0, 0, 1])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 4
        pos_arr_column_1 = 3
        pos_arr_column_2 = 5
        pos_arr_direct_column = 1
        pos_arr_direct_column_1 = 0
        pos_arr_direct_column_2 = 2
    if direction in ['z', 'Z']:
        #number_order_text = 'xyz'
        number_order_text = 'z'
        l_arr_row = 2
        l_arr_column_1 = 0
        l_arr_column_2 = 1
        side_vector = poscar_dict['vec_c']
        side_vector_len = poscar_dict['len_vec_c']
        vec_1 = poscar_dict['vec_a']
        vec_2 = poscar_dict['vec_b']
        unit_vector = funcs.unit_vec(side_vector)
        ortho_vector = np.array([0, 0, 1])
        ortho_vector_1 = np.array([1, 0, 0])
        ortho_vector_2 = np.array([0, 1, 0])
        box_len_ortho = (poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
        cos_angle = unit_vector.dot(ortho_vector)
        pos_arr_column = 5
        pos_arr_column_1 = 3
        pos_arr_column_2 = 4
        pos_arr_direct_column = 2
        pos_arr_direct_column_1 = 0
        pos_arr_direct_column_2 = 1
    poscar_direction_dict['number_order_text'] = number_order_text
    poscar_direction_dict['l_arr_row'] = l_arr_row
    poscar_direction_dict['l_arr_column_1'] = l_arr_column_1
    poscar_direction_dict['l_arr_column_2'] = l_arr_column_2
    poscar_direction_dict['side_vector'] = side_vector
    poscar_direction_dict['side_vector_len'] = side_vector_len
    poscar_direction_dict['vec_1'] = vec_1
    poscar_direction_dict['vec_2'] = vec_2
    poscar_direction_dict['unit_vector'] = unit_vector
    poscar_direction_dict['ortho_vector'] = ortho_vector
    poscar_direction_dict['ortho_vector_1'] = ortho_vector_1
    poscar_direction_dict['ortho_vector_2'] = ortho_vector_2
    poscar_direction_dict['box_len_ortho'] = box_len_ortho
    poscar_direction_dict['cos_angle'] = cos_angle
    poscar_direction_dict['pos_arr_column'] = pos_arr_column
    poscar_direction_dict['pos_arr_column_1'] = pos_arr_column_1
    poscar_direction_dict['pos_arr_column_2'] = pos_arr_column_2
    poscar_direction_dict['pos_arr_direct_column'] = pos_arr_direct_column
    poscar_direction_dict['pos_arr_direct_column_1'] = pos_arr_direct_column_1
    poscar_direction_dict['pos_arr_direct_column_2'] = pos_arr_direct_column_2
    return poscar_direction_dict

def poscar_layer_params(poscar_dict, poscar_direction_dict, criteria = 'auto', delta = 0.05, layer_dist_tolerance = 'auto', radius_style = 'csd', radius_scale_factor = 1.15, write_layer_info = False, suppress_warning = True):
    '''
    Get layer parameters from POSCAR
    delta: This is the spacial resolution (in unit of Angstrom) for the atoms, for example z=0.24 and z=0.25 are considered to be in the same position if delta=0.01. Default: delta = 0.05
    '''
    import numpy as np
    from .. import funcs
    from . import vasp_read
    from .. import periodic_table
    from copy import copy, deepcopy
    import math
    import os
    from . import vasp_build

    periodic_table_dict = periodic_table.periodic_tab()

    original_criteria = criteria
    original_layer_dist_tolerance = layer_dist_tolerance
    #recommended_layer_dist_tolerance = 0.664 
    #################################################
    # Creating layer property dictionary: poscar_layer_dict
    #################################################
    poscar_layer_dict = {}
    poscar_layer_dict['delta'] = delta

    poscar_dict['atom_number_ortho_' + poscar_direction_dict['number_order_text']] = funcs.atom_number_ortho(
        atom_key_arr = poscar_dict['atom_key_arr'],
        pos_arr_cartesian= poscar_dict['pos_arr'][:,3:6],
        delta = delta,
        number_order = poscar_direction_dict['number_order_text']
        )
    layer_dist_tolerance = poscar_layer_dist_tolerance(poscar_dict, radius_style = radius_style, radius_scale_factor = radius_scale_factor)
    # group layers according to different properties and adding tags
    poscar_layer_dict['num_layers'] = 0 
    poscar_layer_dict['layer_coord_list'] = [None]
    poscar_layer_dict['indx_list'] = [None]
    poscar_layer_dict['atoms_indx_list'] =  [list()]
    #bond_status_list is a list for identifying the bonding status for each atomic layer. The bonding status are listed below:
    #bond_status = 0: no bond with the adjacent two layers
    #bond_status = -1: forms bond only with one side (left or bottom) of the two layers
    #bond_status = 1: forms bond only with one side (right or top) of the two layers
    #bond_status = 2: forms bond with both the adjacent two layers
    poscar_layer_dict['bond_status_list'] = [None]
    poscar_layer_dict['layer_bonding_indx_list'] = [None]
    poscar_layer_dict['layer_bonding_atom_name_list'] = [None]
    # periodicity_status_list is a list of indices that dinstinguishes each layer, the layer indices are labeled in a layer-by-layer fashion, the index for the first layer is 0. If the next layer is different from the previous layer, add 1 to the previous layer index. If there is another layer that is the same as the first layer then that layer is indexed to 0, and this means a new periodicity of lattice starts. The example of periodicity_status_list is :
    # [0,0,0,0,0] is a model with five equivalent layers
    # [0,1,2,3,0,1,2,3] is a model with 8 layers with a periodicity of 2, each periodicity conatains 4 layers labeled with 0,1,2,3
    # [0,1,2,1,2,1,2] is a model with 7 layers with only 1 periodicity, the only periodicity contains 7 layers labeled with 0,1,2,1,2,1,2
    poscar_layer_dict['periodicity_status_list'] = [None]
    poscar_layer_dict['vacuum_exist'] = None
    poscar_layer_dict['vacuum_width'] = None
    poscar_layer_dict['crystal_width'] = None

    ######################################################################################################################################
    # Check whether the value of layer_dist_tolerance is appropriate or not.
    # First, reorder the atoms in certain direction, then get the difference of the projected distances of the atom and the last atom
    # Next, find the largest distance between atoms in certain direction and denote it as max_ordered_atom_dist
    # if layer_dist_tolerance is larger than max_ordered_atom_dist, then set layer_dist_tolerance to layer_dist_tolerance = max_ordered_atom_dist - delta_dist, in which delta_dist is a small value. e.g delta_dist = 0.01
    ######################################################################################################################################
    delta_layer_dist = 0.01
    atom_pos_list = [None] * poscar_dict['n_atoms']
    for i in range(1, poscar_dict['n_atoms'] + 1):
        i_atom_indx = int(np.argwhere(poscar_dict['atom_number_ortho_' + poscar_direction_dict['number_order_text']] == i))
        i_atom_pos = poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_column']]
        atom_pos_list[i - 1] = i_atom_pos
    last_atom_pos_list = [0] + atom_pos_list
    del last_atom_pos_list[-1]
    ordered_atom_pos_arr = np.array(atom_pos_list)
    ordered_last_atom_pos_arr = np.array(last_atom_pos_list)
    ordered_atom_dist_arr = ordered_atom_pos_arr - ordered_last_atom_pos_arr
    abs_ordered_atom_dist_arr = np.abs(ordered_atom_dist_arr, out = ordered_atom_dist_arr)
    max_ordered_atom_dist = np.max(abs_ordered_atom_dist_arr)	
    ###############################################################
    #find the reasonable interlayer distance 
    ###############################################################
    if radius_style in ['csd', 'CSD']:
        csd_covalent_radius_list = [periodic_table_dict['csd_covalent_radius'][i_elmt] for i_elmt in poscar_dict['elmt_species_arr']]
        csd_covalent_radius_list.sort()
        smallest_atom_radius = csd_covalent_radius_list[0]
        largest_atom_radius = csd_covalent_radius_list[-1]
        min_atom_bonding_len = smallest_atom_radius * 2 / 100 * radius_scale_factor
        max_atom_bonding_len = largest_atom_radius * 2 / 100 * radius_scale_factor
    low_interlayer_dist = poscar_layer_dict['delta']
    high_interlayer_dist = max_atom_bonding_len 
    normal_interlayer_dist_arr = np.array([x for x in abs_ordered_atom_dist_arr if x > low_interlayer_dist and x < high_interlayer_dist])
    #normal_interlayer_dist_arr = np.array([x for x in normal_interlayer_dist_arr if abs(x - low_interlayer_dist) > delta])    
    #print('normal_interalyer_dist_arr=',normal_interlayer_dist_arr)
    if len(normal_interlayer_dist_arr) != 0:
        poscar_layer_dict['corrected_layer_dist_tolerance'] = np.min(normal_interlayer_dist_arr) - delta_layer_dist
    else:
        poscar_layer_dict['corrected_layer_dist_tolerance'] = max_ordered_atom_dist - delta_layer_dist
    #print('corrected_layer_dist_tolerance=',poscar_layer_dict['corrected_layer_dist_tolerance'])

    #layer_dist_tolerance = poscar_layer_dict['corrected_layer_dist_tolerance']

    # Set criteria
    if layer_dist_tolerance >= max_ordered_atom_dist:
        ##print('WARNING #2103202201 (from vasp_build.exfoliation): layer_dist_tolerance = ' + str(layer_dist_tolerance) + ' is too large for model ' + str(poscar_dict['file_path']) + '. It is automatically set to layer_dist_tolerance = ' + str(poscar_layer_dict['corrected_layer_dist_tolerance']))
        layer_dist_tolerance = poscar_layer_dict['corrected_layer_dist_tolerance']
        if criteria == 'auto':
            criteria = 'periodicity'
        elif criteria == 'bonding':
            if suppress_warning == False:
                print('# WARNING #2103211849 (from vasp_tools): Layers are connected by atom bonds. It is recommended to use the \'periodicity\' criterion to exfoliate the 2D material.')
    elif layer_dist_tolerance < max_ordered_atom_dist:
        layer_dist_tolerance = poscar_layer_dict['corrected_layer_dist_tolerance']
        if criteria == 'auto':
            criteria = 'bonding'
        elif criteria == 'periodicity':
            if suppress_warning == False:
                print('There are atomic layers without interlayer atom bonding. It is recommended to use the \'bonding\' criterion to exfoliate the 2D material.')
    poscar_layer_dict['criteria'] = criteria
    #print('******',original_layer_dist_tolerance, layer_dist_tolerance)
    # Use user-defined or automatically generated layer_dist_tolerance
    if original_layer_dist_tolerance == 'auto':
        pass
    else:
        # use user-defined value 
        layer_dist_tolerance = original_layer_dist_tolerance
    poscar_layer_dict['layer_dist_tolerance'] = layer_dist_tolerance
    ##print('layer_dist_tolerance=================',layer_dist_tolerance)

    ##if layer_dist_tolerance < recommended_layer_dist_tolerance:
    ##    print('WARNING #2012301056 (from vasp_build.exfoliation): layer_dist_tolerance is smaller than ' + str(recommended_layer_dist_tolerance) + '(Default value of layer_dist_tolerance is (1.15 * 2 * (CSD_covalent_radius of He)) + arbitrary_delta = 1.15 * 2 * 0.28 + 0.02 = 0.664 Angstrom), It is recommended to be larger than this value.')
    ##def poscar_layer_status(poscar_dict, poscar_direction_dict, poscar_layer_dict, radius_style = 'csd', radius_scale_factor = 1.15, write_layer_info = False, suppress_warning = True):
    def poscar_layer_status(poscar_dict, poscar_direction_dict, delta = 0.05, radius_style = 'csd', radius_scale_factor = 1.15, write_layer_info = False, suppress_warning = True):
        '''
        Get the status of each layer
        Two modes: bonding, periodicity
        '''
        import numpy as np
        import math
        ###############################################################################
        # Reinitialize the parameters in poscar_layer_dict
        ###############################################################################
        poscar_layer_dict = {}
        poscar_layer_dict['delta'] = delta
        poscar_layer_dict['num_layers'] = 0 
        poscar_layer_dict['layer_coord_list'] = [None]
        poscar_layer_dict['indx_list'] = [None]
        poscar_layer_dict['atoms_indx_list'] =  [list()]
        poscar_layer_dict['atom_name_list'] =  [list()]
        poscar_layer_dict['bond_status_list'] = [None]
        poscar_layer_dict['layer_bonding_indx_list'] = [None]
        poscar_layer_dict['layer_bonding_atom_name_list'] = [None]
        poscar_layer_dict['periodicity_status_list'] = [None]
        poscar_layer_dict['vacuum_exist'] = None
        poscar_layer_dict['vacuum_width'] = None
        poscar_layer_dict['crystal_width'] = None
        ###############################################################################
        # Check bonding and periodicity status of atomic layers
        ###############################################################################
        # temp_pos can be used as a representitive position of each layer
        # indx0 and temp_pos_0 denote the index and position of the atom in the bottom layer
        indx_0 = np.argwhere(poscar_dict['atom_number_ortho_' + poscar_direction_dict['number_order_text']] == 1)
        temp_pos_0 = poscar_dict['pos_arr'][indx_0, poscar_direction_dict['pos_arr_column']]
        temp_pos = temp_pos_0
        i_layer_indx = 0
        bond_breaking = True
        shift_vec = poscar_dict['l_arr'][poscar_direction_dict['l_arr_row'],:]
        # First, find the atoms in each layer according to the parameter 'delta'
        for i in range(1, poscar_dict['n_atoms'] + 1):
            shift_top2bottom = False
            i_atom_indx = int(np.argwhere(poscar_dict['atom_number_ortho_' + poscar_direction_dict['number_order_text']] == i))
            i_atom_pos = poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_column']]
            i_atom_pos_shift = (poscar_dict['pos_arr'][i_atom_indx, 3:6] - shift_vec)[poscar_direction_dict['pos_arr_column'] -3]
            if abs(i_atom_pos_shift - temp_pos_0) <= poscar_layer_dict['delta']:
                # for the atoms near the topmost boundary, check whether it is in the same layer as the bottom layer
                i_layer_indx = 0
                poscar_layer_dict['atoms_indx_list'][i_layer_indx].append(i_atom_indx)
                poscar_layer_dict['atom_name_list'][i_layer_indx].append(poscar_dict['atomname_list'][i_atom_indx])
                shift_top2bottom = True
            elif abs(i_atom_pos - temp_pos) <= poscar_layer_dict['delta']:
                # The i-th atom is in the same layer as the representitive atom
                pass
            else:
                i_layer_indx = i_layer_indx + 1
                poscar_layer_dict['indx_list'].append(None)
                poscar_layer_dict['layer_coord_list'].append(None)
                poscar_layer_dict['atoms_indx_list'].append(list())
                poscar_layer_dict['atom_name_list'].append(list())
                poscar_layer_dict['bond_status_list'].append(None)
                poscar_layer_dict['layer_bonding_indx_list'].append(None)
                poscar_layer_dict['layer_bonding_atom_name_list'].append(None)
                poscar_layer_dict['periodicity_status_list'].append(None)
            poscar_layer_dict['indx_list'][i_layer_indx] = i_layer_indx
            poscar_layer_dict['layer_coord_list'][i_layer_indx] = i_atom_pos
            if shift_top2bottom == False:
                poscar_layer_dict['atoms_indx_list'][i_layer_indx].append(i_atom_indx)
                poscar_layer_dict['atom_name_list'][i_layer_indx].append(poscar_dict['atomname_list'][i_atom_indx])
            poscar_layer_dict['bond_status_list'][i_layer_indx] = 0
            poscar_layer_dict['layer_bonding_indx_list'][i_layer_indx] = [] 
            poscar_layer_dict['layer_bonding_atom_name_list'][i_layer_indx] = [] 
            poscar_layer_dict['periodicity_status_list'][i_layer_indx] = None
            temp_pos = i_atom_pos
        poscar_layer_dict['num_layers'] = len(poscar_layer_dict['indx_list'])
         
        #######################################################################################
        # check bonding status of atoms among the i-th layer with the j-th and the k-th layer
        # where j = i + n, k = i - n, n = 1, 2, 3, 4, 5 and n should be smaller than half the 
        # number of layers
        #######################################################################################
        bottom_indx = 0
        top_indx = poscar_layer_dict['num_layers'] - 1
        n_indx_max = 15
        if n_indx_max > math.floor(poscar_layer_dict['num_layers'] / 2):
            n_indx_max = math.floor(poscar_layer_dict['num_layers'] / 2)
        if n_indx_max < 3:
            n_indx_max = 3
        for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
            ij_atom_bonding_list = []
            ik_atom_bonding_list = []
            for n_indx in range(1, n_indx_max):
                # j and k layers are the adjacent layers for layer i
                ij_layer_bonding = False
                ik_layer_bonding = False
                j_layer_indx_temp = i_layer_indx - n_indx
                k_layer_indx_temp = i_layer_indx + n_indx
                base_j = math.floor(j_layer_indx_temp / poscar_layer_dict['num_layers'])
                base_k = math.floor(k_layer_indx_temp / poscar_layer_dict['num_layers'])
                j_layer_indx = j_layer_indx_temp - base_j * poscar_layer_dict['num_layers']
                k_layer_indx = k_layer_indx_temp - base_k * poscar_layer_dict['num_layers']

                i_layer_atom_indx_list = poscar_layer_dict['atoms_indx_list'][i_layer_indx]
                j_layer_atom_indx_list = poscar_layer_dict['atoms_indx_list'][j_layer_indx]
                k_layer_atom_indx_list = poscar_layer_dict['atoms_indx_list'][k_layer_indx]
                for i_atom_indx in i_layer_atom_indx_list:
                    i_elmt = poscar_dict['atom_species_arr'][i_atom_indx]
                    i_atom_pos = poscar_dict['pos_arr'][i_atom_indx, 3:6]
                    for j_atom_indx in j_layer_atom_indx_list:
                        j_elmt = poscar_dict['atom_species_arr'][j_atom_indx]
                        j_atom_pos = poscar_dict['pos_arr'][j_atom_indx, 3:6] + base_j * shift_vec
                        ##print(n_indx, base_j, poscar_dict['pos_arr'][j_atom_indx, 3:6])
                        ##print(j_atom_pos)
                        bonding = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos, mode = radius_style, factor = radius_scale_factor) 
                        ij_atom_bonding_list.append(bonding)
                        if bonding == True:
                            ij_layer_bonding = True
                            poscar_layer_dict['layer_bonding_atom_name_list'][i_layer_indx].append(poscar_dict['atomname_list'][i_atom_indx] + '-' + poscar_dict['atomname_list'][j_atom_indx])  
                    for k_atom_indx in k_layer_atom_indx_list:
                        k_elmt = poscar_dict['atom_species_arr'][k_atom_indx]
                        k_atom_pos = poscar_dict['pos_arr'][k_atom_indx, 3:6] + base_k * shift_vec
                        bonding = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos, mode = radius_style, factor = radius_scale_factor) 
                        ik_atom_bonding_list.append(bonding)
                        if bonding == True:
                            ik_layer_bonding = True
                            poscar_layer_dict['layer_bonding_atom_name_list'][i_layer_indx].append(poscar_dict['atomname_list'][i_atom_indx] + '-' + poscar_dict['atomname_list'][k_atom_indx])  
                # Determine the layer_bonding_indx_list, this list records the layer indices that form bonds with the i-th layer
                if ij_layer_bonding == True: 
                    poscar_layer_dict['layer_bonding_indx_list'][i_layer_indx].append(j_layer_indx)  
                if ik_layer_bonding == True: 
                    poscar_layer_dict['layer_bonding_indx_list'][i_layer_indx].append(k_layer_indx)

                # Determine the bond_status_list
                if True in ij_atom_bonding_list and True in ik_atom_bonding_list:
                    # The i-th layer form bonds with the adjacent j and k layers
                    poscar_layer_dict['bond_status_list'][i_layer_indx] = 2
                elif True in ij_atom_bonding_list and True not in ik_atom_bonding_list:
                    # The i-th layer does not form bonds with the adjacent k layer, but forms bond with the adjacent j layer
                    poscar_layer_dict['bond_status_list'][i_layer_indx] = 1
                elif True not in ij_atom_bonding_list and True in ik_atom_bonding_list:
                    # The i-th layer does not form bonds with the adjacent j layer, but forms bond with the adjacent k layer
                    poscar_layer_dict['bond_status_list'][i_layer_indx] = -1
                elif True not in ij_atom_bonding_list and True not in ik_atom_bonding_list:
                    # The i-th layer does not form bonds with the adjacent two layers
                    poscar_layer_dict['bond_status_list'][i_layer_indx] = 0
            ##print('i=', i_layer_indx, poscar_layer_dict['layer_bonding_indx_list'][i_layer_indx])
        ########################################################################################################
        # Check the validity of the bond_status list
        #
        # The following bond_status_list is not correct:
        # bond_status_list= [-1, 1, 2, -1, 1, -1, 1, 2, -1, 1, -1, 1, 2, -1, 1, -1, 1, 2, -1, 1] e.g. mp-9789
        # bond_status_list= [-1, -1, -1, 2, 1, -1, 2, 1, 1, 1, -1, -1, -1, 2, 1, -1, 2, 1, 1, 1] e.g. mp-27178
        #
        # We should be careful with the following situations, in case of i-th layer does not form bonds with 
        # the j-th layer, but forms bond with the (j-1) layer,
        # Situation 1:
        # "1" should not be on the left side of "2", and "-1" should not be on the right side of "2"
        # The left side of "2" should only be "-1" or "2"
        # The right side of "2" should only be "1" or "2"
        # if "1,2" is encountered, it will be changed to "2,2"
        # else if "2, -1" is encountered, it will be changed to "2,2"
        # Situation 2:
        # "1" should not be on the left side of "1", and "-1" should not be on the right side of "-1"
        # The left side of "1" should only be "-1" or "2"
        # The right side of "-1" should only be "1" or "2"
        # if "1,1" is encountered, it will be changed to "2,1"
        # else if "-1, -1" is encountered, it will be changed to "-1,2"
        #
        # in the above description, "left" deontes "bottom" and "right" denotes "top"
        ########################################################################################################
        if 0 not in poscar_layer_dict['bond_status_list']:
            for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
                for n_indx in range(2, n_indx_max):
                    j_layer_indx_temp = i_layer_indx - n_indx
                    k_layer_indx_temp = i_layer_indx + n_indx
                    base_j = math.floor(j_layer_indx_temp / poscar_layer_dict['num_layers'])
                    base_k = math.floor(k_layer_indx_temp / poscar_layer_dict['num_layers'])
                    j_layer_indx = j_layer_indx_temp - base_j * poscar_layer_dict['num_layers']
                    k_layer_indx = k_layer_indx_temp - base_k * poscar_layer_dict['num_layers']

                    if j_layer_indx in poscar_layer_dict['layer_bonding_indx_list'][i_layer_indx]: 
                        for temp_n_indx in range(1, n_indx):
                            j_layer_indx_temp1 = i_layer_indx - temp_n_indx
                            base_j1 = math.floor(j_layer_indx_temp1 / poscar_layer_dict['num_layers'])
                            j_layer_indx1 = j_layer_indx_temp1 - base_j1 * poscar_layer_dict['num_layers']
                            poscar_layer_dict['bond_status_list'][j_layer_indx1] = 2
                    if k_layer_indx in poscar_layer_dict['layer_bonding_indx_list'][i_layer_indx]: 
                        for temp_n_indx in range(1, n_indx):
                            k_layer_indx_temp1 = i_layer_indx + temp_n_indx
                            base_k1 = math.floor(k_layer_indx_temp1 / poscar_layer_dict['num_layers'])
                            k_layer_indx1 = k_layer_indx_temp1 - base_k1 * poscar_layer_dict['num_layers']
                            poscar_layer_dict['bond_status_list'][k_layer_indx1] = 2
            ########################################################################################
            # Treating model with  all layers bonded, but with a vacuum layer
            # Bond status cannot be solely used to treat this situation.
            # if the layer separation between the i-th and the k-th layer is larger than
            # the thresold_vacuum_width, the vacuum exists. And the i-th item and the k-th items of 
            # the bond_status_list take the values of 1 and -1
            ########################################################################################
            thresold_vacuum_width = 5.0
            if len(set(poscar_layer_dict['bond_status_list'])) == 1 and len(poscar_layer_dict['bond_status_list']) != 1 and 2 in poscar_layer_dict['bond_status_list']:
                top_layer_indx = poscar_layer_dict['num_layers'] - 1
                bottom_layer_indx = 0
                top_width = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['layer_coord_list'][top_layer_indx]
                bottom_width = poscar_layer_dict['layer_coord_list'][bottom_layer_indx]
                top_plus_bottom_width = top_width + bottom_width
                if top_plus_bottom_width > thresold_vacuum_width:
                    poscar_layer_dict['bond_status_list'][top_layer_indx] = 1
                    poscar_layer_dict['bond_status_list'][bottom_layer_indx] = -1
                # Get layer spearation between the i-th and the k-th (k-th layer is the next layer to i-th layer)
                ik_layer_separation_list = []
                n_indx = 1
                for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
                    k_layer_indx_temp = i_layer_indx + n_indx
                    base_k = math.floor(k_layer_indx_temp / poscar_layer_dict['num_layers'])
                    k_layer_indx = k_layer_indx_temp - base_k * poscar_layer_dict['num_layers']
                    ik_layer_separation = abs(poscar_layer_dict['layer_coord_list'][k_layer_indx] - poscar_layer_dict['layer_coord_list'][i_layer_indx])
                    ik_layer_separation_list.append(ik_layer_separation)

                    if ik_layer_separation > thresold_vacuum_width and i_layer_indx != top_layer_indx:
                        poscar_layer_dict['bond_status_list'][i_layer_indx] = 1
                        poscar_layer_dict['bond_status_list'][k_layer_indx] = -1
                        break
                    if ik_layer_separation > thresold_vacuum_width and i_layer_indx == top_layer_indx:
                        if len(poscar_layer_dict['bond_status_list']) == 2:
                            poscar_layer_dict['bond_status_list'][i_layer_indx] = 1
                            poscar_layer_dict['bond_status_list'][k_layer_indx] = -1
                            break
                        elif len(poscar_layer_dict['bond_status_list']) > 2:
                            # get the average layer separation of all the layers 
                            mean_layer_separation = np.mean(ik_layer_separation_list)
                            if top_plus_bottom_width > mean_layer_separation * 2:
                                poscar_layer_dict['bond_status_list'][i_layer_indx] = 1
                                poscar_layer_dict['bond_status_list'][k_layer_indx] = -1

        ######################################################################
        # Get vacuum width and crystal_width
        ######################################################################
        n_indx = 1
        poscar_direction_dict['pos_arr_column']
        i_layer_atom_pos_list = []
        k_layer_atom_pos_list = []
        poscar_layer_dict['vacuum_exist'] = False
        if len(poscar_layer_dict['bond_status_list']) == 1 and 0 in poscar_layer_dict['bond_status_list']:
            # Dealing with single layer model
            poscar_layer_dict['vacuum_exist'] = True
            poscar_layer_dict['crystal_width'] = 0
            poscar_layer_dict['vacuum_width'] = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['crystal_width']
        elif len(set(poscar_layer_dict['bond_status_list'])) == 1 and 2 in poscar_layer_dict['bond_status_list']:
            poscar_layer_dict['vacuum_exist'] = False
            poscar_layer_dict['crystal_width'] = poscar_direction_dict['box_len_ortho']
            poscar_layer_dict['vacuum_width'] = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['crystal_width']
        else:
            # Deal with multilayer model
            for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
                k_layer_indx_temp = i_layer_indx + n_indx
                base_k = math.floor(k_layer_indx_temp / poscar_layer_dict['num_layers'])
                k_layer_indx = k_layer_indx_temp - base_k * poscar_layer_dict['num_layers']
                if poscar_layer_dict['bond_status_list'][i_layer_indx] == 1 and poscar_layer_dict['bond_status_list'][k_layer_indx] == -1:
                    poscar_layer_dict['vacuum_exist'] = True
                    for i_atom_indx in poscar_layer_dict['atoms_indx_list'][i_layer_indx]:
                        i_layer_atom_pos_list.append(poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_column']])
                    for k_atom_indx in poscar_layer_dict['atoms_indx_list'][k_layer_indx]:
                        k_layer_atom_pos_list.append(poscar_dict['pos_arr'][k_atom_indx, poscar_direction_dict['pos_arr_column']])
                    if np.min(i_layer_atom_pos_list) <= np.min(k_layer_atom_pos_list):
                        poscar_layer_dict['vacuum_width'] = abs(np.min(k_layer_atom_pos_list) - np.max(i_layer_atom_pos_list))
                        poscar_layer_dict['crystal_width'] = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['vacuum_width']
                    elif np.min(i_layer_atom_pos_list) > np.min(k_layer_atom_pos_list):
                        poscar_layer_dict['crystal_width'] = abs(np.max(i_layer_atom_pos_list) - np.min(k_layer_atom_pos_list))
                        poscar_layer_dict['vacuum_width'] = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['crystal_width']
                elif len(poscar_layer_dict['bond_status_list']) > 1 and poscar_layer_dict['bond_status_list'][i_layer_indx] == 0:
                    poscar_layer_dict['vacuum_exist'] = True
                    for i_atom_indx in poscar_layer_dict['atoms_indx_list'][i_layer_indx]:
                        i_layer_atom_pos_list.append(poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_column']])
                    for k_atom_indx in poscar_layer_dict['atoms_indx_list'][k_layer_indx]:
                        k_layer_atom_pos_list.append(poscar_dict['pos_arr'][k_atom_indx, poscar_direction_dict['pos_arr_column']])
                    if i_layer_indx != (poscar_layer_dict['num_layers'] - 1): 
                        poscar_layer_dict['vacuum_width'] = abs(np.min(k_layer_atom_pos_list) - np.max(i_layer_atom_pos_list))
                        poscar_layer_dict['crystal_width'] = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['vacuum_width']
                    elif i_layer_indx == (poscar_layer_dict['num_layers'] - 1): 
                        poscar_layer_dict['crystal_width'] = abs(np.max(i_layer_atom_pos_list) - np.min(k_layer_atom_pos_list))
                        poscar_layer_dict['vacuum_width'] = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['crystal_width']

        #########################################
        # check periodicity status of layers
        #########################################
        #poscar_direction_dict = poscar_direction_params(poscar_dict = poscar_dict, direction = direction)
        atom_dist_tolerance = 0.05 
        periodicity_number = 0
        poscar_layer_dict['periodicity_status_list'][0] = periodicity_number 
        ##for i_atom_indx in range(poscar_dict['n_atoms']):
        ##    print(poscar_dict['atomname_list'][i_atom_indx], poscar_dict['pos_arr'][i_atom_indx,3:6])
        for i_layer_indx in range(0, poscar_layer_dict['num_layers'] - 1):
            if poscar_layer_dict['periodicity_status_list'][i_layer_indx] == None:
                periodicity_number = periodicity_number + 1
                poscar_layer_dict['periodicity_status_list'][i_layer_indx] = periodicity_number 
            original_periodicity_number = periodicity_number
            i_layer_atom_indx_list = poscar_layer_dict['atoms_indx_list'][i_layer_indx] 
            i_layer_coord = poscar_layer_dict['layer_coord_list'][i_layer_indx]
            for j_layer_indx in range(i_layer_indx + 1, poscar_layer_dict['num_layers']):
                if poscar_layer_dict['periodicity_status_list'][j_layer_indx] != None:
                    continue
                j_layer_atom_indx_list = poscar_layer_dict['atoms_indx_list'][j_layer_indx] 
                j_layer_coord = poscar_layer_dict['layer_coord_list'][j_layer_indx]

                # find atoms with equivalent position in certain direction
                # Using Direct instead of Cartesian coordinate to realize this function is a good choice
                identical_atoms_list1 = []
                for i_atom_indx in i_layer_atom_indx_list:
                    identical_atoms_list2 = []
                    i_elmt = poscar_dict['atom_species_arr'][i_atom_indx]
                    i_atom_pos_direct = poscar_dict['pos_arr'][i_atom_indx, 0:3]
                    for j_atom_indx in j_layer_atom_indx_list:
                        j_elmt = poscar_dict['atom_species_arr'][j_atom_indx]
                        j_atom_pos_direct = poscar_dict['pos_arr'][j_atom_indx, 0:3]

                        ortho_layer_direct_dist = poscar_dict['pos_arr'][j_atom_indx, poscar_direction_dict['pos_arr_direct_column']] - poscar_dict['pos_arr'][i_atom_indx, poscar_direction_dict['pos_arr_direct_column']]
                        disp_direct_vec = ortho_layer_direct_dist * poscar_direction_dict['ortho_vector']
                        ij_atom_direct_vec = j_atom_pos_direct - (i_atom_pos_direct + disp_direct_vec)
                        ij_atom_cartesian_vec = np.array([None] * 3)
                        ij_atom_cartesian_vec[0] = np.dot(ij_atom_direct_vec, poscar_dict['l_arr'][:,[0]])
                        ij_atom_cartesian_vec[1] = np.dot(ij_atom_direct_vec, poscar_dict['l_arr'][:,[1]])
                        ij_atom_cartesian_vec[2] = np.dot(ij_atom_direct_vec, poscar_dict['l_arr'][:,[2]])
                        ij_atom_cartesian_dist = np.linalg.norm(ij_atom_cartesian_vec) 

                        if i_elmt == j_elmt and ij_atom_cartesian_dist <= atom_dist_tolerance:
                            identical_atoms_list2.append(True)
                        elif i_elmt != j_elmt and ij_atom_cartesian_dist > atom_dist_tolerance:
                            identical_atoms_list2.append(False)
                    if True in identical_atoms_list2:
                        identical_atoms_list1.append(True)
                    else:
                        identical_atoms_list1.append(False)
                if False in identical_atoms_list1:                 
                    pass
                else:
                    poscar_layer_dict['periodicity_status_list'][j_layer_indx] = original_periodicity_number 
        ##################################
        # print atoms in each layer
        ##################################
        if write_layer_info == True:
            layer_info_str = ''
            layer_info_str = layer_info_str + '# ' + poscar_dict['file_path'] + '\n'
            if poscar_direction_dict['number_order_text'] == 'z':
                layer_info_str = layer_info_str + 'direction: c' + '\n'
            elif poscar_direction_dict['number_order_text'] == 'y':
                layer_info_str = layer_info_str + 'direction: b' + '\n'
            elif poscar_direction_dict['number_order_text'] == 'x':
                layer_info_str = layer_info_str + 'direction: a' + '\n'
            layer_info_str = layer_info_str + ( 
                    funcs.str_format('layer', 7, ' ') + 
                    funcs.str_format('bonding', 9, ' ') + 
                    funcs.str_format('periodicity', 13, ' ') + 
                    funcs.str_format('atoms in layer', 25, ' ') + 
                    'bonded atoms list' + '\n')
            for i_layer_indx in reversed(range(poscar_layer_dict['num_layers'])): 
                i_layer_indx_str = str(i_layer_indx + 1)
                i_bond_status_str = str(poscar_layer_dict['bond_status_list'][i_layer_indx])
                i_periodicity_status_str = str(poscar_layer_dict['periodicity_status_list'][i_layer_indx])
                i_atom_name_list_str = str(poscar_layer_dict['atom_name_list'][i_layer_indx])
                i_bonded_atom_name_list_str = str(poscar_layer_dict['layer_bonding_atom_name_list'][i_layer_indx])
                layer_info_str = layer_info_str + ( 
                    funcs.str_format(i_layer_indx_str, 7, ' ') + 
                    funcs.str_format(i_bond_status_str, 9, ' ') + 
                    funcs.str_format(i_periodicity_status_str, 13, ' ') + 
                    funcs.str_format(i_atom_name_list_str, 25, ' ') +
                    i_bonded_atom_name_list_str + '\n')
            print(layer_info_str)
        return poscar_layer_dict

    temp_poscar_layer_dict = poscar_layer_dict
    poscar_layer_dict = poscar_layer_status(poscar_dict = poscar_dict, poscar_direction_dict = poscar_direction_dict, delta = delta, radius_style = radius_style, radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning = suppress_warning)

    ##########################################################################
    # Check the bond_status_list: check for '0' in the bond_status_list
    # Enlarge the atom bonding length by tuning the radius_scale_factor from 
    # a small value to a larger value (add 0.05 to radius_scale_factor each 
    # time up to 0.50) till there is no "0" in the bond_status_list
    ##########################################################################
    radius_scale_factor_temp = radius_scale_factor
    temp_poscar_layer_dict = poscar_layer_dict
    condition_1 = len(poscar_layer_dict['bond_status_list']) > 1
    condition_2 = 0 in poscar_layer_dict['bond_status_list']

    num_iter = 10
    radius_scale_factor_list = radius_scale_factor + 0.05 * np.array([sum(np.arange(1, x)) for x in np.arange(2, num_iter + 2)])

    if condition_1 and condition_2:
        for radius_scale_factor in radius_scale_factor_list:
            if suppress_warning == False:
                print('# WARNING #2104121511 (from vasp_tools): radius_scale_factor too small. Set radius_scale_factor = ', '{:.4f}'.format(radius_scale_factor))
            poscar_layer_dict = poscar_layer_status(poscar_dict = poscar_dict, poscar_direction_dict = poscar_direction_dict, delta = delta, radius_style = radius_style, radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning= suppress_warning)
            if 0 not in poscar_layer_dict['bond_status_list']:
                poscar_layer_dict['radius_scale_factor'] = radius_scale_factor
                break
        if radius_scale_factor == radius_scale_factor_list[-1] and 0 in poscar_layer_dict['bond_status_list']:
            # There are atoms in a layer too far from other adjacent layers
            poscar_layer_dict = temp_poscar_layer_dict
            radius_scale_factor = radius_scale_factor_temp
    elif 0 in poscar_layer_dict['bond_status_list'] and len(poscar_layer_dict['bond_status_list']) == 1:
        pass
    elif 2 in poscar_layer_dict['bond_status_list'] and len(poscar_layer_dict['bond_status_list']) == 1:
        pass
        ########################################################
        # Check single layer configuration
        ########################################################
        # Build another layer in the next period and check the bonding status of atoms in the two layers
        # align atoms to the bottom of the model 
        #print('ENTERING single layeri===============================')
        poscar_dict_temp = deepcopy(poscar_dict)
        # atom_species_arr, n_atoms, pos_arr, atom_name_list
        ############################################################
        # Expand model in a certain direction by creating supercell
        ############################################################
        fpath, fname = os.path.split(poscar_dict['file_path'])
        output_poscar_file_path_temp = os.path.join(fpath, 'POSCAR_temp')
        if poscar_direction_dict['pos_arr_direct_column'] == 0:
            num_xyz_list = [4, 1, 1]
        elif poscar_direction_dict['pos_arr_direct_column'] == 1:
            num_xyz_list = [1, 4, 1]
        elif poscar_direction_dict['pos_arr_direct_column'] == 2:
            num_xyz_list = [1, 1, 4]
        vasp_build.create_supercell(poscar_file_path = poscar_dict['file_path'],
            num_xyz_list = num_xyz_list,
            output_poscar_file_path = output_poscar_file_path_temp,
            ) 
        poscar_dict_temp = vasp_read.read_poscar(output_poscar_file_path_temp)
        poscar_dict_temp['atom_number_ortho_' + poscar_direction_dict['number_order_text']] = funcs.atom_number_ortho(
            atom_key_arr = poscar_dict_temp['atom_key_arr'],
            pos_arr_cartesian= poscar_dict_temp['pos_arr'][:,3:6],
            delta = delta,
            number_order = poscar_direction_dict['number_order_text']
            )
        poscar_direction_dict_temp = poscar_direction_params(poscar_dict_temp, direction = poscar_direction_dict['number_order_text'])
        poscar_layer_dict_temp = poscar_layer_status(poscar_dict = poscar_dict_temp, poscar_direction_dict = poscar_direction_dict_temp, delta = delta, radius_style = radius_style, radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning= suppress_warning)
        #os.remove(output_poscar_file_path_temp)
        if 0 in poscar_layer_dict_temp['bond_status_list']:
            poscar_layer_dict['vacuum_exist'] = True
            poscar_layer_dict['bond_status_list'][0] = 0
            poscar_layer_dict['crystal_width'] = 0
            poscar_layer_dict['vacuum_width'] = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['crystal_width']
        else:
            poscar_layer_dict['vacuum_exist'] = False
            poscar_layer_dict['bond_status_list'][0] = 2
            poscar_layer_dict['crystal_width'] = poscar_direction_dict['box_len_ortho']
            poscar_layer_dict['vacuum_width'] = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['crystal_width']
        ##print('EXITING single layeri===============================')

    # Check if there is only 1 or -1 in the bond_status_list
    condition_3 = (0 not in poscar_layer_dict['bond_status_list'])
    condition_4 = ((1 in poscar_layer_dict['bond_status_list'] and -1 not in poscar_layer_dict['bond_status_list']) or (-1 in poscar_layer_dict['bond_status_list'] and 1 not in poscar_layer_dict['bond_status_list'])) and condition_3
    radius_scale_factor_list = radius_scale_factor + 0.05 * np.array([sum(np.arange(1, x)) for x in np.arange(2, num_iter + 2)])
    if condition_1 and condition_4:
        for radius_scale_factor in radius_scale_factor_list:
            if suppress_warning == False:
                print('# WARNING #2104121511 (from vasp_tools): radius_scale_factor too small. Set radius_scale_factor = ', '{:.4f}'.format(radius_scale_factor))
            poscar_layer_dict = poscar_layer_status(poscar_dict = poscar_dict, poscar_direction_dict = poscar_direction_dict, delta = delta, radius_style = radius_style, radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning= suppress_warning)
            condition_temp1 = (1 in poscar_layer_dict['bond_status_list'] and -1 in poscar_layer_dict['bond_status_list'])
            condition_temp2 = (1 not in poscar_layer_dict['bond_status_list'] and -1 not in poscar_layer_dict['bond_status_list'])
            condition_temp3 = condition_temp1 or condition_temp2
            condition_temp4 = (1 in poscar_layer_dict['bond_status_list'] and -1 not in poscar_layer_dict['bond_status_list'])
            condition_temp5 = (1 not in poscar_layer_dict['bond_status_list'] and -1 in poscar_layer_dict['bond_status_list'])
            condition_temp6 = condition_temp4 or condition_temp5
            if condition_temp3:
                poscar_layer_dict['radius_scale_factor'] = radius_scale_factor
                break
        if radius_scale_factor == radius_scale_factor_list[-1] and condition_temp6:
            # There are atoms in a layer too far from other adjacent layers
            poscar_layer_dict = temp_poscar_layer_dict
            radius_scale_factor = radius_scale_factor_temp

    poscar_layer_dict['bond_status_list'] = poscar_layer_dict['bond_status_list']
    poscar_layer_dict['periodicity_status_list'] = poscar_layer_dict['periodicity_status_list']
    return poscar_layer_dict
