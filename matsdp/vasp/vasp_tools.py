# -*- coding: utf-8 -*-
def ordered_incar_dict(incar_dict):
    '''
    rearrange the order of the keys of the incar_dict according to vasp_default.incar_default()
    '''  
    args_dict = locals()
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
    args_dict = locals()
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

    ############################################################################
    #WARNING: random wavefunctions but no delay for mixing, default for NELMDL
    ############################################################################

    ##################################################################
    #WARNING in EDDRMM: call to ZHEGV failed, returncode =   6  3 27
    ##################################################################

    ######################################################################################################
    # RNING in EDDRMM: call to ZHEGV failed, returncode =   6  3 31
    # num prob  num prob  num prob  num prob  num prob  num prob  num prob  num prob  num prob  num prob
    ######################################################################################################

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
    # |      The number of bands has been changed from the values supplied          |
    # |      in the INCAR file. This is a result of running the parallel version.   |
    # |      The orbitals not found in the WAVECAR file will be initialized with    |
    # |      random numbers, which is usually adequate. For correlated              |
    # |      calculations, however, you should redo the groundstate calculation.    |
    # |      I found NBANDS    =      172  now  NBANDS  =     192                   |
    # |                                                                             |
    #  -----------------------------------------------------------------------------
    #################################################################################################

    #################################################################################################
    # WARNING: dimensions on CHGCAR file are different
    #################################################################################################
    # This could be caused by a change in lattice constant for the POSCAR and the CHGCAR files

    pass
    return 0

def check_error(job_dir):
    '''
    check ERROR of the vasp job
    it returns a error list
    NOTE: errors are those prompt that cannot be ignored, and the job has to stop due to an error, which is different from warnings
    '''
    args_dict = locals()
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
    # |      VASP internal routines  have requested a change of the k-point set.    |
    # |      Unfortunately this is only possible if NPAR=number of nodes.           |
    # |      Please remove the tag NPAR from the INCAR file and restart the         |
    # |      calculations.                                                          |
    # |                                                                             |
    # |      ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ..., BYE!!! <----       |
    # |                                                                             |
    #  -----------------------------------------------------------------------------
    #################################################################################################



    #################################################################################################
    # ERROR: charge density could not be read from file CHGCAR for ICHARG>10
    #################################################################################################

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
    # |      VASP internal routines  have requested a change of the k-point set.    |
    # |      Unfortunately this is only possible if NPAR=number of nodes.           |
    # |      Please remove the tag NPAR from the INCAR file and restart the         |
    # |      calculations.                                                          |
    # |                                                                             |
    # |      ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ..., BYE!!! <----       |
    # |                                                                             |
    #  -----------------------------------------------------------------------------
    #################################################################################################
    error_kwds_list_dict['VASP internal routines  have requested a change of the k-point set'] = ['VASP internal routines  have requested a change of the k-point set',]


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
    args_dict = locals()
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
        if error_type == 'VASP internal routines  have requested a change of the k-point set':
            print('Error type ' + error_type + ' is detected for the directory: ' + i_dir_path + '\n' + 
                  'This error usually appears in the DFPT calculation (IBRION = 8).' + '\n' + 
                  'SOLUTION: ' + 'n' + 
                  '    (1) This may be caused by specifying NPAR value in the INCAR, Try not setting the NPAR value in the INCAR file. The VASP program then set the NPAR to its default value' + '\n' + 
                  '    (2) When the kpoint changed, it may cause an increase of k-points, which gives rise to a demand of large memory. In this case, please check your memory of the current job.'
                  )
                  
    return 0

def vasp_job_finished(outcar_file_path, suppress_warning = True):
    '''
    check whether the VASP job has been finished or not
    '''
    args_dict = locals()
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

def job_status(job_parent_dir = './', write_abspath = False, write_band_gap = False, write_e_fermi = False, write_formula = False, write_volume = False, write_area = False, sort_by = None, screen_by = None, suppress_warning = False):
    '''
    Check the job status for multiple jobs
    job_parent_dir: This is the parent directory which contains multiple VASP jobs
    '''
    args_dict = locals()
    import os
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from . import vasp_analyze
    from . import vasp_default
    from .. import default_params
    import multiprocessing
    from multiprocessing import Pool
    import time
    import math
    import numpy as np

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)
    incar_default_dict = vasp_default.incar_default() 

    job_parent_dir = os.path.abspath(job_parent_dir)

    if funcs.dir_status(job_parent_dir) == 0:
        quit()
    #############################################
    # initialize parameters
    #############################################
    # write_abspath tag
    if not write_abspath == True and (not screen_by is None or not sort_by is None):
        if suppress_warning == False:
            print('#WARNING #2104251041 (from vasp_tools): screen_by or sort_by tag is not None, the write_abspath is automatically set to write_abspath = True')
        write_abspath = True

    if not sort_by is None:
        if not write_band_gap == True and all(x in sort_by for x in ['band', 'gap']):
            if suppress_warning == False:
                print('#WARNING #2104251426 (from vasp_tools): sort by band gap is designated, the write_band_gap tag is automatically set to write_band_gap = True')
            write_band_gap = True
        if not write_e_fermi == True and all(x in sort_by for x in ['fermi']):
            if suppress_warning == False:
                print('#WARNING #2104251616 (from vasp_tools): sort by e_fermi is designated, the write_e_fermi tag is automatically set to write_e_fermi = True')
            write_e_fermi = True
        if not write_formula == True and all(x in sort_by for x in ['formula']):
            if suppress_warning == False:
                print('#WARNING #2104251617 (from vasp_tools): sort by formula is designated, the write_formula tag is automatically set to write_formula = True')
            write_formula = True
        if not write_volume == True and all(x in sort_by for x in ['volume']):
            if suppress_warning == False:
                print('#WARNING #2104252001 (from vasp_tools): sort by volume is designated, the write_volume tag is automatically set to write_volume = True')
            write_volume = True
        if not write_area == True and all(x in sort_by for x in ['area']):
            if suppress_warning == False:
                print('#WARNING #2104252002 (from vasp_tools): sort by area is designated, the write_area tag is automatically set to write_area = True')
            write_area = True

    if not screen_by is None:
        if not write_band_gap == True and all(x in screen_by for x in ['band', 'gap']):
            if suppress_warning == False:
                print('#WARNING #2104251442 (from vasp_tools): screen by band gap is designated, the write_band_gap tag is automatically set to write_band_gap = True')
            write_band_gap = True
        if not write_e_fermi == True and all(x in screen_by for x in ['fermi']):
            if suppress_warning == False:
                print('#WARNING #2104251617 (from vasp_tools): screen by e_fermi is designated, the write_e_fermi tag is automatically set to write_e_fermi = True')
            write_e_fermi = True
        if not write_formula == True and all(x in screen_by for x in ['formula']):
            if suppress_warning == False:
                print('#WARNING #2104251618 (from vasp_tools): screen by formula is designated, the write_formula tag is automatically set to write_formula = True')
            write_formula = True
        if not write_volume == True and all(x in screen_by for x in ['volume']):
            if suppress_warning == False:
                print('#WARNING #2104252002 (from vasp_tools): screen by volume is designated, the write_volume tag is automatically set to write_volume = True')
            write_volume = True
        if not write_area == True and all(x in screen_by for x in ['area']):
            if suppress_warning == False:
                print('#WARNING #2104252003 (from vasp_tools): screen by area is designated, the write_area tag is automatically set to write_area = True')
            write_area = True

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
    job_status_arr = np.array([None] * len(folder_path_list) * 7)
    job_status_arr.shape = len(folder_path_list), 7
    
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
    if write_band_gap == True:
        temp_str = temp_str +  funcs.str_format('Egap', max_len = 8)
    if write_e_fermi == True:
        temp_str = temp_str +  funcs.str_format('Efermi', max_len = 8)
    if write_formula == True:
        temp_str = temp_str + funcs.str_format('atoms', max_len = 6) + funcs.str_format('formula', max_len = 21)
    if write_volume == True:
        temp_str = temp_str +  funcs.str_format('volume', max_len = 11)
    if write_area == True:
        # area_star is the cross product of the two vectors with the smallest lengths
        temp_str = temp_str +  funcs.str_format('area_ab', max_len = 10) + funcs.str_format('area_bc', max_len = 10) + funcs.str_format('area_ca', max_len = 10) + funcs.str_format('area_star', max_len = 10)
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
        temp_str = temp_str + funcs.str_format('-' * len('Egap'), max_len = 8)
    if write_e_fermi == True:
        temp_str = temp_str + funcs.str_format('-' * len('Efermi'), max_len = 8)
    if write_formula == True:
        temp_str = temp_str + funcs.str_format('-' * len('atoms'), max_len = 6) + funcs.str_format('-' * len('formula'), max_len = 21)
    if write_volume == True:
        temp_str = temp_str + funcs.str_format('-' * len('volume'), max_len = 11)
    if write_area == True:
        temp_str = temp_str + funcs.str_format('-' * len('area_ab'), max_len = 10) + funcs.str_format('-' * len('area_bc'), max_len = 10) + funcs.str_format('-' * len('area_ca'), max_len = 10) + funcs.str_format('-' * len('area_star'), max_len = 10)
    temp_str = temp_str + '\n'

    # check each job
    args = []
    for i_indx in range(0, len(folder_path_list)):
        #primitive_structure = data[i_indx]['cif']
        #material_id_list.append(data[i_indx]['material_id'])
        args.append((folder_path_list, folder_path_str_list, i_indx, job_status_file_path, write_band_gap, write_e_fermi, write_formula, write_volume, write_area))
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

    results_list = [None] * len(folder_path_list)
    if parallel_run == True:
        with Pool(cores) as p:
            results_list = p.starmap(single_job_status, args)
    elif parallel_run == False:
        for i_indx in range(0, len(folder_path_list)):
            # use unpacked arguments in the tuple
            results_list[i_indx] = single_job_status(*args[i_indx])

    ##########################################################
    # Screen the jobs according to certain criteria
    ##########################################################
    job_status_kwd_lowercase_list = ['job_dir', 'status', 'energy_without_entropy', 'toten', 'energy(sigma->0)', 'elapsed_time', 'egap', 'efermi', 'atoms', 'formula', 'volume', 'area_ab', 'area_bc', 'area_ca', 'area_star']
    job_status_kwd_list = ['folder_path_str', 'job_status_str', 'energy_without_entropy', 'TOTEN', 'energy(sigma->0)', 'elapsed_time_hour', 'band_gap', 'e_fermi', 'num_atoms', 'formula', 'volume', 'area_ab', 'area_bc', 'area_ca', 'area_star']

    float_tag_list = ['energy_without_entropy', 'TOTEN', 'energy(sigma->0)', 'elapsed_time_hour', 'band_gap', 'num_atoms', 'e_fermi', 'volume', 'area_ab', 'area_bc', 'area_ca', 'area_star']
    str_tag_list = ['folder_path_str', 'job_status_str', 'formula']

    screen_list = []
    if not screen_by is None:
        if '!=' in screen_by:
            temp_logic_val = False
        else:
            temp_logic_val = True

        for i_indx in range(0, len(folder_path_list)):
            if '=' in screen_by and all(x not in screen_by for x in ['>', '<']):
                screen_key = screen_by.split('=')[0].replace('\t', '    ').strip('\n').strip()
                screen_value = screen_by.split('=')[1].replace('\t', '    ').strip('\n').strip()
                if screen_key in incar_default_dict.keys() and not results_list[i_indx]['incar_dict'] is None:
                    if results_list[i_indx]['incar_dict']['read_status'] == 1 and screen_key in results_list[i_indx]['incar_dict']:
                        results_value = results_list[i_indx]['incar_dict'][screen_key]
                    else:
                        results_value = None
                    if not (isinstance(results_value, int) or isinstance(results_value, float)):
                        results_value = str(results_value).replace('\t', '    ').strip('\n').strip()
                    else:
                        screen_value = float(screen_value)
                    screen_status = (screen_value == results_value)
                    if screen_status == temp_logic_val:
                        screen_list.append(results_list[i_indx])
                if screen_key.lower() in job_status_kwd_lowercase_list:
                    temp_indx = job_status_kwd_lowercase_list.index(screen_key.lower())
                    if job_status_kwd_list[temp_indx] in str_tag_list:
                        results_value = results_list[i_indx][job_status_kwd_list[temp_indx]]
                    elif job_status_kwd_list[temp_indx] in float_tag_list:
                        if results_list[i_indx][job_status_kwd_list[temp_indx]] != '--':
                            results_value = float(results_list[i_indx][job_status_kwd_list[temp_indx]])
                        else:
                            results_value = results_list[i_indx][job_status_kwd_list[temp_indx]]
                        if screen_value != '--':
                            screen_value = float(screen_value)
                    screen_status = (screen_value == results_value)
                    if screen_status == temp_logic_val:
                        screen_list.append(results_list[i_indx])
                if screen_key == 'dir_level':
                    results_value = float(folder_level_list[i_indx])
                    if not (isinstance(results_value, int) or isinstance(results_value, float)):
                        results_value = str(results_value).replace('\t', '    ').strip('\n').strip()
                    else:
                        screen_value = float(screen_value)
                    screen_status = (screen_value == results_value)
                    if screen_status == temp_logic_val:
                        screen_list.append(results_list[i_indx])
            elif '>' in screen_by:
                if '>=' not in screen_by:
                    screen_key = screen_by.split('>')[0].replace('\t', '    ').strip('\n').strip()
                    screen_value = float(screen_by.split('>')[1].replace('\t', '    ').strip('\n').strip())
                elif '>=' in screen_by:
                    screen_key = screen_by.split('>=')[0].replace('\t', '    ').strip('\n').strip()
                    screen_value = float(screen_by.split('>=')[1].replace('\t', '    ').strip('\n').strip())
                if screen_key in incar_default_dict.keys() and not results_list[i_indx]['incar_dict'] is None:
                    if results_list[i_indx]['incar_dict']['read_status'] == 1 and screen_key in results_list[i_indx]['incar_dict']:
                        results_value = results_list[i_indx]['incar_dict'][screen_key]
                    else:
                        results_value = None
                    screen_status = (results_value > screen_value)
                    if '>=' in screen_by:
                        screen_status = (results_value >= screen_value)
                    if screen_status == True:
                        screen_list.append(results_list[i_indx])
                if screen_key.lower() in job_status_kwd_lowercase_list:
                    temp_indx = job_status_kwd_lowercase_list.index(screen_key.lower())
                    if job_status_kwd_list[temp_indx] in str_tag_list:
                        results_value = results_list[i_indx][job_status_kwd_list[temp_indx]]
                    elif job_status_kwd_list[temp_indx] in float_tag_list:
                        if results_list[i_indx][job_status_kwd_list[temp_indx]] != '--':
                            results_value = float(results_list[i_indx][job_status_kwd_list[temp_indx]])
                        else:
                            results_value = results_list[i_indx][job_status_kwd_list[temp_indx]]
                        if screen_value != '--':
                            screen_value = float(screen_value)
                    screen_status = False
                    if type(results_value) == type(screen_value):
                        screen_status = (results_value > screen_value)
                        if '>=' in screen_by:
                            screen_status = (results_value >= screen_value)
                        if screen_status == True:
                            screen_list.append(results_list[i_indx])
                if screen_key == 'dir_level':
                    results_value = float(folder_level_list[i_indx])
                    if not (isinstance(results_value, int) or isinstance(results_value, float)):
                        results_value = str(results_value).replace('\t', '    ').strip('\n').strip()
                    else:
                        screen_value = float(screen_value)
                    screen_status = False
                    if type(results_value) == type(screen_value):
                        screen_status = (results_value > screen_value)
                    if '>=' in screen_by:
                        screen_status = (results_value >= screen_value)
                    if screen_status == True:
                        screen_list.append(results_list[i_indx])
            elif '<' in screen_by:
                if '<=' not in screen_by:
                    screen_key = screen_by.split('<')[0].replace('\t', '    ').strip('\n').strip()
                    screen_value = float(screen_by.split('<')[1].replace('\t', '    ').strip('\n').strip())
                if '<=' in screen_by:
                    screen_key = screen_by.split('<=')[0].replace('\t', '    ').strip('\n').strip()
                    screen_value = float(screen_by.split('<=')[1].replace('\t', '    ').strip('\n').strip())
                if screen_key in incar_default_dict.keys() and not results_list[i_indx]['incar_dict'] is None:
                    if results_list[i_indx]['incar_dict']['read_status'] == 1 and screen_key in results_list[i_indx]['incar_dict']:
                        results_value = results_list[i_indx]['incar_dict'][screen_key]
                    else:
                        results_value = None
                    screen_status = (results_value < screen_value)
                    if '<=' in screen_by:
                        screen_status = (results_value <= screen_value)
                    if screen_status == True:
                        screen_list.append(results_list[i_indx])
                if screen_key.lower() in job_status_kwd_lowercase_list:
                    temp_indx = job_status_kwd_lowercase_list.index(screen_key.lower())
                    if job_status_kwd_list[temp_indx] in str_tag_list:
                        results_value = results_list[i_indx][job_status_kwd_list[temp_indx]]
                    elif job_status_kwd_list[temp_indx] in float_tag_list:
                        if results_list[i_indx][job_status_kwd_list[temp_indx]] != '--':
                            results_value = float(results_list[i_indx][job_status_kwd_list[temp_indx]])
                        else:
                            results_value = results_list[i_indx][job_status_kwd_list[temp_indx]]
                        if screen_value != '--':
                            screen_value = float(screen_value)
                    screen_status = False
                    if type(results_value) == type(screen_value):
                        screen_status = (results_value < screen_value)
                        if '<=' in screen_by:
                            screen_status = (results_value <= screen_value)
                        if screen_status == True:
                            screen_list.append(results_list[i_indx])
                if screen_key == 'dir_level':
                    results_value = float(folder_level_list[i_indx])
                    if not (isinstance(results_value, int) or isinstance(results_value, float)):
                        results_value = str(results_value).replace('\t', '    ').strip('\n').strip()
                    else:
                        screen_value = float(screen_value)
                    screen_status = False
                    if type(results_value) == type(screen_value):
                        screen_status = (results_value < screen_value)
                    if '<=' in screen_by:
                        screen_status = (results_value <= screen_value)
                    if screen_status == True:
                        screen_list.append(results_list[i_indx])

        results_list = screen_list

    ##########################################################
    # Sort the jobs according to certain criteria
    ##########################################################
    if not sort_by is None:
        i_tag = None
        if all(x in sort_by for x in ['job', 'dir']):
            i_tag = 'folder_path_str'
            sorted_arg_arr = np.argsort(results_list['folder_path_str'])
        elif all(x in sort_by for x in ['status']):
            i_tag = 'job_status_str'
            sorted_arg_arr = np.argsort([x['job_status_str'] for x in results_list])
        elif all(x in sort_by for x in ['formula']):
            i_tag = 'formula'
            sorted_arg_arr = np.argsort([x['formula'] for x in results_list])
        elif all(x in sort_by for x in ['energy', 'without', 'entropy']):
            i_tag = 'energy_without_entropy'
        elif all(x in sort_by for x in ['TOTEN']):
            i_tag = 'TOTEN'
        elif all(x in sort_by for x in ['energy', 'sigma', '0']):
            i_tag = 'energy(sigma->0)'
        elif all(x in sort_by for x in ['elapsed', 'time']):
            i_tag = 'elapsed_time_hour'
        elif all(x in sort_by for x in ['band', 'gap']):
            i_tag = 'band_gap'
        elif all(x in sort_by for x in ['atoms']):
            i_tag = 'num_atoms'
        elif all(x in sort_by for x in ['fermi']):
            i_tag = 'e_fermi'
        elif all(x in sort_by for x in ['volume']):
            i_tag = 'volume'
        elif all(x in sort_by for x in ['area_ab']):
            i_tag = 'area_ab'
        elif all(x in sort_by for x in ['area_bc']):
            i_tag = 'area_bc'
        elif all(x in sort_by for x in ['area_ca']):
            i_tag = 'area_ca'
        elif all(x in sort_by for x in ['area_star']):
            i_tag = 'area_star'
        else:
            pass

        if i_tag in float_tag_list:
            # sort by float value
            float_list = [x for x in results_list if x[i_tag] != '--']
            float_value_list = [x[i_tag] for x in results_list if x[i_tag] != '--']
            dash_list = [x for x in results_list if x[i_tag] == '--']
            sorted_arg_arr = np.argsort(float_value_list)
            sort_list = [None] * len(float_list)
            for i_indx in range(len(float_list)):
                sort_list[i_indx] = float_list[sorted_arg_arr[i_indx]]
            sort_list = sort_list + dash_list
            results_list = sort_list
        elif i_tag in str_tag_list:
            # sort by alphabetical order
            sort_list = [None] * len(results_list)
            for i_indx in range(len(results_list)):
                sort_list[i_indx] = results_list[sorted_arg_arr[i_indx]]
            results_list = sort_list

    for i_indx in range(0, len(results_list)):
        temp_str = temp_str + results_list[i_indx]['full_info_str']

    with open(job_status_file_path, 'w') as f:
        f.write(temp_str)
    funcs.write_log(logfile,
        'vasp_tools.job_status(\n' +
        '    job_parent_dir = ' + 'r\'' + job_parent_dir + '\'' + '\n' + 
        '    write_abspath = ' + str(write_abspath) + '\n' + 
        '    write_band_gap = ' + str(write_band_gap) + '\n' + 
        '    write_e_fermi = ' + str(write_e_fermi) + '\n' + 
        '    write_formula = ' + str(write_formula) + '\n' + 
        '    write_volume = ' + str(write_volume) + '\n' + 
        '    write_area = ' + str(write_area) + '\n' + 
        '    sort_by = ' + str(sort_by) + '\n' + 
        '    screen_by = ' + str(screen_by) + '\n' + 
        '    suppress_warning = ' + str(suppress_warning) + '\n' + 
        '    )\n' +
        '#############\n'
        )
    return results_list

def single_job_status(folder_path_list, folder_path_str_list, i_folder_indx, job_status_file_path, write_band_gap, write_e_fermi, write_formula, write_volume, write_area):
    '''
    Get the job status
    '''
    args_dict = locals()
    import os
    from .. import funcs
    from .. import convert
    from . import vasp_read
    from . import vasp_analyze

    single_job_status_dict = {}
    single_job_status_dict['folder_path_str'] = None 
    single_job_status_dict['job_status_str'] = None
    single_job_status_dict['energy_without_entropy'] = None
    single_job_status_dict['TOTEN'] = None
    single_job_status_dict['energy(sigma->0)'] = None
    single_job_status_dict['elapsed_time_hour'] = None
    single_job_status_dict['band_gap'] = None
    single_job_status_dict['e_fermi'] = None
    single_job_status_dict['num_atoms'] = None
    single_job_status_dict['formula'] = None
    single_job_status_dict['volume'] = None
    single_job_status_dict['area_ab'] = None
    single_job_status_dict['area_bc'] = None
    single_job_status_dict['area_ca'] = None
    single_job_status_dict['area_star'] = None
    single_job_status_dict['full_info_str'] = None
    single_job_status_dict['incar_dict'] = None

    max_length_folder_path_str = max([len(x) for x in folder_path_str_list])
#for i_folder_indx in range(0, len(folder_path_list)):
    output_exist = False
    # check job status
    job_status = 'unfinished'
    job_status_str = 'undone'
    poscar_file_path = os.path.join(folder_path_list[i_folder_indx], 'POSCAR')
    poscar_file_status = funcs.file_status(poscar_file_path)
    if funcs.check_job_type(folder_path_list[i_folder_indx])['VASP'] == False:
        job_status = 'not_vasp_job'
        job_status_str = 'not_VASP_job'
        num_atoms = '--'
        formula = '--'
        volume = '--'
        area_ab = '--'
        area_bc = '--'
        area_ca = '--'
        area_star = '--'
    else:
        if poscar_file_status != 1:
            num_atoms = '--'
            formula = '--'
            volume = '--'
            area_ab = '--'
            area_bc = '--'
            area_ca = '--'
            area_star = '--'
        else:
            poscar_dict = vasp_read.read_poscar(poscar_file_path)
            num_atoms = poscar_dict['n_atoms']
            formula = poscar_dict['alphabetical_formula']
            volume = '{:.4f}'.format(poscar_dict['volume'])
            area_ab = '{:.4f}'.format(poscar_dict['area_ab'])
            area_bc = '{:.4f}'.format(poscar_dict['area_bc'])
            area_ca = '{:.4f}'.format(poscar_dict['area_ca'])
            area_star = '{:.4f}'.format(poscar_dict['area_star'])
    energy_without_entropy = None
    temp_str = ''
    incar_dict = None
    for file0 in os.listdir(folder_path_list[i_folder_indx]):
        if 'OUTCAR' == file0:
            output_exist = True
            output_file_path = os.path.join(folder_path_list[i_folder_indx], 'OUTCAR')
            oszicar_file_path = os.path.join(folder_path_list[i_folder_indx], 'OSZICAR')
            kpoints_file_path = os.path.join(folder_path_list[i_folder_indx], 'KPOINTS')
            eigenval_file_path = os.path.join(folder_path_list[i_folder_indx], 'EIGENVAL')
            incar_file_path = os.path.join(folder_path_list[i_folder_indx], 'INCAR')
            if vasp_job_finished(output_file_path) == True:
                job_status = 'finished'
                job_status_str = 'finished'
                outcar_params_dict = vasp_read.read_outcar(output_file_path)
                elapsed_time_hour = convert.time_converter(hour = 0, minute = 0, second = outcar_params_dict['elapsed_time'], unit = 'hour')
                kpoints_dict = vasp_read.read_kpoints(kpoints_file_path)
                oszicar_dict = vasp_read.read_oszicar(oszicar_file_path)
                incar_dict = vasp_read.read_incar(incar_file_path)

                # check if the parameters in the outcar_params_dict are valid
                valid_outcar_par = True
                temp_outcar_par_list = ['file_path', 'NELM', 'energy_without_entropy', 'TOTEN', 'energy(sigma->0)', 'elapsed_time', 'e_fermi']
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
                band_gap = None
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
                if write_e_fermi == True and valid_outcar_par == True:
                    e_fermi = '{:.4f}'.format(outcar_params_dict['e_fermi'])
                else:
                    e_fermi = '--'
                if write_formula == True and poscar_file_status == 1:
                    num_atoms = poscar_dict['n_atoms']
                    formula = poscar_dict['alphabetical_formula']
                if write_volume == True and valid_outcar_par == True and poscar_file_status == 1:
                    volume = '{:.4f}'.format(poscar_dict['volume'])
                else:
                    volume = '--'
                if write_area == True and valid_outcar_par == True and poscar_file_status == 1:
                    area_ab = '{:.4f}'.format(poscar_dict['area_ab'])
                    area_bc = '{:.4f}'.format(poscar_dict['area_bc'])
                    area_ca = '{:.4f}'.format(poscar_dict['area_ca'])
                    area_star = '{:.4f}'.format(poscar_dict['area_star'])
                else:
                    area_ab = '--'
                    area_bc = '--'
                    area_ca = '--'
                    area_star = '--'
            continue
  
    if job_status == 'finished' and valid_outcar_par and valid_oszicar_par:
        single_job_status_dict['folder_path_str'] = folder_path_str_list[i_folder_indx] 
        single_job_status_dict['job_status_str'] = job_status_str
        single_job_status_dict['energy_without_entropy'] = outcar_params_dict['energy_without_entropy']
        single_job_status_dict['TOTEN'] = outcar_params_dict['TOTEN']
        single_job_status_dict['energy(sigma->0)'] = outcar_params_dict['energy(sigma->0)']
        single_job_status_dict['elapsed_time_hour'] = elapsed_time_hour
        single_job_status_dict['band_gap'] = band_gap
        single_job_status_dict['e_fermi'] = e_fermi
        single_job_status_dict['num_atoms'] = num_atoms
        single_job_status_dict['formula'] = formula
        single_job_status_dict['area_ab'] = area_ab
        single_job_status_dict['area_bc'] = area_bc
        single_job_status_dict['area_ca'] = area_ca
        single_job_status_dict['area_star'] = area_star
        single_job_status_dict['volume'] = volume
        
        temp_str = temp_str + (
            funcs.str_format(folder_path_str_list[i_folder_indx], max_len = max_length_folder_path_str) + '   '
            + funcs.str_format(job_status_str, max_len = 13) 
            + funcs.float_format(outcar_params_dict['energy_without_entropy'], max_len = 16, format_str = '{:.4f}')
            + funcs.float_format(outcar_params_dict['TOTEN'], max_len = 13, format_str = '{:.4f}')
            + funcs.float_format(outcar_params_dict['energy(sigma->0)'], max_len = 17, format_str = '{:.4f}')
            + funcs.float_format(elapsed_time_hour, max_len = 14, format_str = '{:.4f}')
            )
        if write_band_gap == True:
            temp_str = temp_str + funcs.str_format(band_gap, max_len = 8)
        if write_e_fermi == True:
            temp_str = temp_str + funcs.str_format(e_fermi, max_len = 8)
        if write_formula == True:
            temp_str = temp_str + funcs.str_format(num_atoms, max_len = 6) + funcs.str_format(formula, max_len = 21)
        if write_volume == True:
            temp_str = temp_str + funcs.str_format(volume, max_len = 11)
        if write_area == True:
            temp_str = temp_str + funcs.str_format(area_ab, max_len = 10) + funcs.str_format(area_bc, max_len = 10) + funcs.str_format(area_ca, max_len = 10) + funcs.str_format(area_star, max_len = 10)
        temp_str = temp_str + ' \n'
    else:
        single_job_status_dict['folder_path_str'] = folder_path_str_list[i_folder_indx] 
        single_job_status_dict['job_status_str'] = job_status_str
        single_job_status_dict['energy_without_entropy'] = '--'
        single_job_status_dict['TOTEN'] = '--'
        single_job_status_dict['energy(sigma->0)'] = '--'
        single_job_status_dict['elapsed_time_hour'] = '--'
        single_job_status_dict['band_gap'] = '--'
        single_job_status_dict['e_fermi'] = '--'
        single_job_status_dict['num_atoms'] = num_atoms
        single_job_status_dict['formula'] = formula
        single_job_status_dict['area_ab'] = area_ab
        single_job_status_dict['area_bc'] = area_bc
        single_job_status_dict['area_ca'] = area_ca
        single_job_status_dict['area_star'] = area_star
        single_job_status_dict['volume'] = volume
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
        if write_e_fermi == True:
            temp_str = temp_str + funcs.str_format('--', max_len = 8)
        if write_formula == True:
            temp_str = temp_str + funcs.str_format(num_atoms, max_len = 6) + funcs.str_format(formula, max_len = 21)
        if write_volume == True:
            temp_str = temp_str + funcs.str_format(volume, max_len = 11)
        if write_area == True:
            temp_str = temp_str + funcs.str_format(area_ab, max_len = 10) + funcs.str_format(area_bc, max_len = 10) + funcs.str_format(area_ca, max_len = 10) + funcs.str_format(area_star, max_len = 10)
        temp_str = temp_str + ' \n'
    single_job_status_dict['full_info_str'] = temp_str
    single_job_status_dict['incar_dict'] = incar_dict
    return single_job_status_dict

def check_file_type(file_path):
    '''
    check the file type
    file_type: string format. It defines the file type, e.g. 'EIGENVAL', 'PROCAR'.
    the separator of each blocks in a line is white space.
    '''
    args_dict = locals()
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
    args_dict = locals()
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

                    volume_poscar = poscar_dict['volume']
                    volume_contcar = contcar_dict['volume']
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

                '{:.4f}'.format(volume_poscar) + ' ' * (13 - len('{:.4f}'.format(volume_poscar))) + 
                '{:.4f}'.format(volume_contcar) + ' ' * (13 - len('{:.4f}'.format(volume_contcar))) + 

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
    ##################################
    # Determine the args string
    ##################################
    log_str = ''
    func_name = 'vasp_tools.get_latt_param'
    args_str = func_name + '(' + '\n'
    for i_arg in args_dict.keys():
        arg_value = args_dict[i_arg]
        if isinstance(arg_value,str):
            arg_value_str = '\'' + arg_value + '\''
        else:
            arg_value_str = str(arg_value)

        if i_arg == 'job_parent_dir':
            arg_value_str = 'r\'' + job_parent_dir + '\''
        args_str += '    ' + i_arg + ' = ' + arg_value_str + ',\n'
    args_str += '    )\n'
    args_str += '################################################\n'
    log_str += args_str
    funcs.write_log(logfile, log_str)
    return temp_str

def convert_coord_system(poscar_file_path, mode = 'direct2cartesian'):
    '''
    Convert from Direct coordinate to Cartesian coordinate
    mode: 'direct2cartesian' or 'cartesian2direct'
    '''
    args_dict = locals()
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
    args_dict = locals()
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
    '''
    Get the tolerance value of the layer distance for a specific POSCAR file
    '''
    args_dict = locals()
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
    args_dict = locals()
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
        box_len_ortho = abs(poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
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
        box_len_ortho = abs(poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
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
        box_len_ortho = abs(poscar_dict['vec_a'] + poscar_dict['vec_b'] + poscar_dict['vec_c']).dot(ortho_vector)
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

def poscar_layer_params(poscar_dict, poscar_direction_dict, criteria = 'auto', delta = 0.05, layer_dist_tolerance = 'auto', radius_style = 'csd', radius_scale_factor = 2.00, write_layer_info = False, suppress_warning = True):
    '''
    Get layer parameters from POSCAR
    delta: This is the spacial resolution (in unit of Angstrom) for the atoms, for example z=0.24 and z=0.25 are considered to be in the same position if delta=0.01. Default: delta = 0.05
    '''
    args_dict = locals()
    import numpy as np
    from .. import funcs
    from . import vasp_read
    from .. import periodic_table
    from copy import copy, deepcopy
    import math
    import os
    from . import vasp_build
    import itertools

    periodic_table_dict = periodic_table.periodic_tab()

    original_criteria = criteria
    original_layer_dist_tolerance = layer_dist_tolerance
    #recommended_layer_dist_tolerance = 0.664 
    #######################################################
    # Creating layer property dictionary: poscar_layer_dict
    #######################################################
    poscar_layer_dict = {}
    poscar_layer_dict['delta'] = delta

    ##import pprint
    ##pprint.pprint(poscar_direction_dict)
    poscar_dict['atom_number_ortho_' + poscar_direction_dict['number_order_text']] = funcs.atom_number_ortho(
        atom_key_arr = poscar_dict['atom_key_arr'],
        pos_arr_cartesian= poscar_dict['pos_arr'][:,3:6],
        delta = delta,
        number_order = poscar_direction_dict['number_order_text']
        )
    layer_dist_tolerance = poscar_layer_dist_tolerance(poscar_dict, radius_style = radius_style, radius_scale_factor = radius_scale_factor)
    # group layers according to different properties and adding tags
    poscar_layer_dict['direction'] = poscar_direction_dict['direction']
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
    poscar_layer_dict['layer_info_str'] = None

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

    def poscar_layer_status(poscar_dict, poscar_direction_dict, delta = 0.05, radius_style = 'csd', radius_scale_factor = 2.00, write_layer_info = False, suppress_warning = True):
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
        poscar_layer_dict['direction'] = poscar_direction_dict['direction']
        poscar_layer_dict['delta'] = delta
        poscar_layer_dict['num_layers'] = 0 
        poscar_layer_dict['layer_coord_list'] = [list()]
        poscar_layer_dict['indx_list'] = [None]
        poscar_layer_dict['atoms_indx_list'] =  [list()]
        poscar_layer_dict['atom_name_list'] =  [list()]
        poscar_layer_dict['bond_status_list'] = [None]
        poscar_layer_dict['layer_bonding_indx_list'] = [None]
        poscar_layer_dict['layer_bonding_atom_name_list'] = [None]
        # periodicity_status_list: Each layer is labled with the periodicity of the layer
        poscar_layer_dict['periodicity_status_list'] = [None]
        # cell_by_bond_list: Each layer is labeled with the bond status of the cell. If all the items in the list are zeros, then this means the structure cannot be analyzed by bonding crietrion and atom layers are all bonded.
        # cell_by_periodicity_list: Each layer is labeled with the periodicity status of the cell
        poscar_layer_dict['cell_by_bond_list'] = [None]
        poscar_layer_dict['cell_by_periodicity_list'] = [None]
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
        s1 = poscar_dict['l_arr'][poscar_direction_dict['l_arr_column_1'],:]
        s2 = poscar_dict['l_arr'][poscar_direction_dict['l_arr_column_2'],:]
        shift_1 = s1
        shift_2 = -s1
        shift_3 = s2
        shift_4 = -s2
        shift_5 = s1 + s2
        shift_6 = s1 - s2
        shift_7 = -s1 + s2
        shift_8 = -s1 - s2
        # First, find the atoms in each layer according to the parameter 'delta'
        poscar_layer_dict['cell_by_bond_list'][0] = 1
        poscar_layer_dict['cell_by_periodicity_list'][0] = 1
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
                poscar_layer_dict['layer_coord_list'].append(list())
                poscar_layer_dict['atoms_indx_list'].append(list())
                poscar_layer_dict['atom_name_list'].append(list())
                poscar_layer_dict['bond_status_list'].append(None)
                poscar_layer_dict['layer_bonding_indx_list'].append(None)
                poscar_layer_dict['layer_bonding_atom_name_list'].append(None)
                poscar_layer_dict['periodicity_status_list'].append(None)
                poscar_layer_dict['cell_by_bond_list'].append(1)
                poscar_layer_dict['cell_by_periodicity_list'].append(1)
            poscar_layer_dict['indx_list'][i_layer_indx] = i_layer_indx
            poscar_layer_dict['layer_coord_list'][i_layer_indx].append(i_atom_pos)
            if shift_top2bottom == False:
                ##poscar_layer_dict['layer_coord_list'][i_layer_indx] = i_atom_pos
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
                        # Due to periodic boundary condition, some atoms may form bond across the boundary. The j atom is shifted to account for this 
                        bonding_0 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos, mode = radius_style, factor = radius_scale_factor) 
                        bonding_1 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos + shift_1, mode = radius_style, factor = radius_scale_factor) 
                        bonding_2 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos + shift_2, mode = radius_style, factor = radius_scale_factor) 
                        bonding_3 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos + shift_3, mode = radius_style, factor = radius_scale_factor) 
                        bonding_4 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos + shift_4, mode = radius_style, factor = radius_scale_factor) 
                        bonding_5 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos + shift_5, mode = radius_style, factor = radius_scale_factor) 
                        bonding_6 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos + shift_6, mode = radius_style, factor = radius_scale_factor) 
                        bonding_7 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos + shift_7, mode = radius_style, factor = radius_scale_factor) 
                        bonding_8 = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos + shift_8, mode = radius_style, factor = radius_scale_factor) 
                        ##bonding = funcs.form_bond(i_elmt, j_elmt, i_atom_pos, j_atom_pos, mode = radius_style, factor = radius_scale_factor) 
                        bonding = any([bonding_0, bonding_1, bonding_2, bonding_3, bonding_4, bonding_5 ,bonding_6, bonding_7 ,bonding_8])
                        ij_atom_bonding_list.append(bonding)
                        if bonding == True:
                            ij_layer_bonding = True
                            poscar_layer_dict['layer_bonding_atom_name_list'][i_layer_indx].append(poscar_dict['atomname_list'][i_atom_indx] + '-' + poscar_dict['atomname_list'][j_atom_indx])  
                    for k_atom_indx in k_layer_atom_indx_list:
                        k_elmt = poscar_dict['atom_species_arr'][k_atom_indx]
                        k_atom_pos = poscar_dict['pos_arr'][k_atom_indx, 3:6] + base_k * shift_vec
                        ##bonding = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos, mode = radius_style, factor = radius_scale_factor) 
                        bonding_0 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos, mode = radius_style, factor = radius_scale_factor) 
                        bonding_1 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos + shift_1, mode = radius_style, factor = radius_scale_factor) 
                        bonding_2 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos + shift_2, mode = radius_style, factor = radius_scale_factor) 
                        bonding_3 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos + shift_3, mode = radius_style, factor = radius_scale_factor) 
                        bonding_4 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos + shift_4, mode = radius_style, factor = radius_scale_factor) 
                        bonding_5 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos + shift_5, mode = radius_style, factor = radius_scale_factor) 
                        bonding_6 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos + shift_6, mode = radius_style, factor = radius_scale_factor) 
                        bonding_7 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos + shift_7, mode = radius_style, factor = radius_scale_factor) 
                        bonding_8 = funcs.form_bond(i_elmt, k_elmt, i_atom_pos, k_atom_pos + shift_8, mode = radius_style, factor = radius_scale_factor) 
                        bonding = any([bonding_0, bonding_1, bonding_2, bonding_3, bonding_4, bonding_5 ,bonding_6, bonding_7 ,bonding_8])
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
        ##print(poscar_layer_dict['bond_status_list'])
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
            # the thresold_vacuum_width, the vacuum exists. And the i-th layer and the k-th layer of 
            # the bond_status_list take the values of 1 and -1
            ########################################################################################
            thresold_vacuum_width = 5.0
            if len(set(poscar_layer_dict['bond_status_list'])) == 1 and len(poscar_layer_dict['bond_status_list']) != 1 and 2 in poscar_layer_dict['bond_status_list']:
                top_layer_indx = poscar_layer_dict['num_layers'] - 1
                bottom_layer_indx = 0
                ##top_width = poscar_direction_dict['box_len_ortho'] - poscar_layer_dict['layer_coord_list'][top_layer_indx]
                ##bottom_width = poscar_layer_dict['layer_coord_list'][bottom_layer_indx]
                top_width = poscar_direction_dict['box_len_ortho'] - np.max([item for sublist in poscar_layer_dict['layer_coord_list'] for item in sublist])
                bottom_width = np.min([item for sublist in poscar_layer_dict['layer_coord_list'] for item in sublist])
                top_plus_bottom_width = top_width + bottom_width
                ##print('box_len_ortho=',poscar_direction_dict['box_len_ortho'],poscar_layer_dict['layer_coord_list'])
                ##print('top=',top_width,'bottom=',bottom_width,'top+bottom=',top_plus_bottom_width)
                if top_plus_bottom_width > thresold_vacuum_width:
                    poscar_layer_dict['bond_status_list'][top_layer_indx] = 1
                    poscar_layer_dict['bond_status_list'][bottom_layer_indx] = -1
                # Get layer spearation between the i-th and the k-th layers (the k-th layer is the next layer to the i-th layer)
                ik_layer_separation_list = []
                n_indx = 1
                for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
                    k_layer_indx_temp = i_layer_indx + n_indx
                    base_k = math.floor(k_layer_indx_temp / poscar_layer_dict['num_layers'])
                    k_layer_indx = k_layer_indx_temp - base_k * poscar_layer_dict['num_layers']
                    # find the minimum separation between the i-th and the k-th layers
                    if i_layer_indx == 0 or k_layer_indx == 0: 
                        temp_separation_list = []
                        i_layer_pos_min = np.min(poscar_layer_dict['layer_coord_list'][i_layer_indx])
                        i_layer_pos_max = np.max(poscar_layer_dict['layer_coord_list'][i_layer_indx])
                        k_layer_pos_min = np.min(poscar_layer_dict['layer_coord_list'][k_layer_indx])
                        k_layer_pos_max = np.max(poscar_layer_dict['layer_coord_list'][k_layer_indx])
                        for i_pos in [i_layer_pos_min, i_layer_pos_max]:
                            for k_pos in [k_layer_pos_min, k_layer_pos_max]:
                                temp_separation_list.append(abs(i_pos - k_pos))
                        ##ik_layer_separation = abs(poscar_layer_dict['layer_coord_list'][k_layer_indx] - poscar_layer_dict['layer_coord_list'][i_layer_indx])
                        dis_1 = abs((i_pos - k_pos) + poscar_direction_dict['box_len_ortho'])
                        dis_2 = abs((i_pos - k_pos) - poscar_direction_dict['box_len_ortho'])
                        dis_3 = abs(-(i_pos - k_pos) + poscar_direction_dict['box_len_ortho'])
                        dis_4 = abs(-(i_pos - k_pos) - poscar_direction_dict['box_len_ortho'])
                        ik_layer_separation = np.min([np.min(temp_separation_list), dis_1, dis_2, dis_3, dis_4])
                    else:
                        ik_layer_separation = abs(np.min(poscar_layer_dict['layer_coord_list'][k_layer_indx]) - np.min(poscar_layer_dict['layer_coord_list'][i_layer_indx]))
                    ik_layer_separation_list.append(ik_layer_separation)

                    #print('ik_layer_separation=',ik_layer_separation, 'k=', k_layer_indx,poscar_layer_dict['layer_coord_list'][k_layer_indx], 'i=',i_layer_indx, poscar_layer_dict['layer_coord_list'][i_layer_indx])
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
        ##if write_layer_info == True:
        ##    print('bond_status_list = ',poscar_layer_dict['bond_status_list'])
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
                    break
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
                    break

        #########################################
        # Get periodicity status of layers
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
                if not poscar_layer_dict['periodicity_status_list'][j_layer_indx] is None:
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
                        ij_atom_cartesian_vec[0:3] = np.dot(ij_atom_direct_vec, poscar_dict['l_arr'])
                        ij_atom_cartesian_dist = np.linalg.norm(ij_atom_cartesian_vec) 

                        if i_elmt == j_elmt and ij_atom_cartesian_dist <= atom_dist_tolerance:
                            identical_atoms_list2.append(True)
                        #elif i_elmt != j_elmt and ij_atom_cartesian_dist > atom_dist_tolerance:
                        elif i_elmt == j_elmt and ij_atom_cartesian_dist > atom_dist_tolerance:
                            identical_atoms_list2.append(False)
                        else:
                            identical_atoms_list2.append(False)
                    if True in identical_atoms_list2:
                        identical_atoms_list1.append(True)
                    else:
                        identical_atoms_list1.append(False)
                if False in identical_atoms_list1:                 
                    pass
                else:
                    poscar_layer_dict['periodicity_status_list'][j_layer_indx] = original_periodicity_number 
            #print('layer ', i_layer_indx + 1, identical_atoms_list2)
            #print(i_layer_indx, poscar_layer_dict['periodicity_status_list'][i_layer_indx])

        ################################################################################################################
        # Get cell periodicity based on bond_status_list and periodicity_status_list, each layer is labeled with the cell periodicity number
        ################################################################################################################
        cell_by_periodicity_number = 0
        # Each layer is labeled by a number , so find out the elements with the value of 0 or 1 in the periodicity status list
        layer_indices_0_list = [i for i, x in enumerate(poscar_layer_dict['periodicity_status_list']) if x == 0]
        layer_indices_1_list = [i for i, x in enumerate(poscar_layer_dict['periodicity_status_list']) if x == 1]
        # characteristic_periodicity_indx denotes the characteristic index for identifying the new layers. When checking the layer index from layer to layer, if the characteristic_periodicithy_indx is met, a new periodicity is identified. 
        # Due to the periodic boundary condition, in some model, the first layer (layer index =0) may be the same layer in the other end of the model, The program starts to screen atoms from the first layer, but in this case the first layer only contains part of the atoms of the actual layer,  thus in this case, in searching for the new layers, better set the second layer (layer index = 1) as the starting layer.
        characteristic_periodicity_indx = 0
        if len(layer_indices_0_list) == 1:
            # this situation corresponds to a model with only one layer
            characteristic_periodicity_indx = 0
        if len(layer_indices_0_list) > 1 and len(layer_indices_1_list) > 1:
            delta_indices_0 = layer_indices_0_list[1] - layer_indices_0_list[0]
            delta_indices_1 = layer_indices_1_list[1] - layer_indices_1_list[0]
            ##if delta_indices_0 > delta_indices_1:
            ##    characteristic_periodicity_indx = 1
            ##else:
            ##    characteristic_periodicity_indx = 0
            characteristic_periodicity_indx = 0
        # get the layer index of the layer with the maximum periodicity status index
        # Only when the maximum periodicity index is exceeded, and when the characteristic_periodicity_indx is met, the new period is admitted.
        max_periodicity_indx = np.max([x for x in poscar_layer_dict['periodicity_status_list'] if x is not None])
        max_p_layer_indx_list = np.argwhere(poscar_layer_dict['periodicity_status_list'] == max_periodicity_indx)
        max_p_layer_indx_list = list(itertools.chain.from_iterable(max_p_layer_indx_list))

        cell_by_bond_number = 0
        cell_by_periodicity_number = 0
        admit_new_period = True
        for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
           # Deal with bond_status_list
           if poscar_layer_dict['bond_status_list'][i_layer_indx] == -1:
               cell_by_bond_number = cell_by_bond_number + 1
           poscar_layer_dict['cell_by_bond_list'][i_layer_indx] = cell_by_bond_number

           # Deal with periodicity_status_list
           if i_layer_indx in max_p_layer_indx_list:
               # once the maximum periodicity index is found, the new period will be admitted when the characteristic_periodicity_indx is met.
               admit_new_period = True
           if (poscar_layer_dict['periodicity_status_list'][i_layer_indx] == characteristic_periodicity_indx) and admit_new_period == True:
               cell_by_periodicity_number = cell_by_periodicity_number + 1
               admit_new_period = False
           poscar_layer_dict['cell_by_periodicity_list'][i_layer_indx] = cell_by_periodicity_number
        # If the layers can be distinguished by bonding criterion, then modify the cell period number with value 0 to a resonable value
        max_b_layer_number = np.max(poscar_layer_dict['cell_by_bond_list'])
        for i_layer_indx in range(0, poscar_layer_dict['num_layers']):
            if poscar_layer_dict['cell_by_bond_list'][i_layer_indx] == 0:
                poscar_layer_dict['cell_by_bond_list'][i_layer_indx] = max_b_layer_number

        def get_layer_info_str(poscar_dict, poscar_layer_dict):

            '''
            Get a string of layer information
            '''
            layer_info_str = ''
            layer_info_str = layer_info_str + '# ' + poscar_dict['file_path'] + '\n'
            #if poscar_direction_dict['number_order_text'] == 'z':
            if poscar_layer_dict['direction'] == 'z':
                layer_info_str = layer_info_str + 'direction: c' + '\n'
            #elif poscar_direction_dict['number_order_text'] == 'y':
            elif poscar_layer_dict['direction'] == 'y':
                layer_info_str = layer_info_str + 'direction: b' + '\n'
            #elif poscar_direction_dict['number_order_text'] == 'x':
            elif poscar_layer_dict['direction'] == 'x':
                layer_info_str = layer_info_str + 'direction: a' + '\n'
            layer_info_str = layer_info_str + ( 
                    funcs.str_format('layer', 7, ' ') + 
                    funcs.str_format('bonding', 9, ' ') + 
                    funcs.str_format('periodicity', 13, ' ') + 
                    funcs.str_format('cell_bond', 11, ' ') + 
                    funcs.str_format('cell_period', 13, ' ') + 
                    funcs.str_format('atoms in layer', 25, ' ') + 
                    'bonded atoms list' + '\n')
            for i_layer_indx in reversed(range(poscar_layer_dict['num_layers'])): 
                i_layer_indx_str = str(i_layer_indx + 1)
                i_bond_status_str = str(poscar_layer_dict['bond_status_list'][i_layer_indx])
                i_periodicity_status_str = str(poscar_layer_dict['periodicity_status_list'][i_layer_indx])
                i_cell_bond_str = str(poscar_layer_dict['cell_by_bond_list'][i_layer_indx])
                i_cell_period_str = str(poscar_layer_dict['cell_by_periodicity_list'][i_layer_indx])
                i_atom_name_list_str = str(poscar_layer_dict['atom_name_list'][i_layer_indx])
                i_bonded_atom_name_list_str = str(poscar_layer_dict['layer_bonding_atom_name_list'][i_layer_indx])
                layer_info_str = layer_info_str + ( 
                    funcs.str_format(i_layer_indx_str, 7, ' ') + 
                    funcs.str_format(i_bond_status_str, 9, ' ') + 
                    funcs.str_format(i_periodicity_status_str, 13, ' ') + 
                    funcs.str_format(i_cell_bond_str, 11, ' ') + 
                    funcs.str_format(i_cell_period_str, 13, ' ') + 
                    funcs.str_format(i_atom_name_list_str, 25, ' ') +
                    i_bonded_atom_name_list_str + '\n')
            return layer_info_str
        poscar_layer_dict['layer_info_str'] = get_layer_info_str(poscar_dict, poscar_layer_dict)
        ##################################
        # print atoms in each layer
        ##################################
        if write_layer_info == True:
            print(poscar_layer_dict['layer_info_str'])
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
    condition_a1 = (poscar_layer_dict['bond_status_list'].count(1) > 1) and (poscar_layer_dict['bond_status_list'].count(-1) > 1)

    num_iter = 10
    radius_scale_factor_list = radius_scale_factor + 0.05 * np.array([sum(np.arange(1, x)) for x in np.arange(2, num_iter + 2)])

    if (condition_1 and condition_2) or condition_a1:
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
    return poscar_layer_dict

def kpoint_rec2car(vec_rec, reciprocal_arr):
    '''
    vec_rec: vector of kpoint in reciprocal coordinate
    vec_car: vector of kpoint in Cartesian coordinate
    reciprocal_arr: an array (3x3) of the reciprocal lattice
    '''
    args_dict = locals()
    import numpy as np
    b_1 = reciprocal_arr[0,:]
    b_2 = reciprocal_arr[1,:]
    b_3 = reciprocal_arr[2,:]
    vec_car = np.array([None]*3)    
    vec_car[0] = (vec_rec[0] * b_1 + vec_rec[1] * b_2  + vec_rec[2] * b_3)[0]
    vec_car[1] = (vec_rec[0] * b_1 + vec_rec[1] * b_2  + vec_rec[2] * b_3)[1]
    vec_car[2] = (vec_rec[0] * b_1 + vec_rec[1] * b_2  + vec_rec[2] * b_3)[2]
    return vec_car

def get_vacuum_status(poscar_file_path):
    '''Check whether vacuum exists in x/y/z directions'''
    from . import vasp_read
    vacuum_in_xyz = [False, False, False]
    for direction in ['x','y','z']:
        poscar_dict = vasp_read.read_poscar(poscar_file_path)
        poscar_direction_dict = poscar_direction_params(poscar_dict, direction = direction)
        delta = 0.05
        radius_scale_factor = 2.0 
        write_layer_info = False
        suppress_warning = True
        poscar_layer_dict = poscar_layer_params(poscar_dict, poscar_direction_dict, criteria = 'bonding', delta = delta, layer_dist_tolerance = 'auto', radius_style = 'csd', radius_scale_factor = radius_scale_factor, write_layer_info = write_layer_info, suppress_warning = suppress_warning)
