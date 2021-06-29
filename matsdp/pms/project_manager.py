# -*- coding: utf-8 -*-
'''
This is the Project Management Tool (PMT)
'''

def gen_submit_script(project_dir, submit_script_path_dict = None, queue_system = 'PBS'):
    '''
    generate submit script for the specified project
    submit_script_path_dict: the keyword 'universal' must exist in submit_script_path_dict, if no submit files are provided for other calculations, the submit file in submit_script_path_dict['universal'] will be used. 
    '''
    args_dict = locals()
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
    project_dir = os.path.abspath(project_dir)
    task_name_list = [x for x in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, x))]
    submit_script_name = 'submit.sh'
    run_submit_script_name = 'run_submit.sh'
    for i_task_indx in range(0, len(task_name_list)):
        i_task_dir = os.path.join(project_dir, task_name_list[i_task_indx])
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

def queue_management(project_dir, max_num_running_jobs = 2, 
                     submit_script_name = 'submit.sh', queue_system = 'PBS', job_kwd = 'vaspjob',
                     task_flow_list = ['test_encut', 'test_kpoints', 'opt', 'scf', 'bs'],
                     submitted_job_dir_list = [], 
                     elmt_potcar_dir = None,
                     apply_to_existed_task = True,
                     submit_script_path_dict = None,
                     num_atoms_testing = 4,
                     use_primitive_cell = True, 
                     kpath_scheme = 'vaspkit', 
                     kpath_params = '''echo 3 ; echo 303 | vaspkit''',
                     debug_mode = False,
                    ):
    '''
    queuing system management tool
    debug_mode: under the debug mode, the jobs won't be submitted, but the programmer can see the job submission command.
    '''
    args_dict = locals()
    import os
    import time
    from .. import funcs
    from .. import default_params
    from ..vasp import vasp_analyze
    from ..vasp import vasp_tools
    from ..vasp import vasp_read
    from . import task_manager

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    projects_dir = os.path.join(output_dir, defaults_dict['projects_dir_name'])

    if elmt_potcar_dir not in [None, 'None', 'none']:
        elmt_potcar_dir = os.path.abspath(elmt_potcar_dir)

    project_dir = os.path.abspath(project_dir)
    project_name = os.path.split(project_dir)[-1]
    task_name_list = [x for x in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, x))]

    if submitted_job_dir_list in [None, 'None', 'none'] or len(submitted_job_dir_list) == 0:
        submitted_job_dir_list = []
    else:
        for i_indx in range(0, len(submitted_job_dir_list)):
            submitted_job_dir_list[i_indx] = os.path.abspath(submitted_job_dir_list[i_indx])

    log_str = ''
    control_sh_file_path = os.path.join(output_dir, 'control_submit.sh')
    funcs.touch(control_sh_file_path)
    run_submit_script_name = 'run_submit.sh'
    job_id_file_name = 'job_id.txt'
    num_submitted_jobs = 0
    for i_task_indx in range(0, len(task_name_list)):
        start_next_task = False
        ##print(task_name_list[i_task_indx])
        i_task_dir = os.path.join(project_dir, task_name_list[i_task_indx])
        opt_dir = os.path.join(i_task_dir, 'opt/')
        scf_dir = os.path.join(i_task_dir, 'scf/') 
        dos_dir = os.path.join(i_task_dir, 'dos/')
        bs_dir = os.path.join(i_task_dir, 'bs/')
        bs_soc_dir = os.path.join(i_task_dir, 'bs_soc/')
        nscf_dir = os.path.join(i_task_dir, 'nscf/')
        test_encut_dir = os.path.join(i_task_dir, 'test_encut/')
        test_kpoints_dir = os.path.join(i_task_dir, 'test_kpoints/')
        for i_job_indx in range(0, len(task_flow_list)):
            i_job_name = task_flow_list[i_job_indx] 
            i_job_dir = os.path.join(i_task_dir, task_flow_list[i_job_indx])
            ##print('    ' + task_flow_list[i_job_indx])
            if 'test' in task_flow_list[i_job_indx]:
                if vasp_read.read_poscar(os.path.join(opt_dir, 'POSCAR'))['n_atoms'] > num_atoms_testing:
                    # for models with larger than num_atoms_testing atoms, we don't test ENCUT and KPOINTS.
                    continue
                # for subjods
                subjob_name_list = [x for x in os.listdir(i_job_dir) if os.path.isdir(os.path.join(i_job_dir, x))]
                for i_subjob_indx in range(0, len(subjob_name_list)):
                    i_subjob_dir = os.path.join(i_job_dir, subjob_name_list[i_subjob_indx])
                    i_job_finished = vasp_tools.vasp_job_finished(os.path.join(i_subjob_dir, 'OUTCAR'), suppress_warning = True)
                    i_job_type = task_manager.check_job_type(i_subjob_dir)
                    ##if i_job_type != 'VASP':
                    ##    continue 
                    if i_job_finished == True:
                        continue
                    else:
                        num_running_jobs = get_num_running_jobs(queue_system, job_kwd)
                        if num_running_jobs < max_num_running_jobs:
                            # Check whether the job has been interrupted (unfinihsed OUTCAR exist and the job is not submitted), if interrupted, then the job will not be submitted again and the user must check the job and submit the job manually.
                            i_outcar_file_path = os.path.join(i_subjob_dir, 'OUTCAR')
                            if os.path.exists(i_outcar_file_path) and os.path.isfile(i_outcar_file_path):
                            #if task_manager.job_submission_status(job_dir = i_subjob_dir, job_kwd = job_kwd, queue_system = queue_system) == False:
                                print('WARNING #2102051224 (from project_manager): The job ' + i_subjob_dir + ' is unfinished or interrupted or running, pleased check it manually.')
                                continue
                            #check whether the job exists in the queuing list
                            submit_script_file_path = os.path.join(i_subjob_dir, submit_script_name)
                            job_id_file_path = os.path.join(i_subjob_dir, job_id_file_name)
                            if os.path.exists(job_id_file_path):
                                with open(job_id_file_path, 'r') as f:
                                    line = f.readlines()
                                    try:
                                        job_id = int(line[0])
                                    except:
                                        job_id = None
                            else:
                                job_id = None
                            if task_manager.job_id_exists(job_id) == True:
                                # if the job is running, then pass
                                pass
                            else:
                                if i_subjob_dir in submitted_job_dir_list:
                                    print('WARNING #20110601 (from project_manager.queue_management): Repetitive submission! Please check your inputs. The directories to be checked is: ' + str(i_subjob_dir))
                                    ##exit()
                                    continue
                                elif vasp_tools.vasp_job_finished(os.path.join(i_subjob_dir, 'OUTCAR'), suppress_warning = True) == True:
                                    pass
                                else:
                                    if (i_job_name in ['test_kpoints'] or i_job_name in ['test_encut']) and vasp_read.read_poscar(os.path.join(opt_dir, 'POSCAR'))['n_atoms'] > num_atoms_testing:
                                        # if the number of atom is larger than num_atoms_testing, then encut and kpionts tests will not be performed.
                                        continue
                                    elif i_job_name in ['test_kpoints'] and test_job_finished(test_encut_dir) == False:
                                        continue
                                    else:
                                        i_job_poscar_file_path_opt = os.path.join(opt_dir, 'POSCAR') 
                                        task_manager.gen_inputs(poscar_file_path_list = [i_job_poscar_file_path_opt], project_name = project_name, task_type = 'VASP', elmt_potcar_dir = elmt_potcar_dir, apply_to_existed_task = apply_to_existed_task, use_primitive_cell = use_primitive_cell, kpath_scheme = kpath_scheme, kpath_params = kpath_params)
                                        task_manager.gen_submit_script(i_task_dir, submit_script_path_dict = submit_script_path_dict, task_flow_list = task_flow_list, queue_system = queue_system)
                                        # check errors
                                        error_type = vasp_tools.check_error(i_subjob_dir)
                                        # handle errors
                                        vasp_tools.error_handler(i_subjob_dir)
                                        if error_type == 'BAD TERMINATION CODE: 9':
                                            continue
                                        # submit job
                                        pwd_dir = os.getcwd()
                                        fpath, fname = os.path.split(submit_script_file_path)
                                        with open(control_sh_file_path, 'w') as f:
                                            f.write('#!/bin/bash\n')
                                            f.write('cd ' + fpath + '\n')
                                            f.write(r'sed -i "s/\r//" ' + submit_script_file_path + '\n')
                                            # redirect the job id to the file job_id_file_path
                                            submit_cmd = 'qsub ' + submit_script_file_path + ' > ' + job_id_file_path + '\n'
                                            f.write(submit_cmd)
                                            current_time = time.time()
                                            formatted_time = time.strftime('## %Y%m%d %H:%M:%S',time.localtime(current_time))
                                            print(formatted_time)
                                            print(submit_cmd)
                                            f.write('cd ' + pwd_dir + '\n')
                                        if debug_mode == False:
                                            os.system('sh ' + control_sh_file_path)
                                        log_str = log_str + '# ' + submit_cmd + '\n'
                                        funcs.write_log(logfile, log_str)
                                        submitted_job_dir_list.append(i_subjob_dir)
                                        num_submitted_jobs = num_submitted_jobs + 1
            else:
                i_job_finished = vasp_tools.vasp_job_finished(os.path.join(i_job_dir, 'OUTCAR'), suppress_warning = True)
                i_job_type = task_manager.check_job_type(i_job_dir)
                #if i_job_type != 'VASP':
                #    ##print('        ' + 'not vasp job')
                #    continue 
                if i_job_finished == True:
                    ##print('        ' + 'finished')
                    continue
                else:
                    num_running_jobs = get_num_running_jobs(queue_system, job_kwd) 
                    if num_running_jobs < max_num_running_jobs:
                        # Check whether the job has been interrupted (unfinihsed OUTCAR exist and the job is not submitted), if interrupted, then the job will not be submitted again and the user must check the job and submit the job manually.
                        i_outcar_file_path = os.path.join(i_job_dir, 'OUTCAR')
                        if os.path.exists(i_outcar_file_path) and os.path.isfile(i_outcar_file_path):
                        ##if task_manager.job_submission_status(job_dir = i_job_dir, job_kwd = job_kwd, queue_system = queue_system) == True:
                            print('WARNING #2102051227 (from project_manager): The job ' + i_job_dir + ' is unfinished or interrupted or running, pleased check it manually.')
                            continue
                        # submit job
                        submit_script_file_path = os.path.join(i_job_dir, submit_script_name)
                        run_submit_script_file_path = os.path.join(i_job_dir, run_submit_script_name)
                        job_id_file_path = os.path.join(i_job_dir, job_id_file_name)
                        #check whether the job exists in the queuing list
                        if os.path.exists(job_id_file_path):
                            with open(job_id_file_path, 'r') as f:
                                line = f.readlines()
                                try:
                                    job_id = int(line[0])
                                except:
                                    job_id = None
                        else:
                            job_id = None
                        if task_manager.job_id_exists(job_id) == True:
                            pass
                        else:
                            if i_job_dir in submitted_job_dir_list:
                                print('WARNING #20110603 (from project_manager.queue_management): Repetitive submission! Please check your inputs. The directories to be checked is: ' + str(i_job_dir))
                                ##exit()
                                continue
                            elif vasp_tools.vasp_job_finished(os.path.join(i_job_dir, 'OUTCAR'), suppress_warning = True) == True:
                                pass
                            else:
                                if i_job_name in ['scf'] and vasp_tools.vasp_job_finished(os.path.join(opt_dir, 'OUTCAR'), suppress_warning = True) == False:
                                    continue
                                elif i_job_name in ['bs', 'bs_soc', 'nscf', 'dos'] and vasp_tools.vasp_job_finished(os.path.join(scf_dir, 'OUTCAR'), suppress_warning = True) == False:
                                    continue
                                elif i_job_name in ['opt'] and vasp_read.read_poscar(os.path.join(opt_dir, 'POSCAR'))['n_atoms'] <= num_atoms_testing and  test_job_finished(test_kpoints_dir) == False:
                                    continue  
                                else:
                                    i_job_poscar_file_path_opt = os.path.join(opt_dir, 'POSCAR')
                                    task_manager.gen_inputs(poscar_file_path_list = [i_job_poscar_file_path_opt], project_name = project_name, task_type = 'VASP', elmt_potcar_dir = elmt_potcar_dir, apply_to_existed_task = apply_to_existed_task, use_primitive_cell = use_primitive_cell, kpath_scheme = kpath_scheme, kpath_params = kpath_params)
                                    task_manager.gen_submit_script(i_task_dir, submit_script_path_dict = submit_script_path_dict, task_flow_list = task_flow_list, queue_system = queue_system)
                                    # check errors
                                    error_type = vasp_tools.check_error(i_job_dir)
                                    # handle errors
                                    vasp_tools.error_handler(i_job_dir)
                                    if error_type == 'BAD TERMINATION CODE: 9':
                                        continue
                                    # submit job
                                    pwd_dir = os.getcwd()
                                    fpath, fname = os.path.split(submit_script_file_path)
                                    with open(control_sh_file_path, 'w') as f:
                                        f.write('#!/bin/bash\n')
                                        f.write('cd ' + fpath + '\n')
                                        f.write(r'sed -i "s/\r//" ' + submit_script_file_path + '\n')
                                        submit_cmd = 'qsub ' + submit_script_file_path + ' > ' + job_id_file_path + '\n'
                                        f.write(submit_cmd)
                                        current_time = time.time()
                                        formatted_time = time.strftime('## %Y%m%d %H:%M:%S',time.localtime(current_time))
                                        print(formatted_time)
                                        print(submit_cmd)
                                        f.write('cd ' + pwd_dir + '\n')
                                    if debug_mode == False:
                                        os.system('sh ' + control_sh_file_path)
                                    log_str = log_str + '# ' + submit_cmd + '\n'
                                    funcs.write_log(logfile, log_str)
                                    submitted_job_dir_list.append(i_job_dir)
                                    num_submitted_jobs = num_submitted_jobs + 1
    return submitted_job_dir_list

# get the number of running jobs according to the jobs keyword (job_kwd)
def get_num_running_jobs(queue_system, job_kwd):
    '''
    get number of running jobs
    Note that the queuing jobs is also counted
    '''
    args_dict = locals()
    import os
    from .. import funcs
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])

    queue_log_file_name = 'matsdp_pms.log'
    queue_log_file_path = os.path.join(output_dir, queue_log_file_name)
    funcs.touch(queue_log_file_path)
    if queue_system == 'PBS':
        #os.system('qstat -r ' + ' | grep ' + job_kwd + ' > ' + queue_log_file_path)   #running jobs
        os.system('(qstat -r ; qstat -i)' + ' | grep ' + job_kwd + ' > ' + queue_log_file_path)   #running jobs and queuing jobs
        with open(queue_log_file_path, 'r') as f:
            line = f.readlines() 
            num_lines = len(line)
    num_running_jobs = num_lines
    return num_running_jobs

def test_job_finished(test_job_dir):
    '''
    whether the test job (test_encut or test_kpoints) has been finishied or not
    '''
    args_dict = locals()
    import os
    from ..vasp import vasp_analyze
    job_finished = True
    test_job_dir = os.path.abspath(test_job_dir)
    test_job_list = [x for x in os.listdir(test_job_dir) if os.path.isdir(os.path.join(test_job_dir, x))]
    for i_test_subjob_indx in range(0, len(test_job_list)):
        i_test_subjob_dir = os.path.join(test_job_dir, test_job_list[i_test_subjob_indx])
        i_outcar_file_path = os.path.join(i_test_subjob_dir, 'OUTCAR')
        if os.path.exists(i_outcar_file_path) and os.path.isfile(i_outcar_file_path): 
            if vasp_tools.vasp_job_finished(i_outcar_file_path, suppress_warning = True) == False:
                job_finished = False
        else:
            job_finished = False
    return job_finished


def job_monitor(project_dir, max_num_running_jobs = 3, 
                submit_script_name = 'submit.sh', queue_system = 'PBS', job_kwd = 'vaspjob',
                task_flow_list = ['test_encut', 'test_kpoints', 'opt', 'scf', 'bs'],
                submitted_job_dir_list = [],
                elmt_potcar_dir = None,
                apply_to_existed_task = True,
                job_check_frequency = None, max_monitor_time = None,
                submit_script_path_dict = None,
                num_atoms_testing = 4, 
                use_primitive_cell = True, 
                kpath_scheme = 'vaspkit', 
                kpath_params = '''echo 3 ; echo 303 | vaspkit''',
                debug_mode = False,
               ):
    '''
    execute code repeately
    units of job_check_frequency and max_monitor_time are all in min.
    ''' 
    args_dict = locals()
    import os
    import time
    from .. import convert
    from ..vasp import vasp_tools
    ##import threading

    max_monitor_time_temp = max_monitor_time
    job_check_frequency = convert.time_converter(hour = 0, minute = job_check_frequency, second = 0, unit = 'seconds')
    max_monitor_time = convert.time_converter(hour = 0, minute = max_monitor_time, second = 0, unit = 'seconds')

    start = time.time()
    ##formatted_time = time.strftime('## %Y%m%d %H:%M:%S',time.localtime(current_time))
    def monitor_timer(submitted_job_dir_list, last_submitted_job_dir_list):
        import random
        import time
        ##print("Please execute your code here")
        last_submitted_job_dir_list = submitted_job_dir_list
        submitted_job_dir_list = queue_management(
            project_dir = project_dir,
            max_num_running_jobs = max_num_running_jobs,
            submit_script_name = submit_script_name, 
            queue_system = queue_system, 
            job_kwd = job_kwd,
            task_flow_list = task_flow_list,
            submitted_job_dir_list = submitted_job_dir_list,
            elmt_potcar_dir = elmt_potcar_dir,
            apply_to_existed_task = apply_to_existed_task,
            submit_script_path_dict = submit_script_path_dict,
            num_atoms_testing = num_atoms_testing, 
            use_primitive_cell = use_primitive_cell, 
            kpath_scheme = kpath_scheme, 
            kpath_params = kpath_params,
            debug_mode = debug_mode,
            )
        ##print('last_submitted: ',last_submitted_job_dir_list)
        ##print('newly_submitted: ',submitted_job_dir_list) 
        end = time.time()
        duration = end - start
        return duration, submitted_job_dir_list, last_submitted_job_dir_list

    duration = 0 
    num_loops = 10000
    last_submitted_job_dir_list = []
    if submitted_job_dir_list in [None, 'None', 'none'] or len(submitted_job_dir_list) == 0:
        submitted_job_dir_list = []
    else:
        for i_indx in range(0, len(submitted_job_dir_list)):
            submitted_job_dir_list[i_indx] = os.path.abspath(submitted_job_dir_list[i_indx])
    while 1:
        #call monitor_timer
        for i_loop in range(0, num_loops):
            duration, submitted_job_dir_list, last_submitted_job_dir_list = monitor_timer(submitted_job_dir_list, last_submitted_job_dir_list)
            if duration > max_monitor_time:
                print('Maximum time (' + str(max_monitor_time_temp) + ' min) reached. Exiting...')
                vasp_tools.job_status(project_dir)
                exit()
            else:
                pass 
            ##print(i_loop, duration)
            time.sleep(job_check_frequency)


##    ######################################
##    # use threading
##    ######################################
##
##    def monitor_timer(submitted_job_dir_list, last_submitted_job_dir_list):
##        import random
##        ##print("Please execute your code here")
##        last_submitted_job_dir_list = submitted_job_dir_list
##        submitted_job_dir_list = queue_management(
##            project_dir = project_dir,
##            max_num_running_jobs = max_num_running_jobs,
##            submit_script_name = submit_script_name, 
##            queue_system = queue_system, 
##            job_kwd = job_kwd,
##            task_flow_list = task_flow_list,
##            )
##        ##submitted_job_dir_list[0] = submitted_job_dir_list[0] + str(random.randint(0,9))  # This line is only used for debugging
##        #build timer repeatly
##        timer = threading.Timer(job_check_frequency, monitor_timer, (submitted_job_dir_list, last_submitted_job_dir_list, ))
##        timer.start()
##        ##print('last_submitted: ',last_submitted_job_dir_list)
##        ##print('newly_submitted: ',submitted_job_dir_list) 
##        for i_job in submitted_job_dir_list:
##            if i_job in last_submitted_job_dir_list:
##                print('ERROR (from job_monitor): Job submission failed. Please check your inputs. The directories to be checked are: ' + str(i_job))
##                timer.cancel()
##         
##    #call monitor_timer
##    stop_indicator = 0
##    last_submitted_job_dir_list = []
##    submitted_job_dir_list = []
##    
##    timer = threading.Timer(0, monitor_timer, (submitted_job_dir_list, last_submitted_job_dir_list, ))
##    timer.start()
##
##    #stop after max_monitor_time
##    time.sleep(max_monitor_time)
##    timer.cancel()

    return 0 

def write_project_summary(project_dir, task_type = 'VASP', plot_fatband = False, band_struct_ylim = [-2, 2]):
    '''write project summary'''
    args_dict = locals()
    import os
    from . import task_manager
    project_dir = os.path.abspath(project_dir)
    task_name_list = [x for x in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, x))]
    for i_task_indx in range(0, len(task_name_list)):
        i_task_path = os.path.join(project_dir, task_name_list[i_task_indx])
        task_manager.write_task_summary(i_task_path, task_type = 'VASP', plot_fatband = plot_fatband, band_struct_ylim = band_struct_ylim)
    return 0
