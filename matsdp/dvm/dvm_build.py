def create_multiple_dvm_jobs(poscar_file_path_dict, origin_atom_name_list, radius, elmt_ind_file_dir):
    '''
    create multiple DVM jobs (the *.incar, *.input, IND.DAT files will be created) based on atom selection (spherical) of the POSCAR files.
    poscar_file_path_dict: Dictionary type. A dictionary which contains the POSCAR file path, the key of the dictionary will be used as part of the DVM job name.
    origin_atom_name_list: the origin atom in the center of the sphere in the atom selection (spherical) of the POSCAR
    radius: the radius of the sphere in the atom selection (spherical) of the POSCAR
    elmt_ind_file_dir: the top directory which contains the IND.DAT files of the elements
    '''
    import os
    from .. import funcs
    from ..vasp import vasp_build
    from . import dvm_write
    from .. import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.getcwd() + '/' + defaults_dict['output_dir_name']
    funcs.mkdir(output_dir)

    job_id_list = list(poscar_file_path_dict.keys())
    for i_job in range(0, len(job_id_list)):
        for i_elmt_name in origin_atom_name_list:
            job_name = job_id_list[i_job] + '_' + i_elmt_name + '_R' + str(radius)
            pos_file_path_str = output_dir + '/atom_selection_sphere/' + job_name + '/' + job_name + '.vasp'
            # create the *.incar file and the *.vasp file
            vasp_build.selection_sphere(
                poscar_file_path = poscar_file_path_dict[job_id_list[i_job]],
                origin_atom_name = i_elmt_name,
                radius = radius,
                include_mirror_atoms = True,
                output_file_name = job_id_list[i_job] + '_' + i_elmt_name + '_R' + str(radius)
                )
            # create the *.input file
            dvm_write.write_input(
                pos_file_path = pos_file_path_str,
                dvm_input_dict = None
                )
            # create the IND.DAT file
            dvm_write.write_ind(
                pos_file_path = pos_file_path_str,
                elmt_ind_file_dir = elmt_ind_file_dir,
                )
    funcs.write_log(
        logfile,
        'dvm_build.create_multiple_dvm_jobs(\n' +
        '    poscar_file_path_dict = ' + str(poscar_file_path_dict) + ',\n' +
        '    origin_atom_name_list = ' + str(origin_atom_name_list) + ',\n' +
        '    radius = ' + str(radius) + ',\n' +
        '    elmt_ind_file_dir = ' + 'r\'' + str(elmt_ind_file_dir) + '\')' + '\n'
        )
    return 0
