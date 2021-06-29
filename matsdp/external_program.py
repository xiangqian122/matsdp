# -*- coding: utf-8 -*-

def mp_query(api_key_str = None, query_output_dir_name = None, criteria='*PS3', properties=['material_id', 'pretty_formula', 'cif'], final = True, elect_struct_type_list = ['bs'], save_bs_data = False, bs_data_format = 'json', save_bs_fig = False, timeout = 120):
    '''
    Query materials by some criteria and properties from Materials Project.
    '''
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    #from pymatgen.ext.matproj import MPRester
    from pymatgen import MPRester
    from pymatgen.io.cif import CifParser
    import os
    import re
    from . import funcs
    from . import default_params
    from pymatgen.io.vasp.outputs import Vasprun
    from pymatgen.electronic_structure.plotter import BSDOSPlotter,BSPlotter,BSPlotterProjected,DosPlotter
    import matplotlib.pyplot as plt
    import multiprocessing
    from multiprocessing import Pool, TimeoutError
    import time
    import platform
    import math

    start_time = time.time()
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    if query_output_dir_name in [None, 'None', 'none']:
        query_output_dir_name = 'mp_query'
    query_output_dir = os.path.join(output_dir, query_output_dir_name)
    funcs.mkdir(query_output_dir)

    # initializes the REST adaptor
    mpr = MPRester(api_key_str)
    
    # This supresses warnings.
    import warnings
    warnings.filterwarnings('ignore')
    # This is a helper function to shorten lists during the 
    # live presentation of this lesson for better readability. 
    # You can ignore it. 
    def shortlist(long_list, n=5):
        print("First {} of {} items:".format(min(n, 5), len(long_list)))
        for item in long_list[0:n]:
            print(item)

    mp_query_dict = {}
    mp_query_dict['material_id_list'] = []
            
    for i_property in ['material_id', 'pretty_formula', 'cif']:
        if i_property not in properties:
            properties.append(i_property)
    #data = mpr.query('*PS3', properties=["material_id", "pretty_formula"])
    data = mpr.query(criteria=criteria, properties=properties)
    #shortlist(data)
    print('#####################################################################################')
    print('criteria = ', str(criteria))
    print('query_output_dir_name = ', query_output_dir_name) 
    print('Number of entries found: ' + str(len(data)))

    file_path_list, file_name_list = funcs.get_files(dir_path = query_output_dir)
    for i_indx in range(0, len(data)):
        #primitive_structure = data[i_indx]['cif']
        material_id = data[i_indx]['material_id']
        pretty_formula = data[i_indx]['pretty_formula']
        poscar_name_1 = 'POSCAR_primitive_'+ pretty_formula + '_' + material_id
        poscar_name_1 = re.sub('[()]', '', poscar_name_1)
        mp_query_dict['material_id_list'].append(material_id)
        if poscar_name_1 in file_name_list:
            continue
        primitive_structure = CifParser.from_string(data[i_indx]['cif']).get_structures(primitive=True)[0]
        primitive_structure.to(filename= os.path.join(query_output_dir, poscar_name_1))
##        # get structure from material_id
##        structure = mpr.get_structure_by_material_id(data[i_indx]['material_id'])
    
    if save_bs_data == True:
        args = []
        for i_indx in range(0, len(data)):
            pretty_formula = data[i_indx]['pretty_formula']
            material_id = data[i_indx]['material_id']
            bs_json_file_name = 'bs_'+ pretty_formula + '_' + material_id + '.json'
            if bs_json_file_name in file_name_list:
                continue
            #primitive_structure = data[i_indx]['cif']
            #material_id_list.append(data[i_indx]['material_id'])
            material_id = data[i_indx]['material_id']
            args.append((api_key_str, query_output_dir_name, material_id, elect_struct_type_list, bs_data_format, save_bs_fig))
        print('number of json files to be retrieved = ', len(args))
        cores = math.ceil(multiprocessing.cpu_count() * 3 / 4)
        # https://note.qidong.name/2018/11/python-multiprocessing/
        # https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.AsyncResult.wait
        with Pool(cores) as p:
            #results = p.starmap(mp_get_elect_struct, args)
            results = p.starmap_async(mp_get_elect_struct, args)
            results.wait(timeout)
        try:
            print(results.get(timeout=timeout))
        except TimeoutError:
            print(TimeoutError.__name__)

    end_time = time.time()
    print('time=',end_time-start_time)
    return mp_query_dict

#api_key_str = ''
#mp_query(api_key_str = api_key_str, criteria={'spacegroup.number': 213, 'band_gap':0}, properties = ['material_id', 'pretty_formula', 'cif'])

def mp_get_elect_struct(api_key_str = None, query_output_dir_name = None, material_id = False, elect_struct_type_list = ['bs'], bs_data_format = 'json', save_fig = True):
    '''
    Get band structure data form Materials Project by materials_id
    elect_struct_type_list: 'bs' or 'dos'
    data_format: 'json', 'formed_json', 'npz'
    '''
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    #from pymatgen.ext.matproj import MPRester
    from pymatgen import MPRester
    from pymatgen.io.cif import CifParser
    import os
    import re
    from . import funcs
    from . import default_params
    from pymatgen.io.vasp.outputs import Vasprun
    from pymatgen.electronic_structure.plotter import BSDOSPlotter,BSPlotter,BSPlotterProjected,DosPlotter
    import matplotlib.pyplot as plt
    import pickle
    import json
    import numpy as np
    import time

    start_time = time.time()
    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    if query_output_dir_name in [None, 'None', 'none']:
        query_output_dir_name = 'mp_query'
    query_output_dir = os.path.join(output_dir, query_output_dir_name)
    funcs.mkdir(query_output_dir)

    # initializes the REST adaptor
    mpr = MPRester(api_key_str)

    data = mpr.query(criteria={'material_id':material_id}, properties=['pretty_formula'])
    pretty_formula = data[0]['pretty_formula']
    mat_label = pretty_formula + '_' + material_id
    ##print('{} starts'.format(mat_label))
    
    # This supresses warnings.
    import warnings
    warnings.filterwarnings('ignore')

    elect_struct_data_dict = {}
    elect_struct_data_dict['bs'] = None
    elect_struct_data_dict['dos'] = None

    bs_found = False
    dos_found = False
    bs_json_file_path = None
    dos_json_file_path = None

    write_json_file = False
    write_formed_json_file = False
    write_binary_file = False
    if bs_data_format == 'json':
        write_json_file = True
    elif bs_data_format == 'formed_json':
        write_formed_json_file = True
    elif bs_data_format == 'npz':
        write_binary_file = True

    if 'bs' in elect_struct_type_list:
        # get band structure from material_id
        try:
            # bs_data is a BandStructureSymmLine object, converting it to dict (.as_dict()) and save it to a file may cause a memory error)
            # https://pymatgen.org/pymatgen.electronic_structure.bandstructure.html
            bs_data = mpr.get_bandstructure_by_material_id(material_id = material_id)
            bs_json_file_path = os.path.join(query_output_dir, 'bs_'+ pretty_formula + '_' + material_id + '.json')
            bs_formed_json_file_path = os.path.join(query_output_dir, 'bs_formed_'+ pretty_formula + '_' + material_id + '.json')
            bs_npz_file_path = os.path.join(query_output_dir, 'bs_'+ pretty_formula + '_' + material_id + '.npz')
            if bs_data is not None:
                bs_found = True
                elect_struct_data_dict['bs'] = bs_data.as_dict()
                if write_json_file == True: 
                    with open(bs_json_file_path, 'w') as f_json:
                        json.dump(bs_data.as_dict(), f_json)
                if write_formed_json_file == True: 
                    with open(bs_formed_json_file_path, 'w') as f_formed_json:
                        data_formed = json.dumps(bs_data.as_dict(), sort_keys = True, indent = 4, separators=(',', ': '))
                        f_formed_json.write(data_formed)
                if write_binary_file == True: 
                    with open(bs_npz_file_path, 'wb') as f_npz:
                        np.savez_compressed(f_npz, bs = bs_data)      
        except:
            pass
    if 'dos' in elect_struct_type_list:
        # get DOS from material_id
        try:
            # get dos from material_id
            dos_json_file_path = os.path.join(query_output_dir, 'dos_'+ pretty_formula + '_' + material_id + '.json')
            dos_formed_json_file_path = os.path.join(query_output_dir, 'dos_formed_'+ pretty_formula + '_' + material_id + '.json')
            dos_npz_file_path = os.path.join(query_output_dir, 'dos_'+ pretty_formula + '_' + material_id + '.npz')
            dos_data = mpr.get_dos_by_material_id(material_id = material_id)
            ##bs_vasprun = Vasprun("vasprun.xml", parse_projected_eigen = True)
            ##dos_data = bs_vasprun.complete_dos 
            if dos_data is not None:
                dos_found = True
                elect_struct_data_dict['dos'] = dos_data.as_dict()
                with open(dos_json_file_path, 'w') as f_json:
                    json.dump(dos_data.as_dict(), f_json)
                if write_formed_json_file == True: 
                    with open(dos_formed_json_file_path, 'w') as f_formed_json:
                        data_formed = json.dumps(dos_data.as_dict(), sort_keys = True, indent = 4, separators=(',', ': '))
                        f_formed_json.write(data_formed)
                if write_binary_file == True: 
                    with open(dos_npz_file_path, 'wb') as f_npz:
                        np.savez_compressed(f_npz, dos=dos_data)      
        except:
            pass
    if bs_found and save_fig:
        dst_bs_fig_file_path = os.path.join(query_output_dir, 'bs_'+ pretty_formula + '_' + material_id + '.png')
        if write_json_file == True: 
            with open(bs_json_file_path, 'r') as f_json:
                bs_data_dict = json.load(f_json)
        if write_formed_json_file == True: 
            with open(bs_formed_json_file_path, 'r') as f_formed_json:
                bs_data_dict = json.load(f_formed_json)
        if write_binary_file == True: 
            with open(bs_npz_file_path, 'rb') as f_npz:
                npdata = np.load(f_npz)
                #bs_data = npdata['bs'].iterms()
                bs_data = npdata['bs'].item()
                npdata.close()
                bs_data_dict = bs_data.as_dict()
        bs_fig_file = mp_plot_bs(bs_data_dict = bs_data_dict, xlim = None, ylim = (-3,3))
        fpath, fname = os.path.split(bs_fig_file)
        funcs.mv(bs_fig_file, dst_bs_fig_file_path)
    end_time = time.time()
    #print(('{} ends. time=' + str(end_time-start_time)).format(mat_label))
    return elect_struct_data_dict

#def mp_plot_bs_from_json(bs_json_file_path, xlim = None, ylim = (-3,3)):
def mp_plot_bs(bs_data_dict, xlim = None, ylim = (-3,3)):
    '''
    plot band structure from json file.
    '''
    import os
    from pymatgen.electronic_structure.plotter import BSDOSPlotter,BSPlotter,BSPlotterProjected,DosPlotter
    from pymatgen.electronic_structure.bandstructure import BandStructure, BandStructureSymmLine
    import matplotlib.pyplot as plt
    import json
    import numpy as np
    import time
    from . import default_params
    from . import funcs

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    log_str = ''
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time.time()))
    plotter = BSPlotter(bs = BandStructureSymmLine.from_dict(bs_data_dict))
    plotter.get_plot(zero_to_efermi = True, vbm_cbm_marker = True, ylim = ylim, smooth = False)
    fig_file = os.path.join(output_dir, 'fig_' + str(formatted_time) + '_band.png')
    plt.savefig(os.path.join(fig_file), img_format='png')

    ##dos_json_file_path = None
    ##if dos_json_file_path not in [None, 'None', 'none']:
    ##    dos_json_file_path = os.path.abspath(dos_json_file_path)
    ##    with open(dos_json_file_path, 'r') as f:
    ##        dos_data = json.load(f)
    ##if dos_json_file_path not in [None, 'None', 'none']:
    ##    # bs with dos
    ##    pbandpdos_fig = BSDOSPlotter(bs_projection = 'elements', dos_projection = 'orbitals', vb_energy_range = 3, fixed_cb_energy = 5)
    ##    pbandpdos_fig.get_plot(bs = bs_data, dos = dos_data)
    ##    plt.savefig(os.path.join(bs_fpath, bs_fname + dos_fname + '.png'), img_format='png')
    return fig_file

def mp_get_structure_by_material_id(api_key_str = None, material_id = 'mp-1234', final = True, conventional_unit_cell = False):
    """
    Get a Structure with a certain material_id.

    Args:
        material_id (str): Materials Project material_id (a string,
            e.g., mp-1234).
        final (bool): Whether to get the final structure, or the initial
            (pre-relaxation) structure. Defaults to True.
        conventional_unit_cell (bool): Whether to get the standard
            conventional unit cell

    Returns:
        Structure object.
    """
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen import MPRester
    import os
    from . import funcs
    from . import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    output_dir = os.path.join(os.getcwd(), defaults_dict['output_dir_name'])
    funcs.mkdir(output_dir)

    # initializes the REST adaptor
    mpr = MPRester(api_key_str)
    structure = mpr.get_structure_by_material_id(mp_id_str)
    sga = SpacegroupAnalyzer(structure)
    conventional_structure = sga.get_conventional_standard_structure()
    ### Get a primitive version of the Structure
    primitive_structure = structure.get_primitive_structure()
    primitive_structure.to(filename = os.path.join(output_dir, 'POSCAR_primitive_' + material_id))
    conventional_structure.to(filename = os.path.join(output_dir, 'POSCAR_conentional' + material_id))
    return structure

def mp_cif2poscar(cif_file_path, poscar_file_path):
    import os
    from pymatgen.io.cif import CifParser
    try:
        cif_file_path = os.path.abspath(cif_file_path)
        poscar_file_path = os.path.abspath(poscar_file_path)
        parser = CifParser(cif_file_path)
        structure = parser.get_structures()[0]
        structure.to(filename=poscar_file_path)
        return 0
    except:
        print('# WARNING #2104252136 (from external_program): cif2poscar failed: ' + str(cif_file_path))
        return 1

#mp_id_str = 'mp-3342'
#get_structure_by_material_id(api_key_str = api_key_str, material_id = mp_id_str, final = True, conventional_unit_cell = False)

