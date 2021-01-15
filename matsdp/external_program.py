# -*- coding: utf-8 -*-

def mp_query(api_key_str = None, query_output_dir_name = None, criteria='*PS3', properties=['material_id', 'pretty_formula', 'cif'], final = True):
    '''
    Query materials by some criteria and properties from Materials Project.
    '''
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    #from pymatgen.ext.matproj import MPRester
    from pymatgen import MPRester
    from pymatgen.io.cif import CifParser
    import os
    from . import funcs
    from . import default_params

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
            
    for i_property in ['material_id', 'pretty_formula', 'cif']:
        if i_property not in properties:
            properties.append(i_property)
    #data = mpr.query('*PS3', properties=["material_id", "pretty_formula"])
    data = mpr.query(criteria=criteria, properties=properties)
    #shortlist(data)
    print('Number of entries found: ' + str(len(data)))

    for i_indx in range(0, len(data)):
        #primitive_structure = data[i_indx]['cif']
        primitive_structure = CifParser.from_string(data[i_indx]['cif']).get_structures(primitive=True)[0]
        primitive_structure.to(filename= os.path.join(query_output_dir, 'POSCAR_primitive_'+ data[i_indx]['pretty_formula'] + '_' + data[i_indx]['material_id']))

##        # get structure from material_id
##        structure = mpr.get_structure_by_material_id(data[i_indx]['material_id'])

##        # get band structure from material_id
##        bs_found = False
##        try:
##            bs = mpr.get_bandstructure_by_material_id(material_id = data[i_indx]['material_id'])
##            if bs is not None:
##                bs_found = True
##        except:
##            pass
##        if bs_found:
##            plotter = BSPlotter(bs)
##            plotter.show()
##            pass

##        # get dos from material_id
##        dos = mpr.get_dos_by_material_id(data[i_indx]['material_id'])

    return 0

#api_key_str = ''
#mp_query(api_key_str = api_key_str, criteria={'spacegroup.number': 213, 'band_gap':0}, properties = ['material_id', 'pretty_formula', 'cif'])

def mp_get_structure_by_material_id(api_key_str = None, material_id = 'mp-1234', final = True, conventional_unit_cell = False):
    """
    Get a Structure wither a certain material_id.

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

#mp_id_str = 'mp-3342'
#get_structure_by_material_id(api_key_str = api_key_str, material_id = mp_id_str, final = True, conventional_unit_cell = False)

