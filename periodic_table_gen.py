def periodic_tab_gen():
    '''
    Part of the information is taken from PerioDict Table from Paul A Freshney 2010
    http://perioDicttableexplorer.com/pc_pte.htm
    The covalent radius is taken from the literature: [Beatriz Cordero et al. Covalent radii revisited. Dalton Transactions, 2008 2832--2838.]
    
    This is only a module for the developers!    
    This module generates periodic_table.py which contains the periodic_table_dict (dictionary of the periodic table)

    This module relies on the periodic_table.csv data file located in the same directory.

    This module relies on the periodic_table_bond.csv data file located in the same directory.
    The bond pair information is partially taken from the VESTA configuration file (style.ini)
    '''
    import os
    import numpy as np
    from matsdp import funcs

    periodic_table_pyfile = os.path.abspath('./matsdp/periodic_table.py')
    periodic_table_dict = {}
    periodic_table_file = os.path.join(os.path.dirname(__file__), 'periodic_table.csv')
    with open(periodic_table_file,'r') as f:
        lines = f.readlines()
        line_num = len(lines)
        col_num = len(funcs.split_line(lines[0],','))
        #Define two basic lists
        property_list = funcs.split_line(lines[0],',')
        elmt_symbol_list = [None] * (line_num - 1)
        for i in range(1,line_num):         
            elmt_symbol_list[i - 1] = str(funcs.split_line(lines[i],',')[1])

        for i in range(0,col_num):
            periodic_table_dict[property_list[i]] = {}
            
        #Atomic Number
        col_indx = property_list.index('atomic_number')
        for i in range(0,line_num - 1):
            periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = int(funcs.split_line(lines[i + 1],',')[col_indx])
        #symbol and name
        col_indx1 = property_list.index('symbol')
        col_indx2 = property_list.index('name')
        for i in range(0,line_num - 1):
                periodic_table_dict[property_list[col_indx1]][elmt_symbol_list[i]] = str(funcs.split_line(lines[i + 1],',')[col_indx1])
                periodic_table_dict[property_list[col_indx2]][elmt_symbol_list[i]] = str(funcs.split_line(lines[i + 1],',')[col_indx2])
        # elmt_color
        col_indx = property_list.index('elmt_color')
        for i in range(0,line_num - 1):
            if funcs.split_line(lines[i + 1],',')[col_indx] != 'None':
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = str(funcs.split_line(lines[i + 1],',')[col_indx])
            else:
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = None
        # recommended POTCAR setting for VASP: PAW
        col_indx = property_list.index('vasppot_paw')
        for i in range(0,line_num - 1):
            if funcs.split_line(lines[i + 1],',')[col_indx] != 'None':
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = '_' + str(funcs.split_line(lines[i + 1],',')[col_indx])
            else:
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = ''
        # recommended POTCAR setting for VASP: GW
        col_indx = property_list.index('vasppot_gw')
        for i in range(0,line_num - 1):
            if funcs.split_line(lines[i + 1],',')[col_indx] != 'GW':
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = '_' + str(funcs.split_line(lines[i + 1],',')[col_indx])
            else:
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = '_GW'
        # valence electron configuration. Also refer to the recommended PAW POTCAR of VASP
        col_indx = property_list.index('valence_elec')
        for i in range(0, line_num -1):
            if funcs.split_line(lines[i + 1],',')[col_indx] != 'None':
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = str(funcs.split_line(lines[i + 1],',')[col_indx])
            else:
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = None
        # number of valence electrons
        col_indx = property_list.index('num_valence_elec')
        for i in range(0, line_num -1):
            if funcs.split_line(lines[i + 1],',')[col_indx] != 'None':
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = int(funcs.split_line(lines[i + 1],',')[col_indx])
            else:
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = None
        ##########################
        # for all float type data
        ###########################
        # atomic_radius, unit: pm
        col_indx = property_list.index('atomic_radius')
        for i in range(0,line_num - 1):
            if funcs.split_line(lines[i + 1],',')[col_indx] != 'None':
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = float(funcs.split_line(lines[i + 1],',')[col_indx])
            else:
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = None
        # metallic_radius, unit: pm
        col_indx = property_list.index('metallic_radius')
        for i in range(0,line_num - 1):
            if funcs.split_line(lines[i + 1],',')[col_indx] != 'None':
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = float(funcs.split_line(lines[i + 1],',')[col_indx])
            else:
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = None
        # csd_covalent_radius, unit: pm
        col_indx = property_list.index('csd_covalent_radius')
        for i in range(0,line_num - 1):
            if funcs.split_line(lines[i + 1],',')[col_indx] != 'None':
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = float(funcs.split_line(lines[i + 1],',')[col_indx])
            else:
                periodic_table_dict[property_list[col_indx]][elmt_symbol_list[i]] = None
    ############################################################
    # for the dos_mode. dos_mode is useful for plotting the DOS
    ############################################################
    periodic_table_dict['dos_mode'] = {}
    for i_elmt in range(0, len(periodic_table_dict['atomic_number'])):
        elmt_number = i_elmt + 1
        elmt_symbol = funcs.get_keys_by_value(periodic_table_dict['atomic_number'], elmt_number)[0]
        if elmt_number >= 1 and elmt_number <= 4:
            periodic_table_dict['dos_mode'][elmt_symbol] = ['s']
        elif elmt_number >= 5 and elmt_number <= 20:
            periodic_table_dict['dos_mode'][elmt_symbol] = ['s','p']
        elif elmt_number >= 21:
            periodic_table_dict['dos_mode'][elmt_symbol] = ['s','p','d']    

    #############################################################
    # write the periodic_table_dict to the periodic_table_pyfile
    #############################################################
    with open(periodic_table_pyfile,'w') as f:
        f.write('# -*- coding: utf-8 -*-' + '\n' +
                'def periodic_tab():' + '\n' +
                '  periodic_table_dict = {}' + '\n' +
                '  periodic_table_dict[\'atomic_number\']=' + str(periodic_table_dict['atomic_number']) + '\n' + 
                '  periodic_table_dict[\'symbol\']=' + str(periodic_table_dict['symbol']) + '\n' +
                '  periodic_table_dict[\'name\']=' + str(periodic_table_dict['name']) + '\n' +
                '  periodic_table_dict[\'elmt_color\']=' + str(periodic_table_dict['elmt_color']) + '\n' +
                '  periodic_table_dict[\'atomic_radius\']=' + str(periodic_table_dict['atomic_radius']) + '\n' +
                '  periodic_table_dict[\'metallic_radius\']=' + str(periodic_table_dict['metallic_radius']) + '\n' +
                '  periodic_table_dict[\'csd_covalent_radius\']=' + str(periodic_table_dict['csd_covalent_radius']) + '\n' + 
                '  periodic_table_dict[\'vasppot_paw\']=' + str(periodic_table_dict['vasppot_paw']) + '\n' + 
                '  periodic_table_dict[\'vasppot_gw\']=' + str(periodic_table_dict['vasppot_gw']) + '\n' + 
                '  periodic_table_dict[\'dos_mode\']=' + str(periodic_table_dict['dos_mode']) + '\n' +
                '  periodic_table_dict[\'valence_elec\']=' + str(periodic_table_dict['valence_elec']) + '\n' +
                '  periodic_table_dict[\'num_valence_elec\']=' + str(periodic_table_dict['num_valence_elec']) + '\n' +
                '  return periodic_table_dict'
                )
    return periodic_table_dict

periodic_tab_gen()
