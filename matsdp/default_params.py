# -*- coding: utf-8 -*-
def default_params():
    '''
    Description:
        It defines the default parameters of the program.
    Args:
        None
    Return:
        defaults_dict
    '''
    defaults_dict = {}
    defaults_dict['version'] = '0.2.1'
    defaults_dict['logfile'] = 'matsdp.log'
    defaults_dict['output_dir_name'] = 'outputs'
    defaults_dict['projects_dir_name'] = 'projects'
    defaults_dict['projects_summary_dir_name'] = 'projects_summary'
    defaults_dict['task_summary_dir_name'] = 'task_summary'
    defaults_dict['test_dir_name'] = 'test'
    defaults_dict['greek_capital_letter_list'] = ['Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Eta', 'Theta', 'Iota', 'Kappa', 'Lambda', 'Mu', 'Nu', 'Xi', 'Omicron', 'Pi', 'Rho', 'Sigma', 'Tau', 'Upsilon', 'Phi', 'Chi', 'Psi', 'Omega']
    defaults_dict['greek_small_letter_list'] = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'omicron', 'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega']
    return defaults_dict
