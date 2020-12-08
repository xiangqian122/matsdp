# -*- coding: utf-8 -*-
##########################################
#Folder and file manipulation
##########################################
def cp(src_filedir,dst_filedir, suppress_warning = False):
    import os
    import shutil
    if not os.path.isfile(src_filedir):
        if suppress_warning == False:
            print('WARNING: funcs error. The file ' + src_filedir + ' does not exist!')
    else:
        fpath, fname = os.path.split(dst_filedir)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        if src_filedir == dst_filedir:
            if suppress_warning == False:
                print('WARNING(from funcs.cp):' + dst_filedir + ' already exist and is skipped.')
        else:
            shutil.copyfile(src_filedir, dst_filedir)
    return 0
        
def file_path_name(file_dir):
    '''Get parent path and file name for the input file directory'''
    import os
    
    parent_path = os.path.dirname(file_dir)  
    filename = os.path.split(file_dir)[-1]
    return parent_path, filename

def grep(kwd, file_dir, suppress_warning = False):
    '''extract the line with the keyword(kwd) specified
    the function is similar to the "grep" command in the linux shell
    '''
    import os
    file_dir = os.path.abspath(file_dir)
    if os.path.exists(file_dir) and os.path.getsize(file_dir) > 0:
        pass
    else:
        if suppress_warning == False:
            print('ERROR (from funcs.grep): the file ' + str(file_dir) + ' does not exist or is empty')
        ##funcs.write_log(logfile, '#' + file_dir + " doesn't exist of is empty")
        ##quit()
    line_str = ''
    line_list = []
    try:
        with open(file_dir,'r') as f:
            lines = f.readlines()
            for line in lines:
                if kwd in line:
                    line_list.append(line)
    except:
        pass
    line_str = ''.join(line_list)
    return line_str

def line_num(input_file, suppress_warning = False):
    'Get number of lines from a file'
    import os
    input_file = os.path.abspath(input_file)
    if os.path.exists(input_file) and os.path.getsize(input_file) > 0:
        pass
    else:
        if suppress_warning == False:
            print('funcs error: the file ' + str(input_file) + ' does not exist or is empty')
        funcs.write_log(logfile, '#' + input_file + " doesn't exist of is empty")
        quit()
    count=0  
    thefile=open(input_file)  
    while True:  
        buffer=thefile.read(1024*8192)  
        if not buffer:  
            break  
        count+=buffer.count('\n')  
    thefile.close()  
    return count

def merge_files(file1,file2):
    '''
    Merge Files
    append contents of file2 to file1
    '''
    import os
    cwd = os.getcwd()
    file1_dir = cwd + '/' + file1
    file2_dir = cwd + '/' + file2
    with open(file1_dir,'a') as file1_obj:
        with open(file2_dir) as file2_obj:
            file2line = file2_obj.readlines()
            file1_obj.write("".join(file2line))
    return 0

def mkdir(path):
    'create a directory'
    import os
    path = path.strip()
    path = path.rstrip("/")
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        #print(path + ' already exists')
        return False

def mv(src_filedir,dst_filedir):
    import os
    import shutil
    if not os.path.isfile(src_filedir):
        print(src_filedir + ' does not exist!')
    else:
        fpath,fname=os.path.split(dst_filedir)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        shutil.move(src_filedir,dst_filedir)
    return 0

def write_log(logfile,output_str):
    '''Write specific string into a log file,
       When the log file exists, then append output to log file.
       If the log file doesn't exist, then create it and write output to log file.
       Note that output_str must be string format'''
    import os
    import time
    from . import default_params
    defaults_dict = default_params.default_params()
    version = defaults_dict['version']

    current_time = time.time()
    formatted_time = time.strftime('## time: %Y%m%d %H:%M:%S',time.localtime(current_time))
    if os.path.exists(logfile) == False or (os.path.exists(logfile) and os.path.getsize(logfile) == 0):
        # write header to the log file
        with open(logfile,'w') as logfile_object:
            logfile_object.write(
                '################log file###############\n' +
                '########version: matsdp-' + version + '##########\n' +
                '# -*- coding: utf-8 -*-' + '\n' +
                '##possible imports:' + '\n' +
                'import matsdp' + '\n' +
                'from matsdp import convert' + '\n' +
                'from matsdp import funcs' + '\n' +
                'from matsdp import vasp' + '\n' +
                'from matsdp import apt' + '\n' +
                'from matsdp import dvm' + '\n' +
                'from matsdp import pms' + '\n' +
                'from matsdp.vasp import vasp_read' + '\n' +
                'from matsdp.vasp import vasp_plot' + '\n' +
                'from matsdp.vasp import vasp_analyze' + '\n' +
                'from matsdp.vasp import vasp_build' + '\n' +
                'from matsdp.vasp import vasp_write' + '\n' +
                'from matsdp.vasp import vasp_tools' + '\n' +
                'from matsdp.vasp import vasp_default' + '\n' +
                'from matsdp.vasp import vasp_help' + '\n' +
                'from matsdp.apt import apt_read' + '\n' +
                'from matsdp.apt import apt_plot' + '\n' +
                'from matsdp.dvm import dvm_read' + '\n' +
                'from matsdp.dvm import dvm_analyze' + '\n' +
                'from matsdp.dvm import dvm_build' + '\n' +
                'from matsdp.dvm import dvm_write' + '\n' +
                'from matsdp.dvm import dvm_default' + '\n' +
                'from matsdp.dvm import dvm_help' + '\n' +
                'from matsdp.pms import project_manager' + '\n' +
                'from matsdp.pms import task_manager' + '\n' +
                '######################################\n')
        with open(logfile,'a') as logfile_object:
            logfile_object.write(
                formatted_time + '\n' +
                output_str +  '\n') 
    else:
        with open(logfile,'a') as logfile_object:
            logfile_object.write(
                formatted_time + '\n' +
                output_str + '\n')
    return 0

def rm_empty_lines_in_file(in_file):
    '''
    remove emtpy lines in a file
    '''
    import os
    in_file = os.path.abspath(in_file)
    fpath,fname = os.path.split(in_file)
    fname_temp = fname + '.temp'
    ftemp = os.path.join(fpath,fname_temp)
    with open(in_file) as f1, open(ftemp, 'w') as f2:
        for line in f1:
            if not line.strip():
                continue
            f2.write(line)
    mv(ftemp,in_file)
    return 0

def get_file_str(file_path):
    
    '''
    get file content, and convert it to string
    '''
    import os
    file_path = os.path.abspath(file_path)
    file_str = ''
    with open(file_path, 'r') as f:
        line = f.readlines()
        for i_line in line:
            file_str = file_str + i_line
    return file_str

def write_file(input_str, dest_file_path, mode = 'w'):
    '''
    write sepcified string to a desiginated file
    input_str: the string which is to be included in the destination file
    dest_file_path: destination file path, the specific string will be written to the destination file
    mode: the mode for writing. The values could be 'w' or 'a' (write or append)
    if the destination file exists, then it will be created.
    '''    
    import os
    input_str = str(input_str)
    dest_file_path = os.path.abspath(dest_file_path)
    workdir, dest_file = file_path_name(dest_file_path)
    mkdir(workdir)
    with open(dest_file_path, mode) as f:
        f.write(input_str)
    return 0

def touch(file_path):
    '''create file'''
    import os
    file_path = os.path.abspath(file_path)
    workdir, file = file_path_name(file_path)
    mkdir(workdir)
    with open(file_path, 'w') as f:
        pass
    return 0

def dir_tree(input_dir):
    '''
    list the directory tree structure
    '''
    import os
    import time
    current_time = time.time()
    formatted_time = time.strftime('%Y%m%d %H:%M:%S',time.localtime(current_time))

    input_dir = os.path.abspath(input_dir)
    if os.path.exists(input_dir) == True and os.path.isdir(input_dir) == True:
        pass
    else:
        print('funcs error: the directory ' + input_dir + ' does not exist. Exiting ...')
        quit()
    dir_tree_str = ''
    for root, dirs, files in os.walk(input_dir):
        dirs.sort()
        level = root.replace(input_dir, '').count(os.sep)
        indent = ' ' * 4 * (level)
        dir_tree_str = dir_tree_str + '{}{}/'.format(indent, os.path.basename(root)) + '\n'
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            dir_tree_str = dir_tree_str + '{}{}'.format(subindent, f) + '\n'
    input_dir_tree_file = os.path.join(input_dir, 'dir_tree.txt')
    with open(input_dir_tree_file, 'w') as f:
        f.write(
            'directory:\n' +
            input_dir + '\n\n' +
            'time:\n' +
            formatted_time + '\n\n' +
            'directory tree structure:\n' +
            dir_tree_str
            )
    return dir_tree_str

def file_status(file_path, suppress_warning = False):
    '''
    check the status of the file: file exist, file not found, empty file etc.
    return value: the status of the file
    '''
    import os
    from . import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    status = 0
    file_path = os.path.abspath(file_path)
    if os.path.exists(file_path) == False:
        # file does not exist
        status = 0
        if suppress_warning == False:
            print('# WARNING #20120305 (from funcs.file_status): The file ' + str(file_path) + ' does not exist')
        ##write_log(logfile, '# WARNING (from funcs.file_status): The file ' + file_path + " does not exist")
    elif os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        # the file exists and is non-empty
        status = 1
    elif os.path.exists(file_path) and os.path.getsize(file_path) == 0:
        # emtpy file exist
        status = 2
        if suppress_warning == False:
            print('# WARNING 20120306 (from funcs.file_status): The file ' + str(file_path) + ' is empty')
        ##write_log(logfile, '# WARNING (from funcs.file_status): The file ' + file_path + " is empty")
    return status

def dir_status(dir_path):
    '''
    check the status of the directory: dir exist, dir not found, empty dir
    return value: the status of the directory
    '''
    import os
    from . import default_params

    defaults_dict = default_params.default_params()
    logfile = defaults_dict['logfile']
    status = 0
    dir_path = os.path.abspath(dir_path)
    if os.path.exists(dir_path) == False:
        # file directory not exist
        status = 0
        print('# WARNING: The directory ' + str(dir_path) + ' does not exist')
        write_log(logfile, '# WARNING: The directory ' + dir_path + " does not exist")
    elif os.path.exists(dir_path) and os.path.isdir(dir_path) ==True:
        # the directory exists
        status = 1
    return status

def get_files(dir_path, extension = None):
    '''
    Get files in a specific directory
    dir_path: the path of the directory.
    extension: get the files with the extension "extension_filter".
    '''
    import os
    file_path_list = []
    file_name_list = []
    for ifile in os.listdir(dir_path):
        file_path = os.path.join(dir_path, ifile)
        if os.path.isfile(file_path):
            if extension not in [None,'None','none']:
                if file_path.endswith(extension):
                    file_path_list.append(file_path)
                    file_name_list.append(ifile)
            else:
                file_path_list.append(file_path)
                file_name_list.append(ifile)
    return file_path_list, file_name_list

def get_dirs(dir_path):
    '''
    get all the directories (include the directory itself and its subdirectories) in dir_path and put them in a list
    '''
    import os
    dir_list = []
    path1, name1 = os.path.split(dir_path)
    path_recorder_list = [path1]
    path_level_list = [None] * 99
    temp_indx = 0
    for root, dirs, files in os.walk(dir_path):
        dirs.sort()
        temp_indx = temp_indx + 1
        level = root.replace(dir_path, '').count(os.sep)
        path_level_list[level] = os.path.basename(root)
        ##if level >= 1:
        ##    path_recorder_list = [path1] + path_level_list[:level]
        path_recorder_list = [path1] + path_level_list[:level]
        path_recorder_list.append(os.path.basename(root))
        i_dir_path = os.path.join('', *path_recorder_list)
        dir_list.append(i_dir_path)
    return dir_list

########################################
#String manipulation
##########################################
def extract_num(input_str):
    '''Extract number from string using regular expression, and put them in a list'''
    import re
    num = re.findall(r"\d+\.?\d*", input_str)
    return num

def extract_alpha_str(input_str):
    '''Extract english letters from string using regular expression, and put them in a list'''
    import re
    alpha_str = [x for x in re.split(r'[^A-Za-z]', input_str) if x != '']
    return alpha_str

def insert_str(old_str,indx,inserted_str):
    '''Insert a specified string(IndsetedStr) into sepecified location of an input string(old_str)
    indx is where you want to insert your string
    inserted_str is the string that you want to insert'''
    new_str = old_str[:indx] + inserted_str + old_str[indx:]
    return new_str

def replace(file_path, pattern, subst):
    '''
    Search and replace a string in a line of a file
    https://stackoverflow.com/questions/39086/search-and-replace-a-line-in-a-file-in-python
    '''
    from tempfile import mkstemp
    from shutil import move
    from os import fdopen, remove
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
    return 0

def replace_file_content(file_path, old_str, new_str):
    'Replace string in file'
    import os
    file_obj = os.path.abspath(file_path)
    f = open(file_obj,'r+')
    all_the_lines=f.readlines()
    f.seek(0)
    f.truncate()
    for line in all_the_lines:
        f.write(line.replace(old_str,new_str))
    f.close()
    return 0
  
def split_line(line, separator = ' '):
    '''split the line using the specified separator, the default separator is white space 
    and put the splitted line string into a list'''
    temp_list = [item.strip() for item in str(line).replace('\t', '    ').strip('\n').strip().split(separator)]
    while '' in temp_list:
        temp_list.remove('')
    return temp_list

def find_kwd_line_indx(kwd, file_path):
    '''
    find the keyword line index in the file (the line index starts from 0)
    the list which contains the line indices will be returned
    '''
    import os
    file_path = os.path.abspath(file_path)
    kwd = str(kwd)
    with open(file_path, 'r') as f:
        lines = f.readlines()
        indx = [x for x in range(len(lines)) if kwd in lines[x]]
    return indx

def str_format(input_str, max_len = 13):
    '''
    get formatted string
    e.g. str_format('aaa', 6) would return 'aaa   '
    '''
    if str(max_len).isdigit() == False:
        print('ERROR: funcs error. Not an integer')
        exit()
    input_str = str(input_str)
    if len(input_str) > max_len:
        max_len = len(input_str) + 1
    formatted_str = input_str + ' ' * (max_len - len(input_str))
    return formatted_str
    
def float_format(input_float, max_len = 13, format_str = '{:.2f}'):
    '''
    get formatted float
    e.g. float_format('3.1415926', 6, '{:.2f}') would return '3.14  '
    '''    
    if str(max_len).isdigit() == False:
        print('ERROR: funcs error. Not an integer')
        exit()
    try :
        float(input_float)
    except ValueError :
        print('ERROR: funcs error. Not a float')
    if len(str(input_float)) > max_len:
        max_len = len(str(input_float)) + 1
    formatted_float = format_str.format(float(input_float)) + ' ' * (max_len - len(format_str.format(float(input_float))))
    return formatted_float

def str_strip(text):
    '''string manipulation: strip'''
    try:
        return text.strip()
    except AttributeError:
        return text

def str2int(text):
    '''
    convert string to integer
    Note: if there are NaN values in the column of integers, the whole column would be designated as float type.
    '''
    try:
        return int(text.strip('" ').strip())
    except AttributeError:
        return text
    except ValueError:
        return None

def str2float(text):
    '''convert string to float'''
    import math
    try:
        return float(text.strip('" ').strip())
    except AttributeError:
        return text
    except ValueError:
        return math.nan

##########################################
#Other functions
##########################################
def convert_bool(var):
    '''
    tackle the problem: bool('False') returns True
    convert var to real boolean type variable
    '''
    if type(var) == str:
        if var == 'False':
            var = False
        elif var == 'True':
            var = True
    if type(var) == bool:
        pass
    return var

def logic_retn_val(input_logic_val, true_return_val, false_return_val):
    '''Return value according to the input logic value
    If input logic value is True then return true_return_val
    If input logic value is False then return false_return_val'''
    if input_logic_val == True:
        return true_return_val
    elif input_logic_val == False:
        return false_return_val


def get_keys_by_value(input_dict, value_to_find):
    '''
    Get a list of keys from the dictionary which has the given value
    '''
    list_of_keys = list()
    list_of_items = input_dict.items()
    for item in list_of_items:
        if item[1] == value_to_find:
            list_of_keys.append(item[0])
    return  list_of_keys

def sigmoid(x, a, b, c, d):
    '''
    fitting using the sigmoid function
    '''
    import numpy as np
    return a / (1. + np.exp(-c * (x - d))) + b

def linear(x, a, b):
    '''
    fitting using the linear function
    '''
    import numpy as np
    return a * x + b

def data_normalize(input_data_value, data_list, colormap_vmin = None, colormap_vmax = None):
    '''
    data normalization
    '''
    if min(data_list) != max(data_list):
        if colormap_vmin == None or colormap_vmax == None:
            colormap_vmin = min(data_list)
            colormap_vmax = max(data_list)
        elif colormap_vmin != None and colormap_vmax != None and colormap_vmin >= colormap_vmax:
            colormap_vmin = min(data_list)
            colormap_vmax = max(data_list)        
        elif colormap_vmin != None and colormap_vmax != None and colormap_vmin < colormap_vmax:
            ##print(type(colormap_vmin), type(min(data_list)))
            if (colormap_vmin <= min(data_list)) and (colormap_vmax >= max(data_list)):
                pass
            elif (colormap_vmin > min(data_list)) and (colormap_vmax >= max(data_list)):
                colormap_vmin = min(data_list)
            elif (colormap_vmin <= min(data_list)) and (colormap_vmax < max(data_list)):
                colormap_vmax = max(data_list)
            elif (colormap_vmin > min(data_list)) and (colormap_vmax < max(data_list)):
                colormap_vmin = min(data_list)
                colormap_vmax = max(data_list)
        normalized_value = (input_data_value - colormap_vmin) / (colormap_vmax - colormap_vmin)
    else:
        real_data_value = min(data_list)
        if colormap_vmin == None or colormap_vmax == None:
            colormap_vmin = real_data_value - 1
            colormap_vmax = real_data_value + 1          
        elif colormap_vmin != None and colormap_vmax != None and colormap_vmin >= colormap_vmax:
            colormap_vmin = real_data_value - 1
            colormap_vmax = real_data_value + 1          
        elif colormap_vmin != None and colormap_vmax != None and colormap_vmin < colormap_vmax:
            if colormap_vmax < real_data_value:
                colormap_vmin = real_data_value - 1
                colormap_vmax = real_data_value + 1          
            elif colormap_vmin > real_data_value:
                colormap_vmin = real_data_value - 1
                colormap_vmax = real_data_value + 1          
            elif colormap_vmin <= real_data_value and colormap_vmax >= real_data_value:
                pass
        normalized_value = (real_data_value - colormap_vmin) / (colormap_vmax - colormap_vmin)
    return normalized_value, colormap_vmin, colormap_vmax

def userdefined_colormap(data_list , colormap_vmin = None, colormap_vmax = None, vmax_color = 'red', vmin_color = 'blue',
                         alignment = 'horizontal', left = None, bottom = None, width = None, height = None):
    '''
    plot colormap by color mixing
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

    color_marker_num = 100
    # get colormap_vmin and colormap_vmax
    colormap_normalized_data, colormap_vmin, colormap_vmax = data_normalize(input_data_value = float(data_list[0]),
                                                                            data_list = data_list.astype(float),
                                                                            colormap_vmin = colormap_vmin,
                                                                            colormap_vmax = colormap_vmax)

    y_dense_arr = np.linspace(colormap_vmin, colormap_vmax, num = color_marker_num,endpoint=True)
    x_dense_arr = [0] * len(y_dense_arr)

    af = plt.gcf()
    ax = plt.gca()
##    fig_number = plt.gcf().number
    # define subplot position: [left, bottom, width, height]
    if left == None or bottom == None or width == None or height == None:
##        subplot_position = [0.1, 0.1, 0.8, 0.8]
        if alignment == 'vertical':
            subplot_position = [0.52, 0.17, 0.7, 0.7]
        elif alignment == 'horizontal':
            subplot_position = [0.08, -0.28, 0.7, 0.7]
    else:
        subplot_position = [left, bottom, width, height]
    af.add_axes(subplot_position)
    if alignment == 'vertical':
        pass
    elif alignment == 'horizontal':
        temp = y_dense_arr.copy()
        y_dense_arr = x_dense_arr.copy()
        x_dense_arr = temp.copy()
    colormap_marker_size = 30
    plt.plot(0, 0, color = 'white', markersize = 0)
    for i in range(0,color_marker_num):
        if alignment == 'vertical':
            colormap_normalized_data, colormap_vmin, colormap_vmax = data_normalize(input_data_value = y_dense_arr[i],
                                                                                    data_list = y_dense_arr,
                                                                                    colormap_vmin = colormap_vmin,
                                                                                    colormap_vmax = colormap_vmax)
        elif alignment == 'horizontal':
            colormap_normalized_data, colormap_vmin, colormap_vmax = data_normalize(input_data_value = x_dense_arr[i],
                                                                                    data_list = x_dense_arr,
                                                                                    colormap_vmin = colormap_vmin,
                                                                                    colormap_vmax = colormap_vmax)
        # higher spectrum of the colormap
        plt.plot(x_dense_arr[i],
                 y_dense_arr[i],
                 marker = 's',
                 markersize = colormap_marker_size,
                 markeredgecolor = 'None',
                 color = vmax_color,
                 alpha = colormap_normalized_data,
                 linestyle = '',
                 )
        # lower spectrum of the colormap
        plt.plot(x_dense_arr[i],
                 y_dense_arr[i],
                 marker = 's',
                 markersize = colormap_marker_size,
                 markeredgecolor = 'None',
                 color = vmin_color,
                 alpha = (1 - colormap_normalized_data),
                 linestyle = '',
                 )
    ax_add = plt.gca()
    ax_add.set(title='')
    ax_add.margins(x=0.0, y=0.0)
    ax_add.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax_add.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax_add.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax_add.yaxis.set_minor_locator(AutoMinorLocator(5))
    colormap_aspect_len = (colormap_vmax - colormap_vmin) * 0.02
    if alignment == 'vertical':        
        plt.xlim(-colormap_aspect_len,colormap_aspect_len)
        plt.ylim(colormap_vmin, colormap_vmax)
        ax_add.yaxis.tick_right()
        plt.xticks([])
    elif alignment == 'horizontal':
        plt.xlim(colormap_vmin, colormap_vmax)
        plt.ylim(-colormap_aspect_len,colormap_aspect_len)
        ax_add.xaxis.tick_bottom()
        plt.yticks([])
    ax_add.set(aspect='equal')
    plt.sca(ax)
##    plt.show()
#############################################################
#axis rotation matrix calculation based on Eulerian angles
#############################################################
def t_euler(euler_angle_type, angle_unit, phi, theta, psi):
    '''
    Description:
        axis rotation matrix calculation based on Eulerian angles
    Args:
        euler_angle_type: string of length 3. It specify the type of rotations based on Eulerian angles. Choices are 'zyz', 'zxz', 'zyx', etc.. Usually the 'zyz' type is used.
            'zyz' : proper Euler angle, y-convention. Performe consecutive rotations at axes counter-clockwisely. z-y-z rotation.
                    First rotate the z axes of atoms by an angle phi, then rotate the intermidiate y axis of atoms by an angle theta, finally rotate the final z axis of atoms by an angle psi
            'zxz' : proper Euler angle, x-convention. Performe consecutive rotations at axes counter-clockwisely. z-x-z rotation.
                    First rotate the z axes of atoms by an angle phi, then rotate the intermidiate x axis of atoms by an angle theta, finally rotate the final z axis of atoms by an angle psi
            'zyx' : Tait-Bryan angles. z-y-x rotation. Performe consecutive rotations at axes counter-clockwisely. z-y-x rotation.
                    First rotate the z axes of atoms by an angle phi, then rotate the intermidiate y axis of atoms by an angle theta, finally rotate the final x axis of atoms by an angle psi
            ......
        angle_unit: String of length 3. If angle_unit = 'deg' then units in degrees, else if angle_unit = 'rad' then units in radians
        phi, theta, psi: the first, second, and third rotation Eulerian angles, units in degrees or radians. The unit depends on the choice of the parameter angle_unit
    Reference:
        Euler Angles: Herbert Goldstein and Charles P. Poole Jr. and John L. Safko,Classical Mechanics (3rd Edition),2001.
    Return:
        Eulerian angle rotation matrix
    '''
    import numpy as np
    if angle_unit == 'deg':
        phi_rad = phi / 180 * np.pi
        theta_rad = theta / 180 * np.pi
        psi_rad = psi / 180 * np.pi
    elif angle_unit == 'rad':
        phi_rad = phi
        theta_rad = theta
        psi_rad = psi      
    rot_axis = [euler_angle_type[0],euler_angle_type[1],euler_angle_type[2]]
    rot_angle = [phi_rad, theta_rad, psi_rad]
    t = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    t_x = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    t_y = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    t_z = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    for i in range(0,3):
        Angle = rot_angle[i]
        axis = rot_axis[i]
        if axis == 'z':
            t_z[0,0] = np.cos(Angle); t_z[0,1] = np.sin(Angle); t_z[0,2] = 0
            t_z[1,0] = -np.sin(Angle); t_z[1,1] = np.cos(Angle);  t_z[1,2] = 0
            t_z[2,0] = 0; t_z[2,1] = 0; t_z[2,2] = 1
            t_rot = t_z
        if axis == 'y':
            t_y[0,0] = np.cos(Angle); t_y[0,1] = 0; t_y[0,2] = -np.sin(Angle)
            t_y[1,0] = 0; t_y[1,1] = 1;  t_y[1,2] = 0
            t_y[2,0] = np.sin(Angle); t_y[2,1] = 0; t_y[2,2] = np.cos(Angle)
            t_rot = t_y
        if axis == 'x':
            t_x[0,0] = 1; t_x[0,1] = 0; t_x[0,2] = 0
            t_x[1,0] = 0; t_x[1,1] = np.cos(Angle); t_x[1,2] = np.sin(Angle)
            t_x[2,0] = 0; t_x[2,1] = -np.sin(Angle); t_x[2,2] = np.cos(Angle)
            t_rot = t_x
        t = np.dot(t_rot,t)
    return t

def lorentzian_lineshape_func(x, x_i, delta):
    '''
    Reference:
        Shan-Ying Wang and Jing-Zhi Yu and Hiroshi Mizuseki and Qiang Sun and Chong-Yu Wang and Yoshiyuki Kawazoe, PRB, 2004, 70
        Peter F. Bernath, Spectra of Atoms and Molecules, Oxford University Press, 2015
    delta: The half width at half maximum (HWHM). For the broadening of the DOS, the delta is the broadening width
    '''
    import numpy as np
    result = delta / np.pi /(np.power((x - x_i),2) + np.power(delta,2))
    return result

def lorentzian_broadening(x_arr, y_arr, delta):
    '''
    Reference:
        Shan-Ying Wang and Jing-Zhi Yu and Hiroshi Mizuseki and Qiang Sun and Chong-Yu Wang and Yoshiyuki Kawazoe, PRB, 2004, 70
        Peter F. Bernath, Spectra of Atoms and Molecules, Oxford University Press, 2015
    delta: The half width at half maximum (HWHM). For the broadening of the DOS, the delta is the broadening width
    x_grid: the spacing of the x points for the broadened data
    the dimension of the x_arr should be the same as that of the y_arr
    for the density fo states: x_arr = energy_arr, y_arr = dos_arr
    If the generated number of x points are larger than the original number of x points and the delta is too small, it may induce fictitious peaks. That is to say, if the original number of x points is too sparse (too sparse points across the spectrum will not make the curve look good), do not expect the broadening scheme to make the curves with sparse points look better. If one wants the sparse curve to look better, the interpolation algorithm can be used instead of directly using broadening scheme. Based on this thought, the generated points of the broadened curve will not exceed the number of points of the original dataset.
    '''
    import numpy as np
    ##x_min = np.min(x_arr)
    ##x_max = np.max(x_arr)
    ##x_grid_num = len(x_arr)
    ##x_arr_broadening = np.linspace(x_min, x_max, x_grid_num)
    x_arr_broadening = x_arr
    y_arr_broadening = np.array([0.0000000] * len(x_arr_broadening), dtype = np.float)
    for i in range(0, len(x_arr_broadening)):
        for j in range(0, len(x_arr)):
            y_arr_broadening[i] = y_arr_broadening[i] + y_arr[j] * lorentzian_lineshape_func(float(x_arr_broadening[i]), float(x_arr[j]), delta)
    return x_arr_broadening, y_arr_broadening

def morse(r, D_e, a, r_e):
    '''
    Morse-like potential
    $V\prime(r)=D_{e}(1-e^{-a(r-r_{e})})^{2}$
    '''
    import numpy as np
    return np.power(D_e * (1 - np.exp(-a * (r - r_e))),2)

def lennard_jones(r, epsilon, r_m):
    '''
    Lennard Jones potential
    $V_{LJ}=\epsilon[(r_m/r)^{12}-(r_m/r)^{6}]$
    '''
    import numpy as np
    return epsilon * (np.power((r_m/r),12)-2*np.power((r_m/r),6))

def uber(a, a_m, l):
    import numpy as np
    '''
    UBER relation
    '''
    astar = (a-a_m)/l
    return np.exp(- np.sqrt(2) * astar) - 2 * np.exp(-astar/np.sqrt(2))

def birch_murnaghan(x_list, energy_list, x_type = 'V', b_0prime = 160):
    '''
    Birch Murnaghan equation of state (EOS)
    $E(V)=E_{0}+\frac{B_{0}V}{B_{0}^{\prime}}(\frac{(V_{0}/V)^{B_{0}^{\prime}}}{B_{0}^{\prime}-1}+1)-\frac{B_{0}V_{0}}{B_{0}^{\prime}-1}$
    x_list: the value of the x axis
    x_type: values can be 'V' or 'a0', which means volume or lattice constant.
    energy_list: energy of the systems
    '''
    pass
    return 0

def e_subst(etot_doped, etot_undoped, composition_str, base_compostion_str, chem_pot_dict):
    '''function of calculating the substitution formation energy
    etot_doped: the model after doping;
    etot_undoped: the model before doping;
    composition_str: a string of the composition of the doped system;
    base_composition_str: a string of the composition of the undoped system;
    chem_pot_dict: a dictionary of the chemical potential of elements.
    '''
    import math
    from . import periodic_table

    periodic_table_dict = periodic_table.periodic_tab()

    base_elmt_list = extract_alpha_str(base_compostion_str)
    elmt_list = extract_alpha_str(composition_str)
    base_elmt_num_list = extract_num(base_compostion_str)
    elmt_num_list = extract_num(composition_str)
    # initialization
    base_elmt_num_dict = {}
    elmt_num_dict = {}
    for i_elmt in periodic_table_dict['symbol'].keys():
        base_elmt_num_dict[i_elmt] = 0
        elmt_num_dict[i_elmt] = 0

    # refresh the dictionary
    for i in range(0, len(base_elmt_list)):
        base_elmt_num_dict[base_elmt_list[i]] = int(base_elmt_num_list[i])
    for i in range(0, len(elmt_list)):
        elmt_num_dict[elmt_list[i]] = int(elmt_num_list[i])

    num_subst_dict = {}
    for i_elmt in periodic_table_dict['symbol'].keys():
        num_subst_dict[i_elmt] = abs(elmt_num_dict[i_elmt] - base_elmt_num_dict[i_elmt])

    solute_elmt_list = list(set(elmt_list) - set(base_elmt_list))

    num_tot = 0
    solute_elmt_dict = {}
    solute_elmt_num_dict = {}
    for i_elmt in solute_elmt_list:
        solute_elmt_dict[i_elmt] = i_elmt
        solute_elmt_num_dict[i_elmt] = elmt_num_dict[i_elmt]
        num_tot += solute_elmt_num_dict[i_elmt]
        
    etot_fin = etot_doped
    etot_ini = etot_undoped
    for i_base_elmt in base_elmt_list:
        etot_fin += num_subst_dict[i_base_elmt] * chem_pot_dict[i_base_elmt]

    for i_solute_elmt in solute_elmt_list:
        etot_ini += num_subst_dict[i_solute_elmt] * chem_pot_dict[i_solute_elmt]

    if num_tot != 0:
        res = (etot_fin - etot_ini) / num_tot
    else:
        res = math.nan

#    try:
#        res = (etot_fin - etot_ini) / num_tot
#    except ZeroDivisionError:
#        res = float('NaN')
#    except RuntimeWarning:
#        # RuntimeWarning: invalid value encountered in double_scalars
#        res = float('NaN')

    return res

# atom numbering according to the order of three orthogonalized directions
def atom_number_ortho(atom_key_arr, pos_arr_cartesian, delta = 0.05, number_order = 'xyz'):
    '''
    - Define atom number without considering axis rotation. In the future, axis rotation can be added.
    - The basic idea is to sorting atom indices independently in three dimensions, then
    grouping the atoms according to their coordinate in each direction (the atoms with 
    similar positions in one direction are grouped in the same group), then
    sorting the atom indices again according to the position and groups.
    - The numbering starts from the origin, i.e. (0,0,0)
    - atom_key_arr: Atom index starts from 0
    - pos_arr_cartesian: It is a two dimensional 'N by 3' array, where N denotes number of atoms, 
    the array columns are x, y, and z Cartesian coordinates.
    - delta: This is the spacial resolution, for example z=0.24 and z=0.25 are considered to be in the same position if delta=0.01.
    - number_order: choices are 'xyz', 'yzx', 'zxy', 'xzy', 'yxz', 'zyx'. Designate the atom number according to the order in each axis directions.
    for example, xyz denotes number starts from 0 in the origin, then the index first increases in x direction, then in y and z.
    '''
    import numpy as np
    n_atoms = len(atom_key_arr)
    # sorted and grouped indices
    indx_arr = np.array([0] * n_atoms * 3)
    indx_arr.shape = n_atoms, 3
    #max_indx_arr denotes the maximum index in each direction
    max_indx_arr = np.array([None] * 3)
    # three directions
    for i in range(0,3):
        running_indx = 0
        sorted_pos_arr_cartesian = np.sort(pos_arr_cartesian[:,i])
        sorted_pos_arr_cartesian_indx = np.argsort(pos_arr_cartesian[:,i])
        temp_pos = np.min(sorted_pos_arr_cartesian)
        for temp_indx in sorted_pos_arr_cartesian_indx:
            if (pos_arr_cartesian[temp_indx, i] - temp_pos) <= delta:
                pass
            else:
                running_indx += 1
            temp_pos = pos_arr_cartesian[temp_indx, i]
            indx_arr[temp_indx, i] = running_indx
        max_indx_arr[i] = np.max(indx_arr[:,i])
    #analyzing the numbering order
    number_order_list = [0, 1, 2]
    max_indx_arr_order = np.array([None] * 3)
    for i in range(0,3):
        axis_direction = number_order.strip()[i]
        if axis_direction == 'x':
            axis_direction_indx = 0
        elif axis_direction == 'y':
            axis_direction_indx = 1
        elif axis_direction == 'z':
            axis_direction_indx = 2
        else:
            print('Error: error in reading the parameter number_order')
        number_order_list[i] = axis_direction_indx
        max_indx_arr_order[i] =  max_indx_arr[number_order_list[i]]
    grid_arr = np.array([None] * n_atoms)
    for temp_indx in range(0, n_atoms):
        grid_arr[temp_indx] = str(indx_arr[temp_indx, number_order_list[0]]) + ',' + str(indx_arr[temp_indx, number_order_list[1]]) + ',' + str(indx_arr[temp_indx, number_order_list[2]])
    #from collections import Counter
    #b = dict(Counter(grid_arr))
    #print ([key for key,value in b.items() if value > 1])
    #print ({key:value for key,value in b.items()if value > 1})
    atom_number_ortho_arr = np.array([None] * n_atoms)
    running_number = 1
    for k in range(0, max_indx_arr_order[2]+1):
        for j in range(0, max_indx_arr_order[1]+1):
            for i in range(0, max_indx_arr_order[0]+1):
                temp_grid_str = str(i) + ',' + str(j) + ',' + str(k)
                temp_indx = np.argwhere(grid_arr == temp_grid_str).squeeze()
                #if np.isfinite(temp_indx):
                if temp_indx.size > 0:
                    atom_number_ortho_arr[int(temp_indx)] = running_number
                    #print(pos_arr_cartesian[int(temp_indx),:])
                    running_number += 1
    return atom_number_ortho_arr

######################
#Notes
######################
'''
At the end of this file, there a need to take the following notes:
    1. NaN value: Because the type of nan in python is float, it would not cause too much errors when dealing with floats. So we use nan instead of None, although sometimes the ``None'' value would be better.
'''
