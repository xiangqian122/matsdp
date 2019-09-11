##########################################
#Folder and file manipulation
##########################################
def cp(src_filedir,dst_filedir):
    import os
    import shutil
    if not os.path.isfile(src_filedir):
        print(src_filedir + ' does not exist!')
    else:
        fpath,fname = os.path.split(dst_filedir)
        if not os.path.exists(fpath):
            os.makedirs(fpath)
        shutil.copyfile(src_filedir,dst_filedir)
        
def file_path_name(file_dir):
    '''Get parent path and file name for the input file directory'''
    import os
    
    parent_path = os.path.dirname(file_dir)  
    filename = os.path.split(file_dir)[-1]
    return parent_path, filename

def grep(kwd,file_dir):
    '''extract the line with the keyword(kwd) specified
    the function the is similar to the "grep" command in linux shell'''
    line_list = []
    with open(file_dir,'r') as f:
        lines = f.readlines()
        for line in lines:
            if kwd in line:
                line_list.append(line)
    line_str = ''.join(line_list)
    return line_str

def line_num(input_file):
    'Get number of lines from a file'
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
    'Merge Files'
    import os
    cwd = os.getcwd()
    file1_dir = cwd + '/' + file1
    file2_dir = cwd + '/' + file2
    with open(file1_dir,'a') as file1_obj:
        with open(file2_dir) as file2_obj:
            file2line = file2_obj.readlines()
            file1_obj.write("".join(file2line))

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

def write_log(logfile,output_str):
    '''Write specific string into a log file,
       When the log file exists, then append output to log file.
       If the log file doesn't exist, then create it and write output to log file.
       Note that output_str must be string format'''
    import os
    if os.path.isfile(logfile) == False:
        with open(logfile,'w') as logfile_object:
            logfile_object.write(output_str +  '\n')
    else:
        with open(logfile,'a') as logfile_object:
            logfile_object.write(output_str + '\n')

##########################################
#String manipulation
##########################################
def extract_num(Str):
    '''Extract number from string using regular expression'''
    import re
    Num = re.findall(r"\d+\.?\d*", Str)
    return Num

def insert_str(old_str,indx,inserted_str):
    '''Insert a specified string(IndsetedStr) into sepecified location of an input string(old_str)
    indx is where you want to insert your string
    inserted_str is the string that you want to insert'''
    new_str = old_str[:indx] + inserted_str + old_str[indx:]
    return new_str

def replace_file_content(filename,oldStr,newStr):
    'Replace string in file'
    import os
    cwd = os.getcwd()
    fileObj = cwd + '/' + filename
    f = open(fileObj,'r+')
    all_the_lines=f.readlines()
    f.seek(0)
    f.truncate()
    for line in all_the_lines:
        f.write(line.replace(oldStr,newStr))
    f.close()
  
def split_line(line, separator = ' '):
    '''split the line using white space 
    and put the splitted line string into a list'''
    temp_list = [item.strip() for item in str(line).replace('\t', '    ').strip('\n').strip().split(separator)]
    while '' in temp_list:
        temp_list.remove('')
    return temp_list

##########################################
#Other functions
##########################################
def logic_retn_val(input_logic_val,true_return_val,false_return_val):
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


