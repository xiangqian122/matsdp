# For any inofrmation please contact dianwuwang@163.com

import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")    
import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.filedialog
import tkinter.messagebox
from matsdp import default_params
from matsdp import funcs
from matsdp import periodic_table
from matsdp.vasp import vasp_read
from matsdp.vasp import vasp_build
from matsdp.vasp import vasp_plot
from matsdp.vasp import vasp_analyze
from matsdp.vasp import vasp_write
from matsdp.apt import apt_plot

defaults_dict = default_params.default_params()
logfile = defaults_dict['logfile']

periodic_table_dict = periodic_table.periodic_tab()
dos_mode = periodic_table_dict['dos_mode']
elmt_color = periodic_table_dict['elmt_color']

# set the maximum recursion depth
sys.setrecursionlimit(100000)

time_start = time.time()
formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time_start))
# log file for the program. The results will be dumped into this log file
funcs.write_log(logfile,'## Job started at ' + formatted_time)

program_name = 'matsdp'
version = '0.1.5'
authors = 'dianwuwang@163.com'

root = tk.Tk()
root.geometry('500x500+5+5')
root.title(program_name)

def on_run_nn_map_button_clicked(event):
    '''
    Description:
        Run nn_map()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## nn_map GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='nn_map')
    child_toplevel.geometry('1030x130+20+25')
    
    nn_map_subtoolkit(child_toplevel)
    return 'break'

def on_run_simple_cna_button_clicked(event):
    '''
    Description:
        simple_cna()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## simple_cna GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='simple_cna')
    child_toplevel.geometry('1030x130+20+25')
    
    simple_cna_subtoolkit(child_toplevel)
    return 'break'

def on_run_subst_button_clicked(event):
    '''
    Description:
        substitution()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## substitution GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='substitution')
    child_toplevel.geometry('1020x130+20+25')
    
    subst_subtoolkit(child_toplevel)
    return 'break'

def on_run_vasp_selection_sphere_button_clicked(event):
    '''
    Description:
        selection_sphere()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## selection_sphere GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='selection_sphere')
    child_toplevel.geometry('1020x130+20+25')
    
    vasp_selection_sphere_subtoolkit(child_toplevel)
    return 'break'

def on_run_plot_dos_button_clicked(event):
    '''
    Description:
        plot_dos()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## plot_dos GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='plot_dos')
    child_toplevel.geometry('1225x600+25+25')
    
    plot_dos_subtoolkit(child_toplevel)
    return 'break'

def on_run_plot_poscar_button_clicked(event):
    '''
    Description:
        plot_poscar()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## plot_poscar GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='plot_poscar')
    child_toplevel.geometry('1080x590+25+25')
    
    plot_poscar_subtoolkit(child_toplevel)
    return 'break'

def on_run_estruct_button_clicked(event):
    '''
    Description:
        estruct()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## estruct GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='estruct')
    child_toplevel.geometry('1020x130+20+25')
    
    estruct_subtoolkit(child_toplevel)
    return 'break'

def on_run_write_poscar_with_force_button_clicked(event):
    '''
    Description:
        write_poscar_with_force()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## write_poscar_with_force GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='write_poscar_with_force')
    child_toplevel.geometry('1020x130+20+25')
    
    write_poscar_with_force_subtoolkit(child_toplevel)
    return 'break'

def on_run_concentration_profile_button_clicked(event):
    '''
    Description:
        plot_proxigram_csv()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## plot concentration profile GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='plot concentration profile')
    child_toplevel.geometry('1190x150+20+25')
    
    concentration_profile_subtoolkit(child_toplevel)
    return 'break'


def on_run_read_py_script_button_clicked(event):
    '''
    Description:
        Run ReadPyScript()
    '''
    import tkinter as tk
    funcs.write_log(logfile,'## read_py_script GUI\n')
    child_toplevel = tk.Toplevel(root)
    tk.Label(child_toplevel, text='read_py_script')
    child_toplevel.geometry('1020x110+20+25')
    
    read_py_script_subtoolkit(child_toplevel)
    return 'break'

##############################
#Read .py script
##############################
class read_py_script_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'read .py or *.log file'
        self.max_num_of_files = 1
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.current_file_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.init_all_input_params()
        self.init_gui()
    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                }
                               for k in range(self.max_num_of_files)]
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)

    def on_run_button_clicked(self):
        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        pyfile = self.get_input_dict()['list_of_file_paths'][0]
        
        if os.path.isfile(pyfile) == True:
            with open(pyfile,'r') as f:
                lines = f.read()
                quote_str1 = '''
                ######################################
                # *.py File
                ######################################
                '''
                quote_str2 = '''
                ######################################
                #Execution result of the *.py file
                ######################################
                '''
                funcs.write_log(logfile,quote_str1)
                funcs.write_log(logfile,lines)
                funcs.write_log(logfile,quote_str2)
                exec(lines)
        else:
            funcs.write_log(logfile, input_settingsfile + " doesn't exist, please check current folder")
        funcs.write_log(logfile,'## read_py_script complete\n')
        return 'break'
    
    def create_run_bar(self):
        run_bar_frame = tk.Frame(self.root, height = 15)
        start_row = self.max_num_of_files + 10
        run_bar_frame.grid(row = start_row, columnspan = 13, sticky = tk.W + tk.E, padx = 15, pady = 10)
        self.run_button = tk.Button(run_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 1, padx = 2)
        self.run_button.configure(background='lightgreen')
    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root)
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)
        for i in range(0, self.max_num_of_files):
            if i == 0:
                tk.Label(left_frame, text = '*.py or *.log file:').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text = 'Find', command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 120, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = i, column =4, padx = 6, pady =4)
            if i == 0:
                filename_str.set('Enter or choose the *.py or *.log file (the file needs to be in the same directory as the current executable).')

    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_run_bar()

##########################
#nn_map
##########################
class nn_map_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'nn_map'
        self.max_num_of_n_shell = 10
        self.max_num_of_files = 1
        self.initial_n_shell = 1
        self.initial_a0 = 3.545
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.current_file_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.init_all_input_params()
        self.init_gui()
    def set_n_shell(self):
        self.get_input_dict()['n_shell'] = int(self.n_shell_value_widget.get())
    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                'a0' : self.initial_a0,
                                'n_shell' : self.initial_n_shell,
                                }
                               for k in range(self.max_num_of_files)]
    def on_n_shell_changed(self):
        self.set_n_shell()
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)

    def on_run_button_clicked(self):
        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        self.get_input_dict()['n_shell'] = int(self.n_shell_value_widget.get())
        self.get_input_dict()['a0'] = float(self.a0_value_widget.get())
        
        poscar_dir = self.get_input_dict()['list_of_file_paths'][0]
        n_shell = self.get_input_dict()['n_shell']
        a0 = self.get_input_dict()['a0']
        
        vasp_analyze.nn_map(poscar_dir, a0, n_shell)
        return 'break'
    
    def create_run_bar(self):
        run_bar_frame = tk.Frame(self.root, height = 15)
        start_row = self.max_num_of_files + 10
        run_bar_frame.grid(row = start_row, columnspan = 13, sticky = tk.W + tk.E, padx = 15, pady = 10)
        self.run_button = tk.Button(run_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 1, padx = 2)
        self.run_button.configure(background='lightgreen')
    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root)
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)
        for i in range(0, self.max_num_of_files):
            tk.Label(left_frame, text = 'POSCAR:').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text='Find',command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 120, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = i, column =4, padx = 6, pady =4)
            filename_str.set('Enter or choose POSCAR path here')
    def create_top_bar(self):
        top_bar_frame = tk.Frame(self.root, height = 25)
        top_bar_frame.grid(row = 0, columnspan = 12, rowspan = 10, padx = 5, pady = 5)

        tk.Label(top_bar_frame, text = 'n_shell=').grid(row = 0, column = 0, sticky = tk.E)
        self.n_shell_value_widget = tk.Spinbox(top_bar_frame, from_ = 1, to = self.max_num_of_n_shell, width = 5, command = self.on_n_shell_changed)
        self.n_shell_value_widget.grid(row = 0, column = 1)
        self.n_shell_value_widget.delete(0,'end')
        self.n_shell_value_widget.insert(0,self.initial_n_shell)

        tk.Label(top_bar_frame, text = 'a0=').grid(row = 0, column =2)
        latt_const_str = tk.StringVar()
        self.a0_value_widget = tk.Entry(top_bar_frame, textvariable = latt_const_str)
        self.a0_value_widget.grid(row = 0, column = 3, padx = 3, pady = 2)
        latt_const_str.set(self.get_input_dict()['a0'])

    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_top_bar()
        self.create_run_bar()

##if __name__ == '__main__':
##    print('main')
##    root = tk.Tk()
##    nn_mapToolkit(root)
##    root.mainloop()

##########################
# simple_cna
##########################
class simple_cna_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'simple_cna'
        self.max_num_of_files = 1
        self.initial_a0 = 3.545
        self.initial_common_neighbor_elmt_list = 'Re, Ni'
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.current_file_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.init_all_input_params()
        self.init_gui()
    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                'a0' : self.initial_a0,
                                'common_neighbor_elmt_list' : self.initial_common_neighbor_elmt_list
                                }
                               for k in range(self.max_num_of_files)]
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)

    def on_run_button_clicked(self):
        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        self.get_input_dict()['a0'] = self.a0_value_widget.get()
        self.get_input_dict()['common_neighbor_elmt_list'] = self.common_neighbor_elmt_list_value_widget.get()
        
        poscar_dir = self.get_input_dict()['list_of_file_paths'][0]
        a0 = float(self.get_input_dict()['a0'])
        if ',' in self.get_input_dict()['common_neighbor_elmt_list']:
            common_neighbor_elmt_list = funcs.split_line(line = self.get_input_dict()['common_neighbor_elmt_list'], separator = ',')
        else:
            common_neighbor_elmt_list = funcs.split_line(line = self.get_input_dict()['common_neighbor_elmt_list'], separator = ' ')
        
        vasp_analyze.simple_cna(poscar_dir, a0, common_neighbor_elmt_list)
        return 'break'
    
    def create_run_bar(self):
        run_bar_frame = tk.Frame(self.root, height = 15)
        start_row = self.max_num_of_files + 10
        run_bar_frame.grid(row = start_row, columnspan = 13, sticky = tk.W + tk.E, padx = 15, pady = 10)
        self.run_button = tk.Button(run_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 1, padx = 2)
        self.run_button.configure(background='lightgreen')
    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root)
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)
        for i in range(0, self.max_num_of_files):
            tk.Label(left_frame, text = 'POSCAR:').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text='Find',command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 120, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = i, column =4, padx = 6, pady =4)
            filename_str.set('Enter or choose POSCAR path here')
    def create_top_bar(self):
        top_bar_frame = tk.Frame(self.root, height = 25)
        top_bar_frame.grid(row = 0, columnspan = 12, rowspan = 10, padx = 5, pady = 5)

        tk.Label(top_bar_frame, text = 'a0=').grid(row = 0, column =0)
        a0_str = tk.StringVar()
        self.a0_value_widget = tk.Entry(top_bar_frame, width = 7, textvariable = a0_str)
        self.a0_value_widget.grid(row = 0, column = 1, padx = 3, pady = 2)
        a0_str.set(self.get_input_dict()['a0'])

        tk.Label(top_bar_frame, text = 'common_neighbor_elmt_list=').grid(row = 0, column =2)
        common_neighbor_elmt_list_str = tk.StringVar()
        self.common_neighbor_elmt_list_value_widget = tk.Entry(top_bar_frame, width = 37, textvariable = common_neighbor_elmt_list_str)
        self.common_neighbor_elmt_list_value_widget.grid(row = 0, column = 3, padx = 3, pady = 2)
        common_neighbor_elmt_list_str.set(self.get_input_dict()['common_neighbor_elmt_list'])

    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_top_bar()
        self.create_run_bar()

###############################
#vasp_build -- substitution
###############################
class subst_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'substitution'
        self.max_num_of_files = 2
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.current_file_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.init_all_input_params()
        self.init_gui()
    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                }
                               for k in range(self.max_num_of_files)]
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)

    def on_run_button_clicked(self):
        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        poscar_dir = self.get_input_dict()['list_of_file_paths'][0]
        substitution_list_file = self.get_input_dict()['list_of_file_paths'][1]
        
        vasp_build.substitution(substitution_list_file, poscar_dir)
        funcs.write_log(logfile,'## Substitution complete\n')
        return 'break'
    
    def create_run_bar(self):
        run_bar_frame = tk.Frame(self.root, height = 15)
        start_row = self.max_num_of_files + 10
        run_bar_frame.grid(row = start_row, columnspan = 13, sticky = tk.W + tk.E, padx = 15, pady = 10)
        self.run_button = tk.Button(run_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 1, padx = 2)
        self.run_button.configure(background='lightgreen')
    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root)
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)
        for i in range(0, self.max_num_of_files):
            if i == 0:
                tk.Label(left_frame, text = 'POSCAR:').grid(row = i, column = 0, padx =2, pady = 4)
            elif i == 1:
                tk.Label(left_frame, text = '.subst file:').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text = 'Find', command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 120, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = i, column =4, padx = 6, pady =4)
            if i == 0:
                filename_str.set('Enter or choose POSCAR path here. The atoms in the POSCAR will be substituted for according to the .subst file')
            elif i == 1:
                filename_str.set('Enter or choose a .subst file here')

    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_run_bar()

##########################
# vasp_selection_sphere
##########################
class vasp_selection_sphere_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'vasp_selection_sphere'
        self.max_num_of_files = 1
        self.initial_origin_atom_name = 'Ni1'
        self.initial_radius = 7.0
        self.initial_output_file_name = 'example'

        self.include_mirror_atoms = False
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.current_file_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.init_all_input_params()
        self.init_gui()


    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                'origin_atom_name' : self.initial_origin_atom_name,
                                'radius' : self.initial_radius,
                                'output_file_name' : self.initial_output_file_name
                                }
                               for k in range(self.max_num_of_files)]
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)

    def on_include_mirror_atoms_changed(self,value):
        self.include_mirror_atoms = value

    def on_run_button_clicked(self):
        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        self.get_input_dict()['origin_atom_name'] = self.origin_atom_name_value_widget.get()
        self.get_input_dict()['radius'] = self.radius_value_widget.get()
        self.get_input_dict()['output_file_name'] = self.output_file_name_value_widget.get()
        
        poscar_dir = self.get_input_dict()['list_of_file_paths'][0]
        origin_atom_name = str(self.get_input_dict()['origin_atom_name'])
        radius = float(self.get_input_dict()['radius'])
        self.include_mirror_atoms = funcs.convert_bool(self.include_mirror_atoms_str.get())
        include_mirror_atoms = self.include_mirror_atoms
        output_file_name = str(self.get_input_dict()['output_file_name'])
        
        vasp_build.selection_sphere(poscar_dir, origin_atom_name, radius, include_mirror_atoms, output_file_name)
        return 'break'
    
    def create_run_bar(self):
        run_bar_frame = tk.Frame(self.root, height = 15)
        start_row = self.max_num_of_files + 10
        run_bar_frame.grid(row = start_row, columnspan = 13, sticky = tk.W + tk.E, padx = 15, pady = 10)
        self.run_button = tk.Button(run_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 1, padx = 2)
        self.run_button.configure(background='lightgreen')
    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root)
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)
        for i in range(0, self.max_num_of_files):
            tk.Label(left_frame, text = 'POSCAR:').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text='Find',command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 120, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = i, column =4, padx = 6, pady =4)
            filename_str.set('Enter or choose POSCAR path here')
    def create_top_bar(self):
        top_bar_frame = tk.Frame(self.root, height = 25)
        top_bar_frame.grid(row = 0, columnspan = 12, rowspan = 10, padx = 5, pady = 5)

        tk.Label(top_bar_frame, text = 'origin_atom_name=').grid(row = 0, column =0)
        origin_atom_name_str = tk.StringVar()
        self.origin_atom_name_value_widget = tk.Entry(top_bar_frame, width = 7, textvariable = origin_atom_name_str)
        self.origin_atom_name_value_widget.grid(row = 0, column = 1, padx = 3, pady = 2)
        origin_atom_name_str.set(self.get_input_dict()['origin_atom_name'])

        tk.Label(top_bar_frame, text = 'radius=').grid(row = 0, column =2)
        radius_str = tk.StringVar()
        self.radius_value_widget = tk.Entry(top_bar_frame, width = 8, textvariable = radius_str)
        self.radius_value_widget.grid(row = 0, column = 3, padx = 3, pady = 2)
        radius_str.set(self.get_input_dict()['radius'])

        tk.Label(top_bar_frame,text = 'include_mirror_atoms:').grid(row = 0, column = 4, padx =2 , pady = 4)
        self.include_mirror_atoms_str = tk.StringVar()
        options = ['True','False']
        self.include_mirror_atoms_option_menu_widget = tk.OptionMenu(top_bar_frame,self.include_mirror_atoms_str,*options, command = self.on_include_mirror_atoms_changed)
        self.include_mirror_atoms_option_menu_widget.grid(row = 0, column =5, padx = 6, pady =4, sticky = tk.W + tk.E)
        self.include_mirror_atoms_str.set(str(self.include_mirror_atoms))

        tk.Label(top_bar_frame, text = 'output_file_name=').grid(row = 0, column =6)
        output_file_name_str = tk.StringVar()
        self.output_file_name_value_widget = tk.Entry(top_bar_frame, width = 45, textvariable = output_file_name_str)
        self.output_file_name_value_widget.grid(row = 0, column = 7, padx = 3, pady = 2)
        output_file_name_str.set(self.get_input_dict()['output_file_name'])

    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_top_bar()
        self.create_run_bar()

###################
#plot_dos
###################
class plot_dos_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'plot_dos'
        self.max_num_of_doscar = 10
        self.max_num_of_atoms = 10
        self.max_num_of_subplots = 20
        self.max_num_of_elmts = 10
        self.initial_num_of_doscar = 1
        self.initial_num_of_atoms = 1
        self.initial_num_of_subplots = 1
        self.initial_num_of_elmts = 1
        self.num_of_doscar = self.initial_num_of_doscar
        self.num_of_atoms = self.initial_num_of_atoms
        self.num_of_subplots = self.initial_num_of_subplots
        self.num_of_elmts = self.initial_num_of_elmts

        #General
        self.fermi_shift_zero = True
        self.fig_format = 'png'
        self.fig_format_str = 'png'
        self.fig_width = 15
        self.fig_height = 10
        self.fig_dpi = 600
        self.mainplot_xlabel = False
        self.mainplot_ylabel = False
        self.user_defined_dos_mode = False
        self.initial_dos_mode = dos_mode
        self.dos_mode = dos_mode
        self.peak_analyzer = False
        self.subplot_share_xaxis = False
        self.subplot_share_yaxis = False
                
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_doscar
        self.current_figindx = 0
        #DOSCAR
        self.load_doscar_entry_widget = [None] * self.max_num_of_doscar
        #Atoms        
        self.atom_doscar_indx_widget = [None] * self.max_num_of_atoms
        self.atom_sysname_entry_widget = [None] * self.max_num_of_atoms
        self.atom_name_entry_widget = [None] * self.max_num_of_atoms
        self.atom_palette_entry_widget = [None] * self.max_num_of_atoms     
        self.atom_subplot_arg_entry_widget = [None] * self.max_num_of_atoms
        #Subplots
        self.subplot_arg_entry_widget = [None] * self.max_num_of_subplots
        self.xlo_entry_widget = [None] * self.max_num_of_subplots
        self.xhi_entry_widget = [None] * self.max_num_of_subplots
        self.ylo_entry_widget = [None] * self.max_num_of_subplots
        self.yhi_entry_widget = [None] * self.max_num_of_subplots
        self.subplot_xtick_check_button_widget = [None] * self.max_num_of_subplots
        self.subplot_ytick_check_button_widget = [None] * self.max_num_of_subplots
        self.subplot_xlabel_check_button_widget = [None] * self.max_num_of_subplots
        self.subplot_ylabel_check_button_widget = [None] * self.max_num_of_subplots
        #Elements: dos_mode
        self.dos_mode_elmtname_entry_widget = [None] * self.max_num_of_elmts
        self.s_dos_mode_check_button_widget = [None] * self.max_num_of_elmts
        self.p_dos_mode_check_button_widget = [None] * self.max_num_of_elmts
        self.d_dos_mode_check_button_widget = [None] * self.max_num_of_elmts
        self.l_dos_mode_check_button_widget = [None] * self.max_num_of_elmts
        self.s_dos_mode = True
        self.p_dos_mode = True
        self.d_dos_mode = True
        self.l_dos_mode = False

        self.init_all_input_params()
        self.init_gui()

    def set_doscar_related_params(self):
        for i in range(0, self.num_of_doscar):
            self.get_doscar_related_input_dict()['ListOfDOSCARPaths'][i] = self.load_doscar_entry_widget[i].get()
        return self.all_doscar_related_params[self.current_figindx]
    def set_atom_related_params(self):
        for i in range(0, self.num_of_atoms):
            self.get_atom_related_input_dict()['ListOfAtomDOSCARindx'][i] = self.atom_doscar_indx_widget[i].get()
            self.get_atom_related_input_dict()['ListOfAtomsysname'][i] = self.atom_sysname_entry_widget[i].get()
            self.get_atom_related_input_dict()['ListOfatom_name'][i] = self.atom_name_entry_widget[i].get()
            self.get_atom_related_input_dict()['ListOfAtomPalette'][i] = self.atom_palette_entry_widget[i].get()
            self.get_atom_related_input_dict()['ListOfAtomSubplotArg'][i] = self.atom_subplot_arg_entry_widget[i].get()
        return self.all_atom_related_params[self.current_figindx]
    def set_subplot_related_params(self):
        for i in range(0, self.num_of_subplots):
            self.get_subplot_related_input_dict()['ListOfSubplotArg'][i] = self.subplot_arg_entry_widget[i].get()
            self.get_subplot_related_input_dict()['ListOfXLo'][i] = self.xlo_entry_widget[i].get()
            self.get_subplot_related_input_dict()['ListOfXHi'][i] = self.xhi_entry_widget[i].get()            
            self.get_subplot_related_input_dict()['ListOfYLo'][i] = self.ylo_entry_widget[i].get()
            self.get_subplot_related_input_dict()['ListOfYHi'][i] = self.yhi_entry_widget[i].get()            
            self.get_subplot_related_input_dict()['ListOfSubplotXTick'][i] = self.subplot_xlabel_check_button_widget[i].var.get()
            self.get_subplot_related_input_dict()['ListOfSubplotYTick'][i] = self.subplot_ylabel_check_button_widget[i].var.get()
            self.get_subplot_related_input_dict()['ListOfSubplotXLabel'][i] = self.subplot_xlabel_check_button_widget[i].var.get()
            self.get_subplot_related_input_dict()['ListOfSubplotYLabel'][i] = self.subplot_ylabel_check_button_widget[i].var.get()
        return self.all_subplot_related_params[self.current_figindx]

    def set_num_of_doscar(self):
        self.num_of_doscar = int(self.num_of_doscar_widget.get())
    def set_num_of_atoms(self):
        self.num_of_atoms = int(self.num_of_atoms_widget.get())
    def set_num_of_subplots(self):
        self.num_of_subplots = int(self.num_of_subplots_widget.get())
    def set_num_of_elmts(self):
        self.num_of_elmts = int(self.num_of_elmts_widget.get())
    def get_doscar_related_input_dict(self):
        return self.all_doscar_related_params[self.current_figindx]
    def get_atom_related_input_dict(self):
        return self.all_atom_related_params[self.current_figindx]
    def get_subplot_related_input_dict(self):
        return self.all_subplot_related_params[self.current_figindx]
    def get_list_of_file_paths(self):
        return self.get_doscar_related_input_dict()['ListOfDOSCARPaths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_doscar_related_params = [{'ListOfDOSCARPaths': [None] * self.max_num_of_doscar,
                                        }
                                       for k in range(self.max_num_of_doscar)]
        self.all_atom_related_params = [{'ListOfAtomDOSCARindx': [1] * self.max_num_of_atoms,
                                      'ListOfAtomsysname': ['None'] * self.max_num_of_atoms,
                                      'ListOfatom_name': ['Ni1'] * self.max_num_of_atoms,
                                      'ListOfAtomPalette': ['black'] * self.max_num_of_atoms,
                                      'ListOfAtomSubplotArg': [111]* self.max_num_of_atoms,
                                      }
                                     for k in range(self.max_num_of_atoms)]
        self.all_subplot_related_params = [{'ListOfSubplotArg': [111]* self.max_num_of_subplots,
                                         'ListOfXLo': ['None'] * self.max_num_of_subplots,
                                         'ListOfXHi': ['None'] * self.max_num_of_subplots,
                                         'ListOfYLo': ['None'] * self.max_num_of_subplots,
                                         'ListOfYHi': ['None'] * self.max_num_of_subplots,
                                         'ListOfSubplotXTick': [True]* self.max_num_of_subplots,
                                         'ListOfSubplotYTick': [True]* self.max_num_of_subplots,
                                         'ListOfSubplotXLabel': [True]* self.max_num_of_subplots,
                                         'ListOfSubplotYLabel': [True]* self.max_num_of_subplots,
                                         }
                                        for k in range(self.max_num_of_subplots)]
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def on_num_of_doscar_changed(self):
        self.set_doscar_related_params()
        self.set_num_of_doscar()
        self.create_left_file_loader()
    def on_num_of_atoms_changed(self):
        self.set_atom_related_params()
        self.set_num_of_atoms()
        self.create_atom_params_panel()
    def on_num_of_subplots_changed(self):
        self.set_subplot_related_params()
        self.set_num_of_subplots()
        self.create_subplot_bar()
    def on_num_of_elmts_changed(self):
        self.set_num_of_elmts()
        self.create_elmt_check_bar()
            
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_doscar_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_doscar_entry_widget[text_widget_indx].insert(0, file_path)

    def on_run_button_clicked(self):
        #DOSCAR
        for i in range(0, self.num_of_doscar):
            self.get_doscar_related_input_dict()['ListOfDOSCARPaths'][i] = self.load_doscar_entry_widget[i].get()
        #Atoms
        for i in range(0, self.num_of_atoms):
            self.get_atom_related_input_dict()['ListOfAtomDOSCARindx'][i] = self.atom_doscar_indx_widget[i].get()
            self.get_atom_related_input_dict()['ListOfAtomsysname'][i] = self.atom_sysname_entry_widget[i].get()
            self.get_atom_related_input_dict()['ListOfatom_name'][i] = self.atom_name_entry_widget[i].get()
            self.get_atom_related_input_dict()['ListOfAtomPalette'][i] = self.atom_palette_entry_widget[i].get()
            self.get_atom_related_input_dict()['ListOfAtomSubplotArg'][i] = self.atom_subplot_arg_entry_widget[i].get()
        #Subplots
        for i in range(0, self.num_of_subplots):
            self.get_subplot_related_input_dict()['ListOfSubplotArg'][i] = self.subplot_arg_entry_widget[i].get()
            self.get_subplot_related_input_dict()['ListOfXLo'][i] = self.xlo_entry_widget[i].get()
            self.get_subplot_related_input_dict()['ListOfXHi'][i] = self.xhi_entry_widget[i].get()
            self.get_subplot_related_input_dict()['ListOfYLo'][i] = self.ylo_entry_widget[i].get()
            self.get_subplot_related_input_dict()['ListOfYHi'][i] = self.yhi_entry_widget[i].get()
            self.get_subplot_related_input_dict()['ListOfSubplotXTick'][i] = self.subplot_xtick_check_button_widget[i].var.get()
            self.get_subplot_related_input_dict()['ListOfSubplotYTick'][i] = self.subplot_ytick_check_button_widget[i].var.get()
            self.get_subplot_related_input_dict()['ListOfSubplotXLabel'][i] = self.subplot_xlabel_check_button_widget[i].var.get()
            self.get_subplot_related_input_dict()['ListOfSubplotYLabel'][i] = self.subplot_ylabel_check_button_widget[i].var.get()
        self.fermi_shift_zero = self.fermi_shift_zero_check_button_widget.var.get()
        self.peak_analyzer = self.peak_analyzer_check_button_widget.var.get()
        self.subplot_share_xaxis = self.subplot_share_xaxis_check_button_widget.var.get()
        self.subplot_share_yaxis = self.subplot_share_yaxis_check_button_widget.var.get()
        self.mainplot_xlabel = self.mainplot_xlabel_check_button_widget.var.get()
        self.mainplot_ylabel = self.mainplot_ylabel_check_button_widget.var.get()
        self.fig_format = self.fig_format_str.get()
        self.fig_width = self.fig_width_entry_widget.get()
        self.fig_height = self.fig_height_entry_widget.get()
        self.fig_dpi = int(self.fig_dpi_entry_widget.get())
        #DOS elements
        self.num_of_elmts = int(self.num_of_elmts_widget.get())
        self.dos_mode = self.initial_dos_mode.copy()
        for i in range(0, self.num_of_elmts):
            self.dos_mode[self.dos_mode_elmtname_entry_widget[i].get()] = []
            if self.s_dos_mode_check_button_widget[i].var.get() == True:
                self.dos_mode[self.dos_mode_elmtname_entry_widget[i].get()].append('s')
            if self.p_dos_mode_check_button_widget[i].var.get() == True:
                self.dos_mode[self.dos_mode_elmtname_entry_widget[i].get()].append('p')                
            if self.d_dos_mode_check_button_widget[i].var.get() == True:
                self.dos_mode[self.dos_mode_elmtname_entry_widget[i].get()].append('d')  
            if self.l_dos_mode_check_button_widget[i].var.get() == True:
                self.dos_mode[self.dos_mode_elmtname_entry_widget[i].get()].append('LDOS')  

        atom_doscar_dir_list = []
        atom_subplot_arg_list = []
        subplot_arg_list = []
        subplot_xtick_list = []
        subplot_ytick_list = []
        subplot_xlabel_list = []
        subplot_ylabel_list = []
        for i in range(0, self.num_of_atoms):
            atom_doscar_dir_list.append(self.get_doscar_related_input_dict()['ListOfDOSCARPaths'][int(self.get_atom_related_input_dict()['ListOfAtomDOSCARindx'][i]) - 1])
            atom_subplot_arg_list.append(int(self.get_atom_related_input_dict()['ListOfAtomSubplotArg'][i]))
        for i in range(0, self.num_of_subplots):
            subplot_arg_list.append(int(self.get_subplot_related_input_dict()['ListOfSubplotArg'][i]))
            subplot_xtick_list.append(self.get_subplot_related_input_dict()['ListOfSubplotXTick'][i])
            subplot_ytick_list.append(self.get_subplot_related_input_dict()['ListOfSubplotYTick'][i])
            subplot_xlabel_list.append(self.get_subplot_related_input_dict()['ListOfSubplotXLabel'][i])
            subplot_ylabel_list.append(self.get_subplot_related_input_dict()['ListOfSubplotYLabel'][i])        
        atom_sysname_list = self.get_atom_related_input_dict()['ListOfAtomsysname'][0:self.num_of_atoms]
        atom_indx_list = self.get_atom_related_input_dict()['ListOfatom_name'][0:self.num_of_atoms]
        atom_palette_list = self.get_atom_related_input_dict()['ListOfAtomPalette'][0:self.num_of_atoms]
        subplot_xlo_list = [ii for ii in self.get_subplot_related_input_dict()['ListOfXLo'][0:self.num_of_subplots]]
        subplot_xhi_list = [ii for ii in self.get_subplot_related_input_dict()['ListOfXHi'][0:self.num_of_subplots]]
        subplot_ylo_list = [ii for ii in self.get_subplot_related_input_dict()['ListOfYLo'][0:self.num_of_subplots]]
        subplot_yhi_list = [ii for ii in self.get_subplot_related_input_dict()['ListOfYHi'][0:self.num_of_subplots]]
        for ii in range(0,self.num_of_subplots):
            if subplot_xlo_list[ii] != 'None':
                subplot_xlo_list[ii] = float(subplot_xlo_list[ii])
            else:
                subplot_xlo_list[ii] = None
            if subplot_xhi_list[ii] != 'None':
                subplot_xhi_list[ii] = float(subplot_xhi_list[ii])
            else:
                subplot_xhi_list[ii] = None
            if subplot_ylo_list[ii] != 'None':
                subplot_ylo_list[ii] = float(subplot_ylo_list[ii])
            else:
                subplot_ylo_list[ii] = None
            if subplot_yhi_list[ii] != 'None':
                subplot_yhi_list[ii] = float(subplot_yhi_list[ii])
            else:
                subplot_yhi_list[ii] = None                            
        mainplot_axis_label_list = [self.mainplot_xlabel, self.mainplot_ylabel]
        fermi_shift_zero = self.fermi_shift_zero
        peak_analyzer = self.peak_analyzer
        subplot_share_xy_list = [self.subplot_share_xaxis, self.subplot_share_yaxis]
        fig_format = self.fig_format
        fig_dpi = self.fig_dpi
        fig_size = [float(self.fig_width), float(self.fig_height)]
        dos_mode = self.dos_mode

        vasp_plot.plot_dos(atom_doscar_dir_list, atom_sysname_list, atom_indx_list, atom_palette_list, atom_subplot_arg_list,
                          subplot_arg_list, subplot_xlo_list, subplot_xhi_list, subplot_ylo_list, subplot_yhi_list,
                          subplot_xtick_list, subplot_ytick_list, subplot_xlabel_list, subplot_ylabel_list, subplot_share_xy_list, mainplot_axis_label_list,
                          dos_mode, fermi_shift_zero, peak_analyzer, fig_format, fig_size, fig_dpi)
        return 'break'
    
    def create_bottom_bar(self):
        bottom_bar_frame = tk.Frame(self.root, height = 25,background = 'gray')
        start_row = self.max_num_of_doscar + 1
        bottom_bar_frame.grid(row = start_row, columnspan = 14, sticky = tk.W + tk.E + tk.N + tk.S, padx = 5, pady = 5)

        fermi_ticked_yes = tk.BooleanVar()
        self.fermi_shift_zero_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'fermi_shift_zero', variable = fermi_ticked_yes)
        self.fermi_shift_zero_check_button_widget.grid(row = start_row, column =0, padx = 6, pady =4)
        fermi_ticked_yes.set(bool(self.fermi_shift_zero))
        self.fermi_shift_zero_check_button_widget.var = fermi_ticked_yes

        tk.Label(bottom_bar_frame,text = 'fig_width:').grid(row = start_row, column = 1, padx =2 , pady = 4)
        fig_width_str = tk.StringVar()
        self.fig_width_entry_widget = tk.Entry(bottom_bar_frame, width = 5, textvariable = fig_width_str)
        self.fig_width_entry_widget.grid(row = start_row, column =2, padx = 6, pady =4)
        fig_width_str.set(str(self.fig_width))

        tk.Label(bottom_bar_frame,text = 'fig_height:').grid(row = start_row, column = 3, padx =2 , pady = 4)
        fig_heightStr = tk.StringVar()
        self.fig_height_entry_widget = tk.Entry(bottom_bar_frame, width = 5, textvariable = fig_heightStr)
        self.fig_height_entry_widget.grid(row = start_row, column =4, padx = 6, pady =4)
        fig_heightStr.set(str(self.fig_height))

        tk.Label(bottom_bar_frame,text = 'fig_dpi:').grid(row = start_row, column = 5, padx =2 , pady = 4)
        fig_dpi_str = tk.StringVar()
        self.fig_dpi_entry_widget = tk.Entry(bottom_bar_frame, width = 5, textvariable = fig_dpi_str)
        self.fig_dpi_entry_widget.grid(row = start_row, column =6, padx = 6, pady =4)
        fig_dpi_str.set(str(self.fig_dpi))
        
        tk.Label(bottom_bar_frame,text = 'fig_format:').grid(row = start_row, column = 7, padx =2 , pady = 4)
        self.fig_format_str = tk.StringVar()
        self.fig_format_option_menu_widget = tk.OptionMenu(bottom_bar_frame,self.fig_format_str,'png', 'eps', 'pdf', 'tif', 'tiff', 'jpg', 'jpeg', 'svg', 'svgz', 'pgf', 'ps', 'raw', 'rgba')
        self.fig_format_option_menu_widget.grid(row = start_row, column =8, padx = 6, pady =4)
        self.fig_format_str.set(str(self.fig_format))        
            
        xlabel_ticked_yes = tk.BooleanVar()
        self.mainplot_xlabel_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'mainplot: xlabel', variable = xlabel_ticked_yes)
        self.mainplot_xlabel_check_button_widget.grid(row = start_row, column =9, padx = 6, pady =4)
        xlabel_ticked_yes.set(bool(self.mainplot_xlabel))
        self.mainplot_xlabel_check_button_widget.var = xlabel_ticked_yes

        ylabel_ticked_yes = tk.BooleanVar()
        self.mainplot_ylabel_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'mainplot: ylabel', variable = ylabel_ticked_yes)
        self.mainplot_ylabel_check_button_widget.grid(row = start_row, column =10, padx = 6, pady =4)
        ylabel_ticked_yes.set(bool(self.mainplot_ylabel))
        self.mainplot_ylabel_check_button_widget.var = ylabel_ticked_yes

        peak_analyzer_ticked_yes = tk.BooleanVar()
        self.peak_analyzer_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'peak_analyzer', variable = peak_analyzer_ticked_yes)
        self.peak_analyzer_check_button_widget.grid(row = start_row, column =11, padx = 6, pady =4)
        peak_analyzer_ticked_yes.set(bool(self.peak_analyzer))
        self.peak_analyzer_check_button_widget.var = peak_analyzer_ticked_yes

        subplot_share_xaxis_ticked_yes = tk.BooleanVar()
        self.subplot_share_xaxis_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'Share x axis', variable = subplot_share_xaxis_ticked_yes)
        self.subplot_share_xaxis_check_button_widget.grid(row = start_row, column =12, padx = 6, pady =4)
        subplot_share_xaxis_ticked_yes.set(bool(self.subplot_share_xaxis))
        self.subplot_share_xaxis_check_button_widget.var = subplot_share_xaxis_ticked_yes

        subplot_share_yaxis_ticked_yes = tk.BooleanVar()
        self.subplot_share_yaxis_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'Share y axis', variable = subplot_share_yaxis_ticked_yes)
        self.subplot_share_yaxis_check_button_widget.grid(row = start_row, column =13, padx = 6, pady =4)
        subplot_share_yaxis_ticked_yes.set(bool(self.subplot_share_yaxis))
        self.subplot_share_yaxis_check_button_widget.var = subplot_share_yaxis_ticked_yes

    def create_elmt_check_bar(self):
        elmts_check_bar_frame = tk.Frame(self.root, height= 25, background = 'light gray')
        start_row = 10
        elmts_check_bar_frame.grid(row = start_row, column = 10, columnspan = 4, sticky = tk.W + tk.E + tk.N + tk.S, padx = 5, pady = 5)
        for i in range(0, self.num_of_elmts):
            tk.Label(elmts_check_bar_frame,text = 'ElmtName:').grid(row = start_row + i + 1, column = 1, padx =2 , pady = 4)
            elmtname_str = tk.StringVar()
            self.dos_mode_elmtname_entry_widget[i] = tk.Entry(elmts_check_bar_frame, width = 10, textvariable = elmtname_str)
            self.dos_mode_elmtname_entry_widget[i].grid(row = start_row + i + 1, column = 1, padx = 6, pady =4)
            elmtname_str.set('Ni')

            s_dos_ticked_yes = tk.BooleanVar()
            self.s_dos_mode_check_button_widget[i] = tk.Checkbutton(elmts_check_bar_frame, text = 's', variable = s_dos_ticked_yes)
            self.s_dos_mode_check_button_widget[i].grid(row = start_row + i + 1, column = 2, padx = 6, pady =4)
            s_dos_ticked_yes.set(bool(self.s_dos_mode))
            self.s_dos_mode_check_button_widget[i].var = s_dos_ticked_yes               
                            
            p_dos_ticked_yes = tk.BooleanVar()
            self.p_dos_mode_check_button_widget[i] = tk.Checkbutton(elmts_check_bar_frame, text = 'p', variable = p_dos_ticked_yes)
            self.p_dos_mode_check_button_widget[i].grid(row = start_row + i + 1, column = 3, padx = 6, pady =4)
            p_dos_ticked_yes.set(bool(self.p_dos_mode))
            self.p_dos_mode_check_button_widget[i].var = p_dos_ticked_yes

            d_dos_ticked_yes = tk.BooleanVar()
            self.d_dos_mode_check_button_widget[i] = tk.Checkbutton(elmts_check_bar_frame, text = 'd', variable = d_dos_ticked_yes)
            self.d_dos_mode_check_button_widget[i].grid(row = start_row + i + 1, column = 4, padx = 6, pady =4)
            d_dos_ticked_yes.set(bool(self.d_dos_mode))
            self.d_dos_mode_check_button_widget[i].var = d_dos_ticked_yes 

            l_dos_ticked_yes = tk.BooleanVar()
            self.l_dos_mode_check_button_widget[i] = tk.Checkbutton(elmts_check_bar_frame, text = 'LDOS', variable = l_dos_ticked_yes)
            self.l_dos_mode_check_button_widget[i].grid(row = start_row + i + 1, column = 5, padx = 6, pady =4)
            l_dos_ticked_yes.set(bool(self.l_dos_mode))
            self.l_dos_mode_check_button_widget[i].var = l_dos_ticked_yes 
    
    def create_atom_params_panel(self):
        atom_params_frame = tk.Frame(self.root, background = 'light blue')
        atom_params_frame.grid(row = 10, column = 3, columnspan = 7, sticky = tk.W + tk.N + tk.S +tk.E, padx = 5, pady = 4)
        tk.Label(atom_params_frame, text = 'DOSCAR No.').grid(row = 0, column = 0, padx =2 , pady = 4)
        for i in range(0, self.num_of_atoms):
            tk.Label(atom_params_frame, text = 'DOSCAR No.').grid(row = i, column = 0, padx =2 , pady = 4)
            self.atom_doscar_indx_widget[i] = tk.Spinbox(atom_params_frame, from_ = 1, to = self.num_of_doscar, width = 3)
            self.atom_doscar_indx_widget[i].grid(row = i, column = 1)
            self.atom_doscar_indx_widget[i].delete(0,'end')
            self.atom_doscar_indx_widget[i].insert(0,str(self.get_atom_related_input_dict()['ListOfAtomDOSCARindx'][i]))

            tk.Label(atom_params_frame, text = 'sysname:').grid(row = i, column =2, padx =2, pady =4)
            AtomsysnameStr = tk.StringVar()
            self.atom_sysname_entry_widget[i] = tk.Entry(atom_params_frame, width = 5, textvariable = AtomsysnameStr)
            self.atom_sysname_entry_widget[i].grid(row = i, column =3, padx = 6, pady =4)
            AtomsysnameStr.set(str(self.get_atom_related_input_dict()['ListOfAtomsysname'][i]))
        
            tk.Label(atom_params_frame,text = 'atomname:').grid(row = i, column = 4, padx =2 , pady = 4)
            atom_nameStr = tk.StringVar()
            self.atom_name_entry_widget[i] = tk.Entry(atom_params_frame, width = 6, textvariable = atom_nameStr)
            self.atom_name_entry_widget[i].grid(row = i, column =5, padx = 6, pady =4)
            atom_nameStr.set(str(self.get_atom_related_input_dict()['ListOfatom_name'][i]))

            tk.Label(atom_params_frame,text = 'color:').grid(row = i, column = 6, padx =2 , pady = 4)
            atom_palette_str = tk.StringVar()
            self.atom_palette_entry_widget[i] = tk.Entry(atom_params_frame, width = 6, textvariable = atom_palette_str)
            self.atom_palette_entry_widget[i].grid(row = i, column =7, padx = 6, pady =4)
            atom_palette_str.set(str(self.get_atom_related_input_dict()['ListOfAtomPalette'][i]))
            
            tk.Label(atom_params_frame, text = 'subplot_arg:').grid(row = i, column = 8, padx =2, pady =4)
            Atomsubplot_arg_str = tk.StringVar()
            self.atom_subplot_arg_entry_widget[i] = tk.Entry(atom_params_frame, width = 4, textvariable = Atomsubplot_arg_str)
            self.atom_subplot_arg_entry_widget[i].grid(row = i, column =9, padx = 6, pady =4)
            Atomsubplot_arg_str.set(str(self.get_atom_related_input_dict()['ListOfAtomSubplotArg'][i]))

    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root, background='light gray')
        left_frame.grid(row = 10, column = 0, columnspan = 3, sticky = tk.E + tk.N + tk.S)
        for i in range(0, self.num_of_doscar):
            tk.Label(left_frame, text = 'DOSCAR No.' + str(i + 1) + ':').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text = 'Find', command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_doscar_entry_widget[i] = tk.Entry(left_frame, width = 15, textvariable= filename_str)
            self.load_doscar_entry_widget[i].grid(row = i, column =2, padx = 6, pady =4)
            filename_str.set(self.get_doscar_related_input_dict()['ListOfDOSCARPaths'][i])

    def create_top_bar(self):
        top_bar_frame = tk.Frame(self.root, height = 25)
        top_bar_frame.grid(row = 0, columnspan = 12, rowspan = 10, padx = 5, pady = 5)

        tk.Label(top_bar_frame, text = 'num_doscar:').grid(row = 0, column = 0, sticky = tk.E)
        self.num_of_doscar_widget = tk.Spinbox(top_bar_frame, from_ = 1, to = self.max_num_of_doscar, width = 3, command = self.on_num_of_doscar_changed)
        self.num_of_doscar_widget.grid(row = 0, column = 1)
        self.num_of_doscar_widget.delete(0,'end')
        self.num_of_doscar_widget.insert(0,self.initial_num_of_doscar)

        tk.Label(top_bar_frame, text = 'num_atoms:').grid(row = 0, column =2)
        self.num_of_atoms_widget = tk.Spinbox(top_bar_frame, from_ = 1, to = self.max_num_of_atoms, width = 3, command = self.on_num_of_atoms_changed)
        self.num_of_atoms_widget.grid(row = 0, column = 3)
        self.num_of_atoms_widget.delete(0,'end')
        self.num_of_atoms_widget.insert(0,self.initial_num_of_atoms)

        tk.Label(top_bar_frame, text = 'num_elments:').grid(row = 0, column = 4)
        self.num_of_elmts_widget = tk.Spinbox(top_bar_frame, from_ = 1, to = self.max_num_of_elmts, width = 3, command = self.on_num_of_elmts_changed)
        self.num_of_elmts_widget.grid(row = 0, column = 5)
        self.num_of_elmts_widget.delete(0,'end')
        self.num_of_elmts_widget.insert(0,self.initial_num_of_elmts)

        tk.Label(top_bar_frame, text = 'num_subplots:').grid(row = 0, column =6)
        self.num_of_subplots_widget = tk.Spinbox(top_bar_frame, from_ = 1, to = self.max_num_of_subplots, width = 3, command = self.on_num_of_subplots_changed)
        self.num_of_subplots_widget.grid(row = 0, column = 7)
        self.num_of_subplots_widget.delete(0,'end')
        self.num_of_subplots_widget.insert(0,self.initial_num_of_atoms)

        self.run_button = tk.Button(top_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = 0, column = 9, padx = 50)
        self.run_button.configure(background='lightgreen')


    def create_subplot_bar(self):
        subplot_bar_frame = tk.Frame(self.root, height = 25, background='light gray')
        start_row = self.max_num_of_doscar + 2
        subplot_bar_frame.grid(row = start_row, columnspan = 14, sticky = tk.W + tk.E + tk.N + tk.S, padx = 5, pady = 5)

        for i in range(0, self.num_of_subplots):
            tk.Label(subplot_bar_frame, text = 'subplot_arg:').grid(row = i, column =0, padx =2, pady =4)
            subplot_arg_str = tk.StringVar()
            self.subplot_arg_entry_widget[i] = tk.Entry(subplot_bar_frame, width = 4, textvariable = subplot_arg_str)
            self.subplot_arg_entry_widget[i].grid(row = i, column =1, padx = 6, pady =4)
            subplot_arg_str.set(str(self.get_subplot_related_input_dict()['ListOfSubplotArg'][i]))

            tk.Label(subplot_bar_frame, text = 'xlo:').grid(row = i, column =2, padx =2, pady =4)
            xlo_str = tk.StringVar()
            self.xlo_entry_widget[i] = tk.Entry(subplot_bar_frame, width = 5, textvariable = xlo_str)
            self.xlo_entry_widget[i].grid(row = i, column =3, padx = 6, pady =4)
            xlo_str.set(str(self.get_subplot_related_input_dict()['ListOfXLo'][i]))

            tk.Label(subplot_bar_frame, text = 'xhi:').grid(row = i, column =4, padx =2, pady =4)
            xhi_str = tk.StringVar()
            self.xhi_entry_widget[i] = tk.Entry(subplot_bar_frame, width = 5, textvariable = xhi_str)
            self.xhi_entry_widget[i].grid(row = i, column =5, padx = 6, pady =4)
            xhi_str.set(str(self.get_subplot_related_input_dict()['ListOfXHi'][i]))

            tk.Label(subplot_bar_frame, text = 'ylo:').grid(row = i, column =6, padx =2, pady =4)
            ylo_str = tk.StringVar()
            self.ylo_entry_widget[i] = tk.Entry(subplot_bar_frame, width = 5, textvariable = ylo_str)
            self.ylo_entry_widget[i].grid(row = i, column =7, padx = 6, pady =4)
            ylo_str.set(str(self.get_subplot_related_input_dict()['ListOfYLo'][i]))

            tk.Label(subplot_bar_frame, text = 'yhi:').grid(row = i, column =8, padx =2, pady =4)
            yhi_str = tk.StringVar()
            self.yhi_entry_widget[i] = tk.Entry(subplot_bar_frame, width = 5, textvariable = yhi_str)
            self.yhi_entry_widget[i].grid(row = i, column =9, padx = 6, pady =4)
            yhi_str.set(str(self.get_subplot_related_input_dict()['ListOfYHi'][i]))

            xtick_ticked_yes = tk.BooleanVar()
            self.subplot_xtick_check_button_widget[i] = tk.Checkbutton(subplot_bar_frame, text = 'subplot: xtick', variable = xtick_ticked_yes)
            self.subplot_xtick_check_button_widget[i].grid(row = i, column = 10, padx = 6, pady =4)
            xtick_ticked_yes.set(bool(self.get_subplot_related_input_dict()['ListOfSubplotXTick'][i]))
            self.subplot_xtick_check_button_widget[i].var = xtick_ticked_yes

            ytick_ticked_yes = tk.BooleanVar()
            self.subplot_ytick_check_button_widget[i] = tk.Checkbutton(subplot_bar_frame, text = 'subplot: ytick', variable = ytick_ticked_yes)
            self.subplot_ytick_check_button_widget[i].grid(row = i, column =11, padx = 6, pady =4)
            ytick_ticked_yes.set(bool(self.get_subplot_related_input_dict()['ListOfSubplotYTick'][i]))
            self.subplot_ytick_check_button_widget[i].var = ytick_ticked_yes

            xlabel_ticked_yes = tk.BooleanVar()
            self.subplot_xlabel_check_button_widget[i] = tk.Checkbutton(subplot_bar_frame, text = 'subplot: xlabel', variable = xlabel_ticked_yes)
            self.subplot_xlabel_check_button_widget[i].grid(row = i, column = 12, padx = 6, pady =4)
            xlabel_ticked_yes.set(bool(self.get_subplot_related_input_dict()['ListOfSubplotXLabel'][i]))
            self.subplot_xlabel_check_button_widget[i].var = xlabel_ticked_yes

            ylabel_ticked_yes = tk.BooleanVar()
            self.subplot_ylabel_check_button_widget[i] = tk.Checkbutton(subplot_bar_frame, text = 'subplot: ylabel', variable = ylabel_ticked_yes)
            self.subplot_ylabel_check_button_widget[i].grid(row = i, column =13, padx = 6, pady =4)
            ylabel_ticked_yes.set(bool(self.get_subplot_related_input_dict()['ListOfSubplotYLabel'][i]))
            self.subplot_ylabel_check_button_widget[i].var = ylabel_ticked_yes

    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_atom_params_panel()
        self.create_top_menu()
        self.create_top_bar()
        self.create_bottom_bar()
        self.create_elmt_check_bar()
        self.create_subplot_bar()
        
##################################
#plot_poscar
##################################
class plot_poscar_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'plot_poscar'
        self.initial_num_of_elmts = 1
        self.max_num_of_files = 1
        self.max_num_of_dir = 1
        self.max_num_of_elmts = 10
        self.num_of_elmts = self.initial_num_of_elmts
        self.initial_elmt_colr = elmt_color
        self.elmt_color = elmt_color

        self.Runplot_poscar = False
        self.Runplot_poscar_for_workdir = False
        self.euler_angle_type = 'zyz'
        self.phi = -3
        self.theta = 4
        self.psi = 0
        self.draw_mirror_atom = True
        self.box_on = True
        self.axis_indicator = True
        self.plot_cell_basis_vector_label = True
        self.plot_atom_label = True
        self.plot_poscar_mode = 'single'
        self.plot_poscar_mode_str = 'single'
        self.poscar_or_contcar = 'POSCAR'
        self.poscar_or_contcar_Str = 'POSCAR'
        self.fig_format = 'png'
        self.fig_format_str = 'png'
        self.fig_dpi = 100
        self.user_defined_elmt_color = False
        # colormap
        self.draw_colormap = False
        self.initial_colormap_column_indx = 1        
        self.colormap_vmin = None,
        self.colormap_vmax = None,
        self.vmin_color = 'blue',
        self.vmax_color = 'red',
        self.colorbar_alignment = 'vertical'
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.all_dir_params = [None] * self.max_num_of_dir
        self.current_file_indx = 0
        self.current_dir_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.load_dir_entry_widget = [None] * self.max_num_of_dir
        self.elmtname_entry_widget = [None] * self.max_num_of_elmts
        self.elmt_color_entry_widget = [None] * self.max_num_of_elmts
        
        self.colormap_column_indx_entry_widget = 1     
        self.colormap_vmin_entry_widget = None
        self.colormap_vmax_entry_widget = None
        self.vmin_color_entry_widget = 'blue'
        self.vmax_color_entry_widget = 'red'

        self.init_all_input_params()
        self.init_gui()
    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_dir_dict(self):
        return self.all_dir_params[self.current_dir_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def get_list_of_dir_paths(self):
        return self.get_dir_dict()['list_of_dir_paths']
    def set_num_of_elmts(self):
        self.num_of_elmts = int(self.num_of_elmts_widget.get())
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def set_dir_path(self, Dirindx, dir_path):
        self.get_list_of_dir_paths()[Dirindx] = dir_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                }
                               for k in range(self.max_num_of_files)]
        self.all_dir_params = [{'list_of_dir_paths': [None] * self.max_num_of_dir,
                                }
                               for k in range(self.max_num_of_dir)]

    def on_user_defined_elmt_color_ticked_yes_changed(self):
        self.user_defined_elmt_color = bool(self.user_defined_elmt_color_check_button_widget.var.get())
        if self.user_defined_elmt_color == True:
            self.create_elmt_num_panel()
            self.create_elmt_color_panel()
        elif self.user_defined_elmt_color == False:
            self.create_colormap_mask_panel_1()
            self.create_colormap_mask_panel_2()            
            self.elmt_color = periodic_table_dict['elmt_color']
     
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def on_open_dir_button_clicked(self, Dirindx):  
        def event_handler():
            dir_path = tkinter.filedialog.askdirectory()
            if not dir_path:
                return
            self.set_dir_path(Dirindx, dir_path)
            self.display_all_dir_paths()
        return event_handler

    def on_num_of_elmts_changed(self):
        self.set_num_of_elmts()
        self.create_elmt_color_panel()

    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_all_dir_paths(self):
        '''
        Description:
            display the directory path in each entry box
        '''            
        for indx, dir_path in enumerate(self.get_list_of_dir_paths()):
            self.display_dir_path(indx, dir_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)
    def display_dir_path(self, text_widget_indx, dir_path):
        '''
        Description:
            display the directory path in the TextWidgetNum-th entry box
        '''
        if dir_path is None:
            return
        self.load_dir_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_dir_entry_widget[text_widget_indx].insert(0, dir_path)
    def on_plot_poscar_mode_changed(self,value):
        self.plot_poscar_mode = value
        if self.plot_poscar_mode == 'single':
            self.create_left_file_loader()
        elif self.plot_poscar_mode == 'multiple':
            self.create_multiple_dir_loader()
    def on_draw_colormap_changed(self,value):
        self.draw_colormap = value
        if self.draw_colormap == 'True':
            self.create_colormap_panel()
            self.create_colormap_mask_panel_1()
            self.create_colormap_mask_panel_2()
        elif self.draw_colormap == 'False':
            self.user_defined_elmt_color = False
            self.create_user_defined_elmt_color_ticked_yes_bar()
            self.on_user_defined_elmt_color_ticked_yes_changed()
    def on_colorbar_alignment_changed(self,value):
        self.colorbar_alignment = value
    def on_run_button_clicked(self):
        periodic_table_dict = periodic_table.periodic_tab()
        elmt_color = periodic_table_dict['elmt_color']

        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        for i in range(0, self.max_num_of_dir):
            self.get_dir_dict()['list_of_dir_paths'][i] = self.load_dir_entry_widget[i].get()
        self.plot_poscar_mode = funcs.convert_bool(self.plot_poscar_mode_str.get())
        plot_poscar_mode = self.plot_poscar_mode
        if plot_poscar_mode == 'single':
            poscar_dir = self.get_input_dict()['list_of_file_paths'][0]
        elif plot_poscar_mode == 'multiple':
            workdir = self.get_dir_dict()['list_of_dir_paths'][0]

        self.euler_angle_type = self.euler_angle_type_entry_widget.get()        
        self.phi = self.phi_entry_widget.get()
        self.theta = self.theta_entry_widget.get()
        self.psi = self.psi_entry_widget.get()
        if self.user_defined_elmt_color == True:
            for i in range(0, self.num_of_elmts):
                self.elmt_color[self.elmtname_entry_widget[i].get()] = self.elmt_color_entry_widget[i].get()
            elmt_color.update(self.elmt_color)
        elif self.user_defined_elmt_color == False:
            elmt_color = periodic_table_dict['elmt_color']
        self.draw_mirror_atom = self.draw_mirror_atom_check_button_widget.var.get()
        self.box_on = self.box_on_check_button_widget.var.get()
        self.axis_indicator = self.axis_indicator_check_button_widget.var.get()
        self.plot_cell_basis_vector_label = self.plot_cell_basis_vector_label_check_button_widget.var.get()
        self.plot_atom_label = self.plot_atom_label_check_button_widget.var.get()
        self.poscar_or_contcar = self.poscar_or_contcar_Str.get()
        self.fig_format = self.fig_format_str.get()
        self.fig_dpi = int(self.fig_dpi_entry_widget.get())
        #colormap
        self.draw_colormap = funcs.convert_bool(self.draw_colormap_str.get())
        self.colormap_column_indx = int(self.colormap_column_indx_entry_widget.get())
        if self.colormap_vmin_entry_widget.get() != None and self.colormap_vmin_entry_widget.get() != 'None':
            self.colormap_vmin = float(self.colormap_vmin_entry_widget.get())
        else:
            self.colormap_vmin = None
        if self.colormap_vmax_entry_widget.get() != None and self.colormap_vmax_entry_widget.get() != 'None':
            self.colormap_vmax = float(self.colormap_vmax_entry_widget.get())
        else:
            self.colormap_vmax = None
        self.vmin_color = self.vmin_color_entry_widget.get()
        self.vmax_color = self.vmax_color_entry_widget.get()
        self.colorbar_alignment = self.colorbar_alignment_str.get()

        euler_angle_type = self.euler_angle_type
        phi = float(self.phi)
        theta = float(self.theta)
        psi = float(self.psi)
        draw_mirror_atom = self.draw_mirror_atom
        box_on = self.box_on
        axis_indicator = self.axis_indicator
        plot_cell_basis_vector_label = self.plot_cell_basis_vector_label
        plot_atom_label = self.plot_atom_label
        poscar_or_contcar = self.poscar_or_contcar
        fig_format = self.fig_format
        fig_dpi = self.fig_dpi
        #colormap
        draw_colormap = self.draw_colormap
        colormap_column_indx = self.colormap_column_indx
        colormap_vmin = self.colormap_vmin
        colormap_vmax = self.colormap_vmax
        vmin_color = self.vmin_color
        vmax_color = self.vmax_color
        colorbar_alignment = self.colorbar_alignment

        if plot_poscar_mode == 'single':
            vasp_plot.plot_poscar(
                poscar_dir, euler_angle_type, phi, theta, psi, elmt_color, draw_mirror_atom, box_on,
                axis_indicator, plot_cell_basis_vector_label, plot_atom_label, fig_format, fig_dpi,
                draw_colormap, colormap_column_indx, colormap_vmin, colormap_vmax, vmin_color, vmax_color, colorbar_alignment)
        elif plot_poscar_mode == 'multiple':
            vasp_plot.plot_poscar_for_workdir(
                workdir, euler_angle_type, phi, theta, psi, elmt_color, draw_mirror_atom, box_on,
                axis_indicator, plot_cell_basis_vector_label, plot_atom_label, poscar_or_contcar, fig_format, fig_dpi,
                draw_colormap, colormap_column_indx, colormap_vmin, colormap_vmax, vmin_color, vmax_color, colorbar_alignment)
        funcs.write_log(logfile,'## plot_poscar complete\n')
        return 'break'
    
    def create_elmt_num_panel(self):
        elmt_num_frame = tk.Frame(self.root, height = 25)
        start_row = 14
        elmt_num_frame.grid(row = start_row, columnspan = 6, rowspan = 1, sticky = tk.W + tk.E + tk.N + tk.S, padx = 5, pady = 5)

        tk.Label(elmt_num_frame, text = 'NumElements:').grid(row = start_row, column = 1)
        self.num_of_elmts_widget = tk.Spinbox(elmt_num_frame, from_ = 1, to = self.max_num_of_elmts, width = 3, command = self.on_num_of_elmts_changed)
        self.num_of_elmts_widget.grid(row = start_row, column = 5)
        self.num_of_elmts_widget.delete(0,'end')
        self.num_of_elmts_widget.insert(0,self.initial_num_of_elmts)
        
    def create_elmt_color_panel(self):
        elmt_color_frame = tk.Frame(self.root, height= 25, background = None)
        start_row = 15
        elmt_color_frame.grid(row = start_row, column = 0, columnspan = 2, sticky = tk.W + tk.E + tk.N + tk.S, padx = 5, pady = 5)
        self.num_of_elmts = int(self.num_of_elmts_widget.get())
        for i in range(0,self.num_of_elmts):
            tk.Label(elmt_color_frame,text = 'element:').grid(row = start_row + i + 1, column = 1, padx =2 , pady = 4)
            elmtname_str = tk.StringVar()
            self.elmtname_entry_widget[i] = tk.Entry(elmt_color_frame, width = 5, textvariable = elmtname_str)
            self.elmtname_entry_widget[i].grid(row = start_row + i + 1, column = 2, padx = 6, pady =4)
            elmtname_str.set('Ni')

            tk.Label(elmt_color_frame,text = 'color:').grid(row = start_row + i + 1, column = 3, padx =2 , pady = 4)
            elmt_color_str = tk.StringVar()
            self.elmt_color_entry_widget[i] = tk.Entry(elmt_color_frame, width = 9, textvariable = elmt_color_str)
            self.elmt_color_entry_widget[i].grid(row = start_row + i + 1, column = 4, padx = 6, pady =4)
            elmt_color_str.set('gray')

    def create_colormap_panel(self):
        colormap_frame = tk.Frame(self.root, height= 25, background = None)
        start_row = 13
        colormap_frame.grid(row = start_row, column = 0, columnspan = 2, sticky = tk.W + tk.E + tk.N + tk.S, padx = 5, pady = 5)

        tk.Label(colormap_frame,text = 'colormap_column_indx:').grid(row = start_row, column = 1, padx =2 , pady = 4)
        colormap_column_indx_str = tk.StringVar()
        self.colormap_column_indx_entry_widget = tk.Entry(colormap_frame, width = 5, textvariable = colormap_column_indx_str)
        self.colormap_column_indx_entry_widget.grid(row = start_row, column = 2, padx = 6, pady =4)
        colormap_column_indx_str.set('1')

        tk.Label(colormap_frame,text = 'colormap_vmin:').grid(row = start_row, column = 3, padx =2 , pady = 4)
        colormap_vmin_str = tk.StringVar()
        self.colormap_vmin_entry_widget = tk.Entry(colormap_frame, width = 5, textvariable = colormap_vmin_str)
        self.colormap_vmin_entry_widget.grid(row = start_row, column = 4, padx = 6, pady =4)
        colormap_vmin_str.set('None')

        tk.Label(colormap_frame,text = 'colormap_vmax:').grid(row = start_row, column = 5, padx =2 , pady = 4)
        colormap_vmax_str = tk.StringVar()
        self.colormap_vmax_entry_widget = tk.Entry(colormap_frame, width = 5, textvariable = colormap_vmax_str)
        self.colormap_vmax_entry_widget.grid(row = start_row, column = 6, padx = 6, pady =4)
        colormap_vmax_str.set('None')

        tk.Label(colormap_frame,text = 'vmin_color:').grid(row = start_row, column = 7, padx =2 , pady = 4)
        vmin_color_str = tk.StringVar()
        self.vmin_color_entry_widget = tk.Entry(colormap_frame, width = 5, textvariable = vmin_color_str)
        self.vmin_color_entry_widget.grid(row = start_row, column = 8, padx = 6, pady =4)
        vmin_color_str.set('blue')

        tk.Label(colormap_frame,text = 'vmax_color:').grid(row = start_row, column = 9, padx =2 , pady = 4)
        vmax_color_str = tk.StringVar()
        self.vmax_color_entry_widget = tk.Entry(colormap_frame, width = 5, textvariable = vmax_color_str)
        self.vmax_color_entry_widget.grid(row = start_row, column = 10, padx = 6, pady =4)
        vmax_color_str.set('red')

        tk.Label(colormap_frame,text = 'colorbar_alignment:').grid(row = start_row, column = 11, padx =2 , pady = 4)
        self.colorbar_alignment_str = tk.StringVar()
        options = ['vertical','horizontal']
        self.colorbar_alignment_option_menu_widget = tk.OptionMenu(colormap_frame,self.colorbar_alignment_str,*options, command = self.on_colorbar_alignment_changed)
        self.colorbar_alignment_option_menu_widget.grid(row = start_row, column =12, padx = 6, pady =4, sticky = tk.W + tk.E)
        self.colorbar_alignment_str.set(str(self.colorbar_alignment))

    def create_colormap_mask_panel_1(self):
        colormap_mask1_frame = tk.Frame(self.root, height= 25, background = None)
        start_row = 14
        colormap_mask1_frame.grid(row = start_row, column = 0, columnspan = 2, sticky = tk.W + tk.E + tk.N + tk.S, padx = 5, pady = 5)

    def create_colormap_mask_panel_2(self):
        colormap_mask2_frame = tk.Frame(self.root, height= 25, background = None)
        start_row = 15
        colormap_mask2_frame.grid(row = start_row, column = 0, columnspan = 2, sticky = tk.W + tk.E + tk.N + tk.S, padx = 5, pady = 5)


    def create_user_defined_elmt_color_ticked_yes_bar(self):
        user_defined_elmt_color_ticked_yes_bar_frame = tk.Frame(self.root, height = 25)
        StartRow = 13
        user_defined_elmt_color_ticked_yes_bar_frame.grid(row = StartRow, columnspan = 12, sticky = tk.W + tk.E, padx = 5, pady = 5)
        
        user_defined_elmt_color_ticked_yes = tk.BooleanVar()
        self.user_defined_elmt_color_check_button_widget = tk.Checkbutton(user_defined_elmt_color_ticked_yes_bar_frame, text = 'user-defined element color',
                                                                   variable = user_defined_elmt_color_ticked_yes, command = self.on_user_defined_elmt_color_ticked_yes_changed)
        self.user_defined_elmt_color_check_button_widget .grid(row = StartRow, column =0, padx = 6, pady =4)
        user_defined_elmt_color_ticked_yes.set(bool(self.user_defined_elmt_color))
        self.user_defined_elmt_color_check_button_widget.var = user_defined_elmt_color_ticked_yes  


    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root, background = 'lightblue')
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)

        for i in range(0, self.max_num_of_files):
            tk.Label(left_frame, text = 'single POSCAR file:').grid(row = 0, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text = 'Find', command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = 0, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 86, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = 0, column =2, padx = 6, pady =4)
            filename_str.set('For mode=single: Enter or choose the POSCAR file')

    def create_multiple_dir_loader(self):
        multiple_dir_frame = tk.Frame(self.root, background = 'gray')
        multiple_dir_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)

        for i in range(0, self.max_num_of_dir):
            tk.Label(multiple_dir_frame, text = 'POSCAR work directory:').grid(row = 1, column = 0, padx =2, pady = 4)
            OpenDirButton = tk.Button(multiple_dir_frame, text = 'Find', command = self.on_open_dir_button_clicked(i))
            OpenDirButton.grid(row = 1, column = 1, padx = 5, pady =4)
            DirNameStr = tk.StringVar()
            self.load_dir_entry_widget[i] = tk.Entry(multiple_dir_frame, width = 86, textvariable= DirNameStr)
            self.load_dir_entry_widget[i].grid(row = 1, column =2, padx = 6, pady =4)
            DirNameStr.set('For mode=multiple: Enter or choose the directory which contains subdirectories with multiple POSCARs')

            tk.Label(multiple_dir_frame,text = 'poscar_or_contcar=').grid(row = 1, column = 3, padx =2 , pady = 4, sticky = tk.W)
            self.poscar_or_contcar_Str = tk.StringVar()
            self.poscar_or_contcar_option_menu_widget = tk.OptionMenu(multiple_dir_frame,self.poscar_or_contcar_Str,'POSCAR','CONTCAR')
            self.poscar_or_contcar_option_menu_widget.grid(row = 1, column =4, padx = 6, pady =4)
            self.poscar_or_contcar_Str.set(str(self.poscar_or_contcar))

    def create_top_bar(self):
        top_bar_frame = tk.Frame(self.root, height = 25)
        start_row = 0
        top_bar_frame.grid(row = start_row, columnspan = 12, sticky = tk.W + tk.E, padx = 5, pady = 5)

        tk.Label(top_bar_frame,text = 'mode (plot single POSCAR or multiple POSCARs?):').grid(row = start_row, column = 0, padx =2 , pady = 4)
        self.plot_poscar_mode_str = tk.StringVar()
        options = ['single','multiple']
        self.plot_poscar_mode_option_menu_widget = tk.OptionMenu(top_bar_frame,self.plot_poscar_mode_str,*options, command = self.on_plot_poscar_mode_changed)
        self.plot_poscar_mode_option_menu_widget.grid(row = start_row, column =1, padx = 6, pady =4, sticky = tk.W + tk.E)
        self.plot_poscar_mode_str.set(str(self.plot_poscar_mode))

        tk.Label(top_bar_frame,text = 'color mapping:').grid(row = start_row, column = 2, padx =2 , pady = 4)
        self.draw_colormap_str = tk.StringVar()
        options = ['True','False']
        self.draw_colormap_option_menu_widget = tk.OptionMenu(top_bar_frame,self.draw_colormap_str,*options, command = self.on_draw_colormap_changed)
        self.draw_colormap_option_menu_widget.grid(row = start_row, column =3, padx = 6, pady =4, sticky = tk.W + tk.E)
        self.draw_colormap_str.set(str(self.draw_colormap))

        self.run_button = tk.Button(top_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 9, padx = 2, sticky = tk.W + tk.E)
        self.run_button.configure(background='lightgreen')


    def create_mid_bar(self):
        mid_bar_frame = tk.Frame(self.root, height = 25)
        start_row = 1
        mid_bar_frame.grid(row = start_row, columnspan = 12, sticky = tk.W + tk.E, padx = 5, pady = 5)

        tk.Label(mid_bar_frame,text = 'euler_angle_type=').grid(row = start_row, column = 0, padx =2 , pady = 4)
        euler_angle_typeStr = tk.StringVar()
        self.euler_angle_type_entry_widget = tk.Entry(mid_bar_frame, width = 5, textvariable = euler_angle_typeStr)
        self.euler_angle_type_entry_widget.grid(row = start_row, column =1, padx = 6, pady =4)
        euler_angle_typeStr.set(str(self.euler_angle_type))

        tk.Label(mid_bar_frame,text = 'phi=').grid(row = start_row, column = 2, padx =2 , pady = 4)
        phiStr = tk.StringVar()
        self.phi_entry_widget = tk.Entry(mid_bar_frame, width = 5, textvariable = phiStr)
        self.phi_entry_widget.grid(row = start_row, column =3, padx = 6, pady =4)
        phiStr.set(str(self.phi))

        tk.Label(mid_bar_frame,text = 'theta=').grid(row = start_row, column = 4, padx =2 , pady = 4)
        thetaStr = tk.StringVar()
        self.theta_entry_widget = tk.Entry(mid_bar_frame, width = 5, textvariable = thetaStr)
        self.theta_entry_widget.grid(row = start_row, column =5, padx = 6, pady =4)
        thetaStr.set(str(self.theta))

        tk.Label(mid_bar_frame,text = 'psi=').grid(row = start_row, column = 6, padx =2 , pady = 4)
        psiStr = tk.StringVar()
        self.psi_entry_widget = tk.Entry(mid_bar_frame, width = 5, textvariable = psiStr)
        self.psi_entry_widget.grid(row = start_row, column =7, padx = 6, pady =4)
        psiStr.set(str(self.psi))

        
    def create_bottom_bar(self):
        bottom_bar_frame = tk.Frame(self.root, height = 25)
        start_row = self.max_num_of_files + 3
        bottom_bar_frame.grid(row = start_row, columnspan = 12, sticky = tk.W + tk.E, padx = 5, pady = 5)
        
        draw_mirror_atom_ticked_yes = tk.BooleanVar()
        self.draw_mirror_atom_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'draw_mirror_atom', variable = draw_mirror_atom_ticked_yes)
        self.draw_mirror_atom_check_button_widget.grid(row = start_row, column =0, padx = 6, pady =4)
        draw_mirror_atom_ticked_yes.set(bool(self.draw_mirror_atom))
        self.draw_mirror_atom_check_button_widget.var = draw_mirror_atom_ticked_yes

        box_on_ticked_yes = tk.BooleanVar()
        self.box_on_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'box_on', variable = box_on_ticked_yes)
        self.box_on_check_button_widget.grid(row = start_row, column =1, padx = 6, pady =4)
        box_on_ticked_yes.set(bool(self.box_on))
        self.box_on_check_button_widget.var = box_on_ticked_yes

        axis_indicator_ticked_yes = tk.BooleanVar()
        self.axis_indicator_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'axis_indicator', variable = axis_indicator_ticked_yes)
        self.axis_indicator_check_button_widget.grid(row = start_row, column =2, padx = 6, pady =4)
        axis_indicator_ticked_yes.set(bool(self.axis_indicator))
        self.axis_indicator_check_button_widget.var = axis_indicator_ticked_yes

        plot_cell_basis_vector_label_ticked_yes = tk.BooleanVar()
        self.plot_cell_basis_vector_label_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'plot_cell_basis_vector_label', variable = plot_cell_basis_vector_label_ticked_yes)
        self.plot_cell_basis_vector_label_check_button_widget.grid(row = start_row, column =3, padx = 6, pady =4)
        plot_cell_basis_vector_label_ticked_yes.set(bool(self.box_on))
        self.plot_cell_basis_vector_label_check_button_widget.var = plot_cell_basis_vector_label_ticked_yes

        plot_atom_label_ticked_yes = tk.BooleanVar()
        self.plot_atom_label_check_button_widget = tk.Checkbutton(bottom_bar_frame, text = 'plot_atom_label', variable = plot_atom_label_ticked_yes)
        self.plot_atom_label_check_button_widget.grid(row = start_row, column =4, padx = 6, pady =4)
        plot_atom_label_ticked_yes.set(bool(self.plot_atom_label))
        self.plot_atom_label_check_button_widget.var = plot_atom_label_ticked_yes
        
        tk.Label(bottom_bar_frame,text = 'fig_format=').grid(row = start_row, column = 5, padx =2 , pady = 4)
        self.fig_format_str = tk.StringVar()
        self.fig_format_option_menu_widget = tk.OptionMenu(bottom_bar_frame,self.fig_format_str,'png', 'eps', 'pdf', 'tif', 'tiff', 'jpg', 'jpeg', 'svg', 'svgz', 'pgf', 'ps', 'raw', 'rgba')
        self.fig_format_option_menu_widget.grid(row = start_row, column =6, padx = 6, pady =4)
        self.fig_format_str.set(str(self.fig_format))

        tk.Label(bottom_bar_frame,text = 'fig_dpi=').grid(row = start_row, column = 7, padx =2 , pady = 4)
        fig_dpi_str = tk.StringVar()
        self.fig_dpi_entry_widget = tk.Entry(bottom_bar_frame, width = 5, textvariable = fig_dpi_str)
        self.fig_dpi_entry_widget.grid(row = start_row, column =8, padx = 6, pady =4)
        fig_dpi_str.set(str(self.fig_dpi))
        
    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_colormap_panel()
        self.create_colormap_mask_panel_1()
        self.create_colormap_mask_panel_2()
        self.create_multiple_dir_loader()
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_top_bar()
        self.create_mid_bar()
        self.create_bottom_bar()
        self.create_user_defined_elmt_color_ticked_yes_bar()

    def display_help_messagebox(event=None):
        tkinter.messagebox.showinfo(
            'Help', '\n Please read the manual')


class estruct_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'estruct'
        self.max_num_of_files = 1
        self.initial_sysname = 'System'
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.current_file_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.init_all_input_params()
        self.init_gui()
    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                'sysname' : self.initial_sysname
                                }
                               for k in range(self.max_num_of_files)]
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)

    def on_run_button_clicked(self):
        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        self.get_input_dict()['sysname'] = self.sysname_value_widget.get()
        
        doscar_dir = self.get_input_dict()['list_of_file_paths'][0]
        sysname = self.get_input_dict()['sysname']
        
        vasp_analyze.estruct(doscar_dir, sysname)
        return 'break'
    
    def create_run_bar(self):
        run_bar_frame = tk.Frame(self.root, height = 15)
        start_row = self.max_num_of_files + 10
        run_bar_frame.grid(row = start_row, columnspan = 13, sticky = tk.W + tk.E, padx = 15, pady = 10)
        self.run_button = tk.Button(run_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 1, padx = 2)
        self.run_button.configure(background='lightgreen')
    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root)
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)
        for i in range(0, self.max_num_of_files):
            tk.Label(left_frame, text = 'DOSCAR:').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text='Find',command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 120, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = i, column =4, padx = 6, pady =4)
            filename_str.set('Enter or choose DOSCAR path here')
    def create_top_bar(self):
        top_bar_frame = tk.Frame(self.root, height = 25)
        top_bar_frame.grid(row = 0, columnspan = 12, rowspan = 10, padx = 5, pady = 5)

        tk.Label(top_bar_frame, text = 'sysname=').grid(row = 0, column =2)
        sysname_str = tk.StringVar()
        self.sysname_value_widget = tk.Entry(top_bar_frame, textvariable = sysname_str)
        self.sysname_value_widget.grid(row = 0, column = 3, padx = 3, pady = 2)
        sysname_str.set(self.get_input_dict()['sysname'])

    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_top_bar()
        self.create_run_bar()

##########################
# write_outcar_with_force
##########################
class write_poscar_with_force_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'write_outcar_with_force'
        self.max_num_of_files = 1
        self.initial_ionic_step = 'last'
        self.initial_output_poscar_file_name = 'None'
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.current_file_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.init_all_input_params()
        self.init_gui()
    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                'ionic_step' : self.initial_ionic_step,
                                'output_poscar_file_name' : self.initial_output_poscar_file_name
                                }
                               for k in range(self.max_num_of_files)]
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)

    def on_run_button_clicked(self):
        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        self.get_input_dict()['ionic_step'] = self.ionic_step_value_widget.get()
        self.get_input_dict()['output_poscar_file_name'] = self.output_poscar_file_name_value_widget.get()
        
        outcar_dir = self.get_input_dict()['list_of_file_paths'][0]
        if self.get_input_dict()['ionic_step'] == 'last' or self.get_input_dict()['ionic_step'] == 'first':
            ionic_step = self.get_input_dict()['ionic_step']
        else:
            ionic_step = int(self.get_input_dict()['ionic_step'])
        if self.get_input_dict()['output_poscar_file_name'] == 'None':
            output_poscar_file_name = None
        else:
            pass
        
        vasp_write.write_poscar_with_force(outcar_dir, ionic_step, output_poscar_file_name)
        return 'break'
    
    def create_run_bar(self):
        run_bar_frame = tk.Frame(self.root, height = 15)
        start_row = self.max_num_of_files + 10
        run_bar_frame.grid(row = start_row, columnspan = 13, sticky = tk.W + tk.E, padx = 15, pady = 10)
        self.run_button = tk.Button(run_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 1, padx = 2)
        self.run_button.configure(background='lightgreen')
    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root)
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)
        for i in range(0, self.max_num_of_files):
            tk.Label(left_frame, text = 'OUTCAR:').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text='Find',command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 120, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = i, column =4, padx = 6, pady =4)
            filename_str.set('Enter or choose OUTCAR path here')
    def create_top_bar(self):
        top_bar_frame = tk.Frame(self.root, height = 25)
        top_bar_frame.grid(row = 0, columnspan = 12, rowspan = 10, padx = 5, pady = 5)

        tk.Label(top_bar_frame, text = 'ionic_step=').grid(row = 0, column =0)
        ionic_step_str = tk.StringVar()
        self.ionic_step_value_widget = tk.Entry(top_bar_frame, width = 7, textvariable = ionic_step_str)
        self.ionic_step_value_widget.grid(row = 0, column = 1, padx = 3, pady = 2)
        ionic_step_str.set(self.get_input_dict()['ionic_step'])

        tk.Label(top_bar_frame, text = 'output_poscar_file_name=').grid(row = 0, column =2)
        output_poscar_file_name_str = tk.StringVar()
        self.output_poscar_file_name_value_widget = tk.Entry(top_bar_frame, width = 37, textvariable = output_poscar_file_name_str)
        self.output_poscar_file_name_value_widget.grid(row = 0, column = 3, padx = 3, pady = 2)
        output_poscar_file_name_str.set(self.get_input_dict()['output_poscar_file_name'])

    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_top_bar()
        self.create_run_bar()

class concentration_profile_subtoolkit:
    def __init__(self, root):
        self.subtoolkitname = 'plot_proxigram_csv'
        self.max_num_of_files = 1
        self.initial_sysname = 'system'
        self.initial_visible_elmt = 'Ni'
        self.interpolation_on = False
        self.fig_format = 'png'
        self.fig_format_str = 'png'
        self.fig_width = 6
        self.fig_height = 5
        self.fig_dpi = 600
        
        self.root = root
        self.root.title(self.subtoolkitname)
        self.all_input_params = [None] * self.max_num_of_files
        self.current_file_indx = 0
        self.load_file_entry_widget = [None] * self.max_num_of_files
        self.init_all_input_params()
        self.init_gui()
    def get_input_dict(self):
        return self.all_input_params[self.current_file_indx]
    def get_list_of_file_paths(self):
        return self.get_input_dict()['list_of_file_paths']
    def set_file_path(self, file_indx, file_path):
        self.get_list_of_file_paths()[file_indx] = file_path
    def init_all_input_params(self):
        self.all_input_params = [{'list_of_file_paths': [None] * self.max_num_of_files,
                                'sysname' : self.initial_sysname,
                                'visible_elmt' : self.initial_visible_elmt
                                }
                               for k in range(self.max_num_of_files)]
    def on_open_file_button_clicked(self, file_indx):  
        def event_handler():
            file_path = tkinter.filedialog.askopenfilename()
            if not file_path:
                return
            self.set_file_path(file_indx, file_path)
            self.display_all_file_paths()
        return event_handler
    def display_all_file_paths(self):
        '''
        Description:
            display the file path in each entry box
        '''            
        for indx, file_path in enumerate(self.get_list_of_file_paths()):
            self.display_file_path(indx, file_path)            
    def display_file_path(self, text_widget_indx, file_path):
        '''
        Description:
            display the file path in the TextWidgetNum-th entry box
        '''
        if file_path is None:
            return
        self.load_file_entry_widget[text_widget_indx].delete(0, tk.END)
        self.load_file_entry_widget[text_widget_indx].insert(0, file_path)

    def on_run_button_clicked(self):
        for i in range(0, self.max_num_of_files):
            self.get_input_dict()['list_of_file_paths'][i] = self.load_file_entry_widget[i].get()
        self.get_input_dict()['sysname'] = self.sysname_value_widget.get()
        self.get_input_dict()['visible_elmt'] = [item.strip(' ') for item in self.visible_elmt_widget.get().strip().split(',')]
        
        proxigram_csv_dir = self.get_input_dict()['list_of_file_paths'][0]
        sysname = self.get_input_dict()['sysname']
        visible_elmt_list = self.get_input_dict()['visible_elmt']
        self.interpolation_on = self.interpolation_on_check_button_widget.var.get()
        interpolation_on = self.interpolation_on
        self.fig_format = self.fig_format_str.get()
        self.fig_width = self.fig_width_entry_widget.get()
        self.fig_height = self.fig_height_entry_widget.get()
        self.fig_dpi = int(self.fig_dpi_entry_widget.get())

        fig_format = self.fig_format
        fig_dpi = self.fig_dpi
        fig_width = float(self.fig_width)
        fig_height = float(self.fig_height)
        
        apt_plot.plot_proxigram_csv(proxigram_csv_dir, sysname, visible_elmt_list, interpolation_on,
                                    fig_width, fig_height, fig_dpi, fig_format)
        return 'break'
    
    def create_run_bar(self):
        run_bar_frame = tk.Frame(self.root, height = 15)
        start_row = self.max_num_of_files + 10
        run_bar_frame.grid(row = start_row, columnspan = 13, sticky = tk.W + tk.E, padx = 15, pady = 10)
        self.run_button = tk.Button(run_bar_frame, text = 'Run', compound = 'left', command = self.on_run_button_clicked)
        self.run_button.grid(row = start_row, column = 1, padx = 2)
        self.run_button.configure(background='lightgreen')
    def create_left_file_loader(self):
        left_frame = tk.Frame(self.root)
        left_frame.grid(row = 10, column = 0, columnspan = 6, sticky = tk.W + tk.E + tk.N + tk.S)
        for i in range(0, self.max_num_of_files):
            tk.Label(left_frame, text = 'proxigram_csv_file:').grid(row = i, column = 0, padx =2, pady = 4)
            open_file_button = tk.Button(left_frame, text='Find',command = self.on_open_file_button_clicked(i))
            open_file_button.grid(row = i, column = 1, padx = 5, pady =4)
            filename_str = tk.StringVar()
            self.load_file_entry_widget[i] = tk.Entry(left_frame, width = 120, textvariable= filename_str)
            self.load_file_entry_widget[i].grid(row = i, column =4, padx = 6, pady =4)
            filename_str.set('Enter or choose proxigram_csv_file path here')
    def create_top_bar(self):
        top_bar_frame = tk.Frame(self.root, height = 25)
        top_bar_frame.grid(row = 0, columnspan = 12, rowspan = 10, padx = 5, pady = 5, sticky = tk.W + tk.E + tk.N + tk.S)

        tk.Label(top_bar_frame, text = 'sysname=').grid(row = 0, column =0)
        sysname_str = tk.StringVar()
        self.sysname_value_widget = tk.Entry(top_bar_frame, textvariable = sysname_str)
        self.sysname_value_widget.grid(row = 0, column = 1, padx = 3, pady = 2)
        sysname_str.set(self.get_input_dict()['sysname'])

        tk.Label(top_bar_frame, text = 'visible element=').grid(row = 0, column =2)
        visible_elmt_str = tk.StringVar()
        self.visible_elmt_widget = tk.Entry(top_bar_frame, width = 32, textvariable = visible_elmt_str)
        self.visible_elmt_widget.grid(row = 0, column = 3, padx = 3, pady = 2)
        visible_elmt_str.set(self.get_input_dict()['visible_elmt'])

        interpolation_on_ticked_yes = tk.BooleanVar()
        self.interpolation_on_check_button_widget = tk.Checkbutton(top_bar_frame, text = 'interpolation', variable = interpolation_on_ticked_yes)
        self.interpolation_on_check_button_widget.grid(row = 0, column =4, padx = 6, pady =4)
        interpolation_on_ticked_yes.set(bool(self.interpolation_on))
        self.interpolation_on_check_button_widget.var = interpolation_on_ticked_yes

        tk.Label(top_bar_frame,text = 'fig_width:').grid(row = 0, column = 5, padx =2 , pady = 4)
        fig_width_str = tk.StringVar()
        self.fig_width_entry_widget = tk.Entry(top_bar_frame, width = 5, textvariable = fig_width_str)
        self.fig_width_entry_widget.grid(row = 0, column =6, padx = 6, pady =4)
        fig_width_str.set(str(self.fig_width))

        tk.Label(top_bar_frame,text = 'fig_height:').grid(row = 0, column = 7, padx =2 , pady = 4)
        fig_heightStr = tk.StringVar()
        self.fig_height_entry_widget = tk.Entry(top_bar_frame, width = 5, textvariable = fig_heightStr)
        self.fig_height_entry_widget.grid(row = 0, column =8, padx = 6, pady =4)
        fig_heightStr.set(str(self.fig_height))

        tk.Label(top_bar_frame,text = 'fig_dpi:').grid(row = 0, column = 9, padx =2 , pady = 4)
        fig_dpi_str = tk.StringVar()
        self.fig_dpi_entry_widget = tk.Entry(top_bar_frame, width = 5, textvariable = fig_dpi_str)
        self.fig_dpi_entry_widget.grid(row = 0, column =10, padx = 6, pady =4)
        fig_dpi_str.set(str(self.fig_dpi))
        
        tk.Label(top_bar_frame,text = 'fig_format:').grid(row = 0, column = 11, padx =2 , pady = 4)
        self.fig_format_str = tk.StringVar()
        self.fig_format_option_menu_widget = tk.OptionMenu(top_bar_frame,self.fig_format_str,'png', 'eps', 'pdf', 'tif', 'tiff', 'jpg', 'jpeg', 'svg', 'svgz', 'pgf', 'ps', 'raw', 'rgba')
        self.fig_format_option_menu_widget.grid(row = 0, column =12, padx = 6, pady =4)
        self.fig_format_str.set(str(self.fig_format))        


    def create_top_menu(self):
        self.menu_bar = tk.Menu(self.root)
        self.file_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.file_menu.add_command(label = 'Exit', accelerator='Alt+F4', command=exit_editor)
        self.menu_bar.add_cascade(label = 'File', menu = self.file_menu)
        self.about_menu = tk.Menu(self.menu_bar, tearoff = 0)
        self.about_menu.add_command(label = 'About', command=display_about_messagebox)
        self.about_menu.add_command(label = 'Help', command=display_help_messagebox)
        self.menu_bar.add_cascade(label = 'About', menu = self.about_menu)
        self.root.config(menu = self.menu_bar)
    def init_gui(self):
        self.create_left_file_loader()
        self.create_top_menu()
        self.create_top_bar()
        self.create_run_bar()


#######################
#Main GUI
#######################
def Openfilename():
    import tkinter.filedialog
    file_path = tkinter.filedialog.askopenfilename()
    if not file_path:
        return
    return file_path

def OpenDir():
    import tkinter.filedialog
    file_dir = tkinter.filedialog.askdirectory()
    if not file_dir:
        return
    return file_dir

def new_file(event=None):
    root.title("Untitled")
    global file_name
    file_name = None
    content_text.delete(1.0, tk.END)

def open_file(event=None):
    input_file_name = tkinter.filedialog.askopenfilename(defaultextension=".txt",
                                                         filetypes=[("All Files", "*.*"), ("Text Documents", "*.txt")])
    if input_file_name:
        global file_name
        file_name = input_file_name
        root.title('{} - {}'.format(os.path.basename(file_name), PROGRAM_NAME))
        content_text.delete(1.0, tk.END)
        with open(file_name) as _file:
            content_text.insert(1.0, _file.read())


def write_to_file(file_name):
    try:
        content = content_text.get(1.0, 'end')
        with open(file_name, 'w') as the_file:
            the_file.write(content)
    except IOError:
        pass

def save_as(event=None):
    input_file_name = tkinter.filedialog.asksaveasfilename(defaultextension=".txt",
                                                           filetypes=[("All Files", "*.*"), ("Text Documents", "*.txt")])
    if input_file_name:
        global file_name
        file_name = input_file_name
        write_to_file(file_name)
        root.title('{} - {}'.format(os.path.basename(file_name), PROGRAM_NAME))
    return "break"


def save(event=None):
    global file_name
    if not file_name:
        save_as()
    else:
        write_to_file(file_name)
    return "break"

def display_about_messagebox(event=None):
    tkinter.messagebox.showinfo(
        'About', program_name + '\n version: ' + str(version) + ' \n' + str(authors))


def display_help_messagebox(event=None):
    tkinter.messagebox.showinfo(
        'Help', '\n Please read the manual')


def exit_editor(event=None):
    if tkinter.messagebox.askokcancel("Quit?", "Really quit?"):
        root.destroy()

menu_bar = tk.Menu(root)

file_menu = tk.Menu(menu_bar, tearoff=0)
##file_menu.add_command(label='New', accelerator='Ctrl+N',
##                      compound='left', underline=0, command=new_file)
##file_menu.add_command(label='Open', accelerator='Ctrl+O',
##                      compound='left', underline=0, command=open_file)
##file_menu.add_command(label='Save', accelerator='Ctrl+S',
##                      compound='left', underline=0, command=save)
##file_menu.add_command(label='Save as', accelerator='Shift+Ctrl+S', command=save_as)
##file_menu.add_separator()
file_menu.add_command(label='Exit', accelerator='Alt+F4', command=exit_editor)
menu_bar.add_cascade(label='File', menu=file_menu)

help_menu = tk.Menu(menu_bar, tearoff=0)
help_menu.add_command(label='Help', command=display_help_messagebox)
menu_bar.add_cascade(label='Help',  menu=help_menu)

about_menu = tk.Menu(menu_bar, tearoff=0)
about_menu.add_command(label='About', command=display_about_messagebox)
menu_bar.add_cascade(label='About',  menu=about_menu)
root.config(menu=menu_bar)

root.protocol('WM_DELETE_WINDOW', exit_editor)

#####################
#substitution button
#####################
subst_button = tk.Button(root, text='VASP: Build model by substitution', anchor = tk.W, compound = 'left', width=37)
subst_button.bind('<Button-1>', on_run_subst_button_clicked)
subst_button.pack()

#########################
#selection_sphere button
#########################
subst_button = tk.Button(root, text='VASP: Build model by selection (sphere)', anchor = tk.W, compound = 'left', width=37)
subst_button.bind('<Button-1>', on_run_vasp_selection_sphere_button_clicked)
subst_button.pack()

#####################
#plot_poscar button
#####################
plot_poscar_button = tk.Button(root, text='VASP: Plot POSCAR (support color mapping)', anchor = tk.W, compound = 'left', width=37)
plot_poscar_button.bind('<Button-1>', on_run_plot_poscar_button_clicked)
plot_poscar_button.pack()

#################
#plot_dos button
#################
plot_dos_button = tk.Button(root, text='VASP: Plot DOS', anchor = tk.W, compound = 'left', width=37)
plot_dos_button.bind('<Button-1>', on_run_plot_dos_button_clicked)
plot_dos_button.pack()

################
#nn_map button
################
nn_map_button = tk.Button(root, text='VASP: NN map', anchor = tk.W, compound = 'left', width=37)
nn_map_button.bind('<Button-1>', on_run_nn_map_button_clicked)
nn_map_button.pack()

################
#simple_cna
################
simple_cna_button = tk.Button(root, text='VASP: Simple CNA', anchor = tk.W, compound = 'left', width=37)
simple_cna_button.bind('<Button-1>', on_run_simple_cna_button_clicked)
simple_cna_button.pack()

################
#estruct button
################
estruct_button = tk.Button(root, text='VASP: Calculate E_struct', anchor = tk.W, compound = 'left', width=37)
estruct_button.bind('<Button-1>', on_run_estruct_button_clicked)
estruct_button.pack()

################################
#write_poscar_with_force button
################################
write_poscar_with_force_button = tk.Button(root, text='VASP: Write POSCAR with force', anchor = tk.W, compound = 'left', width=37)
write_poscar_with_force_button.bind('<Button-1>', on_run_write_poscar_with_force_button_clicked)
write_poscar_with_force_button.pack()

###########################
#plot_proxigram_csv button
###########################
plot_proxigram_csv_button = tk.Button(root, text='APT: Plot concentration profile', anchor = tk.W, compound = 'left', width=37)
plot_proxigram_csv_button.bind('<Button-1>', on_run_concentration_profile_button_clicked)
plot_proxigram_csv_button.pack()

#####################
#read_py_script button
#####################
read_py_script_button = tk.Button(root, text='Read *.py or *.log file', anchor = tk.W, compound = 'left', width=37)
read_py_script_button.bind('<Button-1>', on_run_read_py_script_button_clicked)
read_py_script_button.pack()


root.mainloop()

time_end = time.time()
formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time_end))
funcs.write_log(logfile,'## Job ended at '  + formatted_time + '\n' +
               '## time(Main_GUI) ' + str(time_end-time_start) + ' s')

##if __name__ == '__main__':
##    main()
