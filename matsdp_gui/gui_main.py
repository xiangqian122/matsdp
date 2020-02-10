# -*- coding: utf-8 -*-

# Created by: PyQt5 UI code generator 5.13.0

import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")    
import matplotlib.pyplot as plt

import matsdp
from matsdp import default_params
from matsdp import funcs
from matsdp import periodic_table
from matsdp.vasp import vasp_read
from matsdp.vasp import vasp_build
from matsdp.vasp import vasp_plot
from matsdp.vasp import vasp_analyze
from matsdp.vasp import vasp_write
from matsdp.vasp import vasp_default
from matsdp.vasp import vasp_help
from matsdp.apt import apt_plot
from matsdp.dvm import dvm_read
from matsdp.dvm import dvm_build
from matsdp.dvm import dvm_analyze
from matsdp.dvm import dvm_write
from matsdp.dvm import dvm_default
from matsdp.dvm import dvm_help

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog, QMessageBox

##from matsdp_vasp_build import *
##from gui_vasp_build import Ui_vasp_buildWindow
from gui_vasp.gui_vasp_build import Ui_vasp_buildWindow
from gui_vasp.gui_vasp_plot import Ui_vasp_plotWindow
from gui_vasp.gui_vasp_analyze import Ui_vasp_analyzeWindow
from gui_vasp.gui_vasp_write import Ui_vasp_writeWindow
from gui_apt.gui_apt_plot import Ui_apt_plotWindow
from gui_dvm.gui_dvm_analyze import Ui_dvm_analyzeWindow
from gui_run_py_script import Ui_run_py_scriptWindow

defaults_dict = default_params.default_params()
logfile = defaults_dict['logfile']

periodic_table_dict = periodic_table.periodic_tab()
dos_mode = periodic_table_dict['dos_mode']

# set the maximum recursion depth
sys.setrecursionlimit(100000)

time_start = time.time()
formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time_start))
# log file for the program. The results will be dumped into this log file
funcs.write_log(logfile,'## matsdp(GUI) started at ' + formatted_time)

program_name = 'matsdp'
version = '0.2.0'
authors = 'dianwuwang@163.com'
manual_url = 'https://github.com/dianwdw/matsdp/blob/master/doc/matsdp_manual.pdf'
program_url = 'https://github.com/dianwdw/matsdp'

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(565, 413)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(120, 50, 321, 271))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.vasp_build_button = QtWidgets.QPushButton(self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.vasp_build_button.setFont(font)
        self.vasp_build_button.setObjectName("vasp_build_button")
        self.verticalLayout.addWidget(self.vasp_build_button)
        self.vasp_plot_button = QtWidgets.QPushButton(self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.vasp_plot_button.setFont(font)
        self.vasp_plot_button.setObjectName("vasp_plot_button")
        self.verticalLayout.addWidget(self.vasp_plot_button)
        self.vasp_analyze_button = QtWidgets.QPushButton(self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.vasp_analyze_button.setFont(font)
        self.vasp_analyze_button.setObjectName("vasp_analyze_button")
        self.verticalLayout.addWidget(self.vasp_analyze_button)
        self.vasp_write_button = QtWidgets.QPushButton(self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.vasp_write_button.setFont(font)
        self.vasp_write_button.setObjectName("vasp_write_button")
        self.verticalLayout.addWidget(self.vasp_write_button)
        self.apt_plot_button = QtWidgets.QPushButton(self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.apt_plot_button.setFont(font)
        self.apt_plot_button.setObjectName("apt_plot_button")
        self.verticalLayout.addWidget(self.apt_plot_button)
        self.dvm_analyze_button = QtWidgets.QPushButton(self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.dvm_analyze_button.setFont(font)
        self.dvm_analyze_button.setObjectName("dvm_analyze_button")
        self.verticalLayout.addWidget(self.dvm_analyze_button)
        self.run_py_script_button = QtWidgets.QPushButton(self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(12)
        self.run_py_script_button.setFont(font)
        self.run_py_script_button.setObjectName("run_py_script_button")
        self.verticalLayout.addWidget(self.run_py_script_button)

        #set font size
        pointsize = font.pointSize()
        font.setPixelSize(pointsize * 90.0/50.0)
        app.setFont(font)
        
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 565, 18))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuAbout = QtWidgets.QMenu(self.menubar)
        self.menuAbout.setObjectName("menuAbout")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())
        
        self.action_exit_program = QtWidgets.QAction(MainWindow)
        self.action_exit_program.setObjectName("action_exit_program")
        self.menuFile.addAction(self.action_exit_program)
        self.action_exit_program.triggered.connect(self.on_exit_program)

        self.action_about_program = QtWidgets.QAction(MainWindow)
        self.action_about_program.setObjectName("action_about_program")
        self.menuAbout.addAction(self.action_about_program)
        self.action_about_program.triggered.connect(self.on_about_program)

        self.action_manual = QtWidgets.QAction(MainWindow)
        self.action_manual.setObjectName("action_manual")
        self.menuHelp.addAction(self.action_manual)
        self.action_manual.triggered.connect(self.on_manual)

        self.action_updates = QtWidgets.QAction(MainWindow)
        self.action_updates.setObjectName("action_updates")
        self.menuHelp.addAction(self.action_updates)
        self.action_updates.triggered.connect(self.on_updates)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", program_name))
        self.vasp_build_button.setStatusTip(_translate("MainWindow", "build model"))
        self.vasp_build_button.setText(_translate("MainWindow", "VASP: Build"))
        self.vasp_plot_button.setStatusTip(_translate("MainWindow", "vasp plot"))
        self.vasp_plot_button.setText(_translate("MainWindow", "VASP: Plot"))
        self.vasp_analyze_button.setStatusTip(_translate("MainWindow", "VASP: Analyze"))
        self.vasp_analyze_button.setText(_translate("MainWindow", "VASP: Analyze"))
        self.vasp_write_button.setStatusTip(_translate("MainWindow", "vasp write"))
        self.vasp_write_button.setText(_translate("MainWindow", "VASP: Write"))
        self.apt_plot_button.setStatusTip(_translate("MainWindow", "APT plot"))
        self.apt_plot_button.setText(_translate("MainWindow", "APT: Plot"))
        self.dvm_analyze_button.setStatusTip(_translate("MainWindow", "dvm_analyze"))
        self.dvm_analyze_button.setText(_translate("MainWindow", "DVM: Analyze"))
        self.run_py_script_button.setStatusTip(_translate("MainWindow", "run_py_script"))
        self.run_py_script_button.setText(_translate("MainWindow", "Run *.py or *.log file"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuAbout.setTitle(_translate("MainWindow", "About"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.action_exit_program.setText(_translate("MainWindow", "Exit"))
        self.action_about_program.setText(_translate("MainWindow", "About matsdp"))
        self.action_manual.setText(_translate("MainWindow", "Manual"))
        self.action_updates.setText(_translate("MainWindow", "Check for updates"))

    def on_exit_program(self):
        sys.exit(app.exec_())
        
    def on_about_program(self):
        QMessageBox.about(None, "About matsdp", program_name + '-' + version + '\n' + authors)

    def on_manual(self):
        QtGui.QDesktopServices.openUrl(QtCore.QUrl(manual_url))

    def on_updates(self):
        QtGui.QDesktopServices.openUrl(QtCore.QUrl(program_url))

        
class Controller:

    def __init__(self):
        pass

    def Show_MainWindow(self):

        self.MainWindow = QtWidgets.QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.MainWindow)
        self.ui.vasp_build_button.clicked.connect(self.Show_vasp_buildWindow)
        self.ui.vasp_plot_button.clicked.connect(self.Show_vasp_plotWindow)
        self.ui.vasp_analyze_button.clicked.connect(self.Show_vasp_analyzeWindow)
        self.ui.vasp_write_button.clicked.connect(self.Show_vasp_writeWindow)
        self.ui.apt_plot_button.clicked.connect(self.Show_apt_plotWindow)
        self.ui.dvm_analyze_button.clicked.connect(self.Show_dvm_analyzeWindow)
        self.ui.run_py_script_button.clicked.connect(self.Show_run_py_scriptWindow)

        self.MainWindow.show()
        
    #vasp_build
    def Show_vasp_buildWindow(self):        
        self.vasp_buildWindow = QtWidgets.QMainWindow()
        self.ui = Ui_vasp_buildWindow()
        self.ui.setupUi(self.vasp_buildWindow)
        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_substitution)
        self.ui.actionRun_1.triggered.connect(self.on_run_module_button_clicked_selection_sphere)

        self.vasp_buildWindow.show()

    def on_run_module_button_clicked_substitution(self):
        self.ui.poscar_file_path_s0 = self.ui.lineEdit_s7.text()
        self.ui.substitution_list_file_path_s0 = self.ui.lineEdit_s8.text()

        
        vasp_build.substitution(
            substitution_list_file = self.ui.substitution_list_file_path_s0,
            poscar_file_path = self.ui.poscar_file_path_s0,
            )

    def on_run_module_button_clicked_selection_sphere(self):
        self.ui.origin_atom_name_t0 = self.ui.lineEdit_t1.text()
        self.ui.radius_t0 = float(self.ui.lineEdit_t2.text())
        self.ui.output_file_name_t0 = self.ui.lineEdit_t3.text()
        self.ui.poscar_file_path_t0 = self.ui.lineEdit_t7.text()
        
        if self.ui.checkBox_t4.isChecked():
            self.ui.include_mirror_atoms_t0= True
        else:
            self.ui.include_mirror_atoms_t0 = False

        vasp_build.selection_sphere(
            poscar_file_path = self.ui.poscar_file_path_t0,
            origin_atom_name = self.ui.origin_atom_name_t0,
            radius = self.ui.radius_t0,
            include_mirror_atoms = self.ui.include_mirror_atoms_t0,
            output_file_name = self.ui.output_file_name_t0,
            )

    #vasp_plot
    def Show_vasp_plotWindow(self):        
        self.vasp_plotWindow = QtWidgets.QMainWindow()
        self.ui = Ui_vasp_plotWindow()
        self.ui.setupUi(self.vasp_plotWindow)
        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_plot_poscar)
        self.ui.actionRun_1.triggered.connect(self.on_run_module_button_clicked_plot_doscar)

        self.vasp_plotWindow.show()

    def on_run_module_button_clicked_plot_poscar(self):
        elmt_color = periodic_table_dict['elmt_color'].copy()
        self.ui.num_elmts_a0 = self.ui.num_elmts_a0
        self.ui.plot_poscar_mode = self.ui.plot_poscar_mode
        if self.ui.plot_poscar_mode == 'single':
            self.ui.poscar_file_path = self.ui.lineEdit_a7.text()
        elif self.ui.plot_poscar_mode == 'multiple':
            self.ui.workdir = self.ui.lineEdit_a7.text()

        self.ui.euler_angle_type = self.ui.lineEdit_a1.text()
        self.ui.phi = float(self.ui.lineEdit_a2.text())
        self.ui.theta = float(self.ui.lineEdit_a3.text())
        self.ui.psi = float(self.ui.lineEdit_a4.text())
        self.ui.draw_mirror_atom = funcs.convert_bool(self.ui.checkBox_a0.isChecked())
        self.ui.box_on = funcs.convert_bool(self.ui.checkBox_a1.isChecked())
        self.ui.axis_indicator = funcs.convert_bool(self.ui.checkBox_a2.isChecked())
        self.ui.plot_cell_basis_vector_label = funcs.convert_bool(self.ui.checkBox_a3.isChecked())
##        self.ui.plot_atom_label = funcs.convert_bool(self.ui.checkBox_a4.isChecked())
        self.ui.plot_atom_label = self.ui.lineEdit_a44.text()
        self.ui.fig_format = self.ui.fig_format_a0
        self.ui.fig_dpi = float(self.ui.lineEdit_a6.text())
        self.ui.draw_colormap = funcs.convert_bool(self.ui.draw_colormap)
        self.ui.colormap_column_indx = int(self.ui.lineEdit_a31.text())
        self.ui.colormap_vmin = self.ui.lineEdit_a32.text()
        if self.ui.colormap_vmin not in [None,'None','none']:
            self.ui.colormap_vmin = float(self.ui.colormap_vmin)
        else:
            self.ui.colormap_vmin = None
        self.ui.colormap_vmax = self.ui.lineEdit_a33.text()
        if self.ui.colormap_vmax not in [None,'None','none']:
            self.ui.colormap_vmax = float(self.ui.colormap_vmax)
        else:
            self.ui.colormap_vmax = None
        self.ui.vmin_color = self.ui.lineEdit_a34.text()
        self.ui.vmax_color = self.ui.lineEdit_a35.text()
        self.ui.colorbar_alignment = self.ui.colorbar_alignment

        self.ui.poscar_or_contcar = self.ui.poscar_or_contcar

        self.ui.user_defined_elmt_color = self.ui.user_defined_elmt_color

        if self.ui.user_defined_elmt_color == True:
            for i in range(0, self.ui.num_elmts_a0):
                self.ui.elmt_color_dict_a0[self.ui.lineEdit_a22_list[i].text()] = self.ui.lineEdit_a23_list[i].text()
            elmt_color.update(self.ui.elmt_color_dict_a0)
        elif self.ui.user_defined_elmt_color == False:
            elmt_color = periodic_table_dict['elmt_color'].copy()
        self.ui.elmt_color_dict_a0 = elmt_color
        
        if self.ui.plot_poscar_mode == 'single':
            vasp_plot.plot_poscar(
                poscar_file_path = self.ui.poscar_file_path,
                euler_angle_type = self.ui.euler_angle_type,
                phi = self.ui.phi,
                theta = self.ui.theta,
                psi = self.ui.psi,
                elmt_color = self.ui.elmt_color_dict_a0,
                draw_mirror_atom = self.ui.draw_mirror_atom,
                box_on = self.ui.box_on,
                axis_indicator = self.ui.axis_indicator,
                plot_cell_basis_vector_label = self.ui.plot_cell_basis_vector_label,
                plot_atom_label = self.ui.plot_atom_label,
                fig_format = self.ui.fig_format,
                fig_dpi = self.ui.fig_dpi,
                draw_colormap = self.ui.draw_colormap,
                colormap_column_indx = self.ui.colormap_column_indx,
                colormap_vmin = self.ui.colormap_vmin,
                colormap_vmax = self.ui.colormap_vmax,
                vmin_color = self.ui.vmin_color,
                vmax_color = self.ui.vmax_color,
                colorbar_alignment = self.ui.colorbar_alignment,
                )
        elif self.ui.plot_poscar_mode == 'multiple':
            vasp_plot.plot_poscar_for_workdir(
                workdir = self.ui.workdir,
                euler_angle_type = self.ui.euler_angle_type,
                phi = self.ui.phi,
                theta = self.ui.theta,
                psi = self.ui.psi,
                elmt_color = self.ui.elmt_color_dict_a0,
                draw_mirror_atom = self.ui.draw_mirror_atom,
                box_on = self.ui.box_on,
                axis_indicator = self.ui.axis_indicator,
                plot_cell_basis_vector_label = self.ui.plot_cell_basis_vector_label,
                plot_atom_label = self.ui.plot_atom_label,
                poscar_or_contcar = self.ui.poscar_or_contcar,
                fig_format = self.ui.fig_format,
                fig_dpi = self.ui.fig_dpi,
                draw_colormap = self.ui.draw_colormap,
                colormap_column_indx = self.ui.colormap_column_indx,
                colormap_vmin = self.ui.colormap_vmin,
                colormap_vmax = self.ui.colormap_vmax,
                vmin_color = self.ui.vmin_color,
                vmax_color = self.ui.vmax_color,
                colorbar_alignment = self.ui.colorbar_alignment,
                )

    def on_run_module_button_clicked_plot_doscar(self):
        #doscar
        for i_doscar in range(0, self.ui.num_doscar):
            self.ui.doscar_list[i_doscar] = self.ui.lineEdit_list[i_doscar].text()

        #atom
        for i_atom in range(0, self.ui.num_atoms):
            self.ui.atom_doscar_indx_list[i_atom] = self.ui.spinBox_9_list[i_atom].value() - 1
            self.ui.atom_doscar_file_path_list[i_atom] = self.ui.doscar_list[self.ui.atom_doscar_indx_list[i_atom]]
            self.ui.atom_sysname_list[i_atom] = self.ui.lineEdit_12_list[i_atom].text()
            self.ui.atom_name_list[i_atom] = self.ui.lineEdit_13_list[i_atom].text()
            self.ui.atom_palette_list[i_atom] = self.ui.lineEdit_14_list[i_atom].text()
            self.ui.atom_subplot_arg_list[i_atom] = self.ui.lineEdit_15_list[i_atom].text()
        #subplot
        for i_subplot in range(0, self.ui.num_subplots):
            self.ui.subplot_arg_list[i_subplot] = self.ui.lineEdit_21_list[i_subplot].text()
            self.ui.subplot_xlo_list[i_subplot] = self.ui.lineEdit_22_list[i_subplot].text()
            self.ui.subplot_xhi_list[i_subplot] = self.ui.lineEdit_23_list[i_subplot].text()
            self.ui.subplot_ylo_list[i_subplot] = self.ui.lineEdit_24_list[i_subplot].text()
            self.ui.subplot_yhi_list[i_subplot] = self.ui.lineEdit_25_list[i_subplot].text()
            self.ui.subplot_xtick_list[i_subplot] = self.ui.checkBox_26_list[i_subplot].isChecked()
            self.ui.subplot_ytick_list[i_subplot] = self.ui.checkBox_27_list[i_subplot].isChecked()
            self.ui.subplot_xlabel_list[i_subplot] = self.ui.checkBox_28_list[i_subplot].isChecked()
            self.ui.subplot_ylabel_list[i_subplot] = self.ui.checkBox_29_list[i_subplot].isChecked()
        #element
        for i_elmt in range(0, self.ui.num_elmts):
            i_elmt_name = self.ui.lineEdit_31_list[i_elmt].text()
            self.ui.dos_mode_dict[i_elmt_name] = []
            if self.ui.checkBox_32_list[i_elmt].isChecked() == True:
                self.ui.dos_mode_dict[i_elmt_name].append('s')
            if self.ui.checkBox_33_list[i_elmt].isChecked() == True:
                self.ui.dos_mode_dict[i_elmt_name].append('p')
            if self.ui.checkBox_34_list[i_elmt].isChecked() == True:
                self.ui.dos_mode_dict[i_elmt_name].append('d')
            if self.ui.checkBox_35_list[i_elmt].isChecked() == True:
                self.ui.dos_mode_dict[i_elmt_name].append('LDOS')

        self.ui.atom_doscar_indx_list = self.ui.atom_doscar_indx_list[0:self.ui.num_atoms]
        self.ui.atom_doscar_file_path_list = self.ui.atom_doscar_file_path_list[0:self.ui.num_atoms]
        self.ui.atom_sysname_list = self.ui.atom_sysname_list[0:self.ui.num_atoms]
        self.ui.atom_name_list = self.ui.atom_name_list[0:self.ui.num_atoms]
        self.ui.atom_palette_list = self.ui.atom_palette_list[0:self.ui.num_atoms]   
        self.ui.atom_subplot_arg_list = self.ui.atom_subplot_arg_list[0:self.ui.num_atoms]

        self.ui.subplot_arg_list = self.ui.subplot_arg_list[0:self.ui.num_subplots]
        self.ui.subplot_xlo_list = self.ui.subplot_xlo_list[0:self.ui.num_subplots]
        self.ui.subplot_xhi_list = self.ui.subplot_xhi_list[0:self.ui.num_subplots]
        self.ui.subplot_ylo_list = self.ui.subplot_ylo_list[0:self.ui.num_subplots]
        self.ui.subplot_yhi_list = self.ui.subplot_yhi_list[0:self.ui.num_subplots]
        self.ui.subplot_xtick_list = self.ui.subplot_xtick_list[0:self.ui.num_subplots]
        self.ui.subplot_ytick_list = self.ui.subplot_ytick_list[0:self.ui.num_subplots]
        self.ui.subplot_xlabel_list = self.ui.subplot_xlabel_list[0:self.ui.num_subplots]   
        self.ui.subplot_ylabel_list = self.ui.subplot_ylabel_list[0:self.ui.num_subplots]

        self.ui.dos_mode_dict = self.ui.dos_mode_dict
        #General
        self.ui.fermi_shift_zero = self.ui.checkBox_40.isChecked()
        self.ui.fig_width = float(self.ui.lineEdit_41.text())
        self.ui.fig_height = float(self.ui.lineEdit_42.text())
        self.ui.fig_size = [float(self.ui.fig_width), float(self.ui.fig_height)]
        self.ui.fig_dpi = float(self.ui.lineEdit_43.text())
        self.ui.fig_format = self.ui.fig_format
        self.ui.mainplot_xlabel = self.ui.checkBox_45.isChecked()
        self.ui.mainplot_ylabel = self.ui.checkBox_46.isChecked()
        self.ui.share_x_axis = self.ui.checkBox_47.isChecked()
        self.ui.share_y_axis = self.ui.checkBox_48.isChecked()
        
        self.ui.xtick_direction = self.ui.xtick_direction
        self.ui.ytick_direction = self.ui.ytick_direction
        self.ui.peak_analyzer = funcs.convert_bool(self.ui.peak_analyzer)
        self.ui.peak_analyzer_factor = float(self.ui.lineEdit_52.text())
        self.ui.smoothing = funcs.convert_bool(self.ui.smoothing)
        self.ui.smoothing_factor = float(self.ui.lineEdit_54.text())
        self.ui.line_width = float(self.ui.lineEdit_55.text())
        self.ui.font_size = float(self.ui.lineEdit_56.text())            

        vasp_plot.plot_dos(
            atom_doscar_file_path_list = self.ui.atom_doscar_file_path_list,
            atom_sysname_list = self.ui.atom_sysname_list,
            atom_indx_list = self.ui.atom_name_list,
            atom_palette_list = self.ui.atom_palette_list,
            atom_subplot_arg_list = self.ui.atom_subplot_arg_list,
            subplot_arg_list = self.ui.subplot_arg_list,
            subplot_xlo_list = self.ui.subplot_xlo_list,
            subplot_xhi_list = self.ui.subplot_xhi_list,
            subplot_ylo_list = self.ui.subplot_ylo_list,
            subplot_yhi_list = self.ui.subplot_yhi_list,
            subplot_xtick_list = self.ui.subplot_xtick_list,
            subplot_ytick_list = self.ui.subplot_ytick_list,
            subplot_xlabel_list = self.ui.subplot_xlabel_list,
            subplot_ylabel_list = self.ui.subplot_ylabel_list,
	    subplot_share_xy_list = [self.ui.share_x_axis, self.ui.share_y_axis] , 
	    mainplot_axis_label_list = [self.ui.mainplot_xlabel, self.ui.mainplot_ylabel], 
	    xtick_direction = self.ui.xtick_direction, 
	    ytick_direction = self.ui.ytick_direction,
            dos_mode_dict = self.ui.dos_mode_dict, 
	    fermi_shift_zero = self.ui.fermi_shift_zero, 
	    peak_analyzer = self.ui.peak_analyzer, 
	    peak_analyzer_factor = self.ui.peak_analyzer_factor, 
	    smoothing = self.ui.smoothing, 
	    smoothing_factor = self.ui.smoothing_factor, 
	    line_width = self.ui.line_width, 
	    font_size = self.ui.font_size, 
	    fig_format = self.ui.fig_format, 
	    fig_size = self.ui.fig_size, 
	    fig_dpi = self.ui.fig_dpi,
            )

    #vasp_analyze
    def Show_vasp_analyzeWindow(self):        
        self.vasp_analyzeWindow = QtWidgets.QMainWindow()
        self.ui = Ui_vasp_analyzeWindow()
        self.ui.setupUi(self.vasp_analyzeWindow)
        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_nn_map)
        self.ui.actionRun_1.triggered.connect(self.on_run_module_button_clicked_estruct)

        self.vasp_analyzeWindow.show()


    def on_run_module_button_clicked_nn_map(self):
        # nn_map
        self.ui.nn_mode = self.ui.nn_mode
        self.ui.a0 = float(self.ui.lineEdit_b1.text())
        self.ui.n_shell = int(self.ui.lineEdit_b2.text())
        self.ui.common_neighbor_elmt_list = self.ui.lineEdit_b3.text()
        self.ui.poscar_file_path_b0 = self.ui.lineEdit_b7.text()
        
        if ',' in self.ui.lineEdit_b3.text():
            self.ui.common_neighbor_elmt_list = funcs.split_line(line = self.ui.lineEdit_b3.text(), separator = ',')
        else:
            self.ui.common_neighbor_elmt_list = funcs.split_line(line = self.ui.lineEdit_b3.text(), separator = ' ')

        if self.ui.nn_mode == 'NN map':
            vasp_analyze.nn_map(
                poscar_file_path = self.ui.poscar_file_path_b0,
                a0 = self.ui.a0,
                n_shell = self.ui.n_shell,
                )
        elif self.ui.nn_mode == 'simple CNA':
            vasp_analyze.simple_cna(
                poscar_file_path = self.ui.poscar_file_path_b0,
                a0 = self.ui.a0,
                common_neighbor_elmt_list = self.ui.common_neighbor_elmt_list,
                )

    def on_run_module_button_clicked_estruct(self):
        #estruct
        self.ui.sysname = self.ui.lineEdit_c1.text()
        self.ui.doscar_file_path = self.ui.lineEdit_c7.text()
        
        vasp_analyze.estruct(
            doscar_file_path = self.ui.doscar_file_path,
            sysname = self.ui.sysname,
            )
        

    #vasp_write
    def Show_vasp_writeWindow(self):
        self.vasp_writeWindow = QtWidgets.QMainWindow()
        self.ui = Ui_vasp_writeWindow()
        self.ui.setupUi(self.vasp_writeWindow)

        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_write_poscar_with_force)

        self.vasp_writeWindow.show()

    def on_run_module_button_clicked_write_poscar_with_force(self):
        #write_poscar_with_force
        self.ui.ionic_step = self.ui.lineEdit_d1.text()
        self.ui.output_poscar_file_name = self.ui.lineEdit_d2.text()
        self.ui.outcar_file_path = self.ui.lineEdit_d7.text()

        if self.ui.ionic_step in ['last','first']:
            self.ui.ionic_step = self.ui.ionic_step
        else:
            self.ui.ionic_step = int(self.ui.ionic_step)

        if self.ui.output_poscar_file_name in [None,'None','none']:
            self.ui.output_poscar_file_name = None
        else:
            pass
        
        vasp_write.write_poscar_with_force(
            outcar_file_path = self.ui.outcar_file_path,
            ionic_step = self.ui.ionic_step,
            output_poscar_file_name = self.ui.output_poscar_file_name,
            )

    #apt_plot
    def Show_apt_plotWindow(self):
        self.apt_plotWindow = QtWidgets.QMainWindow()
        self.ui = Ui_apt_plotWindow()
        self.ui.setupUi(self.apt_plotWindow)

        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_apt_plot)
##        self.ui.actionRun_1.tri ggered.connect(self.on_run_module_button_clicked_plot_doscar)

        self.apt_plotWindow.show()

    def on_run_module_button_clicked_apt_plot(self):
        self.ui.sysname_e0 = self.ui.lineEdit_e1.text()
        self.ui.proxigram_csv_file_path = self.ui.lineEdit_e8.text()
        self.ui.visible_elmt_list_e0 = [item.strip(' ') for item in self.ui.lineEdit_e2.text().strip().split(',')]          
        self.ui.interpolation_on = self.ui.interpolation_on
        self.ui.fig_width_e0 = float(self.ui.lineEdit_e4.text())
        self.ui.fig_height_e0 = float(self.ui.lineEdit_e5.text())
        self.ui.fig_dpi_e0 = float(self.ui.lineEdit_e6.text())
        self.ui.fig_format_e0 = self.ui.fig_format_e0            

        apt_plot.plot_proxigram_csv(
            proxigram_csv_file_path = self.ui.proxigram_csv_file_path,
            sysname = self.ui.sysname_e0,
            visible_elmt_list = self.ui.visible_elmt_list_e0,
            interplation_on = self.ui.interpolation_on,
            fig_width = self.ui.fig_width_e0,
            fig_height = self.ui.fig_height_e0,
            fig_dpi = self.ui.fig_dpi_e0,
            fig_format = self.ui.fig_format_e0,
            )

    # dvm_analyze
    def Show_dvm_analyzeWindow(self):
        self.dvm_analyzeWindow = QtWidgets.QMainWindow()
        self.ui = Ui_dvm_analyzeWindow()
        self.ui.setupUi(self.dvm_analyzeWindow)

        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_dvm_analyze)
        
        self.dvm_analyzeWindow.show()

    def on_run_module_button_clicked_dvm_analyze(self):
        self.ui.dvm_get_ie_mode = self.ui.dvm_get_ie_mode
        self.ui.dvm_otput_file_path = self.ui.lineEdit_f7.text()
        self.ui.latt_const_f0 = float(self.ui.lineEdit_f1.text())

        if self.ui.dvm_get_ie_mode == 'IE for 1NN atom pairs':
            dvm_analyze.ie_nn(
                dvm_otput_file_path = self.ui.dvm_otput_file_path,
                a0 = self.ui.latt_const_f0,
                )
        elif self.ui.dvm_get_ie_mode == 'IE for all the atom pairs':
            dvm_write.write_ie(
                dvm_otput_file_path = self.ui.dvm_otput_file_path,
                )



    def Show_run_py_scriptWindow(self):
        self.run_py_scriptWindow = QtWidgets.QMainWindow()
        self.ui = Ui_run_py_scriptWindow()
        self.ui.setupUi(self.run_py_scriptWindow)

        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_run_py_script)
##        self.ui.actionRun_1.triggered.connect(self.on_run_module_button_clicked_plot_doscar)

        self.run_py_scriptWindow.show()

    def on_run_module_button_clicked_run_py_script(self):
        self.ui.script_file = self.ui.lineEdit_g7.text()
        if os.path.isfile(self.ui.script_file) == True:
            with open(self.ui.script_file,'r') as f:
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
            funcs.write_log(logfile, self.ui.script_file + " doesn't exist, please check current folder")
        funcs.write_log(logfile,'## read_py_script complete\n')


##app = QtWidgets.QApplication(sys.argv)
##Controller = Controller()
##Controller.Show_MainWindow()
##
####    MainWindow = QtWidgets.QMainWindow()
####    ui = Ui_MainWindow()
####    ui.setupUi(MainWindow)
####    MainWindow.show()
##sys.exit(app.exec_())

time_end = time.time()
formatted_time = time.strftime('%Y%m%d_%H-%M-%S',time.localtime(time_end))
funcs.write_log(logfile,'## matsdp(GUI) ended at '  + formatted_time + '\n' +
               '## time(' + program_name + ') ' + str(time_end-time_start) + ' s')
        
if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Controller = Controller()
    Controller.Show_MainWindow()

##    MainWindow = QtWidgets.QMainWindow()
##    ui = Ui_MainWindow()
##    ui.setupUi(MainWindow)
##    MainWindow.show()
    sys.exit(app.exec_())
