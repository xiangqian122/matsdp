# -*- coding: utf-8 -*-

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
from matsdp.vasp import vasp_plot

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from functools import partial

defaults_dict = default_params.default_params()
logfile = defaults_dict['logfile']

periodic_table_dict = periodic_table.periodic_tab()
dos_mode = periodic_table_dict['dos_mode']
elmt_color = periodic_table_dict['elmt_color']

class Ui_vasp_plotWindow(object):
    def setupUi(self, vasp_plotWindow):
        vasp_plotWindow.setObjectName("vasp_plotWindow")
        vasp_plotWindow.resize(1200, 800)
        self.centralwidget = QtWidgets.QWidget(vasp_plotWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 1160, 711))
        self.tabWidget.setObjectName("tabWidget")

        ###############################################
        #plot_poscar
        ###############################################

        # plot_poscar
##        self.tab = QtWidgets.QWidget()
        self.tab = QtWidgets.QScrollArea()
        self.tab.setObjectName("tab")
        self.tabWidget.addTab(self.tab, "")
        
        content_widget_a0 = QtWidgets.QWidget()
        self.tab.setWidget(content_widget_a0)
        self.flay_a0 = QtWidgets.QFormLayout(content_widget_a0)
        self.tab.setWidgetResizable(True)
        
        self.widget_a0 = QtWidgets.QWidget(self.tab)
        self.widget_a0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_a0.setObjectName("widget")
        self.flay_a0.addRow(self.widget_a0)

        self.widget_a1 = QtWidgets.QWidget(self.tab)
        self.widget_a1.setGeometry(QtCore.QRect(20, 90, 290, 450))
        self.widget_a1.setObjectName("widget")
        self.flay_a0.addRow(self.widget_a1)

        self.widget_a2 = QtWidgets.QWidget(self.tab)
        self.widget_a2.setGeometry(QtCore.QRect(20, 100, 290, 450))
        self.widget_a2.setObjectName("widget")
        self.flay_a0.addRow(self.widget_a2)

        self.widget_a3 = QtWidgets.QWidget(self.tab)
        self.widget_a3.setGeometry(QtCore.QRect(20, 100, 290, 450))
        self.widget_a3.setObjectName("widget")
        self.flay_a0.addRow(self.widget_a3)
        
        #parameters
        self.plot_poscar_mode = 'single'
        self.poscar_or_contcar = 'POSCAR'
        self.num_elmts_a0 = 1
        self.num_elmts_max_a0 = 30
        self.elmt_color_dict_a0 = periodic_table_dict['elmt_color'].copy()

        self.poscar_file_path_a0 = None
        self.euler_angle_type = 'zyx'
        self.phi = -3
        self.theta = 5
        self.psi = 0
        self.draw_mirror_atom = True
        self.box_on = True
        self.axis_indicator = True
        self.plot_cell_basis_vector_label = False
        self.plot_atom_label = None
        self.fig_format_a0 = 'png'
        self.fig_dpi_a0 = 100
        self.draw_colormap = False
        self.colormap_column_indx = 1
        self.colormap_vmin = None
        self.colormap_vmax = None
        self.vmin_color = 'blue'
        self.vmax_color = 'red'
        self.colorbar_alignment = 'vertical'

        self.user_defined_elmt_color = False
        
        self.workdir = None
        
        self.fig_format_list_a0 = ['png','eps','pdf','tif','tiff','jpg','jepg','svg','svgz','pgf','ps','raw','rgba']

        self.element_name_list = list(periodic_table_dict['name'].keys())

        #Initialize the groups
        self.build_group_a0(vasp_plotWindow)
        self.build_group_a1(vasp_plotWindow)
        self.build_group_a2(vasp_plotWindow)
        self.build_group_a3(vasp_plotWindow)

        ###############################################
        #plot_doscar
        ###############################################

              
        # plot_doscar
        
##        self.tab_2 = QtWidgets.QWidget()
        self.tab_2 = QtWidgets.QScrollArea()
        self.tabWidget.addTab(self.tab_2, "")
        content_widget = QtWidgets.QWidget()
        self.tab_2.setWidget(content_widget)
        self.flay = QtWidgets.QFormLayout(content_widget)
        self.tab_2.setWidgetResizable(True)

        self.tab_2.setObjectName("tab_2")
        self.groupBox = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox.setGeometry(QtCore.QRect(20, 10, 800, 61))
        self.groupBox.setObjectName("groupBox")
        self.flay.addRow(self.groupBox)

        self.widget = QtWidgets.QWidget(self.tab_2)
        self.widget.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget.setObjectName("widget")
        self.flay.addRow(self.widget)
        
        self.widget_2 = QtWidgets.QWidget(self.tab_2)
        self.widget_2.setGeometry(QtCore.QRect(300, 70, 500, 420))
        self.widget_2.setObjectName("widget_2")
        self.flay.addRow(self.widget_2)
        
        self.widget_3 = QtWidgets.QWidget(self.tab_2)
        self.widget_3.setGeometry(QtCore.QRect(20, 500, 700, 420))
        self.widget_3.setObjectName("widget_3")
        self.flay.addRow(self.widget_3)

        self.widget_4 = QtWidgets.QWidget(self.tab_2)
        self.widget_4.setGeometry(QtCore.QRect(700, 500, 400, 420))
        self.widget_4.setObjectName("widget_4")
        self.flay.addRow(self.widget_4)

        self.widget_5 = QtWidgets.QWidget(self.tab_2)
        self.widget_5.setGeometry(QtCore.QRect(700, 500, 400, 420))
        self.widget_5.setObjectName("widget_5")
        self.flay.addRow(self.widget_5)

        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget0")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout0")

        self.max_num_doscar = 30
        self.max_num_atoms = 30
        self.max_num_subplots = 30
        self.max_num_elmts = 30
        self.initial_num_doscar = 1
        self.initial_num_atoms = 1
        self.initial_num_subplots = 1
        self.initial_num_elmts = 1
        self.num_doscar = self.initial_num_doscar
        self.num_atoms = self.initial_num_atoms
        self.num_subplots = self.initial_num_subplots
        self.num_elmts = self.initial_num_elmts

        #DOSCAR
        self.doscar_list = [None] * self.max_num_doscar
        self.temp_doscar_list = [None] * self.max_num_doscar
        #Atoms        
        self.atom_doscar_indx_list = [None] * self.max_num_atoms
        self.atom_doscar_file_path_list = [None] * self.max_num_atoms
        self.atom_sysname_list = [None] * self.max_num_atoms
        self.atom_name_list = [None] * self.max_num_atoms
        self.atom_palette_list = [None] * self.max_num_atoms     
        self.atom_subplot_arg_list = [None] * self.max_num_atoms
        #Subplots
        self.subplot_arg_list = [None] * self.max_num_subplots
        self.subplot_xlo_list = [None] * self.max_num_subplots
        self.subplot_xhi_list = [None] * self.max_num_subplots
        self.subplot_ylo_list = [None] * self.max_num_subplots
        self.subplot_yhi_list = [None] * self.max_num_subplots
        self.subplot_xtick_list = [None] * self.max_num_subplots
        self.subplot_ytick_list = [None] * self.max_num_subplots
        self.subplot_xlabel_list = [None] * self.max_num_subplots
        self.subplot_ylabel_list = [None] * self.max_num_subplots
        #element
        self.dos_mode_dict = dos_mode
##        self.dos_mode_elmtname_list = [None] *self.max_num_elmts
##        self.s_dos_mode_checked_list = [None] *self.max_num_elmts
##        self.p_dos_mode_checked_list = [None] *self.max_num_elmts
##        self.d_dos_mode_checked_list = [None] *self.max_num_elmts
##        self.l_dos_mode_checked_list = [None] *self.max_num_elmts
##        self.s_dos_mode = True
##        self.p_dos_mode = True
##        self.d_dos_mode = True
##        self.l_dos_mode = False
        #General
        self.fermi_shift_zero = True
        self.fig_width = 15
        self.fig_height = 10
        self.fig_size = [self.fig_width, self.fig_height]
        self.fig_dpi = 600
        self.fig_format_list = ['png','eps','pdf','tif','tiff','jpg','jepg','svg','svgz','pgf','ps','raw','rgba']
        self.fig_format = 'png'
        self.mainplot_xlabel = True
        self.mainplot_ylabel = True
        self.share_x_axis = False
        self.share_y_axis = False
        self.xtick_direction = 'out'
        self.ytick_direction = 'out'
        self.peak_analyzer = False
        self.peak_analyzer_factor = 0.02
        self.smoothing = False
        self.smoothing_factor = 0.05
        self.line_width = 2.0
        self.font_size = 18

        self.label = QtWidgets.QLabel(self.groupBox)
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label.setObjectName("label")
        self.label.setText("num_doscar")
        self.label.setStatusTip("num_doscar: Number of DOSCAR files")
        self.spinBox = QtWidgets.QSpinBox(self.groupBox)
        self.gridLayout.addWidget(self.spinBox, 0, 1, 1, 1)
        self.spinBox.setObjectName("spinBox")
        self.spinBox.setRange(1,30)
        self.spinBox.setValue(1)
        self.spinBox.valueChanged.connect(self.build_group3)
        
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.gridLayout.addWidget(self.label_2, 0, 2, 1, 1)
        self.label_2.setObjectName("label_2")
        self.label_2.setText("num_atoms")
        self.label_2.setStatusTip("num_atoms: Number of selected atoms for plotting DOS curves")
        self.spinBox_2 = QtWidgets.QSpinBox(self.groupBox)
        self.gridLayout.addWidget(self.spinBox_2, 0, 3, 1, 1)
        self.spinBox_2.setObjectName("spinBox_2")
        self.spinBox_2.setRange(1,30)
        self.spinBox_2.setValue(1)
##        self.spinBox_2.setMaximum(30)
##        num_atoms = self.spinBox_2.value()
        self.spinBox_2.valueChanged.connect(self.build_group4)
##        self.spinBox_2.valueChanged[int].connect(self.set_num_atoms)
##        vasp_plotWindow.update()


        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.gridLayout.addWidget(self.label_3, 0, 4, 1, 1)
        self.label_3.setObjectName("label_3")
        self.label_3.setText("num_subplots")
        self.label_3.setStatusTip("num_subplots: Number of subplots")
        self.spinBox_3 = QtWidgets.QSpinBox(self.groupBox)
        self.gridLayout.addWidget(self.spinBox_3, 0, 5, 1, 1)
        self.spinBox_3.setObjectName("spinBox_3")
        self.spinBox_3.setRange(1,30)
        self.spinBox_3.setValue(1)
        self.spinBox_3.valueChanged.connect(self.build_group5)
        
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.gridLayout.addWidget(self.label_4, 0, 6, 1, 1)
        self.label_4.setObjectName("label_4")
        self.label_4.setText("num_elements")
        self.label_4.setStatusTip("num_elements: Number of elements")
        self.spinBox_4 = QtWidgets.QSpinBox(self.groupBox)
        self.gridLayout.addWidget(self.spinBox_4, 0, 7, 1, 1)
        self.spinBox_4.setObjectName("spinBox_4")
        self.spinBox_4.setRange(1,30)
        self.spinBox_4.setValue(1)
        self.spinBox_4.valueChanged.connect(self.build_group2)
        
        self.groupBox.setLayout(self.gridLayout)
        self.groupBox.setTitle("control")

        
        #Initialize the groups
        self.build_group3(vasp_plotWindow)
        self.build_group4(vasp_plotWindow)
        self.build_group5(vasp_plotWindow)
        self.build_group2(vasp_plotWindow)
        self.build_group6(vasp_plotWindow)

        
        vasp_plotWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(vasp_plotWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
##        self.menuHelp = QtWidgets.QMenu(self.menubar)
##        self.menuHelp.setObjectName("menuHelp")
        vasp_plotWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(vasp_plotWindow)
        self.statusbar.setObjectName("statusbar")
        vasp_plotWindow.setStatusBar(self.statusbar)
        self.actionRun = QtWidgets.QAction(vasp_plotWindow)
        self.actionRun.setObjectName("actionRun")
        self.action_1 = QtWidgets.QAction(vasp_plotWindow)
        self.action_1.setObjectName("action_1")
        self.action_2 = QtWidgets.QAction(vasp_plotWindow)
        self.action_2.setObjectName("action_2")
##        self.actionHelp = QtWidgets.QAction(vasp_plotWindow)
##        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_1 = QtWidgets.QAction(vasp_plotWindow)
        self.actionRun_1.setObjectName("actionRun_1")
        self.menuRun.addAction(self.actionRun)
        self.menuRun.addAction(self.actionRun_1)
        self.menubar.addAction(self.menuRun.menuAction())
##        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(vasp_plotWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(vasp_plotWindow)

    def retranslateUi(self, vasp_plotWindow):
        _translate = QtCore.QCoreApplication.translate
        vasp_plotWindow.setWindowTitle(_translate("vasp_plotWindow", "vasp_plot"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("vasp_plotWindow", "plot_poscar"))


        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("vasp_plotWindow", "plot_doscar"))
        self.menuRun.setTitle(_translate("vasp_plotWindow", "Run"))
##        self.menuHelp.setTitle(_translate("vasp_plotWindow", "Help"))
        self.actionRun.setText(_translate("vasp_plotWindow", "Run plot_poscar"))
        self.actionRun.setStatusTip(_translate("vasp_plotWindow", "Plot POSCAR"))
        self.action_1.setText(_translate("vasp_plotWindow", "plot_poscar"))
        self.action_2.setText(_translate("vasp_plotWindow", "plot_doscar"))
##        self.actionHelp.setText(_translate("vasp_plotWindow", "Help"))
        self.actionRun_1.setText(_translate("vasp_plotWindow", "Run plot_doscar"))
        self.actionRun_1.setStatusTip(_translate("vasp_plotWindow", "Plot DOSCAR"))

    #####################
    # functions
    #####################

    ## plot_poscar
    def build_group_a0(self,vasp_plotWindow):
        self.widget_a0 = QtWidgets.QWidget(self.tab)
        self.widget_a0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_a0.setObjectName("widget")
        self.flay_a0.removeRow(1)
        self.flay_a0.insertRow(1, self.widget_a0)

        self.groupBox_a0 = QtWidgets.QGroupBox(self.widget_a0)
        self.groupBox_a0.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_a0.setObjectName("groupBox_a0")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_a0)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
            
        self.label_a0 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a0, 0, 0, 1, 1)
        self.label_a0.setObjectName("label_a0")
        self.comboBox_a0 = QtWidgets.QComboBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.comboBox_a0, 0, 1, 1, 1)
        self.comboBox_a0.setObjectName("comboBox")
        self.comboBox_a0.addItem("single")
        self.comboBox_a0.addItem("multiple")
        self.comboBox_a0.activated[str].connect(self.mode_handleActivated)

        self.label_a1 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a1, 1, 0, 1, 1)
        self.label_a1.setObjectName("label_a1")
        self.lineEdit_a1 = QtWidgets.QLineEdit(self.groupBox_a0)
        self.gridLayout.addWidget(self.lineEdit_a1, 1, 1, 1, 1)
        self.lineEdit_a1.setObjectName("lineEdit")

        self.label_a2 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a2, 1, 2, 1, 1)
        self.label_a2.setObjectName("label_a2")
        self.lineEdit_a2 = QtWidgets.QLineEdit(self.groupBox_a0)
        self.gridLayout.addWidget(self.lineEdit_a2, 1, 3, 1, 1)
        self.lineEdit_a2.setObjectName("lineEdit")

        self.label_a3 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a3, 1, 4, 1, 1)
        self.label_a1.setObjectName("label_a3")
        self.lineEdit_a3 = QtWidgets.QLineEdit(self.groupBox_a0)
        self.gridLayout.addWidget(self.lineEdit_a3, 1, 5, 1, 1)
        self.lineEdit_a3.setObjectName("lineEdit")

        self.label_a4 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a4, 1, 6, 1, 1)
        self.label_a4.setObjectName("label_a4")
        self.lineEdit_a4 = QtWidgets.QLineEdit(self.groupBox_a0)
        self.gridLayout.addWidget(self.lineEdit_a4, 1, 7, 1, 1)
        self.lineEdit_a4.setObjectName("lineEdit")

        self.checkBox_a0 = QtWidgets.QCheckBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.checkBox_a0, 2, 0, 1, 1)
        self.checkBox_a0.setObjectName("checkBox_a0")
        self.checkBox_a0.toggle()

        self.checkBox_a1 = QtWidgets.QCheckBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.checkBox_a1, 2, 1, 1, 1)
        self.checkBox_a1.setObjectName("checkBox_a1")
        self.checkBox_a1.toggle()

        self.checkBox_a2 = QtWidgets.QCheckBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.checkBox_a2, 2, 2, 1, 1)
        self.checkBox_a2.setObjectName("checkBox_a2")
        self.checkBox_a2.toggle()

        self.checkBox_a3 = QtWidgets.QCheckBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.checkBox_a3, 2, 3, 1, 1)
        self.checkBox_a3.setObjectName("checkBox_a3")
        self.checkBox_a3.toggle()

        self.label_a44 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a44, 2, 4, 1, 1)
        self.label_a44.setObjectName("label_a44")
        self.lineEdit_a44 = QtWidgets.QLineEdit(self.groupBox_a0)
        self.gridLayout.addWidget(self.lineEdit_a44, 2, 5, 1, 1)
        self.lineEdit_a44.setObjectName("lineEdit")

##        self.checkBox_a4 = QtWidgets.QCheckBox(self.groupBox_a0)
##        self.gridLayout.addWidget(self.checkBox_a4, 2, 4, 1, 1)
##        self.checkBox_a4.setObjectName("checkBox_a4")
##        self.checkBox_a4.toggle()

        self.label_a5 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a5, 3, 0, 1, 1)
        self.label_a5.setObjectName("label_a5")
        self.comboBox_a5 = QtWidgets.QComboBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.comboBox_a5, 3, 1, 1, 1)
        self.comboBox_a5.setObjectName("comboBox")
        self.comboBox_a5.addItems(self.fig_format_list_a0)
        self.comboBox_a5.activated[str].connect(self.fig_format_handleActivated)


        self.label_a6 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a6, 3, 2, 1, 1)
        self.label_a6.setObjectName("label_a6")
        self.lineEdit_a6 = QtWidgets.QLineEdit(self.groupBox_a0)
        self.gridLayout.addWidget(self.lineEdit_a6, 3, 3, 1, 1)
        self.lineEdit_a6.setObjectName("lineEdit")

        self.checkBox_a20 = QtWidgets.QCheckBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.checkBox_a20, 4, 0, 1, 1)
        self.checkBox_a20.setObjectName("checkBox_a20")
        self.checkBox_a20.toggle()
        self.checkBox_a20.toggled.connect(self.build_group_a2)
        self.checkBox_a20.toggled.connect(self.set_enable)
        self.user_defined_elmt_color = self.checkBox_a20.isChecked()

        self.label_a21 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a21, 4, 1, 1, 1)
        self.label_a21.setObjectName("label_a21")
        self.spinBox_a21 = QtWidgets.QSpinBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.spinBox_a21, 4, 2, 1, 1)
        self.spinBox_a21.setObjectName("spinBox_a21")
        self.spinBox_a21.setRange(1,30)
        self.spinBox_a21.setValue(1)
        self.spinBox_a21.valueChanged.connect(self.build_group_a2)
##        self.spinBox_a21.setEnabled(False)

        self.label_a10 = QtWidgets.QLabel(self.groupBox_a0)
        self.gridLayout.addWidget(self.label_a10, 5, 0, 1, 1)
        self.label_a10.setObjectName("label_a10")
        self.comboBox_a10 = QtWidgets.QComboBox(self.groupBox_a0)
        self.gridLayout.addWidget(self.comboBox_a10, 5, 1, 1, 1)
        self.comboBox_a10.setObjectName("comboBox")
        self.comboBox_a10.addItem("False")
        self.comboBox_a10.addItem("True")
        self.comboBox_a10.activated[str].connect(self.build_group_a2)
        self.comboBox_a10.activated[str].connect(self.draw_colormap_handleActivated)




        self.label_a0.setText('mode')
        self.label_a1.setText('euler_angle_type')
        self.label_a2.setText('phi')
        self.label_a3.setText('theta')
        self.label_a4.setText('psi')
        self.lineEdit_a1.setText("zyz")
        self.lineEdit_a1.setStatusTip('string of length 3. It specify the type of rotations based on Eulerian angles. Choices are zyz, zxz, zyx, etc.. Usually the zyz type is used.')
        self.lineEdit_a2.setText("-3")
        self.lineEdit_a2.setStatusTip('phi, theta, psi: float formats. The first, second, and third rotation Eulerian angles, units in degrees')
        self.lineEdit_a3.setText("4")
        self.lineEdit_a3.setStatusTip('phi, theta, psi: float formats. The first, second, and third rotation Eulerian angles, units in degrees')
        self.lineEdit_a4.setText("0")
        self.lineEdit_a4.setStatusTip('phi, theta, psi: float formats. The first, second, and third rotation Eulerian angles, units in degrees')
        self.checkBox_a0.setText("draw_mirror_atom")
        self.checkBox_a0.setStatusTip("draw_mirror_atom: Whether to plot the mirror atoms at the periodic boundary")
        self.checkBox_a1.setText("box_on")
        self.checkBox_a1.setStatusTip("box_on: Whether to plot the box or not")
        self.checkBox_a2.setText("axis_indicator")
        self.checkBox_a2.setStatusTip("axis_indicator: Whethether to plot the axis indicator or not")
        self.checkBox_a3.setText("cell_basis_vector")
        self.checkBox_a3.setStatusTip("cell_basis_vector: Whether to plot the cell basis vector labels( i.e., to label the three basis vectors of the cell as a, b, and c")
        self.label_a44.setText('atom_label')
        self.lineEdit_a44.setText('None')
        self.lineEdit_a44.setStatusTip('atom_label: String value. Values: atom_name, atom_index, atom_species, or None. It plot atom label to each atom')

##        self.checkBox_a4.setText("atom_label")
##        self.checkBox_a4.setStatusTip("atom_label: String value. values: 'atom_name', 'atom_index', 'atom_species', or None/'None'. It plot atom label to each atom")
        self.label_a5.setText('fig_format')
        self.label_a5.setStatusTip('A string that defines output figure format')
        self.label_a6.setText('fig_dpi')            
        self.lineEdit_a6.setText("100")
        self.checkBox_a20.setText('user-defined element color')
        self.label_a21.setText('num_elements')
        self.label_a10.setText('color mapping')

        self.groupBox_a0.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_a0)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(240)
        layout = QtWidgets.QVBoxLayout(self.widget_a0)
        layout.addWidget(scroll)
        
        self.groupBox_a0.setTitle("general settings")

    def build_group_a1(self,vasp_plotWindow):
        self.widget_a1 = QtWidgets.QWidget(self.tab)
        self.widget_a1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_a1.setObjectName("widget")
        self.flay_a0.removeRow(2)
        self.flay_a0.insertRow(2, self.widget_a1)

        self.groupBox_a1 = QtWidgets.QGroupBox(self.widget_a1)
        self.groupBox_a1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_a1.setObjectName("groupBox_a1")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_a1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_a7 = QtWidgets.QLabel(self.groupBox_a1)
        self.gridLayout.addWidget(self.label_a7, 0, 0, 1, 1)
        self.label_a7.setObjectName("label_a7")
        self.pushButton_a7 = QtWidgets.QPushButton(self.groupBox_a1)
        self.gridLayout.addWidget(self.pushButton_a7, 0, 1, 1, 1)
        self.pushButton_a7.setObjectName("pushButton_a7")
        self.pushButton_a7.clicked.connect(self.on_pushButton_clicked_a0)        
        self.lineEdit_a7 = QtWidgets.QLineEdit(self.groupBox_a1)
        self.gridLayout.addWidget(self.lineEdit_a7, 0, 2, 1, 1)
        self.lineEdit_a7.setObjectName("lineEdit_a7")

        self.label_a8 = QtWidgets.QLabel(self.groupBox_a1)
        self.gridLayout.addWidget(self.label_a8, 1, 0, 1, 1)
        self.label_a8.setObjectName("label_a8")
        self.comboBox_a8 = QtWidgets.QComboBox(self.groupBox_a1)
        self.gridLayout.addWidget(self.comboBox_a8, 1, 1, 1, 1)
        self.comboBox_a8.setObjectName("comboBox")
        self.comboBox_a8.addItem('POSCAR')
        self.comboBox_a8.addItem('CONTCAR')
        self.comboBox_a8.activated[str].connect(self.poscar_or_contcar_handleActivated)
        self.comboBox_a8.setEnabled(False)
        self.comboBox_a8.setStatusTip('Plot the POSCAR or CONTCAR?')
        
        self.label_a7.setText('single POSCAR file')
        self.pushButton_a7.setText("...")
        self.pushButton_a7.setStatusTip('Choose the POSCAR file')
        self.lineEdit_a7.setText("None")
        self.lineEdit_a7.setStatusTip('Enter the POSCAR file path (relative path is also accepted)')
        self.label_a8.setText('poscar_or_contcar')
        self.label_a8.setStatusTip('Choose the POSCAR or CONTCAR as the input')

        self.groupBox_a1.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_a1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(100)
        layout = QtWidgets.QVBoxLayout(self.widget_a1)
        layout.addWidget(scroll)
        
        self.groupBox_a1.setTitle("file importer")



    def build_group_a2(self,vasp_plotWindow):
        self.widget_a2 = QtWidgets.QWidget(self.tab)
        self.widget_a2.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_a2.setObjectName("widget_a2")
        self.flay_a0.removeRow(3)
        self.flay_a0.insertRow(3, self.widget_a2)

        self.groupBox_a2 = QtWidgets.QGroupBox(self.widget_a2)
        self.groupBox_a2.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_a2.setObjectName("groupBox_a2")
     
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_a2)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.num_elmts_a0 = self.spinBox_a21.value()

        self.lineEdit_a22_list = [None] * self.num_elmts_a0
        self.lineEdit_a23_list = [None] * self.num_elmts_a0

##        self.comboBox_22_list = [None] * self.num_elmts_a0
        
        for i_elmt in range(0, self.num_elmts_a0):
            self.label_a22 = QtWidgets.QLabel(self.groupBox_a2)
            self.gridLayout.addWidget(self.label_a22, 2 + i_elmt, 0, 1, 1)
            self.label_a22.setObjectName("label_a22")
            self.label_a22.setText('element')

##            self.comboBox_a22 = QtWidgets.QComboBox(self.groupBox_6)
##            self.gridLayout.addWidget(self.comboBox_a22, 1, 7, 1, 1)
##            self.comboBox_a22.setObjectName("comboBox")
##            self.comboBox_a22.addItems(self.element_name_list)
##    ##        self.comboBox_a22.addItem("png")
##            self.comboBox_a22.activated[str].connect(self.fig_format_handleActivated)
            
            self.lineEdit_a22_list[i_elmt] = QtWidgets.QLineEdit(self.groupBox_a2)
            self.gridLayout.addWidget(self.lineEdit_a22_list[i_elmt], 2 + i_elmt, 1, 1, 1)
            self.lineEdit_a22_list[i_elmt].setObjectName("lineEdit_a22")
            self.lineEdit_a22_list[i_elmt].setText('Ni')

            self.label_a23 = QtWidgets.QLabel(self.groupBox_a2)
            self.gridLayout.addWidget(self.label_a23, 2 + i_elmt, 2, 1, 1)
            self.label_a23.setObjectName("label_a23")
            self.label_a23.setText('color')
            
            self.lineEdit_a23_list[i_elmt] = QtWidgets.QLineEdit(self.groupBox_a2)
            self.gridLayout.addWidget(self.lineEdit_a23_list[i_elmt], 2 + i_elmt, 3, 1, 1)
            self.lineEdit_a23_list[i_elmt].setObjectName("lineEdit_a23")
            self.lineEdit_a23_list[i_elmt].setText('gray')

        
        
##        self.pushButton_a8.setText("...")
##        self.lineEdit_a8.setText("None")
##        self.label_a9.setText('poscar_or_contcar')
            

        self.groupBox_a2.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_a2)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(120)
        layout = QtWidgets.QVBoxLayout(self.widget_a2)
        layout.addWidget(scroll)
        
        self.groupBox_a2.setTitle("element color")


    def build_group_a3(self,vasp_plotWindow):
        self.widget_a3 = QtWidgets.QWidget(self.tab)
        self.widget_a3.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_a3.setObjectName("widget_a3")
        self.flay_a0.removeRow(4)
        self.flay_a0.insertRow(4, self.widget_a3)

        self.groupBox_a3 = QtWidgets.QGroupBox(self.widget_a3)
        self.groupBox_a3.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_a3.setObjectName("groupBox_a3")
     
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_a3)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_a31 = QtWidgets.QLabel(self.groupBox_a3)
        self.gridLayout.addWidget(self.label_a31, 0, 0, 1, 1)
        self.label_a31.setObjectName("label_a31")
        self.lineEdit_a31 = QtWidgets.QLineEdit(self.groupBox_a3)
        self.gridLayout.addWidget(self.lineEdit_a31, 0, 1, 1, 1)
        self.lineEdit_a31.setObjectName("lineEdit")

        self.label_a32 = QtWidgets.QLabel(self.groupBox_a3)
        self.gridLayout.addWidget(self.label_a32, 1, 0, 1, 1)
        self.label_a32.setObjectName("label_a32")
        self.lineEdit_a32 = QtWidgets.QLineEdit(self.groupBox_a3)
        self.gridLayout.addWidget(self.lineEdit_a32, 1, 1, 1, 1)
        self.lineEdit_a32.setObjectName("lineEdit")

        self.label_a33 = QtWidgets.QLabel(self.groupBox_a3)
        self.gridLayout.addWidget(self.label_a33, 1, 2, 1, 1)
        self.label_a33.setObjectName("label_a33")
        self.lineEdit_a33 = QtWidgets.QLineEdit(self.groupBox_a3)
        self.gridLayout.addWidget(self.lineEdit_a33, 1, 3, 1, 1)
        self.lineEdit_a33.setObjectName("lineEdit")

        self.label_a34 = QtWidgets.QLabel(self.groupBox_a3)
        self.gridLayout.addWidget(self.label_a34, 1, 4, 1, 1)
        self.label_a34.setObjectName("label_a34")
        self.lineEdit_a34 = QtWidgets.QLineEdit(self.groupBox_a3)
        self.gridLayout.addWidget(self.lineEdit_a34, 1, 5, 1, 1)
        self.lineEdit_a34.setObjectName("lineEdit")

        self.label_a35 = QtWidgets.QLabel(self.groupBox_a3)
        self.gridLayout.addWidget(self.label_a35, 1, 6, 1, 1)
        self.label_a35.setObjectName("label_a35")
        self.lineEdit_a35 = QtWidgets.QLineEdit(self.groupBox_a3)
        self.gridLayout.addWidget(self.lineEdit_a35, 1, 7, 1, 1)
        self.lineEdit_a35.setObjectName("lineEdit")

        self.label_a36 = QtWidgets.QLabel(self.groupBox_a3)
        self.gridLayout.addWidget(self.label_a36, 2, 0, 1, 1)
        self.label_a36.setObjectName("label_a36")
        self.comboBox_a36 = QtWidgets.QComboBox(self.groupBox_a3)
        self.gridLayout.addWidget(self.comboBox_a36, 2, 1, 1, 1)
        self.comboBox_a36.setObjectName("comboBox")
        self.comboBox_a36.addItem('vertical')
        self.comboBox_a36.addItem('horizontal')
        self.comboBox_a36.activated[str].connect(self.colorbar_alignment_handleActivated)

        self.label_a31.setText("colormap_column_indx")
        self.lineEdit_a31.setText("1")
        self.label_a32.setText("colormap_vmin")
        self.lineEdit_a32.setText("None")
        self.lineEdit_a32.setStatusTip('Minimum value in the colorbar. If None is set, the minimum value in the colorbar is set automatically')
        self.label_a33.setText("colormap_vmax")
        self.lineEdit_a33.setText("None")
        self.lineEdit_a33.setStatusTip('Maximum value in the colorbar. If None is set, the maximum value in the colorbar is set automatically')
        self.label_a34.setText("vmin_color")
        self.label_a34.setStatusTip('vmix_color: The color for the minimum value end in the colorbar')
        self.lineEdit_a34.setText("blue")
        self.lineEdit_a34.setStatusTip('The color for the minimum value end in the colorbar')
        self.label_a35.setText("vmax_color")
        self.label_a35.setStatusTip('vmax_color: The color for the maximum value end in the colorbar')
        self.lineEdit_a35.setText("red")
        self.lineEdit_a35.setStatusTip('The color for the maximum value end in the colorbar')
        self.label_a36.setText("colorbar_alignment")
        self.label_a36.setStatusTip('colorbar_alignment: Whether teh colorbar is vertically or horizontally aligned')

        self.label_a31.setEnabled(False)
        self.lineEdit_a31.setEnabled(False)
        self.label_a32.setEnabled(False)
        self.lineEdit_a32.setEnabled(False)
        self.label_a33.setEnabled(False)
        self.lineEdit_a33.setEnabled(False)
        self.label_a34.setEnabled(False)
        self.lineEdit_a34.setEnabled(False)
        self.label_a35.setEnabled(False)
        self.lineEdit_a35.setEnabled(False)
        self.label_a36.setEnabled(False)
        self.comboBox_a36.setEnabled(False)  

        self.groupBox_a3.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_a3)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(150)
        layout = QtWidgets.QVBoxLayout(self.widget_a3)
        layout.addWidget(scroll)
        
        self.groupBox_a3.setTitle("color mapping")

    def mode_handleActivated(self, text):
        self.plot_poscar_mode = text
        if self.plot_poscar_mode == 'single':
            self.label_a7.setText('single POSCAR file')
            self.lineEdit_a7.setText('{}'.format(''))
            self.pushButton_a7.setStatusTip('Choose the POSCAR file')
            self.lineEdit_a7.setStatusTip('Enter the POSCAR file path (relative path is also accepted)')
            self.comboBox_a8.setEnabled(False)
##            self.comboBox_a8.clear()
##            self.comboBox_a8.addItem('POSCAR')
##            self.comboBox_a8.addItem('CONTCAR')
        elif self.plot_poscar_mode == 'multiple':
            self.label_a7.setText('POSCAR work directory')
            self.pushButton_a7.setStatusTip('Choose the directory (folder) which contains the multiple POSCAR/CONTCAR files')
            self.lineEdit_a7.setText('{}'.format(''))
            self.lineEdit_a7.setStatusTip('Enter the directory (folder) which contains multiple VASP calculation jobs (relative path is also accepted)')
            self.comboBox_a8.setEnabled(True)
            self.comboBox_a8.setStatusTip('Plot the POSCAR or CONTCAR?')
##            self.comboBox_a8.clear()
            
    def set_enable(self, text):
        self.user_defined_elmt_color = funcs.convert_bool(text)
        if self.user_defined_elmt_color == True:
            self.spinBox_a21.setEnabled(True)
            self.num_elmts_a0 = self.spinBox_a21.value()
            for i_elmt in range(0, self.num_elmts_a0):
                self.lineEdit_a22_list[i_elmt].setEnabled(True)
                self.lineEdit_a23_list[i_elmt].setEnabled(True)
        elif self.user_defined_elmt_color == False:
            self.spinBox_a21.setEnabled(False)

            self.num_elmts_a0 = self.spinBox_a21.value()
            for i_elmt in range(0, self.num_elmts_a0):
                self.lineEdit_a22_list[i_elmt].setEnabled(False)
                self.lineEdit_a23_list[i_elmt].setEnabled(False)
##                print(self.lineEdit_22_list[i_elmt].isEnabled())


    def poscar_or_contcar_handleActivated(self, text):
        self.poscar_or_contcar = text

    def draw_colormap_handleActivated(self, text):
        self.draw_colormap = funcs.convert_bool(text)

        if self.draw_colormap == True:
            self.label_a31.setEnabled(True)
            self.lineEdit_a31.setEnabled(True)
            self.label_a32.setEnabled(True)
            self.lineEdit_a32.setEnabled(True)
            self.label_a33.setEnabled(True)
            self.lineEdit_a33.setEnabled(True)
            self.label_a34.setEnabled(True)
            self.lineEdit_a34.setEnabled(True)
            self.label_a35.setEnabled(True)
            self.lineEdit_a35.setEnabled(True)
            self.label_a36.setEnabled(True)
            self.comboBox_a36.setEnabled(True)

            self.num_elmts_a0 = self.spinBox_a21.value()
            for i_elmt in range(0, self.num_elmts_a0):
                self.lineEdit_a22_list[i_elmt].setEnabled(False)
                self.lineEdit_a23_list[i_elmt].setEnabled(False)

            
        elif self.draw_colormap == False:
            self.label_a31.setEnabled(False)
            self.lineEdit_a31.setEnabled(False)
            self.label_a32.setEnabled(False)
            self.lineEdit_a32.setEnabled(False)
            self.label_a33.setEnabled(False)
            self.lineEdit_a33.setEnabled(False)
            self.label_a34.setEnabled(False)
            self.lineEdit_a34.setEnabled(False)
            self.label_a35.setEnabled(False)
            self.lineEdit_a35.setEnabled(False)
            self.label_a36.setEnabled(False)
            self.comboBox_a36.setEnabled(False)

            self.num_elmts_a0 = self.spinBox_a21.value()
            for i_elmt in range(0, self.num_elmts_a0):
                self.lineEdit_a22_list[i_elmt].setEnabled(True)
                self.lineEdit_a23_list[i_elmt].setEnabled(True)

            
    def colorbar_alignment_handleActivated(self, text):
        self.colorbar_alignment = text

    def on_pushButton_clicked_a0(self):
        if self.plot_poscar_mode == 'single':
            self.label_a7.setText('single POSCAR file')
            self.open_dialog_box_a0_file()
        elif self.plot_poscar_mode == 'multiple':
            self.label_a7.setText('POSCAR work directory')
            self.open_dialog_box_a0_dir()            
        
    def open_dialog_box_a0_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_a7.setText('{}'.format(file_path))

    def open_dialog_box_a0_dir(self):
        dirname = QFileDialog.getExistingDirectory()
        dir_path = dirname
        self.lineEdit_a7.setText('{}'.format(dir_path))





    ## plot_doscar

    def build_group3(self,vasp_plotWindow):
        self.widget = QtWidgets.QWidget(self.tab_2)
        self.widget.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget.setObjectName("widget")
        self.flay.removeRow(1)
        self.flay.insertRow(1, self.widget)

        self.groupBox_3 = QtWidgets.QGroupBox(self.widget)
        self.groupBox_3.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_3.setObjectName("groupBox_3")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_3)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
            
        self.num_doscar = self.spinBox.value()
        self.pushButton_list = [None] * self.num_doscar
        self.lineEdit_list = [None] * self.num_doscar

        for i_doscar in range(0, self.num_doscar):
            self.label_5 = QtWidgets.QLabel(self.groupBox_3)
            self.gridLayout.addWidget(self.label_5, i_doscar, 0, 1, 1)
            self.label_5.setObjectName("label_5")
            
            self.pushButton_list[i_doscar] = QtWidgets.QPushButton(self.groupBox_3)
            self.gridLayout.addWidget(self.pushButton_list[i_doscar], i_doscar, 1, 1, 1)
            self.pushButton_list[i_doscar].setObjectName("pushButton_" + str(i_doscar))
            
            self.lineEdit_list[i_doscar] = QtWidgets.QLineEdit(self.groupBox_3)
            self.gridLayout.addWidget(self.lineEdit_list[i_doscar], i_doscar, 2, 1, 1)
            self.lineEdit_list[i_doscar].setObjectName("lineEdit")
            
            self.label_5.setText('DOSCAR ' + str(i_doscar + 1))
            self.pushButton_list[i_doscar].setText("...")
            self.pushButton_list[i_doscar].setStatusTip("Choose a DOSCAR file")
            self.lineEdit_list[i_doscar].setText("None")
            self.lineEdit_list[i_doscar].setStatusTip("Enter the DOSCAR file path (relative path is also accepted)")

        for i_doscar in range(0, self.num_doscar):
            self.pushButton_list[i_doscar].clicked.connect(partial(self.on_pushButton_clicked,i_doscar))
##            self.pushButton_list[i_doscar].clicked.connect(self.on_pushButton_clicked)

        self.groupBox_3.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_3)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget)
        layout.addWidget(scroll)

        
        self.groupBox_3.setTitle("DOSCAR")


    def build_group4(self,vasp_plotWindow):
        self.widget_2 = QtWidgets.QWidget(self.tab_2)
        self.widget_2.setGeometry(QtCore.QRect(300, 70, 500, 420))
        self.widget_2.setObjectName("widget_2")
        self.flay.removeRow(2)
        self.flay.insertRow(2, self.widget_2)
        
        self.groupBox_4 = QtWidgets.QGroupBox(self.widget_2)
        self.groupBox_4.setGeometry(QtCore.QRect(310, 80, 361, 141))
        self.groupBox_4.setObjectName("groupBox_4")


        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_4)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        
        self.num_atoms = self.spinBox_2.value()
        self.spinBox_9_list = [None] * self.num_atoms
        self.lineEdit_12_list = [None] * self.num_atoms
        self.lineEdit_13_list = [None] * self.num_atoms
        self.lineEdit_14_list = [None] * self.num_atoms
        self.lineEdit_15_list = [None] * self.num_atoms

        for i_atom in range(0, self.num_atoms):
            self.label_11 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_11, i_atom, 0, 1, 1)
            self.label_11.setObjectName("label_11")
            self.spinBox_9_list[i_atom] = QtWidgets.QSpinBox(self.groupBox_4)
            self.gridLayout.addWidget(self.spinBox_9_list[i_atom], i_atom, 1, 1, 1)
            self.spinBox_9_list[i_atom].setObjectName("spinBox_9")
            self.spinBox_9_list[i_atom].setRange(1,30)
            self.spinBox_9_list[i_atom].setValue(1)
            
            self.label_12 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_12, i_atom, 2, 1, 1)
            self.label_12.setObjectName("label_12")
            self.lineEdit_12_list[i_atom] = QtWidgets.QLineEdit(self.groupBox_4)
            self.gridLayout.addWidget(self.lineEdit_12_list[i_atom], i_atom, 3, 1, 1)
            self.lineEdit_12_list[i_atom].setObjectName("lineEdit_12")
            
            self.label_13 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_13, i_atom, 4, 1, 1)
            self.label_13.setObjectName("label_13")
            self.lineEdit_13_list[i_atom] = QtWidgets.QLineEdit(self.groupBox_4)
            self.gridLayout.addWidget(self.lineEdit_13_list[i_atom], i_atom, 5, 1, 1)
            self.lineEdit_13_list[i_atom].setObjectName("lineEdit_13")

            self.label_14 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_14, i_atom, 6, 1, 1)
            self.label_14.setObjectName("label_14")
            self.lineEdit_14_list[i_atom] = QtWidgets.QLineEdit(self.groupBox_4)
            self.gridLayout.addWidget(self.lineEdit_14_list[i_atom], i_atom, 7, 1, 1)
            self.lineEdit_14_list[i_atom].setObjectName("lineEdit_14")

            self.label_15 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_15, i_atom, 8, 1, 1)
            self.label_15.setObjectName("label_15")
            self.lineEdit_15_list[i_atom] = QtWidgets.QLineEdit(self.groupBox_4)
            self.gridLayout.addWidget(self.lineEdit_15_list[i_atom], i_atom, 9, 1, 1)
            self.lineEdit_15_list[i_atom].setObjectName("lineEdit_15")

            self.label_11.setText("DOSCAR")
            self.label_11.setStatusTip("Choose the corresponding DOSCAR file for the specified atom")
            self.label_12.setText("sysname")
            self.label_12.setStatusTip("User-defined system name, this value will be shown in the legend of the exported figure")
            self.lineEdit_12_list[i_atom].setText("None")
            self.lineEdit_12_list[i_atom].setStatusTip("User-defined system name, this value will be shown in the legend of the exported figure. Default: None")
            self.label_13.setText("atomname")
            self.label_13.setStatusTip('The atom name for the specific DOS curve')
            self.lineEdit_13_list[i_atom].setText("Ni1")
            self.lineEdit_13_list[i_atom].setStatusTip('The atom name for specific DOS curve')
            self.label_14.setText("color")
            self.label_14.setStatusTip('The color of the specific DOS curve')
            self.lineEdit_14_list[i_atom].setText("black")
            self.lineEdit_14_list[i_atom].setStatusTip('The color of the specific DOS curve')
            self.label_15.setText("subplot_arg")
            self.label_15.setStatusTip('The subplot argument for the specific DOS curve, e.g. 111 or 211 or 223 etc.')
            self.lineEdit_15_list[i_atom].setText("111")
            self.lineEdit_15_list[i_atom].setStatusTip('The subplot argument for the specific DOS curve, e.g. 111 or 211 or 223 etc.')
            
        self.groupBox_4.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_4)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget_2)
        layout.addWidget(scroll)

        self.groupBox_4.setTitle("Atom")


    def build_group5(self,vasp_plotWindow):
        self.widget_3 = QtWidgets.QWidget(self.tab_2)
        self.widget_3.setGeometry(QtCore.QRect(20, 500, 700, 420))
        self.widget_3.setObjectName("widget_3")
        self.flay.removeRow(3)
        self.flay.insertRow(3, self.widget_3)

        self.groupBox_5 = QtWidgets.QGroupBox(self.widget_3)
        self.groupBox_5.setGeometry(QtCore.QRect(20, 600, 431, 61))
        self.groupBox_5.setObjectName("groupBox_5")

        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_5)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.num_subplots = self.spinBox_3.value()
        self.lineEdit_21_list = [None] * self.num_subplots
        self.lineEdit_22_list = [None] * self.num_subplots
        self.lineEdit_23_list = [None] * self.num_subplots
        self.lineEdit_24_list = [None] * self.num_subplots
        self.lineEdit_25_list = [None] * self.num_subplots
        self.lineEdit_25_list = [None] * self.num_subplots
        self.checkBox_26_list = [None] * self.num_subplots
        self.checkBox_27_list = [None] * self.num_subplots
        self.checkBox_28_list = [None] * self.num_subplots
        self.checkBox_29_list = [None] * self.num_subplots

        for i_subplot in range(0, self.num_subplots):
            self.label_21 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_21, i_subplot, 0, 1, 1)
            self.label_21.setObjectName("label_21")
            self.lineEdit_21_list[i_subplot] = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_21_list[i_subplot], i_subplot, 1, 1, 1)
            self.lineEdit_21_list[i_subplot].setObjectName("lineEdit_21")

            self.label_22 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_22, i_subplot, 2, 1, 1)
            self.label_22.setObjectName("label_22")
            self.lineEdit_22_list[i_subplot] = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_22_list[i_subplot], i_subplot, 3, 1, 1)
            self.lineEdit_22_list[i_subplot].setObjectName("lineEdit_22")

            self.label_23 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_23, i_subplot, 4, 1, 1)
            self.label_23.setObjectName("label_23")
            self.lineEdit_23_list[i_subplot] = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_23_list[i_subplot], i_subplot, 5, 1, 1)
            self.lineEdit_23_list[i_subplot].setObjectName("lineEdit_23")

            self.label_24 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_24, i_subplot, 6, 1, 1)
            self.label_24.setObjectName("label_24")
            self.lineEdit_24_list[i_subplot] = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_24_list[i_subplot], i_subplot, 7, 1, 1)
            self.lineEdit_24_list[i_subplot].setObjectName("lineEdit_24")

            self.label_25 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_25, i_subplot, 8, 1, 1)
            self.label_25.setObjectName("label_25")
            self.lineEdit_25_list[i_subplot] = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_25_list[i_subplot], i_subplot, 9, 1, 1)
            self.lineEdit_25_list[i_subplot].setObjectName("lineEdit_25")

            self.checkBox_26_list[i_subplot] = QtWidgets.QCheckBox(self.groupBox_5)
            self.gridLayout.addWidget(self.checkBox_26_list[i_subplot], i_subplot, 10, 1, 1)
            self.checkBox_26_list[i_subplot].setObjectName("checkBox_26")
            self.checkBox_26_list[i_subplot].toggle()

            self.checkBox_27_list[i_subplot] = QtWidgets.QCheckBox(self.groupBox_5)
            self.gridLayout.addWidget(self.checkBox_27_list[i_subplot], i_subplot, 11, 1, 1)
            self.checkBox_27_list[i_subplot].setObjectName("checkBox_27")
            self.checkBox_27_list[i_subplot].toggle()

            self.checkBox_28_list[i_subplot] = QtWidgets.QCheckBox(self.groupBox_5)
            self.gridLayout.addWidget(self.checkBox_28_list[i_subplot], i_subplot, 12, 1, 1)
            self.checkBox_28_list[i_subplot].setObjectName("checkBox_28")
            self.checkBox_28_list[i_subplot].toggle()
            
            self.checkBox_29_list[i_subplot] = QtWidgets.QCheckBox(self.groupBox_5)
            self.gridLayout.addWidget(self.checkBox_29_list[i_subplot], i_subplot, 13, 1, 1)
            self.checkBox_29_list[i_subplot].setObjectName("checkBox_29")
            self.checkBox_29_list[i_subplot].toggle()

            self.label_21.setText("subplot_arg")
            self.label_21.setStatusTip('Subplot argument. e.g. 111 or 212 or 334 etc.')
            self.lineEdit_21_list[i_subplot].setText("111")
            self.lineEdit_21_list[i_subplot].setStatusTip('Subplot argument. e.g. 111 or 212 or 334 etc.')
            self.label_22.setText("xlo")
            self.label_22.setStatusTip('The minimum value of the x axis in the specific subplot')
            self.lineEdit_22_list[i_subplot].setText("None")
            self.lineEdit_22_list[i_subplot].setStatusTip('The minimum value of the x axis in the specific subplot. Default: None')
            self.label_23.setText("xhi")
            self.label_23.setStatusTip('The maximum value of the x axis in the specific subplot')
            self.lineEdit_23_list[i_subplot].setText("None")
            self.lineEdit_23_list[i_subplot].setStatusTip('The maximum value of the x axis in the specific subplot. Default: None')
            self.label_24.setText("ylo")
            self.label_24.setStatusTip('The minimum value of the y axis in the specific subplot')
            self.lineEdit_24_list[i_subplot].setText("None")
            self.lineEdit_24_list[i_subplot].setStatusTip('The minimum value of the y axis in the specific subplot. Default: None')
            self.label_25.setText("yhi")
            self.label_25.setStatusTip('The maximum value of the y axis in the specific subplot')
            self.lineEdit_25_list[i_subplot].setText("None")
            self.lineEdit_25_list[i_subplot].setStatusTip('The maximum value of the y axis in the specific subplot. Default: None')

            self.checkBox_26_list[i_subplot].setText("xtick")
            self.checkBox_26_list[i_subplot].setStatusTip('Whether to plot the x axis tick of a specific subplot or not')
            self.checkBox_27_list[i_subplot].setText("ytick")
            self.checkBox_27_list[i_subplot].setStatusTip('Whether to plot the y axis tick of a specific subplot or not')
            self.checkBox_28_list[i_subplot].setText("xlabel")
            self.checkBox_28_list[i_subplot].setStatusTip('Whether to plot the x axis label of a specific subplot or not')
            self.checkBox_29_list[i_subplot].setText("ylabel")
            self.checkBox_29_list[i_subplot].setStatusTip('Whether to plot the y axis label of a specific subplot or not')
            
        self.groupBox_5.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_5)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget_3)
        layout.addWidget(scroll)

        self.groupBox_5.setTitle("Subplot")

    def build_group2(self,vasp_plotWindow):
        '''element'''
        self.widget_4 = QtWidgets.QWidget(self.tab_2)
        self.widget_4.setGeometry(QtCore.QRect(700, 500, 400, 420))
        self.widget_4.setObjectName("widget_4")
        self.flay.removeRow(4)
        self.flay.insertRow(4, self.widget_4)

        self.groupBox_2 = QtWidgets.QGroupBox(self.widget_4)
        self.groupBox_2.setGeometry(QtCore.QRect(710, 20, 261, 101))
        self.groupBox_2.setObjectName("groupBox_2")

        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_2)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.num_elmts = self.spinBox_4.value()

        self.lineEdit_31_list = [None] * self.num_elmts
        self.checkBox_32_list = [None] * self.num_elmts
        self.checkBox_33_list = [None] * self.num_elmts
        self.checkBox_34_list = [None] * self.num_elmts
        self.checkBox_35_list = [None] * self.num_elmts
        
        for i_elmt in range(0, self.num_elmts):
            self.label_31 = QtWidgets.QLabel(self.groupBox_2)
            self.gridLayout.addWidget(self.label_31, i_elmt, 0, 1, 1)
            self.label_31.setObjectName("label_31")
            
            self.lineEdit_31_list[i_elmt] = QtWidgets.QLineEdit(self.groupBox_2)
            self.gridLayout.addWidget(self.lineEdit_31_list[i_elmt], i_elmt, 1, 1, 1)
            self.lineEdit_31_list[i_elmt].setObjectName("lineEdit_31")

            self.checkBox_32_list[i_elmt] = QtWidgets.QCheckBox(self.groupBox_2)
            self.gridLayout.addWidget(self.checkBox_32_list[i_elmt], i_elmt, 2, 1, 1)
            self.checkBox_32_list[i_elmt].setObjectName("checkBox_32")
            self.checkBox_32_list[i_elmt].toggle()

            self.checkBox_33_list[i_elmt] = QtWidgets.QCheckBox(self.groupBox_2)
            self.gridLayout.addWidget(self.checkBox_33_list[i_elmt], i_elmt, 3, 1, 1)
            self.checkBox_33_list[i_elmt].setObjectName("checkBox_33")
            self.checkBox_33_list[i_elmt].toggle()

            self.checkBox_34_list[i_elmt] = QtWidgets.QCheckBox(self.groupBox_2)
            self.gridLayout.addWidget(self.checkBox_34_list[i_elmt], i_elmt, 4, 1, 1)
            self.checkBox_34_list[i_elmt].setObjectName("checkBox_34")
            self.checkBox_34_list[i_elmt].toggle()
            
            self.checkBox_35_list[i_elmt] = QtWidgets.QCheckBox(self.groupBox_2)
            self.gridLayout.addWidget(self.checkBox_35_list[i_elmt], i_elmt, 5, 1, 1)
            self.checkBox_35_list[i_elmt].setObjectName("checkBox_35")

            self.label_31.setText("element")
            self.lineEdit_31_list[i_elmt].setText("Ni")
            self.lineEdit_31_list[i_elmt].setStatusTip('Specific the element name')
            self.checkBox_32_list[i_elmt].setText("s")
            self.checkBox_32_list[i_elmt].setStatusTip('Plot the s-PDOS of the specific element?')
            self.checkBox_33_list[i_elmt].setText("p")
            self.checkBox_33_list[i_elmt].setStatusTip('Plot the p-PDOS of the specific element?')
            self.checkBox_34_list[i_elmt].setText("d")
            self.checkBox_34_list[i_elmt].setStatusTip('Plot the d-PDOS of the specific element?')
            self.checkBox_35_list[i_elmt].setText("LDOS")
            self.checkBox_35_list[i_elmt].setStatusTip('Plot the LDOS of the specific element?')

            
        self.groupBox_2.setLayout(self.gridLayout)
        
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_2)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget_4)
        layout.addWidget(scroll)

        self.groupBox_2.setTitle("Element")

    def build_group6(self,vasp_plotWindow):
        self.widget_5 = QtWidgets.QWidget(self.tab_2)
        self.widget_5.setGeometry(QtCore.QRect(700, 500, 400, 420))
        self.widget_5.setObjectName("widget_5")
        self.flay.addRow(self.widget_5)
        self.flay.removeRow(5)
        self.flay.insertRow(5, self.widget_5)

        self.groupBox_6 = QtWidgets.QGroupBox(self.widget_5)
        self.groupBox_6.setGeometry(QtCore.QRect(710, 30, 481, 101))
        self.groupBox_6.setObjectName("groupBox_6")
        self.flay.addRow(self.groupBox_6)

        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_6)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.checkBox_40 = QtWidgets.QCheckBox(self.groupBox_6)
        self.gridLayout.addWidget(self.checkBox_40, 0, 0, 1, 1)
        self.checkBox_40.setObjectName("checkBox_40")
        self.checkBox_40.toggle()

        self.label_41 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_41, 1, 0, 1, 1)
        self.label_41.setObjectName("label")
        self.lineEdit_41 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_41, 1, 1, 1, 1)
        self.lineEdit_41.setObjectName("lineEdit_41")
        
        self.label_42 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_42, 1, 2, 1, 1)
        self.label_42.setObjectName("label_42")
        self.lineEdit_42 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_42, 1, 3, 1, 1)
        self.lineEdit_42.setObjectName("lineEdit_42")
        
        self.label_43 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_43, 1, 4, 1, 1)
        self.label_43.setObjectName("label_43")
        self.lineEdit_43 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_43, 1, 5, 1, 1)
        self.lineEdit_43.setObjectName("lineEdit_43")

        self.label_44 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_44, 1, 6, 1, 1)
        self.label_44.setObjectName("label_44")
        self.comboBox_44 = QtWidgets.QComboBox(self.groupBox_6)
        self.gridLayout.addWidget(self.comboBox_44, 1, 7, 1, 1)
        self.comboBox_44.setObjectName("comboBox")
        self.comboBox_44.addItems(self.fig_format_list)
##        self.comboBox_44.addItem("png")
        self.comboBox_44.activated[str].connect(self.fig_format_handleActivated)
        
        self.checkBox_45 = QtWidgets.QCheckBox(self.groupBox_6)
        self.gridLayout.addWidget(self.checkBox_45, 2, 0, 1, 1)
        self.checkBox_45.setObjectName("checkBox_45")

        self.checkBox_46 = QtWidgets.QCheckBox(self.groupBox_6)
        self.gridLayout.addWidget(self.checkBox_46, 2, 1, 1, 1)
        self.checkBox_46.setObjectName("checkBox_46")

        self.checkBox_47 = QtWidgets.QCheckBox(self.groupBox_6)
        self.gridLayout.addWidget(self.checkBox_47, 3, 0, 1, 1)
        self.checkBox_47.setObjectName("checkBox_47")

        self.checkBox_48 = QtWidgets.QCheckBox(self.groupBox_6)
        self.gridLayout.addWidget(self.checkBox_48, 3, 1, 1, 1)
        self.checkBox_48.setObjectName("checkBox_48")

        self.label_49 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_49, 4, 0, 1, 1)
        self.label_49.setObjectName("label_49")
        self.comboBox_49 = QtWidgets.QComboBox(self.groupBox_6)
        self.gridLayout.addWidget(self.comboBox_49, 4, 1, 1, 1)
        self.comboBox_49.setObjectName("comboBox")
        self.comboBox_49.addItem("out")
        self.comboBox_49.addItem("in")
        self.comboBox_49.activated[str].connect(self.xtick_direction_handleActivated)

        self.label_50 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_50, 4, 2, 1, 1)
        self.label_50.setObjectName("label_50")
        self.comboBox_50 = QtWidgets.QComboBox(self.groupBox_6)
        self.gridLayout.addWidget(self.comboBox_50, 4, 3, 1, 1)
        self.comboBox_50.setObjectName("comboBox")
        self.comboBox_50.addItem("out")
        self.comboBox_50.addItem("in")
        self.comboBox_50.activated[str].connect(self.ytick_direction_handleActivated)

        self.label_51 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_51, 5, 0, 1, 1)
        self.label_51.setObjectName("label_51")
        self.comboBox_51 = QtWidgets.QComboBox(self.groupBox_6)
        self.gridLayout.addWidget(self.comboBox_51, 5, 1, 1, 1)
        self.comboBox_51.setObjectName("comboBox")
        self.comboBox_51.addItem("False")
        self.comboBox_51.addItem("True")
        self.comboBox_51.activated[str].connect(self.peak_analyzer_handleActivated)

        self.label_52 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_52, 5, 2, 1, 1)
        self.label_52.setObjectName("label_52")
        self.lineEdit_52 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_52, 5, 3, 1, 1)
        self.lineEdit_52.setObjectName("lineEdit_52")

        self.label_53 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_53, 6, 0, 1, 1)
        self.label_53.setObjectName("label_53")
        self.comboBox_53 = QtWidgets.QComboBox(self.groupBox_6)
        self.gridLayout.addWidget(self.comboBox_53, 6, 1, 1, 1)
        self.comboBox_53.setObjectName("comboBox")
        self.comboBox_53.addItem("False")
        self.comboBox_53.addItem("True")
        self.comboBox_53.activated[str].connect(self.smoothing_handleActivated)

        self.label_54 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_54, 6, 2, 1, 1)
        self.label_54.setObjectName("label_54")
        self.lineEdit_54 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_54, 6, 3, 1, 1)
        self.lineEdit_54.setObjectName("lineEdit_54")

        self.label_55 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_55, 7, 0, 1, 1)
        self.label_55.setObjectName("label")
        self.lineEdit_55 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_55, 7, 1, 1, 1)
        self.lineEdit_55.setObjectName("lineEdit_55")
        
        self.label_56 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_56, 7, 2, 1, 1)
        self.label_56.setObjectName("label_56")
        self.lineEdit_56 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_56, 7, 3, 1, 1)
        self.lineEdit_56.setObjectName("lineEdit_56")        

        self.checkBox_40.setText("fermi_shift_zero")
        self.checkBox_40.setStatusTip("fermi_shift_zero: Whether the Fermi energy is shifted to zero or not")
        self.label_41.setText("fig_width")        
        self.lineEdit_41.setText("15")
        self.label_42.setText("fig_height")        
        self.lineEdit_42.setText("10")
        self.label_43.setText("fig_dpi")        
        self.lineEdit_43.setText("600")
        self.label_44.setText("fig_format")        
        self.checkBox_45.setText("mainplot: xlabel")
        self.checkBox_45.setStatusTip("It defines whether the x axis label of the main plot will be shown or not")
        self.checkBox_46.setText("mainplot: ylabel")        
        self.checkBox_46.setStatusTip("It defines whether the y axis label of the main plot will be shown or not")
        self.checkBox_47.setText("share x axis")        
        self.checkBox_47.setStatusTip("share x axis: it defines whether the subplots share the same x axis or not.")        
        self.checkBox_48.setText("share y axis")
        self.checkBox_48.setStatusTip("share y axis: it defines whether the subplots share the same y axis or not.")        
        self.label_49.setText("xtick_direction")
        self.label_49.setStatusTip('The direction of the x axis tick in the plot (pointing inward or pointing outward)')
        self.label_50.setText("ytick_direction") 
        self.label_50.setStatusTip('The direction of the y axis tick in the plot (pointing inward or pointing outward)')
        
        self.label_51.setText("peak_analyzer")
        self.label_51.setStatusTip('Peak finder. If True, the peaks will be labeled out')
        self.label_52.setText("peak_analyzer_factor")        
        self.lineEdit_52.setText("0.02")
        self.lineEdit_52.setStatusTip('The smaller this value, the more fine peaks can be found')

        self.label_53.setText("smoothing")              
        self.label_53.setStatusTip('Lorentian broadening scheme')
        self.label_54.setText("smoothing_factor")
        self.lineEdit_54.setText("0.05")
        self.lineEdit_54.setStatusTip('Float type. This defines the smoothing factor. In the case of the Lorentzian broadening scheme, the smoothing factor is the broadening width')

        self.label_55.setText("line_width")        
        self.label_55.setStatusTip('Float type. Line width. Recommended value 0.5 -- 3.0 (from thin to fat)')
        self.lineEdit_55.setText("2.0")
        self.lineEdit_55.setStatusTip('Float type. Line width. Recommended value 0.5 -- 3.0 (from thin to fat)')
        self.label_56.setText("font_size")        
        self.label_56.setStatusTip('This value designate the font size for the axis label font size and the legend font size. Recommended value is 18')
        self.lineEdit_56.setText("18")
        self.lineEdit_56.setStatusTip('This value designate the font size for the axis label font size and the legend font size. Recommended value is 18')
        
        self.groupBox_6.setLayout(self.gridLayout)
        
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_6)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget_5)
        layout.addWidget(scroll)

        self.groupBox_6.setTitle("General settings")


    def on_pushButton_clicked(self, i_doscar):
##        print('Button {0} is pressed'.format(i_doscar))
        self.open_dialog_box(i_doscar)
        
    def open_dialog_box(self, i_doscar):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]

        self.lineEdit_list[i_doscar].setText('{}'.format(file_path))

    def fig_format_handleActivated(self, text):
        if text in self.fig_format_list:
            self.fig_format = text
            self.fig_format_a0 = text
    def xtick_direction_handleActivated(self, text):
        self.xtick_direction = text
    def ytick_direction_handleActivated(self, text):
        self.ytick_direction = text
    def peak_analyzer_handleActivated(self, text):
        self.peak_analyzer = text
    def smoothing_handleActivated(self, text):
        self.smoothing = text
