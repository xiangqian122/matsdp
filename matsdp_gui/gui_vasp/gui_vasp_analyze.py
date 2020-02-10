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
##from matsdp.vasp import vasp_read
from matsdp.vasp import vasp_analyze

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from functools import partial

defaults_dict = default_params.default_params()
logfile = defaults_dict['logfile']

periodic_table_dict = periodic_table.periodic_tab()
dos_mode = periodic_table_dict['dos_mode']
elmt_color = periodic_table_dict['elmt_color']

class Ui_vasp_analyzeWindow(object):
    def setupUi(self, vasp_analyzeWindow):
        vasp_analyzeWindow.setObjectName("vasp_analyzeWindow")
        vasp_analyzeWindow.resize(1030, 600)
        self.centralwidget = QtWidgets.QWidget(vasp_analyzeWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 990, 511))
        self.tabWidget.setObjectName("tabWidget")

        ###############################################
        #nn_map
        ###############################################

        # nn_map
        self.tab_b1 = QtWidgets.QScrollArea()
        self.tab_b1.setObjectName("tab")
        self.tabWidget.addTab(self.tab_b1, "")
        
        content_widget_b0 = QtWidgets.QWidget()
        self.tab_b1.setWidget(content_widget_b0)
        self.flay_b0 = QtWidgets.QFormLayout(content_widget_b0)
        self.tab_b1.setWidgetResizable(True)
        
        self.widget_b0 = QtWidgets.QWidget(self.tab_b1)
        self.widget_b0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_b0.setObjectName("widget")
        self.flay_b0.addRow(self.widget_b0)

        self.widget_b1 = QtWidgets.QWidget(self.tab_b1)
        self.widget_b1.setGeometry(QtCore.QRect(20, 90, 290, 450))
        self.widget_b1.setObjectName("widget")
        self.flay_b0.addRow(self.widget_b1)
        
        #parameters
        self.nn_mode = 'NN map'
        self.a0 = 3.545
        self.n_shell = 1
        self.poscar_file_path_b0 = None

        #Initialize the groups
        self.build_group_b0(vasp_analyzeWindow)
        self.build_group_b1(vasp_analyzeWindow)

        ###############################################
        # e_struct
        ###############################################

              
        # e_struct
        self.tab_c2 = QtWidgets.QScrollArea()
        self.tabWidget.addTab(self.tab_c2, "")
        content_widget_c2 = QtWidgets.QWidget()
        self.tab_c2.setWidget(content_widget_c2)
        self.flay_c2 = QtWidgets.QFormLayout(content_widget_c2)
        self.tab_c2.setWidgetResizable(True)

        self.tab_c2.setObjectName("tab_c2")
        self.groupBox_c2 = QtWidgets.QGroupBox(self.tab_c2)
        self.groupBox_c2.setGeometry(QtCore.QRect(20, 10, 800, 61))
        self.groupBox_c2.setObjectName("groupBox_c2")
        self.flay_c2.addRow(self.groupBox_c2)

        self.widget_c1 = QtWidgets.QWidget(self.tab_c2)
        self.widget_c1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_c1.setObjectName("widget_c1")
        self.flay_c2.addRow(self.widget_c1)
        
        self.widget_c2 = QtWidgets.QWidget(self.tab_c2)
        self.widget_c2.setGeometry(QtCore.QRect(300, 70, 500, 420))
        self.widget_c2.setObjectName("widget_c2")
        self.flay_c2.addRow(self.widget_c2)

        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_c2)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget0")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout0")
        
        #parameters
        self.sysname = 'System'
        self.doscar_file_path = None

        #Initialize the groups
        self.build_group_c1(vasp_analyzeWindow)
        self.build_group_c2(vasp_analyzeWindow)

        
        vasp_analyzeWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(vasp_analyzeWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
##        self.menuHelp = QtWidgets.QMenu(self.menubar)
##        self.menuHelp.setObjectName("menuHelp")
        vasp_analyzeWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(vasp_analyzeWindow)
        self.statusbar.setObjectName("statusbar")
        vasp_analyzeWindow.setStatusBar(self.statusbar)
        self.actionRun = QtWidgets.QAction(vasp_analyzeWindow)
        self.actionRun.setObjectName("actionRun")
        self.action_1 = QtWidgets.QAction(vasp_analyzeWindow)
        self.action_1.setObjectName("action_1")
        self.action_2 = QtWidgets.QAction(vasp_analyzeWindow)
        self.action_2.setObjectName("action_2")
##        self.actionHelp = QtWidgets.QAction(vasp_analyzeWindow)
##        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_1 = QtWidgets.QAction(vasp_analyzeWindow)
        self.actionRun_1.setObjectName("actionRun_1")
        self.menuRun.addAction(self.actionRun)
        self.menuRun.addAction(self.actionRun_1)
        self.menubar.addAction(self.menuRun.menuAction())
##        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(vasp_analyzeWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(vasp_analyzeWindow)


    def retranslateUi(self, vasp_analyzeWindow):
        _translate = QtCore.QCoreApplication.translate
        vasp_analyzeWindow.setWindowTitle(_translate("vasp_analyzeWindow", "vasp_analyze"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_b1), _translate("vasp_analyzeWindow", "NN map analysis"))

        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_c2), _translate("vasp_analyzeWindow", "E_struct"))
        self.menuRun.setTitle(_translate("vasp_analyzeWindow", "Run"))
##        self.menuHelp.setTitle(_translate("vasp_analyzeWindow", "Help"))
        self.actionRun.setText(_translate("vasp_analyzeWindow", "Run NN map"))
        self.actionRun.setStatusTip(_translate("vasp_analyzeWindow", "Run NN map"))
        self.action_1.setText(_translate("vasp_analyzeWindow", "nn_map"))
        self.action_2.setText(_translate("vasp_analyzeWindow", "e_struct"))
##        self.actionHelp.setText(_translate("vasp_analyzeWindow", "Help"))
        self.actionRun_1.setText(_translate("vasp_analyzeWindow", "Run e_struct"))
        self.actionRun_1.setStatusTip(_translate("vasp_analyzeWindow", "Run e_struct"))

    #####################
    # functions
    #####################

    ## nn map
    
    def build_group_b0(self,vasp_analyzeWindow):
        self.widget_b0 = QtWidgets.QWidget(self.tab_b1)
        self.widget_b0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_b0.setObjectName("widget_b0")
        self.flay_b0.removeRow(1)
        self.flay_b0.insertRow(1, self.widget_b0)

        self.groupBox_b0 = QtWidgets.QGroupBox(self.widget_b0)
        self.groupBox_b0.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_b0.setObjectName("groupBox_b0")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_b0)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
  
        self.label_b0 = QtWidgets.QLabel(self.groupBox_b0)
        self.gridLayout.addWidget(self.label_b0, 0, 0, 1, 1)
        self.label_b0.setObjectName("label_b0")
        self.comboBox_b0 = QtWidgets.QComboBox(self.groupBox_b0)
        self.gridLayout.addWidget(self.comboBox_b0, 0, 1, 1, 1)
        self.comboBox_b0.setObjectName("comboBox")
        self.comboBox_b0.addItem("NN map")
        self.comboBox_b0.addItem("simple CNA")
        self.comboBox_b0.activated[str].connect(self.nn_mode_handleActivated)
 

        self.label_b1 = QtWidgets.QLabel(self.groupBox_b0)
        self.gridLayout.addWidget(self.label_b1, 1, 0, 1, 1)
        self.label_b1.setObjectName("label_b1")
        self.lineEdit_b1 = QtWidgets.QLineEdit(self.groupBox_b0)
        self.gridLayout.addWidget(self.lineEdit_b1, 1, 1, 1, 1)
        self.lineEdit_b1.setObjectName("lineEdit")
 

        self.label_b2 = QtWidgets.QLabel(self.groupBox_b0)
        self.gridLayout.addWidget(self.label_b2, 2, 0, 1, 1)
        self.label_b2.setObjectName("label_a2")
        self.lineEdit_b2 = QtWidgets.QLineEdit(self.groupBox_b0)
        self.gridLayout.addWidget(self.lineEdit_b2, 2, 1, 1, 1)
        self.lineEdit_b2.setObjectName("lineEdit")
  

        self.label_b3 = QtWidgets.QLabel(self.groupBox_b0)
        self.gridLayout.addWidget(self.label_b3, 3, 0, 1, 1)
        self.label_b1.setObjectName("label_b3")
        self.lineEdit_b3 = QtWidgets.QLineEdit(self.groupBox_b0)
        self.gridLayout.addWidget(self.lineEdit_b3, 3, 1, 1, 1)
        self.lineEdit_b3.setObjectName("lineEdit")
        self.label_b3.setEnabled(False)
        self.lineEdit_b3.setEnabled(False)

            
        self.label_b0.setText('mode')
        self.label_b0.setStatusTip('NN map: nearest neighbor map. simple CNA: simple common neighbor analysis')
        self.label_b1.setText('a0')
        self.label_b1.setStatusTip('The lattice constant of the lattice. Currently, it is only available for the FCC lattice')
        self.lineEdit_b1.setText('3.545')
        self.lineEdit_b1.setStatusTip('Enter the lattice constant of the lattice (Angstrom). Currently, it is only available for the FCC lattice')
        self.label_b2.setText('n_shell')
        self.label_b2.setStatusTip('n_shell: The atoms inside nth nearest neighbor shell are included for nearest neighbor analysis')
        self.lineEdit_b2.setText('1')
        self.lineEdit_b2.setStatusTip('n_shell: The atoms inside nth nearest neighbor shell are included for nearest neighbor analysis')
        self.label_b3.setText('common_neighbor_elmt_list')
        self.label_b3.setStatusTip('common_neighbor_elmt_list: Find out the common neighbor atoms among the selected type of elements. e.g. Ni, Al, Re')
        self.lineEdit_b3.setText('Re, Ni')
        self.lineEdit_b3.setStatusTip('common_neighbor_elmt_list: Find out the common neighbor atoms among the selected type of elements. e.g. Ni, Al, Re')

        self.groupBox_b0.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_b0)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(240)
        layout = QtWidgets.QVBoxLayout(self.widget_b0)
        layout.addWidget(scroll)
        
        self.groupBox_b0.setTitle("general settings")


    def build_group_b1(self,vasp_analyzeWindow):
        self.widget_b1 = QtWidgets.QWidget(self.tab_b1)
        self.widget_b1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_b1.setObjectName("widget")
        self.flay_b0.removeRow(2)
        self.flay_b0.insertRow(2, self.widget_b1)

        self.groupBox_b1 = QtWidgets.QGroupBox(self.widget_b1)
        self.groupBox_b1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_b1.setObjectName("groupBox_b1")

        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_b1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")


        self.label_b7 = QtWidgets.QLabel(self.groupBox_b1)
        self.gridLayout.addWidget(self.label_b7, 0, 0, 1, 1)
        self.label_b7.setObjectName("label_b7")
        self.pushButton_b7 = QtWidgets.QPushButton(self.groupBox_b1)
        self.gridLayout.addWidget(self.pushButton_b7, 0, 1, 1, 1)
        self.pushButton_b7.setObjectName("pushButton_b7")
        self.pushButton_b7.clicked.connect(self.on_pushButton_clicked_b0)        
        self.lineEdit_b7 = QtWidgets.QLineEdit(self.groupBox_b1)
        self.gridLayout.addWidget(self.lineEdit_b7, 0, 2, 1, 1)
        self.lineEdit_b7.setObjectName("lineEdit_b7")

        
        self.label_b7.setText('POSCAR file')
        self.pushButton_b7.setText("...")
        self.pushButton_b7.setStatusTip('Choose a POSCAR file')
        self.lineEdit_b7.setText("None")
        self.lineEdit_b7.setStatusTip("Enter the POSCAR file path (relative path is also accepted)")


        self.groupBox_b1.setLayout(self.gridLayout)


        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_b1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(100)
        layout = QtWidgets.QVBoxLayout(self.widget_b1)
        layout.addWidget(scroll)

        
        self.groupBox_b1.setTitle("file importer")


    def nn_mode_handleActivated(self, text):
        self.nn_mode = text
        if self.nn_mode == 'NN map':
            self.label_b2.setEnabled(True)
            self.lineEdit_b2.setEnabled(True)
            self.label_b3.setEnabled(False)
            self.lineEdit_b3.setEnabled(False)
        elif self.nn_mode == 'simple CNA':
            self.label_b2.setEnabled(False)
            self.lineEdit_b2.setEnabled(False)
            self.label_b3.setEnabled(True)
            self.lineEdit_b3.setEnabled(True)

    def on_pushButton_clicked_b0(self):
        self.open_dialog_box_b0_file()

        
    def open_dialog_box_b0_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_b7.setText('{}'.format(file_path))



    ## e_struct

    def build_group_c1(self,vasp_analyzeWindow):
        self.widget_c1 = QtWidgets.QWidget(self.tab_c2)
        self.widget_c1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_c1.setObjectName("widget")
        self.flay_c2.removeRow(1)
        self.flay_c2.insertRow(1, self.widget_c1)

        self.groupBox_c1 = QtWidgets.QGroupBox(self.widget_c1)
        self.groupBox_c1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_c1.setObjectName("groupBox_c1")

        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_c1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_c1 = QtWidgets.QLabel(self.groupBox_c1)
        self.gridLayout.addWidget(self.label_c1, 1, 0, 1, 1)
        self.label_c1.setObjectName("label_c1")
        self.lineEdit_c1 = QtWidgets.QLineEdit(self.groupBox_c1)
        self.gridLayout.addWidget(self.lineEdit_c1, 1, 1, 1, 1)
        self.lineEdit_c1.setObjectName("lineEdit")

        self.label_c1.setText('sysname')
        self.label_c1.setStatusTip('sysname: user defined system name')
        self.lineEdit_c1.setText('System')
        self.lineEdit_c1.setStatusTip('sysname: user defined system name')

        self.groupBox_c1.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_c1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(240)
        layout = QtWidgets.QVBoxLayout(self.widget_c1)
        layout.addWidget(scroll)
        
        self.groupBox_c1.setTitle("general settings")

    def build_group_c2(self,vasp_analyzeWindow):
        self.widget_c2 = QtWidgets.QWidget(self.tab_c2)
        self.widget_c2.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_c2.setObjectName("widget")
        self.flay_c2.removeRow(2)
        self.flay_c2.insertRow(2, self.widget_c2)

        self.groupBox_c2 = QtWidgets.QGroupBox(self.widget_c2)
        self.groupBox_c2.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_c2.setObjectName("groupBox_c2")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_c2)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_c7 = QtWidgets.QLabel(self.groupBox_c2)
        self.gridLayout.addWidget(self.label_c7, 0, 0, 1, 1)
        self.label_c7.setObjectName("label_c7")
        self.pushButton_c7 = QtWidgets.QPushButton(self.groupBox_c2)
        self.gridLayout.addWidget(self.pushButton_c7, 0, 1, 1, 1)
        self.pushButton_c7.setObjectName("pushButton_c7")
        self.pushButton_c7.clicked.connect(self.on_pushButton_clicked_c2)        
        self.lineEdit_c7 = QtWidgets.QLineEdit(self.groupBox_c2)
        self.gridLayout.addWidget(self.lineEdit_c7, 0, 2, 1, 1)
        self.lineEdit_c7.setObjectName("lineEdit_c7")
        
        self.label_c7.setText('DOSCAR file')
        self.pushButton_c7.setText("...")
        self.pushButton_c7.setStatusTip('Choose a DOSCAR file')
        self.lineEdit_c7.setText("None")
        self.lineEdit_c7.setStatusTip('Enter a DOSCAR file path (relative path is also accepted)')

        self.groupBox_c2.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_c2)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(100)
        layout = QtWidgets.QVBoxLayout(self.widget_c2)
        layout.addWidget(scroll)
        
        self.groupBox_c2.setTitle("file importer")

    def on_pushButton_clicked_c2(self):
        self.open_dialog_box_c2_file()
         
    def open_dialog_box_c2_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_c7.setText('{}'.format(file_path))
        



