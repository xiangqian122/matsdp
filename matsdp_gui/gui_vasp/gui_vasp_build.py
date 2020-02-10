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
from matsdp.vasp import vasp_build

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from functools import partial

defaults_dict = default_params.default_params()
logfile = defaults_dict['logfile']

periodic_table_dict = periodic_table.periodic_tab()
dos_mode = periodic_table_dict['dos_mode']
elmt_color = periodic_table_dict['elmt_color']

class Ui_vasp_buildWindow(object):
    def setupUi(self, vasp_buildWindow):
        vasp_buildWindow.setObjectName("vasp_buildWindow")
        vasp_buildWindow.resize(1030, 600)
        self.centralwidget = QtWidgets.QWidget(vasp_buildWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 990, 511))
        self.tabWidget.setObjectName("tabWidget")

        ###############################################
        # substitution
        ###############################################

        # substitution
        self.tab_s1 = QtWidgets.QScrollArea()
        self.tab_s1.setObjectName("tab")
        self.tabWidget.addTab(self.tab_s1, "")
        
        content_widget_s0 = QtWidgets.QWidget()
        self.tab_s1.setWidget(content_widget_s0)
        self.flay_s0 = QtWidgets.QFormLayout(content_widget_s0)
        self.tab_s1.setWidgetResizable(True)
        
        self.widget_s0 = QtWidgets.QWidget(self.tab_s1)
        self.widget_s0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_s0.setObjectName("widget")
        self.flay_s0.addRow(self.widget_s0)

        self.widget_s1 = QtWidgets.QWidget(self.tab_s1)
        self.widget_s1.setGeometry(QtCore.QRect(20, 90, 290, 450))
        self.widget_s1.setObjectName("widget")
        self.flay_s0.addRow(self.widget_s1)
        
        #parameters
        self.poscar_file_path_s0 = None
        self.substitution_list_file_path_s0 = None

        #Initialize the groups
##        self.build_group_s0(vasp_buildWindow)
        self.build_group_s1(vasp_buildWindow)

        ###############################################
        # selection_sphere
        ###############################################
              
        # selection_sphere
        self.tab_t2 = QtWidgets.QScrollArea()
        self.tabWidget.addTab(self.tab_t2, "")
        content_widget_t2 = QtWidgets.QWidget()
        self.tab_t2.setWidget(content_widget_t2)
        self.flay_t2 = QtWidgets.QFormLayout(content_widget_t2)
        self.tab_t2.setWidgetResizable(True)

        self.tab_t2.setObjectName("tab_t2")
        self.groupBox_t2 = QtWidgets.QGroupBox(self.tab_t2)
        self.groupBox_t2.setGeometry(QtCore.QRect(20, 10, 800, 61))
        self.groupBox_t2.setObjectName("groupBox_t2")
        self.flay_t2.addRow(self.groupBox_t2)

        self.widget_t1 = QtWidgets.QWidget(self.tab_t2)
        self.widget_t1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_t1.setObjectName("widget_t1")
        self.flay_t2.addRow(self.widget_t1)
        
        self.widget_t2 = QtWidgets.QWidget(self.tab_t2)
        self.widget_t2.setGeometry(QtCore.QRect(300, 70, 500, 420))
        self.widget_t2.setObjectName("widget_t2")
        self.flay_t2.addRow(self.widget_t2)

        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_t2)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget0")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout0")
        
        #parameters
        self.poscar_file_path_t0 = None
        self.origin_atom_name_t0 = 'Ni1'
        self.radius_t0 = 7.0
        self.output_file_name_t0 = 'example'
        self.include_mirror_atoms_t0 = True

        #Initialize the groups
        self.build_group_t1(vasp_buildWindow)
        self.build_group_t2(vasp_buildWindow)

        
        vasp_buildWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(vasp_buildWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
##        self.menuHelp = QtWidgets.QMenu(self.menubar)
##        self.menuHelp.setObjectName("menuHelp")
        vasp_buildWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(vasp_buildWindow)
        self.statusbar.setObjectName("statusbar")
        vasp_buildWindow.setStatusBar(self.statusbar)
        self.actionRun = QtWidgets.QAction(vasp_buildWindow)
        self.actionRun.setObjectName("actionRun")
        self.action_1 = QtWidgets.QAction(vasp_buildWindow)
        self.action_1.setObjectName("action_1")
        self.action_2 = QtWidgets.QAction(vasp_buildWindow)
        self.action_2.setObjectName("action_2")
##        self.actionHelp = QtWidgets.QAction(vasp_buildWindow)
##        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_1 = QtWidgets.QAction(vasp_buildWindow)
        self.actionRun_1.setObjectName("actionRun_1")
        self.menuRun.addAction(self.actionRun)
        self.menuRun.addAction(self.actionRun_1)
        self.menubar.addAction(self.menuRun.menuAction())
##        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(vasp_buildWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(vasp_buildWindow)


    def retranslateUi(self, vasp_buildWindow):
        _translate = QtCore.QCoreApplication.translate
        vasp_buildWindow.setWindowTitle(_translate("vasp_buildWindow", "vasp_build"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_s1), _translate("vasp_buildWindow", "substitution"))

        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_t2), _translate("vasp_buildWindow", "selection_sphere"))
        self.menuRun.setTitle(_translate("vasp_buildWindow", "Run"))
##        self.menuHelp.setTitle(_translate("vasp_buildWindow", "Help"))
        self.actionRun.setText(_translate("vasp_buildWindow", "Run substitution"))
        self.actionRun.setStatusTip(_translate("vasp_buildWindow", "Run substitution"))
        self.action_1.setText(_translate("vasp_buildWindow", "substitution"))
        self.action_2.setText(_translate("vasp_buildWindow", "selection_sphere"))
##        self.actionHelp.setText(_translate("vasp_buildWindow", "Help"))
        self.actionRun_1.setText(_translate("vasp_buildWindow", "Run selection_sphere"))
        self.actionRun_1.setStatusTip(_translate("vasp_buildWindow", "Run selection_sphere"))

    #####################
    # functions
    #####################

    ## substitution
    
    def build_group_s1(self,vasp_buildWindow):
        self.widget_s1 = QtWidgets.QWidget(self.tab_s1)
        self.widget_s1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_s1.setObjectName("widget")
        self.flay_s0.removeRow(2)
        self.flay_s0.insertRow(2, self.widget_s1)

        self.groupBox_s1 = QtWidgets.QGroupBox(self.widget_s1)
        self.groupBox_s1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_s1.setObjectName("groupBox_s1")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_s1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_s7 = QtWidgets.QLabel(self.groupBox_s1)
        self.gridLayout.addWidget(self.label_s7, 0, 0, 1, 1)
        self.label_s7.setObjectName("label_s7")
        self.pushButton_s7 = QtWidgets.QPushButton(self.groupBox_s1)
        self.gridLayout.addWidget(self.pushButton_s7, 0, 1, 1, 1)
        self.pushButton_s7.setObjectName("pushButton_s7")
        self.pushButton_s7.clicked.connect(self.on_pushButton_clicked_s0)        
        self.lineEdit_s7 = QtWidgets.QLineEdit(self.groupBox_s1)
        self.gridLayout.addWidget(self.lineEdit_s7, 0, 2, 1, 1)
        self.lineEdit_s7.setObjectName("lineEdit_s7")

        self.label_s8 = QtWidgets.QLabel(self.groupBox_s1)
        self.gridLayout.addWidget(self.label_s8, 1, 0, 1, 1)
        self.label_s8.setObjectName("label_s8")
        self.pushButton_s8 = QtWidgets.QPushButton(self.groupBox_s1)
        self.gridLayout.addWidget(self.pushButton_s8, 1, 1, 1, 1)
        self.pushButton_s8.setObjectName("pushButton_s8")
        self.pushButton_s8.clicked.connect(self.on_pushButton_clicked_s8)        
        self.lineEdit_s8 = QtWidgets.QLineEdit(self.groupBox_s1)
        self.gridLayout.addWidget(self.lineEdit_s8, 1, 2, 1, 1)
        self.lineEdit_s8.setObjectName("lineEdit_s8")

        
        self.label_s7.setText('POSCAR file')
        self.pushButton_s7.setText("...")
        self.pushButton_s7.setStatusTip("Choose a POSCAR file")
        self.lineEdit_s7.setText("None")
        self.lineEdit_s7.setStatusTip("Enter the POSCAR file path (relative path is also accepted)")
        self.label_s8.setText('.subst file')
        self.pushButton_s8.setText("...")
        self.pushButton_s8.setStatusTip("Choose a .subst file")
        self.lineEdit_s8.setText("None")
        self.lineEdit_s8.setStatusTip("Enter the .subst file path (relative path is also accepted)")

        self.groupBox_s1.setLayout(self.gridLayout)


        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_s1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(120)
        layout = QtWidgets.QVBoxLayout(self.widget_s1)
        layout.addWidget(scroll)

        
        self.groupBox_s1.setTitle("file importer")

    def on_pushButton_clicked_s0(self):
        self.open_dialog_box_s0_file()

        
    def open_dialog_box_s0_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_s7.setText('{}'.format(file_path))

    def on_pushButton_clicked_s8(self):
        self.open_dialog_box_s8_file()

        
    def open_dialog_box_s8_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_s8.setText('{}'.format(file_path))



    ## selection_sphere

    def build_group_t1(self,vasp_buildWindow):
        self.widget_t1 = QtWidgets.QWidget(self.tab_t2)
        self.widget_t1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_t1.setObjectName("widget")
        self.flay_t2.removeRow(1)
        self.flay_t2.insertRow(1, self.widget_t1)

        self.groupBox_t1 = QtWidgets.QGroupBox(self.widget_t1)
        self.groupBox_t1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_t1.setObjectName("groupBox_t1")

        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_t1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_t1 = QtWidgets.QLabel(self.groupBox_t1)
        self.gridLayout.addWidget(self.label_t1, 0, 0, 1, 1)
        self.label_t1.setObjectName("label_t1")
        self.lineEdit_t1 = QtWidgets.QLineEdit(self.groupBox_t1)
        self.gridLayout.addWidget(self.lineEdit_t1, 0, 1, 1, 1)
        self.lineEdit_t1.setObjectName("lineEdit")

        self.label_t2 = QtWidgets.QLabel(self.groupBox_t1)
        self.gridLayout.addWidget(self.label_t2, 1, 0, 1, 1)
        self.label_t2.setObjectName("label_t2")
        self.lineEdit_t2 = QtWidgets.QLineEdit(self.groupBox_t1)
        self.gridLayout.addWidget(self.lineEdit_t2, 1, 1, 1, 1)
        self.lineEdit_t2.setObjectName("lineEdit")
        
        self.label_t3 = QtWidgets.QLabel(self.groupBox_t1)
        self.gridLayout.addWidget(self.label_t3, 2, 0, 1, 1)
        self.label_t3.setObjectName("label_t3")
        self.lineEdit_t3 = QtWidgets.QLineEdit(self.groupBox_t1)
        self.gridLayout.addWidget(self.lineEdit_t3, 2, 1, 1, 1)
        self.lineEdit_t3.setObjectName("lineEdit")
        
        self.checkBox_t4 = QtWidgets.QCheckBox(self.groupBox_t1)
        self.gridLayout.addWidget(self.checkBox_t4, 3, 0, 1, 1)
        self.checkBox_t4.setObjectName("checkBox_t4")
        self.checkBox_t4.toggle()        
##        self.checkBox_t4.toggled.connect(self.set_value_t4)


        
        self.label_t1.setText('origin atom name')
        self.label_t1.setStatusTip('The name of the atom which is at the center of the sphere')
        self.lineEdit_t1.setText('Ni1')
        self.lineEdit_t1.setStatusTip('Enter the name of the atom which is at the center of the sphere. e.g. Ni1, Al3, Re2')
        self.label_t2.setText('radius')
        self.label_t2.setStatusTip('The radius of the sphere')
        self.lineEdit_t2.setText('7.0')
        self.lineEdit_t2.setStatusTip('Enter the radius of the sphere. Float type.')
        self.label_t3.setText('output file name')
        self.label_t3.setStatusTip('It species the name of the output POSCAR file(without suffix)')
        self.lineEdit_t3.setText('example')
        self.lineEdit_t3.setStatusTip('String value. It species the name of the output POSCAR file(without suffix)')
        self.checkBox_t4.setText("include mirror atoms")
        self.checkBox_t4.setStatusTip('When selecting the atoms inside the sphere, whether to include the mirror atoms or not')
        
        self.groupBox_t1.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_t1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(240)
        layout = QtWidgets.QVBoxLayout(self.widget_t1)
        layout.addWidget(scroll)
        
        self.groupBox_t1.setTitle("general settings")

    def build_group_t2(self,vasp_buildWindow):
        self.widget_t2 = QtWidgets.QWidget(self.tab_t2)
        self.widget_t2.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_t2.setObjectName("widget")
        self.flay_t2.removeRow(2)
        self.flay_t2.insertRow(2, self.widget_t2)

        self.groupBox_t2 = QtWidgets.QGroupBox(self.widget_t2)
        self.groupBox_t2.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_t2.setObjectName("groupBox_t2")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_t2)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_t7 = QtWidgets.QLabel(self.groupBox_t2)
        self.gridLayout.addWidget(self.label_t7, 0, 0, 1, 1)
        self.label_t7.setObjectName("label_t7")
        self.pushButton_t7 = QtWidgets.QPushButton(self.groupBox_t2)
        self.gridLayout.addWidget(self.pushButton_t7, 0, 1, 1, 1)
        self.pushButton_t7.setObjectName("pushButton_t7")
        self.pushButton_t7.clicked.connect(self.on_pushButton_clicked_t2)        
        self.lineEdit_t7 = QtWidgets.QLineEdit(self.groupBox_t2)
        self.gridLayout.addWidget(self.lineEdit_t7, 0, 2, 1, 1)
        self.lineEdit_t7.setObjectName("lineEdit_t7")
        
        self.label_t7.setText('POSCAR file')
        self.pushButton_t7.setText("...")
        self.pushButton_t7.setStatusTip("Choose a POSCAR file")
        self.lineEdit_t7.setText("None")
        self.lineEdit_t7.setStatusTip("Enter the POSCAR file path (relative path is also accepted)")

        self.groupBox_t2.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_t2)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(100)
        layout = QtWidgets.QVBoxLayout(self.widget_t2)
        layout.addWidget(scroll)
        
        self.groupBox_t2.setTitle("file importer")

    def on_pushButton_clicked_t2(self):
        self.open_dialog_box_t2_file()
         
    def open_dialog_box_t2_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_t7.setText('{}'.format(file_path))
        



