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
from matsdp.vasp import vasp_write

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from functools import partial

defaults_dict = default_params.default_params()
logfile = defaults_dict['logfile']

periodic_table_dict = periodic_table.periodic_tab()
dos_mode = periodic_table_dict['dos_mode']
elmt_color = periodic_table_dict['elmt_color']

class Ui_dvm_analyzeWindow(object):
    def setupUi(self, dvm_analyzeWindow):
        dvm_analyzeWindow.setObjectName("dvm_analyzeWindow")
        dvm_analyzeWindow.resize(1030, 600)
        self.centralwidget = QtWidgets.QWidget(dvm_analyzeWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 990, 511))
        self.tabWidget.setObjectName("tabWidget")

        ###############################################
        # vasp_write
        ###############################################

        # write_poscar_with_force
        self.tab_f1 = QtWidgets.QScrollArea()
        self.tab_f1.setObjectName("tab")
        self.tabWidget.addTab(self.tab_f1, "")
        
        content_widget_f0 = QtWidgets.QWidget()
        self.tab_f1.setWidget(content_widget_f0)
        self.flay_f0 = QtWidgets.QFormLayout(content_widget_f0)
        self.tab_f1.setWidgetResizable(True)

        self.widget_f0 = QtWidgets.QWidget(self.tab_f1)
        self.widget_f0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_f0.setObjectName("widget")
        self.flay_f0.addRow(self.widget_f0)

        self.widget_f1 = QtWidgets.QWidget(self.tab_f1)
        self.widget_f1.setGeometry(QtCore.QRect(20, 90, 290, 450))
        self.widget_f1.setObjectName("widget")
        self.flay_f0.addRow(self.widget_f1)
        
        #parameters
        self.dvm_get_ie_mode = 'IE for 1NN atom pairs'
        self.dvm_otput_file_path = None
        self.latt_const_f0 = 3.545
##        self.ionic_step = 'last'
##        self.outcar_file_path = None
##        self.output_poscar_file_name = None

        #Initialize the groups
        self.build_group_f0(dvm_analyzeWindow)
        self.build_group_f1(dvm_analyzeWindow)

        dvm_analyzeWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(dvm_analyzeWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
##        self.menuHelp = QtWidgets.QMenu(self.menubar)
##        self.menuHelp.setObjectName("menuHelp")
        dvm_analyzeWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(dvm_analyzeWindow)
        self.statusbar.setObjectName("statusbar")
        dvm_analyzeWindow.setStatusBar(self.statusbar)
        self.actionRun = QtWidgets.QAction(dvm_analyzeWindow)
        self.actionRun.setObjectName("actionRun")
        self.action_1 = QtWidgets.QAction(dvm_analyzeWindow)
        self.action_1.setObjectName("action_1")
        self.action_2 = QtWidgets.QAction(dvm_analyzeWindow)
        self.action_2.setObjectName("action_2")
##        self.actionHelp = QtWidgets.QAction(dvm_analyzeWindow)
##        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_1 = QtWidgets.QAction(dvm_analyzeWindow)
        self.actionRun_1.setObjectName("actionRun_1")
        self.menuRun.addAction(self.actionRun)
        self.menubar.addAction(self.menuRun.menuAction())
##        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(dvm_analyzeWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(dvm_analyzeWindow)


    def retranslateUi(self, dvm_analyzeWindow):
        _translate = QtCore.QCoreApplication.translate
        dvm_analyzeWindow.setWindowTitle(_translate("dvm_analyzeWindow", "dvm_analyze"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_f1), _translate("dvm_analyzeWindow", "extract IE"))

        self.menuRun.setTitle(_translate("dvm_analyzeWindow", "Run"))
##        self.menuHelp.setTitle(_translate("dvm_analyzeWindow", "Help"))
        self.actionRun.setText(_translate("dvm_analyzeWindow", "extract IE"))
        self.actionRun.setStatusTip(_translate("dvm_analyzeWindow", "extract IE"))
        self.action_1.setText(_translate("dvm_analyzeWindow", "extract IE"))
##        self.actionHelp.setText(_translate("dvm_analyzeWindow", "Help"))

    #####################
    # functions
    #####################

    ## write_poscar_with_force
    
    def build_group_f0(self,dvm_analyzeWindow):
        self.widget_f0 = QtWidgets.QWidget(self.tab_f1)
        self.widget_f0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_f0.setObjectName("widget_f0")
        self.flay_f0.removeRow(1)
        self.flay_f0.insertRow(1, self.widget_f0)

        self.groupBox_f0 = QtWidgets.QGroupBox(self.widget_f0)
        self.groupBox_f0.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_f0.setObjectName("groupBox_f0")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_f0)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_f0 = QtWidgets.QLabel(self.groupBox_f0)
        self.gridLayout.addWidget(self.label_f0, 0, 0, 1, 1)
        self.label_f0.setObjectName("label_f0")
        self.comboBox_f0 = QtWidgets.QComboBox(self.groupBox_f0)
        self.gridLayout.addWidget(self.comboBox_f0, 0, 1, 1, 1)
        self.comboBox_f0.setObjectName("comboBox")
        self.comboBox_f0.addItem("IE for 1NN atom pairs")
        self.comboBox_f0.addItem("IE for all the atom pairs")
        self.comboBox_f0.activated[str].connect(self.dvm_get_ie_mode_handleActivated)

  
        self.label_f1 = QtWidgets.QLabel(self.groupBox_f0)
        self.gridLayout.addWidget(self.label_f1, 1, 0, 1, 1)
        self.label_f1.setObjectName("label_f1")
        self.lineEdit_f1 = QtWidgets.QLineEdit(self.groupBox_f0)
        self.gridLayout.addWidget(self.lineEdit_f1, 1, 1, 1, 1)
        self.lineEdit_f1.setObjectName("lineEdit")

        self.label_f0.setText('mode')            
        self.label_f1.setText('a0')
        self.label_f1.setStatusTip('a0: Lattice constant')
        self.lineEdit_f1.setText('3.545')
        self.lineEdit_f1.setStatusTip('Set the lattice constant of the lattice (Currently, it it only available for the FCC lattice)')

        self.groupBox_f0.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_f0)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(240)
        layout = QtWidgets.QVBoxLayout(self.widget_f0)
        layout.addWidget(scroll)
        
        self.groupBox_f0.setTitle("general settings")


    def build_group_f1(self,dvm_analyzeWindow):
        self.widget_f1 = QtWidgets.QWidget(self.tab_f1)
        self.widget_f1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_f1.setObjectName("widget")
        self.flay_f0.removeRow(2)
        self.flay_f0.insertRow(2, self.widget_f1)

        self.groupBox_f1 = QtWidgets.QGroupBox(self.widget_f1)
        self.groupBox_f1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_f1.setObjectName("groupBox_f1")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_f1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")


        self.label_f7 = QtWidgets.QLabel(self.groupBox_f1)
        self.gridLayout.addWidget(self.label_f7, 0, 0, 1, 1)
        self.label_f7.setObjectName("label_f7")
        self.pushButton_f7 = QtWidgets.QPushButton(self.groupBox_f1)
        self.gridLayout.addWidget(self.pushButton_f7, 0, 1, 1, 1)
        self.pushButton_f7.setObjectName("pushButton_f7")
        self.pushButton_f7.clicked.connect(self.on_pushButton_clicked_f0)        
        self.lineEdit_f7 = QtWidgets.QLineEdit(self.groupBox_f1)
        self.gridLayout.addWidget(self.lineEdit_f7, 0, 2, 1, 1)
        self.lineEdit_f7.setObjectName("lineEdit_f7")
        
        self.label_f7.setText('*.otput file')
        self.pushButton_f7.setText("...")
        self.pushButton_f7.setStatusTip('Choose the .otput file')
        self.lineEdit_f7.setText("None")
        self.lineEdit_f7.setStatusTip('Enter the .otput file path (relative path is also accepted)')

        self.groupBox_f1.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_f1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(100)
        layout = QtWidgets.QVBoxLayout(self.widget_f1)
        layout.addWidget(scroll)
        
        self.groupBox_f1.setTitle("file importer")

    def dvm_get_ie_mode_handleActivated(self, text):
        self.dvm_get_ie_mode = text
        if self.dvm_get_ie_mode == 'IE for 1NN atom pairs':
            self.label_f1.setEnabled(True)
            self.lineEdit_f1.setEnabled(True)
        elif self.dvm_get_ie_mode == 'IE for all the atom pairs':
            self.label_f1.setEnabled(False)
            self.lineEdit_f1.setEnabled(False)

    def on_pushButton_clicked_f0(self):
        self.open_dialog_box_f0_file()

        
    def open_dialog_box_f0_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_f7.setText('{}'.format(file_path))




