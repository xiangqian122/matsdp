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
from matsdp.vasp import vasp_write

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from functools import partial

defaults_dict = default_params.default_params()
logfile = defaults_dict['logfile']

periodic_table_dict = periodic_table.periodic_tab()
dos_mode = periodic_table_dict['dos_mode']
elmt_color = periodic_table_dict['elmt_color']

class Ui_vasp_writeWindow(object):
    def setupUi(self, vasp_writeWindow):
        vasp_writeWindow.setObjectName("vasp_writeWindow")
        vasp_writeWindow.resize(1030, 600)
        self.centralwidget = QtWidgets.QWidget(vasp_writeWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 990, 511))
        self.tabWidget.setObjectName("tabWidget")

        ###############################################
        # vasp_write
        ###############################################

        # write_poscar_with_force
        self.tab_d1 = QtWidgets.QScrollArea()
        self.tab_d1.setObjectName("tab")
        self.tabWidget.addTab(self.tab_d1, "")
        
        content_widget_d0 = QtWidgets.QWidget()
        self.tab_d1.setWidget(content_widget_d0)
        self.flay_d0 = QtWidgets.QFormLayout(content_widget_d0)
        self.tab_d1.setWidgetResizable(True)

        self.widget_d0 = QtWidgets.QWidget(self.tab_d1)
        self.widget_d0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_d0.setObjectName("widget")
        self.flay_d0.addRow(self.widget_d0)

        self.widget_d1 = QtWidgets.QWidget(self.tab_d1)
        self.widget_d1.setGeometry(QtCore.QRect(20, 90, 290, 450))
        self.widget_d1.setObjectName("widget")
        self.flay_d0.addRow(self.widget_d1)
        
        #parameters
        self.ionic_step = 'last'
        self.outcar_file_path = None
        self.output_poscar_file_name = None

        #Initialize the groups
        self.build_group_d0(vasp_writeWindow)
        self.build_group_d1(vasp_writeWindow)

        vasp_writeWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(vasp_writeWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
##        self.menuHelp = QtWidgets.QMenu(self.menubar)
##        self.menuHelp.setObjectName("menuHelp")
        vasp_writeWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(vasp_writeWindow)
        self.statusbar.setObjectName("statusbar")
        vasp_writeWindow.setStatusBar(self.statusbar)
        self.actionRun = QtWidgets.QAction(vasp_writeWindow)
        self.actionRun.setObjectName("actionRun")
        self.action_1 = QtWidgets.QAction(vasp_writeWindow)
        self.action_1.setObjectName("action_1")
        self.action_2 = QtWidgets.QAction(vasp_writeWindow)
        self.action_2.setObjectName("action_2")
##        self.actionHelp = QtWidgets.QAction(vasp_writeWindow)
##        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_1 = QtWidgets.QAction(vasp_writeWindow)
        self.actionRun_1.setObjectName("actionRun_1")
        self.menuRun.addAction(self.actionRun)
        self.menubar.addAction(self.menuRun.menuAction())
##        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(vasp_writeWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(vasp_writeWindow)


    def retranslateUi(self, vasp_writeWindow):
        _translate = QtCore.QCoreApplication.translate
        vasp_writeWindow.setWindowTitle(_translate("vasp_writeWindow", "vasp_write"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_d1), _translate("vasp_writeWindow", "write_poscar_with_force"))

        self.menuRun.setTitle(_translate("vasp_writeWindow", "Run"))
##        self.menuHelp.setTitle(_translate("vasp_writeWindow", "Help"))
        self.actionRun.setText(_translate("vasp_writeWindow", "write_poscar_with_force"))
        self.actionRun.setStatusTip(_translate("vasp_writeWindow", "write POSCAR with force"))
        self.action_1.setText(_translate("vasp_writeWindow", "write_poscar_with_force"))
##        self.actionHelp.setText(_translate("vasp_writeWindow", "Help"))

    #####################
    # functions
    #####################

    ## write_poscar_with_force
    
    def build_group_d0(self,vasp_writeWindow):
        self.widget_d0 = QtWidgets.QWidget(self.tab_d1)
        self.widget_d0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_d0.setObjectName("widget_d0")
        self.flay_d0.removeRow(1)
        self.flay_d0.insertRow(1, self.widget_d0)

        self.groupBox_d0 = QtWidgets.QGroupBox(self.widget_d0)
        self.groupBox_d0.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_d0.setObjectName("groupBox_d0")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_d0)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
  
        self.label_d1 = QtWidgets.QLabel(self.groupBox_d0)
        self.gridLayout.addWidget(self.label_d1, 0, 0, 1, 1)
        self.label_d1.setObjectName("label_d1")
        self.lineEdit_d1 = QtWidgets.QLineEdit(self.groupBox_d0)
        self.gridLayout.addWidget(self.lineEdit_d1, 0, 1, 1, 1)
        self.lineEdit_d1.setObjectName("lineEdit")
 

        self.label_d2 = QtWidgets.QLabel(self.groupBox_d0)
        self.gridLayout.addWidget(self.label_d2, 1, 0, 1, 1)
        self.label_d2.setObjectName("label_d2")
        self.lineEdit_d2 = QtWidgets.QLineEdit(self.groupBox_d0)
        self.gridLayout.addWidget(self.lineEdit_d2, 1, 1, 1, 1)
        self.lineEdit_d2.setObjectName("lineEdit")
            
        self.label_d1.setText('ionic_step')
        self.label_d1.setStatusTip('ionic_step: The specific ionic step. Values: first, last, or integer number')
        self.lineEdit_d1.setText('last')
        self.lineEdit_d1.setStatusTip('ionic_step: The specific ionic step. Values: first, last, or integer number')
        self.label_d2.setText('output_poscar_file_name')
        self.label_d2.setStatusTip('output_poscar_file_name: Specify the user-defined output POSCAR file name')
        self.lineEdit_d2.setText('None')
        self.lineEdit_d2.setStatusTip('Specify the user-defined output POSCAR file name. If value is None, the output file name will be set automatically.')

        self.groupBox_d0.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_d0)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(240)
        layout = QtWidgets.QVBoxLayout(self.widget_d0)
        layout.addWidget(scroll)
        
        self.groupBox_d0.setTitle("general settings")


    def build_group_d1(self,vasp_writeWindow):
        self.widget_d1 = QtWidgets.QWidget(self.tab_d1)
        self.widget_d1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_d1.setObjectName("widget")
        self.flay_d0.removeRow(2)
        self.flay_d0.insertRow(2, self.widget_d1)

        self.groupBox_d1 = QtWidgets.QGroupBox(self.widget_d1)
        self.groupBox_d1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_d1.setObjectName("groupBox_d1")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_d1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")


        self.label_d7 = QtWidgets.QLabel(self.groupBox_d1)
        self.gridLayout.addWidget(self.label_d7, 0, 0, 1, 1)
        self.label_d7.setObjectName("label_d7")
        self.pushButton_d7 = QtWidgets.QPushButton(self.groupBox_d1)
        self.gridLayout.addWidget(self.pushButton_d7, 0, 1, 1, 1)
        self.pushButton_d7.setObjectName("pushButton_d7")
        self.pushButton_d7.clicked.connect(self.on_pushButton_clicked_d0)        
        self.lineEdit_d7 = QtWidgets.QLineEdit(self.groupBox_d1)
        self.gridLayout.addWidget(self.lineEdit_d7, 0, 2, 1, 1)
        self.lineEdit_d7.setObjectName("lineEdit_d7")
        
        self.label_d7.setText('OUTCAR file')
        self.pushButton_d7.setText("...")
        self.pushButton_d7.setStatusTip('Choose a OUTCAR file')
        self.lineEdit_d7.setText("None")
        self.lineEdit_d7.setStatusTip('Enter the OUTCAR file path (relative path is also accepted)')

        self.groupBox_d1.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_d1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(100)
        layout = QtWidgets.QVBoxLayout(self.widget_d1)
        layout.addWidget(scroll)
        
        self.groupBox_d1.setTitle("file importer")

    def on_pushButton_clicked_d0(self):
        self.open_dialog_box_d0_file()

        
    def open_dialog_box_d0_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_d7.setText('{}'.format(file_path))




