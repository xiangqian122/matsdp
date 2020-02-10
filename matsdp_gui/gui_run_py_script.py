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

class Ui_run_py_scriptWindow(object):
    def setupUi(self, run_py_scriptWindow):
        run_py_scriptWindow.setObjectName("run_py_scriptWindow")
        run_py_scriptWindow.resize(1030, 300)
        self.centralwidget = QtWidgets.QWidget(run_py_scriptWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 990, 200))
        self.tabWidget.setObjectName("tabWidget")

        ###############################################
        # run_py_script
        ###############################################

        # run_py_script
        self.tab_g1 = QtWidgets.QScrollArea()
        self.tab_g1.setObjectName("tab")
        self.tabWidget.addTab(self.tab_g1, "")
        
        content_widget_g0 = QtWidgets.QWidget()
        self.tab_g1.setWidget(content_widget_g0)
        self.flay_g0 = QtWidgets.QFormLayout(content_widget_g0)
        self.tab_g1.setWidgetResizable(True)

        self.widget_g0 = QtWidgets.QWidget(self.tab_g1)
        self.widget_g0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_g0.setObjectName("widget")
        self.flay_g0.addRow(self.widget_g0)

        self.widget_g1 = QtWidgets.QWidget(self.tab_g1)
        self.widget_g1.setGeometry(QtCore.QRect(20, 90, 290, 450))
        self.widget_g1.setObjectName("widget")
        self.flay_g0.addRow(self.widget_g1)
        
        #parameters
        self.script_file = None

        #Initialize the groups
##        self.build_group_g0(run_py_scriptWindow)
        self.build_group_g1(run_py_scriptWindow)

        run_py_scriptWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(run_py_scriptWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
##        self.menuHelp = QtWidgets.QMenu(self.menubar)
##        self.menuHelp.setObjectName("menuHelp")
        run_py_scriptWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(run_py_scriptWindow)
        self.statusbar.setObjectName("statusbar")
        run_py_scriptWindow.setStatusBar(self.statusbar)
        self.actionRun = QtWidgets.QAction(run_py_scriptWindow)
        self.actionRun.setObjectName("actionRun")
        self.action_1 = QtWidgets.QAction(run_py_scriptWindow)
        self.action_1.setObjectName("action_1")
        self.action_2 = QtWidgets.QAction(run_py_scriptWindow)
        self.action_2.setObjectName("action_2")
##        self.actionHelp = QtWidgets.QAction(run_py_scriptWindow)
##        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_1 = QtWidgets.QAction(run_py_scriptWindow)
        self.actionRun_1.setObjectName("actionRun_1")
        self.menuRun.addAction(self.actionRun)
        self.menubar.addAction(self.menuRun.menuAction())
##        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(run_py_scriptWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(run_py_scriptWindow)


    def retranslateUi(self, run_py_scriptWindow):
        _translate = QtCore.QCoreApplication.translate
        run_py_scriptWindow.setWindowTitle(_translate("run_py_scriptWindow", "run_py_scriptWindow"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_g1), _translate("run_py_scriptWindow", "Run .py or .log"))

        self.menuRun.setTitle(_translate("run_py_scriptWindow", "Run"))
##        self.menuHelp.setTitle(_translate("run_py_scriptWindow", "Help"))
        self.actionRun.setText(_translate("run_py_scriptWindow", "run .py or .log script"))
        self.actionRun.setStatusTip(_translate("run_py_scriptWindow", "run .py or .log script"))
        self.action_1.setText(_translate("run_py_scriptWindow", "run .py or .log script"))
##        self.actionHelp.setText(_translate("run_py_scriptWindow", "Help"))

    #####################
    # functions
    #####################

    def build_group_g1(self,run_py_scriptWindow):
        self.widget_g1 = QtWidgets.QWidget(self.tab_g1)
        self.widget_g1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_g1.setObjectName("widget")
        self.flay_g0.removeRow(2)
        self.flay_g0.insertRow(2, self.widget_g1)

        self.groupBox_g1 = QtWidgets.QGroupBox(self.widget_g1)
        self.groupBox_g1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_g1.setObjectName("groupBox_g1")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_g1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")


        self.label_g7 = QtWidgets.QLabel(self.groupBox_g1)
        self.gridLayout.addWidget(self.label_g7, 0, 0, 1, 1)
        self.label_g7.setObjectName("label_g7")
        self.pushButton_g7 = QtWidgets.QPushButton(self.groupBox_g1)
        self.gridLayout.addWidget(self.pushButton_g7, 0, 1, 1, 1)
        self.pushButton_g7.setObjectName("pushButton_g7")
        self.pushButton_g7.clicked.connect(self.on_pushButton_clicked_g0)        
        self.lineEdit_g7 = QtWidgets.QLineEdit(self.groupBox_g1)
        self.gridLayout.addWidget(self.lineEdit_g7, 0, 2, 1, 1)
        self.lineEdit_g7.setObjectName("lineEdit_g7")
        
        self.label_g7.setText('*.py or *.log file')
        self.pushButton_g7.setText("...")
        self.pushButton_g7.setStatusTip('Choose a .py or matsdp.log file')
        self.lineEdit_g7.setText("None")
        self.lineEdit_g7.setStatusTip('Enter a .py or matsdp.log file path (relative path is also accepted)')

        self.groupBox_g1.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_g1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(70)
        layout = QtWidgets.QVBoxLayout(self.widget_g1)
        layout.addWidget(scroll)
        
        self.groupBox_g1.setTitle("file importer")

    def on_pushButton_clicked_g0(self):
        self.open_dialog_box_g0_file()

        
    def open_dialog_box_g0_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_g7.setText('{}'.format(file_path))




