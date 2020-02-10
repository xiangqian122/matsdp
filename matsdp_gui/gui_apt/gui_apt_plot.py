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
from matsdp.apt import apt_plot

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog
from functools import partial

defaults_dict = default_params.default_params()
logfile = defaults_dict['logfile']

periodic_table_dict = periodic_table.periodic_tab()
dos_mode = periodic_table_dict['dos_mode']
elmt_color = periodic_table_dict['elmt_color']

class Ui_apt_plotWindow(object):
    def setupUi(self, apt_plotWindow):
        apt_plotWindow.setObjectName("apt_plotWindow")
        apt_plotWindow.resize(1030, 600)
        self.centralwidget = QtWidgets.QWidget(apt_plotWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 990, 511))
        self.tabWidget.setObjectName("tabWidget")

        ###############################################
        # vasp_write
        ###############################################

        # write_poscar_with_force
        self.tab_e1 = QtWidgets.QScrollArea()
        self.tab_e1.setObjectName("tab")
        self.tabWidget.addTab(self.tab_e1, "")
        
        content_widget_e0 = QtWidgets.QWidget()
        self.tab_e1.setWidget(content_widget_e0)
        self.flay_e0 = QtWidgets.QFormLayout(content_widget_e0)
        self.tab_e1.setWidgetResizable(True)

        self.widget_e0 = QtWidgets.QWidget(self.tab_e1)
        self.widget_e0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_e0.setObjectName("widget")
        self.flay_e0.addRow(self.widget_e0)

        self.widget_e1 = QtWidgets.QWidget(self.tab_e1)
        self.widget_e1.setGeometry(QtCore.QRect(20, 90, 290, 450))
        self.widget_e1.setObjectName("widget")
        self.flay_e0.addRow(self.widget_e1)
        
        #parameters
        self.sysname_e0 = 'system'
        self.proxigram_csv_file_path = None
        self.visible_elmt_list_e0 = ['Ni']
        self.interpolation_on = False
        self.fig_width_e0 = 6
        self.fig_height_e0 = 5
        self.fig_dpi_e0 = 600
        self.fig_format_e0 = 'png'
        
        self.fig_format_list_e0 = ['png','eps','pdf','tif','tiff','jpg','jepg','svg','svgz','pgf','ps','raw','rgba']


        #Initialize the groups
        self.build_group_e0(apt_plotWindow)
        self.build_group_e1(apt_plotWindow)

        apt_plotWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(apt_plotWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
##        self.menuHelp = QtWidgets.QMenu(self.menubar)
##        self.menuHelp.setObjectName("menuHelp")
        apt_plotWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(apt_plotWindow)
        self.statusbar.setObjectName("statusbar")
        apt_plotWindow.setStatusBar(self.statusbar)
        self.actionRun = QtWidgets.QAction(apt_plotWindow)
        self.actionRun.setObjectName("actionRun")
        self.action_1 = QtWidgets.QAction(apt_plotWindow)
        self.action_1.setObjectName("action_1")
        self.action_2 = QtWidgets.QAction(apt_plotWindow)
        self.action_2.setObjectName("action_2")
##        self.actionHelp = QtWidgets.QAction(apt_plotWindow)
##        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_1 = QtWidgets.QAction(apt_plotWindow)
        self.actionRun_1.setObjectName("actionRun_1")
        self.menuRun.addAction(self.actionRun)
        self.menubar.addAction(self.menuRun.menuAction())
##        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(apt_plotWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(apt_plotWindow)


    def retranslateUi(self, apt_plotWindow):
        _translate = QtCore.QCoreApplication.translate
        apt_plotWindow.setWindowTitle(_translate("apt_plotWindow", "apt_plot"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_e1), _translate("apt_plotWindow", "plot_csv_proxigram"))

        self.menuRun.setTitle(_translate("apt_plotWindow", "Run"))
##        self.menuHelp.setTitle(_translate("apt_plotWindow", "Help"))
        self.actionRun.setText(_translate("apt_plotWindow", "apt_plot"))
        self.actionRun.setStatusTip(_translate("apt_plotWindow", "apt_plot"))
        self.action_1.setText(_translate("apt_plotWindow", "apt_plot"))
##        self.actionHelp.setText(_translate("apt_plotWindow", "Help"))

    #####################
    # functions
    #####################

    ## apt_plot
    
    def build_group_e0(self,apt_plotWindow):
        self.widget_e0 = QtWidgets.QWidget(self.tab_e1)
        self.widget_e0.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_e0.setObjectName("widget_e0")
        self.flay_e0.removeRow(1)
        self.flay_e0.insertRow(1, self.widget_e0)

        self.groupBox_e0 = QtWidgets.QGroupBox(self.widget_e0)
        self.groupBox_e0.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_e0.setObjectName("groupBox_e0")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_e0)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
  
        self.label_e1 = QtWidgets.QLabel(self.groupBox_e0)
        self.gridLayout.addWidget(self.label_e1, 0, 0, 1, 1)
        self.label_e1.setObjectName("label_e1")
        self.lineEdit_e1 = QtWidgets.QLineEdit(self.groupBox_e0)
        self.gridLayout.addWidget(self.lineEdit_e1, 0, 1, 1, 1)
        self.lineEdit_e1.setObjectName("lineEdit")
 

        self.label_e2 = QtWidgets.QLabel(self.groupBox_e0)
        self.gridLayout.addWidget(self.label_e2, 0, 2, 1, 1)
        self.label_e2.setObjectName("label_e2")
        self.lineEdit_e2 = QtWidgets.QLineEdit(self.groupBox_e0)
        self.gridLayout.addWidget(self.lineEdit_e2, 0, 3, 1, 1)
        self.lineEdit_e2.setObjectName("lineEdit")

        self.checkBox_e3 = QtWidgets.QCheckBox(self.groupBox_e0)
        self.gridLayout.addWidget(self.checkBox_e3, 0, 4, 1, 1)
        self.checkBox_e3.setObjectName("checkBox_e3")
##        self.checkBox_e3.toggle()        
        self.checkBox_e3.toggled.connect(self.set_value_e3)

        self.label_e4 = QtWidgets.QLabel(self.groupBox_e0)
        self.gridLayout.addWidget(self.label_e4, 1, 0, 1, 1)
        self.label_e4.setObjectName("label_e4")
        self.lineEdit_e4 = QtWidgets.QLineEdit(self.groupBox_e0)
        self.gridLayout.addWidget(self.lineEdit_e4, 1, 1, 1, 1)
        self.lineEdit_e4.setObjectName("lineEdit")

        self.label_e5 = QtWidgets.QLabel(self.groupBox_e0)
        self.gridLayout.addWidget(self.label_e5, 1, 2, 1, 1)
        self.label_e5.setObjectName("label_e5")
        self.lineEdit_e5 = QtWidgets.QLineEdit(self.groupBox_e0)
        self.gridLayout.addWidget(self.lineEdit_e5, 1, 3, 1, 1)
        self.lineEdit_e5.setObjectName("lineEdit")

        self.label_e6 = QtWidgets.QLabel(self.groupBox_e0)
        self.gridLayout.addWidget(self.label_e6, 1, 4, 1, 1)
        self.label_e6.setObjectName("label_e6")
        self.lineEdit_e6 = QtWidgets.QLineEdit(self.groupBox_e0)
        self.gridLayout.addWidget(self.lineEdit_e6, 1, 5, 1, 1)
        self.lineEdit_e6.setObjectName("lineEdit")
 
        self.label_e7 = QtWidgets.QLabel(self.groupBox_e0)
        self.gridLayout.addWidget(self.label_e7, 1, 6, 1, 1)
        self.label_e7.setObjectName("label_e7")
        self.comboBox_e7 = QtWidgets.QComboBox(self.groupBox_e0)
        self.gridLayout.addWidget(self.comboBox_e7, 1, 7, 1, 1)
        self.comboBox_e7.setObjectName("comboBox")
        self.comboBox_e7.addItems(self.fig_format_list_e0)
        self.comboBox_e7.activated[str].connect(self.fig_format_e0_handleActivated)

        
        self.label_e1.setText('sysname')
        self.label_e1.setStatusTip('sysname: User-defined system name')
        self.lineEdit_e1.setText('system')
        self.lineEdit_e1.setStatusTip('sysname: User-defined system name')
        self.label_e2.setText('visible element')
        self.label_e2.setStatusTip('Specify which elements are to be plotted')
        self.lineEdit_e2.setText('Ni')
        self.lineEdit_e2.setStatusTip('Specify which elements are to be plotted. e.g. Ni, Re, Al')
        self.checkBox_e3.setText("interpolation")
        self.checkBox_e3.setStatusTip('Whether to interpolate the concentration profile or not')
        self.label_e4.setText('fig_width')
        self.lineEdit_e4.setText('6')
        self.label_e5.setText('fig_height')
        self.lineEdit_e5.setText('5')
        self.label_e6.setText('fig_dpi')
        self.lineEdit_e6.setText('600')
        self.lineEdit_e6.setStatusTip('Figure DPI')
        self.label_e7.setText('fig_format')

        self.groupBox_e0.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_e0)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(240)
        layout = QtWidgets.QVBoxLayout(self.widget_e0)
        layout.addWidget(scroll)
        
        self.groupBox_e0.setTitle("general settings")


    def build_group_e1(self,apt_plotWindow):
        self.widget_e1 = QtWidgets.QWidget(self.tab_e1)
        self.widget_e1.setGeometry(QtCore.QRect(20, 70, 290, 450))
        self.widget_e1.setObjectName("widget")
        self.flay_e0.removeRow(2)
        self.flay_e0.insertRow(2, self.widget_e1)

        self.groupBox_e1 = QtWidgets.QGroupBox(self.widget_e1)
        self.groupBox_e1.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_e1.setObjectName("groupBox_e1")
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_e1)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.label_e8 = QtWidgets.QLabel(self.groupBox_e1)
        self.gridLayout.addWidget(self.label_e8, 0, 0, 1, 1)
        self.label_e8.setObjectName("label_e8")
        self.pushButton_e8 = QtWidgets.QPushButton(self.groupBox_e1)
        self.gridLayout.addWidget(self.pushButton_e8, 0, 1, 1, 1)
        self.pushButton_e8.setObjectName("pushButton_e8")
        self.pushButton_e8.clicked.connect(self.on_pushButton_clicked_e0)        
        self.lineEdit_e8 = QtWidgets.QLineEdit(self.groupBox_e1)
        self.gridLayout.addWidget(self.lineEdit_e8, 0, 2, 1, 1)
        self.lineEdit_e8.setObjectName("lineEdit_e8")
        
        self.label_e8.setText('proxigram_csv file')
        self.pushButton_e8.setText("...")
        self.pushButton_e8.setStatusTip('Choose the proxigram .csv file')
        self.lineEdit_e8.setText("None")
        self.lineEdit_e8.setStatusTip('Enter the proxigram .csv file path (relative path is also accepted)')

        self.groupBox_e1.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_e1)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(100)
        layout = QtWidgets.QVBoxLayout(self.widget_e1)
        layout.addWidget(scroll)
        
        self.groupBox_e1.setTitle("file importer")

    def on_pushButton_clicked_e0(self):
        self.open_dialog_box_e0_file()

        
    def open_dialog_box_e0_file(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_e8.setText('{}'.format(file_path))

    def fig_format_e0_handleActivated(self, text):
        self.fig_format_e0 = text

    def set_value_e3(self):
        self.interpolation_on = self.checkBox_e3.isChecked()

