# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!

import matsdp
from matsdp.vasp import vasp_build

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QFileDialog

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(565, 413)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(70, 50, 321, 271))
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
        self.actionNew = QtWidgets.QAction(MainWindow)
        self.actionNew.setObjectName("actionNew")
        self.menuFile.addAction(self.actionNew)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "matsdp"))
        self.vasp_build_button.setStatusTip(_translate("MainWindow", "build model by substitution"))
        self.vasp_build_button.setText(_translate("MainWindow", "VASP: Build"))
        self.vasp_plot_button.setStatusTip(_translate("MainWindow", "vasp plot"))
        self.vasp_plot_button.setText(_translate("MainWindow", "VASP: Plot"))
        self.vasp_analyze_button.setStatusTip(_translate("MainWindow", "plot POSCAR"))
        self.vasp_analyze_button.setText(_translate("MainWindow", "VASP: Analyze"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuAbout.setTitle(_translate("MainWindow", "About"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.actionNew.setText(_translate("MainWindow", "New"))

class Ui_vasp_buildWindow(object):
    def setupUi(self, vasp_buildWindow):
        vasp_buildWindow.setObjectName("vasp_buildWindow")
        vasp_buildWindow.resize(800, 600)
        self.centralwidget = QtWidgets.QWidget(vasp_buildWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(40, 20, 511, 341))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.layoutWidget = QtWidgets.QWidget(self.tab)
        self.layoutWidget.setGeometry(QtCore.QRect(30, 80, 451, 131))
        self.layoutWidget.setObjectName("layoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_2 = QtWidgets.QLabel(self.layoutWidget)
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout.addWidget(self.label_2)
        self.pushButton = QtWidgets.QPushButton(self.layoutWidget)
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.horizontalLayout.addWidget(self.pushButton)
        self.lineEdit = QtWidgets.QLineEdit(self.layoutWidget)
        self.lineEdit.setObjectName("lineEdit")
        self.horizontalLayout.addWidget(self.lineEdit)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label = QtWidgets.QLabel(self.layoutWidget)
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        self.pushButton_2 = QtWidgets.QPushButton(self.layoutWidget)
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        self.horizontalLayout_2.addWidget(self.pushButton_2)
        self.lineEdit_2 = QtWidgets.QLineEdit(self.layoutWidget)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.horizontalLayout_2.addWidget(self.lineEdit_2)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.checkBox = QtWidgets.QCheckBox(self.tab_2)
        self.checkBox.setGeometry(QtCore.QRect(110, 170, 201, 21))
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.checkBox.setFont(font)
        self.checkBox.setObjectName("checkBox")
        self.layoutWidget1 = QtWidgets.QWidget(self.tab_2)
        self.layoutWidget1.setGeometry(QtCore.QRect(96, 40, 241, 22))
        self.layoutWidget1.setObjectName("layoutWidget1")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.layoutWidget1)
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_3 = QtWidgets.QLabel(self.layoutWidget1)
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_3.addWidget(self.label_3)
        self.pushButton_3 = QtWidgets.QPushButton(self.layoutWidget1)
        self.pushButton_3.setObjectName("pushButton_3")
        self.horizontalLayout_3.addWidget(self.pushButton_3)
        self.lineEdit_3 = QtWidgets.QLineEdit(self.layoutWidget1)
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.horizontalLayout_3.addWidget(self.lineEdit_3)
        self.layoutWidget2 = QtWidgets.QWidget(self.tab_2)
        self.layoutWidget2.setGeometry(QtCore.QRect(100, 70, 241, 22))
        self.layoutWidget2.setObjectName("layoutWidget2")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.layoutWidget2)
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.label_4 = QtWidgets.QLabel(self.layoutWidget2)
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_4.addWidget(self.label_4)
        self.lineEdit_4 = QtWidgets.QLineEdit(self.layoutWidget2)
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.horizontalLayout_4.addWidget(self.lineEdit_4)
        self.layoutWidget3 = QtWidgets.QWidget(self.tab_2)
        self.layoutWidget3.setGeometry(QtCore.QRect(96, 100, 241, 22))
        self.layoutWidget3.setObjectName("layoutWidget3")
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout(self.layoutWidget3)
        self.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_5 = QtWidgets.QLabel(self.layoutWidget3)
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_5.addWidget(self.label_5)
        self.lineEdit_5 = QtWidgets.QLineEdit(self.layoutWidget3)
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.horizontalLayout_5.addWidget(self.lineEdit_5)
        self.layoutWidget4 = QtWidgets.QWidget(self.tab_2)
        self.layoutWidget4.setGeometry(QtCore.QRect(101, 130, 221, 22))
        self.layoutWidget4.setObjectName("layoutWidget4")
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout(self.layoutWidget4)
        self.horizontalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_6 = QtWidgets.QLabel(self.layoutWidget4)
        font = QtGui.QFont()
##        font.setFamily("3ds")
##        font.setPointSize(16)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.horizontalLayout_6.addWidget(self.label_6)
        self.lineEdit_6 = QtWidgets.QLineEdit(self.layoutWidget4)
        self.lineEdit_6.setObjectName("lineEdit_6")
        self.horizontalLayout_6.addWidget(self.lineEdit_6)
        self.tabWidget.addTab(self.tab_2, "")
        vasp_buildWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(vasp_buildWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
        self.menuAbout = QtWidgets.QMenu(self.menubar)
        self.menuAbout.setObjectName("menuAbout")
        vasp_buildWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(vasp_buildWindow)
        self.statusbar.setObjectName("statusbar")
        vasp_buildWindow.setStatusBar(self.statusbar)
        self.actionRun = QtWidgets.QAction(vasp_buildWindow)
        self.actionRun.setObjectName("actionRun")
        self.actionsubstitution = QtWidgets.QAction(vasp_buildWindow)
        self.actionsubstitution.setObjectName("actionsubstitution")
        self.actionsphere_seclection = QtWidgets.QAction(vasp_buildWindow)
        self.actionsphere_seclection.setObjectName("actionsphere_seclection")
        self.actionHelp = QtWidgets.QAction(vasp_buildWindow)
        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_selection_sphere = QtWidgets.QAction(vasp_buildWindow)
        self.actionRun_selection_sphere.setObjectName("actionRun_selection_sphere")
        self.menuRun.addAction(self.actionRun)
        self.menuRun.addAction(self.actionRun_selection_sphere)
        self.menubar.addAction(self.menuRun.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())

        self.retranslateUi(vasp_buildWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(vasp_buildWindow)

    def retranslateUi(self, vasp_buildWindow):
        _translate = QtCore.QCoreApplication.translate
        vasp_buildWindow.setWindowTitle(_translate("vasp_buildWindow", "MainWindow"))
        self.label_2.setText(_translate("vasp_buildWindow", "POSCAR"))
        self.pushButton.setStatusTip(_translate("vasp_buildWindow", "Find the POSCAR file"))
        self.pushButton.setText(_translate("vasp_buildWindow", "..."))
        self.lineEdit.setStatusTip(_translate("vasp_buildWindow", "Enter the POSCAR file path"))
        self.label.setText(_translate("vasp_buildWindow", ".subst"))
        self.pushButton_2.setStatusTip(_translate("vasp_buildWindow", "Find the .subst file"))
        self.pushButton_2.setText(_translate("vasp_buildWindow", "..."))
        self.lineEdit_2.setStatusTip(_translate("vasp_buildWindow", "Enter the .subst file path"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("vasp_buildWindow", "substitution"))
        self.checkBox.setText(_translate("vasp_buildWindow", "include mirror atoms?"))
        self.label_3.setText(_translate("vasp_buildWindow", "POSCAR"))
        self.pushButton_3.setText(_translate("vasp_buildWindow", "..."))
        self.lineEdit_3.setStatusTip(_translate("vasp_buildWindow", "Enter the POSCAR file path"))
        self.label_4.setText(_translate("vasp_buildWindow", "origin atom name"))
        self.label_5.setText(_translate("vasp_buildWindow", "radius"))
        self.label_6.setText(_translate("vasp_buildWindow", "output file name"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("vasp_buildWindow", "selection_sphere"))
        self.menuRun.setTitle(_translate("vasp_buildWindow", "Run"))
        self.menuAbout.setTitle(_translate("vasp_buildWindow", "About"))
        self.actionRun.setText(_translate("vasp_buildWindow", "Run substitution"))
        self.actionRun.setStatusTip(_translate("vasp_buildWindow", "Run the module"))
        self.actionsubstitution.setText(_translate("vasp_buildWindow", "substitution"))
        self.actionsphere_seclection.setText(_translate("vasp_buildWindow", "sphere seclection"))
        self.actionHelp.setText(_translate("vasp_buildWindow", "Help"))
        self.actionRun_selection_sphere.setText(_translate("vasp_buildWindow", "Run selection_sphere"))
        

        self.pushButton.clicked.connect(self.on_pushButton_clicked)
        self.pushButton_2.clicked.connect(self.on_pushButton_2_clicked)
        self.pushButton_3.clicked.connect(self.on_pushButton_3_clicked)

    def set_text(self, text):
        return text
    
    def on_pushButton_clicked(self):
        self.open_dialog_box()

    def on_pushButton_2_clicked(self):
        self.open_dialog_box_2()

    def on_pushButton_3_clicked(self):
        self.open_dialog_box_3()
    
    def open_dialog_box(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit.setText('{}'.format(file_path))

    def open_dialog_box_2(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_2.setText('{}'.format(file_path))

    def open_dialog_box_3(self):
        filename = QFileDialog.getOpenFileName()
        file_path = filename[0]
        self.lineEdit_3.setText('{}'.format(file_path))


class Ui_vasp_plotWindow(object):
    def setupUi(self, vasp_plotWindow):
        vasp_plotWindow.setObjectName("vasp_plotWindow")
        vasp_plotWindow.resize(1140, 800)
        self.centralwidget = QtWidgets.QWidget(vasp_plotWindow)
        self.centralwidget.setObjectName("centralwidget")
        
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 1041, 711))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.tabWidget.addTab(self.tab, "")
##        self.tab_2 = QtWidgets.QWidget()
        self.tab_2 = QtWidgets.QScrollArea()
        content_widget = QtWidgets.QWidget()
        self.tab_2.setWidget(content_widget)
        self.flay = QtWidgets.QFormLayout(content_widget)
        self.tab_2.setWidgetResizable(True)

        self.tab_2.setObjectName("tab_2")
        self.groupBox = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox.setGeometry(QtCore.QRect(20, 10, 631, 61))
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
        
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget0")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout0")

        self.label = QtWidgets.QLabel(self.groupBox)
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label.setObjectName("label")
        self.spinBox = QtWidgets.QSpinBox(self.groupBox)
        self.gridLayout.addWidget(self.spinBox, 0, 1, 1, 1)
        self.spinBox.setObjectName("spinBox")
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.gridLayout.addWidget(self.label_2, 0, 2, 1, 1)
        self.label_2.setObjectName("label_2")
        self.spinBox_2 = QtWidgets.QSpinBox(self.groupBox)
        self.gridLayout.addWidget(self.spinBox_2, 0, 3, 1, 1)
        self.spinBox_2.setObjectName("spinBox_2")        
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.gridLayout.addWidget(self.label_3, 0, 4, 1, 1)
        self.label_3.setObjectName("label_3")
        self.spinBox_3 = QtWidgets.QSpinBox(self.groupBox)
        self.gridLayout.addWidget(self.spinBox_3, 0, 5, 1, 1)
        self.spinBox_3.setObjectName("spinBox_3")
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.gridLayout.addWidget(self.label_4, 0, 6, 1, 1)
        self.label_4.setObjectName("label_4")
        self.spinBox_4 = QtWidgets.QSpinBox(self.groupBox)
        self.gridLayout.addWidget(self.spinBox_4, 0, 7, 1, 1)
        self.spinBox_4.setObjectName("spinBox_4")

        self.groupBox.setLayout(self.gridLayout)
        
        self.groupBox_3 = QtWidgets.QGroupBox(self.widget)
        self.groupBox_3.setGeometry(QtCore.QRect(20, 80, 271, 151))
        self.groupBox_3.setObjectName("groupBox_3")
        self.groupBox_4 = QtWidgets.QGroupBox(self.widget_2)
        self.groupBox_4.setGeometry(QtCore.QRect(310, 80, 361, 141))
        self.groupBox_4.setObjectName("groupBox_4")
        
        num_doscar = 3
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_3)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        for i_doscar in range(0, num_doscar):
            self.label_5 = QtWidgets.QLabel(self.groupBox_3)
            self.gridLayout.addWidget(self.label_5, i_doscar, 0, 1, 1)
            self.label_5.setObjectName("label_5")
            
            self.pushButton = QtWidgets.QPushButton(self.groupBox_3)
            self.gridLayout.addWidget(self.pushButton, i_doscar, 1, 1, 1)
            self.pushButton.setObjectName("pushButton")
            
            self.lineEdit = QtWidgets.QLineEdit(self.groupBox_3)
            self.gridLayout.addWidget(self.lineEdit, i_doscar, 2, 1, 1)
            self.lineEdit.setObjectName("lineEdit")
            
            self.label_5.setText("DOSCAR 1")
            self.pushButton.setText("...")
            self.lineEdit.setText("None")
        self.groupBox_3.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_3)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget)
        layout.addWidget(scroll)

##        self.formLayout = QtWidgets.QFormLayout()
##        print('2')
##        contcar_label_list = []
##        contcar_spin_box_list = []
##        sysname_label_list = []
##        sysname_line_edit_list = []
##        print("a")
##        for i in range(num_atoms):
##            contcar_label_list.append(QtWidgets.QLabel('DOSCAR'))
##            print("b")
####            contcar_spin_box_list.append(QtWidgets.QSpinBox())
####            print("c")
##            sysname_label_list.append(QtWidgets.QLabel('sysname'))
####            sysname_line_edit_list.append(QtWidgets.QLineEdit())
##            self.formLayout.addRow(contcar_label_list[i], sysname_label_list[i])
##            print("d")



        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_4)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        num_atoms = 12
        for i_atom in range(0, num_atoms):
            self.label_11 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_11, i_atom, 0, 1, 1)
            self.label_11.setObjectName("label_11")
            self.spinBox_9 = QtWidgets.QSpinBox(self.groupBox_4)
            self.gridLayout.addWidget(self.spinBox_9, i_atom, 1, 1, 1)
            self.spinBox_9.setObjectName("spinBox_9")
            
            self.label_12 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_12, i_atom, 2, 1, 1)
            self.label_12.setObjectName("label_12")
            self.lineEdit_12 = QtWidgets.QLineEdit(self.groupBox_4)
            self.gridLayout.addWidget(self.lineEdit_12, i_atom, 3, 1, 1)
            self.lineEdit_12.setObjectName("lineEdit_12")
            
            self.label_13 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_13, i_atom, 4, 1, 1)
            self.label_13.setObjectName("label_13")
            self.lineEdit_13 = QtWidgets.QLineEdit(self.groupBox_4)
            self.gridLayout.addWidget(self.lineEdit_13, i_atom, 5, 1, 1)
            self.lineEdit_13.setObjectName("lineEdit_13")

            self.label_14 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_14, i_atom, 6, 1, 1)
            self.label_14.setObjectName("label_14")
            self.lineEdit_14 = QtWidgets.QLineEdit(self.groupBox_4)
            self.gridLayout.addWidget(self.lineEdit_14, i_atom, 7, 1, 1)
            self.lineEdit_14.setObjectName("lineEdit_14")

            self.label_15 = QtWidgets.QLabel(self.groupBox_4)
            self.gridLayout.addWidget(self.label_15, i_atom, 8, 1, 1)
            self.label_15.setObjectName("label_15")
            self.lineEdit_15 = QtWidgets.QLineEdit(self.groupBox_4)
            self.gridLayout.addWidget(self.lineEdit_15, i_atom, 9, 1, 1)
            self.lineEdit_15.setObjectName("lineEdit_15")

            
            self.label_11.setText("DOSCAR")
            self.label_12.setText("sysname")
            self.lineEdit_12.setText("None")
            self.label_13.setText("atomname")
            self.lineEdit_13.setText("Ni1")
            self.label_14.setText("color")
            self.lineEdit_14.setText("black")
            self.label_15.setText("subplot_arg")
            self.lineEdit_15.setText("111")
            
        self.groupBox_4.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_4)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget_2)
        layout.addWidget(scroll)
  
        self.groupBox_5 = QtWidgets.QGroupBox(self.widget_3)
        self.groupBox_5.setGeometry(QtCore.QRect(20, 600, 431, 61))
        self.groupBox_5.setObjectName("groupBox_5")

        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_5)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        num_subplots = 4
        for i_subplot in range(0, num_subplots):
            self.label_21 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_21, i_subplot, 0, 1, 1)
            self.label_21.setObjectName("label_21")
            self.lineEdit_21 = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_21, i_subplot, 1, 1, 1)
            self.lineEdit_21.setObjectName("lineEdit_21")

            self.label_22 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_22, i_subplot, 2, 1, 1)
            self.label_22.setObjectName("label_22")
            self.lineEdit_22 = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_22, i_subplot, 3, 1, 1)
            self.lineEdit_22.setObjectName("lineEdit_22")

            self.label_23 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_23, i_subplot, 4, 1, 1)
            self.label_23.setObjectName("label_23")
            self.lineEdit_23 = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_23, i_subplot, 5, 1, 1)
            self.lineEdit_23.setObjectName("lineEdit_23")

            self.label_24 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_24, i_subplot, 6, 1, 1)
            self.label_24.setObjectName("label_24")
            self.lineEdit_24 = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_24, i_subplot, 7, 1, 1)
            self.lineEdit_24.setObjectName("lineEdit_24")

            self.label_25 = QtWidgets.QLabel(self.groupBox_5)
            self.gridLayout.addWidget(self.label_25, i_subplot, 8, 1, 1)
            self.label_25.setObjectName("label_25")
            self.lineEdit_25 = QtWidgets.QLineEdit(self.groupBox_5)
            self.gridLayout.addWidget(self.lineEdit_25, i_subplot, 9, 1, 1)
            self.lineEdit_25.setObjectName("lineEdit_25")

            self.checkBox_26 = QtWidgets.QCheckBox(self.groupBox_5)
            self.gridLayout.addWidget(self.checkBox_26, i_subplot, 10, 1, 1)
            self.checkBox_26.setObjectName("checkBox_26")

            self.checkBox_27 = QtWidgets.QCheckBox(self.groupBox_5)
            self.gridLayout.addWidget(self.checkBox_27, i_subplot, 11, 1, 1)
            self.checkBox_27.setObjectName("checkBox_27")

            self.checkBox_28 = QtWidgets.QCheckBox(self.groupBox_5)
            self.gridLayout.addWidget(self.checkBox_28, i_subplot, 12, 1, 1)
            self.checkBox_28.setObjectName("checkBox_28")
            
            self.checkBox_29 = QtWidgets.QCheckBox(self.groupBox_5)
            self.gridLayout.addWidget(self.checkBox_29, i_subplot, 13, 1, 1)
            self.checkBox_29.setObjectName("checkBox_29")
            
            self.label_21.setText("subplot_arg")
            self.lineEdit_21.setText("111")
            self.label_22.setText("xlo")
            self.lineEdit_22.setText("None")
            self.label_23.setText("xhi")
            self.lineEdit_23.setText("None")
            self.label_24.setText("ylo")
            self.lineEdit_24.setText("None")
            self.label_25.setText("yhi")
            self.lineEdit_25.setText("None")
            self.checkBox_26.setText("subplot_xtick")
            self.checkBox_27.setText("subplot_ytick")
            self.checkBox_28.setText("subplot_xlabel")
            self.checkBox_29.setText("subplot_xlabel")
            
        self.groupBox_5.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_5)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget_3)
        layout.addWidget(scroll)

        
        self.groupBox_2 = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox_2.setGeometry(QtCore.QRect(710, 20, 261, 101))
        self.groupBox_2.setObjectName("groupBox_2")

        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_2)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        num_elmts = 4
        for i_elmt in range(0, num_elmts):
            self.label_31 = QtWidgets.QLabel(self.groupBox_2)
            self.gridLayout.addWidget(self.label_31, i_elmt, 0, 1, 1)
            self.label_31.setObjectName("label_31")
            
            self.lineEdit_31 = QtWidgets.QLineEdit(self.groupBox_2)
            self.gridLayout.addWidget(self.lineEdit_31, i_elmt, 1, 1, 1)
            self.lineEdit_31.setObjectName("lineEdit_31")

            self.checkBox_32 = QtWidgets.QCheckBox(self.groupBox_2)
            self.gridLayout.addWidget(self.checkBox_32, i_elmt, 2, 1, 1)
            self.checkBox_32.setObjectName("checkBox_32")

            self.checkBox_33 = QtWidgets.QCheckBox(self.groupBox_2)
            self.gridLayout.addWidget(self.checkBox_33, i_elmt, 3, 1, 1)
            self.checkBox_33.setObjectName("checkBox_33")

            self.checkBox_34 = QtWidgets.QCheckBox(self.groupBox_2)
            self.gridLayout.addWidget(self.checkBox_34, i_elmt, 4, 1, 1)
            self.checkBox_34.setObjectName("checkBox_34")
            
            self.checkBox_35 = QtWidgets.QCheckBox(self.groupBox_2)
            self.gridLayout.addWidget(self.checkBox_35, i_elmt, 5, 1, 1)
            self.checkBox_35.setObjectName("checkBox_35")

            self.label_31.setText("element")        
            self.lineEdit_31.setText("Ni")
            self.checkBox_32.setText("s")
            self.checkBox_33.setText("p")
            self.checkBox_34.setText("d")
            self.checkBox_35.setText("LDOS")
            
        self.groupBox_2.setLayout(self.gridLayout)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.groupBox_2)
        scroll.setWidgetResizable(True)
        scroll.setFixedHeight(400)
        layout = QtWidgets.QVBoxLayout(self.widget_4)
        layout.addWidget(scroll)
        
        self.groupBox_6 = QtWidgets.QGroupBox(self.tab_2)
        self.groupBox_6.setGeometry(QtCore.QRect(20, 610, 481, 51))
        self.groupBox_6.setObjectName("groupBox_6")
        self.flay.addRow(self.groupBox_6)

        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_6)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")

        self.checkBox_40 = QtWidgets.QCheckBox(self.groupBox_6)
        self.gridLayout.addWidget(self.checkBox_40, i_elmt, 0, 1, 1)
        self.checkBox_40.setObjectName("checkBox_40")

        self.label_41 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_41, 0, 0, 1, 1)
        self.label_41.setObjectName("label")
        self.lineEdit_41 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_41, i_elmt, 1, 1, 1)
        self.lineEdit_41.setObjectName("lineEdit_41")
        
        self.label_42 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_42, 0, 2, 1, 1)
        self.label_42.setObjectName("label_42")
        self.lineEdit_42 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_42, i_elmt, 2, 1, 1)
        self.lineEdit_42.setObjectName("lineEdit_42")
        
        self.label_43 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_43, 0, 4, 1, 1)
        self.label_43.setObjectName("label_43")
        self.lineEdit_43 = QtWidgets.QLineEdit(self.groupBox_6)
        self.gridLayout.addWidget(self.lineEdit_43, i_elmt, 3, 1, 1)
        self.lineEdit_43.setObjectName("lineEdit_43")
        
        self.label_44 = QtWidgets.QLabel(self.groupBox_6)
        self.gridLayout.addWidget(self.label_44, 0, 6, 1, 1)
        self.label_44.setObjectName("label_44")


        self.groupBox_6.setLayout(self.gridLayout)
        
        self.tabWidget.addTab(self.tab_2, "")

        vasp_plotWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(vasp_plotWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 18))
        self.menubar.setObjectName("menubar")
        self.menuRun = QtWidgets.QMenu(self.menubar)
        self.menuRun.setObjectName("menuRun")
        self.menuAbout = QtWidgets.QMenu(self.menubar)
        self.menuAbout.setObjectName("menuAbout")
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
        self.actionHelp = QtWidgets.QAction(vasp_plotWindow)
        self.actionHelp.setObjectName("actionHelp")
        self.actionRun_1 = QtWidgets.QAction(vasp_plotWindow)
        self.actionRun_1.setObjectName("actionRun_1")
        self.menuRun.addAction(self.actionRun)
        self.menuRun.addAction(self.actionRun_1)
        self.menubar.addAction(self.menuRun.menuAction())
        self.menubar.addAction(self.menuAbout.menuAction())


        self.retranslateUi(vasp_plotWindow)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(vasp_plotWindow)


    def retranslateUi(self, vasp_plotWindow):
        _translate = QtCore.QCoreApplication.translate
        vasp_plotWindow.setWindowTitle(_translate("vasp_plotWindow", "vasp_plotWindow"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("vasp_plotWindow", "plot_poscar"))
        self.groupBox.setTitle(_translate("vasp_plotWindow", "control"))
        self.label_3.setText(_translate("vasp_plotWindow", "num_doscar"))
        self.label.setText(_translate("vasp_plotWindow", "num_atoms"))
        self.label_4.setText(_translate("vasp_plotWindow", "num_elements"))
        self.label_2.setText(_translate("vasp_plotWindow", "num_subplots"))

        self.groupBox_3.setTitle(_translate("vasp_plotWindow", "DOSCAR"))
        self.groupBox_4.setTitle(_translate("vasp_plotWindow", "Atom"))

        self.label_12.setText(_translate("vasp_plotWindow", "sysname"))
        self.groupBox_5.setTitle(_translate("vasp_plotWindow", "Subplot"))
        self.groupBox_2.setTitle(_translate("vasp_plotWindow", "Element"))
        self.groupBox_6.setTitle(_translate("vasp_plotWindow", "General setting"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("vasp_plotWindow", "plot_doscar"))
        print("1")
        self.menuRun.setTitle(_translate("vasp_plotWindow", "Run"))
        self.menuAbout.setTitle(_translate("vasp_plotWindow", "About"))
        self.actionRun.setText(_translate("vasp_plotWindow", "Run plot_poscar"))
        self.actionRun.setStatusTip(_translate("vasp_plotWindow", "Plot POSCAR"))
        self.action_1.setText(_translate("vasp_plotWindow", "plot_poscar"))
        self.action_2.setText(_translate("vasp_plotWindow", "plot_doscar"))
        self.actionHelp.setText(_translate("vasp_plotWindow", "Help"))
        self.actionRun_1.setText(_translate("vasp_plotWindow", "Run plot_doscar"))
        self.actionRun_1.setStatusTip(_translate("vasp_plotWindow", "Plot DOSCAR"))



class Controller:

    def __init__(self):
        pass

    def Show_MainWindow(self):

        self.MainWindow = QtWidgets.QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.MainWindow)
        self.ui.vasp_build_button.clicked.connect(self.Show_vasp_buildWindow)
        self.ui.vasp_plot_button.clicked.connect(self.Show_vasp_plotWindow)

        self.MainWindow.show()

    def Show_vasp_buildWindow(self):        
        self.vasp_buildWindow = QtWidgets.QMainWindow()
        self.ui = Ui_vasp_buildWindow()
        self.ui.setupUi(self.vasp_buildWindow)
        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_substitution)
        self.ui.actionRun_selection_sphere.triggered.connect(self.on_run_module_button_clicked_selection_sphere)

        self.vasp_buildWindow.show()

    def on_run_module_button_clicked_substitution(self):
        poscar_file_path = self.ui.lineEdit.text()
        substitution_list_file = self.ui.lineEdit_2.text()
        vasp_build.substitution(substitution_list_file, poscar_file_path)

    def on_run_module_button_clicked_selection_sphere(self):
        poscar_file_path = self.ui.lineEdit_3.text()
        origin_atom_name = self.ui.lineEdit_4.text()
        radius = float(self.ui.lineEdit_5.text())
        output_file_name = self.ui.lineEdit_6.text()
        if self.ui.checkBox.isChecked():
            include_mirror_atoms = True
        else:
            include_mirror_atoms = False
        
        vasp_build.selection_sphere(poscar_file_path, origin_atom_name, radius, include_mirror_atoms, output_file_name)

    def Show_vasp_plotWindow(self):        
        self.vasp_plotWindow = QtWidgets.QMainWindow()
        self.ui = Ui_vasp_plotWindow()
        self.ui.setupUi(self.vasp_plotWindow)
        self.ui.actionRun.triggered.connect(self.on_run_module_button_clicked_plot_poscar)
        self.ui.actionRun_1.triggered.connect(self.on_run_module_button_clicked_plot_doscar)

        self.vasp_plotWindow.show()

    def on_run_module_button_clicked_plot_poscar(self):
        poscar_file_path = self.ui.lineEdit.text()
        substitution_list_file = self.ui.lineEdit_2.text()
        vasp_build.substitution(substitution_list_file, poscar_file_path)

    def on_run_module_button_clicked_plot_doscar(self):
        poscar_file_path = self.ui.lineEdit_3.text()
        origin_atom_name = self.ui.lineEdit_4.text()
        radius = float(self.ui.lineEdit_5.text())
        output_file_name = self.ui.lineEdit_6.text()
        if self.ui.checkBox.isChecked():
            include_mirror_atoms = True
        else:
            include_mirror_atoms = False

        plot_poscar_mode = 'single'
        if plot_poscar_mode == 'single':
            vasp_plot.plot_poscar(
                poscar_file_path, euler_angle_type, phi, theta, psi, elmt_color, draw_mirror_atom, box_on,
                axis_indicator, plot_cell_basis_vector_label, plot_atom_label, fig_format, fig_dpi,
                draw_colormap, colormap_column_indx, colormap_vmin, colormap_vmax, vmin_color, vmax_color, colorbar_alignment)
        elif plot_poscar_mode == 'multiple':
            vasp_plot.plot_poscar_for_workdir(
                workdir, euler_angle_type, phi, theta, psi, elmt_color, draw_mirror_atom, box_on,
                axis_indicator, plot_cell_basis_vector_label, plot_atom_label, poscar_or_contcar, fig_format, fig_dpi,
                draw_colormap, colormap_column_indx, colormap_vmin, colormap_vmax, vmin_color, vmax_color, colorbar_alignment)

 
        vasp_plot.plot_dos(atom_doscar_file_path_list, atom_sysname_list, atom_indx_list, atom_palette_list, atom_subplot_arg_list,
                          subplot_arg_list, subplot_xlo_list, subplot_xhi_list, subplot_ylo_list, subplot_yhi_list,
                          subplot_xtick_list, subplot_ytick_list, subplot_xlabel_list, subplot_ylabel_list, subplot_share_xy_list, mainplot_axis_label_list,
                          dos_mode, fermi_shift_zero, peak_analyzer, fig_format, fig_size, fig_dpi)

        
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
