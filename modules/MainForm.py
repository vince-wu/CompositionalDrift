# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainForm.ui'
#
# Created by: PyQt5 UI code generator 5.12.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1350, 950)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.tab)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.saveButton = QtWidgets.QPushButton(self.tab)
        self.saveButton.setObjectName("saveButton")
        self.gridLayout.addWidget(self.saveButton, 0, 2, 1, 1)
        self.simulateButton = QtWidgets.QPushButton(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.simulateButton.sizePolicy().hasHeightForWidth())
        self.simulateButton.setSizePolicy(sizePolicy)
        self.simulateButton.setObjectName("simulateButton")
        self.gridLayout.addWidget(self.simulateButton, 0, 0, 1, 1)
        self.exportButton = QtWidgets.QPushButton(self.tab)
        self.exportButton.setObjectName("exportButton")
        self.gridLayout.addWidget(self.exportButton, 0, 1, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 0, 4, 1, 1)
        self.loadButton = QtWidgets.QPushButton(self.tab)
        self.loadButton.setObjectName("loadButton")
        self.gridLayout.addWidget(self.loadButton, 0, 3, 1, 1)
        self.verticalLayout_5.addLayout(self.gridLayout)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.displayView = QtWidgets.QTextBrowser(self.tab)
        self.displayView.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.displayView.setLineWrapMode(QtWidgets.QTextEdit.WidgetWidth)
        self.displayView.setLineWrapColumnOrWidth(100)
        self.displayView.setReadOnly(False)
        self.displayView.setObjectName("displayView")
        self.horizontalLayout_6.addWidget(self.displayView)
        self.scrollArea = QtWidgets.QScrollArea(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(8)
        sizePolicy.setHeightForWidth(self.scrollArea.sizePolicy().hasHeightForWidth())
        self.scrollArea.setSizePolicy(sizePolicy)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents_3 = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_3.setGeometry(QtCore.QRect(0, 0, 391, 606))
        self.scrollAreaWidgetContents_3.setObjectName("scrollAreaWidgetContents_3")
        self.formLayout_6 = QtWidgets.QFormLayout(self.scrollAreaWidgetContents_3)
        self.formLayout_6.setObjectName("formLayout_6")
        self.label_3 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_3.setObjectName("label_3")
        self.formLayout_6.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.numAvDp_DoubleSpinBox = QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents_3)
        self.numAvDp_DoubleSpinBox.setEnabled(False)
        self.numAvDp_DoubleSpinBox.setMaximum(10000000.0)
        self.numAvDp_DoubleSpinBox.setObjectName("numAvDp_DoubleSpinBox")
        self.formLayout_6.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.numAvDp_DoubleSpinBox)
        self.label_4 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_4.setObjectName("label_4")
        self.formLayout_6.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_4)
        self.weightAvDp_DoubleSpinBox = QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents_3)
        self.weightAvDp_DoubleSpinBox.setEnabled(False)
        self.weightAvDp_DoubleSpinBox.setMaximum(10000000.0)
        self.weightAvDp_DoubleSpinBox.setObjectName("weightAvDp_DoubleSpinBox")
        self.formLayout_6.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.weightAvDp_DoubleSpinBox)
        self.label_5 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_5.setObjectName("label_5")
        self.formLayout_6.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_5)
        self.dispersityDoubleSpinBox = QtWidgets.QDoubleSpinBox(self.scrollAreaWidgetContents_3)
        self.dispersityDoubleSpinBox.setEnabled(False)
        self.dispersityDoubleSpinBox.setMaximum(10000000.0)
        self.dispersityDoubleSpinBox.setObjectName("dispersityDoubleSpinBox")
        self.formLayout_6.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.dispersityDoubleSpinBox)
        self.line = QtWidgets.QFrame(self.scrollAreaWidgetContents_3)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.formLayout_6.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.line)
        self.label_6 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_6.setObjectName("label_6")
        self.formLayout_6.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.label_6)
        self.graphComboBox = QtWidgets.QComboBox(self.scrollAreaWidgetContents_3)
        self.graphComboBox.setObjectName("graphComboBox")
        self.graphComboBox.addItem("")
        self.graphComboBox.addItem("")
        self.graphComboBox.addItem("")
        self.graphComboBox.addItem("")
        self.formLayout_6.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.graphComboBox)
        self.label_9 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_9.setObjectName("label_9")
        self.formLayout_6.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.label_9)
        self.runLengthSpinBox = QtWidgets.QSpinBox(self.scrollAreaWidgetContents_3)
        self.runLengthSpinBox.setMinimum(1)
        self.runLengthSpinBox.setObjectName("runLengthSpinBox")
        self.formLayout_6.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.runLengthSpinBox)
        self.label_10 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_10.setObjectName("label_10")
        self.formLayout_6.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.label_10)
        self.rowsSpinBox = QtWidgets.QSpinBox(self.scrollAreaWidgetContents_3)
        self.rowsSpinBox.setMaximum(1000000)
        self.rowsSpinBox.setProperty("value", 10)
        self.rowsSpinBox.setObjectName("rowsSpinBox")
        self.formLayout_6.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.rowsSpinBox)
        self.label_7 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_7.setObjectName("label_7")
        self.formLayout_6.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.label_7)
        self.animateComboBox = QtWidgets.QComboBox(self.scrollAreaWidgetContents_3)
        self.animateComboBox.setObjectName("animateComboBox")
        self.animateComboBox.addItem("")
        self.animateComboBox.addItem("")
        self.formLayout_6.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.animateComboBox)
        self.label_8 = QtWidgets.QLabel(self.scrollAreaWidgetContents_3)
        self.label_8.setObjectName("label_8")
        self.formLayout_6.setWidget(8, QtWidgets.QFormLayout.LabelRole, self.label_8)
        self.legendComboBox = QtWidgets.QComboBox(self.scrollAreaWidgetContents_3)
        self.legendComboBox.setObjectName("legendComboBox")
        self.legendComboBox.addItem("")
        self.legendComboBox.addItem("")
        self.formLayout_6.setWidget(8, QtWidgets.QFormLayout.FieldRole, self.legendComboBox)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents_3)
        self.horizontalLayout_6.addWidget(self.scrollArea)
        self.horizontalLayout_6.setStretch(0, 70)
        self.horizontalLayout_6.setStretch(1, 30)
        self.verticalLayout_5.addLayout(self.horizontalLayout_6)
        self.visualizeWindow = QtWidgets.QGraphicsView(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(8)
        sizePolicy.setHeightForWidth(self.visualizeWindow.sizePolicy().hasHeightForWidth())
        self.visualizeWindow.setSizePolicy(sizePolicy)
        self.visualizeWindow.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.visualizeWindow.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.visualizeWindow.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.visualizeWindow.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContentsOnFirstShow)
        self.visualizeWindow.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.visualizeWindow.setDragMode(QtWidgets.QGraphicsView.RubberBandDrag)
        self.visualizeWindow.setObjectName("visualizeWindow")
        self.verticalLayout_5.addWidget(self.visualizeWindow)
        self.scrollArea_2 = QtWidgets.QScrollArea(self.tab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(2)
        sizePolicy.setHeightForWidth(self.scrollArea_2.sizePolicy().hasHeightForWidth())
        self.scrollArea_2.setSizePolicy(sizePolicy)
        self.scrollArea_2.setWidgetResizable(True)
        self.scrollArea_2.setObjectName("scrollArea_2")
        self.inputLayout = QtWidgets.QWidget()
        self.inputLayout.setGeometry(QtCore.QRect(0, 0, 1314, 113))
        self.inputLayout.setObjectName("inputLayout")
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout(self.inputLayout)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.formLayout_2 = QtWidgets.QFormLayout()
        self.formLayout_2.setObjectName("formLayout_2")
        self.label = QtWidgets.QLabel(self.inputLayout)
        self.label.setObjectName("label")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.systemComboBox = QtWidgets.QComboBox(self.inputLayout)
        self.systemComboBox.setMinimumSize(QtCore.QSize(120, 0))
        self.systemComboBox.setEditable(False)
        self.systemComboBox.setObjectName("systemComboBox")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.systemComboBox.addItem("")
        self.formLayout_2.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.systemComboBox)
        self.label_14 = QtWidgets.QLabel(self.inputLayout)
        self.label_14.setObjectName("label_14")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_14)
        self.modelComboBox = QtWidgets.QComboBox(self.inputLayout)
        self.modelComboBox.setMinimumSize(QtCore.QSize(100, 0))
        self.modelComboBox.setObjectName("modelComboBox")
        self.modelComboBox.addItem("")
        self.modelComboBox.addItem("")
        self.formLayout_2.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.modelComboBox)
        self.label_2 = QtWidgets.QLabel(self.inputLayout)
        self.label_2.setObjectName("label_2")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.holdCheckBox = QtWidgets.QCheckBox(self.inputLayout)
        self.holdCheckBox.setText("")
        self.holdCheckBox.setObjectName("holdCheckBox")
        self.formLayout_2.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.holdCheckBox)
        self.horizontalLayout.addLayout(self.formLayout_2)
        self.formLayout_4 = QtWidgets.QFormLayout()
        self.formLayout_4.setObjectName("formLayout_4")
        self.label_12 = QtWidgets.QLabel(self.inputLayout)
        self.label_12.setObjectName("label_12")
        self.formLayout_4.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_12)
        self.dpSpinBox = QtWidgets.QSpinBox(self.inputLayout)
        self.dpSpinBox.setMinimum(1)
        self.dpSpinBox.setMaximum(1000000)
        self.dpSpinBox.setProperty("value", 100)
        self.dpSpinBox.setObjectName("dpSpinBox")
        self.formLayout_4.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.dpSpinBox)
        self.label_13 = QtWidgets.QLabel(self.inputLayout)
        self.label_13.setObjectName("label_13")
        self.formLayout_4.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_13)
        self.label_15 = QtWidgets.QLabel(self.inputLayout)
        self.label_15.setObjectName("label_15")
        self.formLayout_4.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_15)
        self.poolSpinBox = QtWidgets.QSpinBox(self.inputLayout)
        self.poolSpinBox.setMaximum(1000000000)
        self.poolSpinBox.setStepType(QtWidgets.QAbstractSpinBox.AdaptiveDecimalStepType)
        self.poolSpinBox.setProperty("value", 200000)
        self.poolSpinBox.setObjectName("poolSpinBox")
        self.formLayout_4.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.poolSpinBox)
        self.conversionDoubleSpinBox = QtWidgets.QDoubleSpinBox(self.inputLayout)
        self.conversionDoubleSpinBox.setMinimum(0.0)
        self.conversionDoubleSpinBox.setMaximum(100.0)
        self.conversionDoubleSpinBox.setProperty("value", 100.0)
        self.conversionDoubleSpinBox.setObjectName("conversionDoubleSpinBox")
        self.formLayout_4.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.conversionDoubleSpinBox)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.formLayout_4.setItem(3, QtWidgets.QFormLayout.LabelRole, spacerItem1)
        self.horizontalLayout.addLayout(self.formLayout_4)
        self.monomerRatio_layout = QtWidgets.QFormLayout()
        self.monomerRatio_layout.setObjectName("monomerRatio_layout")
        self.horizontalLayout.addLayout(self.monomerRatio_layout)
        self.horizontalLayout_5.addLayout(self.horizontalLayout)
        self.scrollArea_2.setWidget(self.inputLayout)
        self.verticalLayout_5.addWidget(self.scrollArea_2)
        self.verticalLayout_5.setStretch(1, 10)
        self.verticalLayout_5.setStretch(2, 2)
        self.verticalLayout_4.addLayout(self.verticalLayout_5)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")
        self.verticalLayout.addWidget(self.tabWidget)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.label.setBuddy(self.systemComboBox)
        self.label_14.setBuddy(self.modelComboBox)
        self.label_12.setBuddy(self.dpSpinBox)
        self.label_13.setBuddy(self.conversionDoubleSpinBox)
        self.label_15.setBuddy(self.poolSpinBox)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        self.systemComboBox.setCurrentIndex(2)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Compositional Drift"))
        self.saveButton.setToolTip(_translate("MainWindow", "Save current input values"))
        self.saveButton.setText(_translate("MainWindow", "Save"))
        self.simulateButton.setText(_translate("MainWindow", "Simulate"))
        self.exportButton.setToolTip(_translate("MainWindow", "<html><head/><body><p>Export polymer data as csv. Right click figures to export.</p></body></html>"))
        self.exportButton.setText(_translate("MainWindow", "Export"))
        self.loadButton.setToolTip(_translate("MainWindow", "Load input values from save file"))
        self.loadButton.setText(_translate("MainWindow", "Load"))
        self.label_3.setText(_translate("MainWindow", "Number Average DP"))
        self.label_4.setText(_translate("MainWindow", "Weight Average DP"))
        self.label_5.setText(_translate("MainWindow", "Dispersity Index"))
        self.label_6.setText(_translate("MainWindow", "Graph Type"))
        self.graphComboBox.setItemText(0, _translate("MainWindow", "Monomer Occurrence"))
        self.graphComboBox.setItemText(1, _translate("MainWindow", "Monomer Usage"))
        self.graphComboBox.setItemText(2, _translate("MainWindow", "Run Length"))
        self.graphComboBox.setItemText(3, _translate("MainWindow", "DP Distribution"))
        self.label_9.setToolTip(_translate("MainWindow", "Monomer to plot for run length"))
        self.label_9.setText(_translate("MainWindow", "Run Length Monomer"))
        self.label_10.setToolTip(_translate("MainWindow", "Number of polymers to visualize"))
        self.label_10.setText(_translate("MainWindow", "Rows to Show"))
        self.label_7.setText(_translate("MainWindow", "Animation"))
        self.animateComboBox.setItemText(0, _translate("MainWindow", "On"))
        self.animateComboBox.setItemText(1, _translate("MainWindow", "Off"))
        self.label_8.setText(_translate("MainWindow", "Legend"))
        self.legendComboBox.setItemText(0, _translate("MainWindow", "On"))
        self.legendComboBox.setItemText(1, _translate("MainWindow", "Off"))
        self.label.setToolTip(_translate("MainWindow", "Number of unique monomers in the system"))
        self.label.setText(_translate("MainWindow", "Polymer System "))
        self.systemComboBox.setCurrentText(_translate("MainWindow", "3-Monomer"))
        self.systemComboBox.setItemText(0, _translate("MainWindow", "1-Monomer"))
        self.systemComboBox.setItemText(1, _translate("MainWindow", "2-Monomer"))
        self.systemComboBox.setItemText(2, _translate("MainWindow", "3-Monomer"))
        self.systemComboBox.setItemText(3, _translate("MainWindow", "4-Monomer"))
        self.systemComboBox.setItemText(4, _translate("MainWindow", "5-Monomer"))
        self.systemComboBox.setItemText(5, _translate("MainWindow", "6-Monomer"))
        self.systemComboBox.setItemText(6, _translate("MainWindow", "7-Monomer"))
        self.systemComboBox.setItemText(7, _translate("MainWindow", "8-Monomer"))
        self.systemComboBox.setItemText(8, _translate("MainWindow", "9-Monomer"))
        self.systemComboBox.setItemText(9, _translate("MainWindow", "10-Monomer"))
        self.label_14.setToolTip(_translate("MainWindow", "Model used by simulation"))
        self.label_14.setText(_translate("MainWindow", "Model"))
        self.modelComboBox.setItemText(0, _translate("MainWindow", "Mayo-Lewis"))
        self.modelComboBox.setItemText(1, _translate("MainWindow", "Penultimate"))
        self.label_2.setToolTip(_translate("MainWindow", "If checked, monomers will never be used up"))
        self.label_2.setText(_translate("MainWindow", "Hold Composition"))
        self.label_12.setToolTip(_translate("MainWindow", "Average polymer length"))
        self.label_12.setText(_translate("MainWindow", "Average DP"))
        self.label_13.setToolTip(_translate("MainWindow", "Reaction conversion"))
        self.label_13.setText(_translate("MainWindow", "Percent Conversion"))
        self.label_15.setToolTip(_translate("MainWindow", "Total number of monomers in reaction"))
        self.label_15.setText(_translate("MainWindow", "Monomer Pool Size"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Compositional Drift"))


