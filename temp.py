# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Orbit_Ground_Track.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.setEnabled(True)
        MainWindow.resize(900, 560)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.OGT_label = QtWidgets.QLabel(self.centralwidget)
        self.OGT_label.setGeometry(QtCore.QRect(280, 0, 321, 31))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        font.setPointSize(18)
        font.setUnderline(True)
        self.OGT_label.setFont(font)
        self.OGT_label.setObjectName("OGT_label")
        self.leftGroup = QtWidgets.QGroupBox(self.centralwidget)
        self.leftGroup.setGeometry(QtCore.QRect(0, 400, 450, 130))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV25")
        font.setPointSize(12)
        font.setUnderline(True)
        self.leftGroup.setFont(font)
        self.leftGroup.setAutoFillBackground(True)
        self.leftGroup.setObjectName("leftGroup")
        self.line1_label = QtWidgets.QLabel(self.leftGroup)
        self.line1_label.setGeometry(QtCore.QRect(10, 20, 51, 16))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.line1_label.setFont(font)
        self.line1_label.setObjectName("line1_label")
        self.line2_label = QtWidgets.QLabel(self.leftGroup)
        self.line2_label.setGeometry(QtCore.QRect(10, 40, 51, 16))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.line2_label.setFont(font)
        self.line2_label.setObjectName("line2_label")
        self.line3_label = QtWidgets.QLabel(self.leftGroup)
        self.line3_label.setGeometry(QtCore.QRect(10, 60, 51, 16))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.line3_label.setFont(font)
        self.line3_label.setObjectName("line3_label")
        self.line1_entry = QtWidgets.QLineEdit(self.leftGroup)
        self.line1_entry.setGeometry(QtCore.QRect(70, 20, 371, 20))
        self.line1_entry.setText("")
        self.line1_entry.setObjectName("line1_entry")
        self.line2_entry = QtWidgets.QLineEdit(self.leftGroup)
        self.line2_entry.setGeometry(QtCore.QRect(70, 40, 371, 20))
        self.line2_entry.setObjectName("line2_entry")
        self.line3_entry = QtWidgets.QLineEdit(self.leftGroup)
        self.line3_entry.setGeometry(QtCore.QRect(70, 60, 371, 20))
        self.line3_entry.setObjectName("line3_entry")
        self.readButton = QtWidgets.QPushButton(self.leftGroup)
        self.readButton.setGeometry(QtCore.QRect(380, 90, 61, 23))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.readButton.setFont(font)
        self.readButton.setObjectName("readButton")
        self.rightGroup = QtWidgets.QGroupBox(self.centralwidget)
        self.rightGroup.setGeometry(QtCore.QRect(450, 400, 261, 130))
        self.rightGroup.setAutoFillBackground(True)
        self.rightGroup.setTitle("")
        self.rightGroup.setObjectName("rightGroup")
        self.SatName_label = QtWidgets.QLabel(self.rightGroup)
        self.SatName_label.setGeometry(QtCore.QRect(20, 10, 47, 14))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.SatName_label.setFont(font)
        self.SatName_label.setObjectName("SatName_label")
        self.CatNo_Label = QtWidgets.QLabel(self.rightGroup)
        self.CatNo_Label.setGeometry(QtCore.QRect(20, 30, 91, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.CatNo_Label.setFont(font)
        self.CatNo_Label.setObjectName("CatNo_Label")
        self.Class_label = QtWidgets.QLabel(self.rightGroup)
        self.Class_label.setGeometry(QtCore.QRect(20, 50, 91, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.Class_label.setFont(font)
        self.Class_label.setObjectName("Class_label")
        self.EpochDate_label = QtWidgets.QLabel(self.rightGroup)
        self.EpochDate_label.setGeometry(QtCore.QRect(20, 70, 91, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.EpochDate_label.setFont(font)
        self.EpochDate_label.setObjectName("EpochDate_label")
        self.EpochTime_label = QtWidgets.QLabel(self.rightGroup)
        self.EpochTime_label.setGeometry(QtCore.QRect(20, 90, 81, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.EpochTime_label.setFont(font)
        self.EpochTime_label.setObjectName("EpochTime_label")
        self.verticalRightGroup = QtWidgets.QGroupBox(self.centralwidget)
        self.verticalRightGroup.setGeometry(QtCore.QRect(710, 30, 181, 501))
        self.verticalRightGroup.setAutoFillBackground(True)
        self.verticalRightGroup.setTitle("")
        self.verticalRightGroup.setObjectName("verticalRightGroup")
        self.Prediction_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Prediction_label.setGeometry(QtCore.QRect(10, 310, 151, 21))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        font.setUnderline(True)
        self.Prediction_label.setFont(font)
        self.Prediction_label.setObjectName("Prediction_label")
        self.Prediction_entry_label = QtWidgets.QDateTimeEdit(self.verticalRightGroup)
        self.Prediction_entry_label.setGeometry(QtCore.QRect(10, 330, 171, 21))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.Prediction_entry_label.setFont(font)
        self.Prediction_entry_label.setObjectName("Prediction_entry_label")
        self.PredLat_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.PredLat_label.setGeometry(QtCore.QRect(10, 360, 61, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.PredLat_label.setFont(font)
        self.PredLat_label.setObjectName("PredLat_label")
        self.PredLon_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.PredLon_label.setGeometry(QtCore.QRect(10, 380, 61, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.PredLon_label.setFont(font)
        self.PredLon_label.setObjectName("PredLon_label")
        self.runButton = QtWidgets.QPushButton(self.verticalRightGroup)
        self.runButton.setGeometry(QtCore.QRect(10, 10, 75, 23))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.runButton.setFont(font)
        self.runButton.setObjectName("runButton")
        self.Lat_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Lat_label.setGeometry(QtCore.QRect(10, 60, 51, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Lat_label.setFont(font)
        self.Lat_label.setObjectName("Lat_label")
        self.Lon_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Lon_label.setGeometry(QtCore.QRect(10, 80, 61, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Lon_label.setFont(font)
        self.Lon_label.setObjectName("Lon_label")
        self.Alt_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Alt_label.setGeometry(QtCore.QRect(10, 100, 51, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Alt_label.setFont(font)
        self.Alt_label.setObjectName("Alt_label")
        self.Vel_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Vel_label.setGeometry(QtCore.QRect(10, 120, 51, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Vel_label.setFont(font)
        self.Vel_label.setObjectName("Vel_label")
        self.Period_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Period_label.setGeometry(QtCore.QRect(10, 140, 81, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Period_label.setFont(font)
        self.Period_label.setObjectName("Period_label")
        self.Ecc_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Ecc_label.setGeometry(QtCore.QRect(10, 180, 31, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Ecc_label.setFont(font)
        self.Ecc_label.setObjectName("Ecc_label")
        self.Inc_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Inc_label.setGeometry(QtCore.QRect(10, 160, 31, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Inc_label.setFont(font)
        self.Inc_label.setObjectName("Inc_label")
        self.Semi_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Semi_label.setGeometry(QtCore.QRect(10, 240, 91, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Semi_label.setFont(font)
        self.Semi_label.setObjectName("Semi_label")
        self.Num_Orbit_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.Num_Orbit_label.setGeometry(QtCore.QRect(10, 30, 101, 21))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.Num_Orbit_label.setFont(font)
        self.Num_Orbit_label.setObjectName("Num_Orbit_label")
        self.Num_Orbit_SpinBox = QtWidgets.QSpinBox(self.verticalRightGroup)
        self.Num_Orbit_SpinBox.setGeometry(QtCore.QRect(120, 30, 42, 22))
        self.Num_Orbit_SpinBox.setObjectName("Num_Orbit_SpinBox")
        self.ARGP_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.ARGP_label.setGeometry(QtCore.QRect(10, 200, 41, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.ARGP_label.setFont(font)
        self.ARGP_label.setObjectName("ARGP_label")
        self.RAAN_label = QtWidgets.QLabel(self.verticalRightGroup)
        self.RAAN_label.setGeometry(QtCore.QRect(10, 220, 41, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex_IV50")
        self.RAAN_label.setFont(font)
        self.RAAN_label.setObjectName("RAAN_label")
        self.label_23 = QtWidgets.QLabel(self.verticalRightGroup)
        self.label_23.setGeometry(QtCore.QRect(10, 460, 141, 31))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.label_23.setFont(font)
        self.label_23.setObjectName("label_23")
        self.label = QtWidgets.QLabel(self.verticalRightGroup)
        self.label.setGeometry(QtCore.QRect(10, 260, 71, 16))
        font = QtGui.QFont()
        font.setFamily("Simplex")
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.Map = QtWidgets.QGraphicsView(self.centralwidget)
        self.Map.setGeometry(QtCore.QRect(0, 30, 711, 371))
        self.Map.setObjectName("Map")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 900, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.OGT_label.setText(_translate("MainWindow", "ORBIT GROUND TRACK"))
        self.leftGroup.setTitle(_translate("MainWindow", "TLE"))
        self.line1_label.setText(_translate("MainWindow", "Line 1:"))
        self.line2_label.setText(_translate("MainWindow", "Line 2:"))
        self.line3_label.setText(_translate("MainWindow", "Line 3:"))
        self.readButton.setText(_translate("MainWindow", "READ"))
        self.SatName_label.setText(_translate("MainWindow", "Name:"))
        self.CatNo_Label.setText(_translate("MainWindow", "Catalog no."))
        self.Class_label.setText(_translate("MainWindow", "Classification:"))
        self.EpochDate_label.setText(_translate("MainWindow", "Epoch Date:"))
        self.EpochTime_label.setText(_translate("MainWindow", "Epoch Time:"))
        self.Prediction_label.setText(_translate("MainWindow", "Predicted Pass Search"))
        self.PredLat_label.setText(_translate("MainWindow", "Latitude:"))
        self.PredLon_label.setText(_translate("MainWindow", "Longitude:"))
        self.runButton.setText(_translate("MainWindow", "RUN"))
        self.Lat_label.setText(_translate("MainWindow", "Latitude:"))
        self.Lon_label.setText(_translate("MainWindow", "Longitude:"))
        self.Alt_label.setText(_translate("MainWindow", "Altitude:"))
        self.Vel_label.setText(_translate("MainWindow", "Velocity:"))
        self.Period_label.setText(_translate("MainWindow", "Orbit Period:"))
        self.Ecc_label.setText(_translate("MainWindow", "Ecc:"))
        self.Inc_label.setText(_translate("MainWindow", "Inc:"))
        self.Semi_label.setText(_translate("MainWindow", "Semi-major a:"))
        self.Num_Orbit_label.setText(_translate("MainWindow", "Number of Orbits:"))
        self.ARGP_label.setText(_translate("MainWindow", "ARGP:"))
        self.RAAN_label.setText(_translate("MainWindow", "RAAN:"))
        self.label_23.setText(_translate("MainWindow", "Author: Abinay Brown"))
        self.label.setText(_translate("MainWindow", "UTC time:"))